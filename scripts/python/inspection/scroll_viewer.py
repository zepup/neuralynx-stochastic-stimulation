#!/usr/bin/env python
"""
Interactive viewer for inspecting Neuralynx `.ncs` recordings alongside `.nev` events.

Main purpose:
- open a small set of channels from one recording block
- pan and zoom through raw traces
- mark config loads and stimulation-delivery events from the matching `.nev` file

Smoothly pan/zoom multichannel traces over a whole recording, with vertical dashed
markers at config-file loads (orange) and stimulation-delivery starts (green), each
labelled with the active Condition_N. Built on pyqtgraph, which auto-downsamples and
clips to the visible range, so it stays responsive even over millions of samples --
unlike the notebook's Plotly scrollers, which reship a whole window over the Jupyter
comm channel on every interaction.

Usage:
    python scroll_viewer.py patient
    python scroll_viewer.py folder session_subfolder
    python scroll_viewer.py folder session_subfolder --full

Navigation once open:
    - drag left/right to pan; scroll wheel or right-drag to zoom
    - drag a box (left-drag after clicking the 'A' auto-range off) or use the mouse to zoom x/y
    - the toolbar 'A' button (bottom-left of each plot) auto-ranges (zoom to fit)
    - markers within the visible range are labelled; labels refresh as you pan/zoom
"""
import os
import re
import sys
import glob
import argparse
import datetime
import warnings

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, ".."))
from neuralynx_io.neuralynx_io import load_ncs, load_nev
from shared.project_paths import get_data_dir

import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtCore

COND_RE = re.compile(r'(Condition_\d+)\.cfg', re.IGNORECASE)

# ---- built-in patient example config ------------------------------------------------
# Dataset-specific fields to edit for a new session:
# Example values below are kept as one working session template.
# - suffix
# - channels
# - events_name
# - date
# - start_clock / end_clock
# - title
# - patient ID
# Current example uses a placeholder patient label.
# Note:
# - Neuralynx recordings are commonly capped at about 4 hours per file.
# - One stimulation session can span multiple recording suffixes such as 0029 and 0030.
# - Replace the patient label if needed, but the example suffixes, channels,
#   filenames, and time window are intentionally left as a working example.

PATIENT_DIR = get_data_dir()
PATIENT_CFG = dict(
    directory=PATIENT_DIR,
    suffix="0029",  # working example recording suffix
    channels=["RPHIPP1", "RPULV7", "RPULV8", "RVMPFC2"],  # working example channels
    events_name="Events_0029.nev",  # working example event file
    date="2026-03-15",  # working example session date
    start_clock="13:36:00",  # working example crop start
    end_clock="13:56:00",  # working example crop end
    title="PATIENT_ID / boxA / 0029",
)


# ---- event parsing (same logic as the notebook's parse_config_and_stim_events) ------

def event_string(rec):
    s = rec['EventString']
    if isinstance(s, bytes):
        s = s.decode('latin1', errors='replace')
    return s.rstrip('\x00')


def parse_config_and_stim_events(nev):
    """Return (config_loads, stim_begins), each a list of (timestamp_us, condition_label).

    config_loads counts a load only if its raw '-ProcessConfigurationFile ...' command was
    followed by 'Process Config Started' before the 'ASHB NetCom Command Echo', which drops
    orphaned loads abandoned mid-retry (they never lead to stimulation). stim_begins marks the
    actual artifact onset, labelled with the condition active at the time (None if delivery
    began before any config load appeared in this nev).
    """
    order = np.argsort(nev['records']['TimeStamp'])
    config_loads, stim_begins = [], []
    current_condition = None
    pending_started = False
    for idx in order:
        rec = nev['records'][idx]
        label = event_string(rec)
        ts = int(rec['TimeStamp'])
        if 'ProcessConfigurationFile' in label and 'Echo' not in label:
            pending_started = False
        elif 'Process Config Started' in label:
            pending_started = True
        elif 'ProcessConfigurationFile' in label and 'Echo' in label:
            if pending_started:
                m = COND_RE.search(label)
                current_condition = m.group(1) if m else label
                config_loads.append((ts, current_condition))
            pending_started = False
        elif 'Stimulation delivery begins' in label:
            stim_begins.append((ts, current_condition))
    return config_loads, stim_begins


# ---- data loading -------------------------------------------------------------------

def clock_to_unix_us(clock_str, date):
    dt = datetime.datetime.strptime(f"{date} {clock_str}", "%Y-%m-%d %H:%M:%S")
    dt = dt.replace(tzinfo=datetime.timezone.utc)
    return int(dt.timestamp() * 1e6)


def load_channels(paths, crop=None):
    """Load {name: ncs} for the given {name: path}, optionally cropping to a (start_us, end_us)
    window so we don't hold a multi-hour recording in memory."""
    out = {}
    for name, path in paths.items():
        ncs = load_ncs(path)
        t, d = ncs['time'], ncs['data']
        if crop is not None:
            mask = (t >= crop[0]) & (t <= crop[1])
            t, d = t[mask], d[mask]
        out[name] = dict(time=t, data=d)
    return out


def resolve_config(args):
    if args.mode == "patient":
        c = dict(PATIENT_CFG)
        paths = {ch: os.path.join(c["directory"], f"{ch}_{c['suffix']}.ncs") for ch in c["channels"]}
        events_path = os.path.join(c["directory"], c["events_name"])
        crop = None
        if not args.full and c.get("start_clock"):
            crop = (clock_to_unix_us(c["start_clock"], c["date"]),
                    clock_to_unix_us(c["end_clock"], c["date"]))
        title = c["title"] + ("" if args.full else f"  ({c['start_clock']}-{c['end_clock']})")
        return paths, events_path, crop, title

    folder = os.path.join(PATIENT_DIR, args.folder)
    if not os.path.isdir(folder):
        sys.exit(f"No such data folder: {folder}")
    paths = {}
    for f in sorted(glob.glob(os.path.join(folder, "*.ncs"))):
        stem = os.path.splitext(os.path.basename(f))[0]
        canon = re.sub(r'_\d+$', '', stem)
        paths[canon] = f
    events_path = os.path.join(folder, "Events.nev")
    return paths, events_path, None, f"inspection / {args.folder}"


# ---- viewer -------------------------------------------------------------------------

class ScrollViewer(pg.GraphicsLayoutWidget):
    """Stacked, x-linked channel plots with dynamic in-view event markers."""

    MARKER_COLORS = {"config": (230, 159, 0), "stim": (0, 158, 115)}  # orange / green (color-blind safe)

    def __init__(self, channel_data, config_loads, stim_begins, title, initial_window_s=30,
                 downsample=True):
        super().__init__(show=True, title=title)
        self.channel_data = channel_data
        self.t0 = min(ch["time"][0] for ch in channel_data.values())

        # events as (seconds_since_t0, condition, kind)
        self.events = []
        for ts, cond in config_loads:
            self.events.append(((ts - self.t0) / 1e6, cond, "config"))
        for ts, cond in stim_begins:
            self.events.append(((ts - self.t0) / 1e6, cond, "stim"))
        self.events.sort(key=lambda e: e[0])

        self._marker_items = []  # currently-shown InfiniteLines / labels, cleared each redraw
        self.plots = []

        self.setWindowTitle(f"scroll_viewer — {title}")
        names = list(channel_data.items())
        for row, (name, ch) in enumerate(names):
            p = self.addPlot(row=row, col=0)
            t_s = (ch["time"].astype(np.float64) - self.t0) / 1e6
            data = ch["data"]
            p.plot(t_s, data, pen=pg.mkPen((200, 200, 200), width=1))
            if downsample:
                # 'peak' keeps the min AND max sample within each drawn pixel column, so no spike
                # is ever hidden; zooming in progressively decimates less until it's raw 1:1.
                p.setDownsampling(auto=True, mode="peak")
            p.setClipToView(True)  # only build/draw the visible x-range (fidelity-preserving, just faster)
            # Fix the y-range so scrolling/zooming along x does NOT auto-rescale y. The user can
            # still change y themselves (wheel or drag over the y axis) and it then stays put.
            ylo, yhi = float(np.min(data)), float(np.max(data))
            pad = 0.05 * ((yhi - ylo) or 1.0)
            p.setYRange(ylo - pad, yhi + pad, padding=0)
            p.enableAutoRange(y=False)
            p.setLabel("left", name, units="uV")
            p.showGrid(x=True, y=True, alpha=0.15)
            if row > 0:
                p.setXLink(self.plots[0])
            if row < len(names) - 1:
                p.getAxis("bottom").setStyle(showValues=False)
            self.plots.append(p)
        self.plots[-1].setLabel("bottom", "Time since start", units="s")

        # redraw markers when the (shared) x-range changes; debounce so panning stays smooth
        self._timer = QtCore.QTimer(singleShot=True)
        self._timer.timeout.connect(self._refresh_markers)
        self.plots[0].sigXRangeChanged.connect(lambda *_: self._timer.start(60))

        total_s = max((ch["time"][-1] - self.t0) / 1e6 for ch in channel_data.values())
        self.plots[0].setXRange(0, min(initial_window_s, total_s), padding=0)
        self._refresh_markers()

    def _refresh_markers(self):
        for it, p in self._marker_items:
            p.removeItem(it)
        self._marker_items = []

        (xmin, xmax), _ = self.plots[0].viewRange()
        in_view = [e for e in self.events if xmin <= e[0] <= xmax]
        # guard against label clutter when zoomed way out
        show_labels = len(in_view) <= 40

        top = self.plots[0]
        for t_s, cond, kind in in_view:
            color = self.MARKER_COLORS[kind]
            label_txt = None
            if show_labels:
                cond_label = cond if cond is not None else "(pre-recording)"
                label_txt = f"{cond_label} ({'config load' if kind == 'config' else 'stim begins'})"
            for i, p in enumerate(self.plots):
                line = pg.InfiniteLine(
                    pos=t_s, angle=90, movable=False,
                    pen=pg.mkPen(color, width=1, style=QtCore.Qt.DashLine),
                    label=label_txt if (i == 0 and label_txt) else None,
                    labelOpts=dict(position=0.95, color=color, movable=False,
                                   fill=(0, 0, 0, 120), angle=90),
                )
                p.addItem(line, ignoreBounds=True)
                self._marker_items.append((line, p))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", nargs="?", choices=["patient", "folder"], default="patient",
                    help="load the built-in patient config or inspect a subfolder")
    ap.add_argument("folder", nargs="?", help="subfolder inside NCS_PROJECT_DATA_DIR when mode=folder")
    ap.add_argument("--full", action="store_true", help="load the whole recording (no time crop)")
    ap.add_argument("--window", type=float, default=30, help="initial visible window length in seconds")
    ap.add_argument("--raw", action="store_true",
                    help="disable display downsampling (draw every sample; slower when zoomed far out)")
    args = ap.parse_args()
    if args.mode == "folder" and not args.folder:
        ap.error("folder mode needs a subfolder name")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        paths, events_path, crop, title = resolve_config(args)
        missing = [p for p in paths.values() if not os.path.exists(p)]
        if missing:
            sys.exit("Missing .ncs files:\n  " + "\n  ".join(missing))
        print(f"[{title}] loading {len(paths)} channels: {', '.join(paths)}")
        channel_data = load_channels(paths, crop=crop)
        nev = load_nev(events_path)
        config_loads, stim_begins = parse_config_and_stim_events(nev)
    print(f"[{title}] {len(config_loads)} config loads, {len(stim_begins)} stim-delivery starts")

    app = pg.mkQApp("scroll_viewer")
    pg.setConfigOptions(antialias=False)  # faster with dense data
    viewer = ScrollViewer(channel_data, config_loads, stim_begins, title,
                          initial_window_s=args.window, downsample=not args.raw)
    viewer.resize(1300, 800)
    viewer.show()
    app.exec() if hasattr(app, "exec") else app.exec_()


if __name__ == "__main__":
    main()
