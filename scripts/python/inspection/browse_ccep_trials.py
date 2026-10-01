"""
Interactive trial browser for inspecting stimulation-aligned responses.

Main purpose:
- load `stim_metadata.pkl` together with anatomy labels and `.ncs` files
- browse single stimulation trials for one response region/contact
- step through stimulation events in chronological order
- compare individual trials with the mean response

Session-specific values near the top of this file are left as working examples.
Users should adapt the patient label if needed and can also change channels,
suffixes, filenames, and labels for their own dataset.
"""

import os
import datetime

import warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("MacOSX")   # change to Qt5Agg or TkAgg if needed
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox
import neuralynx_io.neuralynx_io as nlio

warnings.filterwarnings("ignore", category=UserWarning, module="neuralynx_io")

# ============================================================
# SETTINGS — edit these
# ============================================================
BASE_DIR         = "/path/to/PatientData/byPatient/Epilepsy/PATIENT_ID/2026-03-15/stochasticStim"
STIM_PKL         = os.path.join(BASE_DIR, "stim_metadata.pkl")
ANATOMY_XLSX     = os.path.join(BASE_DIR, "PATIENT_ID_2mm.xlsx")
ANAT_CONTACT_COL = "contact"
ANAT_LABEL_COL   = "Unit"
# WM_PKL         = os.path.join(BASE_DIR, "white_matter_avg_data.pkl")
TARGET_INTENSITY   = 10
EPOCH_PRE_MS       = 200
EPOCH_POST_MS      = 300
USE_WM             = False
RECORDING_SUFFIX   = "0029"   # working example: 0029 or 0030
QUIET_VIEW_DECIM   = 20       # decimation for full-recording quiet view

ARTIFACT_SEARCH_MS = [0, 50]



# ============================================================
# Dropdown widget
# ============================================================

class Dropdown:
    """Upward-opening dropdown for matplotlib figures."""

    ITEM_H   = 0.026
    MAX_SHOW = 18

    def __init__(self, fig, rect, label, options, on_change):
        self.fig       = fig
        self._label    = label
        self.on_change = on_change
        self._open     = False
        self._options  = []
        self._selected = ""

        x, y, w, h = rect
        self._rect = rect

        self.ax_btn = fig.add_axes(rect)
        self.btn = Button(self.ax_btn, label, color="#f0f0f0", hovercolor="#d8d8d8")
        self.btn.label.set_fontsize(8)
        self.btn.on_clicked(self._toggle)

        # List axes — positioned above the button, hidden until opened
        self.ax_list = fig.add_axes([x, y + h, w, self.ITEM_H])
        self.ax_list.set_visible(False)
        self.ax_list.set_zorder(20)
        self.ax_list.set_xticks([])
        self.ax_list.set_yticks([])
        for sp in self.ax_list.spines.values():
            sp.set_edgecolor("gray")

        fig.canvas.mpl_connect("button_press_event", self._on_fig_click)
        self.set_options(options or [])

    def set_options(self, options, keep_selection=False):
        prev = self._selected
        self._options = list(options)
        if keep_selection and prev in self._options:
            self._selected = prev
        else:
            self._selected = self._options[0] if self._options else ""
        self._update_btn_label()
        self._rebuild_list()
        self._close()

    def get(self):
        return self._selected

    def _update_btn_label(self):
        self.btn.label.set_text(f"{self._label}: {self._selected or '—'}")

    def _rebuild_list(self):
        self.ax_list.cla()
        self.ax_list.set_xticks([])
        self.ax_list.set_yticks([])
        for sp in self.ax_list.spines.values():
            sp.set_edgecolor("gray")

        n_vis = min(len(self._options), self.MAX_SHOW)
        x, y, w, h = self._rect
        lh = max(n_vis, 1) * self.ITEM_H
        self.ax_list.set_position([x, y + h, w, lh])

        for i, opt in enumerate(self._options[:self.MAX_SHOW]):
            frac_y = 1 - (i + 0.5) / max(n_vis, 1)
            color  = "#a8d4ff" if opt == self._selected else "#fafafa"
            self.ax_list.text(
                0.03, frac_y, opt,
                ha="left", va="center", fontsize=8,
                transform=self.ax_list.transAxes,
                bbox=dict(boxstyle="square,pad=0.1", facecolor=color,
                          edgecolor="none"),
                clip_on=True,
            )

    def _toggle(self, event):
        if self._open:
            self._close()
        else:
            self._open = True
            self.ax_list.set_visible(True)
            self.fig.canvas.draw_idle()

    def _close(self):
        self._open = False
        self.ax_list.set_visible(False)
        self.fig.canvas.draw_idle()

    def _on_fig_click(self, event):
        if not self._open:
            return
        if event.inaxes == self.ax_list:
            try:
                inv = self.ax_list.transAxes.inverted()
                xf, yf = inv.transform((event.x, event.y))
            except Exception:
                return
            if 0 <= xf <= 1 and 0 <= yf <= 1:
                n_vis = min(len(self._options), self.MAX_SHOW)
                idx   = int((1 - yf) * n_vis)
                idx   = max(0, min(n_vis - 1, idx))
                if idx < len(self._options):
                    self._selected = self._options[idx]
                    self._update_btn_label()
                    self._rebuild_list()
                    self._close()
                    if self.on_change:
                        self.on_change(self._selected)
        elif event.inaxes != self.ax_btn:
            self._close()


# ============================================================
# helpers
# ============================================================

def load_stim_metadata(pkl_path):
    df              = pd.read_pickle(pkl_path)
    ts              = df["stim_timestamp_us"].to_numpy(dtype=np.float64)
    electrode1      = df["electrode1"].fillna("").astype(str).to_numpy()
    electrode2      = df["electrode2"].fillna("").astype(str).to_numpy()
    intensity1      = df["intensity1"].to_numpy(dtype=np.float64)
    intensity2      = df["intensity2"].to_numpy(dtype=np.float64)
    run_idx         = df["run_idx"].to_numpy(dtype=int)
    pair_idx        = df["pair_idx"].to_numpy(dtype=int)   # 0-based position within config
    rec_suffix      = df["recording_suffix"].astype(str).str.strip().to_numpy()
    box             = df["box"].fillna("").astype(str).str.strip().to_numpy()
    return ts, electrode1, electrode2, intensity1, intensity2, run_idx, pair_idx, rec_suffix, box



def load_ncs_file(path, avg_wm=None):
    ncs   = nlio.load_ncs(path)
    sig   = np.asarray(ncs["data"], dtype=np.float64).ravel()
    Fs    = float(ncs["sampling_rate"])
    ts0   = float(np.asarray(ncs["timestamp"]).ravel()[0])
    tsPer = 1e6 / Fs

    n_nan = np.sum(np.isnan(sig))
    if n_nan > 0:
        print(f"  NaN samples in {os.path.basename(path)}: {n_nan} ({n_nan/len(sig)*100:.2f}%)")
    sig[np.isnan(sig)] = 0.0

    if avg_wm is not None and len(avg_wm) == len(sig):
        wm = avg_wm.copy(); wm[np.isnan(wm)] = 0.0
        sig -= wm
        sig[np.isnan(sig)] = 0.0

    return sig, Fs, ts0, tsPer


def extract_epochs(sig, Fs, ts0, tsPer, stim_times_us, pre_ms, post_ms):
    pre_n  = int(round(pre_ms  / 1000 * Fs))
    post_n = int(round(post_ms / 1000 * Fs))
    N      = len(sig)
    stim_samp = np.round((stim_times_us - ts0) / tsPer).astype(int)

    t_ms = np.arange(-pre_n, post_n + 1) / Fs * 1000

    trials, valid_idx = [], []
    for i, s0 in enumerate(stim_samp):
        a = s0 - pre_n; b = s0 + post_n + 1
        if a < 0 or b > N: continue
        trials.append(sig[a:b].copy())
        valid_idx.append(i)

    if not trials:
        return np.empty((0, pre_n + post_n + 1)), t_ms, []
    return np.array(trials), t_ms, valid_idx


def find_artifact(trial, t_ms, search_ms):
    mask = (t_ms >= search_ms[0]) & (t_ms <= search_ms[1])
    if not np.any(mask):
        return None, None
    idx     = np.argmax(np.abs(trial[mask]))
    art_t   = t_ms[mask][idx]
    art_amp = trial[mask][idx]
    return art_t, art_amp


# ============================================================
# main
# ============================================================

def main():
    # ── load anatomy ─────────────────────────────────────────
    anat = pd.read_excel(ANATOMY_XLSX)
    anat.columns = anat.columns.str.strip()
    if ANAT_CONTACT_COL not in anat.columns or ANAT_LABEL_COL not in anat.columns:
        raise KeyError(f"Anatomy columns: {anat.columns.tolist()}\n"
                       f"Expected '{ANAT_CONTACT_COL}' and '{ANAT_LABEL_COL}'.")
    anat[ANAT_CONTACT_COL] = anat[ANAT_CONTACT_COL].astype(str).str.strip()
    anat[ANAT_LABEL_COL]   = anat[ANAT_LABEL_COL].astype(str).str.strip()

    contact_to_region = dict(zip(anat[ANAT_CONTACT_COL], anat[ANAT_LABEL_COL]))

    # region → sorted list of contacts (for resp dropdown)
    resp_region_to_contacts = defaultdict(list)
    for _, row in anat.iterrows():
        resp_region_to_contacts[row[ANAT_LABEL_COL]].append(row[ANAT_CONTACT_COL])
    all_resp_regions = sorted(resp_region_to_contacts.keys())

    # ── load stim metadata ───────────────────────────────────
    ts, electrode1, electrode2, intensity1, intensity2, run_idx, pair_idx, rec_suffix, box_col = load_stim_metadata(STIM_PKL)
    config_size = int(pair_idx.max()) + 1   # 40 (pair_idx is 0-based)

    avg_wm = None

    # ── state ────────────────────────────────────────────────
    state = {
        "resp_contact" : resp_region_to_contacts[all_resp_regions[0]][0] if all_resp_regions else "",
        "trial_idx"    : 0,
        "show_mean"    : False,
        "trials"       : None,
        "stim_labels"      : [],   # list of (e1, e2) per trial, in chronological order
        "config_ids"       : [],   # config file index (1-based, chronological) per trial
        "within_config_pos": [],   # ordinal position within that config file
        "n_configs"        : 0,
        "t_ms"         : None,
        "mean_wave"    : None,
        "n_trials"     : 0,
        "tz_offset_s"  : 0,       # seconds to add to fromtimestamp(ts_us/1e6) → hospital wall-clock
    }

    def load_data():
        rc = state["resp_contact"]
        state["trial_idx"] = 0
        state["show_mean"] = False

        sort_idx   = np.argsort(ts)
        stim_times = ts[sort_idx]
        e1_sorted  = electrode1[sort_idx]
        e2_sorted  = electrode2[sort_idx]
        run_sorted = run_idx[sort_idx]
        pos_sorted = pair_idx[sort_idx] + 1   # convert 0-based to 1-based

        ncs_path = os.path.join(BASE_DIR, f"{rc}_{RECORDING_SUFFIX}.ncs")
        if not os.path.isfile(ncs_path):
            print(f"  NCS not found: {ncs_path}")
            _set_empty(); return

        try:
            sig, Fs, ts0_sig, tsPer = load_ncs_file(ncs_path, avg_wm if USE_WM else None)
            trials, t_ms, valid_idx = extract_epochs(sig, Fs, ts0_sig, tsPer,
                                                     stim_times, EPOCH_PRE_MS, EPOCH_POST_MS)
        except Exception as e:
            print(f"  Error loading {rc}: {e}")
            _set_empty(); return

        # Compute timezone offset: TimeCreated (hospital wall-clock) minus fromtimestamp(ts0)
        try:
            with open(ncs_path, "rb") as f:
                raw_hdr_bytes = f.read(16384).decode("utf-8", errors="ignore")
            for line in raw_hdr_bytes.splitlines():
                if line.startswith("-TimeCreated"):
                    tc_str = line.split(None, 1)[1].strip()
                    tc_dt  = datetime.datetime.strptime(tc_str, "%Y/%m/%d %H:%M:%S")
                    state["tz_offset_s"] = (tc_dt - datetime.datetime.fromtimestamp(ts0_sig / 1e6)).total_seconds()
                    break
        except Exception as e:
            print(f"  Warning: could not compute tz offset: {e}")
            state["tz_offset_s"] = 0

        if len(trials) == 0:
            print("  No valid trials.")
            _set_empty(); return

        n_cfgs = len(set(run_sorted[i] for i in valid_idx))

        # stim window boundaries in minutes (for quiet view)
        first_stim_s = int(round((stim_times[valid_idx[0]]  - ts0_sig) / tsPer))
        last_stim_s  = int(round((stim_times[valid_idx[-1]] - ts0_sig) / tsPer))

        state["ncs_path"]          = ncs_path
        state["ts0_sig"]           = ts0_sig          # µs, Neuralynx epoch = Unix epoch
        state["rec_dur_min"]       = len(sig) / Fs / 60
        state["stim_start_min"]    = first_stim_s / Fs / 60
        state["stim_end_min"]      = last_stim_s  / Fs / 60
        state["stim_start_us"]     = ts0_sig + first_stim_s * tsPer
        state["stim_end_us"]       = ts0_sig + last_stim_s  * tsPer
        state["trials"]            = trials
        state["stim_labels"]       = [(e1_sorted[i], e2_sorted[i]) for i in valid_idx]
        state["config_ids"]        = [int(run_sorted[i]) for i in valid_idx]
        state["within_config_pos"] = [int(pos_sorted[i]) for i in valid_idx]
        state["n_configs"]         = n_cfgs
        state["t_ms"]              = t_ms
        state["mean_wave"]         = np.nanmean(trials, axis=0)
        state["n_trials"]          = len(trials)
        print(f"\nResp: {rc} ({contact_to_region.get(rc,'?')}) | {len(trials)} trials across {n_cfgs} configs of {config_size} stims each")

    def _set_empty():
        dummy_t = np.linspace(-EPOCH_PRE_MS, EPOCH_POST_MS, 100)
        state["trials"]      = np.full((1, 100), np.nan)
        state["stim_labels"]       = [("?", "?")]
        state["config_ids"]        = [0]
        state["within_config_pos"] = [0]
        state["n_configs"]         = 0
        state["ncs_path"]          = None
        state["ts0_sig"]           = None
        state["rec_dur_min"]       = None
        state["stim_start_min"]    = None
        state["stim_end_min"]      = None
        state["stim_start_us"]     = None
        state["stim_end_us"]       = None
        state["tz_offset_s"]       = 0
        state["t_ms"]        = dummy_t
        state["mean_wave"]   = np.full(100, np.nan)
        state["n_trials"]    = 0

    # ── build figure ─────────────────────────────────────────
    fig = plt.figure(figsize=(14, 7))
    fig.suptitle("CCEP Single Trial Browser", fontsize=12, fontweight="bold")

    ax = fig.add_axes([0.08, 0.38, 0.88, 0.55])

    # Resp region dropdown
    dd_resp_region = Dropdown(
        fig, [0.08, 0.22, 0.38, 0.07], "Resp region",
        all_resp_regions, on_change=None
    )

    # Resp contact dropdown (depends on resp region)
    dd_resp_contact = Dropdown(
        fig, [0.52, 0.22, 0.40, 0.07], "Resp contact",
        resp_region_to_contacts[all_resp_regions[0]] if all_resp_regions else [],
        on_change=None
    )

    def on_resp_region(region):
        contacts = resp_region_to_contacts.get(region, [])
        dd_resp_contact.set_options(contacts)
        state["resp_contact"] = dd_resp_contact.get()
        load_data(); draw()

    def on_resp_contact(contact):
        state["resp_contact"] = contact
        load_data(); draw()

    dd_resp_region.on_change  = on_resp_region
    dd_resp_contact.on_change = on_resp_contact

    # Navigation buttons
    btn_prev  = Button(fig.add_axes([0.08,  0.05, 0.10, 0.06]), "◀  Prev")
    btn_next  = Button(fig.add_axes([0.19,  0.05, 0.10, 0.06]), "Next  ▶")
    btn_mean  = Button(fig.add_axes([0.33,  0.05, 0.13, 0.06]), "Show Mean")
    ax_lbl    = fig.add_axes([0.50,  0.05, 0.08, 0.06]); ax_lbl.axis("off")
    ax_lbl.text(0.5, 0.5, "Go to trial:", ha="center", va="center", fontsize=9)
    txt_box   = TextBox(fig.add_axes([0.585, 0.05, 0.07, 0.06]), "", initial="1")
    btn_go    = Button(fig.add_axes([0.66,  0.05, 0.06, 0.06]), "Go")
    btn_quiet    = Button(fig.add_axes([0.74,  0.05, 0.11, 0.06]), "Quiet _29")
    btn_quiet_30 = Button(fig.add_axes([0.86,  0.05, 0.12, 0.06]), "Quiet _30")

    # ── draw ─────────────────────────────────────────────────
    def draw():
        ax.cla()
        trials    = state["trials"]
        t_ms      = state["t_ms"]
        mean_wave = state["mean_wave"]
        n_trials  = state["n_trials"]
        tidx      = state["trial_idx"]
        rc        = state["resp_contact"]
        rc_label  = contact_to_region.get(rc, "?")

        finite = trials[np.isfinite(trials)]
        if len(finite) > 0:
            y_min = np.percentile(finite, 1)
            y_max = np.percentile(finite, 99)
            margin = (y_max - y_min) * 0.1 or 5
            y_min -= margin; y_max += margin
        else:
            y_min, y_max = -50, 50

        ax.axvspan(-EPOCH_PRE_MS, 0, color="#e8f0fb", alpha=0.5, zorder=0)
        ax.axvspan(ARTIFACT_SEARCH_MS[0], ARTIFACT_SEARCH_MS[1],
                   color="#fff0cc", alpha=0.4, zorder=0, label="Artifact search")
        ax.axvline(0, color="black", linestyle="--", linewidth=1.5, label="t=0")
        ax.axhline(0, color="grey",  linewidth=0.5)

        if state["show_mean"]:
            ax.plot(t_ms, mean_wave, color="royalblue",
                    linewidth=2, label=f"Mean ({n_trials} trials)")
            art_t, art_amp = find_artifact(mean_wave, t_ms, ARTIFACT_SEARCH_MS)
            if art_t is not None:
                ax.plot(art_t, art_amp, "ro", markersize=8, zorder=5,
                        label=f"Artifact peak: {art_t:.0f} ms ({art_amp:.1f} µV)")
            title = f"Mean ({n_trials} trials)  |  all stims  →  resp: {rc} ({rc_label})  |  {TARGET_INTENSITY} mA"
        else:
            ax.plot(t_ms, mean_wave, color="royalblue", linewidth=1,
                    linestyle="--", alpha=0.4, label="Mean")
            trial = trials[tidx]
            e1, e2 = state["stim_labels"][tidx]
            e1_label = contact_to_region.get(e1, "?")
            e2_label = contact_to_region.get(e2, "?")
            cfg = state["config_ids"][tidx]
            pos = state["within_config_pos"][tidx]
            ax.plot(t_ms, trial, color="#222222", linewidth=1.5,
                    label=f"Trial {tidx+1}")
            art_t, art_amp = find_artifact(trial, t_ms, ARTIFACT_SEARCH_MS)
            if art_t is not None:
                ax.plot(art_t, art_amp, "ro", markersize=9, zorder=5)
                ax.annotate(f"{art_t:.0f} ms\n{art_amp:.1f} µV",
                            xy=(art_t, art_amp),
                            xytext=(art_t + 15, art_amp + (y_max - y_min) * 0.08),
                            fontsize=8, color="red",
                            arrowprops=dict(arrowstyle="->", color="red", lw=1))
            title = (f"Trial {tidx+1} / {n_trials}  |  Config {cfg} / {state['n_configs']}"
                     f", stim {pos} / {config_size}"
                     f"  |  stim: {e1} ({e1_label}) — {e2} ({e2_label})"
                     f"  →  resp: {rc} ({rc_label})")

        ax.set_title(title, fontsize=9, pad=4)
        ax.set_xlabel("Time relative to stim onset (ms)", fontsize=10)
        ax.set_ylabel("Amplitude (µV, raw)", fontsize=10)
        ax.set_xlim(-EPOCH_PRE_MS, EPOCH_POST_MS)
        ax.set_ylim(y_min, y_max)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.25)

        txt_box.set_val(str(tidx + 1))
        fig.canvas.draw_idle()

    # ── callbacks ────────────────────────────────────────────
    def on_prev(event):
        if state["n_trials"] == 0: return
        state["trial_idx"] = max(0, state["trial_idx"] - 1)
        state["show_mean"] = False
        draw()

    def on_next(event):
        if state["n_trials"] == 0: return
        state["trial_idx"] = min(state["n_trials"] - 1, state["trial_idx"] + 1)
        state["show_mean"] = False
        draw()

    def on_mean(event):
        state["show_mean"] = not state["show_mean"]
        btn_mean.label.set_text("Show Trial" if state["show_mean"] else "Show Mean")
        draw()

    def on_go(event):
        try:
            n = int(txt_box.text) - 1
            if 0 <= n < state["n_trials"]:
                state["trial_idx"] = n
                state["show_mean"] = False
                draw()
        except ValueError:
            pass

    def on_key(event):
        if   event.key == "left":  on_prev(None)
        elif event.key == "right": on_next(None)
        elif event.key == "m":     on_mean(None)

    def on_quiet(_event):
        if not state["ncs_path"]:
            print("No NCS file loaded yet.")
            return
        rc       = state["resp_contact"]
        rc_label = contact_to_region.get(rc, "?")
        print(f"Loading full recording for quiet view: {os.path.basename(state['ncs_path'])}")

        ncs     = nlio.load_ncs(state["ncs_path"])
        sig_full = np.asarray(ncs["data"], dtype=np.float64).ravel()
        Fs_full  = float(ncs["sampling_rate"])
        sig_full[np.isnan(sig_full)] = 0.0

        idx   = np.arange(0, len(sig_full), QUIET_VIEW_DECIM)
        t_min = idx / Fs_full / 60
        sig_d = sig_full[idx]
        rail  = 999.0
        sat_pct = np.mean(np.abs(sig_d) >= rail) * 100

        # Calibrated wall-clock time: ts0 + tz_offset_s computed from TimeCreated at load time
        ts0_plotted = float(np.asarray(ncs["timestamp"]).ravel()[0])
        rec_end_us  = ts0_plotted + len(sig_full) / Fs_full * 1e6

        def ts_to_wall(us):
            return (datetime.datetime.fromtimestamp(us / 1e6)
                    + datetime.timedelta(seconds=state["tz_offset_s"]))

        markers = [
            (0,                       ts_to_wall(ts0_plotted).strftime("%H:%M:%S"),             "navy",       "rec start"),
            (state["stim_start_min"], ts_to_wall(state["stim_start_us"]).strftime("%H:%M:%S"),  "darkorange", "stim start"),
            (state["stim_end_min"],   ts_to_wall(state["stim_end_us"]).strftime("%H:%M:%S"),    "darkorange", "stim end"),
            (state["rec_dur_min"],    ts_to_wall(rec_end_us).strftime("%H:%M:%S"),              "navy",       "rec end"),
        ]

        fig2, ax2 = plt.subplots(figsize=(14, 4))
        ax2.plot(t_min, sig_d, color="#222222", linewidth=0.3)
        ax2.axhline( rail, color="red", linestyle="--", linewidth=0.8, label=f"±{rail:.0f} µV rail")
        ax2.axhline(-rail, color="red", linestyle="--", linewidth=0.8)
        if state["stim_start_min"] is not None:
            ax2.axvspan(state["stim_start_min"], state["stim_end_min"],
                        color="orange", alpha=0.15, label="stim period")
        for x_min, label, color, _ in markers:
            ax2.axvline(x_min, color=color, linewidth=1.2, linestyle="--")
            ax2.text(x_min, -0.01, label,
                     color=color, fontsize=11, fontweight="bold",
                     va="top", ha="center", rotation=90, clip_on=False,
                     transform=ax2.get_xaxis_transform())
        ax2.set_xlabel("Time (min)")
        ax2.set_ylabel("Amplitude (µV)")
        ax2.set_title(f"{rc} ({rc_label})  |  {os.path.basename(state['ncs_path'])}"
                      f"  |  {sat_pct:.1f}% at rail  —  zoom to inspect")
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.25)
        fig2.tight_layout()
        fig2.show()

    def on_quiet_30(_event):
        rc       = state["resp_contact"]
        rc_label = contact_to_region.get(rc, "?")
        ncs_path_30 = os.path.join(BASE_DIR, f"{rc}_0030.ncs")
        if not os.path.isfile(ncs_path_30):
            print(f"  NCS not found: {ncs_path_30}")
            return
        print(f"Loading full _0030 recording for quiet view: {os.path.basename(ncs_path_30)}")

        ncs      = nlio.load_ncs(ncs_path_30)
        sig_full = np.asarray(ncs["data"], dtype=np.float64).ravel()
        Fs_full  = float(ncs["sampling_rate"])
        sig_full[np.isnan(sig_full)] = 0.0
        ts0_30 = float(np.asarray(ncs["timestamp"]).ravel()[0])

        # Compute tz_offset from _0030 TimeCreated
        tz_offset_30 = state["tz_offset_s"]   # fallback: same system, same offset
        try:
            with open(ncs_path_30, "rb") as f:
                raw_hdr_30 = f.read(16384).decode("utf-8", errors="ignore")
            for line in raw_hdr_30.splitlines():
                if line.startswith("-TimeCreated"):
                    tc_str = line.split(None, 1)[1].strip()
                    tc_dt  = datetime.datetime.strptime(tc_str, "%Y/%m/%d %H:%M:%S")
                    tz_offset_30 = (tc_dt - datetime.datetime.fromtimestamp(ts0_30 / 1e6)).total_seconds()
                    break
        except Exception as e:
            print(f"  Warning: could not compute _0030 tz offset: {e}")

        def ts_to_wall(us):
            return (datetime.datetime.fromtimestamp(us / 1e6)
                    + datetime.timedelta(seconds=tz_offset_30))

        def us_to_min(us):
            return (us - ts0_30) / 1e6 / 60

        stim_to_wall = ts_to_wall

        rec_dur_min = len(sig_full) / Fs_full / 60
        rec_end_us  = ts0_30 + len(sig_full) / Fs_full * 1e6

        # Split stims for _0030 by box; detect boxB_a vs boxB_b via largest gap
        # rec_suffix values may be "0030" or "30" — match both
        mask_30 = np.isin(rec_suffix, ["0030", "30"])
        ts_30   = ts[mask_30]
        box_30  = box_col[mask_30]
        print(f"  _0030 stim events : {len(ts_30)}")
        print(f"  Unique box values : {np.unique(box_30)}")

        groups = []
        ts_boxa = np.sort(ts_30[box_30 == "boxA"])
        if len(ts_boxa) > 0:
            groups.append(("boxA", ts_boxa.min(), ts_boxa.max(), "steelblue"))

        # Handle "boxB_a"/"boxB_b" directly in metadata, or split "boxB" by largest gap
        if np.any(box_30 == "boxB_a") or np.any(box_30 == "boxB_b"):
            ts_ba = np.sort(ts_30[box_30 == "boxB_a"])
            ts_bb = np.sort(ts_30[box_30 == "boxB_b"])
            if len(ts_ba) > 0: groups.append(("boxB_a", ts_ba.min(), ts_ba.max(), "darkorange"))
            if len(ts_bb) > 0: groups.append(("boxB_b", ts_bb.min(), ts_bb.max(), "purple"))
        else:
            ts_boxb_all = np.sort(ts_30[box_30 == "boxB"])
            if len(ts_boxb_all) > 1:
                split_i   = int(np.argmax(np.diff(ts_boxb_all)))
                ts_boxb_a = ts_boxb_all[:split_i + 1]
                ts_boxb_b = ts_boxb_all[split_i + 1:]
                if len(ts_boxb_a) > 0: groups.append(("boxB_a", ts_boxb_a.min(), ts_boxb_a.max(), "darkorange"))
                if len(ts_boxb_b) > 0: groups.append(("boxB_b", ts_boxb_b.min(), ts_boxb_b.max(), "purple"))
            elif len(ts_boxb_all) > 0:
                groups.append(("boxB", ts_boxb_all.min(), ts_boxb_all.max(), "darkorange"))

        for lbl, us_s, us_e, _ in groups:
            print(f"  {lbl}: {stim_to_wall(us_s).strftime('%H:%M:%S')} – "
                  f"{stim_to_wall(us_e).strftime('%H:%M:%S')}  "
                  f"({(us_e - us_s)/1e6/60:.1f} min)")

        idx   = np.arange(0, len(sig_full), QUIET_VIEW_DECIM)
        t_min = idx / Fs_full / 60
        sig_d = sig_full[idx]
        rail  = 999.0
        sat_pct = np.mean(np.abs(sig_d) >= rail) * 100

        fig2, ax2 = plt.subplots(figsize=(14, 4))
        ax2.plot(t_min, sig_d, color="#222222", linewidth=0.3)
        ax2.axhline( rail, color="red", linestyle="--", linewidth=0.8, label=f"±{rail:.0f} µV rail")
        ax2.axhline(-rail, color="red", linestyle="--", linewidth=0.8)

        for lbl, us_start, us_end, color in groups:
            ax2.axvspan(us_to_min(us_start), us_to_min(us_end),
                        color=color, alpha=0.15, label=f"stim {lbl}")
            for us in [us_start, us_end]:
                ax2.axvline(us_to_min(us), color=color, linewidth=1.2, linestyle="--")
                ax2.text(us_to_min(us), -0.01,
                         stim_to_wall(us).strftime("%H:%M:%S"),
                         color=color, fontsize=11, fontweight="bold",
                         va="top", ha="center", rotation=90, clip_on=False,
                         transform=ax2.get_xaxis_transform())

        for x_min, wall_str in [(0, ts_to_wall(ts0_30).strftime("%H:%M:%S")),
                                 (rec_dur_min, ts_to_wall(rec_end_us).strftime("%H:%M:%S"))]:
            ax2.axvline(x_min, color="navy", linewidth=1.2, linestyle="--")
            ax2.text(x_min, -0.01, wall_str,
                     color="navy", fontsize=11, fontweight="bold",
                     va="top", ha="center", rotation=90, clip_on=False,
                     transform=ax2.get_xaxis_transform())

        ax2.set_xlabel("Time (min)", labelpad=60)
        ax2.set_ylabel("Amplitude (µV)")
        ax2.set_title(f"{rc} ({rc_label})  |  {os.path.basename(ncs_path_30)}"
                      f"  |  {sat_pct:.1f}% at rail  —  boxA / boxB_a / boxB_b")
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.25)
        fig2.tight_layout()
        fig2.show()

    btn_prev.on_clicked(on_prev)
    btn_next.on_clicked(on_next)
    btn_mean.on_clicked(on_mean)
    btn_go.on_clicked(on_go)
    btn_quiet.on_clicked(on_quiet)
    btn_quiet_30.on_clicked(on_quiet_30)
    fig.canvas.mpl_connect("key_press_event", on_key)

    load_data()
    draw()
    plt.show()


if __name__ == "__main__":
    main()
