#!/usr/bin/env python3
"""
Generate static example plots from the included Neuralynx test datasets.

Main purpose:
- create simple overview figures for the repository README/examples
- show raw traces from the included test recordings
- mark event timing from the matching `.nev` files
"""

import os
import struct
from pathlib import Path


HEADER_LENGTH = 16 * 1024
NCS_RECORD_BYTES = 1044
NCS_SAMPLES_PER_RECORD = 512
NCS_STRUCT = struct.Struct("<QIIi512h")
NEV_RECORD_BYTES = 182
NEV_STRUCT = struct.Struct("<hhhQhhhhiiiiiiii128s")


REPO_ROOT = Path(__file__).resolve().parents[3]
EXAMPLE_ROOT = REPO_ROOT / "examples" / "neuralynx_test"
PLOTS_DIR = REPO_ROOT / "examples" / "plots"


def read_ncs(path):
    times = []
    data = []
    with open(path, "rb") as handle:
        handle.seek(HEADER_LENGTH)
        while True:
            raw = handle.read(NCS_RECORD_BYTES)
            if len(raw) != NCS_RECORD_BYTES:
                break
            rec = NCS_STRUCT.unpack(raw)
            ts = rec[0]
            fs = rec[2]
            nvalid = rec[3]
            samples = rec[4:4 + max(0, nvalid)]
            base_time = ts
            for idx, sample in enumerate(samples):
                times.append(base_time + int(idx * (1e6 / fs)))
                data.append(sample)
    if not times:
        return [], []
    t0 = times[0]
    times_s = [(t - t0) / 1e6 for t in times]
    return times_s, data


def read_nev_events(path):
    config_loads = []
    stim_begins = []
    first_ts = None
    with open(path, "rb") as handle:
        handle.seek(HEADER_LENGTH)
        while True:
            raw = handle.read(NEV_RECORD_BYTES)
            if len(raw) != NEV_RECORD_BYTES:
                break
            rec = NEV_STRUCT.unpack(raw)
            ts = rec[3]
            label = rec[-1].decode("latin1", errors="replace").rstrip("\x00")
            if first_ts is None:
                first_ts = ts
            t_rel = (ts - first_ts) / 1e6
            if "ProcessConfigurationFile" in label and "Echo" in label:
                config_loads.append(t_rel)
            if "Stimulation delivery begins" in label:
                stim_begins.append(t_rel)
    return config_loads, stim_begins


def downsample(xs, ys, max_points=3000):
    if len(xs) <= max_points:
        return xs, ys
    step = max(1, len(xs) // max_points)
    return xs[::step], ys[::step]


def scale_points(xs, ys, left, top, width, height):
    x_min = min(xs)
    x_max = max(xs)
    y_min = min(ys)
    y_max = max(ys)
    if x_max == x_min:
        x_max = x_min + 1.0
    if y_max == y_min:
        y_max = y_min + 1.0
    pts = []
    for x, y in zip(xs, ys):
        sx = left + (x - x_min) / (x_max - x_min) * width
        sy = top + height - (y - y_min) / (y_max - y_min) * height
        pts.append(f"{sx:.2f},{sy:.2f}")
    return " ".join(pts), x_min, x_max


def line_x(ts, x_min, x_max, left, width):
    if x_max == x_min:
        return left
    return left + (ts - x_min) / (x_max - x_min) * width


def write_svg(out_path, title, traces, config_loads, stim_begins):
    panel_height = 180
    header_height = 70
    width = 1200
    height = header_height + panel_height * len(traces)
    left = 90
    right = 40
    plot_width = width - left - right
    plot_height = 110

    all_x = []
    for _, xs, _ in traces:
        all_x.extend(xs)
    if all_x:
        global_x_min = min(all_x)
        global_x_max = max(all_x)
    else:
        marker_times = config_loads + stim_begins
        global_x_min = 0.0
        global_x_max = max(marker_times) if marker_times else 1.0

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="28" text-anchor="middle" font-size="20" font-family="Arial">{title}</text>',
        '<text x="90" y="50" font-size="12" font-family="Arial" fill="#444">Orange dashed = config load, green dotted = stimulation begins</text>',
    ]

    for idx, (channel_name, xs, ys) in enumerate(traces):
        panel_top = header_height + idx * panel_height
        xs_ds, ys_ds = downsample(xs, ys)
        pts, _, _ = scale_points(xs_ds, ys_ds, left, panel_top, plot_width, plot_height)

        lines.append(f'<text x="18" y="{panel_top + plot_height/2:.1f}" font-size="12" font-family="Arial">{channel_name}</text>')
        lines.append(f'<rect x="{left}" y="{panel_top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#cccccc"/>')

        for ts in config_loads:
            x = line_x(ts, global_x_min, global_x_max, left, plot_width)
            lines.append(f'<line x1="{x:.2f}" y1="{panel_top}" x2="{x:.2f}" y2="{panel_top + plot_height}" stroke="#d97706" stroke-width="1" stroke-dasharray="5,4"/>')
        for ts in stim_begins:
            x = line_x(ts, global_x_min, global_x_max, left, plot_width)
            lines.append(f'<line x1="{x:.2f}" y1="{panel_top}" x2="{x:.2f}" y2="{panel_top + plot_height}" stroke="#15803d" stroke-width="1" stroke-dasharray="2,4"/>')

        lines.append(f'<polyline fill="none" stroke="black" stroke-width="0.8" points="{pts}"/>')
        lines.append(f'<text x="{left}" y="{panel_top + plot_height + 18}" font-size="11" font-family="Arial" fill="#555">0 s</text>')
        lines.append(f'<text x="{left + plot_width}" y="{panel_top + plot_height + 18}" text-anchor="end" font-size="11" font-family="Arial" fill="#555">{global_x_max:.1f} s</text>')

    lines.append("</svg>")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def build_plot(dataset_dir, channel_stems, title, output_name):
    config_loads, stim_begins = read_nev_events(dataset_dir / "Events.nev")
    traces = []
    for stem in channel_stems:
        matches = sorted(dataset_dir.glob(f"{stem}*.ncs"))
        if not matches:
            continue
        xs, ys = read_ncs(matches[0])
        if xs and ys:
            traces.append((stem, xs, ys))
    if not traces:
        traces.append(("Event timeline only", [0.0, max(config_loads + stim_begins + [1.0])], [0.0, 0.0]))
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = PLOTS_DIR / output_name
    write_svg(out_path, title, traces, config_loads, stim_begins)
    return out_path


def main():
    outputs = [
        build_plot(
            EXAMPLE_ROOT / "1 config loop",
            ["LAMY1", "LPAR1"],
            "Neuralynx test example: 1 config loop",
            "example_1_config_loop.svg",
        ),
    ]
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
