"""
Build project-specific stimulation metadata from event files and config files.

Main purpose:
- parse the active condition from Neuralynx event logs
- read stimulation pairing/intensity information from config files
- create a `.pkl` / `.mat` metadata table used for inspection

Dataset-specific fields near the top of this file are left as editable examples.
Users should adapt them for their own local dataset while avoiding hard-coded public patient identifiers.
"""

import os
import re
import numpy as np
import pandas as pd
import scipy.io
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, ".."))
from shared.project_paths import data_path, get_data_dir

BASE_DIR = get_data_dir()

# Event file in tabular form.
# Example value below is intentionally left as a working session example.
EVENTS_TABLE = data_path("Events_0029.csv")

# Where the condition CFGs are stored.
# Example values below are intentionally left as a working session example.
CFG_BASE = {
    "boxA": data_path("configs", "PATIENT_ID_boxA", "config"),
    "boxB": data_path("configs", "PATIENT_ID_boxB_A", "config"),
}

SUBJECT = "SUBJECT_ID"



# Lines in config files:
# -SetStimulationSequenceEntry "13" 0  B:\mltask\stochasticStimV2\e0043_boxA\wavFiles\160us_1Hz.wav 3 10
CFG_LINE_RE = re.compile(
    r'-SetStimulationSequenceEntry\s+"(?P<contact>\d+)"\s+'
    r'(?P<offset>\d+)\s+(?P<wavpath>\S+)\s+(?P<chan>\d+)\s+(?P<intensity>-?\d+)'
)

# In the event label, find which box (boxA/boxB) and which condition
BOX_RE = re.compile(r'box([ABab])')
COND_RE = re.compile(r'(Condition_\d+)\.cfg')

# From events filename, extract recording suffix: Events_0029.csv -> "0029"
SUFFIX_RE = re.compile(r'Events_(\d+)', re.IGNORECASE)


# helpers

def load_condition_pairs(box: str, condition: str):
    """
    Given 'boxA' or 'boxB' and a condition name 'Condition_1',
    load the corresponding CFG file and return a list of bipole pairs.

    Each pair dict has:
        pair_idx, contact1, contact2, offset_us, intensity1, intensity2, wavpath
    """
    cfg_dir = CFG_BASE[box]
    cfg_path = os.path.join(cfg_dir, f"{condition}.cfg")

    if not os.path.isfile(cfg_path):
        raise FileNotFoundError(f"Condition file not found: {cfg_path}")

    entries = []
    with open(cfg_path, "r") as f:
        for line in f:
            m = CFG_LINE_RE.match(line.strip())
            if not m:
                continue
            gd = m.groupdict()
            entries.append(
                {
                    "contact": int(gd["contact"]),
                    "offset_us": int(gd["offset"]),  # µs relative to run start
                    "wavpath": gd["wavpath"],
                    "chan": int(gd["chan"]),
                    "intensity": int(gd["intensity"]),
                }
            )

    if len(entries) % 2 != 0:
        print(f"WARNING: odd number of sequence entries in {cfg_path}; dropping last one.")
        entries = entries[:-1]

    pairs = []
    for i in range(0, len(entries), 2):
        e1 = entries[i]
        e2 = entries[i + 1]
        pair_idx = i // 2

        # offsets should match within a pair; use e1's offset
        offset_us = e1["offset_us"]

        pairs.append(
            {
                "pair_idx": pair_idx,
                "contact1": e1["contact"],
                "contact2": e2["contact"],
                "offset_us": offset_us,
                "intensity1": e1["intensity"],
                "intensity2": e2["intensity"],
                "wavpath": e1["wavpath"],
            }
        )

    return pairs


def load_events_table(path):
    """
    Load the event table from either CSV or XLSX.
    Real example sessions may use CSV exports; some users may keep XLSX.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        return pd.read_csv(path)
    if ext in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    raise ValueError(f"Unsupported event table format: {path}")


def build_stim_metadata():
    """
    Reads the selected event table, identifies which cfg file is active using
    ASHB NetCom Command Echo lines, and generates a stim entry for each
    bipole pair in each run.

    Trigger anchor: 'Stimulation delivery begins' (not the earlier Echo event).
    Runs where 'Stimulation delivery failed' appears are excluded entirely.

    Returns a DataFrame with one row per bipole pair per run.
    """
    events = load_events_table(EVENTS_TABLE)
    events = events.sort_values("timestamp_us").reset_index(drop=True)

    condition_cache = {}  # (box, condition) -> list of pairs
    stim_rows = []

    current_box = None        # 'boxA' or 'boxB'
    current_condition = None  # 'Condition_1', etc.
    run_idx = 0

    # State for pairing Echo → Delivery begins / failed
    pending_box = None
    pending_condition = None
    pending_trigger_ts = None  # timestamp from ExternalTrigger Echo

    for _, row in events.iterrows():
        label = str(row["label"])
        ts_us = int(row["timestamp_us"])

        # 1) Config selection: which condition is now active?
        if "ASHB NetCom Command Echo" in label and "ProcessConfigurationFile" in label:
            box_match = BOX_RE.search(label)
            cond_match = COND_RE.search(label)

            if box_match and cond_match:
                box_letter = box_match.group(1).upper()
                current_box = f"box{box_letter}"
                current_condition = cond_match.group(1)
                key = (current_box, current_condition)
                if key not in condition_cache:
                    condition_cache[key] = load_condition_pairs(current_box, current_condition)
            else:
                current_box = None
                current_condition = None
            pending_box = None
            pending_condition = None
            pending_trigger_ts = None

        # 2) External trigger echo: arm pending state AND capture timestamp
        elif "ASHB ExternalTrigger Echo: 1 1" in label:
            if current_box is not None and current_condition is not None:
                pending_box = current_box
                pending_condition = current_condition
                pending_trigger_ts = ts_us  # anchor: when the TTL trigger arrived
            else:
                pending_box = None
                pending_condition = None
                pending_trigger_ts = None

        # 3) Delivery failed: discard this run
        elif "Stimulation delivery failed" in label:
            pending_box = None
            pending_condition = None
            pending_trigger_ts = None

        # 4) Delivery begins: confirm delivery; use Echo timestamp as anchor
        elif "Stimulation delivery begins" in label:
            if pending_box is None or pending_condition is None or pending_trigger_ts is None:
                continue

            key = (pending_box, pending_condition)
            pairs = condition_cache.get(key, [])
            if not pairs:
                pending_box = None
                pending_condition = None
                pending_trigger_ts = None
                continue

            run_idx += 1
            trigger_ts = pending_trigger_ts  # ExternalTrigger Echo timestamp

            for p in pairs:
                stim_ts = trigger_ts + p["offset_us"]
                stim_rows.append(
                    {
                        "subject": SUBJECT,
                        "box": pending_box,
                        "condition": pending_condition,
                        "run_idx": run_idx,
                        "pair_idx": p["pair_idx"],
                        "contact1": p["contact1"],
                        "contact2": p["contact2"],
                        "intensity1": p["intensity1"],
                        "intensity2": p["intensity2"],
                        "offset_us": p["offset_us"],
                        "trigger_timestamp_us": trigger_ts,
                        "stim_timestamp_us": stim_ts,
                        "wavpath": p["wavpath"],
                    }
                )

            pending_box = None
            pending_condition = None
            pending_trigger_ts = None

    stim_df = pd.DataFrame(stim_rows)
    return stim_df


def load_contact_maps():
    """
    Load BoxA_contact.csv and BoxB_contact.csv and build maps:
        ("boxA", contact_num) -> electrode_name
        ("boxB", contact_num) -> electrode_name
    Uses C1/C2 -> C1_electrode/C2_electrode.
    """
    contact_map = {}

    for box in ["boxA", "boxB"]:
        csv_path = data_path(f"{box.capitalize()}_contact.csv")
        if not os.path.isfile(csv_path):
            print(f"WARNING: contact file missing: {csv_path}")
            continue

        df = pd.read_csv(csv_path)
        df.columns = df.columns.str.strip()

        # Expect columns: Lead, C1, C2, C1_electrode, C2_electrode
        for _, r in df.iterrows():
            c1 = r.get("C1")
            c2 = r.get("C2")
            e1 = r.get("C1_electrode")
            e2 = r.get("C2_electrode")

            if pd.notna(c1) and pd.notna(e1):
                contact_map[(box, int(c1))] = str(e1)
            if pd.notna(c2) and pd.notna(e2):
                contact_map[(box, int(c2))] = str(e2)

    return contact_map


def infer_recording_suffix():
    """
    Extract recording suffix from EVENTS_XLSX filename, e.g.
    Events_0036.xlsx -> '0036'
    """
    base = os.path.basename(EVENTS_TABLE)
    m = SUFFIX_RE.search(base)
    if not m:
        raise ValueError(f"Could not parse recording suffix from {base}")
    return m.group(1)


def add_electrodes_and_paths(stim_df):
    """
    Adds electrode1, electrode2, ncs_path1, ncs_path2 to stim_df.

    Electrode names come from BoxA_contact.csv / BoxB_contact.csv.
    NCS paths are built as:
        <BASE_DIR>/<electrode>_<suffix>.ncs
    where <suffix> comes from Events_XXXX.xlsx.
    """
    contact_map = load_contact_maps()
    recording_suffix = infer_recording_suffix()
    ncs_dir = BASE_DIR

    electrode1 = []
    electrode2 = []
    ncs1 = []
    ncs2 = []

    for _, row in stim_df.iterrows():
        box = row["box"]  # 'boxA' or 'boxB'
        c1 = int(row["contact1"])
        c2 = int(row["contact2"])

        e1 = contact_map.get((box, c1), None)
        e2 = contact_map.get((box, c2), None)

        electrode1.append(e1)
        electrode2.append(e2)

        if e1 is not None:
            ncs1.append(os.path.join(ncs_dir, f"{e1}_{recording_suffix}.ncs"))
        else:
            ncs1.append(None)

        if e2 is not None:
            ncs2.append(os.path.join(ncs_dir, f"{e2}_{recording_suffix}.ncs"))
        else:
            ncs2.append(None)

    stim_df["electrode1"] = electrode1
    stim_df["electrode2"] = electrode2
    stim_df["ncs_path1"] = ncs1
    stim_df["ncs_path2"] = ncs2

    return stim_df




if __name__ == "__main__":
    # Build raw stim metadata from events + cfgs
    stim_df = build_stim_metadata()

    # Add electrode names + NCS file paths
    stim_df = add_electrodes_and_paths(stim_df)

    print(stim_df.head())
    print("Total stim rows:", len(stim_df))

    out_base = data_path("stim_metadata")
    stim_df.to_pickle(out_base + ".pkl")
    stim_df.to_excel(out_base + ".xlsx", index=False)

    # Save stim_metadata_full.mat for MATLAB (browse_ccep_trials.m)
    mat_dict = {
        "stim_timestamp_us": stim_df["stim_timestamp_us"].to_numpy(dtype=np.float64),
        "trigger_timestamp_us": stim_df["trigger_timestamp_us"].to_numpy(dtype=np.float64),
        "electrode1": np.array(stim_df["electrode1"].fillna("").tolist(), dtype=object),
        "electrode2": np.array(stim_df["electrode2"].fillna("").tolist(), dtype=object),
        "intensity1": stim_df["intensity1"].to_numpy(dtype=np.float64),
        "intensity2": stim_df["intensity2"].to_numpy(dtype=np.float64),
        "contact1":   stim_df["contact1"].to_numpy(dtype=np.float64),
        "contact2":   stim_df["contact2"].to_numpy(dtype=np.float64),
        "run_idx":    stim_df["run_idx"].to_numpy(dtype=np.float64),
        "pair_idx":   stim_df["pair_idx"].to_numpy(dtype=np.float64),
        "offset_us":  stim_df["offset_us"].to_numpy(dtype=np.float64),
    }
    mat_path = data_path("stim_metadata_full.mat")
    scipy.io.savemat(mat_path, mat_dict)

    print("Saved:")
    print("  ", out_base + ".pkl")
    print("  ", out_base + ".xlsx")
    print("  ", mat_path)
