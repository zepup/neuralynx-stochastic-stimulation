# Bipolar Stimulation Data Inspection

This project intends to conduct stochastic bipolar stimulation on patients using a system provided by Neuralynx.

The work belongs to Dr. Wael Asaad's lab.

Stimulation intensity, stimulation frequency, number of repetitions per stimulated pair, and the number of available contacts are defined by the stimulation configuration and contact-planning files associated with each session.

This repository focuses on inspection of Neuralynx stimulation recordings rather than downstream response analysis.

## How The Data Are Structured

The recording workflow usually contains the following file types:

- `.ncs`
  Neuralynx continuous recording files. These are the raw voltage recordings, usually one file per contact per recording block.
- `.nev`
  Neuralynx event files. These document event timing, including configuration loads, triggers, and stimulation-delivery events.
- `.cfg`
  Stimulation configuration files. These contain the contact pairs, stimulation intensities, and timing offsets that are read into Pegasus and executed as stimulation sequences.
- `_2mm.xlsx`
  Usually the anatomy file. This is the contact table used to map contact names and anatomy labels.
- `stimplan/` or planning tables
  These usually list which contacts are considered usable based on imaging constraints. For example, a pair should not be both in white matter and should not span across boxes.
- `stim_metadata.pkl`
  A project-specific metadata table used for inspection. This file is not generated automatically by Neuralynx. It is created by the custom metadata code in this repository.

Not every contact in the anatomy or planning sheets is ultimately stimulated. Final stimulated pairs can be reduced further after impedance testing.
Neuralynx recording files are also commonly capped at about 4 hours per file, so a single stimulation session can span multiple recording suffixes such as `0029` and `0030`.

A compact view of the structure is:

```text
session_folder/
├── Events_XXXX.nev or Events_XXXX.xlsx
├── CONTACT_A_XXXX.ncs
├── CONTACT_B_XXXX.ncs
├── CONTACT_C_XXXX.ncs
├── anatomy_2mm.xlsx
├── stimplan/
├── conditions/
│   ├── boxA/
│   │   └── config/
│   │       ├── Condition_1.cfg
│   │       └── Condition_2.cfg
│   └── boxB/
│       └── config/
├── BoxA_contact.csv
├── BoxB_contact.csv
└── stim_metadata.pkl
```

## What The Config Files Mean

The `.cfg` files define the stimulation sequence that Pegasus runs. In practice, they specify:

- which contacts form a bipolar pair
- stimulation intensity for each contact in the pair
- the within-sequence timing offsets
- the waveform file used by the sequence

These config files are essential for understanding what the recording is supposed to contain, because the `.nev` file documents when things happened, while the `.cfg` file documents what Pegasus was instructed to do.

## Timestamp Alignment And Why `stim_metadata` Exists

One of the main inspection problems is timestamp alignment.

The `.cfg` file does not live on the Neuralynx recording timeline by itself. It defines:

- which bipolar pair should be stimulated
- the intensity for that pair
- the offset of each stimulation entry within the sequence

The Neuralynx recording timeline lives in:

- `.nev`, which records events such as config loads, trigger echoes, and stimulation-delivery messages
- `.ncs`, which contains the raw recorded signal

To inspect data correctly, the stimulation instructions from the config file must be aligned to the recording timestamps from Neuralynx.

In this workflow, `stim_metadata.pkl` is the bridge between those two domains.

The metadata builder uses the ASHB export/event labels in the event file as timing markers to align:

- the active config file
- the sequence offsets inside that config
- the Neuralynx event timestamp
- the final expected stimulation timestamp on the recording timeline

Conceptually, the alignment works like this:

1. The ASHB / Pegasus-related event export in the event file identifies which config was active.
2. The event stream marks when an external trigger and stimulation-delivery event occurred.
3. The `.cfg` file provides the ordered bipolar pairs and their within-sequence offsets.
4. The metadata code combines the event timestamp anchor with the config offsets.
5. The result is `stim_metadata.pkl`, which assigns each expected stimulation pair a recording-aligned timestamp.

This matters because the config file alone is not enough to know where a stimulation should appear in the recording, and the `.nev` file alone is not enough to recover pair identity and sequence structure. The inspection-ready metadata file is what resolves that timestamp-alignment issue.

## Included Example Data

This repository includes a small Neuralynx test dataset as an example under:

- `examples/neuralynx_test/1 config loop/`

### `1 config loop`

This example contains:

- one config file: `Condition_1.cfg`
- one event file: `Events.nev`
- two recording channels: `LAMY1.ncs`, `LPAR1.ncs`

This example is intended to represent one configuration loop. The config file contains one stimulation sequence over Box A contacts, and the loop is repeated 20 times.

This example folder is included so users can inspect a known small dataset before opening a real patient dataset.

A static example plot generated from this folder is included under:

- `examples/plots/example_1_config_loop.svg`

## Repository Structure

The repository is arranged by purpose:

- `scripts/python/neuralynx_io/`
  Minimal Python reader for Neuralynx `.ncs` and `.nev` files.
- `scripts/python/shared/project_paths.py`
  Shared Python path helper for pointing the scripts to one data folder.
- `scripts/python/metadata/build_stim_metadata.py`
  Creates the custom `stim_metadata.pkl` file from `.nev`, `.cfg`, ASHB event markers, and contact-mapping files.
- `scripts/python/inspection/scroll_viewer.py`
  Interactive trace viewer for inspecting recordings and event markers.
- `scripts/python/inspection/browse_ccep_trials.py`
  Interactive browser for stepping through stimulation-aligned single trials using `stim_metadata.pkl`.
- `scripts/matlab/shared/project_paths.m`
  Shared MATLAB path helper.
- `scripts/matlab/io/preprocessMat.m`
  Minimal MATLAB entry point for loading `.ncs` and `.nev` files into one structure.
- `scripts/matlab/io/pullTimeLabels.m`
  Helper for matching requested clock times to Neuralynx recording suffixes.
- `examples/neuralynx_test/`
  Small bench-style example dataset included in the repository.
- `examples/plots/`
  Static SVG overview plots generated from the included Neuralynx test data.

## Dataset-Specific Values To Replace

Before sharing or reusing the repository, users should replace patient- and session-specific values inside the retained scripts.

Examples include:

- `PATIENT_ID` labels in viewer titles or script comments
- recording suffixes such as `0029`
- event filenames such as `Events_0029.nev`
- selected inspection channels
- session date and optional clock-time crop windows

The repository intentionally masks the patient identifier while leaving other session values in place as working examples.

## Minimum Code Needed To Read And Inspect Data

The minimal files needed for inspection are:

- `scripts/python/neuralynx_io/neuralynx_io.py`
- `scripts/python/inspection/scroll_viewer.py`
- `scripts/python/inspection/browse_ccep_trials.py`
- `scripts/python/shared/project_paths.py`

If metadata must be rebuilt, also use:

- `scripts/python/metadata/build_stim_metadata.py`

If a MATLAB loading path is preferred, use:

- `scripts/matlab/io/preprocessMat.m`
- `scripts/matlab/io/pullTimeLabels.m`
- `scripts/matlab/shared/project_paths.m`

## Required Packages

### Python

The retained Python scripts expect:

- `numpy`
- `scipy`
- `pandas`
- `matplotlib`
- `pyqtgraph`

Notes:

- `neuralynx_io.py` works with `numpy`
- `build_stim_metadata.py` requires `pandas` and `scipy`
- `browse_ccep_trials.py` requires `numpy`, `pandas`, and `matplotlib`
- `scroll_viewer.py` requires `pyqtgraph`

### MATLAB

The retained MATLAB scripts expect:

- MATLAB
- FieldTrip, with Neuralynx-compatible readers available through `ft_read_header`, `ft_read_data`, and `ft_read_event`

### MATLAB Path Setup

The MATLAB scripts are set up to read the bundled example data first. The example folder contains only:

- `LAMY1.ncs`
- `LPAR1.ncs`
- `Events.nev`
- `Condition_1.cfg`

It does not require an anatomy spreadsheet. This is intended only to help students and collaborators become familiar with the Neuralynx file format.

Example-data usage:

```matlab
addpath('/path/to/repository/scripts/matlab/io')
setenv('FIELDTRIP_DIR', '/path/to/fieldtrip')

collData = preprocessMat();
```

By default, `project_paths.m` points to:

```text
examples/neuralynx_test/1 config loop
```

Optional real-data override:

```matlab
setenv('NCS_PROJECT_DATA_DIR', '/path/to/PatientData/byPatient/Epilepsy/PATIENT_ID/2026-03-15/stochasticStim')
setenv('NCS_PROJECT_ANATOMY_FILE', '/path/to/PATIENT_ID_2mm.xlsx')

collData = preprocessMat('PATIENT_ID', '2026-03-15', 'HH:MM:SS', 'HH:MM:SS');
```

Optional environment variables for real datasets:

- `NCS_PROJECT_OUTPUT_DIR`: where `.mat` and optional `.h5` outputs are written. Defaults to the data folder.
- `NCS_PROJECT_ANATOMY_FILE`: full path to the anatomy spreadsheet. The example data do not use this.

## Minimum Steps For A Real Shared Dataset

For a real dataset shared outside this repository, the minimum steps depend on whether `stim_metadata.pkl` already exists.

### If `stim_metadata.pkl` is already present in the data folder

This is the minimum needed to start inspection:

1. Put the data outside the repository in one local folder.
2. Set:

```bash
export NCS_PROJECT_DATA_DIR="/path/to/real/data"
```

3. Open the Python viewer on the desired recording block:

```bash
python scripts/python/inspection/scroll_viewer.py patient
```

or point to a subfolder:

```bash
python scripts/python/inspection/scroll_viewer.py folder session_subfolder
```

In this case, users do not need to rebuild metadata before basic inspection if the shared folder already contains the metadata file they need.

In practice, if `stim_metadata.pkl` is already present, that means the config-to-recording timestamp alignment work has already been done.
Users can then either inspect the full recording timeline with `scroll_viewer.py` or inspect stimulation-aligned single trials with `browse_ccep_trials.py`.

### If `stim_metadata.pkl` is not present

Run:

```bash
python scripts/python/metadata/build_stim_metadata.py
```

before inspection.

## Raw Data Policy

Real patient raw data should still remain outside the repository. The included example Neuralynx test data are small bench-style examples for documentation and inspection only.

Do not commit:

- raw patient `.ncs` recordings
- raw patient `.nev` files
- identifying anatomy spreadsheets
- patient-specific planning or metadata files

## Notes For Reuse

- `build_stim_metadata.py` creates metadata by aligning ASHB / Pegasus event markers with config-defined offsets on the Neuralynx recording timeline; Neuralynx does not provide that final inspection-ready table automatically.
- Session-specific values such as stimulation intensity, frequency, repetition count, and usable contacts come from the session's config, anatomy, planning, and impedance-filtered stimulation records.
- The retained scripts contain placeholders such as subject ID, recording suffix, and example channel names that should be edited for each dataset.
- Users should insert their own patient or session label where placeholder text such as `PATIENT_ID` appears; real patient identifiers should not be hard-coded into the repository.
