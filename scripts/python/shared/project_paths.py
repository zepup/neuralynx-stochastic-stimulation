"""
Central path helper for the inspection scripts.

Main purpose:
- let users set the Neuralynx data directory once with an environment variable
- reuse that path across the minimal Python scripts in this repository
"""

import os


def get_data_dir() -> str:
    data_dir = os.environ.get("NCS_PROJECT_DATA_DIR")
    if not data_dir:
        raise RuntimeError(
            "Set the NCS_PROJECT_DATA_DIR environment variable to the folder "
            "containing Neuralynx data, event files, and derived metadata."
        )
    return os.path.abspath(os.path.expanduser(data_dir))


def data_path(*parts: str) -> str:
    return os.path.join(get_data_dir(), *parts)
