"""
Test module for the calibration module of the duet_tools package.
"""
from __future__ import annotations

import pytest
from pathlib import Path

from duet_tools.calibration import (
    DuetRun,
    import_duet,
    assign_targets,
    fueltype_targets,
    calibrate,
    _validate_target_args,
    _do_calibration,
    _maxmin_calibration,
    _meansd_calibration,
    _constant_calibration,
    _moisture_weights_from_density,
    _truncate_at_0,
)

from duet_tools.utils import write_array_to_dat, read_dat_to_array

TEST_DIR = Path(__file__).parent
TMP_DIR = TEST_DIR / "tmp"
TMP_DIR.mkdir(exist_ok=True)
