"""
Test module for the calibration module of the duet_tools package.
"""
from __future__ import annotations

import pytest
import numpy as np
from pathlib import Path

from duet_tools.calibration import (
    DuetRun,
    Targets,
    FuelParameters,
    LandfireQuery,
    import_duet,
    assign_targets,
    assign_fuel_parameters,
    calibrate,
)

from duet_tools.utils import write_array_to_dat, read_dat_to_array

TEST_DIR = Path(__file__).parent
TMP_DIR = TEST_DIR / "tmp"
TMP_DIR.mkdir(exist_ok=True)


class TestDuetRun:
    def test_import_duet(self):
        duet_run = import_duet(directory=TMP_DIR, nx=252, ny=252)
        # test that data types are correct
        assert isinstance(duet_run, DuetRun)
        assert isinstance(duet_run.density, np.ndarray)
        assert isinstance(duet_run.depth, np.ndarray)
        assert duet_run.moisture is None
        # test array shapes
        assert duet_run.density.shape == (2, 252, 252)
        assert duet_run.depth.shape == (2, 252, 252)
        # test that wrong dimensions raise error
        with pytest.raises(ValueError):
            duet_run = import_duet(directory=TMP_DIR, nx=252, ny=252, nz=3)


class TestAssignTargets:
    def test_assign_targets(self):
        maxmin_targets = assign_targets(method="maxmin", max=1.0, min=0.2)
        meansd_targets = assign_targets(method="meansd", mean=0.6, sd=0.03)
        constant_target = assign_targets(method="constant", target=1.0)
        assert isinstance(maxmin_targets, Targets)
        assert isinstance(meansd_targets, Targets)
        assert isinstance(constant_target, Targets)


class TestCalibrate:
    def test_calibrate_maxmin(self):
        # try 1 fueltype and 1 parameter
        duet_run = import_duet(TMP_DIR, 252, 252)
        grass_density = assign_targets(1.0, 0.3, method="maxmin")
        grass_targets = assign_fuel_parameters(density=grass_density)
        calibrated_duet = calibrate(duet_run, grass=grass_targets)
        assert isinstance(calibrated_duet, DuetRun)
