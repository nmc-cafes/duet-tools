"""
Test module for the calibration module of the duet_tools package.
"""
from __future__ import annotations

import pytest
import numpy as np
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
        litter_moisture = assign_targets(0.1, method="constant")
        litter_density = assign_targets(0.6, 0.05, method="meansd")
        litter_depth = assign_targets(0.3, 0.01, method="maxmin")
        # test data types
        assert isinstance(litter_moisture, dict)
        assert isinstance(litter_density, dict)
        assert isinstance(litter_depth, dict)
        # test dict structure
        assert list(litter_moisture.keys()) == ["method", "constant"]
        assert list(litter_moisture.values()) == ["constant", 0.1]
        assert list(litter_density.keys()) == ["method", "mean", "sd"]
        assert list(litter_density.values()) == ["meansd", 0.6, 0.05]
        assert list(litter_depth.keys()) == ["method", "max", "min"]
        assert list(litter_depth.values()) == ["maxmin", 0.3, 0.01]
        # test kwarg instead of arg
        with pytest.raises(TypeError):
            litter_moisture = assign_targets(constant=0.1, method="constant")
        # test wrong number of arguments
        with pytest.raises(ValueError):
            litter_moisture = assign_targets(0.1, 0.2, method="constant")
        with pytest.raises(ValueError):
            litter_moisture = assign_targets(0.1, method="meansd")
        with pytest.raises(ValueError):
            litter_moisture = assign_targets(0.1, 0.2, 0.3, method="maxmin")


class TestFueltypeTargets:
    def test_fueltype_targets(self):
        litter_moisture = assign_targets(0.1, method="constant")
        litter_density = assign_targets(0.6, 0.05, method="meansd")
        litter_depth = assign_targets(0.3, 0.01, method="maxmin")
        litter_targets = fueltype_targets(
            density=litter_density, moisture=litter_moisture, depth=litter_depth
        )
        # test data type
        assert isinstance(litter_targets, dict)
        # test dict structure
        assert list(litter_targets.keys()) == ["density", "moisture", "depth"]
        assert list(litter_targets["density"].keys()) == ["method", "mean", "sd"]
        # test not having targets for some parameters
        litter_targets = fueltype_targets(density=litter_density)
        assert len(litter_targets.keys()) == 1
        litter_targets = fueltype_targets(density=litter_density, depth=litter_depth)
        assert len(litter_targets.keys()) == 2


class TestCalibrate:
    def test_calibrate_maxmin(self):
        # try 1 fueltype and 1 parameter
        duet_run = import_duet(TMP_DIR, 252, 252)
        grass_density = assign_targets(1.0, 0.3, method="maxmin")
        grass_targets = fueltype_targets(density=grass_density)
        calibrated_duet = calibrate(duet_run, grass=grass_targets)
        assert isinstance(calibrated_duet, DuetRun)
