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
    _maxmin_calibration,
    _meansd_calibration,
    _constant_calibration,
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
        assert maxmin_targets.targets == [1.0, 0.2]
        assert meansd_targets.targets == [0.6, 0.03]
        assert constant_target.targets == [1.0]
        # test other attributes
        assert maxmin_targets.args == ["max", "min"]
        assert maxmin_targets.method == "maxmin"
        assert maxmin_targets.calibration_function == _maxmin_calibration

    def test_target_validation(self):
        # test wrong method
        with pytest.raises(ValueError):
            assign_targets(method="minmax", min=0.2, max=1.0)
        # test wrong kwargs
        with pytest.raises(ValueError):
            assign_targets(method="maxmin", max=1.0)
        with pytest.raises(ValueError):
            assign_targets(method="constant", max=1.0, min=0.2)
        with pytest.raises(ValueError):
            assign_targets(method="meansd", max=1.0, min=0.02)
        # test incorrect inputs
        with pytest.raises(ValueError):
            assign_targets(method="maxmin", max=0.2, min=1.0)
        with pytest.warns(UserWarning):
            assign_targets(method="meansd", mean=0.03, sd=0.6)


class TestCalibrate:
    def test_calibrate_maxmin(self):
        # try 1 fueltype and 1 parameter
        duet_run = import_duet(TMP_DIR, 252, 252)
        print(np.max(duet_run.depth[0, :, :]))
        grass_density = assign_targets(method="maxmin", max=1.0, min=0.2)
        grass_targets = assign_fuel_parameters(fuel_type="grass", density=grass_density)
        calibrated_duet = calibrate(duet_run, fuel_type_targets=grass_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.density, np.ndarray)
        assert np.allclose(calibrated_duet.depth, duet_run.depth)
        assert np.allclose(calibrated_duet.density, duet_run.density) == False
        assert np.allclose(calibrated_duet.density[1, :, :], duet_run.density[1, :, :])
        assert np.max(calibrated_duet.density[0, :, :]) == 1.0
        # assert np.min(calibrated_duet.density[0, :, :]) == 0.2 #why does this
        # try density and depth
        grass_depth = assign_targets(method="maxmin", max=1.0, min=0.2)
        grass_targets = assign_fuel_parameters(
            fuel_type="grass", density=grass_density, depth=grass_depth
        )
        # can't calibrate depth with maxmin because there's only one value!
        with pytest.raises(ValueError):
            calibrated_duet = calibrate(duet_run, fuel_type_targets=grass_targets)
        # now do moisture.. it will also raise an error since it doesn't exist in the og duet
        grass_moisture = assign_targets(method="maxmin", max=0.5, min=0.05)
        grass_targets = assign_fuel_parameters(
            fuel_type="grass", density=grass_density, moisture=grass_moisture
        )
        with pytest.raises(ValueError):
            calibrated_duet = calibrate(duet_run, fuel_type_targets=grass_targets)
