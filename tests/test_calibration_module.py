"""
Test module for the calibration module of the duet_tools package.
"""

from __future__ import annotations

import pytest
import numpy as np
from pathlib import Path
import math
import matplotlib.pyplot as plt

from duet_tools.calibration import (
    DuetRun,
    Targets,
    FuelParameter,
    import_duet,
    assign_targets,
    set_fuel_parameter,
    set_density,
    set_moisture,
    set_height,
    calibrate,
    _maxmin_calibration,
)

from duet_tools.utils import read_dat_to_array

TEST_DIR = Path(__file__).parent
TMP_DIR = TEST_DIR / "tmp"
TMP_DIR.mkdir(exist_ok=True)
DATA_DIR = TEST_DIR / "test-data"


class TestDuetRun:
    def test_import_duet_v1(self):
        duet_run = import_duet(
            directory=DATA_DIR / "v1", nx=976, ny=1998, nsp=2, version="v1"
        )
        # test that data types are correct
        assert isinstance(duet_run, DuetRun)
        assert isinstance(duet_run.density, np.ndarray)
        assert isinstance(duet_run.moisture, np.ndarray)
        assert isinstance(duet_run.height, np.ndarray)
        # test array shapes
        assert duet_run.density.shape == (2, 1998, 976)
        assert duet_run.moisture.shape == (2, 1998, 976)
        assert duet_run.height.shape == (2, 1998, 976)
        # test that wrong dimensions raises error
        with pytest.raises(ValueError):
            duet_run = import_duet(
                directory=DATA_DIR / "v1", nx=976, ny=1998, nsp=3, version="v1"
            )
        # test that wrong version number raises error
        with pytest.raises(ValueError):
            duet_run = import_duet(
                directory=DATA_DIR / "v1", nx=976, ny=1998, nsp=2, version="v3"
            )

    def test_import_duet_v2(self):
        duet_run = import_duet(
            directory=DATA_DIR / "v2", nx=333, ny=295, nsp=9, version="v2"
        )
        # test that data types are correct
        assert isinstance(duet_run, DuetRun)
        assert isinstance(duet_run.density, np.ndarray)
        assert isinstance(duet_run.moisture, np.ndarray)
        assert isinstance(duet_run.height, np.ndarray)
        # test array shapes
        assert duet_run.density.shape == (2, 295, 333)
        assert duet_run.moisture.shape == (2, 295, 333)
        assert duet_run.height.shape == (2, 295, 333)
        # test that wrong dimensions raise error
        with pytest.raises(ValueError):
            duet_run = import_duet(directory=DATA_DIR / "v2", nx=333, ny=295, nsp=3)

    def test_to_numpy(self):
        duet_run = import_duet(directory=DATA_DIR / "v2", nx=333, ny=295, nsp=9)
        # test each fuel parameter and type
        grass_density = duet_run.to_numpy("grass", "density")
        litter_density = duet_run.to_numpy("litter", "density")
        grass_moisture = duet_run.to_numpy("grass", "moisture")
        litter_moisture = duet_run.to_numpy("litter", "moisture")
        grass_height = duet_run.to_numpy("grass", "height")
        litter_height = duet_run.to_numpy("litter", "height")
        assert np.array_equal(grass_density, duet_run.density[0, :, :])
        assert np.array_equal(litter_density, duet_run.density[1, :, :])
        assert np.array_equal(grass_moisture, duet_run.moisture[0, :, :])
        assert np.array_equal(litter_moisture, duet_run.moisture[1, :, :])
        assert np.array_equal(grass_height, duet_run.height[0, :, :])
        assert np.array_equal(litter_height, duet_run.height[1, :, :])
        # test integrated and separated
        separated_density = duet_run.to_numpy("separated", "density")
        assert np.array_equal(separated_density, duet_run.density)
        integrated_density = duet_run.to_numpy("integrated", "density")
        integrated_height = duet_run.to_numpy("integrated", "height")
        integrated_moisture = duet_run.to_numpy("integrated", "moisture")
        assert np.array_equal(integrated_density, np.sum(duet_run.density, axis=0))
        assert np.array_equal(integrated_height, np.max(duet_run.height, axis=0))
        weights = _maxmin_calibration(duet_run.density.copy(), max=1.0, min=0)
        weights[weights == 0] = 0.01
        masked = np.ma.masked_array(duet_run.moisture, duet_run.moisture == 0)
        averaged = np.ma.average(masked, axis=0, weights=weights)
        weighted_average_moisture = np.ma.filled(averaged, 0)
        assert np.array_equal(integrated_moisture, weighted_average_moisture)

    def test_to_quicfire(self):
        duet_run = import_duet(directory=DATA_DIR / "v2", nx=333, ny=295, nsp=9)
        duet_run.to_quicfire(TMP_DIR, overwrite=True)
        with pytest.raises(FileExistsError):
            duet_run.to_quicfire(TMP_DIR)
        duet_run.to_quicfire(TMP_DIR, overwrite=True)
        treesrhof = read_dat_to_array(TMP_DIR, "treesrhof.dat", nx=333, ny=295, nz=1)
        treesrhof = treesrhof[0, :, :]
        treesmoist = read_dat_to_array(
            TMP_DIR,
            "treesmoist.dat",
            nx=333,
            ny=295,
            nz=1,
        )
        treesmoist = treesmoist[0, :, :]
        treesfueldepth = read_dat_to_array(
            TMP_DIR,
            "treesfueldepth.dat",
            nx=333,
            ny=295,
            nz=1,
        )
        treesfueldepth = treesfueldepth[0, :, :]
        assert np.array_equal(
            treesrhof, duet_run._integrate("density").astype("float32")
        )
        assert np.array_equal(
            treesfueldepth, duet_run._integrate("height").astype("float32")
        )
        assert np.array_equal(
            treesmoist, duet_run._integrate("moisture").astype("float32")
        )

        files = ["treesrhof.dat", "treesmoist.dat", "treesfueldepth.dat"]
        for file in files:
            path = TMP_DIR / file
            path.unlink()
        duet_run.to_quicfire(TMP_DIR, density=False, moisture=False, height=False)
        for file in files:
            path = TMP_DIR / file
            assert path.exists() == False


class TestAssignTargets:
    def test_assign_targets(self):
        maxmin_targets = assign_targets(method="maxmin", max=1.0, min=0.2)
        meansd_targets = assign_targets(method="meansd", mean=0.6, sd=0.03)
        constant_target = assign_targets(method="constant", value=1.0)
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


class TestSetFuelParameter:
    def test_set_fuel_parameter(self):
        maxmin_targets = assign_targets(method="maxmin", max=1.0, min=0.2)
        meansd_targets = assign_targets(method="meansd", mean=0.6, sd=0.03)
        constant_target = assign_targets(method="constant", value=1.0)
        # test separated functions
        density_targets = set_density(grass=maxmin_targets, litter=meansd_targets)
        moisture_targets = set_moisture(grass=constant_target)
        height_targets = set_height(all=maxmin_targets)
        assert isinstance(density_targets, FuelParameter)
        assert isinstance(moisture_targets, FuelParameter)
        assert isinstance(height_targets, FuelParameter)
        assert density_targets.parameter == "density"
        assert moisture_targets.parameter == "moisture"
        # test generic function
        density_targets = set_fuel_parameter(
            parameter="density", grass=maxmin_targets, litter=meansd_targets
        )
        assert isinstance(density_targets, FuelParameter)
        assert density_targets.parameter == "density"
        assert density_targets.fuel_types == ["grass", "litter"]
        assert density_targets.targets == [maxmin_targets, meansd_targets]
        # test validation
        with pytest.raises(ValueError):
            density_targets = set_fuel_parameter(
                parameter="bulk density", litter=meansd_targets
            )
        with pytest.raises(ValueError):
            density_targets = set_fuel_parameter(
                parameter="density", both=meansd_targets
            )
        with pytest.raises(ValueError):
            density_targets = set_fuel_parameter(
                parameter="density", litter=meansd_targets, all=constant_target
            )


class TestCalibrate:
    def test_calibrate_maxmin(self):
        # try 1 fueltype and 1 parameter
        duet_run = import_duet(DATA_DIR / "v2", 295, 333, 9)
        grass_density = assign_targets(method="maxmin", max=1.0, min=0.2)
        density_targets = set_fuel_parameter(parameter="density", grass=grass_density)
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=density_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.density, np.ndarray)
        assert np.allclose(calibrated_duet.height, duet_run.height)
        assert np.allclose(calibrated_duet.density, duet_run.density) == False
        assert np.allclose(calibrated_duet.density[1, :, :], duet_run.density[1, :, :])
        assert np.max(calibrated_duet.density[0, :, :]) == 1.0
        assert np.min(calibrated_duet.density[0, :, :]) == 0.2
        # try density and height
        grass_height = assign_targets(method="maxmin", max=1.0, min=0.2)
        density_targets = set_fuel_parameter(parameter="density", grass=grass_density)
        height_targets = set_fuel_parameter(parameter="height", grass=grass_height)
        calibrated_duet = calibrate(
            duet_run, fuel_parameter_targets=[density_targets, height_targets]
        )
        # now do moisture
        grass_moisture = assign_targets(method="maxmin", max=0.5, min=0.05)
        moisture_targets = set_fuel_parameter(
            parameter="moisture", grass=grass_moisture
        )
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=moisture_targets)

        # try with two fueltypes
        litter_density = assign_targets(method="maxmin", max=0.1, min=0.01)
        density_targets = set_fuel_parameter(
            parameter="density", litter=litter_density, grass=grass_density
        )
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=density_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.density, np.ndarray)
        assert np.allclose(calibrated_duet.height, duet_run.height)
        assert np.allclose(calibrated_duet.density, duet_run.density) == False
        assert np.max(calibrated_duet.density[0, :, :]) == 1.0
        assert np.max(calibrated_duet.density[1, :, :]) == 0.1

    def test_calibrate_meansd(self):
        # try 1 fueltype and 1 parameter
        duet_run = import_duet(DATA_DIR / "v2", 295, 333, 9)
        grass_density = assign_targets(method="meansd", mean=0.6, sd=0.3)
        density_targets = set_fuel_parameter(parameter="density", grass=grass_density)
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=density_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.density, np.ndarray)
        assert np.allclose(calibrated_duet.height, duet_run.height)
        assert np.allclose(calibrated_duet.density, duet_run.density) == False
        assert np.allclose(calibrated_duet.density[1, :, :], duet_run.density[1, :, :])
        assert math.isclose(
            np.mean(calibrated_duet.density[0, :, :]), np.float32(0.6), abs_tol=10**-6
        )
        assert math.isclose(
            np.std(calibrated_duet.density[0, :, :]), np.float32(0.3), abs_tol=10**-6
        )
        # try density and height
        grass_height = assign_targets(method="meansd", mean=0.5, sd=0.05)
        density_targets = set_fuel_parameter(parameter="density", grass=grass_density)
        height_targets = set_fuel_parameter(parameter="height", grass=grass_height)
        calibrated_duet = calibrate(
            duet_run, fuel_parameter_targets=[density_targets, height_targets]
        )

        # try with two fueltypes
        litter_density = assign_targets(method="meansd", mean=0.05, sd=0.005)
        density_targets = set_fuel_parameter(
            parameter="density", litter=litter_density, grass=grass_density
        )
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=density_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.density, np.ndarray)
        assert np.allclose(calibrated_duet.height, duet_run.height)
        assert np.allclose(calibrated_duet.density, duet_run.density) == False

    def test_constant_calibration(self):
        duet_run = import_duet(DATA_DIR / "v2", 295, 333, 9)
        grass_height = assign_targets(method="constant", value=0.5)
        litter_height = assign_targets(method="constant", value=0.05)
        height_targets = set_fuel_parameter(
            parameter="height", grass=grass_height, litter=litter_height
        )
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=height_targets)
        assert isinstance(calibrated_duet.height, np.ndarray)
        assert math.isclose(
            np.mean(
                calibrated_duet.height[0, :, :][calibrated_duet.height[0, :, :] > 0]
            ),
            np.float32(0.5),
            abs_tol=10**-6,
        )
        assert math.isclose(
            np.mean(
                calibrated_duet.height[1, :, :][calibrated_duet.height[1, :, :] > 0]
            ),
            np.float32(0.05),
            abs_tol=10**-6,
        )

    def test_fueltype_all(self):
        duet_run = import_duet(DATA_DIR / "v2", 295, 333, 9)
        density = assign_targets(method="maxmin", max=2.0, min=0.5)
        grass_height = assign_targets(method="constant", value=0.75)
        litter_height = assign_targets(method="constant", value=0.15)
        density_targets = set_fuel_parameter(parameter="density", all=density)
        height_targets = set_fuel_parameter(
            parameter="height", grass=grass_height, litter=litter_height
        )
        calibrated_duet = calibrate(
            duet_run,
            fuel_parameter_targets=[density_targets, height_targets],
        )
        assert isinstance(calibrated_duet, DuetRun)
        # first assert that the height values are what we set them to be
        assert math.isclose(
            np.mean(
                calibrated_duet.height[0, :, :][calibrated_duet.height[0, :, :] > 0]
            ),
            np.float32(0.75),
            abs_tol=10**-6,
        )
        assert math.isclose(
            np.mean(
                calibrated_duet.height[1, :, :][calibrated_duet.height[1, :, :] > 0]
            ),
            np.float32(0.15),
            abs_tol=10**-6,
        )
        # then assert that a) the density array has two layers, and b) the overall max and min matches
        assert calibrated_duet.density.shape[0] == 2
        assert np.max(np.sum(calibrated_duet.density, axis=0)) == np.float32(2.0)
        assert np.min(np.sum(calibrated_duet.density, axis=0)) == np.float32(0.5)


def plot_array(x, title):
    plt.figure(2)
    plt.set_cmap("viridis")
    plt.imshow(x, origin="lower")
    plt.colorbar()
    plt.title(title, fontsize=18)
    plt.show()
