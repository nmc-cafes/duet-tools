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
    import_duet_manual,
    assign_targets,
    set_fuel_parameter,
    set_loading,
    set_moisture,
    set_depth,
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
        duet_run = import_duet_manual(
            directory=DATA_DIR / "v1",
            loading_grid_name="surface_rhof.dat",
            moisture_grid_name="surface_moist.dat",
            depth_grid_name="surface_depth.dat",
            nx=976,
            ny=1998,
            nsp=2,
            version="v1",
        )
        # test that data types are correct
        assert isinstance(duet_run, DuetRun)
        assert isinstance(duet_run.loading, np.ndarray)
        assert isinstance(duet_run.moisture, np.ndarray)
        assert isinstance(duet_run.depth, np.ndarray)
        # test array shapes
        assert duet_run.loading.shape == (2, 1998, 976)
        assert duet_run.moisture.shape == (2, 1998, 976)
        assert duet_run.depth.shape == (2, 1998, 976)
        # test that wrong dimensions raises error
        with pytest.raises(ValueError):
            duet_run = import_duet_manual(
                directory=DATA_DIR / "v1",
                loading_grid_name="surface_rhof.dat",
                moisture_grid_name="surface_moist.dat",
                depth_grid_name="surface_depth.dat",
                nx=976,
                ny=1998,
                nsp=3,
                version="v1",
            )
        # test that wrong version number raises error
        with pytest.raises(ValueError):
            duet_run = import_duet_manual(
                directory=DATA_DIR / "v1",
                loading_grid_name="surface_rhof.dat",
                moisture_grid_name="surface_moist.dat",
                depth_grid_name="surface_depth.dat",
                nx=976,
                ny=1998,
                nsp=2,
                version="v3",
            )

    def test_import_duet_v2(self):
        duet_run = import_duet(directory=DATA_DIR / "v2")
        assert isinstance(duet_run, DuetRun)
        # test that data types are correct
        assert isinstance(duet_run.loading, np.ndarray)
        assert isinstance(duet_run.moisture, np.ndarray)
        assert isinstance(duet_run.depth, np.ndarray)
        # test array shapes
        assert duet_run.loading.shape == (3, 295, 333)
        assert duet_run.moisture.shape == (3, 295, 333)
        assert duet_run.depth.shape == (3, 295, 333)

    def test_import_duet_v2_manual(self):
        duet_run = import_duet_manual(
            directory=DATA_DIR / "v2",
            loading_grid_name="surface_rhof_layered.dat",
            moisture_grid_name="surface_moist_layered.dat",
            depth_grid_name="surface_depth_layered.dat",
            nx=333,
            ny=295,
            nsp=9,
            version="v2",
        )
        assert isinstance(duet_run, DuetRun)
        # test that data types are correct
        assert isinstance(duet_run.loading, np.ndarray)
        assert isinstance(duet_run.moisture, np.ndarray)
        assert isinstance(duet_run.depth, np.ndarray)
        # test array shapes
        assert duet_run.loading.shape == (3, 295, 333)
        assert duet_run.moisture.shape == (3, 295, 333)
        assert duet_run.depth.shape == (3, 295, 333)
        # test that wrong dimensions raise error
        with pytest.raises(ValueError):
            duet_run = import_duet_manual(
                directory=DATA_DIR / "v2",
                loading_grid_name="surface_rhof.dat",
                moisture_grid_name="surface_moist.dat",
                depth_grid_name="surface_depth.dat",
                nx=333,
                ny=295,
                nsp=3,
                version="v2",
            )

    def test_to_numpy(self):
        duet_run = import_duet(directory=DATA_DIR / "v2")
        # test each fuel parameter and type
        grass_loading = duet_run.to_numpy("grass", "loading")
        coniferous_loading = duet_run.to_numpy("coniferous", "loading")
        deciduous_loading = duet_run.to_numpy("deciduous", "loading")
        litter_loading = duet_run.to_numpy("litter", "loading")
        grass_moisture = duet_run.to_numpy("grass", "moisture")
        coniferous_moisture = duet_run.to_numpy("coniferous", "moisture")
        deciduous_moisture = duet_run.to_numpy("deciduous", "moisture")
        litter_moisture = duet_run.to_numpy("litter", "moisture")
        grass_depth = duet_run.to_numpy("grass", "depth")
        coniferous_depth = duet_run.to_numpy("coniferous", "depth")
        deciduous_depth = duet_run.to_numpy("deciduous", "depth")
        litter_depth = duet_run.to_numpy("litter", "depth")
        assert np.array_equal(grass_loading, duet_run.loading[0, :, :])
        assert np.array_equal(coniferous_loading, duet_run.loading[1, :, :])
        assert np.array_equal(deciduous_loading, duet_run.loading[2, :, :])
        assert np.array_equal(
            litter_loading, np.sum(duet_run.loading[1:, :, :], axis=0)
        )
        assert np.array_equal(grass_moisture, duet_run.moisture[0, :, :])
        assert np.array_equal(coniferous_moisture, duet_run.moisture[1, :, :])
        assert np.array_equal(deciduous_moisture, duet_run.moisture[2, :, :])
        litter_weights = _maxmin_calibration(duet_run.loading[1:, :, :], max=1.0, min=0)
        litter_weights[litter_weights == 0] = 0.01
        litter_masked = np.ma.masked_array(
            duet_run.moisture[1:, :, :], duet_run.moisture[1:, :, :] == 0
        )
        litter_averaged = np.ma.average(litter_masked, axis=0, weights=litter_weights)
        weighted_average_litter_moisture = np.ma.filled(litter_averaged, 0)
        assert np.array_equal(litter_moisture, weighted_average_litter_moisture)
        assert np.array_equal(grass_depth, duet_run.depth[0, :, :])
        assert np.array_equal(coniferous_depth, duet_run.depth[1, :, :])
        assert np.array_equal(deciduous_depth, duet_run.depth[2, :, :])
        assert np.array_equal(litter_depth, np.sum(duet_run.depth[1:, :, :], axis=0))
        # test integrated and separated
        separated_loading = duet_run.to_numpy("separated", "loading")
        assert np.array_equal(separated_loading, duet_run.loading)
        integrated_loading = duet_run.to_numpy("integrated", "loading")
        integrated_depth = duet_run.to_numpy("integrated", "depth")
        integrated_moisture = duet_run.to_numpy("integrated", "moisture")
        assert np.array_equal(integrated_loading, np.sum(duet_run.loading, axis=0))
        assert np.array_equal(integrated_depth, np.max(duet_run.depth, axis=0))
        weights = _maxmin_calibration(duet_run.loading.copy(), max=1.0, min=0)
        weights[weights == 0] = 0.01
        masked = np.ma.masked_array(duet_run.moisture, duet_run.moisture == 0)
        averaged = np.ma.average(masked, axis=0, weights=weights)
        weighted_average_moisture = np.ma.filled(averaged, 0)
        assert np.array_equal(integrated_moisture, weighted_average_moisture)

    def test_to_quicfire(self):
        duet_run = import_duet(directory=DATA_DIR / "v2")
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
            treesrhof, duet_run._integrate_all("loading").astype("float32")
        )
        assert np.array_equal(
            treesfueldepth, duet_run._integrate_all("depth").astype("float32")
        )
        assert np.array_equal(
            treesmoist, duet_run._integrate_all("moisture").astype("float32")
        )

        files = ["treesrhof.dat", "treesmoist.dat", "treesfueldepth.dat"]
        for file in files:
            path = TMP_DIR / file
            path.unlink()
        duet_run.to_quicfire(TMP_DIR, loading=False, moisture=False, depth=False)
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
        loading_targets = set_loading(grass=maxmin_targets, litter=meansd_targets)
        moisture_targets = set_moisture(grass=constant_target)
        depth_targets = set_depth(all=maxmin_targets)
        assert isinstance(loading_targets, FuelParameter)
        assert isinstance(moisture_targets, FuelParameter)
        assert isinstance(depth_targets, FuelParameter)
        assert loading_targets.parameter == "loading"
        assert moisture_targets.parameter == "moisture"
        # test generic function
        loading_targets = set_fuel_parameter(
            parameter="loading", grass=maxmin_targets, litter=meansd_targets
        )
        assert isinstance(loading_targets, FuelParameter)
        assert loading_targets.parameter == "loading"
        assert loading_targets.fuel_types == ["grass", "litter"]
        assert loading_targets.targets == [maxmin_targets, meansd_targets]
        # test validation
        with pytest.raises(ValueError):
            loading_targets = set_fuel_parameter(
                parameter="bulk loading", litter=meansd_targets
            )
        with pytest.raises(ValueError):
            loading_targets = set_fuel_parameter(
                parameter="loading", both=meansd_targets
            )
        with pytest.raises(ValueError):
            loading_targets = set_fuel_parameter(
                parameter="loading", litter=meansd_targets, all=constant_target
            )


class TestCalibrate:
    def test_calibrate_maxmin(self):
        # try 1 fueltype and 1 parameter
        duet_run = import_duet(DATA_DIR / "v2")
        grass_loading = assign_targets(method="maxmin", max=1.0, min=0.2)
        loading_targets = set_fuel_parameter(parameter="loading", grass=grass_loading)
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=loading_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.loading, np.ndarray)
        assert np.allclose(calibrated_duet.depth, duet_run.depth)
        assert np.allclose(calibrated_duet.loading, duet_run.loading) == False
        assert np.allclose(calibrated_duet.loading[1, :, :], duet_run.loading[1, :, :])
        assert np.allclose(calibrated_duet.loading[2, :, :], duet_run.loading[2, :, :])
        assert np.max(calibrated_duet.loading[0, :, :]) == 1.0
        assert np.min(calibrated_duet.loading[0, :, :]) == 0.2
        # try loading and depth
        grass_depth = assign_targets(method="maxmin", max=1.0, min=0.2)
        loading_targets = set_fuel_parameter(parameter="loading", grass=grass_loading)
        depth_targets = set_fuel_parameter(parameter="depth", grass=grass_depth)
        calibrated_duet = calibrate(
            duet_run, fuel_parameter_targets=[loading_targets, depth_targets]
        )
        # now do moisture
        grass_moisture = assign_targets(method="maxmin", max=0.5, min=0.05)
        moisture_targets = set_fuel_parameter(
            parameter="moisture", grass=grass_moisture
        )
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=moisture_targets)

        # try with two fueltypes
        litter_loading = assign_targets(method="maxmin", max=0.1, min=0.01)
        loading_targets = set_fuel_parameter(
            parameter="loading", litter=litter_loading, grass=grass_loading
        )
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=loading_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.loading, np.ndarray)
        assert np.allclose(calibrated_duet.depth, duet_run.depth)
        assert np.allclose(calibrated_duet.loading, duet_run.loading) == False
        assert np.max(calibrated_duet.loading[0, :, :]) == 1.0
        assert np.max(calibrated_duet.loading[1:, :, :]) == 0.1

        # try with three fueltypes
        coniferous_loading = assign_targets(method="maxmin", max=0.1, min=0.01)
        deciduous_loading = assign_targets(method="maxmin", max=0.5, min=0.2)
        loading_targets = set_fuel_parameter(
            parameter="loading",
            coniferous=coniferous_loading,
            deciduous=deciduous_loading,
            grass=grass_loading,
        )
        # there is no deciduous litter, so it will throw an error
        with pytest.raises(ValueError):
            calibrated_duet = calibrate(
                duet_run, fuel_parameter_targets=loading_targets
            )
        # now without deciduous litter targets
        loading_targets = set_fuel_parameter(
            parameter="loading",
            coniferous=coniferous_loading,
            grass=grass_loading,
        )
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=loading_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.loading, np.ndarray)
        assert np.allclose(calibrated_duet.depth, duet_run.depth)
        assert np.allclose(calibrated_duet.loading, duet_run.loading) == False
        assert np.max(calibrated_duet.loading[0, :, :]) == 1.0
        assert np.max(calibrated_duet.loading[1, :, :]) == 0.1
        assert np.max(calibrated_duet.loading[2, :, :]) == 0

    def test_calibrate_meansd(self):
        # try 1 fueltype and 1 parameter
        duet_run = import_duet(DATA_DIR / "v2")
        grass_loading = assign_targets(method="meansd", mean=0.6, sd=0.3)
        loading_targets = set_fuel_parameter(parameter="loading", grass=grass_loading)
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=loading_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.loading, np.ndarray)
        assert np.allclose(calibrated_duet.depth, duet_run.depth)
        assert np.allclose(calibrated_duet.loading, duet_run.loading) == False
        assert np.allclose(
            calibrated_duet.loading[1:, :, :], duet_run.loading[1:, :, :]
        )
        assert math.isclose(
            np.mean(calibrated_duet.loading[0, :, :]), np.float32(0.6), abs_tol=10**-6
        )
        assert math.isclose(
            np.std(calibrated_duet.loading[0, :, :]), np.float32(0.3), abs_tol=10**-6
        )
        # try loading and depth
        grass_depth = assign_targets(method="meansd", mean=0.5, sd=0.05)
        loading_targets = set_fuel_parameter(parameter="loading", grass=grass_loading)
        depth_targets = set_fuel_parameter(parameter="depth", grass=grass_depth)
        calibrated_duet = calibrate(
            duet_run, fuel_parameter_targets=[loading_targets, depth_targets]
        )

        # try with two fueltypes
        litter_loading = assign_targets(method="meansd", mean=0.05, sd=0.005)
        loading_targets = set_fuel_parameter(
            parameter="loading", litter=litter_loading, grass=grass_loading
        )
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=loading_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.loading, np.ndarray)
        assert np.allclose(calibrated_duet.depth, duet_run.depth)
        assert np.allclose(calibrated_duet.loading, duet_run.loading) == False

        # try with three fueltypes
        coniferous_loading = assign_targets(method="meansd", mean=0.1, sd=0.03)
        deciduous_loading = assign_targets(method="meansd", mean=0.5, sd=0.2)
        loading_targets = set_fuel_parameter(
            parameter="loading",
            coniferous=coniferous_loading,
            deciduous=deciduous_loading,
            grass=grass_loading,
        )
        # there is no deciduous litter, so it will throw an error
        with pytest.raises(ValueError):
            calibrated_duet = calibrate(
                duet_run, fuel_parameter_targets=loading_targets
            )
        # now without deciduous litter targets
        loading_targets = set_fuel_parameter(
            parameter="loading",
            coniferous=coniferous_loading,
            grass=grass_loading,
        )
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=loading_targets)
        assert isinstance(calibrated_duet, DuetRun)
        assert isinstance(calibrated_duet.loading, np.ndarray)
        assert np.allclose(calibrated_duet.depth, duet_run.depth)
        assert np.allclose(calibrated_duet.loading, duet_run.loading) == False
        assert math.isclose(
            np.mean(calibrated_duet.loading[0, :, :]), np.float32(0.6), abs_tol=10**-6
        )
        assert math.isclose(
            np.std(calibrated_duet.loading[0, :, :]), np.float32(0.3), abs_tol=10**-6
        )

    def test_constant_calibration(self):
        duet_run = import_duet(DATA_DIR / "v2")
        grass_depth = assign_targets(method="constant", value=0.5)
        litter_depth = assign_targets(method="constant", value=0.05)
        depth_targets = set_fuel_parameter(
            parameter="depth", grass=grass_depth, litter=litter_depth
        )
        calibrated_duet = calibrate(duet_run, fuel_parameter_targets=depth_targets)
        assert isinstance(calibrated_duet.depth, np.ndarray)
        assert math.isclose(
            np.mean(calibrated_duet.depth[0, :, :][calibrated_duet.depth[0, :, :] > 0]),
            np.float32(0.5),
            abs_tol=10**-6,
        )
        assert math.isclose(
            np.mean(
                calibrated_duet.depth[1:, :, :][calibrated_duet.depth[1:, :, :] > 0]
            ),
            np.float32(0.05),
            abs_tol=10**-6,
        )

    def test_fueltype_all(self):
        duet_run = import_duet(DATA_DIR / "v2")
        loading = assign_targets(method="maxmin", max=2.0, min=0.5)
        grass_depth = assign_targets(method="constant", value=0.75)
        litter_depth = assign_targets(method="constant", value=0.15)
        loading_targets = set_fuel_parameter(parameter="loading", all=loading)
        depth_targets = set_fuel_parameter(
            parameter="depth", grass=grass_depth, litter=litter_depth
        )
        calibrated_duet = calibrate(
            duet_run,
            fuel_parameter_targets=[loading_targets, depth_targets],
        )
        assert isinstance(calibrated_duet, DuetRun)
        # first assert that the depth values are what we set them to be
        assert math.isclose(
            np.mean(calibrated_duet.depth[0, :, :][calibrated_duet.depth[0, :, :] > 0]),
            np.float32(0.75),
            abs_tol=10**-6,
        )
        assert math.isclose(
            np.mean(calibrated_duet.depth[1, :, :][calibrated_duet.depth[1, :, :] > 0]),
            np.float32(0.15),
            abs_tol=10**-6,
        )
        # then assert that a) the loading array has three layers, and b) the overall max and min matches
        assert calibrated_duet.loading.shape[0] == 3
        assert np.max(np.sum(calibrated_duet.loading, axis=0)) == np.float32(2.0)
        assert np.min(np.sum(calibrated_duet.loading, axis=0)) == np.float32(0.5)


def plot_array(x, title):
    plt.figure(2)
    plt.set_cmap("viridis")
    plt.imshow(x, origin="lower")
    plt.colorbar()
    plt.title(title, fontsize=18)
    plt.show()
