"""
Test module for the landfire module of the duet_tools package.
"""

from __future__ import annotations

import pytest
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import geojson

from duet_tools.calibration import Targets
from duet_tools.landfire import LandfireQuery, query_landfire, assign_targets_from_sb40

TEST_DIR = Path(__file__).parent
TMP_DIR = TEST_DIR / "tmp"
TMP_DIR.mkdir(exist_ok=True)
DATA_DIR = TEST_DIR / "test-data"


class TestLandfireTargets:
    @classmethod
    def get_geojson(self):
        with open(DATA_DIR / "anderson_butte.geojson") as fid:
            sample_geojson = geojson.load(fid)
        return geojson.Polygon(sample_geojson["features"][0]["geometry"]["coordinates"])

    def test_query_landfire(self):
        sample_aoi = self.get_geojson()
        query = query_landfire(
            area_of_interest=sample_aoi, directory=TMP_DIR, input_epsg=4326
        )
        assert isinstance(query, LandfireQuery)
        assert isinstance(query.fuel_types, np.ndarray)
        assert isinstance(query.density, np.ndarray)
        assert isinstance(query.moisture, np.ndarray)
        assert isinstance(query.height, np.ndarray)

    def test_assign_targets_from_sb40(self):
        sample_aoi = self.get_geojson()
        query = query_landfire(
            area_of_interest=sample_aoi, directory=TMP_DIR, input_epsg=4326
        )
        # test just grass density
        grass_density = assign_targets_from_sb40(query, "grass", "density")
        assert isinstance(grass_density, Targets)
        assert grass_density.method == "maxmin"
        assert grass_density.args == ["max", "min"]
        assert len(grass_density.targets) == 2
        # test fuel and parameter with only one value
        with pytest.warns(UserWarning):
            grass_moisture = assign_targets_from_sb40(query, "grass", "moisture")
        grass_moisture = assign_targets_from_sb40(
            query, "grass", "moisture", method="constant"
        )
        # get the rest of the fuels and params with maxmin
        grass_height = assign_targets_from_sb40(query, "grass", "height")
        litter_density = assign_targets_from_sb40(query, "litter", "density")
        litter_moisture = assign_targets_from_sb40(query, "litter", "density")
        litter_height = assign_targets_from_sb40(query, "litter", "density")
        all_density = assign_targets_from_sb40(query, "all", "density")
        all_moisture = assign_targets_from_sb40(query, "all", "moisture")
        all_height = assign_targets_from_sb40(query, "all", "height")
        # test a couple with meansd
        grass_density_meansd = assign_targets_from_sb40(
            query, "grass", "density", method="meansd"
        )
        all_height = assign_targets_from_sb40(query, "all", "height", "meansd")
        # test wrong inputs
        with pytest.raises(ValueError):
            grass_height = assign_targets_from_sb40(query, "both", "height")
        with pytest.raises(ValueError):
            grass_height = assign_targets_from_sb40(query, "all", "moist")
        with pytest.raises(ValueError):
            grass_height = assign_targets_from_sb40(
                query, "all", "height", method="minmax"
            )


def plot_array(x, title):
    plt.figure(2)
    plt.set_cmap("viridis")
    plt.imshow(x, origin="lower")
    plt.colorbar()
    plt.title(title, fontsize=18)
    plt.show()
