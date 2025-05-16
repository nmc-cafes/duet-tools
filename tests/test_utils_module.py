"""
Test module for the utils module of the duet_tools package.
"""

from __future__ import annotations

import pytest
import geojson
import numpy as np
from pathlib import Path
from scipy.io import FortranFile

from duet_tools.utils import (
    read_dat_to_array,
    write_array_to_dat,
    read_shapefile_to_geojson,
)

TEST_DIR = Path(__file__).parent
TMP_DIR = TEST_DIR / "tmp"
TMP_DIR.mkdir(exist_ok=True)
DATA_DIR = TEST_DIR / "test-data"


def get_test_array(dir: Path, dat: str, shape: tuple) -> np.ndarray:
    with open(Path(dir, dat), "rb") as fin:
        array = FortranFile(fin).read_reals(dtype=np.float32).reshape(shape, order="F")
    return array


class TestUtilFunctions:
    def test_read_dat(self):
        dir = DATA_DIR / "v2"
        dat = "treesrhof.dat"
        arr = read_dat_to_array(dir, dat, nx=333, ny=295, nz=84)
        assert isinstance(arr, np.ndarray)
        assert arr.shape == (84, 295, 333)
        test = get_test_array(dir, dat, (333, 295, 84))
        assert arr.all() == test.T.all()
        # try with invalid shape
        with pytest.raises(ValueError):
            arr = read_dat_to_array(dir, dat, nx=333, ny=294, nz=84)
        # try with invalid order
        with pytest.raises(ValueError):
            arr = read_dat_to_array(dir, dat, nx=333, ny=295, nz=84, order="P")
        # try with nsp instead of nz
        dat = "surface_rhof_layered.dat"
        arr = read_dat_to_array(dir, dat, nx=333, ny=295, nsp=9)
        assert isinstance(arr, np.ndarray)
        assert arr.shape == (9, 295, 333)
        test = get_test_array(dir, dat, (9, 333, 295))
        assert arr.all() == np.moveaxis(test, 1, 2).all()

    def test_write_dat(self):
        dir_in = DATA_DIR / "v2"
        dir_out = TMP_DIR
        dat_in = "treesrhof.dat"
        dat_out = "temp.dat"
        arr1 = get_test_array(dir_in, dat_in, (333, 295, 84)).T
        write_array_to_dat(arr1, dat_out, dir_out)
        arr2 = get_test_array(dir_out, dat_out, (333, 295, 84)).T
        assert arr1.all() == arr2.all()

    def test_shapefile_to_geojson(self):
        # test with .shp
        json = read_shapefile_to_geojson(DATA_DIR / "new_mixedcon_area.shp")
        assert isinstance(json, geojson.Polygon)
        # test with .zip
        json = read_shapefile_to_geojson(DATA_DIR / "new_mixedcon_area.zip")
        assert isinstance(json, geojson.Polygon)
