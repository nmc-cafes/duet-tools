"""
Test module for the inputs module of the duet_tools package.
"""

from __future__ import annotations

import pytest
from pathlib import Path

from duet_tools.inputs import InputFile

TEST_DIR = Path(__file__).parent
TMP_DIR = TEST_DIR / "tmp"
TMP_DIR.mkdir(exist_ok=True)
DATA_DIR = TEST_DIR / "test-data"


class TestInputFile:
    def test_create_input_file(self):
        input_file = InputFile.create(
            nx=234, ny=345, nz=30, duration=5, wind_direction=225
        )
        # test that data types are correct
        assert isinstance(input_file, InputFile)
        assert isinstance(input_file.nx, int)
        assert isinstance(input_file.dx, float)
        assert isinstance(input_file.duration, int)
        assert isinstance(input_file.random_seed, int)

        # test that input values are correct
        assert input_file.ny == 345
        assert input_file.dy == 2.0
        assert input_file.wind_variability == 359

    def test_read_input_file(self):
        input_file = InputFile.from_directory(DATA_DIR)

        assert isinstance(input_file, InputFile)
        assert input_file.nx == 234
        assert input_file.ny == 345
        assert input_file.nz == 30
        assert input_file.duration == 5
        assert input_file.random_seed == 47
        assert input_file.dx == 2.0
        assert input_file.dy == 2.0
        assert input_file.dz == 1.0
        assert input_file.wind_direction == 225
        assert input_file.wind_variability == 359.0

    def test_modify_input_file(self):
        verification = [
            "234  ! number of cells in x direction\n",
            "345  ! number of cells in y direction\n",
            "30  ! number of cells in z direction\n",
            "2.0  ! cell size in x direction (in meters)\n",
            "2.0  ! cell size in y direction (in meters)\n",
            "1.0  ! cell size in z direction (in meters)\n",
            "47  ! random number seed\n",
            "225  ! wind direction (in degrees)\n",
            "359.0  ! wind direction variability (in degrees)\n",
            "5  ! duration of simulation (in years)\n",
        ]
        input_file = InputFile.create(
            nx=235, ny=344, nz=30, duration=5, wind_direction=225, random_seed=47
        )

        input_file.nx = 234
        input_file.ny = 345
        input_file.to_file(TMP_DIR)

        with open(TMP_DIR / "duet.in", "r") as f:
            lines = f.readlines()

        for i in range(len(lines)):
            assert lines[i] == verification[i]

    def test_write_input_file(self):
        verification = [
            "234  ! number of cells in x direction\n",
            "345  ! number of cells in y direction\n",
            "30  ! number of cells in z direction\n",
            "2.0  ! cell size in x direction (in meters)\n",
            "2.0  ! cell size in y direction (in meters)\n",
            "1.0  ! cell size in z direction (in meters)\n",
            "47  ! random number seed\n",
            "225  ! wind direction (in degrees)\n",
            "359.0  ! wind direction variability (in degrees)\n",
            "5  ! duration of simulation (in years)\n",
        ]

        input_file = InputFile.create(
            nx=234, ny=345, nz=30, duration=5, wind_direction=225, random_seed=47
        )
        input_file.to_file(TMP_DIR)

        with open(TMP_DIR / "duet.in", "r") as f:
            lines = f.readlines()

        for i in range(len(lines)):
            assert lines[i] == verification[i]

        # test validation
        for invalid_int in ["234", 234.5, -234]:
            input_file.nx = invalid_int
            with pytest.raises(ValueError):
                input_file.to_file(TMP_DIR)
        for invalid_float in ["2", -2]:
            input_file.dx = invalid_float
            with pytest.raises(ValueError):
                input_file.to_file(TMP_DIR)
        for invalid_deg in [-15, "270", 360, 361]:
            input_file.wind_direction = invalid_deg
            with pytest.raises(ValueError):
                input_file.to_file(TMP_DIR)
