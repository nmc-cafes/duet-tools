"""
DUET inputs module
"""

from __future__ import annotations

# Core Imports
from pathlib import Path
from random import randint


class InputFile:
    """
    Class representing a DUET input file.

    Attributes
    ----------
    nx: int
        Number of cells in the x-direction.
    ny: int
        Number of cells in the y-direction.
    nz: int
        Number of cells in the z-direction.
    dx: float
        Cell size in the x-direction (m)
    dy: float
        Cell size in the y-direction (m)
    dz: float
        Cell size in the z-direction (m)
    random_seed: int
        Random number seed.
    wind_direction: float
        Wind direction (degrees).
    wind_variability: float
        Wind direction variability (degrees).
    duration:
        Duration of simulation (years).
    """

    def __init__(
        self,
        nx: int,
        ny: int,
        nz: int,
        dx: float,
        dy: float,
        dz: float,
        random_seed: int,
        wind_direction: float,
        wind_variability: float,
        duration: int,
    ):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.random_seed = random_seed
        self.wind_direction = wind_direction
        self.wind_variability = wind_variability
        self.duration = duration

    @classmethod
    def create(
        cls,
        nx: int,
        ny: int,
        nz: int,
        duration: int,
        wind_direction: float,
        dx: float = 2.0,
        dy: float = 2.0,
        dz: float = 1.0,
        random_seed: int = None,
        wind_variability: float = 359.0,
    ):
        if random_seed is None:
            random_seed = randint(1, 50000)
        return cls(
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            random_seed,
            wind_direction,
            wind_variability,
            duration,
        )

    def to_file(self, directory: Path | str):
        """
        Writes a DUET input file to the specified path

        Parameters
        ----------
        directory: Path | str
            Directory for writing DUET input file

        Returns
        -------
        None:
            Writes duet.in to directory
        """
        if isinstance(directory, str):
            directory = Path(directory)

        out_path = directory / "duet.in"
        with open(out_path, "w") as f:
            f.write(f"{self.nx}  ! number of cells in x direction\n")
            f.write(f"{self.ny}  ! number of cells in y direction\n")
            f.write(f"{self.nz}  ! number of cells in z direction\n")
            f.write(f"{self.dx}  ! cell size in x direction (in meters)\n")
            f.write(f"{self.dy}  ! cell size in y direction (in meters)\n")
            f.write(f"{self.dz}  ! cell size in z direction (in meters)\n")
            f.write(f"{self.random_seed}  ! random number seed\n")
            f.write(f"{self.wind_direction}  ! wind direction (in degrees)\n")
            f.write(
                f"{self.wind_variability}  ! wind direction variability (in degrees)\n"
            )
            f.write(f"{self.duration}  ! duration of simulation (in years)\n")
