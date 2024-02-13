"""
Utility functions for DUET tools modules
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from scipy.io import FortranFile


def read_dat_to_array(
    directory: str | Path, filename: str, nx: int, ny: int, nz: int, order: str = "F"
) -> np.ndarray:
    if isinstance(directory, str):
        directory = Path(directory)
    with open(Path(directory, filename), "rb") as fin:
        array = (
            FortranFile(fin)
            .read_reals(dtype="float32")
            .reshape((nz, ny, nx), order=order)
        )
    return array


def write_array_to_dat(
    array: np.ndarray,
    dat_name: str,
    output_dir: Path,
    dtype: type = np.float32,
    reshape: bool = True,
) -> None:
    """
    Write a numpy array to a fortran binary file. Array must be cast to the
    appropriate data type before calling this function. If the array is 3D,
    the array will be reshaped from (y, x, z) to (z, y, x) for fortran.
    """
    # Reshape array from (y, x, z) to (z, y, x) (also for fortran)
    if reshape:
        if len(array.shape) == 3:
            array = np.moveaxis(array, 2, 0).astype(dtype)
        else:
            array = array.astype(dtype)
    else:
        array = array.astype(dtype)

    # Write the zarr array to a dat file with scipy FortranFile package
    with FortranFile(Path(output_dir, dat_name), "w") as f:
        f.write_record(array)
