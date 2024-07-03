"""
Utility functions for DUET tools modules
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from scipy.io import FortranFile


def read_dat_to_array(
    directory: str | Path, filename: str, nx: int, ny: int, nz: int, order: str = "C"
) -> np.ndarray:
    """
    Reads a fortran binary file (.dat) to a numpy array

    Parameters
    ----------
    directory: Path | str
        Path to directory of the .dat file.
    filename: str
        Name of the .dat file
    nx : int
        Number of cells in the x-direction
    ny : int
        Number of cells in the y-direction
    nz : int
        Number of cells in the z-direction
    order : str
        Order of the .dat file. Must be one of "C" or "F". Defaults to "C".

    Returns
    -------
        A numpy array with shape (nz, ny, nx).
    """
    if order not in ["C", "F"]:
        raise ValueError('Order must be either "C" or "F".')
    if isinstance(directory, str):
        directory = Path(directory)
    with open(Path(directory, filename), "rb") as fin:
        array = (
            FortranFile(fin)
            .read_reals(dtype="float32")
            .reshape((ny, nx, nz), order=order)
        )
    return np.moveaxis(array, 2, 0)


def write_array_to_dat(
    array: np.ndarray,
    dat_name: str,
    output_dir: Path | str,
    dtype: type = np.float32,
    reshape: bool = True,
) -> None:
    """
    Write a numpy array to a fortran binary file (.dat).

    Parameters
    ----------
    array : np.ndarray
        numpy array to be written to a file
    dat_name : str
        Filename ending with .dat
    output_dir : Path | str
        Directory where file will be written
    dtype : type
        Data type of the array. Defaults to np.float32
    reshape: bool
        Whether to reshape the array to (y,x,z). Defaults to True.
    """
    if isinstance(output_dir, str):
        output_dir = Path(output_dir)
    # Reshape array from (y, x, z) to (z, y, x) (also for fortran)
    if reshape:
        if len(array.shape) == 3:
            array = np.moveaxis(array, 0, 2).astype(dtype)
        else:
            array = array.astype(dtype)
    else:
        array = array.astype(dtype)

    # Write the zarr array to a dat file with scipy FortranFile package
    with FortranFile(Path(output_dir, dat_name), "w") as f:
        f.write_record(array)
