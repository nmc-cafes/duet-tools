"""
Utility functions for DUET tools modules
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from scipy.io import FortranFile


def read_dat_to_array(
    directory: str | Path,
    filename: str,
    nx: int,
    ny: int,
    nz: int = None,
    nsp: int = None,
    order: str = "F",
    dtype: type = np.float32,
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
    nsp: int
        Number of species
    order : str
        Order of the .dat file. Must be one of "C" or "F". Defaults to "F".

    Returns
    -------
        A numpy array with shape (nz, ny, nx).
    """
    if order not in ["C", "F"]:
        raise ValueError('Order must be either "C" or "F".')
    if isinstance(directory, str):
        directory = Path(directory)

    if (nz is None) and (nsp is None):
        shape = (nx, ny)
    elif nz is None:
        shape = (nsp, nx, ny)
    elif nsp is None:
        shape = (nx, ny, nz)
    else:
        shape = (nsp, nx, ny, nz)

    with open(Path(directory, filename), "rb") as fin:
        array = FortranFile(fin).read_reals(dtype=dtype).reshape(shape, order=order)

    if (nz is None) and (nsp is None):
        return np.moveaxis(array, 1, 0)
    elif nz is None:
        return np.moveaxis(array, 2, 1)
    elif nsp is None:
        return np.transpose(array)
    else:
        return np.moveaxis(np.moveaxis(array, 3, 1), 3, 1)


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
        Whether to reshape the array. Array dimensions in duet-tools are either (nz,ny,nx)
        or (nsp,ny,nx) and will be written in row-major order, meaning that with column-major
        order (default for fortran), they will be (nx, ny, nz). Reshaping (nsp, ny, nx) arrays
        to (ny, nx, nsp) will result in the column-major order of (nsp, nx, ny), which is
        expected by DUET and LANL Trees. If True, reshaping will only be applied to 3D arrays.
        Defaults to False.
    """
    if isinstance(output_dir, str):
        output_dir = Path(output_dir)
    # Reshape array from (nsp, ny, nx) to (ny, nx, nsp)
    if reshape:
        if len(array.shape) == 3:
            array = np.moveaxis(array, 0, 2).astype(dtype)
        else:
            array = array.astype(dtype)
    else:
        array = array.astype(dtype)

    # Written in row-major order
    with FortranFile(Path(output_dir, dat_name), "w") as f:
        f.write_record(array)
