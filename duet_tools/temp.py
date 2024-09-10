from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from scipy.io import FortranFile


def plot_array(x, title):
    plt.figure(2)
    plt.set_cmap("viridis")
    plt.imshow(x, origin="lower")
    plt.colorbar()
    plt.title(title, fontsize=18)
    plt.show()


def zxy_order_f_flip(
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
            .reshape((nz, nx, ny), order="F")
        )
    return np.moveaxis(array, 2, 1)


def zyx_order_c(
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
            .reshape((nz, ny, nx), order="C")
        )
    return array


def yxz_order_c_move_2_0(
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
            .reshape((ny, nx, nz), order="C")
        )
    return np.moveaxis(array, 2, 0)


def read_dat_to_array(
    directory: str | Path,
    filename: str,
    nx: int,
    ny: int,
    nz: int = None,
    nsp: int = None,
    order: str = "F",
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
        Order of the .dat file. Must be one of "C" or "F". Defaults to "C".

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
        array = FortranFile(fin).read_reals(dtype="float32").reshape(shape, order=order)

    if (nz is None) and (nsp is None):
        return np.moveaxis(array, 1, 0)
    elif nz is None:
        return np.moveaxis(array, 2, 1)
    elif nsp is None:
        return np.transpose(array)
    else:
        return np.moveaxis(np.moveaxis(array, 3, 1), 3, 1)


path = Path(__file__).parent.parent / "tests" / "test-data" / "v2"

nx = 333
ny = 295
nz = 84
nsp = 9

# surface_rhof = zyx_order_c(path, "surface_rhof_layered.dat", nx, ny, nz=nsp)
# plot_array(surface_rhof[0, :, :], "surface_rhof Previous method")
# surface_rhof = zxy_order_f_flip(path, "surface_rhof_layered.dat", nx, ny, nz=nsp)
# plot_array(surface_rhof[0, :, :], "surface_rhof Previous method")
# surface_rhof = yxz_order_c_move_2_0(path, "surface_rhof_layered.dat", nx, ny, nz=nsp)
# plot_array(np.sum(surface_rhof[1:, :, :], axis=0), "surface_rhof Alternate order=C")

# treesrhof = zyx_order_c(path, "treesrhof.dat", nx, ny, nz)
# plot_array(np.sum(treesrhof, axis=0), "treesrhof summed Previous method")
# treesrhof = zxy_order_f_flip(path, "treesrhof.dat", nx, ny, nz)
# plot_array(np.sum(treesrhof, axis=0), "treesrhof summed Jenna's method")
# treesrhof = yxz_order_c_move_2_0(path, "treesrhof.dat", nx, ny, nz)
# plot_array(np.sum(treesrhof, axis=0), "treesrhof summed Alternate order=C")

surface_rhof = read_dat_to_array(path, "surface_rhof_layered.dat", nx, ny, nsp=nsp)
plot_array(np.sum(surface_rhof[1:, :, :], axis=0), "surface_rhof new routine")
treesrhof = read_dat_to_array(path, "treesrhof.dat", nx, ny, nz=nz)
plot_array(np.sum(treesrhof[1:, :, :], axis=0), "treesrhof new routine")
topo = read_dat_to_array(path, "topo.dat", nx, ny)
plot_array(topo, "topo new routine")
