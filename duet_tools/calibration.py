"""
DUET Tools Calibration module
"""
from __future__ import annotations

# Core Imports
from pathlib import Path
import importlib.resources

# External Imports
import numpy as np

# from scipy.io import FortranFile
import geojson

# import landfire
# from landfire.geospatial import get_bbox_from_polygon
# import zipfile
# import rasterio as rio
# import rasterio.mask
# from rasterio.enums import Resampling
# import pandas as pd
# import re
import warnings

# Internal Imports
from duet_tools.utils import read_dat_to_array, write_array_to_dat


try:  # Python 3.9+
    DATA_PATH = importlib.resources.files("duet_tools").joinpath("data")
except AttributeError:  # Python 3.6-3.8
    from pkg_resources import resource_filename

    DATA_PATH = resource_filename("duet_tools", "data")


# TODO: These shouldn't need to be here but testing won't run without them
def fueltype_targets():
    pass


def _validate_target_args():
    pass


def _moisture_weights_from_density():
    pass


class DuetRun:
    """
    Class containing all arrays for a DUET run.
    """

    def __init__(
        self, density: np.ndarray, depth: np.ndarray, moisture: np.ndarray = None
    ):
        self.density = density
        self.moisture = moisture
        self.depth = depth

    def add_moisture_array(self, moisture_array=np.ndarray) -> None:
        """
        Add an array of moisture values to the DUET run

        Parameters
        ----------
        moisture_array : np.ndarray
            3D array of fuel moisture content values. Must be the same shape as
            density and depth arrays already present in the object. Values must
            be positive where fuel is present.

        Returns
        -------
        None :
            Sets the DuetRun.moisture attribute to the input array
        """
        self._validate_input_moisture(moisture_array)
        self.moisture = moisture_array

    def to_quicfire(self, directory: str | Path) -> None:
        """
        Writes a DuetRun object to QUIC-fire fuel .dat inputs to a directory:
        treesrhof.dat, treesmoist.dat, treesfueldepth.dat

        Parameters
        ----------
        directory : str | Path
            Path to directory for writing QUIC-fire files

        Returns
        -------
        None :
            Writes QUIC-Fire .dat files to the provided directory.
        """
        if isinstance(directory, str):
            directory = Path(directory)
        # TODO: write to_quicfire

    def to_numpy(self, fuel_type: str, fuel_parameter: str) -> np.ndarray:
        """
        Returns a numpy array of the provided fuel type and parameter.

        Parameters
        ----------
        fuel_type : str
            Fuel type of desired array. Must be one of "grass", "litter", "separated",
            or "integrated".
            "separated" : returns a 3D array of shape (2,ny,nx), where the first layer
                is grass, and the second layer is litter.
            "integrated" : returns a vertically-integrated array of both fuel types.
                Array remains 3D, which shape (1,nx,ny). Integration method depends on
                fuel parameter.
        fuel_parameter : str
            Fuel parameter of desired array. Must be one of "density", "moisture", or
            "depth".

        Returns
        -------
        np.ndarray :
            3D array of the provided fuel type and parameter.
        """
        fueltypes_allowed = ["grass", "litter", "separated", "integrated"]
        parameters_allowed = ["density", "moisture", "depth"]
        if fuel_type not in fueltypes_allowed:
            raise ValueError(
                f"Fuel type {fuel_type} not supported. Must be one of {fueltypes_allowed}"
            )
        if fuel_parameter not in parameters_allowed:
            raise ValueError(
                f"Fuel parameter {fuel_parameter} not supported. Must be one of {parameters_allowed}"
            )
        # TODO: write to_numpy

    def _validate_input_moisture(self, moisture: np.ndarray):
        if moisture.shape != self.density.shape:
            raise ValueError(
                f"Input array shape {moisture.shape} must match existing arrays {self.density.shape}."
            )
        if self.density[np.where(moisture == 0)].any() != 0:
            raise ValueError(
                "Value of moisture array cannot be zero where fuel is present"
            )


class Targets:
    """
    Class containing and validating target methods and values for fuel parameters.
    """

    def __init__(self, method: str, args: list(str), targets: list):
        self.method = self._validate_method(method)
        self.args, self.targets = self._validate_target_args(method, args, targets)
        self.calibration_function = self._get_calibration_function(method)

    def _get_calibration_function(self, method):
        if method == "maxmin":
            return _maxmin_calibration
        if method == "meansd":
            return _meansd_calibration
        if method == "constant":
            return _constant_calibration

    def _validate_method(self, method: str):
        methods_allowed = ["maxmin", "meansd", "constant", "sb40"]
        if method not in methods_allowed:
            raise ValueError(
                f"Method {method} not supported. Must be one of {methods_allowed}"
            )
        return method

    def _validate_target_args(self, method: str, args: list(str), targets: list(float)):
        method_dict = {
            "constant": ["target"],
            "maxmin": ["max", "min"],
            "meansd": ["mean", "sd"],
            "sb40": ["landfire"],
        }
        args_allowed = method_dict.get(method)
        if set(args_allowed) != set(args):
            raise ValueError(f"Invalid **kwargs for method {method}. Must be {args}")

        targets_dict = dict(zip(args, targets))
        if method == "maxmin":
            if targets_dict["max"] <= targets_dict["min"]:
                raise ValueError("Target maximum must be greater than target minimum")
        if method == "meansd":
            if targets_dict["mean"] < targets_dict["sd"]:
                warnings.warn(
                    "Target mean is smaller than target sd. Were they input correctly?"
                )
        if method == "sb40":
            if not isinstance(targets_dict["landfire"], LandfireQuery):
                raise ValueError(
                    "Value of landfire **kwarg must be of class LandfireQuery. Please use query_landfire()"
                )

        return args, targets


class FuelParameter:
    """
    Class containing and validating calibration targets for a single fuel parameter. Targets can
    be set for multiple fuel types.
    """

    def __init__(self, parameter: str, fuel_types: list(str), targets: list(Targets)):
        self.parameter = self._validate_fuel_parameter(parameter)
        self.fuel_types = self._validate_fuel_types(fuel_types)
        self.targets = targets

    def _validate_fuel_types(self, fuel_types):
        fueltypes_allowed = ["grass", "litter", "all"]
        for fuel_type in fuel_types:
            if fuel_type not in fueltypes_allowed:
                raise ValueError(
                    f"Method {fuel_type} not supported. Must be one of {fueltypes_allowed}"
                )
        if "all" in fuel_types and len(fuel_types) > 1:
            raise ValueError(
                "When fuel parameter targets are assigned to all fuel types, "
                "no other fuel parameter objects should be provided"
            )
        return fuel_types

    def _validate_fuel_parameter(self, parameter):
        fuel_parameters_allowed = ["density", "moisture", "depth"]
        if parameter not in fuel_parameters_allowed:
            raise ValueError(
                f"Fuel parameter {parameter} not supported. Must be one of {fuel_parameters_allowed}"
            )
        return parameter


class LandfireQuery:
    """
    Class containing the information from a LandFire query, to be passed to assign_targets()
    """


def import_duet(directory: str | Path, nx: int, ny: int, nz: int = 2) -> DuetRun:
    """
    Creates a DuetRun object from DUET output files

    Parameters
    ----------
    directory : str | Path
        Path to directory storing the DUET output files surface_rhof.dat and surface_depth.dat
    nx: int
        Number of DUET domain cells in the x-direction
    ny: int
        Number of DUET domain cells in the y-direction
    nz: int
        Number of layers (fuel types) in the DUET outputs. Defaults to 2 (grass and litter).

    Returns
    -------
    Instance of class DuetRun
    """
    if isinstance(directory, str):
        directory = Path(directory)
    density = read_dat_to_array(
        directory=directory, filename="surface_rhof.dat", nx=nx, ny=ny, nz=nz
    )
    depth = read_dat_to_array(
        directory=directory, filename="surface_depth.dat", nx=nx, ny=ny, nz=nz
    )
    return DuetRun(density=density, depth=depth)


def assign_targets(method: str, **kwargs: float | LandfireQuery) -> Targets:
    """
    Assigns target values and calculation method for exactly one fuel type and parameter

    Parameters
    ----------
    method : str
        Calibration method for the target values provided. Must be one of:
        "constant", "maxmin", "meansd", "sb40".
    **kwargs : str
        Keyword arguments correspond to the calibration method.
        #TODO: kwargs in docstring

    Returns
    -------
    Instance of class Targets
    """
    args = list(kwargs.keys())
    targets = list(kwargs.values())

    return Targets(method=method, args=args, targets=targets)


def set_fuel_parameter(parameter: str, **kwargs: Targets):
    """
    Sets calibration targets for grass, litter, both separately, or all
    fuel types together, for a single fuel parameter.

    Parameters
    ----------
    parameter : str
        Fuel parameter for which to set targets
    grass : Targets | None
        Grass calibration targets. Only the grass layer of the DUET bulk
        density array will be calibrated.
    litter : Targets | None
        Litter calibration targets. Only the litter layer of the DUET bulk
        density array will be calibrated.
    all : Targets | None
        Calibration targets for all (both) fuel types. Both layers of the
        DUET bulk density array will be calibrated together.

    Returns
    -------
    FuelParameter :
        Object representing targets for the given fuel parameter, for each provided fuel type
    """
    parameter = parameter
    fuel_types = list(kwargs.keys())
    targets = list(kwargs.values())

    return FuelParameter(parameter, fuel_types, targets)


def set_density(**kwargs: Targets):
    """
    Sets bulk density calibration targets for grass, litter, both separately, or all
    fuel types together.

    Parameters
    ----------
    grass : Targets | None
        Grass bulk density calibration targets. Only the grass layer of the DUET bulk
        density array will be calibrated.
    litter : Targets | None
        Litter bulk density calibration targets. Only the litter layer of the DUET bulk
        density array will be calibrated.
    all : Targets | None
        Bulk density calibration targets for all (both) fuel types. Both layers of the
        DUET bulk density array will be calibrated together.

    Returns
    -------
    FuelParameter :
        Object representing bulk density targets for each provided fuel type
    """
    parameter = "density"
    fuel_types = list(kwargs.keys())
    targets = list(kwargs.values())

    return FuelParameter(parameter, fuel_types, targets)


def set_moisture(**kwargs: Targets):
    """
    Sets moisture calibration targets for grass, litter, both separately, or all
    fuel types together.

    Parameters
    ----------
    grass : Targets | None
        Grass moisture calibration targets. Only the grass layer of the DUET bulk
        density array will be calibrated.
    litter : Targets | None
        Litter moisture calibration targets. Only the litter layer of the DUET bulk
        density array will be calibrated.
    all : Targets | None
        Moisture calibration targets for all (both) fuel types. Both layers of the
        DUET bulk density array will be calibrated together.

    Returns
    -------
    FuelParameter :
        Object representing moisture targets for each provided fuel type
    """
    parameter = "moisture"
    fuel_types = list(kwargs.keys())
    targets = list(kwargs.values())

    return FuelParameter(parameter, fuel_types, targets)


def set_depth(**kwargs: Targets):
    """
    Sets depth calibration targets for grass, litter, both separately, or all
    fuel types together.

    Parameters
    ----------
    grass : Targets | None
        Grass depth calibration targets. Only the grass layer of the DUET bulk
        density array will be calibrated.
    litter : Targets | None
        Litter depth calibration targets. Only the litter layer of the DUET bulk
        density array will be calibrated.
    all : Targets | None
        Depth calibration targets for all (both) fuel types. Both layers of the
        DUET bulk density array will be calibrated together.

    Returns
    -------
    FuelParameter :
        Object representing depth targets for each provided fuel type
    """
    parameter = "depth"
    fuel_types = list(kwargs.keys())
    targets = list(kwargs.values())

    return FuelParameter(parameter, fuel_types, targets)


def calibrate(
    duet_run: DuetRun, fuel_parameter_targets: list(FuelParameter) | FuelParameter
) -> DuetRun:
    """
    Calibrates the arrays in a DuetRun object using the provided targets and methods for one
    or more fuel types.

    Parameters
    ----------
    duet_run : DuetRun
        The DUET run to calibrate

    fuel_type_targets : FuelParameters | list(FuelParameters)
        FuelParameters object or list of FuelParameters objects for the fuel types
        to be calibrated.

    Returns
    -------
    Instance of class DuetRun with calibrated fuel arrays
    """
    if isinstance(fuel_parameter_targets, FuelParameter):
        fuel_parameter_targets = [fuel_parameter_targets]

    calibrated_duet = _duplicate_duet_run(duet_run)
    for fuelparameter in fuel_parameter_targets:
        fuelparam = fuelparameter.parameter
        for i in range(len(fuelparameter.fuel_types)):
            fueltype = fuelparameter.fuel_types[i]
            array_to_calibrate = _get_array_to_calibrate(duet_run, fueltype, fuelparam)
            calibrated_array = _do_calibration(
                array_to_calibrate, fuelparameter.targets[i]
            )
            calibrated_duet = _add_calibrated_array(
                calibrated_duet, calibrated_array, fueltype, fuelparam
            )
    return calibrated_duet


def get_unit_from_fastfuels(zroot):
    """
    Creates a geojson bounding box of a fastfuels domain.

    Returns
    -------
    geojson
    """
    # TODO: write get_unit_from_fastfuels


def get_unit_from_shapefile(directory: str | Path):
    """
    Reads in a shapefile and returns a geojson bounding box.

    Returns
    -------
    geojson
    """
    # TODO: write get_unit_from_shapefile


def write_numpy_to_quicfire(array: np.ndarray, directory: str | Path, filename: str):
    if isinstance(directory, str):
        directory = Path(directory)
    write_array_to_dat(array=array, dat_name=filename, output_dir=directory)


def _get_array_to_calibrate(duet_run: DuetRun, fueltype: str, fuelparam: str):
    if fuelparam == "density":
        if fueltype == "grass":
            return duet_run.density[0, :, :].copy()
        if fueltype == "litter":
            return duet_run.density[1, :, :].copy()
        return np.sum(duet_run.density, axis=0)
    if fuelparam == "depth":
        if fueltype == "grass":
            return duet_run.depth[0, :, :].copy()
        if fueltype == "litter":
            return duet_run.depth[1, :, :].copy()
        return np.max(duet_run.depth, axis=0)
    if fuelparam == "moisture":
        if duet_run.moisture:
            if fueltype == "grass":
                return duet_run.moisture[0, :, :].copy()
            if fueltype == "litter":
                return duet_run.moisture[1, :, :].copy()
            density_weights = duet_run.density.copy()
            density_weights[density_weights == 0] = 0.01
            return np.average(duet_run.moisture, weights=density_weights, axis=0)
        raise ValueError(
            "No moisture array available to calibrate. Please add moisture"
            "array using DuetRun.add_moisture_array"
        )


def _duplicate_duet_run(duet_run: DuetRun) -> DuetRun:
    new_density = duet_run.density.copy() if duet_run.density is not None else None
    new_moisture = duet_run.moisture.copy() if duet_run.moisture is not None else None
    new_depth = duet_run.depth.copy() if duet_run.depth is not None else None

    new_duet = DuetRun(
        density=new_density,
        moisture=new_moisture,
        depth=new_depth,
    )

    return new_duet


def _do_calibration(array: np.ndarray, target_obj: Targets):
    kwarg_dict = {}
    for i in range(len(target_obj.args)):
        kwarg_dict[target_obj.args[i]] = target_obj.targets[i]
    new_array = target_obj.calibration_function(array, **kwarg_dict)
    return new_array


def _maxmin_calibration(x: np.ndarray, **kwargs: float) -> np.ndarray:
    """
    Scales and shifts values in a numpy array based on an observed range. Does not assume
    data is normally distributed.
    """
    max_val = kwargs["max"]
    min_val = kwargs["min"]
    x1 = x[x > 0]
    if np.max(x1) == np.min(x1):
        raise ValueError(
            "maxmin calibration cannot be used when array has only one positive value. "
            "Please use 'constant' calibration method"
        )
    x2 = (x1 - np.min(x1)) / (np.max(x1) - np.min(x1))
    x3 = x2 * (max_val - min_val)
    x4 = x3 + min_val
    xnew = x.copy()
    xnew[np.where(x > 0)] = x4
    return xnew


def _meansd_calibration(x: np.ndarray, **kwargs: float) -> np.ndarray:
    """
    Scales and shifts values in a numpy array based on an observed mean and standard deviation.
    Assumes data is normally distributed.
    """
    mean_val = kwargs["mean"]
    sd_val = kwargs["sd"]
    x1 = x[x > 0]
    if np.max(x1) == np.min(x1):
        raise ValueError(
            "meansd calibration should not be used when array has only one positive value. "
            "Please use 'constant' calibration method"
        )
    x2 = mean_val + (x1 - np.mean(x1)) * (sd_val / np.std(x1))
    xnew = x.copy()
    xnew[np.where(x > 0)] = x2
    if np.min(xnew) < 0:
        xnew = _truncate_at_0(xnew)
    return xnew


# TODO: add option for fuel vs cell bulk density?
def _constant_calibration(x: np.ndarray, **kwargs: float) -> np.ndarray:
    constant = kwargs["target"]
    arr = x.copy()
    arr[arr > 0] = constant
    return arr


def _add_calibrated_array(
    duet_to_calibrate: DuetRun,
    calibrated_array: np.ndarray,
    fueltype: str,
    fuelparam: str,
) -> DuetRun:
    for param in ["density", "moisture", "depth"]:
        if fuelparam == param:
            if fueltype == "grass":
                duet_to_calibrate.__dict__[param][0, :, :] = calibrated_array
            if fueltype == "litter":
                duet_to_calibrate.__dict__[param][1, :, :] = calibrated_array
            if fueltype == "all":
                duet_to_calibrate.__dict__[param] = _separate_2d_array(
                    calibrated_array, param, duet_to_calibrate.density
                )
    return duet_to_calibrate


def _separate_2d_array(
    calibrated: np.ndarray, param: str, duet_run: DuetRun
) -> np.ndarray:
    separated = np.array([calibrated, calibrated])
    if param == "density":
        weights = duet_run.density.copy()
        weights[0, :, :] = duet_run.density[0, :, :] / np.sum(duet_run.density, axis=0)
        weights[1, :, :] = duet_run.density[1, :, :] / np.sum(duet_run.density, axis=0)
        separated[0, :, :] = calibrated * weights[0, :, :]
        separated[1, :, :] = calibrated * weights[1, :, :]
    if param == "moisture":
        separated[0, :, :][np.where(duet_run.moisture[0, :, :] == 0)] = 0
        separated[1, :, :][np.where(duet_run.moisture[1, :, :] == 0)] = 0
    if param == "depth":
        weights = duet_run.depth.copy()
        weights[0, :, :] = duet_run.depth[0, :, :] / np.max(duet_run.depth, axis=0)
        weights[1, :, :] = duet_run.depth[1, :, :] / np.max(duet_run.depth, axis=0)
        separated[0, :, :] = calibrated * weights[0, :, :]
        separated[1, :, :] = calibrated * weights[1, :, :]
    return separated


def _truncate_at_0(arr: np.ndarray) -> np.ndarray:
    """
    Artificially truncates data to positive values by scaling all values below the median
    to the range (0, mean), effectively "compressing" those values.
    """
    arr2 = arr.copy()
    bottom_half = arr2[arr2 < np.median(arr2)]
    squeezed = (bottom_half - np.min(bottom_half)) / (
        np.max(bottom_half) - np.min(bottom_half)
    ) * (np.median(arr2) - 0) + 0
    arr2[np.where(arr2 < np.median(arr2))] = squeezed
    arr2[np.where(arr == 0)] = 0
    return arr2


# class DuetCalibrator(BaseModel):
#     # TODO: Instead of saving just the most recent array to self.calibrated_array, append to a list or a dict of calibrated arrays, so that you can access a bunch without reading in dat files
#     nx: PositiveInt
#     ny: PositiveInt
#     nz: PositiveInt
#     dx: PositiveInt
#     dy: PositiveInt
#     xmin: float
#     ymin: float
#     xmax: float
#     ymax: float
#     output_dir: Union[str, Path]
#     calibrated: bool = False
#     calibrated_array: Optional[np.ndarray] = None
#     calibrated_fuel_type: list = []
#     calibration_method: list = []
#     saved_files: list = []

#     @computed_field
#     @property
#     def original_duet_array(self) -> np.ndarray:
#         return self._read_original_duet()

#     @computed_field
#     @property
#     def duet_dict(self) -> dict:
#         return self._get_input_array()

#     def calibrate_max_min(
#         self, fuel_type: str | list, max_val: float | list, min_val: float | list
#     ) -> None:
#         """
#         Calibrate the values of the surface bulk density output from DUET by setting the
#         range (maximum and minimum).

#         Parameters
#         ----------

#         fuel_type : str | list
#             Fuel type(s) to calibrate. May be one of "total", "grass", or "litter" given
#             as a string, or both "grass" and "litter" given as a list. When "total", both
#             fuel types are calibrated based on one set of inputs. When "litter" or "grass",
#             the given fuel type is calibrated, the other is left unchanged and added to the
#             calibrated fuel type to produce the final array. When ["grass","litter"] or
#             ["litter","grass"] both fuel types are calibrated based on their respective inputs.

#         max_val : float | list
#             Target maximim value for calibration. If multiple values are given for multiple fuel
#             types, the position of each value in the list must match the position of their
#             corresponding fuel type.

#         min_val : float | list
#             Target minimum value for calibration. If multiple values are given for multiple
#             fuel types, the position of each value in the list must match the position of their
#             corresponding fuel type.

#         Returns
#         -------
#         None
#             Calibrated array of DUET surface bulk density is saved to the output directory.
#             Filename indicates the fuel type and the max/min calibration method, and is
#             incremented if previous calibrations of the same fuel type and method have been
#             conducted.

#         """
#         self._validate_inputs(fuel_type, max_val, min_val)
#         if isinstance(fuel_type, str):
#             fuel_type = [fuel_type]
#         if isinstance(max_val, int) or isinstance(max_val, float):
#             max_val = [max_val]
#         if isinstance(min_val, int) or isinstance(min_val, float):
#             min_val = [min_val]
#         calibrated = {}
#         for f in range(len(fuel_type)):
#             arr = self.duet_dict[fuel_type[f]]
#             calibrated[fuel_type[f]] = self._maxmin_calibration(
#                 arr, max_val[f], min_val[f]
#             )
#         self.calibrated_array = self._combine_fuel_types(calibrated_dict=calibrated)
#         self.calibrated = True
#         self.calibrated_fuel_type.append(fuel_type)
#         self.calibrated_fuel_type = self._flatten(self.calibrated_fuel_type)
#         self.calibration_method.append("maxmin")
#         self.duet_dict = self._get_input_array()

#     def calibrate_mean_sd(
#         self, fuel_type: str | list, mean_val: float | list, sd_val: float | list
#     ) -> None:
#         """
#         Calibrate the values of the surface bulk density output from DUET by setting the
#         center and spread (mean and standard deviation).

#         Parameters
#         ----------

#         fuel_type : str | list
#             Fuel type(s) to calibrate. May be one of "total", "grass", or "litter" given
#             as a string, or both "grass" and "litter" given as a list. When "total", both
#             fuel types are calibrated based on one set of inputs. When "litter" or "grass",
#             the given fuel type is calibrated, the other is left unchanged and added to the
#             calibrated fuel type to produce the final array. When ["grass","litter"] or
#             ["litter","grass"] both fuel types are calibrated based on their respective inputs.

#         mean_val : float | list
#             Target mean value for calibration. If multiple values are given for multiple fuel
#             types, the position of each value in the list must match the position of their
#             corresponding fuel type.

#         sd_val : float | list
#             Target standard deviation for calibration. If multiple values are given for multiple
#             fuel types, the position of each value in the list must match the position of their
#             corresponding fuel type.

#         Returns
#         -------
#         None
#             Calibrated array of DUET surface bulk density is saved to the output directory.
#             Filename indicates the fuel type and the max/min calibration method, and is
#             incremented if previous calibrations of the same fuel type and method have been
#             conducted.

#         """
#         self._validate_inputs(fuel_type, mean_val, sd_val)
#         if isinstance(fuel_type, str):
#             fuel_type = [fuel_type]
#         if isinstance(mean_val, int) or isinstance(mean_val, float):
#             mean_val = [mean_val]
#         if isinstance(sd_val, int) or isinstance(sd_val, float):
#             sd_val = [sd_val]
#         calibrated = {}
#         for f in range(len(fuel_type)):
#             arr = self.duet_dict[fuel_type[f]]
#             calibrated[fuel_type[f]] = self._meansd_calibration(
#                 arr, mean_val[f], sd_val[f]
#             )
#         self.calibrated_array = self._combine_fuel_types(calibrated_dict=calibrated)
#         self.calibrated = True
#         self.calibrated_fuel_type.append(fuel_type)
#         self.calibrated_fuel_type = self._flatten(self.calibrated_fuel_type)
#         self.calibration_method.append("meansd")
#         self.duet_dict = self._get_input_array()

#     def calibrate_with_sb40(self, fuel_type: str | list) -> None:
#         self._validate_inputs(fuel_type)
#         print("Querying LandFire...\n")
#         # Query Landfire and return array of SB40 keys
#         sb40_arr = self._query_landfire()
#         # Import SB40 FBFM parameters table
#         sb40_params_path = DATA_PATH / "sb40_parameters.csv"
#         sb40_params = pd.read_csv(sb40_params_path)
#         # Generate dict of fastfuels bulk density values and apply to Landfire query
#         sb40_dict = self._get_sb40_fuel_params(sb40_params)
#         sb40_ftype, sb40_rhof = self._get_sb40_arrays(sb40_arr, sb40_dict)
#         if isinstance(fuel_type, str):
#             fuel_type = [fuel_type]
#         calibrated = {}
#         for f in fuel_type:
#             if f == "grass":
#                 if 1 in sb40_ftype:
#                     max_val = np.max(sb40_rhof[sb40_ftype == 1])
#                     grass_arr = sb40_rhof[sb40_ftype == 1]
#                     min_val = np.min(grass_arr[grass_arr > 0])
#                 else:
#                     print(
#                         "WARNING: grass fuel not present in sb40. Continuing with litter only."
#                     )
#                     max_val = 0
#                     min_val = 0
#             elif f == "litter":
#                 if -1 in sb40_ftype:
#                     max_val = np.max(sb40_rhof[sb40_ftype == -1])
#                     litter_arr = sb40_rhof[sb40_ftype == -1]
#                     min_val = np.min(litter_arr[litter_arr > 0])
#                 else:
#                     print(
#                         "WARNING: litter fuel not present in sb40. Continuing with grass only."
#                     )
#                     max_val = 0
#                     min_val = 0
#             else:
#                 max_val = np.max(sb40_rhof)
#                 min_val = np.min(sb40_rhof[sb40_rhof > 0])
#             calibrated[f] = self._maxmin_calibration(
#                 self.duet_dict[f], max_val, min_val
#             )
#         self.calibrated_array = self._combine_fuel_types(calibrated_dict=calibrated)
#         self.calibrated = True
#         self.calibrated_fuel_type.append(fuel_type)
#         self.calibrated_fuel_type = self._flatten(self.calibrated_fuel_type)
#         self.calibration_method.append("sb40")
#         self.duet_dict = self._get_input_array()

#     def revert_to_original_duet(self, delete_files: bool = False) -> None:
#         """
#         Ensure that the next calibration will be conducted on the original DUET output and
#         optionally delete all files saved from previous calibrations of the DuetCalibrator
#         instance.

#         Parameters
#         ----------

#         delete_files : bool
#             Whether to delete the previously saved .dat files. Default is False,
#             meaning files will not be deleted and any subsequent calibrations of the
#             same method and fuel type will be saved with incremented filenames.

#         """
#         if delete_files:
#             [
#                 Path(self.output_dir, self.saved_files[file]).unlink()
#                 for file in range(len(self.saved_files))
#                 if Path(self.output_dir, self.saved_files[file]).exists()
#             ]
#             self.saved_files = []
#         self.calibrated = False
#         self.calibrated_array = None
#         self.calibrated_fuel_type = []
#         self.calibration_method = []
#         self.original_duet_array = self._read_original_duet()
#         self.duet_dict = self._get_input_array()

#     def to_file(self) -> None:
#         """
#         Write the most recently calibrated surface fuel array to a .dat file.
#         File will be saved to the output directory of the DuetCalibrator instance.
#         """
#         if self.calibrated:
#             arr_name = self._name_calibrated_file()
#             _write_np_array_to_dat(self.calibrated_array, arr_name, self.output_dir)
#             self.saved_files.append(arr_name)
#         else:
#             raise Exception("Must calibrate array before writing to file.")

#     def replace_quicfire_surface_fuels(self):
#         """
#         Replace surface fuel bulk density in quicfire output
#         (from export_zarr_to_quicfire) with DUET output.

#         Parameters
#         ----------
#         quicfire_dir: Path | str
#             Directory where QUIO-Fire .dat files are located,
#             and to where updated .dat files are written to.

#         Returns
#         -------
#         None
#             Modified bulk density array (treesrhof.dat) is written to the QUIC-Fire directory
#         """
#         with open(Path(self.output_dir, "treesrhof.dat"), "rb") as fin:
#             qf_arr = (
#                 FortranFile(fin)
#                 .read_reals(dtype="float32")
#                 .reshape((self.nz, self.ny, self.nx), order="C")
#             )
#         if self.calibrated:
#             tag = "calibrated"
#             duet_arr = np.add(
#                 self.calibrated_array[0, :, :], self.calibrated_array[1, :, :]
#             )
#         else:
#             tag = "unmodified"
#             duet_arr = np.add(
#                 self.original_duet_array[0, :, :], self.original_duet_array[1, :, :]
#             )
#         qf_arr[0, :, :] = duet_arr
#         _write_np_array_to_dat(
#             qf_arr, "treesrhof.dat", self.output_dir, np.float32, reshape=False
#         )
#         print(
#             "Replaced FastFuels surface fuel layer with {} DUET surface fuels".format(
#                 tag
#             )
#         )

#     def _validate_inputs(self, fuel_type, val1=None, val2=None):
#         # Validate fuel types
#         valid_ftypes = [
#             "litter",
#             "grass",
#             "total",
#             ["litter", "grass"],
#             ["grass", "litter"],
#         ]
#         if fuel_type not in valid_ftypes:
#             raise ValueError(
#                 "Invalid fuel type. Must be one of {}.".format(valid_ftypes)
#             )
#         if fuel_type == "total" and self.calibrated == True:
#             raise ValueError(
#                 "Invalide fuel type: 'total' fuel calibration cannot be applied to a previously calibrated array. Choose a different fuel type or use revert_to_original_duet() before calibrating total fuels."
#             )
#         if fuel_type in self.calibrated_fuel_type:
#             warnings.warn(
#                 "Fuel type '{}' already calibrated. Replacing previous calibrated values.".format(
#                     fuel_type
#                 )
#             )

#         # Validate fuel summary arguments
#         if val1 is not None:
#             if isinstance(fuel_type, list):
#                 if not isinstance(val1, list) or not isinstance(val2, list):
#                     raise TypeError(
#                         "Fuel input values must be provided as a list when fuel_type is a list."
#                     )
#                 if len(fuel_type) != len(val1) | len(val2):
#                     raise ValueError(
#                         "Number of fuel value inputs must match number of fuel types ({} fuel inputs for fuel_type = {}).".format(
#                             len(fuel_type), fuel_type
#                         )
#                     )

#     def _maxmin_calibration(
#         self, x: np.ndarray, max_val: float | int, min_val: float | int
#     ) -> np.ndarray:
#         """
#         Scales and shifts values in a numpy array based on an observed range. Does not assume
#         data is normally distributed.
#         """
#         x1 = x[x > 0]
#         x2 = (x1 - np.min(x1)) / (np.max(x1) - np.min(x1))
#         x3 = x2 * (max_val - min_val)
#         x4 = x3 + min_val
#         xnew = x.copy()
#         xnew[np.where(x > 0)] = x4
#         return xnew

#     def _meansd_calibration(
#         self, x: np.ndarray, mean: float | int, sd: float | int
#     ) -> np.ndarray:
#         """
#         Scales and shifts values in a numpy array based on an observed mean and standard deviation.
#         Assumes data is normally distributed.
#         """
#         x1 = x[x > 0]
#         x2 = mean + (x1 - np.mean(x1)) * (sd / np.std(x1))
#         xnew = x.copy()
#         xnew[np.where(x > 0)] = x2
#         if np.min(xnew) < 0:
#             xnew = self._truncate_at_0(xnew)
#         return xnew

#     def _truncate_at_0(self, arr: np.ndarray) -> np.ndarray:
#         """
#         Artificially truncates data to positive values by scaling all values below the median
#         to the range (0, mean), effectively "compressing" those values.
#         """
#         arr2 = arr.copy()
#         bottom_half = arr2[arr2 < np.median(arr2)]
#         squeezed = (bottom_half - np.min(bottom_half)) / (
#             np.max(bottom_half) - np.min(bottom_half)
#         ) * (np.median(arr2) - 0) + 0
#         arr2[np.where(arr2 < np.median(arr2))] = squeezed
#         arr2[np.where(arr == 0)] = 0
#         return arr2

#     def _query_landfire(self, delete_files: bool = True) -> np.ndarray:
#         """
#         Download a grid of SB40 fuel models from Landfire for the unit and convert to a numpy array

#         Parameters
#         ----------
#         delete_files: bool = True
#             Whether to delete intermediate tif files. Defaults to True

#         Returns
#         -------
#         NumPy Array
#             A numpy array of the SB40 FBFM keys for the site
#         """

#         # Name intermediate files
#         temp = [
#             "landfire_sb40.zip",
#             "landfire_sb40.tif",
#             "sb40_upsampled.tif",
#             "sb40_cropped.tif",
#         ]

#         # Create a bounding box from fuelgrid zarr
#         coords = [
#             self.xmin,
#             self.ymin,
#             self.xmin,
#             self.ymax,
#             self.xmax,
#             self.ymax,
#             self.xmax,
#             self.ymin,
#             self.xmin,
#             self.ymin,
#         ]
#         poly = geojson.Polygon(coordinates=[coords], precision=8)
#         bbox = get_bbox_from_polygon(aoi_polygon=poly, crs=5070)

#         # Download Landfire data to output directory
#         lf = landfire.Landfire(bbox, output_crs="5070")
#         lf.request_data(
#             layers=["200F40_19"], output_path=Path(self.output_dir, "landfire_sb40.zip")
#         )

#         # Exctract tif from compressed download folder and rename
#         with zipfile.ZipFile(Path(self.output_dir, "landfire_sb40.zip")) as zf:
#             extension = ".tif"
#             rename = "landfire_sb40.tif"
#             info = zf.infolist()
#             for file in info:
#                 if file.filename.endswith(extension):
#                     file.filename = rename
#                     zf.extract(file, self.output_dir)

#         # Upsample landfire raster to the quicfire resolution
#         with rio.open(Path(self.output_dir, "landfire_sb40.tif")) as sb:
#             upscale_factor = 30 / self.dx  # lf res/qf res
#             profile = sb.profile.copy()
#             # resample data to target shape
#             data = sb.read(
#                 out_shape=(
#                     sb.count,
#                     int(sb.height * upscale_factor),
#                     int(sb.width * upscale_factor),
#                 ),
#                 resampling=Resampling.nearest,
#             )

#             # scale image transform
#             transform = sb.transform * sb.transform.scale(
#                 (sb.width / data.shape[-1]), (sb.height / data.shape[-2])
#             )
#             profile.update(
#                 {
#                     "height": data.shape[-2],
#                     "width": data.shape[-1],
#                     "transform": transform,
#                 }
#             )
#             with rio.open(
#                 Path(self.output_dir, "sb40_upsampled.tif"), "w", **profile
#             ) as dataset:
#                 dataset.write(data)

#         # Crop the upsampled raster to the unit bounds
#         with rio.open(Path(self.output_dir, "sb40_upsampled.tif"), "r+") as rst:
#             out_image, out_transform = rasterio.mask.mask(rst, [poly], crop=True)
#             out_meta = rst.meta
#             out_meta.update(
#                 {
#                     "driver": "GTiff",
#                     "height": out_image.shape[1],
#                     "width": out_image.shape[2],
#                     "transform": out_transform,
#                 }
#             )
#             with rio.open(
#                 Path(self.output_dir, "sb40_cropped.tif"), "w", **out_meta
#             ) as cropped:
#                 cropped.write(out_image)

#         # Read in the cropped raster as a numpy array
#         with rio.open(Path(self.output_dir, "sb40_cropped.tif")) as rst:
#             arr = rst.read(1)

#         if delete_files:
#             [
#                 Path(self.output_dir, file).unlink()
#                 for file in temp
#                 if Path(self.output_dir, file).exists()
#             ]

#         return arr[arr > 0]

#     # TODO: fix rasterio cropping issue (grr) so that landfire raster is same size as fuelgrid

#     def _get_sb40_fuel_params(self, params: pd.DataFrame = None) -> dict:
#         """
#         Builds a dictionary of SB40 fuel parameter values and converts them to
#         the official FastFuels units

#         Returns:
#             dict: SB40 parameters for each fuel model
#         """
#         # Load in the SB40 fuel parameters
#         if params is not None:
#             sb40_df = params.copy()
#         else:
#             fpath = Path("src", "data_module", "data", "sb40_parameters.csv")
#             with open(fpath) as fin:
#                 sb40_df = pd.read_csv(fin)

#         # Convert tons/ac-ft to kg/m^3
#         sb40_df["1_hr_kg_per_m3"] = sb40_df["1_hr_t_per_ac"] * 0.22417
#         sb40_df["10_hr_kg_per_m3"] = sb40_df["10_hr_t_per_ac"] * 0.22417
#         sb40_df["100_hr_kg_per_m3"] = sb40_df["100_hr_t_per_ac"] * 0.22417
#         sb40_df["live_herb_kg_per_m3"] = sb40_df["live_herb_t_per_ac"] * 0.22417
#         sb40_df["live_woody_kg_per_m3"] = sb40_df["live_woody_t_per_ac"] * 0.22417

#         # Convert inverse feet to meters
#         sb40_df["dead_1_hr_sav_ratio_1_per_m"] = (
#             sb40_df["dead_1_hr_sav_ratio_1_per_ft"] * 3.2808
#         )
#         sb40_df["live_herb_sav_ratio_1_per_m"] = (
#             sb40_df["live_herb_sav_ratio_1_per_ft"] * 3.2808
#         )
#         sb40_df["live_wood_sav_ratio_1_per_m"] = (
#             sb40_df["live_wood_sav_ratio_1_per_ft"] * 3.2808
#         )

#         # Convert percent to ratio
#         sb40_df["dead_fuel_extinction_moisture"] /= 100

#         # Convert feet to meters
#         sb40_df["fuel_bed_depth_m"] = sb40_df["fuel_bed_depth_ft"] * 0.3048

#         # Compute wet loading
#         sb40_df["wet_load"] = sb40_df["1_hr_kg_per_m3"] + sb40_df["live_herb_kg_per_m3"]

#         # Compute a live herb curing factor alpha as a function of wet loading.
#         # This is kind of a B.S. approach raised by Rod on a phone call with
#         # Anthony on 02/28/2023. I don't like this at all, but it is a temporary
#         # Fix for the BP3D team to run some simulations.
#         # low_load_fuel_models = [
#         sb40_df["alpha"] = [0.5 if rho > 1 else 1.0 for rho in sb40_df["wet_load"]]

#         # Compute dry loading
#         sb40_df["dry_herb_load"] = sb40_df["live_herb_kg_per_m3"] * sb40_df["alpha"]
#         sb40_df["dry_load"] = sb40_df["1_hr_kg_per_m3"] + sb40_df["dry_herb_load"]

#         # Compute SAV
#         sb40_df["sav_1hr_ratio"] = sb40_df["1_hr_kg_per_m3"] / sb40_df["dry_load"]
#         sb40_df["sav_1hr"] = (
#             sb40_df["sav_1hr_ratio"] * sb40_df["dead_1_hr_sav_ratio_1_per_m"]
#         )
#         sb40_df["sav_herb_ratio"] = sb40_df["dry_herb_load"] / sb40_df["dry_load"]
#         sb40_df["sav_herb"] = (
#             sb40_df["sav_herb_ratio"] * sb40_df["live_herb_sav_ratio_1_per_m"]
#         )
#         sb40_df["sav"] = sb40_df["sav_1hr"] + sb40_df["sav_herb"]

#         # Convert nan to 0
#         sb40_df.fillna(0, inplace=True)

#         ## NIKO ADDITIONS
#         # Create dictionary for assigning fuel types for DUET calibration
#         duet_dict = {
#             "NB": 0,  # 0 = NEUTRAL, i.e. not predominantly grass or litter
#             "GR": 1,  # 1 = GRASS predominantly
#             "GS": 1,
#             "SH": 1,  # I am considering shrubs as grass
#             "TU": 0,
#             "TL": -1,  # -1 = LITTER predominantly
#             "SB": 0,
#         }

#         # Add column to df with DUET designations
#         pattern = r"[0-9]"  # take out numbers from fbfm_type strings
#         sb40_df["fbfm_cat"] = sb40_df["fbfm_code"].apply(
#             lambda x: re.sub(pattern, "", x)
#         )
#         sb40_df["duet_fuel_type"] = sb40_df["fbfm_cat"].apply(
#             lambda x: duet_dict.get(x)
#         )
#         ## END NIKO ADDITIONS

#         # Build the dictionary with fuel parameters for the Scott and Burgan 40
#         # fire behavior fuel models. Dict format: key ->
#         # [name, loading (tons/ac), SAV (1/ft), ext. MC (percent), bed depth (ft)]
#         # Note: Eventually we want to get rid of this and just use the dataframe.
#         # This is legacy from the old parameter table json.
#         sb40_dict = {}
#         for key in sb40_df["key"]:
#             row = sb40_df[sb40_df["key"] == key]
#             sb40_dict[key] = [
#                 row["fbfm_code"].values[0],
#                 row["dry_load"].values[0],
#                 row["sav"].values[0],
#                 row["dead_fuel_extinction_moisture"].values[0],
#                 row["fuel_bed_depth_m"].values[0],
#                 row["duet_fuel_type"].values[0],
#             ]

#         return sb40_dict

#     def _get_sb40_arrays(self, sb40_keys: np.ndarray, sb40_dict: dict) -> tuple:
#         """
#         Use a dictionary of bulk density and fuel types that correspond to SB40
#         fuel models to assign those values across the study area.

#         Fuel types are as follows:
#         - 1: Predominantly grass. All cells with a GR, GS, or SH designation from SB40.
#         - -1: Predominantly tree litter. All cells with a TL designation from SB40.
#         - 0: Neither predominantly grass or tree litter. All other SB40 designations.

#         Returns:
#             - np.ndarray of fuel types
#             - np.ndarray of bulk density values as calculated by fastfuels
#         """
#         ftype_dict = {key: val[5] for key, val in sb40_dict.items()}
#         ftype_arr = np.vectorize(ftype_dict.get)(sb40_keys)

#         rhof_dict = {key: val[1] for key, val in sb40_dict.items()}
#         rhof_arr = np.vectorize(rhof_dict.get)(sb40_keys)

#         return ftype_arr, rhof_arr

#     def _combine_fuel_types(self, calibrated_dict) -> np.ndarray:
#         calibrated_duet = np.zeros((2, self.ny, self.nx))
#         if len(calibrated_dict.keys()) == 1:
#             ftype = list(calibrated_dict.keys())[0]
#             if ftype == "grass":
#                 calibrated_duet[0, :, :] = calibrated_dict["grass"]
#                 calibrated_duet[1, :, :] = self.duet_dict["litter"]
#             elif ftype == "litter":
#                 calibrated_duet[0, :, :] = self.duet_dict["grass"]
#                 calibrated_duet[1, :, :] = calibrated_dict["litter"]
#             else:
#                 grass_weights = self.duet_dict["grass"] / self.duet_dict["total"]
#                 litter_weights = self.duet_dict["litter"] / self.duet_dict["total"]
#                 calibrated_duet[0, :, :][
#                     np.where(self.original_duet_array[0, :, :] > 0)
#                 ] = (calibrated_dict["total"] * grass_weights)[
#                     np.where(self.original_duet_array[0, :, :] > 0)
#                 ]
#                 calibrated_duet[1, :, :][
#                     np.where(self.original_duet_array[0, :, :] > 0)
#                 ] = (calibrated_dict["total"] * litter_weights)[
#                     np.where(self.original_duet_array[0, :, :] > 0)
#                 ]
#         else:
#             calibrated_duet[0, :, :] = calibrated_dict["grass"]
#             calibrated_duet[1, :, :] = calibrated_dict["litter"]
#         return calibrated_duet

#     def _get_input_array(self):
#         duet_dict = {}
#         if self.calibrated:
#             duet_dict["grass"] = self.calibrated_array[0, :, :]
#             duet_dict["litter"] = self.calibrated_array[1, :, :]
#             duet_dict["total"] = np.add(
#                 self.calibrated_array[0, :, :], self.calibrated_array[1, :, :]
#             )
#         else:
#             duet_dict["grass"] = self.original_duet_array[0, :, :]
#             duet_dict["litter"] = self.original_duet_array[1, :, :]
#             duet_dict["total"] = np.add(
#                 self.original_duet_array[0, :, :], self.original_duet_array[1, :, :]
#             )
#         return duet_dict

#     def _name_calibrated_file(self) -> str:
#         delim = "_"
#         ftype_str = (
#             self.calibrated_fuel_type[0]
#             if len(self.calibrated_fuel_type) == 1
#             else delim.join([str(ele) for ele in self.calibrated_fuel_type])
#         )
#         method_str = (
#             self.calibration_method[0]
#             if len(self.calibration_method) == 1
#             else delim.join([str(ele) for ele in self.calibration_method])
#         )
#         arr_name = delim.join(
#             ["surface_rhof_calibrated", ftype_str, "{}.dat".format(method_str)]
#         )
#         if Path(self.output_dir, arr_name).exists():
#             i = 2
#             while Path(
#                 self.output_dir,
#                 delim.join(
#                     ["surface_rhof_calibrated", ftype_str, method_str, "%s.dat" % i]
#                 ),
#             ).exists():
#                 i += 1
#             arr_name = delim.join(
#                 ["surface_rhof_calibrated", ftype_str, method_str, "%s.dat" % i]
#             )
#         return arr_name

#     def _flatten(self, A):
#         rt = []
#         for i in A:
#             if isinstance(i, list):
#                 rt.extend(self._flatten(i))
#             else:
#                 rt.append(i)
#         return rt


# def _write_np_array_to_dat(
#     array: np.ndarray,
#     dat_name: str,
#     output_dir: Path,
#     dtype: type = np.float32,
#     reshape: bool = True,
# ) -> None:
#     """
#     Write a numpy array to a fortran binary file. Array must be cast to the
#     appropriate data type before calling this function. If the array is 3D,
#     the array will be reshaped from (y, x, z) to (z, y, x) for fortran.
#     """
#     # Reshape array from (y, x, z) to (z, y, x) (also for fortran)
#     if reshape:
#         if len(array.shape) == 3:
#             array = np.moveaxis(array, 2, 0).astype(dtype)
#         else:
#             array = array.astype(dtype)
#     else:
#         array = array.astype(dtype)

#     # Write the zarr array to a dat file with scipy FortranFile package
#     with FortranFile(Path(output_dir, dat_name), "w") as f:
#         f.write_record(array)
