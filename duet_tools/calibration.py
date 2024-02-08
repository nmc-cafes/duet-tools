"""
DUET Tools Calibration module
"""

from __future__ import annotations

# Core Imports
from pathlib import Path
import importlib.resources
import warnings
import re

# External Imports
import numpy as np
import pandas as pd
import geojson
import shapely
import zipfile
from pyproj import Transformer
import landfire
from landfire.geospatial import get_bbox_from_polygon
import rasterio as rio
from rasterio.enums import Resampling
from rasterio.mask import mask

# Internal Imports
from duet_tools.utils import read_dat_to_array, write_array_to_dat


try:  # Python 3.9+
    DATA_PATH = importlib.resources.files("duet_tools").joinpath("data")
except AttributeError:  # Python 3.6-3.8
    from pkg_resources import resource_filename

    DATA_PATH = resource_filename("duet_tools", "data")


class DuetRun:
    """
    Class containing all arrays for a DUET run.
    """

    def __init__(
        self, density: np.ndarray, height: np.ndarray, moisture: np.ndarray = None
    ):
        self.density = density
        self.moisture = moisture
        self.height = height

    def add_moisture_array(self, moisture_array=np.ndarray) -> None:
        """
        Add an array of moisture values to the DUET run

        Parameters
        ----------
        moisture_array : np.ndarray
            3D array of fuel moisture content values. Must be the same shape as
            density and height arrays already present in the object. Values must
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
        written_files = []
        if isinstance(directory, str):
            directory = Path(directory)
        if self.density is not None:
            treesrhof = self._integrate("density")
            write_array_to_dat(treesrhof, "treesrhof.dat", directory)
            written_files.append("treesrhof.dat")
        if self.moisture is not None:
            treesmoist = self._integrate("moisture")
            write_array_to_dat(treesmoist, "treesmoist.dat", directory)
            written_files.append("treesmoist.dat")
        if self.height is not None:
            treesfueldepth = self._integrate("height")
            write_array_to_dat(treesfueldepth, "treesfueldepth.dat", directory)
            written_files.append("treesfueldepth.dat")
        if len(written_files) == 0:
            print("No files were written")
        else:
            print(
                f"QUIC-Fire files {written_files} were written to directory {directory}"
            )

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
            "height".

        Returns
        -------
        np.ndarray :
            Numpy array of the provided fuel type and parameter.
        """
        self._validate_fuel_inputs(fuel_type, fuel_parameter)
        if fuel_type == "separated":
            return self.__dict__[fuel_parameter].copy()
        if fuel_type == "integrated":
            return self._integrate(fuel_parameter)
        if fuel_type == "grass":
            return self.__dict__[fuel_parameter][0, :, :].copy()
        if fuel_type == "litter":
            return self.__dict__[fuel_parameter][1, :, :].copy()

    def _integrate(self, fuel_parameter: str) -> np.ndarray:
        if fuel_parameter == "density":
            return np.sum(self.density, axis=0)
        if fuel_parameter == "moisture":
            return _density_weighted_average(self.moisture, self.density)
        if fuel_parameter == "height":
            return np.max(self.height, axis=0)

    def _validate_input_moisture(self, moisture: np.ndarray):
        if moisture.shape != self.density.shape:
            raise ValueError(
                f"Input array shape {moisture.shape} must match existing arrays {self.density.shape}."
            )
        if self.density[np.where(moisture == 0)].any() != 0:
            raise ValueError(
                "Value of moisture array cannot be zero where fuel is present"
            )

    def _validate_fuel_inputs(self, fuel_type: str, fuel_parameter: str):
        fueltypes_allowed = ["grass", "litter", "separated", "integrated"]
        if fuel_type not in fueltypes_allowed:
            raise ValueError(
                f"Fuel type {fuel_type} not supported. Must be one of {fueltypes_allowed}"
            )
        parameters_allowed = ["density", "moisture", "height"]
        if fuel_parameter not in parameters_allowed:
            raise ValueError(
                f"Fuel parameter {fuel_parameter} not supported. Must be one of {parameters_allowed}"
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
        methods_allowed = ["maxmin", "meansd", "constant"]
        if method not in methods_allowed:
            raise ValueError(
                f"Method {method} not supported. Must be one of {methods_allowed}"
            )
        return method

    def _validate_target_args(self, method: str, args: list(str), targets: list(float)):
        method_dict = {
            "constant": ["value"],
            "maxmin": ["max", "min"],
            "meansd": ["mean", "sd"],
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
        fuel_parameters_allowed = ["density", "moisture", "height"]
        if parameter not in fuel_parameters_allowed:
            raise ValueError(
                f"Fuel parameter {parameter} not supported. Must be one of {fuel_parameters_allowed}"
            )
        return parameter


class LandfireQuery:
    """
    Class containing the information from a LandFire query, to be passed to assign_targets()
    """

    def __init__(
        self,
        fuel_types: np.ndarray,
        density: np.ndarray,
        moisture: np.ndarray,
        height: np.ndarray,
    ):
        self.fuel_types = fuel_types
        self.density = density
        self.moisture = moisture
        self.height = height
        self._validate_arrays_shape()

    def get_targets(self, fuel_type: str, parameter: str, method: str):
        self._validate_get_targets(fuel_type, parameter, method)
        # select fuel parameter
        if parameter == "density":
            param_arr = self.density
        elif parameter == "moisture":
            param_arr = self.moisture
        else:
            param_arr = self.height
        # select fuel type
        if fuel_type == "grass":
            fuel_arr = param_arr[self._get_fueltype_indices(self.fuel_types, 1)]
        elif fuel_type == "litter":
            fuel_arr = param_arr[self._get_fueltype_indices(self.fuel_types, -1)]
        else:
            fuel_arr = param_arr
        # get targets based on method
        fuel_arr = fuel_arr[np.where(fuel_arr > 0)]
        if method == "maxmin":
            if np.max(fuel_arr) == np.min(fuel_arr):
                raise ValueError(
                    f"There is only one value for {fuel_type} {parameter}. "
                    "Please use 'constant' calibration method"
                )
            return Targets(
                method="maxmin",
                args=["max", "min"],
                targets=[np.max(fuel_arr), np.min(fuel_arr)],
            )
        if method == "meansd":
            if np.max(fuel_arr) == np.min(fuel_arr):
                raise ValueError(
                    f"There is only one value for {fuel_type} {parameter}. "
                    "Please use 'constant' calibration method"
                )
            return Targets(
                method="meansd",
                args=["mean", "sd"],
                targets=[np.mean(fuel_arr), np.std(fuel_arr)],
            )
        if method == "constant":
            if np.max(fuel_arr) != np.min(fuel_arr):
                raise ValueError(
                    "Multiple values present in Landfire query. Please use either maxmin "
                    "or meansd calibration method."
                )
            return Targets(
                method="constant", args=["value"], targets=[np.mean(fuel_arr)]
            )

    def _get_fueltype_indices(self, arr: np.ndarray, ft: int):
        ft_dict = {1: "grass", -1: "litter"}
        if ft not in arr.flatten():
            raise ValueError(f"Landfire query does not contain any {ft_dict[ft]} fuel.")
        return np.where(arr == ft)

    def _validate_arrays_shape(self):
        attributes = [
            getattr(self, attr)
            for attr in dir(self)
            if not callable(getattr(self, attr)) and not attr.startswith("__")
        ]
        first_shape = None
        for arr in attributes:
            if isinstance(arr, np.ndarray):
                if first_shape is None:
                    first_shape = arr.shape
                elif arr.shape != first_shape:
                    raise ValueError(
                        "All arrays in LandfireQuery must have the same shape."
                    )

    def _validate_get_targets(self, fuel_type, parameter, method):
        fueltypes_allowed = ["grass", "litter", "all"]
        parameters_alowed = ["density", "moisture", "height"]
        methods_allowed = ["maxmin", "meansd", "constant"]

        if fuel_type not in fueltypes_allowed:
            raise ValueError(
                f"Fuel type {fuel_type} not supported. Must be one of {fueltypes_allowed}."
            )
        if parameter not in parameters_alowed:
            raise ValueError(
                f"Parameter {parameter} not supported. Must be one of {parameters_alowed}."
            )
        if method not in methods_allowed:
            raise ValueError(
                f"Method {method} not supported. Must be one of {methods_allowed}."
            )


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
    height = read_dat_to_array(
        directory=directory, filename="surface_depth.dat", nx=nx, ny=ny, nz=nz
    )
    return DuetRun(density=density, height=height)


def query_landfire(
    area_of_interest: geojson.Polygon | shapely.Polygon,
    directory: str | Path,
    input_epsg: int,
    delete_files: bool = True,
) -> LandfireQuery:
    """
    Creates and submits a LANDFIRE query for a specified area of interest.

    Parameters
    ----------
    area_of_interest : geojson.Polygon | shapely.Polygon
        Area in which to query LANDFIRE data. For best results, dimensions in meters should
        match (nx*dx, ny*dy) of DUET domain.
    directory : Path | str
        Directory where files associated with the LANDFIRE query will be saved.
    input_epsg : int
        EPSG number for CRS of area_of_interest polyong
    delete_files : bool = True
        Whether to delete intermediate files created in the process of querying LANDFIRE data. Defaults to True

    Returns
    -------
    LandfireQuery
    """
    if isinstance(directory, str):
        directory = Path(directory)

    if isinstance(area_of_interest, geojson.Polygon):
        area_of_interest = shapely.Polygon(area_of_interest["coordinates"][0])

    if input_epsg != 4236:
        area_of_interest = _reproject_polygon(
            area_of_interest, input_epsg, target_epsg=4326
        )

    _query_landfire(poly=area_of_interest, directory=directory)
    landfire_arr = _landfire_to_array(directory)

    # Import SB40 FBFM parameters table
    sb40_params_path = DATA_PATH / "sb40_parameters.csv"
    sb40_params = pd.read_csv(sb40_params_path)

    # Generate dict of fastfuels bulk density values and apply to Landfire query
    sb40_dict = _get_sb40_fuel_params(sb40_params)
    sb40_arr = _get_sb40_arrays(landfire_arr, sb40_dict)

    if delete_files:
        _delete_intermediate_files(directory)

    return LandfireQuery(
        fuel_types=sb40_arr[0, :, :],
        density=sb40_arr[1, :, :],
        moisture=sb40_arr[2, :, :],
        height=sb40_arr[3, :, :],
    )


def assign_targets(method: str, **kwargs: float) -> Targets:
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


def assign_targets_from_sb40(
    query: LandfireQuery, fuel_type: str, parameter: str, method: str = "maxmin"
):
    """
    Assign a calibration target and method for a given fuel type and parameter.

    Parameters
    ----------
    query : LandfireQuery
        An object of class LandfireQuery created with query_landfire. Calibration targets
        will be calculated from these values.
    fuel_type : str
        The fuel type to obtain target values for. Must be one of "grass", "litter", or "all".
    parameter : str
        The fuel parameter to obtain target values for. Must be one of "density", "moisture", or "height".
    method : str
        The desired calibration method for the sb40-derived targets. Must be one of "maxmin", "meandsd",
        or "constant". Default is "maxmin". "constant" is only recommended if only one parameter value
        is present for the given fuel type. "meansd" is not recommended since values often do not follow a
        normal distribution.

    Returns
    -------
    Targets :
        A Targets object with values derived from Landfire and SB40 fuel models
    """
    return query.get_targets(fuel_type, parameter, method)


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


def set_height(**kwargs: Targets):
    """
    Sets height calibration targets for grass, litter, both separately, or all
    fuel types together.

    Parameters
    ----------
    grass : Targets | None
        Grass height calibration targets. Only the grass layer of the DUET bulk
        density array will be calibrated.
    litter : Targets | None
        Litter height calibration targets. Only the litter layer of the DUET bulk
        density array will be calibrated.
    all : Targets | None
        Height calibration targets for all (both) fuel types. Both layers of the
        DUET bulk density array will be calibrated together.

    Returns
    -------
    FuelParameter :
        Object representing height targets for each provided fuel type
    """
    parameter = "height"
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
    if fuelparam == "height":
        if fueltype == "grass":
            return duet_run.height[0, :, :].copy()
        if fueltype == "litter":
            return duet_run.height[1, :, :].copy()
        return np.max(duet_run.height, axis=0)
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
    new_height = duet_run.height.copy() if duet_run.height is not None else None

    new_duet = DuetRun(
        density=new_density,
        moisture=new_moisture,
        height=new_height,
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
    value = kwargs["value"]
    arr = x.copy()
    arr[arr > 0] = value
    return arr


def _add_calibrated_array(
    duet_to_calibrate: DuetRun,
    calibrated_array: np.ndarray,
    fueltype: str,
    fuelparam: str,
) -> DuetRun:
    for param in ["density", "moisture", "height"]:
        if fuelparam == param:
            if fueltype == "grass":
                duet_to_calibrate.__dict__[param][0, :, :] = calibrated_array
            if fueltype == "litter":
                duet_to_calibrate.__dict__[param][1, :, :] = calibrated_array
            if fueltype == "all":
                duet_to_calibrate.__dict__[param] = _separate_2d_array(
                    calibrated_array, param, duet_to_calibrate
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
    if param == "height":
        weights = duet_run.height.copy()
        weights[0, :, :] = duet_run.height[0, :, :] / np.max(duet_run.height, axis=0)
        weights[1, :, :] = duet_run.height[1, :, :] / np.max(duet_run.height, axis=0)
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


def _density_weighted_average(moisture: np.ndarray, density: np.ndarray) -> np.ndarray:
    """
    Vertically integrate moisture by a weighted mean, where the weights comd from cell bulk density
    """
    weights = _maxmin_calibration(density, max=1.0, min=0)
    weights[weights == 0] = 0.01
    integrated = np.average(moisture, axis=0, weights=weights)
    return integrated


def _query_landfire(
    poly: shapely.Polygon,
    directory: Path,
) -> None:
    """
    Download a grid of SB40 fuel models from Landfire for the unit and convert to a numpy array
    """

    bbox = get_bbox_from_polygon(aoi_polygon=poly, crs=4326)

    # Download Landfire data to output directory
    lf = landfire.Landfire(bbox, output_crs="5070")
    lf.request_data(
        layers=["200F40_19"], output_path=Path(directory, "landfire_sb40.zip")
    )

    # Exctract tif from compressed download folder and rename
    with zipfile.ZipFile(Path(directory, "landfire_sb40.zip")) as zf:
        extension = ".tif"
        rename = "landfire_sb40.tif"
        info = zf.infolist()
        for file in info:
            if file.filename.endswith(extension):
                file.filename = rename
                zf.extract(file, directory)


def _landfire_to_array(
    directory: Path,
) -> np.ndarray:
    # Upsample landfire raster to the quicfire resolution
    with rio.open(Path(directory, "landfire_sb40.tif")) as sb:
        arr = sb.read(1)

    return arr


def _get_sb40_fuel_params(params: pd.DataFrame) -> dict:
    """
    Builds a dictionary of SB40 fuel parameter values and converts them to
    the official FastFuels units

    Returns:
        dict: SB40 parameters for each fuel model
    """

    # Convert tons/ac-ft to kg/m^3
    params["1_hr_kg_per_m3"] = params["1_hr_t_per_ac"] * 0.22417
    params["10_hr_kg_per_m3"] = params["10_hr_t_per_ac"] * 0.22417
    params["100_hr_kg_per_m3"] = params["100_hr_t_per_ac"] * 0.22417
    params["live_herb_kg_per_m3"] = params["live_herb_t_per_ac"] * 0.22417
    params["live_woody_kg_per_m3"] = params["live_woody_t_per_ac"] * 0.22417

    # Convert inverse feet to meters
    params["dead_1_hr_sav_ratio_1_per_m"] = (
        params["dead_1_hr_sav_ratio_1_per_ft"] * 3.2808
    )
    params["live_herb_sav_ratio_1_per_m"] = (
        params["live_herb_sav_ratio_1_per_ft"] * 3.2808
    )
    params["live_wood_sav_ratio_1_per_m"] = (
        params["live_wood_sav_ratio_1_per_ft"] * 3.2808
    )

    # Convert percent to ratio
    params["dead_fuel_extinction_moisture"] /= 100

    # Convert feet to meters
    params["fuel_bed_depth_m"] = params["fuel_bed_depth_ft"] * 0.3048

    # Compute wet loading
    params["wet_load"] = params["1_hr_kg_per_m3"] + params["live_herb_kg_per_m3"]

    # Compute a live herb curing factor alpha as a function of wet loading.
    # This is kind of a B.S. approach raised by Rod on a phone call with
    # Anthony on 02/28/2023. I don't like this at all, but it is a temporary
    # Fix for the BP3D team to run some simulations.
    # low_load_fuel_models = [
    params["alpha"] = [0.5 if rho > 1 else 1.0 for rho in params["wet_load"]]

    # Compute dry loading
    params["dry_herb_load"] = params["live_herb_kg_per_m3"] * params["alpha"]
    params["dry_load"] = params["1_hr_kg_per_m3"] + params["dry_herb_load"]

    # Compute SAV
    params["sav_1hr_ratio"] = params["1_hr_kg_per_m3"] / params["dry_load"]
    params["sav_1hr"] = params["sav_1hr_ratio"] * params["dead_1_hr_sav_ratio_1_per_m"]
    params["sav_herb_ratio"] = params["dry_herb_load"] / params["dry_load"]
    params["sav_herb"] = (
        params["sav_herb_ratio"] * params["live_herb_sav_ratio_1_per_m"]
    )
    params["sav"] = params["sav_1hr"] + params["sav_herb"]

    # Convert nan to 0
    params.fillna(0, inplace=True)

    # Create dictionary for assigning fuel types for DUET calibration
    duet_dict = {
        "NB": 0,  # 0 = NEUTRAL, i.e. not predominantly grass or litter
        "GR": 1,  # 1 = GRASS predominantly
        "GS": 1,
        "SH": 1,  # I am considering shrubs as grass
        "TU": 0,
        "TL": -1,  # -1 = LITTER predominantly
        "SB": 0,
    }

    # Add column to df with DUET designations
    pattern = r"[0-9]"  # take out numbers from fbfm_type strings
    params["fbfm_cat"] = params["fbfm_code"].apply(lambda x: re.sub(pattern, "", x))
    params["duet_fuel_type"] = params["fbfm_cat"].apply(lambda x: duet_dict.get(x))

    # Build the dictionary with fuel parameters for the Scott and Burgan 40
    # fire behavior fuel models. Dict format: key ->
    # [name, loading (tons/ac), SAV (1/ft), ext. MC (percent), bed depth (ft)]
    # Note: Eventually we want to get rid of this and just use the dataframe.
    # This is legacy from the old parameter table json.
    sb40_dict = {}
    for key in params["key"]:
        row = params[params["key"] == key]
        sb40_dict[key] = [
            row["fbfm_code"].values[0],
            row["dry_load"].values[0],
            row["sav"].values[0],
            row["dead_fuel_extinction_moisture"].values[0],
            row["fuel_bed_depth_m"].values[0],
            row["duet_fuel_type"].values[0],
        ]

    return sb40_dict


def _get_sb40_arrays(sb40_keys: np.ndarray, sb40_dict: dict) -> np.ndarray:
    """
    Use a dictionary of bulk density and fuel types that correspond to SB40
    fuel models to assign those values across the study area.

    Fuel types are as follows:
    - 1: Predominantly grass. All cells with a GR, GS, or SH designation from SB40.
    - -1: Predominantly tree litter. All cells with a TL designation from SB40.
    - 0: Neither predominantly grass or tree litter. All other SB40 designations.

    Returns:
    3D np.ndarray:
    4 layers:
        1. fuel types
        2. bulk density values as calculated by fastfuels
        3. fuel moisture content values as calculated by fastfuels
        4. fuel height values as caluclated by fastfuels
    """
    val_idx = [5, 1, 3, 4]
    fuel_arr = np.zeros((4, sb40_keys.shape[0], sb40_keys.shape[1]))
    for i in range(len(val_idx)):
        layer_dict = {key: val[val_idx[i]] for key, val in sb40_dict.items()}
        layer = np.vectorize(layer_dict.get)(sb40_keys)
        fuel_arr[i, :, :] = layer

    return fuel_arr


def _reproject_polygon(poly: shapely.Polygon, input_epsg: int, target_epsg: int):
    coords = list(poly.exterior.coords)
    transformer = Transformer.from_crs(input_epsg, target_epsg, always_xy=True)
    transformed_coords = transformer.transform(*zip(*coords))
    transformed_poly = shapely.Polygon(list(zip(*transformed_coords)))

    return transformed_poly


def _delete_intermediate_files(directory: Path):
    # Name intermediate files
    temp = [
        "landfire_sb40.zip",
        "landfire_sb40.tif",
        "sb40_upsampled.tif",
        "sb40_cropped.tif",
    ]
    [Path(directory, file).unlink() for file in temp if Path(directory, file).exists()]
