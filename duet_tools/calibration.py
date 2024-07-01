"""
DUET Tools Calibration module
"""

from __future__ import annotations

# Core Imports
from pathlib import Path
import warnings

# External Imports
import numpy as np

# Internal Imports
from duet_tools.utils import read_dat_to_array, write_array_to_dat
from duet_tools.landfire import LandfireQuery


class DuetRun:
    """
    Class containing all arrays for a DUET run.

    Attributes
    ----------
    density : np.ndarray
        3D Array of fuel bulk density (kg/m^3) values in the format exported by DUET:
        Grass bulk denisty in first layer, litter bulk density in second layer.
    moisture : np.ndarray
        3D Array of fuel moisture content (%) values in the format exported by DUET:
        Grass moisture content in first layer, litter moisture content in second layer.
    height : np.ndarray
        3D Array of fuel fuel height (m) values in the format exported by DUET:
        Grass height in first layer, litter depth in second layer.

    """

    def __init__(self, density: np.ndarray, height: np.ndarray, moisture: np.ndarray):
        self.density = density
        self.moisture = moisture
        self.height = height

    def to_quicfire(
        self,
        directory: str | Path,
        density: bool = True,
        moisture: bool = True,
        height: bool = True,
        overwrite: bool = False,
    ) -> None:
        """
        Writes a DuetRun object to QUIC-fire fuel .dat inputs to a directory:
        treesrhof.dat, treesmoist.dat, treesfueldepth.dat

        Parameters
        ----------
        directory : str | Path
            Path to directory for writing QUIC-fire files
        density : bool
            Whether to export the bulk density array. Defaults to True.
        moisture : bool
            Whether to export the moisture content array. Defaults to True.
        height : bool
            Whether to export the fuel height (depth) arra. Defaults to True.
        overwrite : bool
            Whether to overwrite trees*.dat files already present in the directory.
            If files exist, raises an error if True, warning if False. Defaults to False.

        Returns
        -------
        None :
            Writes QUIC-Fire .dat files to the provided directory.
        """
        written_files = []
        if isinstance(directory, str):
            directory = Path(directory)

        files = ["treesrhof.dat", "treesmoist.dat", "treesfueldepth.dat"]
        to_overwrite = []
        for file in files:
            path = directory / file
            if path.exists():
                to_overwrite.append(file)
        if len(to_overwrite) > 0:
            if overwrite:
                warnings.warn(
                    f"File(s) {to_overwrite} already exist(s) and will be overwritten."
                )
            else:
                raise FileExistsError(
                    f"File(s) {to_overwrite} already exist(s) in {directory}. "
                    f"Please set overwrite = True."
                )
        if density:
            if self.density is not None:
                treesrhof = self._integrate("density")
                write_array_to_dat(treesrhof, "treesrhof.dat", directory)
                written_files.append("treesrhof.dat")
        if moisture:
            if self.moisture is not None:
                treesmoist = self._integrate("moisture")
                write_array_to_dat(treesmoist, "treesmoist.dat", directory)
                written_files.append("treesmoist.dat")
        if height:
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
    Should be instantiated using [`assign_targets`](reference.md#duet_tools.calibration.assign_targets).

    Attributes
    ----------
    method : str
        Method by which to calibrate to the target values. Must be one of
        "maxmin", "meansd", or "constant".
    args : list[str]
        Sting(s) to be used as keyword arguments for calibration, which correspond
        to the calibration method. For maxmin calibration, use ["max","min"];
        for meansd calibration, use ["mean","sd"]; for constant calibration, use
        []"value"].
    targets : list
        Calibration targets, which correspond to the elements of Targets.args.
    """

    def __init__(self, method: str, args: list[str], targets: list):
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

    def _validate_target_args(self, method: str, args: list[str], targets: list[float]):
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
    Class containing and validating calibration targets for a single fuel parameter.
    A single Target object can be set for multiple fuel types. Should be instantiated using
    [`set_fuel_parameter`](reference.md#duet_tools.calibration.set_fuel_parameter)

    Attributes
    ----------
    parameter : str
        Fuel parameter for which targets should be set. Must be one of "density", "moisture", or "height".
    fuel_type : str
        Fuel type to which targets should be set. Must be one of "grass", "litter", or "all".
    targets : list[Targets]
        Targets to be set to the provided parameter and fuel type.
    """

    def __init__(self, parameter: str, fuel_types: list[str], targets: list[Targets]):
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
    moisture = read_dat_to_array(
        directory=directory, filename="surface_moist.dat", nx=nx, ny=ny, nz=nz
    )
    return DuetRun(density=density, height=height, moisture=moisture)


def assign_targets(method: str, **kwargs: float) -> Targets:
    """
    Assigns target values and calculation method for exactly one fuel type and parameter

    Parameters
    ----------
    method : str
        Calibration method for the target values provided. Must be one of:
        "constant", "maxmin", "meansd", "sb40".
    **kwargs : float
        Keyword arguments correspond to the calibration method.
        For "maxmin" method, **kwargs keys must be `max` and `min`.
        For "meansd" method, **kwargs keys must be `mean` and `sd`.
        For "constant" method, **kwargs key must be `value`.

    Returns
    -------
    Instance of class Targets
    """
    args = list(kwargs.keys())
    targets = list(kwargs.values())

    return Targets(method=method, args=args, targets=targets)


def assign_targets_from_sb40(
    query: LandfireQuery, fuel_type: str, parameter: str, method: str = "maxmin"
) -> Targets:
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
    query._validate_get_targets(fuel_type, parameter, method)
    # select fuel parameter
    if parameter == "density":
        param_arr = query.density
    elif parameter == "moisture":
        param_arr = query.moisture
    else:
        param_arr = query.height
    # select fuel type
    if fuel_type == "grass":
        fuel_arr = param_arr[query._get_fueltype_indices(query.fuel_types, 1)]
    elif fuel_type == "litter":
        fuel_arr = param_arr[query._get_fueltype_indices(query.fuel_types, -1)]
    else:
        fuel_arr = param_arr
    # get targets based on method
    fuel_arr = fuel_arr[np.where(fuel_arr > 0)]
    if method == "maxmin":
        if np.max(fuel_arr) == np.min(fuel_arr):
            warnings.warn(
                f"There is only one value for {fuel_type} {parameter}. "
                "Setting calibration method to 'constant'",
                UserWarning,
            )
            method = "constant"
            args = ["value"]
            targets = [np.mean(fuel_arr)]
        else:
            method = "maxmin"
            args = ["max", "min"]
            targets = [np.max(fuel_arr), np.min(fuel_arr)]
        return Targets(
            method=method,
            args=args,
            targets=targets,
        )
    if method == "meansd":
        if np.max(fuel_arr) == np.min(fuel_arr):
            warnings.warn(
                f"There is only one value for {fuel_type} {parameter}. "
                "Setting calibration method to 'constant'",
                UserWarning,
            )
            method = "constant"
            args = ["value"]
            targets = [np.mean(fuel_arr)]
        else:
            method = "meansd"
            args = ["mean", "sd"]
            targets = [np.mean(fuel_arr), np.std(fuel_arr)]
        return Targets(
            method=method,
            args=args,
            targets=targets,
        )
    if method == "constant":
        if np.max(fuel_arr) != np.min(fuel_arr):
            raise ValueError(
                "Multiple values present in Landfire query. Please use either maxmin "
                "or meansd calibration method."
            )
        return Targets(method="constant", args=["value"], targets=[np.mean(fuel_arr)])


def set_fuel_parameter(parameter: str, **kwargs: Targets):
    """
    Sets calibration targets for grass, litter, both separately, or all
    fuel types together, for a single fuel parameter.

    Parameters
    ----------
    parameter : str
        Fuel parameter for which to set targets
    **kwargs : Targets
        grass : Targets
            Grass calibration targets. Only the grass layer of the DUET bulk
            density array will be calibrated.
        litter : Targets
            Litter calibration targets. Only the litter layer of the DUET bulk
            density array will be calibrated.
        all : Targets
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
    duet_run: DuetRun, fuel_parameter_targets: list[FuelParameter] | FuelParameter
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
    """
    Writes a numpy array to a QUIC-Fire fuel input (.dat) in the chosen directory.

    Parameters
    ----------
    array : np.ndarray
        The numpy array to be written. Must be 3D.
    directory : str | Path
        The directory where the file will be written.
    filename : str
        The name of the file to be written. Must end in ".dat".

    Returns
    -------
    None
        File is written to disk.
    """
    if isinstance(directory, str):
        directory = Path(directory)
    write_array_to_dat(array=array, dat_name=filename, output_dir=directory)


def _get_array_to_calibrate(duet_run: DuetRun, fueltype: str, fuelparam: str):
    """
    Identifies and returns the layer of the duet outputs should be processed based
    on the fuel parameter and fuel type.
    """
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
        if fueltype == "grass":
            return duet_run.moisture[0, :, :].copy()
        if fueltype == "litter":
            return duet_run.moisture[1, :, :].copy()
        density_weights = duet_run.density.copy()
        density_weights[density_weights == 0] = 0.01
        return np.average(duet_run.moisture, weights=density_weights, axis=0)


def _duplicate_duet_run(duet_run: DuetRun) -> DuetRun:
    """
    Makes a copy of a DuetRun object.
    """
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
    """
    Calibrates an array based on the method and values in a Targets object.
    """
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
    """
    Conducts constant calibration by changing all nonzero values in a
    duet fuel parameter/type to a single given value.
    """
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
    """
    Replaces or creates calibrated array(s) in a DuetRun object.
    """
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
    """
    Separates a combined array into its component fuel types based on the
    fuel parameter.
    """
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
