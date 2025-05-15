"""
DUET Tools Calibration module
"""

from __future__ import annotations

# Core Imports
from pathlib import Path
import warnings

# External Imports
import numpy as np
import pandas as pd

# Internal Imports
from duet_tools.utils import read_dat_to_array, write_array_to_dat

DATA_DIR = Path(__file__).parent / "data"


class DuetRun:
    """
    Class containing all arrays for a DUET run.

    Attributes
    ----------
    density : np.ndarray
        3D Array of fuel bulk density (kg/m^3) values in the format exported by DUET:
        Grass bulk denisty in first layer, litter bulk density for each tree species in
        subsequent layers.
    moisture : np.ndarray
        3D Array of fuel moisture content (%) values in the format exported by DUET:
        Grass moisture content in first layer, litter moisture content for each tree
        species in subsequent layers.
    height : np.ndarray
        3D Array of fuel fuel height (m) values in the format exported by DUET:
        Grass height in first layer, litter depth for each tree species in subsequent
        layers.
    duet_version : str
        DUET version. Must be one of "v1" or "v2".

    """

    def __init__(
        self,
        density: np.ndarray,
        height: np.ndarray,
        moisture: np.ndarray,
        duet_version: str,
    ):
        self.density = density
        self.moisture = moisture
        self.height = height
        self.duet_version = duet_version

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
                treesrhof = self._integrate_all("density")
                write_array_to_dat(treesrhof, "treesrhof.dat", directory)
                written_files.append("treesrhof.dat")
        if moisture:
            if self.moisture is not None:
                treesmoist = self._integrate_all("moisture")
                write_array_to_dat(treesmoist, "treesmoist.dat", directory)
                written_files.append("treesmoist.dat")
        if height:
            if self.height is not None:
                treesfueldepth = self._integrate_all("height")
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
            Fuel type of desired array. Must be one of "integrated, "separated", "grass",
            "litter", "deciduous", or "coniferous.
            "integrated" : returns a vertically-integrated array of all fuel types.
                Array remains 3D, with shape (1,ny,nx). Integration method depends on
                fuel parameter.
            "separated" : returns a 3D array of shape (nlitter,ny,nx), where the first layer
                is grass, and the subsequent layers are litter. If using DUET v1, nlitter = 1;
                if using DUET v2, nlitter is 2 (deciduous and coniferous).
            "grass" : returns a 3D array of the chosen parameter for grass, with shape (1,ny,nx)
            "litter" : returns a 3D array of integrated litter values for all tree species,
            with shape (1,ny,nx).
            "coniferous" : returns a 3D array of litter values for coniferous tree species, with
            shape (1,ny,nx).
            "deciduous" : returns a 3D array of litter values for deciduous tree species, with
            shape (1,ny,nx).

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
            return self._integrate_all(fuel_parameter)
        if fuel_type == "grass":
            return self.__dict__[fuel_parameter][0, :, :].copy()
        if fuel_type == "litter":
            return self._integrate_litter(fuel_parameter)
        if fuel_type == "coniferous":
            return self.__dict__[fuel_parameter][1, :, :].copy()
        if fuel_type == "deciduous":
            return self.__dict__[fuel_parameter][2, :, :].copy()

    def _integrate_all(self, fuel_parameter: str) -> np.ndarray:
        if fuel_parameter == "density":
            return np.sum(self.density, axis=0)
        if fuel_parameter == "moisture":
            return _density_weighted_average(self.moisture, self.density)
        if fuel_parameter == "height":
            return np.max(self.height, axis=0)

    def _integrate_litter(self, fuel_parameter: str) -> np.ndarray:
        if fuel_parameter == "density":
            return np.sum(self.density[1:, :, :], axis=0)
        if fuel_parameter == "moisture":
            return _density_weighted_average(
                self.moisture[1:, :, :], self.density[1:, :, :]
            )
        if fuel_parameter == "height":
            return np.max(self.height[1:, :, :], axis=0)

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
        fueltypes_allowed = {
            "v1": ["grass", "litter", "separated", "integrated"],
            "v2": [
                "grass",
                "litter",
                "separated",
                "integrated",
                "deciduous",
                "coniferous",
            ],
        }
        if fuel_type not in fueltypes_allowed.get(self.duet_version):
            raise ValueError(
                f"Fuel type {fuel_type} not supported for DUET version {self.duet_version}. "
                f"Must be one of {fueltypes_allowed.get(self.duet_version)}"
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
    fuel_types : list[str]
        Fuel type(s) to which targets should be set. May be any of "grass", "litter",
        "coniferous", "deciduous", or "all".
    targets : list[Targets]
        Targets to be set to the provided parameter and fuel types.
    """

    def __init__(self, parameter: str, fuel_types: list[str], targets: list[Targets]):
        self.parameter = self._validate_fuel_parameter(parameter)
        self.fuel_types = self._validate_fuel_types(fuel_types)
        self.targets = targets

    def _validate_fuel_types(self, fuel_types):
        fueltypes_allowed = ["grass", "litter", "coniferous", "deciduous", "all"]
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


def import_duet(
    directory: str | Path, nx: int, ny: int, nsp: int, version: str = "v2"
) -> DuetRun:
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
    nsp: int
        Number of vegetation species (tree species + grass) in the DUET outputs.
        Must be 2 (grass and litter) for DUET v1.
    version: str
        DUET version that produced the outputs. Must be one of ["v1","v2"]. Defaults to "v2".

    Returns
    -------
    Instance of class DuetRun
    """
    supported = ["v1", "v2"]
    if version not in supported:
        raise ValueError(
            f"Version {version} not supported. Please use one of {supported}"
        )
    if isinstance(directory, str):
        directory = Path(directory)

    name_dict = {
        "rhof": {"v1": "surface_rhof.dat", "v2": "surface_rhof_layered.dat"},
        "depth": {"v1": "surface_depth.dat", "v2": "surface_depth_layered.dat"},
        "moist": {"v1": "surface_moist.dat", "v2": "surface_moist_layered.dat"},
    }

    density_nsp = read_dat_to_array(
        directory=directory,
        filename=name_dict["rhof"].get(version),
        nx=nx,
        ny=ny,
        nsp=nsp,
    )
    height_nsp = read_dat_to_array(
        directory=directory,
        filename=name_dict["depth"].get(version),
        nx=nx,
        ny=ny,
        nsp=nsp,
    )
    moisture_nsp = read_dat_to_array(
        directory=directory,
        filename=name_dict["moist"].get(version),
        nx=nx,
        ny=ny,
        nsp=nsp,
    )
    if version == "v1":  # arrays are kept as-is
        density = density_nsp.copy()
        height = height_nsp.copy()
        moisture = moisture_nsp.copy()
    if version == "v2":
        density = np.zeros((3, ny, nx))
        height = np.zeros((3, ny, nx))
        moisture = np.zeros((3, ny, nx))
        density[0, :, :] = density_nsp[0, :, :]
        height[0, :, :] = height_nsp[0, :, :]
        moisture[0, :, :] = moisture_nsp[0, :, :]

        # identify which species are deciduous and coniferous
        groups = _group_litter_species(directory)
        coniferous_indices = [i for i, v in groups.items() if v == "coniferous"]
        deciduous_indices = [i for i, v in groups.items() if v == "deciduous"]

        # for each parameter, coniferous is layer 1, deciduous is layer 2
        density[1, :, :] = (
            density_nsp[coniferous_indices].sum(axis=0)
            if coniferous_indices
            else np.zeros(density.shape[1:])
        )
        density[2, :, :] = (
            density_nsp[deciduous_indices].sum(axis=0)
            if deciduous_indices
            else np.zeros(density.shape[1:])
        )
        height[1, :, :] = (
            height_nsp[coniferous_indices].sum(axis=0)
            if coniferous_indices
            else np.zeros(height.shape[1:])
        )
        height[2, :, :] = (
            height_nsp[deciduous_indices].sum(axis=0)
            if deciduous_indices
            else np.zeros(height.shape[1:])
        )
        moisture[1, :, :] = (
            _density_weighted_average(
                moisture_nsp[coniferous_indices], density_nsp[coniferous_indices]
            )
            if coniferous_indices
            else np.zeros(moisture.shape[1:])
        )
        moisture[2, :, :] = (
            _density_weighted_average(
                moisture_nsp[deciduous_indices], density_nsp[deciduous_indices]
            )
            if deciduous_indices
            else np.zeros(density.shape[1:])
        )

    return DuetRun(
        density=density, height=height, moisture=moisture, duet_version=version
    )


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
            Grass calibration targets. Only the grass layer of the DUET parameter
            array will be calibrated.
        litter : Targets
            Litter calibration targets. Only the litter layer(s) of the DUET
            parameter array will be calibrated. Coniferous and deciduous litter will
            be calibrated together.
        coniferous : Targets
            Coniferous litter calibration targets. Only the coniferous litter layer of
            the DUET parameter array will be calibrated.
        deciduous : Targets
            Deciduous litter calibration targets. Only the deciduous litter layer of
            the DUET parameter array will be calibrated.
        all : Targets
            Calibration targets for all (both) fuel types. All layers of the
            DUET parameter array will be calibrated together.

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
        Grass calibration targets. Only the grass layer of the DUET bulk density
        array will be calibrated.
    litter : Targets | None
        Litter calibration targets. Only the litter layer(s) of the DUET
        bulk density array will be calibrated. Coniferous and deciduous litter will
        be calibrated together.
    coniferous : Targets | None
        Coniferous litter calibration targets. Only the coniferous litter layer of
        the DUET parameter array will be calibrated.
    deciduous : Targets | None
        Deciduous litter calibration targets. Only the deciduous litter layer of
        the DUET bulk density array will be calibrated.
    all : Targets | None
        Calibration targets for all (both) fuel types. All layers of the
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
        Grass calibration targets. Only the grass layer of the DUET moisture
        array will be calibrated.
    litter : Targets | None
        Litter calibration targets. Only the litter layer(s) of the DUET
        moisture array will be calibrated. Coniferous and deciduous litter will
        be calibrated together.
    coniferous : Targets | None
        Coniferous litter calibration targets. Only the coniferous litter layer of
        the DUET parameter array will be calibrated.
    deciduous : Targets | None
        Deciduous litter calibration targets. Only the deciduous litter layer of
        the DUET moisture array will be calibrated.
    all : Targets | None
        Calibration targets for all (both) fuel types. All layers of the
        DUET moisture array will be calibrated together.

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
        Grass calibration targets. Only the grass layer of the DUET height
        array will be calibrated.
    litter : Targets | None
        Litter calibration targets. Only the litter layer(s) of the DUET
        height array will be calibrated. Coniferous and deciduous litter will
        be calibrated together.
    coniferous : Targets | None
        Coniferous litter calibration targets. Only the coniferous litter layer of
        the DUET parameter array will be calibrated.
    deciduous : Targets | None
        Deciduous litter calibration targets. Only the deciduous litter layer of
        the DUET height array will be calibrated.
    all : Targets | None
        Calibration targets for all (both) fuel types. All layers of the
        DUET height array will be calibrated together.

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
            if np.sum(array_to_calibrate) == 0:
                raise ValueError(f"No fuels present for fuel type {fueltype}")
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
        if fueltype == "coniferous":
            return duet_run.density[1, :, :].copy()
        if fueltype == "deciduous":
            return duet_run.density[2, :, :].copy()
        if fueltype == "litter":
            return duet_run.density[1:, :, :].sum(axis=0)
        else:
            return np.sum(duet_run.density, axis=0)
    if fuelparam == "height":
        if fueltype == "grass":
            return duet_run.height[0, :, :].copy()
        if fueltype == "coniferous":
            return duet_run.height[1, :, :].copy()
        if fueltype == "deciduous":
            return duet_run.height[2, :, :].copy()
        if fueltype == "litter":
            return duet_run.height[1:, :, :].sum(axis=0)
        else:
            return np.max(duet_run.height, axis=0)
    if fuelparam == "moisture":
        if fueltype == "grass":
            return duet_run.moisture[0, :, :].copy()
        if fueltype == "coniferous":
            return duet_run.moisture[1, :, :].copy()
        if fueltype == "deciduous":
            return duet_run.moisture[2, :, :].copy()
        if fueltype == "litter":
            return _density_weighted_average(
                duet_run.moisture[1:, :, :], duet_run.density[1:, :, :]
            )
        else:
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
        duet_version=duet_run.duet_version,
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
            if fueltype == "coniferous":
                duet_to_calibrate.__dict__[param][1, :, :] = calibrated_array
            if fueltype == "deciduous":
                duet_to_calibrate.__dict__[param][2, :, :] = calibrated_array
            if fueltype == "litter":
                duet_to_calibrate.__dict__[param][1:, :, :] = _separate_2d_array(
                    calibrated_array, param, fueltype, duet_to_calibrate
                )
            if fueltype == "all":
                duet_to_calibrate.__dict__[param] = _separate_2d_array(
                    calibrated_array, param, fueltype, duet_to_calibrate
                )
    return duet_to_calibrate


def _separate_2d_array(
    calibrated: np.ndarray, param: str, fueltype: str, duet_run: DuetRun
) -> np.ndarray:
    """
    Separates a combined array into its component fuel types based on the
    fuel parameter.
    """
    if fueltype == "all":
        separated = np.zeros(duet_run.density.shape)
        weights = np.zeros(duet_run.density.shape)
    if fueltype == "litter":
        separated = np.zeros(duet_run.density[1:, :, :].shape)
        weights = np.zeros(duet_run.density[1:, :, :].shape)
    if param == "density":
        for s in range(separated.shape[0]):
            density = (
                duet_run.density[1:, :, :] if fueltype == "litter" else duet_run.density
            )
            density_sum = np.sum(density, axis=0)
            weights[s, :, :] = np.where(
                density_sum != 0, density[s, :, :] / density_sum, 0
            )
            separated[s, :, :] = calibrated * weights[s, :, :]
    if param == "moisture":
        for s in range(duet_run.density.shape[0]):
            separated[s, :, :][np.where(duet_run.moisture[s, :, :] == 0)] = 0
    if param == "height":
        for s in range(separated.shape[0]):
            height = (
                duet_run.height[1:, :, :] if fueltype == "litter" else duet_run.height
            )
            if fueltype == "all":
                height_max = np.max(height, axis=0)
                weights[s, :, :] = np.where(
                    height_max != 0, height[s, :, :] / height_max, 0
                )
            if fueltype == "litter":
                height_sum = np.sum(height, axis=0)
                weights[s, :, :] = np.where(
                    height_sum != 0, height[s, :, :] / height_sum, 0
                )
            separated[s, :, :] = calibrated * weights[s, :, :]
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
    masked = np.ma.masked_array(moisture, moisture == 0)
    averaged = np.ma.average(masked, axis=0, weights=weights)
    integrated = np.ma.filled(averaged, 0)
    return integrated


def _group_litter_species(dir: Path) -> dict:
    """
    Returns a dictionary indicating whether the tree species is deciduous
    or coniferous for each layer in the surface_*_layered.dat files.
    """
    ref_species = pd.read_csv(DATA_DIR / "REF_SPECIES.csv")
    ref_species["Group"] = ref_species["MAJOR_SPGRPCD"].apply(_classify_spgrpcd)
    lookup_dict = (
        ref_species.drop_duplicates(subset="SPCD").set_index("SPCD")["Group"].to_dict()
    )

    with open(dir / "surface_species.dat", "r") as dat:
        lines = dat.readlines()

    spcd = [int(line.strip()) for line in lines]
    spcd_dict = {}
    for i in range(len(spcd)):
        if spcd[i] is None:
            raise Exception(
                f"DUET run contains species with FIA code {spcd[i]} "
                f"which does not have a valid major species group code"
            )
        spcd_dict[i + 1] = lookup_dict.get(spcd[i])

    return spcd_dict


def _classify_spgrpcd(code):
    if code in [1, 2]:
        return "coniferous"
    elif code in [3, 4]:
        return "deciduous"
    else:
        return None
