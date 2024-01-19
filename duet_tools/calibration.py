"""
DUET Tools Calibration module
"""
from __future__ import annotations

# Core Imports
from pathlib import Path
import importlib.resources

# from typing import ...

# External Imports
import numpy as np
from scipy.io import FortranFile
import geojson
import landfire
from landfire.geospatial import get_bbox_from_polygon
import zipfile
import rasterio as rio
import rasterio.mask
from rasterio.enums import Resampling
import pandas as pd
import re
import warnings
from pydantic import (
    BaseModel,
    NonNegativeFloat,
    NonNegativeInt,
    PositiveFloat,
    PositiveInt,
    computed_field,
    field_validator,
)

# Internal Imports


try:  # Python 3.9+
    DATA_PATH = importlib.resources.files("duet_tools").joinpath("data")
except AttributeError:  # Python 3.6-3.8
    from pkg_resources import resource_filename

    DATA_PATH = resource_filename("duet_tools", "data")


class DuetCalibrator:
    # TODO: Instead of saving just the most recent array to self.calibrated_array, append to a list or a dict of calibrated arrays, so that you can access a bunch without reading in dat files
    def __init__(self, nx, ny, nz, xmin, ymin, xmax, ymax, dx, dy, output_dir):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dy = dy
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax
        self.output_dir = Path(output_dir)
        self.sb40_data = DATA_PATH / "sb40_parameters.csv"
        self.calibrated = False
        self.calibrated_array = None
        self.calibrated_fuel_type = []
        self.calibration_method = []
        self.saved_files = []
        self.original_duet_array = self._read_original_duet()
        self.duet_dict = self._get_input_array()

    def calibrate_max_min(
        self, fuel_type: str | list, max_val: float | list, min_val: float | list
    ) -> None:
        """
        Calibrate the values of the surface bulk density output from DUET by setting the
        range (maximum and minimum).

        Parameters
        ----------

        fuel_type : str | list
            Fuel type(s) to calibrate. May be one of "total", "grass", or "litter" given
            as a string, or both "grass" and "litter" given as a list. When "total", both
            fuel types are calibrated based on one set of inputs. When "litter" or "grass",
            the given fuel type is calibrated, the other is left unchanged and added to the
            calibrated fuel type to produce the final array. When ["grass","litter"] or
            ["litter","grass"] both fuel types are calibrated based on their respective inputs.

        max_val : float | list
            Target maximim value for calibration. If multiple values are given for multiple fuel
            types, the position of each value in the list must match the position of their
            corresponding fuel type.

        min_val : float | list
            Target minimum value for calibration. If multiple values are given for multiple
            fuel types, the position of each value in the list must match the position of their
            corresponding fuel type.

        Returns
        -------
        None
            Calibrated array of DUET surface bulk density is saved to the output directory.
            Filename indicates the fuel type and the max/min calibration method, and is
            incremented if previous calibrations of the same fuel type and method have been
            conducted.

        """
        self._validate_inputs(fuel_type, max_val, min_val)
        if isinstance(fuel_type, str):
            fuel_type = [fuel_type]
        if isinstance(max_val, int) or isinstance(max_val, float):
            max_val = [max_val]
        if isinstance(min_val, int) or isinstance(min_val, float):
            min_val = [min_val]
        calibrated = {}
        for f in range(len(fuel_type)):
            arr = self.duet_dict[fuel_type[f]]
            calibrated[fuel_type[f]] = self._maxmin_calibration(
                arr, max_val[f], min_val[f]
            )
        self.calibrated_array = self._combine_fuel_types(calibrated_dict=calibrated)
        self.calibrated = True
        self.calibrated_fuel_type.append(fuel_type)
        self.calibrated_fuel_type = self._flatten(self.calibrated_fuel_type)
        self.calibration_method.append("maxmin")
        self.duet_dict = self._get_input_array()

    def calibrate_mean_sd(
        self, fuel_type: str | list, mean_val: float | list, sd_val: float | list
    ) -> None:
        """
        Calibrate the values of the surface bulk density output from DUET by setting the
        center and spread (mean and standard deviation).

        Parameters
        ----------

        fuel_type : str | list
            Fuel type(s) to calibrate. May be one of "total", "grass", or "litter" given
            as a string, or both "grass" and "litter" given as a list. When "total", both
            fuel types are calibrated based on one set of inputs. When "litter" or "grass",
            the given fuel type is calibrated, the other is left unchanged and added to the
            calibrated fuel type to produce the final array. When ["grass","litter"] or
            ["litter","grass"] both fuel types are calibrated based on their respective inputs.

        mean_val : float | list
            Target mean value for calibration. If multiple values are given for multiple fuel
            types, the position of each value in the list must match the position of their
            corresponding fuel type.

        sd_val : float | list
            Target standard deviation for calibration. If multiple values are given for multiple
            fuel types, the position of each value in the list must match the position of their
            corresponding fuel type.

        Returns
        -------
        None
            Calibrated array of DUET surface bulk density is saved to the output directory.
            Filename indicates the fuel type and the max/min calibration method, and is
            incremented if previous calibrations of the same fuel type and method have been
            conducted.

        """
        self._validate_inputs(fuel_type, mean_val, sd_val)
        if isinstance(fuel_type, str):
            fuel_type = [fuel_type]
        if isinstance(mean_val, int) or isinstance(mean_val, float):
            mean_val = [mean_val]
        if isinstance(sd_val, int) or isinstance(sd_val, float):
            sd_val = [sd_val]
        calibrated = {}
        for f in range(len(fuel_type)):
            arr = self.duet_dict[fuel_type[f]]
            calibrated[fuel_type[f]] = self._meansd_calibration(
                arr, mean_val[f], sd_val[f]
            )
        self.calibrated_array = self._combine_fuel_types(calibrated_dict=calibrated)
        self.calibrated = True
        self.calibrated_fuel_type.append(fuel_type)
        self.calibrated_fuel_type = self._flatten(self.calibrated_fuel_type)
        self.calibration_method.append("meansd")
        self.duet_dict = self._get_input_array()

    def calibrate_with_sb40(self, fuel_type: str | list) -> None:
        self._validate_inputs(fuel_type)
        print("Querying LandFire...\n")
        # Query Landfire and return array of SB40 keys
        sb40_arr = self._query_landfire()
        # Import SB40 FBFM parameters table
        sb40_params = pd.read_csv(self.sb40_data)
        # Generate dict of fastfuels bulk density values and apply to Landfire query
        sb40_dict = self._get_sb40_fuel_params(sb40_params)
        sb40_ftype, sb40_rhof = self._get_sb40_arrays(sb40_arr, sb40_dict)
        if isinstance(fuel_type, str):
            fuel_type = [fuel_type]
        calibrated = {}
        for f in fuel_type:
            if f == "grass":
                if 1 in sb40_ftype:
                    max_val = np.max(sb40_rhof[sb40_ftype == 1])
                    grass_arr = sb40_rhof[sb40_ftype == 1]
                    min_val = np.min(grass_arr[grass_arr > 0])
                else:
                    print(
                        "WARNING: grass fuel not present in sb40. Continuing with litter only."
                    )
                    max_val = 0
                    min_val = 0
            elif f == "litter":
                if -1 in sb40_ftype:
                    max_val = np.max(sb40_rhof[sb40_ftype == -1])
                    litter_arr = sb40_rhof[sb40_ftype == -1]
                    min_val = np.min(litter_arr[litter_arr > 0])
                else:
                    print(
                        "WARNING: litter fuel not present in sb40. Continuing with grass only."
                    )
                    max_val = 0
                    min_val = 0
            else:
                max_val = np.max(sb40_rhof)
                min_val = np.min(sb40_rhof[sb40_rhof > 0])
            calibrated[f] = self._maxmin_calibration(
                self.duet_dict[f], max_val, min_val
            )
        self.calibrated_array = self._combine_fuel_types(calibrated_dict=calibrated)
        self.calibrated = True
        self.calibrated_fuel_type.append(fuel_type)
        self.calibrated_fuel_type = self._flatten(self.calibrated_fuel_type)
        self.calibration_method.append("sb40")
        self.duet_dict = self._get_input_array()

    def revert_to_original_duet(self, delete_files: bool = False) -> None:
        """
        Ensure that the next calibration will be conducted on the original DUET output and
        optionally delete all files saved from previous calibrations of the DuetCalibrator
        instance.

        Parameters
        ----------

        delete_files : bool
            Whether to delete the previously saved .dat files. Default is False,
            meaning files will not be deleted and any subsequent calibrations of the
            same method and fuel type will be saved with incremented filenames.

        """
        if delete_files:
            [
                Path(self.output_dir, self.saved_files[file]).unlink()
                for file in range(len(self.saved_files))
                if Path(self.output_dir, self.saved_files[file]).exists()
            ]
            self.saved_files = []
        self.calibrated = False
        self.calibrated_array = None
        self.calibrated_fuel_type = []
        self.calibration_method = []
        self.original_duet_array = self._read_original_duet()
        self.duet_dict = self._get_input_array()

    def to_file(self) -> None:
        """
        Write the most recently calibrated surface fuel array to a .dat file.
        File will be saved to the output directory of the DuetCalibrator instance.
        """
        if self.calibrated:
            arr_name = self._name_calibrated_file()
            _write_np_array_to_dat(self.calibrated_array, arr_name, self.output_dir)
            self.saved_files.append(arr_name)
        else:
            raise Exception("Must calibrate array before writing to file.")

    def replace_quicfire_surface_fuels(self):
        """
        Replace surface fuel bulk density in quicfire output
        (from export_zarr_to_quicfire) with DUET output.

        Parameters
        ----------
        quicfire_dir: Path | str
            Directory where QUIO-Fire .dat files are located,
            and to where updated .dat files are written to.

        Returns
        -------
        None
            Modified bulk density array (treesrhof.dat) is written to the QUIC-Fire directory
        """
        with open(Path(self.output_dir, "treesrhof.dat"), "rb") as fin:
            qf_arr = (
                FortranFile(fin)
                .read_reals(dtype="float32")
                .reshape((self.nz, self.ny, self.nx), order="C")
            )
        if self.calibrated:
            tag = "calibrated"
            duet_arr = np.add(
                self.calibrated_array[0, :, :], self.calibrated_array[1, :, :]
            )
        else:
            tag = "unmodified"
            duet_arr = np.add(
                self.original_duet_array[0, :, :], self.original_duet_array[1, :, :]
            )
        qf_arr[0, :, :] = duet_arr
        _write_np_array_to_dat(
            qf_arr, "treesrhof.dat", self.output_dir, np.float32, reshape=False
        )
        print(
            "Replaced FastFuels surface fuel layer with {} DUET surface fuels".format(
                tag
            )
        )

    def _validate_inputs(self, fuel_type, val1=None, val2=None):
        # Validate fuel types
        valid_ftypes = [
            "litter",
            "grass",
            "total",
            ["litter", "grass"],
            ["grass", "litter"],
        ]
        if fuel_type not in valid_ftypes:
            raise ValueError(
                "Invalid fuel type. Must be one of {}.".format(valid_ftypes)
            )
        if fuel_type == "total" and self.calibrated == True:
            raise ValueError(
                "Invalide fuel type: 'total' fuel calibration cannot be applied to a previously calibrated array. Choose a different fuel type or use revert_to_original_duet() before calibrating total fuels."
            )
        if fuel_type in self.calibrated_fuel_type:
            warnings.warn(
                "Fuel type '{}' already calibrated. Replacing previous calibrated values.".format(
                    fuel_type
                )
            )

        # Validate fuel summary arguments
        if val1 is not None:
            if isinstance(fuel_type, list):
                if not isinstance(val1, list) or not isinstance(val2, list):
                    raise TypeError(
                        "Fuel input values must be provided as a list when fuel_type is a list."
                    )
                if len(fuel_type) != len(val1) | len(val2):
                    raise ValueError(
                        "Number of fuel value inputs must match number of fuel types ({} fuel inputs for fuel_type = {}).".format(
                            len(fuel_type), fuel_type
                        )
                    )

    def _maxmin_calibration(
        self, x: np.array, max_val: float | int, min_val: float | int
    ) -> np.array:
        """
        Scales and shifts values in a numpy array based on an observed range. Does not assume
        data is normally distributed.
        """
        x1 = x[x > 0]
        x2 = (x1 - np.min(x1)) / (np.max(x1) - np.min(x1))
        x3 = x2 * (max_val - min_val)
        x4 = x3 + min_val
        xnew = x.copy()
        xnew[np.where(x > 0)] = x4
        return xnew

    def _meansd_calibration(
        self, x: np.array, mean: float | int, sd: float | int
    ) -> np.array:
        """
        Scales and shifts values in a numpy array based on an observed mean and standard deviation.
        Assumes data is normally distributed.
        """
        x1 = x[x > 0]
        x2 = mean + (x1 - np.mean(x1)) * (sd / np.std(x1))
        xnew = x.copy()
        xnew[np.where(x > 0)] = x2
        if np.min(xnew) < 0:
            xnew = self._truncate_at_0(xnew)
        return xnew

    def _truncate_at_0(self, arr: np.array) -> np.array:
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

    def _query_landfire(self, delete_files: bool = True) -> np.array:
        """
        Download a grid of SB40 fuel models from Landfire for the unit and convert to a numpy array

        Parameters
        ----------
        delete_files: bool = True
            Whether to delete intermediate tif files. Defaults to True

        Returns
        -------
        NumPy Array
            A numpy array of the SB40 FBFM keys for the site
        """

        # Name intermediate files
        temp = [
            "landfire_sb40.zip",
            "landfire_sb40.tif",
            "sb40_upsampled.tif",
            "sb40_cropped.tif",
        ]

        # Create a bounding box from fuelgrid zarr
        coords = [
            self.xmin,
            self.ymin,
            self.xmin,
            self.ymax,
            self.xmax,
            self.ymax,
            self.xmax,
            self.ymin,
            self.xmin,
            self.ymin,
        ]
        poly = geojson.Polygon(coordinates=[coords], precision=8)
        bbox = get_bbox_from_polygon(aoi_polygon=poly, crs=5070)

        # Download Landfire data to output directory
        lf = landfire.Landfire(bbox, output_crs="5070")
        lf.request_data(
            layers=["200F40_19"], output_path=Path(self.output_dir, "landfire_sb40.zip")
        )

        # Exctract tif from compressed download folder and rename
        with zipfile.ZipFile(Path(self.output_dir, "landfire_sb40.zip")) as zf:
            extension = ".tif"
            rename = "landfire_sb40.tif"
            info = zf.infolist()
            for file in info:
                if file.filename.endswith(extension):
                    file.filename = rename
                    zf.extract(file, self.output_dir)

        # Upsample landfire raster to the quicfire resolution
        with rio.open(Path(self.output_dir, "landfire_sb40.tif")) as sb:
            upscale_factor = 30 / self.dx  # lf res/qf res
            profile = sb.profile.copy()
            # resample data to target shape
            data = sb.read(
                out_shape=(
                    sb.count,
                    int(sb.height * upscale_factor),
                    int(sb.width * upscale_factor),
                ),
                resampling=Resampling.nearest,
            )

            # scale image transform
            transform = sb.transform * sb.transform.scale(
                (sb.width / data.shape[-1]), (sb.height / data.shape[-2])
            )
            profile.update(
                {
                    "height": data.shape[-2],
                    "width": data.shape[-1],
                    "transform": transform,
                }
            )
            with rio.open(
                Path(self.output_dir, "sb40_upsampled.tif"), "w", **profile
            ) as dataset:
                dataset.write(data)

        # Crop the upsampled raster to the unit bounds
        with rio.open(Path(self.output_dir, "sb40_upsampled.tif"), "r+") as rst:
            out_image, out_transform = rasterio.mask.mask(rst, [poly], crop=True)
            out_meta = rst.meta
            out_meta.update(
                {
                    "driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform,
                }
            )
            with rio.open(
                Path(self.output_dir, "sb40_cropped.tif"), "w", **out_meta
            ) as cropped:
                cropped.write(out_image)

        # Read in the cropped raster as a numpy array
        with rio.open(Path(self.output_dir, "sb40_cropped.tif")) as rst:
            arr = rst.read(1)

        if delete_files:
            [
                Path(self.output_dir, file).unlink()
                for file in temp
                if Path(self.output_dir, file).exists()
            ]

        return arr[arr > 0]

    # TODO: fix rasterio cropping issue (grr) so that landfire raster is same size as fuelgrid

    def _get_sb40_fuel_params(self, params: pd.DataFrame = None) -> dict:
        """
        Builds a dictionary of SB40 fuel parameter values and converts them to
        the official FastFuels units

        Returns:
            dict: SB40 parameters for each fuel model
        """
        # Load in the SB40 fuel parameters
        if params is not None:
            sb40_df = params.copy()
        else:
            fpath = Path("src", "data_module", "data", "sb40_parameters.csv")
            with open(fpath) as fin:
                sb40_df = pd.read_csv(fin)

        # Convert tons/ac-ft to kg/m^3
        sb40_df["1_hr_kg_per_m3"] = sb40_df["1_hr_t_per_ac"] * 0.22417
        sb40_df["10_hr_kg_per_m3"] = sb40_df["10_hr_t_per_ac"] * 0.22417
        sb40_df["100_hr_kg_per_m3"] = sb40_df["100_hr_t_per_ac"] * 0.22417
        sb40_df["live_herb_kg_per_m3"] = sb40_df["live_herb_t_per_ac"] * 0.22417
        sb40_df["live_woody_kg_per_m3"] = sb40_df["live_woody_t_per_ac"] * 0.22417

        # Convert inverse feet to meters
        sb40_df["dead_1_hr_sav_ratio_1_per_m"] = (
            sb40_df["dead_1_hr_sav_ratio_1_per_ft"] * 3.2808
        )
        sb40_df["live_herb_sav_ratio_1_per_m"] = (
            sb40_df["live_herb_sav_ratio_1_per_ft"] * 3.2808
        )
        sb40_df["live_wood_sav_ratio_1_per_m"] = (
            sb40_df["live_wood_sav_ratio_1_per_ft"] * 3.2808
        )

        # Convert percent to ratio
        sb40_df["dead_fuel_extinction_moisture"] /= 100

        # Convert feet to meters
        sb40_df["fuel_bed_depth_m"] = sb40_df["fuel_bed_depth_ft"] * 0.3048

        # Compute wet loading
        sb40_df["wet_load"] = sb40_df["1_hr_kg_per_m3"] + sb40_df["live_herb_kg_per_m3"]

        # Compute a live herb curing factor alpha as a function of wet loading.
        # This is kind of a B.S. approach raised by Rod on a phone call with
        # Anthony on 02/28/2023. I don't like this at all, but it is a temporary
        # Fix for the BP3D team to run some simulations.
        # low_load_fuel_models = [
        sb40_df["alpha"] = [0.5 if rho > 1 else 1.0 for rho in sb40_df["wet_load"]]

        # Compute dry loading
        sb40_df["dry_herb_load"] = sb40_df["live_herb_kg_per_m3"] * sb40_df["alpha"]
        sb40_df["dry_load"] = sb40_df["1_hr_kg_per_m3"] + sb40_df["dry_herb_load"]

        # Compute SAV
        sb40_df["sav_1hr_ratio"] = sb40_df["1_hr_kg_per_m3"] / sb40_df["dry_load"]
        sb40_df["sav_1hr"] = (
            sb40_df["sav_1hr_ratio"] * sb40_df["dead_1_hr_sav_ratio_1_per_m"]
        )
        sb40_df["sav_herb_ratio"] = sb40_df["dry_herb_load"] / sb40_df["dry_load"]
        sb40_df["sav_herb"] = (
            sb40_df["sav_herb_ratio"] * sb40_df["live_herb_sav_ratio_1_per_m"]
        )
        sb40_df["sav"] = sb40_df["sav_1hr"] + sb40_df["sav_herb"]

        # Convert nan to 0
        sb40_df.fillna(0, inplace=True)

        ## NIKO ADDITIONS
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
        sb40_df["fbfm_cat"] = sb40_df["fbfm_code"].apply(
            lambda x: re.sub(pattern, "", x)
        )
        sb40_df["duet_fuel_type"] = sb40_df["fbfm_cat"].apply(
            lambda x: duet_dict.get(x)
        )
        ## END NIKO ADDITIONS

        # Build the dictionary with fuel parameters for the Scott and Burgan 40
        # fire behavior fuel models. Dict format: key ->
        # [name, loading (tons/ac), SAV (1/ft), ext. MC (percent), bed depth (ft)]
        # Note: Eventually we want to get rid of this and just use the dataframe.
        # This is legacy from the old parameter table json.
        sb40_dict = {}
        for key in sb40_df["key"]:
            row = sb40_df[sb40_df["key"] == key]
            sb40_dict[key] = [
                row["fbfm_code"].values[0],
                row["dry_load"].values[0],
                row["sav"].values[0],
                row["dead_fuel_extinction_moisture"].values[0],
                row["fuel_bed_depth_m"].values[0],
                row["duet_fuel_type"].values[0],
            ]

        return sb40_dict

    def _get_sb40_arrays(self, sb40_keys: np.array, sb40_dict: dict) -> tuple:
        """
        Use a dictionary of bulk density and fuel types that correspond to SB40
        fuel models to assign those values across the study area.

        Fuel types are as follows:
        - 1: Predominantly grass. All cells with a GR, GS, or SH designation from SB40.
        - -1: Predominantly tree litter. All cells with a TL designation from SB40.
        - 0: Neither predominantly grass or tree litter. All other SB40 designations.

        Returns:
            - np.array of fuel types
            - np.array of bulk density values as calculated by fastfuels
        """
        ftype_dict = {key: val[5] for key, val in sb40_dict.items()}
        ftype_arr = np.vectorize(ftype_dict.get)(sb40_keys)

        rhof_dict = {key: val[1] for key, val in sb40_dict.items()}
        rhof_arr = np.vectorize(rhof_dict.get)(sb40_keys)

        return ftype_arr, rhof_arr

    def _combine_fuel_types(self, calibrated_dict) -> np.array:
        calibrated_duet = np.zeros((2, self.ny, self.nx))
        if len(calibrated_dict.keys()) == 1:
            ftype = list(calibrated_dict.keys())[0]
            if ftype == "grass":
                calibrated_duet[0, :, :] = calibrated_dict["grass"]
                calibrated_duet[1, :, :] = self.duet_dict["litter"]
            elif ftype == "litter":
                calibrated_duet[0, :, :] = self.duet_dict["grass"]
                calibrated_duet[1, :, :] = calibrated_dict["litter"]
            else:
                grass_weights = self.duet_dict["grass"] / self.duet_dict["total"]
                litter_weights = self.duet_dict["litter"] / self.duet_dict["total"]
                calibrated_duet[0, :, :][
                    np.where(self.original_duet_array[0, :, :] > 0)
                ] = (calibrated_dict["total"] * grass_weights)[
                    np.where(self.original_duet_array[0, :, :] > 0)
                ]
                calibrated_duet[1, :, :][
                    np.where(self.original_duet_array[0, :, :] > 0)
                ] = (calibrated_dict["total"] * litter_weights)[
                    np.where(self.original_duet_array[0, :, :] > 0)
                ]
        else:
            calibrated_duet[0, :, :] = calibrated_dict["grass"]
            calibrated_duet[1, :, :] = calibrated_dict["litter"]
        return calibrated_duet

    def _read_original_duet(self):
        nx = self.nx
        ny = self.nx
        nz = 2  # number of duet layers, right now grass and litter. Will be number of tree species + 1
        with open(Path(self.output_dir, "surface_rhof.dat"), "rb") as fin:
            duet_rhof = (
                FortranFile(fin)
                .read_reals(dtype="float32")
                .reshape((nz, ny, nx), order="F")
            )
        return duet_rhof

    def _get_input_array(self):
        duet_dict = {}
        if self.calibrated:
            duet_dict["grass"] = self.calibrated_array[0, :, :]
            duet_dict["litter"] = self.calibrated_array[1, :, :]
            duet_dict["total"] = np.add(
                self.calibrated_array[0, :, :], self.calibrated_array[1, :, :]
            )
        else:
            duet_dict["grass"] = self.original_duet_array[0, :, :]
            duet_dict["litter"] = self.original_duet_array[1, :, :]
            duet_dict["total"] = np.add(
                self.original_duet_array[0, :, :], self.original_duet_array[1, :, :]
            )
        return duet_dict

    def _name_calibrated_file(self) -> str:
        delim = "_"
        ftype_str = (
            self.calibrated_fuel_type[0]
            if len(self.calibrated_fuel_type) == 1
            else delim.join([str(ele) for ele in self.calibrated_fuel_type])
        )
        method_str = (
            self.calibration_method[0]
            if len(self.calibration_method) == 1
            else delim.join([str(ele) for ele in self.calibration_method])
        )
        arr_name = delim.join(
            ["surface_rhof_calibrated", ftype_str, "{}.dat".format(method_str)]
        )
        if Path(self.output_dir, arr_name).exists():
            i = 2
            while Path(
                self.output_dir,
                delim.join(
                    ["surface_rhof_calibrated", ftype_str, method_str, "%s.dat" % i]
                ),
            ).exists():
                i += 1
            arr_name = delim.join(
                ["surface_rhof_calibrated", ftype_str, method_str, "%s.dat" % i]
            )
        return arr_name

    def _flatten(self, A):
        rt = []
        for i in A:
            if isinstance(i, list):
                rt.extend(self._flatten(i))
            else:
                rt.append(i)
        return rt


def _write_np_array_to_dat(
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
