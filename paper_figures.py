from __future__ import annotations

from pathlib import Path
from duet_tools import (
    import_duet,
    assign_targets,
    set_fuel_parameter,
    set_loading,
    calibrate,
)
from duet_tools.landfire import (
    LandfireQuery,
    _landfire_to_array,
    _get_sb40_arrays,
    _get_sb40_fuel_params,
    _delete_intermediate_files,
    assign_targets_from_sb40,
)
import numpy as np
import pandas as pd
from shutil import copyfile
import zipfile

DATA_PATH = Path(__file__).parent / "duet_tools" / "data"
duet_path = Path(__file__).parent / "figures-data"

# Calibration figure

# Import DUET outputs
duet_run = import_duet(directory=duet_path)

# Assign targets for each fuel type and fuel parameter
## method options: "maxmin", "meansd", "constant"
grass_loading = assign_targets(method="meansd", mean=0.5, sd=0.25)
deciduous_loading = assign_targets(method="meansd", mean=0.5, sd=0.1)
coniferous_loading = assign_targets(method="maxmin", max=5.0, min=0)

# Bring together fuel types for each parameter
loading_targets = set_fuel_parameter(
    parameter="loading",
    grass=grass_loading,
    deciduous=deciduous_loading,
    coniferous=coniferous_loading,
)

# Calibrate the DUET run
calibrated_duet = calibrate(duet_run=duet_run, fuel_parameter_targets=loading_targets)

# Look at individual numpy arrays
original_loading = duet_run.to_numpy(
    fuel_type="integrated", fuel_parameter="loading"
)  # 2D array
original_grass_loading = duet_run.to_numpy(
    fuel_type="grass", fuel_parameter="loading"
)  # 2D array
original_litter_loading = duet_run.to_numpy(
    fuel_type="litter", fuel_parameter="loading"
)  # 2D array
original_deciduous_loading = duet_run.to_numpy(
    fuel_type="deciduous", fuel_parameter="loading"
)  # 2D array
original_coniferous_loading = duet_run.to_numpy(
    fuel_type="coniferous", fuel_parameter="loading"
)  # 2D array

calibrated_loading = calibrated_duet.to_numpy(
    fuel_type="integrated", fuel_parameter="loading"
)  # 2D array
calibrated_grass_loading = calibrated_duet.to_numpy(
    fuel_type="grass", fuel_parameter="loading"
)  # 2D array
calibrated_litter_loading = calibrated_duet.to_numpy(
    fuel_type="litter", fuel_parameter="loading"
)  # 2D array
calibrated_deciduous_loading = calibrated_duet.to_numpy(
    fuel_type="deciduous", fuel_parameter="loading"
)  # 2D array
calibrated_coniferous_loading = calibrated_duet.to_numpy(
    fuel_type="coniferous", fuel_parameter="loading"
)  # 2D array

save_path = duet_path / "Arrays"
save_path.mkdir(exist_ok=True)

np.savetxt(save_path / "original_loading.txt", original_loading)
np.savetxt(save_path / "original_grass.txt", original_grass_loading)
np.savetxt(save_path / "original_litter.txt", original_litter_loading)
np.savetxt(save_path / "original_deciduous.txt", original_deciduous_loading)
np.savetxt(save_path / "original_coniferous.txt", original_coniferous_loading)

np.savetxt(save_path / "calibrated_loading.txt", calibrated_loading)
np.savetxt(save_path / "calibrated_grass.txt", calibrated_grass_loading)
np.savetxt(save_path / "calibrated_litter.txt", calibrated_litter_loading)
np.savetxt(save_path / "calibrated_deciduous.txt", calibrated_deciduous_loading)
np.savetxt(save_path / "calibrated_coniferous.txt", calibrated_coniferous_loading)

# Landfire figure

# Import sample Landfire data
test_path = Path(__file__).parent / "tests" / "test-data"
test_zip = test_path / "landfire_test_data.zip"
target_zip = duet_path / "landfire_sb40.zip"
copyfile(test_zip, target_zip)

# Exctract tif from compressed download folder and rename
with zipfile.ZipFile(Path(duet_path, "landfire_sb40.zip")) as zf:
    extension = ".tif"
    rename = "landfire_sb40.tif"
    info = zf.infolist()
    for file in info:
        if file.filename.endswith(extension):
            file.filename = rename
            zf.extract(file, duet_path)

landfire_arr = _landfire_to_array(duet_path)

# Import SB40 FBFM parameters table
sb40_params_path = DATA_PATH / "sb40_parameters.csv"
sb40_params = pd.read_csv(sb40_params_path)

# Generate dict of fastfuels bulk density values and apply to Landfire query
sb40_dict = _get_sb40_fuel_params(sb40_params)
sb40_arr = _get_sb40_arrays(landfire_arr, sb40_dict)

_delete_intermediate_files(duet_path)

landfire_query = LandfireQuery(
    fuel_types=sb40_arr[0, :, :],
    loading=sb40_arr[1, :, :],
    moisture=sb40_arr[2, :, :],
    depth=sb40_arr[3, :, :],
)

# Replace the above code in publication listing
# duet_run = import_duet(directory=duet_path)
# landfire_query = query_landfire(
#     area_of_interest=aoi_geojson, year=2019, directory=duet_path, input_epsg=4326
# )

grass_loading_targets = assign_targets_from_sb40(landfire_query, "grass", "loading")
litter_loading_targets = assign_targets_from_sb40(landfire_query, "litter", "loading")
landfire_loading = set_loading(
    grass=grass_loading_targets, litter=litter_loading_targets
)
calibrated_duet_landfire = calibrate(duet_run, landfire_loading)

calibrated_grass_landfire = calibrated_duet_landfire.to_numpy("grass", "loading")
calibrated_litter_landfire = calibrated_duet_landfire.to_numpy("litter", "loading")
calibrated_loading_landfire = calibrated_duet_landfire.to_numpy("integrated", "loading")

np.savetxt(save_path / "landfire_loading.txt", landfire_query.loading)

np.savetxt(save_path / "calibrated_loading_landfire.txt", calibrated_loading_landfire)
np.savetxt(save_path / "calibrated_grass_landfire.txt", calibrated_grass_landfire)
np.savetxt(save_path / "calibrated_litter_landfire.txt", calibrated_litter_landfire)
