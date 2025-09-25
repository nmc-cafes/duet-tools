from __future__ import annotations

from pathlib import Path
from duet_tools import (
    import_duet,
    assign_targets,
    set_fuel_parameter,
    calibrate,
)
import numpy as np

duet_path = Path(__file__).parent.parent / "figures-data"

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
