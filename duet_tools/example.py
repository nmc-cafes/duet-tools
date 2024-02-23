from __future__ import annotations

from pathlib import Path
from duet_tools import (
    import_duet,
    assign_targets,
    set_fuel_parameter,
    calibrate,
)

# Import DUET outputs
duet_path = Path(__file__).parent.parent / "tests" / "tmp"
duet_run = import_duet(directory=duet_path, nx=100, ny=100)

# Assign targets for each fuel type and fuel parameter
## method options: "maxmin", "meansd", "constant"
grass_density = assign_targets(method="maxmin", max=1.0, min=0.1)
litter_density = assign_targets(method="meansd", mean=0.6, sd=0.05)
grass_height = assign_targets(method="constant", value=1.0)
litter_height = assign_targets(method="constant", value=0.1)

# Bring together fuel types for each parameter
density_targets = set_fuel_parameter(
    parameter="density", grass=grass_density, litter=litter_density
)
height_targets = set_fuel_parameter(
    parameter="height", grass=grass_height, litter=litter_height
)

# Calibrate the DUET run
calibrated_duet = calibrate(
    duet_run=duet_run, fuel_parameter_targets=[density_targets, height_targets]
)

# Look at individual numpy arrays
calibrated_litter_density = calibrated_duet.to_numpy(
    fuel_type="litter", fuel_parameter="density"
)  # 2D array
calibrated_height = calibrated_duet.to_numpy(
    fuel_type="separated", fuel_parameter="height"
)  # 3D array, two layers
calibrated_density = calibrated_duet.to_numpy(
    fuel_type="integrated", fuel_parameter="density"
)  # 2D array

# Export to QUIC-Fire .dat files
qf_path = Path(__file__).parent.parent / "tests" / "tmp"
calibrated_duet.to_quicfire(directory=qf_path)
