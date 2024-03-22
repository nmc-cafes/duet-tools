## Calibration

DUET output files can be calibrated to match user-provided or data-derived ranges or distributions of fuel parameters. Please see [calibration](reference.md#duet_tools.calibration) for full documentation

### How to import DUET output files

DUET output files can be read in using the [`import_duet`](reference.md#duet_tools.calibration.import_duet) function.

```python
from duet_tools.calibration import import_duet
```

To import DUET outputs, simply specify the path to their directory, and the x and y dimensions of the domain. The resulting object is of class [`DuetRun`](reference.md#duet_tools.calibration.DuetRun).

```python
duet_path = "path/to/duet/files"
duet_run = import_duet(directory=duet_path,nx=200,ny=200)
```

### How to calibrate DUET outputs to target ranges and/or distributions

A target range for each fuel parameter can be defined using method="maxmin". If instead you want to define a target distribution of values, use method="meansd".
First, make [`Targets`](reference.md#duet_tools.calibration.Targets) objects for each fuel parameter and fuel type you wish to calibrate using [`assign_targets`](reference.md#duet_tools.calibration.assign_targets). Then, set each fuel parameter to the associated target(s) using [`set_fuel_parameter`](reference.md#duet_tools.calibration.set_fuel_parameter). Last, provide a list of the resulting [`FuelParameter`](reference.md#duet_tools.calibration.FuelParameter) objects to the [`calibrate`](reference.md#duet_tools.calibration.calibrate) function, along with the `DuetRun` to calibrate.

```python
from duet_tools.calibrattion import assign_targets, set_fuel_parameter, calibrate

# Assign targets
grass_density = assign_targets(method="maxmin", max=1.0, min=0.1)
litter_density = assign_targets(method="meansd", mean=0.6, sd=0.1)
grass_height = assign_targets(method="constant", value=0.75)
litter_height = assign_targets(method="constant", value=0.2)

# Set fuel parameters
density_targets = set_fuel_parameter(
    parameter="density", grass=grass_density, litter=litter_density
)
height_targets = set_fuel_parameter(
    parameter="height", grass=grass_height, litter=litter_height
)

# Calibrate
calibrated_duet = calibrate(
    duet_run=duet_run, fuel_parameter_targets=[density_targets, grass_targets]
)
```

If you have targets that apply to all fuel types, rather than litter or grass separately, simply use the `all` keyword argument in `set_fuel_parameter`.

```python
all_density = assign_targets(method="maxmin", max= 1.0, min=0.1)
density_targets = set_fuel_parameter(parameter="density", all=all_density)
```

### How to calibrate DUET using LANDFIRE data

### How to create numpy arrays from calibrated DUET outputs

Once a `DuetRun` object is created by the `calibrate` function, its constituent arrays can be exported to numpy ndarrays. Specify the parameter and fuel type.

```python
grass_density_array = calibrated_duet.to_numpy(
    fuel_type="grass", fuel_parameter="density"
) #2D array
litter_height_array = calibrated_duet.to_numpy(
    fuel_type="litter", fuel_parameter="height"
)
integrated_density = calibrated_duet.to_numpy(
    fuel_type="integrated", fuel_parameter="density
) #2D array of grass and litter
separated_density = calibrated_duet.to_numpy(
    fuel_type="separated", fuel_parameter="height
) #3D array where grass is z-layer 0 and litter is z-layer 1
```

### How to export quicfire .dat files from calibrated DUET outputs

The calibrated `DuetRun` object can also be exported as QUIC-Fire fuel input files to a specified directory. The resulting .dat files are 3D, with a single z-layer of integrated fuel types for each parameter. Naming follows the explected QUIC-Fire filenames of "treesrhof.dat", "treesmoist.dat", and/or "treesfueldepth.dat"

```python
quicfire_path = "path/to/quicfire"
calibrated_duet.to_quicfire(directory=quicfire_path)
```