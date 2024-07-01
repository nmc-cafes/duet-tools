## Calibration

DUET output files can be calibrated to match user-provided or data-derived ranges or distributions of fuel parameters. Please see [calibration](reference.md#duet_tools.calibration) for full documentation.

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

- **nx** and **ny** define the number of cells in the x and y direction of the DUET grid.
- **directory** is the path to the DUET output files.

### How to calibrate DUET outputs to target ranges and/or distributions

A target range for each fuel parameter can be defined using method="maxmin". If instead you want to define a target distribution of values, use method="meansd".
First, make [`Targets`](reference.md#duet_tools.calibration.Targets) objects for each fuel parameter and fuel type you wish to calibrate using [`assign_targets`](reference.md#duet_tools.calibration.assign_targets). Then, set each fuel parameter to the associated target(s) using [`set_fuel_parameter`](reference.md#duet_tools.calibration.set_fuel_parameter). Last, provide a list of the resulting [`FuelParameter`](reference.md#duet_tools.calibration.FuelParameter) objects to the [`calibrate`](reference.md#duet_tools.calibration.calibrate) function, along with the `DuetRun` to calibrate.

```python
from duet_tools.calibrattion import assign_targets, set_fuel_parameter, calibrate

grass_density = assign_targets(method="maxmin", max=1.0, min=0.1)
litter_density = assign_targets(method="meansd", mean=0.6, sd=0.1)
grass_height = assign_targets(method="constant", value=0.75)
litter_height = assign_targets(method="constant", value=0.2)
```

- **method** specifies how the calibration will be conducted. When using the `"maxmin"` method, a target range of values should be supplied using the keyword arguments **max** and **min**. To specify a target distribution, set the method to `"meansd"` and use the keyword arguments **mean** and **sd**. To assign the same value to everywhere a fuel type is present, use the `"constant"` calibration method with a keyword argument of **value**.

Once any number of `Targets` objects are created, they are used to set the targets of each desired fuel parameter.

```python
density_targets = set_fuel_parameter(
    parameter="density", grass=grass_density, litter=litter_density
)
height_targets = set_fuel_parameter(
    parameter="height", grass=grass_height, litter=litter_height
)
```
- **parameter** can be one of `"density"`, `"height"`, or `"moisture"`. A `FuelParameter` object represents only one of thes parameters.

- **keyword arguments** specify which fuel type(s) should be set for a given parameter. To set fuel types individually, use either or both of **grass** and **litter**. If you have targets that apply to all fuel types, rather than litter or grass separately, simply use the **all** keyword argument.

```python
all_density = assign_targets(method="maxmin", max= 1.0, min=0.1)
density_targets = set_fuel_parameter(parameter="density", all=all_density)
```

Last, use the calibrate function to return a new `DuetRun` object with calibrated fuel arrays.

```python
# Calibrate
calibrated_duet = calibrate(
    duet_run=duet_run, fuel_parameter_targets=[density_targets, grass_targets]
)
```

- **duet_run** is the `DuetRun` object that will be calibrated.
- **fuel_parameter_targets** is the `FuelParameter` object or list of `FuelParameter` objects that contain calibration targets.

### How to calibrate DUET using LANDFIRE data

When fuel parameter targets are not known for the study area, calibration can be conducted using targets derived from LANDFIRE data (Scott and Burgan 40 Fuel Models; SB40). These data can be queried using the secondary `landfire` module. Please see [landfire](reference.md#duet_tools.landfire) for full documentation.

The first step is to use the [`query_landfire`](reference.md#duet_tools.landfire.query_landfire) function to access fuels data for the specific area of interest. A spatial bounding box must be supplied in the form of either a geojson or shapely polygon. The function returns a [`LandfireQuery`](reference.md#duet_tools.landfire.LandfireQuery) class object.

```python
import geojson
from duet_tools.landfire import query_landfire

geojson_path = "path/to/geojson/file"
landfire_path = "path/where/landfire/files/are/saved"
with open(geojson_path) as fid:
    aoi_geojson = geojson.load(fid)
landfire_query = query_landfire(
    area_of_interest=aoi_geojson,
    directory=landfire_path,
    input_epsg=4326
)
```

- **area_of_interest** may be given as either a geojson polygon or a shapely polygon. It is the spatial bounding box of the burn plot/area of interest.
- **directory** is the path to the directory where landfire-associated files will be download, namely a .zip file that is read in and processed under the hood.
- **input_epsg** is the EPSG code for the coordinate reference system and projects of the area of interest polygon.
- **delete_files** specifies whether or not to delete to files downloaded from the LANDFIRE website. Since the files are usually not needed after the `LandfireQuery` object is returned, it defaults to True.

Once LANDFIRE data is queried, targets can be assigned for whatever fuel parameters and types the user desires using [`assign_targets_from_sb40`](reference.md#duet_tools.calibration.assign_targets_from_sb40). Unlike `assign_targets`, the fuel parameter and fuel type must be specified for targets to be assigned.

```python
litter_density_sb40 = assign_targets_from_sb40(
    query=landfire_query,
    fuel_type="litter",
    parameter="density",
    method="maxmin",
)
```

A `Targets` object is returned with values derived from the SB40 designations in the area of interest. An error will be issued if the given fuel type is not present in the area of interest.

The default calibration method is "maxmin". If the fuel parameter has only a single value for the given fuel type, the chosen calibration method will be changed to "constant" (if not already) and a warning will be issued. If "meansd" is used, a warning will be issued, since specifying a distribution from SB40-derived values is discouraged. If the "constant" calibration method is specified and more than one value for the given fuel parameter and type is present, an error will be issued.

Once a `Targets` object is obtained from `assing_targets_from_sb40` calibration proceeds in the normal way.

```python
density_targets = set_fuel_parameter(
    parameter="density", litter=litter_density_sb40
)
calibrated_duet = calibrate(
    duet_run=duet_run, fuel_parameter_targets=density_targets
)
```

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
    fuel_type="integrated", fuel_parameter="density"
) #2D array of grass and litter
separated_density = calibrated_duet.to_numpy(
    fuel_type="separated", fuel_parameter="height
) #3D array where grass is z-layer 0 and litter is z-layer 1
```

### How to export quicfire .dat files from calibrated DUET outputs

The calibrated `DuetRun` object can also be exported as QUIC-Fire fuel input files to a specified directory. The resulting .dat files are 3D, with a single z-layer of integrated fuel types for each parameter. Naming follows the explected QUIC-Fire filenames of "treesrhof.dat", "treesmoist.dat", and/or "treesfueldepth.dat".

```python
quicfire_path = "path/to/quicfire"
calibrated_duet.to_quicfire(directory=quicfire_path)
```

If a parameter array is not present in the `duet_run` object, it will not be written to a file. The user may also choose not to export any parameter by setting its argement to False.

```python
quicfire_path = "path/to/quicfire"
calibrated_duet.to_quicfire(directory=quicfire_path, moisture=False)
```

If the directory already contains any of those three files, an error will be raised. To overwrite those files, set overwrite to True.

```python
quicfire_path = "path/to/quicfire"
calibrated_duet.to_quicfire(directory=quicfire_path, moisture=False, overwrite=True)
```
