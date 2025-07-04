## Inputs Module

The [`inputs`](reference.md#duet_tools.inputs) module in `duet-tools` offers a convenient way to programmatically write a DUET input file. Other inputs, such as canopy and species grids, must be supplied by the user. A DUET input file takes information about the 3D domain size (**nx, ny, nz**) and resolution (**dx, dy, dz**), the wind information (**wind_direction, wind_variability**), and the simulation parameters (**duration, random_seed**).

## Calibration Module

The [`calibration`](reference.md#duet_tools.calibration) module in `duet-tools` is a central part of the package that allows the user to scale the numeric values in DUET output files to match user-supplied target values, without altering the spatial distributions calculated in DUET. Calibrated DUET files (and DUET files in general) are meant to be used by QUIC-Fire or FIRETEC fire behavior models, and export methods for those models are included.

### Working with DUET

The `duet-tools` package assumes that the user has a license and acces to the DUET program, developed by Los Alamos National Lab. DUET distributes surface fuels based on the locations of tree canopies, the prevailing wind direction, and the time since fire. The user is responsible for parameterizing and running DUET; all calibibration in `duet-tools` is done post-hoc. By default, DUET will output two files called *surface_rhof.dat* and *surface_depth.dat*, which are Fortran-formatted arrays representing the surface fuel loading (bulk density) and surface fuel depth. They are 3D arrays where the first z layer (axis 0) is represents grass loading/moisture/depth, and the next two layers represent coniferous and deciduous tree litter, respectively.

The [`import_duet`](reference.md#duet_tools.calibration.import_duet) function assumes that DUET input and output files are present in the given directory, and that they have not been renamed or altered. Alternatively, the user may import duet arrays manually by providing file names and array sizes for the necessary DUET outputs using [`import_duet_manual`](reference.md#duet_tools.calibration.import_duet_manual). Once DUET output files are imported, the resulting [`DuetRun`](reference.md#duet_tools.calibration.DuetRun) object stores and orgainzes the arrays for all subsequent manipulations.

### Calibration Methods and Targets

The goal of the `calibration` module revolves around scaling the magnitude of DUET values without changing their spatial distribution. This assumes that the user has a target for the *range* or *distribution* of the values for each fuel type and/or parameter. Depending on ther user's knowledge of the target fuels, they may choose the appropriate **method**, specified in the [`assign_targets`](reference.md#duet_tools.calibration.assign_targets) function.

#### Methods
- **maxmin:** If the user has a target *range* of values, they can specify the maximum and minimum, and the existing DUET values will be scaled to that range without changing the underlying distribution.
- **meansd:** If the user has a target *distribution* of values, they can specify the mean and standard deviation. The distribution of DUET values will be scaled by shifting the mean and altering the spread. This method assumes that both the DUET and target values are normally distributed. Resulting values are truncated at 0.
- **constant:** This method may be used to change all nonzero values of a fuel type/parameter to a single specifed value. This affects all nonzero values of the fuel type/parameter, so the spatial distribution is preserved.

Specifying methods and target values in `assign_targets` creates an instance of class [`Targets`](reference.md#duet_tools.calibration.Targets), which stores and validates the provided information.

### Fuel Types and Parameters

Fuels are described in QUIC-Fire and FIRETEC using different *parameters* that affect fire behavior. The main parameters are fuel *loading*, fuel *moisture*, and fuelbed *depth* of the surface fuels. DUET distributes fuel loading and depth. DUET version 1 does not output fuel moisture.

These parameters are described for two *fuel types* in DUET: leaf/needle litter, which is stochastically distributed near tree canopies, and grass, which is more likely to "grow" in more open areas between tree canopies. DUET outputs separate arrays for each fuel parameter, with layers representing grass and litter fuel types.

`duet-tools` can calibrate any combination of fuel type and parameter with methods and targets contained in a `Targets` class object. The `Targets` class is agnostic of both fuel type and fuel parameter, and can thus be assigned to any number of fuel types and parameters. The [`FuelParameter`](reference.md#duet_tools.calibration.FuelParameter) class represents a single fuel parameter (*e.g.* fuel loading), storing and validating all assigned `Targets`.

A `FuelParameter` class instance is created using the [`set_fuel_parameter`](reference.md#duet_tools.calibration.set_fuel_parameter) function. The parameter to be calibrated is specified, and the `Targets` objects are supplied to whichever fuel type they are meant to calibrate. Grass, coniferous litter, and deciduous litter can be calibrated individually, separately by specifying them as separate arguments, or together by using the argument *all* or *litter* (they can optionally be separated back into their component fuel types after calibration).

### Calibration

The [`calibrate`](reference.md#duet_tools.calibration.calibrate) function handles all the calculations for scaling DUET values based on the supplied `FuelParameter` objects and their associated `Targets`. The `calibrate` function is meant to be called only once, after all `FuelParameter` objects have been created with all desired `Targets`. Calibration is performed on a `DuetRun` object, and new `DuetRun` object is returned.

Because a `Targets` object can be assigned to any fuel parameter and type, the user may mix and match methods and target values for whichever fuels they wish to alter. There are no requirements to the number of parameters or types that can or cannot be calibrated.

### Using LANDFIRE Targets

In many situations, fuels data for specific burn units may not be available. LANDFIRE is a national dataset that includes 30m-resoltuion Scott and Burgan 40 Fuel Model (SB40) designations. From these data, values for fuel loading, moisture, and depth can be derived. The `landfire` module offers a method of querying SB40 LANDFIRE data for a specific burn unit and assigning calibration targets based on those designations.

To use LANDFIRE data, the user must use the [`landfire`](reference.md#duet_tools.landfire) module to query LANDFIRE data (see below). Once data has been queried, the resulting [`LandfireQuery`](reference.md#duet_tools.landfire.LandfireQuery) object is passed to the [`assign_targets_from_sb40`](reference.md#duet_tools.landfire.assign_targets_from_sb40) function. Because LANDFIRE data most often does not follow a normal distribution, the calibration method defaults to "maxmin", and "meansd" is not recommended. Unlike `assign_targets`, which does not specify fuel parameter or type, both must be specified in `assign_targets_from_sb40`. However, the resulting `Targets` object is treated like any other, and must be assigned to the correct fuel parameter and type(s) in `set_fuel_parameter`.
**NOTE:** the `landfire` module is only available for Python v3.10 and must be installed using `pip install duet-tools[landfire]`


### DuetRun Class

The [`DuetRun`](reference.md#duet_tools.calibration.DuetRun) class is instantiated using the `import_duet` function. It's three attributes correspond to the 3D fuel parameter arrays that may be calibrated: `loading`, `moisture`, and `depth`.

#### Pre-Calibration

Before a DuetRun class is calibrated, it simply holds the surface_*.dat output arrays from a specified directory. A user can conveniently interface with these arrays outside of duet-tools by using the `to_numpy` function (see below).

#### Post-Calibration

The `calibrate` function returns a new instance of the `DuetRun` class containing calibrated arrays. These arrays can then be exported in a few ways using methods from the `DuetRun` class.

The [`to_numpy`](reference.md#duet_tools.calibration.DuetRun.to_numpy)method  returns a numpy array of a specifed fuel parameter and type. Arrays can be returned as an *integrated* 2D array, where the fuel types are combined for a given parameter, or a *separated* 3D array, where the fuel types occupy different z-layers (axis 0).

The [`to_quicfire`](reference.md#duet_tools.calibration.DuetRun.to_quicfire) method writes Fortran files (.dat) to be used in QUIC-Fire. Fuel types are integrated for each fuel parameter, and exported as a 3D array with 1 z-layer. The filenames are set to match QUIC-Fire's expected file naming system.

## Landfire Module

The [`Landfire`](reference.md#duet_tools.landfire) module is an auxiliary module handling the interfacing and processing of LANDFIRE data. Data is queried using the [`query_landfire`](reference.md#duet_tools.landfire.query_landfire) function by providing spatial data and information for the area of interest. An instance of class [`LandfireQuery`](reference.md#duet_tools.landfire.LandfireQuery) is returned, which can be provided to `assign_targets_from_sb40` in the `calibration` module.

Values for fuel loading, fuel moisture, and surface fuel depth are derived from Scott and Burgan 40 Fuel Model designations, using methods developed for FastFuels (citation). These data are available at a 30x30m resolution for the contiguous United States.

When a fuel type is selected (*e.g.* grass or litter), fuel parameter values are derived from only SB40 Fuel Models that are predominantly comprised of that fuel type. Because DUET does not have a designation for shrub fuels, any SB40 Fuel Model with major shrub components are categorized as grass, since their growth patterns will also follow light availability.

If a user attempts to assign targets from a fuel type that is not present in the area of interest, an error will be given. If there is a single parameter value for a given fuel type in the area of interest, targets will automatically be given the "constant" calibration method, and a warning will be issued. Using the "meansd" calibration method is generally discouraged, since values derived from SB40 most often do not follow a normal distribution.
