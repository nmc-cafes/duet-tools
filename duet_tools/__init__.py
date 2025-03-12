from duet_tools.calibration import (
    DuetRun,
    Targets,
    FuelParameter,
    import_duet,
    assign_targets,
    set_density,
    set_moisture,
    set_height,
    set_fuel_parameter,
    calibrate,
)

from duet_tools.utils import (
    write_array_to_dat,
    read_dat_to_array,
)

from duet_tools.inputs import InputFile


__all__ = [
    "DuetRun",
    "Targets",
    "FuelParameter",
    "import_duet",
    "assign_targets",
    "set_density",
    "set_moisture",
    "set_height",
    "set_fuel_parameter",
    "calibrate",
    "assign_targets_from_sb40",
    "LandfireQuery",
    "query_landfire",
    "write_array_to_dat",
    "read_dat_to_array",
    "InputFile",
]


def _missing_landfire(*args, **kwargs):
    raise ImportError(
        "The 'landfire' module requires additional dependencies. "
        "Please reinstall with: pip install package-name[landfire]"
    )


try:
    from duet_tools.landfire import (
        LandfireQuery,
        query_landfire,
        assign_targets_from_sb40,
    )

    __all__.extend(["LandfireQuery", "query_landfire", "assign_targets_from_sb40"])
except ImportError:
    LandfireQuery = _missing_landfire
    query_landfire = _missing_landfire
    assign_targets_from_sb40 = _missing_landfire
