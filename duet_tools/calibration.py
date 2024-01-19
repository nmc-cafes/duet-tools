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
