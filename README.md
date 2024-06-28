# duet-tools

## Quick-Links

[Documentation](https://nmc-cafes.github.io/duet-tools/) - [PyPi Package](https://pypi.org/project/duet-tools/)

## What is duet-tools?

duet-tools is a Python package that provides a convenient interface for programmatically working with the outputs of the DUET program developed by [McDanold et al. (2023)](https://doi.org/10.1016/j.ecolmodel.2023.110425) at Los Alamos National Lab. Central to the package is the ability to calibrate the values in DUET outputs to match targets supplied by the user or national datasets.

The goals of duet-tools are to:

1. Read in DUET outputs for easy manipulation.
2. Calibrate DUET outputs by scaling the magnitude of the parameter values while retaining the spatial distributions from DUET.
3. Provide a platform for the future development of additional tools centered around DUET

## Installation

quicfire-tools can be installed using `pip`.

### pip

```bash
pip install duet-tools
```


## Issues

If you encounter any issues with the quicfire-tools package, please submit an issue on the duet-tools GitHub
repository [issues page](https://github.com/nmc-cafes/duet-tools/issues).
