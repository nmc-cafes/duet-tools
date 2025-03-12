import requests
from setuptools import find_packages, setup
import sys

installing_landfire = any(arg.startswith("duet-tools[landfire]") for arg in sys.argv)
using_python_312 = sys.version_info >= (3, 12)
if using_python_312 and installing_landfire:
    sys.exit(
        "Error: 'duet-tools[landfire]' is only available for Python <3.12.\n"
        "Use Python 3.10 or 3.11 to install with the 'landfire' extra."
    )


def read_file(fname):
    with open(fname, encoding="utf-8") as fd:
        return fd.read()


def get_version():
    """Get the version number."""
    url = "https://api.github.com/repos/nmc-cafes/duet-tools/releases/latest"
    response = requests.get(url)
    response.raise_for_status()
    version = response.json()["tag_name"]
    return version[1:]  # Remove the leading "v" from the version number


def get_requirements(file_name: str):
    """Get the requirements from the requirements.txt file."""
    requirements = []
    with open(f"requirements/{file_name}.txt", encoding="utf-8") as fd:
        for line in fd:
            line = line.strip()
            if line and not line.startswith("#"):
                requirements.append(line)
    return requirements


NAME = "duet-tools"
DESCRIPTION = "Output management tools for LANL DUET program"
LONG_DESCRIPTION = read_file("README.md")
VERSION = get_version()
LICENSE = "MIT"
URL = "https://github.com/nmc-cafes/duet-tools"
PROJECT_URLS = {"Bug Tracker": f"{URL}/issues"}
INSTALL_REQUIRES = get_requirements("requirements")
EXTRAS_REQUIRE = {"landfire": get_requirements("requirements_landfire")}

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    license=LICENSE,
    url=URL,
    project_urls=PROJECT_URLS,
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering",
    ],
    package_dir={"": "."},
    packages=find_packages(exclude=["docs", "tests"]),
    package_data={"duet_tools": ["data/*"]},
    include_package_data=True,
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    python_requires=">=3.10",
)
