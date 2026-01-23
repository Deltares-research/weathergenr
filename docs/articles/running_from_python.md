# Running from Python

## Overview

This vignette demonstrates how to run the **weathergenr** R package from
a **Python** environment using **rpy2**. The workflow is designed to be
portable across machines and operating systems and avoids hard-coded,
user-specific paths.

The recommended sequence is:

1.  Install R and make it discoverable.
2.  Create a Python environment with `rpy2`.
3.  Configure minimal environment variables (`R_HOME`, `PATH`).
4.  Install and call **weathergenr** from Python.

## Prerequisites

### Install R

Install the latest version of R from CRAN and verify that it is
available from the command line:

``` bash
R --version
```

### Create a Python environment

``` bash
conda create -n weathergenr-py python=3.11 -y
conda activate weathergenr-py
conda install -c conda-forge rpy2 -y
pip install numpy pandas
```

Avoid setting `PYTHONHOME` and `PYTHONPATH` unless debugging a
misconfigured environment.

## Configure environment variables for rpy2

``` python
import os
import shutil
from pathlib import Path

def _windows_guess_r_home():
    candidates = [
        Path("C:/Program Files/R"),
        Path("C:/Program Files (x86)/R"),
    ]
    for base in candidates:
        if base.exists():
            versions = sorted(base.glob("R-*"), reverse=True)
            if versions:
                return str(versions[0])
    return None

r_home = os.environ.get("R_HOME")

if not r_home and os.name == "nt":
    r_home = _windows_guess_r_home()

if r_home:
    os.environ["R_HOME"] = r_home
    r_bin = Path(r_home) / "bin"
    if os.name == "nt":
        r_bin_x64 = r_bin / "x64"
        if r_bin_x64.exists():
            r_bin = r_bin_x64
    os.environ["PATH"] = str(r_bin) + os.pathsep + os.environ.get("PATH", "")
else:
    if shutil.which("R") is None:
        raise RuntimeError("R not found. Install R or set R_HOME.")

print("R_HOME:", os.environ.get("R_HOME", "<not set>"))
print("R on PATH:", shutil.which("R") or "<not found>")
```

## Install and load weathergenr from GitHub

``` python
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

translations = {
    "package.dependencies": "package_dot_dependencies",
    "package_dependencies": "package_uscore_dependencies",
}

base = importr("base")
utils = importr("utils")
utils.chooseCRANmirror(ind = 1)

if not ro.r('requireNamespace("devtools", quietly = TRUE)')[0]:
    utils.install_packages("devtools")

devtools = importr("devtools", robject_translations = translations)
devtools.install_github("tanerumit/weathergenr")

weathergenr = importr("weathergenr", robject_translations = translations)
```

## Run a minimal example

``` python
import os
from pathlib import Path
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

base = importr("base")
weathergenr = importr("weathergenr")

ncfile = base.system_file("extdata", "ntoum_era5_data.nc", package="weathergenr")
ncdata = weathergenr.read_netcdf(ncfile)

outdir = Path(os.environ.get("WEATHERGENR_OUTDIR", "weathergenr-output")).resolve()
outdir.mkdir(parents=True, exist_ok=True)

variables = ro.StrVector(["precip", "temp", "temp_min", "temp_max"])

seed = 123

stochastic_weather = weathergenr.generate_weather(
  obs_data = ncdata[0],
  obs_grid = ncdata[1],
  obs_dates = ncdata[2],
  vars = variables,
  n_years = 20,
  start_year = 2020,
  year_start_month = 1,
  n_realizations = 3,
  warm_var = "precip",
  warm_signif = 0.90,
  warm_pool_size = 10000,
  annual_knn_n = 100,
  wet_q = 0.2,
  extreme_q = 0.8,
  out_dir = str(outdir),
  seed = seed,
  parallel = False,
  verbose = False
)

print("Completed. Output written to:", outdir)
```
