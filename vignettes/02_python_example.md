# Steps to Run Gridded Weather Generator from Phython Environment

**Warning: This is still in progress and may not work**

Prerequisites: 
1. Install latest version of R via https://mirror.lyrahosting.com/CRAN/
2. Create a new python environment with packages rpy2, r-base, r-essentials, pandas, numpy2



```python
# Step 1:
#  set all PATH variables

import os
os.environ['PATH'] = 'C:/Program Files/R/R-4.1.2/bin/x64' + os.pathsep + os.environ.get('PATH', '')
os.environ['PYTHONHOME'] = 'C:/Users/taner/Anaconda3/envs/wegentest'
os.environ['PYTHONPATH'] = 'C:/Users/taner/Anaconda3/envs/wegentest/Lib/site-packages'
os.environ['R_HOME'] = 'C:/Program Files/R/R-4.1.2'
os.environ['R_USER'] = 'C:/Users/taner/Anaconda3/envs/wegentest/Lib/site-packages/rpy2'

# Check if variables are correctly defined for rpy2 
import rpy2.situation
for row in rpy2.situation.iter_info():
    print(row)
```


```python
#Import necessary packages
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

# This is needed for conversion between R and Python syntax
d = {'package.dependencies': 'package_dot_dependencies',
     'package_dependencies': 'package_uscore_dependencies'}

# Load core packages
base = importr('base')
utils = importr('utils')
utils.chooseCRANmirror(ind=1) # select the first mirror in the list
devtools = utils.install_packages('devtools')
devtools = importr('devtools', robject_translations = d)

# Install gridwegen from Github master branch
gridwegen = devtools.install_github("tanerumit/gridwegen")
gridwegen = importr('gridwegen', robject_translations = d)
```


```python
# Path to output files
out_path = "C:/wegentest/ntoum/"
nc_file =  "ntoum.nc"
nc_dimnames = base.list(x = "lon", y = "lat", time = "time")
variables = base.c("precip", "temp", "temp_min", "temp_max")
variable_labels = base.c("Precipitation", "Mean Temperature", "Minimum Temperature", "Maximum Temperature")
variable_units = base.c("mm/day", "DegC", "DegC", "DegC")

# Read-in gridded weather data from netcdf
nc_data_ini = gridwegen.readNetcdf(
    nc_path = base.paste0(out_path, "data/"),
    nc_file = nc_file,
    nc_dimnames = nc_dimnames,
    nc_variables = variables,
    origin_date = base.as_Date("1981-01-01"),
    has_leap_days = True)

nc_data_ini
```


```python
#### NOTES ####

# R to Python syntax changes
# replace .(dots) in the parameters with underscores (e.g., nc.path becomes nc_path)
# replace TRUE with True

```
