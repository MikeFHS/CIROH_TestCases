# CIROH_TestCases
AutoRoute  / FloodSpreader Test Cases to be used in the CIROH Workshop.

AutoRoute samples cross-sections from the DEM data, and then calculates the depth, top-width, and velocity for a given flow at each 
    cross-section using a volume-fill approach.
    
FloodSpreader uses the information from AutoRoute to create a flood inundation map.

# Conda Environment:
    Typically AutoRoute and FloodSpreader are run using Anaconda Prompts (https://www.anaconda.com/download)
    To activate the "autoroute" environment open an Anaconda Command Prompt and type "conda env create -f environment.yml"

# Model and Data Sources:
    autoroute.exe and floodspreader.exe are both from ERDC: https://erdc-library.erdc.dren.mil/jspui/handle/11681/38783
    
    DEM data is from the National Elevation Dataset (1/3 arc second) - https://apps.nationalmap.gov/downloader/
    Land Cover data is from the National Land Cover Database 2011 - https://www.mrlc.gov/data/nlcd-2011-land-cover-conus
    Streamlines are from GeoGLoWS - http://geoglows-v2.s3-website-us-west-2.amazonaws.com/#streams/
    Return period flow data is from GeoGLoWS - http://geoglows-v2-retrospective.s3-website-us-west-2.amazonaws.com/#return-periods/
    FEMA S_FLD_HAZ_AR shapefile is only provided for context.
        The file (from 30067C_20221116) was obtained from https://hazards-fema.maps.arcgis.com/apps/webappviewer/index.html

# Model/Data/Files are AS-IS with No Warranty:
Any scripts, files, etc. provided in this document are for testing purposes only and are provided "as is" and without any warranty of any kind.
