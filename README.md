# Command-line Zonal Stats and Volume Extraction Tool
Computes RGB and Multispectral vegetation indices, cut/fill volumes using digital surface models, and allows extraction of reduced (zonal statistics) values of pixels contained within the spatial domain of polygon geometries.

[Back to home](https://github.com/alexanderhernandez-USDA/Zonal-Stats-QGIS/blob/main/homedir.md)

# Contents
[Workflow](#Workflow)

[Packages and Environment](#Packages/Environment)

[Installation and Setup](#Installation/Setup)

[Usage - Running Zonal Stats](#Running-Zonal-Stats)

[Adding New Indices](#Customizing-or-adding-your-own-indices)

### Workflow
The Zonal Stats tool has the following basic operational workflow:
1. Read input geopackage with geopandas
2. Read input rasters and perform desired calculations (vegetation indices, volume, etc.) using rasterio and scikit-image 
3. Write calculations to temporary rasters using scikit-image, preserve metadata with exiftool
4. Using exactextract python package, get median values (or sum in the case of volume calculation) of calculations within user created polygons or within polygons generated from point geometries 
5. Append exactextract results as new columns to inputted geopackage
6. Write updated geopackage to file using geopandas

### Packages/Environment
Currently, the most recent version of the environment for Zonal Stats is using Python 3.9.  
Zonal Stats uses the following major python packages (not including those part of the default Python installation)
- Geopandas (and therefore Pandas): https://geopandas.org/
- Rasterio: https://rasterio.readthedocs.io/
- Exactextract: https://github.com/isciences/exactextract
- Shapely: https://shapely.readthedocs.io/en/stable/manual.html
- Scikit-image (skimage): https://scikit-image.org/
- Numpy: https://numpy.org/

## Installation/Setup

The first thing you need is to have miniconda already installed on your computer (https://docs.anaconda.com/miniconda/install/).

You will need to create a conda environment to run zonal stats 3. This can be done by using the geospatial3.9_env.yml file in this repository.
The YML file can be accessed directly from:

https://github.com/alexanderhernandez-USDA/Zonal_Stats_3/blob/main/geospatial3.9_env.yml

You can then create the miniconda environment using the following sintax:

```
conda create -f geospatial3.9_env.yml       # This will create an environement called 'geospatial'
```


## Running Zonal Stats
To run zonal_stats, you will need to activate the conda environment created from the geospatial3.9_env.yml.

```
conda activate geospatial       # This will create an environement called 'geospatial'
```

Now that you have created and activated the miniconda environment - you can download the following python script (*.py) and configuration file (*.conf) that contains all of the available vegetation indices to your disk:

https://github.com/alexanderhernandez-USDA/Zonal_Stats_3/blob/main/zonal_stats_3.py
https://github.com/alexanderhernandez-USDA/Zonal_Stats_3/blob/main/indices.conf

Next `cd` or navigate to the `zonal_stats` folder where you saved the Python script. In here, you will find `zonal_stats_3.py`. This program can run any number of index/data extractions as specified by flags, and has two required inputs, an input geopackage and an output geopackage. The command is run like so:
```
python3 zonal_stats_3.py [OPTIONS] <index_flags> <index_requirements> /path/to/input.gpkg /path/to/output.gpkg
```

Please note that at least one index option (-a, -i, -n, -v or -V) and its corresponding input directory must be used in order to run

All non-volumetric calculations and extractions require band order to be specified as well. The default band order is [red,green,blue,rededge,nir], and you can
use the default by using empty square brackets, like [], when specifying band order.

Additionally, images in your input folder need to have a date in their file, delimited by ., _, or -. For example, three vaild date formats would be:

YYYY_MM_DD  
DD.MM.YYYY  
MM-DD-YYYY  

Zonal stats has the following global options:
```
Global Options
    -h      Show this help
    -q      Suppress non-error messages
    -l      Show available vegetation indices
    -u <uid>        Specify a unique ID column name from geopackage/shapefile (verify that it is actually unique, there cannot be any repeats)
    -t <threads>    Specify number of threads for multithreading, by default one thread is used, you can specify up to total thread count - 2 (recommended amount is number of images to be processed)
    -p              Use a geopackage that has points instead of polygons, returns values at those points
    -S <distance>   Use a geopackage that has points instead of polygons, and specify a distance in meters for a square buffer around the point
    -C <distance>   Use a geopackage that has points instead of polygons, and specify a distance in meters for a circular buffer around the point
    -o <directory>  Generate multi-layer tifs containing intermediate calculation/extraction data by date in the specified output directory
```

The options should be specified in the `[OPTIONS]` part of the command. After specifying global options, you have to specify _at least one_ index or extraction to run. The options for index/extractions are as follows:
```
Index / Calculation options : Specify a directory (and band order if not a volume calc) per flag. Multiple index/calculation opetions may be used in a single command
    -i [<indices>] <directory> [<bands>]     Pass a list of vegetation indices to be run, in format [index1,index2,...] (try -l option to see available indices)
    -a <directory> [<bands>]          Run all vegetation indices (except DGCI, VOLUME, and NONE)
    -n <directory> [<bands>]          Run zonal_stats without any indices (i.e. on a thermal image to get median temperatures)
    -v <directory>          Specify a path to a folder to run volume calculation on using a plane average calculation.
    -V <directory> <second_dsm.tif>      Specify a path to a folder to run volume calculation on using a secondary DSM to calculate heights, the second DSM should be specified after the folder name.
```

For the -i option, a complete list of vegetation indices can be found by running `python3 zonal_stats_3.py -l`

Below are some examples for running Zonal stats:
```
python3 zonal_stats_3.py -i [BI,SCI,GLI] flight/rasters/ [red,green,blue,rededge,nir] flight/package.gpkg zonal_stats.gpkg
#Runs with indices BI, SCI, and GLI

python3 zonal_stats_3.py -u pid -a flight/rasters/ [red,green,blue,rededge,nir] flight/package.gpkg zonal_stats.gpkg
#Runs with all indices, unique ID column set to 'pid'

python3 zonal_stats_3.py -a flight/rasters/ [red,green,blue,rededge,nir] flight/package.gpkg zonal_stats.gpkg
#Runs all indices with band order red, green, blue, NIR, RedEdge

python3 zonal_stats_3.py -n flight/thermals/ [swir] flight/rasters/ flight/package.gpkg zonal_stats.gpkg
#Runs zonal stats on the thermals directory getting raw values

python3 zonal_stats_3.py -v flight/dsms/ flight/package.gpkg zonal_stats.gpkg
# Performs volume calculation using a plane average

python3 zonal_stats_3.py -V flight/dsms/ flight/ref_dsm.tif flight/package.gpkg zonal_stats.gpkg
# Performs volume calculation using a reference DSM raster

python3 zonal_stats_3.py -t 12 -v flight/dsms/ flight/package.gpkg zonal_stats.gpkg
# Uses 12 threads to perform volume calculation

python3 zonal_stats_3.py -p -i [BI] flight/rasters/ [red,green,blue,rededge,nir] flight/point_package.gpkg point_stats.gpkg
# Gets values at each point of a point-based geopackage after BI calculation

python3 zonal_stats_3.py -C 1 -i [BI] flight/rasters/ [red,green,blue,rededge,nir] flight/point_package.gpkg circle_stats.gpkg
# Makes circular buffers with a 1 meter radius from point-based geopackage and performs BI calculation and then gets stats for each buffer

python3 zonal_stats_3.py -S 1 -i [BI] flight/rasters/ [red,green,blue,rededge,nir] flight/point_package.gpkg square_stats.gpkg
# Makes square buffers with a 1 meter radius from point-based geopackage and performs BI calculation and then gets stats for each buffer

python3 zonal_stats_3.py -n flight/thermals/ [swir] -i [BI] flight/rasters/ [] flight/package.gpkg output.gpkg
# Get raw values from SWIR images and gets the BI index from images in rasters with default band order

python3 zonal_stats_3.py -o flight/outputs/ -i [BI,SCI,GLI] flight/raster/ [] flight/package.gpkg output.gpkg
# Get BI, SCI and GLI indices and output rasters with calculation data in flight/outputs/
```

# Customizing or adding your own indices
Aside from the built in vegetation indices and calculations, custom indices/calculation can be added via the indices.conf file within the same folder as zonal_stats_3.py. 

These function can be defined in indices.conf like so:
```
[NAME]
    desc: Description
    calc: Calculation
```
Where NAME is the index name (in all caps) and cannot be an existing index name (like SCI, BI, etc). Description is a brief description of the index, and calculation is a single line of valid python operations (+,-,/,\*,%,(),etc.) Note that functions from Python's Math library and the Numpy library are available and can be used here with `math.<function>` for math functions and `np.<numpy_function>` for Numpy functions. Raster bands are defined as variables in these equations like so:

r = red  
b = blue  
g = green  
n = nir  
re = rededge  

For example:
```
[EX]
    desc: Example index that adds red, green and blue bands and divides by 10
    calc: (r+g+b)/10
```
Another example with Numpy functions:
```
[EXNP]
    desc: Example index using a Numpy function
    calc: np.arctan(r*b)
```

Indices defined in indices.conf can be used in zonal stats by their name via the -i flag, like so:
```
python3 zonal_stats_3.py -i [EX] /path/to/images /path/to/input.gpkg /path/to/output.gpkg
```

For more information, run `python3 zonal_stats_3.py -h` and for all vegetation indices run `python3 zonal_stats_3.py -l`
