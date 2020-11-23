# An automated GIS toolbox for watershed delineation with lakes

# An automated GIS toolbox for watershed delineation with lakes -- Overview

Before introducing methods and application procedures of the lake river routing toolbox, an overview about how this toolbox will represent lakes in the routing network is described here. Catchments defined by predefined river system without considering lakes are showed in Figure 1A. Both lake's inlets and outlets are not represented as catchment outlets. With this routing structure (Figure 1A), both lake’s inflow and outflow can’t be explicitly simulated by semi-distributed hydrological models such as SWAT, Raven and HYPE, because the streamflow is only explicitly simulated at each catchment’s outlet in these models. The lake river routing toolbox presented in this paper is developed to solve this problem. 


Lakes within the drainage area (Figure 1A) are divided into two categories by the lake river routing toolbox: 1) lakes that are connected by the predefined river network and 2) lakes that are not connected by the predefined river network (Figure 1A). Connected lake polygons will be used to identify each lake's inlet and outlets and divide catchment polygons in Figure 1A into catchments with connected lakes and catchments without connected lakes (Figure 1B). Each connected lake will be represented by a catchment with connected lake Figure 1B. Both streamflow that routes into the connected lakes via each lake’s inlets and that release from each connected lake’s outlet can be explicitly represented by semi-distributed hydrological models with the routing structure showed in Figure 1B. 


For lakes that are within the drainage area but not connected by the predefined river network, each of them will be represented by a catchment which is generated using each none connected lake’s outlet as a catchment outlet (Figure 1B). The water released from catchment with a none connected lake will directly move into the river channel of downstream catchment if downstream catchment do not have lakes or directly move into the lake of the downstream catchment if downstream catchment has connected lake or none connected lake. For example, water released from Cat_1 will directly move into none connected lake in Cat_2 and water released from Cat_2 will directly move into the river channel in Cat_3 (Figure 1B). The water transfer time from the outlet of the catchment with a none connected lake (Cat_1 or Cat_2) to downstream catchment's lake (Cat_2) or river channel (Cat_3) is neglected.

<figure>
    <p align="center">
    <img src="https://github.com/dustming/RoutingTool/wiki/Figures/Figure1.png" width="100%" height="100%" />
    </p>
    <font size="1">
    <figcaption width="50%"> <b>Figure 1</b>: Lakes in the generated routing network by the lake river routing toolbox. A is the predefined river network and catchment boundary (Catchment boundary (River)) defined by the predefined river network without considering lakes. Figure B, the generated lake river routing structure by this toolbox using predefined river network and lake's polygons.<br>
    </figcaption>
    </font>
</figure>



**Table of Contents**

1. [Install in Windows system](https://github.com/dustming/RoutingTool/wiki/Installation-of-the-toolbox#Install-in-Windows-system)
2. [Install in Ubuntu system](https://github.com/dustming/RoutingTool/wiki/Installation-of-the-toolbox#Install-in-Ubuntu-system) 

## Install in Windows system
### Installation required softwares (QGIS, GRASS and GDAL)
QGIS, GRASS and GDAL needs to be installed and can be dowloaded from [here](https://qgis.org/en/site/forusers/download.html). QGIS 3.x standalone installation is prefered. Normally QGIS 3.x intaller will automatecally intall both Gdal and GRASS 7.x within QGSI folders. Note GRASS and GDAL is not needed for several post processing tools. 

### Install required addtional packages

- Required python packages 

Following python packages are needed: os, sys, numpy, shutli, distutils, tempfile, copy, pandas, scipy, simpledbf, sqlite3 and 'grass_session'. Those packaged need to be installed into the QGIS python. Methods to install python packages in QGIS python can be found [here](https://landscapearchaeology.org/2018/installing-python-packages-in-qgis-3-for-windows/). Option 2 in this link was tested. 

- Required GRASS addons 

Following GRASS addon is needed 'r.accumulate', 'r.clip', 'r.stream.basins', 'r.stream.order', 'r.stream.segment', 'r.stream.snap' and . Procedures to install GRASS addon can be found [here](https://grasswiki.osgeo.org/wiki/AddOns). Note GRASS and GDAL is not needed for several post processing tools.  


### Install Lake river routing toolbox  

- Download the lake river routing toolbox 

Download source code of the lake river routing toolbox (Toolbox_QGIS) from [here](https://github.com/dustming/RoutingTool/tree/master/Toolbox_QGIS). Put the downloaded source code in any prefered folder for example 'C:/Document/Toolbox_QGIS/'. And then add 'C:/Document/Toolbox_QGIS/' into system variables: 'PYTHONPATH' and 'Path'. 

- Modify 'QGIS_Toolbox.bat' in the toolbox source code folder 

An example of the 'QGIS_Toolbox.bat' is shown here, you need modify the path of each variables based on your installation of QGIS. 

```
@echo off

SET OSGEO4W_ROOT=C:\QGIS310
SET GRASS_ROOT=C:\QGIS310\apps\grass\grass78
SET GISBASE=C:\QGIS310\apps\grass\grass78
SET QGISPrefixPath=C:\QGIS310\apps\qgis
SET GRASSBIN=C:\QGIS310\bin\grass78.bat
SET RoutingToolFolder=C:\Users\dustm\Documents\GitHub\RoutingTool\Toolbox_QGIS

@echo off
call %OSGEO4W_ROOT%\bin\o4w_env.bat
call %OSGEO4W_ROOT%\apps\grass\grass78\etc\env.bat
call qt5_env.bat
call py3_env.bat
@echo off

path %PATH%;%OSGEO4W_ROOT%\apps\qgis\bin
path %PATH%;%OSGEO4W_ROOT%\apps\qgis\python\plugins
path %PATH%;%OSGEO4W_ROOT%\apps\Qt5\bin
path %PATH%;%OSGEO4W_ROOT%\apps\Python37\Scripts

path %PATH%;%GRASS_ROOT%\lib
path %PATH%;%GRASS_ROOT%\bin
path %PATH%;%GRASS_ROOT%\script
path %PATH%;%RoutingToolFolder%


set GDAL_FILENAME_IS_UTF8=YES

SET PYTHONHOME=C:\QGIS310\apps\Python37

SET PYTHONPATH=%PYTHONPATH%;%OSGEO4W_ROOT%\apps\qgis\python;%OSGEO4W_ROOT%\apps\qgis\python\plugins;%OSGEO4W_ROOT%\apps\qgis\python\plugins\processing
SET PYTHONPATH=%PYTHONPATH%;%GRASS_ROOT%\etc\python\;%GRASS_ROOT%\etc\python\grass;%GRASS_ROOT%\etc\python\grass\script;%RoutingToolFolder%;%RoutingToolFolder%\Toolbox_Common_Functions

rem Set VSI cache to be used as buffer, see #6448
set VSI_CACHE=TRUE
set VSI_CACHE_SIZE=1000000

set QT_PLUGIN_PATH=%OSGEO4W_ROOT%\apps\qgis\qtplugins;%OSGEO4W_ROOT%\apps\qt5\plugins
set QGIS_PREFIX_PATH=%OSGEO4W_ROOT:\=/%/apps/qgis

cmd.exe    
```

### Validate the installation

First, open a command line window and type

```
QGIS_Toolbox.bat
```
If all system variables are successfully assigned, a new command line window will pop up. Then type

```
python
``` 
And then 

```
from ToolboxClass import LRRT
```
If everything goes on well, 'LRRT' will be successfully loaded.

## Install in Ubuntu system
### Installation required software (QGIS, GRASS, and GDAL)
- Check system default python 3 path
```
which python3 
/user/bin/python3
```
- Install qgis with grass plugin 
Please following the instruction from [here](https://qgis.org/en/site/forusers/alldownloads.html#debian-ubuntu) to install QGIS with GRASS plugin. 
```
sudo apt install gnupg software-properties-common

wget -qO - https://qgis.org/downloads/qgis-2020.gpg.key | sudo gpg --no-default-keyring --keyring gnupg-ring:/etc/apt/trusted.gpg.d/qgis-archive.gpg --import

sudo chmod a+r /etc/apt/trusted.gpg.d/qgis-archive.gpg

sudo add-apt-repository "deb https://qgis.org/debian `lsb_release -c -s` main"

sudo apt update

sudo apt install qgis qgis-plugin-grass

```
- Check QGIS and Grass python environment

In an ubuntu system, the QGIS will be installed with system default python3. we can directly import QGIS using system default python3. 

```
python3 
>>>import qgis
>>>qgis.__file__  ### print qgis module path 
'/usr/lib/python3/dist-packages/qgis/__init__.py'
```
However, GRASS is not installed in the default as python3 site-packages. Need to setup GRASS python environment by following steps.  

a) install GRASS GUI and GRASS development package 

```
sudo apt install grass-gui ### install grass GUI 
sudo apt install grass-dev ### install grass development package  
```

b) load GRASS GUI and find path grass python modules

```
grass ## to load the grass 
### open the interactive python shell within grass GUI.
>>>  import os 
>>>  os.environ['PYTHONPATH']
/usr/lib/grass78/etc/python
>>>  os.environ['GISBASE']
'/usr/lib/grass78'
```

We need to add printed PYTHONPATH and GISBASE into the system environment variable every time we use the toolbox. 

```

export GISBASE='/usr/lib/grass78'
export PYTHONPATH=$PYTHONPATH:'/usr/lib/grass78/etc/python'
python3
>>> import grass.script as grass 
```
Now we know the key variables needed to set up the QGIS and GRASS python environment.

## Install the needed library
This part is the same as the procedure in the windows system. 

## Install toolbox 
Download the toolbox, and define the following system variables before using toolbox

```
export RoutingToolFolder='/home/dustming/Documents/RoutingTool-master/Toolbox_QGIS'
export GISBASE='/usr/lib/grass78'
export QGISPrefixPath='/usr'

export PYTHONPATH=$PYTHONPATH:'/usr/lib/grass78/etc/python'### folder has a grass folder
export PYTHONPATH=$PYTHONPATH:$RoutingToolFolder
export PYTHONPATH=$PYTHONPATH:'/usr/share/qgis/python/plugins' ## folder has db_manager and processing
export PYTHONPATH=$PYTHONPATH:'/usr/share/qgis/python' ## folder has plugin and console 
export PYTHONPATH=$PYTHONPATH:$RoutingToolFolder/Toolbox_Common_Functions
export PATH=$PATH:$RoutingToolFolder
export PATH=$PATH:$RoutingToolFolder/Toolbox_Common_Functions
```
### test run 
Go to tests folder and then 
```
python3 -m pytest test_RoutingNetworkTopologyUpdateToolset_riv
```
