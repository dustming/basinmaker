=====================
Installation overview
=====================

Overview
========

BasinMaker is a python package depends on several existing GIS platforms. So, the installation of BasinMaker includes two steps: 1) setup the python environment for the dependent GIS platforms; and 2) install BasinMaker itself.

Two installation modes (light installation and full installation) are available. The light installation will allow user to use BasinMaker to post process an existing routing product, such as the Basinmaker derived North American lake-river routing product. But it cannot be used to delineate a lake-river routing structure from DEM. The combination of BasinMaker light installation and the North American lake-river routing product could generate lake-river routing structures satisfying many user modeling demands. While the full installation of BasinMaker enables users to delineate a new lake-river routing structure from a user specified DEM.

For light installation (recommended unless users know for sure they need to delineate their watershed from scratch from a DEM), only QGIS or ArcGIS pro is needed. The python environment for both QGIS and ArcGIS pro can be easily compiled within the anaconda environment under different OS systems. The instruction about light installation procedure can be found in :ref:`Light installation`. Note that users can install both the light version and the full version on the same OS system.   

For full installation, both GRASS GIS and QGIS are needed. It is quite a challenge to setup a python environment for QGIS and GRASS together. Here, two procedures are provided for Windows and Ubuntu OS systems, respectively. Two procedures have been tested on several machines. But it is possible that provided procedures are not working on your machine. Please feel free to create an issue on the GitHub or email m43han@uwaterloo.ca, we are happy to make it work on your machine. The instruction for Windows and Ubuntu system can be found in :ref:`Full installation`.
    

Light installation
==================

QGIS with anaconda
------------------

#. Install anaconda

    The installer of anaconda can be installed from `here <https://www.anaconda.com/>`_


#. Create an empty python environment and tjen active it  

    .. code-block::
      
      conda create --name <any_name_for_env>
      conda activate <any_name_for_env>
   
   
#. Install QGIS

    .. code-block:: 

      conda install -c conda-forge qgis
   
   
#. Install BasinMaker

    .. code-block::
    
      git clone https://github.com/dustming/basinmaker.git basinmaker
   
      cd ./basinmaker
   
      python setup.py develop
   
   
#. Install dependent packages 

    .. code-block::
  
      pip install pandas pytest scipy simpledbf netCDF4


ArcGIS pro with anaconda (Windows only)
---------------------------------------


#. Install anaconda

    The installer of anaconda can be installed from `here <https://www.anaconda.com/>`_


#. Create an empty python environment and then active it 

    .. code-block::
    
      conda create --name <any_name_for_env>
      conda activate <any_name_for_env>
   
   
#. Install arcpy and arcgis 

    .. code-block::
    
      conda install -c esri arcpy arcgis
   
   
#. Install BasinMaker 

    .. code-block::
    
      git clone https://github.com/dustming/basinmaker.git basinmaker
   
      cd ./basinmaker
   
      python setup.py develop
   
   
#. Install dependent packages 

    .. code-block::
    
      pip install pandas pytest scipy simpledbf netCDF4


Full installation
==================

QGIS and GRASS in Windows
-------------------------

#. Installation of QGIS and GRASS using OSGEO4W: 
    
    For the Windows system, we can install both GRASS and QGIS within OSGEO4W environment.  
    
    The QGIS installer can be downloaded from `here <https://qgis.org/en/site/forusers/download.html>`_. Please using **OSGeo4W installer**.
    
    The GRASS installer can be found in `here <https://grass.osgeo.org/download/windows/>`_.  Please using **OSGeo4W installer**.
    
    We would suggest to install QGIS and GRASS outside the **C/:Program Files**. Better to install them into a folder path without space in the folder name.

#. Install BasinMaker 

    .. code-block::
    
      git clone https://github.com/dustming/basinmaker.git basinmaker
            
#. Setup GRASS and QGIS python environment

    The python environment for QGIS and GRASS GIS in Windows can be set up by modifying the following basinmaker.bat file. 
    
    * Please change OSGEO4W_ROOT to your OSGEO4W installation folder at line 3.
    * Please change the grass78.* in line 11 and 13 to your GRASS GIS version number.
    * Please double check the paths defined in the basinamker.bat file exists in your machine
    * Copy the basinmaker.bat file into path_to_basinmaker_folder/basinmaker/basinmaker.bat

    .. code-block::
      :linenos:
      
      @echo off
      rem define OSGEO4W_ROOT, change it to your OSGEO4W installation folder
      set OSGEO4W_ROOT=C:\OSGeo4W64
      
      rem setup OSGEO4W environment 
      call "%OSGEO4W_ROOT%\bin\o4w_env.bat"
      call qt5_env.bat
      call py3_env.bat
      
      rem  setup environment variables for GRASS GIS
      set GRASS_ROOT=%OSGEO4W_ROOT%\apps\grass\grass78
      set GISBASE=%GRASS_ROOT%
      set GRASSBIN=%OSGEO4W_ROOT%\bin\grass78.bat
      call "%GRASS_ROOT%\etc\env.bat"
      path %PATH%;%GRASS_ROOT%\lib
      path %PATH%;%GRASS_ROOT%\bin
      path %PATH%;%GRASS_ROOT%\script
      set PYTHONPATH=%GRASS_ROOT%\etc\python;%GRASS_ROOT%\etc\python\grass;%GRASS_ROOT%\etc\python\grass\script;%PYTHONPATH%
      
      rem for qgis 
      path %OSGEO4W_ROOT%\apps\qgis\bin;%PATH%
      set QGIS_PREFIX_PATH=%OSGEO4W_ROOT:\=/%/apps/qgis
      set GDAL_FILENAME_IS_UTF8=YES
      rem Set VSI cache to be used as buffer, see #6448
      set VSI_CACHE=TRUE
      set VSI_CACHE_SIZE=1000000
      set QT_PLUGIN_PATH=%OSGEO4W_ROOT%\apps\qgis\qtplugins;%OSGEO4W_ROOT%\apps\qt5\plugins
      set PYTHONPATH=%OSGEO4W_ROOT%\apps\qgis\python;%OSGEO4W_ROOT%\apps\qgis\python\plugins;%PYTHONPATH%
      
      cd ..
      python setup.py develop 
      
      cmd.exe
    
#. Validate the GRASS and QGIS python environment
     
    * Run the saved basinmaker.bat file in step 3.
    * Try to load following packages

    .. code-block::
       
      >where python    
      >C:\OSGeo4W64\apps\Python37\python.exe
  
      >python
      >>>from qgis.core import *
      >>>import qgis
      >>>from qgis.analysis import QgsNativeAlgorithms
      >>>from qgis.PyQt.QtCore import *
      >>>from qgis import processing
      Application path not initialized
      >>>from processing.core.Processing import Processing
      >>>from processing.tools import dataobjects
      >>>import grass.script as grass
      >>>from grass.script import array as garray
      >>>from grass.script import core as gcore
      >>>import grass.script.setup as gsetup
      >>>from grass.pygrass.modules.shortcuts import general as g
      >>>from grass.pygrass.modules.shortcuts import raster as r
      >>>from grass.pygrass.modules import Module
      
      
#. Install dependent packages

    .. code-block::

      pip install simpledbf grass_session sqlite3 pandas distutils


#. Install GRASS GIS addons

    Following GRASS GIS addons(r.accumulate,r.clip,r.stream.basins and r.stream.snap) needs to be installed. How to install GRASS GIS addon 
    can be found in `here <https://grass.osgeo.org/download/addons/>`_. 
  
#. Run basinmaker.bat everytime you want to use basinmaker python packages


QGIS and GRASS in Ubuntu
------------------------
    
#. Installation of QGIS and GRASS 
    
    For ubuntu system, both QGIS and GRASS GIS can be installed at the same time by installing the QGIS with GRASS addon. 
    The installation procedure is the following comes from `here <https://qgis.org/en/site/forusers/alldownloads.html#debian-ubuntu>`_. 
    
    .. code-block::
    
      $sudo apt install gnupg software-properties-common
      $wget -qO - https://qgis.org/downloads/qgis-2020.gpg.key | sudo gpg --no-default-keyring --keyring gnupg-ring:/etc/apt/trusted.gpg.d/qgis-archive.gpg --import
      $sudo chmod a+r /etc/apt/trusted.gpg.d/qgis-archive.gpg
      $sudo add-apt-repository "deb https://qgis.org/debian `lsb_release -c -s` main"
      $sudo apt update
      $sudo apt install qgis qgis-plugin-grass
      
    * Install GRASS GIS GUI and development packages 
    
    .. code-block::
      
      $sudo apt install grass-gui 
      $sudo apt install grass-dev        

#. Install BasinMaker::

    $git clone https://github.com/dustming/basinmaker.git basinmaker
    $cd basinmaker
    $sudo python3 ./basinmaker/setup.py develop 
    
#. Setup GRASS and QGIS python environment

    The python environment for QGIS and GRASS GIS in Ubuntu can be set up by modifying the basinmaker.sh file 
    
    * Please change the grass78.* in line 1 and 4 to your GRASS GIS version number.
    * Please double check the paths defined in the basinamker.sh file exists in your machine
    * Copy the basinmaker.sh file into path_to_basinmaker_folder/basinmaker/basinmaker.sh

    .. code-block::
      :linenos:
      
      export GISBASE='/usr/lib/grass78'
      export QGISPrefixPath='/usr'
      
      export PYTHONPATH=$PYTHONPATH:'/usr/lib/grass78/etc/python'  ### folder has a grass folder
      export PYTHONPATH=$PYTHONPATH:'/usr/share/qgis/python/plugins' ## folder has db_manager and processing
      export PYTHONPATH=$PYTHONPATH:'/usr/share/qgis/python' ## folder has plugin and console 
      
#. Validate the GRASS and QGIS python environment
     
    * Run the saved basinmaker.sh file in step 3.
    
    .. code-block::

      $source ./basinmaker.sh
    
    * Try to load following packages

    .. code-block::
         
      $python3
      >>>from qgis.core import *
      >>>import qgis
      >>>from qgis.analysis import QgsNativeAlgorithms
      >>>from qgis.PyQt.QtCore import *
      >>>from qgis import processing
      Application path not initialized
      >>>from processing.core.Processing import Processing
      >>>from processing.tools import dataobjects
      >>>import grass.script as grass
      >>>from grass.script import array as garray
      >>>from grass.script import core as gcore
      >>>import grass.script.setup as gsetup
      >>>from grass.pygrass.modules.shortcuts import general as g
      >>>from grass.pygrass.modules.shortcuts import raster as r
      >>>from grass.pygrass.modules import Module

#. Install dependent packages

    .. code-block::

      pip install simpledbf grass_session sqlite3 pandas distutils


#. Install GRASS GIS addons

    Following GRASS GIS addons(r.accumulate,r.clip,r.stream.basins and r.stream.snap) needs to be installed. How to install GRASS GIS addon 
    can be found in `here <https://grass.osgeo.org/download/addons/>`_.     

#. Run basinmaker.sh everytime you want to use basinmaker python packages

