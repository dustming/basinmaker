============
Installation
============



Overview
========

BasinMaker is a python package work on several existing GIS platforms. So, the installation of BasinMaker includes two steps: 1) setup the python environment for the dependent GIS platforms; and 2) install BasinMaker itself.

Two installation modes (light installation and full installation) are available. The light installation will allow user to use BasinMaker to post process an existing routing product, such as the Basinmaker derived North American lake-river routing product. But it cannot be used to delineate a lake-river routing structure from DEM. The combination of BasinMaker light installation and the North American lake-river routing product could generate lake-river routing structures satisfying many user modeling demands. While the full installation of BasinMaker enables users to delineate a new lake-river routing structure from a user specified DEM.

For light installation (recommended unless users know for sure they need to delineate their watershed from scratch from a DEM), only QGIS or ArcGIS pro is needed. The python environment for both QGIS and ArcGIS pro can be easily compiled within the anaconda environment under different OS systems. The instruction about light installation procedure can be found in :ref:`Light installation`. Note that users can install both the light version and the full version on the same OS system.   

For full installation, both GRASS GIS and QGIS are needed. It is quite a challenge to setup a python environment for QGIS and GRASS together. Here, two procedures are provided for Windows and Ubuntu OS systems, respectively. Two procedures have been tested on several machines. But We can't guarantee install procedures work on every machine, even Windows or Ubuntu machines, but if you run into a problem create an issue on the GitHub and time permitting, we will try to help. The instruction for Windows and Ubuntu system can be found in :ref:`Full installation`. If you managed to do the full installation on a different operating system, we would be grateful if you could document and share the detailed installation procedure that was successful (email: m43han@uwaterloo.ca).
    

Light installation
==================

QGIS with anaconda
------------------

#. Install anaconda

    The installer of anaconda can be installed from `here <https://www.anaconda.com/>`_. Note for windows system, please activate the 'Register Anaconda3 as my default python 3.8' 


#. Create an empty python environment and then active it  
    
    For windows system, search and open "Anacoda Prompt" (Windows) to active a conda command line. Then
   
    .. code-block::
      
      conda create --name <any_name_for_env>
      conda activate <any_name_for_env>
   
   
#. Install QGIS

    .. code-block:: 

      conda install -c conda-forge qgis
   
   
#. Install BasinMaker 

    .. code-block::
      
      python -m pip install basinmaker   
   
#. Install dependent packages (This may take a few minutes) 

    .. code-block::
  
      python -m pip install pandas pytest scipy simpledbf netCDF4 jupyter

#. Test validation 
     
    Please download the test data and scripts from `here <https://github.com/dustming/RoutingTool/wiki/Files/test.zip>`_. and unzip it to a folder, the path of this folder will refer as path_test_data in following section. Then

    * PyTables is not installed. No support for HDF output.
    * SQLalchemy is not installed. No support for SQL output.    
    * Warnings
    
    .. code-block::
     
       
      cd path_test_data/test
      python test_light_installation.py
      (... some messages)
      ####################################
      BasinMaker is successfully installed
      ####################################

#. Users must active this conda environment when they wish to use functionalities from BasinMaker.
            

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
      
      python -m pip install basinmaker
   
#. Install dependent packages 

    .. code-block::
    
      python -m pip install pandas pytest scipy simpledbf netCDF4 jupyter


#. Test validation 
     
    Please download the test data and scripts from `here <https://github.com/dustming/RoutingTool/wiki/Files/test.zip>`_. and unzip it to a folder, the path of this folder will refer as path_test_data in following section. Then
    
    Please ignore following output messages 
    
    * PyTables is not installed. No support for HDF output.
    * SQLalchemy is not installed. No support for SQL output.    
    * Warnings

    .. code-block::
     
       
      cd path_test_data/test
      python test_light_installation_arcgis.py
      (... some messages)
      ####################################
      BasinMaker is successfully installed
      ####################################

#. Users must active this conda environment when they wish to use functionalities from BasinMaker.


Full installation
==================

QGIS and GRASS in Windows
-------------------------

#. Installation of QGIS and GRASS using OSGEO4W: 
    
    For the Windows system, we can install both GRASS and QGIS within OSGEO4W environment.
    
    The OSGeo4W is a binary distribution of a broad set of open source geospatial software for Windows environments, including both GRASS GIS and QGIS.  
    
    The OSGeo4W installer can be downloaded from `here <https://qgis.org/en/site/forusers/download.html>`_. Please using OSGeo4W Network Installer (64 bit).
    
    We would suggest to 
    
    * Install QGIS and GRASS outside the **C/:Program Files**. Better to install them into a folder path without space in the folder name.
    * Use ‘Express Desktop Install’ 
    * Choose the default 3 packages
    * Run the downloaded installation file 
                
#. Setup GRASS and QGIS python environment

    The python environment for QGIS and GRASS GIS in Windows can be set up by modifying the following :download:`basinmaker.bat.txt <./_static/basinmaker.bat.txt>`.

    * Please rename 'basinmaker.bat.txt' to 'basinmaker.bat'.    
    * Please change OSGEO4W_ROOT to your OSGEO4W installation folder at line 2.
    * Please change the grass78.* in line 8 and 10 to your GRASS GIS version number.
    * Please double check the paths defined in the basinamker.bat file exists in your machine
    * Save the modified basinmaker.bat to a handy directory.  Run basinmaker.bat every time before using basinmaker.
    
#. Install BasinMaker (do not activate anaconda) 

    .. code-block::
      
      >basinmaker.bat
      Microsoft Windows [Version 10.0.19041.867]
      (c) 2020 Microsoft Corporation. All rights reserved
      >
      >python -m pip install basinmaker
    
#. Validate the GRASS and QGIS python environment
     
    * Please check if the python executable comes from the OSGeo4W64 installation folder
      by typing following commands after run basinmaker.bat. If the output is not 
      similar to the output showed in following output block. Please go back to step 2 and check
      the basinmaker.bat file  

    .. code-block::
       
      >where python    
      C:\OSGeo4W64\apps\Python37\python.exe

    * Check if all dependent QGIS and GRASS libraries can be imported in current python 
      environment by type following commands. 

    .. code-block::
      
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

      pip install simpledbf grass_session


#. Install GRASS GIS addons

    Install following GRASS GIS addons: 
    
    * r.accumulate
    * r.clip
    * r.stream.basins
    * r.stream.snap  
    
    How to install GRASS GIS addon can be found in `open GRASS GIS GUI <https://grass.osgeo.org/grass78/manuals/helptext.html>`_ and `add GRASS GIS addon <https://grass.osgeo.org/download/addons/>`_.  
  
#. Test validation 
     
    * Please download the test data and scripts from `here <https://github.com/dustming/RoutingTool/wiki/Files/test.zip>`_. and unzip it to a folder, the path of this folder will refer as path_test_data in following section. Then
    * run basinmaker.bat

    * PyTables is not installed. No support for HDF output.
    * SQLalchemy is not installed. No support for SQL output.    
    * Warnings
      
    .. code-block::
     
       
      cd path_test_data/test
      python test_full_installation.py
      (... some messages)
      ####################################
      BasinMaker is successfully installed
      ####################################
      
#. Users must run basinmaker.bat every time they wish to use functionalities from BasinMaker.

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

#. Setup GRASS and QGIS python environment

    The python environment for QGIS and GRASS GIS in Ubuntu can be set up by modifying the following :download:`basinmaker.sh <./_static/basinmaker.sh>`.
    
    * Please change the grass78.* in line 2 and 5 to your GRASS GIS version number.
    * Please double check the paths defined in the basinamker.sh file exists in your machine
    * Save the modified basinmaker.sh
    
#. Install BasinMaker 

    .. code-block::
      
      $source ./basinmaker.sh
      $Python3 -m pip install basinmaker
      
#. Validate the GRASS and QGIS python environment
     
    * Check if all dependent QGIS and GRASS libraries can be imported in current python 
      environment by type following commands.

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

      python3 -m pip install simpledbf grass_session

#. Install GRASS GIS addons

    Install following GRASS GIS addons: 
    
    * r.accumulate
    * r.clip
    * r.stream.basins
    * r.stream.snap  
    
    How to install GRASS GIS addon can be found in `open GRASS GIS GUI <https://grass.osgeo.org/grass78/manuals/helptext.html>`_ and `add GRASS GIS addon <https://grass.osgeo.org/download/addons/>`_. 

#. Test validation 
     
    * Please download the test data and scripts from `here <https://github.com/dustming/RoutingTool/wiki/Files/test.zip>`_. and unzip it to a folder, the path of this folder will refer as path_test_data in following section. Then
    * run basinmaker.sh

    * PyTables is not installed. No support for HDF output.
    * SQLalchemy is not installed. No support for SQL output.    
    * Warnings
        
    .. code-block::
     
       
      cd path_test_data/test
      python test_full_installation.py
      (... some messages)
      ####################################
      BasinMaker is successfully installed
      ####################################
          
  
#. Users must run basinmaker.sh every time they wish to use functionalities from BasinMaker.

