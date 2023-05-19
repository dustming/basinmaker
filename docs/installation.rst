===================
Installation v3.0.3
===================


Overview
========

BasinMaker is a python package work on several existing GIS platforms. So, the installation of BasinMaker includes two steps: 1) setup the python environment for the dependent GIS platforms; and 2) install BasinMaker itself.

Two installation modes (light installation and full installation) are available. The light installation will allow user to use BasinMaker to post process an existing routing product, such as the Basinmaker derived North American lake-river routing product. But it cannot be used to delineate a lake-river routing structure from DEM. The combination of BasinMaker light installation and the North American lake-river routing product could generate lake-river routing structures satisfying many user modeling demands. While the full installation of BasinMaker enables users to delineate a new lake-river routing structure from a user specified DEM.

For light installation (recommended unless users know for sure they need to delineate their watershed from scratch from a DEM), only geopandas or ArcGIS pro is needed. The python environment for both geopandas can be easily compiled within the anaconda environment under windows OS environment. The instruction about light installation procedure can be found in :ref:`Light installation`. Note that users can install both the light version and the full version on the same OS system.

Users, who want to use BasinMaker post processing functions without a Windows OS system, could try BasinMaker on Google Colab in :ref:`BasinMaker on Google Colab`

For full installation, BasinMaker can work under two python enviroments: the ArcGIS Pro and GRASS/QGIS GIS python enviroments. The instruction to install BasinMaker for each python enviroments are provided at :ref:`Full installation`.It is quite a challenge to setup a python environment for QGIS and GRASS together. Here, two procedures are provided for Windows OS systems, respectively. The procedure have been tested on several machines. But We can't guarantee install procedures work on every machin, but if you run into a problem create an issue on the GitHub and time permitting, we will try to help. If you managed to do the full installation on a different operating system, we would be grateful if you could document and share the detailed installation procedure that was successful (email: m43han@uwaterloo.ca).


Updating BasinMaker to v3.0
===========================
For existing users of v2.0 who want to update BasinMaker to the version v3.0.0. Users should simply reinstall BasinMaker 3.0 to a new working environment (e.g., called 'Basinmaker3' perhaps).


BasinMaker on Google Colab
==========================

A post-processing example via google colab (no installation on your local machine necessary!) can be found at here `here <https://colab.research.google.com/drive/14OC8l4ZeabOGGi0bL0ZFK1QzTOY8M9yM?usp=sharing>`_. The google colab is an online python notebook dose not require installation. This example will show you how to discretize, simplify, and revise the provided routing product for your purposes.


Light installation
==================

Geopandas with anaconda
-----------------------

#. Install anaconda

    The installer of anaconda can be installed from `here <https://www.anaconda.com/>`_. Note for windows system, please activate the 'Register Anaconda3 as my default python 3.9'

#. Create an empty python environment and then active it

    For windows system, search and open "Anacoda Prompt" (Windows) to active a conda command line. **Users must make sure:**

    * They have the proper privileges to create environment variables (e.g., run Anacoda Prompt as administrator will work)

    * DO NOT USE Anaconda Powershell Prompt

    Then

    .. code-block::

      conda create --name <any_name_for_env> python==3.11.3  
      conda activate <any_name_for_env>

#. Install packages (may take some time)

    .. code-block::

      conda config --append channels conda-forge

    .. code-block::

      conda install gdal

    .. code-block::

      conda install numpy=1.21.6 pandas==1.3.5 rasterstats==0.17.0 geopandas==0.11.0 GDAL==3.5.1 Shapely==1.8.5 rasterio==1.2.10

    .. code-block::

      python -m pip install  pytest  simpledbf netCDF4 joblib jupyter requests wget gdown pandas rasterstats geopandas rasterio scipy

#. Install BasinMaker

    .. code-block::

      python -m pip install basinmaker

#. Test validation

    Please download the test data and scripts from `here <https://github.com/dustming/RoutingTool/wiki/Files/test.zip>`_. and unzip it to a folder, the path of this folder will refer as path_test_data in following section. Then

    Please ignore following output messages

        * PyTables is not installed. No support for HDF output.
        * SQLalchemy is not installed. No support for SQL output.
        * Warnings

    .. code-block::


      cd path_test_data/test
      python test_light_installation_qgis.py
      (... some messages)
      ####################################
      BasinMaker is successfully installed
      ####################################

#. Users must active this conda environment when they wish to use functionalities from BasinMaker.

Full installation
==================

    The BasinMaker watershed delineation mode can be used under both ArcGIS Pro and GRASS GIS environments. We recommend using BasinMaker under the ArcGIS Pro environment. However, installation instructions for both Python environments are provided in the following two sections.


ArcGIS Pro in Windows
-------------------------

#. Install ArcGIS Pro

    BasinMaker has been tested with ArcGIS Pro version 3.0.3. To use the software, please ensure that you have installed this version of ArcGIS Pro. If you need assistance with installing ArcGIS Pro, please contact your IT department for detailed instructions.

#. Setup the python environment for BasinMaker in ArcGIS Pro

    Below are the key steps to create an ArcGIS Python environment. For detailed instructions, please refer to this `link <https://pro.arcgis.com/en/pro-app/2.9/arcpy/get-started/work-with-python-environments.htm>`_ .

    * Open ArcGIS Pro and click on the "Settings/Project" icon in the upper left corner of the ArcGIS Pro window.
    * Click the "Package Manager" tab. And then click the "Manage Environments" button in the upper right corner of the window.
    * Clone the "ArcGIS Pro" environment and name it <any_name_for_env>. The clone process will take a few minutes.
    * Select and active the newly created environment and restart ArcGIS Pro for the changes to take effect

#. Install BasinMaker in ArcGIS Pro

    * Open the ArcGIS Pro Python command prompt. To open the ArcGIS Pro Python command prompt, navigate to the Windows program directory: Programs > ArcGIS > Python Command Prompt.
    * Install BasinMaker and related pacakges using the following command:

    .. code-block::

      > python -m pip install basinmaker


#. Install dependent packages

    .. code-block::

      > python -m pip install pytest simpledbf netCDF4 joblib jupyter requests wget gdown geopandas rasterstats

#. Validate the installation with the package of test files.

    * Please download the `test data and script <https://github.com/dustming/RoutingTool/wiki/Files/test_arcgis_full.zip>`_ and unzip it to a folder, the path of this folder will refer as path_test_data in following section. Then
    * Open the ArcGIS Pro Python command prompt and run the following command:

    .. code-block::

      cd path_test_data
      python test_full_delineation_arcgis.py
      (... some messages)
      ####################################
      BasinMaker is successfully installed
      ####################################

    * Please ignore following output messages

      PyTables is not installed. No support for HDF output.

      SQLalchemy is not installed. No support for SQL output.

      Warnings

#. Users must Open the ArcGIS Pro Python command prompt every time they wish to use functionalities from BasinMaker.


QGIS and GRASS in Windows
-------------------------

#. Installation of QGIS and GRASS using OSGEO4W:

    For the Windows system, we can install both GRASS and QGIS within OSGEO4W environment.

    The OSGeo4W is a binary distribution of a broad set of open source geospatial software for Windows environments, including both GRASS GIS and QGIS.

    The OSGeo4W installer can be downloaded from `here <https://qgis.org/en/site/forusers/download.html>`_.

    Please use the advanced install option and keep the default selection in all pop up pages, except in the 'select package page'.


    In the select package:


    * In the Desktop group, please select 1) grass: GRASS GIS 7.8; 2) qgis: QGIS DESKTOP; 3)qt5_tools:Qt5 tools (development); 4)saga:SAGA(...)


    * In the Libs group, please select 1)python3-geopandas; 2)python3-rtree; 3)python3-rasterstats


    We would suggest to

    * Install QGIS and GRASS outside the **C/:Program Files**. Better to install them into a folder path without space in the folder name.
    * Run the downloaded installation file

#. Setup GRASS and QGIS python environment

    The python environment for QGIS and GRASS GIS in Windows can be set up by modifying the following :download:`basinmaker.bat.txt <./_static/basinmaker.bat.txt>`.

    * Please rename 'basinmaker.bat.txt' to 'basinmaker.bat'.
    * Please change OSGEO4W_ROOT to your OSGEO4W installation folder at line 2.
    * Please change the grass78.* in line 8 and 10 to your GRASS GIS version number.
    * Please double check the paths defined in the basinmaker.bat file exist in your machine
    * Save the modified basinmaker.bat to a handy directory.  Run basinmaker.bat every time before using basinmaker.

#. Install BasinMaker (do not activate anaconda)

    .. code-block::

      >basinmaker.bat
      Microsoft Windows [Version 10.0.19041.867]
      (c) 2020 Microsoft Corporation. All rights reserved
      >
      >python -m pip install basinmaker

#. Validate the GRASS and QGIS python environment

    * Please check if the python executable comes from the OSGeo4W installation folder
      by typing following commands after run basinmaker.bat. If the output is not
      similar to the output showed in following output block. Please go back to step 2 and check
      the basinmaker.bat file

    .. code-block::

      >where python
      C:\OSGeo4W\apps\Python37\python.exe

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
      >>>quit()


#. Install dependent packages

    .. code-block::

      python -m pip install simpledbf grass_session scipy joblib
      python -m pip install --upgrade pip
      python -m pip install geopandas -U
      python -m pip install numpy -U

#. Install GRASS GIS addons
    Install following GRASS GIS addons:

    * r.accumulate
    * r.clip
    * r.stream.basins
    * r.stream.snap

    For new GRASS users, see how to install GRASS GIS addon `here <https://github.com/dustming/RoutingTool/wiki/Files/GRASS_GIS_Addons_Install_Instruction.pdf>`_.

    If you want to learn how to use GRASS for more than BasinMaker, `this site <https://grass.osgeo.org/download/addons/>`_.  may help you.

#. Test validation

    * Please download the test data and scripts from `here <https://github.com/dustming/RoutingTool/wiki/Files/test.zip>`_. and unzip it to a folder, the path of this folder will refer as path_test_data in following section. Then
    * run basinmaker.bat
    * Please ignore following output messages

        PyTables is not installed. No support for HDF output.

        SQLalchemy is not installed. No support for SQL output.

        Warnings

    .. code-block::


      cd path_test_data/test
      python test_full_installation.py
      (... some messages)
      ####################################
      BasinMaker is successfully installed
      ####################################

#. Users must run basinmaker.bat every time they wish to use functionalities from BasinMaker.
