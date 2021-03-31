@echo off
set OSGEO4W_ROOT=C:\OSGeo4W64
call "%OSGEO4W_ROOT%\bin\o4w_env.bat"
call qt5_env.bat
call py3_env.bat

rem  for grass gis 
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
 

cmd.exe
