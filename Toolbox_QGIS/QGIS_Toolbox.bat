@echo off

SET OSGEO4W_ROOT=C:\PROGRA~1\QGIS3~1.12
SET GRASS_ROOT=C:\PROGRA~1\QGIS3~1.12\apps\grass\grass78
SET GISBASE= C:\PROGRA~1\QGIS3~1.12\apps\grass\grass78
SET GRASS_PYTHON=C:\PROGRA~1\QGIS3~1.12\apps\Python37
SET QGISPrefixPath=C:\PROGRA~1\QGIS3~1.12\apps\qgis
SET GRASSBIN=C:\PROGRA~1\QGIS3~1.12\bin\grass78.bat
SET RoutingToolFolder=C:\Users\m43han\Documents\GitHub\RoutingTool\Toolbox_QGIS

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

SET PYTHONHOME=C:\PROGRA~1\QGIS3~1.12\apps\Python37

SET PYTHONPATH=%PYTHONPATH%;%OSGEO4W_ROOT%\apps\qgis\python;%OSGEO4W_ROOT%\apps\qgis\python\plugins;%OSGEO4W_ROOT%\apps\qgis\python\plugins\processing
SET PYTHONPATH=%PYTHONPATH%;%GRASS_ROOT%\etc\python\;%GRASS_ROOT%\etc\python\grass;%GRASS_ROOT%\etc\python\grass\script;%RoutingToolFolder%;%RoutingToolFolder%\Toolbox_Common_Functions

rem Set VSI cache to be used as buffer, see #6448
set VSI_CACHE=TRUE
set VSI_CACHE_SIZE=1000000

set QT_PLUGIN_PATH=%OSGEO4W_ROOT%\apps\qgis\qtplugins;%OSGEO4W_ROOT%\apps\qt5\plugins
set QGIS_PREFIX_PATH=%OSGEO4W_ROOT:\=/%/apps/qgis

cmd.exe
