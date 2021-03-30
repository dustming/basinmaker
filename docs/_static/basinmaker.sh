#!/bin/sh
export GISBASE='/usr/lib/grass78'
export QGISPrefixPath='/usr'

export PYTHONPATH=$PYTHONPATH:'/usr/lib/grass78/etc/python'  ### folder has a grass folder
export PYTHONPATH=$PYTHONPATH:'/usr/share/qgis/python/plugins' ## folder has db_manager and processing
export PYTHONPATH=$PYTHONPATH:'/usr/share/qgis/python' ## folder has plugin and console

