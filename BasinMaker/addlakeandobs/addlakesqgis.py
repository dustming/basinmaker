from processing_functions_raster_array import *
from processing_functions_raster_grass import *
from processing_functions_raster_qgis import *
from processing_functions_vector_grass import *
from processing_functions_vector_qgis import *
from utilities import *
from preprocessing.preprocessinglakeply import preprocessing_lake_polygon
import sqlite3
from addlakeandobs.definelaketypeqgis import define_connected_and_non_connected_lake_type
from addlakeandobs.filterlakesqgis import select_lakes_by_area_r
from addlakeandobs.pourpoints import define_pour_points_with_lakes

def add_lakes_into_existing_watershed_delineation(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    path_lakefile_in,
    lake_attributes,
    threshold_con_lake,
    threshold_non_con_lake,
    alllake = 'all_lakes',
    lake_boundary ='lake_boundary',
    connected_lake = 'connect_lake', 
    non_connected_lake = 'nonconnect_lake',
    str_connected_lake = 'str_connected_lake', 
    sl_connected_lake = 'sl_connected_lake',  
    sl_non_connected_lake = 'sl_nonconnect_lake', 
    sl_lakes = 'selected_lakes' ,
    sl_str_connected_lake = 'str_sl_connected_lake',
    max_memroy = 1024*4,
):

    fdr_arcgis=input_geo_names['fdr_arcgis']
    fdr_grass=input_geo_names['fdr_grass']
    str_r=input_geo_names['str_r']
    str_v=input_geo_names['str_v']
    acc=input_geo_names['acc']
    cat_no_lake=input_geo_names['cat_no_lake']
    mask = input_geo_names['mask']
    dem = input_geo_names['dem']
                    
    #prepropessing lakes inputs 
    preprocessing_lake_polygon(
        path_lakefile_in=path_lakefile_in,
        lake_attributes=lake_attributes,
        mask=mask,
        grassdb=grassdb,
        grass_location=grass_location,
        qgis_prefix_path=qgis_prefix_path,
        gis_platform="qgis",
        lake_name=alllake,
        lake_boundary_name = lake_boundary
    )
    
    import grass.script as grass
    import grass.script.setup as gsetup
    from grass.pygrass.modules import Module
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.script import array as garray
    from grass.script import core as gcore
    from grass_session import Session

    os.environ.update(
        dict(GRASS_COMPRESS_NULLS="1", GRASS_COMPRESSOR="ZSTD", GRASS_VERBOSE="1")
    )
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")
    
    con = sqlite3.connect(os.path.join(
        grassdb, grass_location, "PERMANENT", "sqlite", "sqlite.db"
    ))
    
    write_grass_and_arcgis_fdr_rules(grassdb)
    exp = "%s = int(%s)"%(alllake,alllake)
    grass_raster_r_mapcalc(grass, exp)
    exp = "%s = int(%s)"%(lake_boundary,lake_boundary)
    grass_raster_r_mapcalc(grass, exp)
    

    #Define connected and non connected lakes and
    #identify which str make certain lake have two outlet
    define_connected_and_non_connected_lake_type(
        grass,
        con,
        garray,
        str_r=str_r,
        lake=alllake,
        connected_lake = connected_lake,
        non_connected_lake = non_connected_lake,
        str_connected_lake = str_connected_lake,
    )
    
    select_lakes_by_area_r(
        grass=grass,
        con=con,
        lake_v_path = os.path.join(grassdb,alllake + ".shp"),
        threshold_con_lake = 1,
        threshold_non_con_lake = 0,
        lakes = alllake,
        connected_lake = connected_lake,
        non_connected_lake = non_connected_lake,
        str_connected_lake = str_connected_lake,
        sl_connected_lake = sl_connected_lake,
        sl_non_connected_lake = sl_non_connected_lake,
        sl_lakes = sl_lakes,
        sl_str_connected_lake = sl_str_connected_lake,
    )
    
    define_pour_points_with_lakes(
        grass = grass,
        con = con,
        garray=garray,
        str_r=str_r,
        cat_no_lake = cat_no_lake,
        sl_lakes = sl_lakes,
        sl_connected_lake = sl_connected_lake, 
        sl_str_connected_lake = sl_str_connected_lake,
        acc =acc,
        lake_pourpoints = "lake_pourpoints",
    )    



    PERMANENT.close()
    return
