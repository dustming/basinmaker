import os

import numpy as np
import pandas as pd
import qgis
from qgis.analysis import QgsNativeAlgorithms
from qgis.core import *
from qgis.PyQt.QtCore import *
from joblib import Parallel, delayed
from basinmaker.utilities.utilities import *
import tempfile
from json import load, JSONEncoder
import json
import shutil

def qgis_raster_gdal_warpreproject(processing, Input, TARGET_CRS, Output):
    """Functions reproject raster layer
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    out = processing.run(
        "gdal:warpreproject",
        {
            "INPUT": Input,
            "SOURCE_CRS": None,
            "TARGET_CRS": QgsCoordinateReferenceSystem(TARGET_CRS),
            "RESAMPLING": 0,
            "NODATA": None,
            "TARGET_RESOLUTION": None,
            "OPTIONS": "",
            "DATA_TYPE": 0,
            "TARGET_EXTENT": None,
            "TARGET_EXTENT_CRS": None,
            "MULTITHREADING": False,
            "EXTRA": "",
            "OUTPUT": Output,
        },
    )
    return out

####
def create_geo_jason_file(processing,Input_Polygon_path):


    product_dir = os.path.dirname(Input_Polygon_path)
    Names_in = os.path.basename(Input_Polygon_path).split('_')
    n_charc = len(Names_in)
    version  = Names_in[n_charc - 1][0:4]
    TOLERANCEs = [0.0001,0.0005,0.001,0.005,0.01,0.05]


    head_name_cat = "finalcat_info"
    head_name_riv = "finalcat_info_riv"
    head_name_slake = "sl_connected_lake"
    head_name_nlake = "sl_non_connected_lake"

    Input_file_name = []
    Output_file_name = []                            
    if 'v' in version:
        Input_file_name = [
                           head_name_cat + "_"+version+'.shp',
                           head_name_riv + "_"+version+'.shp',
                           head_name_slake + "_"+version+'.shp',
                           head_name_nlake + "_"+version+'.shp',
                           ]
        Output_file_name = [
                           head_name_cat + "_"+version+'.geojson',
                           head_name_riv + "_"+version+'.geojson',
                           head_name_slake + "_"+version+'.geojson',
                           head_name_nlake + "_"+version+'.geojson',
                           ]
    else:
        Input_file_name = [
                           head_name_cat +'.shp',
                           head_name_riv +'.shp',
                           head_name_slake +'.shp',
                           head_name_nlake +'.shp',
                           ]
        Output_file_name = [
                           head_name_cat +'.geojson',
                           head_name_riv +'.geojson',
                           head_name_slake +'.geojson',
                           head_name_nlake +'.geojson',
                           ]
    created_jason_files = []   
    created_jason_files_lake_riv = [] 
    
    for i  in  range(0,len(Input_file_name)):
        input_path = os.path.join(product_dir,Input_file_name[i])
        output_jason_path = os.path.join(product_dir,Output_file_name[i])
        if not os.path.exists(input_path):
            continue 
        created_jason_files.append(output_jason_path) 
        
        if 'finalcat_info_riv' in Input_file_name[i] or 'connected_lake' in Input_file_name[i]:
            created_jason_files_lake_riv.append(output_jason_path)  
                            
        # reproject to WGS84
        input_wgs_84 = processing.run("native:reprojectlayer", {'INPUT':input_path,
                                      'TARGET_CRS':QgsCoordinateReferenceSystem('EPSG:4326'),
                                       'OUTPUT':"memory:"})['OUTPUT']  

        # input_wgs_84_smth =  processing.run("native:smoothgeometry", {'INPUT':input_wgs_84,
        #                                                         'ITERATIONS':5,'OFFSET':0.5,'MAX_ANGLE':180,
        #                                                         'OUTPUT':"memory:"})['OUTPUT']

        
        if 'finalcat_info' in Input_file_name[i] or "finalcat_info_riv" in Input_file_name[i]:
            input_tojson = processing.run("native:fieldcalculator", {'INPUT':input_wgs_84,
                                                      'FIELD_NAME':'rvhName','FIELD_TYPE':2,'FIELD_LENGTH':20,'FIELD_PRECISION':0,
                                                      'FORMULA':'\'sub\' +to_string(to_int(\"SubId\")) ','OUTPUT':"memory:"})['OUTPUT']
        else:
             input_tojson = input_wgs_84
        
        for TOLERANCE in TOLERANCEs:                               
            input_wgs_84_simplify = processing.run("native:simplifygeometries", {
                                                   'INPUT':input_tojson,
                                                   'METHOD':0,
                                                   'TOLERANCE':TOLERANCE,
                                                   'OUTPUT':"memory:"
                                                   }
                                                   )['OUTPUT']
                                               
            qgis.core.QgsVectorFileWriter.writeAsVectorFormat(input_wgs_84_simplify,output_jason_path, 'utf-8', input_wgs_84_simplify.crs(), 'GeoJson',layerOptions=['COORDINATE_PRECISION=3'])
            
            json_file_size = os.stat(output_jason_path).st_size/1024/1024 #to MB
            if json_file_size <= 100:
                break
    # if len(created_jason_files) > 1 and os.stat(os.path.join(product_dir,Output_file_name[0])).st_size/1024/1024 < 500::
    #     for i in range(0,len(created_jason_files)):
    #         injson = load(open(created_jason_files[i]))
    #         if i == 0:
    #             output_jason = injson
    #         else:
    #             output_jason['features'] += injson['features']
    # 
    #     with open(os.path.join(product_dir,'routing_product.geojson'), 'w', encoding='utf-8') as f:
    #         json.dump(output_jason, f, ensure_ascii=False, indent=4)
    # else:
    #     shutil.copy(created_jason_files[0], os.path.join(product_dir,'routing_product.geojson')) 
        
    if len(created_jason_files_lake_riv) > 1 and os.stat(os.path.join(product_dir,Output_file_name[0])).st_size/1024/1024 < 500:
        for i in range(0,len(created_jason_files_lake_riv)):
            injson2 = load(open(created_jason_files_lake_riv[i]))
            if 'finalcat_info_riv' in created_jason_files_lake_riv[i]:
                new_features = []
                for element in injson2["features"]:
                    if element["properties"]["Lake_Cat"] == 0:    
                        new_features.append(element) 
                injson2["features"] = new_features

            if i == 0:
                output_jason_lake_riv = injson2
            else:
                output_jason_lake_riv['features'] += injson2['features']
                
        with open(os.path.join(product_dir,'routing_product_lake_river.geojson'), 'w', encoding='utf-8') as f:
            json.dump(output_jason_lake_riv, f, ensure_ascii=False, indent=4)
    else:
        shutil.copy(created_jason_files_lake_riv[0], os.path.join(product_dir,'routing_product_lake_river.geojson'))             
    return 
    
##############


def qgis_raster_clip_raster_by_mask(processing, Input, MASK, TARGET_CRS, Output):
    """Functions clip raster by mask
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    out = processing.run(
        "gdal:cliprasterbymasklayer",
        {
            "INPUT": Input,
            "MASK": MASK,
            "SOURCE_CRS": None,
            "TARGET_CRS": QgsCoordinateReferenceSystem(TARGET_CRS),
            "NODATA": None,
            "ALPHA_BAND": False,
            "CROP_TO_CUTLINE": True,
            "KEEP_RESOLUTION": False,
            "SET_RESOLUTION": False,
            "X_RESOLUTION": None,
            "Y_RESOLUTION": None,
            "MULTITHREADING": False,
            "OPTIONS": "",
            "DATA_TYPE": 0,
            "EXTRA": "",
            "OUTPUT": Output,
        },
    )
    return out


def qgis_raster_clip_raster_by_extent(processing, Input, PROJWIN, Output):
    """Functions clip raster by extent
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """
    out = processing.run(
        "gdal:cliprasterbyextent",
        {
            "INPUT": Input,
            "PROJWIN": PROJWIN,
            "NODATA": None,
            "OPTIONS": "",
            "DATA_TYPE": 0,
            "EXTRA": "",
            "OUTPUT": Output,
        },
    )

    return out


def qgis_raster_saga_clip_raster_with_polygon(processing, context, Input, MASK, Output):
    """Functions clip raster by mask,
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    out = processing.run(
        "saga:cliprasterwithpolygon",
        {"INPUT": Input, "POLYGONS": MASK, "OUTPUT": Output},
        context=context,
    )

    return out


def qgis_raster_slope(processing, Input, Output):
    """Functions calculate slope from DEM
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    out = processing.run(
        "qgis:slope", {"INPUT": Input, "Z_FACTOR": 1, "OUTPUT": Output}
    )
    return out


def qgis_raster_aspect(processing, Input, Output):
    """Functions calculate aspect from DEM
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    out = processing.run(
        "qgis:aspect", {"INPUT": Input, "Z_FACTOR": 1, "OUTPUT": Output}
    )
    return out


def qgis_raster_zonal_statistics(processing, INPUT_RASTER, INPUT_VECTOR, COLUMN_PREFIX):
    """Functions calculate zonal statistics between polygon and raster
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    processing.run(
        "qgis:zonalstatistics",
        {
            "INPUT_RASTER": INPUT_RASTER,
            "RASTER_BAND": 1,
            "INPUT_VECTOR": INPUT_VECTOR,
            "COLUMN_PREFIX": COLUMN_PREFIX,
            "STATS": [2],
        },
    )
    return


def qgis_raster_read_raster(processing, INPUT):
    """qgis function read raster into memeory
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = QgsRasterLayer(INPUT, "")
    return out


def qgis_raster_gdal_translate(processing, INPUT, OUTPUT, format="GTiff"):
    """qgis function read raster into memeory
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    params = {"INPUT": INPUT, "format": format, "OUTPUT": OUTPUT}
    out = processing.run("gdal:translate", params)
    return out


def qgis_raster_return_raster_properties(processing, INPUT):
    """qgis function return raster projection and cellsize
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    cellSize = float(INPUT.rasterUnitsPerPixelX())  ### Get Raster cell size
    SpRef_in = INPUT.crs().authid()  ### get Raster spatialReference id
    return cellSize, SpRef_in


def qgis_raster_gdal_polygonize(processing, context, INPUT, OUTPUT):
    """qgis function polygonize rater to polygon
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run(
        "gdal:polygonize",
        {
            "INPUT": INPUT,
            "BAND": 1,
            "FIELD": "DN",
            "EIGHT_CONNECTEDNESS": False,
            "EXTRA": "",
            "OUTPUT": OUTPUT,
        },
        context=context,
    )

    return out


def qgis_raster_gdal_rasterize(
    processing, context, INPUT, Column_nm, cellsize, w, s, e, n, OUTPUT
):
    """qgis gdal function polygonize rater to polygon
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    cmd = (
        "gdal_rasterize -at -of GTiff -a_nodata -9999 -a "
        + Column_nm
        + " -tr  "
        + str(cellsize)
        + "  "
        + str(cellsize)
        + "  -te   "
        + str(w)
        + "   "
        + str(s)
        + "   "
        + str(e)
        + "   "
        + str(n)
        + "   "
        + '"'
        + INPUT
        + '"'
        + "    "
        + '"'
        + OUTPUT
        + '"'
    )
    os.system(cmd)
    return


# os.system('gdal_rasterize -at -of GTiff -a_nodata -9999 -a Hylak_id -tr  '+ str(self.cellSize) + "  " +str(self.cellSize)+'  -te   '+ str(grsregion['w'])+"   " +str(grsregion['s'])+"   " +str(grsregion['e'])+"   " +str(grsregion['n'])+"   " + "\"" +  self.Path_allLakeply +"\""+ "    "+ "\""+self.Path_allLakeRas+"\"")

def Copy_Pddataframe_to_shpfile_by_file(
    Path_shpfile,
    Pddataframe,
    link_col_nm_shp="SubId",
    link_col_nm_df="SubId",
    UpdateColNM=["#"],
    Input_Is_Feature_In_Mem=False,
):

    if Input_Is_Feature_In_Mem:
        layer_cat = Path_shpfile
    else:
        layer_cat = QgsVectorLayer(Path_shpfile, "")
    
    # table to csv 
    Attri_table.to_csv(os.path.join(tempfile.gettempdir(), "attribute.csv"),index = False)
    csv_uri = 'file:///%s?delimiter=,'%(os.path.join(tempfile.gettempdir(), "attribute.csv"))
    csv = QgsVectorLayer(csv_uri, 'attributes', 'delimitedtext')

    
    


def Copy_Pddataframe_to_shpfile(
    Path_shpfile,
    Pddataframe,
    link_col_nm_shp="SubId",
    link_col_nm_df="SubId",
    UpdateColNM=["#"],
    Input_Is_Feature_In_Mem=False,
):
    """Function modify attribute table of Path_shpfile using value from Pddataframe
    Parameters
    ----------
    Path_shpfile                        : shpfile
        Path to the shpfile
    Pddataframe                         : dataframe
        Dataframe constains data that will be used to update attribute table of
        Path_shpfile
    link_col_nm_shp                     : string
        The column name that link Pddataframe attribute table in Path_shpfile
    link_col_nm_df                      : string
        The column name that link Path_shpfile attribute table in Pddataframe
    UpdateColNM                         : list
        It is a list of column name, it is equal to '#', all column value will
        be updated, otherwise only column name in UpdateColNM will be updated.
    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """
    if Input_Is_Feature_In_Mem:
        layer_cat = Path_shpfile
    else:
        layer_cat = QgsVectorLayer(Path_shpfile, "")
        
    Attri_Name = layer_cat.fields().names()
    features = layer_cat.getFeatures()
    with edit(layer_cat):
        for sf in features:
            Atti_Valu = sf.attributes()
            sf_subid = sf[link_col_nm_shp]
            tarinfo = Pddataframe[Pddataframe[link_col_nm_df] == sf_subid]
            
            if UpdateColNM[0] == "#":
                for icolnm in range(0, len(Attri_Name)):  ### copy infomaiton
                    if (
                        Attri_Name[icolnm] == "Obs_NM"
                        or Attri_Name[icolnm] == "SRC_obs"
                        or Attri_Name[icolnm] == "LAND_USE_C"
                        or Attri_Name[icolnm] == "SOIL_PROF"
                        or Attri_Name[icolnm] == "VEG_C"
                    ):  
                        if str(tarinfo[Attri_Name[icolnm]].values[0]) == 'nan':
                            sf[Attri_Name[icolnm]] = qgis.core.NULL
                        else:
                            sf[Attri_Name[icolnm]] = str(
                                tarinfo[Attri_Name[icolnm]].values[0]
                            )
                    elif (
                        Attri_Name[icolnm] == "cat"
                        or Attri_Name[icolnm] == "layer"
                        or Attri_Name[icolnm] == "path"
                    ):
                        continue
                    else:
                        sf[Attri_Name[icolnm]] = float(
                            tarinfo[Attri_Name[icolnm]].values[0]
                        )
            else:
                for icolnm in range(0, len(UpdateColNM)):
                    if (
                        UpdateColNM[icolnm] == "Obs_NM"
                        or UpdateColNM[icolnm] == "SRC_obs"
                        or UpdateColNM[icolnm] == "LAND_USE_C"
                        or UpdateColNM[icolnm] == "SOIL_PROF"
                        or UpdateColNM[icolnm] == "VEG_C"
                    ):
                        sf[UpdateColNM[icolnm]] = str(
                            tarinfo[UpdateColNM[icolnm]].values[0]
                        )
                    elif (
                        UpdateColNM[icolnm] == "cat"
                        or UpdateColNM[icolnm] == "layer"
                        or UpdateColNM[icolnm] == "path"
                    ):
                        continue
                    else:
                        sf[UpdateColNM[icolnm]] = float(
                            tarinfo[UpdateColNM[icolnm]].values[0]
                        )

            layer_cat.updateFeature(sf)
    if Input_Is_Feature_In_Mem:
        return layer_cat
    else:
        del layer_cat
        return

def copy_data_and_dissolve(all_subids,tempfolder,processing,Path_Temp_final_rviply,Path_Temp_final_rvi,
    mapoldnew_info,COLUMN_NAMES_CONSTANT_CLEAN,OutputFolder,Path_Catchment_Polygon,context,Path_final_rviply_input = '#',Path_final_rvi_input = '#'):    
    
    if len(all_subids) > 5000:

        Path_final_rviply = os.path.join(
            tempfolder,
            "temp_finalriv_ply"  + ".shp",
        )
        Path_final_rvi = os.path.join(
            tempfolder,
            "temp_finalriv"  + ".shp",
        )
    
        Copy_Pddataframe_to_shpfile_main(
            Path_shpfile = Path_Temp_final_rviply,
            Pddataframe = mapoldnew_info,
            all_subids = all_subids,
            tempfolder = tempfolder,
            processing = processing,
            output = Path_final_rviply,
            link_col_nm_shp="SubId",
            link_col_nm_df="Old_SubId",
            UpdateColNM=["#"],
        )  
        
        Copy_Pddataframe_to_shpfile_main(
            Path_shpfile = Path_Temp_final_rvi,
            Pddataframe = mapoldnew_info,
            all_subids = all_subids,
            tempfolder = tempfolder,
            processing = processing,  
            output =   Path_final_rvi,    
            link_col_nm_shp="SubId",
            link_col_nm_df="Old_SubId",
            UpdateColNM=["#"],
        )    
                              
    else:
            
        Path_final_rviply = Path_Temp_final_rviply
        Path_final_rvi = Path_Temp_final_rvi
                    
        Copy_Pddataframe_to_shpfile(
            Path_Temp_final_rviply,
            mapoldnew_info,
            link_col_nm_shp="SubId",
            link_col_nm_df="Old_SubId",
            UpdateColNM=["#"],
        )
        Copy_Pddataframe_to_shpfile(
            Path_Temp_final_rvi,
            mapoldnew_info,
            link_col_nm_shp="SubId",
            link_col_nm_df="Old_SubId",
            UpdateColNM=["#"],
        )

    # dissolve shpfile based on new subid
    if len(os.path.basename(Path_Catchment_Polygon).split('_')) == 5:
        Path_final_rviply_out = os.path.join(OutputFolder, "finalcat_info_"+os.path.basename(Path_Catchment_Polygon).split('_')[4])
        Path_final_rvi_out = os.path.join(OutputFolder, "finalcat_info_riv_"+os.path.basename(Path_Catchment_Polygon).split('_')[4])
    else:
        Path_final_rviply_out = os.path.join(OutputFolder, "finalcat_info.shp")
        Path_final_rvi_out = os.path.join(OutputFolder, "finalcat_info_riv.shp")

    if Path_final_rviply_input != '#' and Path_final_rvi_input != '#':
        Path_final_rviply_out = os.path.join(OutputFolder, os.path.basename(Path_final_rviply_input))
        Path_final_rvi_out = os.path.join(OutputFolder, os.path.basename(Path_final_rvi_input))
        
    COLUMN_NAMES_CONSTANT_CLEAN_RIV = COLUMN_NAMES_CONSTANT_CLEAN.copy()
    COLUMN_NAMES_CONSTANT_CLEAN_RIV.remove('centroid_x') 
    COLUMN_NAMES_CONSTANT_CLEAN_RIV.remove('centroid_y')  
    Clean_Attribute_Name(Path_final_rvi, COLUMN_NAMES_CONSTANT_CLEAN_RIV)
    Clean_Attribute_Name(Path_final_rviply, COLUMN_NAMES_CONSTANT_CLEAN)
    
                    
    riv_dissolve = qgis_vector_dissolve(
        processing,
        context,
        INPUT=Path_final_rvi,
        FIELD=["SubId"],
        OUTPUT="memory:",# Path_final_rvi_out,
    )["OUTPUT"]

    processing.run("native:extractbyexpression", {
                    'INPUT':riv_dissolve,
                    'EXPRESSION':' \"Lake_Cat\" > 0  OR  \"RivLength\" > 0',
                    'OUTPUT':Path_final_rvi_out}
                    )

    ply_draft = qgis_vector_dissolve(
        processing,
        context,
        INPUT=Path_final_rviply,
        FIELD=["SubId"],
        OUTPUT="memory:",
    )["OUTPUT"]
    
    # clean attribute table of shpfile
    Clean_Attribute_Name(Path_final_rvi_out, COLUMN_NAMES_CONSTANT_CLEAN)
    ply_draft,tempnum = Clean_Attribute_Name(ply_draft, COLUMN_NAMES_CONSTANT_CLEAN,Input_Is_Feature_In_Mem = True)

    # add centroid to new drived polygons
    
    formular = "x(centroid(transform($geometry,'%s','%s')))" % (
        ply_draft.crs().authid(),
        "EPSG:4326",
    )

    ply_draft = qgis_vector_field_calculator(
        processing=processing,
        context=context,
        FORMULA=formular,
        FIELD_NAME="centroid_x",
        INPUT=ply_draft,
        OUTPUT="memory:",
    )["OUTPUT"]

    formular = "y(centroid(transform($geometry,'%s','%s')))" % (
        ply_draft.crs().authid(),
        "EPSG:4326",
    )

    qgis_vector_field_calculator(
        processing=processing,
        context=context,
        FORMULA=formular,
        FIELD_NAME="centroid_y",
        INPUT=ply_draft,
        OUTPUT=Path_final_rviply_out,
    )
    if 'finalcat_info' in Path_final_rviply_out:
        create_geo_jason_file(processing,Path_final_rviply_out)
    return 

def Copy_Pddataframe_to_shpfile_main(
    Path_shpfile,
    Pddataframe,
    all_subids,
    tempfolder,
    processing,
    output,
    link_col_nm_shp="SubId",
    link_col_nm_df="SubId",
    UpdateColNM=["#"],
    Input_Is_Feature_In_Mem=False,
):

    path_to_temp_files = []
    ncores = 3
    feature_per_file = int(len(all_subids)/ncores)
    
    for i in range(0,ncores):
        num_random = str(np.random.randint(1, 10000 + 1)) +'_'+str(i)
        path_of_i_temp_file = os.path.join(
            tempfolder,
            "temp_finalriv_ply" + num_random + ".shp",
        )
        if i != ncores - 1:
            subids_of_i_temp_file = all_subids[i*feature_per_file:(i+1)*feature_per_file]
        else:
            subids_of_i_temp_file = all_subids[i*feature_per_file:len(all_subids)]
                
        Selectfeatureattributes(
            processing,
            Input=Path_shpfile,
            Output=path_of_i_temp_file,
            Attri_NM="SubId",
            Values=subids_of_i_temp_file,
        )
        path_to_temp_files.append(path_of_i_temp_file)

    Parallel(n_jobs=ncores, verbose=1, backend="threading")(
        delayed(Copy_Pddataframe_to_shpfile)(temp_file,Pddataframe,
        link_col_nm_shp="SubId",link_col_nm_df="Old_SubId",UpdateColNM=["#"]
        ) for temp_file in path_to_temp_files)
    
    processing.run("native:mergevectorlayers", {'LAYERS':path_to_temp_files,'CRS':None,'OUTPUT':output})
    
    return 
    

def Remove_Unselected_Lake_Attribute_In_Finalcatinfo(Path_Finalcatinfo, Conn_Lake_Ids):
    """Functions will set lake id not in Conn_Lake_Ids to -1.2345 in attribute
        table of Path_Finalcatinfo
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    layer_cat = QgsVectorLayer(Path_Finalcatinfo, "")
    Attri_Name = layer_cat.fields().names()
    features = layer_cat.getFeatures()
    with edit(layer_cat):
        for sf in features:
            sf_subid = float(sf["HyLakeId"])

            if sf_subid in Conn_Lake_Ids or float(sf["Lake_Cat"]) == 2:
                continue
            sf["HyLakeId"] = float(0)
            sf["LakeVol"] = float(0)
            sf["LakeArea"] = float(0)
            sf["LakeDepth"] = float(0)
            sf["Laketype"] = float(0)
            sf["Lake_Cat"] = float(0)
            layer_cat.updateFeature(sf)
    del layer_cat
    return


#########
def Add_centroid_to_feature(
    Path_feagure, centroidx_nm="#", centroidy_nm="#", Input_Is_Feature_In_Mem=False
):
    """Functions will add centorid x y to Path_feagure
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """
    if Input_Is_Feature_In_Mem:
        layer_cat = Path_feagure
    else:
        layer_cat = QgsVectorLayer(Path_feagure, "")
    Attri_Name = layer_cat.fields().names()
    features = layer_cat.getFeatures()
    with edit(layer_cat):
        for sf in features:
            centroidxy = sf.geometry().centroid().asPoint()
            sf[centroidx_nm] = centroidxy[0]
            sf[centroidy_nm] = centroidxy[1]
            layer_cat.updateFeature(sf)
    if Input_Is_Feature_In_Mem:
        return layer_cat
    else:
        del layer_cat


def Selectfeatureattributes(processing, Input="#", Output="#", Attri_NM="#", Values=[],Is_str=False):
    """Functions extract features from Input, based on values in column Attri_NM
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """
    if Is_str:
        
        #'  \"STATION_NU\" IN (\'02KB001\',\'02KB003\')'
        
        exp = "\"%s\" IN (" %(Attri_NM)
        exp = exp + "\'%s\'" %(Values[0])
        for i in range(1, len(Values)):
            exp =  exp + " , \'%s\'" %(Values[i])
        exp = exp + ")"
    else:
        ncores = 3
        nvalues_per_group = 1000
        k = 0
        exp = Attri_NM + "  IN  (  " + str(int(Values[0]))
        for i in range(1, len(Values)):
            k = k + 1
            if k != 0:
                exp = exp + " , " + str(int(Values[i]))
            if k == nvalues_per_group and k < len(Values) - 3:
                k = 0
                exp = exp + ")"
                exp = exp + ' OR ' + Attri_NM + "  IN  (  " + str(int(Values[i]))
        
        exp = exp + ")"
    
    processing.run(
        "native:extractbyexpression",
        {"INPUT": Input, "EXPRESSION": exp, "OUTPUT": Output},
    )


def Copyfeature_to_another_shp_by_attribute(
    Source_shp, Target_shp, Col_NM="SubId", Values=[-1], Attributes=[-1]
):

    """Functions that will copy features in Source_shp to Target_shp
    based on attribute values in Values
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    layer_src = QgsVectorLayer(Source_shp, "")
    layer_trg = QgsVectorLayer(Target_shp, "")

    src_features = layer_src.getFeatures()

    Selected_Features = []
    for sf in src_features:
        # centroidxy = sf.geometry().centroid().asPoint()
        Select_value = sf[Col_NM]
        if Select_value in Values:
            src_geometry = sf.geometry()
            attribute = Attributes.loc[Attributes[Col_NM] == Select_value].values
            temp_feature = QgsFeature()
            temp_feature.setGeometry(src_geometry)
            temp_feature.setAttributes(attribute.tolist()[0])
            Selected_Features.append(temp_feature)

    layer_trg.startEditing()
    layer_trg.addFeatures(Selected_Features)
    layer_trg.commitChanges()
    layer_trg.updateExtents()
    del layer_src
    del layer_trg


def Add_New_SubId_To_Subregion_shpfile(
    processing, context, Layer, SubID_info="#", OutputPath="#", Region_ID=1
):
    """Asign new subbasin Id to each subbasins in each subregion
    Parameters
    ----------
    processing                        : qgis object
    context                           : qgis object
    Layer                             : vector layer
        it is the subbasin polygon or polyline of watershed delineation
        result in each subregion
    SubID_info                        : dataframe
        it is a dataframe contains new subbasin id for each subbasin
        in Layer
    OutputPath                        : string
        Path to the output file
    Region_ID                         : integer
        it is the subregion id of layer

    Notes
    -------
        the output will be the same type of input Layer
        stored in OutputPath, a new SubId will be given
    Returns:
    -------
        None
    """

    layer2 = qgis_vector_field_calculator(
        processing=processing,
        context=context,
        FORMULA=str(int(Region_ID)),
        FIELD_NAME="Region_ID",
        INPUT=Layer,
        OUTPUT="memory:",
    )['OUTPUT']
    qgis_vector_field_calculator(
        processing=processing,
        context=context,
        FORMULA="0",
        FIELD_NAME="Use_region",
        INPUT=layer2,
        OUTPUT=OutputPath,
    )
    
    layer_new = QgsVectorLayer(OutputPath, "")

    features = layer_new.getFeatures()
    with edit(layer_new):
        for sf in features:
            cSubId = int(sf["SubId"])
            cDowSubId = int(sf["DowSubId"])
            nSubId = SubID_info.loc[SubID_info["SubId"] == cSubId]["nSubId"]
            if len(SubID_info.loc[SubID_info["SubId"] == cDowSubId]) == 0 or cSubId == cDowSubId:
                nDowSubId = -1
            else:
                nDowSubId = SubID_info.loc[SubID_info["SubId"] == cDowSubId]["nSubId"]
            nSeg_ID = SubID_info.loc[SubID_info["SubId"] == cSubId]["nSeg_ID"]
            sf["SubId"] = int(nSubId)
            sf["DowSubId"] = int(nDowSubId)
            sf["Seg_ID"] = int(nSeg_ID)
            layer_new.updateFeature(sf)
    del layer_new
    return


def qgis_vector_field_calculator(
    processing,
    context,
    FORMULA,
    INPUT,
    OUTPUT,
    FIELD_NAME,
    FIELD_PRECISION=0,
    FIELD_TYPE=0,
    NEW_FIELD=True,
    FIELD_LENGTH=10,
):
    """qgis filed calcuator
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    alg_params = {
        "FIELD_LENGTH": FIELD_LENGTH,
        "FIELD_NAME": FIELD_NAME,
        "FIELD_PRECISION": FIELD_PRECISION,
        "FIELD_TYPE": FIELD_TYPE,
        "FORMULA": FORMULA,
        "INPUT": INPUT,
        "NEW_FIELD": NEW_FIELD,
        "OUTPUT": OUTPUT,
    }

    out = processing.run("qgis:fieldcalculator", alg_params, context=context)
    return out


def qgis_vector_dissolve(
    processing, context, INPUT, FIELD, OUTPUT, USING_GDAL_FUNCTION=False
):
    """qgis dissolve input vector based on values in FIELD list
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    if USING_GDAL_FUNCTION:
        out = processing.run(
            "gdal:dissolve",
            {"INPUT": INPUT, "FIELD": FIELD, "OUTPUT": OUTPUT},
            context=context,
        )
    else:
        out = processing.run(
            "native:dissolve",
            {"INPUT": INPUT, "FIELD": FIELD, "OUTPUT": OUTPUT},
            context=context,
        )
    return out


def qgis_vector_fix_geometries(processing, context, INPUT, OUTPUT):
    """qgis fixgeometries
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run("native:fixgeometries", {"INPUT": INPUT, "OUTPUT": OUTPUT})
    return out


def qgis_vector_merge_vector_layers(processing, context, INPUT_Layer_List, OUTPUT):
    """qgis merge_vector_layers
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run(
        "native:mergevectorlayers", {"LAYERS": INPUT_Layer_List, "OUTPUT": OUTPUT}
    )

    return out


def qgis_vector_return_crs_id(
    processing, context, INPUT_Layer, Input_Is_Feature_In_Mem=True
):
    """qgis return vector layer projection crs id
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    if Input_Is_Feature_In_Mem:
        out = INPUT_Layer.crs().authid()
    else:
        layer = QgsVectorLayer(INPUT_Layer, "")
        out = layer.crs().authid()

    return out


def qgis_vector_add_polygon_attribute_to_points(
    processing, context, INPUT_Layer, POLYGONS, FIELDS, OUTPUT
):
    """qgis return vector layer projection crs id
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run(
        "saga:addpolygonattributestopoints",
        {
            "INPUT": INPUT_Layer,
            "POLYGONS": POLYGONS,
            "FIELDS": FIELDS,
            "OUTPUT": OUTPUT,
        },
    )

    return out


def qgis_vector_extract_by_attribute(
    processing, context, INPUT_Layer, FIELD, OPERATOR, VALUE, OUTPUT
):
    """qgis extract vector by attribute
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run(
        "native:extractbyattribute",
        {
            "INPUT": INPUT_Layer,
            "FIELD": FIELD,
            "OPERATOR": OPERATOR,
            "VALUE": VALUE,
            "OUTPUT": OUTPUT,
        },
    )

    return out


def qgis_vector_add_attributes(processing, context, INPUT_Layer, attribute_list):
    """qgis add attributes to vector
    ----------
    Notes
    -------

    Returns:
    -------
        None,
    """
    INPUT_Layer.dataProvider().addAttributes(attribute_list)
    INPUT_Layer.updateFields()
    INPUT_Layer.commitChanges()
    return INPUT_Layer


def qgis_vector_get_attributes(processing, context, INPUT_Layer, attribute_NM):
    """qgis retrun attribute value of vector
    ----------
    Notes
    -------

    Returns:
    -------
        None,
    """
    if attribute_NM == "count":
        out = INPUT_Layer.featureCount()
        return out
    elif attribute_NM == "field_name":
        out = INPUT_Layer.fields().names()
    elif attribute_NM == "features":
        out = INPUT_Layer.getFeatures()
        return out
    else:
        print("wrong attribute name")
        return -1


####


def qgis_vector_union_two_layers(
    processing, context, INPUT, OVERLAY, OUTPUT, OVERLAY_FIELDS_PREFIX=""
):
    """qgis union two layers
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """

    out = processing.run(
        "native:union",
        {
            "INPUT": INPUT,
            "OVERLAY": OVERLAY,
            "OVERLAY_FIELDS_PREFIX": OVERLAY_FIELDS_PREFIX,
            "OUTPUT": OUTPUT,
        },
        context=context,
    )

    return out


def qgis_vector_reproject_layers(processing, context, INPUT, TARGET_CRS, OUTPUT):
    """qgis function reproject vector layer
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run(
        "native:reprojectlayer",
        {
            "INPUT": INPUT,
            "TARGET_CRS": QgsCoordinateReferenceSystem(TARGET_CRS),
            "OUTPUT": OUTPUT,
        },
    )
    return out


def qgis_vector_buffer(processing, context, INPUT, Buffer_Distance, OUTPUT):
    """qgis function buffer vector layer with a buffer distance
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run(
        "native:buffer",
        {
            "INPUT": INPUT,
            "DISTANCE": Buffer_Distance,
            "SEGMENTS": 5,
            "END_CAP_STYLE": 0,
            "JOIN_STYLE": 0,
            "MITER_LIMIT": 2,
            "DISSOLVE": True,
            "OUTPUT": OUTPUT,
        },
        context=context,
    )
    return out


def qgis_vector_ectract_by_location(processing, context, INPUT, INTERSECT, OUTPUT):
    """qgis function extract vector by location
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run(
        "native:extractbylocation",
        {"INPUT": INPUT, "PREDICATE": [6], "INTERSECT": INTERSECT, "OUTPUT": OUTPUT},
        context=context,
    )
    return out


def qgis_vector_polygon_stro_lines(processing, context, INPUT, OUTPUT):
    """qgis function to obtain polygon boundary lines
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run(
        "native:polygonstolines", {"INPUT": INPUT, "OUTPUT": OUTPUT}, context=context
    )
    return out


def qgis_vector_clip(processing, context, INPUT, OVERLAY, OUTPUT):
    """qgis function reproject vector layer
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    layer_clip = processing.run(
        "native:clip", {"INPUT": INPUT, "OVERLAY": OVERLAY, "OUTPUT": OUTPUT}
    )

    return layer_clip


def qgis_vector_join_attribute_table(
    processing,
    context,
    INPUT1,
    FIELD1,
    INPUT2,
    FIELD2,
    OUTPUT,
    FIELDS_TO_COPY=[],
    METHOD=1,
    DISCARD_NONMATCHING=False,
    PREFIX="",
):
    """qgis function join attibute table
    ----------

    Notes
    -------

    Returns:
    -------
    """

    out = processing.run(
        "native:joinattributestable",
        {
            "INPUT": INPUT1,
            "FIELD": FIELD1,
            "INPUT_2": INPUT2,
            "FIELD_2": FIELD2,
            "FIELDS_TO_COPY": FIELDS_TO_COPY,
            "METHOD": METHOD,
            "DISCARD_NONMATCHING": DISCARD_NONMATCHING,
            "PREFIX": PREFIX,
            "OUTPUT": OUTPUT,
        },
        context=context,
    )
    #          processing.run("native:joinattributestable", {'INPUT':HRU_draft,'FIELD':Sub_ID,'INPUT_2':Path_Subbasin_Ply,'FIELD_2':Sub_ID,
    #                    'FIELDS_TO_COPY':[],'METHOD':1,'DISCARD_NONMATCHING':False,'PREFIX':'',
    #                    'OUTPUT':'memory:'},context = context)['OUTPUT']

    return out


def qgis_vector_create_spatial_index(processing, context, INPUT):
    """qgis function create spatial index
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run("qgis:createspatialindex", {"INPUT": INPUT})
    return out


def qgis_vector_read_vector(processing, context, INPUT):
    """qgis function read vector layer into memeory
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = QgsVectorLayer(INPUT, "")
    return out


def Obtain_Attribute_Table(
    processing, context, input_layer, Input_Is_Feature_In_Mem=True
):
    """Read QGIS vector layer attribute table

    Parameters
    ----------
    processing                        : qgis object
    context                           : qgis object
    vec_layer                         : qgis object
        a vector layer readed or generated by QGIS
    Returns:
    -------
    Attri_Tbl                  : Dataframe
       The attribute table in the input vector layer
    """
    
    if Input_Is_Feature_In_Mem:
        
        output_layer = os.path.join(tempfile.gettempdir(),'basinmaker_layer_'+str(np.random.randint(1, 10000 + 1))+".shp")
        
        _writer = QgsVectorFileWriter.writeAsVectorFormat(
        input_layer,
        output_layer,
        'utf-8',
        driverName='ESRI Shapefile'
        )
        vec_layer = output_layer

    else:
        vec_layer = input_layer
    
    Attri_Tbl = Dbf_To_Dataframe(vec_layer)

    # N_new_features = vec_layer.featureCount()
    # field_names = []
    # for field in vec_layer.fields():
    #     field_names.append(field.name())
    # Attri_Tbl = pd.DataFrame(
    #     data=np.full((N_new_features, len(field_names)), np.nan), columns=field_names
    # )
    # 
    # src_features = vec_layer.getFeatures()
    # i_sf = 0
    # for sf in src_features:
    #     for field in vec_layer.fields():
    #         filednm = field.name()
    #         Attri_Tbl.loc[Attri_Tbl.index[i_sf], filednm] = sf[filednm]
    #     i_sf = i_sf + 1
    return Attri_Tbl


def Clean_Attribute_Name(
    Input, FieldName_List, Input_Is_Feature_In_Mem=False, Col_NM_Max="SubId"
):
    """Function clean feature attribute table, all colnmun not in FieldName_List
        will be removed
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """

    fieldnames = set(FieldName_List)
    if Input_Is_Feature_In_Mem:
        layer_cat = Input
    else:
        layer_cat = QgsVectorLayer(Input, "")

    field_ids = []
    for field in layer_cat.fields():
        if field.name() not in fieldnames:
            field_ids.append(layer_cat.dataProvider().fieldNameIndex(field.name()))
        if field.name() == Col_NM_Max:
            max_subbasin_id = layer_cat.maximumValue(
                layer_cat.dataProvider().fieldNameIndex(field.name())
            )

    layer_cat.dataProvider().deleteAttributes(field_ids)
    layer_cat.updateFields()
    layer_cat.commitChanges()

    if Input_Is_Feature_In_Mem:
        return layer_cat, max_subbasin_id
    else:
        del layer_cat
        return max_subbasin_id

def change_neg_value_to_null_in_attribute_table(processing,Path_input_Layer,Path_output_layer):
        
    col_Name = 'Obs_NM'
    formular = 'CASE \r\n  WHEN \"%s\" =-9999 THEN NULL\r\n  ELSE \"%s\"\r\nEND' %(col_Name,col_Name)
    memout1 = processing.run("qgis:fieldcalculator", {
                   'INPUT':Path_input_Layer,
                   'FIELD_NAME':col_Name,
                   'FIELD_TYPE':2,
                   'FIELD_LENGTH':0,
                   'FIELD_PRECISION':0,
                   'FORMULA':formular,
                   'OUTPUT':"memory:",
                   })["OUTPUT"] 

    col_Name = 'SRC_obs'
    formular = 'CASE \r\n  WHEN \"%s\" =-9999 THEN NULL\r\n  ELSE \"%s\"\r\nEND' %(col_Name,col_Name)
    memout1 = processing.run("qgis:fieldcalculator", {
                   'INPUT':memout1,
                   'FIELD_NAME':col_Name,
                   'FIELD_TYPE':2,
                   'FIELD_LENGTH':0,
                   'FIELD_PRECISION':0,
                   'FORMULA':formular,
                   'OUTPUT':Path_output_layer,
                   })    
    return 
########
