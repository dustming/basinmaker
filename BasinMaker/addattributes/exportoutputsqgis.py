import pandas as pd
import numpy as np
from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *


def export_files_to_output_folder(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    output_riv,
    output_cat,
    output_folder,
):
    cat_riv_info = input_geo_names["cat_riv_info"]
    cat_ply_info = input_geo_names["cat_ply_info"]
    snapped_obs_points = input_geo_names["snapped_obs_points"]
    all_lakes = input_geo_names["all_lakes"]
    
    QgsApplication.setPrefixPath(qgis_prefix_path, True)
    Qgs = QgsApplication([], False)
    Qgs.initQgis()
    from processing.core.Processing import Processing
    from processing.tools import dataobjects
    from qgis import processing

    feedback = QgsProcessingFeedback()
    Processing.initialize()
    QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
    context = dataobjects.createContext()
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    processing.run(
        "native:dissolve",
        {
            "INPUT": os.path.join(grassdb, cat_riv_info + ".shp"),
            "FIELD": ["SubId"],
            "OUTPUT": os.path.join(grassdb, cat_riv_info + "_dis.shp"),
        },
        context=context,
    )
    processing.run(
        "native:dissolve",
        {
            "INPUT": os.path.join(grassdb, cat_ply_info + ".shp"),
            "FIELD": ["SubId"],
            "OUTPUT": os.path.join(grassdb, cat_ply_info + "_dis.shp"),
        },
        context=context,
    )

    subinfo = Dbf_To_Dataframe(os.path.join(grassdb, cat_ply_info + "_dis.shp"))

    layer_cat = QgsVectorLayer(os.path.join(grassdb, cat_ply_info + "_dis.shp"), "")
    # add attribute to layer
    layer_cat = qgis_vector_add_attributes(
        processing,
        context,
        INPUT_Layer=layer_cat,
        attribute_list=[
            QgsField("centroid_x", QVariant.Double),
            QgsField("centroid_y", QVariant.Double),
        ],
    )

    Selectfeatureattributes(
        processing,
        Input=layer_cat,
        Output=os.path.join(output_folder, output_cat + ".shp"),
        Attri_NM="SubId",
        Values=subinfo[subinfo["SubId"] > 0]["SubId"].values,
    )
    Selectfeatureattributes(
        processing,
        Input=os.path.join(grassdb, cat_riv_info + "_dis.shp"),
        Output=os.path.join(output_folder, output_riv + ".shp"),
        Attri_NM="SubId",
        Values=subinfo[subinfo["SubId"] > 0]["SubId"].values,
    )

    Add_centroid_to_feature(
        os.path.join(output_folder, output_cat + ".shp"), "centroid_x", "centroid_y"
    )

    subinfo = Dbf_To_Dataframe(os.path.join(output_folder, output_cat + ".shp"))
        
    if os.path.exists(os.path.join(grassdb,all_lakes+'.shp')):
        cl_lakeids = subinfo.loc[subinfo["IsLake"] == 1]["HyLakeId"].values
        ncl_lakeids = subinfo.loc[subinfo["IsLake"] == 2]["HyLakeId"].values

        if len(cl_lakeids) > 0:
            Selectfeatureattributes(
                processing,
                Input=os.path.join(grassdb,all_lakes+'.shp'),
                Output=os.path.join(output_folder, "sl_connected_lake.shp"),
                Attri_NM="Hylak_id",
                Values=cl_lakeids,
            )
        if len(ncl_lakeids) > 0:
            Selectfeatureattributes(
                processing,
                Input=os.path.join(grassdb,all_lakes+'.shp'),
                Output=os.path.join(output_folder, "sl_non_connected_lake.shp"),
                Attri_NM="Hylak_id",
                Values=ncl_lakeids,
            )
            
    if os.path.exists(os.path.join(grassdb,snapped_obs_points+'.shp')) and len(subinfo.loc[subinfo["IsObs"] > 0]["IsObs"].values) > 0:

        Selectfeatureattributes(
            processing,
            Input=os.path.join(grassdb, snapped_obs_points + ".shp"),
            Output=os.path.join(output_folder, "obs_gauges.shp"),
            Attri_NM="Obs_ID",
            Values=subinfo.loc[subinfo["IsObs"] > 0]["IsObs"].values,
        )

    Clean_Attribute_Name(
        os.path.join(output_folder, output_cat + ".shp"), COLUMN_NAMES_CONSTANT
    )
    Clean_Attribute_Name(
        os.path.join(output_folder, output_riv + ".shp"), COLUMN_NAMES_CONSTANT
    )

    return
