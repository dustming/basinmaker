import pandas as pd
import numpy as np
from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *


def export_files_to_output_folder(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    output_riv,
    output_cat,
    output_folder,
    obs_attributes,
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

    formular = "x(centroid(transform($geometry,'%s','%s')))" % (
        layer_cat.crs().authid(),
        "EPSG:4326",
    )

    layer_cat = qgis_vector_field_calculator(
        processing=processing,
        context=context,
        FORMULA=formular,
        FIELD_NAME="centroid_x",
        INPUT=layer_cat,
        OUTPUT="memory:",
    )["OUTPUT"]
        

    formular = "y(centroid(transform($geometry,'%s','%s')))" % (
        layer_cat.crs().authid(),
        "EPSG:4326",
    )

    layer_cat = qgis_vector_field_calculator(
        processing=processing,
        context=context,
        FORMULA=formular,
        FIELD_NAME="centroid_y",
        INPUT=layer_cat,
        OUTPUT="memory:",
    )["OUTPUT"]


    Selectfeatureattributes(
        processing,
        Input=layer_cat,
        Output=os.path.join(grassdb, output_cat + "withnegv.shp"),
        Attri_NM="SubId",
        Values=subinfo[subinfo["SubId"] > 0]["SubId"].values,
    )
    Selectfeatureattributes(
        processing,
        Input=os.path.join(grassdb, cat_riv_info + "_dis.shp"),
        Output=os.path.join(grassdb, output_riv + "withnegv.shp"),
        Attri_NM="SubId",
        Values=subinfo[subinfo["SubId"] > 0]["SubId"].values,
    )
    
    change_neg_value_to_null_in_attribute_table(
                 processing,
                 os.path.join(grassdb, output_cat + "withnegv.shp"),
                 os.path.join(output_folder, output_cat + ".shp")
    )
    change_neg_value_to_null_in_attribute_table(
                 processing,
                 os.path.join(grassdb, output_riv + "withnegv.shp"),
                 os.path.join(output_folder, output_riv + ".shp")
    )
    
    # Add_centroid_to_feature(
    #     os.path.join(output_folder, output_cat + ".shp"), "centroid_x", "centroid_y"
    # )

    subinfo = Dbf_To_Dataframe(os.path.join(output_folder, output_cat + ".shp"))

    if os.path.exists(os.path.join(grassdb, all_lakes + ".shp")):
        cl_lakeids = subinfo.loc[subinfo["Lake_Cat"] == 1]["HyLakeId"].values
        ncl_lakeids = subinfo.loc[subinfo["Lake_Cat"] == 2]["HyLakeId"].values

        if len(cl_lakeids) > 0:
            Selectfeatureattributes(
                processing,
                Input=os.path.join(grassdb, all_lakes + ".shp"),
                Output=os.path.join(output_folder, "sl_connected_lake.shp"),
                Attri_NM="Hylak_id",
                Values=cl_lakeids,
            )
        if len(ncl_lakeids) > 0:
            Selectfeatureattributes(
                processing,
                Input=os.path.join(grassdb, all_lakes + ".shp"),
                Output=os.path.join(output_folder, "sl_non_connected_lake.shp"),
                Attri_NM="Hylak_id",
                Values=ncl_lakeids,
            )

    if (
        os.path.exists(os.path.join(grassdb, snapped_obs_points + ".shp"))
        and len(subinfo.loc[subinfo["Has_POI"] > 0]["Has_POI"].values) > 0
    ):

        Selectfeatureattributes(
            processing,
            Input=os.path.join(grassdb, snapped_obs_points + ".shp"),
            Output=os.path.join(output_folder, "obs_gauges.shp"),
            Attri_NM=obs_attributes[1],
            Values=subinfo.loc[subinfo["Has_POI"] > 0]["Obs_NM"].values,
            Is_str = True,
        )

    Clean_Attribute_Name(
        os.path.join(output_folder, output_cat + ".shp"), COLUMN_NAMES_CONSTANT_CLEAN
    )
    Clean_Attribute_Name(
        os.path.join(output_folder, output_riv + ".shp"), COLUMN_NAMES_CONSTANT_CLEAN
    )

    return
