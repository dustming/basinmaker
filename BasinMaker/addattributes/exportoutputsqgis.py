from utilities import *
import pandas as pd
import numpy as np
from processing_functions_vector_qgis import *


def export_files_to_output_folder(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_riv,
    input_cat,
    output_riv,
    output_cat,
    input_lake_path,
    selected_basin_ids,
):

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
            "INPUT": os.path.join(grassdb, input_riv + ".shp"),
            "FIELD": ["SubId"],
            "OUTPUT": os.path.join(grassdb, input_riv + "_dis.shp"),
        },
        context=context,
    )
    processing.run(
        "native:dissolve",
        {
            "INPUT": os.path.join(grassdb, input_cat + ".shp"),
            "FIELD": ["SubId"],
            "OUTPUT": os.path.join(grassdb, input_cat + "_dis.shp"),
        },
        context=context,
    )

    subinfo = Dbf_To_Dataframe(os.path.join(grassdb, input_cat + "_dis.shp"))

    if input_lake_path != "#":
        cl_lakeids = subinfo.loc[subinfo["IsLake"] == 1]["HyLakeId"].values
        ncl_lakeids = subinfo.loc[subinfo["IsLake"] == 2]["HyLakeId"].values

        Selectfeatureattributes(
            processing,
            Input=input_lake_path,
            Output=os.path.join(grassdb, "sl_connected_lake.shp"),
            Attri_NM="Hylak_id",
            Values=cl_lakeids,
        )
        Selectfeatureattributes(
            processing,
            Input=input_lake_path,
            Output=os.path.join(grassdb, "sl_non_connected_lake.shp"),
            Attri_NM="Hylak_id",
            Values=ncl_lakeids,
        )

    return
