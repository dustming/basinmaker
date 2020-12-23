from processing_functions_raster_array import *
from processing_functions_raster_grass import *
from processing_functions_raster_qgis import *
from processing_functions_vector_grass import *
from processing_functions_vector_qgis import *
from utilities import *
import os
from preprocessing.reprojectandclipvectorbyplyqgis import (
    reproject_clip_vectors_by_polygon,
)
from preprocessing.rasterizevectorsandloadtodbqgis import (
    rasterize_vectors_and_load_to_db,
)


def preprocessing_input_polygon(
    grassdb,
    grass_location,
    qgis_prefix_path,
    mask,
    path_polygon,
    ply_name,
    attribute_names=[],
    raster_names=[],
):

    reproject_clip_vectors_by_polygon(
        grassdb,
        grass_location,
        qgis_prefix_path,
        mask=os.path.join(grassdb, mask + ".shp"),
        path_polygon=path_polygon,
        ply_name=ply_name,
    )
    for i in range(0, len(attribute_names)):
        rasterize_vectors_and_load_to_db(
            grassdb,
            grass_location,
            qgis_prefix_path,
            mask,
            vector_path=os.path.join(grassdb, ply_name + ".shp"),
            attribue_name=attribute_names[i],
            raster_name=raster_names[i],
        )
