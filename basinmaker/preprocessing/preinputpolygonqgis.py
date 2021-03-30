from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import os
from basinmaker.preprocessing.reprojectandclipvectorbyplyqgis import (
    reproject_clip_vectors_by_polygon,
)
from basinmaker.preprocessing.rasterizevectorsandloadtodbqgis import (
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
