from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import os


def preprocessing_lake_polygon(
    mask,
    path_lakefile_in="#",
    lake_attributes=[],
    grassdb="#",
    grass_location="#",
    qgis_prefix_path="#",
    threshold_con_lake=999999999,
    threshold_non_con_lake=999999999,
    gis_platform="qgis",
    lake_name="alllake",
    lake_boundary_name="lake_boundary",
    lakes_lg_cl_thres="lakes_lg_cl_thres",
    lakes_lg_ncl_thres="lakes_lg_ncl_thres",
):

    if gis_platform == "qgis":
        assert (
            grassdb != "#"
        ), "grass database folder is needed, when gis_platform = qgis"
        assert (
            grass_location != "#"
        ), "grass location name is needed, when gis_platform = qgis"
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from basinmaker.preprocessing.preinputpolygonqgis import preprocessing_input_polygon
        from basinmaker.preprocessing.reprojectandclipvectorbyplyqgis import (
            obtain_polygon_boundary,
        )
        from basinmaker.preprocessing.rasterizevectorsandloadtodbqgis import (
            rasterize_vectors_and_load_to_db,
            obtain_lake_vectors_larger_than_threstholds,
        )

    if path_lakefile_in != "#":
        preprocessing_input_polygon(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            mask=mask,
            path_polygon=path_lakefile_in,
            ply_name=lake_name,
            attribute_names=[lake_attributes[0]],
            raster_names=[lake_name],
        )

        # obtain lake boundary lines
        obtain_polygon_boundary(
            grassdb=grassdb,
            qgis_prefix_path=qgis_prefix_path,
            ply_path=os.path.join(grassdb, lake_name + ".shp"),
            output=os.path.join(grassdb, lake_boundary_name + ".shp"),
        )
        # obtain_lake_vectors_larger_than_threstholds(
        #     grassdb=grassdb,
        #     qgis_prefix_path=qgis_prefix_path,
        #     all_lakes=lake_name,
        #     lake_attributes=lake_attributes,
        #     threshold_con_lake=threshold_con_lake,
        #     threshold_non_con_lake=threshold_non_con_lake,
        #     lakes_lg_cl_thres=lakes_lg_cl_thres,
        #     lakes_lg_ncl_thres=lakes_lg_ncl_thres,
        # )

        rasterize_vectors_and_load_to_db(
            grassdb,
            grass_location,
            qgis_prefix_path,
            mask,
            vector_path=os.path.join(grassdb, lake_boundary_name + ".shp"),
            attribue_name=lake_attributes[0],
            raster_name=lake_boundary_name,
        )

        # rasterize_vectors_and_load_to_db(
        #     grassdb,
        #     grass_location,
        #     qgis_prefix_path,
        #     mask,
        #     vector_path=os.path.join(grassdb, lakes_lg_cl_thres + ".shp"),
        #     attribue_name=lake_attributes[0],
        #     raster_name=lakes_lg_cl_thres,
        # )
        # rasterize_vectors_and_load_to_db(
        #     grassdb,
        #     grass_location,
        #     qgis_prefix_path,
        #     mask,
        #     vector_path=os.path.join(grassdb, lakes_lg_ncl_thres + ".shp"),
        #     attribue_name=lake_attributes[0],
        #     raster_name=lakes_lg_ncl_thres,
        # )
