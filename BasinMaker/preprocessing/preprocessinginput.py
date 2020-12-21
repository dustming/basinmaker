from processing_functions_vector_qgis import qgis_vector_polygon_stro_lines
import os


def preprocessing_inputs(
    mask,
    path_lakefile_in="#",
    lake_attributes=[],
    grassdb="#",
    grass_location="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
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
        from preprocessing.preinputpolygonqgis import preprocessing_input_polygon
        from preprocessing.reprojectandclipvectorbyplyqgis import (
            obtain_polygon_boundary,
        )
        from preprocessing.rasterizevectorsandloadtodbqgis import (
            rasterize_vectors_and_load_to_db,
        )

    if path_lakefile_in != "#":
        preprocessing_input_polygon(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            mask=mask,
            path_polygon=path_lakefile_in,
            ply_name="Lake",
            attribute_names=[lake_attributes[0]],
            raster_names=["alllake"],
        )

        # obtain lake boundary lines
        obtain_polygon_boundary(
            grassdb=grassdb,
            qgis_prefix_path=qgis_prefix_path,
            ply_path=os.path.join(grassdb, "Lake_clipped.shp"),
            output=os.path.join(grassdb, "Lake_clipped_boundary.shp"),
        )
        rasterize_vectors_and_load_to_db(
            grassdb,
            grass_location,
            qgis_prefix_path,
            mask,
            vector_path=os.path.join(grassdb, "Lake_clipped_boundary.shp"),
            attribue_name=lake_attributes[0],
            raster_name="Lake_Bound",
        )
