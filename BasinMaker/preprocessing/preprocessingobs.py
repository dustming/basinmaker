import os


def preprocessing_obs_point(
    mask,
    path_obsin_in="#",
    obs_attributes=[],
    grassdb="#",
    grass_location="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
    obsname="obs",
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
        from processing_functions_vector_grass import grass_raster_v_import

        preprocessing_input_polygon(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            mask=mask,
            path_polygon=path_obsin_in,
            ply_name=obsname,
            attribute_names=[obs_attributes[0]],
            raster_names=[obsname],
        )
