from utilities.utilities import *
from func.pdtable import *


def add_attributes_to_catchments(
    input_geo_names,
    path_bkfwidthdepth="#",
    bkfwd_attributes=[],
    path_landuse="#",
    path_landuse_info="#",
    out_cat_name="catchment_without_merging_lakes",
    out_riv_name="river_without_merging_lakes",
    grassdb="#",
    grass_location="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
    projection="EPSG:3573",
    obs_attributes=[],
    lake_attributes=[],
    outlet_obs_id=1,
    path_sub_reg_outlets_v="#",
    output_folder="#",
):

    columns = COLUMN_NAMES_CONSTANT
    coltypes = COLUMN_TYPES_CONSTANT

    # local geo file names
    cat_ply_info = Internal_Constant_Names["cat_ply_info"]
    cat_riv_info = Internal_Constant_Names["cat_riv_info"]
    outlet_pt_info = Internal_Constant_Names["outlet_pt_info"]

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
        from addattributes.createattributestemplateqgis import (
            create_catchments_attributes_template_table,
        )
        from addattributes.calculatebasicattributesqgis import (
            calculate_basic_attributes,
        )
        from addattributes.addlakeattributesqgis import add_lake_attributes
        from addattributes.joinpandastoattributesqgis import (
            join_pandas_table_to_vector_attributes,
        )
        from addattributes.exportoutputsqgis import export_files_to_output_folder
        from addattributes.addgaugeattributesqgis import add_gauge_attributes
        from addattributes.calfloodmanningnqgis import calculate_flood_plain_manning_n
        from addattributes.calbkfwidthdepthqgis import (
            calculate_bankfull_width_depth_from_polyline,
        )

        attr_template = create_catchments_attributes_template_table(
            grassdb=grassdb,
            grass_location=grass_location,
            columns=columns,
            input_geo_names=input_geo_names,
        )

        attr_basic = calculate_basic_attributes(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            input_geo_names=input_geo_names,
            projection=projection,
            catinfo=attr_template,
            cat_ply_info=cat_ply_info,
            cat_riv_info=cat_riv_info,
            outlet_pt_info=outlet_pt_info,
        )
        input_geo_names["cat_ply_info"] = cat_ply_info
        input_geo_names["cat_riv_info"] = cat_riv_info
        input_geo_names["outlet_pt_info"] = outlet_pt_info

        if len(lake_attributes) > 0:
            attr_lake = add_lake_attributes(
                grassdb=grassdb,
                grass_location=grass_location,
                qgis_prefix_path=qgis_prefix_path,
                input_geo_names=input_geo_names,
                catinfo=attr_basic,
            )
        else:
            attr_lake = attr_basic

        if len(obs_attributes) > 0:
            attr_obs = add_gauge_attributes(
                grassdb=grassdb,
                grass_location=grass_location,
                qgis_prefix_path=qgis_prefix_path,
                input_geo_names=input_geo_names,
                obs_attributes=obs_attributes,
                catinfo=attr_lake,
            )
        else:
            attr_obs = attr_lake

        if outlet_obs_id > 0:
            attr_select = return_interest_catchments_info(
                catinfo=attr_obs,
                outlet_obs_id=outlet_obs_id,
                path_sub_reg_outlets_v=path_sub_reg_outlets_v,
            )
        else:
            attr_select = attr_obs

        if path_landuse != "#":
            attr_landuse = calculate_flood_plain_manning_n(
                grassdb=grassdb,
                grass_location=grass_location,
                qgis_prefix_path=qgis_prefix_path,
                catinfo=attr_select,
                input_geo_names=input_geo_names,
                path_landuse=path_landuse,
                path_landuse_info=path_landuse_info,
            )
        else:
            attr_landuse = attr_select

        attr_da = streamorderanddrainagearea(attr_landuse)

        if path_bkfwidthdepth != "#":
            attr_bkf = calculate_bankfull_width_depth_from_polyline(
                grassdb=grassdb,
                grass_location=grass_location,
                qgis_prefix_path=qgis_prefix_path,
                path_bkfwidthdepth=path_bkfwidthdepth,
                bkfwd_attributes=bkfwd_attributes,
                catinfo=attr_da,
                input_geo_names=input_geo_names,
            )
        else:
            attr_bkf = attr_da

        attr_ncl = update_non_connected_catchment_info(attr_bkf)

        join_pandas_table_to_vector_attributes(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            vector_name=cat_ply_info,
            pd_table=attr_ncl,
            column_types=coltypes,
            columns_names=columns,
        )

        join_pandas_table_to_vector_attributes(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            vector_name=cat_riv_info,
            pd_table=attr_ncl,
            column_types=coltypes,
            columns_names=columns,
        )

        export_files_to_output_folder(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            input_geo_names=input_geo_names,
            output_riv=out_riv_name,
            output_cat=out_cat_name,
            output_folder=output_folder,
        )
