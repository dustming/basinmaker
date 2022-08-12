import os


def combine_catchments_covered_by_the_same_lake_method(
    Routing_Product_Folder = '#',
    qgis_prefix_path="#",
    area_thresthold = 10*30*30/1000/1000,
    gis_platform="qgis",
):

    if gis_platform == "qgis":
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from basinmaker.postprocessing.combine import (
            combine_catchments_covered_by_the_same_lake_qgis,
        )

        combine_catchments_covered_by_the_same_lake_qgis(
            # Path_final_rivply=Path_final_rivply,
            # Path_final_riv=Path_final_riv,
            Routing_Product_Folder = Routing_Product_Folder,
            qgis_prefix_path=qgis_prefix_path,
        )
    if gis_platform == "arcgis":
        from basinmaker.postprocessing.combinearcgis import (
            combine_catchments_covered_by_the_same_lake_arcgis,
        )

        combine_catchments_covered_by_the_same_lake_arcgis(
            Routing_Product_Folder=Routing_Product_Folder,
        )

    if gis_platform == "purepy":
        from basinmaker.postprocessing.combinepurepy import (
            combine_catchments_covered_by_the_same_lake_purepy,
        )

        combine_catchments_covered_by_the_same_lake_purepy(
            Routing_Product_Folder=Routing_Product_Folder,
            area_thresthold = area_thresthold,
        )


def simplify_routing_structure_by_filter_lakes_method(
    Path_final_riv_ply="#",
    Path_final_riv="#",
    Path_Con_Lake_ply="#",
    Path_NonCon_Lake_ply="#",
    Routing_Product_Folder = '#',
    Thres_Area_Conn_Lakes=-1,
    Thres_Area_Non_Conn_Lakes=-1,
    Selected_Lake_List_in=[],
    OutputFolder="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
):

    if gis_platform == "qgis":
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from basinmaker.postprocessing.selectlake import (
            simplify_routing_structure_by_filter_lakes_qgis,
        )

        simplify_routing_structure_by_filter_lakes_qgis(
            Routing_Product_Folder = Routing_Product_Folder,
            Thres_Area_Conn_Lakes=Thres_Area_Conn_Lakes,
            Thres_Area_Non_Conn_Lakes=Thres_Area_Non_Conn_Lakes,
            Selected_Lake_List_in=Selected_Lake_List_in,
            OutputFolder=OutputFolder,
            qgis_prefix_path=qgis_prefix_path,
            gis_platform=gis_platform,
        )
    if gis_platform == "arcgis":
        from basinmaker.postprocessing.selectlakearcgis import (
            simplify_routing_structure_by_filter_lakes_arcgis,
        )

        simplify_routing_structure_by_filter_lakes_arcgis(
            Routing_Product_Folder = Routing_Product_Folder,
            Thres_Area_Conn_Lakes=Thres_Area_Conn_Lakes,
            Thres_Area_Non_Conn_Lakes=Thres_Area_Non_Conn_Lakes,
            Selected_Lake_List_in=Selected_Lake_List_in,
            OutputFolder=OutputFolder,
            qgis_prefix_path=qgis_prefix_path,
            gis_platform=gis_platform,
        )

    if gis_platform == "purepy":
        from basinmaker.postprocessing.selectlakepurepy import (
            simplify_routing_structure_by_filter_lakes_purepy,
        )

        simplify_routing_structure_by_filter_lakes_purepy(
            Routing_Product_Folder = Routing_Product_Folder,
            Thres_Area_Conn_Lakes=Thres_Area_Conn_Lakes,
            Thres_Area_Non_Conn_Lakes=Thres_Area_Non_Conn_Lakes,
            Selected_Lake_List_in=Selected_Lake_List_in,
            OutputFolder=OutputFolder,
            qgis_prefix_path=qgis_prefix_path,
            gis_platform=gis_platform,
        )


def simplify_routing_structure_by_drainage_area_method(
    Routing_Product_Folder='#',
    Area_Min=-1,
    OutputFolder="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
):

    if gis_platform == "qgis":
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from basinmaker.postprocessing.increaseda import (
            simplify_routing_structure_by_drainage_area_qgis,
        )

        simplify_routing_structure_by_drainage_area_qgis(
            Routing_Product_Folder = Routing_Product_Folder,
            Area_Min=Area_Min,
            OutputFolder=OutputFolder,
            qgis_prefix_path=qgis_prefix_path,
        )

    if gis_platform == "arcgis":
        from basinmaker.postprocessing.increasedaarcgis import (
            simplify_routing_structure_by_drainage_area_arcgis,
        )

        simplify_routing_structure_by_drainage_area_arcgis(
            Routing_Product_Folder = Routing_Product_Folder,
            Area_Min=Area_Min,
            OutputFolder=OutputFolder,
        )

    if gis_platform == "purepy":
        from basinmaker.postprocessing.increasedapurepy import (
            simplify_routing_structure_by_drainage_area_purepy,
        )

        simplify_routing_structure_by_drainage_area_purepy(
            Routing_Product_Folder = Routing_Product_Folder,
            Area_Min=Area_Min,
            OutputFolder=OutputFolder,
        )

def select_part_of_routing_product_method(
    Path_Points,
    Gauge_NMS,
    OutputFolder,
    mostdownid,
    mostupid,
    Path_Catchment_Polygon="#",
    Path_River_Polyline="#",
    Path_Con_Lake_ply="#",
    Path_NonCon_Lake_ply="#",
    qgis_prefix_path="#",
    Routing_Product_Folder = '#',
    gis_platform="qgis",
):

    if gis_platform == "qgis":
        from basinmaker.postprocessing.selectprod import (
            Locate_subid_needsbyuser_qgis,
            Select_Routing_product_based_SubId_qgis,
        )

        if Gauge_NMS[0] != "#" or Path_Points != "#":
            subids = Locate_subid_needsbyuser_qgis(
                Path_Points=Path_Points,
                Gauge_NMS=Gauge_NMS,
                Path_products=Path_Catchment_Polygon,
                qgis_prefix_path=qgis_prefix_path,
            )
            if len(subids) <= 0:
                print("subbasin did not found        ")
                return

            subid = subids[0]
            Select_Routing_product_based_SubId_qgis(
                OutputFolder=OutputFolder,
                Path_Catchment_Polygon=Path_Catchment_Polygon,
                Path_River_Polyline=Path_River_Polyline,
                Path_Con_Lake_ply=Path_Con_Lake_ply,
                Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
                mostdownid=subid,
                mostupstreamid=-1,
                qgis_prefix_path=qgis_prefix_path,
            )

        else:
            Select_Routing_product_based_SubId_qgis(
                OutputFolder=OutputFolder,
                # Path_Catchment_Polygon=Path_Catchment_Polygon,
                # Path_River_Polyline=Path_River_Polyline,
                # Path_Con_Lake_ply=Path_Con_Lake_ply,
                # Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
                mostdownid=mostdownid,
                mostupstreamid=mostupid,
                Routing_Product_Folder = Routing_Product_Folder,
                qgis_prefix_path=qgis_prefix_path,
            )

    if gis_platform == "arcgis":
        from basinmaker.postprocessing.selectprodarcgis import (
        Select_Routing_product_based_SubId_arcgis
        )
        Select_Routing_product_based_SubId_arcgis(
            OutputFolder=OutputFolder,
            mostdownid=mostdownid,
            mostupstreamid=mostupid,
            Routing_Product_Folder = Routing_Product_Folder,
        )

    if gis_platform == "purepy":
        from basinmaker.postprocessing.selectprodpurepy import (
        Select_Routing_product_based_SubId_purepy
        )
        Select_Routing_product_based_SubId_purepy(
            OutputFolder=OutputFolder,
            mostdownid=mostdownid,
            mostupstreamid=mostupid,
            Routing_Product_Folder = Routing_Product_Folder,
        )

    return


def generate_hrus_method(
    Path_Subbasin_Ply,
    Landuse_info,
    Soil_info,
    Veg_info,
    Sub_Lake_ID="HyLakeId",
    Sub_ID="SubId",
    Path_Connect_Lake_ply="#",
    Path_Non_Connect_Lake_ply="#",
    Lake_Id="Hylak_id",
    Path_Landuse_Ply="#",
    Landuse_ID="Landuse_ID",
    Path_Soil_Ply="#",
    Soil_ID="Soil_ID",
    Path_Veg_Ply="#",
    Veg_ID="Veg_ID",
    Path_Other_Ply_1="#",
    Other_Ply_ID_1="O_ID_1",
    Path_Other_Ply_2="#",
    Other_Ply_ID_2="O_ID_2",
    DEM="#",
    Inmportance_order = ["Hylak_id","Landuse_ID"],
    min_hru_area_pct_sub = 0.2,
    Project_crs="EPSG:3573",
    OutputFolder="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
    pixel_size = 30,
    area_ratio_thresholds = [0,0,0]
):
    if gis_platform == "qgis":
        from basinmaker.postprocessing.hru import GenerateHRUS_qgis

        GenerateHRUS_qgis(
            Path_Subbasin_Ply=Path_Subbasin_Ply,
            Landuse_info=Landuse_info,
            Soil_info=Soil_info,
            Veg_info=Veg_info,
            Sub_Lake_ID=Sub_Lake_ID,
            Sub_ID=Sub_ID,
            Path_Connect_Lake_ply=Path_Connect_Lake_ply,
            Path_Non_Connect_Lake_ply=Path_Non_Connect_Lake_ply,
            Lake_Id=Lake_Id,
            Path_Landuse_Ply=Path_Landuse_Ply,
            Landuse_ID=Landuse_ID,
            Path_Soil_Ply=Path_Soil_Ply,
            Soil_ID=Soil_ID,
            Path_Veg_Ply=Path_Veg_Ply,
            Veg_ID=Veg_ID,
            Path_Other_Ply_1=Path_Other_Ply_1,
            Other_Ply_ID_1=Other_Ply_ID_1,
            Path_Other_Ply_2=Path_Other_Ply_2,
            Other_Ply_ID_2=Other_Ply_ID_2,
            DEM=DEM,
            Project_crs=Project_crs,
            OutputFolder=OutputFolder,
            qgis_prefix_path=qgis_prefix_path,
            importance_order = Inmportance_order,
            min_hru_pct_sub_area = min_hru_area_pct_sub,
        )

    if gis_platform == "arcgis":
        from basinmaker.postprocessing.hruarcgis import GenerateHRUS_arcgis

        GenerateHRUS_arcgis(
            Path_Subbasin_Ply=Path_Subbasin_Ply,
            Landuse_info=Landuse_info,
            Soil_info=Soil_info,
            Veg_info=Veg_info,
            Sub_Lake_ID=Sub_Lake_ID,
            Sub_ID=Sub_ID,
            Path_Connect_Lake_ply=Path_Connect_Lake_ply,
            Path_Non_Connect_Lake_ply=Path_Non_Connect_Lake_ply,
            Lake_Id=Lake_Id,
            Path_Landuse_Ply=Path_Landuse_Ply,
            Landuse_ID=Landuse_ID,
            Path_Soil_Ply=Path_Soil_Ply,
            Soil_ID=Soil_ID,
            Path_Veg_Ply=Path_Veg_Ply,
            Veg_ID=Veg_ID,
            Path_Other_Ply_1=Path_Other_Ply_1,
            Other_Ply_ID_1=Other_Ply_ID_1,
            Path_Other_Ply_2=Path_Other_Ply_2,
            Other_Ply_ID_2=Other_Ply_ID_2,
            Inmportance_order = Inmportance_order,
            min_hru_area_pct_sub = min_hru_area_pct_sub,
            DEM=DEM,
            Project_crs=Project_crs,
            OutputFolder=OutputFolder,
        )


    if gis_platform == "purepy":
        from basinmaker.postprocessing.hrupurepy import GenerateHRUS_purepy

        GenerateHRUS_purepy(
            Path_Subbasin_Ply=Path_Subbasin_Ply,
            Landuse_info=Landuse_info,
            Soil_info=Soil_info,
            Veg_info=Veg_info,
            Sub_Lake_ID=Sub_Lake_ID,
            Sub_ID=Sub_ID,
            Path_Connect_Lake_ply=Path_Connect_Lake_ply,
            Path_Non_Connect_Lake_ply=Path_Non_Connect_Lake_ply,
            Lake_Id=Lake_Id,
            Path_Landuse_Ply=Path_Landuse_Ply,
            Landuse_ID=Landuse_ID,
            Path_Soil_Ply=Path_Soil_Ply,
            Soil_ID=Soil_ID,
            Path_Veg_Ply=Path_Veg_Ply,
            Veg_ID=Veg_ID,
            Path_Other_Ply_1=Path_Other_Ply_1,
            Other_Ply_ID_1=Other_Ply_ID_1,
            Path_Other_Ply_2=Path_Other_Ply_2,
            Other_Ply_ID_2=Other_Ply_ID_2,
            Inmportance_order = Inmportance_order,
            min_hru_area_pct_sub = min_hru_area_pct_sub,
            DEM=DEM,
            Project_crs=Project_crs,
            OutputFolder=OutputFolder,
            pixel_size = pixel_size,
            area_ratio_thresholds = area_ratio_thresholds,
        )


    return


def obtain_grids_polygon_from_netcdf_file(
    qgis_prefix_path,
    netcdf_path="#",
    output_folder="#",
    coor_x_nm="lon",
    coor_y_nm="lat",
    is_rotated_grid=1,
    r_coor_x_nm="rlon",
    r_coor_y_nm="rlat",
    spatial_ref="EPSG:4326",
    x_add=-360,
    y_add=0,
    gis_platform = "qgis",
):



    if gis_platform == "qgis":
        from basinmaker.postprocessing.gridweight import (
            Generate_Grid_Poly_From_NetCDF_QGIS,
        )
        Generate_Grid_Poly_From_NetCDF_QGIS(
            NetCDF_Path=netcdf_path,
            Output_Folder=output_folder,
            Coor_x_NM=coor_x_nm,
            Coor_y_NM=coor_y_nm,
            Is_Rotated_Grid=is_rotated_grid,
            R_Coor_x_NM=r_coor_x_nm,
            R_Coor_y_NM=r_coor_y_nm,
            SpatialRef=spatial_ref,
            x_add=x_add,
            y_add=y_add,
            qgis_prefix_path=qgis_prefix_path,
        )



def generate_area_weight_of_two_polygons(
    target_polygon_path="#",
    mapping_polygon_path="#",
    col_nm="HRU_ID",
    output_folder="#",
    qgis_prefix_path = '#',
    gis_platform = "qgis"
):

    if gis_platform == "qgis":
        from basinmaker.postprocessing.gridweight import (
            Area_Weighted_Mapping_Between_Two_Polygons_QGIS,
        )
        Area_Weighted_Mapping_Between_Two_Polygons_QGIS(
            Target_Ply_Path=target_polygon_path,
            Mapping_Ply_Path=mapping_polygon_path,
            Col_NM=col_nm,
            Output_Folder=output_folder,
            qgis_prefix_path = qgis_prefix_path,
        )

    if gis_platform == "arcgis":
        from basinmaker.postprocessing.gridweightarcgis import (
            Area_Weighted_Mapping_Between_Two_Polygons_Arcgis,
        )
        Area_Weighted_Mapping_Between_Two_Polygons_Arcgis(
            Target_Ply_Path=target_polygon_path,
            Mapping_Ply_Path=mapping_polygon_path,
            Col_NM=col_nm,
            Output_Folder=output_folder,
        )
