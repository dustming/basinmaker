def combine_catchments_covered_by_the_same_lake(
    OutputFolder="#",
    Path_final_rivply="#",
    Path_final_riv="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
):

    if gis_platform == "qgis":
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from postprocessing.combine import (
            combine_catchments_covered_by_the_same_lake_qgis,
        )

        combine_catchments_covered_by_the_same_lake_qgis(
            OutputFolder=OutputFolder,
            Path_final_rivply=Path_final_rivply,
            Path_final_riv=Path_final_riv,
            qgis_prefix_path=qgis_prefix_path,
        )


def simplify_routing_structure_by_filter_lakes(
    Path_final_riv_ply="#",
    Path_final_riv="#",
    Path_Con_Lake_ply="#",
    Path_NonCon_Lake_ply="#",
    Thres_Area_Conn_Lakes=-1,
    Thres_Area_Non_Conn_Lakes=-1,
    Selection_Method="ByArea",
    Selected_Lake_List_in=[],
    OutputFolder="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
):

    if gis_platform == "qgis":
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from postprocessing.selectlake import (
            simplify_routing_structure_by_filter_lakes_qgis,
        )

        simplify_routing_structure_by_filter_lakes_qgis(
            Path_final_riv_ply=Path_final_riv_ply,
            Path_final_riv=Path_final_riv,
            Path_Con_Lake_ply=Path_Con_Lake_ply,
            Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
            Thres_Area_Conn_Lakes=Thres_Area_Conn_Lakes,
            Thres_Area_Non_Conn_Lakes=Thres_Area_Non_Conn_Lakes,
            Selection_Method=Selection_Method,
            Selected_Lake_List_in=Selected_Lake_List_in,
            OutputFolder=OutputFolder,
            qgis_prefix_path=qgis_prefix_path,
            gis_platform=gis_platform,
        )


def simplify_routing_structure_by_drainage_area(
    Path_final_riv_ply="#",
    Path_final_riv="#",
    Path_Con_Lake_ply="#",
    Path_NonCon_Lake_ply="#",
    Area_Min=-1,
    OutputFolder="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
):

    if gis_platform == "qgis":
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from postprocessing.increaseda import (
            simplify_routing_structure_by_drainage_area_qgis,
        )

        simplify_routing_structure_by_drainage_area_qgis(
            Path_final_riv_ply=Path_final_riv_ply,
            Path_final_riv=Path_final_riv,
            Path_Con_Lake_ply=Path_Con_Lake_ply,
            Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
            Area_Min=Area_Min,
            OutputFolder=OutputFolder,
            qgis_prefix_path=qgis_prefix_path,
        )
