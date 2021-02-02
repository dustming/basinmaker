import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from postprocessing.selectlakearcgis import simplify_routing_structure_by_filter_lakes_arcgis
from postprocessing.combinearcgis import combine_catchments_covered_by_the_same_lake_arcgis

Path_final_riv_ply=sys.argv[1]
Path_final_riv=sys.argv[2]
Path_Con_Lake_ply=sys.argv[3]
Path_NonCon_Lake_ply=sys.argv[4]
Thres_Area_Conn_Lakes=float(sys.argv[5])
Thres_Area_Non_Conn_Lakes=float(sys.argv[6])
OutputFolder =sys.argv[7]

simplify_routing_structure_by_filter_lakes_arcgis(
    Path_final_riv_ply=Path_final_riv_ply,
    Path_final_riv=Path_final_riv,
    Path_Con_Lake_ply=Path_Con_Lake_ply,
    Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
    Thres_Area_Conn_Lakes=Thres_Area_Conn_Lakes,
    Thres_Area_Non_Conn_Lakes=Thres_Area_Non_Conn_Lakes,
    Selection_Method="ByArea",
    Selected_Lake_List_in=[],
    OutputFolder=OutputFolder,
    gis_platform="arcgis",
)

combine_catchments_covered_by_the_same_lake_arcgis(
    OutputFolder = OutputFolder, 
    Path_final_rivply=os.path.join(OutputFolder,"catchment_without_merging_lakes.shp"), 
    Path_final_riv=os.path.join(OutputFolder,"river_without_merging_lakes.shp"),
)
