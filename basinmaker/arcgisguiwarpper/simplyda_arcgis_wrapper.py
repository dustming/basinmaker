import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from postprocessing.increasedaarcgis import simplify_routing_structure_by_drainage_area_arcgis
from postprocessing.combinearcgis import combine_catchments_covered_by_the_same_lake_arcgis

Path_final_riv_ply=sys.argv[1]
Path_final_riv=sys.argv[2]
Path_Con_Lake_ply=sys.argv[3]
Path_NonCon_Lake_ply=sys.argv[4]
Area_Min=int(sys.argv[5])
OutputFolder =sys.argv[6]

simplify_routing_structure_by_drainage_area_arcgis(
    Path_final_riv_ply=Path_final_riv_ply,
    Path_final_riv=Path_final_riv,
    Path_Con_Lake_ply=Path_Con_Lake_ply,
    Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
    Area_Min=Area_Min,
    OutputFolder=OutputFolder,
)

combine_catchments_covered_by_the_same_lake_arcgis(
    OutputFolder = OutputFolder, 
    Path_final_rivply=os.path.join(OutputFolder,"catchment_without_merging_lakes.shp"), 
    Path_final_riv=os.path.join(OutputFolder,"river_without_merging_lakes.shp"),
)