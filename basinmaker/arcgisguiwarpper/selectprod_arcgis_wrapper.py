import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from basinmaker.postprocessing.selectprodarcgis import Select_Routing_product_based_SubId_arcgis
from basinmaker.postprocessing.combinearcgis import combine_catchments_covered_by_the_same_lake_arcgis

OutputFolder =sys.argv[7]
Path_Catchment_Polygon=sys.argv[1]
Path_River_Polyline=sys.argv[2]
Path_Con_Lake_ply=sys.argv[3]
Path_NonCon_Lake_ply=sys.argv[4]
mostdownid=int(sys.argv[5])
mostupstreamid=int(sys.argv[6])

arcpy.AddMessage(mostupstreamid)
arcpy.AddMessage(Path_NonCon_Lake_ply)
arcpy.AddMessage(OutputFolder)

Select_Routing_product_based_SubId_arcgis(
    OutputFolder = OutputFolder,
    Path_Catchment_Polygon=Path_Catchment_Polygon,
    Path_River_Polyline=Path_River_Polyline,
    Path_Con_Lake_ply=Path_Con_Lake_ply,
    Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
    mostdownid=mostdownid,
    mostupstreamid=mostupstreamid,
)


combine_catchments_covered_by_the_same_lake_arcgis(
    OutputFolder = OutputFolder, 
    Path_final_rivply=os.path.join(OutputFolder,"catchment_without_merging_lakes.shp"), 
    Path_final_riv=os.path.join(OutputFolder,"river_without_merging_lakes.shp"),
)