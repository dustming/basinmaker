import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from postprocessing.selectprodarcgis import Select_Routing_product_based_SubId_arcgis

OutputFolder =sys.argv[6]
Path_Catchment_Polygon=sys.argv[1]
Path_River_Polyline=sys.argv[2]
Path_Con_Lake_ply=sys.argv[3]
Path_NonCon_Lake_ply=sys.argv[4]
mostdownid=int(sys.argv[5])
mostupstreamid=int(sys.argv[6])

Select_Routing_product_based_SubId_arcgis(
    OutputFolder = OutputFolder,
    Path_Catchment_Polygon=Path_Catchment_Polygon,
    Path_River_Polyline=Path_River_Polyline,
    Path_Con_Lake_ply=Path_Con_Lake_ply,
    Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
    mostdownid=mostdownid,
    mostupstreamid=mostupstreamid,
)
