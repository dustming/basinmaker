import sys, os
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from postprocessing.selectlakearcgis import simplify_routing_structure_by_filter_lakes_arcgis
from postprocessing.combinearcgis import combine_catchments_covered_by_the_same_lake_arcgis

Routing_Product_Folder = sys.argv[1]
Thres_Area_Conn_Lakes=float(sys.argv[2])
Thres_Area_Non_Conn_Lakes=float(sys.argv[3])
Selected_Lake_List_in = sys.argv[4] 
OutputFolder =sys.argv[5]
gis_platform = 'arcgis'

Selected_Lake_List_in = Selected_Lake_List_in.split(";")
Selected_Lake_List = []
for i in range(0,len(Selected_Lake_List_in)):
    try:
        Selected_Lake_List.append(int(Selected_Lake_List_in[i]))
    except:
        Selected_Lake_List = []

simplify_routing_structure_by_filter_lakes_arcgis(
    Routing_Product_Folder = Routing_Product_Folder,
    Thres_Area_Conn_Lakes=Thres_Area_Conn_Lakes,
    Thres_Area_Non_Conn_Lakes=Thres_Area_Non_Conn_Lakes,
    Selected_Lake_List_in=Selected_Lake_List,
    OutputFolder=OutputFolder,
    qgis_prefix_path='#',
    gis_platform=gis_platform,
)
            
combine_catchments_covered_by_the_same_lake_arcgis(
    Routing_Product_Folder=OutputFolder,
)  
