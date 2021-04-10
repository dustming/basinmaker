import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from postprocessing.increasedaarcgis import simplify_routing_structure_by_drainage_area_arcgis
from postprocessing.combinearcgis import combine_catchments_covered_by_the_same_lake_arcgis

Routing_Product_Folder = sys.argv[1]
Area_Min=float(sys.argv[2])
OutputFolder =sys.argv[3]


simplify_routing_structure_by_drainage_area_arcgis(
    Routing_Product_Folder = Routing_Product_Folder,
    Area_Min=Area_Min,
    OutputFolder=OutputFolder,
)
        

combine_catchments_covered_by_the_same_lake_arcgis(
    Routing_Product_Folder=OutputFolder,
)  