import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from postprocessing.selectprodarcgis import Select_Routing_product_based_SubId_arcgis
from postprocessing.combinearcgis import combine_catchments_covered_by_the_same_lake_arcgis

OutputFolder =sys.argv[4]
mostdownid_in=sys.argv[2]
mostupstreamid_in=sys.argv[3]
Routing_Product_Folder = sys.argv[1]

arcpy.AddMessage(Routing_Product_Folder)

mostdownid_in = mostdownid_in.split(";")
mostdownid = []
for i in range(0,len(mostdownid_in)):
    try:
        mostdownid.append(int(mostdownid_in[i]))
    except:
        mostdownid = []
        arcpy.AddMessage("Invalid input most down stream subbasin IDs  ", mostdownid_in)
        exist()

mostupstreamid_in = mostupstreamid_in.split(";")
        
mostupid = []
for i in range(0,len(mostupstreamid_in)):
    try:
        mostupid.append(int(mostupstreamid_in[i]))
    except:
        mostupid = []
        arcpy.AddMessage("Invalid input most up stream subbasin IDs  ", mostupid)
        exist()

arcpy.AddMessage(mostdownid)
arcpy.AddMessage(mostupid)
        
Select_Routing_product_based_SubId_arcgis(
    OutputFolder=OutputFolder,
    mostdownid=mostdownid,
    mostupstreamid=mostupid,
    Routing_Product_Folder = Routing_Product_Folder,
)
        

combine_catchments_covered_by_the_same_lake_arcgis(
    Routing_Product_Folder=OutputFolder,
)  