import pytest
from ToolboxClass import LRRT

    
    
def test_Locate_subid_needsbyuser():
    

    ###The second version of routing product 
    Routing_Product = './testdata/Routing_product_V2/finalcat_info.shp'
    ###The interested Gauge Name 
    Gauge_NM        = ['02LE024']
    
    Expect_SubIDs = [1024]
    
    RTtool=LRRT()
    subids = RTtool.Locate_subid_needsbyuser(Guage_NMS = Gauge_NM,Path_products=Routing_Product)
    assert len(subids) == len(Expect_SubIDs)
    assert all([a == b for a, b in zip(subids, Expect_SubIDs)])
    
    #### for using point shpfile 
    Pointshpfile = './testdata/Routing_product_V2/obspoint_snap.shp'
    Expect_SubIDs = [519, 629,720,1024,999,884,883,1213,1103,1200,1276,1224,
                     1276, 1276, 1276, 1277, 1286, 1340, 1347, 1349]
    subids = RTtool.Locate_subid_needsbyuser(Path_Points = Pointshpfile, Path_products=Routing_Product)
    assert len(subids) == len(Expect_SubIDs)
    assert all([a == b for a, b in zip(subids, Expect_SubIDs)])