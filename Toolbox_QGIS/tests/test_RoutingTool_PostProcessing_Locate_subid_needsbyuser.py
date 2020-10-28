import pytest
from ToolboxClass import LRRT

    
    
def test_Locate_subid_needsbyuser():
    """test function that will:
    
    Function that used to obtain subbasin ID of certain gauge.
    or subbasin ID of the polygon that includes the given point
    shapefile. 
    
    """  
        

    ###Floder where store the inputs for tests function 
    Routing_Product = './testdata/Routing_product_V2/finalcat_info.shp'
    
    
    #######
    #test for option 1
    #######
    ### Expected generated result 
    Expect_SubIDs = [1024]

    ###Define path of input dataset    
    Gauge_NM        = ['02LE024']
    ###Generate test resuts                     
    RTtool=LRRT()
    subids = RTtool.Locate_subid_needsbyuser(Gauge_NMS = Gauge_NM,Path_products=Routing_Product)
    ### compare the length of result suinds and expected result Expect_SubIDs
    assert len(subids) == len(Expect_SubIDs)
    ### compare the content of result suinds and expected result Expect_SubIDs
    assert all([a == b for a, b in zip(subids, Expect_SubIDs)])

    #######
    #test for option 2
    #######
    ### Expected generated result     
    Expect_SubIDs = [519, 629,720,1024,999,884,883,1213,1103,1200,1276,1224,
                     1276, 1276, 1276, 1277, 1286, 1340, 1347, 1349]                 
    ###Define path of input dataset    
    Pointshpfile = './testdata/Routing_product_V2/obspoint_snap.shp'
    ###Generate test resuts                 
    subids = RTtool.Locate_subid_needsbyuser(Path_Points = Pointshpfile, Path_products=Routing_Product)
    ### compare the length of result suinds and expected result Expect_SubIDs
    assert len(subids) == len(Expect_SubIDs)
    ### compare the content of result suinds and expected result Expect_SubIDs
    assert all([a == b for a, b in zip(subids, Expect_SubIDs)])