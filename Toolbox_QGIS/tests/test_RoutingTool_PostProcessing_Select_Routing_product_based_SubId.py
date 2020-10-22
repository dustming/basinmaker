import pytest
from ToolboxClass import LRRT
import os
import pandas as pd 
from simpledbf import Dbf5
import shutil 

def Dbf_To_Dataframe(file_path):
    tempinfo = Dbf5(file_path[:-3] + "dbf")
    dataframe = tempinfo.to_dataframe()
    return dataframe
        
def test_Locate_subid_needsbyuser():
    

    ###The second version of routing product 
    Routing_Product_Folder = './testdata/Routing_product_V2/'
    
    Expect_Result_Folder = os.path.join('./testdata','02LE024')
    Output_Folder = os.path.join('./testdata','test_out')
    
    ###Read expexted result in pd dataframe 
    Expect_Finalcat_info     = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_info.shp'))
    Expect_Finalcat_riv_info = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_info_riv.shp'))
    
    Expect_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalriv_info_ply.shp'))
    Expect_Finalriv_info     = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalriv_info.shp'))
    
    Expect_Con_Lake_Ply      = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'Con_Lake_Ply.shp'))
    Expect_Non_Con_Lake_Ply  = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'Non_Con_Lake_Ply.shp'))
    
    ###The interested Gauge Name 
    Path_Con_Lake_ply    = os.path.join(Routing_Product_Folder,'Con_Lake_Ply.shp')      ### Connected lake polygons
    Path_NonCon_Lake_ply = os.path.join(Routing_Product_Folder,'Non_Con_Lake_Ply.shp')  ### None connected lake polygons
    ###product that need futher process 
    Path_final_riv_ply   = os.path.join(Routing_Product_Folder,'finalriv_info_ply.shp') ### River polyline 
    Path_final_riv       = os.path.join(Routing_Product_Folder,'finalriv_info.shp')     ### Catchment polygons 
    ###product that do not need need futher process 
    Path_final_cat_ply   = os.path.join(Routing_Product_Folder,'finalcat_info.shp')     ### catchment polygons 
    Path_final_cat_riv   = os.path.join(Routing_Product_Folder,'finalcat_info_riv.shp') ### CRiver polyline 
    
    RTtool=LRRT()
    subids = RTtool.Select_Routing_product_based_SubId(OutputFolder = Output_Folder, 
                                                       Path_final_cat = Path_final_cat_ply,
                                                       Path_final_cat_riv = Path_final_cat_riv,
                                                       Path_Con_Lake_ply = Path_Con_Lake_ply,
                                                       Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
                                                       mostdownid = 1024
                                                       )
                                                       
    subids = RTtool.Select_Routing_product_based_SubId(OutputFolder = Output_Folder, 
                                                       Path_final_riv = Path_final_riv,
                                                       Path_final_riv_ply = Path_final_riv_ply,
                                                       Path_Con_Lake_ply = Path_Con_Lake_ply,
                                                       Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
                                                       mostdownid = 1024
                                                       )

    ###Read  result in pd dataframe 
    Result_Finalcat_info     = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_info.shp'))
    Result_Finalcat_riv_info = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_info_riv.shp'))
    
    Result_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalriv_info_ply.shp'))
    Result_Finalriv_info     = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalriv_info.shp'))
    
    Result_Con_Lake_Ply      = Dbf_To_Dataframe(os.path.join(Output_Folder,'Con_Lake_Ply.shp'))
    Result_Non_Con_Lake_Ply  = Dbf_To_Dataframe(os.path.join(Output_Folder,'Non_Con_Lake_Ply.shp'))
                                    
    assert Result_Finalcat_info.equals(Expect_Finalcat_info)
    assert Result_Finalcat_riv_info.equals(Expect_Finalcat_riv_info)
    assert Result_Finalriv_info_ply.equals(Expect_Finalriv_info_ply)
    assert Result_Finalriv_info.equals(Expect_Finalriv_info)
    assert Result_Con_Lake_Ply.equals(Expect_Con_Lake_Ply)
    assert Result_Non_Con_Lake_Ply.equals(Expect_Non_Con_Lake_Ply)
    
    shutil.rmtree(Output_Folder) 

    
        