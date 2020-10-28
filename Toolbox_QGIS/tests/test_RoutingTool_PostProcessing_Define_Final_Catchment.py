import pytest
from ToolboxClass import LRRT
import os
import pandas as pd 
from simpledbf import Dbf5
import shutil 

def Dbf_To_Dataframe(file_path):
    """Transfer an input dbf file to dataframe
    
    Parameters
    ---------- 
    file_path   : string
    Full path to a shapefile 
    
    Returns:
    -------
    dataframe   : datafame 
    a pandas dataframe of attribute table of input shapefile    
    """
    tempinfo = Dbf5(file_path[:-3] + "dbf")
    dataframe = tempinfo.to_dataframe()
    return dataframe
        
def test_Define_Final_Catchment():
    """test function that will:
    Generate the final lake river routing structure by merging subbasin 
    polygons that are covered by the same lake.
    The input are the catchment polygons and river segements 
    before merging for lakes. The input files can be output of
    any of following functions:
    SelectLakes, Select_Routing_product_based_SubId, 
    Customize_Routing_Topology,RoutingNetworkTopologyUpdateToolset_riv
    The result is the final catchment polygon that ready to be used for 
    hydrological modeling     
    """    
    ###Floder where store the inputs for tests function
    Routing_Product_Folder = './testdata/Simplified_By_DA/'
    
    ###Folder where store the expected resuts            
    Expect_Result_Folder = os.path.join('./testdata','Simplified_By_DA')
    ###Folder where the output will be generated 
    Output_Folder = os.path.join('./testdata','testout10')
    
    
    ###The pathes for all inputs 
    Path_final_riv_ply   = os.path.join(Routing_Product_Folder,'finalriv_info_ply.shp') ### River polyline 
    Path_final_riv       = os.path.join(Routing_Product_Folder,'finalriv_info.shp')     ### Catchment polygons 

    ###Generate test resuts            
    RTtool=LRRT()
    RTtool.Define_Final_Catchment(Path_final_rivply = Path_final_riv_ply,
                                           Path_final_riv = Path_final_riv, 
                                           OutputFolder = Output_Folder)  
                                           
    """Evaluate total number of subbasin, total subbasin area, total river length
       and total lake area 
    N_Cat is the total number of subbasins in the routing network 
    len_Riv is the total river length in the routing network 
    Bas_Area is the total subbasin area in the routing network
    Lake_Area is the total lake area in the routing network     
    """ 
    ### transfer expected  product into pandas dataframe
    Expect_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_info.shp')).sort_values(by=['SubId'])
    ### calcuate expected total number of catchment:Expect_N_Cat
    Expect_N_Cat = len(Expect_Finalriv_info_ply)
    ### calcuate expected total river length :Expect_len_Riv
    Expect_len_Riv = sum(Expect_Finalriv_info_ply['RivLength'])
    ### calcuate expected total basin area :Expect_Bas_Area
    Expect_Bas_Area = sum(Expect_Finalriv_info_ply['BasArea'])
    ### calcuate expected total lake area :Expect_Lake_Area
    Expect_Lake_Area = sum(Expect_Finalcat_info['LakeArea'])  
      
    ### transfer resulted  product into pandas dataframe    
    Result_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_info.shp')).sort_values(by=['SubId'])
    ### calcuate resulted total number of catchment:Result_N_Cat
    Result_N_Cat = len(Result_Finalriv_info_ply)
    ### calcuate resulted total river length :Result_len_Riv
    Result_len_Riv = sum(Result_Finalriv_info_ply['RivLength'])
    ### calcuate resulted total basin area :Result_Bas_Area
    Result_Bas_Area = sum(Result_Finalriv_info_ply['BasArea'])
    ### calcuate resulted total lake area :Result_Lake_Area
    Result_Lake_Area = sum(Result_Finalcat_info['LakeArea'])

    ### compare Expect_N_Cat and Result_N_Cat
    assert Expect_N_Cat == Result_N_Cat
    ### compare Expect_len_Riv and Result_len_Riv
    assert Expect_len_Riv == pytest.approx(Result_len_Riv, 0.1)
    ### compare Expect_Bas_Area and Result_Bas_Area
    assert Expect_Bas_Area == pytest.approx(Result_Bas_Area, 0.1)
    ### compare Expect_Lake_Area and Result_Lake_Area
    assert Expect_Lake_Area == pytest.approx(Result_Lake_Area, 0.1)
    
    shutil.rmtree(Output_Folder) 

    
            