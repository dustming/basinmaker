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
        
def test_Define_Final_Catchment():
    
#################### Test for result in Simplified_By_DA
    ###The second version of routing product 
    Routing_Product_Folder = './testdata/Simplified_By_DA/'
    
    Expect_Result_Folder = os.path.join('./testdata','Simplified_By_DA')
    Output_Folder = os.path.join('./testdata','testout')
    
    
    ###product that need futher process 
    Path_final_riv_ply   = os.path.join(Routing_Product_Folder,'finalriv_info_ply.shp') ### River polyline 
    Path_final_riv       = os.path.join(Routing_Product_Folder,'finalriv_info.shp')     ### Catchment polygons 
     
    RTtool=LRRT()
    subids = RTtool.Define_Final_Catchment(Path_final_rivply = Path_final_riv_ply,
                                           Path_final_riv = Path_final_riv, 
                                           OutputFolder = Output_Folder)                                           
    ###Read  result in pd dataframe 
    Result_Finalcat_info         = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_info.shp'))
    Result_Finalcat_info_riv     = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_info_riv.shp'))

    ###Read expexted result in pd dataframe 
    Expect_Finalcat_info         = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_info.shp'))
    Expect_Finalcat_info_riv     = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_info_riv.shp'))
    
    Expect_N_Cat = len(Expect_Finalcat_info)
    Expect_len_Riv = sum(Expect_Finalcat_info['RivLength'])
    Expect_Bas_Area = sum(Expect_Finalcat_info['BasArea'])
    Expect_Lake_Area = sum(Expect_Finalcat_info['LakeArea'])
    Result_N_Cat = len(Result_Finalcat_info)
    Result_len_Riv = sum(Result_Finalcat_info['RivLength'])
    Result_Bas_Area = sum(Result_Finalcat_info['BasArea'])
    Result_Lake_Area = sum(Result_Finalcat_info['LakeArea'])
    assert Expect_N_Cat == Result_N_Cat
    assert Expect_len_Riv == pytest.approx(Result_len_Riv, 0.1)
    assert Expect_Bas_Area == pytest.approx(Result_Bas_Area, 0.1)
    assert Expect_Lake_Area == pytest.approx(Result_Lake_Area, 0.1)
    shutil.rmtree(Output_Folder) 

# #################### Test for result in Simplified_By_DA
#     ###The second version of routing product 
#     Routing_Product_Folder = './testdata/Simplify_By_Lake_Area/'
# 
#     Expect_Result_Folder = os.path.join('./testdata','Simplify_By_Lake_Area')
#     Output_Folder = os.path.join('./testdata','testout')
# 
# 
#     ###product that need futher process 
#     Path_final_riv_ply   = os.path.join(Routing_Product_Folder,'finalriv_info_ply.shp') ### River polyline 
#     Path_final_riv       = os.path.join(Routing_Product_Folder,'finalriv_info.shp')     ### Catchment polygons 
# 
#     RTtool=LRRT()
#     subids = RTtool.Define_Final_Catchment(Path_final_rivply = Path_final_riv_ply,
#                                            Path_final_riv = Path_final_riv, 
#                                            OutputFolder = Output_Folder)                                           
#     ###Read  result in pd dataframe 
#     Result_Finalcat_info         = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_info.shp'))
#     Result_Finalcat_info_riv     = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_info_riv.shp'))
# 
#     ###Read expexted result in pd dataframe 
#     Expect_Finalcat_info         = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_info.shp'))
#     Expect_Finalcat_info_riv     = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_info_riv.shp'))
# 
#     Expect_N_Cat = len(Expect_Finalcat_info)
#     Expect_len_Riv = sum(Expect_Finalcat_info['RivLength'])
#     Expect_Bas_Area = sum(Expect_Finalcat_info['BasArea'])
#     Expect_Lake_Area = sum(Expect_Finalcat_info['LakeArea'])
#     Result_N_Cat = len(Result_Finalcat_info)
#     Result_len_Riv = sum(Result_Finalcat_info['RivLength'])
#     Result_Bas_Area = sum(Result_Finalcat_info['BasArea'])
#     Result_Lake_Area = sum(Result_Finalcat_info['LakeArea'])
#     assert Expect_N_Cat == Result_N_Cat
#     assert Expect_len_Riv == pytest.approx(Result_len_Riv, 0.1)
#     assert Expect_Bas_Area == pytest.approx(Result_Bas_Area, 0.1)
#     assert Expect_Lake_Area == pytest.approx(Result_Lake_Area, 0.1)
#     shutil.rmtree(Output_Folder) 
    
            