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
        
def test_SelectLakes():
    

    ###The second version of routing product 
    Routing_Product_Folder = './testdata/Simplified_By_DA/'
    
    Expect_Result_Folder = os.path.join('./testdata','Simplify_By_Lake_Area')
    Output_Folder = os.path.join('./testdata','testout')
    
    
    ###The define pathes for inputs  
    Path_Con_Lake_ply    = os.path.join(Routing_Product_Folder,'Con_Lake_Ply.shp')      ### Connected lake polygons
    Path_NonCon_Lake_ply = os.path.join(Routing_Product_Folder,'Non_Con_Lake_Ply.shp')  ### None connected lake polygons
    ###product that need futher process 
    Path_final_riv_ply   = os.path.join(Routing_Product_Folder,'finalriv_info_ply.shp') ### River polyline 
    Path_final_riv       = os.path.join(Routing_Product_Folder,'finalriv_info.shp')     ### Catchment polygons 
    
    ## run extract function 
    Lake_Area_thresthold_Connected_Lake = 10
    Lake_Area_thresthold_NonConnected_Lake = 3  
    
    RTtool=LRRT()
    subids = RTtool.SelectLakes(Path_final_riv_ply = Path_final_riv_ply,
                                Path_final_riv = Path_final_riv, 
                                Path_Con_Lake_ply=Path_Con_Lake_ply,
                                Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
                                Thres_Area_Conn_Lakes = Lake_Area_thresthold_Connected_Lake,
                                Thres_Area_Non_Conn_Lakes = Lake_Area_thresthold_NonConnected_Lake,
                                Selection_Method = 'ByArea',
                                OutputFolder = Output_Folder)
                                                       
    ###Read  result in pd dataframe 
    Result_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalriv_info_ply.shp'))
    Result_Finalriv_info     = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalriv_info.shp'))
    Result_Con_Lake_Ply      = Dbf_To_Dataframe(os.path.join(Output_Folder,'Con_Lake_Ply.shp'))
    Result_Non_Con_Lake_Ply  = Dbf_To_Dataframe(os.path.join(Output_Folder,'Non_Con_Lake_Ply.shp'))

    ###Read expexted result in pd dataframe 
    Expect_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalriv_info_ply.shp'))
    Expect_Finalriv_info     = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalriv_info.shp'))
    Expect_Con_Lake_Ply      = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'Con_Lake_Ply.shp'))
    Expect_Non_Con_Lake_Ply  = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'Non_Con_Lake_Ply.shp'))
    
    Expect_N_Cat = len(Expect_Finalriv_info_ply)
    Expect_len_Riv = sum(Expect_Finalriv_info_ply['RivLength'])
    Expect_Bas_Area = sum(Expect_Finalriv_info_ply['BasArea'])
    Expect_Lake_Area = sum(Expect_Finalriv_info_ply['LakeArea'])

    Result_N_Cat = len(Result_Finalriv_info_ply)
    Result_len_Riv = sum(Result_Finalriv_info_ply['RivLength'])
    Result_Bas_Area = sum(Result_Finalriv_info_ply['BasArea'])
    Result_Lake_Area = sum(Result_Finalriv_info_ply['LakeArea'])

    
    
    assert Expect_N_Cat == Result_N_Cat
    assert Expect_len_Riv == pytest.approx(Result_len_Riv, 0.1)
    assert Expect_Bas_Area == pytest.approx(Result_Bas_Area, 0.1)
    assert Expect_Lake_Area == pytest.approx(Result_Lake_Area, 0.1)
#    assert Result_Finalriv_info.equals(Expect_Finalriv_info)

#    assert Result_Finalriv_info.equals(Expect_Finalriv_info)
    assert Result_Con_Lake_Ply.equals(Expect_Con_Lake_Ply)
    assert Result_Non_Con_Lake_Ply.equals(Expect_Non_Con_Lake_Ply)
    
    shutil.rmtree(Output_Folder) 

        