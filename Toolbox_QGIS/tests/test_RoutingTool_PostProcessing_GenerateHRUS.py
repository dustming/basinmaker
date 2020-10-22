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
        
def test_GenerateHRUS():
    
#################### Test for result in Simplified_By_DA
    ###The second version of routing product 
    Routing_Product_Folder = './testdata/Simplified_By_DA/'
    HRU_Folder = os.path.join('./testdata','HRU')
    Expect_Result_Folder = os.path.join('./testdata','Simplified_By_DA')
    Output_Folder = os.path.join('./testdata','testout')
         
    RTtool=LRRT()
    subids = RTtool.GenerateHRUS(OutputFolder = Output_Folder,
                                 Path_Subbasin_Ply =  os.path.join(Routing_Product_Folder,'finalcat_info.shp'),
         						 Path_Connect_Lake_ply = os.path.join(Routing_Product_Folder,'Con_Lake_Ply.shp'),
								 Path_Non_Connect_Lake_ply = os.path.join(Routing_Product_Folder,'Non_Con_Lake_Ply.shp'),
                                 Path_Landuse_Ply = '#',Landuse_ID = 'Landuse_ID',
                                 Path_Soil_Ply = '#',Soil_ID = 'Soil_ID',
                                 Path_Veg_Ply = '#',Veg_ID = 'Veg_ID',
                                 Path_Other_Ply_1='#', Other_Ply_ID_1='O_ID_1',
                                 Path_Other_Ply_2='#', Other_Ply_ID_2='O_ID_2',
								 Landuse_info=os.path.join(HRU_Folder,'landuse_info.csv'),
								 Soil_info=os.path.join(HRU_Folder,'soil_info.csv'),
								 Veg_info=os.path.join(HRU_Folder,'veg_info.csv'),
								 DEM = '#')                                       
    ###Read  result in pd dataframe 
    Result_Finalhru_info         = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_hru_info.shp'))

    ###Read expexted result in pd dataframe 
    Expect_Finalhru_info         = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_hru_info.shp'))
    
    Expect_N_HRU = len(Expect_Finalhru_info)
    Expect_len_Riv = sum(Expect_Finalhru_info['RivLength'])
    Expect_HRU_Area = sum(Expect_Finalhru_info['HRU_Area'])

    Result_N_HRU = len(Result_Finalhru_info)
    Result_len_Riv = sum(Result_Finalhru_info['RivLength'])
    Result_HRU_Area = sum(Result_Finalhru_info['HRU_Area'])

    
    assert Expect_N_HRU == Result_N_HRU
    assert Expect_len_Riv == pytest.approx(Result_len_Riv, 0.1)
    assert Expect_HRU_Area == pytest.approx(Result_HRU_Area, 0.1)

    shutil.rmtree(Output_Folder) 
    
#################### Test for result in Simplified_By_DA
    ###The second version of routing product 
    Routing_Product_Folder = './testdata/Simplify_By_Lake_Area/'
    HRU_Folder = os.path.join('./testdata','HRU')
    Expect_Result_Folder = os.path.join('./testdata','Simplify_By_Lake_Area')
    Output_Folder = os.path.join('./testdata','testout')
         
    RTtool=LRRT()
    subids = RTtool.GenerateHRUS(OutputFolder = Output_Folder,
                                 Path_Subbasin_Ply =  os.path.join(Routing_Product_Folder,'finalcat_info.shp'),
         						 Path_Connect_Lake_ply = os.path.join(Routing_Product_Folder,'Con_Lake_Ply.shp'),
								 Path_Non_Connect_Lake_ply = os.path.join(Routing_Product_Folder,'Non_Con_Lake_Ply.shp'),
                                 Path_Landuse_Ply = '#',Landuse_ID = 'Landuse_ID',
                                 Path_Soil_Ply = '#',Soil_ID = 'Soil_ID',
                                 Path_Veg_Ply = '#',Veg_ID = 'Veg_ID',
                                 Path_Other_Ply_1='#', Other_Ply_ID_1='O_ID_1',
                                 Path_Other_Ply_2='#', Other_Ply_ID_2='O_ID_2',
								 Landuse_info=os.path.join(HRU_Folder,'landuse_info.csv'),
								 Soil_info=os.path.join(HRU_Folder,'soil_info.csv'),
								 Veg_info=os.path.join(HRU_Folder,'veg_info.csv'),
								 DEM = '#')                                       
    ###Read  result in pd dataframe 
    Result_Finalhru_info         = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_hru_info.shp'))

    ###Read expexted result in pd dataframe 
    Expect_Finalhru_info         = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_hru_info.shp'))
    
    Expect_N_HRU = len(Expect_Finalhru_info)
    Expect_len_Riv = sum(Expect_Finalhru_info['RivLength'])
    Expect_HRU_Area = sum(Expect_Finalhru_info['HRU_Area'])

    Result_N_HRU = len(Result_Finalhru_info)
    Result_len_Riv = sum(Result_Finalhru_info['RivLength'])
    Result_HRU_Area = sum(Result_Finalhru_info['HRU_Area'])

    
    assert Expect_N_HRU == Result_N_HRU
    assert Expect_len_Riv == pytest.approx(Result_len_Riv, 0.1)
    assert Expect_HRU_Area == pytest.approx(Result_HRU_Area, 0.1)

    shutil.rmtree(Output_Folder) 
                