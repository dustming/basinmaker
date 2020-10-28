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
        
def test_GenerateHRUS():
    """test function that will:
    
    be used to overlay: subbasin polygon, lake polygon (optional)
    , Land use polygon (optional), soil type polygon(optional),
    vegetation polygon (optional), and two other user defined polygons
    (optional).
    Non-lake HRU polygons in a subbasin is defined by an unique 
    combination of all user provided datasets.
    A lake HRU polygon is defined the same as the provided lake polygon.
    All value of landuse and Veg polygon covered by lake will
    be changed to 1, indicating it is a covered by lake.
    All value of the soil polygon covered by the lake will be change to
    the soil id of the polygon covered by the lake with largest area.  
    """      
    ###Floder where store the inputs for tests function
    Routing_Product_Folder = './testdata/Simplified_By_DA/'
    HRU_Folder = os.path.join('./testdata','HRU')
    
    ###Folder where store the expected resuts            
    Expect_Result_Folder = os.path.join('./testdata','Simplified_By_DA')
    
    ###Folder where the output will be generated 
    Output_Folder = os.path.join('./testdata','testout10')
    
    ###Generate test resuts            
    RTtool=LRRT()
    RTtool.GenerateHRUS(OutputFolder = Output_Folder,
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
                                                                  
    """Evaluate total number of HRU, total HRU area and total river length
    N_HRU is the total number of HRUs in the HRU polygon 
    len_Riv is the total river length in the routing network 
    HRU_Area is the total HRU area in the routing network     
    """ 
    
    ### transfer expected product into pandas dataframe
    Expect_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_hru_info.shp'))
    ### calcuate expected total number of HRUS:Expect_N_HRU
    Expect_N_HRU = len(Expect_Finalriv_info_ply)
    ### calcuate expected total river length :Expect_len_Riv
    Expect_len_Riv = sum(Expect_Finalriv_info_ply['RivLength'])
    ### calcuate expected total HRU area :Expect_HRU_Area
    Expect_HRU_Area = sum(Expect_Finalriv_info_ply['HRU_Area'])
    
    ### transfer resulted product into pandas dataframe    
    Result_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_hru_info.shp'))
    ### calcuate resulted total number of HRUS:Result_N_HRU
    Result_N_HRU = len(Result_Finalriv_info_ply)
    ### calcuate resulted total river length :Result_len_Riv
    Result_len_Riv = sum(Result_Finalriv_info_ply['RivLength'])
    ### calcuate resulted total HRU area :Result_HRU_Area
    Result_HRU_Area = sum(Result_Finalriv_info_ply['HRU_Area'])
    
    ### compare Expect_N_HRU and Result_N_HRU
    assert Expect_N_HRU == Result_N_HRU
    ### compare Expect_len_Riv and Result_len_Riv
    assert Expect_len_Riv == pytest.approx(Result_len_Riv, 0.1)
    ### compare Expect_HRU_Area and Expect_HRU_Area
    assert Expect_HRU_Area == pytest.approx(Expect_HRU_Area, 0.1)
                                          
    shutil.rmtree(Output_Folder) 
    
                    