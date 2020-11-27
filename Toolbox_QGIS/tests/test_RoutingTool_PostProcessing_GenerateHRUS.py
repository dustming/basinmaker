import pytest
from ToolboxClass import LRRT
import os
import pandas as pd 
from simpledbf import Dbf5
import shutil 
from processing_functions_attribute_table import Evaluate_Two_Dataframes
from utilities import Dbf_To_Dataframe

        
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
    Expect_finalcat_hru_info = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalcat_hru_info.shp')).sort_values(by=['SubId','HRU_Area'])

    ### transfer resulted product into pandas dataframe    
    Result_finalcat_hru_info = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalcat_hru_info.shp')).sort_values(by=['SubId','HRU_Area'])
    
    print(Evaluate_Two_Dataframes(Result_finalcat_hru_info,Expect_finalcat_hru_info,Check_Col_NM = 'HRU_ID'))
    assert Evaluate_Two_Dataframes(Result_finalcat_hru_info,Expect_finalcat_hru_info,Check_Col_NM = 'HRU_ID')
    
                                          
    shutil.rmtree(Output_Folder)    
                    