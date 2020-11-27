import pytest
from ToolboxClass import LRRT
import os
import pandas as pd 
from simpledbf import Dbf5
import shutil 
from processing_functions_attribute_table import Evaluate_Two_Dataframes
from utilities import Dbf_To_Dataframe
        
def test_SelectLakes():
    """test function that will: 
    Function that used to simplify the routing product by user 
    provided lake area thresthold. 
    The input catchment polygons is the routing product before  
    merging for lakes. It is provided with the routing product.
    The result is the simplified catchment polygons. But 
    result from this fuction still not merging catchment 
    covering by the same lake. Thus, The result generated 
    from this tools need further processed by 
    Define_Final_Catchment, or can be further processed by 
    Customize_Routing_Topology 
    """      

    ###Floder where store the inputs for tests function
    Routing_Product_Folder = './testdata/Simplified_By_DA/'
    
    ###Folder where store the expected resuts            
    Expect_Result_Folder = os.path.join('./testdata','Simplify_By_Lake_Area')
    
    ###Folder where the output will be generated 
    Output_Folder = os.path.join('./testdata','testout10')
    
    
    ###The pathes for all inputs 
    Path_Con_Lake_ply    = os.path.join(Routing_Product_Folder,'Con_Lake_Ply.shp')      ### Connected lake polygons
    Path_NonCon_Lake_ply = os.path.join(Routing_Product_Folder,'Non_Con_Lake_Ply.shp')  ### None connected lake polygons
    Path_final_riv_ply   = os.path.join(Routing_Product_Folder,'finalriv_info_ply.shp') ### River polyline 
    Path_final_riv       = os.path.join(Routing_Product_Folder,'finalriv_info.shp')     ### Catchment polygons 
    
    ###Generate test resuts       
    Lake_Area_thresthold_Connected_Lake = 10
    Lake_Area_thresthold_NonConnected_Lake = 3  
    RTtool=LRRT()
    RTtool.SelectLakes(Path_final_riv_ply = Path_final_riv_ply,
                                Path_final_riv = Path_final_riv, 
                                Path_Con_Lake_ply=Path_Con_Lake_ply,
                                Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
                                Thres_Area_Conn_Lakes = Lake_Area_thresthold_Connected_Lake,
                                Thres_Area_Non_Conn_Lakes = Lake_Area_thresthold_NonConnected_Lake,
                                Selection_Method = 'ByArea',
                                OutputFolder = Output_Folder)
    
    """Evaluate result attribute table of result polygon,polylie, and two 
       lake polygons      
    """ 

    ### transfer expected siplified product into pandas dataframe
    Expect_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalriv_info_ply.shp')).sort_values(by=['SubId'])

    ### transfer resulted siplified product into pandas dataframe    
    Result_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalriv_info_ply.shp')).sort_values(by=['SubId'])

    
    assert Evaluate_Two_Dataframes(Result_Finalriv_info_ply,Expect_Finalriv_info_ply,Check_Col_NM = 'SubId')


    """Evaluate lake polygon files 
    Con_Lake_Ply is the lake polygon that connected by river network 
    Non_Con_Lake_Ply is the lake polygon that did not connected by 
    river network 
    """ 
    ### transfer expected siplified connected lake polygon  into pandas dataframe Expect_Con_Lake_Ply
    Expect_Con_Lake_Ply      = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'Con_Lake_Ply.shp')).sort_values(by=['Hylak_id'])
    ### transfer expected siplified non connected lake polygon  into pandas dataframe Expect_Non_Con_Lake_Ply
    Expect_Non_Con_Lake_Ply  = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'Non_Con_Lake_Ply.shp')).sort_values(by=['Hylak_id'])
    
    ### transfer resulted siplified connected lake polygon  into pandas dataframe Result_Con_Lake_Ply
    Result_Con_Lake_Ply      = Dbf_To_Dataframe(os.path.join(Output_Folder,'Con_Lake_Ply.shp')).sort_values(by=['Hylak_id'])
    ### transfer resulted siplified non connected lake polygon  into pandas dataframe Result_Non_Con_Lake_Ply
    Result_Non_Con_Lake_Ply  = Dbf_To_Dataframe(os.path.join(Output_Folder,'Non_Con_Lake_Ply.shp')).sort_values(by=['Hylak_id'])
    
    ### compare two pandas dataframe Expect_Con_Lake_Ply and Result_Con_Lake_Ply
    assert Evaluate_Two_Dataframes(Expect_Con_Lake_Ply,Result_Con_Lake_Ply,Check_Col_NM = 'Hylak_id')
    ### compare two pandas dataframe Expect_Non_Con_Lake_Ply and Result_Non_Con_Lake_Ply
    assert Evaluate_Two_Dataframes(Expect_Non_Con_Lake_Ply,Result_Non_Con_Lake_Ply,Check_Col_NM = 'Hylak_id')
    
    shutil.rmtree(Output_Folder) 
        