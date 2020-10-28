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
    Output_Folder = os.path.join('./testdata','testout')
    
    
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
    
    """Evaluate total number of subbasin, total subbasin area, total river length
       and total lake area. 
    N_Cat is the total number of subbasins in the simplified routing network 
    len_Riv is the total river length in the simplified routing network 
    Bas_Area is the total subbasin area in the simplified routing network
    Lake_Area is the total lake area in the simplified routing network     
    """ 

    ### transfer expected siplified product into pandas dataframe
    Expect_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'finalriv_info_ply.shp')).sort_values(by=['SubId'])
    ### calcuate expected total number of catchment:Expect_N_Cat
    Expect_N_Cat = len(Expect_Finalriv_info_ply)
    ### calcuate expected total river length :Expect_len_Riv
    Expect_len_Riv = sum(Expect_Finalriv_info_ply['RivLength'])
    ### calcuate expected total basin area :Expect_Bas_Area
    Expect_Bas_Area = sum(Expect_Finalriv_info_ply['BasArea'])
    ### calcuate expected total lake area: Expect_Lake_Area
    Expect_Lake_Area = sum(Expect_Finalriv_info_ply['LakeArea'])

    ### transfer resulted siplified product into pandas dataframe    
    Result_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Output_Folder,'finalriv_info_ply.shp')).sort_values(by=['SubId'])
    ### calcuate resulted total number of catchment:Result_N_Cat
    Result_N_Cat = len(Result_Finalriv_info_ply)
    ### calcuate resulted total river length :Result_len_Riv
    Result_len_Riv = sum(Result_Finalriv_info_ply['RivLength'])
    ### calcuate resulted total basin area :Result_Bas_Area
    Result_Bas_Area = sum(Result_Finalriv_info_ply['BasArea'])
    ### calcuate resulted total lake area: Result_Lake_Area
    Result_Lake_Area = sum(Result_Finalriv_info_ply['LakeArea'])
    
    ### compare Expect_N_Cat and Result_N_Cat
    assert Expect_N_Cat == Result_N_Cat
    ### compare Expect_len_Riv and Result_len_Riv
    assert Expect_len_Riv == pytest.approx(Result_len_Riv, 0.1)
    ### compare Expect_Bas_Area and Result_Bas_Area
    assert Expect_Bas_Area == pytest.approx(Result_Bas_Area, 0.1)
    ### compare Expect_Lake_Area and Result_Lake_Area    
    assert Expect_Lake_Area == pytest.approx(Result_Lake_Area, 0.1)

    """Evaluate lake polygon files 
    Con_Lake_Ply is the lake polygon that connected by river network 
    Non_Con_Lake_Ply is the lake polygon that did not connected by 
    river network 
    """ 
    ### transfer expected siplified connected lake polygon  into pandas dataframe Expect_Con_Lake_Ply
    Expect_Con_Lake_Ply      = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'Con_Lake_Ply.shp'))
    ### transfer expected siplified non connected lake polygon  into pandas dataframe Expect_Non_Con_Lake_Ply
    Expect_Non_Con_Lake_Ply  = Dbf_To_Dataframe(os.path.join(Expect_Result_Folder,'Non_Con_Lake_Ply.shp'))
    
    ### transfer resulted siplified connected lake polygon  into pandas dataframe Result_Con_Lake_Ply
    Result_Con_Lake_Ply      = Dbf_To_Dataframe(os.path.join(Output_Folder,'Con_Lake_Ply.shp'))
    ### transfer resulted siplified non connected lake polygon  into pandas dataframe Result_Non_Con_Lake_Ply
    Result_Non_Con_Lake_Ply  = Dbf_To_Dataframe(os.path.join(Output_Folder,'Non_Con_Lake_Ply.shp'))
    
    ### compare two pandas dataframe Expect_Con_Lake_Ply and Result_Con_Lake_Ply
    assert Result_Con_Lake_Ply.equals(Expect_Con_Lake_Ply)
    ### compare two pandas dataframe Expect_Non_Con_Lake_Ply and Result_Non_Con_Lake_Ply
    assert Result_Non_Con_Lake_Ply.equals(Expect_Non_Con_Lake_Ply)
    
    shutil.rmtree(Output_Folder) 
    
        