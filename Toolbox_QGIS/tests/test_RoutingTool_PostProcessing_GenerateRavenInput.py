import pytest
from ToolboxClass import LRRT
import os
import shutil 
import filecmp

def test_GenerateRavenInput(HYDAT_Path):
    """test function that will:
    
    Function that used to generate Raven input files. All output 
    will be stored in folder "<OutputFolder>/RavenInput". 
    
    Parameters
    ----------    
    CA_HYDAT       (string):     Path and filename of previously downloaded 
                                 external database containing streamflow observations, 
                                 e.g. HYDAT for Canada ("Hydat.sqlite3"). 
                                 Read from command line
    """    

    ###Floder where store the inputs for tests function 
    Routing_Product_Folder = './testdata/Simplified_By_DA/'
    
    ###Folder where store the expected resuts                
    Expect_Result_Folder = os.path.join('./testdata','Simplified_By_DA')

    ###Folder where the output will be generated 
    Output_Folder = os.path.join('./testdata','testout')
    
    ###Define path of input dataset    
    CA_HYDAT   = HYDAT_Path
    Path_final_hru_info   = os.path.join(Routing_Product_Folder,'finalcat_hru_info.shp') ### River polyline 

    ###Generate test resuts                 
    RTtool=LRRT()
    RTtool.GenerateRavenInput(Path_final_hru_info = Path_final_hru_info,
                                       OutputFolder = Output_Folder,
                                       Model_Name   = 'Mytest_model',
                                       WriteObsrvt  = True,
                                       Startyear    = 2010,
                                       EndYear      = 2014,
                                       CA_HYDAT     = CA_HYDAT,
                                       WarmUp       = 1) 
                                       
    """Evaluate generated raven input files in result and expected folder
      
    Mytest_model.rvh contains subbasins and HRUs
    channel_properties.rvp  contains definition and parameters for channels
    Lakes.rvh contains definition and parameters of lakes
    02LE013_720.rvt contains downloaded streamflow observation data
    02LE015_519.rvt contains downloaded streamflow observation data
    02LE024_1024.rvt contains downloaded streamflow observation data
    obsinfo.csv contains information of each gauge such as number of missing 
        value and drainage area. 
    """    
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','Mytest_model.rvh'),os.path.join(Expect_Result_Folder,'RavenInput','Mytest_model.rvh'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','channel_properties.rvp'),os.path.join(Expect_Result_Folder,'RavenInput','channel_properties.rvp'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','Lakes.rvh'),os.path.join(Expect_Result_Folder,'RavenInput','Lakes.rvh'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','obs','02LE013_720.rvt'),os.path.join(Expect_Result_Folder,'RavenInput','obs','02LE013_720.rvt'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','obs','02LE015_519.rvt'),os.path.join(Expect_Result_Folder,'RavenInput','obs','02LE015_519.rvt'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','obs','02LE024_1024.rvt'),os.path.join(Expect_Result_Folder,'RavenInput','obs','02LE024_1024.rvt'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','obs','obsinfo.csv'),os.path.join(Expect_Result_Folder,'RavenInput','obs','obsinfo.csv'))                      

    RTtool.Output_Clean()
