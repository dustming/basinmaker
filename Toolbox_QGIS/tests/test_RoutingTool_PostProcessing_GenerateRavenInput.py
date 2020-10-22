import pytest
from ToolboxClass import LRRT
import os
import shutil 
import filecmp

def test_GenerateRavenInput(HYDAT_Path):
    '''
    CA_HYDAT       (string):     Path and filename of previously downloaded 
                                 external database containing streamflow observations, 
                                 e.g. HYDAT for Canada ("Hydat.sqlite3"). 
                                 Read from command line
    '''    
#################### Test for result in Simplified_By_DA
    ###The second version of routing product 
    Routing_Product_Folder = './testdata/Simplified_By_DA/'
    
    Expect_Result_Folder = os.path.join('./testdata','Simplified_By_DA')
    Output_Folder = os.path.join('./testdata','testout')
    CA_HYDAT   = HYDAT_Path
    ###product that need futher process 
    Path_final_hru_info   = os.path.join(Routing_Product_Folder,'finalcat_hru_info.shp') ### River polyline 
     
    RTtool=LRRT()
    subids = RTtool.GenerateRavenInput(Path_final_hru_info = Path_final_hru_info,
                                       OutputFolder = Output_Folder,
                                       Model_Name   = 'Mytest_model',
                                       WriteObsrvt  = True,
                                       Startyear    = 2010,
                                       EndYear      = 2014,
                                       CA_HYDAT     = CA_HYDAT,
                                       WarmUp       = 1) ##require hydat databse  
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','Mytest_model.rvh'),os.path.join(Expect_Result_Folder,'RavenInput','Mytest_model.rvh'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','channel_properties.rvp'),os.path.join(Expect_Result_Folder,'RavenInput','channel_properties.rvp'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','Lakes.rvh'),os.path.join(Expect_Result_Folder,'RavenInput','Lakes.rvh'))                      
    
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','obs','02LE013_720.rvt'),os.path.join(Expect_Result_Folder,'RavenInput','obs','02LE013_720.rvt'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','obs','02LE015_519.rvt'),os.path.join(Expect_Result_Folder,'RavenInput','obs','02LE015_519.rvt'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','obs','02LE024_1024.rvt'),os.path.join(Expect_Result_Folder,'RavenInput','obs','02LE024_1024.rvt'))                      
    assert filecmp.cmp(os.path.join(Output_Folder,'RavenInput','obs','obsinfo.csv'),os.path.join(Expect_Result_Folder,'RavenInput','obs','obsinfo.csv'))                      

    shutil.rmtree(Output_Folder) 