import pytest
from ToolboxClass import LRRT
import os


#
#BasinMaker_Folder = "C:/Users/dustm/Documents/GitHub/RoutingTool"
        
    ###Floder where store the inputs for tests function
Routing_Product_Folder = '../../tests/testdata/Routing_product_V1/'
HRU_Folder = os.path.join('../../tests/testdata','HRU')
    

###Folder where the output will be generated 
Output_Folder = Routing_Product_Folder #os.path.join('../../tests/testdata','test11')
    
RTtool=LRRT()


# Generate HRU 
# Connected lake polygon comes form hydrolake, needs to be extracted first 
# it is not provied in routing product 
#
RTtool.GenerateHRUS(OutputFolder = Output_Folder,
                    Path_Subbasin_Ply =  os.path.join(Routing_Product_Folder,'finalcat_info.shp'),
     				Path_Connect_Lake_ply = os.path.join(Routing_Product_Folder,'Con_Lake_Ply.shp'),
                    Path_Non_Connect_Lake_ply = '#',
                    Path_Landuse_Ply = '#',Landuse_ID = 'Landuse_ID',
                    Path_Soil_Ply = '#',Soil_ID = 'Soil_ID',
                    Path_Veg_Ply = '#',Veg_ID = 'Veg_ID',
                    Path_Other_Ply_1='#', Other_Ply_ID_1='O_ID_1',
                    Path_Other_Ply_2='#', Other_Ply_ID_2='O_ID_2',
                    Landuse_info=os.path.join(HRU_Folder,'landuse_info.csv'),
                    Soil_info=os.path.join(HRU_Folder,'soil_info.csv'),
                    Veg_info=os.path.join(HRU_Folder,'veg_info.csv'),
                    DEM = os.path.join(Routing_Product_Folder,'Liever_DEM.tif')) 
# define raven input, using output of generathru,'os.path.join(Output_Folder,'finalcat_hru_info.shp')'                     
RTtool.GenerateRavenInput(Path_final_hru_info = os.path.join(Output_Folder,'finalcat_hru_info.shp'), 
                          OutputFolder = Output_Folder,
                          Model_Name   = 'Mytest_model',
                          WriteObsrvt  = False,
#                          Startyear    = 2010, 
#                          EndYear      = 2014,
#                          CA_HYDAT     = CA_HYDAT,
                          Old_Product  = True
#                          WarmUp       = 1
                          ) 
                                       
