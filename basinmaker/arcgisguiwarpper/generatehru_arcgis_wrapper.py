import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from postprocessing.hruarcgis import GenerateHRUS_arcgis


Path_Subbasin_Ply = sys.argv[1]
Path_Connect_Lake_ply = sys.argv[2]
Path_Non_Connect_Lake_ply = sys.argv[3]
Path_Landuse_Ply = sys.argv[4]
Path_Soil_Ply = sys.argv[5]
Path_Veg_Ply = sys.argv[6]
Path_Other_Ply_1 = sys.argv[7]
Path_Other_Ply_2 = sys.argv[8]
Landuse_info = sys.argv[9]
Soil_info = sys.argv[10]
Veg_info = sys.argv[11]
DEM = sys.argv[12]
Inmportance_order = sys.argv[13]
min_hru_area_pct_sub = float(sys.argv[14])
Project_crs = sys.argv[15]
OutputFolder = sys.argv[16]

Other_Ply_ID_1="O_ID_1"
Other_Ply_ID_2="O_ID_2"
Landuse_ID = 'Landuse_ID'
Soil_ID = 'Soil_ID'
Veg_ID = 'Veg_ID'

Sub_Lake_ID="HyLakeId"
Sub_ID="SubId"
Lake_Id="Hylak_id"
Landuse_ID="Landuse_ID"
Soil_ID="Soil_ID"
Other_Ply_ID_1="O_ID_1"
Veg_ID="Veg_ID"
Other_Ply_ID_2="O_ID_2"
            
Inmportance_order = Inmportance_order.split(";")

Sub_Lake_ID = 'HyLakeId',
Sub_ID = 'SubId',
Lake_Id = 'Hylak_id',

GenerateHRUS_arcgis(
    Path_Subbasin_Ply=Path_Subbasin_Ply,
    Landuse_info=Landuse_info,
    Soil_info=Soil_info,
    Veg_info=Veg_info,
    Sub_Lake_ID=Sub_Lake_ID,
    Sub_ID=Sub_ID,
    Path_Connect_Lake_ply=Path_Connect_Lake_ply,
    Path_Non_Connect_Lake_ply=Path_Non_Connect_Lake_ply,
    Lake_Id=Lake_Id,
    Path_Landuse_Ply=Path_Landuse_Ply,
    Landuse_ID=Landuse_ID,
    Path_Soil_Ply=Path_Soil_Ply,
    Soil_ID=Soil_ID,
    Path_Veg_Ply=Path_Veg_Ply,
    Veg_ID=Veg_ID,
    Path_Other_Ply_1=Path_Other_Ply_1,
    Other_Ply_ID_1=Other_Ply_ID_1,
    Path_Other_Ply_2=Path_Other_Ply_2,
    Other_Ply_ID_2=Other_Ply_ID_2,
    Inmportance_order = Inmportance_order,
    min_hru_area_pct_sub = min_hru_area_pct_sub,
    DEM=DEM,
    Project_crs=Project_crs,
    OutputFolder=OutputFolder,
)