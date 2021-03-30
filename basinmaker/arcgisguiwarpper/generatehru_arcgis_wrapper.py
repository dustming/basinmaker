import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from basinmaker.postprocessing.hruarcgis import GenerateHRUS_arcgis


Path_Subbasin_Ply = arcpy.GetParameterAsText(0)
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
Landuse_ID = sys.argv[12]
Soil_ID = sys.argv[13]
Veg_ID = sys.argv[14]
Other_Ply_ID_1 = sys.argv[15]
Other_Ply_ID_2 = sys.argv[16]
DEM = sys.argv[17]
Inmportance_order = sys.argv[18]
min_hru_area_pct_sub = float(sys.argv[19])
Project_crs = sys.argv[20]
OutputFolder = sys.argv[21]

arcpy.AddMessage(Path_Subbasin_Ply)
arcpy.AddMessage(Path_Connect_Lake_ply)
arcpy.AddMessage(Inmportance_order)

if Other_Ply_ID_1 == '#':
    Other_Ply_ID_1="O_ID_1"
    
if Path_Other_Ply_2 == '#':
    Other_Ply_ID_2="O_ID_2"
    
Inmportance_order = Inmportance_order.split(";")

Sub_Lake_ID = 'HyLakeId',
Sub_ID = 'SubId',
Lake_Id = 'Hylak_id',

GenerateHRUS_arcgis(
    Path_Subbasin_Ply = Path_Subbasin_Ply,
    Landuse_info = Landuse_info,
    Soil_info = Soil_info,
    Veg_info = Veg_info,
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
    DEM=DEM,
    Project_crs = Project_crs,
    Inmportance_order = Inmportance_order,
    min_hru_area_pct_sub = min_hru_area_pct_sub,
    OutputFolder = OutputFolder,
)