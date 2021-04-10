import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from hymodin.raveninput import GenerateRavenInput

Path_final_hru_info =sys.argv[1]
Startyear=-1
EndYear=-1
CA_HYDAT='#'
Lake_As_Gauge=0
Model_Name = sys.argv[2]
OutputFolder = sys.argv[3]

if Startyear != '#':
    Startyear = int(Startyear)
    DownLoadObsData = True
else:
    Startyear = -1
    DownLoadObsData = False

if EndYear != '#':
    EndYear = int(EndYear)
else:
    EndYear = -1 
    
#Inmportance_order = Inmportance_order.split(";")

GenerateRavenInput(
    Path_final_hru_info=Path_final_hru_info,
    lenThres=1,
    iscalmanningn=1,
    Startyear=Startyear,
    EndYear=EndYear,
    CA_HYDAT=CA_HYDAT,
    WarmUp=0,
    Template_Folder="#",
    Lake_As_Gauge=Lake_As_Gauge,
    WriteObsrvt=False,
    DownLoadObsData=DownLoadObsData,
    Model_Name=Model_Name,
    Old_Product=False,
    SubBasinGroup_NM_Channel=["Allsubbasins"],
    SubBasinGroup_Length_Channel=[-1],
    SubBasinGroup_NM_Lake=["AllLakesubbasins"],
    SubBasinGroup_Area_Lake=[-1],
    OutputFolder=OutputFolder,
    Forcing_Input_File="#",
    aspect_from_gis = "arcgis",
)
