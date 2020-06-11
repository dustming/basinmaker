import time 
from joblib import Parallel, delayed 
import os
import pandas as pd



###Define inputs for toolbox
Datafolder = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Data_Lakeofwoods/"
Outputfolder = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Outputs1/"

###Input Dataset for the whole watersheds
in_wd = Datafolder + "narivs_new.shp"
in_lake = Datafolder +  "HyLake2.shp"
in_obs = Datafolder + "obsfinal.shp"
landuse = Datafolder + "landuse2"
landuseinfo = Datafolder + "landuse.csv"
CA_HYDAT =os.path.join(Datafolder,'Hydat.sqlite3')
Out_Sub_Reg_Dem_Folder = os.path.join(Datafolder,'SubRegion')

####define functtion
def Generaterouting_subregion(i):
    SubReg_info = pd.read_csv(os.path.join(Out_Sub_Reg_Dem_Folder,'Sub_reg_info.csv')) 
    basinid = int(SubReg_info['Sub_Reg_ID'].values[i])            
    ProjectName = SubReg_info['ProjectNM'].values[i]
    Path_Sub_Polygon = os.path.join(Out_Sub_Reg_Dem_Folder,SubReg_info['Ply_Name'].values[i])
    print("Processing for subregion       ", basinid)
    file_object = open(os.path.join(Out_Sub_Reg_Dem_Folder,'Log'), 'a')
	## run for each sub region
    try:
        RTtool=LRRT(WidDep = in_wd,Lakefile = in_lake,Landuse = landuse,Landuseinfo = landuseinfo,obspoint = in_obs,OutputFolder = Outputfolder, ProjectNM = ProjectName,Path_Sub_Reg_Out_Folder = Out_Sub_Reg_Dem_Folder,Is_Sub_Region = 1)
        RTtool.Generatmaskregion(Path_Sub_Polygon = Path_Sub_Polygon)
        RTtool.Generateinputdata()
        RTtool.WatershedDiscretizationToolset()
        RTtool.AutomatedWatershedsandLakesFilterToolset(Thre_Lake_Area_Connect = 5,Thre_Lake_Area_nonConnect = 5)
        RTtool.RoutingNetworkTopologyUpdateToolset_riv('EPSG:3573',Outlet_Obs_ID = basinid + 1000)
        RTtool.Output_Clean(clean = 'False')
        Datafolder_final = os.path.join(Outputfolder,ProjectName)
        RTtool.Define_Final_Catchment(Datafolder = Datafolder_final)
        file_object.close()
    except:
        shutil.rmtree(os.path.join(Outputfolder,ProjectName),ignore_errors=True)
        file_object.write(str(basinid) + '          Failed ' + '\n')
        file_object.close()
        pass

	

SubReg_info_main = pd.read_csv(os.path.join(Out_Sub_Reg_Dem_Folder,'Sub_reg_info.csv'))

file_object = open(os.path.join(Out_Sub_Reg_Dem_Folder,'Log'), 'w')
file_object.write('SubRegion Id,           Status' + '\n')
file_object.close()

arg_instances = range(0,len(SubReg_info_main))

Parallel(n_jobs=4, verbose=1, backend="threading")(delayed(Generaterouting_subregion)(i) for i in arg_instances)
