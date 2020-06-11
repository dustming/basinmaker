import time 
from joblib import Parallel, delayed 
import os
import pandas as pd



####define functtion
def Generaterouting_subregion(i):
	os.system("python QGIS_script_for_Tool_paper.py  " + str(i))
	

###Read Sub Region info
Datafolder = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Data_Lakeofwoods/"
Out_Sub_Reg_Dem_Folder = os.path.join(Datafolder,'SubRegion') 
SubReg_info_main = pd.read_csv(os.path.join(Out_Sub_Reg_Dem_Folder,'Sub_reg_info.csv'))

arg_instances = range(0,len(SubReg_info_main))

#### Run
Parallel(n_jobs=4, verbose=1, backend="threading")(delayed(Generaterouting_subregion)(i) for i in arg_instances)
