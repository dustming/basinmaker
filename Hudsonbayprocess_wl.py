import arcpy
import pandas as pd
import os
na_hyshedply = "C:/Users/m43han/Documents/Routing/Sample/Database/Shapefiles/hybas_na_lev01-12_v1c/hybas_na_lev12_v1c.shp"
na_hydem = "C:/Users/m43han/Documents/Routing/Sample/Database/Rasters/dem15s"
na_dir = "C:/Users/m43han/Documents/Routing/Sample/Database/Rasters/dir15s"
na_acc = "C:/Users/m43han/Documents/Routing/Sample/Database/Rasters/acc15s"
in_wd = "C:/Users/m43han/Documents/Routing/Sample/Database/Shapefiles/Width_Depth/narivs.shp"
na_shddbf = "C:/Users/m43han/Documents/Routing/Sample/Database/Shapefiles/hybas_na_lev01-12_v1c/hybas_na_lev12_v1c.dbf"
in_lakevol = "0"
in_lake = "C:/Users/m43han/Documents/Routing/Sample/Database/Shapefiles/HyLake.shp"
in_obs = "C:/Users/m43han/Documents/Routing/Sample/Database/Obspoints/obs_hu_wgs84_da.shp"
in_output = "C:/Users/m43han/Documents/Routing/Sample/Outputs/"

level= "08"  #"12"
islake="wl" #"nola"
list = pd.read_csv("C:/Users/m43han/Documents/Routing/Code/Routing/Code/Toolbox/basins.csv",sep=",",low_memory=False)
for i in range(0,100):
	in_tarid = int(list.ix[i]['Outid'])
	tttt = in_output+level+islake+str(in_tarid)+'/'+'finalcat_info.shp'
	if list.ix[i]['isna'] == 1:# and not os.path.isfile(tttt):
		in_outputf = in_output+level+islake+str(in_tarid)
		print in_outputf
		if not os.path.exists(in_outputf):
			os.makedirs(in_outputf)
		arcpy.ImportToolbox("C:/Users/m43han/Documents/Routing/Code/Routing/Code/Toolbox/RAVEN Tool Routing Box.tbx")
		arcpy.Generateinputs(str(in_tarid),na_hyshedply,na_hydem,na_dir,na_acc,in_wd,na_shddbf,in_lakevol,in_lake,in_obs,in_outputf)
		in_fishd ="C:/Users/m43han/Documents/Routing/Sample/Database/Shapefiles/hybas_na_lev01-12_v1c/hybas_na_lev"+level+"_v1c.shp"
		arcpy.GenFirstcatchments(in_outputf,in_lakevol,in_fishd)
		arcpy.AddLakeandObs(in_outputf,in_lakevol)
		z_factor = "0.00001395"	
		arcpy.Generatecatchmentinformation(in_outputf,z_factor)
		Fgrid = "C:/Users/m43han/Documents/Routing/Sample/Inputs/forcingncgrid.asc"
		ravenout = in_outputf +"/"+"RavenInput"
		if not os.path.exists(ravenout):
			os.makedirs(ravenout)
		arcpy.GenerateRavenInputs(in_outputf,ravenout,"100",Fgrid)

