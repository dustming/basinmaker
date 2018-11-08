import arcpy
import pandas as pd
import os
na_hyshedply = "C:/Users/m43han/Documents/Routing/Sample/Database/Shapefiles/hybas_ar_lev01-12_v1c/hybas_ar_lev12_v1c.shp"
na_hydem = "C:/Users/m43han/Documents/Routing/Sample/Database/Rasters/dem15s"
na_dir = "C:/Users/m43han/Documents/Routing/Sample/Database/Rasters/dir15s"
na_acc = "C:/Users/m43han/Documents/Routing/Sample/Database/Rasters/acc15s"
in_wd = "C:/Users/m43han/Documents/Routing/Sample/Database/Shapefiles/Width_Depth/narivs_new.shp"
na_shddbf = "C:/Users/m43han/Documents/Routing/Sample/Database/Shapefiles/hybas_ar_lev01-12_v1c/hybas_ar_lev12_v1c.dbf"
in_lake = "C:/Users/m43han/Documents/Routing/Sample/Database/Shapefiles/HyLake.shp"
in_obs = "C:/Users/m43han/Documents/Routing/Sample/Database/Obspoints/obs_ca.shp" #obs_hu_wgs84_da.shp"
in_output = "C:/Users/m43han/Documents/Routing/Sample/Outputs/"
OutHyID2 = "-1"
cellSize = "0.0041666667"
landuse = "C:/Users/m43han/Documents/Routing/Sample/Database/Rasters/landuse"
landuseinfo = "C:/Users/m43han/Documents/Routing/Sample/Database/Rasters/landuseinfo.csv"
InNonConL = "-1"
NonConLThres = "-1"

forcinggrid = "#"
forcingply = "#"
Boundaryply ="#"
missrow = "0"
misscol = "0"
iscalmanningn = "-1"

ilevel = "08"
level= "08"  #"12"
islake="wl" #"nola"
in_lakevol = "0"
list = pd.read_csv("C:/Users/m43han/Documents/Routing/Code/Routing/Code/Toolbox/basins_ca.csv",sep=",",low_memory=False)
for i in range(0,len(list)):
	in_tarid = int(list.ix[i]['Outid'])
	tttt = in_output+ilevel+islake+str(in_tarid)+'/'+'RavenInput/' + 'test.rvh'
	if list.ix[i]['isna'] == 0 and not os.path.isfile(tttt):
		in_outputf = in_output+ilevel+islake+str(in_tarid)
		print in_outputf
		if not os.path.exists(in_outputf):
			os.makedirs(in_outputf)
		arcpy.ImportToolbox("C:/Users/m43han/Documents/Routing/Code/Routing/Code/Toolbox/RAVEN Tool Routing Box.tbx")
		arcpy.Generateinputs(str(in_tarid),na_hyshedply,na_hydem,na_dir,na_acc,in_wd,na_shddbf,in_lakevol,in_lake,in_obs,in_outputf,cellSize,OutHyID2,landuse,landuseinfo)
		in_fishd ="C:/Users/m43han/Documents/Routing/Sample/Database/Shapefiles/hybas_ar_lev01-12_v1c/hybas_ar_lev"+level+"_v1c.shp"
		arcpy.GenFirstcatchments(in_outputf,in_lakevol,in_fishd,cellSize)
		arcpy.AddLakeandObs(in_outputf,in_lakevol,cellSize,InNonConL,NonConLThres)
		arcpy.Generatecatchmentinformation(in_outputf,cellSize)
#		Fgrid = "C:/Users/m43han/Documents/Routing/Sample/Inputs/forcingncgrid.asc"
		ravenout = in_outputf +"/"+"RavenInput"
		if not os.path.exists(ravenout):
			os.makedirs(ravenout)
		arcpy.GenerateRavenInputs(in_outputf,ravenout,"100",forcinggrid,forcingply,cellSize,Boundaryply,missrow,misscol,iscalmanningn)
