##################################################################3
def dbftocsv(filename,outname):
    if filename.endswith('.dbf'):
        print "Converting %s to csv" % filename
        csv_fn = outname
        with open(csv_fn,'wb') as csvfile:
            in_db = dbf.Dbf(filename)
            out_csv = csv.writer(csvfile)
            names = []
            for field in in_db.header.fields:
                names.append(field.name)
            out_csv.writerow(names)
            for rec in in_db:
                out_csv.writerow(rec.fieldData)
            in_db.close()
            print "Done..."
    else:
        print "Filename does not end with .dbf"

def Maphru2force(orank,cat,catinfo,fnrows,fncols,outfolder,forcinggrid,outFolderraven):
    arcpy.RasterToPoint_conversion(outfolder + "finalcat.asc", outfolder + "Finalcat_Point.shp", "VALUE")
    ExtractValuesToPoints(outfolder + "Finalcat_Point.shp", forcinggrid, outfolder + "MapForcing.shp",
                      "NONE", "VALUE_ONLY")
    dbftocsv(outfolder + "MapForcing.dbf",outfolder + "MapForcing.csv")
    Mapforcing = pd.read_csv(outfolder + "MapForcing.csv",sep=",",low_memory=False)
    ogridforc = open(outFolderraven+"GriddedForcings2.txt","w")
    ogridforc.write(":GridWeights" +"\n")
    ogridforc.write("   #      " +"\n")
    ogridforc.write("   # [# HRUs]"+"\n")
    sNhru = len(catinfo)
    ogridforc.write("   :NumberHRUs       "+ str(sNhru) + "\n")
    sNcell = fnrows*fncols
    ogridforc.write("   :NumberGridCells  "+str(sNcell)+"\n")
    ogridforc.write("   #            "+"\n")
    ogridforc.write("   # [HRU ID] [Cell #] [w_kl]"+"\n")
    ncncols = fncols
    ncnrows = fnrows
    tab = '       '
    for i in range(0,len(catinfo.index)):
        catid = int(catinfo.iloc[i]['SUBID'])
        catmapf = Mapforcing.loc[Mapforcing['GRID_CODE'] == catid]
        rankids = catmapf.values[:,2]
        ids = np.unique(rankids)
        sumwt = 0.0
        for j in range(0,len(ids)):
            StrGid = str(int(catid))+tab
            ncrowcol = np.argwhere(orank==ids[j])
            Strcellid = str((ncrowcol[0,0] * ncncols + ncrowcol[0,1]))+tab
            if len(ids) == 1:
                pesr = 1
            else:
                if j < len(ids) - 1:
                    pesr = float(len(rankids[np.argwhere(rankids == ids[j])]))/float(len(rankids))
                    sumwt = sumwt + pesr
#                print j,pesr,sumwt,float(len(rankids[np.argwhere(rankids == ids[j])])),float(len(rankids))
                else:
                    pesr = 1 - sumwt
#                print j,pesr,sumwt
            ogridforc.write("    "+StrGid+Strcellid+str(pesr) + "\n")
    ogridforc.write(":EndGridWeights")
    ogridforc.close()


def Writervhchanl(ocatinfo,outFolder,nrows,ncols,lenThres):
    catinfo = copy.copy(ocatinfo)
#    print int(catinfo.iloc[0]['SUBID']),len(catinfo.index)
    ochn = open(outFolder+"modelchannel.rvp","w")
##################3
    orvh = open(outFolder+"test.rvh","w")
    orvh.write("# --------------------------------------------"+"\n")
    orvh.write("# Raven HRU Input file"+"\n")
    orvh.write("#  lake catchment emulation"+"\n")
    orvh.write("# --------------------------------------------"+"\n")
    orvh.write(":SubBasins"+"\n")
    orvh.write("  :Attributes   NAME  DOWNSTREAM_ID       PROFILE REACH_LENGTH  GAUGED"+"\n")
    orvh.write("  :Units        none           none          none           km    none"+"\n")
    tab = "     "
    for i in range(0,len(catinfo.index)):
        ### Get catchment width and dpeth
        catid = int(catinfo.iloc[i]['NSUBID'])
        temp = catinfo.iloc[i]['RIVLEN']
        if (temp >= lenThres):
            catlen = float(temp)/1000 #### in km
            strRlen = str(catlen)
        else:
            catlen = -9999
            strRlen = 'ZERO-'
        #####################################################3
        Strcat = str(catid)
        StrDid = str(int(catinfo.iloc[i]['NDOWSUBID']))
        pronam = 'Chn_'+ Strcat
        chslope = catinfo.iloc[i]['RIVSLOPE']
        if chslope <= 0.0001:
            chslope = max(catinfo.iloc[i]['BASINSLOPE'],0.0001)
        writechanel(pronam,catinfo.iloc[i]['BKFWIDTH'],catinfo.iloc[i]['BKFDEPTH'],
                    chslope,ochn,max(catinfo.iloc[i]['MEANELEV'],300))
        if catinfo.iloc[i]['ISOBS'] >= 0 :
            Guage = '1'
        else:
            Guage = '0'
        orvh.write("  "+Strcat+tab+'sub'+Strcat+tab+StrDid+tab+pronam+tab+strRlen+tab+Guage+"\n")
    orvh.write(":EndSubBasins"+"\n")
    orvh.write("\n")
##########################################
    orvh.write(":HRUs"+"\n")
    orvh.write("  :Attributes AREA ELEVATION  LATITUDE  LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT"+"\n")
    orvh.write("  :Units       km2         m       deg        deg       none            none       none           none             none            none     deg      deg"+"\n")
    for i in range(0,len(catinfo.index)):
        catid = int(catinfo.iloc[i]['NSUBID'])
        catslope = max(catinfo.iloc[i]['BASINSLOPE'],0.0001)
        catarea2 = float(catinfo.iloc[i]['AREA2'])/1000.00/1000.00
        StrGid =  str(catid)+tab
        StrGidarea = str(catarea2)+tab
        StrGidelev = str(max(catinfo.iloc[i]['MEANELEV'],300))+tab
        lat = str(catinfo.iloc[i]['INSIDE_Y'])+tab
        lon = str(catinfo.iloc[i]['INSIDE_X'])+tab
        LAND_USE_CLASS = 'FOREST'+tab
        VEG_CLASS = 'FOREST'+tab
        SOIL_PROFILE ='SOILPROF'+tab
        AQUIFER_PROFILE ='[NONE]'+tab
        TERRAIN_CLASS ='[NONE]'+tab
        SLOPE = str(catslope)+tab
        ASPECT = '200'+tab
        orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+StrGid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")
    orvh.write(":EndHRUs"+"\n")
    orvh.write(":RedirectToFile TestLake.rvh")
    orvh.close()
    ochn.close()
    return catinfo
##############################

#########################################################
def writechanel(chname,chwd,chdep,chslope,orchnl,elev):
    ### Following SWAT instructions, assume a trapezoidal shape channel, with channel sides has depth and width ratio of 2. zch = 2
    zch = 2
    chslope = max(chslope,0.001)
    sidwd = zch * chdep ###river side width
    tab = "          "
    botwd = chwd - 2*sidwd ### river
    if (botwd < 0):
        botwd = 0.5*chwd
        sidwd = 0.5*0.5*chwd
    mann = "0.035"
    elev = max(elev,300)
    zfld = 4 + elev
    zbot = elev - chdep
    sidwdfp = 4/0.25
    Channame = ":ChannelProfile"+tab+chname+tab
    orchnl.write(Channame+"\n")
    Chanslop = "  :Bedslope"+tab+str(chslope)
    orchnl.write(Chanslop+"\n")
    orchnl.write("  :SurveyPoints"+"\n")
    orchnl.write("    0"+tab+str(zfld)+"\n")
    orchnl.write("    "+str(sidwdfp)+tab+str(elev)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd)+tab+str(elev)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd + sidwd)+tab+str(zbot)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd + sidwd + botwd)+tab+str(zbot)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd + 2*sidwd + botwd)+tab+str(elev)+"\n")
    orchnl.write("    "+str(sidwdfp + 4*chwd + 2*sidwd + botwd)+tab+str(elev)+"\n")
    orchnl.write("    "+str(2*sidwdfp + 4*chwd + 2*sidwd + botwd)+tab+str(zfld)+"\n")
    orchnl.write("  :EndSurveyPoints"+"\n")
    orchnl.write("  :RoughnessZones"+"\n")
    orchnl.write("    0" + tab + mann +"\n")
    orchnl.write("  :EndRoughnessZones"+"\n")
    orchnl.write(":EndChannelProfile"+"\n")
    orchnl.write("\n")
    orchnl.write("##############new channel ##############################\n")
#########################################################################################################33

def writelake(catinfo,outFolderraven,lakevol,WeirCoe):
    f2 = open(outFolderraven+"TestLake.rvh","w")
    tab = '       '
    for i in range(0,len(catinfo.index)):
        if catinfo.iloc[i]['HYLAKEID'] > 0 and catinfo.iloc[i]['LAKEAREA'] >= lakevol:
            lakeid = int(catinfo.iloc[i]['HYLAKEID'])
            catid = catinfo.iloc[i]['NSUBID']
            A = catinfo.iloc[i]['LAKEAREA']*1000*1000
            h0 = catinfo.iloc[i]['LAKEDEPTH']
            Crewd = catinfo.iloc[i]['BKFWIDTH']
#            if slakeinfo.iloc[0]['Wshd_area'] < 6000 and slakeinfo.iloc[0]['Wshd_area'] > 0:
        ######write lake information to file
            f2.write(":Reservoir"+ "   Lake_"+ str(int(lakeid))+ "   ######## " +"\n")
            f2.write("  :SubBasinID  "+str(int(catid))+ "\n")
            f2.write("  :HRUID   "+str(int(catid))+ "\n")
            f2.write("  :Type RESROUTE_STANDARD   "+"\n")
            f2.write("  :WeirCoefficient  "+str(WeirCoe)+ "\n")
            f2.write("  :CrestWidth "+str(Crewd)+ "\n")
            f2.write("  :MaxDepth "+str(h0)+ "\n")
            f2.write("  :LakeArea    "+str(A)+ "\n")
            f2.write(":EndReservoir   "+"\n")
            f2.write("#############################################"+"\n")
            f2.write("###New Lake starts"+"\n")
    f2.close()


import arcpy
import pandas as pd
import os
from dbfpy import dbf
import csv
import copy
import numpy as np
#in_output = "C:/Users/m43han/Documents/Routing/Sample/Outputs/"
in_output = "C:/Users/m43han/Documents/Routing/Sample/Outputs/panca/"

ilevel = "08"
level= "08"  #"12"
islake="wl/" #"nola" "wl"
lenThres = 2000
list = pd.read_csv("C:/Users/m43han/Documents/Routing/Code/Routing/Code/Toolbox/Basins_ca.csv",sep=",",low_memory=False)
arcpy.env.overwriteOutput = True
forcinggrid = "C:/Users/m43han/Documents/Routing/Sample/Inputs/forcingncgrid.asc"
tarshp =in_output+'finalout_'+islake+level+"/"
lakevol = 0
WeirCoe = 0.6
k = 0
if not os.path.exists(tarshp):
	os.makedirs(tarshp)
for i in range(0,len(list)):
	in_tarid = int(list.ix[i]['Outid'])
	in_outputf = in_output+ilevel+islake+str(in_tarid)+'/'+'finalcat_info.shp'
	in_outputf2 = in_output+ilevel+islake+str(in_tarid)+'/RavenInput/'+'test.rvh'
	print os.path.isfile(in_outputf2)
	if os.path.isfile(in_outputf2):
		if len(arcpy.ListFields(in_outputf,"BasinID"))<=0:
			arcpy.AddField_management(in_outputf, "BasinID", "LONG",10,"","", "", "NULLABLE")
		arcpy.CalculateField_management(in_outputf, "BasinID", str(k), "PYTHON")
		if len(arcpy.ListFields(in_outputf,"BasinID2"))<=0:
			arcpy.AddField_management(in_outputf, "BasinID2", "LONG",20,"","", "", "NULLABLE")
		arcpy.CalculateField_management(in_outputf, "BasinID2", str(int(in_tarid)), "PYTHON")
		if k == 0:
			arcpy.Copy_management(in_outputf, tarshp+'/'+'finalcat_info.shp')
            arcpy.Copy_management(in_outputf, tarshp+'/'+str(int(in_tarid))+'_finalcat'+'.shp')
		else:
			arcpy.Append_management(in_outputf,tarshp+'/'+'finalcat_info.shp','TEST')
            arcpy.Copy_management(in_outputf, tarshp+'/'+str(int(in_tarid))+'_finalcat'+'.shp')
		k = k + 1
dbftocsv(tarshp +"/"+ "finalcat_info.dbf",tarshp +"/"+ "finalcat_info.csv")
ncatinfo = pd.read_csv(tarshp +"/"+"finalcat_info.csv",sep=",",low_memory=False)
ncatinfo['NSUBID'] = -9
ncatinfo['NDOWSUBID'] = -9
baseid = 0
obasinid = 0
currcatid = 0
for i in range(len(ncatinfo)):
	if ncatinfo.ix[i]['BASINID'] == obasinid:
		ncatinfo['NSUBID'][i] = ncatinfo.ix[i]['SUBID'] + baseid
		if ncatinfo.ix[i]['DOWSUBID'] > 0:
			ncatinfo['NDOWSUBID'][i] = ncatinfo.ix[i]['DOWSUBID'] + baseid
		else:
			ncatinfo['NDOWSUBID'][i] = -1
		currcatid = ncatinfo.ix[i]['SUBID'] + baseid
	else:
		baseid = currcatid + 2
		ncatinfo['NSUBID'][i] =  ncatinfo.ix[i]['SUBID'] + baseid
		if ncatinfo.ix[i]['DOWSUBID'] > 0:
			ncatinfo['NDOWSUBID'][i] = ncatinfo.ix[i]['DOWSUBID'] + baseid
		else:
			ncatinfo['NDOWSUBID'][i] = -1
		obasinid = ncatinfo.ix[i]['BASINID']
		currcatid = ncatinfo.ix[i]['SUBID'] + baseid
ncatinfo.to_csv(tarshp +"/"+"newfinal_info.csv",sep=",")
if not os.path.exists(tarshp + "/"+'raveninput'+str(lakevol)+str(WeirCoe)+'/'):
    os.makedirs(tarshp + "/"+'raveninput'+str(lakevol)+str(WeirCoe)+'/')
Writervhchanl(ncatinfo,tarshp + "/"+'raveninput'+str(lakevol)+str(WeirCoe)+'/',0,0,lenThres)
writelake(ncatinfo,tarshp + "/"+'raveninput'+str(lakevol)+str(WeirCoe)+'/',lakevol,WeirCoe)



###########################################################################################################
#ogridforc = open(tarshp + "/"+'raveninput'+str(lakevol)+'/'+"GriddedForcings2.txt","w")
#orank = np.loadtxt(forcinggrid,dtype = 'i4',skiprows = 6)
#fncols = int(arcpy.GetRasterProperties_management(forcinggrid, "COLUMNCOUNT").getOutput(0))
#fnrows = int(arcpy.GetRasterProperties_management(forcinggrid, "ROWCOUNT").getOutput(0))

#ogridforc.write(":GridWeights" +"\n")
#ogridforc.write("   #      " +"\n")
#ogridforc.write("   # [# HRUs]"+"\n")
#sNhru = len(ncatinfo)
#ogridforc.write("   :NumberHRUs       "+ str(sNhru) + "\n")
#sNcell = fnrows*fncols
#ogridforc.write("   :NumberGridCells  "+str(sNcell)+"\n")
#ogridforc.write("   #            "+"\n")
#ogridforc.write("   # [HRU ID] [Cell #] [w_kl]"+"\n")
#ncncols = fncols
#ncnrows = fnrows
#tab = "      "

#for i in range(len(ncatinfo)):
#	in_tarid = int(ncatinfo.ix[i]['BASINID2'])
#	in_mapgrid = in_output+level+islake+str(in_tarid)+'/'+'RavenInput/'+'GriddedForcings2.txt'
#	Mapforcing = np.genfromtxt (in_mapgrid,skip_header =7,skip_footer =1)
#	nsubid = ncatinfo.ix[i]['NSUBID']
#	osubid = ncatinfo.ix[i]['SUBID']
#	ids = Mapforcing[Mapforcing[:,0] == osubid,:]
#	print in_tarid,osubid,nsubid
#	for j in range(len(ids)):
#		StrGid = str(int(nsubid)) + tab
#		Strcellid = str(int(ids[j,1])) + tab
#		pesr = ids[j,2]
#		ogridforc.write("    "+StrGid+Strcellid+str(pesr) + "\n")
#ogridforc.write(":EndGridWeights")
#ogridforc.close()
