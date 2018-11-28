


def GenerateGeoClass(geoclass,landuseinfo):
    geoclass['Main_cropid'] = 'NA'
    geoclass['Second_crop_id'] = 'NA'
    geoclass['Crop_rotation'] = 'NA'
    geoclass['Veg_type'] = 'NA'
    geoclass['Spec_Code'] = 'NA'
    geoclass['Tile_depth'] = 0
    geoclass['Stream_depth'] = 0
    geoclass['Num_Soil_Layer'] = 'NA'
    geoclass['Dep_Soil_Layer1'] = 'NA'
    geoclass['Dep_Soil_Layer2'] = 'NA'
    geoclass['Dep_Soil_Layer3'] = 'NA'
    waterid = landuseinfo[landuseinfo['Landuse'] == 'Water']['LandID'].values[0]
    inlakeid = int(max(landuseinfo['LandID_Model'].values)) + 1
    olakeid = int(max(landuseinfo['LandID_Model'].values)) + 2
    for i in range(0,len(geoclass.index)):
        idx = geoclass.index[i]
        landid = geoclass.loc[idx,'LandID']
        ###################################### for olake and ilake
        geoclass.loc[idx,'Spec_Code'] = 0
        if landid == inlakeid or landid == olakeid:
            geoclass.loc[idx,'Main_cropid'] = 0
            geoclass.loc[idx,'Veg_type'] = 3
            if landid == inlakeid:
                geoclass.loc[idx,'Spec_Code'] = 1
            else:
                geoclass.loc[idx,'Spec_Code'] = 2
            geoclass.loc[idx,'LandID'] = waterid
            geoclass.loc[idx,'Second_crop_id'] = 0
            geoclass.loc[idx,'Crop_rotation'] = 0
            geoclass.loc[idx,'Num_Soil_Layer'] = 1
            geoclass.loc[idx,'Dep_Soil_Layer1'] = 1
            geoclass.loc[idx,'Dep_Soil_Layer2'] = 1
            geoclass.loc[idx,'Dep_Soil_Layer3'] = 1
            continue
        ##############################################33
        ilandinfo = landuseinfo[landuseinfo['LandID'] == landid]
        if ilandinfo['Landuse'].values[0] == 'Water' or ilandinfo['Landuse'].values[0] == 'Urban':
            geoclass.loc[idx,'Main_cropid'] = 0
        else:
            geoclass.loc[idx,'Main_cropid'] = landid
        geoclass.loc[idx,'Second_crop_id'] = 0
        if ilandinfo['Landuse'].values[0] != 'Agri':
            geoclass.loc[idx,'Crop_rotation'] = 0
        else:
            geoclass.loc[idx,'Crop_rotation'] = landid
        if  ilandinfo['Landuse'].values[0] == 'Water':
            geoclass.loc[idx,'Veg_type'] = 3
        elif ilandinfo['Landuse'].values[0] == 'Forest':
            geoclass.loc[idx,'Veg_type'] = 2
        else:
            geoclass.loc[idx,'Veg_type'] = 1

        if  ilandinfo['Landuse'].values[0] == 'Water':
            geoclass.loc[idx,'Num_Soil_Layer'] = 1
            geoclass.loc[idx,'Dep_Soil_Layer1'] = 1
            geoclass.loc[idx,'Dep_Soil_Layer2'] = 1
            geoclass.loc[idx,'Dep_Soil_Layer3'] = 1
        else:
            geoclass.loc[idx,'Num_Soil_Layer'] = 3
            geoclass.loc[idx,'Dep_Soil_Layer1'] = 0.2
            geoclass.loc[idx,'Dep_Soil_Layer2'] = 1
            geoclass.loc[idx,'Dep_Soil_Layer3'] = 2
    return geoclass

def GenerateLakeClass(LakeDataclasss,catinfo):
    cats = np.unique(catinfo['SubId'].values)
    idx = 0
    for i in range(0,len(cats)):
        catid = cats[i]
        icat = catinfo[catinfo['SubId'] == catid]
        if icat['IsLake'].values[0] < 0:
            continue
        LakeDataclasss.loc[idx,'lakedataid'] =  icat['HyLakeId'].values[0]
        LakeDataclasss.loc[idx,'lakeid'] =  0
        LakeDataclasss.loc[idx,'ldtype'] =  1
        LakeDataclasss.loc[idx,'lake_depth'] =  icat['LakeDepth'].values[0]
        LakeDataclasss.loc[idx,'area'] =  icat['LakeArea'].values[0]
        LakeDataclasss.loc[idx,'w0ref'] =  icat['LakeDepth'].values[0]
        LakeDataclasss.loc[idx,'rate'] =  1
        LakeDataclasss.loc[idx,'exp'] =  1
        idx = idx +1
#        LakeDataclasss['qprod1'] =  1    #m3/s
#        LakeDataclasss['qprod2'] =  100   #m3/s
#        LakeDataclasss['datum1'] =  0101
#        LakeDataclasss['datum2'] =  1001
    return LakeDataclasss

def GenerateCropclass(landuseinfo,CropDataclass):
    idx = 0
    for i in range(0,len(landuseinfo)):
        if landuseinfo['Landuse'].values[i] == 'Water' or landuseinfo['Landuse'].values[i] == 'Urban':
            continue
        else:
            CropDataclass.loc[idx,'cropid'] = landuseinfo['LandID'].values[i]
            CropDataclass.loc[idx,'reg'] = 1
            CropDataclass.loc[idx,'fn1'] = 9
            CropDataclass.loc[idx,'fp1'] = 9
            CropDataclass.loc[idx,'mn1'] = 9
            CropDataclass.loc[idx,'mp1'] = 9
            CropDataclass.loc[idx,'fday1'] = 9
            CropDataclass.loc[idx,'mday1'] = 9
            CropDataclass.loc[idx,'fdown1'] = 9
            CropDataclass.loc[idx,'mdown1'] = 9
            CropDataclass.loc[idx,'fn2'] = 9
            CropDataclass.loc[idx,'fp2'] = 9
            CropDataclass.loc[idx,'mn2'] = 9
            CropDataclass.loc[idx,'mp2'] = 9
            CropDataclass.loc[idx,'fday2'] = 9
            CropDataclass.loc[idx,'mday2'] = 9
            CropDataclass.loc[idx,'fdown2'] = 9
            CropDataclass.loc[idx,'mdown2'] = 9
            idx = idx + 1
    return CropDataclass



# In[176]:

def Generatepar(Outfolder,soilinfo,landuseinfo):
    opar = open(Outfolder+'par.txt',"w") #sep='\t
    tab = '\t'
##########Soil type dependened parameters
    wcfc = 'wcfc'  + tab ## fieled capicity ~~  -kp33 watercontent
    wcwp = 'wcwp'  + tab ## fieled capicity ~~  -kp1500 watercontent
    wcep = 'wcep'  + tab ## fieled capicity ~~  -kp1500 watercontent
    mperc1 = 'mperc1'  + tab ## fieled capicity ~~  KSAT
    mperc2 = 'mperc2'  + tab ## fieled capicity ~~  KSAT centimeters per hour (cm/h) ---> mm/day
    sfrost = 'sfrosts' + tab
    rrcs1 = 'rrcs1' + tab
    rrcs2 = 'rrcs2' + tab
    trrcs = 'trrcs' + tab
    srrate = 'srrate' + tab
    macrate = 'macrate' + tab
    mactrsm = 'mactrsm' + tab
    mactrinf = 'mactrinf' + tab
    for i in range(0,len(soilinfo)):
        wcfc = wcfc + str(round(soilinfo.ix[i]['KP33']/100.00, 2) - round(soilinfo.ix[i]['KP1500']/100.00, 2)) + tab
        wcwp = wcwp + str(round(soilinfo.ix[i]['KP1500']/100.00, 2)) + tab
        wcep = wcep + str(round(soilinfo.ix[i]['KP0']/100.00, 2) - round(soilinfo.ix[i]['KP33']/100.00, 2)) + tab
        mperc1 = mperc1 + str(round(soilinfo.ix[i]['KSAT']*24*10, 2)) + tab  # cm/h to mm/day
        mperc2 = mperc2 + str(round(soilinfo.ix[i]['KSAT']*24*10, 2)) + tab
        sfrost = sfrost + str(1.00)+tab
        rrcs1 = rrcs1 + str(round(soilinfo.ix[i]['KSAT']*24*10/10, 2)) + tab
        rrcs2 = rrcs2 + str(round(soilinfo.ix[i]['KSAT']*24*10/10, 2)) + tab
        trrcs = trrcs + str(round(soilinfo.ix[i]['KSAT']*24*10/10, 2)) + tab
        srrate = srrate + str(round(1, 2)) + tab
        macrate = macrate + str(round(0, 2)) + tab
        mactrinf = mactrinf + str(round(soilinfo.ix[i]['KSAT']*24*10/10, 2)) + tab
        mactrsm = mactrsm + str(round(0, 2)) + tab

    wcfc = wcfc + '\n'
    wcwp = wcwp + '\n'
    wcep = wcep + '\n'
    mperc1 = mperc1 + '\n'
    mperc2 = mperc2 + '\n'
#        sfrost = sfrost + '\n'
    rrcs1 = rrcs1 + '\n'
    rrcs2 = rrcs2 + '\n'
    trrcs = trrcs + '\n'
    srrate = srrate + '\n'
    macrate = macrate + '\n'
    mactrinf = mactrinf +'\n'
    mactrsm = mactrsm + '\n'

    opar.write(wcfc)
    opar.write(wcwp)
    opar.write(wcep)
    opar.write(mperc1)
    opar.write(mperc2)
    opar.write(rrcs1)
    opar.write(rrcs2)
    opar.write(trrcs)
#    opar.write(sfrost)
    opar.write(srrate)
    opar.write(macrate)
    opar.write(mactrinf)
    opar.write(mactrsm)
############landuse dependent paramters
    cmlt = 'cmlt'  + tab ## melting factor
    ttmp = 'ttmp' + tab
    frost = 'frost' + tab
    srrcs = 'srrcs' + tab
    for i in range(0,len(landuseinfo)+2):
        cmlt = cmlt + str(3.00) + tab
        ttmp = ttmp + str(0.00) + tab
        frost = frost + str(1.00) + tab
        srrcs = srrcs + str(0.8) + tab
    cmlt = cmlt +  '\n'
    ttmp = ttmp +  '\n'
    frost = frost + '\n'
    srrcs = srrcs + '\n'
    opar.write(cmlt)
    opar.write(ttmp)
    opar.write(srrcs)
#    opar.write(frost)

##########General paramters
    rrcs3 = 'rrcs3' + tab + str(0.2)+ '\n'
    ttpi = 'ttpi' +tab + str(1.00)+ '\n'
    ttpd = 'ttpd' + tab + str(1.00)+ '\n'
    rrcscorr = 'rrcscorr' +tab +  str(1.00) + '\n'
    gratk = 'gratk' + tab + str(1.00) + '\n'
    grata = 'grata' + tab + str(0.00) + '\n'
    gratp = 'gratp' + tab + str(1.00) + '\n'
    rivvel = 'rivvel' + tab + str(1.00) + '\n'
    damp = 'damp' + tab + str(0.20) + '\n'
    gldepi = 'gldepi' + tab + str(5.00) + '\n'
    rivvel = 'rivvel' + tab + str(1.00) + '\n'
    opar.write(ttpi)
    opar.write(ttpd)
    opar.write(rrcs3)
    opar.write(rrcscorr)
    opar.write(gratk)
    opar.write(grata)
    opar.write(gratp)
    opar.write(rivvel)
    opar.write(damp)
    opar.write(gldepi)
    opar.write(damp)
    opar.write(damp)
    opar.close()
# In[170]:

# coding: utf-8

# In[159]:


def creatpdstructure(soilinfo,landuseinfo,catinfo):
#### http://www.smhi.net/hype/wiki/doku.php?id=start:hype_file_reference:geodata.txt
    geodata = pd.DataFrame(np.full(maxsubnum,1),columns=['Drop'])
    landclass = pd.DataFrame(np.full(maxsubnum,1),columns=['Drop'])
    LakeDataclasss = pd.DataFrame(np.full(maxsubnum,1),columns=['Drop'])
    CropDataclass = pd.DataFrame(np.full(maxsubnum,1),columns=['Drop'])
    CropDataclass['cropid'] = 'NA'
    CropDataclass['reg'] = 'NA'
    CropDataclass['fn1'] = 'NA'
    CropDataclass['fp1'] = 'NA'
    CropDataclass['mn1'] = 'NA'
    CropDataclass['mp1'] = 'NA'
    CropDataclass['fday1'] = 'NA'
    CropDataclass['mday1'] = 'NA'
    CropDataclass['fdown1'] = 'NA'
    CropDataclass['mdown1'] = 'NA'
    CropDataclass['fn2'] = 'NA'
    CropDataclass['fp2'] = 'NA'
    CropDataclass['mn2'] = 'NA'
    CropDataclass['mp2'] = 'NA'
    CropDataclass['fday2'] = 'NA'
    CropDataclass['mday2'] = 'NA'
    CropDataclass['fdown2'] = 'NA'
    CropDataclass['mdown2'] = 'NA'

    landclass['slc_ID'] = 'NA'
    landclass['slc_nn'] = 'NA'
    landclass['LandID'] = 'NA'
    landclass['LandID_Model'] = 'NA'
    landclass['Soilid'] = 'NA'
    geodata['area'] = 'NA'
    geodata['subid'] = 'NA'
    geodata['maindown'] = 'NA'
    geodata['latitude'] = 'NA'
    geodata['region'] = 'NA'
    geodata['parreg'] = 'NA'
    geodata['wqparreg'] = 'NA'
    geodata['lakeregion'] = 'NA'
    geodata['ilregion'] = 'NA'
    geodata['olregion'] = 'NA'
    geodata['elev_mean'] = 'NA'
    geodata['slope_mean'] = 'NA'
    geodata['lakedataid'] = 'NA'
    geodata['rivlen'] = 'NA'
    LakeDataclasss['lakedataid'] = 'NA'
    LakeDataclasss['lakeid'] = 'NA'
    LakeDataclasss['ldtype'] = 'NA'
    LakeDataclasss['lake_depth'] = 'NA'
    LakeDataclasss['area'] = 'NA'
    LakeDataclasss['w0ref'] = 'NA'
    LakeDataclasss['rate'] = 'NA'
    LakeDataclasss['exp'] = 'NA'
#    LakeDataclasss['deltaw0'] = 'NA'
#    LakeDataclasss['qprod1'] = 'NA'
#    LakeDataclasss['qprod2'] = 'NA'
#    LakeDataclasss['datum1'] = 'NA'
#    LakeDataclasss['datum2'] = 'NA'
#    LakeDataclasss['qamp'] = 'NA'
#    LakeDataclasss['qpha'] = 'NA'
#    LakeDataclasss['regvol'] = 'NA'
#    LakeDataclasss['wamp'] = 'NA'
#    LakeDataclasss['maxQprod'] = 'NA'
#    LakeDataclasss['minflow'] = 'NA'
#    LakeDataclasss['obsflow'] = 'NA'
#    LakeDataclasss['prodpp'] = 'NA'
#    LakeDataclasss['prodsp'] = 'NA'
#    LakeDataclasss['Qmean'] = 'NA'
#    LakeDataclasss['tpmean'] = 'NA'
#    LakeDataclasss['tnmean'] = 'NA'
#    LakeDataclasss['tocmean'] = 'NA'
#    LakeDataclasss['limqprod'] = 'NA'
#    LakeDataclasss['sedon'] = 'NA'
#    LakeDataclasss['sedpp'] = 'NA'
#    LakeDataclasss['sedoc'] = 'NA'
#    LakeDataclasss['wprodn'] = 'NA'
#    LakeDataclasss['wprodp'] = 'NA'
#    LakeDataclasss['wprodc'] = 'NA'
#    LakeDataclasss['denitwl'] = 'NA'
#    LakeDataclasss['deeplake'] = 'NA'
#    LakeDataclasss['fastlake'] = 'NA'
#    LakeDataclasss['t2mix'] = 'NA'
    kk = 1
    idxwater = -9
    for i in range(0,len(landuseinfo['LandID'])):
        if landuseinfo['Landuse'][i] == 'Water':
#            colname = 'slc_Wate
#            geodata[colname] = 0.0
            idxwater = i
        else:
            for j in range(0,len(soilinfo['Soilid'])):
#                colname = 'slc_'+str(landuseinfo['LandID'][i])+'_'+str(soilinfo['Soilid'][j])
                colname = 'slc_' + str(kk)
                geodata[colname] = 0.0
                landclass.loc[kk,'slc_ID'] = kk
                landclass.loc[kk,'slc_nn'] = 'slc_' + str(kk)
                landclass.loc[kk,'LandID'] = int(landuseinfo['LandID'][i])
                landclass.loc[kk,'LandID_Model'] = int(landuseinfo['LandID_Model'][i])
                landclass.loc[kk,'Soilid'] = int(soilinfo['Soilid'][j])
                kk = kk + 1
    colname = 'slc_' + str(kk)
    geodata[colname] = 0.0
    landclass.loc[kk,'slc_ID'] = kk
    landclass.loc[kk,'slc_nn'] = 'slc_' + str(kk)
    landclass.loc[kk,'LandID'] = int(landuseinfo['LandID'][idxwater])
    landclass.loc[kk,'LandID_Model'] = int(landuseinfo['LandID_Model'][idxwater])
    landclass.loc[kk,'Soilid'] = 1
    kk = kk + 1
    colname = 'slc_' + str(kk)
    geodata[colname] = 0.0
    landclass.loc[kk,'slc_ID'] = kk
    landclass.loc[kk,'slc_nn'] = 'slc_' + str(kk)
    landclass.loc[kk,'LandID'] = int(max(landuseinfo['LandID_Model'].values))+1 #### incatchment lake
    landclass.loc[kk,'LandID_Model'] = int(max(landuseinfo['LandID_Model'].values))+1
    landclass.loc[kk,'Soilid'] = 1
    kk = kk + 1
    colname = 'slc_' + str(kk)
    geodata[colname] = 0.0
    landclass.loc[kk,'slc_ID'] = kk
    landclass.loc[kk,'slc_nn'] = 'slc_' + str(kk)
    landclass.loc[kk,'LandID'] = int(max(landuseinfo['LandID_Model'].values))+2   ### outcatchment lake
    landclass.loc[kk,'LandID_Model'] = int(max(landuseinfo['LandID_Model'].values))+2
    landclass.loc[kk,'Soilid'] = 1
#    geodata['grwdown'] = 'NA'
#    geodata['grwolake'] = 'NA'
    geodata['grwolake'] = 'NA'
    geodata['loc_tp'] = 'NA'
    geodata['loc_tn'] = 'NA'
    geodata['loc_vol'] = 'NA'
    geodata['loc_sp'] = 'NA'
    geodata['loc_in'] = 'NA'
    geodata['wetdep_n'] = 'NA'
    geodata['drydep_n1'] = 'NA'
    geodata['drydep_n2'] = 'NA'
    geodata['drydep_n3'] = 'NA'
    geodata['deploadn1'] = 'NA'
    geodata['deploadn2'] = 'NA'
    geodata['deploadn3'] = 'NA'
    geodata['deploadn4'] = 'NA'
    geodata['deploadn5'] = 'NA'
    geodata['deploadn6'] = 'NA'
    geodata['deploadn7'] = 'NA'
    geodata['deploadn8'] = 'NA'
    geodata['deploadn9'] = 'NA'
    geodata['deploadn10'] = 'NA'
    geodata['deploadn11'] = 'NA'
    geodata['deploadn12'] = 'NA'
    geodata['lrwet_area'] = 'NA'
    geodata['mrwet_area'] = 'NA'
    geodata['lrwet_dep'] = 'NA'
    geodata['mrwet_dep'] = 'NA'
    geodata['lrwet_part'] = 'NA'
    geodata['mrwet_part'] = 'NA'
    geodata['buffer'] = 'NA'
    geodata['petmodel'] = 'NA'
    geodata = geodata.drop(['Drop'], axis=1)
    LakeDataclasss = LakeDataclasss.drop(['Drop'], axis=1)
    CropDataclass = CropDataclass.drop(['Drop'], axis=1)
    return geodata,landclass,LakeDataclasss,CropDataclass
####################3
def writegeodata(catid,i,geodata,hrus,landuseinfo,landclass,soilinfo):
    cathrus = hrus[hrus['CATCHMENTS'] == catid]
    icat = catinfo[catinfo['SubId'] == catid]
    geodata.loc[i,'area'] = icat['Area2'].values[0]/1000.0/1000.0
    geodata.loc[i,'subid'] = catid
    geodata.loc[i,'maindown'] = icat['DowSubId'].values[0]
    geodata.loc[i,'latitude'] = icat['INSIDE_Y'].values[0]
    geodata.loc[i,'region'] = 1
    geodata.loc[i,'parreg'] = 1
    geodata.loc[i,'wqparreg'] = 1
    geodata.loc[i,'lakeregion'] = 1
    geodata.loc[i,'ilregion'] = 1
    geodata.loc[i,'olregion'] = 1
    geodata.loc[i,'elev_mean'] = icat['MeanElev'].values[0]
    geodata.loc[i,'slope_mean'] = icat['BasinSlope'].values[0]
    geodata.loc[i,'lakedataid'] = max(0,icat['HyLakeId'].values[0])
    geodata.loc[i,'rivlen'] = max(0,icat['Rivlen'].values[0])
    waterhruarea = 0.0
    otherarea = 0.0
    wetlandarea = 0.0
    forestarea = 0.0
    sumwt = 0.0
    geodata.loc[i,'area'] = sum(cathrus['COUNT'].values)*30*30/1000.0/1000.0
    catarea = sum(cathrus['COUNT'].values)*30*30/1000.0/1000.0
    cathrus = cathrus.sort_values(by=['COUNT'])

##### Deal with lake first   landclass water == water id
    inlakeid = int(max(landuseinfo['LandID_Model'].values)) + 1
    outlakeid = int(max(landuseinfo['LandID_Model'].values)) + 2
    waterid = int(landuseinfo[landuseinfo['Landuse'] == 'Water']['LandID'].values[0])
    if icat['HyLakeId'].values[0] > 0:
        olakecolnam = landclass[landclass['LandID'] == outlakeid]['slc_nn'].values[0]
        geodata.loc[i,olakecolnam] = float(icat['LakeArea'].values[0])/catarea
        sumwt = sumwt + float(icat['LakeArea'].values[0])/catarea

    for j in range(len(cathrus.index)):
        idx = cathrus.index[j]
        soilid = cathrus.ix[idx]['SOILCA']
        landuseid = cathrus.ix[idx]['LANDUSE']
        isoil = soilinfo[soilinfo['Soilid'] == soilid]
        ilanduse = landuseinfo[landuseinfo['LandID'] == landuseid]
        hruarea = cathrus.ix[idx]['COUNT']*cellsize*cellsize/1000.0/1000.0
        if ilanduse['Landuse'].values == "Water":
            waterhruarea = waterhruarea + cathrus.ix[idx]['COUNT']*cellsize*cellsize/1000.0/1000.0
            continue
        elif ilanduse['Landuse'].values == "Forest" or ilanduse['Landuse'].values == "WetLand":
            continue
        else:
            colname = landclass.loc[np.logical_and(landclass['LandID'] == landuseid,landclass['Soilid'] == soilid)]['slc_nn'].values[0]
            geodata.loc[i,colname] =  min(float(hruarea)/catarea,1-sumwt)
            sumwt = sumwt + min(float(hruarea)/catarea,1-sumwt)
    if sumwt < 1:
        for j in range(len(cathrus.index)):
            idx = cathrus.index[j]
            soilid = cathrus.ix[idx]['SOILCA']
            landuseid = cathrus.ix[idx]['LANDUSE']
            isoil = soilinfo[soilinfo['Soilid'] == soilid]
            ilanduse = landuseinfo[landuseinfo['LandID'] == landuseid]
            hruarea = cathrus.ix[idx]['COUNT']*cellsize*cellsize/1000.0/1000.0
            if ilanduse['Landuse'].values != "Forest":
                continue
            else:
                colname = landclass.loc[np.logical_and(landclass['LandID'] == landuseid,landclass['Soilid'] == soilid)]['slc_nn'].values[0]
                geodata.loc[i,colname] =  min(float(hruarea)/catarea,1-sumwt)
                sumwt = sumwt + min(float(hruarea)/catarea,1-sumwt)
    if sumwt < 1:
        for j in range(len(cathrus.index)):
            idx = cathrus.index[j]
            soilid = cathrus.ix[idx]['SOILCA']
            landuseid = cathrus.ix[idx]['LANDUSE']
            isoil = soilinfo[soilinfo['Soilid'] == soilid]
            ilanduse = landuseinfo[landuseinfo['LandID'] == landuseid]
            hruarea = cathrus.ix[idx]['COUNT']*cellsize*cellsize/1000.0/1000.0
            if ilanduse['Landuse'].values != "WetLand":
                continue
            else:
                colname = landclass.loc[np.logical_and(landclass['LandID'] == landuseid,landclass['Soilid'] == soilid)]['slc_nn'].values[0]
                geodata.loc[i,colname] =  min(float(hruarea)/catarea,1-sumwt)
                sumwt = sumwt + min(float(hruarea)/catarea,1-sumwt)

    if waterhruarea > 0 and sumwt < 1:  ###Last consider add water
        waterid = landuseinfo[landuseinfo['Landuse'] == 'Water']['LandID'].values[0]
        colname= landclass[landclass['LandID'] == waterid]['slc_nn'].values[0]
        geodata.loc[i,colname] = 1 - sumwt
    geodata.loc[i,'grwdown'] = 0.0
    geodata.loc[i,'grwolake'] = 0.0
    geodata.loc[i,'loc_tp'] = 0.0
    geodata.loc[i,'loc_tn'] = 0.0
    geodata.loc[i,'loc_vol'] = 0.0
    geodata.loc[i,'loc_sp'] = 0.0
    geodata.loc[i,'loc_in'] = 0.0
    geodata.loc[i,'wetdep_n'] = 0.0
    geodata.loc[i,'drydep_n1'] = 0.0
    geodata.loc[i,'drydep_n2'] = 0.0
    geodata.loc[i,'drydep_n3'] = 0.0
    geodata.loc[i,'deploadn1'] = 0.0
    geodata.loc[i,'deploadn2'] = 0.0
    geodata.loc[i,'deploadn3'] = 0.0
    geodata.loc[i,'deploadn4'] = 0.0
    geodata.loc[i,'deploadn5'] = 0.0
    geodata.loc[i,'deploadn6'] = 0.0
    geodata.loc[i,'deploadn7'] = 0.0
    geodata.loc[i,'deploadn8'] = 0.0
    geodata.loc[i,'deploadn9'] = 0.0
    geodata.loc[i,'deploadn10'] = 0.0
    geodata.loc[i,'deploadn11'] = 0.0
    geodata.loc[i,'deploadn12'] = 0.0
    geodata.loc[i,'lrwet_area'] = 0.0
    geodata.loc[i,'mrwet_area'] = 0.0
    geodata.loc[i,'lrwet_dep'] = 0.0
    geodata.loc[i,'mrwet_dep'] = 0.0
    geodata.loc[i,'lrwet_part'] = 0.0
    geodata.loc[i,'mrwet_part'] = 0.0
    geodata.loc[i,'buffer'] = 0.0
    geodata.loc[i,'petmodel'] = 0.0
    return geodata
##################################3
def Defcat(out,outletid):
    otsheds = np.full((1,1),outletid)
    Shedid = np.full((10000000,1),-99999999999999999)
    psid = 0
    rout = copy.copy(out)
    while len(otsheds) > 0:
        noutshd = np.full((10000000,1),-999999999999999)
        poshdid = 0
        for i in range(0,len(otsheds)):
            Shedid[psid] = otsheds[i]
            psid = psid + 1
            irow = np.argwhere(rout[:,1]==otsheds[i]).astype(int)
            for j in range(0,len(irow)):
                noutshd[poshdid] = rout[irow[j],0]
                poshdid = poshdid + 1
        noutshd = np.unique(noutshd)
        otsheds = noutshd[noutshd>=0]
    Shedid = Shedid[Shedid>=0]
    Shedid,idx = np.unique(Shedid,return_index=True)
    outshed = np.full((len(Shedid),2),-9)
    outshed[:,0] = Shedid
    outshed[:,1] = idx
    outshed2 = outshed[outshed[:,1].argsort()[::-1]]
#    print "!@#!@#@!#@!#@!#@!#@!#"
#    print outshed2
    return outshed2

#########################start the program
from arcpy import env
from arcpy.sa import *
import copy
import sys
import shutil
import os
import csv
import csv
from dbfpy import dbf
import pandas as pd
from shutil import copyfile
import ConfigParser
from simpledbf import Dbf5
import numpy as np
arcpy.env.overwriteOutput = True
####### Required parameters
countthreshold = 50
#WorkFolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Grand River Basin/'
WorkFolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/GrandRiverProject/Dataprocess/'
InputsFolder =  WorkFolder + 'Project/'
cellsize = 30
maxsubnum = 500
###################Read inputs
Ourfolder =WorkFolder +  "HYPEINPUTS/"
if not os.path.exists(Ourfolder):
    os.makedirs(Ourfolder)
tab = '\t'
soilinfo = pd.read_csv(InputsFolder+"Soilinfo.csv",sep=",",low_memory=False)
landuseinfo = pd.read_csv(InputsFolder+"Landuseinfo.csv",sep=",",low_memory=False)
dbf2 = Dbf5(WorkFolder+ "finalcat_info.dbf")
catinfo = dbf2.to_dataframe()
dbf = Dbf5(InputsFolder+ "HRU_COMBINE.dbf")
hrus = dbf.to_dataframe()
routinfo = catinfo[['SubId','DowSubId']].values
########
geodata,landclass,LakeDataclasss,CropDataclass = creatpdstructure(soilinfo,landuseinfo,catinfo)
cats = np.unique(catinfo['SubId'].values)
outlets = np.unique(catinfo[catinfo['DowSubId'] == -1]['SubId'].values)### get all outlet ids
for j in range(0,len(outlets)):
    cats = Defcat(routinfo,outlets[0])  ### return sorted catchment id begin from head watershed
    for i in range(0,len(cats)):
        catid = cats[i,0]
        geodata = writegeodata(catid,i,geodata,hrus,landuseinfo,landclass,soilinfo)
geodata = geodata[geodata['subid'] != 'NA']
landclass = landclass[landclass['slc_ID'] != 'NA']
geodata.dropna(axis='columns')
geodata.to_csv(Ourfolder+'GeoData.txt',sep='\t',index = None)
geoclass = copy.copy(landclass)
geoclass = geoclass.drop(['slc_nn','Drop'], axis=1)
geoclassnew = GenerateGeoClass(geoclass,landuseinfo)
geoclassnew=geoclassnew.drop(['LandID'], axis=1)
geoclassnew.dropna(axis='columns')
geoclassnew.to_csv(Ourfolder+'GeoClass.txt',sep='\t',index = None, header = None)
LakeDataclasss = GenerateLakeClass(LakeDataclasss,catinfo)
CropDataclass = GenerateCropclass(landuseinfo,CropDataclass)
LakeDataclasss.replace(["NaN"], np.nan, inplace = True)
LakeDataclasss = LakeDataclasss[LakeDataclasss['lakedataid'] != 'NA']
LakeDataclasss.dropna(axis='columns')
#print LakeDataclasss
LakeDataclasss.to_csv(Ourfolder+'LakeData.txt',sep='\t',index = None)
CropDataclass = CropDataclass[CropDataclass['cropid'] != 'NA']
CropDataclass.dropna(axis='columns')
CropDataclass.to_csv(Ourfolder+'CropData.txt',sep='\t',index = None)
Generatepar(Ourfolder,soilinfo,landuseinfo)

# In[158]:
