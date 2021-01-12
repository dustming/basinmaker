def Generatesubdomain(
    self,
    Min_Num_Domain=9,
    Max_Num_Domain=13,
    Initaial_Acc=5000,
    Delta_Acc=1000,
    Out_Sub_Reg_Dem_Folder="#",
    ProjectNM="Sub_Reg",
    CheckLakeArea=1,
    Acc_Thresthold_stream=500,
    max_memory=2048,
):

    import grass.script as grass
    import grass.script.setup as gsetup
    from grass.pygrass.modules import Module
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.script import array as garray
    from grass.script import core as gcore
    from grass_session import Session

    if not os.path.exists(Out_Sub_Reg_Dem_Folder):
        os.makedirs(Out_Sub_Reg_Dem_Folder)

    #### Determine Sub subregion without lake
    os.environ.update(
        dict(GRASS_COMPRESS_NULLS="1", GRASS_COMPRESSOR="ZSTD", GRASS_VERBOSE="1")
    )
    PERMANENT = Session()
    PERMANENT.open(
        gisdb=self.grassdb, location=self.grass_location_geo, create_opts=""
    )
    N_Basin = 0
    Acc = Initaial_Acc
    print("##############################Loop for suitable ACC ")
    while N_Basin < Min_Num_Domain or N_Basin > Max_Num_Domain:
        grass.run_command(
            "r.watershed",
            elevation="dem",
            flags="sa",
            basin="testbasin",
            drainage="dir_grass_reg",
            accumulation="acc_grass_reg2",
            threshold=Acc,
            overwrite=True,
        )
        strtemp_array = garray.array(mapname="testbasin")
        N_Basin = np.unique(strtemp_array)
        N_Basin = len(N_Basin[N_Basin > 0])
        print(
            "Number of Subbasin:    ",
            N_Basin,
            "Acc  value:     ",
            Acc,
            "Change of ACC ",
            Delta_Acc,
        )
        if N_Basin > Max_Num_Domain:
            Acc = Acc + Delta_Acc
        if N_Basin < Min_Num_Domain:
            Acc = Acc - Delta_Acc
    PERMANENT.close()

    ##### Determine subregion with lake
    try:
        self.Generateinputdata(Is_divid_region=1)
    except:
        print(
            "Check print infomation, some unknown error occured, may have no influence on result"
        )
        pass

    try:
        self.WatershedDiscretizationToolset(
            accthresold=Acc, Is_divid_region=1, max_memroy=max_memory
        )
    except:
        print(
            "Check print infomation, some unknown error occured, may have no influence on result"
        )
        pass

    try:
        self.AutomatedWatershedsandLakesFilterToolset(
            Thre_Lake_Area_Connect=CheckLakeArea,
            Thre_Lake_Area_nonConnect=-1,
            Is_divid_region=1,
            max_memroy=max_memory,
        )
    except:
        print(
            "Check print infomation, some unknown error occured, may have no influence on result"
        )
        pass

    ####Determin river network for whole watersheds
    PERMANENT = Session()
    PERMANENT.open(
        gisdb=self.grassdb, location=self.grass_location_geo, create_opts=""
    )
    grass.run_command(
        "r.stream.extract",
        elevation="dem",
        accumulation="acc_grass",
        threshold=Acc_Thresthold_stream,
        stream_raster="Sub_Reg_str_grass_r",
        stream_vector="Sub_Reg_str_grass_v",
        overwrite=True,
        memory=max_memory,
    )

    #### export outputs
    grass.run_command(
        "r.pack",
        input="ndir_grass",
        output=self.Path_Sub_reg_grass_dir,
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input="ndir_Arcgis",
        output=self.Path_Sub_reg_arcgis_dir,
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input="acc_grass",
        output=self.Path_Sub_reg_grass_acc,
        overwrite=True,
    )
    grass.run_command(
        "r.pack", input="dem", output=self.Path_Sub_reg_dem, overwrite=True
    )
    grass.run_command(
        "v.pack",
        input="Sub_Reg_str_grass_v",
        output=self.Path_Sub_reg_grass_str_v,
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input="Sub_Reg_str_grass_r",
        output=self.Path_Sub_reg_grass_str_r,
        overwrite=True,
    )
    PERMANENT.close()
    return

def Generatesubdomainmaskandinfo(
    self, Out_Sub_Reg_Dem_Folder="#", ProjectNM="Sub_Reg"
):
    #### generate subbregion outlet points and subregion info table
    QgsApplication.setPrefixPath(self.qgisPP, True)
    Qgs = QgsApplication([], False)
    Qgs.initQgis()
    from processing.core.Processing import Processing
    from processing.tools import dataobjects
    from qgis import processing

    feedback = QgsProcessingFeedback()
    Processing.initialize()
    QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
    context = dataobjects.createContext()
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    import grass.script as grass
    import grass.script.setup as gsetup
    from grass.pygrass.modules import Module
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.script import array as garray
    from grass.script import core as gcore
    from grass_session import Session

    os.environ.update(
        dict(GRASS_COMPRESS_NULLS="1", GRASS_COMPRESSOR="ZSTD", GRASS_VERBOSE="-1")
    )
    PERMANENT = Session()
    PERMANENT.open(
        gisdb=self.grassdb, location=self.grass_location_geo, create_opts=""
    )
    grass.run_command("r.mask", raster="dem", maskcats="*", overwrite=True)
    grass.run_command("r.null", map="finalcat", setnull=-9999)

    strtemp_array = garray.array(mapname="finalcat")
    #        mask          = np.isin(sl_lake, Un_selectedlake_ids)

    dir = garray.array(mapname="ndir_Arcgis")
    acc = garray.array(mapname="acc_grass")
    Cat_outlets = copy.copy(strtemp_array)
    Cat_outlets[:, :] = -9999
    Cat_outlets_Down = copy.copy(strtemp_array)
    Cat_outlets_Down[:, :] = -9999
    ncols = int(strtemp_array.shape[1])
    nrows = int(strtemp_array.shape[0])
    Basins = np.unique(strtemp_array)

    Basins = Basins[Basins > 0]

    subregin_info = pd.DataFrame(
        np.full(len(Basins), -99999), columns=["Sub_Reg_ID"]
    )
    subregin_info["Dow_Sub_Reg_Id"] = -9999
    subregin_info["ProjectNM"] = -9999
    subregin_info["Nun_Grids"] = -9999
    subregin_info["Ply_Name"] = -9999
    subregin_info["Max_ACC"] = -9999
    for i in range(0, len(Basins)):
        basinid = int(Basins[i])
        grass.run_command("r.mask", raster="dem", maskcats="*", overwrite=True)
        exp = (
            "dem_reg_"
            + str(basinid)
            + "= if(finalcat == "
            + str(basinid)
            + ",dem, -9999)"
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)
        grass.run_command(
            "r.null", map="dem_reg_" + str(basinid), setnull=[-9999, 0]
        )

        ####define mask
        grass.run_command(
            "r.mask", raster="dem_reg_" + str(basinid), maskcats="*", overwrite=True
        )
        grass.run_command(
            "r.out.gdal",
            input="MASK",
            output=os.path.join(self.tempfolder, "Mask1.tif"),
            format="GTiff",
            overwrite=True,
        )
        processing.run(
            "gdal:polygonize",
            {
                "INPUT": os.path.join(self.tempfolder, "Mask1.tif"),
                "BAND": 1,
                "FIELD": "DN",
                "EIGHT_CONNECTEDNESS": False,
                "EXTRA": "",
                "OUTPUT": os.path.join(
                    self.tempfolder, "HyMask_region_" + str(basinid) + ".shp"
                ),
            },
        )
        processing.run(
            "gdal:dissolve",
            {
                "INPUT": os.path.join(
                    self.tempfolder, "HyMask_region_" + str(basinid) + ".shp"
                ),
                "FIELD": "DN",
                "OUTPUT": os.path.join(
                    self.tempfolder, "HyMask_region_f" + str(basinid) + ".shp"
                ),
            },
        )
        processing.run(
            "gdal:dissolve",
            {
                "INPUT": os.path.join(
                    self.tempfolder, "HyMask_region_" + str(basinid) + ".shp"
                ),
                "FIELD": "DN",
                "OUTPUT": os.path.join(
                    Out_Sub_Reg_Dem_Folder,
                    "HyMask_region_"
                    + str(int(basinid + self.maximum_obs_id))
                    + "_nobuffer.shp",
                ),
            },
        )
        processing.run(
            "native:buffer",
            {
                "INPUT": os.path.join(
                    self.tempfolder, "HyMask_region_f" + str(basinid) + ".shp"
                ),
                "DISTANCE": 0.005,
                "SEGMENTS": 5,
                "END_CAP_STYLE": 0,
                "JOIN_STYLE": 0,
                "MITER_LIMIT": 2,
                "DISSOLVE": True,
                "OUTPUT": os.path.join(
                    Out_Sub_Reg_Dem_Folder,
                    "HyMask_region_"
                    + str(int(basinid + self.maximum_obs_id))
                    + ".shp",
                ),
            },
        )

    grass.run_command("r.mask", raster="dem", maskcats="*", overwrite=True)

    for i in range(0, len(Basins)):
        basinid = int(Basins[i])
        catmask = strtemp_array == basinid
        catacc = acc[catmask]
        trow, tcol = Getbasinoutlet(basinid, strtemp_array, acc, dir, nrows, ncols)
        k = 1
        ttrow, ttcol = trow, tcol
        dowsubreginid = -1
        ### get the downstream catchment id
        while dowsubreginid < 0 and k < 20:
            nrow, ncol = Nextcell(dir, ttrow, ttcol)
            if nrow < 0 or ncol < 0:
                dowsubreginid = -1
                Cat_outlets[trow, tcol] = int(
                    basinid + self.maximum_obs_id
                )  ### for outlet of watershed, use trow tcol
                break
            elif nrow >= nrows or ncol >= ncols:
                dowsubreginid = -1
                Cat_outlets[trow, tcol] = int(basinid + self.maximum_obs_id)
                break
            elif (
                strtemp_array[nrow, ncol] <= 0
                or strtemp_array[nrow, ncol] == basinid
            ):
                dowsubreginid = -1
                Cat_outlets[trow, tcol] = int(basinid + self.maximum_obs_id)
            else:
                dowsubreginid = strtemp_array[nrow, ncol]
                Cat_outlets[ttrow, ttcol] = int(basinid + self.maximum_obs_id)
                Cat_outlets_Down[nrow, ncol] = int(basinid + self.maximum_obs_id)
            k = k + 1
            ttrow = nrow
            ttcol = ncol

        subregin_info.loc[i, "ProjectNM"] = (
            ProjectNM + "_" + str(int(basinid + self.maximum_obs_id))
        )
        subregin_info.loc[i, "Nun_Grids"] = np.sum(catmask)
        subregin_info.loc[i, "Ply_Name"] = (
            "HyMask_region_" + str(int(basinid + self.maximum_obs_id)) + ".shp"
        )
        subregin_info.loc[i, "Max_ACC"] = np.max(np.unique(catacc))
        subregin_info.loc[i, "Dow_Sub_Reg_Id"] = int(
            dowsubreginid + self.maximum_obs_id
        )
        subregin_info.loc[i, "Sub_Reg_ID"] = int(basinid + self.maximum_obs_id)

    ### remove subregion do not contribute to the outlet
    ## find watershed outlet subregion
    outlet_reg_info = subregin_info.loc[
        subregin_info["Dow_Sub_Reg_Id"] == self.maximum_obs_id - 1
    ]
    outlet_reg_info = outlet_reg_info.sort_values(by="Max_ACC", ascending=False)
    outlet_reg_id = outlet_reg_info["Sub_Reg_ID"].values[0]

    mask = subregin_info["Dow_Sub_Reg_Id"] == self.maximum_obs_id - 1
    mask2 = np.logical_not(subregin_info["Sub_Reg_ID"] == outlet_reg_id)
    del_row_mask = np.logical_and(mask2, mask)

    subregin_info = subregin_info.loc[np.logical_not(del_row_mask), :]
    subregin_info.to_csv(
        os.path.join(Out_Sub_Reg_Dem_Folder, "Sub_reg_info.csv"),
        index=None,
        header=True,
    )

    ####### save outlet of each subregion
    grass.run_command("r.mask", raster="dem", maskcats="*", overwrite=True)
    temparray = garray.array()
    temparray[:, :] = Cat_outlets[:, :]
    temparray.write(mapname="Sub_Reg_Outlets", overwrite=True)
    grass.run_command(
        "r.mapcalc",
        expression="Sub_Reg_Outlets_1 = int(Sub_Reg_Outlets)",
        overwrite=True,
    )
    grass.run_command("r.null", map="Sub_Reg_Outlets_1", setnull=-9999)

    temparray = garray.array()
    temparray[:, :] = Cat_outlets_Down[:, :]
    temparray.write(mapname="Sub_Reg_Outlets_Down", overwrite=True)
    grass.run_command(
        "r.mapcalc",
        expression="Sub_Reg_Outlets_Down_1 = int(Sub_Reg_Outlets_Down)",
        overwrite=True,
    )
    grass.run_command("r.null", map="Sub_Reg_Outlets_Down_1", setnull=-9999)

    grass.run_command(
        "r.to.vect",
        input="Sub_Reg_Outlets_1",
        output="Sub_Reg_Outlets_point",
        type="point",
        overwrite=True,
    )
    grass.run_command(
        "v.out.ogr",
        input="Sub_Reg_Outlets_point",
        output=os.path.join(Out_Sub_Reg_Dem_Folder, "Sub_Reg_Outlets.shp"),
        format="ESRI_Shapefile",
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input="Sub_Reg_Outlets_1",
        output=self.Path_Sub_reg_outlets_r,
        overwrite=True,
    )
    grass.run_command(
        "v.pack",
        input="Sub_Reg_Outlets_point",
        output=self.Path_Sub_reg_outlets_v,
        overwrite=True,
    )

    grass.run_command(
        "r.to.vect",
        input="Sub_Reg_Outlets_Down_1",
        output="Sub_Reg_Outlets_Down_point",
        type="point",
        overwrite=True,
    )
    grass.run_command(
        "v.out.ogr",
        input="Sub_Reg_Outlets_Down_point",
        output=os.path.join(Out_Sub_Reg_Dem_Folder, "Sub_Reg_Outlets_Down.shp"),
        format="ESRI_Shapefile",
        overwrite=True,
    )

    return

############################################################################