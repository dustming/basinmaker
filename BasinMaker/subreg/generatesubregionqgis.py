from func.grassgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
from delineationnolake.watdelineationwithoutlake import (
    watershed_delineation_without_lake,
)
from addlakeandobs.addlakesqgis import (
    add_lakes_into_existing_watershed_delineation,
)
        
def Generatesubdomain(
    input_geo_names,
    grassdb,
    grass_location,
    qgis_prefix_path,
    path_lakefile_in,
    lake_attributes,
    Min_Num_Domain=9,
    Max_Num_Domain=13,
    Initaial_Acc=5000,
    Delta_Acc=1000,
    CheckLakeArea=1,
    fdr_path = '#',
    Acc_Thresthold_stream=500,
    max_memory=2048*3,
    Out_Sub_Reg_Folder="#",
    sub_reg_str_r = 'sub_reg_str_r',
    sub_reg_str_v = 'sub_reg_str_v',
    sub_reg_nfdr_grass = 'sub_reg_nfdr_grass',
    sub_reg_nfdr_arcgis = 'sub_reg_nfdr_arcgis',
    sub_reg_acc = 'sub_reg_acc',
    sub_reg_dem = 'sub_reg_dem',
):

    import grass.script as grass
    import grass.script.setup as gsetup
    from grass.pygrass.modules import Module
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.script import array as garray
    from grass.script import core as gcore
    from grass_session import Session

    if not os.path.exists(Out_Sub_Reg_Folder):
        os.makedirs(Out_Sub_Reg_Folder)
    
    #required inputs 
    dem = input_geo_names['dem']
    mask = input_geo_names['mask']
    
    # define local variable file names   
    fdr_arcgis = Internal_Constant_Names["fdr_arcgis"]
    fdr_grass = Internal_Constant_Names["fdr_grass"]
    nfdr_arcgis = Internal_Constant_Names["nfdr_arcgis"]
    nfdr_grass = Internal_Constant_Names["nfdr_grass"]
    str_r = Internal_Constant_Names["str_r"]
    str_v = Internal_Constant_Names["str_v"]
    acc = Internal_Constant_Names["acc"]
    cat_no_lake = Internal_Constant_Names["cat_no_lake"]
    sl_connected_lake = Internal_Constant_Names["sl_connected_lake"]
    sl_non_connected_lake = Internal_Constant_Names["sl_nonconnect_lake"]
    sl_lakes = Internal_Constant_Names["selected_lakes"]
    catchment_without_merging_lakes = Internal_Constant_Names["catchment_without_merging_lakes"]
    river_without_merging_lakes = Internal_Constant_Names["river_without_merging_lakes"]
    cat_use_default_acc = Internal_Constant_Names["cat_use_default_acc"]
    cat_add_lake = Internal_Constant_Names["cat_add_lake"] 
    pourpoints_with_lakes = Internal_Constant_Names["pourpoints_with_lakes"]
    pourpoints_add_obs = Internal_Constant_Names["pourpoints_add_obs"]
    lake_outflow_pourpoints = Internal_Constant_Names["lake_outflow_pourpoints"]
    
    #### Determine Sub subregion without lake
    os.environ.update(
        dict(GRASS_COMPRESS_NULLS="1", GRASS_COMPRESSOR="ZSTD", GRASS_VERBOSE="1")
    )
    PERMANENT = Session()
    PERMANENT.open(
        gisdb=grassdb, location=grass_location, create_opts=""
    )
    N_Basin = 0
    Acc = Initaial_Acc
    print("##############################Loop for suitable ACC ")
    while N_Basin < Min_Num_Domain or N_Basin > Max_Num_Domain:
        grass.run_command(
            "r.watershed",
            elevation=dem,
            flags="sa",
            basin="testbasin",
            drainage="dir_grass_reg",
            accumulation="acc_grass_reg2",
            threshold=Acc,
            overwrite=True,
        )
        N_Basin, temp = generate_stats_list_from_grass_raster(
            grass, mode=1, input_a='testbasin'
        )
        N_Basin = np.unique(N_Basin)
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
    
    if fdr_path == '#':
        mode = 'usingdem'
    else:
        mode = 'usingfdr'
        
    watershed_delineation_without_lake(
        mode=mode,
        input_geo_names=input_geo_names,
        acc_thresold=Acc,
        fdr_path=fdr_path,
        fdr_arcgis=fdr_arcgis,
        fdr_grass=fdr_grass,
        str_r=str_r,
        str_v=str_v,
        acc=acc,
        cat_no_lake=cat_no_lake,
        max_memroy=max_memory,
        grassdb=grassdb,
        grass_location=grass_location,
        qgis_prefix_path=qgis_prefix_path,
        gis_platform='qgis',
    )
    input_geo_names['fdr_arcgis'] = 'fdr_arcgis'
    input_geo_names['fdr_grass'] = 'fdr_grass'
    input_geo_names['str_r'] = 'str_r'
    input_geo_names['str_v'] = 'str_v'
    input_geo_names['acc'] = 'acc'
    input_geo_names['cat_no_lake'] = 'cat_no_lake'
    
    add_lakes_into_existing_watershed_delineation(
        grassdb=grassdb,
        grass_location=grass_location,
        qgis_prefix_path=qgis_prefix_path,
        input_geo_names=input_geo_names,
        path_lakefile_in=path_lakefile_in,
        lake_attributes=lake_attributes,
        threshold_con_lake=CheckLakeArea,
        threshold_non_con_lake=100000,
        sl_connected_lake=sl_connected_lake,
        sl_non_connected_lake=sl_non_connected_lake,
        sl_lakes=sl_lakes,
        nfdr_arcgis=nfdr_arcgis,
        nfdr_grass=nfdr_grass,
        cat_add_lake=cat_add_lake,
        pourpoints_with_lakes=pourpoints_with_lakes,
        cat_use_default_acc=cat_use_default_acc,
        lake_outflow_pourpoints=lake_outflow_pourpoints,
        max_memroy=max_memory,
    )
    
    ####Determin river network for whole watersheds
    PERMANENT = Session()
    PERMANENT.open(
        gisdb=grassdb, location=grass_location, create_opts=""
    )
    grass.run_command(
        "r.stream.extract",
        elevation=dem,
        accumulation=acc,
        threshold=Acc_Thresthold_stream,
        stream_raster=Sub_Reg_str_grass_r,
        stream_vector=Sub_Reg_str_grass_v,
        overwrite=True,
        memory=max_memory,
    )

    #### export outputs
    grass.run_command(
        "r.pack",
        input=nfdr_grass,
        output=os.path.join(Out_Sub_Reg_Folder,sub_reg_nfdr_grass),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input=nfdr_arcgis,
        output=os.path.join(Out_Sub_Reg_Folder,sub_reg_nfdr_arcgis),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input=acc,
        output=os.path.join(Out_Sub_Reg_Folder,sub_reg_acc),
        overwrite=True,
    )
    grass.run_command(
        "r.pack", input=input_geo_names["dem"], output=os.path.join(Out_Sub_Reg_Folder,sub_reg_dem), overwrite=True
    )
    grass.run_command(
        "v.pack",
        input=sub_reg_str_v,
        output=os.path.join(Out_Sub_Reg_Folder,sub_reg_str_v),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input=sub_reg_str_r,
        output=os.path.join(Out_Sub_Reg_Folder,sub_reg_str_r),
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