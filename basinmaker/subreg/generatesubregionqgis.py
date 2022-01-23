from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
from basinmaker.delineationnolake.watdelineationwithoutlake import (
    watershed_delineation_without_lake,
)
from basinmaker.addlakeandobs.addlakesqgis import (
    add_lakes_into_existing_watershed_delineation,
)
from basinmaker.addattributes.calbkfwidthdepthqgis import (
    calculate_bankfull_width_depth_from_polyline,
)
import tempfile
import sqlite3


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
    fdr_path="#",
    Acc_Thresthold_stream=500,
    max_memory=2048 * 3,
    Out_Sub_Reg_Folder="#",
    sub_reg_str_r="sub_reg_str_r",
    sub_reg_str_v="sub_reg_str_v",
    sub_reg_nfdr_grass="sub_reg_nfdr_grass",
    sub_reg_nfdr_arcgis="sub_reg_nfdr_arcgis",
    sub_reg_acc="sub_reg_acc",
    sub_reg_dem="sub_reg_dem",
    cat_add_lake="cat_add_lake",
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

    # required inputs
    dem = input_geo_names["dem"]
    mask = input_geo_names["mask"]

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
    catchment_without_merging_lakes = Internal_Constant_Names[
        "catchment_without_merging_lakes"
    ]
    river_without_merging_lakes = Internal_Constant_Names["river_without_merging_lakes"]
    cat_use_default_acc = Internal_Constant_Names["cat_use_default_acc"]
    pourpoints_with_lakes = Internal_Constant_Names["pourpoints_with_lakes"]
    pourpoints_add_obs = Internal_Constant_Names["pourpoints_add_obs"]
    lake_outflow_pourpoints = Internal_Constant_Names["lake_outflow_pourpoints"]
    all_lakes = input_geo_names["all_lakes"]
    lake_boundary = input_geo_names["lake_boundary"]

    #### Determine Sub subregion without lake
    os.environ.update(
        dict(GRASS_COMPRESS_NULLS="1", GRASS_COMPRESSOR="ZSTD", GRASS_VERBOSE="1")
    )
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")
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
            grass, mode=1, input_a="testbasin"
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

    if fdr_path == "#":
        mode = "usingdem"
    else:
        mode = "usingfdr"

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
        gis_platform="qgis",
    )
    input_geo_names["fdr_arcgis"] = "fdr_arcgis"
    input_geo_names["fdr_grass"] = "fdr_grass"
    input_geo_names["str_r"] = "str_r"
    input_geo_names["str_v"] = "str_v"
    input_geo_names["acc"] = "acc"
    input_geo_names["cat_no_lake"] = "cat_no_lake"

    add_lakes_into_existing_watershed_delineation(
        grassdb=grassdb,
        grass_location=grass_location,
        qgis_prefix_path=qgis_prefix_path,
        input_geo_names=input_geo_names,
        path_lakefile_in=path_lakefile_in,
        lake_attributes=lake_attributes,
        threshold_con_lake=CheckLakeArea,
        threshold_non_con_lake=100000,
        remove_lake_inlets = True,
        only_included_lake_at_river_interction = True,
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
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")
    grass.run_command(
        "r.stream.extract",
        elevation=dem,
        accumulation=acc,
        threshold=Acc_Thresthold_stream,
        stream_raster=sub_reg_str_r,
        stream_vector=sub_reg_str_v,
        overwrite=True,
        memory=max_memory,
    )

    #### export outputs
    grass.run_command(
        "r.pack",
        input=nfdr_grass,
        output=os.path.join(Out_Sub_Reg_Folder, sub_reg_nfdr_grass + ".pack"),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input=nfdr_arcgis,
        output=os.path.join(Out_Sub_Reg_Folder, sub_reg_nfdr_arcgis + ".pack"),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input=fdr_grass,
        output=os.path.join(Out_Sub_Reg_Folder, 'sub_reg_fdr_grass' + ".pack"),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input=fdr_arcgis,
        output=os.path.join(Out_Sub_Reg_Folder, 'sub_reg_fdr_arcgis' + ".pack"),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input=acc,
        output=os.path.join(Out_Sub_Reg_Folder, sub_reg_acc + ".pack"),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input=input_geo_names["dem"],
        output=os.path.join(Out_Sub_Reg_Folder, sub_reg_dem + ".pack"),
        overwrite=True,
    )
    grass.run_command(
        "v.pack",
        input=sub_reg_str_v,
        output=os.path.join(Out_Sub_Reg_Folder, sub_reg_str_v + ".pack"),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input=sub_reg_str_r,
        output=os.path.join(Out_Sub_Reg_Folder, sub_reg_str_r + ".pack"),
        overwrite=True,
    )

    grass.run_command(
        "r.pack",
        input=all_lakes,
        output=os.path.join(Out_Sub_Reg_Folder, all_lakes + ".pack"),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input=lake_boundary,
        output=os.path.join(Out_Sub_Reg_Folder, lake_boundary + ".pack"),
        overwrite=True,
    )

    PERMANENT.close()
    return


def generatesubdomainmaskandinfo(
    Out_Sub_Reg_Dem_Folder,
    input_geo_names,
    grassdb,
    grass_location,
    qgis_prefix_path,
    path_bkfwidthdepth,
    bkfwd_attributes,
):
    ###
    dem = input_geo_names["dem"]
    cat_add_lake = input_geo_names["cat_add_lake"]
    ndir_Arcgis = input_geo_names["nfdr_arcgis"]
    acc_grass = input_geo_names["acc"]
    str_r = input_geo_names["str_r"]
    outlet_pt_info = "outlet_pt_info"

    maximum_obs_id = 80000
    tempfolder = os.path.join(
        tempfile.gettempdir(),
        "basinmaker_subreg" + str(np.random.randint(1, 10000 + 1)),
    )
    if not os.path.exists(tempfolder):
        os.makedirs(tempfolder)

    k = -1
    c = -1
    if path_bkfwidthdepth != "#":
        k, c = calculate_bankfull_width_depth_from_polyline(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            path_bkfwidthdepth=path_bkfwidthdepth,
            bkfwd_attributes=bkfwd_attributes,
            catinfo=[],
            input_geo_names=input_geo_names,
            k_in=-1,
            c_in=-1,
            return_k_c_only=True,
        )

    #### generate subbregion outlet points and subregion info table
    QgsApplication.setPrefixPath(qgis_prefix_path, True)
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
    con = sqlite3.connect(
        os.path.join(grassdb, grass_location, "PERMANENT", "sqlite", "sqlite.db")
    )

    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")
    grass.run_command("r.mask", raster=dem, maskcats="*", overwrite=True)
    grass.run_command("r.null", map=cat_add_lake, setnull=-9999)

    exp = "%s = if(isnull(%s),null(),%s)" % (
        "river_r",
        str_r,
        cat_add_lake,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    routing_temp = generate_routing_info_of_catchments(
        grass,
        con,
        cat=cat_add_lake,
        acc=acc_grass,
        Name="Final",
        str="river_r",
        garray=garray,
    )

    grass.run_command("g.copy", vector=("Final_OL_v", outlet_pt_info), overwrite=True)
    grass.run_command(
        "g.copy", vector=("Final_IL_v_c", "sub_reg_inlet"), overwrite=True
    )
    Paths_Finalcat_ply = []

    sqlstat = (
        "SELECT SubId, DowSubId,ILSubIdmax,ILSubIdmin,MaxAcc_cat,ILpt_ID FROM %s"
        % (outlet_pt_info)
    )
    outletinfo = pd.read_sql_query(sqlstat, con)
    outletinfo = outletinfo.fillna(-1)
    outletinfo = outletinfo.loc[outletinfo["SubId"] > 0]

    sqlstat = "SELECT ILpt_ID,SubId_I FROM %s" % ("Final_IL_v_c")
    inletinfo = pd.read_sql_query(sqlstat, con)
    inletinfo = inletinfo.fillna(-1)

    # update watershed bankfull k and c first
    outletinfo["k"] = k
    outletinfo["c"] = c

    subregin_info = pd.DataFrame(
        np.full(len(outletinfo), -9999), columns=["Sub_Reg_ID"]
    )
    subregin_info["Dow_Sub_Reg_Id"] = -9999
    subregin_info["ProjectNM"] = -9999
    subregin_info["Nun_Grids"] = -9999
    subregin_info["Ply_Name"] = -9999
    subregin_info["Max_ACC"] = -9999
    subregin_info["ILpt_ID"] = -9999

    for i in range(0, len(outletinfo)):

        basinid = int(outletinfo["SubId"].values[i])

        grass.run_command("r.mask", raster=dem, maskcats="*", overwrite=True)
        exp = "%s = if(%s == %s,%s,null())" % (
            "dem_reg_" + str(basinid),
            cat_add_lake,
            str(basinid),
            dem,
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)
        ####define mask
        grass.run_command(
            "r.mask", raster="dem_reg_" + str(basinid), maskcats="*", overwrite=True
        )
        grass.run_command(
            "r.out.gdal",
            input="MASK",
            output=os.path.join(tempfolder, "Mask1.tif"),
            format="GTiff",
            overwrite=True,
        )
        processing.run(
            "gdal:polygonize",
            {
                "INPUT": os.path.join(tempfolder, "Mask1.tif"),
                "BAND": 1,
                "FIELD": "DN",
                "EIGHT_CONNECTEDNESS": False,
                "EXTRA": "",
                "OUTPUT": os.path.join(
                    tempfolder, "HyMask_region_" + str(basinid) + ".shp"
                ),
            },
        )
        processing.run(
            "gdal:dissolve",
            {
                "INPUT": os.path.join(
                    tempfolder, "HyMask_region_" + str(basinid) + ".shp"
                ),
                "FIELD": "DN",
                "OUTPUT": os.path.join(
                    tempfolder, "HyMask_region_f" + str(basinid) + ".shp"
                ),
            },
        )
        processing.run(
            "gdal:dissolve",
            {
                "INPUT": os.path.join(
                    tempfolder, "HyMask_region_" + str(basinid) + ".shp"
                ),
                "FIELD": "DN",
                "OUTPUT": os.path.join(
                    Out_Sub_Reg_Dem_Folder,
                    "HyMask_region_"
                    + str(int(basinid + maximum_obs_id))
                    + "_nobuffer.shp",
                ),
            },
        )

        Paths_Finalcat_ply.append(
            os.path.join(
                Out_Sub_Reg_Dem_Folder,
                "HyMask_region_" + str(int(basinid + maximum_obs_id)) + "_nobuffer.shp",
            )
        )

        processing.run(
            "native:buffer",
            {
                "INPUT": os.path.join(
                    tempfolder, "HyMask_region_f" + str(basinid) + ".shp"
                ),
                "DISTANCE": 0.005,
                "SEGMENTS": 5,
                "END_CAP_STYLE": 0,
                "JOIN_STYLE": 0,
                "MITER_LIMIT": 2,
                "DISSOLVE": True,
                "OUTPUT": os.path.join(
                    Out_Sub_Reg_Dem_Folder,
                    "HyMask_region_" + str(int(basinid + maximum_obs_id)) + ".shp",
                ),
            },
        )

    qgis_vector_merge_vector_layers(
        processing,
        context,
        INPUT_Layer_List=Paths_Finalcat_ply,
        OUTPUT=os.path.join(Out_Sub_Reg_Dem_Folder, "subregion_ply.shp"),
    )

    grass.run_command("r.mask", raster=dem, maskcats="*", overwrite=True)

    problem_subid = []
    for i in range(0, len(outletinfo)):
        basinid = int(outletinfo["SubId"].values[i])
        dowsubreginid = int(outletinfo["DowSubId"].values[i])
        ILpt_ID = outletinfo["ILpt_ID"].values[i]
        subregin_info.loc[i, "ILpt_ID"] = ILpt_ID
        if len(inletinfo[inletinfo["ILpt_ID"] == ILpt_ID]["SubId_I"]) > 0:
            downsubid_inlet = inletinfo[inletinfo["ILpt_ID"] == ILpt_ID][
                "SubId_I"
            ].values[0]
            if dowsubreginid != downsubid_inlet:
                problem_subid.append(basinid)

        catacc = int(outletinfo["MaxAcc_cat"].values[i])

        subregin_info.loc[i, "ProjectNM"] = (
            "sub_reg" + "_" + str(int(basinid + maximum_obs_id))
        )
        subregin_info.loc[i, "Ply_Name"] = (
            "HyMask_region_" + str(int(basinid + maximum_obs_id)) + ".shp"
        )
        subregin_info.loc[i, "Max_ACC"] = catacc

        if basinid == dowsubreginid:
            subregin_info.loc[i, "Dow_Sub_Reg_Id"] = int(-1 + maximum_obs_id)
        else:
            subregin_info.loc[i, "Dow_Sub_Reg_Id"] = int(dowsubreginid + maximum_obs_id)
        subregin_info.loc[i, "Sub_Reg_ID"] = int(basinid + maximum_obs_id)

    subregin_info["k"] = k
    subregin_info["c"] = c

    ### remove subregion do not contribute to the outlet
    ## find watershed outlet subregion
    #    subregin_info  = subregin_info.loc[subregin_info['Dow_Sub_Reg_Id'] == self.maximum_obs_id-1]
    subregin_info = subregin_info.sort_values(by="Max_ACC", ascending=False)
    outlet_reg_id = subregin_info["Sub_Reg_ID"].values[0]
    routing_info = (
        subregin_info[["Sub_Reg_ID", "Dow_Sub_Reg_Id"]].astype("float").values
    )
    needed_sub_reg_ids = defcat(routing_info, outlet_reg_id)

    mask = subregin_info["Sub_Reg_ID"].isin(needed_sub_reg_ids)
    subregin_info = subregin_info.loc[mask, :]
    #        subregin_info.drop(subregin_info.index[del_row_mask]) ###
    subregin_info.to_csv(
        os.path.join(Out_Sub_Reg_Dem_Folder, "Sub_reg_info.csv"),
        index=None,
        header=True,
    )

    subregin_info.to_csv(
        os.path.join(Out_Sub_Reg_Dem_Folder, "Sub_reg_info.csv"),
        index=None,
        header=True,
    )

    grass.run_command("v.db.addcolumn", map=outlet_pt_info, columns="reg_subid int")
    grass.run_command("v.db.addcolumn", map=outlet_pt_info, columns="reg_dowid int")

    grass.run_command("v.db.addcolumn", map=outlet_pt_info, columns="sub_reg_id int")

    grass.run_command(
        "v.db.update",
        map=outlet_pt_info,
        column="reg_subid",
        qcol="SubId + " + str(maximum_obs_id),
    )

    grass.run_command(
        "v.db.update",
        map=outlet_pt_info,
        column="sub_reg_id",
        qcol="SubId + " + str(maximum_obs_id),
    )

    grass.run_command(
        "v.db.update",
        map=outlet_pt_info,
        column="reg_dowid",
        qcol="DowSubId + " + str(maximum_obs_id),
    )

    grass.run_command("v.db.addcolumn", map="sub_reg_inlet", columns="sub_reg_id int")

    grass.run_command(
        "v.db.update",
        map="sub_reg_inlet",
        column="sub_reg_id",
        qcol="SubId_I + " + str(maximum_obs_id),
    )

    grass.run_command(
        "v.out.ogr",
        input=outlet_pt_info,
        output=os.path.join(Out_Sub_Reg_Dem_Folder, outlet_pt_info + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
    )

    grass.run_command(
        "v.out.ogr",
        input="sub_reg_inlet",
        output=os.path.join(Out_Sub_Reg_Dem_Folder, "sub_reg_inlet" + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
    )

    grass.run_command(
        "v.pack",
        input=outlet_pt_info,
        output=os.path.join(Out_Sub_Reg_Dem_Folder, "Sub_Reg_Outlet_v" + ".pack"),
        overwrite=True,
    )
    grass.run_command(
        "r.pack",
        input="Final_OL",
        output=os.path.join(Out_Sub_Reg_Dem_Folder, "Sub_Reg_Outlet_r" + ".pack"),
        overwrite=True,
    )

    print("following subregion's inlet needs to be checked ")
    print(problem_subid)

    return


############################################################################


def Combine_Sub_Region_Results(
    Sub_Region_info,
    Sub_Region_OutputFolder,
    OutputFolder,
    Is_Final_Result,
    qgis_prefix_path,
    subregion_inlet,
    start_sub_id = 0,
    k = 1,
    c = 1,
):
    """Combine subregion watershed delineation results

    It is a function that will combine watershed delineation results
    in different subregions. This function will assgin new subbasin
    id to each polygon in the combined result and update the
    stream orders in the combined result

    Parameters
    ----------
    Sub_Region_info                   : string
        It is the path to a csv file that contains the subregion
        information, such as subregion id etc. It is the output of
        Generatesubdomainmaskandinfo
    Sub_Region_OutputFolder           : string
        It is the path to the subregion output folder.
    OutputFolder                      : string
        It is the path to a folder to save outputs
    Is_Final_Result                   : bool
       Indicate the function is called to combine subbasin polygon of
       final delineation results or subbasin polygons before merging
       for lakes.

    Notes
    -------
    This function has no return values, instead will generate following
    files.
    when Is_Final_Result is true
    os.path.join(OutputFolder,'finalcat_info.shp')
    os.path.join(OutputFolder,'finalcat_info_riv.shp')
    when Is_Final_Result is False
    os.path.join(OutputFolder,'finalrivply_info.shp')
    os.path.join(OutputFolder,'finalriv_info_riv.shp')
    Returns:
    -------
    None

    Examples
    -------

    """
    maximum_obs_id = 80000
    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)

    tempfolder = os.path.join(
        tempfile.gettempdir(),
        "basinmaker_comsubreg" + str(np.random.randint(1, 10000 + 1)),
    )
    if not os.path.exists(tempfolder):
        os.makedirs(tempfolder)

    QgsApplication.setPrefixPath(qgis_prefix_path, True)
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

    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)

    Paths_Finalcat_ply = []
    Paths_Finalcat_line = []
    Paths_Finalriv_ply = []
    Paths_Finalriv_line = []
    Paths_Con_Lake_ply = []
    Paths_None_Con_Lake_ply = []
    Paths_obs_point = []

    Path_Outlet_Down_point = subregion_inlet
    ### add new attribte
    ### find outlet subregion id
    outlet_subregion_id = Sub_Region_info.loc[
        Sub_Region_info["Dow_Sub_Reg_Id"] == maximum_obs_id - 1
    ]["Sub_Reg_ID"].values[0]
    routing_info = (
        Sub_Region_info[["Sub_Reg_ID", "Dow_Sub_Reg_Id"]].astype("float").values
    )
    Subregion_to_outlet = defcat(routing_info, outlet_subregion_id)

    ###remove subregions do not drainge to outlet subregion
    Sub_Region_info = Sub_Region_info.loc[
        Sub_Region_info["Sub_Reg_ID"].isin(Subregion_to_outlet)
    ].copy()
    routing_info = (
        Sub_Region_info[["Sub_Reg_ID", "Dow_Sub_Reg_Id"]].astype("float").values
    )

    Sub_Region_info["N_Up_SubRegion"] = np.nan
    Sub_Region_info["Outlet_SubId"] = np.nan
    subid_strat_iregion = 1
    seg_id_strat_iregion = 1
    for i in range(0, len(Sub_Region_info)):
        isubregion = Sub_Region_info["Sub_Reg_ID"].values[i]
        ProjectNM = Sub_Region_info["ProjectNM"].values[i]
        SubFolder = os.path.join(Sub_Region_OutputFolder, ProjectNM)
    
        ### define path of the output file in this sub region
        Path_Finalcat_ply = os.path.join(SubFolder, "finalcat_info.shp")
        Path_Finalcat_line = os.path.join(SubFolder, "finalcat_info_riv.shp")
        Path_Finalriv_ply = os.path.join(
            SubFolder, "catchment_without_merging_lakes.shp"
        )
        Path_Finalriv_line = os.path.join(SubFolder, "river_without_merging_lakes.shp")
        Path_Con_Lake_ply = os.path.join(SubFolder, "sl_connected_lake.shp")
        Path_None_Con_Lake_ply = os.path.join(SubFolder, "sl_non_connected_lake.shp")
        Path_obs_point = os.path.join(SubFolder, "obs_gauges.shp")
    
    
    
        ### For each subregion, add new subid to each polygon files,
        ### and append result file in the merge list
        if Is_Final_Result == True:
        ### product do not exist
            if os.path.exists(Path_Finalcat_ply) != 1:
                continue
            SubID_info = (
                Dbf_To_Dataframe(Path_Finalcat_ply)
                .drop_duplicates(subset=["SubId"], keep="first")[
                    ["SubId", "DowSubId", "Seg_ID"]
                ]
                .copy()
            )
            SubID_info = SubID_info.reset_index()
            SubID_info["nSubId"] = SubID_info.index + subid_strat_iregion + start_sub_id
            SubID_info["nSeg_ID"] = SubID_info["Seg_ID"] + seg_id_strat_iregion+ start_sub_id
    
            layer_cat = QgsVectorLayer(Path_Finalcat_ply, "")
            Add_New_SubId_To_Subregion_shpfile(
                processing,
                context,
                layer_cat,
                OutputPath=os.path.join(
                    tempfolder,
                    "finalcat_info_Region_" + str(isubregion) + "addatrri.shp",
                ),
                Region_ID=isubregion,
                SubID_info=SubID_info,
            )
            Paths_Finalcat_ply.append(
                os.path.join(
                    tempfolder,
                    "finalcat_info_Region_" + str(isubregion) + "addatrri.shp",
                )
            )
            del layer_cat
    
            layer_cat = QgsVectorLayer(Path_Finalcat_line, "")
            Add_New_SubId_To_Subregion_shpfile(
                processing,
                context,
                layer_cat,
                OutputPath=os.path.join(
                    tempfolder,
                    "finalcat_info_riv_Region_" + str(isubregion) + "addatrri.shp",
                ),
                Region_ID=isubregion,
                SubID_info=SubID_info,
            )
            Paths_Finalcat_line.append(
                os.path.join(
                    tempfolder,
                    "finalcat_info_riv_Region_" + str(isubregion) + "addatrri.shp",
                )
            )
            del layer_cat
    
        else:
            if os.path.exists(Path_Finalriv_ply) != 1:
                continue
            SubID_info = (
                Dbf_To_Dataframe(Path_Finalriv_ply)
                .drop_duplicates(subset=["SubId"], keep="first")[
                    ["SubId", "DowSubId", "Seg_ID"]
                ]
                .copy()
            )
            SubID_info = SubID_info.reset_index()
            SubID_info["nSubId"] = SubID_info.index + subid_strat_iregion + start_sub_id
            SubID_info["nSeg_ID"] = SubID_info["Seg_ID"] + seg_id_strat_iregion + start_sub_id
            layer_cat = QgsVectorLayer(Path_Finalriv_ply, "")
            Add_New_SubId_To_Subregion_shpfile(
                processing,
                context,
                layer_cat,
                OutputPath=os.path.join(
                    tempfolder,
                    "finalriv_info_ply_Region_" + str(isubregion) + "addatrri.shp",
                ),
                Region_ID=isubregion,
                SubID_info=SubID_info,
            )
            Paths_Finalriv_ply.append(
                os.path.join(
                    tempfolder,
                    "finalriv_info_ply_Region_" + str(isubregion) + "addatrri.shp",
                )
            )
            del layer_cat
    
            layer_cat = QgsVectorLayer(Path_Finalriv_line, "")
            Add_New_SubId_To_Subregion_shpfile(
                processing,
                context,
                layer_cat,
                OutputPath=os.path.join(
                    tempfolder,
                    "finalriv_info_Region_" + str(isubregion) + "addatrri.shp",
                ),
                Region_ID=isubregion,
                SubID_info=SubID_info,
            )
            Paths_Finalriv_line.append(
                os.path.join(
                    tempfolder,
                    "finalriv_info_Region_" + str(isubregion) + "addatrri.shp",
                )
            )
            del layer_cat
    
        if os.path.exists(Path_Con_Lake_ply) == 1 and os.stat(Path_Con_Lake_ply).st_size > 100:
            Paths_Con_Lake_ply.append(Path_Con_Lake_ply)
        if os.path.exists(Path_None_Con_Lake_ply) == 1 and os.stat(Path_None_Con_Lake_ply).st_size > 100:
            Paths_None_Con_Lake_ply.append(Path_None_Con_Lake_ply)
        if os.path.exists(Path_obs_point) == 1 and os.stat(Path_obs_point).st_size > 100:
            print(os.stat(Path_obs_point).st_size)
            print(isubregion)
            Paths_obs_point.append(Path_obs_point)
    
        print(
            "Subregion ID is ",
            isubregion,
            "    the start new subid is    ",
            subid_strat_iregion,
            " The end of subid is ",
            max(SubID_info["nSubId"]),
        )
    
        subid_strat_iregion = max(SubID_info["nSubId"]) + 10 - start_sub_id
        seg_id_strat_iregion = max(SubID_info["nSeg_ID"]) + 10 - start_sub_id
    

    # merge connected lake polygons
    if len(Paths_Con_Lake_ply) > 0 and not os.path.exists(os.path.join(OutputFolder, "sl_connected_lake.shp")):
        qgis_vector_merge_vector_layers(
            processing,
            context,
            INPUT_Layer_List=Paths_Con_Lake_ply,
            OUTPUT=os.path.join(OutputFolder, "sl_connected_lake.shp"),
        )
    
    # merge non connected lake polygon
    if len(Paths_None_Con_Lake_ply) > 0 and not os.path.exists(os.path.join(OutputFolder, "sl_non_connected_lake.shp")):
        qgis_vector_merge_vector_layers(
            processing,
            context,
            INPUT_Layer_List=Paths_None_Con_Lake_ply,
            OUTPUT=os.path.join(OutputFolder, "sl_non_connected_lake.shp"),
        )
    
    # merge observation points
    if len(Paths_obs_point) > 0 and not os.path.exists(os.path.join(OutputFolder, "obs_gauges.shp")):
        qgis_vector_merge_vector_layers(
            processing,
            context,
            INPUT_Layer_List=Paths_obs_point,
            OUTPUT=os.path.join(OutputFolder, "obs_gauges.shp"),
        )

    # merge catchment polygon and polyline layers, and update their attirbutes
    if Is_Final_Result == 1:
        #### Obtain downstream # id:
        qgis_vector_merge_vector_layers(
            processing,
            context,
            INPUT_Layer_List=Paths_Finalcat_ply,
            OUTPUT=os.path.join(tempfolder, "finalcat_info.shp"),
        )
        qgis_vector_merge_vector_layers(
            processing,
            context,
            INPUT_Layer_List=Paths_Finalcat_line,
            OUTPUT=os.path.join(tempfolder, "finalcat_info_riv.shp"),
        )
        processing.run(
            "qgis:joinattributesbylocation",
            {
                "INPUT": Path_Outlet_Down_point,
                "JOIN": os.path.join(tempfolder, "finalcat_info.shp"),
                "PREDICATE": [5],
                "JOIN_FIELDS": [],
                "METHOD": 1,
                "DISCARD_NONMATCHING": True,
                "PREFIX": "",
                "OUTPUT": os.path.join(tempfolder, "Down_Sub_ID.shp"),
            },
            context=context,
        )

        AllCatinfo = (
            Dbf_To_Dataframe(os.path.join(tempfolder, "finalcat_info.shp"))
            .drop_duplicates("SubId", keep="first")
            .copy()
        )

        AllCatinfo['SubId'] = AllCatinfo['SubId'].astype('int32')
        AllCatinfo['DowSubId'] = AllCatinfo['DowSubId'].astype('int32')  
        
        DownCatinfo = (
            Dbf_To_Dataframe(os.path.join(tempfolder, "Down_Sub_ID.shp"))
            .drop_duplicates("ILpt_ID", keep="first")
            .copy()
        )

        AllCatinfo, Sub_Region_info = Connect_SubRegion_Update_DownSubId(
            AllCatinfo, DownCatinfo, Sub_Region_info
        )
        AllCatinfo = Update_DA_Strahler_For_Combined_Result(AllCatinfo, Sub_Region_info,k,c)
        
        AllCatinfo['Use_region'] = 0
        Obs_names = AllCatinfo['Obs_NM'].copy(deep=True)
        Obs_names = Obs_names.fillna('nan')
        Obs_names = Obs_names.astype('string')
        mask = Obs_names == 'nan'
        AllCatinfo.loc[mask,'Has_POI'] = 0
        
        Copy_Pddataframe_to_shpfile(
            os.path.join(tempfolder, "finalcat_info.shp"),
            AllCatinfo,
            link_col_nm_shp="SubId",
            link_col_nm_df="SubId",
            UpdateColNM=["DowSubId", "DrainArea", "Strahler","DA_error","Has_POI","Use_region"],
        )
        Copy_Pddataframe_to_shpfile(
            os.path.join(tempfolder, "finalcat_info_riv.shp"),
            AllCatinfo,
            link_col_nm_shp="SubId",
            link_col_nm_df="SubId",
            UpdateColNM=["DowSubId", "DrainArea", "Strahler","DA_error","Has_POI","Use_region"],
        )
        COLUMN_NAMES_CONSTANT_Local = COLUMN_NAMES_CONSTANT_CLEAN 
        COLUMN_NAMES_CONSTANT_Local.append("Region_ID")
        COLUMN_NAMES_CONSTANT_Local.append("Use_region")
        processing.run(
            "native:dissolve",
            {
                "INPUT": os.path.join(tempfolder, "finalcat_info.shp"),
                "FIELD": ["SubId"],
                "OUTPUT": os.path.join(OutputFolder, "finalcat_info.shp"),
            },
            context=context,
        )

        Clean_Attribute_Name(
            os.path.join(OutputFolder, "finalcat_info.shp"), COLUMN_NAMES_CONSTANT_Local
        )

        processing.run(
            "native:dissolve",
            {
                "INPUT": os.path.join(tempfolder, "finalcat_info_riv.shp"),
                "FIELD": ["SubId"],
                "OUTPUT": os.path.join(OutputFolder, "finalcat_info_riv.shp"),
            },
            context=context,
        )
        Clean_Attribute_Name(
            os.path.join(OutputFolder, "finalcat_info_riv.shp"), COLUMN_NAMES_CONSTANT_Local
        )


    else:
        qgis_vector_merge_vector_layers(
            processing,
            context,
            INPUT_Layer_List=Paths_Finalriv_ply,
            OUTPUT=os.path.join(tempfolder, "finalriv_info_ply.shp"),
        )
        qgis_vector_merge_vector_layers(
            processing,
            context,
            INPUT_Layer_List=Paths_Finalriv_line,
            OUTPUT=os.path.join(tempfolder, "finalriv_info.shp"),
        )

        processing.run(
            "qgis:joinattributesbylocation",
            {
                "INPUT": Path_Outlet_Down_point,
                "JOIN": os.path.join(tempfolder, "finalriv_info_ply.shp"),
                "PREDICATE": [5],
                "JOIN_FIELDS": [],
                "METHOD": 1,
                "DISCARD_NONMATCHING": True,
                "PREFIX": "",
                "OUTPUT": os.path.join(tempfolder, "Down_Sub_ID.shp"),
            },
            context=context,
        )

        AllCatinfo = (
            Dbf_To_Dataframe(os.path.join(tempfolder, "finalriv_info_ply.shp"))
            .drop_duplicates("SubId", keep="first")
            .copy()
        )
        DownCatinfo = (
            Dbf_To_Dataframe(os.path.join(tempfolder, "Down_Sub_ID.shp"))
            .drop_duplicates("ILpt_ID", keep="first")
            .copy()
        )
        AllCatinfo['SubId'] = AllCatinfo['SubId'].astype('int32')
        AllCatinfo['DowSubId'] = AllCatinfo['DowSubId'].astype('int32')  
        
        AllCatinfo, Sub_Region_info = Connect_SubRegion_Update_DownSubId(
            AllCatinfo, DownCatinfo, Sub_Region_info
        )
        AllCatinfo['Use_region'] = 0
        AllCatinfo = Update_DA_Strahler_For_Combined_Result(AllCatinfo, Sub_Region_info,k,c)
        Obs_names = AllCatinfo['Obs_NM'].copy(deep=True)
        Obs_names = Obs_names.fillna('nan')
        Obs_names = Obs_names.astype('string')
        mask = Obs_names == 'nan'
        AllCatinfo.loc[mask,'Has_POI'] = 0
        
        DA_obs = AllCatinfo['DA_Obs'].copy(deep=True)
        DA_obs = DA_obs.fillna(0)
        DA_obs = DA_obs.astype('int32')
        mask = DA_obs <= 0
        AllCatinfo.loc[mask,'DA_error'] = 0

                
        Copy_Pddataframe_to_shpfile(
            os.path.join(tempfolder, "finalriv_info_ply.shp"),
            AllCatinfo,
            link_col_nm_shp="SubId",
            link_col_nm_df="SubId",
            UpdateColNM=["DowSubId", "DrainArea", "Strahler","DA_error","Has_POI","Use_region"],
        )
        Copy_Pddataframe_to_shpfile(
            os.path.join(tempfolder, "finalriv_info.shp"),
            AllCatinfo,
            link_col_nm_shp="SubId",
            link_col_nm_df="SubId",
            UpdateColNM=["DowSubId", "DrainArea", "Strahler","DA_error","Has_POI","Use_region"],
        )

        processing.run(
            "native:dissolve",
            {
                "INPUT": os.path.join(tempfolder, "finalriv_info_ply.shp"),
                "FIELD": ["SubId"],
                "OUTPUT": os.path.join(OutputFolder, "catchment_without_merging_lakes.shp"),
            },
            context=context,
        )
        COLUMN_NAMES_CONSTANT_Local = COLUMN_NAMES_CONSTANT 
        COLUMN_NAMES_CONSTANT_Local.append("Region_ID")
        COLUMN_NAMES_CONSTANT_Local.append("Use_region")
        Clean_Attribute_Name(
            os.path.join(OutputFolder, "catchment_without_merging_lakes.shp"), COLUMN_NAMES_CONSTANT_Local
        )

        processing.run(
            "native:dissolve",
            {
                "INPUT": os.path.join(tempfolder, "finalriv_info.shp"),
                "FIELD": ["SubId"],
                "OUTPUT": os.path.join(OutputFolder, "river_without_merging_lakes.shp"),
            },
            context=context,
        )

        Clean_Attribute_Name(
            os.path.join(OutputFolder, "river_without_merging_lakes.shp"), COLUMN_NAMES_CONSTANT_Local
        )
