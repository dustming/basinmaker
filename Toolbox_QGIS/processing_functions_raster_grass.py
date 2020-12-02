import copy

import pandas as pd
from grass.script import array as garray

from processing_functions_attribute_table import Check_If_Lake_Have_Multi_OutLet


####
def grass_raster_setnull(
    grass, raster_nm, null_values, create_new_raster, new_raster_nm="#"
):

    if create_new_raster:
        grass_raster_copy(grass, raster_nm, new_raster_nm)
        grass.run_command("r.null", map=new_raster_nm, setnull=null_values)
    else:
        grass.run_command("r.null", map=raster_nm, setnull=null_values)


#####


def grass_raster_copy(grass, raster_nm_in, raster_nm_new):
    grass.run_command("g.copy", rast=(raster_nm_in, raster_nm_new), overwrite=True)


#####


def Return_Raster_As_Array_With_Db_Path(grassdb, grass_location, raster_mn):
    """Transfer an rater in grass database into np array
    Parameters
    ----------
    grassdb         : string
    Full path to a grass database
    grass_location  : string
    location name in that grass database
    raster_mn       : string
    raster name

    Returns:
    -------
    Array            : array
    np array of the raster.

    """
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")
    Array = copy.deepcopy(garray.array(mapname=raster_mn))
    PERMANENT.close()
    Array[Array <= 0] = -9999
    return Array


###
def Return_Raster_As_Array_With_garray(garray_f, raster_mn):
    """Transfer an rater in grass database into np array
    Parameters
    ----------
    raster_mn       : string
    raster name

    Returns:
    -------
    Array            : array
    np array of the raster.

    """
    Array = copy.deepcopy(garray_f.array(mapname=raster_mn))
    Array[Array <= 0] = -9999
    return Array


def grass_raster_r_in_gdal(grass, raster_path, output_nm, location="#"):
    """import dem to target location
    Parameters
    ----------

    Returns:
    -------

    """
    if location != "#":
        grass.run_command(
            "r.in.gdal",
            input=raster_path,
            output=output_nm,
            overwrite=True,
            location=location,
        )
    else:
        grass.run_command(
            "r.in.gdal", input=raster_path, output=output_nm, overwrite=True
        )


###


def grass_raster_r_stream_extract(
    grass,
    elevation,
    accumulation,
    threshold,
    stream_raster,
    stream_vector,
    direction,
    memory,
):
    """using provided elevation, accumulation and threshold to derived
        stream_raster,and direction
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command(
        "r.stream.extract",
        elevation=elevation,
        accumulation=accumulation,
        threshold=threshold,
        stream_raster=stream_raster,
        direction=direction,
        stream_vector=stream_vector,
        overwrite=True,
        memory=memory,
    )


###


def grass_raster_r_stream_basins(grass, direction, stream, basins, memory, points="#"):
    """generate catchemnt for each river segment in stream raster
    Parameters
    ----------

    Returns:
    -------

    """
    if points == "#":
        grass.run_command(
            "r.stream.basins",
            direction=direction,
            stream=stream,
            basins=basins,
            overwrite=True,
            memory=memory,
        )
    else:
        grass.run_command(
            "r.stream.basins",
            direction=direction,
            points=points,
            basins=basins,
            overwrite=True,
            memory=memory,
        )


def grass_export_array_as_raster(grass, array_tempate, array, raster_name, folder_path):
    array_tempate[:, :] = array[:, :]
    temparray.write(mapname=raster_name, overwrite=True)
    grass.run_command("r.null", map=raster_name, setnull=[-9999, 0])
    grass.run_command(
        "r.out.gdal",
        input=raster_name,
        output=os.path.join(folder_path, raster_name),
        format="GTiff",
        overwrite=True,
        quiet="Ture",
    )


def grass_raster_delineate_watershed_with_point_array(
    grass, garray, pourpoints, output_nm, direction, max_memroy
):
    """generate catchemnt for each river segment in stream raster
    Parameters
    ----------

    Returns:
    -------

    """
    temparray = garray.array()
    temparray[:, :] = -9999

    temparray[:, :] = pourpoints[:, :]
    grass_raster_create_raster_from_array(raster_nm=output_nm + "_pt", array=temparray)
    grass_raster_setnull(
        grass,
        raster_nm=output_nm + "_pt",
        null_values=[-9999, 0],
        create_new_raster=False,
    )
    exp = "%s = int(%s)" % (output_nm + "_pt", output_nm + "_pt")
    grass_raster_r_mapcalc(grass, expression=exp)
    grass_raster_r_to_vect(
        grass,
        input=output_nm + "_pt",
        output=output_nm + "_pt",
        type="point",
        flags="v",
    )

    # using pourpoints generate updated subbasins
    grass_raster_r_stream_basins(
        grass=grass,
        direction="dir_grass",
        stream="#",
        basins=output_nm,
        memory=max_memroy,
        points=output_nm + "_pt",
    )


###


def grass_raster_r_accumulate(grass, direction, accumulation, flags):
    """calculate flow accumulation from flow direction dataset
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command(
        "r.accumulate",
        direction=direction,
        accumulation=accumulation,
        flags=flags,
        overwrite=True,
    )


###


def grass_raster_r_external(grass, input, output):
    """calculate flow accumulation from flow direction dataset
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command("r.external", input=input, output=output)


###


def grass_raster_r_clip(grass, input, output):
    """clip raster with mask in grass env
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command("r.clip", input=input, output=output, overwrite=True)


###


def grass_raster_r_to_vect(grass, input, output, type, flags="#"):
    """convert grass raster to vector
    Parameters
    ----------

    Returns:
    -------

    """
    if flags == "#":
        grass.run_command(
            "r.to.vect", input=input, output=output, type=type, overwrite=True
        )

    else:
        grass.run_command(
            "r.to.vect",
            input=input,
            output=output,
            type=type,
            flags=flags,
            overwrite=True,
        )


###


def grass_raster_r_mask(grass, raster_nm, vector_nm="#"):
    """define grass working mask for current location
    Parameters
    ----------

    Returns:
    -------

    """
    if vector_nm == "#":
        grass.run_command("r.mask", raster=raster_nm, maskcats="*", overwrite=True)
    else:
        grass.run_command("r.mask", vector=vector_nm, maskcats="*", overwrite=True)


###


def grass_raster_g_region(grass, raster_nm, zoom="#"):
    """define grass working region for current location
    Parameters
    ----------

    Returns:
    -------

    """
    if zoom == "#":
        grass.run_command("g.region", raster=raster_nm)
    else:
        grass.run_command("g.region", zoom=zoom)


###


def grass_raster_r_unpack(grass, input, output):
    """unpack pregenerated grass rasters
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command("r.unpack", input=input, output=output, overwrite=True)


###


def grass_raster_v_unpack(grass, input, output):
    """unpack pregenerated grass vectors
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command("v.unpack", input=input, output=output, overwrite=True)


###


def grass_raster_r_reclass(grass, input, output, rules):
    """reclassify grass raster dataset
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command(
        "r.reclass", input=input, output=output, rules=rules, overwrite=True
    )


###


def grass_raster_r_watershed(grass, elevation, drainage, accumulation, flags):
    """generate watershed from dem,output includes flow accumulation
        and flow direction
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command(
        "r.watershed",
        elevation=elevation,
        drainage=drainage,
        accumulation=accumulation,
        flags=flags,
        overwrite=True,
    )


###


def grass_raster_r_mapcalc(grass, expression):
    """grass map calculator
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command("r.mapcalc", expression=expression, overwrite=True)


###


def grass_raster_create_raster_empty_raster(garray, raster_nm):
    """grass create a raster with -9999 in current location
    Parameters
    ----------

    Returns:
    -------

    """
    temparray = garray.array()
    temparray[:, :] = -9999
    temparray.write(mapname=raster_nm, overwrite=True)


def grass_raster_create_raster_from_array(raster_nm, array):
    """grass create a raster with input array
    Parameters
    ----------

    Returns:
    -------

    """
    array.write(mapname=raster_nm, overwrite=True)


###


def grass_raster_v_to_raster(grass, input, output, column, use="attr"):
    """grass create a raster with -9999 in current location
    Parameters
    ----------

    Returns:
    -------

    """
    if use == "attr":
        grass.run_command(
            "v.to.rast",
            input=input,
            output=output,
            use="attr",
            attribute_column=column,
            overwrite=True,
        )
    else:
        grass.run_command(
            "v.to.rast", input=input, output=output, use="cat", overwrite=True
        )


###


def grass_raster_r_water_outlet(
    grass, input_dir_nm, output_watshed_nm, outlet_coordinates
):
    """grass generate watershed based on outlet points
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command(
        "r.water.outlet",
        input=input_dir_nm,
        output=output_watshed_nm,
        coordinates=outlet_coordinates,
        overwrite=True,
    )


###


def grass_raster_r_out_gdal(grass, input_nm, output, format="GTiff"):
    """grass export raster in grass working enviroment to other folder in
        different format.
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command(
        "r.out.gdal", input=input_nm, output=output, format="GTiff", overwrite=True
    )


###


def grass_raster_r_stream_snap(
    grass, input, output, stream_rast, accumulation, radius, memory
):
    """snap observation point to closet river network
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command(
        "r.stream.snap",
        input=input,
        output=output,
        stream_rast=stream_rast,
        accumulation=accumulation,
        radius=radius,
        overwrite=True,
        quiet="Ture",
        memory=memory,
    )


###


def Generate_Routing_structure(
    grass, con, cat="Net_cat", acc="acc_grass", Name="a1", str="#"
):

    Raster_res = grass.core.region()["nsres"]

    ###find outlet point of each cat
    grass.run_command(
        "r.stats.zonal",
        base=cat,
        cover=acc,
        method="max",
        output=Name + "_maxacc",
        overwrite=True,
    )
    ### Find the grid that equal to the max acc, thus this is the outlet grids
    exp = "'%s' =if(%s == int('%s'),%s,null())" % (
        Name + "_OL",
        acc,
        Name + "_maxacc",
        cat,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    ### change outlet points to point vector
    grass.run_command(
        "r.to.vect",
        input=Name + "_OL",
        output=Name + "_OL_v",
        type="point",
        flags="v",
        overwrite=True,
    )
    grass.run_command("v.db.addcolumn", map=Name + "_OL_v", columns="SubId int")
    grass.run_command("v.db.update", map=Name + "_OL_v", column="SubId", qcol="cat")
    grass.run_command("v.what.rast", map=Name + "_OL_v", raster=acc, column="OL_acc")
    grass.run_command(
        "v.buffer",
        input=Name + "_OL_v",
        output=Name + "_OL_v_bf",
        distance=1.5 * Raster_res,
        flags="t",
        overwrite=True,
    )
    ###
    ###find inlet  point of each cat
    exp = "'%s' =if(isnull(%s),null(),%s)" % (Name + "_acc_riv", str, acc)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    grass.run_command("r.null", map=Name + "_acc_riv", setnull=[-9999, 0])
    grass.run_command(
        "r.stats.zonal",
        base=str,
        cover=Name + "_acc_riv",
        method="min",
        output=Name + "_minacc",
        overwrite=True,
    )
    ### Find the grid that equal to the max acc, thus this is the outlet grids
    exp = "'%s' =if(%s == int('%s'),%s,null())" % (
        Name + "_IL",
        Name + "_acc_riv",
        Name + "_minacc",
        str,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    ### change outlet points to point vector
    grass.run_command(
        "r.to.vect",
        input=Name + "_IL",
        output=Name + "_IL_v",
        type="point",
        flags="v",
        overwrite=True,
    )
    grass.run_command("v.db.addcolumn", map=Name + "_IL_v", columns="IL_SubId int")
    #    grass.run_command('v.db.addcolumn', map=  Name+'_IL_v', columns = "IL_acc int")
    grass.run_command("v.db.update", map=Name + "_IL_v", column="IL_SubId", qcol="cat")
    grass.run_command("v.what.rast", map=Name + "_IL_v", raster=acc, column="IL_acc")
    #    grass.run_command('v.db.update', map=  Name+'_IL_v', column = "IL_acc",qcol = 'IL_acc_db - 1')

    #### add max acc inlet point acc value to Name+'_OL_v_bf'
    grass.run_command(
        "v.vect.stats",
        point=Name + "_IL_v",
        areas=Name + "_OL_v_bf",
        method="max_cat",
        points_column="IL_acc",
        count_column="IL_Pt_CT",
        stats_column="max_cat",
    )
    grass.run_command(
        "v.vect.stats",
        point=Name + "_IL_v",
        areas=Name + "_OL_v_bf",
        method="min_cat",
        points_column="IL_acc",
        count_column="IL_Pt_CT2",
        stats_column="min_cat",
    )
    grass.run_command(
        "v.vect.stats",
        point=Name + "_IL_v",
        areas=Name + "_OL_v_bf",
        method="maximum",
        points_column="IL_acc",
        count_column="IL_Pt_CT3",
        stats_column="ILmaxacc",
    )
    grass.run_command(
        "v.vect.stats",
        point=Name + "_IL_v",
        areas=Name + "_OL_v_bf",
        method="minimum",
        points_column="IL_acc",
        count_column="IL_Pt_CT4",
        stats_column="ILminacc",
    )

    #### add outlet's next inlet subid into Name+'_OL_v_bf'
    grass.run_command(
        "v.db.join",
        map=Name + "_OL_v_bf",
        column="max_cat",
        other_table=Name + "_IL_v",
        other_column="cat",
        subset_columns="IL_SubId",
        overwrite=True,
    )
    grass.run_command("v.db.addcolumn", map=Name + "_OL_v_bf", columns="ILSubIdmax int")
    grass.run_command(
        "v.db.addcolumn", map=Name + "_OL_v_bf", columns="ILAccstrmin int"
    )

    grass.run_command(
        "v.db.update", map=Name + "_OL_v_bf", column="ILSubIdmax", qcol="IL_SubId"
    )
    grass.run_command(
        "v.db.update", map=Name + "_OL_v_bf", column="ILAccstrmin", qcol="ILmaxacc"
    )
    ### add min next inlet subid to Name+'_OL_v_bf'
    grass.run_command(
        "v.db.join",
        map=Name + "_OL_v_bf",
        column="min_cat",
        other_table=Name + "_IL_v",
        other_column="cat",
        subset_columns="IL_SubId",
        overwrite=True,
    )
    grass.run_command("v.db.addcolumn", map=Name + "_OL_v_bf", columns="ILSubIdmin int")
    grass.run_command(
        "v.db.addcolumn", map=Name + "_OL_v_bf", columns="ILAccstrmax int"
    )

    grass.run_command(
        "v.db.update", map=Name + "_OL_v_bf", column="ILSubIdmin", qcol="IL_SubId"
    )
    grass.run_command(
        "v.db.update", map=Name + "_OL_v_bf", column="ILAccstrmin", qcol="ILminacc"
    )

    grass.run_command("v.db.addcolumn", map=Name + "_OL_v_bf", columns="DSubId_str int")
    grass.run_command("v.db.addcolumn", map=Name + "_OL_v_bf", columns="DSubAccstr int")
    grass.run_command(
        "v.db.update", map=Name + "_OL_v_bf", column="DSubId_str", qcol="ILSubIdmax"
    )
    grass.run_command(
        "v.db.update", map=Name + "_OL_v_bf", column="DSubAccstr", qcol="ILmaxacc"
    )

    ### in case cat drainage to itself, drainage to another close subbasin in inlet
    exp = " '%s' = '%s' " % ("SubId", "ILSubIdmax")
    grass.run_command(
        "v.db.update",
        map=Name + "_OL_v_bf",
        column="DSubId_str",
        qcol="ILSubIdmin",
        where="SubId = ILSubIdmax",
    )
    grass.run_command(
        "v.db.update",
        map=Name + "_OL_v_bf",
        column="DSubAccstr",
        qcol="ILminacc",
        where="SubId = ILSubIdmax",
    )

    ### add downsubid_Str to point vector
    grass.run_command(
        "v.db.join",
        map=Name + "_OL_v",
        column="SubId",
        other_table=Name + "_OL_v_bf",
        other_column="SubId",
        subset_columns=[
            "DSubId_str",
            "IL_Pt_CT",
            "ILSubIdmin",
            "ILSubIdmax",
            "DSubAccstr",
        ],
        overwrite=True,
    )

    ### find downn subid via cat, the reason needs to do this, is that for some cat
    ### without connected river, above approach can not find downstream cat id for
    ### these catchment

    ###add max acc in of neighbor grids as each growed grids value
    exp = "'%s'=if(isnull('%s'),null(),1)" % (Name + "_OL1", Name + "_OL")
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    grass.run_command(
        "r.grow",
        input=Name + "_OL1",
        output=Name + "_OL1_G",
        radius=1.5,
        overwrite=True,
    )
    grass.run_command(
        "r.clump", input=Name + "_OL1_G", output=Name + "_OL1_G_Clu", overwrite=True
    )

    ###find inlets of each subbasin
    grass.run_command(
        "r.stats.zonal",
        base=Name + "_OL1_G_Clu",
        cover=acc,
        method="max",
        output=Name + "_OL1_G_Clu_maxacc",
        overwrite=True,
    )
    exp = "'%s'=if(%s == int('%s'),%s,null())" % (
        Name + "_IL",
        acc,
        Name + "_OL1_G_Clu_maxacc",
        cat,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    grass.run_command(
        "r.stats.zonal",
        base=Name + "_OL1_G_Clu",
        cover=Name + "_IL",
        method="max",
        output=Name + "_OL1_G_Clu_IL_SubId",
        overwrite=True,
    )

    grass.run_command(
        "v.what.rast",
        map=Name + "_OL_v",
        raster=Name + "_OL1_G_Clu_IL_SubId",
        column="DSubId_cat",
    )
    grass.run_command(
        "v.what.rast", map=Name + "_OL_v", raster=Name + "_maxacc", column="MaxAcc_cat"
    )
    #    grass.run_command('v.db.join', map=  Name+'_OL_v',column = 'SubId', other_table = Name+'_OL_v_bf',other_column ='SubId', subset_columns = ['DSubId_str','OL_acc','IL_Pt_CT','max_cat','IL_SubId_min'],overwrite = True)

    #### create the final downsubid use downsubid from str if it is exist
    grass.run_command("v.db.addcolumn", map=Name + "_OL_v", columns="DowSubId int")
    grass.run_command(
        "v.db.update",
        map=Name + "_OL_v",
        column="DowSubId",
        where="IL_Pt_CT > 0",
        qcol="DSubId_str",
    )
    grass.run_command(
        "v.db.update",
        map=Name + "_OL_v",
        column="DowSubId",
        where="IL_Pt_CT <= 0",
        qcol="DSubId_cat",
    )
    ### if it move to itself using dowsubid from cat
    grass.run_command(
        "v.db.update",
        map=Name + "_OL_v",
        column="DowSubId",
        where="SubId = DSubId_str",
        qcol="DSubId_cat",
    )
    ### if it move to a cat with smaller acc, use down subid from cat
    grass.run_command(
        "v.db.update",
        map=Name + "_OL_v",
        column="DowSubId",
        where="OL_acc > DSubAccstr",
        qcol="DSubId_cat",
    )
    sqlstat = "SELECT SubId, DowSubId, MaxAcc_cat FROM %s" % (Name + "_OL_v")
    Routing_info = pd.read_sql_query(sqlstat, con)
    return Routing_info


def DefineConnected_Non_Connected_Lakes(
    self,
    grass,
    con,
    garray,
    Routing_info,
    str_r="str_grass_r",
    Lake_r="alllake",
    Remove_Lake_Multi_OL=True,
):

    #### overlay str and lakes
    exp = "Connect_Lake_IDs = if(isnull('%s'),null(),%s)" % (str_r, Lake_r)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    ### obtain connected lake ids
    Connect_Lake_Ids, temp = Generate_stats_list_from_grass_raster(
        grass, mode=1, input_a="Connect_Lake_IDs"
    )

    #### create non connected lake raster
    grass.run_command("g.copy", rast=(Lake_r, "Nonconnect_Lake"), overwrite=True)
    grass.run_command(
        "r.null", map="Nonconnect_Lake", setnull=Connect_Lake_Ids, overwrite=True
    )
    #### create potential connected lake raster
    exp = "Connect_Lake = if(isnull('%s'),%s,null())" % ("Nonconnect_Lake", Lake_r)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    if Remove_Lake_Multi_OL == False:
        return

    ##### obtain lake id and coresponding overlaied str id
    CL_Id, Str_Id = Generate_stats_list_from_grass_raster(
        grass, mode=2, input_a="Connect_Lake", input_b=str_r
    )

    Lakes_WIth_Multi_Outlet, Remove_Str = Check_If_Lake_Have_Multi_OutLet(
        CL_Id, Str_Id, Routing_info
    )
    ### no lake has multi outlet

    if len(Lakes_WIth_Multi_Outlet) > 0:
        #        grass.run_command('r.null',map='Connect_Lake', setnull = Lakes_WIth_Multi_Outlet,overwrite = True)
        print("Following lakes have identified outlet :       ")
        print("     ", Lakes_WIth_Multi_Outlet)
    if len(Remove_Str) > 0:
        #        grass.run_command('r.null',map='str_grass_r', setnull = Remove_Str,overwrite = True)
        #        grass.run_command('r.null',map='cat1', setnull = Remove_Str,overwrite = True)
        print(
            "Following stream have been identified to make each lake has one lake outlet:       "
        )
        print("     ", Remove_Str)
    return Remove_Str


def Generate_stats_list_from_grass_raster(
    grass, mode=1, input_a="P_Connect_Lake", input_b="str_grass_r", method="max"
):

    list_a = []
    list_b = []
    if mode == 1:
        p = grass.pipe_command("r.category", map=input_a, separator=" ", quiet=True)
    elif mode == 2:
        p = grass.pipe_command(
            "r.stats", input=[input_a, input_b], flags="n", quiet=True
        )
    else:
        print("error mode opton did not exist")

    for line in p.stdout:
        line_str = line.decode("utf8").rstrip("\r\n").split(" ")
        list_a.append(int(line_str[0]))
        if mode != 1:
            list_b.append(int(line_str[1]))
    p.wait()

    return list_a, list_b
