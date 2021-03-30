from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import sqlite3


def define_cat_and_riv_without_merge_lake_cats(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    path_lakefile_in,
    catchment_without_merging_lakes="catchment_without_merging_lakes",
    river_without_merging_lakes="river_without_merging_lakes",
    max_memroy=1024 * 4,
):

    fdr_arcgis = input_geo_names["fdr_arcgis"]
    fdr_grass = input_geo_names["fdr_grass"]
    str_r = input_geo_names["str_r"]
    str_v = input_geo_names["str_v"]
    acc = input_geo_names["acc"]
    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]
    nfdr_grass = input_geo_names["nfdr_grass"]
    sl_non_connected_lake = input_geo_names["sl_non_connected_lake"]
    pourpoints_add_obs = input_geo_names["pourpoints_add_obs"]
    lake_outflow_pourpoints = input_geo_names["lake_outflow_pourpoints"]

    import grass.script as grass
    import grass.script.setup as gsetup
    from grass.pygrass.modules import Module
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.script import array as garray
    from grass.script import core as gcore
    from grass_session import Session

    os.environ.update(
        dict(GRASS_COMPRESS_NULLS="1", GRASS_COMPRESSOR="ZSTD", GRASS_VERBOSE="1")
    )
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")

    grass.run_command(
        "r.to.vect",
        input=pourpoints_add_obs,
        output=pourpoints_add_obs,
        type="point",
        overwrite=True,
        flags="v",
    )

    if path_lakefile_in != "#":
        grass.run_command(
            "r.stream.basins",
            direction=nfdr_grass,
            points=pourpoints_add_obs,
            basins="define_cat_and_riv_1",
            overwrite=True,
            memory=max_memroy,
        )
        # find non connected lake catchments
        lakeid, catid = generate_stats_list_from_grass_raster(
            grass,
            mode=2,
            input_a=lake_outflow_pourpoints,
            input_b="define_cat_and_riv_1",
        )
        lakecat_ids = np.column_stack((lakeid, catid))

        # find all non connected cat ids
        non_connected_lakeids, temp = generate_stats_list_from_grass_raster(
            grass,
            mode=1,
            input_a=sl_non_connected_lake,
        )

        non_cl_lake_catid = lakecat_ids[
            np.isin(lakecat_ids[:, 0], np.array(non_connected_lakeids)), 1
        ]

        # find all catids
        all_catids, temp = generate_stats_list_from_grass_raster(
            grass,
            mode=1,
            input_a="define_cat_and_riv_1",
        )
        all_catids = np.array(all_catids)
        not_non_cl_lake_catids = all_catids[
            np.logical_not(np.isin(all_catids, non_cl_lake_catid))
        ]
        # grass.run_command(
        #     "g.copy", rast=("define_cat_and_riv_1", "non_cl_lake_cat"), overwrite=True
        # )
        non_cl_lake_cat_array = garray.array(mapname="define_cat_and_riv_1")
        if (len(not_non_cl_lake_catids) > 0):
            mask = np.isin(non_cl_lake_cat_array, not_non_cl_lake_catids)
            non_cl_lake_cat_array[mask] = -9999
        temparray = garray.array()
        temparray[:, :] = non_cl_lake_cat_array[:, :]
        temparray.write(mapname="non_cl_lake_cat", overwrite=True)
        grass.run_command("r.null", map="non_cl_lake_cat", setnull=[-9999, 0])
        exp = "%s = int(%s)" % (
            "non_cl_lake_cat",
            "non_cl_lake_cat",
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)

        # grass.run_command(
        #     "r.null",
        #     map="non_cl_lake_cat",
        #     setnull=not_non_cl_lake_catids,
        #     overwrite=True,
        # )

        # define str and cat
        grass.run_command(
            "r.cross",
            input=[str_r, "define_cat_and_riv_1"],
            output="nstr_seg_t",
            flags="z",
            overwrite=True,
        )
        grass.run_command(
            "r.mapcalc", expression="nstr_seg = int(nstr_seg_t) + 1", overwrite=True
        )
        grass.run_command(
            "r.stream.basins",
            direction=nfdr_grass,
            stream="nstr_seg",
            basins="define_cat_and_riv_2",
            overwrite=True,
        )

        grass.run_command(
            "r.cross",
            input=["define_cat_and_riv_2", "non_cl_lake_cat"],
            output=catchment_without_merging_lakes,
            overwrite=True,
        )

        exp = "%s = %s + 1" % (
            catchment_without_merging_lakes,
            catchment_without_merging_lakes,
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)

        exp = "%s = if(isnull(int(%s)),null(),%s)" % (
            river_without_merging_lakes,
            str_r,
            catchment_without_merging_lakes,
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    else:
        grass.run_command(
            "r.stream.basins",
            direction=fdr_grass,
            points=pourpoints_add_obs,
            basins="define_cat_and_riv_1",
            overwrite=True,
            memory=max_memroy,
        )

        grass.run_command(
            "r.cross",
            input=[str_r, "define_cat_and_riv_1"],
            output="nstr_seg_t",
            flags="z",
            overwrite=True,
        )
        grass.run_command(
            "r.mapcalc", expression="nstr_seg = int(nstr_seg_t) + 1", overwrite=True
        )

        grass.run_command(
            "r.stream.basins",
            direction=fdr_grass,
            stream="nstr_seg",
            basins=catchment_without_merging_lakes,
            overwrite=True,
        )

        exp = "%s = if(isnull(int(%s)),null(),%s)" % (
            river_without_merging_lakes,
            str_r,
            catchment_without_merging_lakes,
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    PERMANENT.close()
    return
