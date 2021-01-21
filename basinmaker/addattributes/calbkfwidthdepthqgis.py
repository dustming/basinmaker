from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
import sqlite3
from scipy.optimize import curve_fit
from preprocessing.reprojectandclipvectorbyplyqgis import (
    reproject_clip_vectors_by_polygon,
)


def calculate_bankfull_width_depth_from_polyline(
    grassdb,
    grass_location,
    qgis_prefix_path,
    path_bkfwidthdepth,
    bkfwd_attributes,
    catinfo,
    input_geo_names,
    k_in=-1,
    c_in=-1,
    return_k_c_only=False,
):
    mask = input_geo_names["mask"]

    default_slope = 0.000012345
    min_manning_n = 0.01
    max_manning_n = 0.15
    default_bkf_width = 1.2345
    default_bkf_depth = 1.2345
    default_bkf_q = 1.2345
    k = -1
    c = -1
    if path_bkfwidthdepth != "#":
        reproject_clip_vectors_by_polygon(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            mask=os.path.join(grassdb, mask + ".shp"),
            path_polygon=path_bkfwidthdepth,
            ply_name="bkf_width_depth",
        )

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

    con = sqlite3.connect(
        os.path.join(grassdb, grass_location, "PERMANENT", "sqlite", "sqlite.db")
    )

    if path_bkfwidthdepth != "#":
        grass.run_command(
            "v.import",
            input=os.path.join(grassdb, "bkf_width_depth" + ".shp"),
            output="bk_full_wid_depth",
            overwrite=True,
        )

    if k_in < 0 and c_in < 0:
        ### read catchment
        sqlstat = "SELECT %s,%s,%s,%s FROM %s" % (
            bkfwd_attributes[3],
            bkfwd_attributes[2],
            bkfwd_attributes[1],
            bkfwd_attributes[0],
            "bk_full_wid_depth",
        )
        bkf_width_depth = pd.read_sql_query(sqlstat, con)
        bkf_width_depth = bkf_width_depth.fillna(-9999)

        da_q = bkf_width_depth[[bkfwd_attributes[3], bkfwd_attributes[2]]].values

        if len(da_q) > 3:
            k, c = return_k_and_c_in_q_da_relationship(da_q)
        elif len(da_q) > 0 and len(da_q) <= 3:
            k = -1
            c = -1
            default_bkf_width = np.average(bkf_width_depth[bkfwd_attributes[0]])
            default_bkf_depth = np.average(bkf_width_depth[bkfwd_attributes[1]])
            default_bkf_q = np.average(bkf_width_depth[bkfwd_attributes[3]])
        else:
            k = -1
            c = -1
    else:
        k = k_in
        c = c_in

    if return_k_c_only:
        return k, c

    idx = catinfo.index
    for i in range(0, len(idx)):
        idx_i = idx[i]
        da = catinfo.loc[idx_i, "DA"] / 1000 / 1000  # m2 to km2
        catid = catinfo.loc[idx_i, "SubId"]
        if k > 0:
            q = func_Q_DA(da, k, c)
            catinfo.loc[idx_i, "BkfWidth"] = 7.2 * q ** 0.5
            catinfo.loc[idx_i, "BkfDepth"] = 0.27 * q ** 0.3
            catinfo.loc[idx_i, "Q_Mean"] = q
        else:
            catinfo.loc[idx_i, "BkfWidth"] = default_bkf_width
            catinfo.loc[idx_i, "BkfDepth"] = default_bkf_depth
            catinfo.loc[idx_i, "Q_Mean"] = default_bkf_q

    # adjust channel parameters

    catinfo_riv = catinfo.loc[catinfo["IsLake"] < 2]
    Seg_IDS = catinfo_riv["Seg_ID"].values
    Seg_IDS = np.unique(Seg_IDS)

    for iseg in range(0, len(Seg_IDS)):
        i_seg_id = Seg_IDS[iseg]
        i_seg_info = catinfo_riv[catinfo_riv["Seg_ID"] == i_seg_id]
        max_elve_seg = np.max(i_seg_info["Max_DEM"].values)
        min_elve_seg = np.max(i_seg_info["Min_DEM"].values)
        length_seg = np.sum(i_seg_info["RivLength"].values)
        qmean_seg = np.average(i_seg_info["Q_Mean"].values)
        width_seg = np.average(i_seg_info["BkfWidth"].values)
        depth_Seg = np.average(i_seg_info["BkfDepth"].values)
        slope_seg = (max_elve_seg - min_elve_seg) / length_seg
        if slope_seg < 0.000000001:
            slope_seg = default_slope  #### Needs to update later

        n_seg = calculateChannaln(width_seg, depth_Seg, qmean_seg, slope_seg)

        for i in range(0, len(i_seg_info)):
            subid = i_seg_info["SubId"].values[i]
            max_elve_rch = i_seg_info["Max_DEM"].values[i]
            min_elve_rch = i_seg_info["Min_DEM"].values[i]
            length_rch = i_seg_info["RivLength"].values[i]
            qmean_rch = i_seg_info["Q_Mean"].values[i]
            width_rch = i_seg_info["BkfWidth"].values[i]
            depth_rch = i_seg_info["BkfDepth"].values[i]
            slope_rch = (max_elve_seg - min_elve_seg) / length_rch

            if slope_rch < 0.000000001:
                slope_rch = slope_seg

            n_rch = calculateChannaln(width_rch, depth_rch, qmean_rch, slope_rch)

            if n_rch < min_manning_n or n_rch > max_manning_n:
                if n_seg < min_manning_n or n_seg > max_manning_n:
                    if n_rch < min_manning_n:
                        n_rch = min_manning_n
                    else:
                        n_rch = max_manning_n
                else:
                    n_rch = n_seg

            catinfo.loc[catinfo["SubId"] == subid, "RivSlope"] = slope_rch
            catinfo.loc[catinfo["SubId"] == subid, "Ch_n"] = n_rch

    PERMANENT.close()
    return catinfo
