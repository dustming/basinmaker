from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *


def delineate_watershed_no_lake_using_subregion_data(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    subreg_fdr_path,
    subreg_acc_path,
    subreg_str_v_path,
    subreg_str_r_path,
    fdr_arcgis,
    fdr_grass,
    str_r,
    str_v,
    acc,
    cat_no_lake,
    max_memroy,
):
    mask = input_geo_names["mask"]

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

    write_grass_and_arcgis_fdr_rules(grassdb)

    # unpack subregion flow accumulation and flow direction data
    grass_raster_r_unpack(grass, input=subreg_fdr_path, output="dir_Arcgis1")
    grass_raster_r_unpack(grass, input=subreg_acc_path, output="grass_acc1")
    exp = "%s  = int(dir_Arcgis1)" % (fdr_arcgis)
    grass_raster_r_mapcalc(grass, expression=exp)
    # calcuate flow direction in grass format
    grass_raster_r_reclass(
        grass,
        input=fdr_arcgis,
        output=fdr_grass,
        rules=os.path.join(grassdb, "Arcgis2GrassDIR.txt"),
    )
    exp = "%s  = grass_acc1" % (acc)
    grass_raster_r_mapcalc(grass, expression=exp)

    # unpack  river network vector
    grass_raster_v_unpack(grass, input=subreg_str_v_path, output=str_v)
    # unpack  river network raster
    grass_raster_r_unpack(grass, input=subreg_str_r_path, output="str_grass_r1")
    # clip river network with mask
    exp = "%s = int(str_grass_r1)" % (str_r)
    grass_raster_r_mapcalc(grass, expression=exp)

    grass_raster_r_stream_basins(
        grass,
        direction=fdr_grass,
        stream=str_r,
        basins=cat_no_lake,
        memory=max_memroy,
    )

    PERMANENT.close()
    return
