from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
import os


def obtain_polygon_boundary(grassdb, grass_location,qgis_prefix_path, lake_name, output):
    
    
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

    exp = "%s = int(%s)" % (lake_name, lake_name)
    grass_raster_r_mapcalc(grass, exp)

    
    grass.run_command(
        "r.to.vect",
        input=lake_name,
        output=lake_name + "plybd",
        type="area",
        overwrite=True,
        flags="v",
    )
    grass.run_command(
        "v.out.ogr",
        input=lake_name + "plybd",
        output=os.path.join(grassdb, lake_name + "plybd.shp"),
        format="ESRI_Shapefile",
        overwrite=True,
    )
    
    PERMANENT.close()
    
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

    # obtain lake boundary lines
    qgis_vector_polygon_stro_lines(
        processing,
        context,
        INPUT=os.path.join(grassdb, lake_name + "plybd.shp"),
        OUTPUT=output,
    )


def reproject_clip_vectors_by_polygon(
    grassdb, grass_location, qgis_prefix_path, mask, path_polygon, ply_name
):

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

    crs_id = qgis_vector_return_crs_id(
        processing, context, mask, Input_Is_Feature_In_Mem=False
    )

    qgis_vector_reproject_layers(
        processing,
        context,
        INPUT=path_polygon,
        TARGET_CRS=crs_id,
        OUTPUT=os.path.join(grassdb, ply_name + "_project.shp"),
    )
    # lake polygon sometime has error in geometry
    try:
        qgis_vector_ectract_by_location(
            processing,
            context,
            INPUT=os.path.join(grassdb, ply_name + "_project.shp"),
            INTERSECT=mask,
            OUTPUT=os.path.join(grassdb, ply_name + ".shp"),
        )
    except:
        print("Need fix lake boundary geometry to speed up")
        qgis_vector_fix_geometries(
            processing,
            context,
            INPUT=os.path.join(grassdb, ply_name + "_project.shp"),
            OUTPUT=os.path.join(grassdb, ply_name + "_fixgeo.shp"),
        )
        qgis_vector_ectract_by_location(
            processing,
            context,
            INPUT=os.path.join(grassdb, ply_name + "_fixgeo.shp"),
            INTERSECT=mask,
            OUTPUT=os.path.join(grassdb, ply_name + ".shp"),
        )
    Qgs.exit()
