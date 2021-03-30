from basinmaker.func.arcgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *


def delineate_watershed_no_lake_using_fdr(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    fdr_path,
    acc_thresold,
    fdr_arcgis,
    fdr_grass,
    str_r,
    str_v,
    acc,
    cat_no_lake,
    max_memroy,
):
    work_folder = grassdb
    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]

    # read and define arcgis work enviroments
    if not os.path.exists(work_folder):
        os.makedirs(work_folder)
    arcpy.env.workspace = work_folder
     
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")
    cellSize = float(arcpy.GetRasterProperties_management(dem+'.tif', "CELLSIZEX").getOutput(0))
    SptailRef = arcpy.Describe(dem+'.tif').spatialReference
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize    
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    arcpy.env.extent = arcpy.Describe(dem + '.tif').extent
    
    
    outFlowDirection = ExtractByMask(fdr_path, mask + ".shp")
    dirraster = SetNull(Raster(outFlowDirection) < 1, Raster(outFlowDirection))
    dirraster.save(fdr_arcgis + '.tif')

    outFlowAccumulation = FlowAccumulation(dirraster)
    outFlowAccumulation.save(acc + ".tif")

    StreamRaster = SetNull(Raster(outFlowAccumulation) < acc_thresold, Raster(outFlowAccumulation))
    StreamRaster = Con(StreamRaster >= 0, 1, 0)
    StreamRaster.save("str_1.tif")
    
    Strlink = StreamLink(StreamRaster, dirraster)
    Strlink.save(str_r + '.tif')
    
    Catchment = Watershed(dirraster,Strlink)

    Catchment.save(cat_no_lake+'.tif')

    StreamToFeature(Strlink, dirraster, str_v + '.shp',"NO_SIMPLIFY")
    
    # copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"Cat1.prj")
    # StreamToFeature(Strlink, dirraster, "DrainL1","NO_SIMPLIFY")
    # copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"DrainL1.prj")
    # arcpy.RasterToPolygon_conversion("Cat1", "Cattemp.shp", "NO_SIMPLIFY")
    # copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"Cattemp.prj")
    # arcpy.Dissolve_management("Cattemp.shp", "Cat1.shp", ["gridcode"])
    # copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"Cat1.prj")
    # arcpy.AddField_management("Cat1.shp", "COAST", "LONG",10,"","", "", "NULLABLE")
    # arcpy.AddField_management("Cat1.shp", "HYBAS_ID", "LONG",10,"","", "", "NULLABLE")
    # arcpy.AddField_management("Cat1.shp", "NEXT_DOWN", "LONG",10,"","", "", "NULLABLE")
    # fieldList = ["GRID_CODE", "FROM_NODE","TO_NODE"]
    # arcpy.JoinField_management("Cat1.shp", "gridcode", "DrainL1.shp", "GRID_CODE", fieldList)
    # arcpy.CalculateField_management("Cat1.shp", "HYBAS_ID", '!FROM_NODE! * 1', "PYTHON")
    # arcpy.CalculateField_management("Cat1.shp", "NEXT_DOWN", '!TO_NODE! * 1', "PYTHON")

    return
