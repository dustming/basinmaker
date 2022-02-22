from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import pandas as pd
import numpy as np
import tempfile

def Generate_Grid_Poly_From_NetCDF_QGIS(
    NetCDF_Path="#",
    Output_Folder="#",
    Coor_x_NM="lon",
    Coor_y_NM="lat",
    Is_Rotated_Grid=1,
    R_Coor_x_NM="rlon",
    R_Coor_y_NM="rlat",
    SpatialRef="EPSG:4326",
    x_add=-360,
    y_add=0,
    qgis_prefix_path='#',
):

    """Generate Grid polygon from NetCDF file
    Function that used to generate grid polygon from a NetCDF file
    Parameters
    ----------
    NetCDF_Path                       : string
        It is the path of the NetCDF file
    Output_Folder                     : string
        It is the path to a folder to save output polygon shpfiles
    Coor_x_NM                         : string
        It is the variable name for the x coordinates of grids in
        the NetCDF file
    Coor_y_NM                         : string
        It is the variable name for the y coordinates of grids in
        the NetCDF file
    Is_Rotated_Grid                   : Integer
        1 : indicate the grid in NetCDF file is rotated
        -1: indicate the grid in NetCDF file is not rotated
    R_Coor_x_NM                       : string
        It is the variable name for the y coordinates of rotated
        grids in the NetCDF file
    R_Coor_y_NM                       : string
        It is the variable name for the y coordinates of rotated
        grids in the NetCDF file
    SpatialRef                        : string
        It is the coordinates system used in the NetCDF file
    x_add                             : float
        It is offset value for x coodinate
    y_add                             : float
        It is offset value for y coodinate
    Notes
    -------
    Nc_Grids.shp                      : Point shpfile (output)
       It is point in the center of each netCDF Grids
    Gridncply.shp                     : Polygon shpfile (output)
       It is the polygon for each grid in the NetCDF
    Returns:
    -------
       None
    Examples
    -------
    """

    if not os.path.exists(Output_Folder):
        os.makedirs(Output_Folder)
    tempfolder = os.path.join(
        tempfile.gettempdir(), "basinmaker_" + str(np.random.randint(1, 10000 + 1))
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
    from netCDF4 import Dataset

    ncfile = NetCDF_Path
    dsin2 = Dataset(
        ncfile, "r"
    )  # sample structure of in nc file converted from fst

    if Is_Rotated_Grid > 0:
        ncols = len(dsin2.variables[R_Coor_x_NM][:])  ### from 0 to (ncols-1).
        nrows = len(dsin2.variables[R_Coor_y_NM][:])
    else:
        ncols = len(dsin2.variables[Coor_x_NM][:])  ### from 0 to (ncols-1).
        nrows = len(dsin2.variables[Coor_y_NM][:])

    latlonrow = np.full((nrows * ncols, 5), -9999.99999)
    latlonrow = np.full((nrows * ncols, 5), -9999.99999)

    ### Create a point layer, each point will be the nc grids
    cmds = (
        "Point?crs=%s&field=FGID:integer&field=Row:integer&field=Col:integer&field=Gridlon:double&field=Gridlat:double&index=yes"
        % SpatialRef
    )
    Point_Nc_Grid = QgsVectorLayer(cmds, "NC Grid Points", "memory")
    DP_Nc_Point = Point_Nc_Grid.dataProvider()
    Point_Nc_Grid.startEditing()

    ### create polygon layter
    cmds = (
        "Polygon?crs=%s&field=FGID:integer&field=Row:integer&field=Col:integer&field=Gridlon:double&field=Gridlat:double&index=yes"
        % SpatialRef
    )
    Polygon_Nc_Grid = QgsVectorLayer(cmds, "NC Grid polygons", "memory")
    DP_Nc_ply = Polygon_Nc_Grid.dataProvider()
    Polygon_Nc_Grid.startEditing()

    for i in range(0, nrows):
        for j in range(0, ncols):
            k = i * ncols + j

            Point_Fea = QgsFeature()

            if Is_Rotated_Grid < 0:
                latlonrow[k, 0] = k
                latlonrow[k, 1] = i  ### irow
                latlonrow[k, 2] = j  ###col
                latlonrow[k, 3] = dsin2.variables[Coor_x_NM][j] + x_add  ## lon
                latlonrow[k, 4] = dsin2.variables[Coor_y_NM][i]  ## lat

            else:
                latlonrow[k, 0] = k
                latlonrow[k, 1] = i  ### irow
                latlonrow[k, 2] = j  ###col
                latlonrow[k, 3] = dsin2.variables[Coor_x_NM][i, j] + x_add  ## lon
                latlonrow[k, 4] = dsin2.variables[Coor_y_NM][i, j]  ## lat

            #### create point for each grid in Net CDF
            NC_Grid_Point = QgsGeometry.fromPointXY(
                QgsPointXY(latlonrow[k, 3], latlonrow[k, 4])
            )
            Point_Fea.setGeometry(NC_Grid_Point)
            Point_Fea.setAttributes(latlonrow[k, :].tolist())
            DP_Nc_Point.addFeature(Point_Fea)

            ### Create a polygon that the current grid is in the center
            ###  of the polygon

            ### find mid point in row direction
            Polygon_Fea = QgsFeature()
   
            x1 = -1
            x3 = -1 
            x4 = -1
            x2 = -1
            if i != 0 and j != 0:
                if Is_Rotated_Grid > 0:

                    # upper left x1,y1 
                    x1 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j - 1] +
                         dsin2.variables[Coor_x_NM][i-1, j] +
                         dsin2.variables[Coor_x_NM][i-1, j - 1]) / 4 + x_add
                         
                    y1 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j - 1] +
                         dsin2.variables[Coor_y_NM][i-1, j] +
                         dsin2.variables[Coor_y_NM][i-1, j - 1]) / 4
                         
                         
            if i!= 0 and j != ncols -1:
                if Is_Rotated_Grid > 0:

                    x2 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i-1, j] +
                         dsin2.variables[Coor_x_NM][i-1, j + 1] +
                         dsin2.variables[Coor_x_NM][i, j + 1]) / 4 + x_add              
                    y2 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i-1, j] +
                         dsin2.variables[Coor_y_NM][i-1, j + 1] +
                         dsin2.variables[Coor_y_NM][i, j + 1]) / 4 
            
            if i != nrows - 1 and j != ncols -1:
                if Is_Rotated_Grid > 0: 
                    # bot right x3,y3 
                    x3 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j +1] +
                         dsin2.variables[Coor_x_NM][i+1, j + 1] +
                         dsin2.variables[Coor_x_NM][i + 1, j]) / 4 + x_add

                    y3 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j +1] +
                         dsin2.variables[Coor_y_NM][i+1, j + 1] +
                         dsin2.variables[Coor_y_NM][i + 1, j]) / 4
                         
            if i != nrows - 1 and j != 0 :                            
                if Is_Rotated_Grid > 0:
                    x4 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j - 1] +
                         dsin2.variables[Coor_x_NM][i+1, j - 1] +
                         dsin2.variables[Coor_x_NM][i + 1, j]) / 4 + x_add

                    y4 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j - 1] +
                         dsin2.variables[Coor_y_NM][i+1, j - 1] +
                         dsin2.variables[Coor_y_NM][i + 1, j]) / 4
            
            # p4 and p1 and p2    
            if i == 0 and j == 0:
                if Is_Rotated_Grid > 0:
                    # x1 = dsin2.variables[Coor_x_NM][i, j]+ x_add  
                    # y1 = dsin2.variables[Coor_y_NM][i, j]
                    
                    x2_t = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j + 1]) / 2 + x_add              
                    y2_t = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j + 1]) / 2 

                    x2 = x2_t + (x2_t - x3)
                    y2 = y2_t + (y2_t - y3)
                    
                    x4 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i + 1, j]) / 2 + x_add

                    y4 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i + 1, j]) / 2

            if i == 0 and j != 0 and j != ncols -1 :
                if Is_Rotated_Grid > 0:
                    x1 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j - 1]) / 2 + x_add
                         
                    y1 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j - 1]) / 2
                    
                    x2 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j + 1]) / 2 + x_add              
                    y2 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j + 1]) / 2 
                         
            if i == 0 and j ==  ncols -1:
                if Is_Rotated_Grid > 0:
                    x1 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j - 1]) / 2 + x_add
                         
                    y1 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j - 1]) / 2

                    x2 = dsin2.variables[Coor_x_NM][i, j]+ x_add  
                    y2 = dsin2.variables[Coor_y_NM][i, j]                    

                    x3 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i + 1, j]) / 2 + x_add

                    y3 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i + 1, j]) / 2                                                                        
            
            if i == nrows - 1 and j != 0 and j != ncols -1:
                if Is_Rotated_Grid > 0:
                    x3 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j +1]) / 2 + x_add

                    y3 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j +1] ) / 2   
                                                                                        
                    x4 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j - 1]) / 2 + x_add

                    y4 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j - 1]) / 2

            if i == nrows - 1 and j == 0:
                if Is_Rotated_Grid > 0:
            
                    x4 = dsin2.variables[Coor_x_NM][i, j] + x_add
            
                    y4 = dsin2.variables[Coor_y_NM][i, j]
            
                    x1 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i-1, j]) / 2 + x_add
                    y1 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i-1, j]) / 2
            
                    x3 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j +1]) / 2 + x_add
            
                    y3 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j +1]) / 2
                                                                                                                                                              

                         
            if i == nrows - 1 and j == ncols -1:
                if Is_Rotated_Grid > 0:
            
                    x3 = dsin2.variables[Coor_x_NM][i, j] + x_add
            
                    y3 = dsin2.variables[Coor_y_NM][i, j]
            
                    x2 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i-1, j]) / 2 + x_add              
                    y2 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i-1, j]) / 2
            
                    x4 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i, j - 1]) / 2 + x_add
            
                    y4 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i, j - 1]) / 2


            if i != nrows - 1 and i !=0 and j == 0:
                if Is_Rotated_Grid > 0:

                    x1 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i-1, j]) / 2 + x_add
                         
                    y1 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i-1, j]) / 2
                                     
                         
                    x4 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i + 1, j]) / 2 + x_add

                    y4 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i + 1, j]) / 2
                         

            if i != nrows - 1 and i !=0 and j == ncols -1:
                if Is_Rotated_Grid > 0:

                    x2 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i-1, j] ) / 2 + x_add              
                    y2 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i-1, j]) / 2
                                                                                                                                             
                    x3 = (dsin2.variables[Coor_x_NM][i, j] +
                         dsin2.variables[Coor_x_NM][i + 1, j]) / 2 + x_add

                    y3 = (dsin2.variables[Coor_y_NM][i, j] +
                         dsin2.variables[Coor_y_NM][i + 1, j]) / 2                    

            if x1 == -1 or x3 == -1 or x4 == -1 or x2 == -1:
                continue
            Point_1 = QgsPointXY(x1, y1)  ## lower left
            Point_2 = QgsPointXY(x2, y2)
            Point_3 = QgsPointXY(x3, y3)
            Point_4 = QgsPointXY(x4, y4)
            # if i == nrows - 1:
            #     print("#########################################")
            #     print(x1, y1)
            #     print(x2, y2)
            #     print(x3, y3)
            #     print(x4, y4)
            #     print(dsin2.variables[Coor_x_NM][i, j] + x_add,dsin2.variables[Coor_y_NM][i, j])
            #     print("#########################################")
            
            gPolygon = QgsGeometry.fromPolygonXY(
                [[Point_1, Point_2, Point_3, Point_4]]
            )
            Polygon_Fea.setGeometry(gPolygon)
            Polygon_Fea.setAttributes(latlonrow[k, :].tolist())
            DP_Nc_ply.addFeature(Polygon_Fea)

    Point_Nc_Grid.commitChanges()
    Point_Nc_Grid.updateExtents()

    Polygon_Nc_Grid.commitChanges()
    Polygon_Nc_Grid.updateExtents()

    pdlatlonrow = pd.DataFrame(
        latlonrow, columns=["FGID", "Row", "Col", "Gridlon", "Gridlat"]
    )
    pdlatlonrow.to_csv(
        os.path.join(Output_Folder, "Gridcorr.csv"), sep=",", index=False
    )

    QgsVectorFileWriter.writeAsVectorFormat(
        layer=Point_Nc_Grid,
        fileName=os.path.join(Output_Folder, "Nc_Grids.shp"),
        fileEncoding="UTF-8",
        destCRS=QgsCoordinateReferenceSystem(SpatialRef),
        driverName="ESRI Shapefile",
    )
    QgsVectorFileWriter.writeAsVectorFormat(
        layer=Polygon_Nc_Grid,
        fileName=os.path.join(Output_Folder, "Gridncply.shp"),
        fileEncoding="UTF-8",
        destCRS=QgsCoordinateReferenceSystem(SpatialRef),
        driverName="ESRI Shapefile",
    )

def Area_Weighted_Mapping_Between_Two_Polygons_QGIS(
    Target_Ply_Path="#",
    Mapping_Ply_Path="#",
    Col_NM="HRU_ID",
    Output_Folder="#",
    qgis_prefix_path = '#',
):

    """Generate Grid polygon from NetCDF file
    Function that used to generate grid polygon from a NetCDF file
    Parameters
    ----------
    Target_Ply_Path                       : string
        It is the path of one inputs HRU polygon file
    Mapping_Ply_Path                      : string
        It is the path of one inputs grid polygon file
    Output_Folder                     : string
        It is the path to a folder to save output polygon shpfiles
    Notes
    -------
    Overlay_Polygons.shp                 : Polygon shpfile (output)
       It is overlay of two input polygons shpfiles
    GriddedForcings2.txt                 : Text file (output)
       It is the polygon area weighted of each polygon in Mapping_Ply_Path
       to each polygon in Target_Ply_Path
    Returns:
    -------
       None
    Examples
    -------
    """
    
    if not os.path.exists(Output_Folder):
        os.makedirs(Output_Folder)
        
    tempfolder = os.path.join(
        tempfile.gettempdir(), "basinmaker_" + str(np.random.randint(1, 10000 + 1))
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

    Path_finalcat_hru_temp = os.path.join(
        tempfolder,
        str(np.random.randint(1, 10000 + 1)) + "finalcat_freferen.shp",
    )
    Path_finalcat_hru_temp2 = os.path.join(
        tempfolder,
        str(np.random.randint(1, 10000 + 1)) + "finalcat_freferen2.shp",
    )
    Path_finalcat_hru_temp_dissolve = os.path.join(
        tempfolder,
        str(np.random.randint(1, 10000 + 1)) + "finalcat_freferen_dissolve.shp",
    )
    Path_finalcat_hru_temp_dissolve_area = os.path.join(
        Output_Folder, "Overlay_Polygons.shp"
    )
    
    Target_Ply_mem = qgis_vector_fix_geometries(
        processing, context, INPUT=Target_Ply_Path, OUTPUT="memory:"
    )["OUTPUT"]

    Mapping_Ply_mem = qgis_vector_fix_geometries(
        processing, context, INPUT=Mapping_Ply_Path, OUTPUT="memory:"
    )["OUTPUT"]
            
    ### create overlay betweeo two polygon and calcuate area of
    ### each new polygon in the overlay
    qgis_vector_union_two_layers(
        processing=processing,
        context=context,
        INPUT=Target_Ply_mem,
        OVERLAY=Mapping_Ply_mem,
        OVERLAY_FIELDS_PREFIX="Map_",
        OUTPUT=Path_finalcat_hru_temp,
    )["OUTPUT"]

    finalcat_hru_temp2 = processing.run(
        "native:extractbyattribute",
        {
            "INPUT": Path_finalcat_hru_temp,
            "FIELD": "HRU_ID",
            "OPERATOR": 2,
            "VALUE": "0",
            "OUTPUT": "memory:",
        },
    )["OUTPUT"]
    processing.run(
        "native:extractbyattribute",
        {
            "INPUT": finalcat_hru_temp2,
            "FIELD": "Map_FGID",
            "OPERATOR": 2,
            "VALUE": "0",
            "OUTPUT": Path_finalcat_hru_temp2,
        },
    )
        
    processing.run(
        "native:dissolve",
        {
            "INPUT": Path_finalcat_hru_temp2,
            "FIELD": ["HRU_ID", "Map_FGID"],
            "OUTPUT": Path_finalcat_hru_temp_dissolve,
        },
        context=context,
    )
    qgis_vector_field_calculator(
        processing=processing,
        context=context,
        FORMULA="area(transform($geometry, 'EPSG:3161','EPSG:3161'))",
        FIELD_NAME="s_area",
        INPUT=Path_finalcat_hru_temp_dissolve,
        OUTPUT=Path_finalcat_hru_temp_dissolve_area,
        FIELD_PRECISION=3,
    )["OUTPUT"]

    ### calculate the area weight of the mapping polygon to target polygon

    dbf1 = Dbf5(Mapping_Ply_Path[:-3] + "dbf")
    Forcinfo = dbf1.to_dataframe()
    Avafgid = Forcinfo["FGID"].values

    dbf2 = Dbf5(Path_finalcat_hru_temp_dissolve_area[:-3] + "dbf")
    Mapforcing = dbf2.to_dataframe()
    Mapforcing = Mapforcing.loc[Mapforcing[Col_NM] > 0]  ### remove
    Mapforcing = Mapforcing.loc[Mapforcing["Map_FGID"] > 0]  ### remove
    Mapforcing = Mapforcing.loc[Mapforcing["s_area"] > 0.000001]  ### remove
    
    grid_weight_string = create_grid_weight_main(Mapforcing,Forcinfo)
    grid_weight_file_path = os.path.join(Output_Folder, "GriddedForcings2.txt")
    WriteStringToFile(grid_weight_string, grid_weight_file_path, "w")
