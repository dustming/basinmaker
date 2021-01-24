from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
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

            if j != 0:
                if Is_Rotated_Grid < 0:  ## left x
                    x1 = latlonrow[k, 3] - 0.5 * (
                        latlonrow[k, 3]
                        - (dsin2.variables[Coor_x_NM][j - 1] + x_add)
                    )
                else:
                    x1 = latlonrow[k, 3] - 0.5 * (
                        latlonrow[k, 3]
                        - (dsin2.variables[Coor_x_NM][i, j - 1] + x_add)
                    )
            else:
                if Is_Rotated_Grid < 0:
                    x1 = latlonrow[k, 3] - 0.5 * (
                        -latlonrow[k, 3]
                        + (dsin2.variables[Coor_x_NM][j + 1] + x_add)
                    )
                else:
                    x1 = latlonrow[k, 3] - 0.5 * (
                        -latlonrow[k, 3]
                        + (dsin2.variables[Coor_x_NM][i, j + 1] + x_add)
                    )

            if j != ncols - 1:  ###  right x
                if Is_Rotated_Grid < 0:
                    x2 = latlonrow[k, 3] + 0.5 * (
                        -latlonrow[k, 3]
                        + (dsin2.variables[Coor_x_NM][j + 1] + x_add)
                    )
                else:
                    x2 = latlonrow[k, 3] + 0.5 * (
                        -latlonrow[k, 3]
                        + (dsin2.variables[Coor_x_NM][i, j + 1] + x_add)
                    )
            else:
                if Is_Rotated_Grid < 0:
                    x2 = latlonrow[k, 3] + 0.5 * (
                        latlonrow[k, 3]
                        - (dsin2.variables[Coor_x_NM][j - 1] + x_add)
                    )
                else:
                    x2 = latlonrow[k, 3] + 0.5 * (
                        latlonrow[k, 3]
                        - (dsin2.variables[Coor_x_NM][i, j - 1] + x_add)
                    )

            if i != nrows - 1:  ## lower y
                if Is_Rotated_Grid < 0:
                    y1 = latlonrow[k, 4] + 0.5 * (
                        -latlonrow[k, 4] + dsin2.variables[Coor_y_NM][i + 1]
                    )
                else:
                    y1 = latlonrow[k, 4] + 0.5 * (
                        -latlonrow[k, 4] + dsin2.variables[Coor_y_NM][i + 1, j]
                    )
            else:
                if Is_Rotated_Grid < 0:
                    y1 = latlonrow[k, 4] + 0.5 * (
                        latlonrow[k, 4] - dsin2.variables[Coor_y_NM][i - 1]
                    )
                else:
                    y1 = latlonrow[k, 4] + 0.5 * (
                        latlonrow[k, 4] - dsin2.variables[Coor_y_NM][i - 1, j]
                    )

            if i != 0:  ## upper y
                if Is_Rotated_Grid < 0:
                    y2 = latlonrow[k, 4] - 0.5 * (
                        latlonrow[k, 4] - dsin2.variables[Coor_y_NM][i - 1]
                    )
                else:
                    y2 = latlonrow[k, 4] - 0.5 * (
                        latlonrow[k, 4] - dsin2.variables[Coor_y_NM][i - 1, j]
                    )
            else:
                if Is_Rotated_Grid < 0:
                    y2 = latlonrow[k, 4] - 0.5 * (
                        -latlonrow[k, 4] + dsin2.variables[Coor_y_NM][i + 1]
                    )
                else:
                    y2 = latlonrow[k, 4] - 0.5 * (
                        -latlonrow[k, 4] + dsin2.variables[Coor_y_NM][i + 1, j]
                    )

            Point_1 = QgsPointXY(x1, y1)  ## lower left
            Point_2 = QgsPointXY(x1, y2)
            Point_3 = QgsPointXY(x2, y2)
            Point_4 = QgsPointXY(x2, y1)

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

    ### create overlay betweeo two polygon and calcuate area of
    ### each new polygon in the overlay
    qgis_vector_union_two_layers(
        processing=processing,
        context=context,
        INPUT=Target_Ply_Path,
        OVERLAY=Mapping_Ply_Path,
        OVERLAY_FIELDS_PREFIX="Map_",
        OUTPUT=Path_finalcat_hru_temp2,
    )["OUTPUT"]

    processing.run(
        "native:extractbyattribute",
        {
            "INPUT": Path_finalcat_hru_temp,
            "FIELD": "HRU_ID",
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
        FORMULA="area(transform($geometry, 'EPSG:4326','EPSG:3573'))",
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

    ####
    hruids = Mapforcing["HRU_ID"].values
    hruids = np.unique(hruids)
    #    Lakeids = np.unique(Lakeids)
    ogridforc = open(os.path.join(Output_Folder, "GriddedForcings2.txt"), "w")
    ogridforc.write(":GridWeights" + "\n")
    ogridforc.write("   #      " + "\n")
    ogridforc.write("   # [# HRUs]" + "\n")
    sNhru = len(hruids)

    ogridforc.write("   :NumberHRUs       " + str(sNhru) + "\n")
    sNcell = (max(Forcinfo["Row"].values) + 1) * (max(Forcinfo["Col"].values) + 1)
    ogridforc.write("   :NumberGridCells  " + str(sNcell) + "\n")
    ogridforc.write("   #            " + "\n")
    ogridforc.write("   # [HRU ID] [Cell #] [w_kl]" + "\n")

    for i in range(len(hruids)):
        hruid = hruids[i]
        cats = Mapforcing.loc[Mapforcing["HRU_ID"] == hruid]
        cats = cats[cats["Map_FGID"].isin(Avafgid)]

        if len(cats) <= 0:
            cats = Mapforcing.loc[Mapforcing["HRU_ID"] == hruid]
            print("Following Grid has to be inluded:.......")
            print(cats["Map_FGID"])
        tarea = sum(cats["s_area"].values)
        fids = cats["Map_FGID"].values
        fids = np.unique(fids)
        sumwt = 0.0
        for j in range(0, len(fids)):
            scat = cats[cats["Map_FGID"] == fids[j]]
            if j < len(fids) - 1:
                sarea = sum(scat["s_area"].values)
                wt = float(sarea) / float(tarea)
                sumwt = sumwt + wt
            else:
                wt = 1 - sumwt

            if len(scat["Map_Row"].values) > 1:  ## should be 1
                print(
                    str(catid)
                    + "error: 1 hru, 1 grid, produce muti polygon need to be merged "
                )
                Strcellid = (
                    str(
                        int(
                            scat["Map_Row"].values[0]
                            * (max(Forcinfo["Col"].values) + 1 + misscol)
                            + scat["Map_Col"].values[0]
                        )
                    )
                    + "      "
                )
            else:
                Strcellid = (
                    str(
                        int(
                            scat["Map_Row"].values
                            * (max(Forcinfo["Col"].values) + 1)
                            + scat["Map_Col"].values
                        )
                    )
                    + "      "
                )

            ogridforc.write(
                "    "
                + str(int(hruid))
                + "     "
                + Strcellid
                + "      "
                + str(wt)
                + "\n"
            )
    #        arcpy.AddMessage(cats)
    ogridforc.write(":EndGridWeights")
    ogridforc.close()
    ########
    # /* example of calcuate grid index
    #           0    1    2    3    4
    #       0    0    1    2    3    4
    #       1    5    6    7    8    9
    #       2    10    11    12    13    14
    #       3    15    16    17    18    19
    ##  we have 4 rows (0-3) and 5 cols (0-4), the index of each cell
    #   should be calaulated by row*(max(colnums)+1) + colnum.
    #   for example row =2, col=0, index = 2*(4+1)+0 = 10
    #   for example row 3, col 3, index = 3*(4+1)+3 = 18