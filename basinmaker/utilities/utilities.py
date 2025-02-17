import os
import glob,shutil
from simpledbf import Dbf5
import pandas as pd 
# define internal file names


import os

def get_lake_paths(Path_Connect_Lake_ply, Path_Non_Connect_Lake_ply, Path_Subbasin_Ply):
    """
    Checks for lake shapefiles based on the provided inputs.
    
    Parameters:
      Path_Connect_Lake_ply (str): User-provided full path for the "connect lake" shapefile 
                                   or "#" if not provided.
      Path_Non_Connect_Lake_ply (str): User-provided full path for the "non-connect lake" shapefile 
                                       or "#" if not provided.
      Path_Subbasin_Ply (str): Full path to a subbasin shapefile whose base folder will be searched.
    
    Behavior:
      - If both Path_Connect_Lake_ply and Path_Non_Connect_Lake_ply are "#", the function returns them as-is.
      - Otherwise, it obtains the base folder from Path_Subbasin_Ply and looks for:
            • a shapefile (.shp) whose filename contains both "con" and "lake" (for the connect lake)
            • a shapefile (.shp) whose filename contains both "non" and "lake" (for the non-connect lake)
      - For each, if a matching file is found, its full path is returned; if not, "#" is returned.
    
    Returns:
      tuple: (connect_lake_full_path, non_connect_lake_full_path)
    """
    
    # If both provided lake paths are "#", just return them.
    if Path_Connect_Lake_ply != "#" or Path_Non_Connect_Lake_ply != "#":
        return Path_Connect_Lake_ply, Path_Non_Connect_Lake_ply

    # Otherwise, get the base folder from the subbasin file path.
    base_folder = os.path.dirname(Path_Subbasin_Ply)
    
    # Initialize the return values to "#" (not found)
    found_connect = "#"
    found_non = "#"
    
    # Loop through files in the base folder.
    try:
        for filename in os.listdir(base_folder):
            # Consider only shapefiles
            if filename.lower().endswith('.shp'):
                lower_filename = filename.lower()
                # Look for the connect lake file: must contain both "con" and "lake"
                if ("con" in lower_filename and "lake" in lower_filename) and "non" not in lower_filename:
                    found_connect = os.path.join(base_folder, filename)
                # Look for the non-connect lake file: must contain both "non" and "lake"
                if "non" in lower_filename and "lake" in lower_filename:
                    found_non = os.path.join(base_folder, filename)
    except FileNotFoundError:
        # If the base folder doesn't exist, we simply return "#" for both.
        pass

    return found_connect, found_non

# Example usage:
# If the user did not provide the lake shapefile paths, they might pass "#":
# result = get_lake_paths("#", "#", r"C:\Data\Subbasins\my_subbasin.shp")
# In this case, the function will simply return ("#", "#").
#
# Otherwise, if at least one of the lake paths is not "#", the function will search the folder
# where "my_subbasin.shp" is located for shapefiles that match the criteria.



def copy_files_with_extension(src_folder,tar_folder,ext):

    files = glob.iglob(os.path.join(src_folder, ext))
    for file in files:
        if os.path.isfile(file):
            # Note here if src_folder and tar_folder are the same, the copy command will not work!
            t = shutil.copy2(file, tar_folder)

def adjust_col_dtypes(data):
    for col in data.columns:
        if col in column_dtypes:
            data[col] = data[col].astype(column_dtypes[col])
            if column_dtypes[col] == 'str':
                data[col] = data[col].fillna('-1.2345')
            elif column_dtypes[col] == 'float':
                data[col]  = pd.to_numeric(data[col] , errors='coerce')
                data[col] = data[col].fillna(-1.2345)
            elif column_dtypes[col] == 'int':
                data[col]  = pd.to_numeric(data[col] , errors='coerce')
                data[col] = data[col].fillna(-12345)
            else:
                continue  
         
    return data

column_dtypes = {
    'SubId': 'int',
    'DowSubId': 'int',
    'RivSlope': 'float',
    'RivLength': 'float',
    'BasSlope': 'float',
    'BasAspect': 'float',
    'BasArea': 'float',
    'BkfWidth': 'float',
    'BkfDepth': 'float',
    'Lake_Cat': 'int',
    'HyLakeId': 'int',
    'LakeVol': 'float',
    'LakeDepth': 'float',
    'LakeArea': 'int',
    'Laketype': 'int',
    'Has_POI': 'int',
    'MeanElev': 'int',
    'FloodP_n': 'float',
    'Q_Mean': 'float',
    'Ch_n': 'float',
    'DrainArea': 'float',
    'Strahler': 'int',
    'Seg_ID': 'int',
    'Seg_order': 'int',
    'Max_DEM': 'float',
    'Min_DEM': 'float',
    'DA_Obs': 'float',
    'DA_Diff': 'str',
    'Obs_NM': 'str',
    'SRC_obs': 'str',
    'centroid_x': 'float',
    'centroid_y': 'float',
    'DA_Chn_L': 'float',
    'DA_Slope': 'float',
    'DA_Chn_Slp': 'float',
    'outletLat': 'float',
    'outletLng': 'float',
    'k': 'float',
    'c': 'float',
    'Use_region' : 'int',
    'Vol_res' : 'float',
}

Internal_Constant_Names = {
    "cat_add_lake_old_fdr": "cat_add_lake_old_fdr",
    "cat_add_lake": "cat_add_lake",
    "pourpoints_add_obs": "pourpoints_add_obs",
    "all_lakes": "all_lakes",
    "lake_boundary": "lake_boundary",
    "connect_lake": "connect_lake",
    "nonconnect_lake": "nonconnect_lake",
    "str_connected_lake": "str_connected_lake",
    "str_sl_connected_lake": "str_sl_connected_lake",
    "nfdr_arcgis": "nfdr_arcgis",
    "nfdr_grass": "nfdr_grass",
    "obs": "obs",
    "pourpoints_with_lakes": "pourpoints_with_lakes",
    "lake_inflow_pourpoints": "lake_inflow_pourpoints",
    "lake_outflow_pourpoints": "lake_outflow_pourpoints",
    "catchment_pourpoints_outside_lake": "catchment_pourpoints_outside_lake",
    "cat_use_default_acc": "cat_use_default_acc",
    "cat_ply_info": "cat_ply_info",
    "cat_riv_info": "cat_riv_info",
    "outlet_pt_info": "outlet_pt_info",
    "fdr_grass": "fdr_grass",
    "fdr_arcgis": "fdr_arcgis",
    "str_r": "str_r",
    "str_v": "str_v",
    "acc": "acc",
    "cat_no_lake": "cat_no_lake",
    "sl_connected_lake": "sl_connected_lake",
    "sl_nonconnect_lake": "sl_nonconnect_lake",
    "selected_lakes": "selected_lakes",
    "catchment_without_merging_lakes": "catchment_without_merging_lakes",
    "river_without_merging_lakes": "river_without_merging_lakes",
    "snapped_obs_points": "snapped_obs_points",
    "cat_use_default_acc": "cat_use_default_acc",
    "snapped_obs_points": "snapped_obs_points",
}


####
COLUMN_NAMES_CONSTANT = [
    "SubId",
    "DowSubId",
    "RivSlope",
    "RivLength",
    "BasSlope",
    "BasAspect",
    "BasArea",
    "BkfWidth",
    "BkfDepth",
    "Lake_Cat",
    "HyLakeId",
    "LakeVol",
    "LakeDepth",
    "LakeArea",
    "Laketype",
    "Has_POI",
    "MeanElev",
    "FloodP_n",
    "Q_Mean",
    "Ch_n",
    "DrainArea",
    "Strahler",
    "Seg_ID",
    "Seg_order",
    "Max_DEM",
    "Min_DEM",
    "DA_Obs",
    "DA_error",
    "Obs_NM",
    "SRC_obs",
    "centroid_x",
    "centroid_y",
    "DA_Chn_L",
    "DA_Slope",
    "DA_Chn_Slp",
    "outletLat",
    "outletLng",
    "k",
    "c",
]

COLUMN_NAMES_CONSTANT_CLEAN = [
    "SubId",
    "DowSubId",
    "RivSlope",
    "RivLength",
    "BasSlope",
    "BasAspect",
    "BasArea",
    "BkfWidth",
    "BkfDepth",
    "Lake_Cat",
    "HyLakeId",
    "LakeVol",
    "LakeDepth",
    "LakeArea",
    "Laketype",
    "Has_POI",
    "MeanElev",
    "FloodP_n",
    "Q_Mean",
    "Ch_n",
    "DrainArea",
    "Strahler",
    "Seg_ID",
    "Seg_order",
    "Max_DEM",
    "Min_DEM",
    "DA_Obs",
    "DA_error",
    "Obs_NM",
    "SRC_obs",
    "centroid_x",
    "centroid_y",
    "DA_Chn_L",
    "DA_Slope",
    "DA_Chn_Slp",
    "outletLat",
    "outletLng",
    "Has_Gauge",
    "k",
    "c",
]


COLUMN_NAMES_With_NULL_Values = [
    "Lake_Cat",
    "HyLakeId",
    "LakeVol",
    "LakeDepth",
    "LakeArea",
    "Laketype",
    "DA_Obs",
    "DA_error",
    "Obs_NM",
    "SRC_obs",
]

COLUMN_NAMES_CONSTANT_HRU = [
    "SubId",
    "DowSubId",
    "RivSlope",
    "RivLength",
    "BasSlope",
    "BasAspect",
    "BasArea",
    "BkfWidth",
    "BkfDepth",
    "Lake_Cat",
    "HyLakeId",
    "LakeVol",
    "LakeDepth",
    "LakeArea",
    "Laketype",
    "Has_POI",
    "MeanElev",
    "FloodP_n",
    "Q_Mean",
    "Ch_n",
    "DrainArea",
    "Strahler",
    "Seg_ID",
    "Seg_order",
    "Max_DEM",
    "Min_DEM",
    "DA_Obs",
    "DA_error",
    "Obs_NM",
    "SRC_obs",
    "centroid_x",
    "centroid_y",
    "HRU_IsLake",
    "Landuse_ID",
    "Soil_ID",
    "Veg_ID",
    "O_ID_1",
    "O_ID_2",
    "HRU_Area",
    "HRU_ID",
    "LAND_USE_C",
    "VEG_C",
    "SOIL_PROF",
    "HRU_CenX",
    "HRU_CenY",
    "HRU_S_mean",
    "HRU_A_mean",
    "HRU_E_mean",
    "SHAPE",
    "DA_Chn_L",
    "DA_Slope",
    "DA_Chn_Slp",
    "Has_Gauge" ,
]


COLUMN_TYPES_CONSTANT = [
    "Integer",
    "Integer",
    "Real",
    "Real",
    "Real",
    "Real",
    "Real",
    "Real",
    "Real",
    "Integer",
    "Integer",
    "Real",
    "Real",
    "Real",
    "Integer",
    "Integer",
    "Real",
    "Real",
    "Real",
    "Real",
    "Real",
    "Integer",
    "Integer",
    "Integer",
    "Real",
    "Real",
    "Real",
    "Real",
    "Character",
    "Character",
    "Real",
    "Real",
    "Real",
    "Real",
    "Real",
    "Real",
    "Real",
    "Real",
    "Real",
]

DEFALUT_FLOOD_N = 0.09
min_manning_n = 0.025
max_manning_n = 0.15
min_riv_slope = 0.000001
max_riv_slope = 1
min_bkf_width = 0.1
min_bkf_depth = 0.1
min_riv_lenth = 90

def Dbf_To_Dataframe(file_path):
    """Transfer an input dbf file to dataframe

    Parameters
    ----------
    file_path   : string
    Full path to a shapefile

    Returns:
    -------
    dataframe   : datafame
    a pandas dataframe of attribute table of input shapefile
    """
    tempinfo = Dbf5(file_path[:-3] + "dbf")
    dataframe = tempinfo.to_dataframe().copy()
    return dataframe


def WriteStringToFile(Out_String, File_Path, WriteMethod):
    """Write String to a file

    Function that used to write Out_String to a file located at the File_Path.

    Parameters
    ----------
    Out_String            : string
        The string that will be writed to the file located at File_Path
    File_Path             : string
        Path and filename of file that will be modified or created
    WriteMethod           : {'a','w'}
        If WriteMethod = "w", a new file will be created at the File_Path
        If WriteMethod = "a", the Out_String will be added to exist file

    Notes
    ------
        The file located at the File_Path will be modified or created

    Returns
    -------
        None

    Examples
    --------
    >>> from WriteRavenInputs import WriteStringToFile
    >>> Out_String = 'sometest at line 1\n some test at line 2\n some test at line 3\n'
    >>> File_Path  = 'C:/Path_to_the_Flie_with_file_name'
    >>> WriteStringToFile(Out_String = Out_String,File_Path = File_Path,WriteMethod = 'w')

    """

    if os.path.exists(
        File_Path
    ):  ### if file exist, we can either modify or overwrite it
        with open(File_Path, WriteMethod) as f:
            f.write(Out_String)
    else:  ## create a new file anyway, since file did not exist
        with open(File_Path, "w") as f:
            f.write(Out_String)


def write_grass_and_arcgis_fdr_rules(grassdb):

    Strlist = [
        "1 = 128",
        "2 = 64",
        "3 = 32",
        "4 = 16",
        "5 = 8",
        "6 = 4",
        "7 = 2",
        "8 = 1",
        "* = NULL",
    ]
    Str = "\n".join(Strlist)
    WriteStringToFile(Str, os.path.join(grassdb, "Grass2ArcgisDIR.txt"), "w")

    Strlist = [
        "1 = 8",
        "2 = 7",
        "4 = 6",
        "8 = 5",
        "16 = 4",
        "32 = 3",
        "64 = 2",
        "128 = 1",
        "* = NULL",
    ]
    Str = "\n".join(Strlist)
    WriteStringToFile(Str, os.path.join(grassdb, "Arcgis2GrassDIR.txt"), "w")

    return


def write_grass_reclass_rule_from_table(table, output_path):
    Strlist = []
    for i in range(0, len(table)):
        Strlist.append(str(table[i, 0]) + " = " + str(table[i, 1]))
    Strlist.append("* = NULL")
    Str = "\n".join(Strlist)
    WriteStringToFile(Str, output_path, "w")
