from basinmaker.utilities.utilities import Internal_Constant_Names


def divide_domain_into_sub_regions(
    input_geo_names,
    grassdb,
    grass_location,
    qgis_prefix_path,
    path_lakefile_in,
    lake_attributes,
    path_bkfwidthdepth="#",
    bkfwd_attributes="#",
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
    gis_platform="qgis",
):
    cat_add_lake = Internal_Constant_Names["cat_add_lake"]

    if gis_platform == "qgis":
        assert (
            grassdb != "#"
        ), "grass database folder is needed, when gis_platform = qgis"
        assert (
            grass_location != "#"
        ), "grass location name is needed, when gis_platform = qgis"
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from basinmaker.subreg.generatesubregionqgis import (
            Generatesubdomain,
            generatesubdomainmaskandinfo,
        )

        Generatesubdomain(
            input_geo_names=input_geo_names,
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            path_lakefile_in=path_lakefile_in,
            lake_attributes=lake_attributes,
            Min_Num_Domain=Min_Num_Domain,
            Max_Num_Domain=Max_Num_Domain,
            Initaial_Acc=Initaial_Acc,
            Delta_Acc=Delta_Acc,
            CheckLakeArea=CheckLakeArea,
            fdr_path=fdr_path,
            Acc_Thresthold_stream=Acc_Thresthold_stream,
            max_memory=max_memory,
            Out_Sub_Reg_Folder=Out_Sub_Reg_Folder,
            sub_reg_str_r=sub_reg_str_r,
            sub_reg_str_v=sub_reg_str_v,
            sub_reg_nfdr_grass=sub_reg_nfdr_grass,
            sub_reg_nfdr_arcgis=sub_reg_nfdr_arcgis,
            sub_reg_acc=sub_reg_acc,
            sub_reg_dem=sub_reg_dem,
            cat_add_lake=cat_add_lake,
        )

        input_geo_names["cat_add_lake"] = cat_add_lake
        generatesubdomainmaskandinfo(
            Out_Sub_Reg_Dem_Folder=Out_Sub_Reg_Folder,
            input_geo_names=input_geo_names,
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            path_bkfwidthdepth=path_bkfwidthdepth,
            bkfwd_attributes=bkfwd_attributes,
        )


def combine_sub_region(
    path_sub_region_info,
    sub_region_outputfolder,
    outputfolder,
    is_final_result,
    qgis_prefix_path,
    path_subregion_inlet,
    gis_platform="qgis",
    start_sub_id = 0,
    k = 1,
    c = 1,
):
    if gis_platform == "qgis":
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from basinmaker.subreg.generatesubregionqgis import Combine_Sub_Region_Results

        Combine_Sub_Region_Results(
            Sub_Region_info=path_sub_region_info,
            Sub_Region_OutputFolder=sub_region_outputfolder,
            OutputFolder=outputfolder,
            Is_Final_Result=is_final_result,
            qgis_prefix_path=qgis_prefix_path,
            subregion_inlet=path_subregion_inlet,
            start_sub_id = start_sub_id,
            k = k,
            c = c
        )
