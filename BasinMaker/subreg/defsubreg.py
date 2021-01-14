from utilities.utilities import Internal_Constant_Names


def divide_domain_into_sub_regions(
    input_geo_names,
    grassdb,
    grass_location,
    qgis_prefix_path,
    path_lakefile_in,
    lake_attributes,
    Min_Num_Domain=9,
    Max_Num_Domain=13,
    Initaial_Acc=5000,
    Delta_Acc=1000,
    CheckLakeArea=1,
    fdr_path = '#',
    Acc_Thresthold_stream=500,
    max_memory=2048*3,
    Out_Sub_Reg_Folder="#",
    sub_reg_str_r = 'sub_reg_str_r',
    sub_reg_str_v = 'sub_reg_str_v',
    sub_reg_nfdr_grass = 'sub_reg_nfdr_grass',
    sub_reg_nfdr_arcgis = 'sub_reg_nfdr_arcgis',
    sub_reg_acc = 'sub_reg_acc',
    sub_reg_dem = 'sub_reg_dem',
    gis_platform="qgis",
):
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
        from subreg.generatesubregionqgis import (
            Generatesubdomain
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
            fdr_path = fdr_path,
            Acc_Thresthold_stream=Acc_Thresthold_stream,
            max_memory=max_memory,
            Out_Sub_Reg_Folder=Out_Sub_Reg_Folder,
            sub_reg_str_r = sub_reg_str_r,
            sub_reg_str_v = sub_reg_str_v,
            sub_reg_nfdr_grass = sub_reg_nfdr_grass,
            sub_reg_nfdr_arcgis = sub_reg_nfdr_arcgis,
            sub_reg_acc = sub_reg_acc,
            sub_reg_dem = sub_reg_dem,
        )

