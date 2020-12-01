import os
import shutil

import pandas as pd
import pytest
from simpledbf import Dbf5

from processing_functions_attribute_table import Evaluate_Two_Dataframes
from ToolboxClass import LRRT
from utilities import Dbf_To_Dataframe


def test_Define_Final_Catchment():
    """test function that will:
    Generate the final lake river routing structure by merging subbasin
    polygons that are covered by the same lake.
    The input are the catchment polygons and river segements
    before merging for lakes. The input files can be output of
    any of following functions:
    SelectLakes, Select_Routing_product_based_SubId,
    Customize_Routing_Topology,RoutingNetworkTopologyUpdateToolset_riv
    The result is the final catchment polygon that ready to be used for
    hydrological modeling
    """
    ###Floder where store the inputs for tests function
    Routing_Product_Folder = "./testdata/Simplified_By_DA/"

    ###Folder where store the expected resuts
    Expect_Result_Folder = os.path.join("./testdata", "Simplified_By_DA")
    ###Folder where the output will be generated
    Output_Folder = os.path.join("./testdata", "testout10")

    ###The pathes for all inputs
    Path_final_riv_ply = os.path.join(
        Routing_Product_Folder, "finalriv_info_ply.shp"
    )  ### River polyline
    Path_final_riv = os.path.join(
        Routing_Product_Folder, "finalriv_info.shp"
    )  ### Catchment polygons

    ###Generate test resuts
    RTtool = LRRT()
    RTtool.Define_Final_Catchment(
        Path_final_rivply=Path_final_riv_ply,
        Path_final_riv=Path_final_riv,
        OutputFolder=Output_Folder,
    )

    """Evaluate attribute table of two polygons   
    """
    ### transfer expected  product into pandas dataframe
    Expect_Finalcat_info = Dbf_To_Dataframe(
        os.path.join(Expect_Result_Folder, "finalcat_info.shp")
    ).sort_values(by=["SubId"])

    ### transfer resulted  product into pandas dataframe
    Result_Finalcat_info = Dbf_To_Dataframe(
        os.path.join(Output_Folder, "finalcat_info.shp")
    ).sort_values(by=["SubId"])

    assert Evaluate_Two_Dataframes(
        Expect_Finalcat_info, Result_Finalcat_info, Check_Col_NM="SubId"
    )

    shutil.rmtree(Output_Folder)
