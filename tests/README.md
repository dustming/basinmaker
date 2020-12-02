###test GenerateRavenInput need provide path to Hydat.sqlite3###

pytest  test_RoutingTool_PostProcessing_Locate_subid_needsbyuser.py

pytest  test_RoutingTool_PostProcessing_Select_Routing_product_based_SubId.py

pytest  test_RoutingTool_PostProcessing_Customize_Routing_Topology.py

pytest  test_RoutingTool_PostProcessing_SelectLakes.py

pytest  test_RoutingTool_PostProcessing_Define_Final_Catchment.py

pytest -q -s  --HYDAT_Path "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Database/Obs/Hydat.sqlite3"  test_RoutingTool_PostProcessing_GenerateRavenInput.py

