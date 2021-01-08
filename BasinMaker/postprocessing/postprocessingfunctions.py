def combine_catchments_covered_by_the_same_lake(
    OutputFolder="#",
    Path_final_rivply="#",
    Path_final_riv="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
):

    if gis_platform == "qgis":
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from postprocessing.combine import (
            combine_catchments_covered_by_the_same_lake_qgis,
        )

        combine_catchments_covered_by_the_same_lake_qgis(
            OutputFolder=OutputFolder,
            Path_final_rivply=Path_final_rivply,
            Path_final_riv=Path_final_riv,
            qgis_prefix_path=qgis_prefix_path,
        )
