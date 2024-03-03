import pandas as pd
import geopandas
from basinmaker.postprocessing.downloadpd import Download_Routing_Product_For_One_Gauge
import os
import numpy as np
from basinmaker.postprocessing.selectprodpurepy import Select_Routing_product_based_SubId_purepy


def Download_Routing_Product_From_Points_Or_LatLon(product_name, Lat=[-1], Lon=[-1]):

    data = pd.DataFrame(
        {
            'Latitude': Lat,
            'Longitude': Lon,
        }
    )

    data_gpd = geopandas.GeoDataFrame(
        data, geometry=geopandas.points_from_xy(data.Longitude, data.Latitude))
    data_gpd.crs = "EPSG:4326"

    if product_name == 'OLRP':
        version = 'v1-0'
        Drainage_region = geopandas.read_file(
            "https://github.com/dustming/RoutingTool/wiki/Files/OLRRP_drainage_region.geojson")

        data_dr = geopandas.overlay(data_gpd, Drainage_region, how='identity')
        if len(data_dr) >= 1 and data_dr['Region'].values[0] in np.unique(Drainage_region['Region'].values):

            dr = data_dr['Region'].values[0]
            Subid, product_path = Download_Routing_Product_For_One_Gauge(
                gauge_name='#', product_name=product_name, region=dr, subreg='#',version=version)
            product_ply_path = os.path.join(
                product_path, 'catchment_without_merging_lakes_'+version+'.shp')
            product = geopandas.read_file(product_ply_path)
            data_gpd = data_gpd.to_crs(product.crs)
            data_pt = geopandas.overlay(data_gpd, product, how='identity')
            Subid = data_pt['SubId'].values[0]
            print("The needed product locates at:", product_path)
            print("The Subbasin Id of the lat, lon is:", Subid)
            return Subid, product_path

        else:
            print("The point did not overlay with the routing product")
            return -1, -1

    if product_name == 'NALRP':
        version = 'v2-1'
        Drainage_region = geopandas.read_file(
            "https://github.com/dustming/RoutingTool/wiki/Files/NA_drainage_region.geojson")

        data_dr = geopandas.overlay(data_gpd, Drainage_region, how='identity')
        if len(data_dr) >= 1:

            dr = data_dr['Region'].values[0]
            Subid, product_path = Download_Routing_Product_For_One_Gauge(
                gauge_name='#', product_name=product_name, region=dr, subreg='#')
            product_ply_path = os.path.join(
                product_path, 'catchment_without_merging_lakes_'+version+'.shp')
            product = geopandas.read_file(product_ply_path)
            data_gpd = data_gpd.to_crs(product.crs)
            data_pt = geopandas.overlay(data_gpd, product, how='identity')
            Subid = data_pt['SubId'].values[0]
            print("The needed product locates at:", product_path)
            print("The Subbasin Id of the lat, lon is:", Subid)
            return Subid, product_path

        else:
            print("The point did not overlay with the routing product")
            return -1, -1

    return -1, -1


def Extract_Routing_Product(version='v1-0', by='Obs_NM', obs_nm='#', subid=-1, region="#", lat=-1, lon=-1, output_path='#'):
    SubId = -1  # default value
    product_path = '#'
    if by == 'Obs_NM':
        SubId, product_path = Download_Routing_Product_For_One_Gauge(
            gauge_name=obs_nm, product_name='OLLRP', region='#', subreg='#', version=version)
    elif by == 'LatLon':
        data = pd.DataFrame(
            {
                'Latitude': [lat],
                'Longitude': [lon],
            }
        )
        data_gpd = geopandas.GeoDataFrame(
            data, geometry=geopandas.points_from_xy(data.Longitude, data.Latitude))
        data_gpd.crs = "EPSG:4326"

        Drainage_region = geopandas.read_file(
            "https://github.com/dustming/RoutingTool/wiki/Files/OLRRP_drainage_region.geojson")

        data_dr = geopandas.overlay(data_gpd, Drainage_region, how='identity')
        if len(data_dr) >= 1 and data_dr['Region'].values[0] in np.unique(Drainage_region['Region'].values):
            dr = data_dr['Region'].values[0]
            SubId, product_path = Download_Routing_Product_For_One_Gauge(
                gauge_name='#', product_name='OLLRP', region=dr, subreg='#',version=version)
            product_folder_real = find_file(
                product_path, 'catchment_without_merging_lakes_'+version+'.shp')
            product_ply_path = os.path.join(
                product_folder_real, 'catchment_without_merging_lakes_'+version+'.shp')
            product = geopandas.read_file(product_ply_path)
            data_gpd = data_gpd.to_crs(product.crs)
            data_pt = geopandas.overlay(data_gpd, product, how='identity')
            SubId = data_pt['SubId'].values[0]
        else:
            print("The point did not overlay with the routing product")
    elif by == 'SubId':
        temp, product_path = Download_Routing_Product_For_One_Gauge(
            gauge_name="#", product_name='OLLRP', region=region, subreg='#', version=version)
        SubId = subid
    else:
        print("The input by is not supported")

    if SubId == -1 or product_path == '#':
        return

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Search for the file
    foldername = find_file(
        product_path, 'catchment_without_merging_lakes_'+version+'.shp')

    Select_Routing_product_based_SubId_purepy(
        OutputFolder=output_path,
        Routing_Product_Folder=foldername,
        mostdownid=[SubId],
        mostupstreamid=[-1],
    )

    from basinmaker.postprocessing.postprocessingfunctions import (
        combine_catchments_covered_by_the_same_lake_method,
    )

    combine_catchments_covered_by_the_same_lake_method(
        Routing_Product_Folder=output_path,
        qgis_prefix_path='#',
        area_thresthold=3,
        length_thresthold=0,
        gis_platform="purepy",
    )

# Function to search for the target file


def find_file(directory, target_file):
    for foldername, subfolders, filenames in os.walk(directory):
        for filename in filenames:
            if filename == target_file:
                return foldername
