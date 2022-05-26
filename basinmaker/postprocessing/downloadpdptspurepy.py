import pandas as pd
import geopandas
from basinmaker.postprocessing.downloadpd import Download_Routing_Product_For_One_Gauge
import os
import numpy as np 

def Download_Routing_Product_From_Points_Or_LatLon(product_name,Lat = [-1],Lon = [-1]):
    
    data = pd.DataFrame(
    {
     'Latitude': Lat,
     'Longitude': Lon,
     }
     )
    
    data_gpd = geopandas.GeoDataFrame(data, geometry=geopandas.points_from_xy(data.Longitude, data.Latitude))
    data_gpd.crs = "EPSG:4326"
    
    if product_name == 'OLRP':
        version = 'v1-0'
        Drainage_region = geopandas.read_file("https://github.com/dustming/RoutingTool/wiki/Files/OLRRP_drainage_region.geojson")
    
        data_dr = geopandas.overlay(data_gpd, Drainage_region, how='identity')
        if len(data_dr) >= 1 and data_dr['Region'].values[0] in np.unique(Drainage_region['Region'].values):

            dr = data_dr['Region'].values[0]
            Subid,product_path = Download_Routing_Product_For_One_Gauge(gauge_name='#',product_name = product_name,region=dr,subreg = '#')
            product_ply_path = os.path.join(product_path,'catchment_without_merging_lakes_'+version+'.shp')
            product = geopandas.read_file(product_ply_path)
            data_gpd = data_gpd.to_crs(product.crs)
            data_pt = geopandas.overlay(data_gpd, product, how='identity')
            Subid = data_pt['SubId'].values[0]
            print("The needed product locates at:",product_path)
            print("The Subbasin Id of the lat, lon is:",Subid)
            return Subid,product_path 

        else:
            print("The point did not overlay with the routing product")
            return -1,-1

    if product_name == 'NALRP':
        version = 'v2-1'
        Drainage_region = geopandas.read_file("https://github.com/dustming/RoutingTool/wiki/Files/NA_drainage_region.geojson")
    
        data_dr = geopandas.overlay(data_gpd, Drainage_region, how='identity')
        if len(data_dr) >= 1:

            dr = data_dr['Region'].values[0]
            Subid,product_path = Download_Routing_Product_For_One_Gauge(gauge_name='#',product_name = product_name,region=dr,subreg = '#')
            product_ply_path = os.path.join(product_path,'catchment_without_merging_lakes_'+version+'.shp')
            product = geopandas.read_file(product_ply_path)
            data_gpd = data_gpd.to_crs(product.crs)
            data_pt = geopandas.overlay(data_gpd, product, how='identity')
            Subid = data_pt['SubId'].values[0]
            print("The needed product locates at:",product_path)
            print("The Subbasin Id of the lat, lon is:",Subid)
            return Subid,product_path 

        else:
            print("The point did not overlay with the routing product")
            return -1,-1
            
                        
    return -1, -1 