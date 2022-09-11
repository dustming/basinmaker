import os
import geopandas

def inverse_topology(path_to_product_folder,version_number = '',end_subid,start_subid):

    path_subbasin = os.path.join(product_folder,'finalcat_info'+version_number+'.shp')
    path_river = os.path.join(product_folder,'finalcat_info_riv'+version_number+'.shp')

    subbasin = geopandas.read_file(path_subbasin)
    subbasin = subbasin.to_crs("EPSG:4326")
