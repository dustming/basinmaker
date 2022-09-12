import os
import geopandas
from basinmaker.func.pdtable import *
import numpy as np
from basinmaker.postprocessing.postprocessingfunctions import (
    combine_catchments_covered_by_the_same_lake_method,
)

def inverse_topology_subids(datafame,end_subid,start_subid):
    datafame['ndownsubid'] =datafame['DowSubId']
    csubid = end_subid
    up_subid = -1

    false_drainage_csubid = datafame.loc[datafame['DowSubId'] == csubid,'SubId'].values
#    print(len(false_drainage_csubid))
    if len(false_drainage_csubid) > 2:
        cmask = datafame['SubId'] == csubid
    else:
        false_drainage_csubid = false_drainage_csubid[false_drainage_csubid != up_subid]
        mask_subid = np.append(false_drainage_csubid,[csubid])
        cmask = datafame['SubId'].isin(mask_subid)

#    print(csubid,datafame.loc[cmask,'DowSubId'].values[0])

    while csubid != -1:
        datafame.loc[cmask,'ndownsubid'] = up_subid
        up_subid = csubid
        csubid = datafame.loc[datafame['SubId'] == csubid,'DowSubId'].values[0]

        false_drainage_csubid = datafame.loc[datafame['DowSubId'] == csubid,'SubId'].values
        if len(false_drainage_csubid) > 2:
            cmask = datafame['SubId'] == csubid
        else:
            false_drainage_csubid = false_drainage_csubid[false_drainage_csubid != up_subid]
            mask_subid = np.append(false_drainage_csubid,[csubid])
            cmask = datafame['SubId'].isin(mask_subid)

#        print(csubid,up_subid,datafame.loc[cmask,'SubId'].values)

    datafame['DowSubId'] = datafame['ndownsubid']
    datafame = datafame.drop(columns=['ndownsubid'])

    datafame = streamorderanddrainagearea(datafame)

    datafame = update_non_connected_catchment_info(datafame)
    datafame.loc[datafame['RivLength'] == -1.2345,'RivSlope'] = -1.2345
    datafame.loc[datafame['RivLength'] == -1.2345,'FloodP_n'] = -1.2345
    datafame.loc[datafame['RivLength'] == -1.2345,'Max_DEM'] = -1.2345
    datafame.loc[datafame['RivLength'] == -1.2345,'Min_DEM'] = -1.2345
    datafame.loc[datafame['RivLength'] == -1.2345,'Ch_n'] = -1.2345


    return datafame


def inverse_topology(path_to_product_folder,path_to_output_folder,version_number,end_subid,start_subid):

    if not os.path.exists(path_to_output_folder):
        os.mkdir(path_to_output_folder)

    path_subbasin = os.path.join(path_to_product_folder,'catchment_without_merging_lakes_'+version_number+'.shp')
    path_river = os.path.join(path_to_product_folder,'river_without_merging_lakes_'+version_number+'.shp')

    subbasin = geopandas.read_file(path_subbasin)

    river = geopandas.read_file(path_river)
    river = river[['SubId','geometry']]
    subbasin = inverse_topology_subids(subbasin,end_subid,start_subid)
    river = pd.merge(river, subbasin.drop(columns=['geometry']), on = 'SubId', how='inner')

    subbasin.to_file(os.path.join(path_to_output_folder,'catchment_without_merging_lakes_'+version_number+'.shp'))
    river.to_file(os.path.join(path_to_output_folder,'river_without_merging_lakes_'+version_number+'.shp'))

    combine_catchments_covered_by_the_same_lake_method(
        Routing_Product_Folder = path_to_output_folder,
        qgis_prefix_path='#',
        gis_platform='purepy',
    )
