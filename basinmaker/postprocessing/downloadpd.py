
import os
import wget
import pandas as pd
import numpy as np
import gdown
# define function two download routing product for a given gauge
def Download_Routing_Product_For_One_Gauge(gauge_name,product_name,region='#',subreg = '#'):
    if product_name == 'NALRP':
        version = 'v2-1'
        if gauge_name != '#':
            gauge_info = pd.read_csv("https://github.com/dustming/RoutingTool/wiki/Files/obs_gauges_NA_v2-1.csv")
            gauge_info_sl = gauge_info[gauge_info['Obs_NM'] == gauge_name]
            if len(gauge_info_sl) < 1:
                print("The gauge ",gauge_name,"is not included in the routing product ")
                SubId = -1
                product_path = '#'
            else:
                if gauge_info_sl['Use_region'].values[0] < 1:
#                    print(gauge_info_sl)
#                    print(gauge_info_sl['Region'].values[0],gauge_info_sl['Sub_Reg'].values[0])
                    region_id = int(gauge_info_sl['Region'].values[0])
                    subreg_id = int(gauge_info_sl['Sub_Reg'].values[0])
                    url = "http://hydrology.uwaterloo.ca/basinmaker/data/original/drainage_region_%s_%s_v2-1.zip" % (str(region_id).zfill(4),str(subreg_id).zfill(5))
                    wget.download(url)
                    os.system('unzip drainage_region_%s_%s_%s.zip' % (str(region_id).zfill(4),str(subreg_id).zfill(5),version))
                    SubId =  gauge_info_sl['SubId'].values[0]
                    product_name = "drainage_region_%s_%s_%s" % (str(region_id).zfill(4),str(subreg_id).zfill(5),version)
                    product_path = os.path.join(os.getcwd(),product_name)
                    print("The needed product locates at:",product_path)
                    print("The Subbasin Id of the interested gauge is:",SubId)
                else:
#                    print(gauge_info_sl)
#                    print(gauge_info_sl['Region'].values[0],gauge_info_sl['Sub_Reg'].values[0])
                    region_id = int(gauge_info_sl['Region'].values[0])
#                    subreg_id = int(gauge_info_sl['Sub_Reg'].values[0])
                    url = "http://hydrology.uwaterloo.ca/basinmaker/data/original/drainage_region_%s_v2-1.zip" % (str(region_id).zfill(4))
                    wget.download(url)
                    os.system('unzip drainage_region_%s_v2-1.zip' % (str(region_id).zfill(4)))
                    SubId =  gauge_info_sl['SubId'].values[0]
                    product_name = "drainage_region_%s_v2-1" % (str(region_id).zfill(4))
                    product_path = os.path.join(os.getcwd(),product_name)
                    print("The needed product locates at:",product_path)
                    print("The Subbasin Id of the interested gauge is:",SubId)

        elif region != '#' and subreg =='#':
            region_id = region
            subreg_id = '-'
            url = "http://hydrology.uwaterloo.ca/basinmaker/data/original/drainage_region_%s_v2-1.zip" % (str(region_id).zfill(4))
            wget.download(url)
            os.system('unzip drainage_region_%s_v2-1.zip' % (str(region_id).zfill(4)))
            product_name =  "drainage_region_%s_v2-1" % (str(region_id).zfill(4))
            product_path = os.path.join(os.getcwd(),product_name)
            SubId =  -1

        elif region != '#' and subreg !='#':
            print("todo")
        else:
            print("Wrong option provided")


    if product_name == 'OLRP':
        version = 'v1-0'
        gauge_info = pd.read_csv("https://github.com/dustming/RoutingTool/wiki/Files/OIH_gauge_info.csv")
        productinfo = pd.read_csv("https://github.com/dustming/RoutingTool/wiki/Files/OLRRP.csv")
        if gauge_name != '#':
            gauge_info_sl = gauge_info[gauge_info['Obs_NM'] == gauge_name]
            if len(gauge_info_sl) < 1 or gauge_info_sl['SubId'].values[0] < 0:
                print("The gauge ",gauge_name,"is not included in the routing product ")
                SubId = -1
                product_path = '#'
            else:
                if gauge_info_sl['Use_region'].values[0] < 1:
                    region_id = gauge_info_sl['Region'].values[0]
                    subreg_id = gauge_info_sl['Sub_Reg'].values[0]
                    mask1 = productinfo['Region'] == region_id
                    mask2 = productinfo['Sub_Reg'] == subreg_id
                    mask = np.logical_and(mask1,mask2)
                    url_veiw = productinfo.loc[mask,'Download'].values[0]

                    url_veiw = url_veiw.split("/")
                    url = "https://drive.google.com/u/0/uc?id=%s&export=download"%(url_veiw[5])
                    output = 'drainage_region_%s_%s_%s.zip' % (region_id,subreg_id,version)

                    os.system('gdown %s -O %s' % (url_veiw[5],output))

                    os.system('unzip drainage_region_%s_%s_%s.zip' % (region_id,subreg_id,version))
                    SubId =  gauge_info_sl['SubId'].values[0]
                    product_name = 'drainage_region_%s_%s_%s' % (region_id,subreg_id,version)
                    product_path = os.path.join(os.getcwd(),product_name)
                    print("The needed product locates at:",product_path)
                    print("The Subbasin Id of the interested gauge is:",SubId)
                else:
                    region_id = gauge_info_sl['Region'].values[0]
                    subreg_id = '-'
                    mask1 = productinfo['Region'] == region_id
                    mask2 = productinfo['Sub_Reg'] == subreg_id
                    mask = np.logical_and(mask1,mask2)
                    url_veiw = productinfo.loc[mask,'Download'].values[0]
                    url_veiw = url_veiw.split("/")
                    url = "https://drive.google.com/u/0/uc?id=%s&export=download"%(url_veiw[5])
                    output = 'drainage_region_%s_%s.zip' % (region_id,version)

                    os.system('gdown %s -O %s' % (url_veiw[5],output))

                    os.system('unzip drainage_region_%s_%s.zip' % (region_id,version))
                    SubId =  gauge_info_sl['SubId'].values[0]
                    product_name ='drainage_region_%s_%s' % (region_id,version)
                    product_path = os.path.join(os.getcwd(),product_name)
                    print("The needed product locates at:",product_path)
                    print("The Subbasin Id of the interested gauge is:",SubId)

        elif region != '#' and subreg =='#':
            region_id = region
            subreg_id = '-'
            mask1 = productinfo['Region'] == region_id
            mask2 = productinfo['Sub_Reg'] == subreg_id
            mask = np.logical_and(mask1,mask2)
            url_veiw = productinfo.loc[mask,'Download'].values[0]
            url_veiw = url_veiw.split("/")
            url = "https://drive.google.com/u/0/uc?id=%s&export=download"%(url_veiw[5])
            output = 'drainage_region_%s_%s.zip' % (region_id,version)

            os.system('gdown %s -O %s' % (url_veiw[5],output))

            os.system('unzip drainage_region_%s_%s.zip' % (region_id,version))
            SubId =  -1
            product_name ='drainage_region_%s_%s' % (region_id,version)
            product_path = os.path.join(os.getcwd(),product_name)

        elif region != '#' and subreg !='#':
            print("todo")
        else:
            print("Wrong option provided")

    return SubId,product_path
