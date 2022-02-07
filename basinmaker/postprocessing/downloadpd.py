
import os 
import wget 


# define function two download routing product for a given gauge
def Download_Routing_Product_For_One_Gauge(gauge_name,product_name):
    if product_name == 'NA':
        version = 'v2-1'
        gauge_info = pd.read_csv("https://github.com/dustming/RoutingTool/wiki/Files/obs_gauges_NA_v2-1.csv")
        gauge_info_sl = gauge_info[gauge_info['Obs_NM'] == gauge_name]
        if len(gauge_info_sl) < 1:
            print("The gauge ",gauge_name,"is not included in the routing product ")
            SubId = -1 
            product_path = '#'
        else:
            if gauge_info_sl['Use_region'].values[0] < 1:
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
                region_id = int(gauge_info_sl['Region'].values[0])
                subreg_id = int(gauge_info_sl['Sub_Reg'].values[0])
                url = "http://hydrology.uwaterloo.ca/basinmaker/data/original/drainage_region_%s_v2-1.zip" % (str(region_id).zfill(4))
                wget.download(url)
                os.system('unzip drainage_region_%s_v2-1.zip' % (str(region_id).zfill(4)))
                SubId =  gauge_info_sl['SubId'].values[0]
                product_name = "drainage_region_%s_v2-1" % (str(region_id).zfill(4))
                product_path = os.path.join(os.getcwd(),product_name)
                print("The needed product locates at:",product_path)
                print("The Subbasin Id of the interested gauge is:",SubId)
                                
    return SubId,product_path