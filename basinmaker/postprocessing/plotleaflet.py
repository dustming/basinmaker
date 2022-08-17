import os
import geopandas
import matplotlib.pyplot as plt ## only needed to plot figures
from ipywidgets import HTML,Layout,IntSlider, ColorPicker, jslink ## only needed to plot figures
from ipyleaflet import Map, GeoData, basemaps, LayersControl,Popup,Marker,Polygon,Choropleth,WidgetControl## only needed to plot figures
#import leafmap.foliumap  as leafmap

def plot_routing_product_with_ipyleaflet(path_to_product_folder,version_number = ''):
    product_folder = path_to_product_folder
    if version_number != '':
        version_number = '_'+version_number
    path_subbasin = os.path.join(product_folder,'finalcat_info'+version_number+'.shp')
    path_river = os.path.join(product_folder,'finalcat_info_riv'+version_number+'.shp')
    path_cllake = os.path.join(product_folder,'sl_connected_lake'+version_number+'.shp')
    path_ncllake = os.path.join(product_folder,'sl_non_connected_lake'+version_number+'.shp')

    subbasin = geopandas.read_file(path_subbasin)
    subbasin = subbasin.to_crs("EPSG:4326")

    if os.path.exists(path_river):
        river = geopandas.read_file(path_river)
        river = river.to_crs("EPSG:4326")


    if os.path.exists(path_cllake):
        cllake = geopandas.read_file(path_cllake)
        cllake = cllake.to_crs("EPSG:4326")
    if os.path.exists(path_ncllake):
        ncllake = geopandas.read_file(path_ncllake)
        ncllake = ncllake.to_crs("EPSG:4326")

    cx = subbasin['centroid_x'].mean()
    cy = subbasin['centroid_y'].mean()
    m = Map(center=(cy,cx), zoom = 9, basemap= basemaps.Esri.WorldTopoMap)

    subnolake = subbasin[subbasin['Lake_Cat'] == 0].copy(deep=True)
    sub_nolake_map = GeoData(geo_dataframe = subnolake,
                      style={'color': '#6E6E6E', 'fillColor': '#89CD66', 'opacity':1, 'weight':1, 'dashArray':'2', 'fillOpacity':1},
                      name = 'Subbasin without lakes')
    subcllake = subbasin[subbasin['Lake_Cat'] == 1].copy(deep=True)
    if len(subcllake) > 0:
        sub_cllake_map = GeoData(geo_dataframe = subcllake,
                          style={'color': '#6E6E6E', 'fillColor': '#CDE389', 'opacity':1, 'weight':1, 'dashArray':'2', 'fillOpacity':1},
                          name = 'Subbasin with connected lakes')
    subncllake = subbasin[subbasin['Lake_Cat'] == 2].copy(deep=True)
    if len(subncllake) > 0:
        sub_ncllake_map = GeoData(geo_dataframe = subncllake,
                          style={'color': '#6E6E6E', 'fillColor': '#F5F57A', 'opacity':1, 'weight':1, 'dashArray':'2', 'fillOpacity':1},
                          name = 'Subbasin with connected lakes')
    if os.path.exists(path_river):
        river_map = GeoData(geo_dataframe = river,
                          style={'color': '#0070FF', 'fillColor': '#0070FF', 'opacity':1, 'weight':2, 'dashArray':'2', 'fillOpacity':1},
                          name = 'River')
    if os.path.exists(path_cllake):
        cllake_map = GeoData(geo_dataframe = cllake,
                            style={'color': '#6E6E6E', 'fillColor': '#0070FF', 'opacity':1, 'weight':1, 'dashArray':'2', 'fillOpacity':1},
                            name = 'Non connected lakes')
    if os.path.exists(path_ncllake):
        ncllake_map = GeoData(geo_dataframe = ncllake,
                             style={'color': '#6E6E6E', 'fillColor': '#0070FF', 'opacity':0, 'weight':1, 'dashArray':'2', 'fillOpacity':1},
                             name = 'Connected lakes')

    m.add_layer(sub_nolake_map)
    if len(subcllake) > 0:
        m.add_layer(sub_cllake_map)
    if len(subncllake) > 0:
        m.add_layer(sub_ncllake_map)
    if os.path.exists(path_river):
        m.add_layer(river_map)
    if os.path.exists(path_cllake):
        m.add_layer(cllake_map)
    if os.path.exists(path_ncllake):
        m.add_layer(ncllake_map)

    m.add_control(LayersControl())
    m.layout.height="700px"

    return m

# def plot_routing_product_with_leafmap(path_to_product_folder,version_number = ''):
#     product_folder = path_to_product_folder
#     if version_number != '':
#         version_number = '_'+version_number
#     path_subbasin = os.path.join(product_folder,'finalcat_info'+version_number+'.shp')
#     path_river = os.path.join(product_folder,'finalcat_info_riv'+version_number+'.shp')
#     path_cllake = os.path.join(product_folder,'sl_connected_lake'+version_number+'.shp')
#     path_ncllake = os.path.join(product_folder,'sl_non_connected_lake'+version_number+'.shp')
#
#     subbasin = geopandas.read_file(path_subbasin)
#     subbasin = subbasin.to_crs("EPSG:4326")
#
#     if os.path.exists(path_river):
#         river = geopandas.read_file(path_river)
#         river = river.to_crs("EPSG:4326")
#
#
#     if os.path.exists(path_cllake):
#         cllake = geopandas.read_file(path_cllake)
#         cllake = cllake.to_crs("EPSG:4326")
#     if os.path.exists(path_ncllake):
#         ncllake = geopandas.read_file(path_ncllake)
#         ncllake = ncllake.to_crs("EPSG:4326")
#
#     subcol_olrrp = ['SubId','geometry','DowSubId','BasArea','Lake_Cat','DrainArea','DA_Chn_L','DA_Slope','DA_Chn_Slp','BasSlope']
#     subcol_na = ['SubId','geometry','DowSubId','BasArea','Lake_Cat','DrainArea','BasSlope']
#
#     rivcol = ['geometry','RivSlope','RivLength','BkfWidth','BkfDepth','Strahler','Q_Mean','FloodP_n','Ch_n']
#     lake_col =  ['geometry','Hylak_id','Lake_type','Lake_area','Vol_total','Depth_avg']
#
#     if version_number == 'v1-0':
#         subcol = subcol_olrrp
#     else:
#         subcol = subcol_na
#
#     subnolake = subbasin[subbasin['Lake_Cat'] == 0].copy(deep=True)
#
#     m = leafmap.Map(
#         draw_control=False,
#         measure_control=False,
#         fullscreen_control=False,
#         attribution_control=True,
#     )
#
#
#     m.add_gdf(subnolake[subcol], layer_name="Subbasin without lakes",style = {'color': '#6E6E6E', 'fillColor': '#89CD66', 'opacity':1, 'weight':1, 'dashArray':'2', 'fillOpacity':1})
#     labels = ["Subbasin without lakes"]
#     colors = ["#89CD66"]
#
#     subcllake = subbasin[subbasin['Lake_Cat'] == 1].copy(deep=True)
#     if len(subcllake) > 0:
#         style = {'color': '#6E6E6E', 'fillColor': '#CDE389', 'opacity':1, 'weight':1, 'dashArray':'2', 'fillOpacity':1}
#         name = 'Subbasin with connected lakes'
#         m.add_gdf(subcllake[subcol], layer_name=name,style = style)
#         labels = labels + ['Subbasin with connected lakes']
#         colors = colors + ['#CDE389']
#
#
#     subncllake = subbasin[subbasin['Lake_Cat'] == 2].copy(deep=True)
#     if len(subncllake) > 0:
#         style = {'color': '#6E6E6E', 'fillColor': '#F5F57A', 'opacity':1, 'weight':1, 'dashArray':'2', 'fillOpacity':1}
#         name = 'Subbasin with non connected lakes'
#         m.add_gdf(subncllake[subcol], layer_name=name,style = style)
#         labels = labels + ['Subbasin with non connected lakes']
#         colors = colors + ['#F5F57A']
#
#     if os.path.exists(path_river):
#         style = {'color': '#0070FF', 'fillColor': '#0070FF', 'opacity':1, 'weight':2, 'dashArray':'2', 'fillOpacity':1}
#         name = 'River'
#         m.add_gdf(river[rivcol], layer_name=name,style = style)
#         labels = labels + ['River/Lakes']
#         colors = colors + ['#0070FF']
#
#     if os.path.exists(path_cllake):
#         style = {'color': '#6E6E6E', 'fillColor': '#0070FF', 'opacity':1, 'weight':1, 'dashArray':'2', 'fillOpacity':1}
#         name = 'Non connected lakes'
#         m.add_gdf(cllake[lake_col], layer_name=name,style = style)
#
#
#     if os.path.exists(path_ncllake):
#         style = {'color': '#6E6E6E', 'fillColor': '#0070FF', 'opacity':1, 'weight':1, 'dashArray':'2', 'fillOpacity':1}
#         name = 'Cconnected lakes'
#         m.add_gdf(ncllake[lake_col], layer_name=name,style = style)
#
#     m.add_legend(title='Legend', labels=labels, colors=colors)
#     return m
