from qgis.core import *
import qgis
from qgis.analysis import QgsNativeAlgorithms
from qgis.PyQt.QtCore import *


def Copy_Pddataframe_to_shpfile(Path_shpfile,Pddataframe,link_col_nm_shp = 'SubId'
                                ,link_col_nm_df = 'SubId',UpdateColNM = ['#']):
    """ Function modify attribute table of Path_shpfile using value from Pddataframe
    Parameters
    ----------
    Path_shpfile                        : shpfile 
        Path to the shpfile 
    Pddataframe                         : dataframe
        Dataframe constains data that will be used to update attribute table of 
        Path_shpfile
    link_col_nm_shp                     : string 
        The column name that link Pddataframe attribute table in Path_shpfile
    link_col_nm_df                      : string 
        The column name that link Path_shpfile attribute table in Pddataframe
    UpdateColNM                         : list 
        It is a list of column name, it is equal to '#', all column value will 
        be updated, otherwise only column name in UpdateColNM will be updated. 
    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated 
    """
    
    layer_cat=QgsVectorLayer(Path_shpfile,"")
    Attri_Name = layer_cat.fields().names()
    features = layer_cat.getFeatures()
    with edit(layer_cat):
        for sf in features:
            Atti_Valu    = sf.attributes()
            sf_subid     = sf[link_col_nm_shp]
            tarinfo      = Pddataframe[Pddataframe[link_col_nm_df] == sf_subid]

            if UpdateColNM[0] == '#':
                for icolnm in range(0,len(Attri_Name)):     ### copy infomaiton
                    if  Attri_Name[icolnm] == 'Obs_NM' or Attri_Name[icolnm] == 'SRC_obs':
                        sf[Attri_Name[icolnm]] = str(tarinfo[Attri_Name[icolnm]].values[0])
                    elif Attri_Name[icolnm] == 'cat' or  Attri_Name[icolnm] == 'layer' or  Attri_Name[icolnm] == 'path':
                        continue
                    else:
                        sf[Attri_Name[icolnm]] = float(tarinfo[Attri_Name[icolnm]].values[0])
            else:
                for icolnm in range(0,len(UpdateColNM)):
                    sf[UpdateColNM[icolnm]] = float(tarinfo[UpdateColNM[icolnm]].values[0])

            layer_cat.updateFeature(sf)
    del layer_cat


def Remove_Unselected_Lake_Attribute_In_Finalcatinfo(Path_Finalcatinfo,Conn_Lake_Ids):
    """ Functions will set lake id not in Conn_Lake_Ids to -1.2345 in attribute 
        table of Path_Finalcatinfo
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated 
    """
    
    layer_cat=QgsVectorLayer(Path_Finalcatinfo,"")
    Attri_Name = layer_cat.fields().names()
    features = layer_cat.getFeatures()
    with edit(layer_cat):
        for sf in features:
            sf_subid        = float(sf['HyLakeId'])

            if sf_subid in Conn_Lake_Ids or float(sf['IsLake']) == 2:
                continue
            sf['HyLakeId']      = float(-1.2345)
            sf['LakeVol']       = float(-1.2345)
            sf['LakeArea']      = float(-1.2345)
            sf['LakeDepth']     = float(-1.2345)
            sf['Laketype']      = float(-1.2345)
            sf['IsLake']        = float(-1.2345)
            layer_cat.updateFeature(sf)
    del layer_cat
    return
    
    
#########
def Add_centroid_to_feature(Path_feagure,centroidx_nm = '#',centroidy_nm='#'):
    """ Functions will add centorid x y to Path_feagure
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated 
    """
    layer_cat=QgsVectorLayer(Path_feagure,"")
    Attri_Name = layer_cat.fields().names()
    features = layer_cat.getFeatures()
    with edit(layer_cat):
        for sf in features:
            centroidxy = sf.geometry().centroid().asPoint()
            sf[centroidx_nm] = centroidxy[0]
            sf[centroidy_nm] = centroidxy[1]
            layer_cat.updateFeature(sf)
    del layer_cat
    return

def Selectfeatureattributes(processing,Input = '#',Output='#',Attri_NM = '#',Values = []):
    """ Functions extract features from Input, based on values in column Attri_NM
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated 
    """
    exp =Attri_NM + '  IN  (  ' +  str(int(Values[0]))
    for i in range(1,len(Values)):
        exp = exp + " , "+str(int(Values[i]))
    exp = exp + ')'
    processing.run("native:extractbyexpression", {'INPUT':Input,'EXPRESSION':exp,'OUTPUT':Output})
    
        