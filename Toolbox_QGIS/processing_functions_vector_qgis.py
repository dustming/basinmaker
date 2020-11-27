from qgis.core import *
import qgis
from qgis.analysis import QgsNativeAlgorithms
from qgis.PyQt.QtCore import *


def Copy_Pddataframe_to_shpfile(Path_shpfile,Pddataframe,link_col_nm_shp = 'SubId',
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
                    if  Attri_Name[icolnm] == 'Obs_NM' or Attri_Name[icolnm] == 'SRC_obs' or  Attri_Name[icolnm] == 'layer' or  Attri_Name[icolnm] == 'path'  :
                        sf[Attri_Name[icolnm]] = str(tarinfo[Attri_Name[icolnm]].values[0])
                    elif Attri_Name[icolnm] == 'cat':
                        continue
                    else:
                        sf[Attri_Name[icolnm]] = float(tarinfo[Attri_Name[icolnm]].values[0])
            else:
                for icolnm in range(0,len(UpdateColNM)):
                    sf[UpdateColNM[icolnm]] = float(tarinfo[UpdateColNM[icolnm]].values[0])

            layer_cat.updateFeature(sf)
    del layer_cat


    