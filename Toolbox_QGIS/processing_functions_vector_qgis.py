from qgis.core import *
import qgis
from qgis.analysis import QgsNativeAlgorithms
from qgis.PyQt.QtCore import *


def Copy_Pddataframe_to_shpfile(Path_shpfile,Pddataframe,link_col_nm_shp = 'SubId'
                                ,link_col_nm_df = 'SubId',UpdateColNM = ['#'],Input_Is_Feature_In_Mem = False):
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
    if Input_Is_Feature_In_Mem:
        layer_cat = Path_shpfile
    else:
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
    if Input_Is_Feature_In_Mem:
        return layer_cat
    else:
        del layer_cat
        return 

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



def Copyfeature_to_another_shp_by_attribute(Source_shp,Target_shp,Col_NM='SubId',Values=[-1],Attributes = [-1]):

    """ Functions that will copy features in Source_shp to Target_shp 
    based on attribute values in Values
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated 
    """
    
    layer_src=QgsVectorLayer(Source_shp,"")
    layer_trg=QgsVectorLayer(Target_shp,"")

    src_features = layer_src.getFeatures()

    Selected_Features = []
    for sf in src_features:
        #centroidxy = sf.geometry().centroid().asPoint()
        Select_value = sf[Col_NM]
        if Select_value in Values:
            src_geometry =  sf.geometry()
            attribute = Attributes.loc[Attributes[Col_NM] == Select_value].values
            temp_feature=QgsFeature()
            temp_feature.setGeometry(src_geometry)
            temp_feature.setAttributes(attribute.tolist()[0])
            Selected_Features.append(temp_feature)

    layer_trg.startEditing()
    layer_trg.addFeatures(Selected_Features)
    layer_trg.commitChanges()
    layer_trg.updateExtents()
    del layer_src
    del layer_trg
    
        


def Add_New_SubId_To_Subregion_shpfile(processing,context,Layer,SubID_info = '#', OutputPath = '#',Region_ID = 1):
    """ Asign new subbasin Id to each subbasins in each subregion
    Parameters
    ----------
    processing                        : qgis object
    context                           : qgis object
    Layer                             : vector layer
        it is the subbasin polygon or polyline of watershed delineation
        result in each subregion
    SubID_info                        : dataframe
        it is a dataframe contains new subbasin id for each subbasin
        in Layer
    OutputPath                        : string
        Path to the output file
    Region_ID                         : integer
        it is the subregion id of layer

    Notes
    -------
        the output will be the same type of input Layer
        stored in OutputPath, a new SubId will be given
    Returns:
    -------
        None
    """
            
    qgis_vector_field_calculator(processing = processing, context = context,FORMULA =str(int(Region_ID)),FIELD_NAME = 'Region_ID',INPUT =Layer,OUTPUT =OutputPath)
    

    layer_new=QgsVectorLayer(OutputPath,"")

    features = layer_new.getFeatures()
    with edit(layer_new):
        for sf in features:
            cSubId = int(sf['SubId'])
            cDowSubId = int(sf['DowSubId'])
            nSubId = SubID_info.loc[SubID_info['SubId'] == cSubId]['nSubId']
            if len(SubID_info.loc[SubID_info['SubId'] == cDowSubId]) == 0:
                nDowSubId = -1
            else:
                nDowSubId = SubID_info.loc[SubID_info['SubId'] == cDowSubId]['nSubId']
            nSeg_ID = SubID_info.loc[SubID_info['SubId'] == cSubId]['nSeg_ID']
            sf['SubId'] = int(nSubId)
            sf['DowSubId'] = int(nDowSubId)
            sf['Seg_ID'] =   int(nSeg_ID)
            layer_new.updateFeature(sf)
    del layer_new
    return
    
def qgis_vector_field_calculator(processing,context,FORMULA,INPUT,OUTPUT,FIELD_NAME,
                                 FIELD_PRECISION = 0,FIELD_TYPE = 0,NEW_FIELD = True,
                                 FIELD_LENGTH = 10):
    """ qgis filed calcuator 
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """    
    alg_params = {
        'FIELD_LENGTH': FIELD_LENGTH,
        'FIELD_NAME': FIELD_NAME,
        'FIELD_PRECISION': FIELD_PRECISION,
        'FIELD_TYPE': FIELD_TYPE,
        'FORMULA':FORMULA, 
        'INPUT': INPUT,
        'NEW_FIELD': NEW_FIELD,
        'OUTPUT':OUTPUT
        }

    out = processing.run('qgis:fieldcalculator', alg_params, context=context)
    return out    



def qgis_vector_dissolve(processing,context,INPUT,FIELD,OUTPUT):
    """ qgis dissolve input vector based on values in FIELD list
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """    
    out = processing.run("native:dissolve", {'INPUT':INPUT,'FIELD':FIELD,'OUTPUT':OUTPUT},context = context)
    return out 


def qgis_vector_fix_geometries(processing,context,INPUT,OUTPUT):
    """ qgis fixgeometries
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """    
    out = processing.run("native:fixgeometries", {'INPUT':INPUT,'OUTPUT':OUTPUT})
    return out 



def qgis_vector_merge_vector_layers(processing,context,INPUT_Layer_List,OUTPUT):
    """ qgis merge_vector_layers
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """    
    out = processing.run("native:mergevectorlayers", {'LAYERS':INPUT_Layer_List,'OUTPUT':OUTPUT})
    
    return out 
    
    
def qgis_vector_return_crs_id(processing,context,INPUT_Layer,Input_Is_Feature_In_Mem = True):
    """ qgis return vector layer projection crs id 
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """
    if Input_Is_Feature_In_Mem:   
        out =INPUT_Layer.crs().authid()
    else:
        layer = QgsVectorLayer(INPUT_Layer, "")
        out =layer.crs().authid()
        
    return out 


def qgis_vector_extract_by_attribute(processing,context,INPUT_Layer,FIELD,OPERATOR,VALUE,OUTPUT):
    """ qgis extract vector by attribute 
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """
    out  = processing.run("native:extractbyattribute", {'INPUT':INPUT_Layer,'FIELD':FIELD,'OPERATOR':OPERATOR,'VALUE':VALUE,'OUTPUT':OUTPUT})
        
    return out 


def qgis_vector_add_attributes(processing,context,INPUT_Layer,attribute_list):
    """ qgis add attributes to vector  
    ----------
    Notes
    -------

    Returns:
    -------
        None, 
    """    
    INPUT_Layer.dataProvider().addAttributes(attribute_list)
    INPUT_Layer.updateFields()
    INPUT_Layer.commitChanges()
    return INPUT_Layer


def qgis_vector_get_attributes(processing,context,INPUT_Layer,attribute_NM):
    """ qgis retrun attribute value of vector  
    ----------
    Notes
    -------

    Returns:
    -------
        None, 
    """    
    if attribute_NM == 'count':
        out = INPUT_Layer.featureCount()
        return out 
    elif attribute_NM == 'field_name':
        out = INPUT_Layer.fields().names()
    elif attribute_NM == 'features':
        out = INPUT_Layer.getFeatures()
        return out 
    else:
        print("wrong attribute name")
        return -1         
#### 
    
        

def qgis_vector_union_two_layers(processing,context,INPUT,OVERLAY,OUTPUT,OVERLAY_FIELDS_PREFIX = ''):
    """ qgis union two layers 
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """
    
    out = processing.run("native:union", {'INPUT':INPUT,'OVERLAY':OVERLAY,'OVERLAY_FIELDS_PREFIX':OVERLAY_FIELDS_PREFIX,'OUTPUT':OUTPUT},context = context)
            
    return out 
    
def qgis_vector_reproject_layers(processing,context,INPUT,TARGET_CRS,OUTPUT):
    """ qgis function reproject vector layer 
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """  
    out = processing.run("native:reprojectlayer", {'INPUT':INPUT,'TARGET_CRS':QgsCoordinateReferenceSystem(TARGET_CRS),'OUTPUT':OUTPUT})
    return out 


def qgis_vector_clip(processing,context,INPUT,OVERLAY,OUTPUT):
    """ qgis function reproject vector layer 
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """  
    layer_clip = processing.run("native:clip", {'INPUT':INPUT,'OVERLAY':OVERLAY,'OUTPUT':OUTPUT})
    
    return out
    
    






def qgis_vector_create_spatial_index(processing,context,INPUT):
    """ qgis function create spatial index
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """  
    out = processing.run("qgis:createspatialindex", {'INPUT':INPUT})
    return out 
    

def Clean_Attribute_Name(Input,FieldName_List,Input_Is_Feature_In_Mem = False,Col_NM_Max ='SubId'):
    """ Function clean feature attribute table, all colnmun not in FieldName_List
        will be removed 
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """  
    
    fieldnames = set(FieldName_List)
    if Input_Is_Feature_In_Mem:
       layer_cat = Input
    else:
        layer_cat  =QgsVectorLayer(Input, "")
        
    field_ids  = []
    for field in layer_cat.fields():
        if field.name() not in fieldnames:
            field_ids.append(layer_cat.dataProvider().fieldNameIndex(field.name()))
        if field.name() == Col_NM_Max:
            max_subbasin_id = layer_cat.maximumValue(layer_cat.dataProvider().fieldNameIndex(field.name()))
            
    layer_cat.dataProvider().deleteAttributes(field_ids)
    layer_cat.updateFields()
    layer_cat.commitChanges()
    
    if Input_Is_Feature_In_Mem:
       return layer_cat,max_subbasin_id
    else:
       del layer_cat 
       return max_subbasin_id

########
                    