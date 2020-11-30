import copy 
from grass.script import array as garray
####
def grass_raster_setnull(grass,raster_nm,null_values,create_new_raster,new_raster_nm):
    
    if create_new_raster:
        grass_raster_copy(grass,raster_nm,new_raster_nm)
        grass.run_command('r.null', map=new_raster_nm,setnull=null_values)
    else:
        grass.run_command('r.null', map=raster_nm,setnull=null_values)

#####

def grass_raster_copy(grass,raster_nm_in,raster_nm_new):
    grass.run_command('g.copy',rast = (raster_nm_in,raster_nm_new),overwrite = True)
    
#####

def Return_Raster_As_Array_With_Db_Path(grassdb,grass_location,raster_mn):
    """Transfer an rater in grass database into np array
    Parameters
    ---------- 
    grassdb         : string
    Full path to a grass database 
    grass_location  : string
    location name in that grass database   
    raster_mn       : string
    raster name 
        
    Returns:
    -------
    Array            : array  
    np array of the raster. 
       
    """    
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location,create_opts='')
    Array = copy.deepcopy(garray.array(mapname=raster_mn))
    PERMANENT.close()
    Array[Array <= 0] = -9999
    return Array
    
###
def Return_Raster_As_Array_With_garray(garray_f,raster_mn):
    """Transfer an rater in grass database into np array
    Parameters
    ----------    
    raster_mn       : string
    raster name 
        
    Returns:
    -------
    Array            : array  
    np array of the raster. 
       
    """    
    Array = copy.deepcopy(garray_f.array(mapname=raster_mn))
    Array[Array <= 0] = -9999
    return Array    
    
    
def grass_raster_r_in_gdal(grass,raster_path,output_nm,location):
    """ import dem to target location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command("r.in.gdal", input = raster_path, output = output_nm, overwrite = True,location =location)
###    

def grass_raster_r_mask(grass,raster_nm):
    """ define grass working mask for current location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command('r.mask', raster=raster_nm, maskcats = '*',overwrite = True)
    
### 

def grass_raster_g_region(grass,raster_nm):
    """ define grass working region for current location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command('g.region', raster=raster_nm)
    
### 
