import copy 
from grass.script import array as garray
####
def grass_raster_setnull(grass,raster_nm,null_values,create_new_raster,new_raster_nm = '#'):
    
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
    
    
def grass_raster_r_in_gdal(grass,raster_path,output_nm,location = '#'):
    """ import dem to target location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    if location != '#':
        grass.run_command("r.in.gdal", input = raster_path, output = output_nm, overwrite = True,location =location)
    else:
        grass.run_command("r.in.gdal", input = raster_path, output = output_nm, overwrite = True)
###    




def grass_raster_r_accumulate(grass,direction,accumulation,flags):
    """ calculate flow accumulation from flow direction dataset 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command('r.accumulate',direction = direction, accumulation = accumulation, flags = flags,overwrite = True)

###



def grass_raster_r_external(grass,input,output):
    """ calculate flow accumulation from flow direction dataset 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command("r.external", input = input, output = output)

###

def grass_raster_r_clip(grass,input,output):
    """ clip raster with mask in grass env 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command("r.clip", input = input, output = output, overwrite = True)

###


def grass_raster_r_mask(grass,raster_nm,vector_nm = '#'):
    """ define grass working mask for current location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    if vector_nm == '#':
        grass.run_command('r.mask', raster=raster_nm, maskcats = '*',overwrite = True)
    else:
        grass.run_command('r.mask', vector=vector_nm, maskcats = '*',overwrite = True)
### 

def grass_raster_g_region(grass,raster_nm,zoom = '#'):
    """ define grass working region for current location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    if zoom == '#':
        grass.run_command('g.region', raster=raster_nm)
    else:
        grass.run_command('g.region', zoom=zoom)
### 

def grass_raster_r_unpack(grass,input,output):
    """ define grass working region for current location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command('r.unpack', input = input, output = output,overwrite = True)
    
### 



def grass_raster_r_reclass(grass,input,output,rules):
    """ reclassify grass raster dataset  
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command('r.reclass', input=input,output = output,rules =rules,overwrite = True)
    
### 


def grass_raster_r_watershed(grass,elevation,drainage,accumulation,flags):
    """ generate watershed from dem,output includes flow accumulation
        and flow direction  
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command('r.watershed',elevation = elevation, drainage = drainage,accumulation = accumulation,flags = flags, overwrite = True)
    
### 


def grass_raster_r_mapcalc(grass,expression):
    """ grass map calculator   
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command('r.mapcalc',expression = expression,overwrite = True)
    
### 

def grass_raster_create_raster_empty_raster(garray,raster_nm):
    """ grass create a raster with -9999 in current location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    temparray = garray.array()
    temparray[:,:] = -9999
    temparray.write(mapname=raster_nm, overwrite=True)
    
### 



def grass_raster_v_to_raster(grass,input,output,column):
    """ grass create a raster with -9999 in current location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command('v.to.rast',input = input,output = output,use = 'attr',attribute_column = column,overwrite = True)
    
### 








def grass_raster_r_water_outlet(grass,input_dir_nm,output_watshed_nm,outlet_coordinates):
    """ grass generate watershed based on outlet points    
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command('r.water.outlet',input = input_dir_nm, output = output_watshed_nm, coordinates  = outlet_coordinates,overwrite = True)
    
### 


def grass_raster_r_out_gdal(grass,input_nm,output,format= 'GTiff'):
    """ grass export raster in grass working enviroment to other folder in
        different format.  
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command('r.out.gdal', input = input_nm,output =output,format= 'GTiff',overwrite = True)
    
### 




