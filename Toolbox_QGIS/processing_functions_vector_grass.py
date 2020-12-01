


def grass_raster_v_in_org(grass,input_path,output_vector_nm,location):
    """ grass load vector into target grass location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command("v.in.ogr", input = input_path,output = output_vector_nm, overwrite = True,location =location)
### 



def grass_raster_v_import(grass,input_path,output_vector_nm):
    """ grass load vector into target grass location 
    Parameters
    ----------    
        
    Returns:
    -------
       
    """      
    grass.run_command("v.import", input = input_path, output = output_vector_nm, overwrite = True)
### 




