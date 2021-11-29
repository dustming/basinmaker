from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import sqlite3


def modify_str_r_to_add_sub_in_head_str(
    grass,
    garray,
    con,
    str_r,
    str_v,
    connected_lake,
    cat_no_lake,
    fdr_grass,
    max_memroy,  
):

    
    
    #### find str in str_r that did not over lay with the lakes 
    
    #### overlay str and lakes
    exp = "%s = if(isnull(int(%s)),null(),%s)" % ('str_with_lake', connected_lake, str_r)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    ### obtain connected lake ids
    str_id_with_lakes, temp = generate_stats_list_from_grass_raster(
        grass, mode=1, input_a='str_with_lake'
    )

    grass_raster_setnull_array(input = str_r,output = 'str_without_lake',values = str_id_with_lakes,grass = grass)
    
    # str_without_lakes = garray.array(mapname=str_r)
    # if (len(str_id_with_lakes) > 0):
    #     mask = np.isin(str_without_lakes, str_id_with_lakes)
    #     str_without_lakes[mask] = -9999
    # 
    # temparray = garray.array()
    # temparray[:, :] = str_without_lakes[:, :]
    # temparray.write(mapname='str_without_lake', overwrite=True)
    # grass.run_command("r.null", map='str_without_lake', setnull=[-9999, 0])
    # exp = "%s = int(%s)" % (
    #     'str_without_lake',
    #     'str_without_lake',
    # )
    # grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    
    
    
    ### add new outlet only to 'str_without_lake'
    
    grass.run_command("v.extract", input=str_v, layer = 1, type = 'point',output = 'str_p',where = 'type_code < 1',overwrite=True)

    sqlstat = "SELECT cat FROM %s" % str_v
    All_ids = pd.read_sql_query(sqlstat, con)
    max_subid = np.max(All_ids['cat'].values)
    
    grass.run_command(
        "v.db.addcolumn", map='str_p', columns="newcat int"
    )
    qcol = "cat + %s" %(str(int(max_subid)))
    
    grass.run_command(
        "v.db.update", map='str_p', column="newcat", qcol=qcol
    )
    grass.run_command(
        "v.to.rast", input='str_p', output='str_p', attribute_column="newcat",use="attr", overwrite=True
    )

    # remove points that 
    exp = "%s = if(isnull(%s),null(),%s)" % (
       'str_p',
       'str_without_lake',
       'str_p',
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    

    grass.run_command(
        "r.grow",
        input='str_p',
        output='str_p_gr',
        radius=3,
        overwrite=True,
    )
    
    exp = "%s = if(isnull(%s),null(),%s)" % (
       'str_t',
        str_r,
       'str_p_gr',
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    exp = "%s = if(isnull(%s),%s,%s)" % (
        str_r,
       'str_t',
        str_r,
       'str_t',
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
                    
    exp = "%s = int(%s)" % (
        str_r,
        str_r,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    grass_raster_r_stream_basins(
        grass,
        direction=fdr_grass,
        stream=str_r,
        basins=cat_no_lake,
        memory=max_memroy,
    )
    
    return 
        
    


