from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import sqlite3


def modify_str_r_to_add_sub_in_head_str(
    grass,
    grassdb,
    grass_location,
    str_r,
    str_v,
    max_memroy,
):

    con = sqlite3.connect(
        os.path.join(grassdb, grass_location, "PERMANENT", "sqlite", "sqlite.db")
    )

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

