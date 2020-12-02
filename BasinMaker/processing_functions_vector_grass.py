def grass_raster_v_in_org(grass, input_path, output_vector_nm, location):
    """grass load vector into target grass location
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command(
        "v.in.ogr",
        input=input_path,
        output=output_vector_nm,
        overwrite=True,
        location=location,
    )


###


def grass_raster_v_import(grass, input_path, output_vector_nm):
    """grass load vector into target grass location
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command(
        "v.import", input=input_path, output=output_vector_nm, overwrite=True
    )


###


def grass_raster_v_db_join(grass, map, column, other_table, other_column):
    """grass join vector database
    Parameters
    ----------

    Returns:
    -------

    """
    grass.run_command(
        "v.db.join",
        map=map,
        column=column,
        other_table=other_table,
        other_column=other_column,
        overwrite=True,
    )


###
