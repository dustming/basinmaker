def Dbf_To_Dataframe(file_path):
    """Transfer an input dbf file to dataframe

    Parameters
    ----------
    file_path   : string
    Full path to a shapefile

    Returns:
    -------
    dataframe   : datafame
    a pandas dataframe of attribute table of input shapefile
    """
    tempinfo = Dbf5(file_path[:-3] + "dbf")
    dataframe = tempinfo.to_dataframe().copy()
    return dataframe
