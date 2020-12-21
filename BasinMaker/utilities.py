import os

from simpledbf import Dbf5


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


def WriteStringToFile(Out_String, File_Path, WriteMethod):
    """Write String to a file

    Function that used to write Out_String to a file located at the File_Path.

    Parameters
    ----------
    Out_String            : string
        The string that will be writed to the file located at File_Path
    File_Path             : string
        Path and filename of file that will be modified or created
    WriteMethod           : {'a','w'}
        If WriteMethod = "w", a new file will be created at the File_Path
        If WriteMethod = "a", the Out_String will be added to exist file

    Notes
    ------
        The file located at the File_Path will be modified or created

    Returns
    -------
        None

    Examples
    --------
    >>> from WriteRavenInputs import WriteStringToFile
    >>> Out_String = 'sometest at line 1\n some test at line 2\n some test at line 3\n'
    >>> File_Path  = 'C:/Path_to_the_Flie_with_file_name'
    >>> WriteStringToFile(Out_String = Out_String,File_Path = File_Path,WriteMethod = 'w')

    """

    if os.path.exists(
        File_Path
    ):  ### if file exist, we can either modify or overwrite it
        with open(File_Path, WriteMethod) as f:
            f.write(Out_String)
    else:  ## create a new file anyway, since file did not exist
        with open(File_Path, "w") as f:
            f.write(Out_String)


def write_grass_and_arcgis_fdr_rules(grassdb):

    Strlist = [
        "1 = 128",
        "2 = 64",
        "3 = 32",
        "4 = 16",
        "5 = 8",
        "6 = 4",
        "7 = 2",
        "8 = 1",
        "* = NULL",
    ]
    Str = "\n".join(Strlist)
    WriteStringToFile(Str, os.path.join(grassdb, "Grass2ArcgisDIR.txt"), "w")

    Strlist = [
        "1 = 8",
        "2 = 7",
        "4 = 6",
        "8 = 5",
        "16 = 4",
        "32 = 3",
        "64 = 2",
        "128	= 1",
        "* = NUL",
    ]
    Str = "\n".join(Strlist)
    WriteStringToFile(Str, os.path.join(grassdb, "Arcgis2GrassDIR.txt"), "w")

    return
