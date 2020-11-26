
def Is_Point_Close_To_Id_In_Raster(prow,pcol,nrows,ncols,id,raster_array):
    """Check if the point is around grids with value equal to Id in raster_array

    Parameters
    ----------

    Returns:
    -------
    noout        : logical
       True  : close to grids with value euqal to id 
       False : not close to grids with value to id 
    """    
    Isclose  = False 
    n_grids_eq_id = 0
    if prow != 0 and prow != nrows -1 and pcol != 0 and pcol != ncols-1:
        if raster_array[prow-1,pcol+1] == id:
            Isclose=True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow-1,pcol-1] == id:
            Isclose=True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow-1,pcol] == id:
            Isclose=True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow,pcol+1] == id:
            Isclose=True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow,pcol-1] == id:
            Isclose=True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow+1,pcol-1] == id:
            Isclose=True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow+1,pcol+1] == id:
            Isclose=True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow+1,pcol] == id:
            Isclose=True
            n_grids_eq_id = n_grids_eq_id + 1
            
    return Isclose,n_grids_eq_id