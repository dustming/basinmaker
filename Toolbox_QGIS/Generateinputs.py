
###################################################################3
def Defcat(out,outletid):
    otsheds = np.full((1,1),outletid)
    Shedid = np.full((10000000,1),-99999999999999999)
    psid = 0
    rout = copy.copy(out)
    while len(otsheds) > 0:
        noutshd = np.full((10000000,1),-999999999999999)
        poshdid = 0
        for i in range(0,len(otsheds)):
            Shedid[psid] = otsheds[i]
            psid = psid + 1
            irow = np.argwhere(rout[:,1]==otsheds[i]).astype(int)
            for j in range(0,len(irow)):
                noutshd[poshdid] = rout[irow[j],0]
                poshdid = poshdid + 1
        noutshd = np.unique(noutshd)
        otsheds = noutshd[noutshd>=0]
    Shedid = np.unique(Shedid)
    Shedid = Shedid[Shedid>=0]
    return Shedid

def writeraster(w_filname,nraster,dataset):
    orvh = open(w_filname,"w")
    ncols = arcpy.GetRasterProperties_management(dataset, "COLUMNCOUNT")
    nrows = arcpy.GetRasterProperties_management(dataset, "ROWCOUNT")
    xllcorner = arcpy.GetRasterProperties_management(dataset, "LEFT")
    yllcorner = arcpy.GetRasterProperties_management(dataset, "BOTTOM")
    orvh.write("ncols      "+str(ncols) + "\n")
    orvh.write("nrows      "+ str(nrows) + "\n")
    orvh.write("xllcorner    "+str(xllcorner) + "\n")
    orvh.write("yllcorner    "+str(yllcorner) + "\n")
    orvh.write("cellsize     "+str(cellSize) + "\n")
    orvh.write("NODATA_value  -9999" + "\n")
    orvh.close()
    f_handle = open(w_filname, 'a')
    np.savetxt(f_handle,nraster,fmt='%i')
    f_handle.close()



def Generateinputdata_hydrosheds(hyshddem,hyshddir,hyshdacc,hyshdply,WidDep,Lakefile,Landuse,Landuseinfo,obspoint,OutHyID,OutHyID2,OutputFolder):

### function to generate inputs that is needed for lake-river routing network delineatation 
### Inputs    
#hyshddem = sys.argv[1]    ## the path for hydroshed dem data 
#hyshddir = sys.argv[2]    ## the path for hydroshed flow direction data 
#hyshdacc = sys.argv[3]    ## the path for hydroshed flow accumulation data 
#hyshdply = sys.argv[4]    ## the path for hydroshed level 12 catchments 
#WidDep = sys.argv[5]      ## the path for channel bankfull width and depth Polyline
#Lakefile = sys.argv[6]    ## the path for hydrolake polygon 
#Landuse = sys.argv[7]     ## the path for landuse data 
#Landuseinfo = sys.argv[8] ## the path for a csv table, which describe the flood plan mannning's coefficient of each land use type 
#obspoint = sys.argv[9]    ## the path for the observation point shp file 
#OutHyID = int(sys.argv[10]) ## the most downstream catchment hydroshed ID 
#OutHyID2 = int(sys.argv[11]) ## the most upstream catchment hydroshed ID
#OutputFolder = sys.argv[12] + "/"  ## the path for output folder to save processed files 

    hyinfocsv = hyshdply[:-3] + "dbf"
    VolThreshold = 0
    print(hyinfocsv)
    