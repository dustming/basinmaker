# Copyright 2018-2018 Ming Han - ming.han(at)uwaterloo.ca
#
# License
# This file is part of Ming Han's personal code library.
#
# Ming Han's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ming Han's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Ming Han's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#
###########################################################################################
#
#This is a python script to use An automated ArcGIS toolbox for watershed delineation with lakes v1.tbx
#This script show how to use this tool box to generate lake-river routing sturcture with non HydroSHEDS dataset.
#
###########################################################################################

import numpy as np
import pandas as pd
import os

####define required inputs
#pythonfile = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Script for woring with HydroSHEDS DEM.py'
#accfile = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/accs.csv'
#Workfolder = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/"


pythonfile = sys.argv[1]
accfile = sys.argv[2]
Workfolder = sys.argv[3] + '/'

os.chdir(Workfolder)####set working folder
accs = pd.read_csv(accfile,sep=',')### read accs

##### generate .bat file to run python script with arcgis py 
#ofile2 = open(Workfolder+"runaccpys.bat", "w")
#ofile2.write('cd' + '      ' + Workfolder +'\n')
#ofile2.write(pythonexc + '      ' + 'acc_py_temp.py'+'\n')

#ofile2.close()

#####begin loop for each acc in the list

### generate 
for i in range(0,len(accs)):
    ##### generate python script for ith acc value
    os.chdir(Workfolder)
    iacc = accs['Accs'][i]
    ofile = open(Workfolder+"acc_py_temp.py", "w")
    ifile = open(pythonfile, "r")
    for line in ifile:
        if 'accthres =' not in line:
            ofile.write(line)
        else:
            ofile.write('accthres =' + "\""+str(iacc)+"\""+'\n')
    ifile.close()
    ofile.close()
    arcpy.AddMessage("Current acc is    " + str(iacc) )
    execfile('acc_py_temp.py') 
    # run python script for ith acc value
#    os.system('cd' + '      ' + Workfolder +'\n')
#    os.system(pythonexc + '      ' + 'acc_py_temp.py'+'\n')    
#    os.system("runaccpys.bat")
    
#### extract subbasins for each point of interest of each acc value.