import pytest
import sys
from WriteRavenInputs import DownloadStreamflowdata_CA
import numpy as np

    
    
def test_DownloadStreamflowdata_CA(HYDAT_Path):
    
    '''
    CA_HYDAT       (string):     Path and filename of previously downloaded 
                                 external database containing streamflow observations, 
                                 e.g. HYDAT for Canada ("Hydat.sqlite3"). 
                                 Read from command line
    '''
    ###Define
    Station_NM = '05PC019'
    StartYear  = 2010
    EndYear    = 2011
    CA_HYDAT   = HYDAT_Path
    
    Expect_Annual_Discharge = 155933.6000137329
    
    flowdata,obs_DA,obtaindata = DownloadStreamflowdata_CA(Station_NM,CA_HYDAT,StartYear,EndYear)
    calculated_Annual_Discharge = np.nansum(flowdata['Flow'].values)
    
    assert Expect_Annual_Discharge == pytest.approx(calculated_Annual_Discharge, 0.1)
