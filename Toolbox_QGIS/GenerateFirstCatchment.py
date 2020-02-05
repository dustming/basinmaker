def WatershedDiscretizationToolset(OutputFolder,Thresholdacc):
    import os 
    print(os.environ['HOME'])
    
    
    dem = OutputFolder + '/' + 'dem.tif'
    processing.run("grass7:r.watershed", {'elevation':'dem.tif','depression':None,'flow':None,'disturbed_land':None,'blocking':None,'threshold':Thresholdacc,'max_slope_length':None,'convergence':5,'memory':1024,
    '-s':True,'-m':False,'-4':False,'-a':False,'-b':False,'accumulation': OutputFolder + '/' + 'acc2.tif',
    'drainage':OutputFolder + '/' + 'dir2.tif','basin':OutputFolder + '/' + 'cat1.tif','stream':OutputFolder + '/' + 'str.tif'})
    
