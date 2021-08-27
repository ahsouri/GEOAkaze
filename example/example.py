from GEOAkaze import GEOAkaze
import os.path
import matplotlib.pyplot as plt
import numpy as np

slave_f = []

slave_f.append('/Users/asouri/Documents/Methane_SAT_OSSEs/Main/GEOAkaze/GEOAkaze/data/MethaneAIR_rf03/MethaneAIR_L1B_O2_20210803T155219_20210803T155249_20210804T101539.nc')
#slave_dir.append('/Users/asouri/Documents/Methane_SAT_OSSEs/Main/GEOAkaze/GEOAkaze/data/MethaneAIR_rf03/MethaneAIR_L1B_O2_20210728T183212_20210728T183242_20210730T100819.nc')

master_f = []
master_f.append('/Volumes/My Passport/S5P/Landsat_Permian_rad.nc')

gkobj = GEOAkaze(slave_f, master_f, 0.0002, 0, 2, 10000, is_histeq=True, is_destriping=False, bandindex_slave=1, w1=0, w2=160)

gkobj.readslave()
#plt.imshow(gkobj.slave,aspect='auto')
#plt.show()

gkobj.readmaster()
#plt.imshow(gkobj.master,aspect='auto')
#plt.show()

#gkobj.append_master()

gkobj.akaze()

print(gkobj.slope_lon)
print(gkobj.intercept_lon)
print(gkobj.slope_lat)
print(gkobj.intercept_lat)
print(gkobj.success)

gkobj.savetokmz('./MethaneAIR_L1B_O2_20210803T155219.kmz')

gkobj.write_to_nc('./test.nc')

gkobj.savetotxt('./test.txt')