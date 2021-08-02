from GEOAkaze import GEOAkaze
import os.path
import matplotlib.pyplot as plt
import numpy as np



slave_dir = '/Users/asouri/Documents/Methane_SAT_OSSEs/Main/GEOAkaze/GEOAkaze/data/MethaneAIR_rf03/'
master_f = '/Volumes/My Passport 1/S5P/msi_temp/band_11_jp2s/'

gkobj = GEOAkaze(slave_dir, master_f, 0.0001, 0, 3, 10000, is_histeq=True, bandindex=1, w1=0, w2=160)

gkobj.readslave()
plt.imshow(gkobj.slave,aspect='auto')
plt.show()
        

gkobj.readmaster()
plt.imshow(gkobj.master,aspect='auto')
plt.show()


gkobj.akaze()
print(gkobj.slope_lon)
print(gkobj.intercept_lon)
print(gkobj.slope_lat)
print(gkobj.intercept_lat)
print(gkobj.success)