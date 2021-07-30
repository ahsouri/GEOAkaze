from GEOAkaze import GEOAkaze
import os.path
import matplotlib.pyplot as plt


slave_dir = '/Users/asouri/Documents/Methane_SAT_OSSEs/Main/GEOAkaze/GEOAkaze/data/MethaneAIR_rf03/'
master_f = '/Volumes/My Passport 1/S5P/Landsat_rf03_rad.nc'
gkobj = GEOAkaze(slave_dir, master_f, 0.0001, 0, 2, 10000)

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