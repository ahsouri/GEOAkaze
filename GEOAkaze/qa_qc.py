# UTF-8
# QA/QC for geoakaze outputs
# Amir Souri (ahsouri@cfa.harvard.edu;ahsouri@gmail.com)
import os.path
import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

class qa_geoakaze(object):

    def __init__(self,geoakaze_nc_fld:str,geoakaze_kmz_fld:str,L1b_folder:str,output_pdf:str,
                 temp_fld:str):
            '''
            Initializing QA/QC with
            ARGS: 
                geoakaze_nc_fld (str): the full path of geoakaze diag files
                geoakaze_kmz_fld (str): the full path of geoakaze kmz files
                L1b_folder (str): the full path of L1b files
                output_pdf (str): the full path of the report in pdf
                temp_fld (str): a temp folder for intermediary files
            '''

            self.geoakaze_nc_fld = geoakaze_nc_fld
            self.geoakaze_kmz_fld = geoakaze_kmz_fld
            self.L1b_folder = L1b_folder
            self.output_pdf = output_pdf
            self.temp_fld = temp_fld
    
    def read_netcdf(self,filename,var):
        ''' 
        Read nc format from a file without a group
        ARGS:
            filename (char): the name of file
            var (char): the target variable
        OUT:
            var (float)
        '''

        nc_f = filename
        nc_fid = Dataset(nc_f, 'r')
        var = nc_fid.variables[var][:]
        nc_fid.close()
        return np.squeeze(var)

    def gray_scale(self):
        '''
        plotting gray_scale images with information from geoakaze_nc_fld
        '''
        geoakaze_diag_fnames = sorted(glob.glob(self.geoakaze_nc_fld + '/*.nc'))
        # loop over files
        for fname in geoakaze_diag_fnames:
            mair_gscale = self.read_netcdf(fname,"slave_gray")
            mair_lat    = self.read_netcdf(fname,"lats_new")
            mair_lon    = self.read_netcdf(fname,"lons_new")
            success     = self.read_netcdf(fname,"success")
            # plate projection at the desired box
            pc = ccrs.PlateCarree()
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(1, 1, 1, projection = pc)
            ax.set_extent([np.nanmin(mair_lon.flatten()), np.nanmax(mair_lon.flatten()),
                           np.nanmin(mair_lat.flatten()), np.nanmax(mair_lat.flatten())], crs = pc)
            # plotting mair
            ax.imshow(mair_gscale,origin='upper',
               extent = [np.nanmin(mair_lon.flatten()), np.nanmax(mair_lon.flatten()),
                           np.nanmin(mair_lat.flatten()), np.nanmax(mair_lat.flatten())],
               interpolation='none')
            # plotting costlines
            ax.coastlines(resolution='50m', color='black', linewidth = 2)
            ax.add_feature(ccrs.cartopy.feature.STATES)
            # plotting title
            fcolor = "green" if success == 1 else "red"
            plt.title('MAIR Grayscale for ' + str(os.path.basename(fname)), loc='left', color=fcolor, fontweight='bold', fontsize=12)
            fig.savefig(self.temp_fld + "/" + os.path.basename(fname) + "_grayscale.png", format='png', dpi=300)
            plt.close()
    def histogram(self):
        '''
        plotting a histogram of shifts
        '''
        geoakaze_diag_fnames = sorted(glob.glob(self.geoakaze_nc_fld + '/*.nc'))
        # loop over files
        hist_x = []
        for fname in geoakaze_diag_fnames:
            mair_lat        = self.read_netcdf(fname,"lats_new")
            mair_lon        = self.read_netcdf(fname,"lons_new")
            mair_lat_old    = self.read_netcdf(fname,"lats_old")
            mair_lon_old    = self.read_netcdf(fname,"lons_old")
            success         = self.read_netcdf(fname,"success")

            if success == 1:
               diff = (mair_lat-mair_lat_old)**2 + (mair_lon-mair_lon_old)**2 
               diff = np.sqrt(diff)
               diff = np.nanmean(diff.flatten())
               hist_x.append(diff*110*1000) # shift in meter
        #plotting       
        fig = plt.figure(figsize=(8, 8))
        plt.hist(hist_x, density=True, bins=30)  
        plt.ylabel('Probability')
        plt.xlabel('Error [m]')
        fig.savefig(self.temp_fld + "/histogram.png", format='png', dpi=300)

    def trajectory(self):

        l1_files = sorted(glob.glob(self.L1b_folder + "/*.nc"))

        lat_c = []
        lon_c = []
        for filename in l1_files:
            try:
               nc_f = filename
               nc_fid = Dataset(nc_f, 'r')
               lat = nc_fid.groups["Geolocation"].variables["Latitude"][:]
               lat_c.append(np.nanmean(lat,axis=1))
               lon = nc_fid.groups["Geolocation"].variables["Longitude"][:]
               lon_c.append(np.nanmean(lon,axis=1))
               nc_fid.close()
            except:
               print("file unreadable")
        fig = plt.figure(figsize=(8, 8))
        plt.scatter(lon_c, lat_c)
        fig.savefig(self.temp_fld + "/trajectory.png", format='png', dpi=300)