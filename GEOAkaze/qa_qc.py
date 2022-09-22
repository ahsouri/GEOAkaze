# UTF-8
# QA/QC for geoakaze outputs
# Amir Souri (ahsouri@cfa.harvard.edu;ahsouri@gmail.com)
import os.path
import glob
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.ticker import FormatStrFormatter
from fpdf import FPDF

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

    def gray_scale(self,random_selection_n=None):
        '''
        plotting gray_scale images with information from geoakaze_nc_fld
        '''
        geoakaze_diag_fnames = sorted(glob.glob(self.geoakaze_nc_fld + '/*.nc'))
        # perform the random sampling if it's needed.
        if (random_selection_n is not None):
            ind = np.random.randint(len(geoakaze_diag_fnames),size=random_selection_n)
            geoakaze_diag_fnames = geoakaze_diag_fnames [ind]
        # loop over files
        for fname in geoakaze_diag_fnames:
            mair_gscale = self.read_netcdf(fname,"slave_gray")
            mst_gscale  = self.read_netcdf(fname,"master_gray")
            mair_lat    = self.read_netcdf(fname,"lats_new")
            mair_lon    = self.read_netcdf(fname,"lons_new")
            success     = self.read_netcdf(fname,"success")
            # plate projection at the desired box
            pc = ccrs.PlateCarree()
            fig = plt.figure(figsize=(8, 8))
            ax = plt.axes(projection = pc)
            ax.set_extent([np.nanmin(mair_lon.flatten()), np.nanmax(mair_lon.flatten()),
                           np.nanmin(mair_lat.flatten()), np.nanmax(mair_lat.flatten())], crs = pc)
            
            # plotting mair
            ax.imshow(mair_gscale,origin='lower',
               extent = [np.nanmin(mair_lon.flatten()), np.nanmax(mair_lon.flatten()),
                           np.nanmin(mair_lat.flatten()), np.nanmax(mair_lat.flatten())],
               interpolation='nearest',aspect='auto')
            
            # plotting msi
            ax.imshow(mair_gscale,origin='lower',
               extent = [np.nanmin(mair_lon.flatten()), np.nanmax(mair_lon.flatten()),
                           np.nanmin(mair_lat.flatten()), np.nanmax(mair_lat.flatten())],
               interpolation='nearest',aspect='auto',alpha=0.4)            
            # plotting costlines
            ax.coastlines(resolution='50m', color='black', linewidth = 2)
            ax.add_feature(ccrs.cartopy.feature.STATES)

            # fixing tickers
            x_ticks = np.arange(np.nanmin(mair_lon.flatten()), np.nanmax(mair_lon.flatten()), 0.03)
            x_labels = np.linspace(np.nanmin(mair_lon.flatten()), np.nanmax(mair_lon.flatten()), np.size(x_ticks))
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_labels, fontsize = 13)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

            y_ticks = np.arange(np.nanmin(mair_lat.flatten()), np.nanmax(mair_lat.flatten()), 0.03)
            y_labels = np.linspace(np.nanmin(mair_lat.flatten()), np.nanmax(mair_lat.flatten()), np.size(y_ticks))
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_labels, fontsize = 13)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

            plt.xlabel('Lon',fontsize = 18)
            plt.ylabel('Lat',fontsize = 18)

            # plotting title and saving
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
        fig = plt.figure(figsize=(8, 8))

        for filename in l1_files:
            try:
               nc_f = filename
               nc_fid = Dataset(nc_f, 'r')
               lat = nc_fid.groups["Geolocation"].variables["Latitude"][:]
               lon = nc_fid.groups["Geolocation"].variables["Longitude"][:]
               nc_fid.close()
               plt.scatter(np.nanmean(lon,axis=1), np.nanmean(lat,axis=1), s=2, color="blue")
            except:
               print("the file is not readable")
        fig.savefig(self.temp_fld + "/trajectory.png", format='png', dpi=300)

    def topdf(self):

        def header(pdfobj,title):
            # Arial bold 15
            pdfobj.set_font('Arial', 'B', 18)
            # Calculate width of title and position
            w = pdfobj.get_string_width(title) + 6
            pdfobj.set_x((210 - w) / 2)
            pdfobj.set_fill_color(255, 255, 255)
            pdfobj.set_text_color(0, 0, 0)
            # Thickness of frame (1 mm)
            pdfobj.set_line_width(1)
            # Title
            pdfobj.cell(w, 9, title, 1, 1, 'C', 1)
            # Line break
            pdfobj.ln(10)

        pdf = FPDF()
        pdf.add_page()
        header(pdf,"KMZ files")
        pdf.output('tuto3.pdf', 'F')