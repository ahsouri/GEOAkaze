# UTF-8
# Apply an akaze algorithm on a satellite image with resepect to
# a reference image to rectify the geolocation errors in the first image
# Amir Souri (ahsouri@cfa.harvard.edu;ahsouri@gmail.com)
# July 2021


class GEOAkaze(object):

    def __init__(self,slavefile,masterfile,gridsize,typesat_slave,typesat_master,dist_thr,is_histeq=True,bandindex=1,w1=None,w2=None):
            import os.path
            import glob
            '''
            Initializing GEOAkaze with the primary inputs
            ARGS: 
                slavefile (char): the name or the folder for the target satellite
                masterfile (char): the name or the folder for the reference image
                gridsize (float): grid size of mosaicing
                is_histeq (bool): 
            '''        
            # check if the slavefile is directory or a folder
            
            if os.path.isdir(os.path.abspath(slavefile)):
                # we need to make a mosaic
                self.is_slave_mosaic = True
                self.slave_bundle = sorted(glob.glob(slavefile + '/*.nc'))
            else:
                self.is_slave_mosaic = False
                self.slave_file = os.path.abspath(slavefile)
            if os.path.isdir(os.path.abspath(masterfile)):
                 # we need to make a mosaic
                self.is_master_mosaic = True
                self.master_bundle = sorted(glob.glob(masterfile + '/*'))
            else:
                self.is_master_mosaic = False
                self.master_bundle = os.path.abspath(masterfile)

            self.gridsize = gridsize
            self.is_histeq = is_histeq 
            self.typesat_slave = typesat_slave
            self.typesat_master = typesat_master
            self.bandindex = bandindex
            self.w1 = w1
            self.w2 = w2
            self.dist_thr = dist_thr

    def read_netcdf(self,filename,var):
        ''' 
        Read nc format from a file without a group
        ARGS:
            filename (char): the name of file
            var (char): the target variable
        OUT:
            variable
        '''
        from netCDF4 import Dataset
        import numpy as np
        nc_f = filename
        nc_fid = Dataset(nc_f, 'r')
        var = nc_fid.variables[var][:]
        nc_fid.close()
        return np.squeeze(var)

    def read_group_nc(self,filename,num_groups,group,var):
        ''' 
        Read nc format from a file with up to 3 subgroups
        ARGS:
            filename (char): the name of file
            num_groups (int): number of groups in the file
            group [num_groups] (list char): the name of group
            var (char): the target variable
        OUT:
            variable
        '''
        from netCDF4 import Dataset
        import numpy as np

        nc_f = filename
        nc_fid = Dataset(nc_f, 'r')
        if num_groups   == 1:
           out = nc_fid.groups[group].variables[var][:]
        elif num_groups == 2:
           out = nc_fid.groups[group[0]].groups[group[1]].variables[var][:]
        elif num_groups == 3:
           out = nc_fid.groups[group[0]].groups[group[1]].groups[group[2]].variables[var][:]
        nc_fid.close()
        return np.squeeze(out) 

    def read_rad(self,fname,typesat):
        '''
        Read the intensity for differrent files/satellites
        ARGS:
            fname (char): the name of the file
            typesat = 0: MethaneAIR
                      1: MethaneSAT_OSSE(nc) (not implemented yet)
                      2: Landsat(nc)
                      3: MSI(jp2)
                      4: MSI(nc)
            bandindex (int): the index of band (e.g., =1 for O2)
            w1,w2 (int): the range of wavelength indices for averaging
        OUT:
            radiance, latitude, longitude
        '''
        import numpy as np
        import cv2
        if typesat == 0:
           rad = self.read_group_nc(fname,1,'Band' + str(self.bandindex),'Radiance')[:]
           lat = self.read_group_nc(fname,1,'Geolocation','Latitude')[:]
           lon = self.read_group_nc(fname,1,'Geolocation','Longitude')[:]
           rad [rad <= 0] = np.nan
           if not (self.w1 is None): #w1 and w2 should be set or none of them
               rad = np.nanmean(rad[self.w1:self.w2,:,:],axis=0)
           else:
               rad = np.nanmean(rad[:,:,:],axis=0)
        elif typesat == 1:
            print('MethaneSAT_OSSE has not been not implemented yet!')
            exit()
        elif typesat == 2:
            rad = self.read_netcdf(fname,'Landsat')
            lat = self.read_netcdf(fname,'Lat')
            lon = self.read_netcdf(fname,'Lon')
        elif self.typesat == 4:
            rad = self.read_netcdf(fname,'MSI_clim')
            lat = self.read_netcdf(fname,'lat')
            lon = self.read_netcdf(fname,'lon')
        
        return rad,lat,lon

    def readslave(self):
        '''
        Read the slave (target) image for different satellite
        OUT:
            radiance, latitude, longitude
        '''
        import numpy as np
        import cv2

        date_slave = []
        if self.typesat_slave == 0:
            # read the data
            rad  = []
            lats = []
            lons = []
            for fname in self.slave_bundle:
                print(fname)
                date_tmp = fname.split("_")
                date_tmp = date_tmp[-3]
                date_tmp = date_tmp.split("T")
                date_tmp = date_tmp[0]
                date_slave.append(float(date_tmp))
                r,la,lo = self.read_rad(fname,self.typesat_slave)
                rad.append(r)
                lats.append(la)
                lons.append(lo)
            date_slave = np.array(date_slave)
            self.yyyymmdd = np.median(date_slave)
            # make a mosaic
            mosaic = self.mosaicing(rad,lats,lons)
            # normalizing
            self.slave = cv2.normalize(mosaic,np.zeros(mosaic.shape, np.double),1.0,0.0,cv2.NORM_MINMAX)
        elif self.typesat_slave == 2 or self.typesat_slave == 3: #landsat or MSI
            r,la,lo = self.read_rad(self.slave_bundle,self.typesat_slave)
            self.slave = cv2.normalize(r,np.zeros(r.shape, np.double),1.0,0.0,cv2.NORM_MINMAX)
        if self.is_histeq:
           clahe = cv2.createCLAHE(clipLimit = 2.0, tileGridSize = (20,20))
           self.slave = clahe.apply(np.uint8(self.slave*255))
        else:
           self.slave = np.uint8(self.slave*255)

        
    def readmaster(self): 
        '''
        Read the master (reference) image for different satellite
        OUT:
            radiance, latitude, longitude
        '''
        import numpy as np
        import cv2
        if self.typesat_master == 0:
            if self.is_master_mosaic:
               # read the data
               rad  = []
               lats = []
               lons = []
               for fname in self.master_bundle:
                   r,la,lo = rad.append(self.read_rad(fname,typesat,bandindex,w1,w2))
                   rad.append(r)
                   lats.append(lats)
                   lons.append(lons)
               # make a mosaic
               mosaic = self.mosaicing(rad,lats,lons)
               # normalizing
               self.master = cv2.normalize(mosaic,np.zeros(mosaic.shape, np.double),1.0,0.0,cv2.NORM_MINMAX)
            else:
                r,la,lo = self.read_rad(self.master_bundle,self.typesat_master)
                self.master = cv2.normalize(r,np.zeros(r.shape, np.double),1.0,0.0,cv2.NORM_MINMAX)
        elif self.typesat_master == 2 or self.typesat_master == 4: #landsat
            r,la,lo = self.read_rad(self.master_bundle,self.typesat_master)
            r = self.cutter(r,la,lo)
            self.master = cv2.normalize(r,np.zeros(r.shape, np.double),1.0,0.0,cv2.NORM_MINMAX)
        elif self.typesat_master == 3: #MSI jp2
            r,la,lo  = self.read_MSI(self.master_bundle)
            r = self.cutter(r,la,lo)
            self.master = cv2.normalize(r,np.zeros(r.shape, np.double),1.0,0.0,cv2.NORM_MINMAX)
        if self.is_histeq:
           clahe = cv2.createCLAHE(clipLimit =2.0, tileGridSize=(20,20))
           self.master = clahe.apply(np.uint8(self.master*255))
        else:
           self.master = np.uint8(self.master*255)

    def mosaicing(self,rads,lats,lons):
        '''
        Merge (mosaic) several images together based on 
        their latitude/longitude. The final box is made of
        min/max of laitude and longitude among all the data
        ARGS:
            rads (list, floats): list of radiance arrays
            lons, lats (list, floats): list of longitude/latitude arrays
        OUT:
            mosaic, gridded_lat, gridded_lon
        ''' 
        import numpy as np     
        from scipy.interpolate import griddata 
        from scipy import interpolate   
        from scipy.interpolate import RegularGridInterpolator
        from scipy.spatial import Delaunay
        import matplotlib.pyplot as plt
        from scipy.interpolate import LinearNDInterpolator
        # first making a mesh
        max_lat = []
        min_lat = []
        max_lon = []
        min_lon = []
        for i in range(len(rads)):
            min_lat.append(np.nanmin(lats[i]))
            max_lat.append(np.nanmax(lats[i]))
            min_lon.append(np.nanmin(lons[i]))
            max_lon.append(np.nanmax(lons[i]))
        min_lat = np.nanmin(min_lat)
        max_lat = np.nanmax(max_lat)
        min_lon = np.nanmin(min_lon)
        max_lon = np.nanmax(max_lon)
        lon = np.arange(min_lon,max_lon,self.gridsize)
        lat = np.arange(min_lat,max_lat,self.gridsize)
        self.lons_grid,self.lats_grid = np.meshgrid(lon,lat)
        full_moasic = np.zeros((np.shape(self.lons_grid)[0],np.shape(self.lons_grid)[1],len(rads)))

        for i in range(len(rads)):
            points = np.zeros((np.size(lons[i]),2))
            points[:,0] = np.array(lons[i]).flatten()
            points[:,1] = np.array(lats[i]).flatten()
            tri = Delaunay(points)
            interpolator = LinearNDInterpolator(tri,rads[i].flatten())
            full_moasic[:,:,i] = interpolator(self.lons_grid, self.lats_grid)
        
        # averaging
        full_moasic[full_moasic<=0] = np.nan
        self.masksleave = np.isnan(full_moasic)
        mosaic = np.nanmean(full_moasic,axis=2)

        return mosaic

    def cutter(self,rad,lat,lon):
        '''
        subset the large msi/landsat data based on min/max lons/lats
        ARGS:
           rad(float) : radiance
           lat(float) : latitude
           lon(float) : longitude
        OUT:
            mosaic, gridded_lat, gridded_lon
        ''' 
        import numpy as np
        from scipy.interpolate import griddata 
        lon_range = np.array([min(self.lons_grid.flatten()),max(self.lons_grid.flatten())])
        lat_range = np.array([min(self.lats_grid.flatten()),max(self.lats_grid.flatten())])
        

        mask_lon = (lon >= lon_range[0]) & (lon <= lon_range[1])
        mask_lat = (lat >= lat_range[0]) & (lat <= lat_range[1])
        
        rad = rad [ mask_lon & mask_lat ]
        lat = lat [ mask_lon & mask_lat ]
        lon = lon [ mask_lon & mask_lat ]

        points = np.zeros((np.size(lat),2))
        points[:,0] = lon.flatten()
        points[:,1] = lat.flatten()

        rad = griddata(points, rad.flatten(), (self.lons_grid, self.lats_grid), method='linear')
        rad[self.masksleave.squeeze()] = np.nan
        
        return rad

    def akaze(self):
        '''
        AKAZE algorihtm : Pablo F Alcantarilla, Jesús Nuevo, and Adrien Bartoli. 
               Fast explicit diffusion for accelerated features in nonlinear scale spaces. Trans. 
               Pattern Anal. Machine Intell, 34(7):1281–1298, 2011.
        OUT:
            slope_lon,slope_lat,intercept_lon,intercept_lat (float) : correction factors
            success (0 or 1): 0->failed, 1->succeed
        ''' 
        import cv2
        import numpy as np
        from scipy import stats

        akaze_mod = cv2.AKAZE_create()
        keypoints_1, descriptors_1 = akaze_mod.detectAndCompute(self.master,None)
        keypoints_2, descriptors_2 = akaze_mod.detectAndCompute(self.slave,None)
        bf = cv2.BFMatcher(cv2.DescriptorMatcher_BRUTEFORCE_HAMMING, crossCheck=True)
        matches = bf.match(descriptors_1,descriptors_2)
        matches = sorted(matches, key = lambda x:x.distance)
        master_matched,slave_matched = self.find_matched_i_j(matches,keypoints_1,keypoints_2,self.dist_thr)
        lat_1 = []
        lon_1 = []
        lat_2 = []
        lon_2 = []
        for i in range(np.shape(master_matched)[0]):
            lat_1.append(self.lats_grid[int(np.round(master_matched[i,1])),int(np.round(master_matched[i,0]))])
            lon_1.append(self.lons_grid[int(np.round(master_matched[i,1])),int(np.round(master_matched[i,0]))])
            lat_2.append(self.lats_grid[int(np.round(slave_matched[i,1])),int(np.round(slave_matched[i,0]))])
            lon_2.append(self.lons_grid[int(np.round(slave_matched[i,1])),int(np.round(slave_matched[i,0]))])

        lat_1 = np.array(lat_1)
        lat_2 = np.array(lat_2)
        lon_1 = np.array(lon_1)
        lon_2 = np.array(lon_2)

        pts1 = np.zeros((len(master_matched),2))
        pts2 = np.zeros((len(master_matched),2))

        pts1[:,0] = lon_1
        pts1[:,1] = lat_1
        pts2[:,0] = lon_2
        pts2[:,1] = lat_2

        print('number of matched points: ' + str(len(master_matched)))

        self.nmatched = len(master_matched)

        self.matched_points_length = len(master_matched)
        data = np.column_stack([lat_1, lat_2])
        good_lat1, good_lat2 = self.robust_inliner(data)
        self.slope_lat, self.intercept_lat, r_value1, p_value, std_err = stats.linregress(good_lat1,good_lat2)
    
        data = np.column_stack([lon_1, lon_2])
        good_lon1, good_lon2 = self.robust_inliner(data)
        self.slope_lon, self.intercept_lon, r_value2, p_value, std_err = stats.linregress(good_lon1,good_lon2)
        
        # this part will be replaced by a cross-validation method
        if (abs(self.slope_lat)<0.9 or abs(self.slope_lat)>1.1 or abs(self.slope_lon)<0.9 or abs(self.slope_lon)>1.1 or
              r_value2<0.95 or r_value1<0.95):
           self.success = 0
        else:
           self.success = 1

    def find_matched_i_j(self,matches_var,keypoints1,keypoints2,dist_thr):
        '''
         A converter to transfer akaze objects to indices
        '''
        import numpy as np
        # Initialize lists
        list_kp1 = []
        list_kp2 = []
        # For each match...
        for mat in matches_var:
        # Get the matching keypoints for each of the images
           if mat.distance>dist_thr:
              continue
           img1_idx = mat.queryIdx
           img2_idx = mat.trainIdx
           # x - columns
           # y - rows
           # Get the coordinates
           (x1, y1) = keypoints1[img1_idx].pt
           (x2, y2) = keypoints2[img2_idx].pt
           # Append to each list
           list_kp1.append((x1, y1))
           list_kp2.append((x2, y2))
        list_kp1 = np.array(list_kp1)
        list_kp2 = np.array(list_kp2)
        return list_kp1,list_kp2

    def robust_inliner(self,data):
        '''
        RANSAC algorithm: https://en.wikipedia.org/wiki/Random_sample_consensus
        ARGS:
           data array [x,y] (float)
        OUT:
           inliners
        ''' 
        # Fit line using all data
        from skimage.measure import LineModelND, ransac
        import numpy as np
        model = LineModelND()
        model.estimate(data)

        # Robustly fit linear model with RANSAC algorithm
        try:
            model_robust, inliers = ransac(data, LineModelND, min_samples=10, residual_threshold=0.0005,
                                    max_trials=100000)
        except:
            print('ransac cannot find outliers, failed!')
            self.success = 0
            exit()
        outliers = inliers == False
        # Predict data of estimated models
        line_x = np.arange(-360, 360)
        line_y = model.predict_y(line_x)
        line_y_robust = model_robust.predict_y(line_x)
        # Compare estimated coefficients
        doplot = False
        file_plot = './ransac_test.png'
        if doplot == True:
           fig, ax = plt.subplots()
           ax.plot(data[inliers, 0], data[inliers, 1], '.b', alpha=0.6,
             label='Inlier data')
           ax.plot(data[outliers, 0], data[outliers, 1], '.r', alpha=0.6,
             label='Outlier data')
           ax.plot(line_x, line_y, '-k', label='Line model from all data')
           ax.plot(line_x, line_y_robust, '-b', label='Robust line model')
           ax.legend(loc='lower left')
           plt.xlim(np.min(data[:,0])-0.01,np.max(data[:,0])+0.01)
           plt.ylim(np.min(data[:,1])-0.01,np.max(data[:,1])+0.01)
           plt.show()
           fig.savefig(file_plot + '.png',dpi=300)
           plt.close(fig)

        return data[inliers, 0],data[inliers, 1]

    def read_MSI(self,msifname):
        '''
        MSI reader
        ARGS:
           msifname (list, str): full address of jp2 files
        OUT:
           msi_gray (float, array): the grayscale image of MSI
           lat_msi, lon_msi (floats): longitudes/latitudes of MSI
        ''' 
        # importing libraries
        from numpy import dtype
        import numpy as np
        import rasterio
        import utm
        from rasterio.merge import merge

        within_box = []
        msi_date = []

        for fname in msifname:
            src = rasterio.open(fname,driver='JP2OpenJPEG')
            zones = (int(str(src.crs)[-2::]))
            out_trans = src.transform
            # check the boundaries
            width = src.width
            height = src.height

            temp =  out_trans * (0,0)
            corner1 = np.array(utm.to_latlon(temp[0],temp[1],int(zones),'T'))
            temp =  out_trans * (width,height)
            corner4 = np.array(utm.to_latlon(temp[0],temp[1],int(zones),'T') )
            if (np.max(np.max(self.lons_grid)) < corner4[1] and np.min(np.min(self.lons_grid)) > corner1[1] and  
                np.max(np.max(self.lats_grid)) < corner1[0] and np.min(np.min(self.lats_grid)) > corner4[0]):
                    within_box.append(fname)
                    date_tmp = fname.split("_")
                    date_tmp = date_tmp[-2]
                    date_tmp = date_tmp.split("T")
                    date_tmp = float(date_tmp[0])
                    msi_date.append(date_tmp)
        
        if not msi_date:
            print('No MSI files being relevant to the targeted location/time were found, please fetch more MSI data')
            self.success = 0
            exit()

        dist_date = np.abs(np.array(msi_date) - float(self.yyyymmdd))
        index_chosen_one = np.argmin(dist_date)
        # now read the most relevant picture
        
        src = rasterio.open(msifname[index_chosen_one],driver='JP2OpenJPEG')
        zones = (int(str(src.crs)[-2::]))
        out_trans = src.transform
        msi_img = src.read(1)
        E_msi = np.zeros_like(msi_img)*np.nan
        N_msi = np.zeros_like(msi_img)*np.nan
        for i in range(np.shape(E_msi)[0]):
            for j in range(np.shape(E_msi)[1]):
                temp = out_trans * (i,j)
                E_msi[i,j] = temp[0]
                N_msi[i,j] = temp[1]

        E_msi = np.float32(E_msi)
        N_msi = np.float32(N_msi)
        temp = np.array(utm.to_latlon(E_msi.flatten(),N_msi.flatten(),int(zones),'T'))
        temp2 = np.reshape(temp,(2,np.shape(msi_img)[0],np.shape(msi_img)[1]))
        lat_msi = np.squeeze(temp2[0,:,:])
        lon_msi = np.squeeze(temp2[1,:,:])
        msi_gray = np.array(msi_img, dtype='uint16').astype('float32')

        return np.transpose(msi_gray),lat_msi,lon_msi
