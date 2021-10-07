# UTF-8
# Apply the akaze algorithm on a satellite image with resepect to
# a reference image to rectify the geolocation errors in the first image
# Amir Souri (ahsouri@cfa.harvard.edu;ahsouri@gmail.com)
# July 2021


class GEOAkaze(object):

    def __init__(self,slavefile,masterfile,gridsize,typesat_slave,typesat_master,dist_thr,
                   msi_clim_fld=None,is_histeq=True,is_destriping=False,bandindex_slave=1,
                   bandindex_master=None,w1=None,w2=None,w3=None,w4=None):
            
            import os.path
            import glob
            '''
            Initializing GEOAkaze with the primary inputs
            ARGS: 
                slavefile (char): the name or the folder for the target satellite
                masterfile (char): the name or the folder for the reference image
                gridsize (float): grid size of mosaicing in degree unit
                is_histeq (bool): applying an adaptive histogram equalization
                typesat_slave and typesat_master (int):     
                      0: MethaneAIR
                      1: MethaneSAT_OSSE(nc) 
                      2: Landsat(nc)
                      3: MSI(jp2)
                      4: MSI(nc)

                dist_thr (int): a threshold used for filtering bad matches
                is_destriping (bool): whether to remove strips 
                bandindex_slave and bandindex_master (int): the index for reading bands in 
                                the netcdf file (1 = Band1)
                w1,w2 (int): boundaries for wavelength index of radiance to be averaged. (slave) 
                w3,w4 (int): boundaries for wavelength index of radiance to be averaged. (master)              
            '''        

            # check if the slavefile is a folder or a file
            
            if os.path.isdir(os.path.abspath(slavefile[0])):
                # we need to make a mosaic
                self.is_slave_mosaic = True
                self.slave_bundle = sorted(glob.glob(slavefile[0] + '/*.nc'))
            else:
                self.slave_bundle = []
                if len(slavefile) > 1:
                   self.is_slave_mosaic = True
                   for fname in slavefile:
                       self.slave_bundle.append(os.path.abspath(fname))
                else:
                    self.is_slave_mosaic = False
                    self.slave_bundle = os.path.abspath(slavefile[0])
            
            # check if the masterfile is a folder or a file

            if os.path.isdir(os.path.abspath(masterfile[0])):
                # we need to make a mosaic
                self.is_master_mosaic = True
                self.master_bundle = sorted(glob.glob(masterfile[0] + '/*'))
            else:
                self.is_master_mosaic = False
                self.master_bundle = []
                if len(masterfile) > 1:
                   self.is_master_mosaic = True
                   for fname in masterfile:
                       self.master_bundle.append(os.path.abspath(fname))  
                else:   
                   self.master_bundle = os.path.abspath(masterfile[0])

            self.gridsize = gridsize
            self.is_histeq = is_histeq 
            self.typesat_slave = typesat_slave
            self.typesat_master = typesat_master
            self.bandindex_slave = bandindex_slave
            self.bandindex_master = bandindex_master
            self.w1 = w1
            self.w2 = w2
            self.w3 = w3
            self.w4 = w4
            self.dist_thr = dist_thr
            self.is_destriping = is_destriping
            self.intercept_lat = 0.0
            self.slope_lat = 1.0
            self.intercept_lon = 0.0
            self.slope_lon = 1.0
            self.success = 0
            self.msi_clim_fld = msi_clim_fld        

    def read_netcdf(self,filename,var):
        ''' 
        Read nc format from a file without a group
        ARGS:
            filename (char): the name of file
            var (char): the target variable
        OUT:
            var (float)
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
            var (float)
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

    def read_rad(self,fname,typesat,bandindex=None,w1=None,w2=None):
        '''
        Read the intensity for differrent files/satellites
        ARGS:
            fname (char): the name of the file
            typesat = 0: MethaneAIR
                      1: MethaneSAT_OSSE(nc)
                      2: Landsat(nc)
                      3: MSI(jp2)
                      4: MSI(nc)

            bandindex (int): the index of band (e.g., =1 for O2)
            w1,w2 (int): the range of wavelength indices for averaging
        OUT:
            radiance, latitude, longitude
        '''
        import numpy as np
        
        if typesat == 0 or typesat == 1:
           rad = self.read_group_nc(fname,1,'Band' + str(bandindex),'Radiance')[:]

           # get flags used for labeling ortho settings
           av_used = self.read_group_nc(fname,1,'SupportingData','AvionicsUsed')
           ak_used= self.read_group_nc(fname,1,'SupportingData','AkazeUsed')
           op_used = self.read_group_nc(fname,1,'SupportingData','OptimizedUsed')

           if (op_used == 0 and ak_used == 0 and av_used == 1):
              lat = self.read_group_nc(fname,1,'Geolocation','Latitude')[:]
              lon = self.read_group_nc(fname,1,'Geolocation','Longitude')[:]
           else:
              lat = self.read_group_nc(fname,1,'SupportingData','AvionicsLatitude')[:]
              lon = self.read_group_nc(fname,1,'SupportingData','AvionicsLongitude')[:]

           rad [rad <= 0] = np.nan
           if not (w1 is None): #w1 and w2 should be set or none of them
               rad = np.nanmean(rad[w1:w2,:,:],axis=0)
           else:
               rad = np.nanmean(rad[:,:,:],axis=0)
        elif typesat == 2:
            rad = self.read_netcdf(fname,'Landsat')
            lat = self.read_netcdf(fname,'Lat')
            lon = self.read_netcdf(fname,'Lon')
        elif typesat == 4:
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

        if self.typesat_slave == 0:
            # read the data
            date_slave = []
            if self.is_slave_mosaic:
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
                   r,la,lo = self.read_rad(fname,self.typesat_slave,self.bandindex_slave,self.w1,self.w2)
                   if self.is_destriping: la = self.destriping(la)
                   rad.append(r)
                   lats.append(la)
                   lons.append(lo)

               date_slave = np.array(date_slave)
               self.yyyymmdd = np.median(date_slave)
               # make a mosaic
               mosaic,self.lats_grid,self.lons_grid,self.maskslave = self.mosaicing(rad,lats,lons)
            else:
                fname = self.slave_bundle
                print(fname)
                date_tmp = fname.split("_")
                date_tmp = date_tmp[-3]
                date_tmp = date_tmp.split("T")
                date_tmp = date_tmp[0]
                date_slave.append(float(date_tmp))
                r,la,lo = self.read_rad(fname,self.typesat_slave,self.bandindex_slave,self.w1,self.w2)
                if self.is_destriping: la = self.destriping(la)
                date_slave = np.array(date_slave)
                self.yyyymmdd = np.median(date_slave)
                # make a mosaic
                mosaic,self.lats_grid,self.lons_grid,self.maskslave = self.mosaicing(r,la,lo)
           

        elif self.typesat_slave == 2 or self.typesat_slave == 3: #landsat or MSI
            r,la,lo = self.read_rad(self.slave_bundle,self.typesat_slave)

        # normalizing
        self.slave = cv2.normalize(mosaic,np.zeros(mosaic.shape, np.double),0.0,1.0,cv2.NORM_MINMAX)
        
        if self.is_histeq:
           clahe = cv2.createCLAHE(clipLimit = 2.0, tileGridSize = (20,20))
           self.slave = clahe.apply(np.uint8(self.slave*255))
        else:
           self.slave = np.uint8(self.slave*255)
        
        # we will need it to append master img to L1 file
        if self.is_slave_mosaic:
           self.slavelat = lats[0]
           self.slavelon = lons[0]
        else:
           self.slavelat = la
           self.slavelon = lo

        self.rawslave =  r

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
               rad  = []
               lats = []
               lons = []
               for fname in self.master_bundle:
                   print(fname)
                   r,la,lo = self.read_rad(fname,self.typesat_master,self.bandindex_master,self.w3,self.w4)
                   if self.is_destriping: la = self.destriping(la)
                   rad.append(r)
                   lats.append(la)
                   lons.append(lo)

               r,lats,lons, _ = self.mosaicing(rad,lats,lons)
               r = self.cutter(r,lats,lons)
               self.rawmaster = r
            else:
                fname = self.master_bundle
                print(fname)
                r,la,lo = self.read_rad(fname,self.typesat_master,self.bandindex_master,self.w3,self.w4)
                if self.is_destriping: la = self.destriping(la)
                # make a mosaic
                r,lats,lons,_ = self.mosaicing(r,la,lo)
                r = self.cutter(r,lats,lons)
                self.rawmaster = r

        elif self.typesat_master == 2 or self.typesat_master == 4: #landsat
            r,la,lo = self.read_rad(self.master_bundle,self.typesat_master)
            r = self.cutter(r,la,lo)
            self.rawmaster = r
        elif self.typesat_master == 3: #MSI jp2
            r,la,lo  = self.read_MSI(self.master_bundle)
            r = self.cutter(r,la,lo)
            self.rawmaster = r/10000.0
        elif self.typesat_master == 5:
            r,la,lo  = self.read_rad(self.master_bundle,self.typesat_master)
            r = self.cutter(r,la,lo)
            self.rawmaster = r        
        
         # normalizing
        self.master = cv2.normalize(r,np.zeros(r.shape, np.double),0.0,1.0,cv2.NORM_MINMAX) 

        if self.is_histeq:
           clahe = cv2.createCLAHE(clipLimit = 2.0, tileGridSize=(20,20))
           self.master = clahe.apply(np.uint8(self.master*255))
        else:
           self.master = np.uint8(self.master*255)

    def mosaicing(self,rads,lats,lons):
        '''
        Merge (mosaic) several images together based on 
        their latitude/longitude. The final box is made of
        min/max of laitude and longitude based on all data
        ARGS:
            rads (list, floats): list of radiance arrays
            lons, lats (list, floats): list of longitude/latitude arrays
        OUT:
            mosaic, gridded_lat, gridded_lon
        ''' 
        import numpy as np      
        from scipy.spatial import Delaunay
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

        lons_grid,lats_grid = np.meshgrid(lon,lat)

        check_list = isinstance(rads, list)

        if check_list:
           full_moasic = np.zeros((np.shape(lons_grid)[0],np.shape(lons_grid)[1],len(rads)))
        
           for i in range(len(rads)):
               points = np.zeros((np.size(lons[i]),2))
               points[:,0] = np.array(lons[i]).flatten()
               points[:,1] = np.array(lats[i]).flatten()
               tri = Delaunay(points)
               interpolator = LinearNDInterpolator(tri,rads[i].flatten())
               full_moasic[:,:,i] = interpolator(lons_grid,lats_grid)
            # averaging
           full_moasic[full_moasic<=0] = np.nan
           mosaic = np.nanmean(full_moasic,axis=2)
           maskslave = np.isnan(mosaic)

        else:
            points = np.zeros((np.size(lons),2))
            points[:,0] = np.array(lons).flatten()
            points[:,1] = np.array(lats).flatten()
            tri = Delaunay(points)
            interpolator = LinearNDInterpolator(tri,rads.flatten())
            mosaic = interpolator(lons_grid, lats_grid)
            mosaic[mosaic<=0] = np.nan
            maskslave = np.isnan(mosaic)

        return mosaic,lats_grid,lons_grid,maskslave

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
        rad[self.maskslave] = np.nan
        
        return rad

    def akaze(self):
        '''
        AKAZE algorihtm : Pablo F Alcantarilla, Jesús Nuevo, and Adrien Bartoli. 
               Fast explicit diffusion for accelerated features in nonlinear scale spaces. Trans. 
               Pattern Anal. Machine Intell, 34(7):1281–1298, 2011.
        OUT:
            slope_lon,slope_lat,intercept_lon,intercept_lat (float) : correction factors
            success (0 or 1): 0->failed, 1->succeeded
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

        print('potential number of matched points: ' + str(len(master_matched)))

        self.nmatched = len(master_matched)
        self.matched_points_length = len(master_matched)

        data = np.column_stack([lat_1, lat_2])
        
        good_lat1, good_lat2 = self.robust_inliner(data)
        self.slope_lat, self.intercept_lat, self.r_value1, p_value, std_err = stats.linregress(good_lat1,good_lat2)
    
        data = np.column_stack([lon_1, lon_2])

        good_lon1, good_lon2 = self.robust_inliner(data)
        self.slope_lon, self.intercept_lon, self.r_value2, p_value, std_err = stats.linregress(good_lon1,good_lon2)
        
        # this part will be replaced by a cross-validation method
        if (abs(self.slope_lat)<0.9 or abs(self.slope_lat)>1.1 or abs(self.slope_lon)<0.9 or abs(self.slope_lon)>1.1 or
              self.r_value2<0.95 or self.r_value1<0.95):
           self.success = 0
        else:
           self.success = 1

    def find_matched_i_j(self,matches_var,keypoints1,keypoints2,dist_thr):
        '''
         A converter to transform the akaze objects to indices
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
           inliners [x,y] (float)
        ''' 
        # Fit line using all data
        from skimage.measure import LineModelND, ransac
        import numpy as np
        import matplotlib.pyplot as plt
        
        model = LineModelND()
        try:
           model.estimate(data)
        except:
           print('not enough matched points to work with')
           self.success = 0
        # Robustly fit linear model with RANSAC algorithm
        try:
            model_robust, inliers = ransac(data, LineModelND, min_samples=10, residual_threshold=0.0005,
                                    max_trials=100000)
        except:
            print('ransac cannot find outliers, failed!')
            self.success = 0
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
        from shapely.geometry import Polygon

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
            temp =  out_trans * (height,width)
            corner4 = np.array(utm.to_latlon(temp[0],temp[1],int(zones),'T') )

            p_master = Polygon([(corner1[1],corner4[0]), (corner4[1],corner4[0]), (corner4[1],corner1[0]), 
                             (corner1[1],corner1[0]), (corner1[1],corner4[0])])

            p_slave = Polygon([(np.min(np.min(self.lons_grid)),np.min(np.min(self.lats_grid))), 
                          (np.max(np.max(self.lons_grid)),np.min(np.min(self.lats_grid))),
                          (np.max(np.max(self.lons_grid)),np.max(np.max(self.lats_grid))), 
                          (np.min(np.min(self.lons_grid)),np.max(np.max(self.lats_grid))),
                          (np.min(np.min(self.lons_grid)),np.min(np.min(self.lats_grid)))])

            dist_cent = (np.array(p_master.centroid.coords)) - (np.array(p_slave.centroid.coords))
            dist_cent = dist_cent**2
            dist_cent = np.sum(dist_cent)
            dist_cent = np.sqrt(dist_cent)
            
            if  (p_master.contains(p_slave)) and (dist_cent<0.45):
                    within_box.append(fname)
                    date_tmp = fname.split("_")
                    date_tmp = date_tmp[-2]
                    date_tmp = date_tmp.split("T")
                    date_tmp = float(date_tmp[0])
                    msi_date.append(date_tmp)
        
        if not within_box:
            print('No MSI files being relevant to the targeted location/time were found, please fetch more MSI data')
            print('trying the climatological files now!')
            msi_gray,lat_msi,lon_msi = self.read_gee_tiff()
            return msi_gray,lat_msi,lon_msi

        dist_date = np.abs(np.array(msi_date) - float(self.yyyymmdd))
        index_chosen_one = np.argmin(dist_date)
        # now read the most relevant picture
        print('The chosen MSI is ' +  within_box[index_chosen_one])

        src = rasterio.open(within_box[index_chosen_one],driver='JP2OpenJPEG')
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

    def read_gee_tiff(self):
        '''
        MSI reader
        Default ARGS:
           geefname (str): msi_clim_fld
        OUT:
           msi_gray (float, array): the grayscale image of MSI
           lat_msi, lon_msi (floats): longitudes/latitudes of MSI
        ''' 
        # importing libraries
        from numpy import dtype
        import glob
        import numpy as np
        import rasterio
        from rasterio.merge import merge
        from shapely.geometry import Polygon
        import matplotlib.pyplot as plt
        
        within_box = []
        intersect_box = []
        geefname = sorted(glob.glob(self.msi_clim_fld + '/*.tif'))

        for fname in geefname:

            try:
                src = rasterio.open(fname,crs='EPSG:3857')
            except:
                continue
            out_trans = src.transform
            # check the boundaries
            width = src.width
            height = src.height

            corner1 =  out_trans * (0,0)     
            corner4 =  out_trans * (height,width)

            p_master = Polygon([(corner1[0],corner4[1]), (corner4[0],corner4[1]), (corner4[0],corner1[1]), 
                             (corner1[0],corner1[1]), (corner1[0],corner4[1])])

            p_slave = Polygon([(np.min(np.min(self.lons_grid)),np.min(np.min(self.lats_grid))), 
                          (np.max(np.max(self.lons_grid)),np.min(np.min(self.lats_grid))),
                          (np.max(np.max(self.lons_grid)),np.max(np.max(self.lats_grid))), 
                          (np.min(np.min(self.lons_grid)),np.max(np.max(self.lats_grid))),
                          (np.min(np.min(self.lons_grid)),np.min(np.min(self.lats_grid)))])
            
            if (p_master.contains(p_slave)):
                    within_box.append(fname)
            elif (p_master.intersects(p_slave)):
                    intersect_box.append(fname)
        
        if ((not within_box) and (not intersect_box)):
            print('The climatology MSI data do not cover this area')
            self.success = 0
            msi_gray = self.slave * 0.0
            lat_msi = self.lats_grid
            lon_msi = self.lons_grid
            return msi_gray,lat_msi,lon_msi

        # now read the most relevant picture
        # if there is no one single master to fully enclose the slave
        if (not within_box) and (intersect_box):
            src_appended = []
            for int_box in range(len(intersect_box)):
                src = rasterio.open(intersect_box[int_box],crs='EPSG:3857')
                src_appended.append(src)
            msi_img, out_trans = rasterio.merge(src_appended)
            print('Two tiles are chosen from the clim ' +  intersect_box)
        # if there is at least one master to fully enclose the slave
        elif within_box:           
            print('The chosen clim MSI is ' +  within_box[0])
            src = rasterio.open(within_box[0],crs='EPSG:3857')
            out_trans = src.transform
            msi_img = src.read(1)

        lat_msi = np.zeros_like(msi_img)*np.nan
        lon_msi = np.zeros_like(msi_img)*np.nan
        for i in range(np.shape(lon_msi)[0]):
            for j in range(np.shape(lon_msi)[1]):
                temp = out_trans * (j,i)
                lon_msi[i,j] = temp[0] 
                lat_msi[i,j] = temp[1]

        lat_msi = np.float32(lat_msi)
        lon_msi = np.float32(lon_msi)
        msi_gray = np.array(msi_img)
                
        return msi_gray,lat_msi,lon_msi


    def write_to_nc(self,output_file):
        ''' 
        Write the final results to a netcdf (for presentation purposes)
        ARGS:
            output_file (char): the name of file to be outputted
        '''
        from netCDF4 import Dataset
        import numpy as np
        from numpy import dtype

        ncfile = Dataset(output_file,'w')
        # create the x and y dimensions.
        ncfile.createDimension('x',np.shape(self.slave)[0])
        ncfile.createDimension('y',np.shape(self.slave)[1])
        ncfile.createDimension('z',1)
        
        data1 = ncfile.createVariable('master_gray',dtype('uint8').char,('x','y'))
        data1[:,:] = self.master
        
        data2 = ncfile.createVariable('slave_gray',dtype('uint8').char,('x','y'))
        data2[:,:] = self.slave
        
        data3 = ncfile.createVariable('lats_old',dtype('float64').char,('x','y'))
        data3[:,:] = self.lats_grid
        
        data4 = ncfile.createVariable('lons_old',dtype('float64').char,('x','y'))
        data4[:,:] = self.lons_grid
        
        data5 = ncfile.createVariable('lats_new',dtype('float64').char,('x','y'))
        data5[:,:] = (self.lats_grid-self.intercept_lat)/self.slope_lat
        
        data6 = ncfile.createVariable('lons_new',dtype('float64').char,('x','y'))
        data6[:,:] = (self.lons_grid-self.intercept_lon)/self.slope_lon
        
        data7 = ncfile.createVariable('success','u1',('z'))
        data7[:] = self.success
        
        ncfile.close()

    def destriping(self,lat):
        ''' 
        in case of the presence of strips in latitude , we can use this
        method (sobel filter + 1D spline interpolation)
        ARGS:
            lat (array, float): latitude with strips
        OUT:
            lat_destriped (array,float): destriped latitude
        '''
        import cv2
        import numpy as np
        from scipy.interpolate import UnivariateSpline


        sobely = cv2.Sobel(lat,cv2.CV_64F,0,1,ksize=5)
        abs_sobel = np.absolute(sobely)
        
        mask = np.zeros_like(lat)
        mask[ abs_sobel>1.2*np.mean(abs_sobel.flatten())] = 1.0
        

        lat [ mask != 0 ] = np.nan
        lat_destriped = np.zeros_like(lat)

        for j in range(0,np.shape(lat)[1]):  
            i = np.arange(0,np.shape(lat)[0])
            lat_line = lat[:,j]
            i_masked = i[~np.isnan(lat_line)]
            lat_masked = lat_line[~np.isnan(lat_line)]

            spl = UnivariateSpline(i_masked, lat_masked)
            spl.set_smoothing_factor(0.5)
            
            lat_destriped[:,j] = spl(i)
    
        return lat_destriped

    def append_master(self):  
        ''' 
          append master image to the L1 data
        ''' 
        from netCDF4 import Dataset
        import numpy as np
        from scipy.interpolate import griddata 

        if(isinstance(self.slave_bundle,list) == True):
            ncfile = Dataset(self.slave_bundle[0],'a',format="NETCDF4")
        else:
            ncfile = Dataset(self.slave_bundle,'a',format="NETCDF4")

        try:
           ncgroup = ncfile.createGroup('SupportingData')
           data = ncgroup.createVariable('Reference_IMG',np.float64,('y','x'))  
        except:
            # already is there
           data = ncfile.groups['SupportingData'].variables['Reference_IMG']

        points = np.zeros((np.size(self.lats_grid),2))
        points[:,0] = self.lons_grid.flatten()
        points[:,1] = self.lats_grid.flatten()

        img_master = griddata(points, self.rawmaster.flatten(), (self.slavelon, self.slavelat), method='nearest')
 
        data[:,:] = img_master

        ncfile.close()

    def savetokmz(self,fname):
        ''' 
        saving the mosaic of slave to a kmz file
        ARGS:
            fname (char): the kmz file
        '''

        from .make_kml import make_kmz
 
        if self.success == 1:
           lats_grid_corrected = (self.lats_grid-self.intercept_lat)/self.slope_lat
           lons_grid_corrected = (self.lons_grid-self.intercept_lon)/self.slope_lon
        else:
           lats_grid_corrected = self.lats_grid
           lons_grid_corrected = self.lons_grid
    
        make_kmz(lons_grid_corrected,lats_grid_corrected,self.slave,fname)
    
    def savetotxt(self,fname):
        ''' 
        saving the correction factors to a txt file
        ARGS:
            fname (char): the prefix part of  the txt file.
        '''
        import os.path

        filename = str(fname) + '_correction_factors_akaze.txt'
        
        if os.path.isfile(filename):
           file1 = open(filename, "a")
           L =  str(self.slope_lon) +',' + str(self.slope_lat) + ',' +str(self.intercept_lon) + \
                    ',' + str(self.intercept_lat) + ',' + str(self.r_value1) +',' + \
                        str(self.r_value2) + ',' + str(self.success)
           file1.writelines(L)
        else:
           L1 = 'file_bundle,slope_lon,slope_lat,intercept_lon,intercept_lat,rvalue_lon,rvalue_lat,success'
           L2 = str(self.slope_lon) +',' + str(self.slope_lat) + ',' +str(self.intercept_lon) + \
                    ',' + str(self.intercept_lat) + ',' + str(self.r_value1) +',' + \
                        str(self.r_value2) + ',' + str(self.success)

           file1 = open(filename, "w")
           #file1.writelines(L1)
           file1.writelines(L2)

    def hammer(self,slave_f,master_f1=None,master_f2=None,factor1=None,factor2=None):
        ''' 
        fixing the failed case (slave_f) using previous/next 
        or both successful cases
        ARGS:
            fname (slave_f): slave file path
            fname (master_f1): master file path (previous)
            fname (master_f2): master file path (next)
        '''
        from scipy import stats
        import numpy as np
        import matplotlib.pyplot as plt

        # read the slave and master
        _ ,lat_sl,lon_sl = self.read_rad(slave_f,0,1)
        if not (master_f1 is None):
           master_rad_1,lat_m1,lo_m1 = self.read_rad(master_f1,0,1)
           factors = np.loadtxt(factor1, delimiter=',')
           print(factors)
           print(np.shape(factors))
           lat_m1 = (lat_m1 - factors[3])/factors[1]
           lon_m1 = (lo_m1 - factors[2])/factors[0] 
        if not (master_f2 is None):
           master_rad_2,lat_m2,lo_m2 = self.read_rad(master_f2,0,1)
           factors = np.loadtxt(factor2, delimiter=',')
           lat_m2 = (lat_m1 - factors[3])/factors[1]
           lon_m2 = (lon_m1 - factors[2])/factors[0] 

        if ~(master_f1 is None) and ~(master_f2 is None): #only previous master is supplied
            
            #find the indices of non-nan gray scales
            saw_first_nan = False
            for i in range(0,np.shape(master_rad_1)[1]):
                if ~np.isnan(master_rad_1[-1,i]):
                    ind1 = i
                    saw_first_nan = True
                if (saw_first_nan) and np.isnan(master_rad_1[-1,i]):
                    ind2 = i - 1    
                
            pts1_m1 = np.zeros((ind2-ind1+1,2))
            pts2_m1 = np.zeros((ind2-ind1+1,2))

            pts1_m1[:,0] = lon_m1[-1,ind1:ind2+1]
            pts1_m1[:,1] = lat_m1[-1,ind1:ind2+1]
            pts2_m1[:,0] = lon_sl[0,ind1:ind2+1]
            pts2_m1[:,1] = lat_sl[0,ind1:ind2+1]

                        #find the indices of non-nan gray scales
            saw_first_nan = False
            for i in range(0,np.shape(master_rad_2)[1]):
                if ~np.isnan(master_rad_2[-1,i]):
                    ind1 = i
                    saw_first_nan = True
                if (saw_first_nan) and np.isnan(master_rad_2[-1,i]):
                    ind2 = i - 1
                    
                
            pts1_m2 = np.zeros((ind2-ind1+1,2))
            pts2_m2 = np.zeros((ind2-ind1+1,2))

            pts1_m2[:,0] = lon_m2[0,ind1:ind2+1]
            pts1_m2[:,1] = lat_m2[0,ind1:ind2+1]
            pts2_m2[:,0] = lon_sl[-1,ind1:ind2+1]
            pts2_m2[:,1] = lat_sl[-1,ind1:ind2+1]

            data_master = np.concatenate([pts1_m1, pts1_m2])
            print(np.shape(data_master))
            data_slave = np.concatenate([pts2_m1, pts2_m2])
            print(np.shape(data_slave))
            self.slope_lat, self.intercept_lat, self.r_value1, \
                 p_value, std_err = stats.linregress(data_master[:,1],
                                                     data_slave[:,1])
            self.slope_lon, self.intercept_lon, self.r_value2, \
                 p_value, std_err = stats.linregress(data_master[:,0],
                                                     data_slave[:,0])

            self.success = 1
        

        

        
        

