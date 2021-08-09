from GEOAkaze import GEOAkaze
import os.path
import matplotlib.pyplot as plt
import numpy as np
import glob
from joblib import Parallel, delayed

def wrapper_akaze(fname):

    slave_dir = []
    slave_dir.append(fname)


    master_f = '/scratch/sao_atmos/ahsouri/Retrieval/MethaneAIR/bundles/RF03_04/B11'

    gkobj = GEOAkaze(slave_dir, master_f, 0.0001, 0, 3, 10000, is_histeq=True, is_destriping=False, bandindex_slave=1, w1=1, w2=160)

    gkobj.readslave()

    gkobj.readmaster()

    gkobj.akaze()

    date_tmp = fname.split("_")
    date_tmp = date_tmp[-3]
    gkobj.write_to_nc('/scratch/sao_atmos/ahsouri/Retrieval/MethaneAIR/bundles/RF06_output/O2_output_' + date_tmp )


slave_dir = []
dir1 = '/scratch/sao_atmos/econway/RF06/O2Avionics/'
files = sorted(glob.glob(dir1 + '/*.nc'))
results = Parallel(n_jobs=24)(delayed(wrapper_akaze)(k) for k in files)


