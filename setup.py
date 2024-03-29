from setuptools import setup,find_packages
from os.path import splitext
from os.path import basename
from glob import glob


with open('README.md') as f:
    readme = f.read()

setup(name='GEOAkaze',
      version='0.2.5',
      description='Geolocation correction using Akaze',
      long_description=readme,
      long_description_content_type='text/markdown',
      author='Amir Souri',
      author_email='ahsouri@gmail.com',
      license='MIT',
      packages=['GEOAkaze'],
      install_requires=[
          'numpy','matplotlib','scipy','netCDF4','opencv-python','rasterio','utm','shapely','simplekml','palettable','scikit-image',
          'fpdf', 'Cartopy'],
      zip_safe=False)
