#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 12:14:08 2024

@author: andreamaccarinelli
"""




from matplotlib.gridspec import GridSpec
from astropy.io import fits
import numpy as np
import ipywidgets as widgets
import plotly.graph_objects as go
from IPython.display import display, HTML
from astropy.visualization import astropy_mpl_style
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
from astropy.wcs import WCS
#import pyregion
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
import random
from astropy.visualization.wcsaxes import WCSAxes
#from pyregion.mpl_helper import properties_func_default
from random import randrange
from astropy.nddata.utils import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.patches import Circle, Rectangle
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy.coordinates import SkyOffsetFrame, ICRS
#from regions import CircleSkyRegion, CircleAnnulusSkyRegion
#from scipy.ndimage.filters import gaussian_filter
import scipy
from matplotlib.patches import Ellipse, Circle
from os.path import exists
import os.path
from os import path
import pandas as pd
import glob
from astropy.cosmology import Planck15 as cosmo
from astropy.table import QTable
from astropy.table import Column
#from regions import Regions
import numpy as np
#from astroML.correlation import two_point
# from astroML.utils import pickle_results  # Update this line
#from astroML.datasets import fetch_sdss_specgals
#from astroML.correlation import bootstrap_two_point_angular
#from matplotlib.patches import PathPatch
#from matplotlib.path import Path



fract=[]
fract=np.zeros((1000))
x = [1,2,3,4,5,6,7,8]

for k in range(1000):
    
    fract[k] = len(np.where( (y > str(x) | (x> 0.04 ))[0])/len(tot))
    fract.append(len(np.where( (y > str(x)) | x> 0.04 )[0])/len(tot))
    
fractTOT=np.mean(fract)
fractERR=np.std(fract)