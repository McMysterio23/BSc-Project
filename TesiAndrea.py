#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 17:42:46 2024

@author: andreatravascio
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
"""
from PyAstronomy import pyasl
from astropy.wcs import WCS
import pyregion
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
import random
from astropy.visualization.wcsaxes import WCSAxes
from pyregion.mpl_helper import properties_func_default
from random import randrange
from astropy.nddata.utils import Cutout2D
"""
from astropy import units as u
from astropy.coordinates import SkyCoord
"""
from matplotlib.patches import Circle, Rectangle
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy.coordinates import SkyOffsetFrame, ICRS
from regions import CircleSkyRegion, CircleAnnulusSkyRegion
"""

#from scipy.ndimage.filters import gaussian_filter
#import scipy
#from matplotlib.patches import Ellipse, Circle
#from os.path import exists
import os.path
#from os import path
#import pandas as pd
#import glob
#from astropy.cosmology import Planck15 as cosmo
#from astropy.table import QTable
#from astropy.table import Column
#from regions import Regions
#import numpy as np
#from astroML.correlation import two_point
#from astroML.utils import pickle_results  # Update this line
#from astroML.datasets import fetch_sdss_specgals
#from astroML.correlation import bootstrap_two_point_angular
#from matplotlib.patches import PathPatch
#from matplotlib.path import Path

#%% Selection BCG

#raBCG,decBCG = np.loadtxt("/Users/andreamaccarinelli/Desktop/SDSS/DATA_BCG.txt",usecols=[1,2],unpack=True,dtype=float)

sampleC4_file = "/Users/andreamaccarinelli/Desktop/SDSS/sampleC4.fits"
fit2file = fits.open(sampleC4_file)
data = fit2file[1].data
raBCG = data['RAdeg']
decBCG = data['DEdeg']

file="/Users/andreamaccarinelli/Desktop/SDSS/gal_info_dr7_v5_2.fits"
obj=fits.open(file)
data=obj[1].data
ra,dec=data['RA'],data['DEC']



indiciBCG=[]
#SELEZIONE BCG
for k in range(len(raBCG)):
    i = np.where((abs(ra-raBCG[k]) < 2/3600) & (abs(dec-decBCG[k]) < 2/3600))[0]
    if len(i) == 1:
        i=i[0]
        indiciBCG.append(i)
    else:
        if len(i) > 1:
            tt = np.sqrt((ra[i]-raBCG[k])**2 + (dec[i]-decBCG[k])**2)
            i=i[np.where(tt == np.min(tt))[0][0]]
            indiciBCG.append(i)
    
       
fileout="/Users/andreamaccarinelli/Desktop/BSc-Project/IndiciBCG.txt"
fmt = "%i"  # Specify the format string
data = np.column_stack((np.array(indiciBCG).astype(int),np.array(indiciBCG).astype(int)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


#%% Remove Index duplicates

file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/duplicates.txt"

def TxtRaws(file):
    dati = []
    # Apri il file e leggi le righe
    with open(file, 'r') as file:
        for riga in file:
            # Ignora le righe vuote o commenti
            if not riga.strip():
                continue
            # Dividi la riga in colonne e converte i valori in interi
            valori = [int(valore) for valore in riga.split()]
            dati.append(valori)

    # Converti la lista in un array numpy
    dati = np.array(dati)
    return dati
dupl=TxtRaws(file)


file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/IndiciBCG.txt"
ind = np.loadtxt(file,usecols=[0],unpack=True,dtype=int)


IndBCG=[]
indnoBCG=[]
indiciDupl=[]

for i in range(len(dupl)):
    if dupl[i][1] > 0 and (dupl[i][0] in indiciDupl) == False:
        for t in range(len(dupl[i])):
            indiciDupl.append(dupl[i][t])
        if dupl[i][0] in ind:
            IndBCG.append(dupl[i][0])
        else:
            indnoBCG.append(dupl[i][0])
    if dupl[i][1] < 0 and (dupl[i][0] in indiciDupl) == False:
        if dupl[i][0] in ind:
            IndBCG.append(dupl[i][0])
        else:
            indnoBCG.append(dupl[i][0])




fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Indici_BCG_nodupl.txt"
fmt = "%i"  # Specify the format string
ddd = np.column_stack((IndBCG,IndBCG))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, ddd, fmt=fmt)


fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Indici_noBCG_nodupl.txt"
fmt = "%i"  # Specify the format string
ddd = np.column_stack((indnoBCG,indnoBCG))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, ddd, fmt=fmt)







#%% Create file txt for only BCG and no BCG


file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/gal_info_dr7_v5_2.fits"
fitfile=fits.open(file)
data=fitfile[1].data
ra,dec,z,zerr,zwarn,vdisp,vdisperr,ebv=data['RA'],data['DEC'],data['Z'],data['Z_ERR'],data['Z_WARNING'],data['V_DISP'],data['V_DISP_ERR'],data['E_BV_SFD']

file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/gal_indx_dr7_v5_2.fits"
fitfile=fits.open(file)
data=fitfile[1].data
Zsun=data['BEST_MODEL_Z']


file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Indici_BCG_nodupl.txt"
ind = np.loadtxt(file,usecols=[0],unpack=True,dtype=int)
file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Indici_noBCG_nodupl.txt"
indno = np.loadtxt(file,usecols=[0],unpack=True,dtype=int)

raBCG,decBCG,zBCG,zerrBCG,zwarnBCG,vdispBCG,vdisperrBCG,ebvBCG,ZsunBCG=np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind)))
k2=0
for k in ind:
    raBCG[k2] = ra[k]
    decBCG[k2] = dec[k]
    zBCG[k2] = z[k]
    zerrBCG[k2] = zerr[k]
    zwarnBCG[k2] = zwarn[k]
    vdispBCG[k2] = vdisp[k]
    vdisperrBCG[k2] = vdisperr[k]
    ebvBCG[k2] = ebv[k]
    ZsunBCG[k2] = Zsun[k]
    k2 += 1
    




fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/InfoBCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(raBCG).astype(float),np.array(decBCG).astype(float),
                        np.array(zBCG).astype(float),np.array(zerrBCG).astype(float),
                        np.array(zwarnBCG).astype(float),np.array(vdispBCG).astype(float),
                        np.array(vdisperrBCG).astype(float),np.array(ebvBCG).astype(float),np.array(ZsunBCG).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)
    
raBCG,decBCG,zBCG,zerrBCG,zwarnBCG,vdispBCG,vdisperrBCG,ebvBCG,ZsunBCG=np.zeros((len(ra)-len(ind))),np.zeros((len(ra)-len(ind))),np.zeros((len(ra)-len(ind))),np.zeros((len(ra)-len(ind))),np.zeros((len(ra)-len(ind))),np.zeros((len(ra)-len(ind))),np.zeros((len(ra)-len(ind))),np.zeros((len(ra)-len(ind))),np.zeros((len(ra)-len(ind)))
k2=0
for k in indno:
    raBCG[k2] = ra[k]
    decBCG[k2] = dec[k]
    zBCG[k2] = z[k]
    zerrBCG[k2] = zerr[k]
    zwarnBCG[k2] = zwarn[k]
    vdispBCG[k2] = vdisp[k]
    vdisperrBCG[k2] = vdisperr[k]
    ebvBCG[k2] = ebv[k]
    ZsunBCG[k2] = Zsun[k]
    k2 += 1
        
    


fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/InfoNoBCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(raBCG).astype(float),np.array(decBCG).astype(float),
                        np.array(zBCG).astype(float),np.array(zerrBCG).astype(float),
                        np.array(zwarnBCG).astype(float),np.array(vdispBCG).astype(float),
                        np.array(vdisperrBCG).astype(float),np.array(ebvBCG).astype(float),np.array(ZsunBCG).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)
    

#%% Create file txt for only BCG and no BCG (line)


file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/gal_line_dr7_v5_2.fits"
fitfile=fits.open(file)
data=fitfile[1].data

sigB,esigB,sigF,esigF=data['SIGMA_BALMER'],data['SIGMA_BALMER_ERR'],data['SIGMA_FORBIDDEN'],data['SIGMA_FORBIDDEN_ERR']
vB,evB,vF,evF=data['V_OFF_BALMER'],data['V_OFF_BALMER_ERR'],data['V_OFF_FORBIDDEN'],data['V_OFF_FORBIDDEN_ERR']

#       OII_3726                 OII_3729                NEIII_3869              H_DELTA                  H_GAMMA               OIII_4363              OIII_4959                OIII_5007                HEI_5876                  OI_6300                 H_ALPHA                NII_6584                   SII_6717              SII_6731                ARIII7135
lines=["OII_3726","OII_3729","NEIII_3869","H_DELTA","H_GAMMA","OIII_4363","OIII_4959","OIII_5007","HEI_5876",
       "OI_6300","H_BETA","H_ALPHA","NII_6584","SII_6717","SII_6731","ARIII7135"]
FluxLines=np.zeros((len(lines),len(sigB)))
eFluxLines=np.zeros((len(lines),len(sigB)))
for i in range(len(lines)):
    FluxLines[i,:]=data[lines[i]+'_FLUX']
    eFluxLines[i,:]=data[lines[i]+'_FLUX_ERR']


file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Indici_BCG_nodupl.txt"
ind = np.loadtxt(file,usecols=[0],unpack=True,dtype=int)
file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Indici_noBCG_nodupl.txt"
indno = np.loadtxt(file,usecols=[0],unpack=True,dtype=int)

sigB2=np.zeros((len(ind)))
esigB2=np.zeros((len(ind)))
sigF2=np.zeros((len(ind)))
esigF2=np.zeros((len(ind)))
vB2=np.zeros((len(ind)))
evB2=np.zeros((len(ind)))
vF2=np.zeros((len(ind)))
evF2=np.zeros((len(ind)))
FluxLines2=np.zeros((len(lines),len(ind)))
eFluxLines2=np.zeros((len(lines),len(ind)))


k2=0
for k in ind:
    sigB2[k2] = sigB[k]
    esigB2[k2] = esigB[k]
    sigF2[k2] = sigF[k]
    esigF2[k2] = esigF[k]
    vB2[k2] = vB[k]
    evB2[k2] = evB[k]
    vF2[k2] = vF[k]
    evF2[k2] = evF[k]
    FluxLines2[:,k2] = FluxLines[:,k]*1E-17
    eFluxLines2[:,k2] = eFluxLines[:,k]*1E-17
    k2 += 1



fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/SigmaLines_BCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(sigB2).astype(float),np.array(esigB2).astype(float),np.array(sigF2).astype(float),np.array(esigF2).astype(float),
                        np.array(vB2).astype(float),np.array(evB2).astype(float),np.array(vF2).astype(float),np.array(evF2).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)
    
    


fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/FluxeLines_BCG.txt"
fmt = "%.5e"  # Specify the format string

data = np.column_stack((np.array(FluxLines2[0,:]).astype(float),np.array(eFluxLines2[0,:]).astype(float),
                        np.array(FluxLines2[1,:]).astype(float),np.array(eFluxLines2[1,:]).astype(float),
                        np.array(FluxLines2[2,:]).astype(float),np.array(eFluxLines2[2,:]).astype(float),
                        np.array(FluxLines2[3,:]).astype(float),np.array(eFluxLines2[3,:]).astype(float),
                        np.array(FluxLines2[4,:]).astype(float),np.array(eFluxLines2[4,:]).astype(float),
                        np.array(FluxLines2[5,:]).astype(float),np.array(eFluxLines2[5,:]).astype(float),
                        np.array(FluxLines2[6,:]).astype(float),np.array(eFluxLines2[6,:]).astype(float),
                        np.array(FluxLines2[7,:]).astype(float),np.array(eFluxLines2[7,:]).astype(float),
                        np.array(FluxLines2[8,:]).astype(float),np.array(eFluxLines2[8,:]).astype(float),
                        np.array(FluxLines2[9,:]).astype(float),np.array(eFluxLines2[9,:]).astype(float),
                        np.array(FluxLines2[10,:]).astype(float),np.array(eFluxLines2[10,:]).astype(float),
                        np.array(FluxLines2[11,:]).astype(float),np.array(eFluxLines2[11,:]).astype(float),
                        np.array(FluxLines2[12,:]).astype(float),np.array(eFluxLines2[12,:]).astype(float),
                        np.array(FluxLines2[13,:]).astype(float),np.array(eFluxLines2[13,:]).astype(float),
                        np.array(FluxLines2[14,:]).astype(float),np.array(eFluxLines2[14,:]).astype(float),
                        np.array(FluxLines2[15,:]).astype(float),np.array(eFluxLines2[15,:]).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)

##### NO BCG ######

sigB2=np.zeros((len(sigB)-len(ind)))
esigB2=np.zeros((len(sigB)-len(ind)))
sigF2=np.zeros((len(sigB)-len(ind)))
esigF2=np.zeros((len(sigB)-len(ind)))
vB2=np.zeros((len(sigB)-len(ind)))
evB2=np.zeros((len(sigB)-len(ind)))
vF2=np.zeros((len(sigB)-len(ind)))
evF2=np.zeros((len(sigB)-len(ind)))
FluxLines2=np.zeros((len(lines),len(sigB)-len(ind)))
eFluxLines2=np.zeros((len(lines),len(sigB)-len(ind)))


k2=0
for k in indno:
    sigB2[k2] = sigB[k]
    esigB2[k2] = esigB[k]
    sigF2[k2] = sigF[k]
    esigF2[k2] = esigF[k]
    vB2[k2] = vB[k]
    evB2[k2] = evB[k]
    vF2[k2] = vF[k]
    evF2[k2] = evF[k]
    FluxLines2[:,k2] = FluxLines[:,k]*1E-17
    eFluxLines2[:,k2] = eFluxLines[:,k]*1E-17
    k2 += 1


fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/SigmaLines_noBCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(sigB2).astype(float),np.array(esigB2).astype(float),np.array(sigF2).astype(float),np.array(esigF2).astype(float),
                        np.array(vB2).astype(float),np.array(evB2).astype(float),np.array(vF2).astype(float),np.array(evF2).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)

    


fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/FluxeLines_noBCG.txt"
fmt = "%.5e"  # Specify the format string

data = np.column_stack((np.array(FluxLines2[0,:]).astype(float),np.array(eFluxLines2[0,:]).astype(float),
                        np.array(FluxLines2[1,:]).astype(float),np.array(eFluxLines2[1,:]).astype(float),
                        np.array(FluxLines2[2,:]).astype(float),np.array(eFluxLines2[2,:]).astype(float),
                        np.array(FluxLines2[3,:]).astype(float),np.array(eFluxLines2[3,:]).astype(float),
                        np.array(FluxLines2[4,:]).astype(float),np.array(eFluxLines2[4,:]).astype(float),
                        np.array(FluxLines2[5,:]).astype(float),np.array(eFluxLines2[5,:]).astype(float),
                        np.array(FluxLines2[6,:]).astype(float),np.array(eFluxLines2[6,:]).astype(float),
                        np.array(FluxLines2[7,:]).astype(float),np.array(eFluxLines2[7,:]).astype(float),
                        np.array(FluxLines2[8,:]).astype(float),np.array(eFluxLines2[8,:]).astype(float),
                        np.array(FluxLines2[9,:]).astype(float),np.array(eFluxLines2[9,:]).astype(float),
                        np.array(FluxLines2[10,:]).astype(float),np.array(eFluxLines2[10,:]).astype(float),
                        np.array(FluxLines2[11,:]).astype(float),np.array(eFluxLines2[11,:]).astype(float),
                        np.array(FluxLines2[12,:]).astype(float),np.array(eFluxLines2[12,:]).astype(float),
                        np.array(FluxLines2[13,:]).astype(float),np.array(eFluxLines2[13,:]).astype(float),
                        np.array(FluxLines2[14,:]).astype(float),np.array(eFluxLines2[14,:]).astype(float),
                        np.array(FluxLines2[15,:]).astype(float),np.array(eFluxLines2[15,:]).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


#%% Create file txt for only BCG and no BCG (sSFR & Mass)


file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/gal_totsfr_dr7_v5_2.fits"
fitfile=fits.open(file)
data=fitfile[1].data
avg,med,p16,p84=data['AVG'],data['MEDIAN'],data['P16'],data['P84']

file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/gal_totspecsfr_dr7_v5_2.fits"
fitfile=fits.open(file)
data=fitfile[1].data
savg,smed,sp16,sp84=data['AVG'],data['MEDIAN'],data['P16'],data['P84']


file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/totlgm_dr7_v5_2.fits"
fitfile=fits.open(file)
data=fitfile[1].data
mavg,mmed,mp16,mp84=data['AVG'],data['MEDIAN'],data['P16'],data['P84']


Ntot=len(avg)

file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Indici_BCG_nodupl.txt"
ind = np.loadtxt(file,usecols=[0],unpack=True,dtype=int)
file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Indici_noBCG_nodupl.txt"
indno = np.loadtxt(file,usecols=[0],unpack=True,dtype=int)

avg2,med2,p162,p842=np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind)))
savg2,smed2,sp162,sp842=np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind)))
mavg2,mmed2,mp162,mp842=np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind))),np.zeros((len(ind)))

k2=0
for k in ind:
    avg2[k2] = avg[k]
    med2[k2] = med[k]
    p162[k2] = p16[k]
    p842[k2] = p84[k]
    savg2[k2] = savg[k]
    smed2[k2] = smed[k]
    sp162[k2] = sp16[k]
    sp842[k2] = sp84[k]
    mavg2[k2] = mavg[k]
    mmed2[k2] = mmed[k]
    mp162[k2] = mp16[k]
    mp842[k2] = mp84[k]
    k2 += 1



fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/MassSFR_BCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(mavg2).astype(float),np.array(mmed2).astype(float),
                        np.array(mp162).astype(float),np.array(mp842).astype(float),
                        np.array(avg2).astype(float),np.array(med2).astype(float),
                        np.array(p162).astype(float),np.array(p842).astype(float),
                        np.array(savg2).astype(float),np.array(smed2).astype(float),
                        np.array(sp162).astype(float),np.array(sp842).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


    
avg2,med2,p162,p842=np.zeros((Ntot-len(ind))),np.zeros((Ntot-len(ind))),np.zeros((Ntot-len(ind))),np.zeros((Ntot-len(ind)))
savg2,smed2,sp162,sp842=np.zeros((Ntot-len(ind))),np.zeros((Ntot-len(ind))),np.zeros((Ntot-len(ind))),np.zeros((Ntot-len(ind)))
mavg2,mmed2,mp162,mp842=np.zeros((Ntot-len(ind))),np.zeros((Ntot-len(ind))),np.zeros((Ntot-len(ind))),np.zeros((Ntot-len(ind)))


k2=0
for k in indno: 
    avg2[k2] = avg[k]
    med2[k2] = med[k]
    p162[k2] = p16[k]
    p842[k2] = p84[k]
    savg2[k2] = savg[k]
    smed2[k2] = smed[k]
    sp162[k2] = sp16[k]
    sp842[k2] = sp84[k]
    mavg2[k2] = mavg[k]
    mmed2[k2] = mmed[k]
    mp162[k2] = mp16[k]
    mp842[k2] = mp84[k]
    k2 += 1



fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/MassSFR_noBCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(mavg2).astype(float),np.array(mmed2).astype(float),
                        np.array(mp162).astype(float),np.array(mp842).astype(float),
                        np.array(avg2).astype(float),np.array(med2).astype(float),
                        np.array(p162).astype(float),np.array(p842).astype(float),
                        np.array(savg2).astype(float),np.array(smed2).astype(float),
                        np.array(sp162).astype(float),np.array(sp842).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)



#%% CALL DATA

Sampl="no" #"no" "y"
 

def calldata(Sampl):
    #INFO BCG
    if Sampl == 'no':
        f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/InfoNoBCG.txt"
    else:
        f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/InfoBCG.txt"
    RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun=np.loadtxt(f,usecols=[0,1,2,3,5,6,7,8],unpack=True,dtype=float)



    #Sigma Balmer and Forbidden Lines
    if Sampl == 'no':
        f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/SigmaLines_noBCG.txt"
    else:
        f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/SigmaLines_BCG.txt"
    SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB=np.loadtxt(f,usecols=[0,1,2,3],unpack=True,dtype=float)
    VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB=np.loadtxt(f,usecols=[4,5,6,7],unpack=True,dtype=float)


    #FLUXES LINES
    if Sampl == 'no':
        f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/FluxeLines_noBCG.txt"
    else:
        f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/FluxeLines_BCG.txt"
    OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869= np.loadtxt(f,usecols=[0,1,2,3,4,5],unpack=True,dtype=float)
    H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959= np.loadtxt(f,usecols=[6,7,8,9,10,11,12,13],unpack=True,dtype=float)
    OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300= np.loadtxt(f,usecols=[14,15,16,17,18,19],unpack=True,dtype=float)
    H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584= np.loadtxt(f,usecols=[20,21,22,23,24,25],unpack=True,dtype=float)
    SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135= np.loadtxt(f,usecols=[26,27,28,29,30,31],unpack=True,dtype=float)

    #Derived Prop. Mass SFR sSFR
    if Sampl == 'no':
        f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/MassSFR_noBCG.txt"
    else:
        f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/MassSFR_BCG.txt"
    Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2=np.loadtxt(f,usecols=[0,2,3,4,6,7,8,10,11],unpack=True,dtype=float)
    return RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2

RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl)


#%% Function PLOTS

from matplotlib.gridspec import GridSpec

def Histograms(dati_x, dati_y, colori, labels, assi_labels, bins=None, Normalized="y"):
    if bins is None:
        bins = [30, 30]

    # Creazione del plot principale
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(4, 4, figure=fig)

    # Plot principale
    ax_main = fig.add_subplot(gs[1:4, 0:3])
    for i in range(len(dati_x)):
        ax_main.scatter(dati_x[i], dati_y[i], color=colori[i], label=labels[i])
    # Istogramma sopra
    ax_top = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    for i in range(len(dati_x)):
        ax_top.hist(dati_x[i], bins=bins[0], color=colori[i], alpha=0.7, edgecolor='black', density=Normalized=='y')
    # ax_top.set_title('Istogramma su X')

    # Istogramma a destra
    ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax_main)
    for i in range(len(dati_y)):
        ax_right.hist(dati_y[i], bins=bins[1], orientation='horizontal', color=colori[i], alpha=0.7, edgecolor='black', density=Normalized=='y')
    # ax_right.set_title('Istogramma su Y')

    # Rimuovi etichette degli assi del subplot principale
    # ax_top.set_xticks([])
    # ax_right.set_yticks([])

    # Visualizza il plot principale
    ax_main.legend()
    ax_main.set_xlabel(assi_labels[0])
    ax_main.set_ylabel(assi_labels[1])
    plt.show()


def PlotScat(x,y,ex=None,ey=None,xlim=None,ylim=None,colore="black",simbolo="o",labels=["X","Y"],Positives=["yes","yes"],overplot=False):
    if (xlim is None) == True:
        pass
    else:
        if np.isnan(xlim[0]) == False and np.isnan(xlim[1]) == False:
            indexes=np.where((x >= xlim[0]) & (x <= xlim[1]))[0]
        else:
            if np.isnan(xlim[0]) == False:
                indexes=np.where(x <= xlim[1])[0]
            if np.isnan(xlim[1]) == False:
                indexes=np.where(x <= xlim[1])[0]
        x = x[indexes]
        y = y[indexes]
        if (ex is None) == False:
            ex = ex[indexes]
        if (ey is None) == False:
            ey = ey[indexes]
    if (ylim is None) == True:
        pass
    else:
        if np.isnan(ylim[0]) == False and np.isnan(ylim[1]) == False:
            indexes=np.where((y >= ylim[0]) & (y <= ylim[1]))[0]
        else:
            if np.isnan(ylim[0]) == False:
                indexes=np.where(y <= ylim[1])[0]
            if np.isnan(ylim[1]) == False:
                indexes=np.where(y <= ylim[1])[0]
        x = x[indexes]
        y = y[indexes]
        if (ex is None) == False:
            ex = ex[indexes]
        if (ey is None) == False:
            ey = ey[indexes]
    #Remove Large Errors
    if Positives[0] == "yes":
        if (ex is None) == False:
            indexes=np.where(x-ex > 0)[0]
            x = x[indexes]
            y = y[indexes]
            ex = ex[indexes]
            if (ey is None) == False:
                ey = ey[indexes]
    if Positives[1] == "yes":
        if (ey is None) == False:
            indexes=np.where(y-ey > 0)[0]
            x = x[indexes]
            y = y[indexes]
            if (ex is None) == False:
                ex = ex[indexes]
            ey = ey[indexes]
    if overplot == False:
        fig,ax=plt.subplots()
    plt.scatter(x,y,color=colore,linestyle='None', marker=simbolo)
    if (ex is None) == False and (ey is None) == False:
        plt.errorbar(x,y,xerr=ex,yerr=ey,color=colore,linestyle='None')
    if (ex is None) == False and (ey is None) == True:
        plt.errorbar(x,y,xerr=ex,color=colore,linestyle='None')
    if (ex is None) == True and (ey is None) == False:
        plt.errorbar(x,y,yerr=ey,color=colore,linestyle='None')
    plt.xlabel(labels[0],fontsize=16)
    plt.ylabel(labels[1],fontsize=16)
    plt.tick_params(axis='both', labelsize=16)
    return len(x)



def ErrLogRatio(num,den,err_num=None,err_den=None,Niter=1000):
    errors=np.zeros((len(num)))
    for i in range(len(num)):
        if (err_num[i] is None) == False:
            rndnum = np.random.normal(num[i], err_num[i], Niter)
            rndnum = np.clip(rndnum, num[i]- err_num[i], num[i]+ err_num[i])

        if (err_den[i] is None) == False:
            rndden = np.random.normal(den[i], err_den[i], Niter)
            rndden = np.clip(rndden, den[i]- err_den[i], den[i]+ err_den[i])
        logval=np.log10(rndnum/rndden)
        errors[i] = np.std(logval)
    return errors


#%% BPT DIAGRAMS 

Sampl="no"
RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl)


"""
OIIIHb = 0.61 / (NIIHa - 0.05) + 1.3     #(Kauffmann+03 line)
OIIIHb = 0.61 / (NIIHa - 0.47) + 1.19    #(Kewley+01 line)

OIIIHb = 0.72 / (SIIHa - 0.32) + 1.30    #(main AGN line)
OIIIHb = 1.89*SIIHa + 0.76               #(LINER/Sy2 line)

OIIIHb = 0.73 / (OIHa + 0.59) + 1.33     #(main AGN line)
OIIIHb = 1.18*OIHa + 1.30                 #(LINER/Sy2 line)
"""

#REMOVE FALSE VALUES
indicitot=np.arange(len(OIII_5007))
indexes1=np.where((OIII_5007 > 0) & (H_ALPHA > 0) & (H_BETA > 0) & (NII_6584 > 0) & (OIII_5007 > eOIII_5007) & (H_ALPHA > eH_ALPHA) & (H_BETA > eH_BETA) & (NII_6584 > eNII_6584))[0]
indexes2=np.where((OIII_5007 > 0) & (H_ALPHA > 0) & (H_BETA > 0) & (SII_6717 > 0) & (SII_6731 > 0) & (OIII_5007 > eOIII_5007) & (H_ALPHA > eH_ALPHA) & (H_BETA > eH_BETA) & (SII_6717 > eSII_6717) & (SII_6731 > eSII_6731))[0]
indexes3=np.where((OIII_5007 > 0) & (H_ALPHA > 0) & (H_BETA > 0) & (OI_6300 > 0) & (OIII_5007 > eOIII_5007) & (H_ALPHA > eH_ALPHA) & (H_BETA > eH_BETA) & (OI_6300 > eOI_6300))[0]

logOIIIHb1 = np.log10(OIII_5007[indexes1]/H_BETA[indexes1])
elogOIIIHb1 = ErrLogRatio(OIII_5007[indexes1],H_BETA[indexes1],err_num=eOIII_5007[indexes1],err_den=eH_BETA[indexes1])
logNIIHa1 = np.log10(NII_6584[indexes1]/H_ALPHA[indexes1])
elogNIIHa1 = ErrLogRatio(NII_6584[indexes1],H_ALPHA[indexes1],err_num=eNII_6584[indexes1],err_den=eH_ALPHA[indexes1])
i1= indicitot[indexes1]

logOIIIHb2 =np.log10(OIII_5007[indexes2]/H_BETA[indexes2])
elogOIIIHb2 = ErrLogRatio(OIII_5007[indexes2],H_BETA[indexes2],err_num=eOIII_5007[indexes2],err_den=eH_BETA[indexes2])
SII=SII_6717[indexes2]+SII_6731[indexes2]
eSII=eSII_6717[indexes2]+eSII_6731[indexes2]
logSIIHa2 = np.log10((SII)/H_ALPHA[indexes2])
elogSIIHa2 = ErrLogRatio(SII,H_ALPHA[indexes2],err_num=eSII,err_den=eH_ALPHA[indexes2])
i2= indicitot[indexes2]

logOIIIHb3 = np.log10(OIII_5007[indexes3]/H_BETA[indexes3])
elogOIIIHb3 = ErrLogRatio(OIII_5007[indexes3],H_BETA[indexes3],err_num=eOIII_5007[indexes3],err_den=eH_BETA[indexes3])
logOIHa3 = np.log10(OI_6300[indexes3]/H_ALPHA[indexes3])
elogOIHa3 = ErrLogRatio(OI_6300[indexes3],H_ALPHA[indexes3],err_num=eOI_6300[indexes3],err_den=eH_ALPHA[indexes3])
i3= indicitot[indexes3]

def BPTd():
    xx1=np.arange(-3,3,0.01)
    xx2=np.arange(-3,3,0.01)
    xx3=np.arange(-3,3,0.01)
    yy1_1 = (0.61/(xx1[xx1 < 0.05] - 0.05)) + 1.3
    yy1_2 = (0.61 / (xx1[xx1 < 0.47] - 0.47)) + 1.19

    yy2_1 = (0.72 / (xx2[xx2 < 0.32] - 0.32)) + 1.30    #(main AGN line)
    yy2_2 = 1.89*xx2[xx2 > -0.33] + 0.76         

    yy3_1 = (0.73 / (xx3[xx3 < -0.59] + 0.59)) + 1.33     #(main AGN line)
    yy3_2 = 1.18*xx3[xx3 > -1.1257] +  1.30 
    return xx1,xx2,xx3,yy1_1,yy1_2,yy2_1,yy2_2,yy3_1,yy3_2


def PBPT(n=1):
    xx1,xx2,xx3,yy1_1,yy1_2,yy2_1,yy2_2,yy3_1,yy3_2=BPTd()
    if n == 1:
        plt.plot(xx1[xx1 < 0.05],yy1_1,color='black')
        plt.plot(xx1[xx1 < 0.47],yy1_2,color='black')
        plt.xlim((-2,3))
        plt.ylim((-2,3))
    if n == 2:
        plt.plot(xx2[xx2 < 0.32],yy2_1,color='black')
        plt.plot(xx2[xx2 > -0.33],yy2_2,color='black')
        plt.xlim((-2,3))
        plt.ylim((-2,3))
    if n == 3:
        plt.plot(xx3[xx3 < -0.59],yy3_1,color='black')
        plt.plot(xx3[xx3 > -1.1257],yy3_2,color='black')
        plt.xlim((-2,3))
        plt.ylim((-2,3))
    plt.subplots_adjust(top=0.8, bottom=0.2, left=0.2, right=0.8, hspace=0.2, wspace=0.2)
    return 0
    

PlotScat(logNIIHa1,logOIIIHb1,ex=elogNIIHa1,ey=elogOIIIHb1,xlim=None,ylim=None,colore="red",simbolo="o",labels=["$log([NII]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"])
PBPT(n=1)

PlotScat(logSIIHa2,logOIIIHb2,ex=elogSIIHa2,ey=elogOIIIHb2,xlim=None,ylim=None,colore="blue",simbolo="o",labels=["$log([SII]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"])
PBPT(n=2)

PlotScat(logOIHa3,logOIIIHb3,ex=elogOIHa3,ey=elogOIIIHb3,xlim=None,ylim=None,colore="green",simbolo="o",labels=["$log([OI]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"])
PBPT(n=3)



if Sampl == 'no':
    fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-NII_gal.txt"
else:
    fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-NII.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(i1).astype(float),np.array(logNIIHa1).astype(float),np.array(elogNIIHa1).astype(float),
                        np.array(logOIIIHb1).astype(float),np.array(elogOIIIHb1).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)
    
if Sampl == 'no':
    fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-SII_gal.txt"
else:
    fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-SII.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(i2).astype(float),np.array(logSIIHa2).astype(float),np.array(elogSIIHa2).astype(float),
                        np.array(logOIIIHb2).astype(float),np.array(elogOIIIHb2).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)

if Sampl == 'no':
    fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-OI_gal.txt"
else:
    fileout="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-OI.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(i3).astype(float),np.array(logOIHa3).astype(float),np.array(elogOIHa3).astype(float),
                        np.array(logOIIIHb3).astype(float),np.array(elogOIIIHb3).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)



#%% SubSample Radiative, Shock and SF


Sampl="no"
RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl)



def SaveType(i,fileout,arrays):
    val=np.zeros((len(i),len(arrays)))
    kk=0
    for k in i:
        for t in range(len(arrays)):
            val[kk,t]=arrays[t][k]
        kk += 1
    if os.path.exists(fileout) == False:
        np.savetxt(fileout, val, delimiter='\t', header='\t'.join(map(str, range(len(arrays)))), comments='')
    return fileout
    

indicitot=np.arange(len(OIII_5007))
if Sampl == 'no':
    f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-NII_gal.txt"
else:
    f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-NII.txt"
i1,x1,ex1,y1,ey1=np.loadtxt(f,usecols=[0,1,2,3,4],unpack=True,dtype=float)

if Sampl == 'no':
    f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-SII_gal.txt"
else:
    f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-SII.txt"
i2,x2,ex2,y2,ey2=np.loadtxt(f,usecols=[0,1,2,3,4],unpack=True,dtype=float)

if Sampl == 'no':
    f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-OI_gal.txt"
else:
    f="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/BPT-OI.txt"
i3,x3,ex3,y3,ey3=np.loadtxt(f,usecols=[0,1,2,3,4],unpack=True,dtype=float)


#BPT NII
iAGN=[]
icomp=[]
ihii1=[]
xAGN=[]
exAGN=[]
yAGN=[]
eyAGN=[]
xcomp=[]
excomp=[]
ycomp=[]
eycomp=[]
xhii1=[]
exhii1=[]
yhii1=[]
eyhii1=[]

for k in range(len(i1)):
    if (y1[k] >= (0.61 / (x1[k] - 0.47)) + 1.19) or x1[k] >= 0.04: #RED
        iAGN.append(int(i1[k]))
        xAGN.append(x1[k])
        yAGN.append(y1[k])
        exAGN.append(ex1[k])
        eyAGN.append(ey1[k])
    if (y1[k] < (0.61 / (x1[k] - 0.47)) + 1.19) and (y1[k] >= (0.61/(x1[k] - 0.05)) + 1.3): #BLUE
        icomp.append(int(i1[k]))
        xcomp.append(x1[k])
        ycomp.append(y1[k])
        excomp.append(ex1[k])
        eycomp.append(ey1[k])
    if y1[k] < (0.61/(x1[k] - 0.05)) + 1.3 and x1[k] < 0.04: # GREEN
        ihii1.append(int(i1[k]))
        xhii1.append(x1[k])
        yhii1.append(y1[k])
        exhii1.append(ex1[k])
        eyhii1.append(ey1[k])


PlotScat(xAGN,yAGN,ex=exAGN,ey=eyAGN,xlim=None,ylim=None,colore="red",simbolo="o",labels=["$log([NII]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"])
PlotScat(xcomp,ycomp,ex=excomp,ey=eycomp,xlim=None,ylim=None,colore="blue",simbolo="o",labels=["$log([NII]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"],overplot=True)
PlotScat(xhii1,yhii1,ex=exhii1,ey=eyhii1,xlim=None,ylim=None,colore="green",simbolo="o",labels=["$log([NII]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"],overplot=True)
PBPT(n=1)


irad=[]
ishock=[]
ihii2=[]
xrad=[]
exrad=[]
yrad=[]
eyrad=[]
xshock=[]
exshock=[]
yshock=[]
eyshock=[]
xhii2=[]
exhii2=[]
yhii2=[]
eyhii2=[]
for k in range(len(i2)):
    if y2[k] >= (0.72 / (x2[k] - 0.32)) + 1.30 or x2[k] > 0.29:
        if y2[k] >= 1.89*x2[k] + 0.76:
            irad.append(int(i2[k]))
            xrad.append(x2[k])
            yrad.append(y2[k])
            exrad.append(ex2[k])
            eyrad.append(ey2[k])
        if y2[k] < 1.89*x2[k] + 0.76:
            ishock.append(int(i2[k]))
            xshock.append(x2[k])
            yshock.append(y2[k])
            exshock.append(ex2[k])
            eyshock.append(ey2[k])
    else:
        ihii2.append(int(i2[k]))
        xhii2.append(x2[k])
        yhii2.append(y2[k])
        exhii2.append(ex2[k])
        eyhii2.append(ey2[k])

PlotScat(xrad,yrad,ex=exrad,ey=eyrad,xlim=None,ylim=None,colore="red",simbolo="o",labels=["$log([SII]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"])
PlotScat(xshock,yshock,ex=exshock,ey=eyshock,xlim=None,ylim=None,colore="blue",simbolo="o",labels=["$log([SII]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"],overplot=True)
PlotScat(xhii2,yhii2,ex=exhii2,ey=eyhii2,xlim=None,ylim=None,colore="green",simbolo="o",labels=["$log([SII]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"],overplot=True)
PBPT(n=2)


irad3=[]
ishock3=[]
ihii3=[]
xrad3=[]
exrad3=[]
yrad3=[]
eyrad3=[]
xshock3=[]
exshock3=[]
yshock3=[]
eyshock3=[]
xhii3=[]
exhii3=[]
yhii3=[]
eyhii3=[]
for k in range(len(i3)):
    if y3[k] >= (0.73 / (x3[k] + 0.59)) + 1.33 or x3[k] > -0.6: 
        if y3[k] >= 1.18*x3[k] + 1.3:  
            irad3.append(int(i3[k]))
            xrad3.append(x3[k])
            yrad3.append(y3[k])
            exrad3.append(ex3[k])
            eyrad3.append(ey3[k])
        if y3[k] < 1.18*x3[k] + 1.3:
            ishock3.append(int(i3[k]))
            xshock3.append(x3[k])
            yshock3.append(y3[k])
            exshock3.append(ex3[k])
            eyshock3.append(ey3[k])
    else:
        ihii3.append(int(i3[k]))
        xhii3.append(x3[k])
        yhii3.append(y3[k])
        exhii3.append(ex3[k])
        eyhii3.append(ey3[k])

PlotScat(xrad3,yrad3,ex=exrad3,ey=eyrad3,xlim=None,ylim=None,colore="red",simbolo="o",labels=["$log([OI]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"])
PlotScat(xshock3,yshock3,ex=exshock3,ey=eyshock3,xlim=None,ylim=None,colore="blue",simbolo="o",labels=["$log([OI]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"],overplot=True)
PlotScat(xhii3,yhii3,ex=exhii3,ey=eyhii3,xlim=None,ylim=None,colore="green",simbolo="o",labels=["$log([OI]/H \\alpha])$","$log([OIII]/H \\beta])$"],Positives=["no","no"],overplot=True)
PBPT(n=3)



#Save subsample properties


arrays=[RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2]
if Sampl == 'no':
    SaveType(iAGN,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_AGN_gal.txt",arrays)
    SaveType(icomp,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_Comp_gal.txt",arrays)
    SaveType(ihii1,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_HII1_gal.txt",arrays)
    SaveType(irad,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_RAD_gal.txt",arrays)
    SaveType(ishock,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_SHOCK_gal.txt",arrays)
    SaveType(ihii2,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_HII2_gal.txt",arrays)
    SaveType(irad3,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_RAD3_gal.txt",arrays)
    SaveType(ishock3,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_SHOCK3_gal.txt",arrays)
    SaveType(ihii3,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_HII3_gal.txt",arrays)
else:
    SaveType(iAGN,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_AGN.txt",arrays)
    SaveType(icomp,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_Comp.txt",arrays)
    SaveType(ihii1,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_HII1.txt",arrays)
    SaveType(irad,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_RAD.txt",arrays)
    SaveType(ishock,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_SHOCK.txt",arrays)
    SaveType(ihii2,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_HII2.txt",arrays)
    SaveType(irad3,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_RAD3.txt",arrays)
    SaveType(ishock3,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_SHOCK3.txt",arrays)
    SaveType(ihii3,"/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_HII3.txt",arrays)

    

#%% PLOTS BPT subsamples

Sampl="no"
#RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl)


if Sampl == 'no':
    f1="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_AGN_gal.txt"
else:
    f1="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_AGN.txt"
RA1,DEC1,Z1,eZ1,SIG1,eSIG1,EBV1,Zsun1=np.loadtxt(f1,usecols=[0,1,2,3,4,5,6,7],unpack=True,dtype=float)
SIGMA_BAL1,eSIGMA_BAL1,SIGMA_FORB1,eSIGMA_FORB1=np.loadtxt(f1,usecols=[8,9,10,11],unpack=True,dtype=float)
VOFF_BAL1,eVOFF_BAL1,VOFF_FORB1,eVOFF_FORB1=np.loadtxt(f1,usecols=[12,13,14,15],unpack=True,dtype=float)
Mass1,eMass11,eMass12,SFR1,eSFR11,eSFR12,sSFR1,esSFR11,esSFR12=np.loadtxt(f1,usecols=[16,17,18,19,20,21,22,23,24],unpack=True,dtype=float)

if Sampl == 'no':
    f2="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_Comp_gal.txt"
else:
    f2="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_Comp.txt"
RA2,DEC2,Z2,eZ2,SIG2,eSIG2,EBV2,Zsun2=np.loadtxt(f2,usecols=[0,1,2,3,4,5,6,7],unpack=True,dtype=float)
SIGMA_BAL2,eSIGMA_BAL2,SIGMA_FORB2,eSIGMA_FORB2=np.loadtxt(f2,usecols=[8,9,10,11],unpack=True,dtype=float)
VOFF_BAL2,eVOFF_BAL2,VOFF_FORB2,eVOFF_FORB2=np.loadtxt(f2,usecols=[12,13,14,15],unpack=True,dtype=float)
Mass2,eMass21,eMass22,SFR2,eSFR21,eSFR22,sSFR2,esSFR21,esSFR22=np.loadtxt(f2,usecols=[16,17,18,19,20,21,22,23,24],unpack=True,dtype=float)

if Sampl == 'no':
    f3="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_HII1_gal.txt"
else:
    f3="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_HII1.txt"
RA3,DEC3,Z3,eZ3,SIG3,eSIG3,EBV3,Zsun3=np.loadtxt(f3,usecols=[0,1,2,3,4,5,6,7],unpack=True,dtype=float)
SIGMA_BAL3,eSIGMA_BAL3,SIGMA_FORB3,eSIGMA_FORB3=np.loadtxt(f3,usecols=[8,9,10,11],unpack=True,dtype=float)
VOFF_BAL3,eVOFF_BAL3,VOFF_FORB3,eVOFF_FORB3=np.loadtxt(f3,usecols=[12,13,14,15],unpack=True,dtype=float)
Mass3,eMass31,eMass32,SFR3,eSFR31,eSFR32,sSFR3,esSFR31,esSFR32=np.loadtxt(f3,usecols=[16,17,18,19,20,21,22,23,24],unpack=True,dtype=float)

#%% Overdensities BCG


f1="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_AGN_gal.txt"
f2="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_Comp_gal.txt"
f3="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_HII1_gal.txt"

RA01,DEC01,Z01,eZ01=np.loadtxt(f1,usecols=[0,1,2,3],unpack=True,dtype=float)
RA02,DEC02,Z02,eZ02=np.loadtxt(f2,usecols=[0,1,2,3],unpack=True,dtype=float)
RA03,DEC03,Z03,eZ03=np.loadtxt(f3,usecols=[0,1,2,3],unpack=True,dtype=float)


f1="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_AGN.txt"
f2="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_Comp.txt"
f3="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/SDSS/outputs/Prop_HII1.txt"

RA1,DEC1,Z1,eZ1=np.loadtxt(f1,usecols=[0,1,2,3],unpack=True,dtype=float)
RA2,DEC2,Z2,eZ2=np.loadtxt(f2,usecols=[0,1,2,3],unpack=True,dtype=float)
RA3,DEC3,Z3,eZ3=np.loadtxt(f3,usecols=[0,1,2,3],unpack=True,dtype=float)




RABCG=np.concatenate((RA1,RA2,RA3))
DECBCG=np.concatenate((DEC1,DEC2,DEC3))

RAnoBCG=np.concatenate((RA01,RA02,RA03))
DECnoBCG=np.concatenate((DEC01,DEC02,DEC03))

RAALL=np.concatenate((RABCG,RAnoBCG))
DECALL=np.concatenate((DECBCG,DECnoBCG))

def DistDeg(RA1,DEC1,RABCG,DECBCG,RAALL,DECALL):
    distGal_AGN=np.zeros((len(RA1),len(RAALL)))
    distBCG_AGN=np.zeros((len(RA1),len(RABCG)))
    c1 = SkyCoord(np.asarray(RAALL)*u.deg,np.asarray(DECALL)*u.deg, frame='icrs')
    cBCG= SkyCoord(np.asarray(RABCG)*u.deg,np.asarray(DECBCG)*u.deg, frame='icrs')
    for i in range(len(RA1)):
        c2 = SkyCoord(np.asarray(RA1[i])*u.deg,np.asarray(DEC1[i])*u.deg, frame='icrs')
        distGal_AGN[i,:]= c1.separation(c2).arcsecond/3600
        distBCG_AGN[i,:] = cBCG.separation(c2).arcsecond/3600
    return distGal_AGN, distBCG_AGN

distGal_AGN,distBCG_AGN=DistDeg(RA1,DEC1,RABCG,DECBCG,RAnoBCG,DECnoBCG)

distGal_Comp,distBCG_Comp=DistDeg(RA2,DEC2,RABCG,DECBCG,RAnoBCG,DECnoBCG)

distGal_HII,distBCG_HII=DistDeg(RA3,DEC3,RABCG,DECBCG,RAnoBCG,DECnoBCG)

def compute_density(ra_bcg, dec_bcg, ra_gal, dec_gal, radius_of_interest):
    # Crea oggetti SkyCoord per BCG e galassie
    bcg_coords = SkyCoord(ra_bcg*u.deg, dec_bcg*u.deg, frame='icrs')
    gal_coords = SkyCoord(ra_gal*u.deg, dec_gal*u.deg, frame='icrs')

    # Calcola separazione in gradi tra BCG e galassie
    separations = bcg_coords.separation(gal_coords).degree

    # Seleziona solo le galassie all'interno del raggio di interesse
    gal_within_radius = gal_coords[separations < radius_of_interest]

    # Calcola densità come numero di galassie entro il raggio diviso l'area del cerchio
    density = len(gal_within_radius) / (np.pi * radius_of_interest**2)

    return density

radius_of_interest = 0.01  # Imposta il raggio di interesse in gradi
bcg_density = compute_density(RABCG, DECBCG, RABCG, DECBCG, radius_of_interest)
no_bcg_density = compute_density(RAnoBCG, DECnoBCG, RAnoBCG, DECnoBCG, radius_of_interest)




fig,ax=plt.subplots()
plt.hist(np.nanmean(distBCG_AGN,axis=0)[np.nanmean(distBCG_AGN,axis=0) < 80], bins=20,color='red',alpha=0.6,density=True)
plt.hist(np.nanmean(distBCG_Comp,axis=0)[np.nanmean(distBCG_Comp,axis=0) < 80], bins=20,color='blue',alpha=0.6,density=True)
plt.hist(np.nanmean(distBCG_HII,axis=0)[np.nanmean(distBCG_HII,axis=0) < 80], bins=20,color='green',alpha=0.6,density=True)

fig,ax=plt.subplots()
plt.hist(np.nanmean(distGal_HII,axis=0)[np.nanmean(distGal_HII,axis=0) < 80], bins=20,color='green',alpha=0.6,density=True)
plt.hist(np.nanmean(distGal_AGN,axis=0)[np.nanmean(distGal_AGN,axis=0) < 80], bins=20,color='red',alpha=0.6,density=True)
plt.hist(np.nanmean(distGal_Comp,axis=0)[np.nanmean(distGal_Comp,axis=0) < 80], bins=20,color='blue',alpha=0.6,density=True)

fig,ax=plt.subplots()
plt.hist(np.nanmean(distGal_AGN,axis=0), bins=20,color='red',alpha=0.4,density=True)
plt.hist(np.nanmean(distBCG_AGN,axis=0), bins=20,color='black',alpha=0.6,density=True)

fig,ax=plt.subplots()
plt.hist(np.nanmean(distGal_Comp,axis=0), bins=20,color='blue',alpha=0.4,density=True)
plt.hist(np.nanmean(distBCG_Comp,axis=0), bins=20,color='black',alpha=0.6,density=True)

fig,ax=plt.subplots()
plt.hist(np.nanmean(distGal_HII,axis=0), bins=20,color='green',alpha=0.4,density=True)
plt.hist(np.nanmean(distBCG_Comp,axis=0), bins=20,color='black',alpha=0.6,density=True)







#%% 

"""
Come la presenza di un AGN modifica le proprietà delle BCG e delle galassie non classificate come BCG?
"""

i1,i2,i3=np.where(Zsun1 < 1)[0],np.where(Zsun2 < 1)[0],np.where(Zsun3 < 1)[0]
dati_x=[Zsun1[i1],Zsun2[i2],Zsun3[i3]]
dati_y=[SIG1[i1],SIG2[i2],SIG3[i3]]
colori=["red","blue","green"]
labels=["AGN","Composite","HII"]
assilabels=["$Z/Z_{\odot}$","$\\sigma~[km/s]$"]
Histograms(dati_x, dati_y, colori,labels,assilabels,bins=[20,20])




PlotScat(Mass1,SIGMA_FORB1,ey=eSIGMA_FORB1,xlim=None,ylim=[0.001,499],colore="red",simbolo="o")
PlotScat(Mass2,SIGMA_FORB2,ey=eSIGMA_FORB2,xlim=None,ylim=[0.001,499],colore="blue",simbolo="o",overplot=True)
PlotScat(Mass3,SIGMA_FORB3,ey=eSIGMA_FORB3,xlim=None,ylim=[0.001,499],colore="green",simbolo="o",labels=["$log(M)~[M_{\\odot}]$","$\\sigma_{Forbidden}~[km/s]$"],overplot=True)


PlotScat(Mass1,SIGMA_BAL1,ey=eSIGMA_BAL1,xlim=None,ylim=[0.001,499],colore="red",simbolo="o")
PlotScat(Mass2,SIGMA_BAL2,ey=eSIGMA_BAL2,xlim=None,ylim=[0.001,499],colore="blue",simbolo="o",overplot=True)
PlotScat(Mass3,SIGMA_BAL3,ey=eSIGMA_BAL3,xlim=None,ylim=[0.001,499],colore="green",simbolo="o",labels=["$log(M)~[M_{\\odot}]$","$\\sigma_{Balmer}~[km/s]$"],overplot=True)

PlotScat(Mass1,VOFF_BAL1,ey=eVOFF_BAL1,xlim=None,ylim=[0.001,299],colore="red",simbolo="o")
PlotScat(Mass2,VOFF_BAL2,ey=eVOFF_BAL2,xlim=None,ylim=[0.001,299],colore="blue",simbolo="o",overplot=True)
PlotScat(Mass3,VOFF_BAL3,ey=eVOFF_BAL3,xlim=None,ylim=[0.001,299],colore="green",simbolo="o",labels=["$log(M)~[M_{\\odot}]$","$v_{Balmer}~[km/s]$"],overplot=True)

PlotScat(sSFR1,VOFF_BAL1,ey=eVOFF_BAL1,xlim=None,ylim=[0.001,299],colore="red",simbolo="o")
PlotScat(sSFR2,VOFF_BAL2,ey=eVOFF_BAL2,xlim=None,ylim=[0.001,299],colore="blue",simbolo="o",overplot=True)
PlotScat(sSFR3,VOFF_BAL3,ey=eVOFF_BAL3,xlim=None,ylim=[0.001,299],colore="green",simbolo="o",labels=["$log(sSFR)~[M_{\\odot}~yr^{-1}]$","$v_{Balmer}~[km/s]$"],overplot=True)

PlotScat(SFR1,VOFF_BAL1,ey=eVOFF_BAL1,xlim=None,ylim=[0.001,299],colore="red",simbolo="o")
PlotScat(SFR2,VOFF_BAL2,ey=eVOFF_BAL2,xlim=None,ylim=[0.001,299],colore="blue",simbolo="o",overplot=True)
PlotScat(SFR3,VOFF_BAL3,ey=eVOFF_BAL3,xlim=None,ylim=[0.001,299],colore="green",simbolo="o",labels=["$log(SFR)~[M_{\\odot}~yr^{-1}]$","$v_{Balmer}~[km/s]$"],overplot=True)


PlotScat(sSFR1,EBV1,ylim=[0,2],colore="red",simbolo="o")
PlotScat(sSFR2,EBV2,ylim=[0,2],colore="blue",simbolo="o",overplot=True)
PlotScat(sSFR3,EBV3,ylim=[0,2],colore="green",simbolo="o",labels=["$log(SFR)~[M_{\\odot}~yr^{-1}]$","$E(B-V)$"],overplot=True)


i1,i2,i3=np.where((VOFF_BAL1 > 0.001) & (VOFF_BAL1 < 299))[0],np.where((VOFF_BAL2 > 0.001) & (VOFF_BAL2 < 299))[0],np.where((VOFF_BAL3 > 0.001) & (VOFF_BAL3 < 299))[0]
dati_x=[sSFR1[i1],sSFR2[i2],sSFR3[i3]]
dati_y=[VOFF_BAL1[i1],VOFF_BAL2[i2],VOFF_BAL3[i3]]
colori=["red","blue","green"]
labels=["AGN","Composite","HII"]
assilabels=["$log(sSFR)~[M_{\\odot}~yr^{-1}]$","$v_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori,labels,assilabels,bins=[100,20])

i1,i2,i3=np.where((SIGMA_BAL1 > 0.001) & (SIGMA_BAL1 < 499))[0],np.where((SIGMA_BAL2 > 0.001) & (SIGMA_BAL2 < 499))[0],np.where((SIGMA_BAL3 > 0.001) & (SIGMA_BAL3 < 499))[0]
dati_x=[sSFR1[i1],sSFR2[i2],sSFR3[i3]]
dati_y=[SIGMA_BAL1[i1],SIGMA_BAL2[i2],SIGMA_BAL3[i3]]
colori=["red","blue","green"]
labels=["AGN","Composite","HII"]
assilabels=["$log(sSFR)~[M_{\\odot}~yr^{-1}]$","$\\sigma_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori,labels,assilabels,bins=[40,20])

i1,i2,i3=np.where((SIGMA_BAL1 > 0.001) & (SIGMA_BAL1 < 499))[0],np.where((SIGMA_BAL2 > 0.001) & (SIGMA_BAL2 < 499))[0],np.where((SIGMA_BAL3 > 0.001) & (SIGMA_BAL3 < 499))[0]
dati_x=[Mass1[i1],Mass2[i2],Mass3[i3]]
dati_y=[SIGMA_BAL1[i1],SIGMA_BAL2[i2],SIGMA_BAL3[i3]]
colori=["red","blue","green"]
labels=["AGN","Composite","HII"]
assilabels=["$log(M)~[M_{\\odot}]$","$\\sigma_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori,labels,assilabels,bins=[40,20])


i1,i2,i3=np.where((VOFF_FORB1 > 0.001) & (VOFF_FORB1 < 299))[0],np.where((VOFF_FORB2 > 0.001) & (VOFF_FORB2 < 299))[0],np.where((VOFF_FORB3 > 0.001) & (VOFF_FORB3 < 299))[0]
dati_x=[sSFR1[i1],sSFR2[i2],sSFR3[i3]]
dati_y=[VOFF_FORB1[i1],VOFF_FORB2[i2],VOFF_FORB3[i3]]
colori=["red","blue","green"]
labels=["AGN","Composite","HII"]
assilabels=["$log(sSFR)~[M_{\\odot}~yr^{-1}]$","$v_{Forbidden}~[km/s]$"]
Histograms(dati_x, dati_y, colori,labels,assilabels,bins=[40,20])




#%%
#Metallicity - Sigma
PlotScat(Zsun,SIG,ey=eSIG,xlim=None,ylim=None,colore="black",simbolo="o",labels=["$Z/Z_{\\odot}$","$\\sigma~[km/s]$"])

#Metallicity - Sigma(FORBIDDEN)
PlotScat(SIGMA_BAL,SIGMA_FORB,ex=eSIGMA_BAL,ey=eSIGMA_FORB,xlim=[0.001,499],ylim=[0.001,499],colore="black",simbolo="o",labels=["$\\sigma_{BALMER}~[km/s]$","$\\sigma_{FORBIDDEN}~[km/s]$"])
plt.plot([0,500],[0,500],color='red')





