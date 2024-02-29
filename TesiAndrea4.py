#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 17:42:46 2024

@author: andreatravascio
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

# %% Selection BCG (OK)

"""
Selezione degli indici delle BCG nel campione totale DR7
"""

file = "/Users/andreamaccarinelli/Desktop/SDSS/sampleC4.fits"
fitfile = fits.open(file)
data = fitfile[1].data
raBCG, decBCG, sigClus, Ngal = data['RAdeg'], data['DEdeg'], data['sigma'], data['Ngal']


file = "/Users/andreamaccarinelli/Desktop/SDSS/gal_info_dr7_v5_2.fits"
fitfile = fits.open(file)
data = fitfile[1].data
ra, dec = data['RA'], data['DEC']


indiciBCG = []
sigmaCl = []
Ngalax = []
# SELEZIONE BCG
for k in range(len(raBCG)):
    i = np.where((abs(ra-raBCG[k]) < 2/3600) &
                 (abs(dec-decBCG[k]) < 2/3600))[0]
    if len(i) == 1:
        i = i[0]
        indiciBCG.append(i)
        sigmaCl.append(sigClus[k])
        Ngalax.append(Ngal[k])
    else:
        if len(i) > 1:
            tt = np.sqrt((ra[i]-raBCG[k])**2 + (dec[i]-decBCG[k])**2)
            i = i[np.where(tt == np.min(tt))[0][0]]
            indiciBCG.append(i)
            sigmaCl.append(sigClus[k])
            Ngalax.append(Ngal[k])

"""    
fileout="/Users/andreamaccarinelli/Desktop/BSc-Project/IndiciBCG.txt"
fmt = ["%i","%f","%i"]  # Specify the format string
data = np.column_stack((np.array(indiciBCG).astype(int),np.array(sigmaCl).astype(float),np.array(Ngalax).astype(int)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)
"""

fileout = "/Users/andreamaccarinelli/Desktop/BSc-Project/IndiciBCG.txt"
fmt = ["%i", "%f", "%i"]  # Specifica la stringa di formato
data = np.column_stack((np.array(indiciBCG).astype(int), np.array(
    sigmaCl).astype(float), np.array(Ngalax).astype(int)))

print("Array dati:", data)
print("Stringa di formato:", fmt)

if not os.path.exists(fileout):
    np.savetxt(fileout, data, fmt=fmt)

# %% Remove Index duplicates (OK)

"""
Rimozione dei duplicati e produzione di file con indici delle BCG e non BCG non duplicati
"""

SamePos = "yes"

def GalCloseBCGFIXED(RA, DEC):
    """
    Fa le stesse cose che faceva la versione precedente questa volta ritagliando
    anche in funzione dello scatterplot degli oggetti analizzati nel paper che ha prodotto
    radioTAB.txt Versione perfezionata un ultima volta da sottoporre ad Andrea...
    """
    
    i = np.where(((RA > -2) & (RA < 53) & (DEC > -10) & (DEC < 4)) | ((RA > 310) & (RA < 360) & (DEC > -10) & (DEC < 4)) | ((RA > 120) & (RA < 250) & (DEC > -5) & (DEC < 6))
                 | ((RA > 111) & (RA < 267) & (DEC > 48) & (DEC < 65)) | ((RA > 111) & (RA < 150) & (DEC > 27) & (DEC < 48)) | ((RA > 222) & (RA < 260) & (DEC > 40) & (DEC < 48))|
                 (RA > 252) & (RA < 267) & (DEC > 25.5) & (DEC < 40))[0]
    return RA[i], DEC[i], i
    

file = "/Users/andreamaccarinelli/Desktop/SDSS/duplicates.txt"

"""
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
"""

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
    #dati = np.array(dati)
    return dati


file = "/Users/andreamaccarinelli/Desktop/SDSS/duplicates.txt"
dupl = []
dupl = TxtRaws(file)


file = "/Users/andreamaccarinelli/Desktop/myOutputs3/IndiciBCG.txt"
ind, NumGalax = np.loadtxt(file, usecols=[0, 2], unpack=True, dtype=int)
SigCluster = np.loadtxt(file, usecols=[1], unpack=True, dtype=float)

IndBCG = []
Ngal = []
SBCG = []

indnoBCG = []
indiciDupl = []



for i in range(len(dupl)):
    if dupl[i][1] > 0 and (dupl[i][0] in indiciDupl) == False:
        for t in range(len(dupl[i])):
            indiciDupl.append(dupl[i][t])
        if dupl[i][0] in ind:
            IndBCG.append(dupl[i][0])
            kii = np.where(dupl[i][0] == ind)[0][0]
            Ngal.append(NumGalax[kii])
            SBCG.append(SigCluster[kii])
        else:
            indnoBCG.append(dupl[i][0])
    if dupl[i][1] < 0 and (dupl[i][0] in indiciDupl) == False:
        if dupl[i][0] in ind:
            IndBCG.append(dupl[i][0])
            kii = np.where(dupl[i][0] == ind)[0][0]
            Ngal.append(NumGalax[kii])
            SBCG.append(SigCluster[kii])
        else:
            indnoBCG.append(dupl[i][0])



fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/Indici_BCG_nodupl.txt"
fmt = ["%i", "%i", "%f"]  # Specify the format string
ddd = np.column_stack((IndBCG, Ngal, SBCG))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, ddd, fmt=fmt)


fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/Indici_noBCG_nodupl.txt"
fmt = "%i"  # Specify the format string
ddd = np.column_stack((indnoBCG, indnoBCG))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, ddd, fmt=fmt)


# %% Create file txt for only BCG and no BCG (R)

"""
Creare file INFO per le BCG e non BCG
"""
"e' stata cambiata la posizione dei file in output, in modo da non fare più"
"alcun upload dei dati stampati su github !"
"La path corrente dei risultati è la seguente :"

"/Users/andreamaccarinelli/Desktop/myOutputs3"

#SamePos = "yes"
SamePos = "no"

file = "/Users/andreamaccarinelli/Desktop/SDSS/gal_info_dr7_v5_2.fits"
fitfile = fits.open(file)
data = fitfile[1].data
ra, dec, z, zerr, zwarn, vdisp, vdisperr, ebv = data['RA'], data['DEC'], data['Z'], data[
    'Z_ERR'], data['Z_WARNING'], data['V_DISP'], data['V_DISP_ERR'], data['E_BV_SFD']

file = "/Users/andreamaccarinelli/Desktop/SDSS/gal_indx_dr7_v5_2.fits"
fitfile = fits.open(file)
data = fitfile[1].data
Zsun = data['BEST_MODEL_Z']


file = "/Users/andreamaccarinelli/Desktop/myOutputs3/Indici_BCG_nodupl.txt"
ind, Ngal = np.loadtxt(file, usecols=[0, 1], unpack=True, dtype=int)
SigCl = np.loadtxt(file, usecols=[2], unpack=True, dtype=float)
file = "/Users/andreamaccarinelli/Desktop/myOutputs3/Indici_noBCG_nodupl.txt"
indno = np.loadtxt(file, usecols=[0], unpack=True, dtype=int)

raBCG, decBCG, zBCG, zerrBCG, zwarnBCG, vdispBCG, vdisperrBCG, ebvBCG, ZsunBCG, SigmaCluster, NumeroGal = np.zeros((len(ind))), np.zeros((len(ind))), np.zeros((len(ind))), np.zeros(
    (len(ind))), np.zeros((len(ind))), np.zeros((len(ind))), np.zeros((len(ind))), np.zeros((len(ind))), np.zeros((len(ind))), np.zeros((len(ind))), np.zeros((len(ind)))
k2 = 0
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
    kk = np.where(k == ind)[0][0]
    SigmaCluster[k2] = SigCl[kk]
    NumeroGal[k2] = Ngal[kk]
    k2 += 1


fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoBCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(raBCG).astype(float), np.array(decBCG).astype(float),
                        np.array(zBCG).astype(float), np.array(
                            zerrBCG).astype(float),
                        np.array(zwarnBCG).astype(float), np.array(
                            vdispBCG).astype(float),
                        np.array(vdisperrBCG).astype(float), np.array(ebvBCG).astype(float), np.array(ZsunBCG).astype(float), np.array(SigmaCluster).astype(float), np.array(NumeroGal).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)

raBCG, decBCG, zBCG, zerrBCG, zwarnBCG, vdispBCG, vdisperrBCG, ebvBCG, ZsunBCG = np.zeros(len(indno)), np.zeros(len(indno)), np.zeros(len(
    indno)), np.zeros(len(indno)), np.zeros(len(indno)), np.zeros(len(indno)), np.zeros(len(indno)), np.zeros(len(indno)), np.zeros(len(indno))

k2 = 0
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


if SamePos == "yes":
    raBCG, decBCG, i_sampleSamePos = GalCloseBCGFIXED(raBCG, decBCG)
    zBCG = zBCG[i_sampleSamePos]
    zerrBCG = zerrBCG[i_sampleSamePos]
    zwarnBCG = zwarnBCG[i_sampleSamePos]
    vdispBCG = vdispBCG[i_sampleSamePos]
    vdisperrBCG = vdisperrBCG[i_sampleSamePos]
    ebvBCG = ebvBCG[i_sampleSamePos]
    ZsunBCG = ZsunBCG[i_sampleSamePos]


if SamePos == "yes":
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoNoBCG_SamePos.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoNoBCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(raBCG).astype(float), np.array(decBCG).astype(float),
                        np.array(zBCG).astype(float), np.array(
                            zerrBCG).astype(float),
                        np.array(zwarnBCG).astype(float), np.array(
                            vdispBCG).astype(float),
                        np.array(vdisperrBCG).astype(float), np.array(ebvBCG).astype(float), np.array(ZsunBCG).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


# %% Radio Emission? (R)

"""
Creare nuovi info file solo per capire cosa è radio loud e cosa no
"""
# S1p4 = 1.4GHz integrated flux density of the source
raR, decR, zRad, S1p4, sizeR = np.loadtxt(
    "/Users/andreamaccarinelli/Desktop/SDSS/RADIOtab.txt", usecols=[3, 4, 5, 6, 10], unpack=True, dtype=float)
radioact = np.loadtxt("/Users/andreamaccarinelli/Desktop/SDSS/RADIOtab.txt",
                      usecols=[11], unpack=True, dtype=str)


f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoBCG.txt"
ra, dec = np.loadtxt(f, usecols=[0, 1], unpack=True, dtype=float)

#print("Valori di latitudine:", dec.shape)

SamePos = "no"

S14 = []
zR = []
SR = []
Radio = []
raggi = []
c1 = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
for k in range(len(raR)):
    c2 = SkyCoord(raR[k]*u.deg, decR[k]*u.deg, frame='icrs')
    dist = c2.separation(c1).to(u.arcsec).value
    if np.min(dist) <= 5:
        S14.append(S1p4[k])
        zR.append(zRad[k])
        SR.append(sizeR[k])
        if radioact[k] == "RLQ":        #Flaggati con il numero uno se sono RadioLoud
            Radio.append(1)
        if radioact[k] == "SFG":        #Flaggati con lo zero se sono Star Forming Galaxy
            Radio.append(0)
    else:
        S14.append(0)
        zR.append(0)
        SR.append(0)
        Radio.append(0)


fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/infoRadioBCG.txt"
fmt = ["%f", "%f", "%f", "%i"]  # Specify the format string
data = np.column_stack((np.array(S14).astype(float), np.array(zR).astype(float),
                        np.array(SR).astype(float), np.array(Radio).astype(int)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


SamePos = "yes"
#SamePos = "no"

if SamePos == "yes":
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoNoBCG_SamePos.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoNoBCG.txt"
ra, dec = np.loadtxt(f, usecols=[0, 1], unpack=True, dtype=float)


S14 = []
zR = []
SR = []
Radio = []
raggi = []

nuova_RA = np.delete(ra, 143189)
nuova_Dec = np.delete(dec, 143189)
c1 = SkyCoord(nuova_RA*u.deg, nuova_Dec*u.deg, frame='icrs')
for k in range(len(raR)):
    c2 = SkyCoord(raR[k]*u.deg, decR[k]*u.deg, frame='icrs')
    dist = c2.separation(c1).to(u.arcsec).value
    if np.min(dist) <= 5:
        S14.append(S1p4[k])
        zR.append(zRad[k])
        SR.append(sizeR[k])
        if radioact[k] == "RLQ":
            Radio.append(1)
        if radioact[k] == "SFG":
            Radio.append(0)
    else:
        S14.append(0)
        zR.append(0)
        SR.append(0)
        Radio.append(0)

if SamePos == "yes":
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoRadionoBCG_SamePos.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoRadioNoBCG.txt"
fmt = ["%f", "%f", "%f", "%i"]  # Specify the format string
data = np.column_stack((np.array(S14).astype(float), np.array(zR).astype(float),
                        np.array(SR).astype(float), np.array(Radio).astype(int)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


# %% Create file txt for only BCG and no BCG (line) (R)



SamePos = "no"
file = "/Users/andreamaccarinelli/Desktop/SDSS/gal_line_dr7_v5_2.fits"
fitfile = fits.open(file)
data = fitfile[1].data

#Riempimento degli array con i dati relativi alle righe della serie di Balmer !!
sigB, esigB, sigF, esigF = data['SIGMA_BALMER'], data['SIGMA_BALMER_ERR'], data['SIGMA_FORBIDDEN'], data['SIGMA_FORBIDDEN_ERR']
vB, evB, vF, evF = data['V_OFF_BALMER'], data['V_OFF_BALMER_ERR'], data['V_OFF_FORBIDDEN'], data['V_OFF_FORBIDDEN_ERR']


#Riempimento degli array con i dati relativi alle altre linee spettrali 
#       OII_3726                 OII_3729                NEIII_3869              H_DELTA                  H_GAMMA               OIII_4363              OIII_4959                OIII_5007                HEI_5876                  OI_6300                 H_ALPHA                NII_6584                   SII_6717              SII_6731                ARIII7135
lines = ["OII_3726", "OII_3729", "NEIII_3869", "H_DELTA", "H_GAMMA", "OIII_4363", "OIII_4959", "OIII_5007", "HEI_5876",
         "OI_6300", "H_BETA", "H_ALPHA", "NII_6584", "SII_6717", "SII_6731", "ARIII7135"]
FluxLines = np.zeros((len(lines), len(sigB)))
eFluxLines = np.zeros((len(lines), len(sigB)))
for i in range(len(lines)):
    FluxLines[i, :] = data[lines[i]+'_FLUX']
    eFluxLines[i, :] = data[lines[i]+'_FLUX_ERR']


file = "/Users/andreamaccarinelli/Desktop/myOutputs3/Indici_BCG_nodupl.txt"
ind = np.loadtxt(file, usecols=[0], unpack=True, dtype=int)
file = "/Users/andreamaccarinelli/Desktop/myOutputs3/Indici_noBCG_nodupl.txt"
indno = np.loadtxt(file, usecols=[0], unpack=True, dtype=int)

sigB2 = np.zeros((len(ind)))
esigB2 = np.zeros((len(ind)))
sigF2 = np.zeros((len(ind)))
esigF2 = np.zeros((len(ind)))
vB2 = np.zeros((len(ind)))
evB2 = np.zeros((len(ind)))
vF2 = np.zeros((len(ind)))
evF2 = np.zeros((len(ind)))
FluxLines2 = np.zeros((len(lines), len(ind)))
eFluxLines2 = np.zeros((len(lines), len(ind)))


k2 = 0
for k in ind:
    sigB2[k2] = sigB[k]
    esigB2[k2] = esigB[k]
    sigF2[k2] = sigF[k]
    esigF2[k2] = esigF[k]
    vB2[k2] = vB[k]
    evB2[k2] = evB[k]
    vF2[k2] = vF[k]
    evF2[k2] = evF[k]
    FluxLines2[:, k2] = FluxLines[:, k]*1E-17
    eFluxLines2[:, k2] = eFluxLines[:, k]*1E-17
    k2 += 1


fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/SigmaLines_BCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(sigB2).astype(float), np.array(esigB2).astype(float), np.array(sigF2).astype(float), np.array(esigF2).astype(float),
                        np.array(vB2).astype(float), np.array(evB2).astype(float), np.array(vF2).astype(float), np.array(evF2).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/FluxeLines_BCG.txt"
fmt = "%.5e"  # Specify the format string

data = np.column_stack((np.array(FluxLines2[0, :]).astype(float), np.array(eFluxLines2[0, :]).astype(float),
                        np.array(FluxLines2[1, :]).astype(float), np.array(
                            eFluxLines2[1, :]).astype(float),
                        np.array(FluxLines2[2, :]).astype(float), np.array(
                            eFluxLines2[2, :]).astype(float),
                        np.array(FluxLines2[3, :]).astype(float), np.array(
                            eFluxLines2[3, :]).astype(float),
                        np.array(FluxLines2[4, :]).astype(float), np.array(
                            eFluxLines2[4, :]).astype(float),
                        np.array(FluxLines2[5, :]).astype(float), np.array(
                            eFluxLines2[5, :]).astype(float),
                        np.array(FluxLines2[6, :]).astype(float), np.array(
                            eFluxLines2[6, :]).astype(float),
                        np.array(FluxLines2[7, :]).astype(float), np.array(
                            eFluxLines2[7, :]).astype(float),
                        np.array(FluxLines2[8, :]).astype(float), np.array(
                            eFluxLines2[8, :]).astype(float),
                        np.array(FluxLines2[9, :]).astype(float), np.array(
                            eFluxLines2[9, :]).astype(float),
                        np.array(FluxLines2[10, :]).astype(float), np.array(
                            eFluxLines2[10, :]).astype(float),
                        np.array(FluxLines2[11, :]).astype(float), np.array(
                            eFluxLines2[11, :]).astype(float),
                        np.array(FluxLines2[12, :]).astype(float), np.array(
                            eFluxLines2[12, :]).astype(float),
                        np.array(FluxLines2[13, :]).astype(float), np.array(
                            eFluxLines2[13, :]).astype(float),
                        np.array(FluxLines2[14, :]).astype(float), np.array(
                            eFluxLines2[14, :]).astype(float),
                        np.array(FluxLines2[15, :]).astype(float), np.array(eFluxLines2[15, :]).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)

##### NO BCG ######

sigB2 = np.zeros((len(sigB)-len(ind)))
esigB2 = np.zeros((len(sigB)-len(ind)))
sigF2 = np.zeros((len(sigB)-len(ind)))
esigF2 = np.zeros((len(sigB)-len(ind)))
vB2 = np.zeros((len(sigB)-len(ind)))
evB2 = np.zeros((len(sigB)-len(ind)))
vF2 = np.zeros((len(sigB)-len(ind)))
evF2 = np.zeros((len(sigB)-len(ind)))
FluxLines2 = np.zeros((len(lines), len(sigB)-len(ind)))
eFluxLines2 = np.zeros((len(lines), len(sigB)-len(ind)))


k2 = 0
for k in indno:
    sigB2[k2] = sigB[k]
    esigB2[k2] = esigB[k]
    sigF2[k2] = sigF[k]
    esigF2[k2] = esigF[k]
    vB2[k2] = vB[k]
    evB2[k2] = evB[k]
    vF2[k2] = vF[k]
    evF2[k2] = evF[k]
    FluxLines2[:, k2] = FluxLines[:, k]*1E-17
    eFluxLines2[:, k2] = eFluxLines[:, k]*1E-17
    k2 += 1


if SamePos == "yes":
    sigB2 = sigB2[i_sampleSamePos]
    esigB2 = esigB2[i_sampleSamePos]
    sigF2 = sigF2[i_sampleSamePos]
    esigF2 = esigF2[i_sampleSamePos]
    vB2 = vB2[i_sampleSamePos]
    evB2 = evB2[i_sampleSamePos]
    vF2 = vF2[i_sampleSamePos]
    evF2 = evF2[i_sampleSamePos]
    FluxLines2 = FluxLines2[:, i_sampleSamePos]
    eFluxLines2 = eFluxLines2[:, i_sampleSamePos]


if SamePos == "yes":
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/SigmaLines_noBCG_SamePos.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/SigmaLines_noBCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(sigB2).astype(float), np.array(esigB2).astype(float), np.array(sigF2).astype(float), np.array(esigF2).astype(float),
                        np.array(vB2).astype(float), np.array(evB2).astype(float), np.array(vF2).astype(float), np.array(evF2).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


if SamePos == "yes":
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/FluxeLines_noBCG_SamePos.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/FluxeLines_noBCG.txt"
fmt = "%.5e"  # Specify the format string

data = np.column_stack((np.array(FluxLines2[0, :]).astype(float), np.array(eFluxLines2[0, :]).astype(float),
                        np.array(FluxLines2[1, :]).astype(float), np.array(
                            eFluxLines2[1, :]).astype(float),
                        np.array(FluxLines2[2, :]).astype(float), np.array(
                            eFluxLines2[2, :]).astype(float),
                        np.array(FluxLines2[3, :]).astype(float), np.array(
                            eFluxLines2[3, :]).astype(float),
                        np.array(FluxLines2[4, :]).astype(float), np.array(
                            eFluxLines2[4, :]).astype(float),
                        np.array(FluxLines2[5, :]).astype(float), np.array(
                            eFluxLines2[5, :]).astype(float),
                        np.array(FluxLines2[6, :]).astype(float), np.array(
                            eFluxLines2[6, :]).astype(float),
                        np.array(FluxLines2[7, :]).astype(float), np.array(
                            eFluxLines2[7, :]).astype(float),
                        np.array(FluxLines2[8, :]).astype(float), np.array(
                            eFluxLines2[8, :]).astype(float),
                        np.array(FluxLines2[9, :]).astype(float), np.array(
                            eFluxLines2[9, :]).astype(float),
                        np.array(FluxLines2[10, :]).astype(float), np.array(
                            eFluxLines2[10, :]).astype(float),
                        np.array(FluxLines2[11, :]).astype(float), np.array(
                            eFluxLines2[11, :]).astype(float),
                        np.array(FluxLines2[12, :]).astype(float), np.array(
                            eFluxLines2[12, :]).astype(float),
                        np.array(FluxLines2[13, :]).astype(float), np.array(
                            eFluxLines2[13, :]).astype(float),
                        np.array(FluxLines2[14, :]).astype(float), np.array(
                            eFluxLines2[14, :]).astype(float),
                        np.array(FluxLines2[15, :]).astype(float), np.array(eFluxLines2[15, :]).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


# %% Create file txt for only BCG and no BCG (sSFR & Mass) (R)

#SamePos = "yes"
SamePos = "no"
file = "/Users/andreamaccarinelli/Desktop/SDSS/gal_totsfr_dr7_v5_2.fits"
fitfile = fits.open(file)
data = fitfile[1].data
avg, med, p16, p84 = data['AVG'], data['MEDIAN'], data['P16'], data['P84']

file = "/Users/andreamaccarinelli/Desktop/SDSS/gal_totspecsfr_dr7_v5_2.fits"
fitfile = fits.open(file)
data = fitfile[1].data
savg, smed, sp16, sp84 = data['AVG'], data['MEDIAN'], data['P16'], data['P84']


file = "/Users/andreamaccarinelli/Desktop/SDSS/totlgm_dr7_v5_2.fits"
fitfile = fits.open(file)
data = fitfile[1].data
mavg, mmed, mp16, mp84 = data['AVG'], data['MEDIAN'], data['P16'], data['P84']


Ntot = len(avg)

file = "/Users/andreamaccarinelli/Desktop/myOutputs3/Indici_BCG_nodupl.txt"
ind = np.loadtxt(file, usecols=[0], unpack=True, dtype=int)
file = "/Users/andreamaccarinelli/Desktop/myOutputs3/Indici_noBCG_nodupl.txt"
indno = np.loadtxt(file, usecols=[0], unpack=True, dtype=int)

avg2, med2, p162, p842 = np.zeros((len(ind))), np.zeros(
    (len(ind))), np.zeros((len(ind))), np.zeros((len(ind)))
savg2, smed2, sp162, sp842 = np.zeros((len(ind))), np.zeros(
    (len(ind))), np.zeros((len(ind))), np.zeros((len(ind)))
mavg2, mmed2, mp162, mp842 = np.zeros((len(ind))), np.zeros(
    (len(ind))), np.zeros((len(ind))), np.zeros((len(ind)))

k2 = 0
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


fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/MassSFR_BCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(mavg2).astype(float), np.array(mmed2).astype(float),
                        np.array(mp162).astype(float), np.array(
                            mp842).astype(float),
                        np.array(avg2).astype(float), np.array(
                            med2).astype(float),
                        np.array(p162).astype(float), np.array(
                            p842).astype(float),
                        np.array(savg2).astype(float), np.array(
                            smed2).astype(float),
                        np.array(sp162).astype(float), np.array(sp842).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


avg2, med2, p162, p842 = np.zeros((Ntot-len(ind))), np.zeros(
    (Ntot-len(ind))), np.zeros((Ntot-len(ind))), np.zeros((Ntot-len(ind)))
savg2, smed2, sp162, sp842 = np.zeros((Ntot-len(ind))), np.zeros(
    (Ntot-len(ind))), np.zeros((Ntot-len(ind))), np.zeros((Ntot-len(ind)))
mavg2, mmed2, mp162, mp842 = np.zeros((Ntot-len(ind))), np.zeros(
    (Ntot-len(ind))), np.zeros((Ntot-len(ind))), np.zeros((Ntot-len(ind)))


k2 = 0
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


if SamePos == "yes":
    avg2 = avg2[i_sampleSamePos]
    med2 = med2[i_sampleSamePos]
    p162 = p162[i_sampleSamePos]
    p842 = p842[i_sampleSamePos]
    savg2 = savg2[i_sampleSamePos]
    smed2 = smed2[i_sampleSamePos]
    sp162 = sp162[i_sampleSamePos]
    sp842 = sp842[i_sampleSamePos]
    mavg2 = mavg2[i_sampleSamePos]
    mmed2 = mmed2[i_sampleSamePos]
    mp162 = mp162[i_sampleSamePos]
    mp842 = mp842[i_sampleSamePos]


if SamePos == "yes":
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/MassSFR_noBCG_SamePos.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/MassSFR_noBCG.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(mavg2).astype(float), np.array(mmed2).astype(float),
                        np.array(mp162).astype(float), np.array(
                            mp842).astype(float),
                        np.array(avg2).astype(float), np.array(
                            med2).astype(float),
                        np.array(p162).astype(float), np.array(
                            p842).astype(float),
                        np.array(savg2).astype(float), np.array(
                            smed2).astype(float),
                        np.array(sp162).astype(float), np.array(sp842).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


# %% CALL DATA (R)

Sampl = "no"  # "no" "y"
#SamePos = "yes"
SamePos = "no"

def calldata(Sampl='y', SamePos='no', indices = None):
    # INFO BCG
    if Sampl == 'no':
        if SamePos == "yes":
            f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoNoBCG_SamePos.txt"
        else:
            f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoNoBCG.txt"
        RA, DEC, Z, eZ, SIG, eSIG, EBV, Zsun = np.loadtxt(
            f, usecols=[0, 1, 2, 3, 5, 6, 7, 8], unpack=True, dtype=float)
        SIGCLUSTER = np.zeros((len(RA)))
        NUMGAL = np.zeros((len(RA)))
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoBCG.txt"
        RA, DEC, Z, eZ, SIG, eSIG, EBV, Zsun, SIGCLUSTER, NUMGAL = np.loadtxt(
            f, usecols=[0, 1, 2, 3, 5, 6, 7, 8, 9, 10], unpack=True, dtype=float)

    # Sigma Balmer and Forbidden Lines
    if Sampl == 'no':
        if SamePos == "yes":
            f = "/Users/andreamaccarinelli/Desktop/myOutputs3/SigmaLines_noBCG_SamePos.txt"
        else:
            f = "/Users/andreamaccarinelli/Desktop/myOutputs3/SigmaLines_noBCG.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/SigmaLines_BCG.txt"
    SIGMA_BAL, eSIGMA_BAL, SIGMA_FORB, eSIGMA_FORB = np.loadtxt(
        f, usecols=[0, 1, 2, 3], unpack=True, dtype=float)
    VOFF_BAL, eVOFF_BAL, VOFF_FORB, eVOFF_FORB = np.loadtxt(
        f, usecols=[4, 5, 6, 7], unpack=True, dtype=float)

    # FLUXES LINES
    if Sampl == 'no':
        if SamePos == "yes":
            f = "/Users/andreamaccarinelli/Desktop/myOutputs3/FluxeLines_noBCG_SamePos.txt"
        else:
            f = "/Users/andreamaccarinelli/Desktop/myOutputs3/FluxeLines_noBCG.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/FluxeLines_BCG.txt"
    OII_3726, eOII_3726, OII_3729, eOII_3729, NEIII_3869, eNEIII_3869 = np.loadtxt(
        f, usecols=[0, 1, 2, 3, 4, 5], unpack=True, dtype=float)
    H_DELTA, eH_DELTA, H_GAMMA, eH_GAMMA, OIII_4363, eOIII_4363, OIII_4959, eOIII_4959 = np.loadtxt(
        f, usecols=[6, 7, 8, 9, 10, 11, 12, 13], unpack=True, dtype=float)
    OIII_5007, eOIII_5007, HEI_5876, eHEI_5876, OI_6300, eOI_6300 = np.loadtxt(
        f, usecols=[14, 15, 16, 17, 18, 19], unpack=True, dtype=float)
    H_BETA, eH_BETA, H_ALPHA, eH_ALPHA, NII_6584, eNII_6584 = np.loadtxt(
        f, usecols=[20, 21, 22, 23, 24, 25], unpack=True, dtype=float)
    SII_6717, eSII_6717, SII_6731, eSII_6731, ARIII7135, eARIII7135 = np.loadtxt(
        f, usecols=[26, 27, 28, 29, 30, 31], unpack=True, dtype=float)

    # Derived Prop. Mass SFR sSFR
    if Sampl == 'no':
        if SamePos == "yes":
            f = "/Users/andreamaccarinelli/Desktop/myOutputs3/MassSFR_noBCG_SamePos.txt"
        else:
            f = "/Users/andreamaccarinelli/Desktop/myOutputs3/MassSFR_noBCG.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/MassSFR_BCG.txt"
    Mass, eMass1, eMass2, SFR, eSFR1, eSFR2, sSFR, esSFR1, esSFR2 = np.loadtxt(
        f, usecols=[0, 2, 3, 4, 6, 7, 8, 10, 11], unpack=True, dtype=float)
    
    # Utilizza gli indici specificati o tutti gli indici se indices è None
    if indices is not None:
        RA = RA[indices]
        DEC = DEC[indices]
        Z = Z[indices]
        eZ = eZ[indices]
        SIG = SIG[indices]
        eSIG = eSIG[indices]
        EBV = EBV[indices]
        Zsun = Zsun[indices]
        SIGCLUSTER = SIGCLUSTER[indices]
        NUMGAL = NUMGAL[indices]
        SIGMA_BAL = SIGMA_BAL[indices]
        eSIGMA_BAL = eSIGMA_BAL[indices]
        SIGMA_FORB = SIGMA_FORB[indices]
        eSIGMA_FORB = eSIGMA_FORB[indices]
        VOFF_BAL = VOFF_BAL[indices]
        eVOFF_BAL = eVOFF_BAL[indices]
        VOFF_FORB = VOFF_FORB[indices]
        eVOFF_FORB = eVOFF_FORB[indices]
        OII_3726 = OII_3726[indices]
        eOII_3726 = eOII_3726[indices]
        OII_3729 = OII_3729[indices]
        eOII_3729 = eOII_3729[indices]
        NEIII_3869 = NEIII_3869[indices]
        eNEIII_3869 = eNEIII_3869[indices]
        H_DELTA = H_DELTA[indices]
        eH_DELTA = eH_DELTA[indices]
        H_GAMMA = H_GAMMA[indices]
        eH_GAMMA = eH_GAMMA[indices]
        OIII_4363 = OIII_4363[indices]
        eOIII_4363 = eOIII_4363[indices]
        OIII_4959 = OIII_4959[indices]
        eOIII_4959 = eOIII_4959[indices]
        OIII_5007 = OIII_5007[indices]
        eOIII_5007 = eOIII_5007[indices]
        HEI_5876 = HEI_5876[indices]
        eHEI_5876 = eHEI_5876[indices]
        OI_6300 = OI_6300[indices]
        eOI_6300 = eOI_6300[indices]
        H_BETA = H_BETA[indices]
        eH_BETA = eH_BETA[indices]
        H_ALPHA = H_ALPHA[indices]
        eH_ALPHA = eH_ALPHA[indices]
        NII_6584 = NII_6584[indices]
        eNII_6584 = eNII_6584[indices]
        SII_6717 = SII_6717[indices]
        eSII_6717 = eSII_6717[indices]
        SII_6731 = SII_6731[indices]
        eSII_6731 = eSII_6731[indices]
        ARIII7135 = ARIII7135[indices]
        eARIII7135 = eARIII7135[indices]
        Mass = Mass[indices]
        eMass1 = eMass1[indices]
        eMass2 = eMass2[indices]
        SFR = SFR[indices]
        eSFR1 = eSFR1[indices]
        eSFR2 = eSFR2[indices]
        sSFR = sSFR[indices]
        esSFR1 = esSFR1[indices]
        esSFR2 = esSFR2[indices]
       
    
    
    
    return RA, DEC, Z, eZ, SIG, eSIG, EBV, Zsun, SIGCLUSTER, NUMGAL, SIGMA_BAL, eSIGMA_BAL, SIGMA_FORB, eSIGMA_FORB, VOFF_BAL, eVOFF_BAL, VOFF_FORB, eVOFF_FORB, OII_3726, eOII_3726, OII_3729, eOII_3729, NEIII_3869, eNEIII_3869, H_DELTA, eH_DELTA, H_GAMMA, eH_GAMMA, OIII_4363, eOIII_4363, OIII_4959, eOIII_4959, OIII_5007, eOIII_5007, HEI_5876, eHEI_5876, OI_6300, eOI_6300, H_BETA, eH_BETA, H_ALPHA, eH_ALPHA, NII_6584, eNII_6584, SII_6717, eSII_6717, SII_6731, eSII_6731, ARIII7135, eARIII7135, Mass, eMass1, eMass2, SFR, eSFR1, eSFR2, sSFR, esSFR1, esSFR2


RA, DEC, Z, eZ, SIG, eSIG, EBV, Zsun, SIGCLUSTER, NUMGAL, SIGMA_BAL, eSIGMA_BAL, SIGMA_FORB, eSIGMA_FORB, VOFF_BAL, eVOFF_BAL, VOFF_FORB, eVOFF_FORB, OII_3726, eOII_3726, OII_3729, eOII_3729, NEIII_3869, eNEIII_3869, H_DELTA, eH_DELTA, H_GAMMA, eH_GAMMA, OIII_4363, eOIII_4363, OIII_4959, eOIII_4959, OIII_5007, eOIII_5007, HEI_5876, eHEI_5876, OI_6300, eOI_6300, H_BETA, eH_BETA, H_ALPHA, eH_ALPHA, NII_6584, eNII_6584, SII_6717, eSII_6717, SII_6731, eSII_6731, ARIII7135, eARIII7135, Mass, eMass1, eMass2, SFR, eSFR1, eSFR2, sSFR, esSFR1, esSFR2 = calldata(
    Sampl=Sampl, SamePos=SamePos)


# %% Function PLOTS (R)


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
        ax_top.hist(dati_x[i], bins=bins[0], color=colori[i],
                    alpha=0.7, edgecolor='black', density=Normalized == 'y')
    # ax_top.set_title('Istogramma su X')

    # Istogramma a destra
    ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax_main)
    for i in range(len(dati_y)):
        ax_right.hist(dati_y[i], bins=bins[1], orientation='horizontal',
                      color=colori[i], alpha=0.7, edgecolor='black', density=Normalized == 'y')
    # ax_right.set_title('Istogramma su Y')

    # Rimuovi etichette degli assi del subplot principale
    # ax_top.set_xticks([])
    # ax_right.set_yticks([])

    # Visualizza il plot principale
    ax_main.legend()
    ax_main.set_xlabel(assi_labels[0])
    ax_main.set_ylabel(assi_labels[1])
    plt.show()


def PlotScat(x, y, ex=None, ey=None, xlim=None, ylim=None, colore="black", simbolo="o", labels=["X", "Y"], Positives=["yes", "yes"], overplot=False):
    if (xlim is None) == True:
        pass
    else:
        if np.isnan(xlim[0]) == False and np.isnan(xlim[1]) == False:
            indexes = np.where((x >= xlim[0]) & (x <= xlim[1]))[0]
        else:
            if np.isnan(xlim[0]) == False:
                indexes = np.where(x <= xlim[1])[0]
            if np.isnan(xlim[1]) == False:
                indexes = np.where(x <= xlim[1])[0]
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
            indexes = np.where((y >= ylim[0]) & (y <= ylim[1]))[0]
        else:
            if np.isnan(ylim[0]) == False:
                indexes = np.where(y <= ylim[1])[0]
            if np.isnan(ylim[1]) == False:
                indexes = np.where(y <= ylim[1])[0]
        x = x[indexes]
        y = y[indexes]
        if (ex is None) == False:
            ex = ex[indexes]
        if (ey is None) == False:
            ey = ey[indexes]
    # Remove Large Errors
    if Positives[0] == "yes":
        if (ex is None) == False:
            indexes = np.where(x-ex > 0)[0]
            x = x[indexes]
            y = y[indexes]
            ex = ex[indexes]
            if (ey is None) == False:
                ey = ey[indexes]
    if Positives[1] == "yes":
        if (ey is None) == False:
            indexes = np.where(y-ey > 0)[0]
            x = x[indexes]
            y = y[indexes]
            if (ex is None) == False:
                ex = ex[indexes]
            ey = ey[indexes]
    if overplot == False:
        fig, ax = plt.subplots()
    plt.scatter(x, y, color=colore, linestyle='None', marker=simbolo)
    if (ex is None) == False and (ey is None) == False:
        plt.errorbar(x, y, xerr=ex, yerr=ey, color=colore, linestyle='None')
    if (ex is None) == False and (ey is None) == True:
        plt.errorbar(x, y, xerr=ex, color=colore, linestyle='None')
    if (ex is None) == True and (ey is None) == False:
        plt.errorbar(x, y, yerr=ey, color=colore, linestyle='None')

    plt.xlabel(labels[0], fontsize=16)
    plt.ylabel(labels[1], fontsize=16)
    plt.tick_params(axis='both', labelsize=16)
    return len(x)

"""
def add_error_cross(x, y, ex, ey):
    x_error_percentile = np.percentile(ex, 75)
    y_error_percentile = np.percentile(ey, 75)

    plt.errorbar(x, y, xerr=ex, yerr=ey, linestyle='None', marker='o', color='blue')
    plt.errorbar(x_error_percentile, y_error_percentile, marker='x', color='red', markersize=10, label='75th Percentile Error')
    plt.legend()
    plt.show()
"""

def ErrLogRatio(num, den, err_num=None, err_den=None, Niter=1000):
    errors = np.zeros((len(num)))
    for i in range(len(num)):
        if (err_num[i] is None) == False:
            rndnum = np.random.normal(num[i], err_num[i], Niter)
            rndnum = np.clip(rndnum, num[i] - err_num[i], num[i] + err_num[i])

        if (err_den[i] is None) == False:
            rndden = np.random.normal(den[i], err_den[i], Niter)
            rndden = np.clip(rndden, den[i] - err_den[i], den[i] + err_den[i])
        logval = np.log10(rndnum/rndden)
        errors[i] = np.std(logval)
    return errors



def scatter_plot(x, y, xlim=None, ylim=None, color="black", marker="o", labels=["X", "Y"], overplot=False):
    if xlim is not None:
        x_mask = np.logical_and(x >= xlim[0], x <= xlim[1])
        x = x[x_mask]
        y = y[x_mask]

    if ylim is not None:
        y_mask = np.logical_and(y >= ylim[0], y <= ylim[1])
        x = x[y_mask]
        y = y[y_mask]

    if not overplot:
        fig, ax = plt.subplots()

    plt.scatter(x, y,s= 5, color=color, marker=marker, label="Data",)

    plt.xlabel(labels[0], fontsize=16)
    plt.ylabel(labels[1], fontsize=16)
    plt.tick_params(axis='both', labelsize=16)

    return len(x)
"""
def add_error_cross(x, y, xerr, yerr, xpos = 1, ypos = -1):
    
    Parameters
    ----------
    x : Array delle ascisse
    y : Array delle ordinate
    xerr : Array degli errori sulle x
    yerr : Array degli errori sulle y
    xpos : Ascissa del centro del mirino
        DESCRIPTION. The default is 1.
    ypos : Ordinata del centro del mirino
        DESCRIPTION. The default is -1.

    Returns
    -------
    None.

    
    x_cross = np.median(xerr)  # 75 percentile per gli errori su x
    y_cross = np.median(yerr)  # 75 percentile per gli errori su y
    plt.errorbar(xpos, ypos, xerr=np.percentile(x_cross, 75), yerr=np.percentile(y_cross, 75),
                 linestyle='None', marker='x', color='red', markersize=10, label='Incertezza al 75%')

"""

def add_error_cross(x, y, xerr, yerr, xpos=1, ypos=-1):
    """
    Parameters
    ----------
    x : Array delle ascisse
    y : Array delle ordinate
    xerr : Array degli errori sulle x
    yerr : Array degli errori sulle y
    xpos : Ascissa del centro del mirino
        DESCRIPTION. The default is 1.
    ypos : Ordinata del centro del mirino
        DESCRIPTION. The default is -1.

    Returns
    -------
    None.

    """
    x_cross = np.median(xerr)  # 75 percentile per gli errori su x
    y_cross = np.median(yerr)  # 75 percentile per gli errori su y

    #plt.errorbar(xpos, ypos, xerr=np.percentile(x_cross, 75), yerr=np.percentile(y_cross, 75),
                 #linestyle='None', marker='x', color='red', markersize=10, label='Incertezza al 75%', capsize=0)
                 
    plt.errorbar(1, -1, xerr=np.percentile(x_cross, 75), yerr=np.percentile(y_cross, 75))

import matplotlib.patches as patches

def add_error_box(x, y, xerr, yerr, xpos=1, ypos=-1):
    """
    Parameters
    ----------
    x : Array delle ascisse
    y : Array delle ordinate
    xerr : Array degli errori sulle x
    yerr : Array degli errori sulle y
    xpos : Ascissa del centro del mirino
        DESCRIPTION. The default is 1.
    ypos : Ordinata del centro del mirino
        DESCRIPTION. The default is -1.

    Returns
    -------
    None.

    """
    x_cross = np.median(xerr)  # 75 percentile per gli errori su x
    y_cross = np.median(yerr)  # 75 percentile per gli errori su y

    # Creare una Rectangle con le dimensioni degli errori sulle x e y
    error_box = patches.Rectangle((xpos - np.percentile(x_cross, 75), ypos - np.percentile(y_cross, 75)),
                                  2 * np.percentile(x_cross, 75), 2 * np.percentile(y_cross, 75),
                                  edgecolor='white', linewidth=2, facecolor='none', label='Incertezza al 75%')

    # Aggiungere la Rectangle al grafico
    plt.gca().add_patch(error_box)

    # Aggiungere linee rosse passanti per i punti medi dei lati del rettangolo
    x_left, y_mid = xpos - np.percentile(x_cross, 75), ypos
    x_right, y_mid = xpos + np.percentile(x_cross, 75), ypos
    plt.plot([x_left, x_right], [y_mid, y_mid], color='red', linestyle='--', linewidth=1.5)

    x_mid, y_top = xpos, ypos + np.percentile(y_cross, 75)
    x_mid, y_bottom = xpos, ypos - np.percentile(y_cross, 75)
    plt.plot([x_mid, x_mid], [y_top, y_bottom], color='red', linestyle='--', linewidth=1.5)



# %% BPT DIAGRAMS (R)

Sampl = "no"
#RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl, SamePos)

RA, DEC, Z, eZ, SIG, eSIG, EBV, Zsun, SIGCLUSTER, NUMGAL, SIGMA_BAL, eSIGMA_BAL, SIGMA_FORB, eSIGMA_FORB, VOFF_BAL, eVOFF_BAL, VOFF_FORB, eVOFF_FORB, OII_3726, eOII_3726, OII_3729, eOII_3729, NEIII_3869, eNEIII_3869, H_DELTA, eH_DELTA, H_GAMMA, eH_GAMMA, OIII_4363, eOIII_4363, OIII_4959, eOIII_4959, OIII_5007, eOIII_5007, HEI_5876, eHEI_5876, OI_6300, eOI_6300, H_BETA, eH_BETA, H_ALPHA, eH_ALPHA, NII_6584, eNII_6584, SII_6717, eSII_6717, SII_6731, eSII_6731, ARIII7135, eARIII7135, Mass, eMass1, eMass2, SFR, eSFR1, eSFR2, sSFR, esSFR1, esSFR2 = calldata(Sampl, SamePos)


#SamePos = "yes"
SamePos = "yes"

"""
OIIIHb = 0.61 / (NIIHa - 0.05) + 1.3     #(Kauffmann+03 line)
OIIIHb = 0.61 / (NIIHa - 0.47) + 1.19    #(Kewley+01 line)

OIIIHb = 0.72 / (SIIHa - 0.32) + 1.30    #(main AGN line)
OIIIHb = 1.89*SIIHa + 0.76               #(LINER/Sy2 line)

OIIIHb = 0.73 / (OIHa + 0.59) + 1.33     #(main AGN line)
OIIIHb = 1.18*OIHa + 1.30                 #(LINER/Sy2 line)
"""

# REMOVE FALSE VALUES
indicitot = np.arange(len(OIII_5007))
indexes1 = np.where((OIII_5007 > 0) & (H_ALPHA > 0) & (H_BETA > 0) & (NII_6584 > 0) & (
    OIII_5007 > eOIII_5007) & (H_ALPHA > eH_ALPHA) & (H_BETA > eH_BETA) & (NII_6584 > eNII_6584))[0]
indexes2 = np.where((OIII_5007 > 0) & (H_ALPHA > 0) & (H_BETA > 0) & (SII_6717 > 0) & (SII_6731 > 0) & (
    OIII_5007 > eOIII_5007) & (H_ALPHA > eH_ALPHA) & (H_BETA > eH_BETA) & (SII_6717 > eSII_6717) & (SII_6731 > eSII_6731))[0]
indexes3 = np.where((OIII_5007 > 0) & (H_ALPHA > 0) & (H_BETA > 0) & (OI_6300 > 0) & (
    OIII_5007 > eOIII_5007) & (H_ALPHA > eH_ALPHA) & (H_BETA > eH_BETA) & (OI_6300 > eOI_6300))[0]

logOIIIHb1 = np.log10(OIII_5007[indexes1]/H_BETA[indexes1])
elogOIIIHb1 = ErrLogRatio(OIII_5007[indexes1], H_BETA[indexes1],
                          err_num=eOIII_5007[indexes1], err_den=eH_BETA[indexes1])
logNIIHa1 = np.log10(NII_6584[indexes1]/H_ALPHA[indexes1])
elogNIIHa1 = ErrLogRatio(NII_6584[indexes1], H_ALPHA[indexes1],
                         err_num=eNII_6584[indexes1], err_den=eH_ALPHA[indexes1])


i1 = indicitot[indexes1]   #Indici dell'array con tutte le galassie su cui è pox fare BPTtype1 



logOIIIHb2 = np.log10(OIII_5007[indexes2]/H_BETA[indexes2])
elogOIIIHb2 = ErrLogRatio(OIII_5007[indexes2], H_BETA[indexes2],
                          err_num=eOIII_5007[indexes2], err_den=eH_BETA[indexes2])
SII = SII_6731[indexes2]
eSII = eSII_6731[indexes2]
logSIIHa2 = np.log10((SII)/H_ALPHA[indexes2])
elogSIIHa2 = ErrLogRatio(
    SII, H_ALPHA[indexes2], err_num=eSII, err_den=eH_ALPHA[indexes2])


i2 = indicitot[indexes2]    #Indici dell'array con tutte le galassie su cui è pox fare BPTtype2 

logOIIIHb3 = np.log10(OIII_5007[indexes3]/H_BETA[indexes3])
elogOIIIHb3 = ErrLogRatio(OIII_5007[indexes3], H_BETA[indexes3],
                          err_num=eOIII_5007[indexes3], err_den=eH_BETA[indexes3])
logOIHa3 = np.log10(OI_6300[indexes3]/H_ALPHA[indexes3])
elogOIHa3 = ErrLogRatio(OI_6300[indexes3], H_ALPHA[indexes3],
                        err_num=eOI_6300[indexes3], err_den=eH_ALPHA[indexes3])


i3 = indicitot[indexes3]   #Indici dell'array con tutte le galassie su cui è pox fare BPTtype3 


def BPTd():
    xx1 = np.arange(-3, 3, 0.01)
    xx2 = np.arange(-3, 3, 0.01)
    xx3 = np.arange(-3, 3, 0.01)
    yy1_1 = (0.61/(xx1[xx1 < 0.05] - 0.05)) + 1.3
    yy1_2 = (0.61 / (xx1[xx1 < 0.47] - 0.47)) + 1.19

    yy2_1 = (0.72 / (xx2[xx2 < 0.32] - 0.32)) + 1.30  # (main AGN line)
    yy2_2 = 1.89*xx2[xx2 > -0.33] + 0.76

    yy3_1 = (0.73 / (xx3[xx3 < -0.59] + 0.59)) + 1.33  # (main AGN line)
    yy3_2 = 1.18*xx3[xx3 > -1.1257] + 1.30
    return xx1, xx2, xx3, yy1_1, yy1_2, yy2_1, yy2_2, yy3_1, yy3_2


def PBPT(n=1):
    xx1, xx2, xx3, yy1_1, yy1_2, yy2_1, yy2_2, yy3_1, yy3_2 = BPTd()
    if n == 1:
        plt.plot(xx1[xx1 < 0.05], yy1_1, color='black')
        plt.plot(xx1[xx1 < 0.47], yy1_2, color='black')
        plt.xlim((-2, 3))
        plt.ylim((-2, 3))
    if n == 2:
        plt.plot(xx2[xx2 < 0.32], yy2_1, color='black')
        plt.plot(xx2[xx2 > -0.33], yy2_2, color='black')
        plt.xlim((-2, 3))
        plt.ylim((-2, 3))
    if n == 3:
        plt.plot(xx3[xx3 < -0.59], yy3_1, color='black')
        plt.plot(xx3[xx3 > -1.1257], yy3_2, color='black')
        plt.xlim((-2, 3))
        plt.ylim((-2, 3))
    plt.subplots_adjust(top=0.8, bottom=0.2, left=0.2,
                        right=0.8, hspace=0.2, wspace=0.2)
    return 0



#DA fare girare uno per volta nella console una volta che la cella è stata fatta girare !!!
"""
#NII
add_error_box(logNIIHa1, logOIIIHb1, elogNIIHa1, elogOIIIHb1, 1, -1)
scatter_plot(logNIIHa1, logOIIIHb1, overplot=True, color = "blue", labels=["$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"])
PBPT(n=1)
plt.text(1, -0.7, 'Error at 75%', ha = 'center')
plt.show()

#SII
add_error_box(logSIIHa2, logOIIIHb2, elogSIIHa2, elogOIIIHb2)
scatter_plot(logSIIHa2, logOIIIHb2, overplot=True, color = "blue", labels=["$log([SII]/H \\alpha])$", "$log([OIII]/H \\beta])$"])
PBPT(n=2)
plt.text(1, -0.7, 'Error at 75%', ha = 'center')
plt.show()

#OI
add_error_box(logOIHa3, logOIIIHb3, elogOIHa3, elogOIIIHb3)
scatter_plot(logOIHa3, logOIIIHb3, overplot=True, color = "blue", labels=["$log([OI]/H \\alpha])$", "$log([OIII]/H \\beta])$"])
PBPT(n=3)
plt.text(1, -0.7, 'Error at 75%', ha = 'center')
plt.show()
"""


"""
PlotScat(logNIIHa1, logOIIIHb1, ex=elogNIIHa1, ey=elogOIIIHb1, xlim=None, ylim=None, colore="red",
         simbolo="o", labels=["$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PBPT(n=1)


PlotScat(logSIIHa2, logOIIIHb2, ex=elogSIIHa2, ey=elogOIIIHb2, xlim=None, ylim=None, colore="blue",
         simbolo="o", labels=["$log([SII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PBPT(n=2)



PlotScat(logOIHa3, logOIIIHb3, ex=elogOIHa3, ey=elogOIIIHb3, xlim=None, ylim=None, colore="green",
         simbolo="o", labels=["$log([OI]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PBPT(n=3)



"""

if Sampl == 'no':
    if SamePos == "yes":
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal_SamePos.txt"
    else:
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(i1).astype(float), np.array(logNIIHa1).astype(float), np.array(elogNIIHa1).astype(float),
                        np.array(logOIIIHb1).astype(float), np.array(elogOIIIHb1).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)

if Sampl == 'no':
    if SamePos == "yes":
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal_SamePos.txt"
    else:
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(i2).astype(float), np.array(logSIIHa2).astype(float), np.array(elogSIIHa2).astype(float),
                        np.array(logOIIIHb2).astype(float), np.array(elogOIIIHb2).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)

if Sampl == 'no':
    if SamePos == "yes":
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI_gal_SamePos.txt"
    else:
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI_gal.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(i3).astype(float), np.array(logOIHa3).astype(float), np.array(elogOIHa3).astype(float),
                        np.array(logOIIIHb3).astype(float), np.array(elogOIIIHb3).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)


# %% SubSample Radiative, Shock and SF  (R)


#Sampl="no"
#RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl)
Sampl = "no"
#RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl)

#SamePos = "yes"
SamePos = "yes"

def SaveType(i, fileout, arrays):
    val = np.zeros((len(i), len(arrays)))
    kk = 0
    for k in i:
        for t in range(len(arrays)):
            val[kk, t] = arrays[t][k]
        kk += 1
    if os.path.exists(fileout) == False:
        np.savetxt(fileout, val, delimiter='\t', header='\t'.join(
            map(str, range(len(arrays)))), comments='')
    return fileout


indicitot = np.arange(len(OIII_5007))
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII.txt"
i1, x1, ex1, y1, ey1 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII.txt"
i2, x2, ex2, y2, ey2 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI_gal_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI.txt"
i3, x3, ex3, y3, ey3 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)


# BPT NII (Quelli che chiamo in altre parti di codice come Type1)
iAGN = []
icomp = []
ihii1 = []
xAGN = []
exAGN = []
yAGN = []
eyAGN = []
xcomp = []
excomp = []
ycomp = []
eycomp = []
xhii1 = []
exhii1 = []
yhii1 = []
eyhii1 = []

for k in range(len(i1)):
    if (y1[k] >= (0.61 / (x1[k] - 0.47)) + 1.19) or x1[k] >= 0.04:  # RED
        iAGN.append(int(i1[k]))
        xAGN.append(x1[k])
        yAGN.append(y1[k])
        exAGN.append(ex1[k])
        eyAGN.append(ey1[k])
    if (y1[k] < (0.61 / (x1[k] - 0.47)) + 1.19) and (y1[k] >= (0.61/(x1[k] - 0.05)) + 1.3):  # BLUE
        icomp.append(int(i1[k]))
        xcomp.append(x1[k])
        ycomp.append(y1[k])
        excomp.append(ex1[k])
        eycomp.append(ey1[k])
    if y1[k] < (0.61/(x1[k] - 0.05)) + 1.3 and x1[k] < 0.04:  # GREEN
        ihii1.append(int(i1[k]))
        xhii1.append(x1[k])
        yhii1.append(y1[k])
        exhii1.append(ex1[k])
        eyhii1.append(ey1[k])

"""
fract=[]
fract=np.zeros((1000))
for k in range(1000):
    
    fract[k] = len(np.where( (y > str(x)) | x> 0.04 )[0])/len(tot))
    fract.append(len(np.where( (y > str(x)) | x> 0.04 )[0])/len(tot))
    
fractTOT=np.mean(fract)
fractERR=np.std(fract)
"""


PlotScat(xAGN, yAGN, ex=exAGN, ey=eyAGN, xlim=None, ylim=None, colore="red", simbolo="o",
         labels=["$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PlotScat(xcomp, ycomp, ex=excomp, ey=eycomp, xlim=None, ylim=None, colore="blue", simbolo="o", labels=[
         "$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PlotScat(xhii1, yhii1, ex=exhii1, ey=eyhii1, xlim=None, ylim=None, colore="green", simbolo="o", labels=[
         "$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PBPT(n=1)


"""
Versione da mettere nella tesi : ( BPT-NII)
add_error_box(logNIIHa1, logOIIIHb1, elogNIIHa1, elogOIIIHb1, 1, -1)
scatter_plot(xAGN, yAGN, color = 'red', marker = 'o', labels = ["$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], overplot = True)
scatter_plot(xcomp, ycomp, color = 'blue', marker = 'o', labels = ["$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], overplot = True)
scatter_plot(xhii1, yhii1, color = 'green', marker = 'o', labels = ["$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], overplot = True)
PBPT(n = 1)
plt.text(1, -0.7, 'Error at 75%', ha = 'center')
plt.show()
"""

#BPT SII (Quelli che in altre parti di codice chiamo Type2)

irad = []
ishock = []
ihii2 = []
xrad = []
exrad = []
yrad = []
eyrad = []
xshock = []
exshock = []
yshock = []
eyshock = []
xhii2 = []
exhii2 = []
yhii2 = []
eyhii2 = []
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

PlotScat(xrad, yrad, ex=exrad, ey=eyrad, xlim=None, ylim=None, colore="red", simbolo="o",
         labels=["$log([SII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PlotScat(xshock, yshock, ex=exshock, ey=eyshock, xlim=None, ylim=None, colore="blue", simbolo="o", labels=[
         "$log([SII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PlotScat(xhii2, yhii2, ex=exhii2, ey=eyhii2, xlim=None, ylim=None, colore="green", simbolo="o", labels=[
         "$log([SII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PBPT(n=2)


#BPT OI (Quelli che in altre parti di codice chiamo Type3)

irad3 = []
ishock3 = []
ihii3 = []
xrad3 = []
exrad3 = []
yrad3 = []
eyrad3 = []
xshock3 = []
exshock3 = []
yshock3 = []
eyshock3 = []
xhii3 = []
exhii3 = []
yhii3 = []
eyhii3 = []
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

PlotScat(xrad3, yrad3, ex=exrad3, ey=eyrad3, xlim=None, ylim=None, colore="red", simbolo="o",
         labels=["$log([OI]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PlotScat(xshock3, yshock3, ex=exshock3, ey=eyshock3, xlim=None, ylim=None, colore="blue", simbolo="o",
         labels=["$log([OI]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PlotScat(xhii3, yhii3, ex=exhii3, ey=eyhii3, xlim=None, ylim=None, colore="green", simbolo="o", labels=[
         "$log([OI]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PBPT(n=3)


# Save subsample properties


arrays = [RA, DEC, Z, eZ, SIG, eSIG, EBV, Zsun, SIGCLUSTER, NUMGAL, SIGMA_BAL, eSIGMA_BAL, SIGMA_FORB, eSIGMA_FORB,
          VOFF_BAL, eVOFF_BAL, VOFF_FORB, eVOFF_FORB, Mass, eMass1, eMass2, SFR, eSFR1, eSFR2, sSFR, esSFR1, esSFR2]
if Sampl == 'no':
    if SamePos == "yes":
        SaveType(
            iAGN, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal_SamePos.txt", arrays)
        SaveType(
            icomp, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal_SamePos.txt", arrays)
        SaveType(
            ihii1, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal_SamePos.txt", arrays)
        SaveType(
            irad, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD_gal_SamePos.txt", arrays)
        SaveType(
            ishock, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK_gal_SamePos.txt", arrays)
        SaveType(
            ihii2, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2_gal_SamePos.txt", arrays)
        SaveType(
            irad3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD3_gal_SamePos.txt", arrays)
        SaveType(ishock3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK3_gal_SamePos.txt", arrays)
        SaveType(
            ihii3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII3_gal_SamePos.txt", arrays)
    else:
        SaveType(
            iAGN, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal.txt", arrays)
        SaveType(
            icomp, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal.txt", arrays)
        SaveType(
            ihii1, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal.txt", arrays)
        SaveType(
            irad, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD_gal.txt", arrays)
        SaveType(
            ishock, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK_gal.txt", arrays)
        SaveType(
            ihii2, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2_gal.txt", arrays)
        SaveType(
            irad3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD3_gal.txt", arrays)
        SaveType(
            ishock3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK3_gal.txt", arrays)
        SaveType(
            ihii3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII3_gal.txt", arrays)
else:
    SaveType(
        iAGN, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN.txt", arrays)
    SaveType(icomp, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp.txt", arrays)
    SaveType(ihii1, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1.txt", arrays)
    SaveType(
        irad, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD.txt", arrays)
    SaveType(ishock, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK.txt", arrays)
    SaveType(ihii2, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2.txt", arrays)
    SaveType(irad3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD3.txt", arrays)
    SaveType(ishock3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK3.txt", arrays)
    SaveType(ihii3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII3.txt", arrays)


# %% PLOTS BPT subsamples (R)

Sampl = "no"
#RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl)
SamePos = "yes"

if Sampl == 'no':
    if SamePos == "yes":
        f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal_SamePos.txt"
    else:
        f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal.txt"
else:
    f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN.txt"
RA1, DEC1, Z1, eZ1, SIG1, eSIG1, EBV1, Zsun1, SIGCLUSTER, NUMGAL = np.loadtxt(
    f1, usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], unpack=True, dtype=float)
SIGMA_BAL1, eSIGMA_BAL1, SIGMA_FORB1, eSIGMA_FORB1 = np.loadtxt(
    f1, usecols=[8, 9, 10, 11], unpack=True, dtype=float)
VOFF_BAL1, eVOFF_BAL1, VOFF_FORB1, eVOFF_FORB1 = np.loadtxt(
    f1, usecols=[12, 13, 14, 15], unpack=True, dtype=float)
Mass1, eMass11, eMass12, SFR1, eSFR11, eSFR12, sSFR1, esSFR11, esSFR12 = np.loadtxt(
    f1, usecols=[16, 17, 18, 19, 20, 21, 22, 23, 24], unpack=True, dtype=float)

if Sampl == 'no':
    if SamePos == "yes":
        f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal_SamePos.txt"
    else:
        f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal.txt"
else:
    f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp.txt"
RA2, DEC2, Z2, eZ2, SIG2, eSIG2, EBV2, Zsun2 = np.loadtxt(
    f2, usecols=[0, 1, 2, 3, 4, 5, 6, 7], unpack=True, dtype=float)
SIGMA_BAL2, eSIGMA_BAL2, SIGMA_FORB2, eSIGMA_FORB2 = np.loadtxt(
    f2, usecols=[8, 9, 10, 11], unpack=True, dtype=float)
VOFF_BAL2, eVOFF_BAL2, VOFF_FORB2, eVOFF_FORB2 = np.loadtxt(
    f2, usecols=[12, 13, 14, 15], unpack=True, dtype=float)
Mass2, eMass21, eMass22, SFR2, eSFR21, eSFR22, sSFR2, esSFR21, esSFR22 = np.loadtxt(
    f2, usecols=[16, 17, 18, 19, 20, 21, 22, 23, 24], unpack=True, dtype=float)

if Sampl == 'no':
    if SamePos == "yes":
        f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal_SamePos.txt"
    else:
        f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal.txt"
else:
    f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1.txt"
RA3, DEC3, Z3, eZ3, SIG3, eSIG3, EBV3, Zsun3 = np.loadtxt(
    f3, usecols=[0, 1, 2, 3, 4, 5, 6, 7], unpack=True, dtype=float)
SIGMA_BAL3, eSIGMA_BAL3, SIGMA_FORB3, eSIGMA_FORB3 = np.loadtxt(
    f3, usecols=[8, 9, 10, 11], unpack=True, dtype=float)
VOFF_BAL3, eVOFF_BAL3, VOFF_FORB3, eVOFF_FORB3 = np.loadtxt(
    f3, usecols=[12, 13, 14, 15], unpack=True, dtype=float)
Mass3, eMass31, eMass32, SFR3, eSFR31, eSFR32, sSFR3, esSFR31, esSFR32 = np.loadtxt(
    f3, usecols=[16, 17, 18, 19, 20, 21, 22, 23, 24], unpack=True, dtype=float)

# %% Overdensities BCG


f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal_SamePos.txt"
f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal_SamePos.txt"
f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal_SamePos.txt"

RA01, DEC01, Z01, eZ01 = np.loadtxt(
    f1, usecols=[0, 1, 2, 3], unpack=True, dtype=float)
RA02, DEC02, Z02, eZ02 = np.loadtxt(
    f2, usecols=[0, 1, 2, 3], unpack=True, dtype=float)
RA03, DEC03, Z03, eZ03 = np.loadtxt(
    f3, usecols=[0, 1, 2, 3], unpack=True, dtype=float)


f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN.txt"
f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp.txt"
f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1.txt"

RA1, DEC1, Z1, eZ1 = np.loadtxt(
    f1, usecols=[0, 1, 2, 3], unpack=True, dtype=float)
RA2, DEC2, Z2, eZ2 = np.loadtxt(
    f2, usecols=[0, 1, 2, 3], unpack=True, dtype=float)
RA3, DEC3, Z3, eZ3 = np.loadtxt(
    f3, usecols=[0, 1, 2, 3], unpack=True, dtype=float)


RABCG = np.concatenate((RA1, RA2, RA3))
DECBCG = np.concatenate((DEC1, DEC2, DEC3))

RAnoBCG = np.concatenate((RA01, RA02, RA03))
DECnoBCG = np.concatenate((DEC01, DEC02, DEC03))

RAALL = np.concatenate((RABCG, RAnoBCG))
DECALL = np.concatenate((DECBCG, DECnoBCG))

"""
Per capire come si è arrivati a ricavare il subsample SamePos ritagliando lo spazio
in cui comparivano le BCG !!

PlotScat(RA01, DEC01, colore="yellow", simbolo="o")
PlotScat(RA1, DEC1, colore = "red", simbolo="o", overplot="True")
PlotScat(RA2, DEC2, colore="blue", simbolo="o", overplot= "True")
PlotScat(RA3, DEC3, colore="green", simbolo="o", overplot="True")
PlotScat(raR, decR, colore="purple", simbolo="o", overplot="True")

"""

def DistDeg(RA1, DEC1, RAALL, DECALL):
    distGal_AGN = np.zeros((len(RA1), len(RAALL)))
    c1 = SkyCoord(np.asarray(RAALL)*u.deg,
                  np.asarray(DECALL)*u.deg, frame='icrs')
    for i in range(len(RA1)):
        c2 = SkyCoord(np.asarray(RA1[i])*u.deg,
                      np.asarray(DEC1[i])*u.deg, frame='icrs')
        distGal_AGN[i, :] = c1.separation(c2).arcsecond/3600
    return distGal_AGN


dBCGall_AGN = DistDeg(RA1, DEC1, RAALL, DECALL)
dBCGall_Comp = DistDeg(RA2, DEC2, RAALL, DECALL)
dBCGall_HII = DistDeg(RA3, DEC3, RAALL, DECALL)

dGALall_AGN = DistDeg(RA01, DEC01, RAALL, DECALL)
dGALall_Comp = DistDeg(RA02, DEC02, RAALL, DECALL)
dGALall_HII = DistDeg(RA03, DEC03, RAALL, DECALL)


def TCCF(distances_BCG_galaxy, Dmax, label, colore):
    # Definisci i bin radiali
    bins = np.linspace(
        np.min(distances_BCG_galaxy[np.where(distances_BCG_galaxy > 0)]), Dmax, 30)

    # Calcola il numero medio di coppie di BCG e galassie in ciascun bin
    DD, _ = np.histogram(distances_BCG_galaxy, bins=bins)

    # Calcola il numero medio di coppie casuali in ciascun bin (esempio: generato casualmente)
    RR = np.random.uniform(0, 10, 10000)  # 10000 coppie casuali
    RR, _ = np.histogram(RR, bins=bins)

    # Calcola la two-point cross-correlation function
    xi = (DD / RR - 1) / len(distances_BCG_galaxy[:, 0])

    # Visualizza la 2PCF
    plt.plot(bins[:-1], xi, marker='o', label=label, color=colore)
    plt.xlabel('Distanza (r)')
    plt.ylabel('Two-Point Cross-Correlation Function (xi)')
    plt.title('Two-Point Cross-Correlation Function BCG-Galassia')
    plt.legend()
    plt.show()

# Dmax=1
# TCCF(dBCGall_AGN,Dmax,"AGN BCG","red")
# TCCF(dBCGall_Comp,Dmax,"Composite BCG","blue")
# TCCF(dBCGall_HII,Dmax,"HII BCG","green")

# Dmax=1
# TCCF(dGALall_AGN,Dmax,"AGN noBCG","purple")
# TCCF(dGALall_Comp,Dmax,"Composite noBCG","yellow")
# TCCF(dGALall_HII,Dmax,"HII noBCG","black")

# %% Selection RA DEC no BCG


def GalCloseBCG(RA, DEC):
    i = np.where(((RA > 0) & (RA < 50) & (DEC > -10) & (DEC < 15)) | ((RA > 310) & (RA < 360) & (DEC > -10) & (DEC < 15)) | ((RA > 120) & (RA < 250) & (DEC > -5) & (DEC < 12))
                 | ((RA > 111) & (RA < 267) & (DEC > 48) & (DEC < 65)) | ((RA > 111) & (RA < 150) & (DEC > 29) & (DEC < 48)) | ((RA > 231) & (RA < 267) & (DEC > 29) & (DEC < 48)))[0]
    return RA[i], DEC[i], i

def GalCloseBCGFIXED(RA, DEC):
    """
    Fa le stesse cose che faceva la versione precedente questa volta ritagliando
    anche in funzione dello scatterplot degli oggetti analizzati nel paper che ha prodotto
    radioTAB.txt Versione perfezionata un ultima volta da sottoporre ad Andrea...
    """
    
    i = np.where(((RA > -2) & (RA < 53) & (DEC > -10) & (DEC < 4)) | ((RA > 310) & (RA < 360) & (DEC > -10) & (DEC < 4)) | ((RA > 120) & (RA < 250) & (DEC > -5) & (DEC < 6))
                 | ((RA > 111) & (RA < 267) & (DEC > 48) & (DEC < 65)) | ((RA > 111) & (RA < 150) & (DEC > 27) & (DEC < 48)) | ((RA > 222) & (RA < 260) & (DEC > 40) & (DEC < 48))|
                 (RA > 252) & (RA < 267) & (DEC > 25.5) & (DEC < 40))[0]
    return RA[i], DEC[i], i
    

raR, decR = np.loadtxt(
    "/Users/andreamaccarinelli/Desktop/SDSS/RADIOtab.txt", usecols=[3, 4], unpack=True, dtype=float)

f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal_SamePos.txt"
f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal_SamePos.txt"
f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal_SamePos.txt"

RA01, DEC01, Z01, eZ01 = np.loadtxt(
    f1, usecols=[0, 1, 2, 3], unpack=True, dtype=float)
RA02, DEC02, Z02, eZ02 = np.loadtxt(
    f2, usecols=[0, 1, 2, 3], unpack=True, dtype=float)
RA03, DEC03, Z03, eZ03 = np.loadtxt(
    f3, usecols=[0, 1, 2, 3], unpack=True, dtype=float)


f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN.txt"
f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp.txt"
f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1.txt"

RA1, DEC1, Z1, eZ1 = np.loadtxt(
    f1, usecols=[0, 1, 2, 3], unpack=True, dtype=float)
RA2, DEC2, Z2, eZ2 = np.loadtxt(
    f2, usecols=[0, 1, 2, 3], unpack=True, dtype=float)
RA3, DEC3, Z3, eZ3 = np.loadtxt(
    f3, usecols=[0, 1, 2, 3], unpack=True, dtype=float)

"""
#Con la funzione non fixata
RA01, DEC01, i01 = GalCloseBCGFIXED(RA01, DEC01)
RA02, DEC02, i02 = GalCloseBCGFIXED(RA02, DEC02)
RA03, DEC03, i03 = GalCloseBCGFIXED(RA03, DEC03)
"""
"""
#con la funzione si ritaglio FIXATA !!
RA01, DEC01, i01 = GalCloseBCGFIXED(RA01, DEC01)
RA02, DEC02, i02 = GalCloseBCGFIXED(RA02, DEC02)
RA03, DEC03, i03 = GalCloseBCGFIXED(RA03, DEC03)
"""

RABCG = np.concatenate((RA1, RA2, RA3))
DECBCG = np.concatenate((DEC1, DEC2, DEC3))

RAnoBCG = np.concatenate((RA01, RA02, RA03))
DECnoBCG = np.concatenate((DEC01, DEC02, DEC03))

RAALL = np.concatenate((RABCG, RAnoBCG))
DECALL = np.concatenate((DECBCG, DECnoBCG))


#PlotScat Esplorativo

PlotScat(RAnoBCG, DECnoBCG, colore="yellow", simbolo="o", labels = ["RA", "DEC"])
PlotScat(RABCG, DECBCG, colore = "red", simbolo="o", overplot="True", labels = ["RA", "DEC"])
#PlotScat(RA2, DEC2, colore="blue", simbolo="o", overplot= "True")
#PlotScat(RA3, DEC3, colore="green", simbolo="o", overplot="True")
#PlotScat(raR, decR, colore="purple", simbolo="o", overplot="True", labels = ["RA", "DEC"])


"""
#SOLO BCG + paper radioloud
PlotScat(RABCG, DECBCG, colore = "red")
PlotScat(raR, decR, colore="purple", simbolo="o", overplot="True", labels = ["RA", "DEC"])
"""


#noBCGSamePosFixed
#PlotScat(RAnoBCG, DECnoBCG, colore = "blue", overplot="True")


#CON TUTTO
#PlotScat(RAALL, DECALL, colore = "blue", overplot="True")



# %% Fraction AGN (Generali ) (R)

# Fraction AGN in z in the same positions

f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal_SamePos.txt"
Zn1, eZn1 = np.loadtxt(f1, usecols=[2, 3], unpack=True, dtype=float)
f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN.txt"
Z1, eZ1 = np.loadtxt(f1, usecols=[2, 3], unpack=True, dtype=float)


f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal_SamePos.txt"
Zn2, eZn2 = np.loadtxt(f2, usecols=[2, 3], unpack=True, dtype=float)
f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp.txt"
Z2, eZ2 = np.loadtxt(f2, usecols=[2, 3], unpack=True, dtype=float)


f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal_SamePos.txt"
Zn3, eZn3 = np.loadtxt(f3, usecols=[2, 3], unpack=True, dtype=float)
f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1.txt"
Z3, eZ3 = np.loadtxt(f3, usecols=[2, 3], unpack=True, dtype=float)

Ztot = np.concatenate((Z1, Z2, Z3))

zrange = np.arange(np.min(Ztot[Ztot > 0.0001]), np.max(Ztot[Ztot < 2]), 0.05)
zrange = np.logspace(np.log10(np.min(Ztot[Ztot > 0.0001])), np.log10(
    np.max(Ztot[Ztot < 2])), num=12)

FrBCG = np.zeros((len(zrange)-1))
FrnoBCG = np.zeros((len(zrange)-1))
for t in range(len(zrange)-1):
    NAGN = len(np.where((Z1 > zrange[t]) & (Z1 < zrange[t+1]))[0])
    NnoAGN1 = len(np.where((Z2 > zrange[t]) & (Z2 < zrange[t+1]))[0])
    NnoAGN2 = len(np.where((Z3 > zrange[t]) & (Z3 < zrange[t+1]))[0])
    NnoAGN = NnoAGN1+NnoAGN2+NAGN
    if NnoAGN == 0:
        FrBCG[t] = np.nan
    else:
        FrBCG[t] = NAGN/NnoAGN
    NAGN = len(np.where((Zn1 > zrange[t]) & (Zn1 < zrange[t+1]))[0])
    NnoAGN1 = len(np.where((Zn2 > zrange[t]) & (Zn2 < zrange[t+1]))[0])
    NnoAGN2 = len(np.where((Zn3 > zrange[t]) & (Zn3 < zrange[t+1]))[0])
    NnoAGN = NnoAGN1+NnoAGN2+NAGN
    if NnoAGN == 0:
        FrnoBCG[t] = np.nan
    else:
        FrnoBCG[t] = NAGN/NnoAGN

plt.scatter(zrange[1:], FrBCG, color='red', label="BCG", marker="o")
plt.scatter(zrange[1:], FrnoBCG, color='blue', label="no BCG", marker="o")
plt.xlabel("$z$", fontsize=16)
plt.ylabel("$f_{AGN}$", fontsize=16)
plt.tick_params(axis='both', labelsize=16)
plt.subplots_adjust(top=0.850, bottom=0.2, left=0.2,
                    right=0.850, hspace=0.2, wspace=0.2)
plt.legend()


# %% Funzioni necessarie per le celle sottostanti !


def weighted_std(data, weights):
    mean = np.average(data, weights=weights)
    variance = np.average((data - mean)**2, weights=weights)
    return np.sqrt(variance)

# %% Fraction OpticalAGN Type1 esclusivo !!

# Fraction ( NumOpticalAGN ) / ( NumOpticalAGN + Comp + HII1 ) 

"""
Esempio di BootStrap algoritm 1D
fract=[]
fract=np.zeros((1000))
for k in range(1000):
    
    fract[k] = len(np.where( (y > str(x)) | x> 0.04 )[0])/len(tot))
    fract.append(len(np.where( (y > str(x)) | x> 0.04 )[0])/len(tot))
    
fractTOT=np.mean(fract)
fractERR=np.std(fract)
"""


#Lettura dei dati 
Sampl = 'yes'
SamePos = "yes"

if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII.txt"

i1, x1, ex1, y1, ey1 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)


#variazione gaussiana valori_random = np.random.normal( media, devs ) 

Modalita = False

if (Modalita == True ):
    
    DEVST_ASCISSE = weighted_std(x1, 1/ex1**2)
    DEVST_ORDINATE = weighted_std(y1, 1/ey1**2)
    
    #CHECK VISIVO SUI VALORI OTTENUTI
    print ("LA devst pesata delle x è :", DEVST_ASCISSE, "mentre quella sulle y è :", DEVST_ORDINATE, "\n")
    
    
    # Implementazione del Bootstrap Algorithm
    fraction = np.zeros(1000)
    for k in range(1000):
        xVAR = np.random.normal(x1, DEVST_ASCISSE)
        yVAR = np.random.normal(y1, DEVST_ORDINATE)
        fraction[k] = len(np.where((yVAR >= (0.61 / (xVAR - 0.47)) + 1.19) | (xVAR >= 0.04))[0]) / len(i1)
    
    fractTOT = np.mean(fraction)
    fractERR = np.std(fraction)

else :
    # Implementazione del Bootstrap Algorithm
    fraction = np.zeros(1000)
    for k in range(1000):
        # Generazione di punti casuali con errori associati
        xVAR = np.random.normal(x1, ex1)
        yVAR = np.random.normal(y1, ey1)
    
        # Calcolo della frazione in base alla relazione specificata
        fraction[k] = len(np.where((yVAR >= (0.61 / (xVAR - 0.47)) + 1.19) | (xVAR >= 0.04))[0]) / len(i1)

    fractTOT = np.mean(fraction)
    fractERR = np.std(fraction)



if(Sampl == "yes"):
    dicitura = "e BCG"
else:
    if(SamePos == "yes"):
        dicitura = "o spazio intorno alle BCG"
    else:
        dicitura = "insieme delle galassie generico "

testo = "La fraction di OpticalAGN nell" + dicitura
print(testo,"è", fractTOT, "cui è associato un errore di ", fractERR)


# %% Fraction OpticalAGN Type2 esclusivo !!


# Fraction ( RADIATIVE ) / ( RADIATIVE + Shock + HII2 ) 

#Lettura dei dati 
Sampl = 'yes'
SamePos = "yes"

if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII.txt"

i2, x2, ex2, y2, ey2 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

Modalita = False

if(Modalita == True):
    DEVST_ASCISSE = weighted_std(x2, 1/ex2**2)
    DEVST_ORDINATE = weighted_std(y2, 1/ey2**2)

    #CHECK VISIVO SUI VALORI OTTENUTI
    print ("LA devst pesata delle x è :", DEVST_ASCISSE, "mentre quella sulle y è :", DEVST_ORDINATE, "\n")

    # Implementazione del Bootstrap Algorithm
    fraction = np.zeros(1000)
    for k in range(1000):
        xVAR = np.random.normal(x2, DEVST_ASCISSE)
        yVAR = np.random.normal(y2, DEVST_ORDINATE)
        fraction[k] = len(np.where((yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | ( xVAR > 0.29) )[0]) / len(i2)

    fractTOT = np.mean(fraction)
    fractERR = np.std(fraction)
    
    
else:
    
    
    # Implementazione del Bootstrap Algorithm
    fraction = np.zeros(1000)
    for k in range(1000):
        # Generazione di punti casuali con errori associati
        xVAR = np.random.normal(x2, ex2)
        yVAR = np.random.normal(y2, ey2)
    
        # Calcolo della frazione in base alla relazione specificata
        fraction[k] = len(np.where((yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29))[0]) / len(i2)

    fractTOT = np.mean(fraction)
    fractERR = np.std(fraction)


if(Sampl == "yes"):
    dicitura = "e BCG"
else:
    if(SamePos == "yes"):
        dicitura = "o spazio intorno alle BCG"
    else:
        dicitura = "insieme delle galassie generico "

testo = "La fraction di AGN nell" + dicitura
print(testo,"è", fractTOT, "cui è associato un errore di ", fractERR)



# %% Fraction OpticalAGN ==> intersec (Type1, Type2) !! (Versione che dovrebbe essere corretta ! )

#Lettura dei dati 
Sampl = 'no'
SamePos = "yes"

#leggo i file del Type1
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal_SamePos.txt"
        f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal_SamePos.txt"
        f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal_SamePos.txt"
        f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal.txt"
        f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal.txt"
        f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal.txt"
        f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII.txt"
    f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN.txt"
    f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp.txt"
    f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1.txt"

i1, x1, ex1, y1, ey1 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

#leggo i file del type2
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal_SamePos.txt"
        f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD_gal_SamePos.txt"
        f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK_gal_SamePos.txt"
        f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2_gal_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal.txt"
        f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD_gal.txt"
        f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK_gal.txt"
        f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII.txt"
    f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD.txt"
    f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK.txt"
    f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2.txt"

i2, x2, ex2, y2, ey2 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)



# Calcola A ∩ B
intersection_AB = np.intersect1d(i1, i2)


#print("Abbiamo a disposizione con il corrente campione di calcolo",len(result), "elementi\n")
print("L'intersezione è composta da un numero di ", len(intersection_AB), "Elementi")

# Usa np.where per ottenere gli indici comuni
elementi_comuni_i1 = i1[np.where(np.isin(i1, intersection_AB))]
elementi_comuni_i2 = i2[np.where(np.isin(i2, intersection_AB))]

#ricavo i relativi elementi x,y facenti capo a questo indice di elementi in comune
x1_selected = x1[np.where(np.isin(np.arange(len(x1)), elementi_comuni_i1))]
y1_selected = y1[np.where(np.isin(np.arange(len(y1)), elementi_comuni_i1))]
x2_selected = x2[np.where(np.isin(np.arange(len(x2)), elementi_comuni_i2))]
y2_selected = y2[np.where(np.isin(np.arange(len(y2)), elementi_comuni_i2))]
ex1_selected = ex1[np.where(np.isin(np.arange(len(x1)), elementi_comuni_i1))]
ex2_selected = ex2[np.where(np.isin(np.arange(len(x2)), elementi_comuni_i2))]
ey1_selected = ey1[np.where(np.isin(np.arange(len(y1)), elementi_comuni_i1))]
ey2_selected = ey2[np.where(np.isin(np.arange(len(x2)), elementi_comuni_i2))]

#Contributo dal type1
fraction1 = np.zeros(5000)
fraction2 = np.zeros(5000)
fraction3 = np.zeros(5000)
NAGN = np.zeros(5000)
for k in range(5000):
    xVAR = np.random.normal(x1_selected, ex1_selected)
    yVAR = np.random.normal(y1_selected, ey1_selected)
    fraction1[k] = len(np.where((yVAR >= (0.61 / (xVAR - 0.47)) + 1.19) | (xVAR >= 0.04))[0]) / len(intersection_AB) #f_1a
    NAGN[k] = len(np.where((yVAR >= (0.61 / (xVAR - 0.47)) + 1.19) | (xVAR >= 0.04))[0])

for k in range(5000):
    xVAR = np.random.normal(x1_selected, ex1_selected)
    yVAR = np.random.normal(y1_selected, ey1_selected)
    fraction2[k] = len(np.where(  yVAR < (0.61 / (xVAR - 0.47)) + 1.19) & (yVAR >= (0.61/(xVAR - 0.05)) + 1.3  )[0])/len(intersection_AB) #f_1b

for k in range(5000):
    xVAR = np.random.normal(x1_selected, ex1_selected)
    yVAR = np.random.normal(y1_selected, ey1_selected)
    fraction3[k] = len(np.where( ( yVAR < (0.61/(xVAR - 0.05)) + 1.3) & (xVAR < 0.04 ) )[0])/len(intersection_AB) #f_1c

fractTOT_1 = np.mean(fraction1)
fractTOT_2 = np.mean(fraction2)
fractTOT_3 = np.mean(fraction3)
fractERR_1 = np.std(fraction1)
fractERR_2 = np.std(fraction2)
fractERR_3 = np.std(fraction3)
sum1 = fractTOT_1 + fractTOT_2 + fractTOT_3
esum1 = (fractERR_1 + fractERR_2 + fractERR_3)

#Contributo dal type2
fraction4 = np.zeros(5000)
for k in range(5000):
    xVAR = np.random.normal(x2_selected, ex2_selected)
    yVAR = np.random.normal(y2_selected, ey2_selected)
    fraction4[k] = len(np.where((yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29))[0]) / len(x2_selected)

# Correggi il nome della variabile usata per calcolare la media e la deviazione standard
fractTOT_4 = np.mean(fraction4)
fractERR_4 = np.std(fraction4)

#Sospetta frazione complessiva 

Frazione_Complex = fractTOT_1 + fractTOT_2
errFrazione_Complex = fractERR_1 + fractERR_2

#mediapesata = ((fractTOT_1 * len(elementi_comuni_i1)) + (fractTOT_2 * len(elementi_comuni_i2)))  / (len(elementi_comuni_i1) + len(elementi_comuni_i2))

print("Ho ottenuto i seguenti risultati : \n")
print("Identificazione con BPT-NII \n", fractTOT_1, fractERR_1, "\n")
print("f1a =", fractTOT_1, "+-", fractERR_1, "\n")
print("f1b =", fractTOT_2, "+-", fractERR_2, "\n")
print("f1c =", fractTOT_3, "+-", fractERR_3, "\n")
print("Normalizzazione :", sum1, esum1 )


print("Identificazione con BPT-SII", fractTOT_4, fractERR_4, "\n")


# %% Versione Funzionante !!!

# Lettura dei dati
Sampl = 'no'
SamePos = "yes"

# Leggo i file del Type1
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal_SamePos.txt"
        
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal.txt"
        
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII.txt"
    

i1, x1, ex1, y1, ey1 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

# Leggo i file del type2
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal_SamePos.txt"
        f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD_gal_SamePos.txt"
        f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK_gal_SamePos.txt"
        f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2_gal_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal.txt"
        f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD_gal.txt"
        f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK_gal.txt"
        f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII.txt"
    f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD.txt"
    f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK.txt"
    f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2.txt"

i2, x2, ex2, y2, ey2 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

# Calcola A ∩ B
intersection_AB = np.intersect1d(i1, i2)

# Stampa il numero di elementi nell'intersezione
print("L'intersezione è composta da un numero di", len(intersection_AB), "Elementi")

# Usa np.where per ottenere gli indici comuni
elementi_comuni_i1 = np.where(np.isin(i1, intersection_AB))[0]
elementi_comuni_i2 = np.where(np.isin(i2, intersection_AB))[0]

# Ricavo i relativi elementi x, y facenti capo a questo indice di elementi in comune
x1_selected = x1[elementi_comuni_i1]
y1_selected = y1[elementi_comuni_i1]
x2_selected = x2[elementi_comuni_i2]
y2_selected = y2[elementi_comuni_i2]
ex1_selected = ex1[elementi_comuni_i1]
ex2_selected = ex2[elementi_comuni_i2]
ey1_selected = ey1[elementi_comuni_i1]
ey2_selected = ey2[elementi_comuni_i2]

# Contributo dal type1
fraction1 = np.zeros(5000)
fraction2 = np.zeros(5000)
fraction3 = np.zeros(5000)
NumOgg1  = np.zeros(5000)
NumOgg2 = np.zeros(5000)
NumOgg3 = np.zeros(5000)

for k in range(5000):
    xVAR = np.random.normal(x1_selected, ex1_selected)
    yVAR = np.random.normal(y1_selected, ey1_selected)
    condition_a = (yVAR >= (0.61 / (xVAR - 0.47)) + 1.19) | (xVAR >= 0.04)
    condition_b = (yVAR < (0.61 / (xVAR - 0.47)) + 1.19) & (yVAR >= (0.61 / (xVAR - 0.05)) + 1.3)
    condition_c = (yVAR < (0.61 / (xVAR - 0.05)) + 1.3) & (xVAR < 0.04)
    fraction1[k] = len(np.where(condition_a)[0]) / len(intersection_AB)
    fraction2[k] = len(np.where(condition_b)[0]) / len(intersection_AB)
    fraction3[k] = len(np.where(condition_c)[0]) / len(intersection_AB)
    NumOgg1[k] = len(np.where(condition_a)[0]) 
    NumOgg2[k] = len(np.where(condition_b)[0])
    NumOgg3[k] = len(np.where(condition_c)[0]) 
    
    

fractTOT_1 = np.mean(fraction1)  #f_1a
fractTOT_2 = np.mean(fraction2)  #f_1b
fractTOT_3 = np.mean(fraction3)  #f_1c
fractERR_1 = np.std(fraction1)
fractERR_2 = np.std(fraction2)
fractERR_3 = np.std(fraction3)
sum1 = fractTOT_1 + fractTOT_2 + fractTOT_3
esum1 = fractERR_1 + fractERR_2 + fractERR_3
N1A = int(np.mean(NumOgg1))
N1B = int(np.mean(NumOgg2))
N1C = int(np.mean(NumOgg3))
eN1A = int(np.std(NumOgg1))
eN1B = int(np.std(NumOgg2))
eN1C = int(np.std(NumOgg3))


# Contributo dal type2
fraction4 = np.zeros(5000)
fraction5 = np.zeros(5000)
fraction6 = np.zeros(5000)
NumOgg4  = np.zeros(5000)
NumOgg5 = np.zeros(5000)
NumOgg6 = np.zeros(5000)
for k in range(5000):
    
    #Creazione delle distribuzioni
    xVAR = np.random.normal(x2_selected, ex2_selected)
    yVAR = np.random.normal(y2_selected, ey2_selected)
    
    #Condizioni di collocamento dei punti
    condition_d = (yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29)
    condition_d = condition_d & (yVAR >= 1.89 * xVAR + 0.76)
    condition_e = (yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29)
    condition_e = condition_e & (yVAR < 1.89 * xVAR + 0.76)
    condition_f = ~((yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29))

    #Popolamento degli array di fraction
    fraction4[k] = len(np.where(condition_d)[0]) / len(x2_selected)
    fraction5[k] = len(np.where(condition_e)[0]) / len(x2_selected)
    fraction6[k] = len(np.where(condition_f)[0]) / len(x2_selected)
    
    #Popolamento degli array di conteggio
    NumOgg4[k] = len(np.where(condition_d)[0]) 
    NumOgg5[k] = len(np.where(condition_e)[0])
    NumOgg6[k] = len(np.where(condition_f)[0]) 

fractTOT_4 = np.mean(fraction4) #f_2a
fractERR_4 = np.std(fraction4)

fractTOT_5 = np.mean(fraction5) #f_2b
fractERR_5 = np.std(fraction5)

fractTOT_6 = np.mean(fraction6) #f_2c
fractERR_6 = np.std(fraction6)

sum2 = fractTOT_4 + fractTOT_5 + fractTOT_6
esum2 = fractERR_4 + fractERR_5 + fractERR_6

N2A = int(np.mean(NumOgg4))
N2B = int(np.mean(NumOgg5))
N2C = int(np.mean(NumOgg6))
eN2A = int(np.std(NumOgg4))
eN2B = int(np.std(NumOgg5))
eN2C = int(np.std(NumOgg6))


# Stampa i risultati
print("Ho ottenuto i seguenti risultati:\n")
print("Identificazione con BPT-NII\n")
print("AGN =", N1A, "+-", eN1A, "\n")
print("AGN =", fractTOT_1, "+-", fractERR_1, "\n")
print("COMPOSITE =", N1B, "+-", eN1B, "\n")
print("COMPOSITE =", fractTOT_2, "+-", fractERR_2, "\n")
print("SF Galaxies =", N1C, "+-", eN1C, "\n")
print("SF Galaxies =", fractTOT_3, "+-", fractERR_3, "\n")
print("Check Normalizzazione:", sum1,"+-" ,  esum1)



# Verifica l'unione delle condizioni
all_points_covered = np.logical_or.reduce([condition_a, condition_b, condition_c])

# Stampa il numero di punti coperti e il totale dei punti
print("Numero di punti coperti:", np.sum(all_points_covered))
print("Totale dei punti:", len(all_points_covered),"\n\n\n")



print("Identificazione con BPT-SII\n")
print("RADIATIVE =", N2A, "+-", eN2A, "\n")
print("RADIATIVE =", fractTOT_4, "+-", fractERR_4, "\n")
print("SHOCK =", N2B, "+-", eN2B, "\n")
print("SHOCK =", fractTOT_5, "+-", fractERR_5, "\n")
print("SF Galaxies =", N2C, "+-", eN2C, "\n")
print("SF Galaxies =", fractTOT_6, "+-", fractERR_6, "\n\n")
print("Check Normalizzazione:", sum2, "+-",esum2)

# Verifica l'unione delle condizioni
all_points_covered = np.logical_or.reduce([condition_d, condition_e, condition_f])

# Stampa il numero di punti coperti e il totale dei punti
print("Numero di punti coperti:", np.sum(all_points_covered))
print("Totale dei punti:", len(all_points_covered))


# %% SENZA INTERSEZIONE

# Lettura dei dati
Sampl = 'yes'
SamePos = "yes"     # Far Girare SOLAMENTE CON SAMEPOS = yes (Manca l'implementazione corretta del percorso !!)

# Leggo i file del Type1
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal_SamePos.txt"
        
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal.txt"
        
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII.txt"
    

i1, x1, ex1, y1, ey1 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

# Leggo i file del type2
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal_SamePos.txt"
        
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal.txt"
        
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII.txt"
    

i2, x2, ex2, y2, ey2 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)


print("Sull'NII abbiamo", len(i1), "Elementi, mantre nel SII abbiamo", len(i2), "Elementi")

# Contributo dal type1
fraction1 = np.zeros(5000)
fraction2 = np.zeros(5000)
fraction3 = np.zeros(5000)
NumOgg1  = np.zeros(5000)
NumOgg2 = np.zeros(5000)
NumOgg3 = np.zeros(5000)

for k in range(5000):
    xVAR = np.random.normal(x1, ex1)
    yVAR = np.random.normal(y1, ey1)
    condition_a = (yVAR >= (0.61 / (xVAR - 0.47)) + 1.19) | (xVAR >= 0.04)
    condition_b = (yVAR < (0.61 / (xVAR - 0.47)) + 1.19) & (yVAR >= (0.61 / (xVAR - 0.05)) + 1.3)
    condition_c = (yVAR < (0.61 / (xVAR - 0.05)) + 1.3) & (xVAR < 0.04)
    fraction1[k] = len(np.where(condition_a)[0]) / len(i1)
    fraction2[k] = len(np.where(condition_b)[0]) / len(i1)
    fraction3[k] = len(np.where(condition_c)[0]) / len(i1)
    NumOgg1[k] = len(np.where(condition_a)[0]) 
    NumOgg2[k] = len(np.where(condition_b)[0])
    NumOgg3[k] = len(np.where(condition_c)[0]) 
    
    

fractTOT_1 = np.mean(fraction1)  #f_1a
fractTOT_2 = np.mean(fraction2)  #f_1b
fractTOT_3 = np.mean(fraction3)  #f_1c
fractERR_1 = np.std(fraction1)
fractERR_2 = np.std(fraction2)
fractERR_3 = np.std(fraction3)
sum1 = fractTOT_1 + fractTOT_2 + fractTOT_3
esum1 = fractERR_1 + fractERR_2 + fractERR_3
N1A = int(np.mean(NumOgg1))
N1B = int(np.mean(NumOgg2))
N1C = int(np.mean(NumOgg3))
eN1A = int(np.std(NumOgg1))
eN1B = int(np.std(NumOgg2))
eN1C = int(np.std(NumOgg3))

# Contributo dal type2
fraction4 = np.zeros(5000)
fraction5 = np.zeros(5000)
fraction6 = np.zeros(5000)
NumOgg4  = np.zeros(5000)
NumOgg5 = np.zeros(5000)
NumOgg6 = np.zeros(5000)
for k in range(5000):
    
    #Creazione delle distribuzioni
    xVAR = np.random.normal(x2, ex2)
    yVAR = np.random.normal(y2, ey2)
    
    #Condizioni di collocamento dei punti
    condition_d = (yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29)
    condition_d = condition_d & (yVAR >= 1.89 * xVAR + 0.76)
    condition_e = (yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29)
    condition_e = condition_e & (yVAR < 1.89 * xVAR + 0.76)
    condition_f = ~((yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29))

    #Popolamento degli array di fraction
    fraction4[k] = len(np.where(condition_d)[0]) / len(i2)
    fraction5[k] = len(np.where(condition_e)[0]) / len(i2)
    fraction6[k] = len(np.where(condition_f)[0]) / len(i2)
    
    #Popolamento degli array di conteggio
    NumOgg4[k] = len(np.where(condition_d)[0]) 
    NumOgg5[k] = len(np.where(condition_e)[0])
    NumOgg6[k] = len(np.where(condition_f)[0]) 

fractTOT_4 = np.mean(fraction4) #f_2a
fractERR_4 = np.std(fraction4)

fractTOT_5 = np.mean(fraction5) #f_2b
fractERR_5 = np.std(fraction5)

fractTOT_6 = np.mean(fraction6) #f_2c
fractERR_6 = np.std(fraction6)

sum2 = fractTOT_4 + fractTOT_5 + fractTOT_6
esum2 = fractERR_4 + fractERR_5 + fractERR_6

N2A = int(np.mean(NumOgg4))
N2B = int(np.mean(NumOgg5))
N2C = int(np.mean(NumOgg6))
eN2A = int(np.std(NumOgg4))
eN2B = int(np.std(NumOgg5))
eN2C = int(np.std(NumOgg6))


# Stampa i risultati
print("Ho ottenuto i seguenti risultati, Senza considerare l'insieme intersezione !!!':\n")
print("Identificazione con BPT-NII\n")
print("AGN =", N1A, "+-", eN1A, "\n")
print("AGN =", fractTOT_1, "+-", fractERR_1, "\n")
print("COMPOSITE =", N1B, "+-", eN1B, "\n")
print("COMPOSITE =", fractTOT_2, "+-", fractERR_2, "\n")
print("SF Galaxies =", N1C, "+-", eN1C, "\n")
print("SF Galaxies =", fractTOT_3, "+-", fractERR_3, "\n")
print("Check Normalizzazione:", sum1,"+-" ,  esum1)



# Verifica l'unione delle condizioni
all_points_covered = np.logical_or.reduce([condition_a, condition_b, condition_c])

# Stampa il numero di punti coperti e il totale dei punti
print("Numero di punti coperti:", np.sum(all_points_covered))
print("Totale dei punti:", len(all_points_covered),"\n\n\n")



print("Identificazione con BPT-SII\n")
print("RADIATIVE =", N2A, "+-", eN2A, "\n")
print("RADIATIVE =", fractTOT_4, "+-", fractERR_4, "\n")
print("SHOCK =", N2B, "+-", eN2B, "\n")
print("SHOCK =", fractTOT_5, "+-", fractERR_5, "\n")
print("SF Galaxies =", N2C, "+-", eN2C, "\n")
print("SF Galaxies =", fractTOT_6, "+-", fractERR_6, "\n\n")
print("Check Normalizzazione:", sum2, "+-",esum2)

# Verifica l'unione delle condizioni
all_points_covered = np.logical_or.reduce([condition_d, condition_e, condition_f])

# Stampa il numero di punti coperti e il totale dei punti
print("Numero di punti coperti:", np.sum(all_points_covered))
print("Totale dei punti:", len(all_points_covered))



# %% Percentuale di RadioLoud maggiore o minore nelle BCG ?

f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/infoRadioBCG.txt"
f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoRadionoBCG_SamePos.txt"



sort, radioact1 = np.loadtxt(f1,usecols=[1,3], unpack=True, dtype=float)
sort2, radioact2 = np.loadtxt(f2,usecols=[1,3], unpack=True, dtype=float)


#Cerco le posizioni in cui c'è una BCG
Bool_esito = np.nonzero(sort)[0]
Bool_esito2 = np.nonzero(sort2)[0]

#print(len(Bool_esito))
RLtemp_1 = radioact1[np.where(np.isin(np.arange(len(radioact1)), Bool_esito))]
RLtemp_2 = radioact2[np.where(np.isin(np.arange(len(radioact2)), Bool_esito2))]

RL_1 = np.count_nonzero(RLtemp_1 == 1)
RL_2 = np.count_nonzero(RLtemp_2 == 1)

#Calcolo dei rapporti
rapporto1, rapporto2 = RL_1/len(Bool_esito) , RL_2 / len(Bool_esito2)
print(rapporto1 *100, rapporto2*100)
#Osservo una percentuale di RLQ molto maggiore nel caso delle BCG 98% mentre nel SamePos è 78.8%

#Inverto ora la direzione del crossmatch per ricavare gli indici RadioLoud nel SDSS ricavato !


#Leggo le Informazioni relative al paper delle RadioEmitters
raR, decR = np.loadtxt("/Users/andreamaccarinelli/Desktop/SDSS/RADIOtab.txt", usecols=[3, 4], unpack=True, dtype=float)

#Leggo le coordinate dai file delle BCG e delle noBCG_SamePos
f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoBCG.txt"
f4 = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoNoBCG_SamePos.txt"
raB, decB = np.loadtxt(f3, usecols=[0, 1], unpack=True, dtype=float)
raS, decS = np.loadtxt(f4, usecols=[0, 1], unpack=True, dtype=float)

#Ottengo gli indici delle posizioni in cui gli oggetti sono studiati e sono RadioLoud 
index1 = np.where((radioact1 != 0) & (sort != 0) )[0]
index2 = np.where((radioact2 != 0) & (sort2 != 0))[0]


#Ottengo le selezioni di coordinate secondo il file del paper
RA_1 = raR[index1]
DEC_1 = decR[index1]
RA_2 = raR[index2]
DEC_2 = decR[index2]

#Per le BCG RadioLoud
c1 = SkyCoord(RA_1*u.deg, DEC_1*u.deg, frame='icrs')
RadioLoud = []
for k in range(len(raB)):
    c2 = SkyCoord(raB[k]*u.deg, decB[k]*u.deg, frame='icrs')
    dist = c2.separation(c1).to(u.arcsec).value
    if np.min(dist) <= 5:
        RadioLoud.append(1)
    else:
        RadioLoud.append(0)

# Leggi tutte le colonne esistenti dal file InfoBCG.txt
existing_data = np.loadtxt(f3, usecols=(0,1,2,3,4,5,6,7,8,9,10), dtype=float)
# Aggiungi la nuova colonna RadioLoud a existing_data
updated_data = np.column_stack((existing_data, RadioLoud))
# Sovrascrivi il file InfoBCG.txt con i nuovi dati
np.savetxt(f3, updated_data, fmt='%f', comments='')

#Per le noBCG_SamePos RadioLoud
c1 = SkyCoord(RA_2*u.deg, DEC_2*u.deg, frame='icrs')
RadioLoud = []
for k in range(len(raS)):
    c2 = SkyCoord(raS[k]*u.deg, decS[k]*u.deg, frame='icrs')
    dist = c2.separation(c1).to(u.arcsec).value
    if np.min(dist) <= 5:
        RadioLoud.append(1)
    else:
        RadioLoud.append(0)

# Leggi tutte le colonne esistenti dal file InfoBCG.txt
existing_data = np.loadtxt(f4, usecols=(0,1,2,3,4,5,6,7,8), dtype=float)
# Aggiungi la nuova colonna RadioLoud a existing_data
updated_data = np.column_stack((existing_data, RadioLoud))
# Sovrascrivi il file InfoBCG.txt con i nuovi dati
np.savetxt(f4, updated_data, fmt='%f', comments='')

# %% Creazione dei Plot riferiti agli oggetti RadioLoud (Collaudato e Funziona !!!)

#Step Zero: Lettura dei dati "Non includo nella trattazione il caso SamePos = "no" 

#Sampl = "yes"
Sampl = "no"
SamePos= "yes"

if (Sampl == "no"):
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoNoBCG_SamePos.txt"
    status = np.loadtxt(f, usecols=[9], unpack=True, dtype=float)
    
else :
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoBCG.txt"
    status = np.loadtxt(f, usecols=[11], unpack=True, dtype=float)
    

#Ricerca degli indici in cui compaiono solo RadioLoud 
index = np.where(status != 0)[0]

#Continua da qui, devi richiamare i dati e poi seguire paro paro quanto facciamo qualche cella sopra
#per creare i BPT diagrams !!!


#calldata(Sampl, "yes", index)
RA, DEC, Z, eZ, SIG, eSIG, EBV, Zsun, SIGCLUSTER, NUMGAL, SIGMA_BAL, eSIGMA_BAL, SIGMA_FORB, eSIGMA_FORB, VOFF_BAL, eVOFF_BAL, VOFF_FORB, eVOFF_FORB, OII_3726, eOII_3726, OII_3729, eOII_3729, NEIII_3869, eNEIII_3869, H_DELTA, eH_DELTA, H_GAMMA, eH_GAMMA, OIII_4363, eOIII_4363, OIII_4959, eOIII_4959, OIII_5007, eOIII_5007, HEI_5876, eHEI_5876, OI_6300, eOI_6300, H_BETA, eH_BETA, H_ALPHA, eH_ALPHA, NII_6584, eNII_6584, SII_6717, eSII_6717, SII_6731, eSII_6731, ARIII7135, eARIII7135, Mass, eMass1, eMass2, SFR, eSFR1, eSFR2, sSFR, esSFR1, esSFR2 = calldata(Sampl, "yes", index)



# REMOVE FALSE VALUES
indicitot = np.arange(len(OIII_5007))
indexes1 = np.where((OIII_5007 > 0) & (H_ALPHA > 0) & (H_BETA > 0) & (NII_6584 > 0) & (
    OIII_5007 > eOIII_5007) & (H_ALPHA > eH_ALPHA) & (H_BETA > eH_BETA) & (NII_6584 > eNII_6584))[0]
indexes2 = np.where((OIII_5007 > 0) & (H_ALPHA > 0) & (H_BETA > 0) & (SII_6717 > 0) & (SII_6731 > 0) & (
    OIII_5007 > eOIII_5007) & (H_ALPHA > eH_ALPHA) & (H_BETA > eH_BETA) & (SII_6717 > eSII_6717) & (SII_6731 > eSII_6731))[0]
indexes3 = np.where((OIII_5007 > 0) & (H_ALPHA > 0) & (H_BETA > 0) & (OI_6300 > 0) & (
    OIII_5007 > eOIII_5007) & (H_ALPHA > eH_ALPHA) & (H_BETA > eH_BETA) & (OI_6300 > eOI_6300))[0]

logOIIIHb1 = np.log10(OIII_5007[indexes1]/H_BETA[indexes1])
elogOIIIHb1 = ErrLogRatio(OIII_5007[indexes1], H_BETA[indexes1],
                          err_num=eOIII_5007[indexes1], err_den=eH_BETA[indexes1])
logNIIHa1 = np.log10(NII_6584[indexes1]/H_ALPHA[indexes1])
elogNIIHa1 = ErrLogRatio(NII_6584[indexes1], H_ALPHA[indexes1],
                         err_num=eNII_6584[indexes1], err_den=eH_ALPHA[indexes1])


i1 = indicitot[indexes1]   #Indici dell'array con tutte le galassie su cui è pox fare BPTtype1 



logOIIIHb2 = np.log10(OIII_5007[indexes2]/H_BETA[indexes2])
elogOIIIHb2 = ErrLogRatio(OIII_5007[indexes2], H_BETA[indexes2],
                          err_num=eOIII_5007[indexes2], err_den=eH_BETA[indexes2])
SII = SII_6731[indexes2]
eSII = eSII_6731[indexes2]
logSIIHa2 = np.log10((SII)/H_ALPHA[indexes2])
elogSIIHa2 = ErrLogRatio(
    SII, H_ALPHA[indexes2], err_num=eSII, err_den=eH_ALPHA[indexes2])


i2 = indicitot[indexes2]    #Indici dell'array con tutte le galassie su cui è pox fare BPTtype2 

logOIIIHb3 = np.log10(OIII_5007[indexes3]/H_BETA[indexes3])
elogOIIIHb3 = ErrLogRatio(OIII_5007[indexes3], H_BETA[indexes3],
                          err_num=eOIII_5007[indexes3], err_den=eH_BETA[indexes3])
logOIHa3 = np.log10(OI_6300[indexes3]/H_ALPHA[indexes3])
elogOIHa3 = ErrLogRatio(OI_6300[indexes3], H_ALPHA[indexes3],
                        err_num=eOI_6300[indexes3], err_den=eH_ALPHA[indexes3])


i3 = indicitot[indexes3]   #Indici dell'array con tutte le galassie su cui è pox fare BPTtype3 



"""
#NII
add_error_box(logNIIHa1, logOIIIHb1, elogNIIHa1, elogOIIIHb1, 1, -1)
scatter_plot(logNIIHa1, logOIIIHb1, overplot=True, color = "blue", labels=["$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"])
PBPT(n=1)
plt.text(1, -0.7, 'Error at 75%', ha = 'center')
plt.show()

#SII
add_error_box(logSIIHa2, logOIIIHb2, elogSIIHa2, elogOIIIHb2)
scatter_plot(logSIIHa2, logOIIIHb2, overplot=True, color = "blue", labels=["$log([SII]/H \\alpha])$", "$log([OIII]/H \\beta])$"])
PBPT(n=2)
plt.text(1, -0.7, 'Error at 75%', ha = 'center')
plt.show()

#OI
add_error_box(logOIHa3, logOIIIHb3, elogOIHa3, elogOIIIHb3)
scatter_plot(logOIHa3, logOIIIHb3, overplot=True, color = "blue", labels=["$log([OI]/H \\alpha])$", "$log([OIII]/H \\beta])$"])
PBPT(n=3)
plt.text(1, -0.7, 'Error at 75%', ha = 'center')
plt.show()
"""



"""
PlotScat(logNIIHa1, logOIIIHb1, ex=elogNIIHa1, ey=elogOIIIHb1, xlim=None, ylim=None, colore="red",
         simbolo="o", labels=["$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PBPT(n=1)

PlotScat(logSIIHa2, logOIIIHb2, ex=elogSIIHa2, ey=elogOIIIHb2, xlim=None, ylim=None, colore="blue",
         simbolo="o", labels=["$log([SII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PBPT(n=2)

PlotScat(logOIHa3, logOIIIHb3, ex=elogOIHa3, ey=elogOIIIHb3, xlim=None, ylim=None, colore="green",
         simbolo="o", labels=["$log([OI]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PBPT(n=3)
"""

if Sampl == 'no':
    if SamePos == "yes":
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal_RadioLoud_SamePos.txt"
    else:
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_RadioLoud_gal.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII-RadioLoud.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(i1).astype(float), np.array(logNIIHa1).astype(float), np.array(elogNIIHa1).astype(float),
                        np.array(logOIIIHb1).astype(float), np.array(elogOIIIHb1).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)

if Sampl == 'no':
    if SamePos == "yes":
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal_RadioLoud_SamePos.txt"
    else:
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_RadioLoud_gal.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII-RadioLoud.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(i2).astype(float), np.array(logSIIHa2).astype(float), np.array(elogSIIHa2).astype(float),
                        np.array(logOIIIHb2).astype(float), np.array(elogOIIIHb2).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)

if Sampl == 'no':
    if SamePos == "yes":
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI_gal_RadioLoud_SamePos.txt"
    else:
        fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI_RadioLoud_gal.txt"
else:
    fileout = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI-RadioLoud.txt"
fmt = "%f"  # Specify the format string
data = np.column_stack((np.array(i3).astype(float), np.array(logOIHa3).astype(float), np.array(elogOIHa3).astype(float),
                        np.array(logOIIIHb3).astype(float), np.array(elogOIIIHb3).astype(float)))
if os.path.exists(fileout) == False:
    np.savetxt(fileout, data, fmt=fmt)

# %% Plot BPT SUBSAMPLES RADIO LOUD Radiative, Shock and SF  (R)


#Sampl="no"
#RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl)
Sampl = "no"
#RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl)

#SamePos = "yes"
SamePos = "yes"

def SaveType(i, fileout, arrays):
    """
    Versione Precedente che dava un messaggio di errore quando fatta girare
    con Sampl = 'no'
    
    Parameters
    ----------
    i : TYPE
        DESCRIPTION.
    fileout : TYPE
        DESCRIPTION.
    arrays : TYPE
        DESCRIPTION.

    Returns
    -------
    fileout : TYPE
        DESCRIPTION.

    """
    val = np.zeros((len(i), len(arrays)))
    kk = 0
    for k in i:
        for t in range(len(arrays)):
            val[kk, t] = arrays[t][k]
        kk += 1
    if os.path.exists(fileout) == False:
        np.savetxt(fileout, val, delimiter='\t', header='\t'.join(
            map(str, range(len(arrays)))), comments='')
    return fileout




indicitot = np.arange(len(OIII_5007))
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal_RadioLoud_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_RadioLoud_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII-RadioLoud.txt"
i1, x1, ex1, y1, ey1 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)


if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal_RadioLoud_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_RadioLoud_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII-RadioLoud.txt"
i2, x2, ex2, y2, ey2 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI_gal_RadioLoud_SamePos.txt"
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI_RadioLoud_gal.txt"
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-OI-RadioLoud.txt"
i3, x3, ex3, y3, ey3 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)


# BPT NII (Quelli che chiamo in altre parti di codice come Type1)
iAGN = []
icomp = []
ihii1 = []
xAGN = []
exAGN = []
yAGN = []
eyAGN = []
xcomp = []
excomp = []
ycomp = []
eycomp = []
xhii1 = []
exhii1 = []
yhii1 = []
eyhii1 = []

for k in range(len(i1)):
    if (y1[k] >= (0.61 / (x1[k] - 0.47)) + 1.19) or x1[k] >= 0.04:  # RED
        iAGN.append(int(i1[k]))
        xAGN.append(x1[k])
        yAGN.append(y1[k])
        exAGN.append(ex1[k])
        eyAGN.append(ey1[k])
    if (y1[k] < (0.61 / (x1[k] - 0.47)) + 1.19) and (y1[k] >= (0.61/(x1[k] - 0.05)) + 1.3):  # BLUE
        icomp.append(int(i1[k]))
        xcomp.append(x1[k])
        ycomp.append(y1[k])
        excomp.append(ex1[k])
        eycomp.append(ey1[k])
    if y1[k] < (0.61/(x1[k] - 0.05)) + 1.3 and x1[k] < 0.04:  # GREEN
        ihii1.append(int(i1[k]))
        xhii1.append(x1[k])
        yhii1.append(y1[k])
        exhii1.append(ex1[k])
        eyhii1.append(ey1[k])

"""
fract=[]
fract=np.zeros((1000))
for k in range(1000):
    
    fract[k] = len(np.where( (y > str(x)) | x> 0.04 )[0])/len(tot))
    fract.append(len(np.where( (y > str(x)) | x> 0.04 )[0])/len(tot))
    
fractTOT=np.mean(fract)
fractERR=np.std(fract)
"""


PlotScat(xAGN, yAGN, ex=exAGN, ey=eyAGN, xlim=None, ylim=None, colore="red", simbolo="o",
         labels=["$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PlotScat(xcomp, ycomp, ex=excomp, ey=eycomp, xlim=None, ylim=None, colore="blue", simbolo="o", labels=[
         "$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PlotScat(xhii1, yhii1, ex=exhii1, ey=eyhii1, xlim=None, ylim=None, colore="green", simbolo="o", labels=[
         "$log([NII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PBPT(n=1)


#BPT SII (Quelli che in altre parti di codice chiamo Type2)

irad = []
ishock = []
ihii2 = []
xrad = []
exrad = []
yrad = []
eyrad = []
xshock = []
exshock = []
yshock = []
eyshock = []
xhii2 = []
exhii2 = []
yhii2 = []
eyhii2 = []
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

PlotScat(xrad, yrad, ex=exrad, ey=eyrad, xlim=None, ylim=None, colore="red", simbolo="o",
         labels=["$log([SII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PlotScat(xshock, yshock, ex=exshock, ey=eyshock, xlim=None, ylim=None, colore="blue", simbolo="o", labels=[
         "$log([SII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PlotScat(xhii2, yhii2, ex=exhii2, ey=eyhii2, xlim=None, ylim=None, colore="green", simbolo="o", labels=[
         "$log([SII]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PBPT(n=2)


#BPT OI (Quelli che in altre parti di codice chiamo Type3)

irad3 = []
ishock3 = []
ihii3 = []
xrad3 = []
exrad3 = []
yrad3 = []
eyrad3 = []
xshock3 = []
exshock3 = []
yshock3 = []
eyshock3 = []
xhii3 = []
exhii3 = []
yhii3 = []
eyhii3 = []
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

PlotScat(xrad3, yrad3, ex=exrad3, ey=eyrad3, xlim=None, ylim=None, colore="red", simbolo="o",
         labels=["$log([OI]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"])
PlotScat(xshock3, yshock3, ex=exshock3, ey=eyshock3, xlim=None, ylim=None, colore="blue", simbolo="o",
         labels=["$log([OI]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PlotScat(xhii3, yhii3, ex=exhii3, ey=eyhii3, xlim=None, ylim=None, colore="green", simbolo="o", labels=[
         "$log([OI]/H \\alpha])$", "$log([OIII]/H \\beta])$"], Positives=["no", "no"], overplot=True)
PBPT(n=3)


# Save subsample properties


arrays = [RA, DEC, Z, eZ, SIG, eSIG, EBV, Zsun, SIGCLUSTER, NUMGAL, SIGMA_BAL, eSIGMA_BAL, SIGMA_FORB, eSIGMA_FORB,
          VOFF_BAL, eVOFF_BAL, VOFF_FORB, eVOFF_FORB, Mass, eMass1, eMass2, SFR, eSFR1, eSFR2, sSFR, esSFR1, esSFR2]
if Sampl == 'no':
    if SamePos == "yes":
        SaveType(
            iAGN, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal_RadioLoud_SamePos.txt", arrays)
        SaveType(
            icomp, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal_RadioLoud_SamePos.txt", arrays)
        SaveType(
            ihii1, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal_RadioLoud_SamePos.txt", arrays)
        SaveType(
            irad, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD_gal_RadioLoud_SamePos.txt", arrays)
        SaveType(
            ishock, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK_gal_RadioLoud_SamePos.txt", arrays)
        SaveType(
            ihii2, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2_gal_RadioLoud_SamePos.txt", arrays)
        SaveType(
            irad3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD3_gal_RadioLoud_SamePos.txt", arrays)
        SaveType(ishock3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK3_gal_RadioLoud_SamePos.txt", arrays)
        SaveType(
            ihii3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII3_gal_RadioLoud_SamePos.txt", arrays)
    else:
        SaveType(
            iAGN, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal_RadioLoud.txt", arrays)
        SaveType(
            icomp, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal_RadioLoud.txt", arrays)
        SaveType(
            ihii1, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal_RadioLoud.txt", arrays)
        SaveType(
            irad, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD_gal_RadioLoud.txt", arrays)
        SaveType(
            ishock, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK_gal_RadioLoud.txt", arrays)
        SaveType(
            ihii2, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2_gal_RadioLoud.txt", arrays)
        SaveType(
            irad3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD3_gal_RadioLoud.txt", arrays)
        SaveType(
            ishock3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK3_gal_RadioLoud.txt", arrays)
        SaveType(
            ihii3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII3_gal_RadioLoud.txt", arrays)
else:
    SaveType(
        iAGN, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN-RadioLoud.txt", arrays)
    SaveType(icomp, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp-RadioLoud.txt", arrays)
    SaveType(ihii1, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1-RadioLoud.txt", arrays)
    SaveType(
        irad, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD-RadioLoud.txt", arrays)
    SaveType(ishock, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK-RadioLoud.txt", arrays)
    SaveType(ihii2, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2-RadioLoud.txt", arrays)
    SaveType(irad3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD3-RadioLoud.txt", arrays)
    SaveType(ishock3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK3-RadioLoud.txt", arrays)
    SaveType(ihii3, "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII3-RadioLoud.txt", arrays)


# %% Determinazione delle Fraction RadioLoud ( SULL'INTERSEZIONE )

# Lettura dei dati
Sampl = 'no'
SamePos = "yes"     # Far Girare SOLAMENTE CON SAMEPOS = yes (Manca l'implementazione corretta del percorso !!)

# Leggo i file del Type1
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal_RadioLoud_SamePos.txt"
        
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal.txt"
        
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII-RadioLoud.txt"
    

i1, x1, ex1, y1, ey1 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

# Leggo i file del type2
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal_RadioLoud_SamePos.txt"
        
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal.txt"
        
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII-RadioLoud.txt"
    

i2, x2, ex2, y2, ey2 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

# Calcola A ∩ B
intersection_AB = np.intersect1d(i1, i2)

# Stampa il numero di elementi nell'intersezione
print("L'intersezione è composta da un numero di", len(intersection_AB), "Elementi")

# Usa np.where per ottenere gli indici comuni
elementi_comuni_i1 = np.where(np.isin(i1, intersection_AB))[0]
elementi_comuni_i2 = np.where(np.isin(i2, intersection_AB))[0]

# Ricavo i relativi elementi x, y facenti capo a questo indice di elementi in comune
x1_selected = x1[elementi_comuni_i1]
y1_selected = y1[elementi_comuni_i1]
x2_selected = x2[elementi_comuni_i2]
y2_selected = y2[elementi_comuni_i2]
ex1_selected = ex1[elementi_comuni_i1]
ex2_selected = ex2[elementi_comuni_i2]
ey1_selected = ey1[elementi_comuni_i1]
ey2_selected = ey2[elementi_comuni_i2]

# Contributo dal type1
fraction1 = np.zeros(5000)
fraction2 = np.zeros(5000)
fraction3 = np.zeros(5000)
NumOgg1  = np.zeros(5000)
NumOgg2 = np.zeros(5000)
NumOgg3 = np.zeros(5000)

for k in range(5000):
    xVAR = np.random.normal(x1_selected, ex1_selected)
    yVAR = np.random.normal(y1_selected, ey1_selected)
    condition_a = (yVAR >= (0.61 / (xVAR - 0.47)) + 1.19) | (xVAR >= 0.04)
    condition_b = (yVAR < (0.61 / (xVAR - 0.47)) + 1.19) & (yVAR >= (0.61 / (xVAR - 0.05)) + 1.3)
    condition_c = (yVAR < (0.61 / (xVAR - 0.05)) + 1.3) & (xVAR < 0.04)
    fraction1[k] = len(np.where(condition_a)[0]) / len(intersection_AB)
    fraction2[k] = len(np.where(condition_b)[0]) / len(intersection_AB)
    fraction3[k] = len(np.where(condition_c)[0]) / len(intersection_AB)
    NumOgg1[k] = len(np.where(condition_a)[0]) 
    NumOgg2[k] = len(np.where(condition_b)[0])
    NumOgg3[k] = len(np.where(condition_c)[0]) 
    
    

fractTOT_1 = np.mean(fraction1)  #f_1a
fractTOT_2 = np.mean(fraction2)  #f_1b
fractTOT_3 = np.mean(fraction3)  #f_1c
fractERR_1 = np.std(fraction1)
fractERR_2 = np.std(fraction2)
fractERR_3 = np.std(fraction3)
sum1 = fractTOT_1 + fractTOT_2 + fractTOT_3
esum1 = fractERR_1 + fractERR_2 + fractERR_3
N1A = int(np.mean(NumOgg1))
N1B = int(np.mean(NumOgg2))
N1C = int(np.mean(NumOgg3))
eN1A = int(np.std(NumOgg1))
eN1B = int(np.std(NumOgg2))
eN1C = int(np.std(NumOgg3))


# Contributo dal type2
fraction4 = np.zeros(5000)
fraction5 = np.zeros(5000)
fraction6 = np.zeros(5000)
NumOgg4  = np.zeros(5000)
NumOgg5 = np.zeros(5000)
NumOgg6 = np.zeros(5000)
for k in range(5000):
    
    #Creazione delle distribuzioni
    xVAR = np.random.normal(x2_selected, ex2_selected)
    yVAR = np.random.normal(y2_selected, ey2_selected)
    
    #Condizioni di collocamento dei punti
    condition_d = (yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29)
    condition_d = condition_d & (yVAR >= 1.89 * xVAR + 0.76)
    condition_e = (yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29)
    condition_e = condition_e & (yVAR < 1.89 * xVAR + 0.76)
    condition_f = ~((yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29))

    #Popolamento degli array di fraction
    fraction4[k] = len(np.where(condition_d)[0]) / len(x2_selected)
    fraction5[k] = len(np.where(condition_e)[0]) / len(x2_selected)
    fraction6[k] = len(np.where(condition_f)[0]) / len(x2_selected)
    
    #Popolamento degli array di conteggio
    NumOgg4[k] = len(np.where(condition_d)[0]) 
    NumOgg5[k] = len(np.where(condition_e)[0])
    NumOgg6[k] = len(np.where(condition_f)[0]) 

fractTOT_4 = np.mean(fraction4) #f_2a
fractERR_4 = np.std(fraction4)

fractTOT_5 = np.mean(fraction5) #f_2b
fractERR_5 = np.std(fraction5)

fractTOT_6 = np.mean(fraction6) #f_2c
fractERR_6 = np.std(fraction6)

sum2 = fractTOT_4 + fractTOT_5 + fractTOT_6
esum2 = fractERR_4 + fractERR_5 + fractERR_6

N2A = int(np.mean(NumOgg4))
N2B = int(np.mean(NumOgg5))
N2C = int(np.mean(NumOgg6))
eN2A = int(np.std(NumOgg4))
eN2B = int(np.std(NumOgg5))
eN2C = int(np.std(NumOgg6))


# Stampa i risultati
print("Ho ottenuto i seguenti risultati:\n")
print("Identificazione con BPT-NII\n")
print("AGN =", N1A, "+-", eN1A, "\n")
print("AGN =", fractTOT_1, "+-", fractERR_1, "\n")
print("COMPOSITE =", N1B, "+-", eN1B, "\n")
print("COMPOSITE =", fractTOT_2, "+-", fractERR_2, "\n")
print("SF Galaxies =", N1C, "+-", eN1C, "\n")
print("SF Galaxies =", fractTOT_3, "+-", fractERR_3, "\n")
print("Check Normalizzazione:", sum1,"+-" ,  esum1)



# Verifica l'unione delle condizioni
all_points_covered = np.logical_or.reduce([condition_a, condition_b, condition_c])

# Stampa il numero di punti coperti e il totale dei punti
print("Numero di punti coperti:", np.sum(all_points_covered))
print("Totale dei punti:", len(all_points_covered),"\n\n\n")



print("Identificazione con BPT-SII\n")
print("RADIATIVE =", N2A, "+-", eN2A, "\n")
print("RADIATIVE =", fractTOT_4, "+-", fractERR_4, "\n")
print("SHOCK =", N2B, "+-", eN2B, "\n")
print("SHOCK =", fractTOT_5, "+-", fractERR_5, "\n")
print("SF Galaxies =", N2C, "+-", eN2C, "\n")
print("SF Galaxies =", fractTOT_6, "+-", fractERR_6, "\n\n")
print("Check Normalizzazione:", sum2, "+-",esum2)

# Verifica l'unione delle condizioni
all_points_covered = np.logical_or.reduce([condition_d, condition_e, condition_f])

# Stampa il numero di punti coperti e il totale dei punti
print("Numero di punti coperti:", np.sum(all_points_covered))
print("Totale dei punti:", len(all_points_covered))

# %% SENZA INTERSEZIONE !!!

# Lettura dei dati
Sampl = 'no'
SamePos = "yes"     # Far Girare SOLAMENTE CON SAMEPOS = yes (Manca l'implementazione corretta del percorso !!)

# Leggo i file del Type1
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal_RadioLoud_SamePos.txt"
        
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII_gal.txt"
        
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-NII-RadioLoud.txt"
    

i1, x1, ex1, y1, ey1 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)

# Leggo i file del type2
if Sampl == 'no':
    if SamePos == "yes":
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal_RadioLoud_SamePos.txt"
        
    else:
        f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII_gal.txt"
        
else:
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/BPT-SII-RadioLoud.txt"
    

i2, x2, ex2, y2, ey2 = np.loadtxt(
    f, usecols=[0, 1, 2, 3, 4], unpack=True, dtype=float)


print("Sull'NII abbiamo", len(i1), "Elementi, mantre nel SII abbiamo", len(i2), "Elementi")

# Contributo dal type1
fraction1 = np.zeros(5000)
fraction2 = np.zeros(5000)
fraction3 = np.zeros(5000)
NumOgg1  = np.zeros(5000)
NumOgg2 = np.zeros(5000)
NumOgg3 = np.zeros(5000)

for k in range(5000):
    xVAR = np.random.normal(x1, ex1)
    yVAR = np.random.normal(y1, ey1)
    condition_a = (yVAR >= (0.61 / (xVAR - 0.47)) + 1.19) | (xVAR >= 0.04)
    condition_b = (yVAR < (0.61 / (xVAR - 0.47)) + 1.19) & (yVAR >= (0.61 / (xVAR - 0.05)) + 1.3)
    condition_c = (yVAR < (0.61 / (xVAR - 0.05)) + 1.3) & (xVAR < 0.04)
    fraction1[k] = len(np.where(condition_a)[0]) / len(i1)
    fraction2[k] = len(np.where(condition_b)[0]) / len(i1)
    fraction3[k] = len(np.where(condition_c)[0]) / len(i1)
    NumOgg1[k] = len(np.where(condition_a)[0]) 
    NumOgg2[k] = len(np.where(condition_b)[0])
    NumOgg3[k] = len(np.where(condition_c)[0]) 
    
    

fractTOT_1 = np.mean(fraction1)  #f_1a
fractTOT_2 = np.mean(fraction2)  #f_1b
fractTOT_3 = np.mean(fraction3)  #f_1c
fractERR_1 = np.std(fraction1)
fractERR_2 = np.std(fraction2)
fractERR_3 = np.std(fraction3)
sum1 = fractTOT_1 + fractTOT_2 + fractTOT_3
esum1 = fractERR_1 + fractERR_2 + fractERR_3
N1A = int(np.mean(NumOgg1))
N1B = int(np.mean(NumOgg2))
N1C = int(np.mean(NumOgg3))
eN1A = int(np.std(NumOgg1))
eN1B = int(np.std(NumOgg2))
eN1C = int(np.std(NumOgg3))

# Contributo dal type2
fraction4 = np.zeros(5000)
fraction5 = np.zeros(5000)
fraction6 = np.zeros(5000)
NumOgg4  = np.zeros(5000)
NumOgg5 = np.zeros(5000)
NumOgg6 = np.zeros(5000)
for k in range(5000):
    
    #Creazione delle distribuzioni
    xVAR = np.random.normal(x2, ex2)
    yVAR = np.random.normal(y2, ey2)
    
    #Condizioni di collocamento dei punti
    condition_d = (yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29)
    condition_d = condition_d & (yVAR >= 1.89 * xVAR + 0.76)
    condition_e = (yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29)
    condition_e = condition_e & (yVAR < 1.89 * xVAR + 0.76)
    condition_f = ~((yVAR >= (0.72 / (xVAR - 0.32)) + 1.30) | (xVAR > 0.29))

    #Popolamento degli array di fraction
    fraction4[k] = len(np.where(condition_d)[0]) / len(i2)
    fraction5[k] = len(np.where(condition_e)[0]) / len(i2)
    fraction6[k] = len(np.where(condition_f)[0]) / len(i2)
    
    #Popolamento degli array di conteggio
    NumOgg4[k] = len(np.where(condition_d)[0]) 
    NumOgg5[k] = len(np.where(condition_e)[0])
    NumOgg6[k] = len(np.where(condition_f)[0]) 

fractTOT_4 = np.mean(fraction4) #f_2a
fractERR_4 = np.std(fraction4)

fractTOT_5 = np.mean(fraction5) #f_2b
fractERR_5 = np.std(fraction5)

fractTOT_6 = np.mean(fraction6) #f_2c
fractERR_6 = np.std(fraction6)

sum2 = fractTOT_4 + fractTOT_5 + fractTOT_6
esum2 = fractERR_4 + fractERR_5 + fractERR_6

N2A = int(np.mean(NumOgg4))
N2B = int(np.mean(NumOgg5))
N2C = int(np.mean(NumOgg6))
eN2A = int(np.std(NumOgg4))
eN2B = int(np.std(NumOgg5))
eN2C = int(np.std(NumOgg6))


# Stampa i risultati
print("Ho ottenuto i seguenti risultati, Senza considerare l'insieme intersezione !!!':\n")
print("Identificazione con BPT-NII\n")
print("AGN =", N1A, "+-", eN1A, "\n")
print("AGN =", fractTOT_1, "+-", fractERR_1, "\n")
print("COMPOSITE =", N1B, "+-", eN1B, "\n")
print("COMPOSITE =", fractTOT_2, "+-", fractERR_2, "\n")
print("SF Galaxies =", N1C, "+-", eN1C, "\n")
print("SF Galaxies =", fractTOT_3, "+-", fractERR_3, "\n")
print("Check Normalizzazione:", sum1,"+-" ,  esum1)



# Verifica l'unione delle condizioni
all_points_covered = np.logical_or.reduce([condition_a, condition_b, condition_c])

# Stampa il numero di punti coperti e il totale dei punti
print("Numero di punti coperti:", np.sum(all_points_covered))
print("Totale dei punti:", len(all_points_covered),"\n\n\n")



print("Identificazione con BPT-SII\n")
print("RADIATIVE =", N2A, "+-", eN2A, "\n")
print("RADIATIVE =", fractTOT_4, "+-", fractERR_4, "\n")
print("SHOCK =", N2B, "+-", eN2B, "\n")
print("SHOCK =", fractTOT_5, "+-", fractERR_5, "\n")
print("SF Galaxies =", N2C, "+-", eN2C, "\n")
print("SF Galaxies =", fractTOT_6, "+-", fractERR_6, "\n\n")
print("Check Normalizzazione:", sum2, "+-",esum2)

# Verifica l'unione delle condizioni
all_points_covered = np.logical_or.reduce([condition_d, condition_e, condition_f])

# Stampa il numero di punti coperti e il totale dei punti
print("Numero di punti coperti:", np.sum(all_points_covered))
print("Totale dei punti:", len(all_points_covered))

# %% Dynamic galaxy versus mass per BPT subsamples

f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN_gal_SamePos.txt"
nSIG1, neSIG1, nSIGMA_BAL1, neSIGMA_BAL1, nSIGMA_FORB1, neSIGMA_FORB1, nM1, nSFR1, nsSFR1 = np.loadtxt(
    f1, usecols=[4, 5, 8, 9, 10, 11, 16, 19, 22], unpack=True, dtype=float)
f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_AGN.txt"
SIG1, eSIG1, SIGMA_BAL1, eSIGMA_BAL1, SIGMA_FORB1, eSIGMA_FORB1, M1, SFR1, sSFR1 = np.loadtxt(
    f1, usecols=[4, 5, 8, 9, 10, 11, 16, 19, 22], unpack=True, dtype=float)


f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp_gal_SamePos.txt"
nSIG2, neSIG2, nSIGMA_BAL2, neSIGMA_BAL2, nSIGMA_FORB2, neSIGMA_FORB2, nM2, nSFR2, nsSFR2 = np.loadtxt(
    f2, usecols=[4, 5, 8, 9, 10, 11, 16, 19, 22], unpack=True, dtype=float)
f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_Comp.txt"
SIG2, eSIG2, SIGMA_BAL2, eSIGMA_BAL2, SIGMA_FORB2, eSIGMA_FORB2, M2, SFR2, sSFR2 = np.loadtxt(
    f2, usecols=[4, 5, 8, 9, 10, 11, 16, 19, 22], unpack=True, dtype=float)


f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1_gal_SamePos.txt"
nSIG3, neSIG3, nSIGMA_BAL3, neSIGMA_BAL3, nSIGMA_FORB3, neSIGMA_FORB3, nM3, nSFR3, nsSFR3 = np.loadtxt(
    f3, usecols=[4, 5, 8, 9, 10, 11, 16, 19, 22], unpack=True, dtype=float)
f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII1.txt"
SIG3, eSIG3, SIGMA_BAL3, eSIGMA_BAL3, SIGMA_FORB3, eSIGMA_FORB3, M3, SFR3, sSFR3 = np.loadtxt(
    f3, usecols=[4, 5, 8, 9, 10, 11, 16, 19, 22], unpack=True, dtype=float)


def HistTot(dati_x, dati_y, colori, labels, assilabels, bins=[30, 30], limsx=None, limsy=None):
    if (limsx is None) == False:
        indicix = []
        for k in range(len(dati_x)):
            if np.isnan(limsx[k][0]) == False:
                i1 = np.where(dati_x[k][:] >= limsx[k][0])[0]
            else:
                i1 = []
            if np.isnan(limsx[k][1]) == False:
                i2 = np.where(dati_x[k][:] <= limsx[k][1])[0]
            else:
                i2 = []
            indicix.append(np.intersect1d(i1, i2))
        indicix = np.asarray(indicix)
    if (limsy is None) == False:
        indiciy = []
        for k in range(len(dati_y)):
            if np.isnan(limsy[k][0]) == False:
                i1 = np.where(dati_y[k][:] >= limsy[k][0])[0]
            else:
                i1 = []
            if np.isnan(limsy[k][1]) == False:
                i2 = np.where(dati_y[k][:] <= limsy[k][1])[0]
            else:
                i2 = []
            indiciy.append(np.intersect1d(i1, i2))
        indiciy = np.asarray(indiciy)
    indici = []
    if (limsx is None) == False and (limsy is None) == False:
        for k in range(len(dati_x)):
            indici.append(np.union1d(indicix[k], indiciy[k]))
    if (limsx is None) == False and (limsy is None) == True:
        for k in range(len(dati_x)):
            indici.append(indicix[k])
    if (limsx is None) == True and (limsy is None) == False:
        for k in range(len(dati_x)):
            indici.append(indiciy[k])
    if (limsx is None) == True and (limsy is None) == True:
        indici.append(np.arange(0, len(dati_x[k]), 1))
    indici = np.asarray(indici)
    datix = []
    datiy = []
    for k in range(len(dati_y)):
        datix.append(dati_x[k][indici[k]])
        datiy.append(dati_y[k][indici[k]])
    datix = np.asarray(datix)
    datiy = np.asarray(datiy)
    Histograms(datix, datiy, colori, labels, assilabels, bins=[40, 20])
    return 0


limsy = [[0.001, 499], [0.001, 499], [0.001, 499]]
dati_x = [sSFR1, sSFR2, sSFR3]
dati_y = [SIGMA_BAL1, SIGMA_BAL2, SIGMA_BAL3]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(sSFR)~[yr^{-1}]$", "$\\sigma_{Balmer}~[km/s]$"]
HistTot(dati_x, dati_y, colori, labels, assilabels,
        bins=[30, 30], limsx=None, limsy=limsy)


# limsx=[[-16,-6.5],[-16,-6.5],[-16,-6.5]]
limsy = [[0.001, 499], [0.001, 499], [0.001, 499]]
dati_x = [nsSFR1, nsSFR2, nsSFR3]
dati_y = [nSIGMA_BAL1, nSIGMA_BAL2, nSIGMA_BAL3]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(M)~[M_{\\odot}]$", "$\\sigma_{Balmer}~[km/s]$"]
HistTot(dati_x, dati_y, colori, labels, assilabels,
        bins=[100, 30], limsx=None, limsy=limsy)


limsy = [[0.001, 299], [0.001, 299], [0.001, 299]]
dati_x = [sSFR1, sSFR2, sSFR3]
dati_y = [VOFF_BAL1, VOFF_BAL2, VOFF_BAL3]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(sSFR)~[yr^{-1}]$", "$v_{off,~Balmer}~[km/s]$"]
HistTot(dati_x, dati_y, colori, labels, assilabels,
        bins=[30, 30], limsx=None, limsy=limsy)


limsy = [[0.001, 499], [0.001, 499], [0.001, 499]]
dati_x = [sSFR1, sSFR2, sSFR3]
dati_y = [SIGMA_BAL1, SIGMA_BAL2, SIGMA_BAL3]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(sSFR)~[yr^{-1}]$", "$\\sigma_{Balmer}~[km/s]$"]
HistTot(dati_x, dati_y, colori, labels, assilabels,
        bins=[30, 30], limsx=None, limsy=limsy)


limsy = [[0.001, 299], [0.001, 299], [0.001, 299]]
dati_x = [M1, M2, M3]
dati_y = [VOFF_BAL1, VOFF_BAL2, VOFF_BAL3]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(M)~[M_{\\odot}]$", "$v_{off,~Balmer}~[km/s]$"]
HistTot(dati_x, dati_y, colori, labels, assilabels,
        bins=[30, 30], limsx=None, limsy=limsy)


limsy = [[0.001, 499], [0.001, 499], [0.001, 499]]
dati_x = [M1, M2, M3]
dati_y = [SIGMA_BAL1, SIGMA_BAL2, SIGMA_BAL3]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(M)~[M_{\\odot}]$", "$\\sigma_{Balmer}~[km/s]$"]
HistTot(dati_x, dati_y, colori, labels, assilabels,
        bins=[30, 30], limsx=None, limsy=limsy)


# %%

"""
Come la presenza di un AGN modifica le proprietà delle BCG e delle galassie non classificate come BCG?
"""

"""
Metallicità in funzione della sigma per AGN Comp e HII
"""

i1, i2, i3 = np.where(Zsun1 < 1)[0], np.where(
    Zsun2 < 1)[0], np.where(Zsun3 < 1)[0]
dati_x = [Zsun1[i1], Zsun2[i2], Zsun3[i3]]
dati_y = [SIG1[i1], SIG2[i2], SIG3[i3]]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$Z/Z_{\odot}$", "$\\sigma~[km/s]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[20, 20])


"""
sSFR in funzione della sigma delle Balmer line per AGN Comp e HII
"""

i1, i2, i3 = np.where((SIGMA_BAL1 > 0.001) & (SIGMA_BAL1 < 499))[0], np.where(
    (SIGMA_BAL2 > 0.001) & (SIGMA_BAL2 < 499))[0], np.where((SIGMA_BAL3 > 0.001) & (SIGMA_BAL3 < 499))[0]
dati_x = [sSFR1[i1], sSFR2[i2], sSFR3[i3]]
dati_y = [SIGMA_BAL1[i1], SIGMA_BAL2[i2], SIGMA_BAL3[i3]]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(sSFR)~[M_{\\odot}~yr^{-1}]$", "$\\sigma_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[40, 20])


"""
MAIN SEQUENCE per AGN Comp e HII
"""

#i1,i2,i3=np.where((SIGMA_BAL1 > 0.001) & (SIGMA_BAL1 < 499))[0],np.where((SIGMA_BAL2 > 0.001) & (SIGMA_BAL2 < 499))[0],np.where((SIGMA_BAL3 > 0.001) & (SIGMA_BAL3 < 499))[0]
# dati_x=[Mass1[i1],Mass2[i2],Mass3[i3]]
# dati_y=[SFR1[i1],SFR2[i2],SFR3[i3]]
dati_x = [Mass1, Mass2, Mass3]
dati_y = [SFR1, SFR2, SFR3]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(M)~[M_{\\odot}]$", "$\log(SFR)~[M_{\\odot}~yr^{-1}]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[40, 40])


"""
Metallicità in funzione della Voffset per AGN Comp e HII
"""
i1, i2, i3 = np.where((VOFF_BAL1 > 0.001) & (VOFF_BAL1 < 299))[0], np.where(
    (VOFF_BAL2 > 0.001) & (VOFF_BAL2 < 299))[0], np.where((VOFF_BAL3 > 0.001) & (VOFF_BAL3 < 299))[0]
dati_x = [Zsun1[i1], Zsun2[i2], Zsun3[i3]]
dati_y = [np.log10(VOFF_BAL1[i1]), np.log10(
    VOFF_BAL2[i2]), np.log10(VOFF_BAL3[i3])]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$Z/Z_{\odot}$", "$v_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[20, 20])


"""
Metallicità in funzione della Voffset per AGN Comp e HII
"""
i1, i2, i3 = np.where((SIGMA_BAL1 > 0.001) & (SIGMA_BAL1 < 499) & (VOFF_BAL1 > 0.001) & (VOFF_BAL1 < 299))[0], np.where((SIGMA_BAL2 > 0.001) & (
    SIGMA_BAL2 < 499) & (VOFF_BAL2 > 0.001) & (VOFF_BAL2 < 299))[0], np.where((SIGMA_BAL3 > 0.001) & (SIGMA_BAL3 < 499) & (VOFF_BAL3 > 0.001) & (VOFF_BAL3 < 299))

dati_x = [Zsun1[i1], Zsun2[i2], Zsun3[i3]]
dati_y = [2*VOFF_BAL1[i1]+SIGMA_BAL1[i1], 2*VOFF_BAL2[i2] +
          SIGMA_BAL2[i2], 2*VOFF_BAL3[i3]+SIGMA_BAL3[i3]]

colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$Z/Z_{\odot}$", "$v_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[20, 20])


# %%
PlotScat(Mass1, SIGMA_FORB1, ey=eSIGMA_FORB1, xlim=None,
         ylim=[0.001, 499], colore="red", simbolo="o")
PlotScat(Mass2, SIGMA_FORB2, ey=eSIGMA_FORB2, xlim=None, ylim=[
         0.001, 499], colore="blue", simbolo="o", overplot=True)
PlotScat(Mass3, SIGMA_FORB3, ey=eSIGMA_FORB3, xlim=None, ylim=[0.001, 499], colore="green", simbolo="o", labels=[
         "$log(M)~[M_{\\odot}]$", "$\\sigma_{Forbidden}~[km/s]$"], overplot=True)


PlotScat(Mass1, SIGMA_BAL1, ey=eSIGMA_BAL1, xlim=None,
         ylim=[0.001, 499], colore="red", simbolo="o")
PlotScat(Mass2, SIGMA_BAL2, ey=eSIGMA_BAL2, xlim=None, ylim=[
         0.001, 499], colore="blue", simbolo="o", overplot=True)
PlotScat(Mass3, SIGMA_BAL3, ey=eSIGMA_BAL3, xlim=None, ylim=[0.001, 499], colore="green", simbolo="o", labels=[
         "$log(M)~[M_{\\odot}]$", "$\\sigma_{Balmer}~[km/s]$"], overplot=True)

PlotScat(Mass1, VOFF_BAL1, ey=eVOFF_BAL1, xlim=None,
         ylim=[0.001, 299], colore="red", simbolo="o")
PlotScat(Mass2, VOFF_BAL2, ey=eVOFF_BAL2, xlim=None, ylim=[
         0.001, 299], colore="blue", simbolo="o", overplot=True)
PlotScat(Mass3, VOFF_BAL3, ey=eVOFF_BAL3, xlim=None, ylim=[0.001, 299], colore="green", simbolo="o", labels=[
         "$log(M)~[M_{\\odot}]$", "$v_{Balmer}~[km/s]$"], overplot=True)

PlotScat(sSFR1, VOFF_BAL1, ey=eVOFF_BAL1, xlim=None,
         ylim=[0.001, 299], colore="red", simbolo="o")
PlotScat(sSFR2, VOFF_BAL2, ey=eVOFF_BAL2, xlim=None, ylim=[
         0.001, 299], colore="blue", simbolo="o", overplot=True)
PlotScat(sSFR3, VOFF_BAL3, ey=eVOFF_BAL3, xlim=None, ylim=[0.001, 299], colore="green", simbolo="o", labels=[
         "$log(sSFR)~[M_{\\odot}~yr^{-1}]$", "$v_{Balmer}~[km/s]$"], overplot=True)

PlotScat(SFR1, VOFF_BAL1, ey=eVOFF_BAL1, xlim=None,
         ylim=[0.001, 299], colore="red", simbolo="o")
PlotScat(SFR2, VOFF_BAL2, ey=eVOFF_BAL2, xlim=None, ylim=[
         0.001, 299], colore="blue", simbolo="o", overplot=True)
PlotScat(SFR3, VOFF_BAL3, ey=eVOFF_BAL3, xlim=None, ylim=[0.001, 299], colore="green", simbolo="o", labels=[
         "$log(SFR)~[M_{\\odot}~yr^{-1}]$", "$v_{Balmer}~[km/s]$"], overplot=True)


PlotScat(sSFR1, EBV1, ylim=[0, 2], colore="red", simbolo="o")
PlotScat(sSFR2, EBV2, ylim=[0, 2], colore="blue", simbolo="o", overplot=True)
PlotScat(sSFR3, EBV3, ylim=[0, 2], colore="green", simbolo="o", labels=[
         "$log(SFR)~[M_{\\odot}~yr^{-1}]$", "$E(B-V)$"], overplot=True)


i1, i2, i3 = np.where((VOFF_BAL1 > 0.001) & (VOFF_BAL1 < 299))[0], np.where(
    (VOFF_BAL2 > 0.001) & (VOFF_BAL2 < 299))[0], np.where((VOFF_BAL3 > 0.001) & (VOFF_BAL3 < 299))[0]
dati_x = [sSFR1[i1], sSFR2[i2], sSFR3[i3]]
dati_y = [VOFF_BAL1[i1], VOFF_BAL2[i2], VOFF_BAL3[i3]]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(sSFR)~[M_{\\odot}~yr^{-1}]$", "$v_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[100, 20])

i1, i2, i3 = np.where((SIGMA_BAL1 > 0.001) & (SIGMA_BAL1 < 499))[0], np.where(
    (SIGMA_BAL2 > 0.001) & (SIGMA_BAL2 < 499))[0], np.where((SIGMA_BAL3 > 0.001) & (SIGMA_BAL3 < 499))[0]
dati_x = [sSFR1[i1], sSFR2[i2], sSFR3[i3]]
dati_y = [SIGMA_BAL1[i1], SIGMA_BAL2[i2], SIGMA_BAL3[i3]]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(sSFR)~[M_{\\odot}~yr^{-1}]$", "$\\sigma_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[40, 20])

i1, i2, i3 = np.where((SIGMA_BAL1 > 0.001) & (SIGMA_BAL1 < 499))[0], np.where(
    (SIGMA_BAL2 > 0.001) & (SIGMA_BAL2 < 499))[0], np.where((SIGMA_BAL3 > 0.001) & (SIGMA_BAL3 < 499))[0]
dati_x = [Mass1[i1], Mass2[i2], Mass3[i3]]
dati_y = [SIGMA_BAL1[i1], SIGMA_BAL2[i2], SIGMA_BAL3[i3]]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(M)~[M_{\\odot}]$", "$\\sigma_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[40, 20])


i1, i2, i3 = np.where((VOFF_FORB1 > 0.001) & (VOFF_FORB1 < 299))[0], np.where(
    (VOFF_FORB2 > 0.001) & (VOFF_FORB2 < 299))[0], np.where((VOFF_FORB3 > 0.001) & (VOFF_FORB3 < 299))[0]
dati_x = [sSFR1[i1], sSFR2[i2], sSFR3[i3]]
dati_y = [VOFF_FORB1[i1], VOFF_FORB2[i2], VOFF_FORB3[i3]]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(sSFR)~[M_{\\odot}~yr^{-1}]$", "$v_{Forbidden}~[km/s]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[40, 20])


# %%
#Metallicity - Sigma
PlotScat(Zsun, SIG, ey=eSIG, xlim=None, ylim=None, colore="black",
         simbolo="o", labels=["$Z/Z_{\\odot}$", "$\\sigma~[km/s]$"])

#Metallicity - Sigma(FORBIDDEN)
PlotScat(SIGMA_BAL, SIGMA_FORB, ex=eSIGMA_BAL, ey=eSIGMA_FORB, xlim=[0.001, 499], ylim=[
         0.001, 499], colore="black", simbolo="o", labels=["$\\sigma_{BALMER}~[km/s]$", "$\\sigma_{FORBIDDEN}~[km/s]$"])
plt.plot([0, 500], [0, 500], color='red')


# %%  Radiative Shock & HII

Sampl = "y"
#RA,DEC,Z,eZ,SIG,eSIG,EBV,Zsun,SIGMA_BAL,eSIGMA_BAL,SIGMA_FORB,eSIGMA_FORB,VOFF_BAL,eVOFF_BAL,VOFF_FORB,eVOFF_FORB,OII_3726,eOII_3726,OII_3729,eOII_3729,NEIII_3869,eNEIII_3869,H_DELTA,eH_DELTA,H_GAMMA,eH_GAMMA,OIII_4363,eOIII_4363,OIII_4959,eOIII_4959,OIII_5007,eOIII_5007,HEI_5876,eHEI_5876,OI_6300,eOI_6300,H_BETA,eH_BETA,H_ALPHA,eH_ALPHA,NII_6584,eNII_6584,SII_6717,eSII_6717,SII_6731,eSII_6731,ARIII7135,eARIII7135,Mass,eMass1,eMass2,SFR,eSFR1,eSFR2,sSFR,esSFR1,esSFR2= calldata(Sampl)


if Sampl == 'no':
    f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD_gal.txt"
else:
    f1 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_RAD.txt"
RA1, DEC1, Z1, eZ1, SIG1, eSIG1, EBV1, Zsun1 = np.loadtxt(
    f1, usecols=[0, 1, 2, 3, 4, 5, 6, 7], unpack=True, dtype=float)
SIGMA_BAL1, eSIGMA_BAL1, SIGMA_FORB1, eSIGMA_FORB1 = np.loadtxt(
    f1, usecols=[8, 9, 10, 11], unpack=True, dtype=float)
VOFF_BAL1, eVOFF_BAL1, VOFF_FORB1, eVOFF_FORB1 = np.loadtxt(
    f1, usecols=[12, 13, 14, 15], unpack=True, dtype=float)
Mass1, eMass11, eMass12, SFR1, eSFR11, eSFR12, sSFR1, esSFR11, esSFR12 = np.loadtxt(
    f1, usecols=[16, 17, 18, 19, 20, 21, 22, 23, 24], unpack=True, dtype=float)

if Sampl == 'no':
    f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK_gal.txt"
else:
    f2 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_SHOCK.txt"
RA2, DEC2, Z2, eZ2, SIG2, eSIG2, EBV2, Zsun2 = np.loadtxt(
    f2, usecols=[0, 1, 2, 3, 4, 5, 6, 7], unpack=True, dtype=float)
SIGMA_BAL2, eSIGMA_BAL2, SIGMA_FORB2, eSIGMA_FORB2 = np.loadtxt(
    f2, usecols=[8, 9, 10, 11], unpack=True, dtype=float)
VOFF_BAL2, eVOFF_BAL2, VOFF_FORB2, eVOFF_FORB2 = np.loadtxt(
    f2, usecols=[12, 13, 14, 15], unpack=True, dtype=float)
Mass2, eMass21, eMass22, SFR2, eSFR21, eSFR22, sSFR2, esSFR21, esSFR22 = np.loadtxt(
    f2, usecols=[16, 17, 18, 19, 20, 21, 22, 23, 24], unpack=True, dtype=float)

if Sampl == 'no':
    f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2_gal.txt"
else:
    f3 = "/Users/andreamaccarinelli/Desktop/myOutputs3/Prop_HII2.txt"
RA3, DEC3, Z3, eZ3, SIG3, eSIG3, EBV3, Zsun3 = np.loadtxt(
    f3, usecols=[0, 1, 2, 3, 4, 5, 6, 7], unpack=True, dtype=float)
SIGMA_BAL3, eSIGMA_BAL3, SIGMA_FORB3, eSIGMA_FORB3 = np.loadtxt(
    f3, usecols=[8, 9, 10, 11], unpack=True, dtype=float)
VOFF_BAL3, eVOFF_BAL3, VOFF_FORB3, eVOFF_FORB3 = np.loadtxt(
    f3, usecols=[12, 13, 14, 15], unpack=True, dtype=float)
Mass3, eMass31, eMass32, SFR3, eSFR31, eSFR32, sSFR3, esSFR31, esSFR32 = np.loadtxt(
    f3, usecols=[16, 17, 18, 19, 20, 21, 22, 23, 24], unpack=True, dtype=float)


# %%

"""
sSFR in funzione della sigma delle Balmer line per AGN Comp e HII
"""

i1, i2, i3 = np.where((SIGMA_BAL1 > 0.001) & (SIGMA_BAL1 < 499))[0], np.where(
    (SIGMA_BAL2 > 0.001) & (SIGMA_BAL2 < 499))[0], np.where((SIGMA_BAL3 > 0.001) & (SIGMA_BAL3 < 499))[0]
dati_x = [sSFR1[i1], sSFR2[i2], sSFR3[i3]]
dati_y = [SIGMA_BAL1[i1], SIGMA_BAL2[i2], SIGMA_BAL3[i3]]
colori = ["red", "blue", "green"]
labels = ["AGN", "Jet", "HII"]
assilabels = ["$log(sSFR)~[M_{\\odot}~yr^{-1}]$", "$\\sigma_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[40, 20])


i1, i2, i3 = np.where((SIGMA_BAL1 > 0.001) & (SIGMA_BAL1 < 499))[0], np.where(
    (SIGMA_BAL2 > 0.001) & (SIGMA_BAL2 < 499))[0], np.where((SIGMA_BAL3 > 0.001) & (SIGMA_BAL3 < 499))[0]
dati_x = [Mass1[i1], Mass2[i2], Mass3[i3]]
dati_y = [SIGMA_BAL1[i1], SIGMA_BAL2[i2], SIGMA_BAL3[i3]]
colori = ["red", "blue", "green"]
labels = ["AGN", "Composite", "HII"]
assilabels = ["$log(M)~[M_{\\odot}]$", "$\\sigma_{Balmer}~[km/s]$"]
Histograms(dati_x, dati_y, colori, labels, assilabels, bins=[40, 20])
