#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 10:08:36 2024

@author: andreamaccarinelli
"""
# %% Creazione immagine insieme delle fractions !
#Script di prova per creare l'immagine della regione in cui calcoliamo le fraction AGN
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def PlotScat(x, y, ex=None, ey=None, xlim=None, ylim=None, colore="black", simbolo="o", labels=["X", "Y"], Positives=["yes", "yes"], overplot=False, alpha=1.0, label = "dato", s = 1):
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
    plt.scatter(x, y, color=colore, linestyle='None', marker=simbolo, alpha=alpha, label = label, s = s)
    if (ex is None) == False and (ey is None) == False:
        plt.errorbar(x, y, xerr=ex, yerr=ey, color=colore, linestyle='None', alpha=alpha)
    if (ex is None) == False and (ey is None) == True:
        plt.errorbar(x, y, xerr=ex, color=colore, linestyle='None', alpha=alpha)
    if (ex is None) == True and (ey is None) == False:
        plt.errorbar(x, y, yerr=ey, color=colore, linestyle='None', alpha=alpha)

    plt.xlabel(labels[0], fontsize=16)
    plt.ylabel(labels[1], fontsize=16)
    plt.tick_params(axis='both', labelsize=16)
    return len(x)

f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoNoBCG.txt"
RAdr7, DECdr7 = np.loadtxt(f, usecols=[0, 1], unpack=True, dtype=float)

f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoBCG.txt"
RAc4, DECc4 = np.loadtxt(f, usecols=[0, 1], unpack=True, dtype=float)

raR, decR = np.loadtxt(
    "/Users/andreamaccarinelli/Desktop/SDSS/RADIOtab.txt", usecols=[3, 4], unpack=True, dtype=float)


# ScatterPlot di tutto SDSS privato delle BCG

PlotScat(RAdr7, DECdr7, xlim = [0, 360], ylim = [-30, +90], colore = 'grey', labels = ['RA', 'DEC'], alpha = 0.01, label = 'SDSS/DR7', s = 0.08, overplot=True)

#ScatterPlot delle BCG
PlotScat(RAc4, DECc4, xlim = [0, 360], ylim = [-30, +90], colore = 'red', labels = ['RA', 'DEC'], alpha = 1, overplot=True, label='BCG', s = 0.5)


# Creazione del plot di dispersione
plt.scatter(raR, decR, color='orange', marker='o', label='RadioEmitters', s = 0.1)





# Supponiamo che tu abbia già un grafico creato e che le regioni siano definite.
# Assicurati di avere le variabili raR e decR definite in precedenza.

# Lista delle condizioni delle regioni

region_conditions = [
    ((raR > -2) & (raR < 53) & (decR > -12) & (decR < 4)),
    ((raR > 310) & (raR < 360) & (decR > -12) & (decR < 4)),
    ((raR > 120) & (raR < 250) & (decR > -5) & (decR < 6)),
    ((raR > 138) & (raR < 267) & (decR > 48) & (decR < 65)),
    ((raR > 130) & (raR < 138) & (decR > 48) & (decR < 60)),
    ((raR > 111) & (raR < 150) & (decR > 40) & (decR < 48)),
    ((raR > 112) & (raR < 136) & (decR > 27) & (decR < 40)),
    ((raR > 222) & (raR < 260) & (decR > 40) & (decR < 48)),
    ((raR > 252) & (raR < 267) & (decR > 25.5) & (decR < 40))
]



# Imposta il colore del bordo delle regioni a nero
edge_color = 'black'
# Imposta lo stile di linea per il bordo delle regioni (tratteggiato)
edge_linestyle = 'dashed'

# Itera attraverso le condizioni delle regioni e colora solo i bordi
for region_condition in region_conditions:
    ra_region = raR[region_condition]
    dec_region = decR[region_condition]

    # Trova i valori limite della regione
    ra_min, ra_max = ra_region.min(), ra_region.max()
    dec_min, dec_max = dec_region.min(), dec_region.max()

    # Plotta il bordo della regione con stile di linea tratteggiato
    plt.plot([ra_min, ra_max, ra_max, ra_min, ra_min], [dec_min, dec_min, dec_max, dec_max, dec_min],
             color=edge_color, linestyle=edge_linestyle)

# Impostazioni del plot
plt.xlabel('RA [deg]')
plt.ylabel('DEC [deg]')
plt.xlim((-2, 363))
plt.ylim((-12, +72))
#plt.title('Bordi delle Regioni')
# Legenda personalizzata
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=5, label='RadioEmitters'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=5, label='BCG'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='grey', markersize=2, label='SDSS/DR7'),
    
]

# Aggiungi la legenda personalizzata al plot
legend = plt.legend(handles=legend_elements, loc='upper left', markerscale=3)
plt.gca().add_artist(legend)


# Mostra il grafico
plt.show()

# %% Realizzazione delle fractions di oggetti RadioLoud 


#Step Zero: Lettura dei dati "Non includo nella trattazione il caso ( SamePos = "no" )

Sampl = "no"
#Sampl = "no"
SamePos= "yes"

if (Sampl == "no"):
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoNoBCG_SamePos.txt"
    ra, dec, status = np.loadtxt(f, usecols=[0,1,9], unpack=True, dtype=float)
    
else :
    f = "/Users/andreamaccarinelli/Desktop/myOutputs3/InfoBCG.txt"
    ra, dec, status = np.loadtxt(f, usecols=[0,1,11], unpack=True, dtype=float)
    

"""
Coordinate dei risultati presentati nell'abstract
region_conditions = [
    ((ra > -2) & (ra < 53) & (dec > -12) & (dec < 4)),
    ((ra > 310) & (ra < 360) & (dec > -12) & (dec < 4)),
    ((ra > 120) & (ra < 250) & (dec > -5) & (dec < 6)),
    ((ra > 130) & (ra < 267) & (dec > 48) & (dec < 65)),
    ((ra > 111) & (ra < 150) & (dec > 40) & (dec < 48)),
    ((ra > 112) & (ra < 136) & (dec > 27) & (dec < 40)),
    ((ra > 222) & (ra < 260) & (dec > 40) & (dec < 48)),
    ((ra > 252) & (ra < 267) & (dec > 25.5) & (dec < 40))
]
"""
region_conditions = [
    ((ra > -2)  & (ra < 53) & (dec > -12) & (dec < 4)),
    ((ra > 310) & (ra < 360) & (dec > -12) & (dec < 4)),
    ((ra > 120) & (ra < 250) & (dec > -5) & (dec < 6)),
    ((ra > 138) & (ra < 267) & (dec > 48) & (dec < 65)),
    ((ra > 130) & (ra < 138) & (dec > 48) & (dec < 60)),
    ((ra > 111) & (ra < 150) & (dec > 40) & (dec < 48)),
    ((ra > 112) & (ra < 136) & (dec > 27) & (dec < 40)),
    ((ra > 222) & (ra < 260) & (dec > 40) & (dec < 48)),
    ((ra > 252) & (ra < 267) & (dec > 25.5) & (dec < 40))
]


# Combina le condizioni utilizzando l'operatore logico OR
combined_condition = np.any(region_conditions, axis=0)

# Ricerca degli indici in cui compaiono solo RadioLoud
indexRLFIXX = np.where((status != 0) & combined_condition)[0]

# Ricerca di tutti gli oggetti nel ritaglio finale
indexFIXX = np.where(combined_condition)[0]

#calcolo del rapporto
frac = len(indexRLFIXX) / len(indexFIXX)
percent = frac * 100

if (Sampl == "no"):
    testo = " galassie non BCG"
else :
    testo = " BCG"

print("Nei ritagli finali",testo," ci sono ", len(indexRLFIXX), "Radioloud")
print("Su un totale di ", len(indexFIXX), testo)
print("Risultante in una fraction di RadioLoud corrispondente a ", frac, " corrispondente al ", percent, "%")


