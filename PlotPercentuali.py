#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 21:01:00 2024

@author: andreamaccarinelli
"""

import matplotlib.pyplot as plt
import numpy as np

# Dati delle percentuali e dei nomi delle popolazioni
percentuali_BCG = [58, 34.7, 7.3, 6.6, 67.6, 25.9]  # Percentuali per il campione BCG
percentuali_nonBCG = [12.04, 17.49, 70.47, 4.22, 9.03, 86.75]  # Percentuali per il campione nonBCG
percentualiVITALEetAl = []
errori_percentuali_BCG = [1.8, 2.2, 0.2, 1.7, 2.5, 2.1]  # Errori sulle percentuali per il campione BCG
errori_percentuali_nonBCG = [0.03, 0.05, 0.03, 0.03, 0.04, 0.04]  # Errori sulle percentuali per il campione nonBCG
errorivitale = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]

nomi_popolazioni = ['AGN', 'Composite', 'HII1', 'Seyfert', 'LINER', 'HII2']

# Simboli per i due campioni
simbolo_BCG = '*'
simbolo_nonBCG = 'o'
simbolo_VITALEetAl = 'D'

# Creazione del plot
plt.figure(figsize=(8, 6))

# Lista di coordinate x
x = np.arange(len(nomi_popolazioni))

# Disegna i simboli e gli errori per ogni percentuale del campione BCG
plt.errorbar(x, percentuali_BCG, yerr=errori_percentuali_BCG, fmt=simbolo_BCG, label='BCG', capsize = 5)

# Disegna i simboli e gli errori per ogni percentuale del campione nonBCG
plt.errorbar(x, percentuali_nonBCG, yerr=errori_percentuali_nonBCG, fmt=simbolo_nonBCG, label='nonBCG', capsize = 5)

# Impostazioni dell'asse y
plt.ylabel('Percentuale')

# Disegna solo le linee orizzontali per le altezze specificate sull'asse y
#plt.hlines(y=[25, 50, 75, 100], xmin=0, xmax=len(nomi_popolazioni)-1, color='black')


# Impostazioni dell'asse x
plt.xticks(x, nomi_popolazioni)
plt.xlabel('Popolazione')

# Aggiungi una legenda
plt.legend(loc='upper left')

# Rimuovi le linee di griglia
plt.grid(False)


# Mostra il plot
plt.show()