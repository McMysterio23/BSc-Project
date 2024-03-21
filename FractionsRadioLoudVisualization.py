#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 10:45:35 2024

@author: andreamaccarinelli
"""

import matplotlib.pyplot as plt
import numpy as np

BCG_RL = 12
nonBCG_RL = 0.6

simbolo_BCG = '*'
simbolo_nonBCG = 'o'



# Nuove percentuali per la specie "Radio Loud" per BCG e nonBCG
percentuale_RadioLoud_BCG = 12
percentuale_RadioLoud_nonBCG = 0.6

# Nuovi simboli per le due specie
simbolo_RadioLoud = 'x'

# Nuovi errori per le due specie
errori_RadioLoud_BCG = 0  # Esempio, inserire l'errore reale
errori_RadioLoud_nonBCG = 0  # Esempio, inserire l'errore reale

# Creazione del plot
plt.figure(figsize=(8, 6))

# Lista di coordinate x
x = np.array([0])  

# Disegna i simboli e gli errori per la percentuale di "Radio Loud" BCG
plt.errorbar(x[0], percentuale_RadioLoud_BCG, yerr=errori_RadioLoud_BCG, fmt=simbolo_BCG, label='Radio Loud BCG', capsize=5, markersize = 13)

# Disegna i simboli e gli errori per la percentuale di "Radio Loud" nonBCG
plt.errorbar(x[0], percentuale_RadioLoud_nonBCG, yerr=errori_RadioLoud_nonBCG, fmt=simbolo_nonBCG, label='Radio Loud nonBCG', capsize=5, markersize = 13)

# Impostazioni dell'asse y
plt.ylabel('Fractions [%]')
plt.yticks(fontsize = 14)

# Impostazioni dell'asse x
plt.xticks(x, ['Radio Loud'], fontsize = 14)
#plt.xlabel('Radio Loud')

# Aggiungi una legenda
#plt.legend(loc='upper left')

# Rimuovi le linee di griglia
plt.grid(False)

# Mostra il plot
plt.show()
