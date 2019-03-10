import pandas as pd
import numpy as np

hllc = 'resultat_HLLC.txt'
hll  = 'resultat_HLL.txt'

# lire resultat hllc
df_hllc = pd.read_fwf(hllc,names=["X", "density", "velocity","energy","pression"])
for ii in df_hllc.columns:
    df_hllc[ii] = pd.to_numeric(df_hllc[ii])

# lire resultat hll
df_hll = pd.read_fwf(hll,names=["X", "density", "velocity","energy","pression"])
for ii in df_hll.columns:
    df_hll[ii] = pd.to_numeric(df_hll[ii])

for ii in df_hll.drop('X',axis=1).columns:
    print('erreur du', str(ii), '= ', np.linalg.norm(df_hll[ii]-df_hllc[ii],ord=np.inf))
