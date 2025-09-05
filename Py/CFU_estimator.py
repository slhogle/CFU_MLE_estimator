import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import optimize

if (args_count := len(sys.argv)) > 2:
    print(f"One argument expected, got {args_count - 1}")
    raise SystemExit(2)
elif args_count < 2:
    print("You must specify the target CFU csv file")
    raise SystemExit(2)

target_f = Path(sys.argv[1])

if not target_f.is_file():
    print("The target file doesn't exist")
    raise SystemExit(1)

### Find the maximum likelihood estimator for CFUs (MPN method)
### samples - a numpy array of colony counts
### dil - a numpy array of dilutions
### V - the Volume
### N - Max number of Colonies
def findMLE(samples, dil, V, N):
    f = lambda x: np.sum(samples*dil/(N*(1-np.exp(-x*dil*V/N))))-np.sum(dil)
    if any(N<samples):
        return [np.inf, np.inf]
    sol = optimize.root_scalar(f, bracket=[0, 100000000], method='brentq')
    r_mle=sol.root
    p0=np.exp(-r_mle*dil*V/N)
    invVar=np.sum(dil*dil*V*V*samples*p0/(N*N*(1-p0)**2))
    estVar=1/invVar
    r_mle_std=np.sqrt(estVar)
    return [r_mle,r_mle_std]

### Find the Poisson estimator for CFUs
### samples - a numpy array of colony counts
### dil - a numpy array of dilutions
### V - the Volume
def findNaivePoisson(samples, dil, V):
    r_p=np.sum(samples)/(np.sum(dil)*V)
    r_p_std=np.sqrt(r_p*r_p/np.sum(samples))
    return [r_p, r_p_std]

### Find the Poisson estimator with a cutoff
### samples - a numpy array of colony counts
### dil - a numpy array of dilutions
### V - the Volume
### N - Cutoff above which there are crowding effects
def findPoissonCutoff(samples, dil, V, N):
    mask=samples<N
    r_p=np.sum(samples[mask])/(np.sum(dil[mask])*V)
    r_p_std=np.sqrt(r_p*r_p/np.sum(samples[mask]))
    return [r_p, r_p_std]

df = pd.read_csv(target_f, header = None)

print(findMLE(df[0], df[1], 1, 500))
print(findPoissonCutoff(df[0], df[1], 1, 100))