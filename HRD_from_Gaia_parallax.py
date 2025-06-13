from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import copy
import math
from astropy.modeling import models
from specutils.fitting import fit_lines
from specutils import SpectralRegion
from specutils.fitting import estimate_line_parameters
from specutils.manipulation import extract_region
from specutils.spectra import Spectrum1D
from astropy import units as u
from astropy.modeling.fitting import LinearLSQFitter, SimplexLSQFitter, SLSQPLSQFitter, LevMarLSQFitter
from specutils import Spectrum1D, SpectralRegion, fitting, analysis, manipulation
from scipy import integrate, optimize, interpolate
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
from astropy.modeling import models, polynomial
import pandas as pd
import csv
import json

%matplotlib inline
%pylab inline

# Open and check the .json file to see what type of columns we have
with open('1697456835684O-result.json') as file:
  contents = json.loads(file.read())
  
metadata = contents["metadata"]
data = contents["data"]
print(len(data))
for i in range(len(metadata)):
    print("{}. {}: {}, unit: {}".format(i, metadata[i]["name"], metadata[i]["description"], metadata[i]["unit"]))

# Extract columns into variables
def extractColumn(data: list, columnIndex: int):
    return np.array(list(map(lambda x: x[columnIndex], data)))

data = list(filter(lambda x: x[1] and x[2] and x[3] and x[4], data))

# HRD from Gaia parallax values
T = extractColumn(data, 3) # Effective temperature
m = extractColumn(data, 2) # Apparent magnitude
d = extractColumn(data, 4) # Distance in parsec
p = extractColumn(data, 1) # Parallax

# Calculate distance and luminosity 
d1 = (1/p)*1000                # Distance from parallax
M1 = m - 5 * np.log10(d1) + 5  # Distance modulus
L0 = 3.0128*10**28             # Sun's luminosity in Watt
M0 = 4.83                      # Absolute magnitude of the Sun

L_L0 = []

L = L0*(10**(M0-M1/2.512))     # Luminosity of the stars
L_L0 = np.array(np.log10(L/L0))

# Data visualization
plt.figure(figsize = (14, 14))
plt.style.use('dark_background')

plt.ylabel('Absolute magnitude',fontsize=20)
plt.xticks(fontsize=13,color='black')
plt.yticks(fontsize=13)

plt.scatter(T,M1,s=1,c=T,cmap='RdYlBu')

ax = plt.gca()
ax.invert_xaxis()
ax.invert_yaxis()

ax2 = ax.twinx()
ax2.scatter(T,L_L0,s=1,c=T,cmap='RdYlBu')
ax2.set_ylabel('Luminosity [$\log{(L/L{\odot})}$]',fontsize=20)
ax2.set_yticklabels([-1.5,-1,0,1,2,3,4,5],fontsize=13)

ax3 = ax.twiny()
ax3.set_xticks((1,3,5,6,7,8,9))
ax3.set_xticklabels(['','','','','','',''],fontsize=16)
ax3.set_xlabel('Spectral types',fontsize=20)

plt.text(0.5,-2.8,'O',fontsize=20)
plt.text(1.95,-2.8,'B',fontsize=20)
plt.text(3.95,-2.8,'A',fontsize=20)
plt.text(5.5,-2.8,'F',fontsize=20)
plt.text(6.5,-2.8,'G',fontsize=20)
plt.text(7.5,-2.8,'K',fontsize=20)
plt.text(8.5,-2.8,'M',fontsize=20)

plt.text(7.5,1.7,'Instability strip',fontsize=20)
plt.text(5.6,8.8,'Subdwarves',fontsize=20)

clb = plt.colorbar(aspect=30, orientation='horizontal',pad=0.0)
clb.ax.invert_xaxis()
clb.set_label('$T_{eff}~~ [K]$',fontsize=20)
clb.ax.tick_params(labelsize=13)

show()