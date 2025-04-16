from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from spectral_cube import SpectralCube


# Open the files FITS
Cube1 = SpectralCube.read("N159_CII_map.fits")
Cube2 = SpectralCube.read("N159_CI809_map.fits")

#print(Cube1)
#print(Cube2)

#Define the number of N (column density)
N = 5 #for now is aleatory

#Obtain the Moment 0 
###Moment 0 correspond to intensity 

M0_Cube1 = Cube1.moment(order=0).value
M0_Cube2 = Cube2.moment(order=0).value

#print(M0_Cube1)
#print(M0_Cube2)

#Normalized the intensities
Cube1_Norm = M0_Cube1/N
Cube2_Norm = M0_Cube2/N

#print(Cube1_Norm)
#print(Cube2_Norm)

#Generate the values for each axis

x = Cube1_Norm.flatten()
y = Cube2_Norm.flatten()

#Filter the nan values
valid = np.isfinite(x) & np.isfinite(y)

x = x[valid]
y = y[valid]

#Create the plots
plt.figure(figsize=(8, 6))
plt.scatter(x, y, s=1, alpha=0.5, color='blue')
plt.xlabel("Integral Normalizada (Cube1)")
plt.ylabel("Integral Normalizada (Cube2)")
plt.title("Scatterplot de Integrales Normalizadas")
plt.grid(True)
plt.show()