from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Abrir archivo FITS
hdu = fits.open("member.uid___A001_X3621_X1a6a.N159-13CII_sci.spw26.cube.I.pbcor.fits")
spectrum = hdu[0].data['CTYPE3']  # Ejemplo: columna con intensidades
freq = hdu[0].data['RESTFRQ']  # Columna de frecuencias
hdu.close()

# Graficar espectro
plt.plot(freq, spectrum)
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Flujo (Jy)")
plt.show()