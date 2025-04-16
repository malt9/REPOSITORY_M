from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Abrir el archivo FITS
with fits.open('N159_CII_map.fit')