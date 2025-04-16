from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.ndimage import gaussian_filter
from spectral_cube import SpectralCube
import astropy.units as u

# Rango de velocidades en km/s
vel_min = 232 * u.km / u.s
vel_max = 241 * u.km / u.s

# Abrir el archivo FITS principal (intensidad)
cube = SpectralCube.read("member.uid___A001_X3621_X1a6a.N159-13CII_sci.spw26.cube.I.pbcor.fits")
cube = cube.with_spectral_unit(u.km / u.s)  # Convertir el eje espectral a unidades de velocidad
subcube = cube.spectral_slab(vel_min, vel_max)  # Seleccionar el rango de velocidades
data = subcube.unmasked_data[:].value  # Extraer los datos como un array numpy
header = cube.header

# Abrir el archivo FITS de contornos
cube_contorno = SpectralCube.read("N159_CII_map.fits")
cube_contorno = cube_contorno.with_spectral_unit(u.km / u.s)  # Convertir el eje espectral a unidades de velocidad
subcube_contorno = cube_contorno.spectral_slab(vel_min, vel_max)  # Seleccionar el rango de velocidades
data_contorno = subcube_contorno.unmasked_data[:].value  # Extraer los datos como un array numpy
header_contorno = cube_contorno.header
wcs_contorno = WCS(header_contorno)

# Validar que ambos cubos tienen el mismo número de canales
if data.shape[0] != data_contorno.shape[0]:
    raise ValueError("El número de canales en los dos cubos no coincide.")

# Reducir el WCS de contornos a 2D (RA y DEC)
try:
    wcs_contorno_2d = wcs_contorno.slice([0])  # Ajusta según la estructura del archivo FITS
except Exception as e:
    raise ValueError(f"Error al reducir el WCS de contornos a 2D: {e}")

# Crear un WCS 2D para las coordenadas espaciales (A partir del cubo de los datos originales CF+)
wcs = WCS(header, naxis=2)

# Crear una figura con subplots compartiendo ejes
fig, axes = plt.subplots(3, 3, figsize=(12, 12), subplot_kw={'projection': wcs})

# Iterar sobre los frames (canales) dentro del rango de velocidades
frames = range(data.shape[0])  # Número de canales en el rango de velocidad
for idx, frame in enumerate(frames):
    row, col = divmod(idx, 3)  # Obtener fila y columna
    ax = axes[row, col]

    # Seleccionar el slice de datos del cubo principal
    data_slice = data[frame, :, :]
    im = ax.imshow(data_slice, origin="lower", cmap=plt.cm.pink, aspect="auto")

    # Seleccionar el slice correspondiente del archivo de contornos
    data_contorno_frame = data_contorno[frame, :, :]

    # Aplicar suavizado al slice de contornos
    data_contorno_smooth = gaussian_filter(data_contorno_frame, sigma=3)

    # Validar que los datos suavizados no estén vacíos
    if np.all(np.isnan(data_contorno_smooth)):
        raise ValueError(f"El slice de contornos para el frame {frame} contiene solo valores NaN.")

    # Calcular niveles de contorno con escala de potencia
    contorno_min = np.nanmin(data_contorno_smooth)
    contorno_max = np.nanmax(data_contorno_smooth)
    if contorno_min == contorno_max:
        raise ValueError(f"Los valores mínimo y máximo de los contornos son iguales para el frame {frame}.")
    levels = contorno_min + (contorno_max - contorno_min) * (np.linspace(0, 1, 11) ** 2)

    # Añadir contornos suavizados del frame actual
    ax.contour(data_contorno_smooth, transform=ax.get_transform(wcs_contorno_2d),
               levels=levels, colors="black", alpha=0.5)

    # Configurar etiquetas solo en los bordes externos
    if col == 0:
        ax.coords[1].set_axislabel("DEC (J2000)", fontsize=10)
    else:
        ax.coords[1].set_ticklabel_visible(False)

    if row == 2:
        ax.coords[0].set_axislabel("RA (J2000)", fontsize=10)
    else:
        ax.coords[0].set_ticklabel_visible(False)

    # Crear un "handle" falso para la leyenda (sin elementos gráficos)
    custom_legend = mpatches.Patch(color='white', label=f"Frame {frame + 1}")
    ax.legend(handles=[custom_legend], loc='upper right', fontsize=8)

# Ocultar subplots vacíos si hay menos de 9 frames
for idx in range(len(frames), 9):
    row, col = divmod(idx, 3)
    axes[row, col].axis("off")

# Ajustar el diseño y añadir una barra de color
fig.tight_layout()
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # Posición de la barra de color
fig.colorbar(im, cax=cbar_ax, label="Intensity")

# Guardar la figura
#plt.savefig("multiplot_with_velocity_range.png", dpi=300)
plt.show()