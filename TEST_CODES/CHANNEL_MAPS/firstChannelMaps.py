from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Abrir el archivo FITS
with fits.open("member.uid___A001_X3621_X1a6a.N159-13CII_sci.spw26.cube.I.pbcor.fits") as hdul:
    data = hdul[0].data
    header = hdul[0].header
    frames = range(1019, 1028, 1)  # Rango de frames para graficar

    #Contornos
    contornos = fits.open("N159_CII_map.fits")
    data_contorno = contornos[0].data
    header_contorno = contornos[0].header
    wcs_contorno = WCS(header_contorno)   

    # Crear un WCS 2D para las coordenadas espaciales (A partir del cubo de los datos originales CF+)
    wcs = WCS(header, naxis=2)

    # Crear una figura con subplots compartiendo ejes
    fig, axes = plt.subplots(3, 3, figsize=(12, 12), subplot_kw={'projection': wcs})
    # fig.subplots_adjust(hspace=0.15, wspace=0.4)  # Ajustar espacio entre subplots

    for idx, frame in enumerate(frames):
        row, col = divmod(idx, 3)  # Obtener fila y columna
        ax = axes[row, col]

        # Seleccionar el slice de datos
        data_slice = data[0, frame, :, :]
        ax.imshow(data_slice, origin="lower", cmap=plt.cm.viridis)


        ax.contour(data_contorno, transform=ax.get_transform(wcs_contorno),
           levels=[109.41,126.053,135.822,142.762,148.149,152.553,156.277,159.504,162.351,164.898,167.202], colors="white", alpha=0.5)

        # Configurar etiquetas solo en los bordes externos
        if col == 0:
            ax.coords[1].set_axislabel("DEC (J2000)", fontsize=10)
        else:
            ax.coords[1].set_ticklabel_visible(False)

        if row == 2:
            ax.coords[0].set_axislabel("RA (J2000)", fontsize=10)
        else:
            ax.coords[0].set_ticklabel_visible(False)

        # Títulos para cada frame
        # ax.legend(f"Frame {frame + 1}", fontsize=10)

        # Crear un "handle" falso para la leyenda (sin elementos gráficos)
        custom_legend = mpatches.Patch(color='white', label=f"Frame {frame + 1}")

        # Añadir la leyenda con el texto personalizado
        ax.legend(handles=[custom_legend], loc='upper right', fontsize=8)

        #Eliminar ticks arriba y a la derecha
        ax.tick_params(axis= "x", bottom="off")
        ax.tick_params(axis="y", left="off")

    # Ocultar subplots vacíos si hay menos de 9 frames
    for idx in range(len(frames), 9):
        row, col = divmod(idx, 3)
        axes[row, col].axis("off")
    plt.show()
    # Guardar y mostrar la figura
    fig.tight_layout()
    plt.savefig("plots.png", dpi=300)
    plt.show()
        

