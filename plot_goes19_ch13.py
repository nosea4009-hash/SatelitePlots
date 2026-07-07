#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_goes19_ch13.py
====================

Procesa un archivo NetCDF (.nc) crudo del satelite GOES-19 (GOES-East),
producto ABI-L2-CMIP (Cloud and Moisture Imagery), Canal 13 (Infrarrojo
de onda larga, "Clean IR", ~10.3 um), y genera un plot con estilo
"MetPy clasico":

    - Mapa dentro de un recuadro blanco fino (spines finos, sin ejes de
      lat/lon visibles).
    - Colormap de infrarrojo tipo GOES (paleta "ir_rgbv" de MetPy: la
      clasica escala violeta -> azul -> verde -> amarillo -> rojo que se
      usa para resaltar tapas de nubes frias) mostrada como colorbar al
      costado del plot.
    - Titulo arriba en fuente Arial con el formato:
          GOES-East CH13
          YYYY-MM-DD HH (hora local Argentina, ART = UTC-3)

No se grafican estaciones automaticas / station plots (ASOS/METAR):
este script SOLO procesa y grafica la imagen cruda del satelite.

Uso
---
    python plot_goes19_ch13.py archivo.nc
    python plot_goes19_ch13.py archivo.nc --out mi_imagen.png
    python plot_goes19_ch13.py archivo.nc --extent -70 -50 -45 -20

Ver el README.md del repositorio para instrucciones de instalacion del
entorno de Conda.
"""

import argparse
import sys
from datetime import datetime, timedelta

import matplotlib

matplotlib.use("Agg")  # backend no interactivo, seguro para correr sin display

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib import font_manager
from metpy.plots import add_timestamp  # noqa: F401  (no se usa directamente, pero deja el import documentado)
from metpy.plots.ctables import registry

# ---------------------------------------------------------------------------
# Configuracion general
# ---------------------------------------------------------------------------

# Desplazamiento de Argentina respecto a UTC (ART = UTC-3, sin horario de
# verano desde 2009).
ART_OFFSET = timedelta(hours=-3)

# Rango de temperatura de brillo (en grados Celsius) sobre el que se
# construye la paleta de colores infrarroja. -80 C ~ +40 C es el rango
# estandar usado en la mayoria de los productos de infrarrojo GOES
# (resalta bien tanto la superficie/mar como los topes de nubes frias).
IR_TEMP_MIN_C = -80.0
IR_TEMP_MAX_C = 40.0

# Nombre de la colortable de MetPy para infrarrojo GOES Canal 13.
# "ir_rgbv" es la paleta clasica violeta-azul-verde-amarillo-rojo usada
# para resaltar el tope de nubes frias en el canal IR.
IR_COLORTABLE = "ir_rgbv"


def _configure_fonts():
    """Intenta usar Arial: si no esta disponible en el sistema, MetPy/
    Matplotlib caeran automaticamente a una fuente sans-serif similar
    (ej. DejaVu Sans / Liberation Sans), sin romper la ejecucion.
    """
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [
        "Arial",
        "Liberation Sans",
        "DejaVu Sans",
        "sans-serif",
    ]

    # Si el usuario tiene Arial.ttf instalada pero fontconfig no la detecto,
    # se puede registrar manualmente aqui (opcional).
    for font_path in font_manager.findSystemFonts():
        try:
            name = font_manager.FontProperties(fname=font_path).get_name()
        except Exception:
            continue
        if name.lower() == "arial":
            font_manager.fontManager.addfont(font_path)
            break


def _find_cmi_variable(ds):
    """Localiza la variable de temperatura de brillo (CMI) dentro del
    dataset. Los archivos ABI-L2-CMIP usan 'CMI'; los archivos
    Multichannel (MCMIP) usan 'CMI_C13'. Se soportan ambos.
    """
    for name in ("CMI", "CMI_C13"):
        if name in ds.variables:
            return name
    raise KeyError(
        "No se encontro la variable de datos del Canal 13 (se esperaba "
        "'CMI' o 'CMI_C13') en el archivo NetCDF. Verifica que el archivo "
        "sea un producto ABI-L2-CMIP o ABI-L2-MCMIP de GOES."
    )


def _get_band_number(ds, cmi_name):
    """Obtiene el numero de banda/canal del archivo. Si no se puede leer,
    asume 13 (este script esta pensado para el canal IR 13)."""
    try:
        if "band_id" in ds.variables:
            band = int(np.asarray(ds["band_id"].values).flatten()[0])
            return band
    except Exception:
        pass
    return 13


def _get_scan_start_datetime(ds):
    """Extrae el instante de inicio del escaneo (UTC) desde los atributos
    estandar de los archivos GOES-R (time_coverage_start), con un
    fallback a la variable de tiempo 't' (punto medio del escaneo)."""
    if "time_coverage_start" in ds.attrs:
        raw = ds.attrs["time_coverage_start"]
        # Formato tipico: '2025-08-22T13:54:00.0Z'
        raw = raw.replace("Z", "")
        for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S"):
            try:
                return datetime.strptime(raw, fmt)
            except ValueError:
                continue

    if "t" in ds.variables:
        t_val = ds["t"].values
        return pd_to_datetime_fallback(t_val)

    raise ValueError(
        "No se pudo determinar la fecha/hora de escaneo del archivo. "
        "Se esperaba el atributo 'time_coverage_start' o la variable 't'."
    )


def pd_to_datetime_fallback(t_val):
    """Convierte el valor de la variable 't' (datetime64 o numero) a un
    objeto datetime estandar sin depender obligatoriamente de pandas."""
    try:
        import pandas as pd

        return pd.to_datetime(t_val).to_pydatetime()
    except Exception:
        # Ultimo recurso: asumir que es datetime64[ns]
        return np.datetime64(t_val, "s").astype(datetime)


def _utc_to_art(dt_utc):
    """Convierte un datetime UTC (naive) a hora local Argentina (ART, UTC-3)."""
    return dt_utc + ART_OFFSET


def load_ch13(nc_path):
    """Abre el archivo NetCDF crudo de GOES-19 y devuelve:
        - data: xarray.DataArray de temperatura de brillo (K), con
          metadatos CF parseados por MetPy (incluye la proyeccion).
        - band: numero de canal (int).
        - scan_start_utc: datetime del inicio del escaneo (UTC).
    """
    ds = xr.open_dataset(nc_path)

    cmi_name = _find_cmi_variable(ds)
    band = _get_band_number(ds, cmi_name)
    scan_start_utc = _get_scan_start_datetime(ds)

    # metpy.parse_cf adjunta la informacion de proyeccion (goes_imager_projection)
    # a la variable de datos, para poder graficarla luego con Cartopy.
    data = ds.metpy.parse_cf(cmi_name)

    return data, band, scan_start_utc


def make_plot(data, band, scan_start_utc, extent=None, dpi=150,
              figsize=(9, 9)):
    """Genera la figura con estilo MetPy clasico: mapa satelital dentro de
    un recuadro blanco fino, colormap IR con colorbar al costado, y
    titulos superiores en Arial (GOES-East CHxx / fecha en hora ART)."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    _configure_fonts()

    # Proyeccion geoestacionaria nativa del dato (via el accessor de MetPy).
    proj = data.metpy.cartopy_crs
    x = data["x"]
    y = data["y"]

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    # Datos crudos vienen en Kelvin (temperatura de brillo); se convierten
    # a Celsius, que es la convencion habitual para la escala de color IR.
    data_celsius = data.values - 273.15

    # --- Colormap infrarrojo GOES Canal 13 (paleta clasica ir_rgbv) -----
    ir_norm, ir_cmap = registry.get_with_range(
        IR_COLORTABLE, IR_TEMP_MIN_C, IR_TEMP_MAX_C
    )

    im = ax.imshow(
        data_celsius,
        extent=(x.min(), x.max(), y.min(), y.max()),
        origin="upper",
        cmap=ir_cmap,
        norm=ir_norm,
        transform=proj,
        interpolation="nearest",
    )

    # Limites geograficos opcionales (lon_min, lon_max, lat_min, lat_max)
    if extent is not None:
        ax.set_extent(extent, crs=ccrs.PlateCarree())

    # --- Costas, paises y provincias/estados (estilo MetPy clasico) ----
    ax.add_feature(cfeature.COASTLINE, edgecolor="white", linewidth=0.75)
    ax.add_feature(cfeature.BORDERS, edgecolor="white", linewidth=0.75)
    ax.add_feature(
        cfeature.STATES.with_scale("50m"),
        edgecolor="white",
        linewidth=0.4,
    )

    # --- Recuadro blanco fino tipo MetPy clasico ------------------------
    # Los "spines" del axes actuan como el marco fino blanco/negro que
    # encierra el plot (comportamiento por defecto de MapPanel de MetPy).
    for spine in ax.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(0.75)

    ax.set_xticks([])
    ax.set_yticks([])

    # --- Colorbar (paleta del producto satelital) al costado -----------
    cbar = fig.colorbar(
        im,
        ax=ax,
        orientation="vertical",
        pad=0.02,
        fraction=0.046,
        extend="both",
    )
    cbar.set_label(u"Temperatura de brillo (\u00b0C)", fontsize=10, family=plt.rcParams["font.sans-serif"][0])
    cbar.ax.tick_params(labelsize=8)

    # --- Titulos superiores (Arial): "GOES-East CHxx" + fecha en ART ----
    scan_start_art = _utc_to_art(scan_start_utc)
    title_line1 = f"GOES-East CH{band:02d}"
    title_line2 = scan_start_art.strftime("%Y-%m-%d %H") + " ART"

    ax.set_title(
        title_line1,
        loc="center",
        fontsize=16,
        fontweight="bold",
        family=plt.rcParams["font.sans-serif"][0],
        pad=28,
    )
    # Segunda linea (fecha), un poco mas chica, arriba del plot tambien.
    ax.text(
        0.5,
        1.02,
        title_line2,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=12,
        family=plt.rcParams["font.sans-serif"][0],
    )

    fig.tight_layout()
    return fig, ax


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Procesa datos crudos GOES-19 (.nc, ABI-L2-CMIP Canal 13) y "
            "genera un plot con estilo MetPy clasico (sin station plots)."
        )
    )
    parser.add_argument(
        "archivo_nc",
        help="Ruta al archivo NetCDF (.nc) crudo de GOES-19, Canal 13.",
    )
    parser.add_argument(
        "--out",
        "-o",
        default=None,
        help=(
            "Ruta de salida para la imagen PNG generada. Por defecto usa "
            "el mismo nombre del .nc con extension .png."
        ),
    )
    parser.add_argument(
        "--extent",
        nargs=4,
        type=float,
        metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"),
        default=None,
        help=(
            "Recorte geografico opcional en grados (PlateCarree): "
            "lon_min lon_max lat_min lat_max. Si se omite, se grafica "
            "el dominio completo del archivo."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Resolucion (DPI) de la imagen de salida. Default: 150.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Ademas de guardar el PNG, intenta mostrar la figura en pantalla.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    try:
        data, band, scan_start_utc = load_ch13(args.archivo_nc)
    except Exception as exc:
        print(f"ERROR al leer '{args.archivo_nc}': {exc}", file=sys.stderr)
        return 1

    fig, ax = make_plot(
        data,
        band,
        scan_start_utc,
        extent=args.extent,
        dpi=args.dpi,
    )

    out_path = args.out
    if out_path is None:
        out_path = args.archivo_nc.rsplit(".", 1)[0] + ".png"

    fig.savefig(out_path, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    print(f"Imagen guardada en: {out_path}")

    if args.show:
        plt.show()

    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
