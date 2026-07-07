#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
goes_ir_plot.py
================

Procesa un archivo crudo de satelite GOES-19 (GOES-East) en formato NetCDF (.nc)
correspondiente al canal infrarrojo Banda 13 (10.3 um, "Clean IR Longwave Window")
y genera un plot en el estilo "MetPy clasico":

    - Proyeccion geoestacionaria real (tomada de los metadatos del propio archivo).
    - El mapa queda encerrado en un recuadro (spine) fino, sobre fondo blanco.
    - Colormap / colorbar del producto satelital al costado del plot.
    - Titulo superior en fuente Arial con el formato:  "GOES-East CH<numero>"
    - Texto inferior (marca de tiempo "estilo MetPy", via metpy.plots.add_timestamp)
      con la fecha/hora LOCAL de Argentina (ART, UTC-3) en formato "YYYY-MM-DD HH".

No se dibujan station plots / observaciones de superficie (todavia no disponibles).

Formatos de archivo soportados
-------------------------------
1) ABI-L2-CMIP (Cloud and Moisture Imagery, producto ya calibrado a temperatura
   de brillo). Variable de datos: ``CMI`` (unidades Kelvin). Este es el formato
   recomendado/mas comun para descargar desde AWS (bucket ``noaa-goes19``,
   producto ``ABI-L2-CMIPF`` / ``ABI-L2-CMIPC`` / ``ABI-L2-CMIPM``, canal 13).

2) ABI-L1b-Rad (Radiancias crudas, sin calibrar). Variable de datos: ``Rad``.
   El script convierte automaticamente la radiancia a temperatura de brillo
   usando la funcion de Planck inversa con los coeficientes que trae el propio
   archivo (``planck_fk1``, ``planck_fk2``, ``planck_bc1``, ``planck_bc2``).

Uso basico
----------
    python goes_ir_plot.py --file OR_ABI-L2-CMIPF-M6C13_G19_sYYYYJJJHHMMSSZ.nc

    python goes_ir_plot.py -f archivo.nc -o salida.png --vmin -90 --vmax 50

Ver ``python goes_ir_plot.py --help`` para todas las opciones.

Nota importante sobre el estilo visual
---------------------------------------
No tuve acceso directo a la imagen de referencia que adjuntaste en el chat (no
llego al entorno de este agente), asi que la paleta de colores exacta y el
rango de temperaturas fueron elegidos con valores tipicos usados en looping de
satelite IR (paleta ``ir_drgb`` de MetPy, rango -90 a 50 C). Todo esto es
facilmente configurable desde la seccion "CONFIGURACION POR DEFECTO" mas abajo
o por linea de comandos (--colormap, --vmin, --vmax) para que puedas ajustarlo
hasta que coincida exactamente con lo que buscas.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

try:
    from zoneinfo import ZoneInfo
except ImportError:  # Python < 3.9 (no deberia pasar en Miniconda moderno)
    from backports.zoneinfo import ZoneInfo  # type: ignore

import cartopy.feature as cfeature
from metpy.plots import add_timestamp, colortables

# ---------------------------------------------------------------------------
# CONFIGURACION POR DEFECTO (ajustar aca si no se usan los flags de linea de
# comandos). Estos valores controlan el "look" final del plot.
# ---------------------------------------------------------------------------
DEFAULT_COLORMAP = "ir_drgb"      # Paleta de MetPy para infrarrojo clasico.
DEFAULT_VMIN_C = -90.0            # Limite frio de la escala, en grados Celsius.
DEFAULT_VMAX_C = 50.0             # Limite calido de la escala, en grados Celsius.
DEFAULT_FIGSIZE = (10, 10)        # Tamano de la figura en pulgadas.
DEFAULT_DPI = 150                 # Resolucion de salida.
DRAW_MAP_LAYERS = True            # Dibujar costas / fronteras / provincias.
FONT_FAMILY = "Arial"             # Fuente solicitada para el titulo.
ARGENTINA_TZ = ZoneInfo("America/Argentina/Buenos_Aires")  # ART, UTC-3 fijo.

# Asegura que matplotlib intente usar Arial y, si no esta instalada en el
# sistema, caiga de forma prolija a una fuente sans-serif equivalente en vez
# de tirar un error.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [FONT_FAMILY, "Liberation Sans", "DejaVu Sans"]


# ---------------------------------------------------------------------------
# Utilidades de lectura / calibracion
# ---------------------------------------------------------------------------
def radiance_to_brightness_temperature(
    rad: xr.DataArray,
    fk1: float,
    fk2: float,
    bc1: float,
    bc2: float,
) -> xr.DataArray:
    """Convierte radiancia ABI cruda (Rad) a temperatura de brillo (Kelvin).

    Formula estandar de NOAA/NESDIS (funcion de Planck invertida) usada para
    los productos ABI-L1b-Rad. Los coeficientes ``fk1``, ``fk2``, ``bc1`` y
    ``bc2`` vienen incluidos como variables escalares en el propio archivo
    NetCDF.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        bt = (fk2 / np.log((fk1 / rad) + 1.0) - bc1) / bc2
    bt.attrs["units"] = "K"
    bt.attrs["long_name"] = "Brightness Temperature"
    return bt


def load_band13(nc_path: Path) -> tuple[xr.DataArray, int, datetime]:
    """Abre el archivo .nc y devuelve (DataArray de temperatura de brillo en K
    con metadatos CF/MetPy ya parseados, numero de banda, datetime UTC del
    inicio del escaneo).
    """
    ds = xr.open_dataset(nc_path)

    if "CMI" in ds.data_vars:
        # Producto Nivel 2 (ABI-L2-CMIP*) ya calibrado a temperatura de brillo.
        var_name = "CMI"
    elif "Rad" in ds.data_vars:
        # Producto Nivel 1b (ABI-L1b-Rad), radiancia cruda: hay que calibrar.
        var_name = "Rad"
    else:
        raise ValueError(
            "El archivo no contiene ni la variable 'CMI' (ABI-L2-CMIP) ni "
            "'Rad' (ABI-L1b-Rad). Verifica que sea un archivo GOES ABI valido."
        )

    # metpy.parse_cf() interpreta la proyeccion geoestacionaria (goes_imager_projection)
    # y deja lista la variable con .metpy.cartopy_crs, .x, .y en las unidades correctas.
    dat = ds.metpy.parse_cf(var_name)

    if var_name == "Rad":
        fk1 = float(ds["planck_fk1"].values)
        fk2 = float(ds["planck_fk2"].values)
        bc1 = float(ds["planck_bc1"].values)
        bc2 = float(ds["planck_bc2"].values)
        # Nota: la coordenada auxiliar "metpy_crs" (creada por parse_cf) se
        # conserva automaticamente en operaciones aritmeticas de xarray, por
        # lo que `dat` sigue teniendo la proyeccion correcta despues de esto.
        dat = radiance_to_brightness_temperature(dat, fk1, fk2, bc1, bc2)

    # Numero de banda (canal). Para archivos de un solo canal viene como escalar.
    if "band_id" in ds.variables:
        band = int(np.atleast_1d(ds["band_id"].values)[0])
    else:
        warnings.warn(
            "No se encontro la variable 'band_id' en el archivo; se asume "
            "banda 13 por defecto."
        )
        band = 13

    if band != 13:
        warnings.warn(
            f"El archivo corresponde a la banda {band}, no a la banda 13. "
            "El script fue pensado para IR banda 13; se continua igual."
        )

    # Hora de inicio del escaneo (UTC), viene en el atributo global del archivo.
    time_str = ds.attrs.get("time_coverage_start")
    if time_str is None:
        raise ValueError(
            "El archivo no tiene el atributo global 'time_coverage_start'."
        )
    scan_start_utc = _parse_iso_time(time_str)

    return dat, band, scan_start_utc


def _parse_iso_time(time_str: str) -> datetime:
    """Parsea el timestamp ISO-8601 (con 'Z') que usa GOES en sus metadatos.

    Los archivos GOES suelen traer fracciones de segundo con una cantidad de
    digitos variable (p. ej. ``.3Z`` o ``.300000Z``). ``datetime.fromisoformat``
    en Python < 3.11 solo acepta 0, 3 o 6 digitos de microsegundos, asi que se
    normaliza manualmente para evitar errores con archivos reales.
    """
    cleaned = time_str.strip().replace("Z", "")
    if "." in cleaned:
        base, frac = cleaned.split(".", 1)
        frac = (frac + "000000")[:6]  # pad/truncar a 6 digitos (microsegundos)
        cleaned = f"{base}.{frac}"
    dt = datetime.fromisoformat(cleaned)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_goes_band13(
    dat: xr.DataArray,
    band: int,
    scan_start_utc: datetime,
    output_path: Path | None,
    colormap: str = DEFAULT_COLORMAP,
    vmin_c: float = DEFAULT_VMIN_C,
    vmax_c: float = DEFAULT_VMAX_C,
    figsize: tuple[float, float] = DEFAULT_FIGSIZE,
    dpi: int = DEFAULT_DPI,
    draw_map_layers: bool = DRAW_MAP_LAYERS,
    extent: tuple[float, float, float, float] | None = None,
    show: bool = False,
) -> plt.Figure:
    """Genera el plot en estilo MetPy clasico y opcionalmente lo guarda/muestra."""

    crs = dat.metpy.cartopy_crs
    x = dat.x
    y = dat.y

    # Temperatura de brillo en Celsius para que el rango de la paleta sea mas
    # intuitivo (los "ir_*" de MetPy son escalas lineales genericas: se les
    # puede pasar cualquier unidad mientras se sea consistente).
    data_c = dat - 273.15

    norm, cmap = colortables.get_with_range(colormap, vmin_c, vmax_c)

    # --- Figura base: fondo blanco, tal como en los ejemplos "clasicos" de MetPy ---
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor="white")
    ax = fig.add_subplot(1, 1, 1, projection=crs)
    ax.set_facecolor("white")

    im = ax.imshow(
        data_c,
        extent=(x.min(), x.max(), y.min(), y.max()),
        origin="upper",
        cmap=cmap,
        norm=norm,
        transform=crs,
        interpolation="nearest",
    )

    if extent is not None:
        import cartopy.crs as ccrs

        ax.set_extent(extent, crs=ccrs.PlateCarree())

    if draw_map_layers:
        ax.coastlines(resolution="50m", color="black", linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor="black")
        ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.25, edgecolor="gray")

    # --- Recuadro fino alrededor del plot (spine geoestacionario de cartopy) ---
    for spine in ax.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(0.8)

    # --- Colorbar del producto satelital, al costado ---
    cbar = fig.colorbar(im, ax=ax, orientation="vertical", pad=0.02, fraction=0.045)
    cbar.set_label("Temperatura de brillo (\u00b0C)", fontsize=11, fontfamily=FONT_FAMILY)
    cbar.ax.tick_params(labelsize=9)

    # --- Titulo superior: "GOES-East CH<numero>" ---
    ax.set_title(
        f"GOES-East CH{band}",
        loc="center",
        fontsize=16,
        fontweight="bold",
        fontfamily=FONT_FAMILY,
        pad=10,
    )

    # --- Fecha/hora local Argentina (ART) abajo, con la utilidad "clasica" de MetPy ---
    scan_start_art = scan_start_utc.astimezone(ARGENTINA_TZ).replace(tzinfo=None)
    add_timestamp(
        ax,
        time=scan_start_art,
        pretext="",
        time_format="%Y-%m-%d %H",
        x=0.5,
        y=-0.04,
        ha="center",
        fontsize=12,
        fontfamily=FONT_FAMILY,
    )

    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=dpi, facecolor="white", bbox_inches="tight")
        print(f"Plot guardado en: {output_path}")

    if show:
        plt.show()

    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Procesa un archivo NetCDF crudo de GOES-19 (Banda 13, IR) y genera "
            "un plot en estilo MetPy clasico."
        )
    )
    parser.add_argument(
        "-f", "--file", required=True, type=Path, help="Ruta al archivo .nc de entrada."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Ruta del archivo de salida (png/jpg). Si no se especifica, se usa "
        "el mismo nombre que el .nc de entrada con extension .png.",
    )
    parser.add_argument(
        "--colormap",
        default=DEFAULT_COLORMAP,
        help=f"Nombre de la colortable de MetPy a usar (default: {DEFAULT_COLORMAP}). "
        "Otras opciones utiles: ir_rgbv, ir_bd, ir_tpc, ir_tv1, WVCIMSS.",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=DEFAULT_VMIN_C,
        help=f"Limite frio de la escala en grados Celsius (default: {DEFAULT_VMIN_C}).",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=DEFAULT_VMAX_C,
        help=f"Limite calido de la escala en grados Celsius (default: {DEFAULT_VMAX_C}).",
    )
    parser.add_argument(
        "--dpi", type=int, default=DEFAULT_DPI, help=f"DPI de salida (default: {DEFAULT_DPI})."
    )
    parser.add_argument(
        "--no-map-layers",
        action="store_true",
        help="No dibujar costas/fronteras/provincias sobre el plot.",
    )
    parser.add_argument(
        "--extent",
        nargs=4,
        type=float,
        metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"),
        default=None,
        help="Recorta el plot a un area lon/lat especifica, ej: --extent -75 -50 -45 -20",
    )
    parser.add_argument(
        "--show", action="store_true", help="Mostrar el plot en una ventana interactiva."
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    nc_path: Path = args.file
    if not nc_path.exists():
        print(f"ERROR: no se encontro el archivo '{nc_path}'.", file=sys.stderr)
        return 1

    output_path = args.output
    if output_path is None:
        output_path = nc_path.with_suffix(".png")

    dat, band, scan_start_utc = load_band13(nc_path)

    extent = tuple(args.extent) if args.extent else None

    plot_goes_band13(
        dat=dat,
        band=band,
        scan_start_utc=scan_start_utc,
        output_path=output_path,
        colormap=args.colormap,
        vmin_c=args.vmin,
        vmax_c=args.vmax,
        dpi=args.dpi,
        draw_map_layers=not args.no_map_layers,
        extent=extent,
        show=args.show,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
