# SatelitePlots

Plots de python de datos crudos satelitales (GOES-19 / GOES-East).

## `goes_ir_plot.py`

Procesa un archivo NetCDF (`.nc`) crudo de GOES-19, canal infrarrojo **Banda 13**
(10.3 µm, "Clean IR Longwave Window"), y genera un plot en **estilo MetPy clásico**:

- Mapa en proyección geoestacionaria real (tomada de los metadatos del propio archivo).
- Recuadro blanco fino alrededor del plot.
- Colormap/colorbar del producto satelital al costado.
- Título superior en fuente **Arial** con el formato `GOES-East CH<numero>`.
- Fecha/hora abajo en formato `YYYY-MM-DD HH`, convertida a **hora local de
  Argentina (ART, UTC-3)**.
- **No** incluye station plots / observaciones de superficie (todavía no disponibles).

Soporta tanto archivos ya calibrados (`ABI-L2-CMIP*`, variable `CMI` en Kelvin)
como archivos de radiancia cruda sin calibrar (`ABI-L1b-Rad`, variable `Rad`),
convirtiendo automáticamente esta última a temperatura de brillo con la
función de Planck inversa y los coeficientes que trae el propio archivo.

### Instalación del entorno (Miniconda 3 + VS Code)

1. Abrí una terminal (Anaconda Prompt / terminal integrada de VS Code) en la
   carpeta del repositorio clonado.
2. Creá el entorno a partir del archivo `environment.yml` incluido:

   ```bash
   conda env create -f environment.yml
   ```

3. Activá el entorno:

   ```bash
   conda activate sateliteplots
   ```

4. En VS Code, seleccioná este entorno como interprete de Python
   (`Ctrl+Shift+P` → *Python: Select Interpreter* → elegir `sateliteplots`).

Si preferís no usar el archivo `environment.yml`, podés crear el entorno a
mano:

```bash
conda create -n sateliteplots python=3.11 -y
conda activate sateliteplots
conda install -c conda-forge metpy xarray netcdf4 cartopy matplotlib numpy -y
```

### Uso

```bash
python goes_ir_plot.py --file OR_ABI-L2-CMIPF-M6C13_G19_sYYYYJJJHHMMSSZ.nc
```

Esto genera un archivo `.png` con el mismo nombre que el `.nc` de entrada.

Opciones principales:

| Flag | Descripción | Default |
|------|-------------|---------|
| `-f, --file` | Archivo `.nc` de entrada (requerido) | — |
| `-o, --output` | Ruta del `.png` de salida | `<nombre_input>.png` |
| `--colormap` | Colortable de MetPy (`ir_drgb`, `ir_rgbv`, `ir_bd`, `ir_tpc`, `ir_tv1`, `WVCIMSS`, etc.) | `ir_drgb` |
| `--vmin` / `--vmax` | Rango de temperatura de brillo en °C para la escala de colores | `-90` / `50` |
| `--extent LON_MIN LON_MAX LAT_MIN LAT_MAX` | Recorta el plot a una región lon/lat | Extensión completa del archivo |
| `--no-map-layers` | No dibuja costas/fronteras/provincias | (dibuja por defecto) |
| `--show` | Muestra el plot en una ventana interactiva | — |

Ejemplo recortando a Argentina y usando otra paleta:

```bash
python goes_ir_plot.py -f mi_archivo.nc --extent -75 -50 -55 -20 --colormap ir_rgbv
```

### Nota sobre el estilo visual

La paleta de colores exacta y el rango de temperaturas (`--colormap`, `--vmin`,
`--vmax`) se dejaron configurables porque no se pudo verificar pixel a pixel
contra la imagen de referencia original. Si el resultado no coincide 100% con
lo esperado, ajustá estos parámetros (o los valores por defecto directamente
en el script, sección `CONFIGURACION POR DEFECTO`) hasta lograr el resultado
deseado.

### Origen de los datos

Los archivos GOES-19 (`G19`) en formato `ABI-L1b-Rad` o `ABI-L2-CMIP` (Banda 13)
se pueden descargar, por ejemplo, desde el bucket público de AWS S3
`noaa-goes19` (sin necesidad de credenciales), o desde el
[NOAA Big Data Program](https://www.noaa.gov/information-technology/open-data-dissemination).
