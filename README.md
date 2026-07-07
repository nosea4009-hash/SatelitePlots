# SatelitePlots

Plots de Python de datos crudos satelitales (GOES-19 / GOES-East).

## Script principal: `plot_goes19_ch13.py`

Procesa un archivo NetCDF (`.nc`) crudo del satelite GOES-19 (soporta
productos **ABI-L1b Radiancias**, **ABI-L2-CMIP** y **ABI-L2-MCMIP**),
Canal 13 (infrarrojo de onda larga, "Clean IR", ~10.3 um), y genera una
imagen con estilo **MetPy clasico**:

- Mapa satelital dentro de un recuadro blanco fino, con costas, paises y
  provincias/estados dibujados en blanco.
- Colormap infrarroja personalizada para el Canal 13: escala de grises
  para temperaturas calidas (superficie, cielo despejado) y una escala
  de colores de realce (cian -> azul -> verde -> amarillo -> naranja ->
  rojo -> negro -> blanco) para las temperaturas frias de los topes de
  nube, con su colorbar al costado del plot, en grados Celsius.
- **Region por defecto: Provincia de Cordoba (Argentina) y alrededores**
  (se puede cambiar con `--extent` o desactivar el recorte con
  `--full-disk` para ver el dominio completo del archivo).
- Titulo arriba del plot en fuente Arial, en dos lineas:
  - `GOES-East CH13`
  - `YYYY-MM-DD HH` (hora **local Argentina / ART**, UTC-3)

Este script **no** grafica estaciones automaticas / station plots
(ASOS, METAR, etc.) — solo procesa y grafica la imagen cruda del satelite.

---

## 1) Instalar el entorno con Miniconda 3

Si todavia no tenes Miniconda instalado, descargalo desde:
https://docs.conda.io/en/latest/miniconda.html

Una vez instalado Miniconda (o Anaconda), abri la terminal integrada de
VSCode (o el "Anaconda Prompt" en Windows) y parate en la carpeta del
repositorio clonado:

```bash
cd ruta/a/SatelitePlots
```

### Opcion A (recomendada): crear el entorno desde `environment.yml`

Este repositorio incluye un archivo `environment.yml` con todas las
dependencias necesarias. Para crear el entorno automaticamente:

```bash
conda env create -f environment.yml
```

Esto crea un entorno llamado **`sateliteplots`** con Python 3.11,
MetPy, xarray, netCDF4, Cartopy, Matplotlib, NumPy y Pandas.

### Opcion B: crear el entorno manualmente

Si preferis instalar los paquetes a mano:

```bash
conda create --name sateliteplots -c conda-forge python=3.11 metpy xarray netcdf4 cartopy matplotlib numpy pandas
```

---

## 2) Activar el entorno

```bash
conda activate sateliteplots
```

(En Linux/Mac tambien funciona `source activate sateliteplots` en
versiones antiguas de conda).

Para verificar que el entorno quedo activo, deberias ver
`(sateliteplots)` al principio de la linea de tu terminal.

### Seleccionar el entorno en VSCode

1. Abri la paleta de comandos (`Ctrl+Shift+P` / `Cmd+Shift+P`).
2. Buscar **"Python: Select Interpreter"**.
3. Elegir el interprete correspondiente al entorno `sateliteplots`
   (deberia aparecer como `Python 3.11.x ('sateliteplots')`).

---

## 3) Ejecutar el script

Con el entorno activado, y con un archivo `.nc` de GOES-19 Canal 13
descargado (por ejemplo desde el bucket publico de AWS
`noaa-goes19` o desde tu estacion receptora GRB/GNC-A):

```bash
python plot_goes19_ch13.py ruta/al/archivo.nc
```

Esto genera una imagen PNG con el mismo nombre del `.nc` (por ejemplo
`archivo.png`) en la misma carpeta.

### Opciones disponibles

```bash
python plot_goes19_ch13.py ruta/al/archivo.nc --out salida.png
python plot_goes19_ch13.py ruta/al/archivo.nc --extent -70 -50 -45 -20
python plot_goes19_ch13.py ruta/al/archivo.nc --full-disk
python plot_goes19_ch13.py ruta/al/archivo.nc --dpi 200
python plot_goes19_ch13.py ruta/al/archivo.nc --show
```

- `--out / -o`: ruta de salida del PNG.
- `--extent LON_MIN LON_MAX LAT_MIN LAT_MAX`: recorta el plot a una
  region geografica especifica (en grados, PlateCarree). Si se omite,
  se usa la region por defecto (**Cordoba y alrededores**,
  aproximadamente 67°O-60°O / 36.5°S-28°S).
- `--full-disk`: grafica el dominio completo del archivo (Full Disk /
  CONUS / Mesoscale, segun el producto), sin recortar a Cordoba. Se
  ignora si tambien se pasa `--extent`.
- `--dpi`: resolucion de salida (por defecto 150).
- `--show`: ademas de guardar el PNG, intenta abrir una ventana con
  la figura.

### Formatos de archivo soportados

El script detecta automaticamente el tipo de producto GOES-R y la
variable de datos del Canal 13:

- Archivos **ABI-L1b-Rad** (Radiancias crudas), variable `Rad`: se
  convierte automaticamente a temperatura de brillo usando la formula
  de calibracion Planck (con los coeficientes que trae el propio
  archivo).
- Archivos **ABI-L2-CMIP** (Cloud and Moisture Imagery de un solo canal),
  variable `CMI`.
- Archivos **ABI-L2-MCMIP** (Multichannel), variable `CMI_C13`.

---

## Notas sobre la fuente Arial

El script intenta usar la fuente **Arial** para los titulos. Si Arial
no esta instalada en tu sistema, Matplotlib usara automaticamente una
fuente sans-serif similar (Liberation Sans / DejaVu Sans) sin que el
script falle. Para asegurar Arial exacta en Windows, generalmente ya
viene preinstalada con el sistema operativo.
