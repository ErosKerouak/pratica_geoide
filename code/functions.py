import numpy as np
import xarray as xr
import os
import rasterio
from rasterio.transform import Affine
from scipy.interpolate import RegularGridInterpolator


def dataarray_to_geotiff(da, out_path, nodata=np.nan, compress="LZW"):
    """
    Exporta um xarray.DataArray 2D (lat/lon) para GeoTIFF em EPSG:4326.

    Parâmetros
    ----------
    da : xr.DataArray
        Deve ter dims ('y', 'x') em graus.
    out_path : str
        Caminho do arquivo GeoTIFF de saída.
    nodata : float
        Valor para marcar NoData (default = np.nan).
    compress : str
        Algoritmo de compressão GTiff (ex: "LZW", "DEFLATE").

    """
    from rasterio.transform import from_origin
    y = da.y.values
    x = da.x.values
    dy = float(np.abs(np.diff(y)).mean())
    dx = float(np.abs(np.diff(x)).mean())

    # Rasterio espera origem no canto superior esquerdo
    transform = from_origin(x.min() - dx/2, y.max() + dy/2, dx, dy)

    arr = np.flipud(da.values.astype("float32"))  # flip Y para "top-down"

    with rasterio.open(
        out_path,
        "w",
        driver="GTiff",
        height=arr.shape[0],
        width=arr.shape[1],
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
        nodata=nodata,
        compress=compress,
        tiled=True,
        predictor=3,
    ) as dst:
        dst.write(arr, 1)

    print(f"[OK] Exportado GeoTIFF: {out_path}")


def geotiff_to_dataarray(geotiff_path: str, band: int = 1) -> xr.DataArray:
    """
    Lê um GeoTIFF (single- or multi-band) como um xarray.DataArray.

    Parâmetros
    ----------
    geotiff_path : str
        Caminho do arquivo GeoTIFF.
    band : int, opcional (default=1)
        Índice da banda a ler (1-based, como no rasterio).

    Retorna
    -------
    xr.DataArray
        DataArray com dims ('y','x') e coords 1D 'y' (lat) e 'x' (lon), quando não há rotação.
        Se houver rotação no transform, retorna coords 2D ('y','x') com arrays 'x' e 'y' bidimensionais.
    """
    if not os.path.exists(geotiff_path):
        raise FileNotFoundError(f"File not found: {geotiff_path}")

    with rasterio.open(geotiff_path) as src:
        arr = src.read(band, masked=True)          # masked array (respeita nodata)
        transform: Affine = src.transform
        crs = src.crs
        nodata = src.nodata
        height, width = src.height, src.width

    # Converte nodata -> NaN (mantém float)
    data = np.asarray(arr.filled(np.nan), dtype=float)

    # Componentes do affine: [a, b, c; d, e, f]
    a, b, c, d, e, f = transform.a, transform.b, transform.c, transform.d, transform.e, transform.f

    # Caso comum (sem rotação): b == 0 e d == 0  -> coords 1D
    if np.isclose(b, 0.0) and np.isclose(d, 0.0):
        # centros de pixel: origem é canto sup-esq (c,f)
        dx, dy = a, e  # dy costuma ser negativo em GeoTIFF norte-acima
        x = c + (np.arange(width) + 0.5) * dx
        y = f + (np.arange(height) + 0.5) * dy

        da = xr.DataArray(
            data,
            dims=("y", "x"),
            coords={"y": y, "x": x},
            name="band1",
            attrs={
                "crs": crs.to_string() if crs is not None else None,
                "transform": (a, b, c, d, e, f),
                "nodata": nodata,
            },
        )

    else:
        # Transformação com rotação -> coordenadas dependem de linha/coluna
        rows = np.arange(height)
        cols = np.arange(width)
        # grade de índices
        C, R = np.meshgrid(cols, rows)
        # centros de pixel: (col+0.5, row+0.5)
        X = a*(C+0.5) + b*(R+0.5) + c
        Y = d*(C+0.5) + e*(R+0.5) + f

        da = xr.DataArray(
            data,
            dims=("y", "x"),
            coords={"y": (("y","x"), Y), "x": (("y","x"), X)},
            name="band1",
            attrs={
                "crs": crs.to_string() if crs is not None else None,
                "transform": (a, b, c, d, e, f),
                "nodata": nodata,
            },
        )

    # Garante eixo y crescente (sul→norte) — se vier decrescente, inverte dados e coords
    if np.ndim(da.coords["y"]) == 1 and da.y.size > 1 and da.y[0] > da.y[-1]:
        da = da.isel(y=slice(None, None, -1))
    elif np.ndim(da.coords["y"]) == 2:
        # 2D: ordena linhas se necessário
        if da.y[0,0] > da.y[-1,0]:
            da = da.isel(y=slice(None, None, -1))

    return da
