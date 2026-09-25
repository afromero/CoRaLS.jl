# South-polar terrain-map preparation

CoRaLS does not version the multi-gigabyte source GeoTIFFs or their generated
map products.  This directory contains the complete, supported preparation
workflow: download the two inputs below, run the command for the map model you
want, then open the notebook or Julia session normally.  Generated outputs go
under `data/` and are ignored by Git.

## Required input TIFFs

A 5 meter-per-pixel map of surface slopes and DEM for 87-90 S can be found here: https://pgda.gsfc.nasa.gov/products/81
CoRaLS has options to randomly sample surface slopes from a Gaussian distribution or Rayleigh distribution. The .tif 
files from the above source can be sampled instead. However, they are large files that may be too intense to carry 
in RAM while running. Instead, we introduce functions to downscale the resolution of these maps in a way that makes
them possible to run even on a local machine. 

| Input file | Contents | Used by |
| --- | --- | --- |
| `ldsm_87s_5mpp.tif` | Direct slope-angle map, with 5/40/80 m overview pages | 80 m, 40 m, and 5 m LDSM models |
| `ldem_87s_5mpp.tif` | 5 m elevation DEM | 5 m DEM-normal model |

The scripts do not assume a particular download directory.  Point `--source`
at the downloaded file.  They require Python 3 with `numpy` and `tifffile`;
the DEM-normal extractor also uses `matplotlib` to create its diagnostic plot.

```sh
python3 -m pip install numpy tifffile matplotlib

LDSM_TIFF=/path/to/ldsm_87s_5mpp.tif
LDEM_TIFF=/path/to/ldem_87s_5mpp.tif
```

## Build a map

Run these commands from the repository root.  Only build the product you plan
to use.

| Julia model | Command | Local output |
| --- | --- | --- |
| `SouthPolarSlopeMap()` | `python3 scripts/terrain/extract_80m_slope_map.py --source "$LDSM_TIFF" --no-geotiff` | `data/ldsm_87s_slope_80mpp.npz` |
| `SouthPolarSlopeMap40m()` | `python3 scripts/terrain/extract_40m_slope_map.py --source "$LDSM_TIFF" --no-geotiff` | `data/ldsm_87s_slope_40mpp.npz` |
| `ShardedSouthPolarSlopeMap()` | `python3 scripts/terrain/extract_5m_south87_quantized_shards.py --source "$LDSM_TIFF"` | `data/ldsm_87s_slope_5mpp_1deg_south87_shards/` |
| `ShardedSouthPolarDEMNormalMap()` | `python3 scripts/terrain/extract_5m_dem_normal_shards.py --source "$LDEM_TIFF"` | `data/ldem_87s_5mpp_normal_1deg_2deg_south87_shards/` |

The 5 m products cover 87--90 degrees south.  The LDSM shards store 1 degree
slope magnitudes and choose a random tilt azimuth at simulation time.  The DEM
shards store 1 degree slope magnitude and 2 degree map-derived normal direction.

Approximate compressed storage is 0.4 GB for the 5 m LDSM product and 1.4 GB
for the DEM-normal product.  Start either sharded model with
`cache_shards=32`: that uses about 128 MiB decoded for LDSM or 256 MiB for the
DEM-normal map.  Increase the cache only if memory permits.

The DEM command also writes its slope-versus-normal-direction histogram to
`data/terrain_diagnostics/`; it is diagnostic output and is ignored by Git.

## Optional diagnostics

After building an 80 m or 40 m map, the helpers in `scripts/terrain/diagnostics/`
can sample its stored slope values or measure fixed-point quantization error:

```sh
python3 scripts/terrain/diagnostics/sample_80m_slope_map.py
python3 scripts/terrain/diagnostics/quantize_slope_map.py \
  --source data/ldsm_87s_slope_40mpp.npz
```

Their figures and summaries are written under `data/terrain_diagnostics/` and
are also ignored by Git.
