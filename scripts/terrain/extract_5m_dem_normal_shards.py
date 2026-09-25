#!/usr/bin/env python3
"""Create bounded-memory, aspect-aware surface-normal shards from the 5 m DEM.

Each valid output pixel stores a 1 degree slope magnitude and the azimuth of
the *outward-normal tilt* in the south-polar stereographic x/y grid.  The
azimuth is quantized to 2 degrees.  Together those values reproduce a local
surface normal, rather than the slope-only maps' randomized azimuth.

The source DEM is much too large to expand at once.  This script decodes only
the TIFF tiles needed for one 2048 x 2048 output shard plus a one-pixel halo.
It writes a regular NPZ shard grid over the square containing 87--90 degrees
south, and accumulates a slope-versus-azimuth histogram while doing so.
"""

from __future__ import annotations

import argparse
import math
from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = REPOSITORY_ROOT / "data" / "ldem_87s_5mpp_normal_1deg_2deg_south87_shards"
DEFAULT_DIAGNOSTICS = REPOSITORY_ROOT / "data" / "terrain_diagnostics"
DEFAULT_FIGURE = DEFAULT_DIAGNOSTICS / "dem_5m_slope_aspect_histogram.png"
DEFAULT_HISTOGRAM = DEFAULT_DIAGNOSTICS / "dem_5m_slope_aspect_histogram.npz"

SOURCE_PAGE = 0
TILE_SIZE_PX = 256
SHARD_SIZE_PX = 2048
MIN_LATITUDE_DEG = -87.0
RADIUS_M = 1_737_400.0
SLOPE_STEP_DEG = 1.0
ASPECT_STEP_DEG = 2.0
ASPECT_BINS = round(360.0 / ASPECT_STEP_DEG)
NODATA_CODE = np.uint16(65535)


def polar_cap_crop(
    height: int, width: int, pixel_size_m: float, origin_x_m: float, origin_y_m: float
) -> tuple[int, int, int, int]:
    """Return a TIFF-tile-aligned square containing the configured cap."""
    cap_radius_m = 2.0 * RADIUS_M * math.tan(
        math.radians(90.0 + MIN_LATITUDE_DEG) / 2.0
    )
    center_col = (0.0 - origin_x_m) / pixel_size_m
    center_row = (origin_y_m - 0.0) / pixel_size_m
    radius_px = math.ceil(cap_radius_m / pixel_size_m)
    row_start = max(0, math.floor((center_row - radius_px) / TILE_SIZE_PX) * TILE_SIZE_PX)
    col_start = max(0, math.floor((center_col - radius_px) / TILE_SIZE_PX) * TILE_SIZE_PX)
    row_stop = min(height, math.ceil((center_row + radius_px) / TILE_SIZE_PX) * TILE_SIZE_PX)
    col_stop = min(width, math.ceil((center_col + radius_px) / TILE_SIZE_PX) * TILE_SIZE_PX)
    return int(row_start), int(row_stop), int(col_start), int(col_stop)


class TileReader:
    """Decode TIFF tiles on demand, retaining a small LRU cache."""

    def __init__(self, tif: tifffile.TiffFile, page: tifffile.TiffPage, max_tiles: int = 128):
        self.tif = tif
        self.page = page
        self.max_tiles = max_tiles
        self.tiles_per_row = math.ceil(page.imagewidth / page.tilewidth)
        self.cache: OrderedDict[tuple[int, int], np.ndarray] = OrderedDict()

    def tile(self, tile_row: int, tile_col: int) -> np.ndarray:
        key = (tile_row, tile_col)
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        index = tile_row * self.tiles_per_row + tile_col
        self.tif.filehandle.seek(self.page.dataoffsets[index])
        compressed = self.tif.filehandle.read(self.page.databytecounts[index])
        decoded, _, _ = self.page.decode(compressed, index)
        tile = np.asarray(decoded).squeeze().astype(np.float32, copy=False)
        self.cache[key] = tile
        if len(self.cache) > self.max_tiles:
            self.cache.popitem(last=False)
        return tile

    def window(self, row_start: int, row_stop: int, col_start: int, col_stop: int) -> np.ndarray:
        """Read a rectangular DEM window without materializing the full raster."""
        if row_start < 0 or col_start < 0 or row_stop > self.page.imagelength or col_stop > self.page.imagewidth:
            raise ValueError("requested DEM window is outside the source raster")
        output = np.empty((row_stop - row_start, col_stop - col_start), dtype=np.float32)
        first_tile_row = row_start // self.page.tilelength
        last_tile_row = (row_stop - 1) // self.page.tilelength
        first_tile_col = col_start // self.page.tilewidth
        last_tile_col = (col_stop - 1) // self.page.tilewidth
        for tile_row in range(first_tile_row, last_tile_row + 1):
            tile_top = tile_row * self.page.tilelength
            tile_bottom = min(tile_top + self.page.tilelength, self.page.imagelength)
            source_top, source_bottom = max(row_start, tile_top), min(row_stop, tile_bottom)
            for tile_col in range(first_tile_col, last_tile_col + 1):
                tile_left = tile_col * self.page.tilewidth
                tile_right = min(tile_left + self.page.tilewidth, self.page.imagewidth)
                source_left, source_right = max(col_start, tile_left), min(col_stop, tile_right)
                tile = self.tile(tile_row, tile_col)
                output[
                    source_top - row_start:source_bottom - row_start,
                    source_left - col_start:source_right - col_start,
                ] = tile[
                    source_top - tile_top:source_bottom - tile_top,
                    source_left - tile_left:source_right - tile_left,
                ]
        return output


def encode_normal_shard(elevation_with_halo: np.ndarray, pixel_size_m: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Central-difference DEM gradients, packed normal codes, and valid mask.

    Array rows increase toward decreasing projected y.  Hence the y derivative
    has the opposite array-row sign.  `azimuth` is the tilt direction of the
    outward normal (downhill), measured counter-clockwise from projected +x.
    """
    center = elevation_with_halo[1:-1, 1:-1]
    left, right = elevation_with_halo[1:-1, :-2], elevation_with_halo[1:-1, 2:]
    above, below = elevation_with_halo[:-2, 1:-1], elevation_with_halo[2:, 1:-1]
    valid = np.isfinite(center) & np.isfinite(left) & np.isfinite(right) & np.isfinite(above) & np.isfinite(below)
    dz_dx = (right - left) / (2.0 * pixel_size_m)
    dz_dy = (above - below) / (2.0 * pixel_size_m)
    slope_deg = np.degrees(np.arctan(np.hypot(dz_dx, dz_dy)))
    normal_azimuth_deg = np.degrees(np.arctan2(-dz_dy, -dz_dx)) % 360.0

    # Quantize finite pixels only.  Converting NaNs to integer emits warnings
    # and is unnecessary because those entries are written as NODATA below.
    slope_code = np.zeros(center.shape, dtype=np.int32)
    aspect_code = np.zeros(center.shape, dtype=np.int32)
    slope_code[valid] = np.rint(slope_deg[valid] / SLOPE_STEP_DEG).astype(np.int32)
    aspect_code[valid] = (
        np.rint(normal_azimuth_deg[valid] / ASPECT_STEP_DEG).astype(np.int32) % ASPECT_BINS
    )
    if np.any(slope_code[valid] < 0) or np.any(slope_code[valid] > 90):
        raise ValueError("DEM gradient produced a slope outside 0--90 degrees")
    packed = np.full(center.shape, NODATA_CODE, dtype=np.uint16)
    packed[valid] = (slope_code[valid] * ASPECT_BINS + aspect_code[valid]).astype(np.uint16)
    return packed, slope_code, aspect_code


def shard_path(directory: Path, shard_row: int, shard_col: int) -> Path:
    return directory / f"shard_r{shard_row:03d}_c{shard_col:03d}.npz"


def write_histogram_figure(histogram: np.ndarray, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(10, 5), constrained_layout=True)
    image = axis.imshow(
        np.log10(histogram.T + 1.0), origin="lower", aspect="auto",
        extent=(0.0, 91.0, 0.0, 360.0), interpolation="nearest", cmap="magma",
    )
    axis.set_xlabel("DEM slope magnitude [degrees, 1° bins]")
    axis.set_ylabel("outward-normal tilt azimuth in projected grid [degrees, 2° bins]")
    axis.set_title("5 m south-polar DEM: slope magnitude versus normal direction")
    colorbar = fig.colorbar(image, ax=axis)
    colorbar.set_label(r"$\log_{10}(N + 1)$")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="source 5 m LDEM BigTIFF")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--histogram", type=Path, default=DEFAULT_HISTOGRAM)
    parser.add_argument("--overwrite", action="store_true", help="replace existing shard files")
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    if not args.source.is_file():
        raise FileNotFoundError(args.source)
    if args.output.exists() and any(args.output.iterdir()) and not args.overwrite:
        raise FileExistsError(f"{args.output} already has files; use --overwrite to replace them")
    args.output.mkdir(parents=True, exist_ok=True)
    if args.overwrite:
        for path in args.output.glob("shard_r*_c*.npz"):
            path.unlink()

    histogram = np.zeros((91, ASPECT_BINS), dtype=np.int64)
    with tifffile.TiffFile(args.source) as tif:
        page = tif.pages[SOURCE_PAGE]
        if page.shape != (40000, 40000) or not page.is_tiled or page.samplesperpixel != 1:
            raise RuntimeError(f"expected a single-band tiled 40000x40000 source page, got {page.shape}")
        if (page.tilelength, page.tilewidth) != (TILE_SIZE_PX, TILE_SIZE_PX):
            raise RuntimeError(f"expected {TILE_SIZE_PX}x{TILE_SIZE_PX} source tiles")
        pixel_size_m = float(page.tags[33550].value[0])
        tie_point = page.tags[33922].value
        origin_x_m, origin_y_m = float(tie_point[3]), float(tie_point[4])
        row_start, row_stop, col_start, col_stop = polar_cap_crop(
            page.imagelength, page.imagewidth, pixel_size_m, origin_x_m, origin_y_m
        )
        crop_height, crop_width = row_stop - row_start, col_stop - col_start
        reader = TileReader(tif, page)
        n_shard_rows = math.ceil(crop_height / SHARD_SIZE_PX)
        n_shard_cols = math.ceil(crop_width / SHARD_SIZE_PX)

        for shard_row in range(n_shard_rows):
            core_row_start = row_start + shard_row * SHARD_SIZE_PX
            core_row_stop = min(core_row_start + SHARD_SIZE_PX, row_stop)
            for shard_col in range(n_shard_cols):
                core_col_start = col_start + shard_col * SHARD_SIZE_PX
                core_col_stop = min(core_col_start + SHARD_SIZE_PX, col_stop)
                elevation = reader.window(
                    core_row_start - 1, core_row_stop + 1, core_col_start - 1, core_col_stop + 1
                )
                packed, slope_code, aspect_code = encode_normal_shard(elevation, pixel_size_m)
                np.savez_compressed(shard_path(args.output, shard_row, shard_col), normal_code=packed)
                valid = packed != NODATA_CODE
                np.add.at(histogram, (slope_code[valid], aspect_code[valid]), 1)
            print(f"wrote DEM-normal shard row {shard_row + 1}/{n_shard_rows}", flush=True)

    np.savez(
        args.output / "manifest.npz",
        shape=np.array((crop_height, crop_width), dtype=np.int64),
        shard_size_px=np.int64(SHARD_SIZE_PX),
        pixel_size_m=np.float64(pixel_size_m),
        origin_x_m=np.float64(origin_x_m + col_start * pixel_size_m),
        origin_y_m=np.float64(origin_y_m - row_start * pixel_size_m),
        radius_m=np.float64(RADIUS_M),
        min_latitude_deg=np.float64(MIN_LATITUDE_DEG),
        slope_step_deg=np.float64(SLOPE_STEP_DEG),
        aspect_step_deg=np.float64(ASPECT_STEP_DEG),
        aspect_bins=np.int64(ASPECT_BINS),
        nodata_code=NODATA_CODE,
        source_page=np.int64(SOURCE_PAGE),
        source_tile_size_px=np.int64(TILE_SIZE_PX),
    )
    args.histogram.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.histogram, count=histogram, slope_edges_deg=np.arange(92), aspect_edges_deg=np.arange(ASPECT_BINS + 1) * ASPECT_STEP_DEG)
    write_histogram_figure(histogram, args.figure)
    print(f"wrote {args.output / 'manifest.npz'}")
    print(f"wrote {args.histogram} and {args.figure}")
    print(f"crop={crop_height}x{crop_width}; packed raw payload is {crop_height * crop_width * 2 / 2**30:.2f} GiB")


if __name__ == "__main__":
    main()
