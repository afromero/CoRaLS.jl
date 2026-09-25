#!/usr/bin/env python3
"""Create bounded-memory 1 degree UInt8 NPZ shards for the 5 m LDSM map.

The 5 m source page is 6.1 GiB as Float32 and cannot be loaded wholesale on
this machine. Tiles are decoded one at a time, quantized to nearest degree,
and accumulated only for one 2048-pixel shard row at a time. The output crops
to the square containing the 87--90 degree south cap and uses a manifest plus
regular row/column shard names for random lookup.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import tifffile


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = REPOSITORY_ROOT / "data" / "ldsm_87s_slope_5mpp_1deg_south87_shards"
SOURCE_PAGE = 0
TILE_SIZE_PX = 256
SHARD_SIZE_PX = 2048
MIN_LATITUDE_DEG = -87.0
RADIUS_M = 1_737_400.0
NODATA_CODE = np.uint8(255)


def polar_cap_crop(
    height: int, width: int, pixel_size_m: float, origin_x_m: float, origin_y_m: float
) -> tuple[int, int, int, int]:
    """Return a tile-aligned bounding square containing the configured cap."""
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


def encode_tile(tile: np.ndarray) -> np.ndarray:
    """Round finite degree values to UInt8 codes and preserve missing pixels."""
    encoded = np.full(tile.shape, NODATA_CODE, dtype=np.uint8)
    valid = np.isfinite(tile)
    values = np.rint(tile[valid])
    if values.size and (values.min() < 0 or values.max() >= NODATA_CODE):
        raise ValueError(f"1 degree UInt8 cannot represent source values {values.min()}..{values.max()}")
    encoded[valid] = values.astype(np.uint8)
    return encoded


def shard_path(directory: Path, shard_row: int, shard_col: int) -> Path:
    return directory / f"shard_r{shard_row:03d}_c{shard_col:03d}.npz"


def flush_shard_row(output: Path, shard_row: int, shards: list[np.ndarray]) -> None:
    for shard_col, codes in enumerate(shards):
        np.savez_compressed(shard_path(output, shard_row, shard_col), slope_code=codes)


def new_shard_row(crop_height: int, crop_width: int, shard_row: int) -> list[np.ndarray]:
    shard_height = min(SHARD_SIZE_PX, crop_height - shard_row * SHARD_SIZE_PX)
    n_shard_cols = math.ceil(crop_width / SHARD_SIZE_PX)
    return [
        np.full(
            (shard_height, min(SHARD_SIZE_PX, crop_width - col * SHARD_SIZE_PX)),
            NODATA_CODE,
            dtype=np.uint8,
        )
        for col in range(n_shard_cols)
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="source 5 m LDSM BigTIFF")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not args.source.is_file():
        raise FileNotFoundError(args.source)
    args.output.mkdir(parents=True, exist_ok=True)

    with tifffile.TiffFile(args.source) as tif:
        page = tif.pages[SOURCE_PAGE]
        if page.shape != (40000, 40000) or not page.is_tiled:
            raise RuntimeError(f"expected tiled 40000x40000 source page, got {page.shape}")
        pixel_size_m = float(page.tags[33550].value[0])
        tie_point = page.tags[33922].value
        origin_x_m, origin_y_m = float(tie_point[3]), float(tie_point[4])
        if (page.tilelength, page.tilewidth) != (TILE_SIZE_PX, TILE_SIZE_PX):
            raise RuntimeError(f"expected {TILE_SIZE_PX}x{TILE_SIZE_PX} source tiles")

        row_start, row_stop, col_start, col_stop = polar_cap_crop(
            page.imagelength, page.imagewidth, pixel_size_m, origin_x_m, origin_y_m
        )
        crop_height, crop_width = row_stop - row_start, col_stop - col_start
        n_shard_cols = math.ceil(crop_width / SHARD_SIZE_PX)
        current_shard_row: int | None = None
        shards: list[np.ndarray] = []

        for data, index, _ in page.segments(maxworkers=1):
            if data is None:
                continue
            source_row, source_col = index[2], index[3]
            if not (row_start <= source_row < row_stop and col_start <= source_col < col_stop):
                continue
            shard_row = (source_row - row_start) // SHARD_SIZE_PX
            if current_shard_row is None:
                current_shard_row = shard_row
                shards = new_shard_row(crop_height, crop_width, shard_row)
            elif shard_row != current_shard_row:
                flush_shard_row(args.output, current_shard_row, shards)
                print(f"wrote shard row {current_shard_row + 1}")
                current_shard_row = shard_row
                shards = new_shard_row(crop_height, crop_width, shard_row)

            tile = encode_tile(data[0, :, :, 0])
            local_row = source_row - row_start
            local_col = source_col - col_start
            shard_col = local_col // SHARD_SIZE_PX
            shard_local_row = local_row % SHARD_SIZE_PX
            shard_local_col = local_col % SHARD_SIZE_PX
            height = min(TILE_SIZE_PX, crop_height - local_row)
            width = min(TILE_SIZE_PX, crop_width - local_col)
            shards[shard_col][
                shard_local_row:shard_local_row + height,
                shard_local_col:shard_local_col + width,
            ] = tile[:height, :width]

        if current_shard_row is None:
            raise RuntimeError("no source tiles intersected the requested cap")
        flush_shard_row(args.output, current_shard_row, shards)
        print(f"wrote shard row {current_shard_row + 1}")

    np.savez(
        args.output / "manifest.npz",
        shape=np.array((crop_height, crop_width), dtype=np.int64),
        shard_size_px=np.int64(SHARD_SIZE_PX),
        pixel_size_m=np.float64(pixel_size_m),
        origin_x_m=np.float64(origin_x_m + col_start * pixel_size_m),
        origin_y_m=np.float64(origin_y_m - row_start * pixel_size_m),
        radius_m=np.float64(RADIUS_M),
        min_latitude_deg=np.float64(MIN_LATITUDE_DEG),
        slope_step_deg=np.float64(1.0),
        nodata_code=NODATA_CODE,
        source_page=np.int64(SOURCE_PAGE),
        source_tile_size_px=np.int64(TILE_SIZE_PX),
    )
    print(f"wrote {args.output / 'manifest.npz'}")
    print(
        f"crop={crop_height}x{crop_width}, shard grid="
        f"{math.ceil(crop_height / SHARD_SIZE_PX)}x{n_shard_cols}, "
        f"working set is one shard row"
    )


if __name__ == "__main__":
    main()
