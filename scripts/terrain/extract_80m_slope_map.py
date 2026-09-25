#!/usr/bin/env python3
"""Extract the native 80 m/pixel LDSM slope-angle map.

Although its filename ends in ``5mpp``, the source LDSM BigTIFF stores slope
angles in degrees, not elevations. Overview level 4 is its native 80 m/pixel
slope map. This script writes that page as an independent, single-band
GeoTIFF so notebook work need not open or decode the 5 GB source file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import tifffile


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = REPOSITORY_ROOT / "data" / "ldsm_87s_slope_80mpp.tif"
DEFAULT_NPZ_OUTPUT = REPOSITORY_ROOT / "data" / "ldsm_87s_slope_80mpp.npz"
OVERVIEW_LEVEL = 4


def geotiff_extratags(base_page: tifffile.TiffPage, pixel_size_m: float):
    """Copy the projection tags while updating the pixel spacing."""
    tags = base_page.tags
    extratags = [
        # ModelPixelScaleTag: X, Y, Z pixel scale in projected metres.
        (33550, "d", 3, (pixel_size_m, pixel_size_m, 0.0), False),
    ]
    for code in (33922, 34735, 34736, 34737):
        tag = tags.get(code)
        if tag is not None:
            extratags.append((code, tag.dtype, tag.count, tag.value, False))
    return extratags


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="source 5 m LDSM BigTIFF")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-geotiff", action="store_true",
                        help="write only the NPZ needed by CoRaLS")
    parser.add_argument("--npz-output", type=Path, default=DEFAULT_NPZ_OUTPUT,
                        help="NPZ file for the Julia slope-map loader")
    args = parser.parse_args()

    if not args.source.is_file():
        raise FileNotFoundError(args.source)

    with tifffile.TiffFile(args.source) as tif:
        base_page = tif.pages[0]
        if len(tif.pages) <= OVERVIEW_LEVEL:
            raise RuntimeError(f"{args.source} has no overview level {OVERVIEW_LEVEL}")
        overview = tif.pages[OVERVIEW_LEVEL]
        base_pixel_size_m = float(base_page.tags[33550].value[0])
        pixel_size_m = base_pixel_size_m * (base_page.imagewidth / overview.imagewidth)
        tie_point = base_page.tags[33922].value
        origin_x_m, origin_y_m = float(tie_point[3]), float(tie_point[4])
        slope_deg = overview.asarray()
        extratags = geotiff_extratags(base_page, pixel_size_m)

    expected_shape = (2500, 2500)
    if slope_deg.shape != expected_shape or pixel_size_m != 80.0:
        raise RuntimeError(
            f"expected the 80 m overview to be {expected_shape}; got "
            f"{slope_deg.shape} at {pixel_size_m:g} m/pixel"
        )

    args.npz_output.parent.mkdir(parents=True, exist_ok=True)
    if not args.no_geotiff:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        tifffile.imwrite(
            args.output,
            slope_deg,
            dtype=slope_deg.dtype,
            compression="deflate",
            tile=(256, 256),
            metadata=None,
            extratags=extratags,
        )
    np.savez_compressed(
        args.npz_output,
        slope_deg=slope_deg,
        pixel_size_m=np.float64(pixel_size_m),
        origin_x_m=np.float64(origin_x_m),
        origin_y_m=np.float64(origin_y_m),
        radius_m=np.float64(1_737_400.0),
        source_page=np.int64(OVERVIEW_LEVEL),
    )
    if not args.no_geotiff:
        print(f"wrote {args.output}")
    print(f"wrote {args.npz_output}")
    print(f"shape={slope_deg.shape}, spacing={pixel_size_m:g} m/pixel, dtype={slope_deg.dtype}, values=slope degrees")


if __name__ == "__main__":
    main()
