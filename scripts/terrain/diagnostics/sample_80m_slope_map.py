#!/usr/bin/env python3
"""Sample stored slope angles directly from the 80 m LDSM slope map.

Copy ``sample_slope_angles`` into a Python notebook after generating
``ldsm_87s_slope_80mpp.npz`` with ``extract_80m_slope_map.py``. It samples
finite pixel values uniformly; it does not calculate a gradient from a DEM.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SLOPE_MAP = REPOSITORY_ROOT / "data" / "ldsm_87s_slope_80mpp.npz"


def read_slope_map(path: str | Path = DEFAULT_SLOPE_MAP) -> tuple[np.ndarray, float]:
    """Return stored slope angle in degrees and projected pixel spacing."""
    with np.load(path) as table:
        slope_deg = table["slope_deg"]
        pixel_size_m = float(table["pixel_size_m"])
    if slope_deg.ndim != 2:
        raise ValueError(f"expected one slope-angle band, got shape {slope_deg.shape}")
    return slope_deg, pixel_size_m


def sample_slope_angles(
    slope_map_deg: np.ndarray,
    n: int,
    *,
    seed: int | None = None,
    isotropic_azimuth: bool = True,
) -> dict[str, np.ndarray]:
    """Uniformly sample finite stored slope angles.

    The LDSM page holds slope magnitude only, not a downhill/aspect direction.
    When ``isotropic_azimuth`` is true (the default), ``slope_azimuth_rad`` is
    therefore an independent uniform draw suitable for CoRaLS's existing
    random-surface convention. It is not a map-derived aspect.
    """
    if n < 1:
        raise ValueError("n must be positive")
    candidates = np.flatnonzero(np.isfinite(slope_map_deg))
    if len(candidates) == 0:
        raise ValueError("slope map has no finite pixels")

    rng = np.random.default_rng(seed)
    indices = rng.choice(candidates, size=n, replace=True)
    row, col = np.unravel_index(indices, slope_map_deg.shape)
    samples = {
        "row": row,
        "col": col,
        "slope_deg": slope_map_deg[row, col],
    }
    if isotropic_azimuth:
        samples["slope_azimuth_rad"] = rng.uniform(0.0, 2.0 * np.pi, size=n)
    return samples


if __name__ == "__main__":
    slope_map_deg, pixel_size_m = read_slope_map()
    samples = sample_slope_angles(slope_map_deg, n=10_000, seed=20260924)
    print(f"Slope-map shape: {slope_map_deg.shape}; pixel spacing: {pixel_size_m:g} m")
    print(f"mean slope: {samples['slope_deg'].mean():.3f} deg")
    print(f"95th percentile: {np.percentile(samples['slope_deg'], 95):.3f} deg")
