#!/usr/bin/env python3
"""Measure fixed-point compression error for an LDSM slope-map NPZ.

The source NPZ stores Float32 slope angles. Rounding while retaining Float32
does not reduce the array's raw size; this diagnostic also evaluates quantized
UInt16 and UInt8 encodings with an explicit scale and a nodata sentinel.
Candidate maps are written only to temporary files for size measurements.
"""

from __future__ import annotations

import argparse
import csv
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = REPOSITORY_ROOT / "data" / "ldsm_87s_slope_40mpp.npz"
DEFAULT_OUTPUT_DIR = REPOSITORY_ROOT / "data" / "terrain_diagnostics" / "quantization"


def compressed_size(payload: dict[str, np.ndarray]) -> int:
    """Return NPZ size without leaving a candidate map in the repository."""
    with tempfile.NamedTemporaryFile(suffix=".npz") as handle:
        np.savez_compressed(handle.name, **payload)
        return Path(handle.name).stat().st_size


def fixed_point(
    slope_deg: np.ndarray, valid: np.ndarray, step_deg: float, dtype: np.dtype
) -> tuple[np.ndarray, np.ndarray]:
    """Encode finite slope values and reconstruct them in degrees."""
    info = np.iinfo(dtype)
    nodata_code = info.max
    codes = np.full(slope_deg.shape, nodata_code, dtype=dtype)
    encoded = np.rint(slope_deg[valid] / step_deg)
    if encoded.max() > nodata_code - 1:
        raise ValueError(f"{step_deg:g}° step does not fit {dtype}")
    codes[valid] = encoded.astype(dtype)

    restored = np.full(slope_deg.shape, np.nan, dtype=np.float32)
    restored[valid] = codes[valid].astype(np.float32) * step_deg
    return codes, restored


def candidate_summary(name: str, residual: np.ndarray, compressed_bytes: int) -> dict:
    return {
        "encoding": name,
        "compressed_mib": compressed_bytes / 2**20,
        "mean_residual_deg": float(np.mean(residual)),
        "rms_residual_deg": float(np.sqrt(np.mean(residual**2))),
        "max_abs_residual_deg": float(np.max(np.abs(residual))),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    table = np.load(args.source)
    slope_deg = table["slope_deg"]
    if slope_deg.dtype != np.float32 or slope_deg.ndim != 2:
        raise ValueError(f"expected a two-dimensional Float32 map, got {slope_deg.dtype} {slope_deg.shape}")
    valid = np.isfinite(slope_deg)
    metadata = {key: table[key] for key in table.files if key != "slope_deg"}
    pixel_size_m = float(table["pixel_size_m"])
    map_name = f"{pixel_size_m:g}m_slope"

    # The Float32 rounded result tests compression alone. Fixed-point candidates
    # are the encodings that actually reduce uncompressed storage.
    rounded = slope_deg.copy()
    rounded[valid] = np.round(rounded[valid], decimals=1)
    rounded_one_degree = slope_deg.copy()
    rounded_one_degree[valid] = np.round(rounded_one_degree[valid])
    candidates: list[tuple[str, np.ndarray, np.ndarray, dict[str, np.ndarray]]] = [
        (
            "Float32 rounded to 0.1°",
            rounded,
            rounded[valid] - slope_deg[valid],
            {"slope_deg": rounded, **metadata},
        ),
        (
            "Float32 rounded to 1°",
            rounded_one_degree,
            rounded_one_degree[valid] - slope_deg[valid],
            {"slope_deg": rounded_one_degree, **metadata},
        ),
    ]
    for name, step_deg, dtype in (
        ("UInt16 at 0.01°", 0.01, np.dtype("uint16")),
        ("UInt16 at 0.1°", 0.1, np.dtype("uint16")),
        ("UInt8 at 0.5°", 0.5, np.dtype("uint8")),
        ("UInt8 at 1°", 1.0, np.dtype("uint8")),
    ):
        codes, restored = fixed_point(slope_deg, valid, step_deg, dtype)
        candidates.append((
            name,
            restored,
            restored[valid] - slope_deg[valid],
            {
                "slope_code": codes,
                "slope_step_deg": np.float64(step_deg),
                "nodata_code": np.array(np.iinfo(dtype).max, dtype=dtype),
                **metadata,
            },
        ))

    rows = [{
        "encoding": "Float32 original",
        "compressed_mib": args.source.stat().st_size / 2**20,
        "mean_residual_deg": 0.0,
        "rms_residual_deg": 0.0,
        "max_abs_residual_deg": 0.0,
    }]
    residuals = {}
    restored_maps = {}
    for name, restored, residual, payload in candidates:
        rows.append(candidate_summary(name, residual, compressed_size(payload)))
        residuals[name] = residual
        restored_maps[name] = restored

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / f"{map_name}_quantization_summary.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(2, 4, figsize=(18, 7), constrained_layout=True)
    plot_names = (
        "UInt16 at 0.01°",
        "UInt16 at 0.1°",
        "UInt8 at 0.5°",
        "UInt8 at 1°",
    )
    for column, name in enumerate(plot_names):
        residual = residuals[name]
        limit = np.max(np.abs(residual))
        axes[0, column].hist(residual, bins=100, density=True, color="C0")
        axes[0, column].set(title=name, xlabel="Residual [deg]", ylabel="Density")
        axes[0, column].grid(alpha=0.3, linestyle=":")

        stride = max(1, int(np.ceil(max(slope_deg.shape) / 1000)))
        image = (restored_maps[name] - slope_deg)[::stride, ::stride]
        im = axes[1, column].imshow(image, cmap="coolwarm", vmin=-limit, vmax=limit)
        axes[1, column].set(title=f"{name} residual map", xticks=[], yticks=[])
        fig.colorbar(im, ax=axes[1, column], label="Residual [deg]")

    figure_path = args.output_dir / f"{map_name}_quantization_residuals.png"
    fig.savefig(figure_path, dpi=180)

    # A separate histogram figure is easier to inspect than the compact
    # histograms above. Quantization residuals should be nearly uniform over
    # half a quantization step, rather than Gaussian-distributed.
    histogram_fig, histogram_axes = plt.subplots(
        1, len(plot_names), figsize=(18, 3.8), constrained_layout=True
    )
    for axis, name in zip(histogram_axes, plot_names):
        axis.hist(residuals[name], bins=200, density=True, color="C0")
        axis.set(title=name, xlabel="Residual [deg]", ylabel="Density")
        axis.grid(alpha=0.3, linestyle=":")
    histogram_path = args.output_dir / f"{map_name}_quantization_residual_histograms.png"
    histogram_fig.savefig(histogram_path, dpi=180)

    print(f"wrote {csv_path}")
    print(f"wrote {figure_path}")
    print(f"wrote {histogram_path}")
    for row in rows:
        print(
            f"{row['encoding']}: {row['compressed_mib']:.2f} MiB; "
            f"RMS={row['rms_residual_deg']:.6f}°; "
            f"max |residual|={row['max_abs_residual_deg']:.6f}°"
        )


if __name__ == "__main__":
    main()
