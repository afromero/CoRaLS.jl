# Plot Visual-Diff Workflow

This directory includes a helper script that keeps Git plot diffs focused on visible changes.

## Script

- `promote_visual_changes.jl`

The script compares each changed PNG in `test/plots/test_plots` against its `HEAD` version and computes RMSE from pixel values.

- If RMSE is at or below the threshold, it restores the exact `HEAD` bytes (removing Git noise from metadata/compression-only changes).
- If RMSE is above threshold, it keeps the new file as a real visual change.

## Usage

From the repository root:

```bash
julia --project=. test/plots/promote_visual_changes.jl
```

Optional threshold override:

```bash
julia --project=. test/plots/promote_visual_changes.jl --threshold 1e-5
```

## Typical flow

1. Run plotting tests to regenerate figures.
2. Run `promote_visual_changes.jl`.
3. Stage only meaningful visual diffs.
