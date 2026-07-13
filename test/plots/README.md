# Plot Visual-Diff Workflow

This directory includes a helper script that keeps Git plot diffs focused on visible changes.

## Script

- `visual_diff_test_plots.jl`

The script compares each changed PNG in `test/plots/test_plots` against its `HEAD` version and computes RMSE from pixel values.

- If RMSE is at or below the threshold, it restores the exact `HEAD` bytes (removing Git noise from metadata/compression-only changes).
- If RMSE is above threshold, it keeps the new file as a real visual change.
- For each changed image, it exports both files for manual inspection:
	- `tmp/plot_visual_compare_<timestamp>/head/...`
	- `tmp/plot_visual_compare_<timestamp>/current/...`
	The summary prints the compare directory path.

## Usage

From the repository root:

```bash
julia --project=. test/plots/visual_diff_test_plots.jl
```

Optional threshold override:

```bash
julia --project=. test/plots/visual_diff_test_plots.jl --threshold 1e-5
```

## Typical flow

1. Run the test suite:

```bash
julia --project=. test/runtests.jl
```

2. Visual diff cleanup runs automatically at the end of `test/runtests.jl`.

3. Stage only meaningful visual diffs.

To run only the visual diff step manually:

```bash
julia --project=. test/plots/visual_diff_test_plots.jl
```
