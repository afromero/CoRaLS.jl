# CoRaLS Slurm rate sweeps

`generic_array.sh` submits array jobs which call `generic_acceptance.jl`.  The
shell script chooses the variable for each task, and the Julia script prints a
CSV-like rate table to the task's standard output.

Before the first submission on OSC, put this branch on the cluster, select
Julia 1.11 (the exact OSC module name may vary), instantiate the project, and
create the log directories from the project root:

```bash
module load julia/1.11
cd "$HOME/BeattyLab/CoRaLS.jl"
julia --project=. -e 'using Pkg; Pkg.instantiate()'
mkdir -p slurm/out slurm/err
```

The `slurm/out` and `slurm/err` directories must exist before `sbatch` is
called, because Slurm opens those log paths before the job script starts. Submit
from the project root as in the commands below.

## Slope-model array

The following submits four otherwise-identical jobs.  Array task 1--4 select
the corresponding entries in `SLOPE_MODELS` in order:

```bash
cd "$HOME/BeattyLab/CoRaLS.jl"
sbatch --array=1-4 \
  --export=ALL,CORALS_PROJECT_DIR="$PWD",VAR=SLOPE,SLOPE_MODELS=no_slope:gaussian_7p6:rayleigh_5p37:data_5m,ALT=10,ENERGY=1,ICE=5,ANT=4,TRIG=4,ANG=-90,FREQ1=300,TEXP=7 \
  slurm/generic_array.sh
```

`SLOPE_MODELS` uses colons, rather than commas, so it is safe inside Slurm's
comma-separated `--export` argument.  The canonical choices are:

- `no_slope` — perfectly smooth surface.
- `gaussian_0` — legacy smooth Gaussian model.
- `gaussian_7p6` — half-normal Gaussian polar slope with 7.6° sigma.
- `rayleigh_5p37` — Rayleigh polar slope with 5.37° scale.
- `data_5m` — empirical 5 m/pixel distribution.

`data_5m` is read from the small repository file
`data/south_polar_5m_slope_distribution.csv`; no terrain raster is loaded by
the jobs.  It represents an area-weighted distribution over valid 87--90°S
terrain in the source analysis.  It is not PSR-only and does not model spatial
correlations between events.  Each sampled slope receives an isotropic azimuth.

`ALT` retains the legacy convention in `generic_acceptance.jl`: its supplied
value is multiplied by 5 km.  Thus `ALT=10` means a 50 km trigger altitude.
The fixed spacecraft position remains the existing `FixedPlatform(-80, 0,
50km)` configuration.

## Other existing array sweeps

Set `VAR=ALT` or `VAR=ICE` to have the task ID replace that variable.  These
older sweeps still use `gaussian_0` unless a fixed `SLOPE_MODEL` is also passed
in the environment.  For example:

```bash
sbatch --array=1-20 \
  --export=ALL,CORALS_PROJECT_DIR="$PWD",VAR=ALT,ENERGY=1,ICE=5,ANT=4,TRIG=4,ANG=-90,FREQ1=300,TEXP=7,SLOPE_MODEL=data_5m \
  slurm/generic_array.sh
```

Each `generic_acceptance.jl` output row now includes the requested slope-model
name, so results from a slope sweep can be concatenated or filtered without
recovering the task-to-model mapping from submission metadata.

Rate sweeps intentionally do not save event payloads.  This avoids a shared
JLD2 output path being overwritten by concurrent array tasks.  If events are
needed for a deliberate single run, set `CORALS_SAVE_EVENTS=true` and provide a
unique `CORALS_SAVEFILE=/path/to/file.jld2` for that job.

For a direct, non-array invocation, the final slope-model argument is optional;
omitting it preserves the prior `gaussian_0` behavior:

```bash
julia --project=. slurm/generic_acceptance.jl 10 1 5 4 4 -90 300 7 data_5m
```
