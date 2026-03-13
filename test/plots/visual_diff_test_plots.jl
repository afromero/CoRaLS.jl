using PyPlot
using Statistics
using Printf
using Dates

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PLOTS_SUBDIR = "test/plots/test_plots"

"""
    changed_plot_paths(repo_root)::Vector{String}

Return absolute paths of changed PNG files under test/plots/test_plots.
"""
function changed_plot_paths(repo_root::AbstractString)::Vector{String}
    cmd = `git -C $repo_root diff --name-only -- $PLOTS_SUBDIR`
    output = read(cmd, String)
    relpaths = filter(!isempty, split(chomp(output), '\n'))
    png_relpaths = filter(path -> endswith(lowercase(path), ".png"), relpaths)
    return [normpath(joinpath(repo_root, relpath)) for relpath in png_relpaths]
end

"""
    write_head_version(repo_root, relpath, outpath)::Bool

Write the HEAD version of relpath to outpath. Returns false if the file does not
exist in HEAD (for example, a new untracked file).
"""
function write_head_version(repo_root::AbstractString, relpath::AbstractString, outpath::AbstractString)::Bool
    open(outpath, "w") do io
        try
            run(pipeline(`git -C $repo_root show HEAD:$relpath`, stdout=io, stderr=devnull))
        catch
            return false
        end
    end
    return true
end

"""
    image_rmse(path_a, path_b)::Float64

Compute root-mean-square pixel difference between two images loaded via PyPlot.
Returns Inf if the image sizes differ.
"""
function image_rmse(path_a::AbstractString, path_b::AbstractString)::Float64
    img_a = PyPlot.imread(path_a)
    img_b = PyPlot.imread(path_b)

    size(img_a) == size(img_b) || return Inf

    diff = Float64.(img_a) .- Float64.(img_b)
    return sqrt(mean(diff .^ 2))
end

"""
    restore_if_visually_same(path, threshold)::Tuple{Bool, Float64}

Compare path against HEAD:path. If RMSE <= threshold, restore exact HEAD bytes
and return (true, rmse). Otherwise return (false, rmse).
"""
function restore_if_visually_same(path::AbstractString, threshold::Float64)::Tuple{Bool, Float64}
    rel_plot_path = Base.Filesystem.relpath(path, REPO_ROOT)
    tmp_dir = mktempdir()
    tmp_head = joinpath(tmp_dir, "head.png")

    has_head = write_head_version(REPO_ROOT, rel_plot_path, tmp_head)
    if !has_head
        rm(tmp_dir; recursive=true, force=true)
        return false, Inf
    end

    rmse = image_rmse(path, tmp_head)
    if rmse <= threshold
        cp(tmp_head, path; force=true)
        rm(tmp_dir; recursive=true, force=true)
        return true, rmse
    end

    rm(tmp_dir; recursive=true, force=true)
    return false, rmse
end

function parse_threshold(args::Vector{String})::Float64
    if isempty(args)
        return 1e-6
    end

    if length(args) == 2 && args[1] == "--threshold"
        try
            return parse(Float64, args[2])
        catch
            error("Invalid threshold value: $(args[2])")
        end
    end

    error("Usage: julia --project=. test/plots/visual_diff_test_plots.jl [--threshold <float>]")
end

"""
    create_compare_dir(repo_root)::String

Create a run-specific temporary folder under repo_root/tmp for visual comparison
artifacts.
"""
function create_compare_dir(repo_root::AbstractString)::String
    ts = Dates.format(now(), "yyyymmdd_HHMMSS")
    compare_dir = joinpath(repo_root, "tmp", "plot_visual_compare_$(ts)")
    mkpath(joinpath(compare_dir, "head"))
    mkpath(joinpath(compare_dir, "current"))
    return compare_dir
end

"""
    export_compare_pair(path, compare_dir)::Bool

Write current and HEAD versions of a plot into compare_dir preserving relative
path layout. Returns true if HEAD version exists.
"""
function export_compare_pair(path::AbstractString, compare_dir::AbstractString)::Bool
    rel = Base.Filesystem.relpath(path, REPO_ROOT)
    current_out = joinpath(compare_dir, "current", rel)
    head_out = joinpath(compare_dir, "head", rel)

    mkpath(dirname(current_out))
    mkpath(dirname(head_out))

    cp(path, current_out; force=true)
    return write_head_version(REPO_ROOT, rel, head_out)
end

function main(args::Vector{String})
    threshold = parse_threshold(args)
    changed = changed_plot_paths(REPO_ROOT)

    if isempty(changed)
        println("No changed plot PNGs found under $(PLOTS_SUBDIR).")
        return
    end

    restored = String[]
    kept = String[]

    println("Checking $(length(changed)) changed plot PNG(s) with RMSE threshold $(threshold)...")

    compare_dir = ""
    exported_pairs = 0

    for path in changed
        was_restored, rmse = restore_if_visually_same(path, threshold)
        rel = Base.Filesystem.relpath(path, REPO_ROOT)

        if was_restored
            push!(restored, rel)
            @printf("RESTORED  %s  (rmse=%.3e)\n", rel, rmse)
        else
            push!(kept, rel)
            if isempty(compare_dir)
                compare_dir = create_compare_dir(REPO_ROOT)
            end
            has_head = export_compare_pair(path, compare_dir)
            exported_pairs += has_head ? 1 : 0
            if isfinite(rmse)
                @printf("CHANGED   %s  (rmse=%.3e)\n", rel, rmse)
            else
                println("CHANGED   $(rel)  (new file or different image shape)")
            end
        end
    end

    println("\nSummary:")
    println("  restored (no visible diff): $(length(restored))")
    println("  kept (visible diff):        $(length(kept))")
    if !isempty(compare_dir)
        println("  visual compare dir:         $(compare_dir)")
        println("  exported HEAD/current pairs: $(exported_pairs)")
    end
end

main(ARGS)
