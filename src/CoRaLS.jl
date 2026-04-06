module CoRaLS
using Reexport
using Unitful: g, cm, m, km, sr, μV, V, eV, GeV, EeV, Hz, MHz, W, K, NoUnits

# Export the CoRaLS Readme.md as its docstring
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) CoRaLS

# exports from geometry.jl
export Rmoon, random_point_on_cap, horizon_angle, random_direction, spherical_cap_area, random_point_on_sphere
export WholeMoonRegion, CircularRegion, PolarRegion, CustomRegion, SouthPolePSR, NorthPolePSR, AllPSR, create_region, is_in_region, region_area, parse_orbit
export FixedPlatform, CircularOrbit, SampledPositions, create_spacecraft, get_position
export random_point_in_aoi, aoi_to_spherical_bounds, is_in_aoi, is_in_region
export WholeMoonMare, WholeMoonHighlands, always_true

# exports from spectrum.jl
export auger_spectrum_2021, auger_spectrum_2020, auger_spectrum, sample_auger, yr, cr_spectrum, sample_power_law

# exports from utils.jl
export rician, polarization_angle, retrigger

# exports from fresnel.jl
export FarFieldDivergence, NearFieldDivergence, MixedFieldDivergence

# exports from ice.jl
export NoIceRoughness, GaussianIceRoughness

# exports from sky.jl
export sky_temperature

# exports from regolith.jl
export regolith_density, regolith_index
export ConstantIndex, SurfaceDeepIndex, StrangwayIndex, StrangwayIndexCB, CE3Index, ConstantDensity, StrangwayDensity, StrangwayDensityCB, CE3LPRDensity
export CE4Index, LunarSourceBookIndex, LSB_DivinerIndex,  CE4LPRDensity_Dong2020, DivinerRadiDensity, LunarSourceBookDensity, LSB_Diviner_Density
#regolith_density(CE3LPRDensity(), )

# exports from surface.jl
export NoSlope, GaussianSlope, NoRoughness, GaussianRoughness, RayleighSlope

# exports from cosmic_ray.jl
export estimate_Xmax, charge_excess_profile

# exports from efield.jl
export regolith_field, JAM, FORTE, GaisserHillasProfile, GaussianProfile, ARW

# exports from detector.jl
export LPDA, ANITA

# exports from simulate.jl
export throw_cosmicray, ScalarGeometry, VectorGeometry

# exports from trigger.jl
export magnitude_trigger, gaussian_trigger, rician_trigger, trigger_all

# exports from acceptance.jl
export acceptance, differential_spectrum, trials_passed, old_acceptance, save_acceptance, load_acceptance, merge_acceptance

# exports from plots.jl
export plot_differential_spectrum, plot_incident_angles, plot_polarization_angle, plot_offaxis_angle
export plot_acceptance, plot_rate_experiment, plot_event_outcomes

# import all of our code
include("constants.jl")
include("utils.jl")
include("geometry.jl")
include("spectrum.jl")
include("cosmic_ray.jl")
include("fresnel.jl")
include("ice.jl")
include("sky.jl")
include("regolith.jl")
include("efield.jl")
include("trigger.jl")
include("simulate.jl")
include("detector.jl")
include("acceptance.jl")
include("surface.jl")
include("plots.jl")

end
