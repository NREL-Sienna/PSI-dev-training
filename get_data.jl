using PowerSystems
const PSY = PowerSystems
using PowerSimulations
const PSI = PowerSimulations
using Dates
using TimeSeries
using Random
rng = MersenneTwister(1234);
using Distributions
using Gurobi
using Logging
using PowerGraphics
using JSON
using DataFrames
using Interpolations
plotlyjs()

solver = optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 1e-3)

# Script to download ata
include(joinpath(pwd(), "utils.jl"))

rts_dir = download("https://github.com/GridMod/RTS-GMLC", pwd(), "master")
rts_src_dir = joinpath(rts_dir, "RTS_Data", "SourceData")
rts_siip_dir = joinpath(rts_dir, "RTS_Data", "FormattedData", "SIIP");

rawsys = PSY.PowerSystemTableData(
    rts_src_dir,
    100.0,
    joinpath(rts_siip_dir, "user_descriptors.yaml"),
    timeseries_metadata_file = joinpath(rts_siip_dir, "timeseries_pointers.json"),
    generator_mapping_file = joinpath(rts_siip_dir, "generator_mapping.yaml"),
)

sys_DA = System(rawsys; forecast_resolution = Dates.Hour(1))
sys_RT = System(rawsys; forecast_resolution = Dates.Minute(5))

#to_json(sys_DA, "RTS_1hr_sys_agc.json", force = true)
