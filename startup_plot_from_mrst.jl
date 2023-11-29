# General
using Revise
using DelimitedFiles
using StatsBase
using Base.Iterators
using MAT

# Plots
using GLMakie
GLMakie.activate!()
using Colors

# Mine
#lib_dir = "/Users/lluis/Documents/julia/jl-tools/"
lib_dir = "C:\\Users\\lsalo\\julia\\jl-tools\\"
includet(joinpath(lib_dir, "modules", "polyconnectivity.jl"))
using .PolyConnectivity
includet(joinpath(lib_dir, "modules", "plotFromMRST.jl"))
using .PlotFromMRST