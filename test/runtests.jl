using Test
using FFTResampling
using GeneralizedGenerated
using Base.Cartesian
using FFTW
using Random
Random.seed!(42)
include("utils.jl")
include("sinc_interpolation.jl")
