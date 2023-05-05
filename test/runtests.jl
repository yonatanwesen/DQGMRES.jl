cd(@__DIR__)
using Pkg
Pkg.activate("..")

#using Test

include("test_utils.jl")
include("../src/DQGMRES.jl")

#@testset "DQGMRES.jl" begin
    # Write your tests here.
FC=Float64
k = 10
A, b = symmetric_definite(FC=Float64)
(x, stats) = dqgmres(A, b,k)
r = b - A * x

A, b = nonsymmetric_indefinite(FC=FC)
(x, stats) = dqgmres(A, b,k)
r = b - A * x

f(r, θ) = -3.0 * cos(θ)
g(r, θ) = 0.0
A, b = polar_poisson(50, 50, f, g)
(x, stats) = dqgmres(A, b,50)
r = b - A * x

A, b, M = square_preconditioned(FC=FC)
(x, stats) = dqgmres(A, b,10,M)
r = b- A*x

A, b, N = square_preconditioned(FC=FC)
(x, stats) = dqgmres(A, b,k,N=N)
r = b - A * x

#end
