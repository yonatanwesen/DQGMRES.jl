cd(@__DIR__)
using Pkg
Pkg.activate("..")

#using Test
include("test_utils.jl")
include("../src/DQGMRES.jl")

#@testset "DQGMRES.jl" begin
    # Write your tests here.
using Krylov
using Plots
using BenchmarkTools
f(r, θ) = -3.0 * cos(θ)
g(r, θ) = 0.0

iters = []
kry_iters = []
gmres_iters =[]

ns = 4:8:500
k = 20
for n in ns
    A,b = symmetric_definite(n)
    x,stats = my_dqgmres(A,b,k)
    x_k,stats_k = Krylov.cg_lanczos(A,b) #Krylov.dqgmres(A,b;memory =k)
    x_g,stats_g = Krylov.gmres(A,b;memory = n)
    push!(iters,stats)
    push!(kry_iters,stats_k.niter)
    push!(gmres_iters,stats_g.niter)
end

plt = plot(ns,iters,label="my_dqgmres",margin=5Plots.mm,marker=:auto)
plot!(plt,ns,kry_iters,label = "Krylov.cg_lanczos",marker =:auto)
plot!(plt,ns,gmres_iters,label = "Krylov.gmres",marker =:auto)

xlabel!("n: size of matrix")
ylabel!("time in seconds")


ylo