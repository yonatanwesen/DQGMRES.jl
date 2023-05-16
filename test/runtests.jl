cd(@__DIR__)
using Pkg
Pkg.activate("..")


using Test

include("test_utils.jl")
include("../src/DQGMRES.jl")

#@testset "DQGMRES.jl" begin
    # Write your tests here.
@testset "DQGMRES" begin
   test_tol = 1e-6
   for elty in (Float64,)
        k = 20
        A, b = symmetric_definite(FC = elty)
        (x, stats) = my_dqgmres(A, b,k)
        r = b - A * x
        resd = norm(r)/norm(b)
        @test (resd ≤test_tol)
        
        A, b = nonsymmetric_indefinite(FC=elty)
        (x, stats) = my_dqgmres(A, b,k)
        r = b - A * x
        resd = norm(r)/norm(b)
        @test (resd ≤ test_tol)

        f(r, θ) = -3.0 * cos(θ)
        g(r, θ) = 0.0
        A, b = polar_poisson(10,10, f, g)
        (x, stats) = my_dqgmres(A, b,k)
        r = b - A * x
        resd = norm(r)/norm(b)
        @test (resd ≤ test_tol)

        A, b, M = square_preconditioned(FC=elty)
        (x, stats) = my_dqgmres(A, b,k,M=M)
        r = b- A*x
        resd = norm(r)/norm(b)
        @test (resd ≤ test_tol)

        A, b, N = square_preconditioned(FC=elty)
        (x, stats) = my_dqgmres(A, b,k,N=N)
        r = b - A * x
        resd = norm(r)/norm(b)
        @test (resd ≤ test_tol)

        A,b,M,N = two_preconditioners(FC=elty)
        (x, stats) = my_dqgmres(A, b,k,M=M,N=N)
        r = b - A * x
        resd = norm(r)/norm(b)
        @test (resd ≤ test_tol)
       
   end 
end


#end
