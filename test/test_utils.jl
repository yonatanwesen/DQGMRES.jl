#copied from https://github.com/JuliaSmoothOptimizers/Krylov.jl/blob/main/test/test_utils.jl
using SparseArrays
eye(n::Int; FC=Float64) = sparse(one(FC) * I, n, n)

function symmetric_definite(n :: Int=10; FC=Float64)
    α = FC <: Complex ? FC(im) : one(FC)
    A = spdiagm(-1 => α * ones(FC, n-1), 0 => 4 * ones(FC, n), 1 => conj(α) * ones(FC, n-1))
    b = A * FC[1:n;]
    return A, b
end

# Symmetric and indefinite systems.
function symmetric_indefinite(n :: Int=10; FC=Float64)
    α = FC <: Complex ? FC(im) : one(FC)
    A = spdiagm(-1 => α * ones(FC, n-1), 0 => ones(FC, n), 1 => conj(α) * ones(FC, n-1))
    b = A * FC[1:n;]
    return A, b
  end
  
  # Nonsymmetric and positive definite systems.
  function nonsymmetric_definite(n :: Int=10; FC=Float64)
    if FC <: Complex
      A = [i == j ? n * one(FC) : FC(im) * one(FC) for i=1:n, j=1:n]
    else
      A = [i == j ? n * one(FC) : i < j ? one(FC) : -one(FC) for i=1:n, j=1:n]
    end
    b = A * FC[1:n;]
    return A, b
  end
  
  # Nonsymmetric and indefinite systems.
  function nonsymmetric_indefinite(n :: Int=10; FC=Float64)
    if FC <: Complex
      A = [i == j ? n * (-one(FC))^(i*j) : FC(im) * one(FC) for i=1:n, j=1:n]
    else
      A = [i == j ? n * (-one(FC))^(i*j) : i < j ? one(FC) : -one(FC) for i=1:n, j=1:n]
    end
    b = A * FC[1:n;]
    return A, b
  end

  # Model Poisson equation in polar coordinates
function polar_poisson(n, m, f, g; R=1.0)
    Δr = 2 * R / (2*n + 1)
    r = [(i - 1/2) * Δr for i = 1 : n+1]
  
    Δθ = 2 * π / m
    θ = [(j - 1) * Δθ for j = 1 : m+1]
  
    λ = [1 / (2 * (k - 1/2)) for k = 1 : n]
    β = [1 / ((k - 1/2)^2 * Δθ^2) for k = 1 : n]
  
    D = spdiagm(0 => β)
    T = spdiagm(-1 => 1.0 .- λ[2:n], 0 => -2.0 * ones(n), 1 => 1.0 .+ λ[1:n-1])
  
    A = spzeros(n * m, n * m)
    for k = 1 : m
      A[1+(k-1)*n : k*n, 1+(k-1)*n : k*n] = T - 2*D
      if k ≤ m-1
        A[1+k*n : (k+1)*n, 1+(k-1)*n : k*n] = D
        A[1+(k-1)*n : k*n, 1+k*n : (k+1)*n] = D
      end
    end
    A[1+(m-1)*n : m*n, 1 : n] = D
    A[1 : n, 1+(m-1)*n : m*n] = D
  
    b = zeros(n * m)
    for i = 1 : n
      for j = 1 : m
        b[i + n*(j-1)] = Δr * Δr * f(r[i], θ[j])
        if i == n
          b[i + n*(j-1)] -= (1.0 + λ[n]) * g(R, θ[j])
        end
      end
    end
  
    return A, b
  end

  function square_preconditioned(n :: Int=10; FC=Float64)
    A   = ones(FC, n, n) + (n-1) * eye(n)
    b   = 10 * FC[1:n;]
    M⁻¹ = FC(1/n) * eye(n)
    return A, b, M⁻¹
  end