module DQGMRES

# Write your package code here.
using LinearAlgebra

"""
#Input arguments

A : a square matrix of size mxn
b : a vector of dimension n
k : the number of previous vectors to orthogonalize against

Q, H = Arnoldi(A,b,k)

Q : unitary matrix containing orthonormal columns basis for Krylov space of dimension m x (k+1)
H : Household matrix of size k+1 x k


"""

function Arnoldi(A,b,k)
    m,m = size(A)

    Q = Matrix(undef,m,k+1)
    H = zeros(Float64,k+1,k)

    num_iterations = 10000

    Q[:,1] = b/norm

    for n in 1:num_iterations
        v = A * view(Q,:,n)

        start_idx = max(1,n-k+1)
        for  j in start_idx:n
            H[j,n] = view(Q,:,j) * v
            v  = v - H[j,n]* view(Q,:,j)
        end
        H[n+1,n] = norm(2,v)
        Q[:,n+1] = v ./ H[n+1,n]
    end

    return Q,H

end

function dqgmres(A,b)
      m,n = size(A)
      if n != m
        throw(DimensionMismatch("A has to be a square matrix!"))
      end
      x0 = zeros(Float64,n)
      r0 = b -A*x0
      γ_0

end


end
