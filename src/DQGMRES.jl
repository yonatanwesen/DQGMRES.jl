#module DQGMRES

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
#export dqgmres


function Arnoldi(A,b,k)
    m,m = size(A)

    Q = Matrix(undef,m,k+1)
    H = zeros(Float64,k+1,k)

    num_iterations = 10000

    Q[:,1] = b/norm(b,2)

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

"""
#Input arguments

A : a square matrix of size mxn
b : a vector of dimension n
k : the number of previous vectors to orthogonalize against
M : Right Preconditioner
N : Left Preconditioner

x, num_iterations = dqgmres(A,b,k;M,N)

x : solution for system Ax = b
num_iterations : num_iterations


"""
const FloatOrComplex{T} = Union{T, Complex{T}} where T <: AbstractFloat


function my_dqgmres(A,b::AbstractVector{T},k; M=I,N = I) where T <: AbstractFloat
      m,n = size(A)
      if n != m
        throw(DimensionMismatch("A has to be a square matrix!"))
      end

      if length(b) != n
        throw(DimensionMismatch("the length of b doesn't match the input of A"))
      end
      #setting up workspace
      k_c = min(k,m)
      P = zeros(T,m,k_c)
      Q = zeros(T,m,k_c)
      H = zeros(T,k_c+1)
      c = zeros(T,k_c)
      s = zeros(T,k_c)

      #iteration counter
      count=0
      cmax = 2*n
      abstol = √eps(T)
      reltol = √eps(T)
    
      #check if the left and right preconditioning are identity or not
      NisI = (N==I)
      MisI = (M==I)

      
      x = zeros(T,n)
      r0 = MisI ? b - A*x : M*(b-A*x)
      r_norm= norm(r0,2)
      tol = abstol + r_norm*reltol
      γ_0 = r_norm

      
      Q[:,1] .= r0 ./ γ_0
      #w = b

      while (r_norm > tol && count < cmax)
         #update the iteration counter
         count = count + 1
         # pos for v_i in circular stacks
         j = mod(count-1,k_c) + 1 #index corresponding to v_j,p_j in circular stacks P and V
         j_next = mod(count,k_c) + 1 # index corresponding to v_j+1

         #do the incomplete Arnoldi procedure
         w = NisI ? A*Q[:,j] : A*N*Q[:,j]
         w = MisI ? w : M*w
         for i = max(1,count-k_c+1) : count
            ipos = mod(i-1,k_c) +1 #pos corrosponding to v_i 
            #h_ij pos in H
            ij = count - i +1 
            q_i = view(Q,:,ipos)
            H[ij] = w'*q_i
            
            w = w - H[ij]*q_i  
            #println(w) 
         end

         #compute h_j+1,j
         h_jn = norm(w,2)
         if h_jn != 0
            Q[:,j_next] .= w ./ h_jn
            #=if any(isnan.(V))
                println("got nans in v after dividing by h_j+1,j")
            end=#
         end
         #v_i'*v_j = 0  for |i-j| < k so we don't want to use r_(j-1-k),(j-1) when we compute r_(j-k,j) for j≥k+1
         if count ≥ k_c+2
            H[k_c+1] = 0.0
         end

         #update the QR factorzation of H_k
         for i = max(1,count - k_c): count-1
            ipos = mod(i-1,k_c) +1 
            ij = count - i 
            ij_n  = ij + 1
            htmp = c[ipos]*H[ij_n] + s[ipos]*H[ij]
            H[ij] = s[ipos]*H[ij_n] - c[ipos]*H[ij]
            H[ij_n] = htmp
         end
         (c[j],s[j],H[1]) = LinearAlgebra.givensAlgorithm(H[1],h_jn)
         γ_n = s[j]*γ_0
         γ_0 = c[j]*γ_0
         p_j = view(P,:,j)
         #println(p_j)
         for i = max(1,count - k_c): count-1
            ipos = mod(i-1,k_c) + 1
            ij = count - i + 1
            if ipos == j
                p_j .= p_j .* (-H[ij]) 
            else
                p_j .= p_j - H[ij]*view(P,:,ipos)
            end
         end
         #println("afterwards")
         #println(p_j)
         p_j .= NisI ? p_j + view(Q,:,j) : p_j + N*view(Q,:,j)
         p_j .= p_j ./ H[1] 
         x = x + γ_0*p_j
         r_norm = abs(γ_n)
         γ_0 = γ_n
        end
    return x,count


          
end


#end
