
export A_mul_B!

function A_mul_B!( alpha::Float64,
                   A::SparseMatrixCSC{Float64,Ti},
                   x::Array{Float64},
                   beta::Float64,
                   y::Array{Float64},
                   nthreads::Int64=0 ) where Ti
# Real:  y = beta*y  +  alpha * A*x 

   if nthreads == 0
      #Base.A_mul_B!( alpha, A, x, beta, y )
	  mul!(y,A,x, alpha, beta)
      return
   elseif nthreads < 1
      throw(ArgumentError("nthreads < 1"))
   end

	n,m = size(A)
    nvec = size(x,2)

   if size(x,1) != m || size(y,1) != n 
      throw(DimensionMismatch("length(x) != m || length(y) != n"))
   elseif size(y,2) != nvec
      throw(DimensionMismatch("length(y,2) != nvec"))
   end
   
   
	if Ti == Int32
      p = ccall( (:a_mul_b_rr_32_, spmatveclib),
		       Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}),
             nthreads, nvec, m, n,    alpha,   beta,              A.nzval,      A.rowval,   A.colptr,   x,   y)
   elseif Ti == Int64
      p = ccall( (:a_mul_b_rr_64_, spmatveclib),
		       Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}),
             nthreads, nvec, m, n,    alpha,   beta,              A.nzval,      A.rowval,   A.colptr,   x,   y) 
   else
      error("Unsupported sparse matrix indexing integer type $Ti")
   end
   return p
end  # function A_mul_B!

#------------------------------------------------------------------------------

function A_mul_B!( alpha::ComplexF64,
                   A::SparseMatrixCSC{Float64,Ti},
                   x::Array{ComplexF64},
                   beta::ComplexF64,
                   y::Array{ComplexF64},
                   nthreads::Int64=0 ) where Ti
# Real, Complex A:  y = beta*y  +  alpha * A*x 

   if nthreads == 0
      mul!(y,A,x, alpha, beta) # Base.A_mul_B!( alpha, complex(A), x, beta, y )
      return
   elseif nthreads < 1
      throw(ArgumentError("nthreads < 1"))
   end

	n,m  = size(A)
   nvec = size(x,2)

   if size(x,1) != m || size(y,1) != n 
      throw(DimensionMismatch("length(x) != m || length(y) != n"))
   elseif size(y,2) != nvec
      throw(DimensionMismatch("length(y,2) != nvec"))
   end
   
	if Ti == Int32
	   p = ccall( (:a_mul_b_rc_32_, spmatveclib),
		 Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{ComplexF64}, Ref{ComplexF64}),
                   nthreads, nvec, m, n,     alpha,   beta,              A.nzval,      A.rowval,   A.colptr,   x,  y)
   elseif Ti == Int64
      p = ccall( (:a_mul_b_rc_64_, spmatveclib),
		 Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}),
                   nthreads, nvec, m, n,     alpha,   beta,              A.nzval,      A.rowval,   A.colptr,   x,  y)
   else
      error("Unsupported sparse matrix indexing integer type $Ti")
   end
   return p
end  # function A_mul_B!

#------------------------------------------------------------------------------

function A_mul_B!( alpha::ComplexF64,
                   A::SparseMatrixCSC{ComplexF64,Ti},
                   x::Array{ComplexF64},
                   beta::ComplexF64,
                   y::Array{ComplexF64},
                   nthreads::Int64=0 ) where Ti
# Complex A:  y = beta*y  +  alpha * A*x 

   if nthreads == 0
      mul!(y,A,x, alpha, beta) # Base.A_mul_B!( alpha, A, x, beta, y )
      return
   elseif nthreads < 1
      throw(ArgumentError("nthreads < 1"))
   end

	n,m  = size(A)
   nvec = size(x,2)

   if size(x,1) != m || size(y,1) != n 
      throw(DimensionMismatch("length(x) != m || length(y) != n"))
   elseif size(y,2) != nvec
      throw(DimensionMismatch("length(y,2) != nvec"))
   end
   
	if Ti == Int32
	   p = ccall( (:a_mul_b_cc_32_, spmatveclib),
		   Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{Int32}, Ref{Int32}, Ref{ComplexF64}, Ref{ComplexF64}),
                   nthreads, nvec, m, n,     alpha,   beta,            A.nzval,      A.rowval,   A.colptr,  x,  y)
   elseif Ti == Int64
      p = ccall( (:a_mul_b_cc_64_, spmatveclib),
		 Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}),
                   nthreads, nvec, m, n,     alpha,   beta,            A.nzval,      A.rowval,   A.colptr,  x,  y)
   else
      error("Unsupported sparse matrix indexing integer type $Ti")
   end
   return p
end  # function A_mul_B!
