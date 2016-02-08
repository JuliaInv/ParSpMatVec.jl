
export A_mul_B!

function A_mul_B!( alpha::Float64,
                   A::SparseMatrixCSC{Float64,Int},
                   x::Array{Float64},
                   beta::Float64,
                   y::Array{Float64},
                   nthreads::Int64=0 )
# Real:  y = beta*y  +  alpha * A*x 

   if nthreads == 0
      Base.A_mul_B!( alpha, A, x, beta, y )
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
   
   
	p  = ccall( (:a_mul_b_rr_, spmatveclib),
		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}),
                &nthreads, &nvec, &m, &n,    &alpha,   &beta,              A.nzval,      A.rowval,   A.colptr,   x,   y);
   
end  # function A_mul_B!

#------------------------------------------------------------------------------

function A_mul_B!( alpha::Complex128,
                   A::SparseMatrixCSC{Float64,Int},
                   x::Array{Complex128},
                   beta::Complex128,
                   y::Array{Complex128},
                   nthreads::Int64=0 )
# Real, Complex A:  y = beta*y  +  alpha * A*x 

   if nthreads == 0
      Base.A_mul_B!( alpha, complex(A), x, beta, y )
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
   
	p  = ccall( (:a_mul_b_rc_, spmatveclib),
		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}),
                   &nthreads, &nvec, &m, &n,     &alpha,   &beta,              A.nzval,      A.rowval,   A.colptr,   convert(Ptr{Complex128}, pointer(x)),  convert(Ptr{Complex128}, pointer(y)));
   
end  # function A_mul_B!

#------------------------------------------------------------------------------

function A_mul_B!( alpha::Complex128,
                   A::SparseMatrixCSC{Complex128,Int},
                   x::Array{Complex128},
                   beta::Complex128,
                   y::Array{Complex128},
                   nthreads::Int64=0 )
# Complex A:  y = beta*y  +  alpha * A*x 

   if nthreads == 0
      Base.A_mul_B!( alpha, A, x, beta, y )
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
   
	p  = ccall( (:a_mul_b_cc_, spmatveclib),
		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}, Ptr{Complex128}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}),
                   &nthreads, &nvec, &m, &n,     &alpha,   &beta,              convert(Ptr{Complex128}, pointer(A.nzval)),      A.rowval,   A.colptr,   convert(Ptr{Complex128}, pointer(x)),  convert(Ptr{Complex128}, pointer(y)));
   
end  # function A_mul_B!
