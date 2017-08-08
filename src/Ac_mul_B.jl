
export Ac_mul_B!

function Ac_mul_B!( alpha::Float64,
                    A::SparseMatrixCSC{Float64,Int},
                    x::Array{Float64},
                    beta::Float64,
                    y::Array{Float64},
                    nthreads::Int64=0 )
# Real:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      Base.Ac_mul_B!( alpha, A, x, beta, y )
      return
   elseif nthreads < 1
      throw(ArgumentError("nthreads < 1"))
   end

	n,m  = size(A)
   nvec = size(x,2)

   if size(x,1) != n || size(y,1) != m 
      throw(DimensionMismatch("size(x) != n || size(y) != m"))
   elseif size(y,2) != nvec
      throw(DimensionMismatch("length(y,2) != nvec"))
   end
   
	p  = ccall( (:ac_mul_b_rr_, spmatveclib),
		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}),
                &nthreads, &nvec, &m, &n,    &alpha,   &beta,              A.nzval,      A.rowval,   A.colptr,   x,   y);
   
end  # function Ac_mul_B!

#------------------------------------------------------------------------------

function Ac_mul_B!( alpha::Complex128,
                    A::SparseMatrixCSC{Float64,Int},
                    x::Array{Complex128},
                    beta::Complex128,
                    y::Array{Complex128},
                    nthreads::Int64=0 )
# Real, Complex A:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      Base.Ac_mul_B!( alpha, A, x, beta, y )
      return
   elseif nthreads < 1
      throw(ArgumentError("nthreads < 1"))
   end

	n,m  = size(A)
   nvec = size(x,2)

   if size(x,1) != n || size(y,1) != m 
      throw(DimensionMismatch("length(x) != n || length(y) != m"))
   elseif size(y,2) != nvec
      throw(DimensionMismatch("length(y,2) != nvec"))
   end
   
	p  = ccall( (:ac_mul_b_rc_, spmatveclib),
		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}),
                   &nthreads, &nvec, &m, &n,     &alpha,   &beta,              A.nzval,      A.rowval,   A.colptr,   convert(Ptr{Complex128}, pointer(x)),  convert(Ptr{Complex128}, pointer(y)));
   
end  # function Ac_mul_B!

#------------------------------------------------------------------------------

function Ac_mul_B!( alpha::Complex128,
                    A::SparseMatrixCSC{Complex128,Int},
                    x::Array{Complex128},
                    beta::Complex128,
                    y::Array{Complex128},
                    nthreads::Int64=0 )
# Complex:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      Base.Ac_mul_B!( alpha, A, x, beta, y )
      return
   elseif nthreads < 1
      throw(ArgumentError("nthreads < 1"))
   end

   n,m  = size(A)
   nvec = size(x,2)
   

   if size(x,1) != n || size(y,1) != m 
      throw(DimensionMismatch("length(x) != n || length(y) != m"))
   elseif size(y,2) != nvec
      throw(DimensionMismatch("length(y,2) != nvec"))
   end
   
	p  = ccall( (:ac_mul_b_cc_, spmatveclib),
		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}, Ptr{Complex128}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}),
                   &nthreads, &nvec, &m, &n,     &alpha,   &beta,              convert(Ptr{Complex128}, pointer(A.nzval)),      A.rowval,   A.colptr,   convert(Ptr{Complex128}, pointer(x)),  convert(Ptr{Complex128}, pointer(y)));
   
end  # function Ac_mul_B!


function Ac_mul_B!( alpha::Complex64,
                    A::SparseMatrixCSC{Complex64,Int32},
                    x::Array{Complex64},
                    beta::Complex64,
                    y::Array{Complex64},
                    nthreads::Int64=0 )
# Complex:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      Base.Ac_mul_B!( alpha, A, x, beta, y )
      return
   elseif nthreads < 1
      throw(ArgumentError("nthreads < 1"))
   end

   n,m  = size(A)
   nvec = size(x,2)

   if size(x,1) != n || size(y,1) != m 
      throw(DimensionMismatch("length(x) != n || length(y) != m"))
   elseif size(y,2) != nvec
      throw(DimensionMismatch("length(y,2) != nvec"))
   end
   
	p  = ccall( (:ac_mul_b_cc_short_, spmatveclib),
		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex64}, Ptr{Complex64}, Ptr{Complex64}, Ptr{Int32}, Ptr{Int32}, Ptr{Complex64}, Ptr{Complex64}),
                   &nthreads, &nvec, &m, &n,     &alpha,   &beta,              convert(Ptr{Complex64}, pointer(A.nzval)),      A.rowval,   A.colptr,   convert(Ptr{Complex64}, pointer(x)),  convert(Ptr{Complex64}, pointer(y)));
   
end  # function Ac_mul_B!

function Ac_mul_B!( alpha::Complex64,
                    A::SparseMatrixCSC{Float32,Int32},
                    x::Array{Complex64},
                    beta::Complex64,
                    y::Array{Complex64},
                    nthreads::Int64=0 )
# Complex:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      Base.Ac_mul_B!( alpha, A, x, beta, y )
      return
   elseif nthreads < 1
      throw(ArgumentError("nthreads < 1"))
   end

   n,m  = size(A)
   nvec = size(x,2)

   if size(x,1) != n || size(y,1) != m 
      throw(DimensionMismatch("length(x) != n || length(y) != m"))
   elseif size(y,2) != nvec
      throw(DimensionMismatch("length(y,2) != nvec"))
   end
   
	p  = ccall( (:ac_mul_b_rc_short_, spmatveclib),
		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex64}, Ptr{Complex64}, Ptr{Float32}, Ptr{Int32}, Ptr{Int32}, Ptr{Complex64}, Ptr{Complex64}),
                   &nthreads, &nvec, &m, &n,     &alpha,   &beta,  A.nzval,      A.rowval,   A.colptr,   convert(Ptr{Complex64}, pointer(x)),  convert(Ptr{Complex64}, pointer(y)));
   
end  # function Ac_mul_B!