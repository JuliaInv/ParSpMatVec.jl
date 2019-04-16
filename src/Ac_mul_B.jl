
export Ac_mul_B!

function Ac_mul_B!( alpha::Float64,
                    A::SparseMatrixCSC{Float64,Ti},
                    x::Array{Float64},
                    beta::Float64,
                    y::Array{Float64},
                    nthreads::Int64=0 ) where Ti
# Real:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      #Ac_mul_B!( alpha, A, x, beta, y )
	  mul!(y,adjoint(A),x, alpha, beta)
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
   
	if Ti == Int32
      p = ccall( (:ac_mul_b_rr_32_, spmatveclib),
      Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}),
              nthreads, nvec, m, n,    alpha,   beta,              A.nzval,      A.rowval,   A.colptr,   x,   y)
   elseif Ti == Int64
      p = ccall( (:ac_mul_b_rr_64_, spmatveclib),
      Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}),
              nthreads, nvec, m, n,    alpha,   beta,              A.nzval,      A.rowval,   A.colptr,   x,   y)
   else
      error("Unsupported sparse matrix indexing integer type $Ti")
   end
   return p
end  # function Ac_mul_B!

#------------------------------------------------------------------------------

function Ac_mul_B!( alpha::ComplexF64,
                    A::SparseMatrixCSC{Float64,Ti},
                    x::Array{ComplexF64},
                    beta::ComplexF64,
                    y::Array{ComplexF64},
                    nthreads::Int64=0 ) where Ti
# Real, Complex A:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      mul!(y,adjoint(A),x, alpha, beta)
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
   
	if Ti == Int32
	   p = ccall( (:ac_mul_b_rc_32_, spmatveclib),
		 Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{ComplexF64}, Ref{ComplexF64}),
                   nthreads, nvec, m, n,     alpha,   beta,          A.nzval,      A.rowval,   A.colptr,   x,  y)
   elseif Ti == Int64
      p = ccall( (:ac_mul_b_rc_64_, spmatveclib),
		 Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}),
                   nthreads, nvec, m, n,     alpha,   beta,          A.nzval,      A.rowval,   A.colptr,   x,  y)
   else
      error("Unsupported sparse matrix indexing integer type $Ti")
   end
   return p
end  # function Ac_mul_B!

#------------------------------------------------------------------------------

function Ac_mul_B!( alpha::ComplexF64,
                    A::SparseMatrixCSC{ComplexF64,Ti},
                    x::Array{ComplexF64},
                    beta::ComplexF64,
                    y::Array{ComplexF64},
                    nthreads::Int64=0 ) where Ti
# Complex:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      mul!(y,adjoint(A),x, alpha, beta) #Base.Ac_mul_B!( alpha, A, x, beta, y )
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
   
	if Ti == Int32
      p = ccall( (:ac_mul_b_cc_32_, spmatveclib),
		 Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{Int32}, Ref{Int32}, Ref{ComplexF64}, Ref{ComplexF64}),
                   nthreads, nvec, m, n,     alpha,   beta,           A.nzval,      A.rowval,   A.colptr,   x,  y)
   elseif Ti == Int64
      p = ccall( (:ac_mul_b_cc_64_, spmatveclib),
		 Nothing, ( Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64}),
                   nthreads, nvec, m, n,     alpha,   beta,           A.nzval,      A.rowval,   A.colptr,   x,  y)
   else
      error("Unsupported sparse matrix indexing integer type $Ti")
   end
   return p
end  # function Ac_mul_B!


function Ac_mul_B!( alpha::ComplexF32,
                    A::SparseMatrixCSC{ComplexF32,Ti},
                    x::Array{ComplexF32},
                    beta::ComplexF32,
                    y::Array{ComplexF32},
                    nthreads::Int64=0 ) where Ti
# Complex:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      mul!(y,adjoint(A),x, alpha, beta) #Base.Ac_mul_B!( alpha, A, x, beta, y )
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
   
   if Ti == Int32
	   p  = ccall( (:ac_mul_b_cc_short_32_, spmatveclib),
		    Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{Int32}, Ptr{Int32}, Ptr{ComplexF32}, Ptr{ComplexF32}),
                   Ref(nthreads), Ref(nvec), Ref(m), Ref(n),     Ref(alpha),   Ref(beta),              convert(Ptr{ComplexF32}, pointer(A.nzval)),      A.rowval,   A.colptr,   convert(Ptr{ComplexF32}, pointer(x)),  convert(Ptr{ComplexF32}, pointer(y)));
   elseif Ti == Int64
      p  = ccall( (:ac_mul_b_cc_short_64_, spmatveclib),
		    Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF32}, Ptr{ComplexF32}),
                   Ref(nthreads), Ref(nvec), Ref(m), Ref(n),     Ref(alpha),   Ref(beta),              convert(Ptr{ComplexF32}, pointer(A.nzval)),      A.rowval,   A.colptr,   convert(Ptr{ComplexF32}, pointer(x)),  convert(Ptr{ComplexF32}, pointer(y)));
   else
      error("Unsupported sparse matrix indexing integer type $Ti")
   end
   return p 
end  # function Ac_mul_B!

function Ac_mul_B!( alpha::ComplexF32,
                    A::SparseMatrixCSC{Float32,Ti},
                    x::Array{ComplexF32},
                    beta::ComplexF32,
                    y::Array{ComplexF32},
                    nthreads::Int64=0 ) where Ti
# Complex:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      mul!(y,adjoint(A),x, alpha, beta) #Base.Ac_mul_B!( alpha, A, x, beta, y )
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
   
   if Ti == Int32
	   p  = ccall( (:ac_mul_b_rc_short_32_, spmatveclib),
		    Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{Float32}, Ptr{Int32}, Ptr{Int32}, Ptr{ComplexF32}, Ptr{ComplexF32}),
                   Ref(nthreads), Ref(nvec), Ref(m), Ref(n),     Ref(alpha),   Ref(beta),  A.nzval,      A.rowval,   A.colptr,   convert(Ptr{ComplexF32}, pointer(x)),  convert(Ptr{ComplexF32}, pointer(y)));
   elseif Ti == Int64
      p  = ccall( (:ac_mul_b_rc_short_64_, spmatveclib),
		    Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{Float32}, Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF32}, Ptr{ComplexF32}),
                   Ref(nthreads), Ref(nvec), Ref(m), Ref(n),     Ref(alpha),   Ref(beta),  A.nzval,      A.rowval,   A.colptr,   convert(Ptr{ComplexF32}, pointer(x)),  convert(Ptr{ComplexF32}, pointer(y)));
   else
      error("Unsupported sparse matrix indexing integer type $Ti")
   end
   return p
end  # function Ac_mul_B!


function Ac_mul_B!( alpha::ComplexF64,
                    A::SparseMatrixCSC{ComplexF32,Ti},
                    x::Array{ComplexF64},
                    beta::ComplexF64,
                    y::Array{ComplexF64},
                    nthreads::Int64=0 ) where Ti
# Complex:  y = beta*y  +  alpha * A'*x 

   if nthreads == 0
      mul!(y,adjoint(A),x, alpha, beta) #Base.Ac_mul_B!( alpha, A, x, beta, y )
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
   
   if Ti == Int32
	   p  = ccall( (:ac_mul_b_cc_mixed_32_, spmatveclib),
		    Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF32}, Ptr{Int32}, Ptr{Int32}, Ptr{ComplexF64}, Ptr{ComplexF64}),
                   Ref(nthreads), Ref(nvec), Ref(m), Ref(n),     Ref(alpha),   Ref(beta),              convert(Ptr{ComplexF32}, pointer(A.nzval)),      A.rowval,   A.colptr,   convert(Ptr{ComplexF64}, pointer(x)),  convert(Ptr{ComplexF64}, pointer(y)));
   elseif Ti == Int64
      p  = ccall( (:ac_mul_b_cc_mixed_64_, spmatveclib),
		    Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF32}, Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF64}, Ptr{ComplexF64}),
                   Ref(nthreads), Ref(nvec), Ref(m), Ref(n),     Ref(alpha),   Ref(beta),              convert(Ptr{ComplexF32}, pointer(A.nzval)),      A.rowval,   A.colptr,   convert(Ptr{ComplexF64}, pointer(x)),  convert(Ptr{ComplexF64}, pointer(y)));
   else
      error("Unsupported sparse matrix indexing integer type $Ti")
   end
   return p
end  # function Ac_mul_B!

