
   

subroutine ac_mul_b_rr( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rr
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rr_':: ac_mul_b_rr

! y = beta*y  + alpha * A'*x

#undef CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

real(kind=8),intent(in):: alpha, beta
real(kind=8),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)
real(kind=8),intent(in):: x(m,nvec)
real(kind=8),intent(inout):: y(n,nvec)

integer ivec, i, j1,j2, j
real(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_rr

!------------------------------------------------------------------------

subroutine ac_mul_b_rc( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rc
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rc_':: ac_mul_b_rc

! y = beta*y  + alpha * A'*x

#undef CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=8),intent(in):: alpha, beta
real(kind=8),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)
complex(kind=8),intent(in):: x(m,nvec)
complex(kind=8),intent(inout):: y(n,nvec)

integer ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_rc

!------------------------------------------------------------------------

subroutine ac_mul_b_cc( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_cc
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_cc_':: ac_mul_b_cc

! y = beta*y  + alpha * A'*x

#define CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=8),intent(in):: alpha, beta
complex(kind=8),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)
complex(kind=8),intent(in):: x(m,nvec)
complex(kind=8),intent(inout):: y(n,nvec)

integer ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_cc
