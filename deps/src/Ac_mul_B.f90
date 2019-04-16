
   

subroutine ac_mul_b_rr_32( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rr_32
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rr_32_':: ac_mul_b_rr_32

! y = beta*y  + alpha * A'*x

#undef CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

real(kind=8),intent(in):: alpha, beta
real(kind=8),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)
real(kind=8),intent(in):: x(m,nvec)
real(kind=8),intent(inout):: y(n,nvec)

integer ivec, i, j1,j2, j
real(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_rr_32

subroutine ac_mul_b_rr_64( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rr_64
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rr_':: ac_mul_b_rr_64

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

integer(kind=8) ivec, i, j1,j2, j
real(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_rr_64

!------------------------------------------------------------------------

subroutine ac_mul_b_rc_32( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rc_32
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rc_32_':: ac_mul_b_rc_32

! y = beta*y  + alpha * A'*x

#undef CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=8),intent(in):: alpha, beta
real(kind=8),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)
complex(kind=8),intent(in):: x(m,nvec)
complex(kind=8),intent(inout):: y(n,nvec)

integer ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_rc_32

subroutine ac_mul_b_rc_64( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rc_64
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rc_64_':: ac_mul_b_rc_64

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

integer(kind=8) ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_rc_64

!------------------------------------------------------------------------

subroutine ac_mul_b_cc_32( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_cc_32
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_cc_32_':: ac_mul_b_cc_32

! y = beta*y  + alpha * A'*x

#define CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=8),intent(in):: alpha, beta
complex(kind=8),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)
complex(kind=8),intent(in):: x(m,nvec)
complex(kind=8),intent(inout):: y(n,nvec)

integer ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_cc_32

subroutine ac_mul_b_cc_64( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_cc_64
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_cc_64_':: ac_mul_b_cc_64

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

integer(kind=8) ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_cc_64

subroutine ac_mul_b_cc_short_64( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_cc_short_64
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_cc_short_64_':: ac_mul_b_cc_short_64

! y = beta*y  + alpha * A'*x

#define CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=4),intent(in):: alpha, beta
complex(kind=4),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)
complex(kind=4),intent(in):: x(m,nvec)
complex(kind=4),intent(inout):: y(n,nvec)



integer(kind=8) ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_cc_short_64

subroutine ac_mul_b_cc_short_32( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_cc_short_32
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_cc_short_32_':: ac_mul_b_cc_short_32

! y = beta*y  + alpha * A'*x

#define CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=4),intent(in):: alpha, beta
complex(kind=4),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)
complex(kind=4),intent(in):: x(m,nvec)
complex(kind=4),intent(inout):: y(n,nvec)



integer(kind=8) ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_cc_short_32

subroutine ac_mul_b_rc_short_64( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rc_short_64
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rc_short_64_':: ac_mul_b_rc_short_64

! y = beta*y  + alpha * A'*x

#undef CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=4),intent(in):: alpha, beta
real(kind=4),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)
complex(kind=4),intent(in):: x(m,nvec)
complex(kind=4),intent(inout):: y(n,nvec)



integer(kind=8) ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_rc_short_64

subroutine ac_mul_b_rc_short_32( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rc_short_32
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rc_short_32_':: ac_mul_b_rc_short_32

! y = beta*y  + alpha * A'*x

#undef CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=4),intent(in):: alpha, beta
real(kind=4),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)
complex(kind=4),intent(in):: x(m,nvec)
complex(kind=4),intent(inout):: y(n,nvec)



integer(kind=8) ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_rc_short_32

subroutine ac_mul_b_cc_mixed_64( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_cc_mixed_64
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_cc_mixed_64_':: ac_mul_b_cc_mixed_64

! y = beta*y  + alpha * A'*x

#define CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=8),intent(in):: alpha
complex(kind=8),intent(in):: beta
complex(kind=4),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)
complex(kind=8),intent(in):: x(m,nvec)
complex(kind=8),intent(inout):: y(n,nvec)



integer(kind=8) ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_cc_mixed_64

subroutine ac_mul_b_cc_mixed_32( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_cc_mixed_32
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_cc_mixed_32_':: ac_mul_b_cc_mixed_32

! y = beta*y  + alpha * A'*x

#define CMPLXA
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec ! # of vectors
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=8),intent(in):: alpha
complex(kind=8),intent(in):: beta
complex(kind=4),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)
complex(kind=8),intent(in):: x(m,nvec)
complex(kind=8),intent(inout):: y(n,nvec)



integer(kind=8) ivec, i, j1,j2, j
complex(kind=8) t

#include "Ac_mul_B.fi"

return
end subroutine ac_mul_b_cc_mixed_32