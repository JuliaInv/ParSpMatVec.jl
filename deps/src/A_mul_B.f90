
subroutine a_mul_b_rr_32( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_rr_32
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_rr_32_':: a_mul_b_rr_32

! y = beta*y  +  alpha * A*x

use omp_lib
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

real(kind=8),intent(in):: alpha, beta
real(kind=8),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)
real(kind=8),intent(in):: x(n,nvec)
real(kind=8),intent(inout):: y(m,nvec)

integer ivec, i, j1,j2, j, jaj, mythread, mm, jm
real(kind=8) xi
real(kind=8),allocatable:: yt(:)

include "A_mul_B.fi"

return
end subroutine a_mul_b_rr_32

subroutine a_mul_b_rr_64( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_rr_64
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_rr_64_':: a_mul_b_rr_64

! y = beta*y  +  alpha * A*x

use omp_lib
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

real(kind=8),intent(in):: alpha, beta
real(kind=8),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)
real(kind=8),intent(in):: x(n,nvec)
real(kind=8),intent(inout):: y(m,nvec)

integer ivec, i, j1,j2, j, jaj, mythread, mm, jm
real(kind=8) xi
real(kind=8),allocatable:: yt(:)

include "A_mul_B.fi"

return
end subroutine a_mul_b_rr_64

!--------------------------------------------------------------------

subroutine a_mul_b_rc_32( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_rc_32
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_rc_32_':: a_mul_b_rc_32

! y = beta*y  +  alpha * A*x

use omp_lib
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=8),intent(in):: alpha, beta
real(kind=8),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)
complex(kind=8),intent(in):: x(n,nvec)
complex(kind=8),intent(inout):: y(m,nvec)

integer ivec, i, j1,j2, j, jaj, mythread, mm, jm
complex(kind=8) xi
complex(kind=8),allocatable:: yt(:)

include "A_mul_B.fi"

return
end subroutine a_mul_b_rc_32

subroutine a_mul_b_rc_64( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_rc_64
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_rc_64_':: a_mul_b_rc_64

! y = beta*y  +  alpha * A*x

use omp_lib
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=8),intent(in):: alpha, beta
real(kind=8),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)
complex(kind=8),intent(in):: x(n,nvec)
complex(kind=8),intent(inout):: y(m,nvec)

integer ivec, i, j1,j2, j, jaj, mythread, mm, jm
complex(kind=8) xi
complex(kind=8),allocatable:: yt(:)

include "A_mul_B.fi"

return
end subroutine a_mul_b_rc_64

!--------------------------------------------------------------------

subroutine a_mul_b_cc_32( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_cc_32
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_cc_32_':: a_mul_b_cc_32

! y = beta*y  +  alpha * A*x

use omp_lib
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=8),intent(in):: alpha, beta
complex(kind=8),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)
complex(kind=8),intent(in):: x(n,nvec)
complex(kind=8),intent(inout):: y(m,nvec)

integer ivec, i, j1,j2, j, jaj, mythread, mm, jm
complex(kind=8) xi
complex(kind=8),allocatable:: yt(:)

include "A_mul_B.fi"

return
end subroutine a_mul_b_cc_32

subroutine a_mul_b_cc_64( nthreads, nvec, n, m, alpha, beta, A, jA, iA, x, y )
!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_cc_64
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_cc_64_':: a_mul_b_cc_64

! y = beta*y  +  alpha * A*x

use omp_lib
implicit none

integer(kind=8),intent(in):: nthreads
integer(kind=8),intent(in):: nvec
integer(kind=8),intent(in):: n  ! # of columns in A
integer(kind=8),intent(in):: m  ! # of rows in A

complex(kind=8),intent(in):: alpha, beta
complex(kind=8),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)
complex(kind=8),intent(in):: x(n,nvec)
complex(kind=8),intent(inout):: y(m,nvec)

integer ivec, i, j1,j2, j, jaj, mythread, mm, jm
complex(kind=8) xi
complex(kind=8),allocatable:: yt(:)

include "A_mul_B.fi"

return
end subroutine a_mul_b_cc_64
