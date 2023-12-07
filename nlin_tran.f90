! nlin_tran.f90

module nlin_tran
#include "macros.h"
  use fftw3
  use params

  implicit none

  real(dl), parameter :: quad_param = 1._dl

  !interface transform_nl
  !   subroutine quad_fld(f, r_param)
  !     import
  !     real(C_DOUBLE), pointer :: f(:,:,:)
  !     real(dl), intent(in) :: r_param
  !   end subroutine quad_fld
  !end interface transform_nl
  
contains

  ! to do: figure out if I should subtract mean
  subroutine quad_fld(f, r_param)
    real(C_DOUBLE), pointer :: f(:,:,:)
    real(dl), intent(in) :: r_param

    f = r_param * f**2
    !f = f - sum(f)/nvol
  end subroutine quad_fld
  
end module nlin_tran
