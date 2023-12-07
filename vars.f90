! vars.f90

module vars
#include "macros.h"
  use fftw3
  use params

  implicit none

  ! field variables
  real(dl), allocatable :: fld(:,:,:)     ! Make this an allocatable array

  ! fft variables
  real(C_DOUBLE), pointer :: f(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)
  type(C_PTR) :: planf, planb
  
contains

  subroutine init_arrays()
    allocate(fld(nx,ny,nz))
    call allocate_3d_array(nx, ny, nz, f, fk)
    planf = fftw_plan_dft_r2c_3d(nz, ny, nx, f, Fk, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
    planb = fftw_plan_dft_c2r_3d(nz, ny, nx, Fk, f, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
  end subroutine init_arrays
  
end module vars
