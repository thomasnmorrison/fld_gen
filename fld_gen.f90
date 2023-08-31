! fld_gen.f90

! to do: option to read spec
! to do: declare variable for nonlinear transfer

program fld_gen
#include "macros.h"
  use params
  use fftw3
  use lin_tran
  use nlin_tran
  use spec_init
  use corr_mod
  use omp_lib
  
  implicit none

  integer :: seed
  
  ! field variables
  real(dl), allocatable :: fld(:,:,:)     ! Make this an allocatable array

  ! spec and transfer variables
  integer, parameter :: kos_spec = kos, kos_tk = 2**0
  real(dl), dimension(nn*kos_spec) :: tk             ! transfer function
  real(dl), dimension(nn*kos_tk) :: spec
  
  ! fft variables
  real(C_DOUBLE), pointer :: f(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)
  type(C_PTR) :: planf, planb

  integer :: terror

  ! allocate arrays and plan fftw
  terror = fftw_init_threads()
  print*,"Error code is ",terror
  terror = omp_get_max_threads()
  print*,"Num Threads = ",terror
  call fftw_plan_with_nthreads(omp_get_max_threads())
  call init_arrays()
  
  ! initialize I/O
  call init_output(seed)

  ! init spec
  select case(SPECOPT)
  case(1)
     call init_corr(y_bg0, alpha_init, phi_init, alpha_fin, h_fac)  ! calculate corr
     spec = corr(:,2,2)  ! get spec for \langle |\chi|^2 \rangle
  case(default)
     print*, "Error"
     exit
  end select

  ! initialize Gaussian field
  call lat_init_spec(corr, bowler_hat_cubic, kos, nkos-1, k_filter, seed, f1, FK, planb, corr_norm)  
  call write_fld(n_file_g, f)  ! output Gaussian field

  ! nonlin transfer
  call quad_fld(f, quad_param)
  call write_fld(n_file_ng, f)  ! output nonlin transfered field
  
  ! lin transfer
  call read_transfer(tk, nn*kos_tk)
  call lin_transfer(tk, kos_tk, kcut, f, fk, planf, planb)
  call write_fld(n_file_ngt, f)  ! output lin transfered field
  
  ! output settings
  call write_header(n_file_temp)
  
contains

  subroutine init_arrays()
    allocate(fld(nx,ny,nz))
    call allocate_3d_array(nx, ny, nz, f, fk)
    planf = fftw_plan_dft_r2c_3d(nz, ny, nx, laplace, Fk, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
    planb = fftw_plan_dft_c2r_3d(nz, ny, nx, Fk, laplace, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
  end subroutine init_arrays

  subroutine write_header(nfile)
    integer, intent(in) :: nfile

    open(unit=nfile, file=header_f)
    write(nfile, ) "SPECOPT: ", SPEC_OPT
    select case(SPEC_OPT)
    case(1)
       write(nfile, ) "y_bg0: ",
       write(nfile, ) "h_fac: ",
       write(nfile, ) "alpha_init: ",
       write(nfile, ) "alpha_fin: ",
       write(nfile, ) "phi_init: ",
    end select
    
    write(nfile, ) "POTOPT: ", POTOPT
    select case(POTOPT)
    case(5)

    end select
    
    write(nfile, ) "spec_f: ", spec_f
    write(nfile, ) "tran_f: ", tran_f
    write(nfile, ) "kos_spec: ", kos_spec
    write(nfile, ) "kos_tk: ", kos_tk
    close(unit=nfile)
  end subroutine write_header
  
end program fld_gen
