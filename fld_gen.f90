! fld_gen.f90

! to do: option to read spec
! to do: declare variable for nonlinear transfer

program fld_gen
#include "macros.h"
#include "pert_macros.h"
  use params
  use vars
  use fftw3
  use lin_tran
  use nlin_tran
  use spec_init
  use corr_mod
  use io_mod
  use omp_lib
  
  implicit none

  integer :: seed
  
  ! spec and transfer variables
  integer, parameter :: kos_spec = kos, kos_tk = 2**0
  real(dl), dimension(nn*kos_tk) :: tk             ! transfer function
  real(dl), dimension(nn*kos_spec) :: spec
  
  ! Lattice initialization variables
  real(dl), parameter, dimension(SYS_DIM_BG) :: y_bg0 = (/0._dl,0._dl,11._dl,0._dl,-0.8129_dl,0._dl/) ! ICs of bg
  real(dl), parameter :: h_fac = 5._dl       ! ln(k/aH) at which to to init modes
  real(dl), parameter :: alpha_init = 0._dl  ! \alpha when \phi = phi_init
  real(dl), parameter :: phi_init = phi_p + phi_w
  real(dl), parameter :: alpha_fin = 0._dl   ! \alpha when corr is set
  real(dl), dimension(2) :: k_filter = dk*(/dble(nn-1)/2._dl, 2._dl*dble(nn-1)/3._dl/)  ! k's for the filter function

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
  case default
     stop "ERROR: SPECOPT INVALID"
  end select

  ! initialize Gaussian field
  call lat_init_spec(spec, bowler_hat_cubic, kos, nn-1, k_filter, seed, f, FK, planb, corr_norm)  
  call write_fld(n_file_g, f)  ! output Gaussian field

  ! nonlin transfer
  call quad_fld(f, quad_param)
  call write_fld(n_file_ng, f)  ! output nonlin transfered field
  
  ! lin transfer
  call read_transfer(tk, nn*kos_tk)
  call lin_transfer((/nx,ny,nz/), tk, kos_tk, int(k_filter(2)), f, fk, planf, planb)
  call write_fld(n_file_ngt, f)  ! output lin transfered field
  
  ! output settings
  call write_header(n_file_temp)
  
contains
  
  subroutine write_header(nfile)
    integer, intent(in) :: nfile

    open(unit=nfile, file=header_f)
    write(nfile, '(A, 2X, I5, 2X)') "SPECOPT: ", SPECOPT
    select case(SPECOPT)
    case(1)
       write(nfile, '(A, 2X, 30(ES22.15, 2X))') "y_bg0: ", y_bg0
       write(nfile, '(A, 2X, 1(ES22.15, 2X))') "h_fac: ", h_fac
       write(nfile, '(A, 2X, 1(ES22.15, 2X))') "alpha_init: ", alpha_init
       write(nfile, '(A, 2X, 1(ES22.15, 2X))') "alpha_fin: ", alpha_fin
       write(nfile, '(A, 2X, 1(ES22.15, 2X))') "phi_init: ", phi_init
    end select
    
    write(nfile, '(A, 2X, I5, 2X)') "POTOPT: ", POTOPT
    select case(POTOPT)
    case(5)
       
    end select
    
    write(nfile, '(2(A, 2X))') "spec_f: ", spec_f
    write(nfile, '(2(A, 2X))') "tran_f: ", tran_f
    write(nfile, '(A, 2X, I5, 2X)') "kos_spec: ", kos_spec
    write(nfile, '(A, 2X, I5, 2X)') "kos_tk: ", kos_tk
    close(unit=nfile)
  end subroutine write_header
  
end program fld_gen
