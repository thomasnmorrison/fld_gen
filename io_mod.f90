! io_mod.f90

module io_mod
#include "macros.h"
  use, intrinsic :: iso_c_binding
  use params
  use vars

  implicit none

  character(len=*), parameter :: run_id = ""
  character(len=*), parameter :: spec_f = ""
  character(len=*), parameter :: tran_f = ""
  character(len=*), parameter :: header_f = ""
  
  integer, parameter :: n_file_temp = 99
  integer, parameter :: n_file_g = 98     ! Gaussian field
  integer, parameter :: n_file_ng = 97    ! nonGaussian field
  integer, parameter :: n_file_ngt = 96   ! transferred nonGaussian field
  
contains

  ! Subroutine to initialize putput files
  subroutine init_output(seed_in)
    integer, intent(in) :: seed_in
    character(len=64) :: ri_char

    write(ri_char,'(A, I3.3)') run_id, seed_in
    open(unit=n_file_g, file="gfld"//trim(ri_char)//".out", form="unformatted",access="stream")
    open(unit=n_file_ng, file="ngfld"//trim(ri_char)//".out", form="unformatted",access="stream")
    open(unit=n_file_ngt, file="ngtfld"//trim(ri_char)//".out", form="unformatted",access="stream")
  end subroutine init_output

  ! to do: check how spectra are output
  ! to do: put in read format
  subroutine read_transfer(tk, nrow)
    real(dl) :: tk(:)
    integer :: nrow
    real(dl) :: junk
    integer :: i
    
    open(unit=n_file_temp, file=tran_f, status="old")
    ! read transfer function
    do i=1, nrow
!       read(n_file_temp, '(1(ES22.15, 2X))') tk(i)
       read(n_file_temp, *) tk(i)
    end do
    close(unit=n_file_temp)
  end subroutine read_transfer

  ! Subroutine to write a field
  subroutine write_fld(nfile, f)
    integer, intent(in) :: nfile
    real(C_DOUBLE), pointer :: f(:,:,:)
    
    write(nfile) f(IRANGE)
  end subroutine write_fld
  
end module io_mod
