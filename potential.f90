! potential.f90

! Module for potential functions

! to do: move the $\Delta V$ stuff to another module with conditional compilation
! to do: add a macro for if $\Delta V$ is added
! to do: add a macro for nfld
! to do: speed test select, if, sign

module potential
#include "macros.h"
  use params

  implicit none

  real(dl), parameter, dimension(nfld) :: m2 = (/1._dl,1._dl/)

  real(dl), parameter :: phi_p = 8.49953_dl
  real(dl), parameter :: phi_w = 0.12547_dl
  real(dl), parameter :: lambda_chi = 2.57010e5
  real(dl), parameter :: vev = 0.07_dl
  
contains

  function bg_v(f) result(v)
    real(dl), intent(in) :: f(:)
    real(dl) :: v

    v = 0.5_dl*sum(m2*f**2) + bg_delta_v(f)
  end function bg_v

  function bg_dv(f,ind) result(dv)
    real(dl), intent(in) :: f(:)
    integer, intent(in) :: ind
    real(dl) :: dv

    dv = m2(ind)*f(ind) + bg_delta_dv(f,ind)
  end function bg_dv

  function bg_ddv(f,ind1,ind2) result(ddv)
    real(dl), intent(in) :: f(:)
    integer, intent(in) :: ind1, ind2
    real(dl) :: ddv

    if (ind1 == ind2) then
       ddv = m2(ind1) + bg_delta_ddv(f,ind1,ind2)
    else
       ddv = bg_delta_ddv(f,ind1,ind2)
    end if
  end function bg_ddv

  function bg_ddv_mat(f) result(ddv)
    real(dl), intent(in) :: f(:)
    real(dl), dimension(nfld,nfld) :: ddv

    integer :: i,j
    do j=1,nfld
       do i=1,nfld
          ddv(i,j) = bg_ddv(f,i,j)
       end do
    end do
  end function bg_ddv_mat
    
  function bg_delta_v(f) result(v)
    real(dl), intent(in) :: f(:)
    real(dl) :: v

    v = (sign(0.5_dl,(phi_p+phi_w-f(1))) + sign(0.5_dl,(f(1)-(phi_p-phi_w)))) &
         * 0.25_dl*lambda_chi*(((f(1)-phi_p)/phi_w)**2 - 1._dl)**2 &
         * ((f(2)**2 - vev**2)**2 - vev**4)
  end function bg_delta_v

  function bg_delta_dv(f,ind) result(dv)
    real(dl), intent(in) :: f(:)
    integer, intent(in) :: ind
    real(dl) :: dv

    select case(ind)
    case(1)
       dv = (sign(0.5_dl,(phi_p+phi_w-f(1))) + sign(0.5_dl,(f(1)-(phi_p-phi_w)))) &
            * lambda_chi*((f(1)-phi_p)/phi_w)*(((f(1)-phi_p)/phi_w)**2 - 1._dl)/phi_w &
            * ((f(2)**2 - vev**2)**2 - vev**4)
    case(2)
       dv = + (sign(0.5_dl,(phi_p+phi_w-f(1))) + sign(0.5_dl,(f(1)-(phi_p-phi_w)))) &
            * lambda_chi*(((f(1)-phi_p)/phi_w)**2 - 1._dl)**2 &
            * f(2)*(f(2)**2 - vev**2)
    end select
   end function bg_delta_dv

   function bg_delta_ddv(f,ind1,ind2) result(ddv)
    real(dl), intent(in) :: f(:)
    integer, intent(in) :: ind1, ind2
    real(dl) :: ddv

    integer :: i,j

    i = min(ind1,ind2); j = max(ind1,ind2)
    if (i == 1 .and. j == 1) then
       ddv = (sign(0.5_dl,(phi_p+phi_w-f(1))) + sign(0.5_dl,(f(1)-(phi_p-phi_w)))) &
            * lambda_chi*((3._dl*(f(1)-phi_p)/phi_w)**2 - 1._dl)/phi_w**2 &
            * ((f(2)**2 - vev**2)**2 - vev**4)
    else if (i == 1 .and. j == 2) then
       ddv = (sign(0.5_dl,(phi_p+phi_w-f(1))) + sign(0.5_dl,(f(1)-(phi_p-phi_w)))) &
            * lambda_chi*((f(1)-phi_p)/phi_w)*(((f(1)-phi_p)/phi_w)**2 - 1._dl)/phi_w &
            * f(2)*(f(2)**2 - vev**2)
    else !if (i == 2 .and. j == 2) then
       ddv = (sign(0.5_dl,(phi_p+phi_w-f(1))) + sign(0.5_dl,(f(1)-(phi_p-phi_w)))) &
            * lambda_chi*(((f(1)-phi_p)/phi_w)**2 - 1._dl)**2 &
            * (3._dl*f(2)**2 - vev**2)
    end if
  end function bg_delta_ddv
  
end module potential
