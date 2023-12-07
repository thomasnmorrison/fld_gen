! lin_tran.f90

module lin_tran
#include "macros.h"
  use fftw3 
  use params

  implicit none

contains

  ! Subroutine to perform a linear transformation on a field
  subroutine lin_transfer(nsize, tk, kos, kcut, f, fk, planf, planb)
    integer, dimension(1:3), intent(in) :: nsize
    real(dl), intent(in) :: tk(:)                 ! transfer function values at oversampled k
    integer, intent(in) :: kos                    ! over sample ratio of k for tk
    integer, intent(in) :: kcut                   ! 
    real(C_DOUBLE), pointer :: f(:,:,:)           ! field pointer
    complex(C_DOUBLE_COMPLEX), pointer :: fk(:,:,:)  ! fft pointer
    type(C_PTR) :: planb, planf                   ! fftw plans

    integer :: n1, n2, n3
     integer :: nn1, nn2, nn3
    integer :: i,j,k,ii,jj,kk,l
    real(dl) :: rad, t_interp

    n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
    nn1 = n1/2+1; nn2=n2/2+1; nn3=n3/2+1

    call fftw_execute_dft_r2c(planf, f, fk)
    
    ! loop over wave numbers
!$OMP PARALLEL DO PRIVATE(ii,jj,kk,l,t_interp,rad) FIRSTPRIVATE(n1,n2,n3,nn3,nn2,nn1,kos)
    do k=1,nz; if (k>nn3) then; kk = k - (n3+1); else; kk=k-1; endif
       do j=1,ny; if (j>nn2) then; jj = j - (n2+1); else; jj=j-1; endif 
          do i=1,nn1; ii=i-1
             rad = sqrt(dble(ii**2)+dble(jj**2)+dble(kk**2))  ! find radius
             l = floor(rad*kos)                               ! get radius index
             if (l == 0 .or. rad .gt. dble(kcut)) then        ! zero the zero/cut mode
                t_interp = 0._dl
             else
                t_interp = (1._dl + l - rad*kos)*tk(l) + (rad*kos - l)*tk(l+1)
             end if
             fk(LATIND) = t_interp*fk(LATIND)
          end do
       end do
    end do
!$OMP END PARALLEL DO

    call fftw_execute_dft_c2r(planb, fk, f)
    f = f/nvol
  end subroutine lin_transfer
  
end module lin_tran
