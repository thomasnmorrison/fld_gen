! spec_init.f90

module spec_init
#include "macros.h"
  use fftw3
  use params
  use grv_mod

  implicit none

  interface
     function filt_template(k_in, param_r) result(f)
       import
       real(dl) :: k_in
       real(dl) :: param_r(:)
       real(dl) :: f
     end function filt_template
  end interface
  
contains

  ! Subroutine to initialize a single field from a supplied spectrum
  subroutine lat_init_spec(spec, filt, kos, kcut, kfilt, seed_in, f, fk, planb, norm)
    real(dl), intent(in) :: spec(:)      ! spectrum to initialize
    procedure(filt_template) :: filt     ! filter function
    integer, intent(in) :: kos           ! oversamping ratio of k in corr relative to dk on the lattice
    integer, intent(in) :: kcut          ! maximum k index of corr, zero all above
    real(dl), intent(in) :: kfilt(:)     ! filter parameters
    integer, intent(in) :: seed_in       ! seed for rng
    real(C_DOUBLE), pointer :: f(:,:,:)              ! pointer field
    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)  ! pointer for fft
    type(C_PTR) :: planb                             ! fft c2r plan
    real(dl), intent(in) :: norm

    integer :: i, j, k, n
    integer :: ii, jj, kk, l
    real(dl) :: rad
    real(dl) :: spec_interp
    complex(C_DOUBLE_COMPLEX), dimension(1) :: grv

    call init_rng(seed_in)  ! initialize rng

    do k=1,nz; if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
       do j=1,ny; if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif 
          do i=1,nnx; ii=i-1
             rad = sqrt(dble(ii**2)+dble(jj**2)+dble(kk**2))  ! find radius
             l = floor(rad*kos)  ! get radius index
             if (l == 0 .or. l .gt. kcut*kos) then
                fk(LATIND) = (0._dl, 0._dl)    ! zero the zero/cut mode and cycle the loop
                cycle
             else
                spec_interp = (1._dl + l - rad*kos)*spec(l) + (rad*kos - l)*spec(l+1)   ! interpolate spec
             endif
             grv(:) = get_grv_complex(1)
             fk(LATIND) = spec_interp*grv(1)*filt(dk*rad, kfilt)  ! initialize mode
          end do
       end do
    end do

    ! Loops to set the ii=0 components to be complex conjugate
    do k=1,nz; if (k.ne.1) then; kk = nz+2-k; else; kk=k; endif
       do j=2,nny; jj = ny+2-j
          fk(1,jj,kk) = conjg(fk(1,j,k))
       end do
    end do

    do k=2,nnz; kk = nz+2-k
       fk(1,1,kk) = conjg(fk(1,1,k))
    end do

    fk = fk*norm/nvol  ! apply normalization
    call fftw_execute_dft_c2r(planb, fk, f)
  end subroutine lat_init_spec

    ! Function for a top_hat filter
  function top_hat(k_in, param_r) result(f)
    real(dl) :: k_in
    real(dl) :: param_r(:)
    real(dl) :: f

    f = 0.5_dl + sign(0.5_dl, param_r(1)-k_in)
  end function top_hat

  ! Function for a bowler hat filter with a cubic tail
  ! n.b. must have param_r(1) .lt. param_r(2)
  function bowler_hat_cubic(k_in, param_r) result(f)
    real(dl) :: k_in
    real(dl) :: param_r(:)
    real(dl) :: f

    real(dl) :: temp
    
    temp = (2._dl*k_in - (param_r(1) + param_r(2))) / ((param_r(2) - param_r(1)))
    f = 0.5_dl + sign(0.5_dl, param_r(1)-k_in) &
         + (sign(0.5_dl, k_in-param_r(1)) + sign(0.5_dl, param_r(2)-k_in)) &
         * (0.5_dl + 0.25_dl*(temp**3 - 3._dl*temp))
  end function bowler_hat_cubic
  
end module spec_init
