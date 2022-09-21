! Modifications by Jesse Buffaloe - June 2022

MODULE Meridional_Int_Module

     USE Constants_Module
     USE GW_Chimera3D_Module

     CONTAINS

!----------------------------------------------------------------------------
! 
!  'integmerid' computes the meridional integral of the calculation of N2m
!

     SUBROUTINE integmerid(nx, ny, x_cf, x_ef, dx, y_cf, y_ef, dy, &
                           RadialZonesPerShell, NumberOfRadialShells, &
                           phi, rho0, vel0, n2merid_r)

     IMPLICIT NONE

     INTEGER, INTENT(IN)              :: nx, ny, RadialZonesPerShell, NumberOfRadialShells
     REAL(KIND=double), INTENT(in)    :: phi
     REAL(KIND=double), INTENT(in)    :: dx, dy

     REAL(KIND=double), DIMENSION(nx,ny), INTENT(in)     :: rho0  ! density on a meridional slice
     REAL(KIND=double), DIMENSION(nx,ny,3), INTENT(in)   :: vel0  ! velocity on a meridional slice

     REAL(KIND=double), DIMENSION(nx), INTENT(in) :: x_cf ! x-direction cell centers and sizes
     REAL(KIND=double), DIMENSION(ny), INTENT(in) :: y_cf ! y-direction cell centers and sizes

     REAL(KIND=double), DIMENSION(nx+1), INTENT(in) :: x_ef ! x-direction cell faces
     REAL(KIND=double), DIMENSION(ny+1), INTENT(in) :: y_ef ! y-direction cell faces

     COMPLEX(KIND=double), DIMENSION(5,NumberOfRadialShells), INTENT(out) :: n2merid_r

     !-----------------------------------------------------------------------
     ! Local Variables
     !-----------------------------------------------------------------------

     INTEGER                          :: i, j, m, nshell, ier
     REAL(KIND=double)                :: Sint_real, Sint_imag

     REAL(KIND=double), DIMENSION(:), ALLOCATABLE        :: yin_real, yin_imag
     REAL(KIND=double), DIMENSION(:,:), ALLOCATABLE      :: cspl_real, cspl_imag
     COMPLEX(KIND=double), DIMENSION(:,:,:), ALLOCATABLE :: n2mintg_org  

     ! Allocate memory
     ALLOCATE(n2mintg_org(nx,ny,5))

     ! Need to allocate these here - Jesse.

     ALLOCATE (dx_cf(nx))
     ALLOCATE (dy_cf(ny))

     ! Compute geometry temporarily (do this in reader instead). - Jesse

     DO i = 1, nx

       dx_cf(i) = dx

       !PRINT*, 'i = ', i
       !PRINT*, 'x center: ', x_cf(i)
       !PRINT*, 'Left face: ', x_ef(i)

     END DO

     DO j = 1, ny

       dy_cf(j) = dy

     END DO

     !-----------------------------------------------------------------------
     !
     ! Compute the N2m integrands
     !

     !PRINT*, 'Beginning n2m loop.'

     DO j = 1, ny
       DO i = 1, nx

              !PRINT*, 'i = ', i
              !PRINT*, 'j = ', j

              !PRINT*, 'x_cf = ', x_cf(i)
              !PRINT*, 'y_cf = ', y_cf(j)

              !PRINT*, 'D = ', rho0(i,j)
             !PRINT*, 'V1_L = ', vel0(i,j,1)
             !PRINT*, 'V2_L = ', vel0(i,j,2)
             !PRINT*, 'V3_L = ', vel0(i,j,3)

              ! m = 2
              n2mintg_org(i,j,1) = x_cf(i)**3*rho0(i,j)  &
                   &  *( vel0(i,j,1)*(1.0d0-dcos(2.0d0*y_cf(j)))*dsin(y_cf(j))  &
                   &     + vel0(i,j,2)*dsin(2.0d0*y_cf(j))*dsin(y_cf(j))  &
                   &     -ui*vel0(i,j,3)*(1.0d0-dcos(2.0d0*y_cf(j))) )*exp(-2.0d0*ui*phi)

              ! m = 1
              n2mintg_org(i,j,2) = - x_cf(i)**3*rho0(i,j)  &
                   &  *( vel0(i,j,1)*dsin(2.0d0*y_cf(j))*dsin(y_cf(j))  &
                   &     +vel0(i,j,2)*dcos(2.0d0*y_cf(j))*dsin(y_cf(j))  &
                   &     -5.0d-1*ui*vel0(i,j,3)*dsin(2.0d0*y_cf(j)) )*exp(-ui*phi)

              ! m = 0
              n2mintg_org(i,j,3) = x_cf(i)**3*rho0(i,j)  &
                   &  *( vel0(i,j,1)*(3.0d0*dcos(y_cf(j))**2-1.0d0)  &
                   &     -3.0d0*vel0(i,j,2)*dcos(y_cf(j))*dsin(y_cf(j)) ) * dsin(y_cf(j))

              ! m = -1
              n2mintg_org(i,j,4) = x_cf(i)**3*rho0(i,j)  &
                   &  *( vel0(i,j,1)*dsin(2.0d0*y_cf(j))*dsin(y_cf(j))  &
                   &     +vel0(i,j,2)*dcos(2.0d0*y_cf(j))*dsin(y_cf(j))  &
                   &     +5.0d-1*ui*vel0(i,j,3)*dsin(2.0d0*y_cf(j)) )*exp(ui*phi)

              ! m = -2
              n2mintg_org(i,j,5) = x_cf(i)**3*rho0(i,j)  &
                   &  *( vel0(i,j,1)*(1.0d0-dcos(2.0d0*y_cf(j)))*dsin(y_cf(j))  &
                   &     + vel0(i,j,2)*dsin(2.0d0*y_cf(j))*dsin(y_cf(j))  &
                   &     +ui*vel0(i,j,3)*(1.0d0-dcos(2.0d0*y_cf(j))) )*exp(2.0d0*ui*phi)

             !PRINT*, 'n2mintg  2: ', n2mintg_org(i,j,1)
             !PRINT*, 'n2mintg  1: ', n2mintg_org(i,j,2)
             !PRINT*, 'n2mintg  0: ', n2mintg_org(i,j,3)
             !PRINT*, 'n2mintg -1: ', n2mintg_org(i,j,4)
             !PRINT*, 'n2mintg -2: ', n2mintg_org(i,j,5)

       END DO
     END DO

     !-----------------------------------------------------------------------
     !
     ! Compute the quadrature with a spline (splq)
     !

     ALLOCATE(cspl_real(ny,3), cspl_imag(ny,3))
     ALLOCATE(yin_real(ny), yin_imag(ny))

     n2merid_r = zero

     DO m = 1, 5

       DO nshell = 1, NumberOfRadialShells

         DO i = (nshell-1)*RadialZonesPerShell+1, nshell*RadialZonesPerShell

           ! Real part

           yin_real(:) = REAL(n2mintg_org(i,:,m))

           CALL splc(y_cf,ny,yin_real,df,iopt,cspl_real,ny,ier)
           CALL splq(y_cf,ny,yin_real,cspl_real,ny,y_ef(1),y_ef(ny),Sint_real,ier)

           ! Imaginary part

           yin_imag(:) = AIMAG(n2mintg_org(i,:,m))

           CALL splc(y_cf,ny,yin_imag,df,iopt,cspl_imag,ny,ier)
           CALL splq(y_cf,ny,yin_imag,cspl_imag,ny,y_ef(1),y_ef(ny),Sint_imag,ier)

           n2merid_r(m,nshell) = n2merid_r(m,nshell) + (Sint_real + ui * Sint_imag) * dx_cf(i)

         END DO ! i

       END DO ! nshell

     END DO

     !-----------------------------------------------------------------------
     !
     ! Add the dimensional factors
     !

     DO m = 1, 5
       n2merid_r(m,:) = Nc(m) * n2merid_r(m,:)
     END DO

     ! Clean up

     DEALLOCATE(dx_cf, dy_cf, n2mintg_org, cspl_real, cspl_imag, yin_real, yin_imag)

     RETURN
     END SUBROUTINE integmerid

!----------------------------------------------------------------------------

END MODULE Meridional_Int_Module
