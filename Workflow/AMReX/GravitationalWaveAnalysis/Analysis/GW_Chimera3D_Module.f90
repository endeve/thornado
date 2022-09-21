MODULE GW_Chimera3D_Module

     USE Constants_Module

     INTEGER            :: nx, ny, nz, nshells, ntime
     INTEGER            :: RadialZonesPerShell, NumberOfRadialShells
     INTEGER            :: k_hyperslab, k_slab_dim, nz_hyperslabs
     INTEGER            :: j_offset, k_offset, k_hyperslab_offset
     INTEGER            :: nz_hyperslab_width
     INTEGER            :: n_proc_y, n_proc_z, ngrid
     INTEGER            :: NumberOfCycles, nprint
     INTEGER            :: mpierror, comm, mpi_size, mpi_rank
     INTEGER            :: intrst, intgw  ! Number of Chimera cycles between checkpointing dumps
                                          ! and between GW data saves

     LOGICAL            :: io_initialized = .FALSE.

     INTEGER                             :: NumberOfDirs
     CHARACTER(LEN=256)                  :: model          ! Model name
     CHARACTER(LEN=256)                  :: path           ! Path to data to be analyzed
     CHARACTER(LEN=256), DIMENSION (100) :: path_list      ! List of directories 

     CHARACTER(LEN=256), DIMENSION(300)  :: FileName   ! Names of files in path directory
     CHARACTER(LEN=9), DIMENSION(300)    :: Cycles     ! Cycles in path directory

     REAL(KIND=double)  :: time, TimeOfBounce, dt, dEdt, E_total

     !-----------------------------------------------------------------------
     !            \\\\\ ARRAYS COLLECTING DATA FOR GW EXTRACTION /////
     !
     !  gw_rho_c    : density: zone average [g cm^{-3}]
     !
     !  gw_u_c      : zone centered average velocity x direction [cm s^{-1}]
     !
     !  gw_v_c      : zone centered average velocity y direction [cm s^{-1}]
     !
     !  gw_w_c      : zone centered average velocity z direction [cm s^{-1}]
     !-----------------------------------------------------------------------

     REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: gw_time
     REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)  :: gw_rho_c
     REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)  :: gw_u_c
     REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)  :: gw_v_c
     REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)  :: gw_w_c

     !-----------------------------------------------------------------------
     !            \\\\\ COORDINATES - VALUES AT CYCLE END /////
     !
     !  x_ef  : x grid zone face locations at cycle end [cm]
     !
     !  y_ef  : y grid zone face locations at cycle end
     !
     !  z_ef  : z grid zone face locations at cycle end
     !
     !  dx_cf : x_ef(i+1) - x_ef(i) [cm]
     !
     !  dy_cf : y_ef(i+1) - y_ef(i)
     !
     !  dz_cf : z_ef(i+1) - z_ef(i)
     !
     !  x_cf  : x grid zone midpoint locations at cycle end [cm]
     !
     !  x_vf  : x grid zone volume midpoint locations at cycle end [cm]
     !
     !  y_cf  : y grid zone midpoint locations at cycle end
     !
     !  z_cf  : z grid zone midpoint locations at cycle end
     !-----------------------------------------------------------------------

      REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: x_ef
      REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: y_ef
      REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: z_ef
      REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: dx_cf
      REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: dy_cf
      REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: dz_cf
      REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: x_cf
      REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: x_vf
      REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: y_cf
      REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: z_cf

     !-----------------------------------------------------------------------
     !            \\\\\ Variables used for the GW extraction /////
     !-----------------------------------------------------------------------

     !-------------------------- Gauss-Legendre ---------------------------

     REAL(KIND=double), DIMENSION(:), ALLOCATABLE        :: thg   ! theta at collocation
     REAL(KIND=double), DIMENSION(:), ALLOCATABLE        :: wth   ! theta weight 
     REAL(KIND=double), DIMENSION(:), ALLOCATABLE        :: sn    ! sin(thg)
     REAL(KIND=double), DIMENSION (:,:), ALLOCATABLE     :: Pleg  ! Legendre at theta
     REAL(KIND=double), DIMENSION(:), ALLOCATABLE        :: phg   ! phi at collocation
     REAL(KIND=double), DIMENSION(:), ALLOCATABLE        :: wph   ! phi weight

     !-------------------------- Hydro arrays -----------------------------

     REAL(KIND=double), DIMENSION(:,:), ALLOCATABLE      :: Rho2D ! 2D array needed for integration
     REAL(KIND=double), DIMENSION(:,:,:), ALLOCATABLE    :: Vel2D ! 2D array of 3-velocity

     REAL(KIND=double), DIMENSION(:), ALLOCATABLE        :: rshell ! average radius of each radial shell
     COMPLEX(KIND=double), DIMENSION(:,:), ALLOCATABLE   :: n2merid_r 
     COMPLEX(KIND=double), DIMENSION(5)                  :: A2m_loc, A2m_tot, d3Idt3_2m_loc, d3Idt3_2m_tot
     COMPLEX(KIND=double), DIMENSION(:,:), ALLOCATABLE   :: A2m_loc_r, d3Idt3_2m_loc_r
     COMPLEX(KIND=double), DIMENSION(5,2)                :: f2m  ! basis tensor f2m(m=2 to -2 ; W or X) 
     COMPLEX(KIND=double), DIMENSION(:,:,:), ALLOCATABLE :: N2m_r     ! N2m in each radial shell  

END MODULE GW_Chimera3D_Module
