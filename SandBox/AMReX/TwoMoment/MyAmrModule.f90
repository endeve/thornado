MODULE MyAmrModule

  ! --- AMReX Modules ---
  USE amrex_fort_module,      ONLY: &
    AR => amrex_real, &
    amrex_spacedim
  USE amrex_base_module,      ONLY: &
    amrex_init, &
    amrex_initialized, &
    amrex_parallel_ioprocessor
  USE amrex_amr_module,       ONLY: &
    amrex_amrcore_init, &
    amrex_amrcore_initialized, &
    amrex_is_all_periodic, &
    amrex_spacedim
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_boxarray_module,  ONLY: &
    amrex_boxarray
  USE amrex_distromap_module, ONLY: &
    amrex_distromap
  USE amrex_geometry_module,  ONLY: &
    amrex_geometry

  ! --- thornado Modules ---
  USE ProgramHeaderModule,   ONLY: &
    InitializeProgramHeader, nDimsX
  ! --- Local Modules ---
  USE MyAmrDataModule, ONLY: &
    InitializeDataAMReX, &
    FinalizeDataAMReX


  ! --- thornado ---
  REAL(AR)                       :: t_end, t_wrt, dt_wrt, t_chk, dt_chk
  REAL(AR),          ALLOCATABLE :: t(:), dt(:)
  REAL(AR),          ALLOCATABLE :: V_0(:)
  REAL(AR)                       :: CFL
  INTEGER                        :: nNodes, nStages
  INTEGER                        :: nE, swE, bcE
  INTEGER                        :: iCycleD, iCycleW, iCycleChk, iRestart, nSpecies
  INTEGER,           ALLOCATABLE :: nX(:), swX(:), bcX(:)
  REAL(AR),          ALLOCATABLE :: xL(:), xR(:)
  REAL(AR)                       :: eL, eR, zoomE
  REAL(AR)                       :: D_0, Chi, Sigma
  CHARACTER(LEN=:),  ALLOCATABLE :: ProgramName
  CHARACTER(LEN=:),  ALLOCATABLE :: Scheme
  CHARACTER(LEN=32), SAVE        :: CoordSys
  LOGICAL,           SAVE        :: UsePhysicalUnits
  LOGICAL,           SAVE        :: DEBUG



  ! --- AMReX  ---
  INTEGER                                    :: MaxLevel, nLevels, coord_sys
  INTEGER                                    :: MaxGridSizeX1
  INTEGER                                    :: MaxGridSizeX2
  INTEGER                                    :: MaxGridSizeX3
  INTEGER                                    :: BlockingFactorX1
  INTEGER                                    :: BlockingFactorX2
  INTEGER                                    :: BlockingFactorX3
  INTEGER                                    :: MaxGridSizeX(3)
  INTEGER,               ALLOCATABLE, SAVE   :: StepNo(:)
  TYPE(amrex_boxarray),  ALLOCATABLE, PUBLIC :: BA(:)
  TYPE(amrex_distromap), ALLOCATABLE, PUBLIC :: DM(:)
  TYPE(amrex_geometry),  ALLOCATABLE, PUBLIC :: GEOM(:)

  ! --- Equation Of State ---
  REAL(AR)                      :: Gamma_IDEAL
  CHARACTER(LEN=:), ALLOCATABLE :: EquationOfState
  CHARACTER(LEN=:), ALLOCATABLE :: EosTableName

  ! --- Positivity limiter ---
  LOGICAL  :: UsePositivityLimiter
  REAL(AR) :: Min_1, Min_2

CONTAINS

  SUBROUTINE MyAmrInit

    REAL(AR), PARAMETER :: Zero = 0.0_AR
    REAL(AR), PARAMETER :: One  = 1.0_AR
    REAL(AR), PARAMETER :: Two  = 2.0_AR

    TYPE(amrex_parmparse) :: PP

    IF( .NOT. amrex_initialized() ) &
      CALL amrex_init()

    IF( .NOT. amrex_amrcore_initialized() ) &
      CALL amrex_amrcore_init()

    DEBUG = .FALSE.
    CALL amrex_parmparse_build( PP )
      CALL PP % query( 'DEBUG', DEBUG )
    CALL amrex_parmparse_destroy( PP )

    UsePhysicalUnits = .FALSE.
    ! --- thornado paramaters thornado.* ---
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get   ( 'dt_wrt',           dt_wrt )
      CALL PP % get   ( 'dt_chk',           dt_chk )
      CALL PP % get   ( 't_end',            t_end )
      CALL PP % get   ( 'nNodes',           nNodes )
      CALL PP % get   ( 'CFL',              CFL )
      CALL PP % get   ( 'ProgramName',      ProgramName )
      CALL PP % get   ( 'Scheme',           Scheme )
      CALL PP % getarr( 'bcX',              bcX )
      CALL PP % getarr( 'swX',              swX ) 
      CALL PP % getarr( 'V_0',              V_0 )
      CALL PP % query ( 'nE', nE )
      CALL PP % get   ( 'swE',              swE )
      CALL PP % get   ( 'bcE',              bcE )
      CALL PP % get   ( 'eL',  eL )
      CALL PP % get   ( 'eR',  eR )  
      CALL PP % get   ( 'D_0',  D_0 )
      CALL PP % get   ( 'Chi',  Chi )
      CALL PP % get   ( 'Sigma',  Sigma )
      CALL PP % get   ( 'zoomE',  zoomE )
      CALL PP % get   ( 'nSpecies',        nSpecies )
      CALL PP % get   ( 'iCycleD',          iCycleD )
      CALL PP % get   ( 'iCycleW',          iCycleW )
      CALL PP % get   ( 'iCycleChk',        iCycleChk )
      CALL PP % get   ( 'iRestart',         iRestart )
      CALL PP % query ( 'UsePhysicalUnits', UsePhysicalUnits )
    CALL amrex_parmparse_destroy( PP )
          


    ! --- Parameters geometry.* ---
    CALL amrex_parmparse_build( PP, 'geometry' )
      CALL PP % get   ( 'coord_sys',  coord_sys )
      CALL PP % getarr( 'prob_lo',    xL )
      CALL PP % getarr( 'prob_hi',    xR )
    CALL amrex_parmparse_destroy( PP )
    IF     ( coord_sys .EQ. 0 )THEN
      CoordSys = 'CARTESIAN'
    ELSE IF( coord_sys .EQ. 1 )THEN
      CoordSys = 'CYLINDRICAL'
    ELSE IF( coord_sys .EQ. 2 )THEN
      CoordSys = 'SPHERICAL'
    ELSE
      STOP 'Invalid choice for coord_sys'
    END IF


    ! --- Parameters amr.* ---
    MaxGridSizeX1    = 1
    MaxGridSizeX2    = 1
    MaxGridSizeX3    = 1
    BlockingFactorX1 = 1
    BlockingFactorX2 = 1
    BlockingFactorX3 = 1
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr( 'n_cell',            nX )
      CALL PP % query ( 'max_grid_size_x',   MaxGridSizeX1 )
      CALL PP % query ( 'max_grid_size_y',   MaxGridSizeX2 )
      CALL PP % query ( 'max_grid_size_z',   MaxGridSizeX3 )
      CALL PP % query ( 'blocking_factor_x', BlockingFactorX1 )
      CALL PP % query ( 'blocking_factor_y', BlockingFactorX2 )
      CALL PP % query ( 'blocking_factor_z', BlockingFactorX3 )
      CALL PP % get   ( 'max_level',         MaxLevel )
    CALL amrex_parmparse_destroy( PP )

    ! --- Equation of state parameters EoS.* ---
    Gamma_IDEAL     = 5.0_AR / 3.0_AR
    EquationOfState = 'IDEAL'
    EosTableName    = ''
    CALL amrex_parmparse_build( PP, 'EoS' )
      CALL PP % query( 'Gamma',           Gamma_IDEAL )
      CALL PP % query( 'EquationOfState', EquationOfState )
      CALL PP % query( 'EosTableName',    EosTableName    )
    CALL amrex_parmparse_destroy( PP )

    ! --- Positivitiy limiter parameters PL.* ---
    UsePositivityLimiter = .TRUE.
    Min_1                = 1.0e-12_AR
    Min_2                = 1.0e-12_AR
    CALL amrex_parmparse_build( PP, 'PL' )
      CALL PP % query( 'UsePositivityLimiter', UsePositivityLimiter )
      CALL PP % query( 'Min_1'               , Min_1                )
      CALL PP % query( 'Min_2'               , Min_2                )
    CALL amrex_parmparse_destroy( PP )

    MaxGridSizeX = [ MaxGridSizeX1, MaxGridSizeX2, MaxGridSizeX3 ]

    nLevels = MaxLevel + 1

    CALL InitializeProgramHeader &
           ( ProgramName_Option = TRIM( ProgramName ), &
             nNodes_Option = nNodes, nX_Option = nX, swX_Option = swX, bcX_Option = bcX, &
             xL_Option = xL, xR_Option = xR, &
             nE_Option = nE, swE_Option = swE, bcE_Option = bcE, eL_option = eL, eR_option = eR, &
             Verbose_Option = amrex_parallel_ioprocessor() )

    IF( nDimsX .NE. amrex_spacedim )THEN
      WRITE(*,'(A)') 'ERROR'
      WRITE(*,'(A)') '-----'
      WRITE(*,'(A)') 'thornado nDimsX different from AMReX amrex_spacedim.'
      WRITE(*,'(A)') 'Check DIM parameter in GNUmakefile. Stopping...'
      STOP
    END IF

    ALLOCATE( StepNo(0:nLevels-1) )
    StepNo = 0

    ALLOCATE( dt(0:nLevels-1) )
    dt = -100.0e0_AR

    ALLOCATE( t(0:nLevels-1) )
    t = 0.0e0_AR


    CALL InitializeDataAMReX( nLevels )    

  END SUBROUTINE MyAmrInit


  SUBROUTINE MyAmrFinalize

    CALL FinalizeDataAMReX( nLevels )

    DEALLOCATE( GEOM )
    DEALLOCATE( DM   )
    DEALLOCATE( BA   )

    DEALLOCATE( t      )
    DEALLOCATE( dt     )
    DEALLOCATE( StepNo )

  END SUBROUTINE MyAmrFinalize



END MODULE MyAmrModule
