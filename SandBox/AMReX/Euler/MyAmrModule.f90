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
  USE ProgramHeaderModule,  ONLY: &
    InitializeProgramHeader, nDOFX, nDimsX
  USE FluidFieldsModule,    ONLY: &
    nCF, nPF, nAF
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE UnitsModule,          ONLY: &
    ActivateUnitsDisplay, &
    DescribeUnitsDisplay, &
    UnitsDisplay

  ! --- Local Modules ---
  USE MyAmrDataModule, ONLY: &
    InitializeDataAMReX, &
    FinalizeDataAMReX

  IMPLICIT NONE

  ! --- thornado ---
  REAL(AR)                       :: t_end, t_wrt, dt_wrt, t_chk, dt_chk
  REAL(AR),          ALLOCATABLE :: t(:), dt(:)
  REAL(AR)                       :: CFL
  INTEGER                        :: nNodes, nStages
  INTEGER                        :: iCycleD, iCycleW, iCycleChk, iRestart
  INTEGER,           ALLOCATABLE :: nX(:), swX(:), bcX(:)
  REAL(AR),          ALLOCATABLE :: xL(:), xR(:)
  CHARACTER(LEN=:),  ALLOCATABLE :: ProgramName
  CHARACTER(LEN=32), SAVE        :: CoordSys
  LOGICAL,           SAVE        :: UsePhysicalUnits
  LOGICAL,           SAVE        :: DEBUG

  ! --- Slope limiter ---
  LOGICAL  :: UseSlopeLimiter
  LOGICAL  :: UseCharacteristicLimiting
  LOGICAL  :: UseTroubledCellIndicator
  REAL(AR) :: SlopeTolerance
  REAL(AR) :: BetaTVD, BetaTVB
  REAL(AR) :: LimiterThresholdParameter
  LOGICAL  :: UseConservativeCorrection

  ! --- Positivity limiter ---
  LOGICAL  :: UsePositivityLimiter
  REAL(AR) :: Min_1, Min_2
  REAL(AR) :: Max_1, Max_2

  ! --- Equation Of State ---
  REAL(AR)                      :: Gamma_IDEAL
  CHARACTER(LEN=:), ALLOCATABLE :: EquationOfState
  CHARACTER(LEN=:), ALLOCATABLE :: EosTableName

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
      CALL PP % get   ( 'nStages',          nStages )
      CALL PP % get   ( 'CFL',              CFL )
      CALL PP % get   ( 'ProgramName',      ProgramName )
      CALL PP % getarr( 'bcX',              bcX )
      CALL PP % getarr( 'swX',              swX )
      CALL PP % get   ( 'iCycleD',          iCycleD )
      CALL PP % get   ( 'iCycleW',          iCycleW )
      CALL PP % get   ( 'iCycleChk',        iCycleChk )
      CALL PP % get   ( 'iRestart',         iRestart )
      CALL PP % query ( 'UsePhysicalUnits', UsePhysicalUnits )
    CALL amrex_parmparse_destroy( PP )

    IF( UsePhysicalUnits ) &
      CALL ActivateUnitsDisplay
    IF( amrex_parallel_ioprocessor() ) &
      CALL DescribeUnitsDisplay

    IF( iCycleW .GT. 0 .AND. dt_wrt .GT. Zero )THEN
      WRITE(*,'(A)') 'iCycleW and dt_wrt cannot both be greater than zero.'
      WRITE(*,'(A)') 'Stopping...'
      STOP
    END IF
    IF( iCycleChk .GT. 0 .AND. dt_chk .GT. Zero )THEN
      WRITE(*,'(A)') 'iCycleChk and dt_chk cannot both be greater than zero.'
      WRITE(*,'(A)') 'Stopping...'
      STOP
    END IF

    CFL = CFL / ( amrex_spacedim * ( Two * nNodes - One ) )

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

    IF( UsePhysicalUnits )THEN

      t_end  = t_end  * UnitsDisplay % TimeUnit
      dt_wrt = dt_wrt * UnitsDisplay % TimeUnit
      dt_chk = dt_chk * UnitsDisplay % TimeUnit

      IF     ( CoordSys .EQ. 'CARTESIAN'   )THEN

        xL(1:3) = xL(1:3) * UnitsDisplay % LengthUnit
        xR(1:3) = xR(1:3) * UnitsDisplay % LengthUnit

      ELSE IF( CoordSys .EQ. 'CYLINDRICAL' )THEN

        xL(1:2) = xL(1:2) * UnitsDisplay % LengthUnit
        xR(1:2) = xR(1:2) * UnitsDisplay % LengthUnit

      ELSE IF( CoordSys .EQ. 'SPHERICAL'   )THEN

        xL(1)   = xL(1) * UnitsDisplay % LengthUnit
        xR(1)   = xR(1) * UnitsDisplay % LengthUnit

      ELSE

        STOP 'Invalid choice for CoordSys'

      END IF

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
    MaxGridSizeX = [ MaxGridSizeX1, MaxGridSizeX2, MaxGridSizeX3 ]
    nLevels = MaxLevel + 1

    ! --- Slope limiter parameters SL.* ---
    UseSlopeLimiter           = .TRUE.
    UseCharacteristicLimiting = .TRUE.
    UseTroubledCellIndicator  = .TRUE.
    SlopeTolerance            = 1.0e-6_AR
    BetaTVD                   = 1.75_AR
    BetaTVB                   = Zero
    LimiterThresholdParameter = 0.03_AR
    UseConservativeCorrection = .TRUE.
    CALL amrex_parmparse_build( PP, 'SL' )
      CALL PP % query( 'UseSlopeLimiter',           UseSlopeLimiter )
      CALL PP % query( 'UseCharacteristicLimiting', UseCharacteristicLimiting )
      CALL PP % query( 'UseTroubledCellIndicator',  UseTroubledCellIndicator )
      CALL PP % query( 'SlopeTolerance',            SlopeTolerance )
      CALL PP % query( 'BetaTVD',                   BetaTVD )
      CALL PP % query( 'BetaTVB',                   BetaTVB )
      CALL PP % query( 'LimiterThresholdParameter', LimiterThresholdParameter )
      CALL PP % query( 'UseConservativeCorrection', UseConservativeCorrection )
    CALL amrex_parmparse_destroy( PP )

    ! --- Positivitiy limiter parameters PL.* ---
    UsePositivityLimiter = .TRUE.
    Min_1                = 1.0e-12_AR
    Min_2                = 1.0e-12_AR
    CALL amrex_parmparse_build( PP, 'PL' )
      CALL PP % query( 'UsePositivityLimiter', UsePositivityLimiter )
      CALL PP % query( 'Min_1',                Min_1 )
      CALL PP % query( 'Min_2',                Min_2 )
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

    CALL InitializeProgramHeader &
           ( ProgramName_Option = TRIM( ProgramName ), &
             nNodes_Option = nNodes, nX_Option = nX, swX_Option = swX, &
             xL_Option = xL, xR_Option = xR, bcX_Option = bcX, &
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
