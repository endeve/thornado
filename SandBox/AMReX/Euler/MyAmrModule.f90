MODULE MyAmrModule

  ! --- AMReX Modules ---
  USE amrex_fort_module,      ONLY: &
    amrex_real, amrex_spacedim
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
    Millisecond, Kilometer

  ! --- Local Modules ---
  USE MyAmrDataModule

  IMPLICIT NONE

  REAL(amrex_real)                    :: t_end, Gamma_IDEAL, CFL
  REAL(amrex_real)                    :: t_wrt, dt_wrt, t_chk, dt_chk
  INTEGER                             :: nNodes, nStages, nLevels, coord_sys
  INTEGER                             :: iCycleD, iCycleW, iCycleChk, iRestart
  INTEGER,          ALLOCATABLE       :: MaxGridSize(:), nX(:), swX(:), bcX(:)
  REAL(amrex_real), ALLOCATABLE       :: xL(:), xR(:), dt(:), t(:)
  CHARACTER(LEN=:), ALLOCATABLE       :: ProgramName
  INTEGER,          ALLOCATABLE, SAVE :: StepNo(:)
  CHARACTER(LEN=32),             SAVE :: Coordsys
  LOGICAL,                       SAVE :: DEBUG, UsePhysicalUnits

  ! --- Slope limiter ---
  LOGICAL          :: UseSlopeLimiter
  LOGICAL          :: UseCharacteristicLimiting
  LOGICAL          :: UseTroubledCellIndicator
  REAL(amrex_real) :: SlopeTolerance
  REAL(amrex_real) :: BetaTVD, BetaTVB
  REAL(amrex_real) :: LimiterThresholdParameter
  LOGICAL          :: UseConservativeCorrection

  ! --- Positivity limiter ---
  LOGICAL          :: UsePositivityLimiter
  REAL(amrex_real) :: Min_1, Min_2, Min_3
  REAL(amrex_real) :: Max_1, Max_2, Max_3

  ! --- Equation Of State ---
  CHARACTER(LEN=:), ALLOCATABLE :: EquationOfState
  CHARACTER(LEN=:), ALLOCATABLE :: EosTableName

  ! --- AMReX Geometry arrays ---
  TYPE(amrex_boxarray),  ALLOCATABLE, PUBLIC :: BA(:)
  TYPE(amrex_distromap), ALLOCATABLE, PUBLIC :: DM(:)
  TYPE(amrex_geometry),  ALLOCATABLE, PUBLIC :: GEOM(:)

CONTAINS


  SUBROUTINE MyAmrInit

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
    Gamma_IDEAL = 5.0_amrex_real / 3.0_amrex_real
    ! --- thornado paramaters thornado.* ---
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get   ( 'dt_wrt',           dt_wrt )
      CALL PP % get   ( 'dt_chk',           dt_chk )
      CALL PP % get   ( 't_end',            t_end )
      CALL PP % get   ( 'nNodes',           nNodes )
      CALL PP % get   ( 'nStages',          nStages )
      CALL PP % get   ( 'CFL',              CFL )
      CALL PP % get   ( 'ProgramName',      ProgramName )
      CALL PP % query ( 'Gamma',            Gamma_IDEAL )
      CALL PP % getarr( 'bcX',              bcX )
      CALL PP % getarr( 'swX',              swX )
      CALL PP % get   ( 'iCycleD',          iCycleD )
      CALL PP % get   ( 'iCycleW',          iCycleW )
      CALL PP % get   ( 'iCycleChk',        iCycleChk )
      CALL PP % get   ( 'iRestart',         iRestart )
      CALL PP % query ( 'UsePhysicalUnits', UsePhysicalUnits )
    CALL amrex_parmparse_destroy( PP )
    IF( iCycleW .GT. 0 .AND. dt_wrt .GT. 0.0_amrex_real )THEN
      WRITE(*,'(A)') 'iCycleW and dt_wrt cannot both be greater than zero.'
      WRITE(*,'(A)') 'Stopping...'
      STOP
    END IF
    IF( iCycleChk .GT. 0 .AND. dt_chk .GT. 0.0_amrex_real )THEN
      WRITE(*,'(A)') 'iCycleChk and dt_chk cannot both be greater than zero.'
      WRITE(*,'(A)') 'Stopping...'
      STOP
    END IF

    CFL = CFL &
            / ( amrex_spacedim * ( 2.0_amrex_real * nNodes - 1.0_amrex_real ) )

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

    IF( ProgramName .EQ. 'StandingAccretionShock_Relativistic' )THEN

      t_end  = t_end  * Millisecond
      dt_wrt = dt_wrt * Millisecond
      dt_chk = dt_chk * Millisecond
      xL(1)  = xL(1)  * Kilometer
      xR(1)  = xR(1)  * Kilometer

    END IF

    ! --- Parameters amr.* ---
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr( 'n_cell',        nX )
      CALL PP % getarr( 'max_grid_size', MaxGridSize )
      CALL PP % get   ( 'max_level',     nLevels )
    CALL amrex_parmparse_destroy( PP )

    ! --- Slope limiter parameters SL.* ---
    CALL amrex_parmparse_build( PP, 'SL' )
      CALL PP % get( 'UseSlopeLimiter',           UseSlopeLimiter )
      CALL PP % get( 'UseCharacteristicLimiting', UseCharacteristicLimiting )
      CALL PP % get( 'UseTroubledCellIndicator',  UseTroubledCellIndicator )
      CALL PP % get( 'SlopeTolerance',            SlopeTolerance )
      CALL PP % get( 'BetaTVD',                   BetaTVD )
      CALL PP % get( 'BetaTVB',                   BetaTVB )
      CALL PP % get( 'LimiterThresholdParameter', LimiterThresholdParameter )
      CALL PP % get( 'UseConservativeCorrection', UseConservativeCorrection )
    CALL amrex_parmparse_destroy( PP )

    ! --- Positivitiy limiter parameters PL.* ---
    Min_3 = 0.0_amrex_real
    Max_1 = 0.0_amrex_real
    Max_2 = 0.0_amrex_real
    Max_3 = 0.0_amrex_real
    CALL amrex_parmparse_build( PP, 'PL' )
      CALL PP % get( 'UsePositivityLimiter', UsePositivityLimiter )
      CALL PP % get( 'Min_1',                Min_1 )
      CALL PP % get( 'Min_2',                Min_2 )
      CALL PP % query( 'Min_3',              Min_3 )
      CALL PP % query( 'Max_1',              Max_1 )
      CALL PP % query( 'Max_2',              Min_2 )
      CALL PP % query( 'Max_3',              Max_3 )
    CALL amrex_parmparse_destroy( PP )

    ! --- Equation of state parameters EoS.* ---
    EquationOfState = 'IDEAL'
    EosTableName    = ''
    CALL amrex_parmparse_build( PP, 'EoS' )
      CALL PP % query( 'EquationOfState', EquationOfState )
      CALL PP % query( 'EosTableName',    EosTableName    )
    CALL amrex_parmparse_destroy( PP )

    IF( EquationOfState .EQ. 'TABLE' )THEN

      t_end  = t_end  * Millisecond
      dt_wrt = dt_wrt * Millisecond
      dt_chk = dt_chk * Millisecond

      IF     ( CoordSys .EQ. 'CARTESIAN'   )THEN
        xL(1:3) = xL(1:3) * Kilometer
        xR(1:3) = xR(1:3) * Kilometer
      ELSE IF( CoordSys .EQ. 'CYLINDRICAL' )THEN
        xL(1:2) = xL(1:2) * Kilometer
        xR(1:2) = xR(1:2) * Kilometer
      ELSE IF( CoordSys .EQ. 'SPHERICAL'   )THEN
        xL(1)   = xL(1)   * Kilometer
        xR(1)   = xR(1)   * Kilometer
      ELSE
        STOP 'Invalid choice for CoordSys'
      END IF

    END IF

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

    ALLOCATE( StepNo(0:nLevels) )
    StepNo = 0

    ALLOCATE( dt(0:nLevels) )
    dt = -100.0e0_amrex_real

    ALLOCATE( t(0:nLevels) )
    t = 0.0e0_amrex_real

    CALL InitializeDataAMReX( nLevels )

  END SUBROUTINE MyAmrInit


  SUBROUTINE MyAmrFinalize

    CALL FinalizeDataAMReX( nLevels )

    DEALLOCATE( GEOM )
    DEALLOCATE( DM )
    DEALLOCATE( BA )

    DEALLOCATE( t )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

  END SUBROUTINE MyAmrFinalize


END MODULE MyAmrModule
