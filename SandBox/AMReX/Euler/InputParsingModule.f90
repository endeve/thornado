MODULE InputParsingModule

  ! --- AMReX modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_init_module, ONLY: &
    amrex_init, &
    amrex_initialized
  USE amrex_amr_module, ONLY: &
    amrex_amrcore_init, &
    amrex_amrcore_initialized

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    InitializeProgramHeader
  USE UtilitiesModule, ONLY: &
    thornado_abort

  ! --- Local modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two

  IMPLICIT NONE

  ! -- thornado ---

  CHARACTER(:), ALLOCATABLE :: ProgramName
  INTEGER, ALLOCATABLE  :: swX(:)
  INTEGER, ALLOCATABLE  :: bcX(:)
  INTEGER               :: nNodes
  INTEGER               :: nStages
  REAL(DP)              :: dt_wrt, dt_chk
  INTEGER               :: iCycleW, iCycleChk, iCycleD
  REAL(DP)              :: t_end
  REAL(DP)              :: Gamma_IDEAL
  REAL(DP)              :: CFL
  CHARACTER(:), ALLOCATABLE :: EquationOfState
  CHARACTER(:), ALLOCATABLE :: EosTableName

  ! --- Boundary Conditions ---

  INTEGER, ALLOCATABLE, PUBLIC, SAVE :: lo_bc(:,:)
  INTEGER, ALLOCATABLE, PUBLIC, SAVE :: hi_bc(:,:)
  INTEGER, ALLOCATABLE, PUBLIC, SAVE :: lo_bc_uCF(:,:)
  INTEGER, ALLOCATABLE, PUBLIC, SAVE :: hi_bc_uCF(:,:)

  ! --- Slope Limiter ---

  LOGICAL                       :: UseSlopeLimiter
  CHARACTER(LEN=:), ALLOCATABLE :: SlopeLimiterMethod
  REAL(DP)                      :: BetaTVD, BetaTVB, SlopeTolerance
  LOGICAL                       :: UseCharacteristicLimiting
  LOGICAL                       :: UseTroubledCellIndicator
  REAL(DP)                      :: LimiterThresholdParameter
  LOGICAL                       :: UseConservativeCorrection

  ! --- Positivity Limiter ---

  LOGICAL  :: UsePositivityLimiter
  REAL(DP) :: Min_1, Min_2
  REAL(DP) :: Max_1, Max_2

  ! --- geometry ---

  INTEGER               :: CoordSys
  REAL(DP), ALLOCATABLE :: xL(:), xR(:)

  ! --- amr ---

  INTEGER, ALLOCATABLE :: nX(:)
  INTEGER :: MaxGridSizeX1
  INTEGER :: MaxGridSizeX2
  INTEGER :: MaxGridSizeX3
  INTEGER :: BlockingFactorX1
  INTEGER :: BlockingFactorX2
  INTEGER :: BlockingFactorX3
  INTEGER :: MaxGridSizeX(3)
  INTEGER :: MaxLevel
  INTEGER :: nLevels
  LOGICAL :: UseTiling
  LOGICAL :: do_reflux
  INTEGER, ALLOCATABLE :: RefinementRatio(:)
  INTEGER, ALLOCATABLE :: StepNo(:)

  REAL(DP), ALLOCATABLE :: dt   (:)
  REAL(DP), ALLOCATABLE :: t_old(:)
  REAL(DP), ALLOCATABLE :: t_new(:)
  CHARACTER(:), ALLOCATABLE :: PlotFileBaseName
  INTEGER :: iOS_CPP(3)

CONTAINS


  SUBROUTINE InitializeParameters

    TYPE(amrex_parmparse) :: PP

    IF( .NOT. amrex_initialized() ) &
      CALL amrex_init()

    IF( .NOT. amrex_amrcore_initialized() ) &
      CALL amrex_amrcore_init()

    ! --- thornado Parameters thornado.* ---

    EquationOfState  = 'IDEAL'
    Gamma_IDEAL      = 4.0_DP / 3.0_DP
    EosTableName     = ''
    PlotFileBaseName = 'plt'
    iCycleD          = 10
    iCycleW          = -1
    iCycleChk        = -1
    dt_wrt           = -1.0_DP
    dt_chk           = -1.0_DP
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get   ( 'ProgramName', ProgramName )
      CALL PP % get   ( 'nNodes', nNodes )
      CALL PP % get   ( 'nStages', nStages )
      CALL PP % getarr( 'swX', swX )
      CALL PP % getarr( 'bcX', bcX )
      CALL PP % get   ( 't_end', t_end )
      CALL PP % get   ( 'CFL', CFL )
      CALL PP % query ( 'EquationOfState', EquationOfState )
      CALL PP % query ( 'Gamma_IDEAL', Gamma_IDEAL )
      CALL PP % query ( 'iCycleD', iCycleD )
      CALL PP % query ( 'EosTableName', EosTableName )
      CALL PP % query ( 'PlotFileBaseName', PlotFileBaseName )
      CALL PP % query ( 'iCycleW', iCycleW )
      CALL PP % query ( 'iCycleChk', iCycleChk )
      CALL PP % query ( 'dt_wrt', dt_wrt )
      CALL PP % query ( 'dt_chk', dt_chk )
    CALL amrex_parmparse_destroy( PP )

    IF( iCycleW * dt_wrt .GT. Zero )THEN

      WRITE(*,'(A)') 'iCycleW and dt_wrt cannot both be greater than zero.'
      WRITE(*,'(A)') 'Stopping...'
      CALL thornado_abort

    END IF

    IF( iCycleChk * dt_chk .GT. Zero )THEN

      WRITE(*,'(A)') 'iCycleChk and dt_chk cannot both be greater than zero.'
      WRITE(*,'(A)') 'Stopping...'
      CALL thornado_abort

    END IF

    CFL = CFL / ( DBLE( amrex_spacedim ) * ( Two * DBLE( nNodes ) - One ) )

    ! --- Slope Limiter Parameters SL.* ---

    UseSlopeLimiter           = .TRUE.
    SlopeLimiterMethod        = 'TVD'
    BetaTVD                   = 1.75_DP
    BetaTVB                   = Zero
    SlopeTolerance            = 1.0e-6_DP
    UseCharacteristicLimiting = .TRUE.
    UseTroubledCellIndicator  = .TRUE.
    LimiterThresholdParameter = 0.03_DP
    UseConservativeCorrection = .TRUE.
    CALL amrex_parmparse_build( PP, 'SL' )
      CALL PP % query( 'UseSlopeLimiter'          , UseSlopeLimiter           )
      CALL PP % query( 'SlopeLimiterMethod'       , SlopeLimiterMethod        )
      CALL PP % query( 'BetaTVD'                  , BetaTVD                   )
      CALL PP % query( 'BetaTVB'                  , BetaTVB                   )
      CALL PP % query( 'SlopeTolerance'           , SlopeTolerance            )
      CALL PP % query( 'UseCharacteristicLimiting', UseCharacteristicLimiting )
      CALL PP % query( 'UseTroubledCellIndicator' , UseTroubledCellIndicator  )
      CALL PP % query( 'LimiterThresholdParameter', LimiterThresholdParameter )
      CALL PP % query( 'UseConservativeCorrection', UseConservativeCorrection )
    CALL amrex_parmparse_destroy( PP )

    ! --- Positivity Limiter Parameters PL.* ---

    UsePositivityLimiter = .TRUE.
    Min_1                = 1.0e-12_DP
    Min_2                = 1.0e-12_DP
    CALL amrex_parmparse_build( PP, 'PL' )
      CALL PP % query( 'UsePositivityLimiter', UsePositivityLimiter )
      CALL PP % query( 'Min_1'               , Min_1                )
      CALL PP % query( 'Min_2'               , Min_2                )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters geometry.* ---

    CALL amrex_parmparse_build( PP, 'geometry' )
      CALL PP % get   ( 'coord_sys', CoordSys )
      CALL PP % getarr( 'prob_lo'  , xL       )
      CALL PP % getarr( 'prob_hi'  , xR       )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters amr.* ---

    MaxGridSizeX1    = 1
    MaxGridSizeX2    = 1
    MaxGridSizeX3    = 1
    BlockingFactorX1 = 1
    BlockingFactorX2 = 1
    BlockingFactorX3 = 1
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr( 'n_cell'           , nX               )
      CALL PP % query ( 'max_grid_size_x'  , MaxGridSizeX1    )
      CALL PP % query ( 'max_grid_size_y'  , MaxGridSizeX2    )
      CALL PP % query ( 'max_grid_size_z'  , MaxGridSizeX3    )
      CALL PP % query ( 'blocking_factor_x', BlockingFactorX1 )
      CALL PP % query ( 'blocking_factor_y', BlockingFactorX2 )
      CALL PP % query ( 'blocking_factor_z', BlockingFactorX3 )
      CALL PP % get   ( 'max_level'        , MaxLevel         )
      CALL PP % get   ( 'UseTiling'        , UseTiling        )
      CALL PP % get   ( 'do_reflux'        , do_reflux        )
      CALL PP % getarr( 'ref_ratio'        , RefinementRatio  )
    CALL amrex_parmparse_destroy( PP )

    MaxGridSizeX = [ MaxGridSizeX1, MaxGridSizeX2, MaxGridSizeX3 ]
    nLevels = MaxLevel + 1

    CALL InitializeProgramHeader &
           ( ProgramName_Option = TRIM( ProgramName ), &
             nNodes_Option      = nNodes,              &
             nX_Option          = nX,                  &
             swX_Option         = swX,                 &
             xL_Option          = xL,                  &
             xR_Option          = xR,                  &
             bcX_Option         = bcX,                 &
             Verbose_Option     = amrex_parallel_ioprocessor() )

    iOS_CPP = 0

    iOS_CPP(1) = 1

    IF( amrex_spacedim .GT. 1 ) iOS_CPP(2) = 1
    IF( amrex_spacedim .GT. 2 ) iOS_CPP(3) = 1

  END SUBROUTINE InitializeParameters

END MODULE InputParsingModule
