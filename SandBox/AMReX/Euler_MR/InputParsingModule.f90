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
    InitializeProgramHeader, &
    nDimsX
  USE UnitsModule, ONLY: &
    ActivateUnitsDisplay, &
    DescribeUnitsDisplay, &
    UnitsDisplay
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem

  ! --- Local modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two
  USE MF_Euler_ErrorModule, ONLY: &
    DescribeError_Euler_MF

  IMPLICIT NONE

  ! -- thornado ---

  CHARACTER(:), ALLOCATABLE :: ProgramName
  INTEGER     , ALLOCATABLE :: swX(:)
  INTEGER     , ALLOCATABLE :: bcX(:)
  INTEGER                   :: nNodes
  INTEGER                   :: nStages
  REAL(DP)                  :: t_wrt, t_chk, dt_wrt, dt_chk
  INTEGER                   :: iCycleW, iCycleChk, iCycleD, iRestart
  REAL(DP)                  :: t_end
  REAL(DP)                  :: CFL
  LOGICAL     , SAVE        :: UsePhysicalUnits

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

  ! --- Equation of State ---

  CHARACTER(:), ALLOCATABLE :: EquationOfState
  CHARACTER(:), ALLOCATABLE :: EosTableName
  REAL(DP)                  :: Gamma_IDEAL

  ! --- geometry ---

  INTEGER               :: coord_sys
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
  INTEGER :: BlockingFactor(3)
  INTEGER :: MaxLevel
  INTEGER :: nLevels
  LOGICAL :: UseTiling
  LOGICAL :: do_reflux
  INTEGER , ALLOCATABLE :: RefinementRatio(:)
  INTEGER , ALLOCATABLE :: StepNo(:)
  INTEGER , ALLOCATABLE :: nRefinementBuffer(:)
  REAL(DP), ALLOCATABLE :: TagCriteria(:)

  REAL(DP), ALLOCATABLE :: dt   (:)
  REAL(DP), ALLOCATABLE :: t_old(:)
  REAL(DP), ALLOCATABLE :: t_new(:)
  CHARACTER(:), ALLOCATABLE :: PlotFileBaseName
  INTEGER :: iOS_CPP(3)

  LOGICAL                   :: WriteNodalData
  CHARACTER(:), ALLOCATABLE :: NodalDataFileName

CONTAINS


  SUBROUTINE InitializeParameters

    TYPE(amrex_parmparse) :: PP

    IF( .NOT. amrex_initialized() ) &
      CALL amrex_init()

    IF( .NOT. amrex_amrcore_initialized() ) &
      CALL amrex_amrcore_init()

    WriteNodalData    = .FALSE.
    NodalDataFileName = ''
    CALL amrex_parmparse_build( PP, 'debug' )
      CALL PP % query( 'WriteNodalData', WriteNodalData )
      CALL PP % query( 'NodalDataFileName', NodalDataFileName )
    CALL amrex_parmparse_destroy( PP )

    ! --- thornado Parameters thornado.* ---

    UsePhysicalUnits = .FALSE.
    PlotFileBaseName = 'plt'
    iCycleD          = 10
    iCycleW          = -1
    iCycleChk        = -1
    iRestart         = -1
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
      CALL PP % query ( 'iCycleD', iCycleD )
      CALL PP % query ( 'PlotFileBaseName', PlotFileBaseName )
      CALL PP % query ( 'iCycleW', iCycleW )
      CALL PP % query ( 'iCycleChk', iCycleChk )
      CALL PP % query ( 'iRestart', iRestart )
      CALL PP % query ( 'dt_wrt', dt_wrt )
      CALL PP % query ( 'dt_chk', dt_chk )
      CALL PP % query ( 'UsePhysicalUnits', UsePhysicalUnits      )
    CALL amrex_parmparse_destroy( PP )

    IF( iCycleW * dt_wrt .GT. Zero )THEN

      CALL DescribeError_Euler_MF( 101 )

    END IF

    IF( iCycleChk * dt_chk .GT. Zero )THEN

      CALL DescribeError_Euler_MF( 102 )

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

    ! --- Equation of State Parameters EoS.* ---

    EquationOfState = 'IDEAL'
    Gamma_IDEAL     = 4.0_DP / 3.0_DP
    EosTableName    = ''
    CALL amrex_parmparse_build( PP, 'EoS' )
      CALL PP % query ( 'EquationOfState', EquationOfState )
      CALL PP % query ( 'Gamma_IDEAL', Gamma_IDEAL )
      CALL PP % query ( 'EosTableName', EosTableName )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters geometry.* ---

    CALL amrex_parmparse_build( PP, 'geometry' )
      CALL PP % get   ( 'coord_sys', coord_sys )
      CALL PP % getarr( 'prob_lo'  , xL       )
      CALL PP % getarr( 'prob_hi'  , xR       )
    CALL amrex_parmparse_destroy( PP )

    IF     ( coord_sys .EQ. 0 )THEN

      CoordinateSystem = 'CARTESIAN'

    ELSE IF( coord_sys .EQ. 1 )THEN

      CoordinateSystem = 'CYLINDRICAL'

    ELSE IF( coord_sys .EQ. 2 )THEN

      CoordinateSystem = 'SPHERICAL'

    ELSE

      CALL DescribeError_Euler_MF &
             ( 103, Message_Option = 'Invalid CoordSys:', &
                      Int_Option = [ coord_sys ] )

    END IF

    IF( UsePhysicalUnits )THEN

      CALL ActivateUnitsDisplay &
             ( CoordinateSystem_Option = TRIM( CoordinateSystem ) )

      t_end  = t_end  * UnitsDisplay % TimeUnit
      dt_wrt = dt_wrt * UnitsDisplay % TimeUnit
      dt_chk = dt_chk * UnitsDisplay % TimeUnit

      xL(1) = xL(1) * UnitsDisplay % LengthX1Unit
      xR(1) = xR(1) * UnitsDisplay % LengthX1Unit
      xL(2) = xL(2) * UnitsDisplay % LengthX2Unit
      xR(2) = xR(2) * UnitsDisplay % LengthX2Unit
      xL(3) = xL(3) * UnitsDisplay % LengthX3Unit
      xR(3) = xR(3) * UnitsDisplay % LengthX3Unit

    END IF

    ! --- Parameters amr.* ---

    IF( amrex_spacedim .EQ. 1 )THEN

      MaxGridSizeX1 = 128
      MaxGridSizeX2 = 1
      MaxGridSizeX3 = 1
      BlockingFactorX1 = 8
      BlockingFactorX2 = 1
      BlockingFactorX3 = 1

    ELSE IF( amrex_spacedim .EQ. 2 )THEN

      MaxGridSizeX1 = 128
      MaxGridSizeX2 = 128
      MaxGridSizeX3 = 1
      BlockingFactorX1 = 8
      BlockingFactorX2 = 8
      BlockingFactorX3 = 1

    ELSE

      MaxGridSizeX1 = 32
      MaxGridSizeX2 = 32
      MaxGridSizeX3 = 32
      BlockingFactorX1 = 8
      BlockingFactorX2 = 8
      BlockingFactorX3 = 8

    END IF
    UseTiling = 0
    do_reflux = 0
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr  ( 'n_cell'           , nX                )
      CALL PP % query   ( 'max_grid_size_x'  , MaxGridSizeX1     )
      CALL PP % query   ( 'max_grid_size_y'  , MaxGridSizeX2     )
      CALL PP % query   ( 'max_grid_size_z'  , MaxGridSizeX3     )
      CALL PP % query   ( 'blocking_factor_x', BlockingFactorX1  )
      CALL PP % query   ( 'blocking_factor_y', BlockingFactorX2  )
      CALL PP % query   ( 'blocking_factor_z', BlockingFactorX3  )
      CALL PP % get     ( 'max_level'        , MaxLevel          )
      CALL PP % query   ( 'UseTiling'        , UseTiling         )
      CALL PP % query   ( 'do_reflux'        , do_reflux         )
      CALL PP % getarr  ( 'ref_ratio'        , RefinementRatio   )
      IF( MaxLevel .GT. 0 )THEN
        CALL PP % getarr( 'TagCriteria'      , TagCriteria       )
        CALL PP % getarr( 'n_error_buf'      , nRefinementBuffer )
      END IF
    CALL amrex_parmparse_destroy( PP )

    MaxGridSizeX   = [ MaxGridSizeX1   , MaxGridSizeX2   , MaxGridSizeX3    ]
    BlockingFactor = [ BlockingFactorX1, BlockingFactorX2, BlockingFactorX3 ]

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

    IF( amrex_parallel_ioprocessor() ) &
      CALL DescribeUnitsDisplay

    IF( nDimsX .NE. amrex_spacedim )THEN

      CALL DescribeError_Euler_MF( 104, &
                                   Int_Option = [ nDimsX, amrex_spacedim ] )

    END IF

    iOS_CPP = 0

    iOS_CPP(1) = 1

    IF( amrex_spacedim .GT. 1 ) iOS_CPP(2) = 1
    IF( amrex_spacedim .GT. 2 ) iOS_CPP(3) = 1

  END SUBROUTINE InitializeParameters

END MODULE InputParsingModule
