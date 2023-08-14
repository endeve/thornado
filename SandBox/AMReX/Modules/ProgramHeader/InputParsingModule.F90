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

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    InitializeProgramHeader, &
    bcZ, &
    nDimsX
  USE UnitsModule, ONLY: &
    ActivateUnitsDisplay, &
    UnitsDisplay, &
    SolarMass, &
    Centimeter
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem
  USE RadiationFieldsModule, ONLY: &
    SetNumberOfSpecies

  ! --- Local modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF

  IMPLICIT NONE

  ! -- thornado ---

  CHARACTER(:), ALLOCATABLE :: ProgramName
  INTEGER     , ALLOCATABLE :: swX(:)
  INTEGER     , ALLOCATABLE :: bcX(:)
  INTEGER                   :: nNodes
  REAL(DP)                  :: t_wrt, t_chk, dt_wrt, dt_chk, dt_rel
  INTEGER                   :: iCycleW, iCycleChk, iCycleD, iRestart, iReGrid
  REAL(DP)                  :: t_end
  LOGICAL     , SAVE        :: UsePhysicalUnits
  LOGICAL     , SAVE        :: DEBUG
  LOGICAL     , SAVE        :: SolveGravity_NR

  ! --- TimeStepping ---

  CHARACTER(:), ALLOCATABLE :: Scheme
  INTEGER                   :: nStages
  REAL(DP)                  :: CFL

  ! --- Transport ---

  INTEGER     , ALLOCATABLE :: bcZ_TwoMoment(:)
  INTEGER  :: nE, nSpecies, swE, bcE
  REAL(DP) :: eL, eR, zoomE

  ! --- Slope Limiter ---

  LOGICAL  :: UseSlopeLimiter_TwoMoment
  REAL(DP) :: BetaTVD_TwoMoment

  ! --- Positivity Limiter ---

  LOGICAL  :: UsePositivityLimiter_TwoMoment
  REAL(DP) :: Min_1_TwoMoment, Min_2_TwoMoment

  ! --- Equation of State ---

  CHARACTER(:), ALLOCATABLE :: EquationOfState
  CHARACTER(:), ALLOCATABLE :: EosTableName
  REAL(DP)                  :: Gamma_IDEAL

  ! --- Opacity Tables ---

  CHARACTER(:), ALLOCATABLE :: OpacityTableName_AbEm
  CHARACTER(:), ALLOCATABLE :: OpacityTableName_Iso
  CHARACTER(:), ALLOCATABLE :: OpacityTableName_NES
  CHARACTER(:), ALLOCATABLE :: OpacityTableName_Pair
  CHARACTER(:), ALLOCATABLE :: OpacityTableName_Brem

  ! --- Non-Linear Solver Parameters ---
  INTEGER  ::  M_outer
  INTEGER  ::  MaxIter_outer
  REAL(DP) ::  Rtol_outer
  INTEGER  ::  M_inner
  INTEGER  ::  MaxIter_inner
  REAL(DP) ::  Rtol_inner
  LOGICAL  ::  Include_NES
  LOGICAL  ::  Include_Pair
  LOGICAL  ::  Include_Brem
  LOGICAL  ::  Include_LinCorr
  REAL(DP), ALLOCATABLE ::  wMatterRHS(:)

  ! --- geometry ---

  INTEGER               :: coord_sys
  REAL(DP), ALLOCATABLE :: xL(:), xR(:)

  ! --- amr ---

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
  INTEGER :: nMaxLevels
  LOGICAL :: UseTiling
  LOGICAL :: UseFluxCorrection_Euler
  LOGICAL :: UseFluxCorrection_TwoMoment
  LOGICAL :: UseAMR
  INTEGER , ALLOCATABLE :: nX(:)
  INTEGER , ALLOCATABLE :: RefinementRatio(:)
  INTEGER , ALLOCATABLE :: StepNo(:)
  INTEGER , ALLOCATABLE :: nRefinementBuffer(:)
  REAL(DP), ALLOCATABLE :: TagCriteria(:)

  REAL(DP), ALLOCATABLE :: dt   (:), dt_TM(:)
  REAL(DP), ALLOCATABLE :: t_old(:)
  REAL(DP), ALLOCATABLE :: t_new(:)
  CHARACTER(:), ALLOCATABLE :: PlotFileNameRoot
  INTEGER :: iOS_CPP(3)

  LOGICAL                   :: WriteNodalData
  CHARACTER(:), ALLOCATABLE :: NodalDataFileName

real(dp)::mass,r0,kt,mu0,e0
real(dp)::d_0,chi,sigma
character(:),allocatable::direction

CONTAINS


  SUBROUTINE InitializeParameters

    TYPE(amrex_parmparse) :: PP

#if defined( THORNADO_AMREX_GIT_HASH )

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)

      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git hash:   ', THORNADO_AMREX_GIT_HASH
      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git date:   ', THORNADO_AMREX_GIT_DATE
      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git branch: ', THORNADO_AMREX_GIT_BRANCH
      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git url:    ', THORNADO_AMREX_GIT_URL

    END IF

#endif

mass=zero
r0=zero
e0=zero
mu0=zero
kt=zero
d_0=zero
chi=zero
sigma=zero
call amrex_parmparse_build( pp, 'ST' )
  call pp % query( 'mass',mass )
  call pp % query( 'r0',r0 )
  call pp % query( 'mu0',mu0 )
  call pp % query( 'e0',e0 )
  call pp % query( 'kt',kt )
  call pp % query( 'd_0',d_0 )
  call pp % query( 'chi',chi )
  call pp % query( 'sigma',sigma )
call amrex_parmparse_destroy( pp )
direction=''
call amrex_parmparse_build( pp, 'thornado' )
  call pp % query( 'direction',direction )
call amrex_parmparse_destroy( pp )

    ! --- debug Parameters debug.* ---

    DEBUG             = .FALSE.
    WriteNodalData    = .FALSE.
    NodalDataFileName = ''
    CALL amrex_parmparse_build( PP, 'debug' )
      CALL PP % query( 'DEBUG', &
                        DEBUG )
      CALL PP % query( 'WriteNodalData', &
                        WriteNodalData )
      CALL PP % query( 'NodalDataFileName', &
                        NodalDataFileName )
    CALL amrex_parmparse_destroy( PP )

    ! --- thornado Parameters thornado.* ---

    UsePhysicalUnits = .FALSE.
    PlotFileNameRoot = 'plt'
    iCycleD          = 10
    iCycleW          = -1
    iCycleChk        = -1
    iRestart         = -1
    dt_wrt           = -1.0_DP
    dt_chk           = -1.0_DP
    dt_rel           = 0.0_DP
    iReGrid          = 1
    SolveGravity_NR  = .FALSE.
    Scheme           = ''
    nE               = 1
    nSpecies         = 1
    swE              = 0
    bcE              = 0
    bcZ_TwoMoment    = [ 0, 0, 0, 0 ]
    eL               = 0.0_DP
    eR               = 1.0_DP
    zoomE            = 1.0_DP
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get   ( 'ProgramName', &
                         ProgramName )
      CALL PP % get   ( 'nNodes', &
                         nNodes )
      CALL PP % get   ( 'nStages', &
                         nStages )
      CALL PP % query ( 'Scheme', &
                         Scheme )
      CALL PP % getarr( 'swX', &
                         swX )
      CALL PP % getarr( 'bcX', &
                         bcX )
      CALL PP % queryarr( 'bcZ_TwoMoment', &
                           bcZ_TwoMoment )
      CALL PP % get   ( 't_end', &
                         t_end )
      CALL PP % get   ( 'CFL', &
                         CFL )
      CALL PP % query ( 'iCycleD', &
                         iCycleD )
      CALL PP % query ( 'PlotFileNameRoot', &
                         PlotFileNameRoot )
      CALL PP % query ( 'iCycleW', &
                         iCycleW )
      CALL PP % query ( 'iCycleChk', &
                         iCycleChk )
      CALL PP % query ( 'iRestart', &
                         iRestart )
      CALL PP % query ( 'dt_wrt', &
                         dt_wrt )
      CALL PP % query ( 'dt_chk', &
                         dt_chk )
      CALL PP % query ( 'dt_rel', &
                         dt_rel )
      CALL PP % query ( 'UsePhysicalUnits', &
                         UsePhysicalUnits )
      CALL PP % query ( 'iReGrid', &
                         iReGrid )
      CALL PP % query ( 'SolveGravity_NR', &
                         SolveGravity_NR )
      CALL PP % query ( 'nE', &
                         nE )
      CALL PP % query ( 'nSpecies', &
                         nSpecies )
      CALL PP % query ( 'swE', &
                         swE )
      CALL PP % query ( 'bcE', &
                         bcE )
      CALL PP % query ( 'eL', &
                         eL )
      CALL PP % query ( 'eR', &
                         eR )
      CALL PP % query ( 'zoomE', &
                         zoomE )
    CALL amrex_parmparse_destroy( PP )

    CALL SetNumberOfSpecies( nSpecies, Verbose_Option = amrex_parallel_ioprocessor() )

    IF( iCycleW * dt_wrt .GT. Zero ) &
      CALL DescribeError_MF &
             ( 101, Int_Option = [ iCycleW ], Real_Option = [ dt_wrt ] )

    IF( iCycleChk * dt_chk .GT. Zero ) &
      CALL DescribeError_MF &
             ( 102, Int_Option = [ iCycleChk ], Real_Option = [ dt_chk ] )

    CFL = CFL / ( DBLE( amrex_spacedim ) * ( Two * DBLE( nNodes ) - One ) )

    ! --- Slope Limiter Parameters SL.* ---

    UseSlopeLimiter_TwoMoment = .TRUE.
    BetaTVD_TwoMoment         = 1.75_DP
    CALL amrex_parmparse_build( PP, 'SL' )
      CALL PP % query( 'UseSlopeLimiter_TwoMoment', &
                        UseSlopeLimiter_TwoMoment )
      CALL PP % query( 'BetaTVD_TwoMoment', &
                        BetaTVD_TwoMoment )
    CALL amrex_parmparse_destroy( PP )

    ! --- Positivity Limiter Parameters PL.* ---

    UsePositivityLimiter_TwoMoment = .TRUE.
    Min_1_TwoMoment                = 1.0e-12_DP
    Min_2_TwoMoment                = 1.0e-12_DP
    CALL amrex_parmparse_build( PP, 'PL' )
      CALL PP % query( 'UsePositivityLimiter_TwoMoment', &
                        UsePositivityLimiter_TwoMoment )
      CALL PP % query( 'Min_1_TwoMoment', &
                        Min_1_TwoMoment )
      CALL PP % query( 'Min_2_TwoMoment', &
                        Min_2_TwoMoment )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters geometry.* ---

    CALL amrex_parmparse_build( PP, 'geometry' )
      CALL PP % get   ( 'coord_sys', &
                         coord_sys )
      CALL PP % getarr( 'prob_lo', &
                         xL )
      CALL PP % getarr( 'prob_hi', &
                         xR )
    CALL amrex_parmparse_destroy( PP )

    IF     ( coord_sys .EQ. 0 )THEN

      CoordinateSystem = 'CARTESIAN'

    ELSE IF( coord_sys .EQ. 1 )THEN

      CoordinateSystem = 'CYLINDRICAL'

    ELSE IF( coord_sys .EQ. 2 )THEN

      CoordinateSystem = 'SPHERICAL'

    ELSE

      CALL DescribeError_MF( 103, Int_Option = [ coord_sys ] )

    END IF

    IF( UsePhysicalUnits )THEN

      CALL ActivateUnitsDisplay &
             ( CoordinateSystem_Option = TRIM( CoordinateSystem ) )

      t_end  = t_end  * UnitsDisplay % TimeUnit
      dt_wrt = dt_wrt * UnitsDisplay % TimeUnit
      dt_chk = dt_chk * UnitsDisplay % TimeUnit
      dt_rel = dt_rel * UnitsDisplay % TimeUnit

      xL(1) = xL(1) * UnitsDisplay % LengthX1Unit
      xR(1) = xR(1) * UnitsDisplay % LengthX1Unit
      xL(2) = xL(2) * UnitsDisplay % LengthX2Unit
      xR(2) = xR(2) * UnitsDisplay % LengthX2Unit
      xL(3) = xL(3) * UnitsDisplay % LengthX3Unit
      xR(3) = xR(3) * UnitsDisplay % LengthX3Unit

      eL = eL * UnitsDisplay % EnergyUnit
      eR = eR * UnitsDisplay % EnergyUnit

      Chi = Chi * ( 1.0_DP / Centimeter )

      Mass = Mass * SolarMass
      E0 = E0 * UnitsDisplay % EnergyUnit
      mu0 = mu0 * UnitsDisplay % EnergyUnit
      kT = kT * UnitsDisplay % EnergyUnit
      R0 = R0 * UnitsDisplay % LengthX1Unit



    END IF

    ! --- Equation of State Parameters EoS.* ---

    EquationOfState = 'IDEAL'
    Gamma_IDEAL     = 4.0_DP / 3.0_DP
    EosTableName    = ''
    CALL amrex_parmparse_build( PP, 'EoS' )
      CALL PP % query ( 'EquationOfState', &
                         EquationOfState )
      CALL PP % query ( 'Gamma_IDEAL', &
                         Gamma_IDEAL )
      CALL PP % query ( 'EosTableName', &
                         EosTableName )
    CALL amrex_parmparse_destroy( PP )

    ! --- Opacity table parameters OP.* ---

    OpacityTableName_AbEm = ''
    OpacityTableName_Iso  = ''
    OpacityTableName_NES  = ''
    OpacityTableName_Pair = ''
    OpacityTableName_Brem = ''
    CALL amrex_parmparse_build( PP, 'OP' )
      CALL PP % query( 'OpacityTableName_AbEm', &
                        OpacityTableName_AbEm )
      CALL PP % query( 'OpacityTableName_Iso', &
                        OpacityTableName_Iso )
      CALL PP % query( 'OpacityTableName_NES', &
                        OpacityTableName_NES )
      CALL PP % query( 'OpacityTableName_Pair', &
                        OpacityTableName_Pair )
      CALL PP % query( 'OpacityTableName_Brem', &
                        OpacityTableName_Brem )
    CALL amrex_parmparse_destroy( PP )

    ! --- Non-Linear Solver parameters NL.* ---

    M_outer         = 2
    MaxIter_outer   = 100
    Rtol_outer      = 1.0d-8
    M_inner         = 2
    MaxIter_inner   = 100
    Rtol_inner      = 1.0d-8
    Include_NES     = .FALSE.
    Include_Pair    = .FALSE.
    Include_Brem    = .FALSE.
    Include_LinCorr = .FALSE.
    wMatterRHS      = [ 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP ]
    CALL amrex_parmparse_build( PP, 'NL' )
      CALL PP % query( 'M_outer', &
                        M_outer )
      CALL PP % query( 'MaxIter_outer', &
                        MaxIter_outer )
      CALL PP % query( 'Rtol_outer', &
                        Rtol_outer )
      CALL PP % query( 'M_inner', &
                        M_inner )
      CALL PP % query( 'MaxIter_inner', &
                        MaxIter_inner )
      CALL PP % query( 'Rtol_inner', &
                        Rtol_inner )
      CALL PP % query( 'Include_NES', &
                        Include_NES )
      CALL PP % query( 'Include_Pair', &
                        Include_Pair )
      CALL PP % query( 'Include_Brem', &
                        Include_Brem )
      CALL PP % query( 'Include_LinCorr', &
                        Include_LinCorr )
      CALL PP % queryarr( 'wMatterRHS', &
                           wMatterRHS )
    CALL amrex_parmparse_destroy( PP )

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
    UseAMR                      = .FALSE.
    UseFluxCorrection_Euler     = .FALSE.
    UseFluxCorrection_TwoMoment = .FALSE.
    UseTiling                   = .FALSE.
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr  ( 'n_cell', &
                           nX )
      CALL PP % query   ( 'max_grid_size_x', &
                           MaxGridSizeX1 )
      CALL PP % query   ( 'max_grid_size_y', &
                           MaxGridSizeX2 )
      CALL PP % query   ( 'max_grid_size_z', &
                           MaxGridSizeX3 )
      CALL PP % query   ( 'blocking_factor_x', &
                           BlockingFactorX1 )
      CALL PP % query   ( 'blocking_factor_y', &
                           BlockingFactorX2 )
      CALL PP % query   ( 'blocking_factor_z', &
                           BlockingFactorX3 )
      CALL PP % get     ( 'max_level', &
                           MaxLevel )
      IF( MaxLevel .GT. 0 )THEN
        CALL PP % query ( 'UseAMR', &
                           UseAMR )
        CALL PP % query ( 'UseFluxCorrection_Euler', &
                           UseFluxCorrection_Euler )
        CALL PP % query ( 'UseFluxCorrection_TwoMoment', &
                           UseFluxCorrection_TwoMoment )
        CALL PP % getarr( 'TagCriteria', &
                           TagCriteria )
        CALL PP % getarr( 'n_error_buf', &
                           nRefinementBuffer )
      END IF
      CALL PP % getarr  ( 'ref_ratio', &
                           RefinementRatio )
      CALL PP % query   ( 'UseTiling', &
                           UseTiling )
    CALL amrex_parmparse_destroy( PP )

    MaxGridSizeX   = [ MaxGridSizeX1   , MaxGridSizeX2   , MaxGridSizeX3    ]
    BlockingFactor = [ BlockingFactorX1, BlockingFactorX2, BlockingFactorX3 ]

    nMaxLevels = MaxLevel + 1
    nLevels    = nMaxLevels

    CALL InitializeProgramHeader &
           ( ProgramName_Option = TRIM( ProgramName ), &
             nNodes_Option      = nNodes, &
             nX_Option          = nX, &
             nE_Option          = nE, &
             swX_Option         = swX, &
             swE_Option         = swE, &
             xL_Option          = xL, &
             xR_Option          = xR, &
             eL_option          = eL, &
             eR_option          = eR, &
             bcX_Option         = bcX, &
             bcE_Option         = bcE, &
             Verbose_Option     = amrex_parallel_ioprocessor() )

    bcZ = bcZ_TwoMoment

    IF( nDimsX .NE. amrex_spacedim ) &
      CALL DescribeError_MF &
             ( 104, Int_Option = [ nDimsX, amrex_spacedim ] )

    iOS_CPP = 0

    iOS_CPP(1) = 1

    IF( amrex_spacedim .GT. 1 ) iOS_CPP(2) = 1
    IF( amrex_spacedim .GT. 2 ) iOS_CPP(3) = 1

  END SUBROUTINE InitializeParameters


  SUBROUTINE DescribeProgramHeader_AMReX

    CHARACTER(32) :: RFMT, IFMT

    WRITE(RFMT,'(A,I2.2,A)') '(4x,A26,1x,', nMaxLevels, 'ES11.3E3)'
    WRITE(IFMT,'(A,I2.2,A)') '(4x,A26,1x,', nMaxLevels, 'I3.2)'

    IF( .NOT. ALLOCATED( TagCriteria ) )THEN
      ALLOCATE( TagCriteria(nMaxLevels) )
      TagCriteria = 0.0_DP
    END IF

    IF( .NOT. ALLOCATED( nRefinementBuffer ) )THEN
      ALLOCATE( nRefinementBuffer(nMaxLevels) )
      nRefinementBuffer = 1
    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(4x,A)')            'INFO: AMReX'
      WRITE(*,'(4x,A)')            '-----------'
      WRITE(*,'(4x,A26,1x,3I4.3)')   'MaxGridSizeX:', &
                                      MaxGridSizeX
      WRITE(*,'(4x,A26,1x,3I4.3)')   'BlockingFactor:', &
                                      BlockingFactor
      WRITE(*,'(4x,A26,1x,I2.2)')    'nMaxLevels:', &
                                      nMaxLevels
      WRITE(*,'(4x,A26,1x,I2.2)')    'iReGrid:', &
                                      iReGrid
      WRITE(*,'(4x,A26,1x,L)')       'UseFluxCorrection_Euler:', &
                                      UseFluxCorrection_Euler
      WRITE(*,'(4x,A26,1x,L)')       'UseTiling:', &
                                      UseTiling
      WRITE(*,'(4x,A26,1x,L)')       'UseAMR:', &
                                      UseAMR
      WRITE(*,TRIM(RFMT))            'TagCriteria:', &
                                      TagCriteria
      WRITE(*,TRIM(IFMT))            'nRefinementBuffer:', &
                                      nRefinementBuffer
      WRITE(*,*)

    END IF

  END SUBROUTINE DescribeProgramHeader_AMReX


subroutine myamrfinalize
end subroutine myamrfinalize

END MODULE InputParsingModule
