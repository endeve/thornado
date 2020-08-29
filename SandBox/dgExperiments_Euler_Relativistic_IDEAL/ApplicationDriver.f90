PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP,   &
    Zero, &
    One,  &
    Two,  &
    Pi,   &
    TwoPi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE ProgramHeaderModule, ONLY: &
    iX_B0,  &
    iX_B1,  &
    iX_E0,  &
    iX_E1,  &
    nDimsX, &
    nDOFX
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE InitializationModule_Relativistic, ONLY: &
    InitializeFields_Relativistic
  USE Euler_SlopeLimiterModule_Relativistic_IDEAL, ONLY: &
    InitializeSlopeLimiter_Euler_Relativistic_IDEAL, &
    FinalizeSlopeLimiter_Euler_Relativistic_IDEAL,   &
    ApplySlopeLimiter_Euler_Relativistic_IDEAL
  USE Euler_PositivityLimiterModule_Relativistic_IDEAL, ONLY: &
    InitializePositivityLimiter_Euler_Relativistic_IDEAL, &
    FinalizePositivityLimiter_Euler_Relativistic_IDEAL,   &
    ApplyPositivityLimiter_Euler_Relativistic_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeFromConserved_Euler_Relativistic, &
    ComputeTimeStep_Euler_Relativistic
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF,  &
    WriteAccretionShockDiagnosticsHDF
  USE FluidFieldsModule, ONLY: &
    uCF, &
    uPF, &
    uAF, &
    uDF, &
    iPF_D
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GravitySolutionModule_CFA_Poseidon, ONLY: &
    InitializeGravitySolver_CFA_Poseidon, &
    FinalizeGravitySolver_CFA_Poseidon,   &
    SolveGravity_CFA_Poseidon
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK,   &
    UpdateFluid_SSPRK
  USE UnitsModule, ONLY: &
    Kilometer,   &
    SolarMass,   &
    Second,      &
    Millisecond, &
    Centimeter,  &
    Gram,        &
    Erg,         &
    UnitsDisplay
  USE Euler_TallyModule_Relativistic_IDEAL, ONLY: &
    InitializeTally_Euler_Relativistic_IDEAL, &
    FinalizeTally_Euler_Relativistic_IDEAL,   &
    ComputeTally_Euler_Relativistic_IDEAL
  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler,            &
    InitializeTimers_Euler,  &
    FinalizeTimers_Euler,    &
    TimersStart_Euler,       &
    TimersStop_Euler,        &
    Timer_Euler_InputOutput, &
    Timer_Euler_Initialize,  &
    Timer_Euler_Finalize
  USE AccretionShockDiagnosticsModule, ONLY: &
    ComputeAccretionShockDiagnostics
  USE Poseidon_UtilitiesModule, ONLY: &
    ComputeSourceTerms_Poseidon

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: AdvectionProfile
  CHARACTER(32) :: RiemannProblemName
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: wrt
  LOGICAL       :: OPTIMIZE = .FALSE.
  LOGICAL       :: SuppressTally = .FALSE.
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UseCharacteristicLimiting
  LOGICAL       :: UseTroubledCellIndicator
  CHARACTER(4)  :: SlopeLimiterMethod
  LOGICAL       :: UsePositivityLimiter
  LOGICAL       :: SelfGravity
  LOGICAL       :: UseConservativeCorrection
  INTEGER       :: iCycle, iCycleD, iCycleW
  INTEGER       :: nX(3), bcX(3), swX(3), nNodes
  INTEGER       :: nStagesSSPRK
  INTEGER       :: RestartFileNumber
  REAL(DP)      :: SlopeTolerance
  REAL(DP)      :: Min_1, Min_2
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, CFL
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: LimiterThresholdParameter
  REAL(DP)      :: Mass = Zero
  REAL(DP)      :: ZoomX(3)

  ! --- Sedov--Taylor blast wave ---
  REAL(DP) :: Eblast
  INTEGER  :: nDetCells
  REAL(DP) :: Vmax, LorentzFactor

  ! --- Standing accretion shock ---
  REAL(DP) :: MassPNS, RadiusPNS, ShockRadius, &
              AccretionRate, PolytropicConstant
  LOGICAL  :: ApplyPerturbation
  INTEGER  :: PerturbationOrder
  REAL(DP) :: PerturbationAmplitude, &
              rPerturbationInner, rPerturbationOuter
  REAL(DP) :: Power(0:2)

  ! --- Yahil Collapse ---
  REAL(DP) :: CentralDensity, CentralPressure, CoreRadius, CollapseTime

  LOGICAL  :: WriteGF = .TRUE., WriteFF = .TRUE.
  LOGICAL  :: ActivateUnits = .FALSE.
  REAL(DP) :: Timer_Evolution

  REAL(DP), ALLOCATABLE :: U_Poseidon(:,:,:,:,:)

  TimeIt_Euler = .TRUE.
  CALL InitializeTimers_Euler
  CALL TimersStart_Euler( Timer_Euler_Initialize )

  ProgramName = 'Advection'
!  ProgramName = 'Advection2D'
!  ProgramName = 'RiemannProblem'
!  ProgramName = 'RiemannProblem2D'
!  ProgramName = 'RiemannProblemSpherical'
!  ProgramName = 'SedovTaylorBlastWave'
!  ProgramName = 'KelvinHelmholtzInstability'
!  ProgramName = 'StandingAccretionShock'
!  ProgramName = 'StaticTOV'
!  ProgramName = 'YahilCollapse'

  swX               = [ 0, 0, 0 ]
  RestartFileNumber = -1
  t                 = 0.0_DP
  ZoomX             = 1.0_DP
  SelfGravity       = .FALSE.

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'Advection' )

      AdvectionProfile = 'SineWave'

      Gamma = 5.0_DP / 3.0_DP
      t_end = 10.0_DP
      bcX = [ 1, 0, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 64, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'Advection2D' )

      AdvectionProfile = 'SineWaveX1X2'

      Gamma = 5.0_DP / 3.0_DP
      t_end = 10.0_DP
      bcX = [ 1, 1, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 32, 32, 1 ]
      swX = [ 1, 1, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ One / SQRT( Two ), One / SQRT( Two ), 1.0_DP ]

    CASE( 'RiemannProblem' )

      RiemannProblemName = 'Sod'

      SELECT CASE ( TRIM( RiemannProblemName ) )

        CASE( 'Sod' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.2d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'IsolatedShock' )

          Gamma = 4.0_DP / 3.0_DP
          t_end = 2.0e1_DP
          bcX   = [ 2, 0, 0 ]

        CASE( 'IsolatedContact' )
          Gamma = 4.0_DP / 3.0_DP
          t_end = 2.0e1_DP
          bcX   = [ 2, 0, 0 ]

        CASE( 'MBProblem1' )
          Gamma = 4.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'MBProblem4' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'PerturbedShockTube' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.35d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'ShockReflection' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.75d0
          bcX   = [ 23, 0, 0 ]

        CASE DEFAULT

          WRITE(*,*)
          WRITE(*,'(A21,A)') 'Invalid RiemannProblemName: ', RiemannProblemName
          WRITE(*,'(A)')     'Valid choices:'
          WRITE(*,'(A)')     '  Sod'
          WRITE(*,'(A)')     '  IsolatedShock'
          WRITE(*,'(A)')     '  IsolatedContact'
          WRITE(*,'(A)')     '  MBProblem1'
          WRITE(*,'(A)')     '  MBProblem4'
          WRITE(*,'(A)')     '  PerturbedShockTube'
          WRITE(*,'(A)')     '  ShockReflection'
          WRITE(*,'(A)')     'Stopping...'
          STOP

      END SELECT

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 128, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'RiemannProblem2D' )

      RiemannProblemName = 'IsolatedShock'

      SELECT CASE ( TRIM( RiemannProblemName ) )

        CASE( 'DzB2002' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 2, 0 ]

        CASE( 'IsolatedShock' )

          Gamma = 4.0_DP / 3.0_DP
          t_end = 25.0_DP
          bcX   = [ 2, 2, 0 ]

      END SELECT

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 64, 64, 1 ]
      swX = [ 1, 1, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'RiemannProblemSpherical' )

      RiemannProblemName = 'SphericalSod'

      CoordinateSystem = 'SPHERICAL'

      Gamma = 5.0_DP / 3.0_DP
      t_end = 5.0d-1
      bcX = [ 2, 0, 0 ]

      nX  = [ 256, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 2.0_DP, Pi, TwoPi ]

      WriteGF = .TRUE.

    CASE( 'SedovTaylorBlastWave' )

      nDetCells = 1
      Eblast    = 1.0d-3

      CoordinateSystem = 'SPHERICAL'

      Gamma = 4.0_DP / 3.0_DP
      t_end = 1.0d0
      bcX = [ 3, 0, 0 ]

      nX  = [ 256, 1, 1 ]
      swX = [ 1, 0, 0 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.2_DP, Pi, TwoPi ]

      WriteGF = .TRUE.

    CASE( 'KelvinHelmholtzInstability' )

       CoordinateSystem = 'CARTESIAN'

       Gamma = 4.0d0 / 3.0d0
       t_end = 1.0d-1
       bcX = [ 1, 1, 0 ]

       nX = [ 16, 32, 1 ]
      swX = [ 1, 1, 0 ]
       xL = [ -0.5d0, -1.0d0, 0.0d0 ]
       xR = [  0.5d0,  1.0d0, 1.0d0 ]

    CASE( 'StandingAccretionShock' )

      Gamma = 4.0e0_DP / 3.0e0_DP
      t_end = 3.0d2 * Millisecond
      bcX = [ 11, 0, 0 ]

      MassPNS            = 1.4_DP    * SolarMass
      RadiusPNS          = 40.0_DP   * Kilometer
      ShockRadius        = 180.0_DP  * Kilometer
      AccretionRate      = 0.3_DP    * ( SolarMass / Second )
      PolytropicConstant = 2.0e14_DP * ( Erg / Centimeter**3 &
                                         / ( Gram / Centimeter**3 )**( Gamma ) )
      ApplyPerturbation     = .TRUE.
      PerturbationOrder     = 0
      PerturbationAmplitude = 0.04_DP
      rPerturbationInner    = 260.0_DP * Kilometer
      rPerturbationOuter    = 280.0_DP * Kilometer

      nX  = [ 960, 1, 1 ]
      swX = [ 1, 1, 0 ]
      xL  = [ RadiusPNS, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0e3_DP * Kilometer, Pi, TwoPi ]

      CoordinateSystem = 'SPHERICAL'

      WriteGF = .TRUE.

      ActivateUnits = .TRUE.

      Mass = MassPNS

    CASE( 'StaticTOV' )

       SelfGravity = .TRUE.

       CoordinateSystem = 'SPHERICAL'

       Gamma = 2.0_DP
       t_end = 1.0e1_DP * Millisecond
       bcX   = [ 30, 0, 0 ]

       nX = [ 128                 , 1   , 1     ]
      swX = [ 1                   , 0   , 0     ]
       xL = [ Zero                , Zero, Zero  ]
       xR = [ 1.0e1_DP * Kilometer,  Pi , TwoPi ]

      WriteGF = .TRUE.

      ActivateUnits = .TRUE.

    CASE( 'YahilCollapse' )

      SelfGravity = .TRUE.

      CoordinateSystem = 'SPHERICAL'

      CentralDensity  = 7.0e9_DP  * ( Gram / Centimeter**3 )
      CentralPressure = 6.0e27_DP * ( Erg  / Centimeter**3 )
      CoreRadius      = 1.0e5_DP  * Kilometer
      CollapseTime    = 1.50e2_DP * Millisecond

      Gamma = 1.30_DP
      t_end = CollapseTime - 0.5_DP * Millisecond
      bcX = [ 30, 0, 0 ]

      nX    = [ 256                 , 1     , 1      ]
      swX   = [ 1                   , 0     , 0      ]
      xL    = [ Zero                , Zero  , Zero   ]
      xR    = [ CoreRadius          , Pi    , TwoPi  ]
      ZoomX = [ 1.032034864238313_DP, 1.0_DP, 1.0_DP ]

      WriteGF = .TRUE.

      ActivateUnits = .TRUE.

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
      WRITE(*,'(A)')     'Valid choices:'
      WRITE(*,'(A)')     '  Advection'
      WRITE(*,'(A)')     '  Advection2D'
      WRITE(*,'(A)')     '  RiemannProblem'
      WRITE(*,'(A)')     '  RiemannProblem2D'
      WRITE(*,'(A)')     '  RiemannProblemSpherical'
      WRITE(*,'(A)')     '  SedovTaylorBlastWave'
      WRITE(*,'(A)')     '  KelvinHelmholtzInstability'
      WRITE(*,'(A)')     '  StandingAccretionShock'
      WRITE(*,'(A)')     '  StaticTOV'
      WRITE(*,'(A)')     '  YahilCollapse'
      WRITE(*,'(A)')     'Stopping...'
      STOP

  END SELECT

  ! --- DG ---

  nNodes = 1
  IF( .NOT. nNodes .LE. 4 ) &
    STOP 'nNodes must be less than or equal to four.'

  ! --- Time Stepping ---

  nStagesSSPRK = 1
  IF( .NOT. nStagesSSPRK .LE. 3 ) &
    STOP 'nStagesSSPRK must be less than or equal to three.'

  CFL = 0.5_DP ! Cockburn & Shu, (2001), JSC, 16, 173

  ! --- Slope Limiter ---

  UseSlopeLimiter           = .TRUE.
  SlopeLimiterMethod        = 'TVD'
  BetaTVD                   = 1.75d0
  BetaTVB                   = 0.0d0
  SlopeTolerance            = 1.0d-6
  UseCharacteristicLimiting = .TRUE.
  UseTroubledCellIndicator  = .TRUE.
  LimiterThresholdParameter = 0.015_DP
  UseConservativeCorrection = .TRUE.

  ! --- Positivity Limiter ---

  UsePositivityLimiter = .TRUE.
  Min_1                = 1.0d-13
  Min_2                = 1.0d-13

  ! === End of User Input ===

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           nX_Option &
             = nX, &
           swX_Option &
             = swX, &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           ZoomX_Option &
             = ZoomX, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = TRIM( CoordinateSystem ), &
           ActivateUnits_Option &
             = ActivateUnits, &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  IF( SelfGravity )THEN

    ALLOCATE( U_Poseidon(1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3),1:6) )

    CALL InitializeGravitySolver_CFA_Poseidon

  ELSE

    CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass_Option = Mass )

  END IF

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeSlopeLimiter_Euler_Relativistic_IDEAL &
         ( UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           SlopeLimiterMethod_Option &
             = TRIM( SlopeLimiterMethod ), &
           BetaTVD_Option &
             = BetaTVD, &
           BetaTVB_Option &
             = BetaTVB, &
           SlopeTolerance_Option &
             = SlopeTolerance, &
           UseCharacteristicLimiting_Option &
             = UseCharacteristicLimiting, &
           UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator, &
           LimiterThresholdParameter_Option &
             = LimiterThresholdParameter, &
           UseConservativeCorrection_Option &
             = UseConservativeCorrection )

  CALL InitializePositivityLimiter_Euler_Relativistic_IDEAL &
         ( UsePositivityLimiter_Option = UsePositivityLimiter, &
           Verbose_Option = .TRUE., &
           Min_1_Option = Min_1, &
           Min_2_Option = Min_2 )

  CALL InitializeFluid_SSPRK( nStages = nStagesSSPRK )
  WRITE(*,*)
  WRITE(*,'(A6,A,ES11.3E3)') '', 'CFL: ', CFL

  CALL InitializeFields_Relativistic &
         ( AdvectionProfile_Option &
             = TRIM( AdvectionProfile ), &
           RiemannProblemName_Option &
             = TRIM( RiemannProblemName ), &
           nDetCells_Option             = nDetCells, &
           Eblast_Option                = Eblast, &
           MassPNS_Option               = MassPNS, &
           ShockRadius_Option           = ShockRadius, &
           AccretionRate_Option         = AccretionRate, &
           PolytropicConstant_Option    = PolytropicConstant, &
           ApplyPerturbation_Option     = ApplyPerturbation, &
           PerturbationOrder_Option     = PerturbationOrder, &
           PerturbationAmplitude_Option = PerturbationAmplitude, &
           rPerturbationInner_Option    = rPerturbationInner, &
           rPerturbationOuter_Option    = rPerturbationOuter, &
           CentralDensity_Option        = CentralDensity, &
           CentralPressure_Option       = CentralPressure, &
           CoreRadius_Option            = CoreRadius, &
           CollapseTime_Option          = CollapseTime )

  IF( RestartFileNumber .GE. 0 )THEN

    CALL ReadFieldsHDF &
           ( RestartFileNumber, t, &
             ReadFF_Option = .TRUE., ReadGF_Option = .TRUE. )

  END IF

  iCycleD = 10
!!$  iCycleW = 1; dt_wrt = -1.0d0
  dt_wrt = 1.0d-2 * ( t_end - t ); iCycleW = -1

  IF( dt_wrt .GT. Zero .AND. iCycleW .GT. 0 ) &
    STOP 'dt_wrt and iCycleW cannot both be present'

  CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

  IF( SelfGravity )THEN

    CALL ComputeSourceTerms_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, U_Poseidon )

    CALL SolveGravity_CFA_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, U_Poseidon )

  END IF

  CALL TimersStop_Euler( Timer_Euler_Initialize )

  IF( .NOT. OPTIMIZE .AND. RestartFileNumber .LT. 0 )THEN

    CALL TimersStart_Euler( Timer_Euler_InputOutput )

    CALL ComputeFromConserved_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

    CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

    CALL TimersStop_Euler( Timer_Euler_InputOutput )

  END IF

  CALL TimersStart_Euler( Timer_Euler_Initialize )

  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'Begin evolution'
  WRITE(*,'(A2,A)') '', '---------------'
  WRITE(*,*)

  t_wrt = t + dt_wrt
  wrt   = .FALSE.

  CALL InitializeTally_Euler_Relativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF, uCF, SuppressTally_Option = SuppressTally )

  CALL TimersStop_Euler( Timer_Euler_Initialize )

  iCycle = 0
  Timer_Evolution = MPI_WTIME()
  DO WHILE( t .LT. t_end )

    iCycle = iCycle + 1

    CALL ComputeTimeStep_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             uGF, uCF, &
             CFL / ( nDimsX * ( Two * DBLE( nNodes ) - One ) ), &
             dt )

    IF( t + dt .LT. t_end )THEN

      t = t + dt

    ELSE

      dt = t_end - t
      t  = t_end

    END IF

    CALL TimersStart_Euler( Timer_Euler_InputOutput )

    IF( MOD( iCycle, iCycleD ) .EQ. 0 )THEN

      WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
        'Cycle: ', iCycle, ' t = ', t / UnitsDisplay % TimeUnit, &
        TRIM( UnitsDisplay % TimeLabel ), &
        ' dt = ', dt /  UnitsDisplay % TimeUnit, &
        TRIM( UnitsDisplay % TimeLabel )

    END IF

    CALL TimersStop_Euler( Timer_Euler_InputOutput )

    IF( SelfGravity )THEN

      CALL UpdateFluid_SSPRK &
             ( t, dt, uGF, uCF, uDF, &
               ComputeIncrement_Euler_DG_Explicit, &
               SolveGravity_CFA_Poseidon )

    ELSE

      CALL UpdateFluid_SSPRK &
             ( t, dt, uGF, uCF, uDF, &
               ComputeIncrement_Euler_DG_Explicit )

    END IF

    IF( .NOT. OPTIMIZE )THEN

      CALL TimersStart_Euler( Timer_Euler_InputOutput )

      IF( iCycleW .GT. 0 )THEN

        IF( MOD( iCycle, iCycleW ) .EQ. 0 ) &
          wrt = .TRUE.

      ELSE

        IF( t + dt .GT. t_wrt )THEN

          t_wrt = t_wrt + dt_wrt
          wrt   = .TRUE.

        END IF

      END IF

      CALL TimersStop_Euler( Timer_Euler_InputOutput )

      IF( wrt )THEN

        CALL TimersStart_Euler( Timer_Euler_InputOutput )

        CALL ComputeFromConserved_Euler_Relativistic &
               ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

        CALL WriteFieldsHDF &
               ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

        CALL TimersStop_Euler( Timer_Euler_InputOutput )

        CALL ComputeTally_Euler_Relativistic_IDEAL &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )

        wrt = .FALSE.

      END IF
    END IF

    IF( TRIM( ProgramName ) .EQ. 'StandingAccretionShock' )THEN

      CALL ComputeFromConserved_Euler_Relativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      CALL ComputeAccretionShockDiagnostics &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, Power )

      CALL WriteAccretionShockDiagnosticsHDF( t, Power )

    END IF

    IF( TRIM( ProgramName ) == 'YahilCollapse' )THEN

      CALL ComputeFromConserved_Euler_Relativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      IF( ANY( uPF(:,:,:,:,iPF_D) .GT. 1.0e15_DP * Gram / Centimeter**3 ) ) EXIT

    END IF

  END DO

  Timer_Evolution = MPI_WTIME() - Timer_Evolution
  WRITE(*,*)
  WRITE(*,'(A,ES13.6E3,A)') 'Total evolution time: ', Timer_Evolution, ' s'
  WRITE(*,*)

  IF( .NOT. OPTIMIZE )THEN

    CALL TimersStart_Euler( Timer_Euler_InputOutput )

    CALL ComputeFromConserved_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

    IF( SelfGravity )THEN

      CALL ComputeSourceTerms_Poseidon &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, U_Poseidon )

      CALL SolveGravity_CFA_Poseidon &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, U_Poseidon )

    END IF

    CALL WriteFieldsHDF &
           ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

    CALL TimersStop_Euler( Timer_Euler_InputOutput )

  END IF

  CALL ComputeTally_Euler_Relativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )

  CALL TimersStart_Euler( Timer_Euler_Finalize )

  CALL FinalizeTally_Euler_Relativistic_IDEAL

  CALL FinalizePositivityLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  IF( SelfGravity )THEN

    DEALLOCATE( U_Poseidon )

    CALL FinalizeGravitySolver_CFA_Poseidon

  END IF

  CALL FinalizeEquationOfState

  CALL FinalizeProgram

  CALL TimersStop_Euler( Timer_Euler_Finalize )

  CALL FinalizeTimers_Euler

  WRITE(*,*)
  WRITE(*,'(2x,A)') 'git info'
  WRITE(*,'(2x,A)') '--------'
  WRITE(*,*)
  WRITE(*,'(2x,A)') 'git branch:'
  CALL EXECUTE_COMMAND_LINE( 'git branch' )
  WRITE(*,*)
  WRITE(*,'(2x,A)') 'git describe --tags:'
  CALL EXECUTE_COMMAND_LINE( 'git describe --tags' )
  WRITE(*,*)
  WRITE(*,'(2x,A)') 'git rev-parse HEAD:'
  CALL EXECUTE_COMMAND_LINE( 'git rev-parse HEAD' )
  WRITE(*,*)
  WRITE(*,'(2x,A)') 'date:'
  CALL EXECUTE_COMMAND_LINE( 'date' )
  WRITE(*,*)

END PROGRAM ApplicationDriver
