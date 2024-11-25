PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two, &
    Pi, &
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
    iX_B0, &
    iX_B1, &
    iX_E0, &
    iX_E1, &
    nDimsX
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE InitializationModule_Relativistic, ONLY: &
    InitializeFields_Relativistic
  USE Euler_SlopeLimiterModule_Relativistic_IDEAL, ONLY: &
    InitializeSlopeLimiter_Euler_Relativistic_IDEAL, &
    FinalizeSlopeLimiter_Euler_Relativistic_IDEAL, &
    ApplySlopeLimiter_Euler_Relativistic_IDEAL
  USE Euler_PositivityLimiterModule_Relativistic_IDEAL, ONLY: &
    InitializePositivityLimiter_Euler_Relativistic_IDEAL, &
    FinalizePositivityLimiter_Euler_Relativistic_IDEAL, &
    ApplyPositivityLimiter_Euler_Relativistic_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeFromConserved_Euler_Relativistic, &
    ComputeTimeStep_Euler_Relativistic
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF, &
    WriteAccretionShockDiagnosticsHDF
  USE FluidFieldsModule, ONLY: &
    uCF, &
    uPF, &
    uAF, &
    uDF
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK,   &
    UpdateFluid_SSPRK
  USE UnitsModule, ONLY: &
    Kilometer, &
    SolarMass, &
    Second, &
    Millisecond, &
    Centimeter, &
    Gram, &
    Erg, &
    UnitsDisplay
  USE Euler_TallyModule_Relativistic, ONLY: &
    InitializeTally_Euler_Relativistic, &
    FinalizeTally_Euler_Relativistic, &
    ComputeTally_Euler_Relativistic
  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler, &
    InitializeTimers_Euler, &
    FinalizeTimers_Euler, &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_InputOutput, &
    Timer_Euler_Initialize, &
    Timer_Euler_Finalize
  USE AccretionShockUtilitiesModule, ONLY: &
    WriteInitialConditionsToFile, &
    ComputeAccretionShockDiagnostics

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: wrt
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UseCharacteristicLimiting
  LOGICAL       :: UseTroubledCellIndicator
  LOGICAL       :: SuppressTally
  CHARACTER(4)  :: SlopeLimiterMethod
  LOGICAL       :: UsePositivityLimiter
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
  LOGICAL       :: WriteGF = .TRUE., WriteFF = .TRUE.
  REAL(DP)      :: Timer_Evolution

  ! --- Standing accretion shock ---

  REAL(DP) :: MassPNS, RadiusPNS, ShockRadius, &
              AccretionRate, PolytropicConstant
  LOGICAL  :: ApplyPerturbation
  INTEGER  :: PerturbationOrder
  REAL(DP) :: PerturbationAmplitude, &
              rPerturbationInner, rPerturbationOuter
  REAL(DP) :: Power(0:2)

  LOGICAL  :: ActivateUnits = .TRUE.
  LOGICAL  :: WriteAccretionShockDiagnostics = .TRUE.

  CHARACTER(32) :: Parameters

  ! --- Relaxation Parameters ---

  LOGICAL            :: InitializeFluidFromFile          = .FALSE.
  LOGICAL            :: WriteInitialConditionsToFile_SAS = .FALSE.
  CHARACTER(LEN=128) :: InitialConditionsFileName

  WRITE(Parameters,'(A)') 'MPNS1.4_Mdot0.3_Rs180'

  WRITE(InitialConditionsFileName,'(A)') &
    'InitialConditions_' // TRIM( Parameters )

  SuppressTally = .FALSE.

  TimeIt_Euler = .TRUE.
  CALL InitializeTimers_Euler
  CALL TimersStart_Euler( Timer_Euler_Initialize )

  ProgramName = 'StandingAccretionShock'

  RestartFileNumber = -1
  t                 = 0.0_DP
  ZoomX             = 1.0_DP

  Gamma = 4.0e0_DP / 3.0e0_DP
  t_end = 5.0e1_DP * Millisecond
  bcX = [ 100, 3, 0 ]

  MassPNS            = 1.4_DP    * SolarMass
  RadiusPNS          = 40.0_DP   * Kilometer
  ShockRadius        = 180.0_DP  * Kilometer
  AccretionRate      = 0.3_DP    * ( SolarMass / Second )
  PolytropicConstant = 2.0e14_DP * ( Erg / Centimeter**3 &
                                     / ( Gram / Centimeter**3 )**( Gamma ) )
  ApplyPerturbation     = .FALSE.
  PerturbationOrder     = 0
  PerturbationAmplitude = 0.04_DP
  rPerturbationInner    = 260.0_DP * Kilometer
  rPerturbationOuter    = 280.0_DP * Kilometer

  nX  = [ 960, 32, 1 ]
  swX = [ 1, 1, 0 ]
  xL  = [ RadiusPNS, 0.0_DP, 0.0_DP ]
  xR  = [ 1.0e3_DP * Kilometer, Pi, TwoPi ]

  CoordinateSystem = 'SPHERICAL'

  Mass = MassPNS

  ! --- DG ---

  nNodes = 3
  IF( .NOT. nNodes .LE. 4 ) &
    STOP 'nNodes must be less than or equal to four.'

  ! --- Time Stepping ---

  nStagesSSPRK = 3
  IF( .NOT. nStagesSSPRK .LE. 3 ) &
    STOP 'nStagesSSPRK must be less than or equal to three.'

  CFL = 0.5_DP ! Cockburn & Shu, (2001), JSC, 16, 173

  ! --- Slope Limiter ---

  UseSlopeLimiter           = .TRUE.
  SlopeLimiterMethod        = 'TVD'
  BetaTVD                   = 1.75d0
  BetaTVB                   = 0.0e0_DP
  SlopeTolerance            = 1.0e-6_DP
  UseCharacteristicLimiting = .TRUE.
  UseTroubledCellIndicator  = .TRUE.
  LimiterThresholdParameter = 0.015_DP
  UseConservativeCorrection = .TRUE.

  ! --- Positivity Limiter ---

  UsePositivityLimiter = .TRUE.
  Min_1                = 1.0e-13_DP
  Min_2                = 1.0e-13_DP

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

  CALL ComputeGeometryX &
       ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass_Option = Mass )

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

  uCF = Zero ! Without this, crashes when copying data in TimeStepper
  uDF = Zero ! Without this, crashes in IO

  CALL InitializeFields_Relativistic &
         ( MassPNS_Option                   = MassPNS, &
           ShockRadius_Option               = ShockRadius, &
           AccretionRate_Option             = AccretionRate, &
           PolytropicConstant_Option        = PolytropicConstant, &
           ApplyPerturbation_Option         = ApplyPerturbation, &
           PerturbationOrder_Option         = PerturbationOrder, &
           PerturbationAmplitude_Option     = PerturbationAmplitude, &
           rPerturbationInner_Option        = rPerturbationInner, &
           rPerturbationOuter_Option        = rPerturbationOuter, &
           InitializeFluidFromFile_Option   = InitializeFluidFromFile, &
           InitialConditionsFileName_Option = InitialConditionsFileName )

  IF( RestartFileNumber .LT. 0 )THEN

    CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

    CALL ComputeFromConserved_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( uCF, uGF )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST       ( uCF, uGF )
#endif

    CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  ELSE

    CALL ReadFieldsHDF &
           ( RestartFileNumber, t, &
             ReadFF_Option = .TRUE., ReadGF_Option = .TRUE. )

  END IF

  iCycleD = 10
!!$  iCycleW = 1; dt_wrt = -1.0d0
!!$  dt_wrt = 1.0d-2 * ( t_end - t ); iCycleW = -1
  dt_wrt = 1.0_DP * Millisecond; iCycleW = -1

  IF( dt_wrt .GT. Zero .AND. iCycleW .GT. 0 ) &
    STOP 'dt_wrt and iCycleW cannot both be present'

  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'Begin evolution'
  WRITE(*,'(A2,A)') '', '---------------'
  WRITE(*,*)

  t_wrt = t + dt_wrt
  wrt   = .FALSE.

  CALL InitializeTally_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, &
           SuppressTally_Option = SuppressTally )

  CALL ComputeTally_Euler_Relativistic &
       ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
         SetInitialValues_Option = .TRUE., Verbose_Option = .FALSE. )

  CALL TimersStop_Euler( Timer_Euler_Initialize )

  IF( nDimsX .GT. 1 .AND. WriteAccretionShockDiagnostics )THEN

    OPEN( 100, FILE = TRIM( Parameters ) &
                        // '_PowerInLegendreModes.dat' )
    WRITE(100,'(A)') 'Time [ms], Power(0), Power(1), Power(2)'

  END IF

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

    IF( MOD( iCycle, iCycleD ) .EQ. 0 )THEN

      WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
        'Cycle: ', iCycle, ' t = ', t / UnitsDisplay % TimeUnit, &
        TRIM( UnitsDisplay % TimeLabel ), &
        ' dt = ', dt /  UnitsDisplay % TimeUnit, &
        TRIM( UnitsDisplay % TimeLabel )

    END IF

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, uDF, &
             ComputeIncrement_Euler_DG_Explicit )

    IF( iCycleW .GT. 0 )THEN

      IF( MOD( iCycle, iCycleW ) .EQ. 0 ) &
        wrt = .TRUE.

    ELSE

      IF( t + dt .GT. t_wrt )THEN

        t_wrt = t_wrt + dt_wrt
        wrt   = .TRUE.

      END IF

    END IF

    IF( wrt )THEN

      CALL TimersStart_Euler( Timer_Euler_InputOutput )

      CALL ComputeFromConserved_Euler_Relativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( uCF, uGF )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST       ( uCF, uGF )
#endif

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

      CALL ComputeTally_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
             Verbose_Option = .FALSE. )

      wrt = .FALSE.

      CALL TimersStop_Euler( Timer_Euler_InputOutput )

    END IF

    IF( nDimsX .GT. 1 .AND. WriteAccretionShockDiagnostics )THEN

      CALL ComputeFromConserved_Euler_Relativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      CALL ComputeAccretionShockDiagnostics &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, Power )

      WRITE(100,'(ES24.16E3,1x,ES24.16E3,1x,ES24.16E3,1x,ES24.16E3)') &
        t / Millisecond, Power(0), Power(1), Power(2)
!!$      CALL WriteAccretionShockDiagnosticsHDF( t, Power )

    END IF

  END DO

  IF( nDimsX .GT. 1 .AND. WriteAccretionShockDiagnostics ) &
    CLOSE( 100 )

  Timer_Evolution = MPI_WTIME() - Timer_Evolution
  WRITE(*,*)
  WRITE(*,'(A,ES13.6E3,A)') 'Total evolution time: ', Timer_Evolution, ' s'
  WRITE(*,*)

  CALL TimersStart_Euler( Timer_Euler_Finalize )

  CALL ComputeFromConserved_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( uCF, uGF )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST       ( uCF, uGF )
#endif

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  CALL ComputeTally_Euler_Relativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t )

  IF( WriteInitialConditionsToFile_SAS )THEN

    CALL WriteInitialConditionsToFile &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, InitialConditionsFileName )

  END IF

  CALL FinalizeTally_Euler_Relativistic

  CALL FinalizePositivityLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

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
