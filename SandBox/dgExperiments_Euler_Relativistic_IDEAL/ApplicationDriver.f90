PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Zero, One, Two, Pi, TwoPi
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
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nDimsX, nDOFX
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
    ReadFieldsHDF
  USE FluidFieldsModule, ONLY: &
    nCF, nPF, nAF, &
    uCF, uPF, uAF, &
    uDF
  USE GeometryFieldsModule, ONLY: &
    nGF, uGF
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE UnitsModule, ONLY: &
    Kilometer, SolarMass, Second, Millisecond
  USE Euler_TallyModule_Relativistic_IDEAL, ONLY: &
    InitializeTally_Euler_Relativistic_IDEAL, &
    FinalizeTally_Euler_Relativistic_IDEAL, &
    ComputeTally_Euler_Relativistic_IDEAL
  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler, &
    InitializeTimers_Euler, FinalizeTimers_Euler, &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_InputOutput, &
    Timer_Euler_Initialize, &
    Timer_Euler_Finalize

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: AdvectionProfile
  CHARACTER(32) :: RiemannProblemName
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: wrt
  LOGICAL       :: OPTIMIZE = .FALSE.
  LOGICAL       :: SuppressTally = .TRUE.
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UseCharacteristicLimiting
  LOGICAL       :: UseTroubledCellIndicator
  CHARACTER(4)  :: SlopeLimiterMethod
  LOGICAL       :: UsePositivityLimiter
  LOGICAL       :: UseConservativeCorrection
  INTEGER       :: iCycle, iCycleD, iCycleW = 0
  INTEGER       :: nX(3), bcX(3), swX(3), nNodes
  INTEGER       :: nStagesSSPRK
  INTEGER       :: RestartFileNumber
  REAL(DP)      :: SlopeTolerance
  REAL(DP)      :: Min_1, Min_2
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, CFL
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: LimiterThresholdParameter

  ! --- Sedov--Taylor blast wave ---
  REAL(DP) :: Eblast
  INTEGER  :: nDetCells
  REAL(DP) :: Vmax, LorentzFactor

  ! --- Standing accretion shock ---
  REAL(DP), ALLOCATABLE :: FluidFieldParameters(:)
  REAL(DP)              :: MassPNS, RadiusPNS, ShockRadius, &
                           AccretionRate, MachNumber
  LOGICAL               :: ApplyPerturbation
  INTEGER               :: PerturbationOrder
  REAL(DP)              :: PerturbationAmplitude, &
                           rPerturbationInner, rPerturbationOuter

  LOGICAL  :: WriteGF = .FALSE., WriteFF = .TRUE.
  LOGICAL  :: ActivateUnits = .FALSE.
  REAL(DP) :: Timer_Evolution

  RestartFileNumber = -1

  t = 0.0_DP

  TimeIt_Euler = .TRUE.
  CALL InitializeTimers_Euler
  CALL TimersStart_Euler( Timer_Euler_Initialize )

!  ProgramName = 'Advection'
!  ProgramName = 'Advection2D'
  ProgramName = 'RiemannProblem'
!  ProgramName = 'RiemannProblem2D'
!  ProgramName = 'RiemannProblemSpherical'
!  ProgramName = 'SedovTaylorBlastWave'
!  ProgramName = 'KelvinHelmholtzInstability'
!  ProgramName = 'StandingAccretionShock'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'Advection' )

      AdvectionProfile = 'TopHat'

      Gamma = 5.0_DP / 3.0_DP
      t_end = 10.0_DP
      bcX = [ 1, 0, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX = [ 64, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'Advection2D' )

      AdvectionProfile = 'SineWaveX1'

      Gamma = 5.0_DP / 3.0_DP
      t_end = 10.0_DP
      bcX = [ 1, 1, 0 ]

      CoordinateSystem = 'CARTESIAN'

      nX = [ 32, 32, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

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
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'RiemannProblemSpherical' )

      RiemannProblemName = 'SphericalSod'

      CoordinateSystem = 'SPHERICAL'

      Gamma = 5.0_DP / 3.0_DP

      nX = [ 128, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 2.0_DP, Pi, TwoPi ]

      bcX = [ 2, 0, 0 ]

      t_end = 5.0d-1

      WriteGF = .TRUE.

    CASE( 'SedovTaylorBlastWave' )

      nDetCells = 1
      Eblast    = 1.0d-3

      CoordinateSystem = 'SPHERICAL'

      Gamma = 4.0_DP / 3.0_DP

      nX = [ 256, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.2_DP, Pi, TwoPi ]

      bcX = [ 3, 0, 0 ]

      t_end = 1.0d0

      WriteGF = .TRUE.

    CASE( 'KelvinHelmholtzInstability' )

       CoordinateSystem = 'CARTESIAN'

       Gamma = 4.0d0 / 3.0d0

       nX = [ 16, 32, 1 ]
       xL = [ -0.5d0, -1.0d0, 0.0d0 ]
       xR = [  0.5d0,  1.0d0, 1.0d0 ]

       bcX = [ 1, 1, 0 ]

       t_end = 1.0d-1

    CASE( 'StandingAccretionShock' )

      CoordinateSystem = 'SPHERICAL'

      MassPNS       = 1.4_DP * SolarMass
      RadiusPNS     = 40.0_DP * Kilometer
      ShockRadius   = 180.0_DP * Kilometer
      AccretionRate = 0.3_DP * SolarMass / Second
      MachNumber    = 10.0_DP

      ApplyPerturbation     = .TRUE.
      PerturbationOrder     = 1
      PerturbationAmplitude = 0.1_DP
      rPerturbationInner    = 260.0_DP * Kilometer
      rPerturbationOuter    = 280.0_DP * Kilometer

      Gamma = 4.0d0 / 3.0d0

      nX = [ 128, 16, 1 ]
      xL = [ RadiusPNS, 0.0_DP, 0.0_DP ]
      xR = [ Two * ShockRadius, Pi, TwoPi ]

      bcX = [ 11, 0, 0 ]

      t_end = 3.0d2 * Millisecond

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
      WRITE(*,'(A)')     'Stopping...'
      STOP

  END SELECT

  nNodes = 3
  IF( .NOT. nNodes .LE. 4 ) &
    STOP 'nNodes must be less than or equal to four.'

  BetaTVD = 1.75d0
  BetaTVB = 0.0d0

  UseSlopeLimiter           = .TRUE.
  SlopeTolerance            = 1.0d-6
  UseCharacteristicLimiting = .TRUE.

  UseTroubledCellIndicator  = .TRUE.

  SlopeLimiterMethod        = 'TVD'

  LimiterThresholdParameter = 0.01_DP

  UseConservativeCorrection = .TRUE.

  UsePositivityLimiter = .TRUE.
  Min_1 = 1.0d-13
  Min_2 = 1.0d-13

  nStagesSSPRK = nNodes
  IF( .NOT. nStagesSSPRK .LE. 3 ) &
    STOP 'nStagesSSPRK must be less than or equal to three.'

  ! --- Cockburn & Shu, (2001), JSC, 16, 173 ---
  CFL = 0.5_DP

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
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
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass_Option = MassPNS )

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
           MachNumber_Option            = MachNumber, &
           ApplyPerturbation_Option     = ApplyPerturbation, &
           PerturbationOrder_Option     = PerturbationOrder, &
           PerturbationAmplitude_Option = PerturbationAmplitude, &
           rPerturbationInner_Option    = rPerturbationInner, &
           rPerturbationOuter_Option    = rPerturbationOuter )

  IF( RestartFileNumber .GE. 0 )THEN

    CALL ReadFieldsHDF( RestartFileNumber, t, ReadFF_Option = .TRUE. )

  END IF

  iCycleD = 10
!!$  iCycleW = 10; dt_wrt = -1.0d0
  dt_wrt = 1.0d-2 * ( t_end - t ); iCycleW = -1

  IF( dt_wrt .GT. Zero .AND. iCycleW .GT. 0 ) &
    STOP 'dt_wrt and iCycleW cannot both be present'

  CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

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
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           SuppressTally_Option = SuppressTally )

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

      IF( ProgramName .EQ. 'StandingAccretionShock' )THEN

        WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A4,A5,ES13.6E3,A3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t / Millisecond, ' ms ', &
          'dt = ', dt / Millisecond, ' ms'

      ELSE

        WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A1,A5,ES13.6E3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

      END IF

    END IF

    CALL TimersStop_Euler( Timer_Euler_InputOutput )

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, uDF, ComputeIncrement_Euler_DG_Explicit )

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
             ( iX_B0, iX_E0, &
               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

        wrt = .FALSE.

      END IF
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

    CALL WriteFieldsHDF &
           ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

    CALL TimersStop_Euler( Timer_Euler_InputOutput )

  END IF

  CALL ComputeTally_Euler_Relativistic_IDEAL &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

  CALL TimersStart_Euler( Timer_Euler_Finalize )

  CALL FinalizeTally_Euler_Relativistic_IDEAL

  CALL FinalizePositivityLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeEquationOfState

  CALL FinalizeProgram

  CALL TimersStop_Euler( Timer_Euler_Finalize )

  CALL FinalizeTimers_Euler

END PROGRAM ApplicationDriver
