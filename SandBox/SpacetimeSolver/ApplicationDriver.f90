PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, One, Two, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nDimsX, nDOFX
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF
  USE GeometryFieldsModule, ONLY: &
    nGF, uGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE FluidFieldsModule, ONLY: &
    nCF, nPF, nAF, &
    uCF, uPF, uAF, &
    uDF
  USE InitializationModule, ONLY: &
    InitializeFields
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE Euler_SlopeLimiterModule_NonRelativistic_IDEAL, ONLY: &
    InitializeSlopeLimiter_Euler_NonRelativistic_IDEAL, &
    FinalizeSlopeLimiter_Euler_NonRelativistic_IDEAL, &
    ApplySlopeLimiter_Euler_NonRelativistic_IDEAL
  USE Euler_PositivityLimiterModule_NonRelativistic_IDEAL, ONLY: &
    InitializePositivityLimiter_Euler_NonRelativistic_IDEAL, &
    FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL, &
    ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeFromConserved_Euler_NonRelativistic, &
    ComputeTimeStep_Euler_NonRelativistic
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit
  USE Euler_TallyModule_NonRelativistic, ONLY: &
    InitializeTally_Euler_NonRelativistic, &
    FinalizeTally_Euler_NonRelativistic, &
    ComputeTally_Euler_NonRelativistic
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
  CHARACTER(32) :: Direction
  CHARACTER(32) :: RiemannProblemName
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: wrt
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UseCharacteristicLimiting
  LOGICAL       :: UseTroubledCellIndicator
  LOGICAL       :: UseConservativeCorrection
  INTEGER       :: iCycle, iCycleD
  INTEGER       :: nX(3), bcX(3), nNodes
  INTEGER       :: RestartFileNumber
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, wTime
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: LimiterThresholdParameter
  REAL(DP)      :: Eblast

  TimeIt_Euler = .TRUE.
  CALL InitializeTimers_Euler
  CALL TimersStart_Euler( Timer_Euler_Initialize )

  CoordinateSystem = 'CARTESIAN'

  ProgramName = 'RiemannProblem'

  RestartFileNumber = -1

  t = 0.0_DP

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'Advection' )

      AdvectionProfile = 'TopHat'

      Direction = 'XY'

      Gamma = 1.4_DP

      nX = [ 64, 64, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 0 ]

      nNodes = 2

      BetaTVD = 1.80_DP
      BetaTVB = 0.0d+00

      iCycleD = 1
      t_end   = SQRT( 2.0d-0 )
      dt_wrt  = 1.0d-2 * t_end

    CASE( 'RiemannProblem' )

      RiemannProblemName = 'Sod'

      Gamma = 1.4_DP

      nX = [ 256, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 2, 0, 0 ]

      nNodes = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.03d0

      iCycleD = 10
      t_end   = 2.0d-1
      dt_wrt  = 1.0d-2 * t_end

   CASE( 'RiemannProblemSpherical' )

      CoordinateSystem = 'SPHERICAL'

      Gamma = 1.4_DP

      nX = [ 128, 16, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 2.0_DP, Pi,     TwoPi  ]

      bcX = [ 3, 3, 0 ]

      nNodes = 2

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 0.03_DP

      iCycleD = 1
      t_end   = 5.0d-1
      dt_wrt  = 2.5d-2

   CASE( 'SphericalSedov' )

      Eblast = 1.0d0

      CoordinateSystem = 'SPHERICAL'

      Gamma = 1.4_DP

      nX = [ 256, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.2_DP, Pi,     TwoPi  ]

      bcX = [ 2, 0, 0 ]

      nNodes = 3

      BetaTVD = 1.75d+00
      BetaTVB = 0.00d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.015_DP

      iCycleD = 1
      t_end   = 1.0d+0
      dt_wrt  = 5.0d-2

    CASE( 'IsentropicVortex' )

      Gamma = 1.4_DP

      nX = [ 50, 50, 1 ]
      xL = [ - 5.0_DP, - 5.0_DP, 0.0_DP ]
      xR = [ + 5.0_DP, + 5.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 0 ]

      nNodes = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .FALSE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.03_DP

      iCycleD = 10
      t_end   = 10.0_DP
      dt_wrt  = 10.0_DP

    CASE( 'KelvinHelmholtz' )

      Gamma = 5.0_DP / 3.0_DP

      nX = [ 256, 256, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 0 ]

      nNodes = 3

      BetaTVD = 1.50_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .FALSE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.03_DP

      iCycleD = 10
      t_end   = 3.00_DP
      dt_wrt  = 0.15_DP

    CASE( 'RayleighTaylor' )

      Gamma = 1.4_DP

      nX = [ 16, 48, 1 ]
      xL = [ - 0.25_DP, + 0.25_DP, 0.0_DP ]
      xR = [ - 0.75_DP, + 0.75_DP, 1.0_DP ]

      bcX = [ 1, 3, 0 ]

      nNodes = 3

      BetaTVD = 2.00_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 0.03_DP

      UseConservativeCorrection = .TRUE.

      iCycleD = 10
      t_end   = 8.5_DP
      dt_wrt  = 0.1_DP

    CASE( 'Implosion' )

      Gamma = 1.4_DP

      nX = [ 256, 256, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 0.3_DP, 0.3_DP, 1.0_DP ]

      bcX = [ 3, 3, 0 ]

      nNodes = 3

      BetaTVD = 1.50_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.30_DP

      iCycleD = 10
      t_end   = 2.500_DP
      dt_wrt  = 0.045_DP

    CASE( 'Explosion' )

      Gamma = 1.4_DP

      nX = [ 256, 256, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.5_DP, 1.5_DP, 1.0_DP ]

      bcX = [ 31, 31, 0 ]

      nNodes = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.05_DP

      iCycleD = 10
      t_end   = 3.20_DP
      dt_wrt  = 0.32_DP

  END SELECT

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
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeSlopeLimiter_Euler_NonRelativistic_IDEAL &
         ( BetaTVD_Option = BetaTVD, &
           BetaTVB_Option = BetaTVB, &
           SlopeTolerance_Option &
             = 1.0d-6, &
           UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           UseCharacteristicLimiting_Option &
             = UseCharacteristicLimiting, &
           UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator, &
           LimiterThresholdParameter_Option &
             = LimiterThresholdParameter )

  CALL InitializePositivityLimiter_Euler_NonRelativistic_IDEAL &
         ( UsePositivityLimiter_Option = .TRUE., &
           Verbose_Option = .TRUE., &
           Min_1_Option = 1.0d-12, &
           Min_2_Option = 1.0d-12 )

  CALL InitializeFields &
         ( AdvectionProfile_Option &
             = TRIM( AdvectionProfile ), &
           Direction_Option &
             = TRIM( Direction ), &
           RiemannProblemName_Option &
             = TRIM( RiemannProblemName ), &
           SedovEnergy_Option = Eblast )

  IF( RestartFileNumber .GE. 0 )THEN

    CALL ReadFieldsHDF( RestartFileNumber, t, ReadFF_Option = .TRUE. )

  END IF

  CALL ApplySlopeLimiter_Euler_NonRelativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

  CALL TimersStart_Euler( Timer_Euler_InputOutput )

  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL TimersStop_Euler( Timer_Euler_InputOutput )

  CALL InitializeFluid_SSPRK( nStages = 3 )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t_wrt = t + dt_wrt
  wrt   = .FALSE.

  CALL InitializeTally_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )

  CALL ComputeTally_Euler_NonRelativistic &
       ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
         SetInitialValues_Option = .TRUE., Verbose_Option = .TRUE. )

  CALL TimersStop_Euler( Timer_Euler_Initialize )

  iCycle = 0
  DO WHILE ( t < t_end )

    iCycle = iCycle + 1

    CALL ComputeTimeStep_Euler_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, &
             CFL = 0.5_DP / ( nDimsX * ( Two * DBLE( nNodes ) - One ) ), &
             TimeStep = dt )

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    CALL TimersStart_Euler( Timer_Euler_InputOutput )
    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A1,A5,ES13.6E3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF
    CALL TimersStop_Euler( Timer_Euler_InputOutput )

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, uDF, ComputeIncrement_Euler_DG_Explicit )

    t = t + dt

    CALL TimersStart_Euler( Timer_Euler_InputOutput )
    IF( wrt )THEN

      CALL ComputeFromConserved_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

      CALL ComputeTally_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
               Verbose_Option = .TRUE. )

      wrt = .FALSE.

    END IF
    CALL TimersStop_Euler( Timer_Euler_InputOutput )

  END DO

  CALL TimersStart_Euler( Timer_Euler_InputOutput )
  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )
  CALL TimersStop_Euler( Timer_Euler_InputOutput )

  CALL ComputeTally_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, Time = t, &
           Verbose_Option = .TRUE. )

  CALL TimersStart_Euler( Timer_Euler_Finalize )
  CALL FinalizeTally_Euler_NonRelativistic_IDEAL

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  CALL FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL

  CALL FinalizeSlopeLimiter_Euler_NonRelativistic_IDEAL

  CALL FinalizeEquationOfState

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

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
