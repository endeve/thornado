PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, One, Two, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1
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
    WriteFieldsHDF
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF
  USE InitializationModule, ONLY: &
    InitializeFields
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE SlopeLimiterModule_Euler, ONLY: &
    InitializeSlopeLimiter_Euler, &
    FinalizeSlopeLimiter_Euler, &
    ApplySlopeLimiter_Euler
  USE PositivityLimiterModule_Euler, ONLY: &
    InitializePositivityLimiter_Euler, &
    FinalizePositivityLimiter_Euler, &
    ApplyPositivityLimiter_Euler
  USE EulerEquationsUtilitiesModule_Beta, ONLY: &
    ComputeFromConserved, &
    ComputeTimeStep
  USE dgDiscretizationModule_Euler, ONLY: &
    ComputeIncrement_Euler_DG_Explicit
  USE Euler_TallyModule, ONLY: &
    InitializeTally_Euler, &
    FinalizeTally_Euler, &
    ComputeTally_Euler

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
  INTEGER       :: iCycle, iCycleD
  INTEGER       :: nX(3), bcX(3), nNodes
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, wTime
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: LimiterThresholdParameter

  CoordinateSystem = 'CARTESIAN'

  ProgramName = 'KelvinHelmholtz'

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

      nX = [ 100, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 2, 0, 0 ]

      nNodes = 3

      BetaTVD = 2.00_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 0.03_DP

      iCycleD = 1
      t_end   = 2.0d-1
      dt_wrt  = 1.0d-2

   CASE( 'RiemannProblemSpherical' )

      CoordinateSystem = 'SPHERICAL'

      Gamma = 1.4_DP

      nX = [ 128, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 2.0_DP, Pi,     TwoPi  ]

      bcX = [ 3, 3, 0 ]

      nNodes = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 0.03_DP

      iCycleD = 1
      t_end   = 5.0d-1
      dt_wrt  = 2.5d-2

    CASE( 'KelvinHelmholtz' )

      Gamma = 5.0_DP / 3.0_DP

      nX = [ 128, 128, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 0 ]

      nNodes = 3

      BetaTVD = 2.00_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.03_DP

      iCycleD = 10
      t_end   = 1.50_DP
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

      iCycleD = 10
      t_end   = 8.5_DP
      dt_wrt  = 0.1_DP

    CASE( 'Implosion' )

      Gamma = 1.4_DP
      
      nX = [ 128, 128, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 0.3_DP, 0.3_DP, 1.0_DP ]

      bcX = [ 3, 3, 0 ]

      nNodes = 3

      BetaTVD = 2.00_DP
      BetaTVB = 0.0d+02

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.03_DP

      iCycleD = 10
      t_end   = 2.500_DP
      dt_wrt  = 0.045_DP

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

  CALL InitializeSlopeLimiter_Euler &
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

  CALL InitializePositivityLimiter_Euler &
         ( Min_1_Option = 1.0d-12, &
           Min_2_Option = 1.0d-12, &
           UsePositivityLimiter_Option = .TRUE. )

  CALL InitializeFields &
         ( AdvectionProfile_Option &
             = TRIM( AdvectionProfile ), &
           Direction_Option &
             = TRIM( Direction ), &
           RiemannProblemName_Option &
             = TRIM( RiemannProblemName ) )

  CALL ApplySlopeLimiter_Euler &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
           uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) )

  CALL ApplyPositivityLimiter_Euler &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
           uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) )

  CALL ComputeFromConserved &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  CALL WriteFieldsHDF &
         ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL InitializeFluid_SSPRK( nStages = 3 )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t     = 0.0_DP
  t_wrt = dt_wrt
  wrt   = .FALSE.

  CALL InitializeTally_Euler &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  iCycle = 0
  DO WHILE ( t < t_end )

    iCycle = iCycle + 1

    CALL ComputeTimeStep &
           ( iX_B0, iX_E0, &
             uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             CFL = 0.3_DP / ( Two * DBLE( nNodes - 1 ) + One ), TimeStep = dt )

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A1,A5,ES13.6E3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, ComputeIncrement_Euler_DG_Explicit )

    t = t + dt

    IF( wrt )THEN

      CALL ComputeFromConserved &
             ( iX_B0, iX_E0, &
               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

      CALL ComputeTally_Euler &
           ( iX_B0, iX_E0, &
             uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

      wrt = .FALSE.

    END IF

  END DO

  CALL ComputeFromConserved &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL ComputeTally_Euler &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

  CALL FinalizeTally_Euler

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  CALL FinalizePositivityLimiter_Euler

  CALL FinalizeSlopeLimiter_Euler

  CALL FinalizeEquationOfState

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
