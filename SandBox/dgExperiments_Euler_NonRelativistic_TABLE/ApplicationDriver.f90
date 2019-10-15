PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, One, Two, Pi, TwoPi
  USE UnitsModule, ONLY: &
    Millisecond, Microsecond, &
    Kilometer
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
  USE EquationOfStateModule_TABLE, ONLY: &
    MinD, MaxD, MinT, MaxT, MinY, MaxY
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF
  USE InitializationModule, ONLY: &
    InitializeFields
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
    InitializeSlopeLimiter_Euler_NonRelativistic_TABLE, &
    FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE, &
    ApplySlopeLimiter_Euler_NonRelativistic_TABLE
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    InitializePositivityLimiter_Euler_NonRelativistic_TABLE, &
    FinalizePositivityLimiter_Euler_NonRelativistic_TABLE, &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeFromConserved_Euler_NonRelativistic, &
    ComputeTimeStep_Euler_NonRelativistic
  USE Euler_dgDiscretizationModule, ONLY: &
    Euler_ComputeIncrement_DG_Explicit
  USE Euler_TallyModule_NonRelativistic_TABLE, ONLY: &
    InitializeTally_Euler_NonRelativistic_TABLE, &
    FinalizeTally_Euler_NonRelativistic_TABLE, &
    ComputeTally_Euler_NonRelativistic_TABLE

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32)  :: ProgramName
  CHARACTER(32)  :: RiemannProblemName
  CHARACTER(32)  :: CoordinateSystem
  CHARACTER(128) :: EosTableName
  LOGICAL        :: wrt
  LOGICAL        :: UseSlopeLimiter
  LOGICAL        :: UseCharacteristicLimiting
  LOGICAL        :: UseTroubledCellIndicator
  LOGICAL        :: UsePositivityLimiter
  INTEGER        :: iCycle, iCycleD
  INTEGER        :: nX(3), bcX(3), nNodes, nStages
  REAL(DP)       :: t, dt, t_end, dt_wrt, t_wrt, wTime
  REAL(DP)       :: xL(3), xR(3)
  REAL(DP)       :: BetaTVD, BetaTVB
  REAL(DP)       :: LimiterThresholdParameter

  ProgramName = 'RiemannProblem'

  EosTableName = 'wl-EOS-SFHo-15-25-50-noBCK.h5'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'RiemannProblem' )

      RiemannProblemName = 'Sod'

      CoordinateSystem = 'CARTESIAN'

      nX = [ 100, 1, 1 ]
      xL = [ - 5.0_DP,   0.0_DP, 0.0_DP ] * Kilometer
      xR = [ + 5.0_DP, + 1.0_DP, 1.0_DP ] * Kilometer

      bcX = [ 2, 0, 0 ]

      nNodes  = 3
      nStages = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 1.0d-2
      UsePositivityLimiter      = .TRUE.

      iCycleD = 10
      t_end   = 2.5d-2 * Millisecond
      dt_wrt  = 2.5d-4 * Millisecond

   CASE( 'RiemannProblemSpherical' )

      RiemannProblemName = 'SphericalSod'

      CoordinateSystem = 'SPHERICAL'

      nX = [ 128, 16, 1 ]
      xL = [ 1.0d-3 * Kilometer, 0.0_DP, 0.0_DP ]
      xR = [ 2.0_DP * Kilometer, Pi,     4.0_DP ]

      bcX = [ 3, 3, 0 ]

      nNodes = 3
      nStages = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 0.03_DP

      UsePositivityLimiter      = .TRUE.

      iCycleD = 10
      t_end   = 5.0d-1 * Millisecond
      dt_wrt  = 2.5d-2 * Millisecond

    CASE( 'Jet' )

      CoordinateSystem = 'CARTESIAN'

      nX = [ 100, 100, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer

      bcX = [ 2, 2, 0 ]

      nNodes  = 3
      nStages = 3

      BetaTVD = 1.75_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .FALSE.
      LimiterThresholdParameter = 1.5d-0
      UsePositivityLimiter      = .TRUE.

      iCycleD = 10
      t_end   = 2.5d-1 * Millisecond
      dt_wrt  = 2.5d-6 * Millisecond

    CASE( 'Implosion' )

      CoordinateSystem = 'CARTESIAN'

      nX = [ 64, 64, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR = [ 0.3_DP, 0.3_DP, 1.0_DP ] * Kilometer

      bcX = [ 3, 3, 0 ]

      nNodes  = 3
      nStages = 3

      BetaTVD = 1.50_DP
      BetaTVB = 0.0d+00

      UseSlopeLimiter           = .TRUE.
      UseCharacteristicLimiting = .TRUE.
      UsePositivityLimiter      = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.30_DP

      iCycleD = 1
      t_end   = 2.500_DP * Millisecond
      dt_wrt  = 0.045_DP * Millisecond

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
      WRITE(*,'(A)')     'Valid choices:'
      WRITE(*,'(A)')     '  RiemannProblem'
      WRITE(*,'(A)')     '  RiemannProblemSpherical'
      WRITE(*,'(A)')     '  Jet'
      WRITE(*,'(A)')     '  Implosion'
      WRITE(*,'(A)')     'Stopping...'
      STOP

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
           ActivateUnits_Option &
             = .TRUE., &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = TRIM( EosTableName ) )

  CALL InitializeSlopeLimiter_Euler_NonRelativistic_TABLE &
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

  CALL InitializePositivityLimiter_Euler_NonRelativistic_TABLE &
         ( UsePositivityLimiter_Option = UsePositivityLimiter, &
           Verbose_Option = .TRUE., &
           Min_1_Option = ( One + EPSILON(One) ) * MinD, &
           Min_2_Option = ( One + EPSILON(One) ) * MinT, &
           Min_3_Option = ( One + EPSILON(One) ) * MinY, &
           Max_1_Option = ( One - EPSILON(One) ) * MaxD, &
           Max_2_Option = ( One - EPSILON(One) ) * MaxT, &
           Max_3_Option = ( One - EPSILON(One) ) * MaxY )

  CALL InitializeFields &
         ( RiemannProblemName_Option &
             = TRIM( RiemannProblemName ) )

  CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
           uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) )

  CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
           uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) )

  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
           uCF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
           uPF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
           uAF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:) )

  CALL WriteFieldsHDF &
         ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL InitializeFluid_SSPRK( nStages )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t     = 0.0_DP * Millisecond
  t_wrt = dt_wrt
  wrt   = .FALSE.

  CALL InitializeTally_Euler_NonRelativistic_TABLE &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  iCycle = 0
  DO WHILE ( t < t_end )

    iCycle = iCycle + 1

    CALL ComputeTimeStep_Euler_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             uGF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
             uCF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
             CFL = 0.5_DP / ( Two * DBLE( nNodes - 1 ) + One ), TimeStep = dt )

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A4,A5,ES13.6E3,A3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t / Millisecond, ' ms ', 'dt = ', dt / Millisecond, ' ms'

    END IF

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, Euler_ComputeIncrement_DG_Explicit )

    t = t + dt

    IF( wrt )THEN

      CALL ComputeFromConserved_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, &
               uGF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
               uCF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
               uPF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
               uAF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:) )

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

      CALL ComputeTally_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, &
             uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

      wrt = .FALSE.

    END IF

  END DO

  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
           uCF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
           uPF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:), &
           uAF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:) )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL ComputeTally_Euler_NonRelativistic_TABLE &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

  CALL FinalizeTally_Euler_NonRelativistic_TABLE

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  CALL FinalizePositivityLimiter_Euler_NonRelativistic_TABLE

  CALL FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE

  CALL FinalizeEquationOfState

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
