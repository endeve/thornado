PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, One, Two
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
    iX_B0, iX_B1, iX_E0, iX_E1
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE InitializationModule_GR, ONLY: &
    InitializeFields_GR
  USE SlopeLimiterModule_Euler_GR, ONLY: &
    InitializeSlopeLimiter_Euler_GR, &
    FinalizeSlopeLimiter_Euler_GR, &
    ApplySlopeLimiter_Euler_GR
  USE PositivityLimiterModule_Euler_GR, ONLY: &
    InitializePositivityLimiter_Euler_GR, &
    FinalizePositivityLimiter_Euler_GR, &
    ApplyPositivityLimiter_Euler_GR
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputeFromConserved_GR, &
    ComputeTimeStep_GR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE dgDiscretizationModule_Euler_GR, ONLY: &
    ComputeIncrement_Euler_GR_DG_Explicit
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK

  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: RiemannProblemName
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: wrt
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UseCharacteristicLimiting
  LOGICAL       :: UseTroubledCellIndicator
  LOGICAL       :: UsePositivityLimiter
  INTEGER       :: iCycle, iCycleD
  INTEGER       :: nX(3), bcX(3), nNodes
  INTEGER       :: nStagesSSPRK
  REAL(DP)      :: SlopeTolerance
  REAL(DP)      :: Min_1, Min_2
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, wTime, CFL
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: LimiterThresholdParameter
  REAL(DP)      :: nRel
  INTEGER       :: nDetCells

  ProgramName = 'RiemannProblem'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'RiemannProblem' )

      RiemannProblemName = 'PerturbedShockTube'

      CoordinateSystem = 'CARTESIAN'

      Gamma = 5.0_DP / 3.0_DP

      nX = [ 100, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 2, 0, 0 ]

      nNodes = 3

      BetaTVD = 2.0_DP
      BetaTVB = 0.0_DP

      UseSlopeLimiter           = .TRUE.
      SlopeTolerance            = 1.0d-6
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.1_DP

      UsePositivityLimiter = .TRUE.
      Min_1 = 1.0d-16
      Min_2 = 1.0d-16

      iCycleD = 100
      t_end   = 3.5d-1
      dt_wrt  = 1.0d-2 * t_end

      nStagesSSPRK = nNodes
      CFL          = 0.1_DP

    CASE( 'SedovBlastWave' )

      nDetCells = 1
      nRel      = 1.0d3

      CoordinateSystem = 'SPHERICAL'

      Gamma = 4.0_DP / 3.0_DP

      nX = [ 128, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 2, 0, 0 ]

      nNodes = 1

      BetaTVD = 2.0_DP
      BetaTVB = 0.0_DP

      UseSlopeLimiter           = .TRUE.
      SlopeTolerance            = 1.0d-6
      UseCharacteristicLimiting = .TRUE.

      UseTroubledCellIndicator  = .TRUE.
      LimiterThresholdParameter = 0.1_DP

      UsePositivityLimiter = .TRUE.
      Min_1 = 1.0d-12
      Min_2 = 1.0d-12

      iCycleD = 1
      t_end   = 0.1_DP
      dt_wrt  = 1.0d-2 * t_end

      nStagesSSPRK = nNodes
      CFL          = 0.1_DP

  END SELECT

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 0, 0 ], &
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

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeFields_GR &
         ( RiemannProblemName_Option &
             = TRIM( RiemannProblemName ) )

  CALL WriteFieldsHDF &
       ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL InitializeFluid_SSPRK( nStages = nStagesSSPRK )

  CALL InitializeSlopeLimiter_Euler_GR &
         ( BetaTVD_Option = BetaTVD, &
           BetaTVB_Option = BetaTVB, &
           SlopeTolerance_Option &
             = SlopeTolerance, &
           UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           UseCharacteristicLimiting_Option &
             = UseCharacteristicLimiting, &
           UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator, &
           LimiterThresholdParameter_Option &
             = LimiterThresholdParameter )

  CALL InitializePositivityLimiter_Euler_GR &
         ( Min_1_Option = Min_1, &
           Min_2_Option = Min_2, &
           UsePositivityLimiter_Option = UsePositivityLimiter )

  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'Evolving Fields...'
  WRITE(*,*)

  t     = 0.0_DP
  t_wrt = dt_wrt
  wrt   = .FALSE.

  iCycle = 0
  DO WHILE( t < t_end )

    iCycle = iCycle + 1

    CALL ComputeTimeStep_GR &
           ( iX_B0, iX_E0, &
             uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             CFL = CFL, TimeStep = dt )

    IF( t + dt .LT. t_end )THEN
      t = t + dt
    ELSE
      dt = t_end - t
      t  = t_end
    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A1,A5,ES13.6E3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, ComputeIncrement_Euler_GR_DG_Explicit )

    CALL ComputeFromConserved_GR &
           ( iX_B0, iX_E0, &
             uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

    IF( t + dt .LT. t_wrt )THEN
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.
    END IF

    IF( wrt )THEN

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

      wrt = .FALSE.

    END IF

  END DO

  CALL ComputeFromConserved_GR &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL FinalizePositivityLimiter_Euler_GR

  CALL FinalizeSlopeLimiter_Euler_GR

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX 

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeEquationOfState

  CALL FinalizeProgram
   
END PROGRAM ApplicationDriver
