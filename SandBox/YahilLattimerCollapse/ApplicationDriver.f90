PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Half, One, Two, &
    Pi, TwoPi, FourPi
  USE UnitsModule, ONLY: &
    Centimeter, &
    Kilometer, &
    Millisecond, &
    Gram, Erg
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
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF
  USE SlopeLimiterModule_Euler, ONLY: &
    InitializeSlopeLimiter_Euler, &
    FinalizeSlopeLimiter_Euler
  USE PositivityLimiterModule_Euler, ONLY: &
    InitializePositivityLimiter_Euler, &
    FinalizePositivityLimiter_Euler
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE InitializationModule, ONLY: &
    InitializeFields
  USE EulerEquationsUtilitiesModule_Beta, ONLY: &
    ComputeFromConserved, &
    ComputeTimeStep
  USE dgDiscretizationModule_Euler, ONLY: &
    ComputeIncrement_Euler_DG_Explicit

  IMPLICIT NONE

  CHARACTER(64), PARAMETER :: FileName = 'YahilHomologousCollapse_Gm_130.dat'
  REAL(DP),      PARAMETER :: Gamma           = 1.3_DP
  REAL(DP),      PARAMETER :: CollapseTime    = 200_DP * Millisecond
  REAL(DP),      PARAMETER :: CentralDensity  = 7.0d9 * Gram / Centimeter**3
  REAL(DP),      PARAMETER :: CentralPressure = 6.0d27 * Erg / Centimeter**3
  REAL(DP),      PARAMETER :: CoreRadius      = 1.0d4 * Kilometer

  INTEGER             :: iCycle, iCycleD, iCycleW
  INTEGER             :: nX(3), nNodes
  REAL(DP)            :: t, dt, t_end

  nX = [ 128, 1, 1 ]

  nNodes = 3

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'YahilLattimerCollapse', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 11, 3, 1 ], &
           xL_Option &
             = [ 0.0d0,      0.0d0, 0.0d0 ], &
           xR_Option &
             = [ CoreRadius, Pi,    TwoPi ], &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           ActivateUnits_Option &
             = .TRUE., &
           BasicInitialization_Option &
             = .TRUE. )

  t_end   = 1.0d+2
  iCycleD = 10
  iCycleW = 1000

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL InitializeSlopeLimiter_Euler &
         ( BetaTVD_Option = 1.15_DP, &
           BetaTVB_Option = 0.0_DP, &
           SlopeTolerance_Option = 1.0d-2, &
           UseSlopeLimiter_Option = .FALSE., &
           UseTroubledCellIndicator_Option = .FALSE., &
           LimiterThresholdParameter_Option = 0.12_DP )

  CALL InitializePositivityLimiter_Euler &
         ( Min_1_Option = 1.0d-12, &
           Min_2_Option = 1.0d-12, &
           UsePositivityLimiter_Option = .FALSE. )

  CALL InitializeFields &
         ( FileName, Gamma, CollapseTime, CentralDensity, CentralPressure )

  CALL WriteFieldsHDF &
         ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

!!$  ! --- Evolve ---
!!$
!!$  t = 0.0_DP
!!$
!!$  iCycle = 0
!!$  DO WHILE ( t < t_end )
!!$
!!$    iCycle = iCycle + 1
!!$
!!$    CALL ComputeTimeStep &
!!$           ( iX_B0, iX_E0, &
!!$             uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
!!$             uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
!!$             CFL = 0.3_DP / ( Two * DBLE( nNodes - 1 ) + One ), TimeStep = dt )
!!$
!!$    IF( t + dt > t_end )THEN
!!$
!!$      dt = t_end - t
!!$
!!$    END IF
!!$
!!$    IF( MOD( iCycle, iCycleD ) == 0 )THEN
!!$
!!$      WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A1,A5,ES13.6E3)') &
!!$          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt
!!$
!!$    END IF
!!$
!!$    CALL UpdateFluid_SSPRK &
!!$           ( t, dt, uGF, uCF, ComputeIncrement_Euler_DG_Explicit )
!!$
!!$    t = t + dt
!!$
!!$    IF( MOD( iCycle, iCycleW ) == 0 )THEN
!!$
!!$      CALL ComputeFromConserved &
!!$             ( iX_B0, iX_E0, &
!!$               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
!!$               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
!!$               uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
!!$               uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )
!!$
!!$      CALL WriteFieldsHDF &
!!$             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )
!!$
!!$    END IF
!!$
!!$  END DO

  ! --- Finalize ---

  CALL FinalizePositivityLimiter_Euler

  CALL FinalizeSlopeLimiter_Euler

  CALL FinalizeEquationOfState

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
