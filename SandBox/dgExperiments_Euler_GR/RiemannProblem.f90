PROGRAM RiemannProblem

  USE KindModule, ONLY: &
    DP, Pi, TwoPi
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
  USE PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter, &
    FinalizePositivityLimiter
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE InitializationModule_GR, ONLY: &
    InitializeFields_RiemannProblem
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE SlopeLimiterModule_Euler_GR, ONLY: &
    InitializeSlopeLimiter_Euler_GR, &
    FinalizeSlopeLimiter_Euler_GR
  USE dgDiscretizationModule_Euler_GR, ONLY: &
    ComputeIncrement_Euler_GR_DG_Explicit
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputeFromConserved
  USE RiemannProblemInitializer, ONLY: &
    RiemannProblemChoice
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState

  IMPLICIT NONE

  INTEGER  :: iCycle, iCycleD, iCycleW, K, bcX(3)
  REAL(DP) :: t, dt, t_end, xL(3), xR(3), x_D, CFL, Gamma, c = 1.0_DP
  REAL(DP) :: D_L, V_L(3), P_L, D_R, V_R(3), P_R
  CHARACTER( LEN = 11 ) :: CS

  REAL(DP)             :: LT
  CHARACTER( LEN = 4 ) :: arg
  INTEGER              :: argv(2), nNodes, i
  LOGICAL              :: ConvergenceRate = .FALSE.

  CALL RiemannProblemChoice &
         ( D_L, V_L, P_L, D_R, V_R, P_R, &
             xL, xR, x_D, K, t, t_end, CFL, Gamma, bcX, CS, iRP = 10 )

  IF ( ConvergenceRate ) THEN
    DO i = 1, IARGC()
      CALL GETARG( i, arg )
      READ( arg, * ) argv(i)
    END DO
    nNodes = argv(1)
    K      = argv(2)
  ELSE
    nNodes = 3
    LT     = 0.1_DP
  END IF
  
  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'RiemannProblem', &
           nX_Option &
             = [ K, 1, 1 ], &
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
             = TRIM( CS ), &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  dt      = CFL * ( xR(1) - xL(1) ) / ( c * K )
  iCycleD = 1
  iCycleW = 1

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeFields_RiemannProblem &
         ( D_L = D_L, V_L = V_L, P_L = P_L, &
           D_R = D_R, V_R = V_R, P_R = P_R, &
           X_D_Option = x_D )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL InitializeFluid_SSPRK( nStages = 3 )

  CALL InitializeSlopeLimiter_Euler_GR &
         ( BetaTVD_Option = 2.0_DP, &
           BetaTVB_Option = 0.0_DP, &
           SlopeTolerance_Option = 1.0d-6, &
           UseSlopeLimiter_Option = .TRUE., &
           UseCharacteristicLimiting_Option = .FALSE., &
           UseTroubledCellIndicator_Option = .TRUE., &
           LimiterThresholdParameter_Option = LT )

  CALL InitializePositivityLimiter &
         ( Min_1_Option = 1.0d-16 , Min_2_Option = 1.0d-16, &
           UsePositivityLimiter_Option = .TRUE. )

  iCycle = 0

  DO WHILE ( t < t_end )

    IF( t + dt < t_end )THEN
      t = t + dt
    ELSE
      dt = t_end - t
      t  = t_end
    END IF

    iCycle = iCycle + 1

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, ComputeIncrement_Euler_GR_DG_Explicit )

    ! --- Update primitive fluid variables, pressure, and sound speed ---
    CALL ComputeFromConserved( iX_B0, iX_E0, uGF, uCF, uPF, uAF )

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

    END IF

  END DO

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL FinalizePositivityLimiter

  CALL FinalizeSlopeLimiter_Euler_GR

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeEquationOfState

  CALL FinalizeProgram

END PROGRAM RiemannProblem
