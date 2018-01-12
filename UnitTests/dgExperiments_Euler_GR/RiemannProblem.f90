PROGRAM RiemannProblem

  USE KindModule, ONLY: &
    DP
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
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule_Beta, ONLY: &
    ComputeGeometryX
  USE FluidFieldsModule, ONLY: &
    uCF, rhsCF, uPF, uAF, iAF_P, iAF_Cs
  USE InputOutputModule, ONLY: &
    WriteFields1D
  USE InitializationModule_GR, ONLY: &
    InitializeFields_RiemannProblem
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE dgDiscretizationModule_Euler_GR, ONLY: &
     ComputeIncrement_Euler_GR_DG_Explicit
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputeFromConserved


  IMPLICIT NONE

  INTEGER  :: iCycle, iCycleD, iCycleW, nX = 512
  REAL(DP) :: t, dt, t_end, c = 1.0_DP, L = 1.0_DP, CFL = 0.5_DP

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'RiemannProblem', &
           nX_Option &
             = [ nX, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 2, 0, 0 ], &
           xL_Option &
             = [ 0.0d0, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ L, 1.0d0, 1.0d0 ], &
           nNodes_Option &
             = 1, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           EquationOfState_Option &
             = 'IDEAL', &
           FluidRiemannSolver_Option & ! --- Hard-coded
             = 'HLLC', &
           Gamma_IDEAL_Option &
             = 4.0_DP / 3.0_DP, &
           Opacity_Option &
             = 'IDEAL', &
           nStages_SSP_RK_Option &
             = 1 )

  t       = 0.0_DP
  dt      = CFL * L / ( c * nX )
  t_end   = 0.2_DP
  iCycleD = 10
  iCycleW = 1

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeFields_RiemannProblem &
         ( D_L = 1.000_DP, V_L = [ 0.0_DP, 0.0_DP, 0.0_DP], P_L = 1.0d+0, &
           D_R = 0.125_DP, V_R = [ 0.0_DP, 0.0_DP, 0.0_DP], P_R = 1.0d-1, &
           X_D_Option = 0.5_DP )

  CALL WriteFields1D( t, .TRUE., .TRUE. )

  CALL InitializeFluid_SSPRK( nStages = 1 )

  iCycle = 0

  DO WHILE ( t < t_end )

    iCycle = iCycle + 1

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    CALL UpdateFluid_SSPRK &
         ( t, dt, uGF, uCF, ComputeIncrement_Euler_GR_DG_Explicit )

    ! --- Update primitive fluid variables, pressure, and sound speed
    CALL ComputeFromConserved( iX_B1, iX_E1, uGF, uCF, uPF, uAF )

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL WriteFields1D( t, .TRUE., .TRUE. )

    END IF

  END DO

  CALL WriteFields1D( t, .TRUE., .TRUE. )

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeProgram

END PROGRAM RiemannProblem
