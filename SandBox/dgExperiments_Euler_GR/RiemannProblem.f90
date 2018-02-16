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
  USE SlopeLimiterModule_Euler_GR, ONLY: &
    InitializeSlopeLimiter, &
    FinalizeSlopeLimiter
  USE dgDiscretizationModule_Euler_GR, ONLY: &
    ComputeIncrement_Euler_GR_DG_Explicit
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputeFromConserved
  USE RiemannProblemInitializer, ONLY: &
    RiemannProblemChoice


  IMPLICIT NONE

  INTEGER  :: iCycle, iCycleD, iCycleW, K
  REAL(DP) :: t, dt, t_end, xR, x_D, CFL, Gamma, c = 1.0_DP
  REAL(DP) :: D_L, V_L(3), P_L, D_R, V_R(3), P_R

!  OPEN(  10 , FILE = 'nIter.dat' )

  CALL RiemannProblemChoice &
         ( D_L, V_L, P_L, D_R, V_R, P_R, &
           xR, x_D, K, t, t_end, CFL, Gamma, iRP = 0 )

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'RiemannProblem', &
           nX_Option &
             = [ K, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 1, 0, 0 ], &
           xL_Option &
             = [ 0.0d0, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ xR, 1.0d0, 1.0d0 ], &
           nNodes_Option &
             = 2, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           EquationOfState_Option &
             = 'IDEAL', &
           FluidRiemannSolver_Option & ! --- Hard-coded
             = 'HLLC', &
           Gamma_IDEAL_Option &
             = Gamma, &
           Opacity_Option &
             = 'IDEAL', &
           nStages_SSP_RK_Option &
             = 1 )

  dt      = CFL * xR / ( c * K )
  iCycleD = 100 * t_end
  iCycleW = 10  * t_end

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  CALL InitializeFields_RiemannProblem &
         ( D_L = D_L, V_L = V_L, P_L = P_L, &
           D_R = D_R, V_R = V_R, P_R = P_R, &
           X_D_Option = x_D )

  CALL WriteFields1D( t, .TRUE., .TRUE. )

  CALL InitializeFluid_SSPRK( nStages = 2 )

!  CALL InitializeSlopeLimiter

  iCycle = 0

  DO WHILE ( t < t_end )

    WRITE(*,*) 't:' , t

    iCycle = iCycle + 1

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    CALL UpdateFluid_SSPRK &
         ( t, dt, uGF, uCF, ComputeIncrement_Euler_GR_DG_Explicit )

    ! --- Update primitive fluid variables, pressure, and sound speed
    !CALL ComputeFromConserved( iX_B0, iX_E0, uGF, uCF, uPF, uAF )

!    WRITE(*,*) 'Finished first time step, stopping...'
!    STOP

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL WriteFields1D( t, .TRUE., .TRUE. )

    END IF

  END DO

  CALL WriteFields1D( t, .TRUE., .TRUE. )

!  CALL FinalizeSlopeLimiter

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeProgram

!  CLOSE(10)

END PROGRAM RiemannProblem
