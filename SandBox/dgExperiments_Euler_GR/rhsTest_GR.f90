PROGRAM rhsTest_GR

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
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule_Beta, ONLY: &
    ComputeGeometryX
  USE FluidFieldsModule, ONLY: &
    uCF, rhsCF
  USE InitializationModule_GR, ONLY: &
    InitializeFields_GR
  USE dgDiscretizationModule_Euler_GR, ONLY: &
    ComputeIncrement_Euler_GR_DG_Explicit

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'rhsTest_GR', &
           nX_Option &
             = [ 32, 16, 32 ], &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 1, 1, 1 ], &
           xL_Option &
             = [ 0.2d0, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 1.2d0,    Pi, TwoPi ], &
           nNodes_Option &
             = 3, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           EquationOfState_Option &
             = 'IDEAL', &
           Gamma_IDEAL_Option &
             = 4.0_DP / 3.0_DP, &
           Opacity_Option &
             = 'IDEAL', &
           nStages_SSP_RK_Option &
             = 1 )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass_Option = 0.05_DP )

  CALL InitializeFields_GR

  CALL ComputeIncrement_Euler_GR_DG_Explicit &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, rhsCF )

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeProgram

END PROGRAM rhsTest_GR
