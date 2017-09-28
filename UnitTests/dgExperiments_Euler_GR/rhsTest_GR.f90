PROGRAM rhsTest_GR

  USE KindModule, ONLY: &
    DP, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE InitializationModule_GR, ONLY: &
    InitializeFields_GR
  USE dgDiscretizationModule_Euler_GR, ONLY: &
    ComputeRHS_Euler_GR

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: mpierr

  CALL MPI_INIT( mpierr )

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'rhsTest_GR', &
           nX_Option &
             = [ 32, 1, 1 ], &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 1, 1, 1 ], &
           xL_Option &
             = [ 0.0d0, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 1.0d0, 1.0d0, 1.0d0 ], &
           nNodes_Option &
             = 4, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           EquationOfState_Option &
             = 'IDEAL', &
           Opacity_Option &
             = 'IDEAL' )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX( [ 32, 1, 1 ], [ 4, 4, 4 ], [ 1, 1, 1 ] )

  CALL InitializeFields

  CALL ComputeRHS_Euler_GR

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeProgram

  CALL MPI_FINALIZE( mpierr )

END PROGRAM rhsTest_GR
