PROGRAM rhsTest

  USE KindModule, ONLY: &
    DP, Pi, TwoPi
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
  USE InitializationModule, ONLY: &
    InitializeFields
  USE dgDiscretizationModule_Euler, ONLY: &
    ComputeRHS_Euler

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: mpierr

  CALL MPI_INIT( mpierr )

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'rhsTest', &
           nX_Option &
             = [ 32, 32, 32 ], &
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

  CALL ComputeGeometryX

  CALL InitializeFields

  CALL ComputeRHS_Euler

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeProgram

  CALL MPI_FINALIZE( mpierr )

END PROGRAM rhsTest
