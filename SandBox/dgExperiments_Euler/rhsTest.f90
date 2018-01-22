PROGRAM rhsTest

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
  USE InitializationModule, ONLY: &
    InitializeFields
  USE dgDiscretizationModule_Euler, ONLY: &
    ComputeIncrement_Euler_DG_Explicit

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(DP) :: wTime

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'rhsTest', &
           nX_Option &
             = [ 32, 16, 32 ], &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = [ 1, 1, 1 ], &
           xL_Option &
             = [ 0.0d0, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 1.0d0, Pi,    TwoPi ], &
           nNodes_Option &
             = 4, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           EquationOfState_Option &
             = 'IDEAL', &
           nStages_SSP_RK_Option &
             = 1 )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  wTime = MPI_WTIME( )
  CALL ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, uGF )
  wTime = MPI_WTIME( ) - wTime
  print*
  print*, "ComputeGeometryX = ", wTime
  print*

  wTime = MPI_WTIME( )
  CALL InitializeFields
  wTime = MPI_WTIME( ) - wTime
  print*
  print*, "InitializeFields = ", wTime
  print*

  wTime = MPI_WTIME( )
  CALL ComputeIncrement_Euler_DG_Explicit &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, rhsCF )
  wTime = MPI_WTIME( ) - wTime
  print*
  print*, "ComputeRHS_Euler = ", wTime
  print*

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeProgram

END PROGRAM rhsTest
