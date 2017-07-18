PROGRAM HomogeneousSphere1D

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE GravitySolutionModule, ONLY: &
    SolveGravity
  USE InputOutputModule, ONLY: &
    WriteFields1D
  USE GravityProblemsInitializationModule, ONLY: &
    InitializeHomogeneousSphere

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'HomogeneousSphere1D', &
           nX_Option &
             = [ 32, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 3, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 2.0_DP, Pi,     4.0_DP ], &
           nNodes_Option &
             = 3, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           GravitySolver_Option &
             = 'Newtonian_Poseidon' )

  CALL InitializeHomogeneousSphere &
         ( SphereRadius_Option = 1.0_DP )

  CALL SolveGravity

  CALL WriteFields1D &
         ( Time = 0.0_DP, WriteGeometryFields_Option = .TRUE. )

  CALL FinalizeProgram

END PROGRAM HomogeneousSphere1D
