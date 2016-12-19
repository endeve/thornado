PROGRAM HomogeneousSphere1D

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Centimeter, &
    Kilometer, &
    Gram, &
    Second, &
    Microsecond, &
    Erg
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
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
             = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
           nNodes_Option &
             = 2, &
           CoordinateSystem_Option &
             = 'SPHERICAL' )

  CALL InitializeHomogeneousSphere

  CALL FinalizeProgram

END PROGRAM HomogeneousSphere1D
