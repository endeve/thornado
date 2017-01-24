PROGRAM ComputeEigensystem_NuclearEOS

  USE KindModule, ONLY: &
    DP
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, & 
    ComputePressure_TABLE, & 
    ComputeSpecificInternalEnergy_TABLE
  USE UnitsModule, ONLY: &
    Gram, Centimeter, Kelvin

  IMPLICIT NONE

  REAL(DP), DIMENSION(1) :: D, T, Y, P, E

  D = [ 1.00d12 * Gram / Centimeter**3 ]
  T = [ 14.5d3 * Kelvin ]
  Y = [ 0.4 ]

  CALL ComputePressure_TABLE(D, T, Y, P)
  !CALL ComputeSpecificInternalEnergy_TABLE(D, T, Y, E)
  
END PROGRAM ComputeEigensystem_NuclearEOS
