PROGRAM GaussianSphericalWave1D

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE StreamingRadiationInitializationModule, ONLY: &
    InitializeGaussianSphericalWave1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'GaussianSphericalWave1D', &
           nX_Option &
             = [ 32, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 10, 0, 0 ], &
           xL_Option &
             = [  0.2_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 10.2_DP, 1.0_DP, 1.0_DP ], &
           nE_Option &
             = 1, &
           eL_Option &
             = 0.0_DP, &
           eR_Option &
             = 1.0_DP, &
           nNodes_Option &
             = 1, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           RadiationSolver_Option &
             = 'M1_DG', &
           EvolveRadiation_Option &
             = .TRUE., &
           nStagesSSPRK_Option &
             = 1 )

  CALL InitializeGaussianSphericalWave1D

  CALL EvolveFields &
         ( t_begin = 0.0d+0, t_end = 7.0d+0, dt_write = 5.0d-1, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM GaussianSphericalWave1D
