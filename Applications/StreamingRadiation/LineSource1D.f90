PROGRAM LineSource1D

  USE KindModule, ONLY: &
    DP, TwoPi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE StreamingRadiationInitializationModule, ONLY: &
    InitializeLineSource1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'LineSource1D', &
           nX_Option &
             = [ 1024, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 2, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 1.5_DP, TwoPi,  1.0_DP ], &
           nE_Option &
             = 1, &
           eL_Option &
             = 0.0_DP, &
           eR_Option &
             = 1.0_DP, &
           nNodes_Option &
             = 1, &
           CoordinateSystem_Option &
             = 'CYLINDRICAL', &
           RadiationSolver_Option &
             = 'M1_DG', &
           EvolveRadiation_Option &
             = .TRUE., &
           nStages_SSP_RK_Option &
             = 2 )

  CALL InitializeLineSource1D( Sigma_Option = 0.005_DP )

  CALL EvolveFields &
         ( t_begin = 0.0d+0, t_end = 1.0d+0, dt_write = 1.0d-1, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM LineSource1D
