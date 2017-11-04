PROGRAM SineWaveDamping1D

  USE KindModule, ONLY: &
    DP
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE OpacityModule_IDEAL, ONLY: &
    InitializeOpacities_IDEAL
  USE InitializationModule, ONLY: &
    InitializeSineWaveDamping1D
  USE ErrorAnalysisModule, ONLY: &
    InitializeErrorAnalysis, &
    FinalizeErrorAnalysis
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    IMEX

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'SineWaveDamping1D', &
           nX_Option &
             = [ 32, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 1, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
           nE_Option &
             = 1, &
           eL_Option &
             = 0.0d0, &
           eR_Option &
             = 1.0d0, &
           nNodes_Option &
             = 3, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           EquationOfState_Option &
             = 'IDEAL', &
           Opacity_Option &
             = 'IDEAL', &
           FluidRadiationCoupling_Option &
             = 'ConstantOpacities', &
           EvolveFluid_Option &
             = .FALSE., &
           RadiationSolver_Option &
             = 'M1_DG', &
           RadiationRiemannSolver_Option &
             = 'HLL', &
           EvolveRadiation_Option &
             = .TRUE., &
           ApplySlopeLimiter_Option &
             = .FALSE., &
           BetaTVB_Option &
             = 0.0d0, &
           BetaTVD_Option &
             = 1.8d0, &
           ApplyPositivityLimiter_Option &
             = .FALSE., &
           IMEX_Scheme_Option &
             = 'IMEX_RKCB2' )

  CALL InitializeOpacities_IDEAL &
         ( Eta_Option = 0.0d0, &
           Chi_Option = 1.0d1, &
           Sig_Option = 0.0d0 )

  CALL InitializeSineWaveDamping1D

  CALL InitializeErrorAnalysis( Time = 1.0d+0 )

  CALL EvolveFields &
         ( t_begin  = 0.0d+0, &
           t_end    = 1.0d+0, &
           dt_write = 1.0d-1, &
           UpdateFields = IMEX )

  CALL FinalizeErrorAnalysis

  CALL FinalizeProgram

END PROGRAM SineWaveDamping1D
