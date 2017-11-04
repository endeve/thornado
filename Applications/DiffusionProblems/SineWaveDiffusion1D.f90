PROGRAM SineWaveDiffusion1D

  USE KindModule, ONLY: &
    DP
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE OpacityModule_IDEAL, ONLY: &
    InitializeOpacities_IDEAL
  USE InitializationModule, ONLY: &
    InitializeSineWaveDiffusion1D
  USE ErrorAnalysisModule, ONLY: &
    InitializeErrorAnalysis, &
    FinalizeErrorAnalysis
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SI_RK, &
    IMEX

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'SineWaveDiffusion1D', &
           nX_Option &
             = [ 8, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 1, 0, 0 ], &
           xL_Option &
             = [ - 3.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ + 3.0_DP, 1.0_DP, 1.0_DP ], &
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
             = 'ElasticScattering', &
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
           nStages_SI_RK_Option &
             = 2, &
           IMEX_Scheme_Option &
             = 'IMEX_RKCB2' )

  CALL InitializeOpacities_IDEAL &
         ( Eta_Option = 0.0d0, &
           Chi_Option = 0.0d0, &
           Sig_Option = 1.0d2 )

  CALL InitializeSineWaveDiffusion1D

  CALL InitializeErrorAnalysis( Time = 2.5d+2 )

  CALL EvolveFields &
         ( t_begin  = 0.0d+0, &
           t_end    = 2.5d+2, &
           dt_write = 5.0d+1, &
           UpdateFields = IMEX )

  CALL FinalizeErrorAnalysis

  CALL FinalizeProgram

END PROGRAM SineWaveDiffusion1D
