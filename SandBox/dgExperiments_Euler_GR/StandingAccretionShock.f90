PROGRAM StandingAccretionShock

  USE KindModule, ONLY: &
    DP, Pi, TwoPi, Two
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE SlopeLimiterModule_Euler_GR, ONLY: &
    InitializeSlopeLimiter, &
    FinalizeSlopeLimiter
  USE PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter, &
    FinalizePositivityLimiter
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule_Beta, ONLY: &
    ComputeGeometryX
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE DataFileReader, ONLY: &
    ReadData, ReadParameters

  
  IMPLICIT NONE

  REAL(DP), ALLOCATABLE  :: FluidFieldData(:,:), FluidFieldParameters(:)

  CALL ReadParameters &
         ( '../StandingAccretionShock_Parameters.dat', FluidFieldParameters )
  CALL ReadData &
         ( '../StandingAccretionShock_Data.dat', FluidFieldData )

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'StandingAccretionShock', &
           nX_Option &
             = [ 64, 32, 1 ], &
           swX_Option &
             = [ 1, 1, 0 ], &
           bcX_Option &
             = [ 1, 1, 0 ], &
           xL_Option &
             = [ FluidFieldParameters(3), 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ Two * FluidFieldParameters(4), Pi,    TwoPi ], &
           nNodes_Option &
             = 2, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = 4.0_DP / 3.0_DP )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass_Option = 0.1_DP )

  CALL WriteFieldsHDF &
         ( 0.0_DP, WriteGF_Option = .TRUE., WriteFF_Option = .TRUE. )

  CALL InitializeFluid_SSPRK( nStages = 3 )

  CALL InitializeSlopeLimiter &
         ( BetaTVD_Option = 2.00_DP, &
           UseSlopeLimiter_Option = .TRUE. , &
           UseTroubledCellIndicator_Option = .TRUE. )

  CALL InitializePositivityLimiter &
         ( Min_1_Option = 1.0d-16 , Min_2_Option = 1.0d-16, &
           UsePositivityLimiter_Option = .TRUE. )

  ! --- 

  DEALLOCATE( FluidFieldParameters, FluidFieldData )

  CALL FinalizePositivityLimiter

  CALL FinalizeSlopeLimiter

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeEquationOfState

  CALL FinalizeProgram

END PROGRAM StandingAccretionShock
