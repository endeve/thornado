PROGRAM ScalarWaveTest

  USE KindModule, ONLY: &
    DP
  USE ScalarWave_ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE InitializationModule, ONLY: &
    InitializeFields
  USE ScalarWave_InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF  
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeScalarWave_SSPRK, &
    FinalizeScalarWave_SSPRK, &
    UpdateScalarWave_SSPRK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option = 'SpacetimeTest', &
           nX_Option = [ 16, 1, 1 ], swX_Option = [ 1, 0, 0 ], nNodes_Option = 2 )

  Call InitializeFields

  CALL WriteFieldsHDF( 0.0_DP, WriteSF_Option=.true.)

  CALL FinalizeProgram

END PROGRAM ScalarWaveTest
