PROGRAM ProgramInitializationTest

  USE KindModule, ONLY: &
    DP
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option = 'ProgramInitializationTest', &
           nX_Option = [ 16, 1, 1 ], nE_Option = 10, nNodes_Option = 2 )

  CALL FinalizeProgram

END PROGRAM ProgramInitializationTest
