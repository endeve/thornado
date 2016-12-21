MODULE GravitySolutionModule_Newtonian_Poseidon

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_Newtonian_Poseidon
  PUBLIC :: FinalizeGravitySolver_Newtonian_Poseidon
  PUBLIC :: SolveGravity_Newtonian_Poseidon

CONTAINS


  SUBROUTINE InitializeGravitySolver_Newtonian_Poseidon

  END SUBROUTINE InitializeGravitySolver_Newtonian_Poseidon


  SUBROUTINE FinalizeGravitySolver_Newtonian_Poseidon

  END SUBROUTINE FinalizeGravitySolver_Newtonian_Poseidon


  SUBROUTINE SolveGravity_Newtonian_Poseidon

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'SolveGravity_Newtonian_Poseidon'
    WRITE(*,*)

  END SUBROUTINE SolveGravity_Newtonian_Poseidon


END MODULE GravitySolutionModule_Newtonian_Poseidon
