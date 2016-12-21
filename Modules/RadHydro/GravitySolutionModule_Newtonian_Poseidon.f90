MODULE GravitySolutionModule_Newtonian_Poseidon

  USE KindModule, ONLY: &
    DP, Pi
  USE GravitySolutionUtilitiesModule, ONLY: &
    ComputeTotalBaryonMass

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

    REAL(DP) :: BaryonMass

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'SolveGravity_Newtonian_Poseidon'
    WRITE(*,*)

    CALL ComputeTotalBaryonMass( BaryonMass )

  END SUBROUTINE SolveGravity_Newtonian_Poseidon


END MODULE GravitySolutionModule_Newtonian_Poseidon
