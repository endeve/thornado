MODULE GravitySolutionModule

  USE KindModule, ONLY: &
    DP
  USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
    InitializeGravitySolver_Newtonian_Poseidon, &
    FinalizeGravitySolver_Newtonian_Poseidon, &
    SolveGravity_Newtonian_Poseidon

  IMPLICIT NONE
  PRIVATE

  CHARACTER(32) :: GravitySolver

  PROCEDURE (GS), POINTER, PUBLIC :: &
    SolveGravity => NULL()

  INTERFACE
    SUBROUTINE GS
    END SUBROUTINE GS
  END INTERFACE

  PUBLIC :: InitializeGravitySolver
  PUBLIC :: FinalizeGravitySolver

CONTAINS


  SUBROUTINE InitializeGravitySolver( GravitySolver_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: GravitySolver_Option

    GravitySolver = 'Dummy'
    IF( PRESENT( GravitySolver_Option ) ) &
      GravitySolver = GravitySolver_Option

    SELECT CASE ( TRIM( GravitySolver ) )

      CASE ( 'Newtonian_Poseidon' )

        CALL InitializeGravitySolver_Newtonian_Poseidon
        SolveGravity &
          => SolveGravity_Newtonian_Poseidon

      CASE DEFAULT

        SolveGravity &
          => SolveGravity_Dummy

    END SELECT

  END SUBROUTINE InitializeGravitySolver


  SUBROUTINE FinalizeGravitySolver

    SELECT CASE ( TRIM( GravitySolver ) )

      CASE ( 'Newtonian_Poseidon' )

        CALL FinalizeGravitySolver_Newtonian_Poseidon

      CASE DEFAULT

    END SELECT

    NULLIFY( SolveGravity )

  END SUBROUTINE FinalizeGravitySolver


  SUBROUTINE SolveGravity_Dummy

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'GravitySolutionModule: SolveGravity_Dummy'
    WRITE(*,*)

  END SUBROUTINE SolveGravity_Dummy


END MODULE GravitySolutionModule
