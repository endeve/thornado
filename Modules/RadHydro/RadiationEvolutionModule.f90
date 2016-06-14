MODULE RadiationEvolutionModule

  USE MomentEquationsSolutionModule_M1_DG, ONLY: &
    ComputeRHS_M1_DG

  IMPLICIT NONE
  PRIVATE

  CHARACTER(32) :: RadiationSolver = 'Dummy'

  PROCEDURE (ComputeRHS), POINTER, PUBLIC :: &
    ComputeRHS_Radiation => NULL()

  INTERFACE
    SUBROUTINE ComputeRHS( iX_Begin, iX_End )
      INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    END SUBROUTINE ComputeRHS
  END INTERFACE

  PUBLIC :: InitializeRadiationEvolution
  PUBLIC :: FinalizeRadiationEvolution

CONTAINS


  SUBROUTINE InitializeRadiationEvolution( RadiationSolver_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RadiationSolver_Option

    IF( PRESENT( RadiationSolver_Option ) )THEN
      RadiationSolver = TRIM( RadiationSolver_Option )
    END IF

    SELECT CASE ( TRIM( RadiationSolver ) )
      CASE ( 'M1_DG' )
        ComputeRHS_Radiation &
          => ComputeRHS_M1_DG
      CASE DEFAULT
        ComputeRHS_Radiation &
          => ComputeRHS_Dummy
    END SELECT

  END SUBROUTINE InitializeRadiationEvolution


  SUBROUTINE FinalizeRadiationEvolution

    NULLIFY( ComputeRHS_Radiation )

  END SUBROUTINE FinalizeRadiationEvolution


  SUBROUTINE ComputeRHS_Dummy( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'RadiationEvolutionModule: ComputeRHS_Dummy'
    WRITE(*,*)

  END SUBROUTINE ComputeRHS_Dummy


END MODULE RadiationEvolutionModule
