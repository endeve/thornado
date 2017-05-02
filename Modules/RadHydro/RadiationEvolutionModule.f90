MODULE RadiationEvolutionModule

  USE KindModule, ONLY: &
    DP
  USE MomentEquationsSolutionModule_M1_DG, ONLY: &
    ComputeRHS_M1_DG
  USE MomentEquationsLimiterModule_DG, ONLY: &
    InitializeLimiters_M1_DG
  USE MomentEquationsSlopeLimiterModule_DG, ONLY: &
    ApplySlopeLimiter_M1_DG
  USE MomentEquationsPositivityLimiterModule_DG, ONLY: &
    ApplyPositivityLimiter_M1_DG

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

  PROCEDURE (ApplyLimiter), POINTER, PUBLIC :: &
    ApplySlopeLimiter_Radiation      => NULL(), &
    ApplyPositivityLimiter_Radiation => NULL()

  INTERFACE
    SUBROUTINE ApplyLimiter
    END SUBROUTINE ApplyLimiter
  END INTERFACE

  PUBLIC :: InitializeRadiationEvolution
  PUBLIC :: FinalizeRadiationEvolution

CONTAINS


  SUBROUTINE InitializeRadiationEvolution &
               ( RadiationSolver_Option, ApplySlopeLimiter_Option, &
                 BetaTVB_Option, BetaTVD_Option, ApplyPositivityLimiter_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RadiationSolver_Option
    LOGICAL,          INTENT(in), OPTIONAL :: ApplySlopeLimiter_Option
    REAL(DP),         INTENT(in), OPTIONAL :: BetaTVB_Option
    REAL(DP),         INTENT(in), OPTIONAL :: BetaTVD_Option
    LOGICAL,          INTENT(in), OPTIONAL :: ApplyPositivityLimiter_Option

    IF( PRESENT( RadiationSolver_Option ) )THEN
      RadiationSolver = TRIM( RadiationSolver_Option )
    END IF

    SELECT CASE ( TRIM( RadiationSolver ) )

      CASE ( 'M1_DG' )

        ComputeRHS_Radiation &
          => ComputeRHS_M1_DG
        ApplySlopeLimiter_Radiation &
          => ApplySlopeLimiter_M1_DG
        ApplyPositivityLimiter_Radiation &
          => ApplyPositivityLimiter_M1_DG

        CALL InitializeLimiters_M1_DG &
               ( ApplySlopeLimiter_Option &
                   = ApplySlopeLimiter_Option, &
                 BetaTVB_Option = BetaTVB_Option, &
                 BetaTVD_Option = BetaTVD_Option, &
                 ApplyPositivityLimiter_Option &
                   = ApplyPositivityLimiter_Option )

      CASE DEFAULT

        ComputeRHS_Radiation &
          => ComputeRHS_Dummy
        ApplySlopeLimiter_Radiation &
          => ApplyLimiter_Dummy
        ApplyPositivityLimiter_Radiation &
          => ApplyLimiter_Dummy

    END SELECT

  END SUBROUTINE InitializeRadiationEvolution


  SUBROUTINE FinalizeRadiationEvolution

    NULLIFY( ComputeRHS_Radiation )

  END SUBROUTINE FinalizeRadiationEvolution


  SUBROUTINE ComputeRHS_Dummy( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

  END SUBROUTINE ComputeRHS_Dummy


  SUBROUTINE ApplyLimiter_Dummy

  END SUBROUTINE ApplyLimiter_Dummy


END MODULE RadiationEvolutionModule
