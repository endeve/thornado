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
    ComputeRHS_Radiation => NULL(), &
    ComputeExplicitIncrement_Radiation => NULL()

  INTERFACE
    SUBROUTINE ComputeRHS( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )
      USE KindModule, ONLY: DP
      INTEGER, INTENT(in)  :: &
        iX_B0(3), iX_B1(3), iX_E0(3), iX_E1(3)
      REAL(DP), INTENT(in) :: &
        U (1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
      REAL(DP), INTENT(out) :: &
        dU(1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,1:)
    END SUBROUTINE ComputeRHS
  END INTERFACE

  PROCEDURE (ApplyLimiter_Radiation), POINTER, PUBLIC :: &
    ApplySlopeLimiter_Radiation      => NULL(), &
    ApplyPositivityLimiter_Radiation => NULL()

  INTERFACE
    SUBROUTINE ApplyLimiter_Radiation( iX_B0, iX_E0, iX_B1, iX_E1, U )
      USE KindModule, ONLY: DP
      INTEGER,  INTENT(in)    :: &
        iX_B0(3), iX_B1(3), iX_E0(3), iX_E1(3)
      REAL(DP), INTENT(inout) :: &
        U(1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
    END SUBROUTINE ApplyLimiter_Radiation
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
        ComputeExplicitIncrement_Radiation &
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
        ComputeExplicitIncrement_Radiation &
          => ComputeRHS_Dummy
        ApplySlopeLimiter_Radiation &
          => ApplyLimiter_Dummy
        ApplyPositivityLimiter_Radiation &
          => ApplyLimiter_Dummy

    END SELECT

  END SUBROUTINE InitializeRadiationEvolution


  SUBROUTINE FinalizeRadiationEvolution

    NULLIFY( ComputeRHS_Radiation )
    NULLIFY( ComputeExplicitIncrement_Radiation )

  END SUBROUTINE FinalizeRadiationEvolution


  SUBROUTINE ComputeRHS_Dummy( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      U (1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
    REAL(DP), INTENT(out) :: &
      dU(1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,1:)

  END SUBROUTINE ComputeRHS_Dummy


  SUBROUTINE ApplyLimiter_Dummy( iX_B0, iX_E0, iX_B1, iX_E1, U )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      U(1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)

  END SUBROUTINE ApplyLimiter_Dummy


END MODULE RadiationEvolutionModule
