MODULE FluidEvolutionModule

  USE KindModule, ONLY: &
    DP
  USE EquationOfStateModule, ONLY: &
    EquationOfState
  USE EulerEquationsSolutionModule_DG, ONLY: &
    ComputeRHS_Euler_DG
  USE EulerEquationsLimiterModule_DG, ONLY: &
    InitializeLimiters_Euler_DG
  USE EulerEquationsSlopeLimiterModule_DG_IDEAL, ONLY: &
    ApplySlopeLimiter_Euler_DG_IDEAL
  USE EulerEquationsPositivityLimiterModule_DG_IDEAL, ONLY: &
    ApplyPositivityLimiter_Euler_DG_IDEAL
  USE EulerEquationsSlopeLimiterModule_DG_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_DG_TABLE
  USE EulerEquationsPositivityLimiterModule_DG_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_DG_TABLE
  USE EulerEquationsSolutionModule_DG_GR, ONLY: &
    ComputeRHS_Euler_DG_GR

  IMPLICIT NONE
  PRIVATE

  CHARACTER(32) :: FluidSolver = 'Dummy'

  PROCEDURE (ComputeRHS), POINTER, PUBLIC :: &
    ComputeRHS_Fluid => NULL()

  INTERFACE
    SUBROUTINE ComputeRHS( iX_Begin, iX_End )
      INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    END SUBROUTINE ComputeRHS
  END INTERFACE

  PROCEDURE (ApplyLimiter_Fluid), POINTER, PUBLIC :: &
    ApplySlopeLimiter_Fluid      => NULL(), &
    ApplyPositivityLimiter_Fluid => NULL()

  INTERFACE
    SUBROUTINE ApplyLimiter_Fluid
    END SUBROUTINE ApplyLimiter_Fluid
  END INTERFACE

  PUBLIC :: InitializeFluidEvolution
  PUBLIC :: FinalizeFluidEvolution

CONTAINS


  SUBROUTINE InitializeFluidEvolution &
               ( FluidSolver_Option, &
                 ApplySlopeLimiter_Option, BetaTVB_Option, BetaTVD_Option, &
                 ApplyPositivityLimiter_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: FluidSolver_Option
    LOGICAL,          INTENT(in), OPTIONAL :: ApplySlopeLimiter_Option
    REAL(DP),         INTENT(in), OPTIONAL :: BetaTVB_Option
    REAL(DP),         INTENT(in), OPTIONAL :: BetaTVD_Option
    LOGICAL,          INTENT(in), OPTIONAL :: ApplyPositivityLimiter_Option

    IF( PRESENT( FluidSolver_Option ) )THEN
      FluidSolver = FluidSolver_Option
    END IF

    SELECT CASE ( TRIM( FluidSolver ) )

      CASE( 'Euler_DG' )

        ComputeRHS_Fluid &
          => ComputeRHS_Euler_DG

        SELECT CASE ( TRIM( EquationOfState ) )

          CASE( 'IDEAL' )

            ApplySlopeLimiter_Fluid &
              => ApplySlopeLimiter_Euler_DG_IDEAL
            ApplyPositivityLimiter_Fluid &
              => ApplyPositivityLimiter_Euler_DG_IDEAL

          CASE( 'TABLE' )

            ApplySlopeLimiter_Fluid &
              => ApplySlopeLimiter_Euler_DG_TABLE
            ApplyPositivityLimiter_Fluid &
              => ApplyPositivityLimiter_Euler_DG_TABLE

        END SELECT

        CALL InitializeLimiters_Euler_DG &
               ( ApplySlopeLimiter_Option &
                   = ApplySlopeLimiter_Option, &
                 BetaTVB_Option = BetaTVB_Option, &
                 BetaTVD_Option = BetaTVD_Option, &
                 ApplyPositivityLimiter_Option &
                   = ApplyPositivityLimiter_Option )

      CASE ( 'Euler_DG_GR' )

        ComputeRHS_Fluid &
          => ComputeRHS_Euler_DG_GR
        ApplySlopeLimiter_Fluid &
          => ApplyLimiter_Dummy
        ApplyPositivityLimiter_Fluid &
          => ApplyLimiter_Dummy

      CASE DEFAULT

        ComputeRHS_Fluid &
          => ComputeRHS_Dummy
        ApplySlopeLimiter_Fluid &
          => ApplyLimiter_Dummy
        ApplyPositivityLimiter_Fluid &
          => ApplyLimiter_Dummy

    END SELECT

  END SUBROUTINE InitializeFluidEvolution


  SUBROUTINE FinalizeFluidEvolution

    NULLIFY( ComputeRHS_Fluid )
    NULLIFY( ApplySlopeLimiter_Fluid )
    NULLIFY( ApplyPositivityLimiter_Fluid )

  END SUBROUTINE FinalizeFluidEvolution


  SUBROUTINE ComputeRHS_Dummy( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    RETURN

  END SUBROUTINE ComputeRHS_Dummy


  SUBROUTINE ApplyLimiter_Dummy

    RETURN

  END SUBROUTINE ApplyLimiter_Dummy


END MODULE FluidEvolutionModule
