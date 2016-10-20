MODULE FluidEvolutionModule

  USE KindModule, ONLY: &
    DP
  USE EulerEquationsSolutionModule_DG, ONLY: &
    ComputeRHS_Euler_DG
  USE EulerEquationsLimiterModule_DG, ONLY: &
    InitializeLimiters_Euler_DG, &
    ApplySlopeLimiter_Euler_DG, &
    ApplyPositivityLimiter_Euler_DG

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

  PROCEDURE (ApplyLimiter), POINTER, PUBLIC :: &
    ApplySlopeLimiter_Fluid      => NULL(), &
    ApplyPositivityLimiter_Fluid => NULL()

  INTERFACE
    SUBROUTINE ApplyLimiter
    END SUBROUTINE ApplyLimiter
  END INTERFACE

  PUBLIC :: InitializeFluidEvolution
  PUBLIC :: FinalizeFluidEvolution

CONTAINS


  SUBROUTINE InitializeFluidEvolution( FluidSolver_Option, BetaTVB_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: FluidSolver_Option
    REAL(DP),         INTENT(in), OPTIONAL :: BetaTVB_Option

    IF( PRESENT( FluidSolver_Option ) )THEN
      FluidSolver = FluidSolver_Option
    END IF

    SELECT CASE ( TRIM( FluidSolver ) )
      CASE( 'Euler_DG' )
        ComputeRHS_Fluid &
          => ComputeRHS_Euler_DG
        ApplySlopeLimiter_Fluid &
          => ApplySlopeLimiter_Euler_DG
        ApplyPositivityLimiter_Fluid &
          => ApplyPositivityLimiter_Euler_DG
        CALL InitializeLimiters_Euler_DG &
               ( BetaTVB_Option = BetaTVB_Option )
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

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'FluidEvolutionModule: ComputeRHS_Dummy'
    WRITE(*,*)

  END SUBROUTINE ComputeRHS_Dummy


  SUBROUTINE ApplyLimiter_Dummy

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'FluidEvolutionModule: ApplyLimiter_Dummy'
    WRITE(*,*)

  END SUBROUTINE ApplyLimiter_Dummy


END MODULE FluidEvolutionModule
