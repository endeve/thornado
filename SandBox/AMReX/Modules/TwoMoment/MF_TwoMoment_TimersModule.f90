MODULE MF_TwoMoment_TimersModule

  ! --- Thornado Modules ---

  USE TwoMoment_TimersModule, ONLY: &
    InitializeTimers, &
    FinalizeTimers

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTimers_AMReX_TwoMoment
  PUBLIC :: FinalizeTimers_AMReX_TwoMoment

  LOGICAL, PUBLIC :: TimeIt_AMReX_TwoMoment = .TRUE.

CONTAINS


  SUBROUTINE InitializeTimers_AMReX_TwoMoment

    IF( .NOT. TimeIt_AMReX_TwoMoment ) RETURN

    CALL InitializeTimers

  END SUBROUTINE InitializeTimers_AMReX_TwoMoment


  SUBROUTINE FinalizeTimers_AMReX_TwoMoment

    IF( .NOT. TimeIt_AMReX_TwoMoment ) RETURN

    CALL FinalizeTimers

  END SUBROUTINE FinalizeTimers_AMReX_TwoMoment


END MODULE MF_TwoMoment_TimersModule
