MODULE TwoMoment_TallyModule_FMC

  USE KindModule, ONLY: &
    DP, Zero
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nCM

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally
  PUBLIC :: ComputeTally
  PUBLIC :: IncrementOffGridTally_TwoMoment
  PUBLIC :: FinalizeTally

CONTAINS


  SUBROUTINE InitializeTally

  END SUBROUTINE InitializeTally


  SUBROUTINE ComputeTally( SetInitialValues_Option )

    LOGICAL, INTENT(in), OPTIONAL :: SetInitialValues_Option

    CALL ComputeTally_TwoMoment

  END SUBROUTINE ComputeTally


  SUBROUTINE ComputeTally_TwoMoment

  END SUBROUTINE ComputeTally_TwoMoment


  SUBROUTINE IncrementOffGridTally_TwoMoment( dM )

    REAL(DP), INTENT(in) :: dM(nCM)

  END SUBROUTINE IncrementOffGridTally_TwoMoment


  SUBROUTINE FinalizeTally

  END SUBROUTINE FinalizeTally


END MODULE TwoMoment_TallyModule_FMC
