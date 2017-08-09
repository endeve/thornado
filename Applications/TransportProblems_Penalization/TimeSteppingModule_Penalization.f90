MODULE TimeSteppingModule_Penalization

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: EvolveFields

CONTAINS


  SUBROUTINE EvolveFields( t_begin, t_end, dt_write, dt_fixed_Option )

    REAL(DP), INTENT(in) :: t_begin, t_end, dt_write
    REAL(DP), INTENT(in), OPTIONAL :: dt_fixed_Option

  END SUBROUTINE EvolveFields


END MODULE TimeSteppingModule_Penalization
