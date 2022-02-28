MODULE AverageDownModule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab
  USE amrex_amrcore_module, ONLY: &
    amrex_max_level, &
    amrex_get_finest_level
  USE amrex_amr_module, ONLY: &
    amrex_geom, &
    amrex_ref_ratio
  USE amrex_multifabutil_module, ONLY: &
    amrex_average_down_dg

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: AverageDown
  PUBLIC :: AverageDownTo


CONTAINS


  SUBROUTINE AverageDown( MF )

    TYPE(amrex_multifab), INTENT(inout) :: MF(0:amrex_max_level)

    INTEGER :: iLevel, FinestLevel

    FinestLevel = amrex_get_finest_level()

    DO iLevel = FinestLevel-1, 0, -1

      CALL AverageDownTo( iLevel, MF )

    END DO

  END SUBROUTINE AverageDown


  SUBROUTINE AverageDownTo( CoarseLevel, MF )

    INTEGER,              INTENT(IN)    :: CoarseLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:amrex_max_level)

    INTEGER :: nComp

    nComp = MF(CoarseLevel) % nComp()

    CALL amrex_average_down_dg &
           ( MF        (CoarseLevel+1), MF        (CoarseLevel), &
             amrex_geom(CoarseLevel+1), amrex_geom(CoarseLevel), &
             1, nComp, amrex_ref_ratio(CoarseLevel))

  END SUBROUTINE AverageDownTo


END MODULE AverageDownModule
