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
    amrex_average_down_dg_order2

  ! --- Local Modules ---


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: AverageDown
  PUBLIC :: AverageDownTo


CONTAINS


  SUBROUTINE AverageDown( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:amrex_max_level)

    INTEGER :: iLevel, FinestLevel, nCompCF, nCompGF

    FinestLevel = amrex_get_finest_level()

    DO iLevel = FinestLevel-1, 0, -1

       nCompGF = MF_uGF(iLevel+1) % nComp()
       nCompCF = MF_uCF(iLevel+1) % nComp()

       CALL amrex_average_down_dg_order2 &
              ( MF_uGF    (iLevel+1), MF_uGF    (iLevel), &
                amrex_geom(iLevel+1), amrex_geom(iLevel), &
                1, nCompGF, amrex_ref_ratio(iLevel))

       CALL amrex_average_down_dg_order2 &
              ( MF_uCF    (iLevel+1), MF_uCF    (iLevel), &
                amrex_geom(iLevel+1), amrex_geom(iLevel), &
                1, nCompCF, amrex_ref_ratio(iLevel))

    END DO

  END SUBROUTINE AverageDown


  SUBROUTINE AverageDownTo( CoarseLevel, MF_uCF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:amrex_max_level)

    INTEGER, INTENT(IN) :: CoarseLevel

    CALL amrex_average_down_dg_order2 &
           ( MF_uCF    (CoarseLevel+1), MF_uCF    (CoarseLevel), &
             amrex_geom(CoarseLevel+1), amrex_geom(CoarseLevel), &
             1, 1, amrex_ref_ratio(CoarseLevel))

  END SUBROUTINE AverageDownTo


END MODULE AverageDownModule
