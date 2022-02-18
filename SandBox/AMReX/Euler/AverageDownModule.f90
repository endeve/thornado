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
    amrex_average_down_dg_order1, &
    amrex_average_down_dg_order2, &
    amrex_average_down_dg_order3

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodes

  ! --- Local Modules ---

  USE MF_Euler_ErrorModule, ONLY: &
    DescribeError_Euler_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: AverageDown
  PUBLIC :: AverageDownTo


CONTAINS


  SUBROUTINE AverageDown( MF )

    TYPE(amrex_multifab), INTENT(inout) :: MF(0:amrex_max_level)

    INTEGER :: iLevel, FinestLevel, nComp

    FinestLevel = amrex_get_finest_level()

    DO iLevel = FinestLevel-1, 0, -1

       nComp = MF(iLevel+1) % nComp()

       IF( nNodes .EQ. 1 )THEN

         CALL amrex_average_down_dg_order1 &
                ( MF        (iLevel+1), MF        (iLevel), &
                  amrex_geom(iLevel+1), amrex_geom(iLevel), &
                  1, nComp, amrex_ref_ratio(iLevel))

       ELSE IF( nNodes .EQ. 2 )THEN

         CALL amrex_average_down_dg_order2 &
                ( MF        (iLevel+1), MF        (iLevel), &
                  amrex_geom(iLevel+1), amrex_geom(iLevel), &
                  1, nComp, amrex_ref_ratio(iLevel))

       ELSE IF( nNodes .EQ. 3 )THEN

         CALL amrex_average_down_dg_order3 &
                ( MF        (iLevel+1), MF        (iLevel), &
                  amrex_geom(iLevel+1), amrex_geom(iLevel), &
                  1, nComp, amrex_ref_ratio(iLevel))

       ELSE

        CALL DescribeError_Euler_MF &
               ( 03, Message_Option = '  SUBROUTINE: AverageDown', &
                     Int_Option = [ nNodes ] )

       END IF

    END DO

  END SUBROUTINE AverageDown


  SUBROUTINE AverageDownTo( CoarseLevel, MF )

    INTEGER,              INTENT(IN)    :: CoarseLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:amrex_max_level)

    IF( nNodes .EQ. 1 )THEN

      CALL amrex_average_down_dg_order1 &
             ( MF        (CoarseLevel+1), MF        (CoarseLevel), &
               amrex_geom(CoarseLevel+1), amrex_geom(CoarseLevel), &
               1, 1, amrex_ref_ratio(CoarseLevel))

    ELSE IF( nNodes .EQ. 2 )THEN

      CALL amrex_average_down_dg_order2 &
             ( MF        (CoarseLevel+1), MF        (CoarseLevel), &
               amrex_geom(CoarseLevel+1), amrex_geom(CoarseLevel), &
               1, 1, amrex_ref_ratio(CoarseLevel))

    ELSE IF( nNodes .EQ. 3 )THEN

      CALL amrex_average_down_dg_order3 &
             ( MF        (CoarseLevel+1), MF        (CoarseLevel), &
               amrex_geom(CoarseLevel+1), amrex_geom(CoarseLevel), &
               1, 1, amrex_ref_ratio(CoarseLevel))

    ELSE

      CALL DescribeError_Euler_MF &
             ( 03, Message_Option = '  SUBROUTINE: AverageDownTo', &
                   Int_Option = [ nNodes ] )

    END IF

  END SUBROUTINE AverageDownTo


END MODULE AverageDownModule
