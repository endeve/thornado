MODULE AverageDownModule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_get_finest_level
  USE amrex_amr_module, ONLY: &
    amrex_geom, &
    amrex_ref_ratio
  USE amrex_multifabutil_module, ONLY: &
    amrex_average_down_dg
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator, &
    amrex_parallel_ioprocessor

  USE MF_UtilitiesModule, ONLY: &
    MultiplyWithMetric
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm

  USE InputParsingModule, ONLY: &
    swX, &
    DEBUG

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: AverageDown
  PUBLIC :: AverageDownTo

CONTAINS


  SUBROUTINE AverageDown( MF_uGF, MF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER :: iLevel, FinestLevel

    FinestLevel = amrex_get_finest_level()

    DO iLevel = FinestLevel-1, 0, -1

      CALL AverageDownTo( iLevel, MF_uGF, MF )

    END DO

  END SUBROUTINE AverageDown


  SUBROUTINE AverageDownTo( CoarseLevel, MF_uGF, MF )

    INTEGER,              INTENT(IN)    :: CoarseLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    TYPE(amrex_multifab) :: SqrtGm(CoarseLevel:CoarseLevel+1)

    INTEGER :: nComp, nF, iErr

    nComp = MF(CoarseLevel) % nComp()

    nF = MF(0) % nComp() / nDOFX

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,'(4x,A,I3.3)') &
          'CALL AverageDownTo, CoarseLevel: ', CoarseLevel

      END IF

    END IF

    CALL amrex_multifab_build &
           ( SqrtGm(CoarseLevel  ), MF_uGF(CoarseLevel  ) % BA, &
                                    MF_uGF(CoarseLevel  ) % DM, nDOFX, swX )
    CALL SqrtGm(CoarseLevel) % COPY &
           ( MF_uGF(CoarseLevel  ), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

    CALL amrex_multifab_build &
           ( SqrtGm(CoarseLevel+1), MF_uGF(CoarseLevel+1) % BA, &
                                    MF_uGF(CoarseLevel+1) % DM, nDOFX, swX )
    CALL SqrtGm(CoarseLevel+1) % COPY &
           ( MF_uGF(CoarseLevel+1), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

    CALL MultiplyWithMetric( SqrtGm(CoarseLevel  ), MF(CoarseLevel  ), nF, +1 )
    CALL MultiplyWithMetric( SqrtGm(CoarseLevel+1), MF(CoarseLevel+1), nF, +1 )

    CALL amrex_average_down_dg &
           ( MF        (CoarseLevel+1), MF        (CoarseLevel), &
             amrex_geom(CoarseLevel+1), amrex_geom(CoarseLevel), &
             1, nComp, amrex_ref_ratio(CoarseLevel))

    CALL MultiplyWithMetric( SqrtGm(CoarseLevel+1), MF(CoarseLevel+1), nF, -1 )
    CALL MultiplyWithMetric( SqrtGm(CoarseLevel  ), MF(CoarseLevel  ), nF, -1 )

    CALL amrex_multifab_destroy( SqrtGm(CoarseLevel+1) )
    CALL amrex_multifab_destroy( SqrtGm(CoarseLevel  ) )

  END SUBROUTINE AverageDownTo


END MODULE AverageDownModule
