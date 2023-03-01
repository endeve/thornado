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
#if defined( THORNADO_USE_MESHREFINEMENT )
  USE amrex_multifabutil_module, ONLY: &
    amrex_average_down_dg
#endif
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

  INTERFACE AverageDown
    MODULE PROCEDURE AverageDown_Geometry
    MODULE PROCEDURE AverageDown_Fluid
  END INTERFACE AverageDown

  INTERFACE AverageDownTo
    MODULE PROCEDURE AverageDownTo_Geometry
    MODULE PROCEDURE AverageDownTo_Fluid
  END INTERFACE AverageDownTo

  PUBLIC :: AverageDown
  PUBLIC :: AverageDownTo

CONTAINS


  SUBROUTINE AverageDown_Geometry( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    INTEGER :: iLevel, FinestLevel

    FinestLevel = amrex_get_finest_level()

    DO iLevel = FinestLevel-1, 0, -1

      CALL AverageDownTo_Geometry( iLevel, MF_uGF )

    END DO

  END SUBROUTINE AverageDown_Geometry


  SUBROUTINE AverageDown_Fluid( MF_uGF, MF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER :: iLevel, FinestLevel

    FinestLevel = amrex_get_finest_level()

    DO iLevel = FinestLevel-1, 0, -1

      CALL AverageDownTo_Fluid( iLevel, MF_uGF, MF )

    END DO

  END SUBROUTINE AverageDown_Fluid


  SUBROUTINE AverageDownTo_Geometry( CoarseLevel, MF_uGF )

    INTEGER             , INTENT(in)    :: CoarseLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    INTEGER :: nComp, nF, iErr

    nComp = MF_uGF(CoarseLevel) % nComp()
    nF    = nComp / nDOFX

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,'(4x,A,I3.3)') &
          'CALL AverageDownTo_Geometry, CoarseLevel: ', CoarseLevel

      END IF

    END IF

    CALL MultiplyWithMetric &
           ( CoarseLevel+1, MF_uGF, nF, +1, swXX_Option = [ 0, 0, 0 ] )

#if defined( THORNADO_USE_MESHREFINEMENT )

    CALL amrex_average_down_dg &
           ( MF_uGF    (CoarseLevel+1), MF_uGF    (CoarseLevel), &
             amrex_geom(CoarseLevel+1), amrex_geom(CoarseLevel), &
             1, nComp, amrex_ref_ratio(CoarseLevel))

#endif

    CALL MultiplyWithMetric &
           ( CoarseLevel+1, MF_uGF, nF, -1, swXX_Option = [ 0, 0, 0 ] )

    CALL MultiplyWithMetric &
           ( CoarseLevel  , MF_uGF, nF, -1, swXX_Option = [ 0, 0, 0 ] )

  END SUBROUTINE AverageDownTo_Geometry


  SUBROUTINE AverageDownTo_Fluid( CoarseLevel, MF_uGF, MF )

    INTEGER             , INTENT(in)    :: CoarseLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    TYPE(amrex_multifab) :: SqrtGm(CoarseLevel:CoarseLevel+1)

    INTEGER :: nComp, nF, iErr

    nComp = MF(CoarseLevel) % nComp()
    nF    =  nComp / nDOFX

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,'(4x,A,I3.3)') &
          'CALL AverageDownTo_Fluid, CoarseLevel: ', CoarseLevel

      END IF

    END IF

    CALL amrex_multifab_build &
           ( SqrtGm(CoarseLevel+1), MF_uGF(CoarseLevel+1) % BA, &
                                    MF_uGF(CoarseLevel+1) % DM, nDOFX, swX )

    CALL SqrtGm(CoarseLevel+1) % COPY &
           ( MF_uGF(CoarseLevel+1), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

    CALL MultiplyWithMetric &
           ( CoarseLevel+1, SqrtGm(CoarseLevel+1), MF, nF, +1, &
             swXX_Option = [ 0, 0, 0 ] )

#if defined( THORNADO_USE_MESHREFINEMENT )

    CALL amrex_average_down_dg &
           ( MF        (CoarseLevel+1), MF        (CoarseLevel), &
             amrex_geom(CoarseLevel+1), amrex_geom(CoarseLevel), &
             1, nComp, amrex_ref_ratio(CoarseLevel))

#endif

    CALL MultiplyWithMetric &
           ( CoarseLevel+1, SqrtGm(CoarseLevel+1), MF, nF, -1, &
             swXX_Option = [ 0, 0, 0 ] )

    CALL amrex_multifab_destroy( SqrtGm(CoarseLevel+1) )

    CALL amrex_multifab_build &
           ( SqrtGm(CoarseLevel), MF_uGF(CoarseLevel) % BA, &
                                  MF_uGF(CoarseLevel) % DM, nDOFX, swX )

    CALL SqrtGm(CoarseLevel) % COPY &
           ( MF_uGF(CoarseLevel), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

    CALL MultiplyWithMetric &
           ( CoarseLevel, SqrtGm(CoarseLevel), MF, nF, -1, &
             swXX_Option = [ 0, 0, 0 ] )

    CALL amrex_multifab_destroy( SqrtGm(CoarseLevel) )

  END SUBROUTINE AverageDownTo_Fluid


END MODULE AverageDownModule
