MODULE AverageDownModule_MHD

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_get_finest_level
  USE amrex_amr_module, ONLY: &
    amrex_ref_ratio
#if defined( THORNADO_USE_MESHREFINEMENT )
  USE thornado_amrex_multifabutil_module, ONLY: &
    amrex_average_down_dg_conservative, &
    amrex_average_down_dg_pointwise, &
    amrex_average_down_cg
#endif
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator, &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    swX
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm
  USE MHD_MeshRefinementModule, ONLY: &
    nFine, &
    pG2L_c, &
    pL2G_c, &
    pF2C_c, &
    vpFineToCoarseProjectionMatrix

  ! --- Local Modules ---

  USE MF_GeometryModule, ONLY: &
    UpdateSpatialMetric_MF
  USE InputParsingModule, ONLY: &
    DEBUG
  USE MF_TimersModule_MHD, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_AverageDown

  IMPLICIT NONE
  PRIVATE

  INTERFACE AverageDown
    MODULE PROCEDURE AverageDown_PointWise
    MODULE PROCEDURE AverageDown_Conservative
  END INTERFACE AverageDown

  INTERFACE AverageDownTo
    MODULE PROCEDURE AverageDownTo_PointWise
    MODULE PROCEDURE AverageDownTo_Conservative
  END INTERFACE AverageDownTo

  PUBLIC :: AverageDown
  PUBLIC :: AverageDownTo

CONTAINS


  SUBROUTINE AverageDown_PointWise( MF, UpdateSpatialMetric_Option )

    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)
    LOGICAL             , INTENT(in), OPTIONAL :: &
      UpdateSpatialMetric_Option

    INTEGER :: iLevel, FinestLevel
    LOGICAL :: UpdateSpatialMetric

    UpdateSpatialMetric = .FALSE.
    IF( PRESENT( UpdateSpatialMetric_Option ) ) &
      UpdateSpatialMetric = UpdateSpatialMetric_Option

    FinestLevel = amrex_get_finest_level()

    DO iLevel = FinestLevel-1, 0, -1

      CALL AverageDownTo_PointWise &
             ( iLevel, MF, UpdateSpatialMetric_Option = UpdateSpatialMetric )

    END DO

  END SUBROUTINE AverageDown_PointWise


  SUBROUTINE AverageDown_Conservative( MF_uGF, MF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER :: iLevel, FinestLevel

    FinestLevel = amrex_get_finest_level()

    DO iLevel = FinestLevel-1, 0, -1

      CALL AverageDownTo_Conservative( iLevel, MF_uGF, MF )

    END DO

  END SUBROUTINE AverageDown_Conservative


  SUBROUTINE AverageDownTo_PointWise &
    ( CoarseLevel, MF, UpdateSpatialMetric_Option )

    INTEGER             , INTENT(in)    :: CoarseLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)
    LOGICAL             , INTENT(in), OPTIONAL :: &
      UpdateSpatialMetric_Option

    INTEGER :: iErr
    LOGICAL :: UpdateSpatialMetric

    CALL TimersStart_AMReX( Timer_AMReX_AverageDown )

#if defined( THORNADO_USE_MESHREFINEMENT )

    UpdateSpatialMetric = .FALSE.
    IF( PRESENT( UpdateSpatialMetric_Option ) ) &
      UpdateSpatialMetric = UpdateSpatialMetric_Option

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,'(4x,A,I3.3)') &
          'CALL AverageDownTo_PointWise, CoarseLevel: ', CoarseLevel

      END IF

    END IF

    CALL amrex_average_down_dg_pointwise &
           ( MF(CoarseLevel+1), MF(CoarseLevel), &
             MF(CoarseLevel) % nComp(), amrex_ref_ratio(CoarseLevel), &
             nDOFX, nFine, vpFineToCoarseProjectionMatrix )

    IF( UpdateSpatialMetric )THEN

      CALL UpdateSpatialMetric_MF( CoarseLevel, MF(CoarseLevel) )

    END IF

#endif

    CALL TimersStop_AMReX( Timer_AMReX_AverageDown )

  END SUBROUTINE AverageDownTo_PointWise


  SUBROUTINE AverageDownTo_Conservative( CoarseLevel, MF_uGF, MF )

    INTEGER             , INTENT(in)    :: CoarseLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    TYPE(amrex_multifab) :: SqrtGm(CoarseLevel:CoarseLevel+1)

    INTEGER :: iErr

    CALL TimersStart_AMReX( Timer_AMReX_AverageDown )

#if defined( THORNADO_USE_MESHREFINEMENT )

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,'(4x,A,I3.3)') &
          'CALL AverageDownTo_Conservative, CoarseLevel: ', CoarseLevel

      END IF

    END IF

    CALL amrex_multifab_build &
           ( SqrtGm(CoarseLevel+1), MF(CoarseLevel+1) % BA, &
                                    MF(CoarseLevel+1) % DM, nDOFX, swX )

    CALL SqrtGm(CoarseLevel+1) % COPY &
           ( MF_uGF(CoarseLevel+1), nDOFX*(iGF_SqrtGm-1)+1, 1, nDOFX, swX )

    CALL amrex_multifab_build &
           ( SqrtGm(CoarseLevel  ), MF(CoarseLevel  ) % BA, &
                                    MF(CoarseLevel  ) % DM, nDOFX, swX )

    CALL SqrtGm(CoarseLevel  ) % COPY &
           ( MF_uGF(CoarseLevel  ), nDOFX*(iGF_SqrtGm-1)+1, 1, nDOFX, swX )

    CALL amrex_average_down_dg_conservative &
           ( MF    (CoarseLevel+1), MF    (CoarseLevel), &
             SqrtGm(CoarseLevel+1), SqrtGm(CoarseLevel), &
             MF(CoarseLevel) % nComp(), amrex_ref_ratio(CoarseLevel), &
             nDOFX, nFine, vpFineToCoarseProjectionMatrix )

    CALL amrex_multifab_destroy( SqrtGm(CoarseLevel+1) )

    CALL amrex_multifab_destroy( SqrtGm(CoarseLevel) )

#endif

    CALL TimersStop_AMReX( Timer_AMReX_AverageDown )

  END SUBROUTINE AverageDownTo_Conservative


END MODULE AverageDownModule_MHD
