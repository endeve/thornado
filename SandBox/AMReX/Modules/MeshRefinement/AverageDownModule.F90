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

  ! --- thornado Modules ---


  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm
  USE Euler_MeshRefinementModule, ONLY: &
    pProjectionMatrix_T_c, &
    pWeightsX_q_c, &
    VolumeRatio, &
    nFine
  USE InputParsingModule, ONLY: &
    nLevels, &
    swX, &
    DEBUG
  USE MF_TimersModule, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_AverageDown

  ! --- Local Modules ---

  USE MF_UtilitiesModule, ONLY: &
    MultiplyWithMetric
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler_MF

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

  INTEGER, PARAMETER :: swXX(3) = [ 0, 0, 0 ]

CONTAINS


  SUBROUTINE AverageDown_Geometry( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    INTEGER :: iLevel, FinestLevel

    FinestLevel = amrex_get_finest_level()

    DO iLevel = FinestLevel-1, 0, -1

      CALL AverageDownTo_Geometry( iLevel, MF_uGF )

    END DO

  END SUBROUTINE AverageDown_Geometry


  SUBROUTINE AverageDown_Fluid &
    ( MF_uGF, MF, MF_uDF, ApplyPositivityLimiter_Option )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)
    TYPE(amrex_multifab), INTENT(inout), OPTIONAL :: &
      MF_uDF(0:)
    LOGICAL             , INTENT(in)   , OPTIONAL :: &
      ApplyPositivityLimiter_Option

    INTEGER :: iLevel, FinestLevel

    LOGICAL :: ApplyPositivityLimiter

    ApplyPositivityLimiter = .FALSE.
    IF( PRESENT( ApplyPositivityLimiter_Option ) ) &
      ApplyPositivityLimiter = ApplyPositivityLimiter_Option

    FinestLevel = amrex_get_finest_level()

    DO iLevel = FinestLevel-1, 0, -1

      IF( ApplyPositivityLimiter )THEN

        CALL AverageDownTo_Fluid &
               ( iLevel, MF_uGF, MF, &
                 MF_uDF, &
                 ApplyPositivityLimiter_Option = ApplyPositivityLimiter )

      ELSE

        CALL AverageDownTo_Fluid( iLevel, MF_uGF, MF )

      END IF

    END DO

  END SUBROUTINE AverageDown_Fluid


  SUBROUTINE AverageDownTo_Geometry( CoarseLevel, MF_uGF )

    INTEGER             , INTENT(in)    :: CoarseLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    INTEGER :: nComp, nF, iErr

    CALL TimersStart_AMReX( Timer_AMReX_AverageDown )

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
           ( CoarseLevel+1, MF_uGF, nF, +1, swXX_Option = swXX )

#if defined( THORNADO_USE_MESHREFINEMENT )

    CALL amrex_average_down_dg &
           ( MF_uGF    (CoarseLevel+1), MF_uGF    (CoarseLevel), &
             amrex_geom(CoarseLevel+1), amrex_geom(CoarseLevel), &
             1, nComp, amrex_ref_ratio(CoarseLevel), &
            nDOFX, nFine, VolumeRatio, pProjectionMatrix_T_c, pWeightsX_q_c )

#endif

    CALL MultiplyWithMetric &
           ( CoarseLevel+1, MF_uGF, nF, -1, swXX_Option = swXX )

    CALL MultiplyWithMetric &
           ( CoarseLevel  , MF_uGF, nF, -1, swXX_Option = swXX )

    CALL TimersStop_AMReX( Timer_AMReX_AverageDown )

  END SUBROUTINE AverageDownTo_Geometry


  SUBROUTINE AverageDownTo_Fluid &
    ( CoarseLevel, MF_uGF, MF, MF_uDF, ApplyPositivityLimiter_Option )

    INTEGER             , INTENT(in)    :: CoarseLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)
    TYPE(amrex_multifab), INTENT(inout), OPTIONAL :: &
      MF_uDF(0:)
    LOGICAL             , INTENT(in)   , OPTIONAL :: &
      ApplyPositivityLimiter_Option

    TYPE(amrex_multifab) :: SqrtGm(CoarseLevel:CoarseLevel+1)

    INTEGER :: nComp, nF, iErr

    LOGICAL :: ApplyPositivityLimiter

    CALL TimersStart_AMReX( Timer_AMReX_AverageDown )

    ApplyPositivityLimiter = .FALSE.
    IF( PRESENT( ApplyPositivityLimiter_Option ) ) &
      ApplyPositivityLimiter = ApplyPositivityLimiter_Option

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
             swXX_Option = swXX )

#if defined( THORNADO_USE_MESHREFINEMENT )

    CALL amrex_average_down_dg &
           ( MF        (CoarseLevel+1), MF        (CoarseLevel), &
             amrex_geom(CoarseLevel+1), amrex_geom(CoarseLevel), &
             1, nComp, amrex_ref_ratio(CoarseLevel), &
             nDOFX, nFine, VolumeRatio, pProjectionMatrix_T_c, pWeightsX_q_c )

#endif

    CALL MultiplyWithMetric &
           ( CoarseLevel+1, SqrtGm(CoarseLevel+1), MF, nF, -1, &
             swXX_Option = swXX )

    CALL amrex_multifab_destroy( SqrtGm(CoarseLevel+1) )

    CALL amrex_multifab_build &
           ( SqrtGm(CoarseLevel), MF_uGF(CoarseLevel) % BA, &
                                  MF_uGF(CoarseLevel) % DM, nDOFX, swX )

    CALL SqrtGm(CoarseLevel) % COPY &
           ( MF_uGF(CoarseLevel), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

    CALL MultiplyWithMetric &
           ( CoarseLevel, SqrtGm(CoarseLevel), MF, nF, -1, &
             swXX_Option = swXX )

    CALL amrex_multifab_destroy( SqrtGm(CoarseLevel) )

    IF( ApplyPositivityLimiter ) &
      CALL ApplyPositivityLimiter_Euler_MF &
             ( CoarseLevel, &
               MF_uGF(CoarseLevel), MF(CoarseLevel), MF_uDF(CoarseLevel) )

    CALL TimersStop_AMReX( Timer_AMReX_AverageDown )

  END SUBROUTINE AverageDownTo_Fluid


END MODULE AverageDownModule
