MODULE MF_MetricInitializationModule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  USE MF_MetricInitializationModule_XCFC_Poseidon_MHD, ONLY: &
    InitializeMetric_XCFC_MF_Poseidon, &
    InitializeMetricFromCheckpoint_XCFC_MF_Poseidon

#endif

  USE MF_KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMetric_MF
  PUBLIC :: InitializeMetricFromCheckpoint_MF

CONTAINS


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE InitializeMetric_MF &
    ( MF_uGF, MF_uCM, MF_uCR, MF_uPM, MF_uAM, TOLERANCE_option )

#else

  SUBROUTINE InitializeMetric_MF &
    ( MF_uGF, MF_uCM, MF_uPM, MF_uAM, TOLERANCE_Option )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
#endif
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPM(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uAM(0:)
    REAL(DP)            , INTENT(in), OPTIONAL :: TOLERANCE_Option

    REAL(DP) :: TOLERANCE

    TOLERANCE = 1.0e-13_DP
    IF( PRESENT( TOLERANCE_Option ) ) &
      TOLERANCE = TOLERANCE_Option

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

#ifndef THORNADO_NOTRANSPORT

    CALL InitializeMetric_XCFC_MF_Poseidon &
           ( MF_uGF, MF_uCM, MF_uCR, MF_uPM, MF_uAM )

#else

    CALL InitializeMetric_XCFC_MF_Poseidon &
           ( MF_uGF, MF_uCM, MF_uPM, MF_uAM, TOLERANCE_Option = TOLERANCE )

#endif

#endif

  END SUBROUTINE InitializeMetric_MF


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE InitializeMetricFromCheckpoint_MF( MF_uGF, MF_uCM, MF_uCR )

#else

  SUBROUTINE InitializeMetricFromCheckpoint_MF( MF_uGF, MF_uCM )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCM(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR(0:)
#endif

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

#ifndef THORNADO_NOTRANSPORT

    CALL InitializeMetricFromCheckpoint_XCFC_MF_Poseidon &
           ( MF_uGF, MF_uCM, MF_uCR )

#else

    CALL InitializeMetricFromCheckpoint_XCFC_MF_Poseidon &
           ( MF_uGF, MF_uCM )

#endif

#endif

  END SUBROUTINE InitializeMetricFromCheckpoint_MF


END MODULE MF_MetricInitializationModule
