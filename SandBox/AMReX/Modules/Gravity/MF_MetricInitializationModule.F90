MODULE MF_MetricInitializationModule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  USE MF_GravitySolutionModule, ONLY: &
    EvolveGravity

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  USE MF_MetricInitializationModule_XCFC_Poseidon, ONLY: &
    InitializeMetric_XCFC_MF_Poseidon, &
    InitializeMetricFromCheckpoint_XCFC_MF_Poseidon

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMetric_MF
  PUBLIC :: InitializeMetricFromCheckpoint_MF

CONTAINS


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE InitializeMetric_MF &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uPF, MF_uAF )

#else

  SUBROUTINE InitializeMetric_MF &
    ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
#endif
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPF(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uAF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    IF( .NOT. EvolveGravity ) RETURN

#ifndef THORNADO_NOTRANSPORT

    CALL InitializeMetric_XCFC_MF_Poseidon &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uPF, MF_uAF )

#else

    CALL InitializeMetric_XCFC_MF_Poseidon &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

#endif

#endif

  END SUBROUTINE InitializeMetric_MF


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE InitializeMetricFromCheckpoint_MF( MF_uGF, MF_uCF, MF_uCR )

#else

  SUBROUTINE InitializeMetricFromCheckpoint_MF( MF_uGF, MF_uCF )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR(0:)
#endif

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    IF( .NOT. EvolveGravity ) RETURN

#ifndef THORNADO_NOTRANSPORT

    CALL InitializeMetricFromCheckpoint_XCFC_MF_Poseidon &
           ( MF_uGF, MF_uCF, MF_uCR )

#else

    CALL InitializeMetricFromCheckpoint_XCFC_MF_Poseidon &
           ( MF_uGF, MF_uCF )

#endif

#endif

  END SUBROUTINE InitializeMetricFromCheckpoint_MF


END MODULE MF_MetricInitializationModule
