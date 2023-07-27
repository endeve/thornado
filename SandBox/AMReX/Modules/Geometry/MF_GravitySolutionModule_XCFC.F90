MODULE MF_GravitySolutionModule_XCFC

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- Local Modules ---

  USE InputParsingModule, ONLY: &
    nLevels
  USE MF_GravitySolutionModule_XCFC_Poseidon, ONLY: &
    InitializeGravitySolver_XCFC_MF_Poseidon, &
    FinalizeGravitySolver_XCFC_MF_Poseidon, &
    ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF_Poseidon, &
    ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF_Poseidon, &
    ComputeConformalFactor_MF_Poseidon, &
    ComputePressureTensorTrace_XCFC_Euler_MF_Poseidon, &
    ComputePressureTensorTrace_XCFC_TwoMoment_MF_Poseidon, &
    ComputeGeometry_MF_Poseidon, &
    InitializeMetric_Euler_MF_Poseidon, &
    InitializeMetric_TwoMoment_MF_Poseidon

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_XCFC_MF
  PUBLIC :: FinalizeGravitySolver_XCFC_MF
  PUBLIC :: ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF
  PUBLIC :: ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF
  PUBLIC :: ComputeConformalFactor_MF
  PUBLIC :: ComputePressureTensorTrace_XCFC_Euler_MF
  PUBLIC :: ComputePressureTensorTrace_XCFC_TwoMoment_MF
  PUBLIC :: ComputeGeometry_MF
  PUBLIC :: InitializeMetric_Euler_MF
  PUBLIC :: InitializeMetric_TwoMoment_MF

CONTAINS


  SUBROUTINE InitializeGravitySolver_XCFC_MF &
    ( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(inout), OPTIONAL :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)   , OPTIONAL :: MF_uCF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL InitializeGravitySolver_XCFC_MF_Poseidon &
           ( MF_uGF, MF_uCF )

#endif

  END SUBROUTINE InitializeGravitySolver_XCFC_MF


  SUBROUTINE FinalizeGravitySolver_XCFC_MF

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL FinalizeGravitySolver_XCFC_MF_Poseidon

#endif

  END SUBROUTINE FinalizeGravitySolver_XCFC_MF


  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF &
    ( MF_uGF, MF_uCF, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF_Poseidon &
           ( MF_uGF, MF_uCF, MF_uGS )

#endif

  END SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF


  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF_Poseidon &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

#endif

  END SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF


  SUBROUTINE ComputeConformalFactor_MF( MF_uGS, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:nLevels-1) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputeConformalFactor_MF_Poseidon( MF_uGS, MF_uGF )

#endif

  END SUBROUTINE ComputeConformalFactor_MF


  SUBROUTINE ComputePressureTensorTrace_XCFC_Euler_MF &
    ( MF_uGF, MF_uCF, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputePressureTensorTrace_XCFC_Euler_MF_Poseidon &
           ( MF_uGF, MF_uCF, MF_uGS )

#endif

  END SUBROUTINE ComputePressureTensorTrace_XCFC_Euler_MF


  SUBROUTINE ComputePressureTensorTrace_XCFC_TwoMoment_MF &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputePressureTensorTrace_XCFC_TwoMoment_MF_Poseidon &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

#endif

  END SUBROUTINE ComputePressureTensorTrace_XCFC_TwoMoment_MF


  SUBROUTINE ComputeGeometry_MF( MF_uGS, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:nLevels-1) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputeGeometry_MF_Poseidon( MF_uGS, MF_uGF )

#endif

  END SUBROUTINE ComputeGeometry_MF


  SUBROUTINE InitializeMetric_Euler_MF( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uAF(0:nLevels-1)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL InitializeMetric_Euler_MF_Poseidon &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

#endif

  END SUBROUTINE InitializeMetric_Euler_MF


  SUBROUTINE InitializeMetric_TwoMoment_MF &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uAF(0:nLevels-1)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL InitializeMetric_TwoMoment_MF_Poseidon &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uPF, MF_uAF )

#endif

  END SUBROUTINE InitializeMetric_TwoMoment_MF


END MODULE MF_GravitySolutionModule_XCFC
