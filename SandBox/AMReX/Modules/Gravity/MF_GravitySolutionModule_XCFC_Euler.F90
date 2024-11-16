MODULE MF_GravitySolutionModule_XCFC_Euler

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- Local Modules ---

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  USE MF_GravitySolutionModule_XCFC_Poseidon_Euler, ONLY: &
    InitializeGravitySolver_XCFC_MF_Poseidon, &
    FinalizeGravitySolver_XCFC_MF_Poseidon, &
    ComputeConformalFactor_XCFC_MF_Poseidon, &
    ComputeLapseShiftCurvature_XCFC_MF_Poseidon

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_XCFC_MF
  PUBLIC :: FinalizeGravitySolver_XCFC_MF
  PUBLIC :: ComputeConformalFactor_XCFC_MF
  PUBLIC :: ComputeLapseShiftCurvature_XCFC_MF

CONTAINS


  SUBROUTINE InitializeGravitySolver_XCFC_MF( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL InitializeGravitySolver_XCFC_MF_Poseidon &
           ( Verbose_Option = Verbose )

#endif

  END SUBROUTINE InitializeGravitySolver_XCFC_MF


  SUBROUTINE FinalizeGravitySolver_XCFC_MF

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL FinalizeGravitySolver_XCFC_MF_Poseidon

#endif

  END SUBROUTINE FinalizeGravitySolver_XCFC_MF


  SUBROUTINE ComputeConformalFactor_XCFC_MF( MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputeConformalFactor_XCFC_MF_Poseidon( MF_uGS, MF_uMF )

#endif

  END SUBROUTINE ComputeConformalFactor_XCFC_MF


  SUBROUTINE ComputeLapseShiftCurvature_XCFC_MF( MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputeLapseShiftCurvature_XCFC_MF_Poseidon( MF_uGS, MF_uMF )

#endif

  END SUBROUTINE ComputeLapseShiftCurvature_XCFC_MF


END MODULE MF_GravitySolutionModule_XCFC_Euler
