MODULE MF_GravitySolutionModule_Newtonian

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- Local Modules ---

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

  USE MF_GravitySolutionModule_Newtonian_Poseidon_MHD, ONLY: &
    InitializeGravitySolver_Newtonian_MF_Poseidon, &
    FinalizeGravitySolver_Newtonian_MF_Poseidon, &
    ComputeGravitationalPotential_Newtonian_MF_Poseidon

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_Newtonian_MF
  PUBLIC :: FinalizeGravitySolver_Newtonian_MF
  PUBLIC :: ComputeGravitationalPotential_Newtonian_MF

CONTAINS


  SUBROUTINE InitializeGravitySolver_Newtonian_MF( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL InitializeGravitySolver_Newtonian_MF_Poseidon &
           ( Verbose_Option = Verbose )

#endif

  END SUBROUTINE InitializeGravitySolver_Newtonian_MF


  SUBROUTINE FinalizeGravitySolver_Newtonian_MF

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    CALL FinalizeGravitySolver_Newtonian_MF_Poseidon

#endif

  END SUBROUTINE FinalizeGravitySolver_Newtonian_MF


  SUBROUTINE ComputeGravitationalPotential_Newtonian_MF( MF_uCF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    CALL ComputeGravitationalPotential_Newtonian_MF_Poseidon( MF_uCF, MF_uGF )

#endif

  END SUBROUTINE ComputeGravitationalPotential_Newtonian_MF


END MODULE MF_GravitySolutionModule_Newtonian
