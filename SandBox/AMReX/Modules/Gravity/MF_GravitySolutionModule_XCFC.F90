MODULE MF_GravitySolutionModule_XCFC

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- Local Modules ---

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  USE MF_GravitySolutionModule_XCFC_Poseidon, ONLY: &
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

  LOGICAL, PUBLIC :: EvolveGravity

CONTAINS


  SUBROUTINE InitializeGravitySolver_XCFC_MF( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    TYPE(amrex_parmparse) :: PP

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    ! --- Initialized to flat space in amrex_init_from_scratch ---

    EvolveGravity = .TRUE.
    CALL amrex_parmparse_build( PP, 'GS' )
      CALL PP % query( 'EvolveGravity', EvolveGravity )
    CALL amrex_parmparse_destroy( PP )

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A)')   'INFO: Gravity Solver'
      WRITE(*,'(4x,A)')   '--------------------'
      WRITE(*,*)
      WRITE(*,'(6x,A18,L)') 'EvolveGravity: ', EvolveGravity

    END IF

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL InitializeGravitySolver_XCFC_MF_Poseidon

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


END MODULE MF_GravitySolutionModule_XCFC
