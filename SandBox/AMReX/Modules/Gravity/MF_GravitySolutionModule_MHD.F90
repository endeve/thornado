MODULE MF_GravitySolutionModule_MHD

  ! --- AMReX Modules ---

  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- Local Modules ---

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  USE MF_GravitySolutionModule_XCFC_MHD, ONLY: &
    InitializeGravitySolver_XCFC_MF, &
    FinalizeGravitySolver_XCFC_MF

#elif GRAVITY_SOLVER_POSEIDON_NEWTONIAN

  USE MF_GravitySolutionModule_Newtonian_MHD, ONLY: &
    InitializeGravitySolver_Newtonian_MF, &
    FinalizeGravitySolver_Newtonian_MF

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_MF
  PUBLIC :: FinalizeGravitySolver_MF

  LOGICAL, PUBLIC :: EvolveGravity

CONTAINS


  SUBROUTINE InitializeGravitySolver_MF( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    TYPE(amrex_parmparse) :: PP

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    EvolveGravity = .FALSE.
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

    ! --- Initialized to flat space in amrex_init_from_scratch ---

    CALL InitializeGravitySolver_XCFC_MF &
           ( Verbose_Option = Verbose )

#elif GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    CALL InitializeGravitySolver_Newtonian_MF &
           ( Verbose_Option = Verbose )

#endif

  END SUBROUTINE InitializeGravitySolver_MF


  SUBROUTINE FinalizeGravitySolver_MF

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL FinalizeGravitySolver_XCFC_MF

#elif GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    CALL FinalizeGravitySolver_Newtonian_MF

#endif

  END SUBROUTINE FinalizeGravitySolver_MF


END MODULE MF_GravitySolutionModule_MHD
