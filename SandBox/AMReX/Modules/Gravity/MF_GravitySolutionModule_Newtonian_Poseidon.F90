MODULE MF_GravitySolutionModule_Newtonian_Poseidon

  ! --- AMReX Modules ---

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodesX
  USE MeshModule, ONLY: &
    MeshX

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    Half

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

  ! --- Poseidon Modules ---

  USE Poseidon_Interface_Initialization, ONLY: &
    Initialize_Poseidon
  USE Poseidon_Interface_Close, ONLY: &
    Poseidon_Close

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_Newtonian_MF_Poseidon
  PUBLIC :: FinalizeGravitySolver_Newtonian_MF_Poseidon

CONTAINS


  SUBROUTINE InitializeGravitySolver_Newtonian_MF_Poseidon( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A)') &
        'INFO: Gravity Solver (Poseidon, Newtonian)'
      WRITE(*,'(4x,A)') &
        '------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'Only implemented for 1D spherical symmetry.'
      WRITE(*,*)

    END IF

    CALL Initialize_Poseidon &
           ( Source_NQ                    = nNodesX,          &
             Source_xL                    = [ -Half, +Half ], &
             Source_RQ_xlocs              = MeshX(1) % Nodes, &
             Source_TQ_xlocs              = MeshX(2) % Nodes, &
             Source_PQ_xlocs              = MeshX(3) % Nodes, &
             Source_Units                 = 'G',              &
             Source_Radial_Boundary_Units = 'km',             &
             Flat_Guess_Option            = .TRUE.,           &
             Verbose_Option               = .FALSE.,          &
             Print_Setup_Option           = .TRUE.,           &
             Print_Results_Option         = .FALSE. )

#endif

  END SUBROUTINE InitializeGravitySolver_Newtonian_MF_Poseidon


  SUBROUTINE FinalizeGravitySolver_Newtonian_MF_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    CALL Poseidon_Close()

#endif

  END SUBROUTINE FinalizeGravitySolver_Newtonian_MF_Poseidon


END MODULE MF_GravitySolutionModule_Newtonian_Poseidon
