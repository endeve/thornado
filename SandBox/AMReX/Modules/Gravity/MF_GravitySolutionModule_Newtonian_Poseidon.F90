MODULE MF_GravitySolutionModule_Newtonian_Poseidon

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nDOFX
  USE MeshModule, ONLY: &
    MeshX
  USE UnitsModule, ONLY: &
    Kilometer, &
    Kilogram, &
    Gram, &
    Meter, &
    GravitationalConstant, &
    Centimeter
  USE FluidFieldsModule, ONLY: &
    iCF_D
  USE GeometryFieldsModule, ONLY: &
    iGF_Phi_N

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    Half, &
    Zero, &
    DP, &
    Pi

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

  ! --- Poseidon Modules ---

  USE Poseidon_Interface_Initialization, ONLY: &
    Initialize_Poseidon
  USE Poseidon_Interface_Boundary_Conditions, ONLY : &
    Poseidon_Set_Uniform_Boundary_Conditions
  USE Poseidon_Interface_Source_Input, ONLY: &
    Poseidon_Input_Sources
  USE Poseidon_Interface_Run, ONLY: &
    Poseidon_Run
  USE Poseidon_Interface_Return_Routines, ONLY: &
    Poseidon_Return_Newtonian_Potential
  USE Poseidon_Interface_Close, ONLY: &
    Poseidon_Close

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_Newtonian_MF_Poseidon
  PUBLIC :: FinalizeGravitySolver_Newtonian_MF_Poseidon
  PUBLIC :: ComputeGravitationalPotential_Newtonian_MF_Poseidon

CONTAINS


! ----------- NEW CHANGES -----!!!!

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
	     Newtonian_Mode_Option        = .TRUE. ,          &
             Verbose_Option               = .TRUE.,           &
             Print_Setup_Option           = .TRUE.,           &
             Print_Results_Option         = .FALSE. )

#endif

  END SUBROUTINE InitializeGravitySolver_Newtonian_MF_Poseidon


  SUBROUTINE ComputeGravitationalPotential_Newtonian_MF_Poseidon( MF_uCF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    TYPE(amrex_parmparse) :: PP
    REAL(DP) :: Phi_r, rho_c, r_c, SolarRadius, R_star, Domain
    REAL(DP) :: Part1, Part2
    INTEGER :: iLevel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    !TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab)                 :: MF_uGS(0:SIZE(MF_uCF)-1)
    TYPE(amrex_multifab)                 :: MF_uMF(0:SIZE(MF_uGF)-1)

    ! Build the multifab



    DO iLevel = 0, SIZE(MF_uCF)-1

      CALL amrex_multifab_build &
              ( MF_uMF(iLevel), MF_uCF(iLevel) % BA, &
                                MF_uCF(iLevel) % DM, nDOFX, 0 )

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), MF_uCF(iLevel) % BA, &
                               MF_uCF(iLevel) % DM, nDOFX, 0 )

      CALL MF_uGS(iLevel) % COPY &
             ( MF_uCF(iLevel), nDOFX*(iCF_D-1)+1, 1, nDOFX, 0 )

    END DO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    CALL Poseidon_Input_Sources( MF_uGS )

    !-- Inner Boundary --!

    CALL Poseidon_Set_Uniform_Boundary_Conditions &
  	   ( BC_Location_Input = 'I', &
    	     BC_Type_Input     = 'N', &
    	     BC_Value_Input    = Zero )

    !-- Outer Boundary --!

    !-- Computing BC at R= 2*R --!
    ! Constants !


    CALL amrex_parmparse_build( PP, 'inputs' )
      CALL PP % get( 'rho_c'      , rho_c )
      CALL PP % get( 'R_star'     , R_star )
      CALL PP % get( 'SolarRadius', SolarRadius )
      CALL PP % get( 'r_c'        ,r_c )
    CALL amrex_parmparse_destroy( PP )

    rho_c       = rho_c * ( Gram / Centimeter**3 )
    SolarRadius = SolarRadius * Kilometer
    R_star      = R_star * SolarRadius
    r_c         = r_c * SolarRadius
    Domain           = 2.0_DP * R_star

    ! -- Evaluating Phi_r at 2* R_star --!
    Part1       = -4.0_DP * Pi * ( r_c**3 / Domain )
    Part2       = (R_star / r_c ) - atan( R_star / r_c )
    Phi_r       = Part1 * Part2

    CALL Poseidon_Set_Uniform_Boundary_Conditions &
           ( BC_Location_Input = 'O', &
    	     BC_Type_Input     = 'D', &
    	     BC_Value_Input    = Phi_r )

   CALL Poseidon_Run ()

   CALL Poseidon_Return_Newtonian_Potential( MF_uMF )


   DO iLevel = 0, SIZE(MF_uGF)-1

     CALL MF_uGF(iLevel) % COPY &
            ( MF_uMF(iLevel), 1, nDOFX*(iGF_Phi_N-1)+1, nDOFX, 0 )

     CALL amrex_multifab_destroy( MF_uGS(iLevel) )

     CALL amrex_multifab_destroy( MF_uMF(iLevel) )

   END DO


#endif

  END SUBROUTINE ComputeGravitationalPotential_Newtonian_MF_Poseidon

!! ----------- END ------!!!!

  SUBROUTINE FinalizeGravitySolver_Newtonian_MF_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTONIAN

    CALL Poseidon_Close()

#endif

  END SUBROUTINE FinalizeGravitySolver_Newtonian_MF_Poseidon


END MODULE MF_GravitySolutionModule_Newtonian_Poseidon
