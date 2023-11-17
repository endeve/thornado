MODULE GravitySolutionModule_XCFC_Poseidon

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half
  USE ProgramHeaderModule, ONLY: &
    nX, &
    nNodesX, &
    nNodes, &
    xL, &
    xR
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryComputationModule, ONLY: &
    LapseFunction, &
    ConformalFactor
  USE XCFC_UtilitiesModule, ONLY: &
    iGS_E, &
    iGS_S1, &
    iGS_S2, &
    iGS_S3, &
    iGS_S, &
    iMF_Psi, &
    iMF_Alpha, &
    iMF_Beta_1, &
    iMF_Beta_2, &
    iMF_Beta_3, &
    iMF_K_dd_11, &
    iMF_K_dd_12, &
    iMF_K_dd_13, &
    iMF_K_dd_22, &
    iMF_K_dd_23, &
    iMF_K_dd_33, &
    ComputeGravitationalMass

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  ! --- Poseidon Modules ---

  USE Poseidon_Interface_Initialization, ONLY: &
    Initialize_Poseidon
  USE Poseidon_Interface_Boundary_Conditions, ONLY : &
    Poseidon_Set_Uniform_Boundary_Conditions
  USE Poseidon_Interface_Source_Input, ONLY: &
    Poseidon_Input_Sources_Part1, &
    Poseidon_Input_Sources_Part2
  USE Poseidon_Interface_Run, ONLY: &
    Poseidon_XCFC_Run_Part1, &
    Poseidon_XCFC_Run_Part2
  USE Poseidon_Interface_Return_Routines, ONLY: &
    Poseidon_Return_Conformal_Factor, &
    Poseidon_Return_Lapse_Function, &
    Poseidon_Return_Shift_Vector, &
    Poseidon_Return_Extrinsic_Curvature
  USE Poseidon_Interface_Close, ONLY: &
    Poseidon_Close

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_XCFC_Poseidon
  PUBLIC :: FinalizeGravitySolver_XCFC_Poseidon
  PUBLIC :: ComputeConformalFactor_XCFC_Poseidon
  PUBLIC :: ComputeLapseShiftCurvature_XCFC_Poseidon

CONTAINS


  SUBROUTINE InitializeGravitySolver_XCFC_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    WRITE(*,*)
    WRITE(*,'(A)') &
      '    INFO: Gravity Solver (Poseidon, XCFC)'
    WRITE(*,'(A)') &
      '    -------------------------------------'
    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Only implemented for 1D spherical symmetry.'
    WRITE(*,*)

    CALL Initialize_Poseidon &
         ( Dimensions_Option            = 3,                         &
           FEM_Degree_Option            = MAX( 1, nNodes - 1 ),      &
           L_Limit_Option               = 0,                         &
           Source_NE                    = nX,                        &
           Domain_Edge_Option           = [ xL(1), xR(1) ],          &
           Source_NQ                    = nNodesX,                   &
           Source_xL                    = [ -Half, +Half ],          &
           Source_RQ_xlocs              = MeshX(1) % Nodes,          &
           Source_TQ_xlocs              = MeshX(2) % Nodes,          &
           Source_PQ_xlocs              = MeshX(3) % Nodes,          &
           Source_Units                 = 'G',                       &
           Source_Radial_Boundary_Units = ' m',                      &
           Source_DR_Option             = MeshX(1) % Width(1:nX(1)), &
           Source_DT_Option             = MeshX(2) % Width(1:nX(2)), &
           Source_DP_Option             = MeshX(3) % Width(1:nX(3)), &
           Max_Iterations_Option        = 20,                        &
           Flat_Guess_Option            = .TRUE.,                    &
           Print_Setup_Option           = .TRUE.,                    &
           Convergence_Criteria_Option  = 1.0e-08_DP,                &
           Verbose_Option               = .FALSE.,                   &
           Print_Results_Option         = .FALSE. )

#endif

  END SUBROUTINE InitializeGravitySolver_XCFC_Poseidon


  SUBROUTINE FinalizeGravitySolver_XCFC_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL Poseidon_Close()

#endif

  END SUBROUTINE FinalizeGravitySolver_XCFC_Poseidon


  SUBROUTINE ComputeConformalFactor_XCFC_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GS, M )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      GS(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:) ! This contains psi^6 * { E, S_i }
    REAL(DP), INTENT(inout) :: &
      M (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    REAL(DP)         :: GravitationalMass, Psi_xR, AlphaPsi_xR, Beta_u_xR(3)
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)

    ! --- Set Boundary Values ---

    GravitationalMass = Zero

    CALL ComputeGravitationalMass &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GS, GravitationalMass )

    ! --- Approximate outer boundary with isotropic expressions ---

    Psi_xR       = ConformalFactor( xR(1), GravitationalMass )
    AlphaPsi_xR  = LapseFunction  ( xR(1), GravitationalMass ) * Psi_xR
    Beta_u_xR(1) = Zero
    Beta_u_xR(2) = Zero
    Beta_u_xR(3) = Zero

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    INNER_BC_TYPES = [ 'N', 'N', 'N', 'N', 'N' ] ! Neumann
    OUTER_BC_TYPES = [ 'D', 'D', 'D', 'D', 'D' ] ! Dirichlet

    INNER_BC_VALUES &
      = [ Zero  , Zero       , Zero        , Zero        , Zero ]
    OUTER_BC_VALUES &
      = [ Psi_xR, AlphaPsi_xR, Beta_u_xR(1), Beta_u_xR(2), Beta_u_xR(3) ]

    CALL Poseidon_Set_Uniform_Boundary_Conditions &
           ( 'I', INNER_BC_TYPES, INNER_BC_VALUES )
    CALL Poseidon_Set_Uniform_Boundary_Conditions &
           ( 'O', OUTER_BC_TYPES, OUTER_BC_VALUES )

    ! --- Set matter sources with current conformal factor ---

    CALL Poseidon_Input_Sources_Part1 &
           ( Input_E  = GS(:,:,:,:,iGS_E), &
             Input_Si = GS(:,:,:,:,iGS_S1:iGS_S3) )

    ! --- Compute conformal factor ---

    CALL Poseidon_XCFC_Run_Part1()

    CALL Poseidon_Return_Conformal_Factor &
         ( Return_ConFactor = M(:,:,:,:,iMF_Psi) )

#endif

  END SUBROUTINE ComputeConformalFactor_XCFC_Poseidon


  SUBROUTINE ComputeLapseShiftCurvature_XCFC_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GS, M )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      GS(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(inout)    :: &
      M (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ! --- Set matter sources with updated conformal factor ---

    CALL Poseidon_Input_Sources_Part2 &
           ( Input_E  = GS(:,:,:,:,iGS_E),  &
             Input_Si = GS(:,:,:,:,iGS_S1:iGS_S3), &
             Input_S  = GS(:,:,:,:,iGS_S)  )

    ! --- Compute lapse and shift ---

    CALL Poseidon_XCFC_Run_Part2()

    CALL Poseidon_Return_Lapse_Function &
           ( Return_Lapse = M(:,:,:,:,iMF_Alpha) )

    CALL Poseidon_Return_Shift_Vector &
           ( Return_Shift = M(:,:,:,:,iMF_Beta_1:iMF_Beta_3) )

    CALL Poseidon_Return_Extrinsic_Curvature &
           ( Return_Kij = M(:,:,:,:,iMF_K_dd_11:iMF_K_dd_33) )

#endif

  END SUBROUTINE ComputeLapseShiftCurvature_XCFC_Poseidon


END MODULE GravitySolutionModule_XCFC_Poseidon
