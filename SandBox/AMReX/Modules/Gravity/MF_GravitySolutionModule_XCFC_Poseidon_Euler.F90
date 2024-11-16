MODULE MF_GravitySolutionModule_XCFC_Poseidon_Euler

  ! --- AMReX Modules ---

  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_multifab_module, ONLY: &
    amrex_multifab

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    xR, &
    swX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryComputationModule, ONLY: &
    LapseFunction, &
    ConformalFactor
  USE XCFC_UtilitiesModule, ONLY: &
    swX_GS

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Half
  USE AverageDownModule_Euler, ONLY: &
    AverageDown
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF
  USE MF_XCFC_UtilitiesModule, ONLY: &
    UpdateConformalFactorAndMetric_XCFC_MF, &
    ComputeGravitationalMass_MF, &
    ComputeConformalFactorSourcesAndMg_XCFC_MF, &
    ComputePressureTensorTrace_XCFC_MF

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
    Poseidon_Return_ALL
  USE Poseidon_Interface_Close, ONLY: &
    Poseidon_Close

#endif

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: FillGhostCells

  PUBLIC :: InitializeGravitySolver_XCFC_MF_Poseidon
  PUBLIC :: FinalizeGravitySolver_XCFC_MF_Poseidon
  PUBLIC :: ComputeConformalFactor_XCFC_MF_Poseidon
  PUBLIC :: ComputeLapseShiftCurvature_XCFC_MF_Poseidon

CONTAINS


  SUBROUTINE InitializeGravitySolver_XCFC_MF_Poseidon( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    TYPE(amrex_parmparse) :: PP

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    FillGhostCells = .FALSE.
    CALL amrex_parmparse_build( PP, 'poseidon' )
      CALL PP % query( 'FillGhostCells', FillGhostCells )
    CALL amrex_parmparse_destroy( PP )

    swX_GS = 0
    IF( FillGhostCells ) swX_GS = swX

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A)') &
        'INFO: Gravity Solver (Poseidon, XCFC)'
      WRITE(*,'(4x,A)') &
        '-------------------------------------'
      WRITE(*,*)
      WRITE(*,'(6x,A,L)') 'FillGhostCells: ', FillGhostCells
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

  END SUBROUTINE InitializeGravitySolver_XCFC_MF_Poseidon


  SUBROUTINE FinalizeGravitySolver_XCFC_MF_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL Poseidon_Close()

#endif

  END SUBROUTINE FinalizeGravitySolver_XCFC_MF_Poseidon


  SUBROUTINE ComputeConformalFactor_XCFC_MF_Poseidon( MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    REAL(DP)         :: GravitationalMass, Psi_xR, AlphaPsi_xR, Beta_u_xR(3)
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)

    ! --- Set Boundary Values ---

    CALL ComputeGravitationalMass_MF( MF_uGS, GravitationalMass )

    ! --- Approximate outer boundary with isotropic expressions ---

    Psi_xR       = ConformalFactor( xR(1), GravitationalMass )
    AlphaPsi_xR  = LapseFunction  ( xR(1), GravitationalMass ) * Psi_xR
    Beta_u_xR(1) = Zero
    Beta_u_xR(2) = Zero
    Beta_u_xR(3) = Zero

    INNER_BC_TYPES = [ 'N', 'N', 'N', 'N', 'N' ] ! Neumann
    OUTER_BC_TYPES = [ 'D', 'D', 'D', 'D', 'D' ] ! Dirichlet

    INNER_BC_VALUES &
      = [ Zero  , Zero       , Zero        , Zero        , Zero ]
    OUTER_BC_VALUES &
      = [ Psi_xR, AlphaPsi_xR, Beta_u_xR(1), Beta_u_xR(2), Beta_u_xR(3) ]

    CALL Poseidon_Set_Uniform_Boundary_Conditions &
           ( 'I', INNER_BC_TYPES, INNER_BC_VALUES )
    CALL Poseidon_Set_Uniform_Boundary_Conditions &
           ( 'O', OUTER_BC_TYPES, OUTER_BC_VALUES)

    ! --- Set XCFC sources with current conformal factor ---
    CALL Poseidon_Input_Sources_Part1( MF_uGS )

    ! --- Compute conformal factor ---

    CALL Poseidon_XCFC_Run_Part1()

    CALL Poseidon_Return_Conformal_Factor &
           ( MF_uMF, FillGhostCells_Option = FillGhostCells )

#endif

  END SUBROUTINE ComputeConformalFactor_XCFC_MF_Poseidon


  SUBROUTINE ComputeLapseShiftCurvature_XCFC_MF_Poseidon( MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ! --- Set gravity sources with updated conformal factor ---

    CALL Poseidon_Input_Sources_Part2( MF_uGS )

    ! --- Compute lapse and shift ---

    CALL Poseidon_XCFC_Run_Part2()

    CALL Poseidon_Return_ALL &
           ( MF_uMF, FillGhostCells_Option = FillGhostCells )

#endif

  END SUBROUTINE ComputeLapseShiftCurvature_XCFC_MF_Poseidon


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE ComputeConformalFactor &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS, MF_uMF )

#else

  SUBROUTINE ComputeConformalFactor &
    ( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
#endif
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

#ifndef THORNADO_NOTRANSPORT

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

#else

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uGS )

#endif

    CALL ComputeConformalFactor_XCFC_MF_Poseidon( MF_uGS, MF_uMF )

    CALL UpdateConformalFactorAndMetric_XCFC_MF( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF, UpdateSpatialMetric_Option = .TRUE. )

#endif

  END SUBROUTINE ComputeConformalFactor


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE ComputeLapseShiftCurvature &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS, MF_uMF )

#else

  SUBROUTINE ComputeLapseShiftCurvature &
    ( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
#endif
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

#ifndef THORNADO_NOTRANSPORT

    CALL ComputePressureTensorTrace_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

#else

    CALL ComputePressureTensorTrace_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uGS )

#endif

    CALL ComputeLapseShiftCurvature_XCFC_MF_Poseidon &
           ( MF_uGS, MF_uMF )

    CALL AverageDown( MF_uGF, UpdateSpatialMetric_Option = .TRUE. )

    CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )

#endif

  END SUBROUTINE ComputeLapseShiftCurvature


END MODULE MF_GravitySolutionModule_XCFC_Poseidon_Euler
