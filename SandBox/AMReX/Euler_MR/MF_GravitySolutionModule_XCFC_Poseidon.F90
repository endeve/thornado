MODULE MF_GravitySolutionModule_XCFC_Poseidon

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE GravitySolutionModule_CFA_Poseidon, ONLY: &
    UpdateConformalFactorAndMetric

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    swX
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X
!!$  USE TimersModule_AMReX_Euler, ONLY: &
!!$    TimersStart_AMReX_Euler, &
!!$    TimersStop_AMReX_Euler,  &
!!$    Timer_AMReX_Euler_Allocate



#ifdef GRAVITY_SOLVER_POSEIDON_CFA

  ! --- Poseidon Modules ---

  USE Initialization_AMReX,    ONLY: &
    Initialize_Poseidon_with_AMReX
  USE Poseidon_Main_Module, ONLY: &
    Poseidon_Close, &
    Poseidon_CFA_Set_Uniform_Boundary_Conditions
  USE Source_Input_AMReX, ONLY: &
    Poseidon_XCFC_Input_Sources1_AMReX, &
    Poseidon_XCFC_Input_Sources_AMReX
  USE Poseidon_XCFC_Interface_Module, ONLY: &
    Poseidon_XCFC_Run_Part1, &
    Poseidon_XCFC_Run_Part2
  USE Return_Functions_AMReX, ONLY: &
    Poseidon_Return_ConFactor_AMReX, &
    Poseidon_Return_ALL_AMReX
  USE Initial_Guess_Module, ONLY: &
    Poseidon_Init_FlatGuess

#endif


  IMPLICIT NONE
  PRIVATE

PUBLIC :: InitializeGravitySolver_XCFC_Poseidon
PUBLIC :: FinalizeGravitySolver_XCFC_Poseidon
PUBLIC :: ComputeConformalFactor_Poseidon
PUBLIC :: ComputeGeometry_Poseidon

  PUBLIC :: UpdateConformalFactorAndMetric_MF

CONTAINS


  SUBROUTINE InitializeGravitySolver_XCFC_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CALL ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    WRITE(*,*)
    WRITE(*,'(A)') &
      '    INFO: Gravity Solver (Poseidon, XCFC)'
    WRITE(*,'(A)') &
      '    ------------------------------------'
    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Only implemented for 1D spherical symmetry.'
    WRITE(*,*)

    CALL Initialize_Poseidon_with_AMReX &
        (   FEM_Degree_Option                   = MAX( 1, nNodes - 1 ), &
            L_Limit_Option                      = 0,                    &
            Units_Option                        = 'G',                  &
            Domain_Edge_Option                  = [ xL(1), xR(1) ],     &
            Coarse_NE_Option                    = nX,                   &
            NQ_Option                           = [ nNodes, 1, 1 ],     &
            Convergence_Criteria_Option         = 1.0e-08_DP,           &
            AMReX_Max_Levels_Option             = nlevels-1,            &
            AMReX_Max_Grid_Size_Option          = MaxGridSizeX,         &   !   ***
            Verbose_Option                      = .FALSE.,              &
            Print_Setup_Option                  = .TRUE.                )


#endif

  END SUBROUTINE InitializeGravitySolver_XCFC_Poseidon


  SUBROUTINE FinalizeGravitySolver_XCFC_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    CALL Poseidon_Close()

#endif

  END SUBROUTINE FinalizeGravitySolver_XCFC_Poseidon


  SUBROUTINE ComputeConformalFactor_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, E, Si, Mg, G )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      E (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):)   ! This is psi^6 * E
    REAL(DP), INTENT(in)    :: &
      Si(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:) ! This is psi^6 * S_i
    REAL(DP), INTENT(in)    :: &
      Mg(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):)
    REAL(DP), INTENT(inout) :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP)         :: Psi_BC, AlphaPsi_BC
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)
    REAL(DP)         :: Tmp_ConFact(nDOFX,nX(1),nX(2),nX(3))

    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    ! --- Set Boundary Values ---

    CALL ComputeGravitationalMass &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mg )

    Psi_BC      = ConformalFactor( xR(1), GravitationalMass )
    AlphaPsi_BC = LapseFunction  ( xR(1), GravitationalMass ) * Psi_BC

    INNER_BC_TYPES = [ "N", "N", "N", "N", "N" ] ! Neumann
    OUTER_BC_TYPES = [ "D", "D", "D", "D", "D" ] ! Dirichlet

    INNER_BC_VALUES = [ Zero  , Zero       , Zero, Zero, Zero ]
    OUTER_BC_VALUES = [ Psi_BC, AlphaPsi_BC, Zero, Zero, Zero ]

    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions &
           ( "I", INNER_BC_TYPES, INNER_BC_VALUES )
    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions &
           ( "O", OUTER_BC_TYPES, OUTER_BC_VALUES)

    CALL Poseidon_Init_FlatGuess() ! Possibly move this to init call

    ! --- Set matter sources with current conformal factor ---
    CALL Poseidon_XCFC_Input_Sources1_AMReX &
           ( MF_Src_Input = MF_Source,            & ! ***
             Input_NQ     = nNodesX,              &
             Input_R_Quad = MeshX(1) % Nodes,     &
             Input_T_Quad = MeshX(2) % Nodes,     &
             Input_P_Quad = MeshX(3) % Nodes,     &
             Input_xL     = [ -Half, +Half ]      )

    ! --- Compute conformal factor ---

    CALL Poseidon_XCFC_Run_Part1()

    CALL Poseidon_Return_ConFactor_AMReX &
         ( NQ               = nNodesX,          &
           RQ_Input         = MeshX(1) % Nodes, &
           TQ_Input         = MeshX(2) % Nodes, &
           PQ_Input         = MeshX(3) % Nodes, &
           Left_Limit       = -Half,            &
           Right_Limit      = +Half,            &
           nLevels          = nLevels,          &   ! ***
           MF_Results       = MF_Results )          ! ***


    CALL UpdateConformalFactorAndMetric &
           ( iX_B0, iX_E0, iX_B1, iX_E1, Tmp_ConFact, G )

#endif

    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE ComputeConformalFactor_Poseidon


  SUBROUTINE ComputeGeometry_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, E, S, Si, G )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      E (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):), &
      S (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):), &
      Si(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(inout)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP) :: Tmp_Lapse             (nDOFX,nX(1),nX(2),nX(3)), &
                Tmp_Shift             (nDOFX,nX(1),nX(2),nX(3),1:3), &
                Tmp_ExtrinsicCurvature(nDOFX,nX(1),nX(2),nX(3),1:6)

    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    ! --- Set matter sources with updated conformal factor ---

    CALL Poseidon_XCFC_Input_Sources_AMReX &
           ( MF_Src_Input = MF_Source,            & ! ***
             Input_NQ     = nNodesX,              &
             Input_R_Quad = MeshX(1) % Nodes,     &
             Input_T_Quad = MeshX(2) % Nodes,     &
             Input_P_Quad = MeshX(3) % Nodes,     &
             Input_xL     = [ -Half, +Half ]      )


    ! --- Compute lapse and shift ---

    CALL Poseidon_XCFC_Run_Part2()



    CALL Poseidon_Return_ALL_AMReX &
         ( NQ               = nNodesX,          &
           RQ_Input         = MeshX(1) % Nodes, &
           TQ_Input         = MeshX(2) % Nodes, &
           PQ_Input         = MeshX(3) % Nodes, &
           Left_Limit       = -Half,            &
           Right_Limit      = +Half,            &
           nLevels          = nLevels,          &
           MF_Results       = MF_Results        )


    ! --- Copy data from Poseidon arrays to thornado arrays ---

    CALL ComputeGeometryFromPoseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             Tmp_Lapse, Tmp_Shift, Tmp_ExtrinsicCurvature, G )

    CALL SetBoundaryConditions &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G )

#endif

    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE ComputeGeometry_Poseidon




  SUBROUTINE UpdateConformalFactorAndMetric_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    ! 1: psi, 2: alpha, 3-5: beta, 6-11: K_ij
    INTEGER, PARAMETER :: nMF = 11

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iLo_G(4), iLo_M(4)

    REAL(DP), CONTIGUOUS, POINTER :: uMF (:,:,:,:)
    REAL(DP), ALLOCATABLE         :: M   (:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G   (:,:,:,:,:)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        iLo_M = LBOUND( uMF )
        iLo_G = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

!!$        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        ALLOCATE( M(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nMF) )

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

!!$        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_X &
               ( nMF, iX_B0, iX_E0, iLo_M, iX_B0, iX_E0, uMF, M )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_G, iX_B1, iX_E1, uGF, G )

        CALL UpdateConformalFactorAndMetric &
               ( iX_B0, iX_E0, iX_B1, iX_E1, M(:,:,:,:,1), G )

        CALL thornado2amrex_X &
               ( nGF, iX_B1, iX_E1, iLo_G, iX_B1, iX_E1, uGF, G )

!!$        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( G )
        DEALLOCATE( M )

!!$        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE UpdateConformalFactorAndMetric_MF


END MODULE MF_GravitySolutionModule_XCFC_Poseidon
