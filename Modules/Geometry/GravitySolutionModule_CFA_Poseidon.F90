MODULE GravitySolutionModule_CFA_Poseidon

  USE KindModule, ONLY: &
    DP, &
    Pi, &
    FourPi, &
    Zero, &
    Half, &
    One, &
    Two
  USE ProgramHeaderModule, ONLY: &
    nX, &
    nNodesX, &
    nDOFX, &
    nNodes, &
    xL, &
    xR
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    NodeNumberTableX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    iGF_Phi_N, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Psi, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_K_dd_11, &
    iGF_K_dd_12, &
    iGF_K_dd_13, &
    iGF_K_dd_22, &
    iGF_K_dd_23, &
    iGF_K_dd_33
  USE GeometryComputationModule, ONLY: &
    LapseFunction, &
    ConformalFactor, &
    ComputeGeometryX_FromScaleFactors, &
    ComputeGeometryX
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler,  &
    Timer_GravitySolver

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

  ! --- Poseidon Modules ---

  USE Initialization_Poseidon, ONLY: &
    Initialize_Poseidon
  USE Poseidon_Main_Module, ONLY: &
    Poseidon_Close, &
    Poseidon_CFA_Set_Uniform_Boundary_Conditions
  USE Source_Input_Module, ONLY: &
    Poseidon_XCFC_Input_Sources1, &
    Poseidon_XCFC_Input_Sources2
  USE Poseidon_XCFC_Interface_Module, ONLY: &
    Poseidon_XCFC_Run_Part1, &
    Poseidon_XCFC_Run_Part2, &
    Poseidon_Return_ConFactor, &
    Poseidon_Return_Lapse, &
    Poseidon_Return_Shift, &
    Poseidon_Return_ExtrinsicCurvature
  USE Initial_Guess_Module, ONLY: &
    Poseidon_Init_FlatGuess

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_CFA_Poseidon
  PUBLIC :: FinalizeGravitySolver_CFA_Poseidon
  PUBLIC :: ComputeConformalFactor_Poseidon
  PUBLIC :: ComputeGeometry_Poseidon

  REAL(DP) :: GravitationalMass


CONTAINS


  SUBROUTINE InitializeGravitySolver_CFA_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CALL ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    WRITE(*,*)
    WRITE(*,'(A)') &
      '    INFO: Gravity Solver (Poseidon, CFA)'
    WRITE(*,'(A)') &
      '    ------------------------------------'
    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Only implemented for 1D spherical symmetry.'
    WRITE(*,*)

    CALL Initialize_Poseidon &
         ( Units_Option                = 'G',                       &
           Dimensions_Option           = 3,                         &
           FEM_Degree_Option           = MAX( 1, nNodes - 1 ),      &
           L_Limit_Option              = 0,                         &
           Domain_Edge_Option          = [ xL(1), xR(1) ],          &
           NE_Option                   = nX,                        &
           NQ_Option                   = [ nNodes, 1, 1 ],          &
           dr_Option                   = MeshX(1) % Width(1:nX(1)), &
           dt_Option                   = MeshX(2) % Width(1:nX(2)), &
           dp_Option                   = MeshX(3) % Width(1:nX(3)), &
           Method_Flag_Option          = 3,                         &
           Print_Setup_Option          = .TRUE.,                    &
           Convergence_Criteria_Option = 1.0e-08_DP,                &
           Verbose_Option              = .FALSE. )

#endif

  END SUBROUTINE InitializeGravitySolver_CFA_Poseidon


  SUBROUTINE FinalizeGravitySolver_CFA_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    CALL Poseidon_Close()

#endif

  END SUBROUTINE FinalizeGravitySolver_CFA_Poseidon


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

    CALL Poseidon_XCFC_Input_Sources1 &
           ( Local_E      = E,                    &
             Local_Si     = Si,                   &
             Local_RE_Dim = nX(1),                &
             Local_TE_Dim = nX(2),                &
             Local_PE_Dim = nX(3),                &
             Local_RQ_Dim = nNodesX(1),           &
             Local_TQ_Dim = nNodesX(2),           &
             Local_PQ_Dim = nNodesX(3),           &
             Input_R_Quad = MeshX(1) % Nodes,     &
             Input_T_Quad = MeshX(2) % Nodes,     &
             Input_P_Quad = MeshX(3) % Nodes,     &
             Left_Limit   = -Half,                &
             Right_Limit  = +Half )

    ! --- Compute conformal factor ---

    CALL Poseidon_XCFC_Run_Part1()

    CALL Poseidon_Return_ConFactor &
         ( NE               = nX,               &
           NQ               = nNodesX,          &
           RQ_Input         = MeshX(1) % Nodes, &
           TQ_Input         = MeshX(2) % Nodes, &
           PQ_Input         = MeshX(3) % Nodes, &
           Left_Limit       = -Half,            &
           Right_Limit      = +Half,            &
           Return_ConFactor = Tmp_ConFact )

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

    CALL Poseidon_XCFC_Input_Sources1 &
           ( Local_E      = E,                    &
             Local_Si     = Si,                   &
             Local_RE_Dim = nX(1),                &
             Local_TE_Dim = nX(2),                &
             Local_PE_Dim = nX(3),                &
             Local_RQ_Dim = nNodesX(1),           &
             Local_TQ_Dim = nNodesX(2),           &
             Local_PQ_Dim = nNodesX(3),           &
             Input_R_Quad = MeshX(1) % Nodes,     &
             Input_T_Quad = MeshX(2) % Nodes,     &
             Input_P_Quad = MeshX(3) % Nodes,     &
             Left_Limit   = -Half,                &
             Right_Limit  = +Half )

    CALL Poseidon_XCFC_Input_Sources2 &
           ( Local_S      = S,                    &
             Local_RE_Dim = nX(1),                &
             Local_TE_Dim = nX(2),                &
             Local_PE_Dim = nX(3),                &
             Local_RQ_Dim = nNodesX(1),           &
             Local_TQ_Dim = nNodesX(2),           &
             Local_PQ_Dim = nNodesX(3),           &
             Input_R_Quad = MeshX(1) % Nodes,     &
             Input_T_Quad = MeshX(2) % Nodes,     &
             Input_P_Quad = MeshX(3) % Nodes,     &
             Left_Limit   = -Half,                &
             Right_Limit  = +Half )

    ! --- Compute lapse and shift ---

    CALL Poseidon_XCFC_Run_Part2()

    CALL Poseidon_Return_Lapse &
         ( NE           = nX,               &
           NQ           = nNodesX,          &
           RQ_Input     = MeshX(1) % Nodes, &
           TQ_Input     = MeshX(2) % Nodes, &
           PQ_Input     = MeshX(3) % Nodes, &
           Left_Limit   = -Half,            &
           Right_Limit  = +Half,            &
           Return_Lapse = Tmp_Lapse )

    CALL Poseidon_Return_Shift &
         ( NE           = nX,               &
           NQ           = nNodesX,          &
           RQ_Input     = MeshX(1) % Nodes, &
           TQ_Input     = MeshX(2) % Nodes, &
           PQ_Input     = MeshX(3) % Nodes, &
           Left_Limit   = -Half,            &
           Right_Limit  = +Half,            &
           Return_Shift = Tmp_Shift )

    CALL Poseidon_Return_ExtrinsicCurvature &
         ( NE          = nX,               &
           NQ          = nNodesX,          &
           RQ_Input    = MeshX(1) % Nodes, &
           TQ_Input    = MeshX(2) % Nodes, &
           PQ_Input    = MeshX(3) % Nodes, &
           Left_Limit  = -Half,            &
           Right_Limit = +Half,            &
           Return_Kij  = Tmp_ExtrinsicCurvature )

    ! --- Copy data from Poseidon arrays to thornado arrays ---

    CALL ComputeGeometryFromPoseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             Tmp_Lapse, Tmp_Shift, Tmp_ExtrinsicCurvature, G )

    CALL SetBoundaryConditions &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G )

#endif

    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE ComputeGeometry_Poseidon


  ! --- PRIVATE Subroutines ---


  SUBROUTINE UpdateConformalFactorAndMetric &
    ( iX_B0, iX_E0, iX_B1, iX_E1, Psi, G )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: Psi(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):)
    REAL(DP), INTENT(inout) :: G  (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iNX, iNX1, iNX2
    REAL(DP) :: X1, X2

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNX = 1, nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

        G(iNX,iX1,iX2,iX3,iGF_Psi) = Psi    (iNX,iX1,iX2,iX3)
        G(iNX,iX1,iX2,iX3,iGF_h_1) = Psi(iNX,iX1,iX2,iX3)**2
        G(iNX,iX1,iX2,iX3,iGF_h_2) = Psi(iNX,iX1,iX2,iX3)**2 * X1
        G(iNX,iX1,iX2,iX3,iGF_h_3) = Psi(iNX,iX1,iX2,iX3)**2 * X1 * SIN( X2 )

      END DO ! iNX

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

    END DO ! iX1
    END DO ! iX2
    END DO ! iX3

  END SUBROUTINE UpdateConformalFactorAndMetric


  SUBROUTINE ComputeGeometryFromPoseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, Alpha, Beta_u, K_dd, G )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: Alpha (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):)
    REAL(DP), INTENT(in)    :: Beta_u(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in)    :: K_dd  (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(inout) :: G     (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iNX, iNX1, iNX2
    REAL(DP) :: X1, X2

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNX = 1, nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

        G(iNX,iX1,iX2,iX3,iGF_Alpha) = Alpha(iNX,iX1,iX2,iX3)

        G(iNX,iX1,iX2,iX3,iGF_Beta_1) = Beta_u(iNX,iX1,iX2,iX3,1)
        G(iNX,iX1,iX2,iX3,iGF_Beta_2) = Beta_u(iNX,iX1,iX2,iX3,2)
        G(iNX,iX1,iX2,iX3,iGF_Beta_3) = Beta_u(iNX,iX1,iX2,iX3,3)

        G(iNX,iX1,iX2,iX3,iGF_h_1) &
          = G(iNX,iX1,iX2,iX3,iGF_Psi)**2
        G(iNX,iX1,iX2,iX3,iGF_h_2) &
          = G(iNX,iX1,iX2,iX3,iGF_Psi)**2 * X1
        G(iNX,iX1,iX2,iX3,iGF_h_3) &
          = G(iNX,iX1,iX2,iX3,iGF_Psi)**2 * X1 * SIN( X2 )

        G(iNX,iX1,iX2,iX3,iGF_K_dd_11) = K_dd(iNX,iX1,iX2,iX3,1)
        G(iNX,iX1,iX2,iX3,iGF_K_dd_12) = K_dd(iNX,iX1,iX2,iX3,2)
        G(iNX,iX1,iX2,iX3,iGF_K_dd_13) = K_dd(iNX,iX1,iX2,iX3,3)
        G(iNX,iX1,iX2,iX3,iGF_K_dd_22) = K_dd(iNX,iX1,iX2,iX3,4)
        G(iNX,iX1,iX2,iX3,iGF_K_dd_23) = K_dd(iNX,iX1,iX2,iX3,5)
        G(iNX,iX1,iX2,iX3,iGF_K_dd_33) = K_dd(iNX,iX1,iX2,iX3,6)

      END DO ! iNX

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

    END DO ! iX1
    END DO ! iX2
    END DO ! iX3

  END SUBROUTINE ComputeGeometryFromPoseidon


  SUBROUTINE SetBoundaryConditions &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CALL SetBoundaryConditions_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G )

  END SUBROUTINE SetBoundaryConditions


  SUBROUTINE SetBoundaryConditions_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX2, iX3
    INTEGER  :: iNX1, jNX1, iNX2, iNX3
    INTEGER  :: iNX, jNX
    REAL(DP) :: X1, X2

    X2 = Half * Pi

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)

      DO iNX3 = 1, nNodesX(3)
      DO iNX2 = 1, nNodesX(2)
      DO iNX1 = 1, nNodesX(1)

        ! --- Inner Boundary: Reflecting ---

        jNX1 = ( nNodesX(1) - iNX1 ) + 1

        iNX = NodeNumberX( iNX1, iNX2, iNX3 )
        jNX = NodeNumberX( jNX1, iNX2, iNX3 )

        G(iNX,0,iX2,iX3,iGF_Alpha) &
          = G(jNX,1,iX2,iX3,iGF_Alpha)

        G(iNX,0,iX2,iX3,iGF_Psi) &
          = G(jNX,1,iX2,iX3,iGF_Psi)

        G(iNX,0,iX2,iX3,iGF_Beta_1) &
          = -G(jNX,1,iX2,iX3,iGF_Beta_1)

        G(iNX,0,iX2,iX3,iGF_h_1) &
          = G(jNX,1,iX2,iX3,iGF_h_1)

        G(iNX,0,iX2,iX3,iGF_h_2) &
          = G(jNX,1,iX2,iX3,iGF_h_2)

        G(iNX,0,iX2,iX3,iGF_h_3) &
          = G(jNX,1,iX2,iX3,iGF_h_3)

        ! --- Outer Boundary: Dirichlet ---

        X1 = NodeCoordinate( MeshX(1), nX(1)+1, iNX1 )

        G(iNX,nX(1)+1,iX2,iX3,iGF_Alpha) &
          = LapseFunction( X1, GravitationalMass )

        G(iNX,nX(1)+1,iX2,iX3,iGF_Psi) &
          = ConformalFactor( X1, GravitationalMass )

        G(iNX,nX(1)+1,iX2,iX3,iGF_Beta_1) &
          = Zero

        G(iNX,nX(1)+1,iX2,iX3,iGF_h_1) &
          = G(iNX,nX(1)+1,iX2,iX3,iGF_Psi)**2

        G(iNX,nX(1)+1,iX2,iX3,iGF_h_2) &
          = G(iNX,nX(1)+1,iX2,iX3,iGF_Psi)**2 * X1

        G(iNX,nX(1)+1,iX2,iX3,iGF_h_3) &
          = G(iNX,nX(1)+1,iX2,iX3,iGF_Psi)**2 * X1 * SIN( X2 )

      END DO
      END DO
      END DO

      CALL ComputeGeometryX_FromScaleFactors( G(:,0      ,iX2,iX3,:) )
      CALL ComputeGeometryX_FromScaleFactors( G(:,nX(1)+1,iX2,iX3,:) )

    END DO
    END DO

  END SUBROUTINE SetBoundaryConditions_X1


  SUBROUTINE ComputeGravitationalMass &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mg )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)  :: &
      Mg(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):)

    INTEGER  :: iX1, iX2, iX3
    REAL(DP) :: d3X

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width(1:nX(1)), &
        dX2 => MeshX(2) % Width(1:nX(2)), &
        dX3 => MeshX(3) % Width(1:nX(3)) )

    ! --- Assuming 1D spherical symmetry ---

    GravitationalMass = Zero

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)

      d3X = Two / Pi * dX1(iX1) * dX2(iX2) * dX3(iX3)

      GravitationalMass &
        = GravitationalMass + d3X                  &
            * SUM( WeightsX_q                      &
                     * Mg(:,iX1,iX2,iX3)           &
                     * G (:,iX1,iX2,iX3,iGF_Alpha) &
                     * G (:,iX1,iX2,iX3,iGF_SqrtGm) )

    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeGravitationalMass


END MODULE GravitySolutionModule_CFA_Poseidon
