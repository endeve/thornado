MODULE GravitySolutionModule_CFA_Poseidon

  USE KindModule, ONLY: &
    DP, &
    Pi, &
    FourPi, &
    Zero, &
    Half, &
    One, &
    Two, &
    Three
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
    NodeNumberX, &
    thornado_abort
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
    iGF_Gm_dd_33
  USE GeometryComputationModule, ONLY: &
    LapseFunction, &
    ConformalFactor, &
    ComputeGeometryX_FromScaleFactors
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler,  &
    Timer_GravitySolver

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

  ! --- Poseidon Modules ---

  USE Initialization_Poseidon, ONLY: &
    Initialize_Poseidon
  USE Poseidon_Main_Module, ONLY: &
    Poseidon_Run, &
    Poseidon_Close, &
    Poseidon_CFA_Set_Uniform_Boundary_Conditions
  USE Source_Input_Module, ONLY: &
    Poseidon_Input_Sources
  USE Variables_Functions, ONLY: &
    Calc_1D_CFA_Values
  USE Initial_Guess_Module, ONLY: &
    Poseidon_Init_FlatGuess

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_CFA_Poseidon
  PUBLIC :: FinalizeGravitySolver_CFA_Poseidon
  PUBLIC :: SolveGravity_CFA_Poseidon


CONTAINS


  SUBROUTINE InitializeGravitySolver_CFA_Poseidon

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
         ( Units_Option       = 'G',                       &
           Dimensions_Option  = 3,                         &
           FEM_Degree_Option  = MAX( 1, nNodes - 1 ),      &
           L_Limit_Option     = 0,                         &
           Domain_Edge_Option = [ xL(1), xR(1) ],          &
           NE_Option          = nX,                        &
           NQ_Option          = [ nNodes, 1, 1 ],          &
           dr_Option          = MeshX(1) % Width(1:nX(1)), &
           dt_Option          = MeshX(2) % Width(1:nX(2)), &
           dp_Option          = MeshX(3) % Width(1:nX(3)), &
           Method_Flag_Option = 2,                         &
           Print_Setup_Option = .TRUE.,                    &
           Verbose_Option     = .FALSE.)

#endif

  END SUBROUTINE InitializeGravitySolver_CFA_Poseidon


  SUBROUTINE FinalizeGravitySolver_CFA_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    CALL Poseidon_Close()

#endif

  END SUBROUTINE FinalizeGravitySolver_CFA_Poseidon


  SUBROUTINE SolveGravity_CFA_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Sources )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout)    :: &
      G      (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      Sources(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    REAL(DP)         :: Psi_BC, AlphaPsi_BC, GravitationalMass
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)

    REAL(DP) :: Tmp_ConFact (nDOFX,iX_B0(1):iX_E0(1), &
                                   iX_B0(2):iX_E0(2), &
                                   iX_B0(3):iX_E0(3))
    REAL(DP) :: Tmp_Lapse   (nDOFX,iX_B0(1):iX_E0(1), &
                                   iX_B0(2):iX_E0(2), &
                                   iX_B0(3):iX_E0(3))
    REAL(DP) :: Tmp_Shift_X1(nDOFX,iX_B0(1):iX_E0(1), &
                                   iX_B0(2):iX_E0(2), &
                                   iX_B0(3):iX_E0(3))

    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    ! --- Set Boundary Values ---

    CALL ComputeGravitationalMass &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, Sources, GravitationalMass )

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

    ! --- Set matter sources ---

    CALL Poseidon_Input_Sources &
           ( 0, 0, 0,                             &
             Local_E      = Sources(:,:,:,:,1),   &
             Local_S      = Sources(:,:,:,:,2),   &
             Local_Si     = Sources(:,:,:,:,3:5), &
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

    CALL Poseidon_Run()

    CALL Calc_1D_CFA_Values &
           ( Num_RE_Input  = nX(1),            &
             Num_RQ_Input  = nNodesX(1),       &
             RQ_Input      = MeshX(1) % Nodes, &
             Left_Limit    = -Half,            &
             Right_Limit   = +Half,            &
             CFA_Lapse     = Tmp_Lapse,        &
             CFA_ConFactor = Tmp_ConFact,      &
             CFA_Shift     = Tmp_Shift_X1 )

    ! --- Copy data from Poseidon arrays to thornado arrays ---

    CALL ComputeGeometryFromPoseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             Tmp_Lapse, Tmp_ConFact, Tmp_Shift_X1, G )

    CALL SetBoundaryConditions &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, GravitationalMass )

#endif

    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE SolveGravity_CFA_Poseidon


  ! --- PRIVATE Subroutines ---


  SUBROUTINE ComputeGeometryFromPoseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, Alpha, Psi, Beta_X1, G )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: Alpha  (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):)
    REAL(DP), INTENT(in)    :: Psi    (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):)
    REAL(DP), INTENT(in)    :: Beta_X1(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):)
    REAL(DP), INTENT(inout) :: G      (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

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

        G(iNX,iX1,iX2,iX3,iGF_Alpha)  = Alpha  (iNX,iX1,iX2,iX3)
        G(iNX,iX1,iX2,iX3,iGF_Psi)    = Psi    (iNX,iX1,iX2,iX3)
        G(iNX,iX1,iX2,iX3,iGF_Beta_1) = Beta_X1(iNX,iX1,iX2,iX3)
        G(iNX,iX1,iX2,iX3,iGF_Beta_2) = Zero
        G(iNX,iX1,iX2,iX3,iGF_Beta_3) = Zero
        G(iNX,iX1,iX2,iX3,iGF_h_1)    = Psi(iNX,iX1,iX2,iX3)**2
        G(iNX,iX1,iX2,iX3,iGF_h_2)    = Psi(iNX,iX1,iX2,iX3)**2 * X1
        G(iNX,iX1,iX2,iX3,iGF_h_3)    = Psi(iNX,iX1,iX2,iX3)**2 * X1 * SIN( X2 )

      END DO ! iNX

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

    END DO ! iX1
    END DO ! iX2
    END DO ! iX3

  END SUBROUTINE ComputeGeometryFromPoseidon


  SUBROUTINE SetBoundaryConditions &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      Mass

    CALL SetBoundaryConditions_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

  END SUBROUTINE SetBoundaryConditions


  SUBROUTINE SetBoundaryConditions_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      Mass

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
          = LapseFunction( X1, Mass )

        G(iNX,nX(1)+1,iX2,iX3,iGF_Psi) &
          = ConformalFactor( X1, Mass )

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
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Sources, Mass )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G      (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)  :: &
      Sources(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(out) :: &
      Mass

    INTEGER  :: iX1, iX2, iX3
    REAL(DP) :: d3X

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width(1:nX(1)), &
        dX2 => MeshX(2) % Width(1:nX(2)), &
        dX3 => MeshX(3) % Width(1:nX(3)) )

    ! --- Assuming 1D spherical symmetry ---

    Mass = Zero

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)

      d3X = Two / Pi * dX1(iX1) * dX2(iX2) * dX3(iX3)

      Mass &
        = Mass + d3X                                     &
            * SUM( WeightsX_q                            &
                     * Sources(:,iX1,iX2,iX3,6)          &
                     * G      (:,iX1,iX2,iX3,iGF_Alpha)  &
                     * G      (:,iX1,iX2,iX3,iGF_SqrtGm) )

    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeGravitationalMass


END MODULE GravitySolutionModule_CFA_Poseidon
