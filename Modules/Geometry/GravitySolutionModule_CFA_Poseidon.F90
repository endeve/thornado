MODULE GravitySolutionModule_CFA_Poseidon

  USE KindModule, ONLY: &
    DP,     &
    Pi,     &
    FourPi, &
    Zero,   &
    Half,   &
    One,    &
    Two
  USE ProgramHeaderModule, ONLY: &
    nX,      &
    nNodesX, &
    nDOFX,   &
    nNodes,  &
    xL,      &
    xR
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    NodeNumberTableX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryComputationModule, ONLY: &
    LapseFunction,   &
    ConformalFactor, &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModule, ONLY: &
    iGF_Phi_N,    &
    iGF_SqrtGm,   &
    iGF_Alpha,    &
    iGF_Psi,      &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3,   &
    iGF_h_1,      &
    iGF_h_2,      &
    iGF_h_3,      &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE UnitsModule, ONLY: &
    SolarMass,centimeter,gram,second,erg,kilometer
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler,  &
    Timer_GravitySolver

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

  ! --- Poseidon Modules --------------------

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

  USE FP_Initial_Guess_Module, ONLY: &
    Init_FP_Guess_Flat

  ! -----------------------------------------

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
           Verbose_Option     = .TRUE. )

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
    REAL(DP)         :: Tmp_Lapse  (nDOFX,nX(1),nX(2),nX(3)), &
                        Tmp_ConFact(nDOFX,nX(1),nX(2),nX(3)), &
                        Tmp_Shift  (nDOFX,nX(1),nX(2),nX(3))

    INTEGER  :: iX1, iX2, iX3, iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3

    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    ! Set Source Values !
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

    ! --- Set Boundary Values ---

    CALL ComputeGravitationalMass &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, Sources, GravitationalMass )

    Psi_BC      = ConformalFactor( xR(1), GravitationalMass )
    AlphaPsi_BC = LapseFunction  ( xR(1), GravitationalMass ) * Psi_BC

    INNER_BC_TYPES = [ "N", "N", "N", "N", "N" ]
    OUTER_BC_TYPES = [ "D", "D", "D", "D", "D" ]

    INNER_BC_VALUES = [ Zero  , Zero       , Zero, Zero, Zero ]
    OUTER_BC_VALUES = [ Psi_BC, AlphaPsi_BC, Zero, Zero, Zero ]

    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions &
           ( "I", INNER_BC_TYPES, INNER_BC_VALUES )
    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions &
           ( "O", OUTER_BC_TYPES, OUTER_BC_VALUES)

    CALL Init_FP_Guess_Flat() ! Possibly move this to init call

    CALL Poseidon_Run()

    CALL Calc_1D_CFA_Values &
           ( Num_RE_Input  = nX(1),            &
             Num_RQ_Input  = nNodesX(1),       &
             RQ_Input      = MeshX(1) % Nodes, &
             Left_Limit    = -Half,            &
             Right_Limit   = +Half,            &
             CFA_Lapse     = Tmp_Lapse,        &
             CFA_ConFactor = Tmp_ConFact,      &
             CFA_Shift     = Tmp_Shift         )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

        G(iNodeX,iX1,iX2,iX3,iGF_Alpha)  &
          = Tmp_Lapse  (iNodeX,iX1,iX2,iX3)
        G(iNodeX,iX1,iX2,iX3,iGF_Psi)    &
          = Tmp_ConFact(iNodeX,iX1,iX2,iX3)
        G(iNodeX,iX1,iX2,iX3,iGF_Beta_1) &
          = Tmp_Shift  (iNodeX,iX1,iX2,iX3)
        G(iNodeX,iX1,iX2,iX3,iGF_Beta_2) &
          = Zero
        G(iNodeX,iX1,iX2,iX3,iGF_Beta_3) &
          = Zero
        G(iNodeX,iX1,iX2,iX3,iGF_h_1)    &
          = Tmp_ConFact(iNodeX,iX1,iX2,iX3)**2
        G(iNodeX,iX1,iX2,iX3,iGF_h_2)    &
          = Tmp_ConFact(iNodeX,iX1,iX2,iX3)**2 * X1
        G(iNodeX,iX1,iX2,iX3,iGF_h_3)    &
          = Tmp_ConFact(iNodeX,iX1,iX2,iX3)**2 * X1 * SIN( X2 )

      END DO ! iNodeX

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

    END DO ! iX1
    END DO ! iX2
    END DO ! iX3

    CALL SetBoundaryConditions &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, GravitationalMass )

#endif

    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE SolveGravity_CFA_Poseidon


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
    INTEGER  :: iNodeX1, jNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iNodeX, jNodeX
    REAL(DP) :: X1, X2

    X2 = Half * Pi

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        ! --- Inner Boundary: Reflecting ---

        jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

        iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
        jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

        G(iNodeX,0,iX2,iX3,iGF_Alpha) &
          = G(jNodeX,1,iX2,iX3,iGF_Alpha)

        G(iNodeX,0,iX2,iX3,iGF_Psi) &
          = G(jNodeX,1,iX2,iX3,iGF_Psi)

        G(iNodeX,0,iX2,iX3,iGF_Beta_1) &
          = -G(jNodeX,1,iX2,iX3,iGF_Beta_1)

        G(iNodeX,0,iX2,iX3,iGF_h_1) &
          = G(jNodeX,1,iX2,iX3,iGF_h_1)

        G(iNodeX,0,iX2,iX3,iGF_h_2) &
          = G(jNodeX,1,iX2,iX3,iGF_h_2)

        G(iNodeX,0,iX2,iX3,iGF_h_3) &
          = G(jNodeX,1,iX2,iX3,iGF_h_3)

        ! --- Outer Boundary: Dirichlet ---

        X1 = NodeCoordinate( MeshX(1), nX(1)+1, iNodeX1 )

        G(iNodeX,nX(1)+1,iX2,iX3,iGF_Alpha) &
          = LapseFunction( X1, Mass )

        G(iNodeX,nX(1)+1,iX2,iX3,iGF_Psi) &
          = ConformalFactor( X1, Mass )

        G(iNodeX,nX(1)+1,iX2,iX3,iGF_Beta_1) &
          = Zero

        G(iNodeX,nX(1)+1,iX2,iX3,iGF_h_1) &
          = G(iNodeX,nX(1)+1,iX2,iX3,iGF_Psi)**2

        G(iNodeX,nX(1)+1,iX2,iX3,iGF_h_2) &
          = G(iNodeX,nX(1)+1,iX2,iX3,iGF_Psi)**2 * X1

        G(iNodeX,nX(1)+1,iX2,iX3,iGF_h_3) &
          = G(iNodeX,nX(1)+1,iX2,iX3,iGF_Psi)**2 * X1 * SIN( X2 )

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
