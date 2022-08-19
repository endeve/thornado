MODULE GravitySolutionModule_Newtonian_Poseidon

  USE KindModule, ONLY: &
    DP, Zero, Half, Pi, FourPi
  USE UnitsModule, ONLY: &
    Kilogram
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, &
    nNodes, xL, xR
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    iGF_Phi_N, iGF_SqrtGm

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTON

  ! --- Poseidon Modules --------------------

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
  
  ! -----------------------------------------

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_Newtonian_Poseidon
  PUBLIC :: FinalizeGravitySolver_Newtonian_Poseidon
  PUBLIC :: SolveGravity_Newtonian_Poseidon

CONTAINS


  SUBROUTINE InitializeGravitySolver_Newtonian_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTON

    CALL Initialize_Poseidon &
           ( FEM_Degree_Option            = MAX( 1, nNodes - 1 ),      &
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
             Newtonian_Mode_Option        = .TRUE.,                    &
             Verbose_Option               = .FALSE. )

#endif

  END SUBROUTINE InitializeGravitySolver_Newtonian_Poseidon


  SUBROUTINE FinalizeGravitySolver_Newtonian_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTON

    CALL Poseidon_Close

#endif

  END SUBROUTINE FinalizeGravitySolver_Newtonian_Poseidon


  SUBROUTINE SolveGravity_Newtonian_Poseidon &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):)

    REAL(DP) :: BaryonMass

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTON

    CALL ComputeTotalBaryonMass &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, D, BaryonMass )

    CALL Poseidon_Input_Sources                     &
           ( Input_Rho = D(:,iX_B0(1):iX_E0(1),     &
                             iX_B0(2):iX_E0(2),     &
                             iX_B0(3):iX_E0(3)) )

    CALL Poseidon_Set_Uniform_Boundary_Conditions   &
           ( BC_Location_Input = 'I',               &
             BC_Type_Input     = 'N',               &
             BC_Value_Input    = Zero )

    CALL Poseidon_Set_Uniform_Boundary_Conditions   &
           ( BC_Location_Input = 'O',               &
             BC_Type_Input     = 'D',               &
             BC_Value_Input    = - BaryonMass / xR(1) )

    CALL Poseidon_Run

    CALL Poseidon_Return_Newtonian_Potential   &
           ( Return_Potential = G(:,1:nX(1),1:nX(2),1:nX(3),iGF_Phi_N) )

    CALL SetBoundaryConditions &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, BaryonMass )

#endif

  END SUBROUTINE SolveGravity_Newtonian_Poseidon


  SUBROUTINE ComputeTotalBaryonMass &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, D, Mass )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)  :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):)
    REAL(DP), INTENT(out) :: &
      Mass

    INTEGER :: iX1, iX2, iX3

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width(1:nX(1)) )

    ! --- Assuming 1D spherical symmetry ---

    Mass = Zero
    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)

      Mass &
        = Mass &
            + FourPi * dX1(iX1) &
                * SUM( WeightsX_q(:) * D(:,iX1,iX2,iX3) &
                         * G(:,iX1,iX2,iX3,iGF_SqrtGm) )

    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTotalBaryonMass


  SUBROUTINE SetBoundaryConditions &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, BaryonMass )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      BaryonMass

    CALL SetBoundaryConditions_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, BaryonMass )

  END SUBROUTINE SetBoundaryConditions


  SUBROUTINE SetBoundaryConditions_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, BaryonMass )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      BaryonMass

    INTEGER  :: iX2, iX3
    INTEGER  :: iNodeX1, jNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iNodeX, jNodeX
    REAL(DP) :: X1

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        ! --- Inner Boundary: Reflecting ---

        jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

        iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
        jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

        G(iNodeX,0,iX2,iX3,iGF_Phi_N) &
          = G(jNodeX,1,iX2,iX3,iGF_Phi_N)

        ! --- Outer Boundary: Dirichlet ---

        X1 = NodeCoordinate( MeshX(1), nX(1)+1, iNodeX1 )

        G(iNodeX,nX(1)+1,iX2,iX3,iGF_Phi_N) &
          = - BaryonMass / X1

      END DO
      END DO
      END DO

    END DO
    END DO

  END SUBROUTINE SetBoundaryConditions_X1


END MODULE GravitySolutionModule_Newtonian_Poseidon
