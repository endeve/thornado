MODULE GravitySolutionModule_Newtonian_Poseidon

  USE KindModule, ONLY: &
    DP, Pi, FourPi
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

  USE Poseidon_Main_Module, ONLY: &
    Poseidon_Initialize, &
    Poseidon_Newtonian_Source_Input, &
    Poseidon_Set_Uniform_Boundary_Condition, &
    Poseidon_Run,  &
    Poseidon_Newtonian_Potential_Output, &
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

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'InitializeGravitySolver_Newtonian_Poseidon'
    WRITE(*,'(A4,A)') '', 'Only implemented for 1D spherical symmetry.'
    WRITE(*,*)

    CALL Poseidon_Initialize &
         ( FEM_Degree_Input          &
             = MAX( 1, nNodes - 1 ), &
           L_Limit_Input = 0,        &
           Inner_Radius = xL(1),     &
           Outer_Radius = xR(1),     &
           R_Elements_Input = nX(1), &
           T_Elements_Input = nX(2), &
           P_Elements_Input = nX(3), &
           Input_Delta_R_Vector      &
             = MeshX(1) % Width(1:nX(1)) )

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

    CALL Poseidon_Newtonian_Source_Input      &
           ( Left_Limit   = - 0.5_DP,         &
             Right_Limit  = + 0.5_DP,         &
             Num_Nodes    = nNodes,           &
             Input_R_Quad = MeshX(1) % Nodes, &
             Rho = D(:,iX_B0(1):iX_E0(1),     &
                       iX_B0(2):iX_E0(2),     &
                       iX_B0(3):iX_E0(3)) )

    CALL Poseidon_Set_Uniform_Boundary_Condition &
           ( BC_Location_Input = 'I',            &
             BC_Type_Input     = 'N',            &
             BC_Value_Input    = 0.0_DP )

    CALL Poseidon_Set_Uniform_Boundary_Condition &
           ( BC_Location_Input = 'O',            &
             BC_Type_Input     = 'D',            &
             BC_Value_Input    = - BaryonMass / xR(1) )

    CALL Poseidon_Run

    CALL Poseidon_Newtonian_Potential_Output   &
           ( Left_Limit    = - 0.5_DP,         &
             Right_Limit   = + 0.5_DP,         &
             Num_Nodes     = nNodes,           &
             Output_R_Quad = MeshX(1) % Nodes, &
             Potential = G(:,1:nX(1),1:nX(2),1:nX(3),iGF_Phi_N) )

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

    Mass = 0.0_DP
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
