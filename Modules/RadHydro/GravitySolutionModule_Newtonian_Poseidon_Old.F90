MODULE GravitySolutionModule_Newtonian_Poseidon_Old

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, &
    nNodes, xL, xR
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Phi_N
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D
  USE GravitySolutionUtilitiesModule, ONLY: &
    ComputeTotalBaryonMass

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


  SUBROUTINE SolveGravity_Newtonian_Poseidon

    REAL(DP) :: BaryonMass

    CALL ComputeTotalBaryonMass( BaryonMass )

#ifdef GRAVITY_SOLVER_POSEIDON_NEWTON

    CALL Poseidon_Newtonian_Source_Input      &
           ( Left_Limit   = - 0.5_DP,         &
             Right_Limit  = + 0.5_DP,         &
             Num_Nodes    = nNodes,           &
             Input_R_Quad = MeshX(1) % Nodes, &
             Rho = uCF(:,1:nX(1),1:nX(2),1:nX(3),iCF_D) )

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
             Potential = uGF(:,1:nX(1),1:nX(2),1:nX(3),iGF_Phi_N) )

#endif

    CALL SetBoundaryConditions( BaryonMass )

  END SUBROUTINE SolveGravity_Newtonian_Poseidon


  SUBROUTINE SetBoundaryConditions( BaryonMass )

    REAL(DP), INTENT(in) :: BaryonMass

    CALL SetBoundaryConditions_X1( BaryonMass )

  END SUBROUTINE SetBoundaryConditions


  SUBROUTINE SetBoundaryConditions_X1( BaryonMass )

    REAL(DP), INTENT(in) :: BaryonMass

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

              uGF(iNodeX,0,iX2,iX3,iGF_Phi_N) &
                = uGF(jNodeX,1,iX2,iX3,iGF_Phi_N)

              ! --- Outer Boundary: Dirichlet ---

              X1 = NodeCoordinate( MeshX(1), nX(1)+1, iNodeX1 )

              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Phi_N) &
                = - BaryonMass / X1

            END DO
          END DO
        END DO

      END DO
    END DO

  END SUBROUTINE SetBoundaryConditions_X1


END MODULE GravitySolutionModule_Newtonian_Poseidon_Old
