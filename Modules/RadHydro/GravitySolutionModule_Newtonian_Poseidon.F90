MODULE GravitySolutionModule_Newtonian_Poseidon

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramHeaderModule, ONLY: &
    nX, nNodes, xL, xR
  USE MeshModule, ONLY: &
    MeshX
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

  END SUBROUTINE SolveGravity_Newtonian_Poseidon


END MODULE GravitySolutionModule_Newtonian_Poseidon
