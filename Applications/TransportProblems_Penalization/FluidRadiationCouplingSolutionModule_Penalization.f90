MODULE FluidRadiationCouplingSolutionModule_Penalization

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE, &
    nDOF
  USE FluidFieldsModule, ONLY: &
    nPF, uPF, iPF_D, &
    nAF, uAF, iAF_T, iAF_Ye
  USE RadiationFieldsModule, ONLY: &
    nPR, uPR, iPR_D
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeNodesX, &
    InitializeWeights, &
    MapForward_FluidField, &
    MapBackward_FluidField, &
    MapForward_RadiationField, &
    MapBackward_RadiationField

  IMPLICIT NONE
  PRIVATE

  INTEGER                                 :: nNodesX_G
  INTEGER                                 :: nNodesE_G
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: D_N, T_N, Y_N
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: X_N
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: J_N, H1_N, H2_N, H3_N

  PUBLIC :: ComputeRHS_Penalization

CONTAINS


  SUBROUTINE ComputeRHS_Penalization( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    CALL InitializeFluidRadiationCoupling

    CALL FinalizeFluidRadiationCoupling

  END SUBROUTINE ComputeRHS_Penalization


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

    ALLOCATE( X_N(nNodesX_G,3) )
    CALL InitializeNodesX( X_N )

    ALLOCATE( D_N(nNodesX_G) )
    CALL MapForward_FluidField &
           ( uPF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iPF_D), &
             D_N(1:nNodesX_G) )

    ALLOCATE( T_N(nNodesX_G) )
    ALLOCATE( Y_N(nNodesX_G) )

    ALLOCATE( J_N(nNodesE_G, nNodesX_G) )
    CALL MapForward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_D,1), &
             J_N(1:nNodesE_G,1:nNodesX_G) )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( X_N, D_N, T_N, Y_N )
    DEALLOCATE( J_N )

  END SUBROUTINE FinalizeFluidRadiationCoupling


END MODULE FluidRadiationCouplingSolutionModule_Penalization
