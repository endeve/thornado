MODULE TwoMoment_DiscretizationModule_Streaming_FMC

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, &
    SqrtTiny, Pi
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply, &
    EigenvaluesSymmetric3
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, nDOFX_X2, nDOFX_X3
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, LX_X1_Dn, LX_X1_Up
  USE ReferenceElementModule, ONLY: &
    nDOF_E, &
    nDOF_X1, nDOF_X2, nDOF_X3, &
    Weights_q, &
    Weights_E, &
    Weights_X1
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_E_Dn,  L_E_Up, &
    dLdE_q, &
    L_X1_Dn, L_X1_Up, &
    dLdX1_q
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE, iGE_Ep1, iGE_Ep2, iGE_Ep3
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_h_1, iGF_h_2, iGF_h_3, &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nPF, iPF_V1, iPF_V2, iPF_V3
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nSpecies, &
    nCM, iCM_E, iCM_F1, iCM_F2, iCM_F3, &
    nPM, iPM_J, iPM_H1, iPM_H2, iPM_H3
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE TwoMoment_UtilitiesModule_FMC, ONLY: &
    ComputeConserved_TwoMoment_FMC, &
    ComputePrimitive_TwoMoment_Richardson_FMC, &
    Flux_X1, &
    Flux_E, &
    NumericalFlux_LLF, &
    ComputeWeakDerivatives_X1, &
    FaceVelocity_X1

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit

  ! --- Geometry Fields ---

  REAL(DP), TARGET, DIMENSION(:,:,:,:,:), ALLOCATABLE :: &
    uGF_K, uGF_F

  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K, &
    Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F

  ! --- Primitive Fluid Fields ---

  REAL(DP), TARGET, DIMENSION(:,:,:,:,:), ALLOCATABLE :: &
    uPF_K, uPF_L, uPF_R

  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    uV1_K, uV2_K, uV3_K, &
    uV1_L, uV2_L, uV3_L, &
    uV1_R, uV2_R, uV3_R

  ! --- Conserved Two-Moment Fields ---

  REAL(DP), TARGET, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: &
    uCM_K, uCM_L, uCM_R

  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    E_K, F_d_1_K, F_d_2_K, F_d_3_K, &
    E_L, F_d_1_L, F_d_2_L, F_d_3_L, &
    E_R, F_d_1_R, F_d_2_R, F_d_3_R

  ! --- Primitive Two-Moment Fields

  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    J_K, H_u_1_K, H_u_2_K, H_u_3_K, &
    J_L, H_u_1_L, H_u_2_L, H_u_3_L, &
    J_R, H_u_1_R, H_u_2_R, H_u_3_R

  ! --- Primitive Fluid Velocities ---

  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    V_u_1_F, V_u_2_F, V_u_3_F, &
    V_u_1_K, V_u_2_K, V_u_3_K

  ! --- Fluxes ---

  REAL(DP), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: &
    NumericalFlux, Flux_q

  ! --- Increment ---
  REAL(DP), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: &
    dU_Z

  ! --- Eigenvalues ---

  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: &
    Alpha

  ! --- Derivatives ---

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: &
    dV_u_dX1, dV_u_dX2, dV_u_dX3, &
    dV_d_dX1, dV_d_dX2, dV_d_dX3, &
    dGm_dd_dX1, dGm_dd_dX2, dGm_dd_dX3

  ! --- ? ---

  INTEGER,  DIMENSION(:), ALLOCATABLE :: &
    nIterations_L, nIterations_R, nIterations_K

  INTEGER,  DIMENSION(:), ALLOCATABLE :: &
    PositionIndexZ_F, PositionIndexZ_K

  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: &
    IndexTableZ_F, IndexTableZ_K

  INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
  INTEGER :: nZ   (4), nK_X , nK_Z , nNodesX_K, nNodesZ_K
  INTEGER :: nZ_E (4), nE_X , nE_Z , nNodesZ_E
  INTEGER :: nZ_X1(4), nX1_X, nX1_Z, nNodesX_X1, nNodesZ_X1
  INTEGER :: nZ_X2(4), nX2_X, nX2_Z, nNodesX_X2, nNodesZ_X2
  INTEGER :: nZ_X3(4), nX3_X, nX3_Z, nNodesX_X3, nNodesZ_X3

CONTAINS

  SUBROUTINE ComputeIncrement_TwoMoment_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_M, dU_M )

    ! --- Input/Output variables ---
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: GE(1:nDOFE, iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: GX(1:nDOFX, &
                               iZ_B1(2):iZ_E1(2), &
                               iZ_B1(3):iZ_E1(3), &
                               iZ_B1(4):iZ_E1(4), &
                               1:nGF)
    REAL(DP), INTENT(inout) :: U_F (1:nDOFX, &
                                    iZ_B1(2):iZ_E1(2), &
                                    iZ_B1(3):iZ_E1(3), &
                                    iZ_B1(4):iZ_E1(4), &
                                    1:nPF)
    REAL(DP), INTENT(inout) :: U_M(1:nDOFZ, &
                                   iZ_B1(1):iZ_E1(1), &
                                   iZ_B1(2):iZ_E1(2), &
                                   iZ_B1(3):iZ_E1(3), &
                                   iZ_B1(4):iZ_E1(4), &
                                   1:nCM, &
                                   1:nSpecies)
    REAL(DP), INTENT(out) :: dU_M(1:nDOFZ, &
                                  iZ_B1(1):iZ_E1(1), &
                                  iZ_B1(2):iZ_E1(2), &
                                  iZ_B1(3):iZ_E1(3), &
                                  iZ_B1(4):iZ_E1(4), &
                                  1:nCM, &
                                  1:nSpecies)

    ! --- Local variables ---
    INTEGER :: iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iCM, iS

    Write(*,*)
    print *,'ComputeIncrement_TwoMoment_Explicit'

    ASSOCIATE ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
                dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    CALL ApplyBoundaryConditions_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U_F )

    CALL ApplyBoundaryConditions_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_M )

    CALL InitializeIncrement_TwoMoment_Explicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    DO iS  = 1, nSpecies
    DO iCM = 1, nCM
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)
  
        DO iNodeZ = 1, nDOFZ
  
          dU_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) = Zero
  
        END DO
  
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL ComputeIncrement_Divergence_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_M, dU_M )
           
    CALL ComputeIncrement_ObserverCorrections &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_M, dU_M )

    ! --- Multiply Inverse Mass Matrix ---

    DO iS  = 1, nSpecies
    DO iCM = 1, nCM
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        dU_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) &
          = dU_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) &
              / ( Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
                    * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm) &
                    * dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_TwoMoment_Explicit

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit


  SUBROUTINE InitializeIncrement_TwoMoment_Explicit( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
    
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)

    nZ         = iZ_E0 - iZ_B0 + 1           ! Number of Elements per Phase Space Dimension
    nK_X       = PRODUCT( nZ(2:4) )          ! Number of Elements in Position Space
    nK_Z       = nSpecies * nZ(1) * nK_X     ! Number of Elements in Phase Space
    nNodesX_K  = nDOFX * nK_X                ! Number of Nodes in Position Space
    nNodesZ_K  = nDOFZ * nK_Z                ! Number of Nodes in Phase Space

    nZ_E  = nZ + [1,0,0,0]                   ! Number of E  Faces per Phase Space Dimension
    nZ_X1 = nZ + [0,1,0,0]                   ! Number of X1 Faces per Phase Space Dimension
    nZ_X2 = nZ + [0,0,1,0]                   ! Number of X2 Faces per Phase Space Dimension
    nZ_X3 = nZ + [0,0,0,1]                   ! Number of X3 Faces per Phase Space Dimension

    nE_X  = PRODUCT( nZ_E (2:4) )            ! Number of E  Faces in Position Space
    nX1_X = PRODUCT( nZ_X1(2:4) )            ! Number of X1 Faces in Position Space
    nX2_X = PRODUCT( nZ_X2(2:4) )            ! Number of X2 Faces in Position Space
    nX3_X = PRODUCT( nZ_X3(2:4) )            ! Number of X3 Faces in Position Space

    nE_Z  = nSpecies * nZ_E (1) * nE_X       ! Number of E  Faces in Phase Space
    nX1_Z = nSpecies * nZ_X1(1) * nX1_X      ! Number of X1 Faces in Phase Space
    nX2_Z = nSpecies * nZ_X2(1) * nX2_X      ! Number of X2 Faces in Phase Space
    nX3_Z = nSpecies * nZ_X3(1) * nX3_X      ! Number of X3 Faces in Phase Space

    nNodesX_X1 = nDOFX_X1 * nX1_X            ! Number of Nodes on X1 Faces in Position space
    nNodesX_X2 = nDOFX_X2 * nX2_X            ! Number of Nodes on X2 Faces in Position space
    nNodesX_X3 = nDOFX_X3 * nX3_X            ! Number of Nodes on X3 Faces in Position space

    nNodesZ_E  = nDOF_E  * nE_Z              ! Number of Nodes on X1 Faces in Phase space
    nNodesZ_X1 = nDOF_X1 * nX1_Z             ! Number of Nodes on X1 Faces in Phase space
    nNodesZ_X2 = nDOF_X2 * nX2_Z             ! Number of Nodes on X1 Faces in Phase space
    nNodesZ_X3 = nDOF_X3 * nX3_Z             ! Number of Nodes on X1 Faces in Phase space

    ALLOCATE( V_u_1_K(nNodesX_K) )
    ALLOCATE( V_u_2_K(nNodesX_K) )
    ALLOCATE( V_u_3_K(nNodesX_K) )

    ALLOCATE( J_K   (nNodesZ_K) )
    ALLOCATE( H_u_1_K(nNodesZ_K) )
    ALLOCATE( H_u_2_K(nNodesZ_K) )
    ALLOCATE( H_u_3_K(nNodesZ_K) )

  END SUBROUTINE InitializeIncrement_TwoMoment_Explicit


  SUBROUTINE FinalizeIncrement_TwoMoment_Explicit

    DEALLOCATE( V_u_1_K, V_u_2_K, V_u_3_K )
    DEALLOCATE( J_K, H_u_1_K, H_u_2_K, H_u_3_K )

  END SUBROUTINE FinalizeIncrement_TwoMoment_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_M, dU_M )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---
    ! --- Input/Output variables ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nPF)
    REAL(DP), INTENT(in)    :: &
      U_M (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCM, &
           1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_M(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCM, &
           1:nSpecies)

    ! --- Geometry Fields ---
    REAL(DP) :: &
      uGF_K(nDOFX, &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            iZ_B0(2)-1:iZ_E0(2)+1, &
            nGF)
    REAL(DP) :: &
      uGF_F(nDOFX_X1, &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            iZ_B0(2):iZ_E0(2)+1, &
            nGF)

    ! --- Primitive Fluid Fields ---

    REAL(DP) :: &
      uPF_K(nDOFX, &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            iZ_B0(2)-1:iZ_E0(2)+1, &
            nPF)
    REAL(DP) :: &
      uPF_L(nDOFX_X1, &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            iZ_B0(2)  :iZ_E0(2)+1, &
            nPF)
    REAL(DP) :: &
      uPF_R(nDOFX_X1, &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            iZ_B0(2)  :iZ_E0(2)+1, &
            nPF)

    ! --- Conserved Two-Moment Fields ---

    REAL(DP) :: &
      uCM_K(nDOFZ, &
            iZ_B0(1)  :iZ_E0(1)  , &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            nSpecies, &
            iZ_B0(2)-1:iZ_E0(2)+1, &
            nCM)
    REAL(DP) :: &
      uCM_L(nDOF_X1, &
            iZ_B0(1)  :iZ_E0(1)  , &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            nSpecies, &
            iZ_B0(2)  :iZ_E0(2)+1, &
            nCM)
    REAL(DP) :: &
      uCM_R(nDOF_X1, &
            iZ_B0(1)  :iZ_E0(1)  , &
            iZ_B0(3)  :iZ_E0(3)  , &
            iZ_B0(4)  :iZ_E0(4)  , &
            nSpecies, &
            iZ_B0(2)  :iZ_E0(2)+1, &
            nCM)

    INTEGER  :: iNodeZ, iNodeE, iNodeX, iNodeZ_X1, iNodeX_X1
    INTEGER  :: iNodeZ_X, nNodesZ_X, iX
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCM, iS, iGF, iPF
    INTEGER  :: iX_F, iZ_F, iX_K, iZ_K
    INTEGER  :: iZP_B0(4), iZP_E0(4)

    REAL(DP) :: Flux_L(nCM), uCM_X1_L(nCM)
    REAL(DP) :: Flux_R(nCM), uCM_X1_R(nCM)
    REAL(DP) :: Flux_K(nCM)
    REAL(DP) :: V_u_1_L, V_u_2_L, V_u_3_L
    REAL(DP) :: V_u_1_R, V_u_2_R, V_u_3_R

    print *, "ComputeIncrement_Divergence_X1"

    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN

    ! --- Permuted Phase Space Limits ---

    iZP_B0(1) = iZ_B0(1) ; iZP_E0(1) = iZ_E0(1)
    iZP_B0(2) = iZ_B0(3) ; iZP_E0(2) = iZ_E0(3)
    iZP_B0(3) = iZ_B0(4) ; iZP_E0(3) = iZ_E0(4)
    iZP_B0(4) = iZ_B0(2) ; iZP_E0(4) = iZ_E0(2)

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

    CALL InitializeIncrement_Divergence_X &
           ( iZP_B0, iZP_E0, nDOFX_X1, nDOF_X1, &
             uGF_K, uGF_F, uPF_K, uPF_L, uPF_R, uCM_K, uCM_L, uCM_R )

    ! --- Permute Geometry Fields ---
    
    DO iGF = 1, nGF
    DO iZ2 = iZ_B0(2)-1, iZ_E0(2)+1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        uGF_K(iNodeX,iZ3,iZ4,iZ2,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

    DO iGF = iGF_h_1, iGF_h_3

      ! --- Face States (Average of Left and Right States) ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               uGF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF), nDOFX, Zero, &
               uGF_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               uGF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX, Half, &
               uGF_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

    END DO

    ! --- Recompute Geometry from Scale Factors ---

    DO iZ2  = iZ_B0(2), iZ_E0(2)+1
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
          = MAX( uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_1)**2, SqrtTiny )
        uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
          = MAX( uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_2)**2, SqrtTiny )
        uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) &
          = MAX( uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_3)**2, SqrtTiny )
        uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_SqrtGm) &
          = SQRT(   uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
                  * uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
                  * uGF_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ! --- Permute Fluid Fields

    DO iPF = 1, nPF
    DO iZ2 = iZ_B0(2)-1, iZ_E0(2)+1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        uPF_K(iNodeX,iZ3,iZ4,iZ2,iPF) = U_F(iNodeX,iZ2,iZ3,iZ4,iPF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iPF = iPF_V1, iPF_V3

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               uPF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iPF), nDOFX, Zero, &
               uPF_L(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iPF), nDOFX_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               uPF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iPF), nDOFX, Zero, &
               uPF_R(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iPF), nDOFX_X1 )

    END DO

    ! --- Compute Face Velocity Components ---

    DO iX_F = 1, nNodesX_X1

      ! --- Left State --- Is this correct?

      V_u_1_L = uV1_L(iX_F)
      V_u_2_L = uV2_L(iX_F)
      V_u_3_L = uV3_L(iX_F)

      ! --- Right State --- Is this correct?

      V_u_1_R = uV1_R(iX_F)
      V_u_2_R = uV2_R(iX_F)
      V_u_3_R = uV3_R(iX_F)

      CALL FaceVelocity_X1 &
             ( V_u_1_L, V_u_2_L, V_u_3_L, &
               V_u_1_R, V_u_2_R, V_u_3_R, &
               V_u_1_F(iX_F), V_u_2_F(iX_F), V_u_3_F(iX_F) )
      
      ! print*,"iX_F    = ", iX_F
      ! print*,"V_u_1_F = ", V_u_1_F(iX_F)
      ! print*,"V_u_2_F = ", V_u_2_F(iX_F)
      ! print*,"V_u_3_F = ", V_u_3_F(iX_F)

    END DO

    ! --- Permute Two-Moment Fields ---
    
    DO iCM = 1, nCM
    DO iZ2 = iZ_B0(2)-1, iZ_E0(2)+1
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        uCM_K(iNodeZ,iZ1,iZ3,iZ4,iS,iZ2,iCM) &
          = U_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Two-Moment Fields ---

    DO iCM = 1, nCM

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_X1, nX1_Z, nDOFZ, One, L_X1_Up, nDOF_X1, &
               uCM_K(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)-1,iCM), nDOFZ, Zero, &
               uCM_L(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ,iCM), nDOF_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_X1, nX1_Z, nDOFZ, One, L_X1_Dn, nDOF_X1, &
               uCM_K(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ,iCM), nDOFZ, Zero, &
               uCM_R(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ,iCM), nDOF_X1 )

    END DO

    ! print *, "uCM_L = ", uCM_L(:,:,:,:,:,:,1)
    ! print *, "uCM_R = ", uCM_R(:,:,:,:,:,:,1)

    ! --- Numerical Flux ---

    ! nNodesZ_X = SIZE( E_L, 1 )
    ! print *, "nNodesZ_X  = ", nNodesZ_X
    ! print *, "nNodesZ_X1 = ", nNodesZ_X1


    DO iZ_F = 1, nNodesZ_X1

      iX = PositionIndexZ_F(iZ_F)

      ! --- Left State Primitive ---

      CALL ComputePrimitive_TwoMoment_Richardson_FMC &
              ( E_L(iZ_F), &
                F_d_1_L(iZ_F), &
                F_d_2_L(iZ_F), &
                F_d_3_L(iZ_F), &
                J_L(iZ_F), &
                H_u_1_L(iZ_F), &
                H_u_2_L(iZ_F), &
                H_u_3_L(iZ_F), &
                V_u_1_F(iX), &
                V_u_2_F(iX), &
                V_u_3_F(iX), &
                Gm_dd_11_F(iX), &
                Gm_dd_22_F(iX), &
                Gm_dd_33_F(iX) )

      ! --- Right State Primitive ---

      CALL ComputePrimitive_TwoMoment_Richardson_FMC &
              ( E_R(iZ_F), &
                F_d_1_R(iZ_F), &
                F_d_2_R(iZ_F), &
                F_d_3_R(iZ_F), &
                J_R(iZ_F), &
                H_u_1_R(iZ_F), &
                H_u_2_R(iZ_F), &
                H_u_3_R(iZ_F), &
                V_u_1_F(iX), &
                V_u_2_F(iX), &
                V_u_3_F(iX), &
                Gm_dd_11_F(iX), &
                Gm_dd_22_F(iX), &
                Gm_dd_33_F(iX) )

      ! print *, "iX             = ", iX, " iZ_F = ", iZ_F
      ! print *, "E_L(iZ_F)     = ", E_L(iZ_F)
      ! print *, "E_R(iZ_F)     = ", E_R(iZ_F)
      ! print *, "F_d_1_R(iNodeZ_X) = ", F_d_1_R(iNodeZ_X)
      ! print *, "F_d_2_R(iNodeZ_X) = ", F_d_2_R(iNodeZ_X)
      ! print *, "F_d_3_R(iNodeZ_X) = ", F_d_3_R(iNodeZ_X)
      ! print *, "J_L(iZ_F)     = ", J_L(iZ_F)
      ! print *, "J_R(iZ_F)     = ", J_R(iZ_F)
      ! print *, "H_u_1_R(iNodeZ_X) = ", H_u_1_R(iNodeZ_X)
      ! print *, "h_u_2_R(iNodeZ_X) = ", H_u_2_R(iNodeZ_X)
      ! print *, "H_u_3_R(iNodeZ_X) = ", H_u_3_R(iNodeZ_X)
      ! print *, "V_u_1_F(iX)       = ", V_u_1_F(iX)
      ! print *, "V_u_3_F(iX)       = ", V_u_2_F(iX)
      ! print *, "V_u_2_F(iX)       = ", V_u_3_F(iX)
      ! print *, "Gm_dd_11_F(iX)    = ", Gm_dd_11_F(iX)
      ! print *, "Gm_dd_22_F(iX)    = ", Gm_dd_22_F(iX)
      ! print *, "Gm_dd_33_F(iX)    = ", Gm_dd_33_F(iX)

    END DO

    DO iZ_F = 1, nNodesZ_X1

      iX_F = PositionIndexZ_F(iZ_F)

      iNodeE    = IndexTableZ_F(1,iZ_F)
      iNodeX_X1 = IndexTableZ_F(2,iZ_F)
      iZ1       = IndexTableZ_F(3,iZ_F)
      iZ3       = IndexTableZ_F(4,iZ_F)
      iZ4       = IndexTableZ_F(5,iZ_F)
      iS        = IndexTableZ_F(6,iZ_F)
      iZ2       = IndexTableZ_F(7,iZ_F)

      iNodeZ_X1 = iNodeE + ( iNodeX_X1 - 1 ) * nDOFE

      ! --- Left State Flux ---

      Flux_L &
        = Flux_X1 &
            ( J_L (iZ_F), H_u_1_L(iZ_F), &
              H_u_2_L(iZ_F), H_u_3_L(iZ_F), &
              V_u_1_F(iX_F), V_u_2_F(iX_F), V_u_3_F(iX_F), &
              Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      CALL ComputeConserved_TwoMoment_FMC &
             ( J_L (iZ_F), H_u_1_L(iZ_F), &
               H_u_2_L(iZ_F), H_u_3_L(iZ_F), &
               uCM_X1_L(iCM_E), uCM_X1_L(iCM_F1), &
               uCM_X1_L(iCM_F2), uCM_X1_L(iCM_F3), &
               V_u_1_F(iX_F), V_u_2_F(iX_F), V_u_3_F(iX_F), &
               Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      ! --- Right State Flux ---

      Flux_R &
        = Flux_X1 &
            ( J_R (iZ_F), H_u_1_R(iZ_F), &
              H_u_2_R(iZ_F), H_u_3_R(iZ_F), &
              V_u_1_F(iX_F), V_u_2_F(iX_F), V_u_3_F(iX_F), &
              Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      CALL ComputeConserved_TwoMoment_FMC &
             ( J_R (iZ_F), H_u_1_R(iZ_F), &
               H_u_2_R(iZ_F), H_u_3_R(iZ_F), &
               uCM_X1_R(iCM_E), uCM_X1_R(iCM_F1), &
               uCM_X1_R(iCM_F2), uCM_X1_R(iCM_F3), &
               V_u_1_F(iX_F), V_u_2_F(iX_F), V_u_3_F(iX_F), &
               Gm_dd_11_F(iX_F), Gm_dd_22_F(iX_F), Gm_dd_33_F(iX_F) )

      ! --- Numerical Flux ---
      ! print*,"iZ_F        = ", iZ_F, "iX_F = ", iX_F
      ! print*,"uCM_X1_L    = ", E_L(iZ_F), F_d_1_L(iZ_F), F_d_2_L(iZ_F), F_d_3_L(iZ_F)
      ! print*,"uCM_X1_R    = ", E_R(iZ_F), F_d_1_R(iZ_F), F_d_2_R(iZ_F), F_d_3_R(iZ_F)
      ! print*,"Flux_L      = ", Flux_L
      ! print*,"Flux_R      = ", Flux_R
      ! print*,"SqrtGm      = ", SqrtGm_F(iX_F)

      DO iCM = 1, nCM

        NumericalFlux(iNodeZ_X1,iCM,iZ1,iZ3,iZ4,iS,iZ2) &
          = NumericalFlux_LLF &
              ( uCM_X1_L(iCM), uCM_X1_R(iCM), Flux_L(iCM), Flux_R(iCM), One )

        NumericalFlux(iNodeZ_X1,iCM,iZ1,iZ3,iZ4,iS,iZ2) &
          = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_X1(iNodeZ_X1) * GE(iNodeE,iZ1,iGE_Ep2) &
              * SqrtGm_F(iX_F) &
              * NumericalFlux(iNodeZ_X1,iCM,iZ1,iZ3,iZ4,iS,iZ2)

      END DO

      ! NumericalFlux(iNodeZ_X1,iCM_E,iZ1,iZ3,iZ4,iS,iZ2) &
      !     = NumericalFlux_LLF &
      !         ( E_L(iZ_F), E_R(iZ_F), Flux_L(iCM_E), Flux_R(iCM_E), One )

      ! NumericalFlux(iNodeZ_X1,iCM_E,iZ1,iZ3,iZ4,iS,iZ2) &
      !     = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
      !         * Weights_X1(iNodeZ_X1) * GE(iNodeE,iZ1,iGE_Ep2) &
      !         * SqrtGm_F(iX_F) &
      !         * NumericalFlux(iNodeZ_X1,iCM_E,iZ1,iZ3,iZ4,iS,iZ2)

      ! NumericalFlux(iNodeZ_X1,iCM_F1,iZ1,iZ3,iZ4,iS,iZ2) &
      !     = NumericalFlux_LLF &
      !         ( F_d_1_L(iZ_F), F_d_1_R(iZ_F), Flux_L(iCM_F1), Flux_R(iCM_F1), One )

      ! NumericalFlux(iNodeZ_X1,iCM_F1,iZ1,iZ3,iZ4,iS,iZ2) &
      !     = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
      !         * Weights_X1(iNodeZ_X1) * GE(iNodeE,iZ1,iGE_Ep2) &
      !         * SqrtGm_F(iX_F) &
      !         * NumericalFlux(iNodeZ_X1,iCM_F1,iZ1,iZ3,iZ4,iS,iZ2)

      ! NumericalFlux(iNodeZ_X1,iCM_F2,iZ1,iZ3,iZ4,iS,iZ2) &
      !     = NumericalFlux_LLF &
      !         ( F_d_2_L(iZ_F), F_d_2_R(iZ_F), Flux_L(iCM_F2), Flux_R(iCM_F2), One )

      ! NumericalFlux(iNodeZ_X1,iCM_F2,iZ1,iZ3,iZ4,iS,iZ2) &
      !     = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
      !         * Weights_X1(iNodeZ_X1) * GE(iNodeE,iZ1,iGE_Ep2) &
      !         * SqrtGm_F(iX_F) &
      !         * NumericalFlux(iNodeZ_X1,iCM_F2,iZ1,iZ3,iZ4,iS,iZ2)

      ! NumericalFlux(iNodeZ_X1,iCM_F3,iZ1,iZ3,iZ4,iS,iZ2) &
      !     = NumericalFlux_LLF &
      !         ( F_d_3_L(iZ_F), F_d_3_R(iZ_F), Flux_L(iCM_F3), Flux_R(iCM_F3), One )

      ! NumericalFlux(iNodeZ_X1,iCM_F3,iZ1,iZ3,iZ4,iS,iZ2) &
      !     = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
      !         * Weights_X1(iNodeZ_X1) * GE(iNodeE,iZ1,iGE_Ep2) &
      !         * SqrtGm_F(iX_F) &
      !         * NumericalFlux(iNodeZ_X1,iCM_F3,iZ1,iZ3,iZ4,iS,iZ2)

    END DO

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---
    ! print *, "dU_Z(:,iCM_E,:,:,:,:,:)  = ", dU_Z(:,iCM_E,:,:,:,:,:)
    
    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCM, nDOF_X1, + One, L_X1_Dn, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), &
             nDOF_X1, Zero, dU_Z, nDOFZ )

    ! print *, "dU_Z(:,iCM_E,:,:,:,:,:)  = ", dU_Z(:,iCM_E,:,:,:,:,:)

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCM, nDOF_X1, - One, L_X1_Up, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)+1), &
             nDOF_X1, One,  dU_Z, nDOFZ )

    ! print *, "dU_Z(:,iCM_E,:,:,:,:,:)  = ", dU_Z(:,iCM_E,:,:,:,:,:)

    !--------------------
    ! --- Volume Term ---
    !--------------------

    ! --- Compute Primitive Fluid in Spatial Elements ---

    DO iX_K = 1, nNodesX_K

      V_u_1_K(iX_K) = uV1_K(iX_K)
      V_u_2_K(iX_K) = uV2_K(iX_K)
      V_u_3_K(iX_K) = uV3_K(iX_K)

    END DO

    DO iZ_K = 1, nNodesZ_K

      iX_K = PositionIndexZ_K(iZ_K)

      CALL ComputePrimitive_TwoMoment_Richardson_FMC &
            ( E_K(iZ_K), &
              F_d_1_K(iZ_K), &
              F_d_2_K(iZ_K), &
              F_d_3_K(iZ_K), &
              J_K(iZ_K), &
              H_u_1_K(iZ_K), &
              H_u_2_K(iZ_K), &
              H_u_3_K(iZ_K), &
              V_u_1_K(iX_K), &
              V_u_2_K(iX_K), &
              V_u_3_K(iX_K), &
              Gm_dd_11_K(iX_K), &
              Gm_dd_22_K(iX_K), &
              Gm_dd_33_K(iX_K) )

    END DO
    
    DO iZ_K = 1, nNodesZ_K

      iX_K = PositionIndexZ_K(iZ_K)

      iNodeE = IndexTableZ_K(1,iZ_K)
      iNodeX = IndexTableZ_K(2,iZ_K)
      iZ1    = IndexTableZ_K(3,iZ_K)
      iZ3    = IndexTableZ_K(4,iZ_K)
      iZ4    = IndexTableZ_K(5,iZ_K)
      iS     = IndexTableZ_K(6,iZ_K)
      iZ2    = IndexTableZ_K(7,iZ_K)

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      Flux_K &
        = Flux_X1 &
            ( J_K (iZ_K), H_u_1_K(iZ_K), &
              H_u_2_K(iZ_K), H_u_3_K(iZ_K), &
              V_u_1_K(iX_K), V_u_2_K(iX_K), V_u_3_K(iX_K), &
              Gm_dd_11_K(iX_K), Gm_dd_22_K(iX_K), Gm_dd_33_K(iX_K) )

      DO iCM = 1, nCM

        Flux_q(iNodeZ,iCM,iZ1,iZ3,iZ4,iS,iZ2) &
          = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * SqrtGm_K(iX_K) &
              * Flux_K(iCM)

      END DO

    END DO

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCM, nDOFZ, One, dLdX1_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_Z, nDOFZ )

    ! print *, "dU_Z(:,iCM_E,:,:,:,:,:)  = ", dU_Z(:,iCM_E,:,:,:,:,:)
    ! print *, "dU_Z(:,iCM_F1,:,:,:,:,:) = ", dU_Z(:,iCM_F1,:,:,:,:,:)
    ! print *, "dE - dF1                 = ", dU_Z(:,iCM_E,:,:,:,:,:) - dU_Z(:,iCM_F1,:,:,:,:,:)

    DO iS  = 1, nSpecies
    DO iCM = 1, nCM
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) &
          = dU_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) &
            + dU_Z(iNodeZ,iCM,iZ1,iZ3,iZ4,iS,iZ2)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence_X

    END ASSOCIATE ! dZ1, etc

  END SUBROUTINE ComputeIncrement_Divergence_X1

  SUBROUTINE ComputeIncrement_ObserverCorrections &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_M, dU_M)

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nPF)
    REAL(DP), INTENT(in)    :: &
      U_M (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCM, &
           1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_M(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCM, &
           1:nSpecies)

    INTEGER  :: iNodeZ, iNodeE, iNodeX, iNodeZ_E
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCM, iS, iGF, iPF
    INTEGER  :: iX, iX_F, iZ_F, iX_K, iZ_K
    INTEGER  :: iZP_B0(4), iZP_E0(4)

    REAL(DP) :: EdgeEnergyCubed, Beta
    REAL(DP) :: A, B_i(1:3), C_ij(1:3,1:3), Lambda(3), vMag, vMagSq, W
    REAL(DP) :: Flux_K(nCM), dFlux_K(nPM)
    REAL(DP) :: Flux_L(nCM), uPM_L(nPM)
    REAL(DP) :: Flux_R(nCM), uPM_R(nPM)
    REAL(DP) :: H_u_0_L, H_u_0_R
    REAL(DP) :: SUM_N, SUM_G1, SUM_G2, SUM_G3

    ! --- Permuted Phase Space Limits ---

    iZP_B0(1) = iZ_B0(2) ; iZP_E0(1) = iZ_E0(2)
    iZP_B0(2) = iZ_B0(3) ; iZP_E0(2) = iZ_E0(3)
    iZP_B0(3) = iZ_B0(4) ; iZP_E0(3) = iZ_E0(4)
    iZP_B0(4) = iZ_B0(1) ; iZP_E0(4) = iZ_E0(1)

    ASSOCIATE &
      ( xZ1 => MeshE    % Center, &
        dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

    CALL InitializeIncrement_ObserverCorrections( iZ_B0, iZ_E0 )

    ! --- Calculate Weak Derivatives ---

    CALL ComputeWeakDerivatives_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
             dV_d_dX1 )

    ! dV_u_dX1 = Zero
    ! dV_d_dX1 = Zero
    dV_u_dX2 = Zero
    dV_d_dX2 = Zero
    dV_u_dX3 = Zero
    dV_d_dX3 = Zero

    ! --- Permute Geometry Fields ---

    DO iGF = 1, nGF
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        uGF_K(iNodeX,iZ2,iZ3,iZ4,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Permute Fluid Fields ---

    DO iPF = 1, nPF
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        uPF_K(iNodeX,iZ2,iZ3,iZ4,iPF) = U_F(iNodeX,iZ2,iZ3,iZ4,iPF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Permute Two-Moment Fields ---

    DO iCM = 1, nCM
    DO iZ1 = iZ_B0(1)-1, iZ_E0(1)+1
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeZ = 1, nDOFZ

        uCM_K(iNodeZ,iZ2,iZ3,iZ4,iS,iZ1,iCM) &
          = U_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Two-Moment Fields ---

    DO iCM = 1, nCM

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_E, nE_Z, nDOFZ, One, L_E_Up, nDOF_E, &
               uCM_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)-1,iCM), nDOFZ, Zero, &
               uCM_L(1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ,iCM), nDOF_E )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOF_E, nE_Z, nDOFZ, One, L_E_Dn, nDOF_E, &
               uCM_K(1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ,iCM), nDOFZ, Zero, &
               uCM_R(1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ,iCM), nDOF_E )

    END DO

    ! --- Eigenvalues ---

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        ! --- Quadratic Form Matrix ---
        ! Assuming velocity is constant with respect to time

        A = Zero

        B_i = Half * [ dV_d_dX1(iNodeX,0,iZ2,iZ3,iZ4), &
                       dV_d_dX2(iNodeX,0,iZ2,iZ3,iZ4), &
                       dV_d_dX3(iNodeX,0,iZ2,iZ3,iZ4) ]

        C_ij(:,1) = Half * [ Two * dV_d_dX1(iNodeX,1,iZ2,iZ3,iZ4), &
                                   dV_d_dX2(iNodeX,1,iZ2,iZ3,iZ4)  &
                                 + dV_d_dX1(iNodeX,2,iZ2,iZ3,iZ4), &
                                   dV_d_dX3(iNodeX,1,iZ2,iZ3,iZ4)  &
                                 + dV_d_dX1(iNodeX,3,iZ2,iZ3,iZ4) ]

        C_ij(:,2) = Half * [       dV_d_dX1(iNodeX,2,iZ2,iZ3,iZ4)  &
                                 + dV_d_dX2(iNodeX,1,iZ2,iZ3,iZ4), &
                             Two * dV_d_dX2(iNodeX,2,iZ2,iZ3,iZ4), &
                                   dV_d_dX3(iNodeX,2,iZ2,iZ3,iZ4)  &
                                 + dV_d_dX2(iNodeX,3,iZ2,iZ3,iZ4) ]

        C_ij(:,3) = Half * [       dV_d_dX1(iNodeX,3,iZ2,iZ3,iZ4)  &
                                 + dV_d_dX3(iNodeX,1,iZ2,iZ3,iZ4), &
                                   dV_d_dX2(iNodeX,3,iZ2,iZ3,iZ4)  &
                                 + dV_d_dX3(iNodeX,2,iZ2,iZ3,iZ4), &
                             Two * dV_d_dX3(iNodeX,3,iZ2,iZ3,iZ4) ]

        ! A(:,0) = Half * [ Zero, &
        !                   dV_d_dX1(iNodeX,0,iZ2,iZ3,iZ4), &
        !                   dV_d_dX2(iNodeX,0,iZ2,iZ3,iZ4), &
        !                   dV_d_dX3(iNodeX,0,iZ2,iZ3,iZ4) ]

        ! A(:,1) = Half * [       dV_d_dX1(iNodeX,0,iZ2,iZ3,iZ4), &
        !                   Two * dV_d_dX1(iNodeX,1,iZ2,iZ3,iZ4), &
        !                         dV_d_dX2(iNodeX,1,iZ2,iZ3,iZ4)  &
        !                       + dV_d_dX1(iNodeX,2,iZ2,iZ3,iZ4), &
        !                         dV_d_dX3(iNodeX,1,iZ2,iZ3,iZ4)  &
        !                       + dV_d_dX1(iNodeX,3,iZ2,iZ3,iZ4) ]
        ! A(:,2) = Half * [       dV_d_dX2(iNodeX,0,iZ2,iZ3,iZ4), &
        !                         dV_d_dX1(iNodeX,2,iZ2,iZ3,iZ4)  &
        !                       + dV_d_dX2(iNodeX,1,iZ2,iZ3,iZ4), &
        !                   Two * dV_d_dX2(iNodeX,2,iZ2,iZ3,iZ4), &
        !                         dV_d_dX3(iNodeX,2,iZ2,iZ3,iZ4)  &
        !                       + dV_d_dX2(iNodeX,3,iZ2,iZ3,iZ4) ]
        ! A(:,3) = Half * [       dV_d_dX3(iNodeX,0,iZ2,iZ3,iZ4), &
        !                         dV_d_dX1(iNodeX,3,iZ2,iZ3,iZ4)  &
        !                       + dV_d_dX3(iNodeX,1,iZ2,iZ3,iZ4), &
        !                         dV_d_dX2(iNodeX,3,iZ2,iZ3,iZ4)  &
        !                       + dV_d_dX3(iNodeX,2,iZ2,iZ3,iZ4), &
        !                   Two * dV_d_dX3(iNodeX,3,iZ2,iZ3,iZ4) ]

        

        CALL EigenvaluesSymmetric3( C_ij, Lambda )

        vMag = SQRT( uPF_K(iNodeX,iZ2,iZ3,iZ4,iPF_V1)**2 &
                     + uPF_K(iNodeX,iZ2,iZ3,iZ4,iPF_V2)**2 &
                     + uPF_K(iNodeX,iZ2,iZ3,iZ4,iPF_V3)**2 )

        W = One / SQRT( One - vMag**2 )

        Alpha(iNodeX,iZ2,iZ3,iZ4) = MAXVAL( ABS( Lambda ) ) &
                                    + Two * SQRT ( B_i(1)**2 &
                                                 + B_i(2)**2 &
                                                 + B_i(3)**2 ) &
                                    + ABS( A )
        Alpha(iNodeX,iZ2,iZ3,iZ4) = Alpha(iNodeX,iZ2,iZ3,iZ4) &
                                    * W**2 * ( One + vMag )**2

      END DO

    END DO
    END DO
    END DO

    DO iZ_F = 1, nNodesZ_E

      iX = PositionIndexZ_F(iZ_F)

      ! --- Left State Primitive ---

      CALL ComputePrimitive_TwoMoment_Richardson_FMC &
              ( E_L(iZ_F), &
                F_d_1_L(iZ_F), &
                F_d_2_L(iZ_F), &
                F_d_3_L(iZ_F), &
                J_L(iZ_F), &
                H_u_1_L(iZ_F), &
                H_u_2_L(iZ_F), &
                H_u_3_L(iZ_F), &
                uV1_K(iX), &
                uV2_K(iX), &
                uV3_K(iX), &
                Gm_dd_11_K(iX), &
                Gm_dd_22_K(iX), &
                Gm_dd_33_K(iX) )

      ! --- Right State Primitive ---

      CALL ComputePrimitive_TwoMoment_Richardson_FMC &
              ( E_R(iZ_F), &
                F_d_1_R(iZ_F), &
                F_d_2_R(iZ_F), &
                F_d_3_R(iZ_F), &
                J_R(iZ_F), &
                H_u_1_R(iZ_F), &
                H_u_2_R(iZ_F), &
                H_u_3_R(iZ_F), &
                uV1_K(iX), &
                uV2_K(iX), &
                uV3_K(iX), &
                Gm_dd_11_K(iX), &
                Gm_dd_22_K(iX), &
                Gm_dd_33_K(iX) )

    END DO

    ! --- Numerical Flux ---

    DO iZ_F = 1, nNodesZ_E

      iX_F = PositionIndexZ_F(iZ_F)

      iNodeZ_E  = IndexTableZ_F(2,iZ_F) ! = iNodeX
      iZ2       = IndexTableZ_F(3,iZ_F)
      iZ3       = IndexTableZ_F(4,iZ_F)
      iZ4       = IndexTableZ_F(5,iZ_F)
      iS        = IndexTableZ_F(6,iZ_F)
      iZ1       = IndexTableZ_F(7,iZ_F)

      ! --- Left State Flux ---

      Flux_L &
        = Flux_E( J_L (iZ_F), H_u_1_L(iZ_F), H_u_2_L(iZ_F), H_u_3_L(iZ_F), &
                  uV1_K(iX_F), uV2_K(iX_F), uV3_K(iX_F), &
                  dV_d_dX1(iNodeZ_E,0,iZ2,iZ3,iZ4), &
                  dV_d_dX1(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_d_dX1(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_d_dX1(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeZ_E,0,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeZ_E,0,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  Gm_dd_11_K(iX_F), Gm_dd_22_K(iX_F), Gm_dd_33_K(iX_F) )

      vMagSq = uV1_K(iX_F)**2 + uV2_K(iX_F)**2 + uV3_K(iX_F)**2

      W = One / SQRT( One - vMagSq )

      H_u_0_L = uV1_K(iX_F) * H_u_1_L(iZ_F) &
                + uV2_K(iX_F) * H_u_2_L(iZ_F) &
                + uV3_K(iX_F) * H_u_3_L(iZ_F)

      uPM_L = [ W * J_L(iZ_F) + H_u_0_L, &
                W * J_L(iZ_F) * uV1_K(iX_F) + H_u_1_L(iZ_F), &
                W * J_L(iZ_F) * uV2_K(iX_F) + H_u_2_L(iZ_F), &
                W * J_L(iZ_F) * uV3_K(iX_F) + H_u_3_L(iZ_F) ]

      ! --- Right State Flux ---

      Flux_R &
        = Flux_E( J_R (iZ_F), H_u_1_R(iZ_F), H_u_2_R(iZ_F), H_u_3_R(iZ_F), &
                  uV1_K(iX_F), uV2_K(iX_F), uV3_K(iX_F), &
                  dV_d_dX1(iNodeZ_E,0,iZ2,iZ3,iZ4), &
                  dV_d_dX1(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_d_dX1(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_d_dX1(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeZ_E,0,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeZ_E,0,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeZ_E,1,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeZ_E,2,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeZ_E,3,iZ2,iZ3,iZ4), &
                  Gm_dd_11_K(iX_F), Gm_dd_22_K(iX_F), Gm_dd_33_K(iX_F) )

      H_u_0_R = uV1_K(iX_F) * H_u_1_R(iZ_F) &
                + uV2_K(iX_F) * H_u_2_R(iZ_F) &
                + uV3_K(iX_F) * H_u_3_R(iZ_F)

      uPM_R = [ W * J_R(iZ_F) + H_u_0_R, &
                W * J_R(iZ_F) * uV1_K(iX_F) + H_u_1_R(iZ_F), &
                W * J_R(iZ_F) * uV2_K(iX_F) + H_u_2_R(iZ_F), &
                W * J_R(iZ_F) * uV3_K(iX_F) + H_u_3_R(iZ_F) ]

      ! --- Numerical Flux (Local Lax-Friedrichs) ---

      EdgeEnergyCubed = ( xZ1(iZ1) - Half * dZ1(iZ1) )**3

      DO iCM = 1, nCM

        NumericalFlux(iNodeZ_E,iCM,iZ2,iZ3,iZ4,iS,iZ1) &
          = NumericalFlux_LLF &
              ( uPM_L(iCM), uPM_R(iCM), Flux_L(iCM), Flux_R(iCM), &
                Alpha(iNodeZ_E,iZ2,iZ3,iZ4) )

        NumericalFlux(iNodeZ_E,iCM,iZ2,iZ3,iZ4,iS,iZ1) &
          = dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_E(iNodeZ_E) * EdgeEnergyCubed &
              * SqrtGm_K(iX_F) &
              * NumericalFlux(iNodeZ_E,iCM,iZ2,iZ3,iZ4,iS,iZ1)

      END DO

    END DO

    ! --- Surface Contributions ---

    ! --- Contributions from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCM, nDOF_E, + One, L_E_Dn, nDOF_E, &
             NumericalFlux(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ), &
             nDOF_E, Zero, dU_Z, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCM, nDOF_E, - One, L_E_Up, nDOF_E, &
             NumericalFlux(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)+1), &
             nDOF_E, One,  dU_Z, nDOFZ )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    DO iZ_K = 1, nNodesZ_K

      iX_K = PositionIndexZ_K(iZ_K)

      CALL ComputePrimitive_TwoMoment_Richardson_FMC &
            ( E_K(iZ_K), &
              F_d_1_K(iZ_K), &
              F_d_2_K(iZ_K), &
              F_d_3_K(iZ_K), &
              J_K(iZ_K), &
              H_u_1_K(iZ_K), &
              H_u_2_K(iZ_K), &
              H_u_3_K(iZ_K), &
              uV1_K(iX_K), &
              uV2_K(iX_K), &
              uV3_K(iX_K), &
              Gm_dd_11_K(iX_K), &
              Gm_dd_22_K(iX_K), &
              Gm_dd_33_K(iX_K) )

    END DO

    DO iZ_K = 1, nNodesZ_K

      iX_K = PositionIndexZ_K(iZ_K)

      iNodeE = IndexTableZ_K(1,iZ_K)
      iNodeX = IndexTableZ_K(2,iZ_K)
      iZ2    = IndexTableZ_K(3,iZ_K)
      iZ3    = IndexTableZ_K(4,iZ_K)
      iZ4    = IndexTableZ_K(5,iZ_K)
      iS     = IndexTableZ_K(6,iZ_K)
      iZ1    = IndexTableZ_K(7,iZ_K)

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      Flux_K &
        = Flux_E( J_K(iZ_K), H_u_1_K(iZ_K), H_u_2_K(iZ_K), H_u_3_K(iZ_K), &
                  uV1_K(iX_K), uV2_K(iX_K), uV3_K(iX_K), &
                  dV_d_dX1(iNodeX,0,iZ2,iZ3,iZ4), &
                  dV_d_dX1(iNodeX,1,iZ2,iZ3,iZ4), &
                  dV_d_dX1(iNodeX,2,iZ2,iZ3,iZ4), &
                  dV_d_dX1(iNodeX,3,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeX,0,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeX,1,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeX,2,iZ2,iZ3,iZ4), &
                  dV_d_dX2(iNodeX,3,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeX,0,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeX,1,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeX,2,iZ2,iZ3,iZ4), &
                  dV_d_dX3(iNodeX,3,iZ2,iZ3,iZ4), &
                  Gm_dd_11_K(iX_K), Gm_dd_22_K(iX_K), Gm_dd_33_K(iX_K) )

      DO iCM = 1, nCM

        Flux_q(iNodeZ,iCM,iZ2,iZ3,iZ4,iS,iZ1) &
          = dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4)  &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep3) &
              * SqrtGm_K(iX_K) &
              * Flux_K(iCM)

      END DO

    END DO

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nK_Z*nCM, nDOFZ, One, dLdE_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_Z, nDOFZ )

    DO iS  = 1, nSpecies
    DO iCM = 1, nCM
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) &
          = dU_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) &
            + dU_Z(iNodeZ,iCM,iZ2,iZ3,iZ4,iS,iZ1)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_ObserverCorrections

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_ObserverCorrections

  SUBROUTINE InitializeIncrement_Divergence_X &
    ( iZP_B0, iZP_E0, nDOFX_X, nDOFZ_X, &
      uGF_K, uGF_F, uPF_K, uPF_L, uPF_R, uCM_K, uCM_L, uCM_R )

    INTEGER, INTENT(in) :: iZP_B0(4), iZP_E0(4) ! Permuted limits
    INTEGER, INTENT(in) :: nDOFX_X ! nDOFX_X1, ...
    INTEGER, INTENT(in) :: nDOFZ_X ! nDOFZ_X1, ...

    ! --- Geometry Fields ---
    REAL(DP), TARGET, INTENT(in) :: &
      uGF_K (1:nDOFX, &
             iZP_B0(2)  :iZP_E0(2)  , &
             iZP_B0(3)  :iZP_E0(3)  , &
             iZP_B0(4)-1:iZP_E0(4)+1, &
            1:nGF), &
      uGF_F (1:nDOFX_X, &
             iZP_B0(2)  :iZP_E0(2)  , &
             iZP_B0(3)  :iZP_E0(3)  , &
             iZP_B0(4)  :iZP_E0(4)+1, &
             1:nGF)

    ! --- Primitive Fluid Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      uPF_K(1:nDOFX, &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            iZP_B0(4)-1:iZP_E0(4)+1, &
            1:nPF), &
      uPF_L(1:nDOFX_X, &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            iZP_B0(4)  :iZP_E0(4)+1, &
            1:nPF), &
      uPF_R(1:nDOFX_X, &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            iZP_B0(4)  :iZP_E0(4)+1, &
            1:nPF)

    ! --- Conserved Two-Moment Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      uCM_K(1:nDOFZ, &
            iZP_B0(1)  :iZP_E0(1)  , &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            1:nSpecies, &
            iZP_B0(4)-1:iZP_E0(4)+1, &
            1:nCM), &
      uCM_L(1:nDOFZ_X, &
            iZP_B0(1)  :iZP_E0(1)  , &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            1:nSpecies, &
            iZP_B0(4)  :iZP_E0(4)+1, &
            1:nCM), &
      uCM_R(1:nDOFZ_X, &
            iZP_B0(1)  :iZP_E0(1)  , &
            iZP_B0(2)  :iZP_E0(2)  , &
            iZP_B0(3)  :iZP_E0(3)  , &
            1:nSpecies, &
            iZP_B0(4)  :iZP_E0(4)+1, &
            1:nCM)

    INTEGER :: iZP1, iZP2, iZP3, iZP4, iS, iNodeE
    INTEGER :: iX_K, iZ_K, iNodeX, iNodeZ
    INTEGER :: iX_F, iZ_F, iNodeX_X, iNodeZ_X
    INTEGER :: nZP(4), nZP_X(4), nX_X, nX_Z, nNodesX_X, nNodesZ_X

    nZP = iZP_E0 - iZP_B0 + 1
    nZP_X = nZP + [0,0,0,1]
    nX_X = PRODUCT( nZP_X(2:4) )
    nX_Z = nSpecies * nZP_X(1) * nX_X
    nNodesX_X = nDOFX_X * nX_X
    nNodesZ_X = nDOFZ_X * nX_Z

    ! --- Geometry Fields ---

    Gm_dd_11_K(1:nNodesX_K) => uGF_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Gm_dd_11)
    Gm_dd_22_K(1:nNodesX_K) => uGF_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Gm_dd_22)
    Gm_dd_33_K(1:nNodesX_K) => uGF_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_Gm_dd_33)
    SqrtGm_K  (1:nNodesX_K) => uGF_K(:,:,:,iZP_B0(4):iZP_E0(4),iGF_SqrtGm  )

    Gm_dd_11_F(1:nNodesX_X) => uGF_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_Gm_dd_11)
    Gm_dd_22_F(1:nNodesX_X) => uGF_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_Gm_dd_22)
    Gm_dd_33_F(1:nNodesX_X) => uGF_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_Gm_dd_33)
    SqrtGm_F  (1:nNodesX_X) => uGF_F(:,:,:,iZP_B0(4):iZP_E0(4)+1,iGF_SqrtGm  )

    ! --- Primitive Fluid Fields

    uV1_K(1:nNodesX_K) => uPF_K(:,:,:,iZP_B0(4):iZP_E0(4),iPF_V1)
    uV2_K(1:nNodesX_K) => uPF_K(:,:,:,iZP_B0(4):iZP_E0(4),iPF_V2)
    uV3_K(1:nNodesX_K) => uPF_K(:,:,:,iZP_B0(4):iZP_E0(4),iPF_V3)

    uV1_L(1:nNodesX_X) => uPF_L(:,:,:,iZP_B0(4):iZP_E0(4)+1,iPF_V1)
    uV2_L(1:nNodesX_X) => uPF_L(:,:,:,iZP_B0(4):iZP_E0(4)+1,iPF_V2)
    uV3_L(1:nNodesX_X) => uPF_L(:,:,:,iZP_B0(4):iZP_E0(4)+1,iPF_V3)

    uV1_R(1:nNodesX_X) => uPF_R(:,:,:,iZP_B0(4):iZP_E0(4)+1,iPF_V1)
    uV2_R(1:nNodesX_X) => uPF_R(:,:,:,iZP_B0(4):iZP_E0(4)+1,iPF_V2)
    uV3_R(1:nNodesX_X) => uPF_R(:,:,:,iZP_B0(4):iZP_E0(4)+1,iPF_V3)

    ! --- Conserved Two-Moment Fields

    E_K    (1:nNodesZ_K) => uCM_K(:,:,:,:,:,iZP_B0(4):iZP_E0(4),iCM_E )
    F_d_1_K(1:nNodesZ_K) => uCM_K(:,:,:,:,:,iZP_B0(4):iZP_E0(4),iCM_F1)
    F_d_2_K(1:nNodesZ_K) => uCM_K(:,:,:,:,:,iZP_B0(4):iZP_E0(4),iCM_F2)
    F_d_3_K(1:nNodesZ_K) => uCM_K(:,:,:,:,:,iZP_B0(4):iZP_E0(4),iCM_F3)

    E_L    (1:nNodesZ_X) => uCM_L(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCM_E )
    F_d_1_L(1:nNodesZ_X) => uCM_L(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCM_F1)
    F_d_2_L(1:nNodesZ_X) => uCM_L(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCM_F2)
    F_d_3_L(1:nNodesZ_X) => uCM_L(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCM_F3)

    E_R    (1:nNodesZ_X) => uCM_R(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCM_E )
    F_d_1_R(1:nNodesZ_X) => uCM_R(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCM_F1)
    F_d_2_R(1:nNodesZ_X) => uCM_R(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCM_F2)
    F_d_3_R(1:nNodesZ_X) => uCM_R(:,:,:,:,:,iZP_B0(4):iZP_E0(4)+1,iCM_F3)

    ! --- Primitive Two-Moment Fields ---

    ALLOCATE( J_L(nNodesZ_X) )
    ALLOCATE( H_u_1_L(nNodesZ_X) )
    ALLOCATE( H_u_2_L(nNodesZ_X) )
    ALLOCATE( H_u_3_L(nNodesZ_X) )

    ALLOCATE( J_R(nNodesZ_X) )
    ALLOCATE( H_u_1_R(nNodesZ_X) )
    ALLOCATE( H_u_2_R(nNodesZ_X) )
    ALLOCATE( H_u_3_R(nNodesZ_X) )

    ! --- Primitive Fluid Velocities ---

    ALLOCATE( V_u_1_F(nNodesX_X) )
    ALLOCATE( V_u_2_F(nNodesX_X) )
    ALLOCATE( V_u_3_F(nNodesX_X) )

    ! --- Fluxes ---

    ALLOCATE( NumericalFlux (nDOFZ_X,nCM, &
                             iZP_B0(1)  :iZP_E0(1)  , &
                             iZP_B0(2)  :iZP_E0(2)  , &
                             iZP_B0(3)  :iZP_E0(3)  , &
                             nSpecies, &
                             iZP_B0(4)  :iZP_E0(4)+1) )
    ALLOCATE( Flux_q (nDOFZ,nCM, &
                      iZP_B0(1)  :iZP_E0(1)  , &
                      iZP_B0(2)  :iZP_E0(2)  , &
                      iZP_B0(3)  :iZP_E0(3)  , &
                      nSpecies, &
                      iZP_B0(4)  :iZP_E0(4)) )

    ! --- Increment ---
    
    ALLOCATE( dU_Z          (nDOFZ,nCM, &
                             iZP_B0(1)  :iZP_E0(1)  , &
                             iZP_B0(2)  :iZP_E0(2)  , &
                             iZP_B0(3)  :iZP_E0(3)  , &
                             nSpecies, &
                             iZP_B0(4)  :iZP_E0(4)  ) )

    ! --- ? ---

    ALLOCATE( PositionIndexZ_F(nNodesZ_X) )
    ALLOCATE( PositionIndexZ_K(nNodesZ_K) )

    ALLOCATE( IndexTableZ_F(7,nNodesZ_X) )
    ALLOCATE( IndexTableZ_K(7,nNodesZ_K) )

    DO iZP4 = iZP_B0(4), iZP_E0(4)+1
    DO iS  = 1, nSpecies
    DO iZP3 = iZP_B0(3), iZP_E0(3)
    DO iZP2 = iZP_B0(2), iZP_E0(2)
    DO iZP1 = iZP_B0(1), iZP_E0(1)
    DO iNodeX_X = 1, nDOFX_X
    DO iNodeE = 1, nDOFE

      iNodeZ_X = iNodeE &
                 + ( iNodeX_X - 1 ) * nDOFE

      iZ_F = iNodeZ_X &
             + ( iZP1 - iZP_B0(1) ) * nDOFZ_X &
             + ( iZP2 - iZP_B0(2) ) * nDOFZ_X * nZP_X(1) &
             + ( iZP3 - iZP_B0(3) ) * nDOFZ_X * nZP_X(1) * nZP_X(2) &
             + ( iS  - 1          ) * nDOFZ_X * nZP_X(1) * nZP_X(2) * nZP_X(3) &
             + ( iZP4 - iZP_B0(4) ) * nDOFZ_X * nZP_X(1) * nZP_X(2) * nZP_X(3) * nSpecies

      iX_F = iNodeX_X &
             + ( iZP2 - iZP_B0(2) ) * nDOFX_X &
             + ( iZP3 - iZP_B0(3) ) * nDOFX_X * nZP_X(2) &
             + ( iZP4 - iZP_B0(4) ) * nDOFX_X * nZP_X(2) * nZP_X(3)

      PositionIndexZ_F(iZ_F) = iX_F

      IndexTableZ_F(:,iZ_F) = [ iNodeE, iNodeX_X, iZP1, iZP2, iZP3, iS, iZP4 ]

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    DO iZP4 = iZP_B0(4), iZP_E0(4)
    DO iS  = 1, nSpecies
    DO iZP3 = iZP_B0(3), iZP_E0(3)
    DO iZP2 = iZP_B0(2), iZP_E0(2)
    DO iZP1 = iZP_B0(1), iZP_E0(1)
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iNodeZ = iNodeE &
               + ( iNodeX - 1 ) * nDOFE

      iZ_K = iNodeZ &
             + ( iZP1 - iZP_B0(1) ) * nDOFZ &
             + ( iZP2 - iZP_B0(2) ) * nDOFZ * nZP(1) &
             + ( iZP3 - iZP_B0(3) ) * nDOFZ * nZP(1) * nZP(2) &
             + ( iS  - 1          ) * nDOFZ * nZP(1) * nZP(2) * nZP(3) &
             + ( iZP4 - iZP_B0(4) ) * nDOFZ * nZP(1) * nZP(2) * nZP(3) * nSpecies

      iX_K = iNodeX &
             + ( iZP2 - iZP_B0(2) ) * nDOFX &
             + ( iZP3 - iZP_B0(3) ) * nDOFX * nZP(2) &
             + ( iZP4 - iZP_B0(4) ) * nDOFX * nZP(2) * nZP(3)

      PositionIndexZ_K(iZ_K) = iX_K

      IndexTableZ_K(:,iZ_K) = [ iNodeE, iNodeX, iZP1, iZP2, iZP3, iS, iZP4 ]

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    RETURN
  END SUBROUTINE InitializeIncrement_Divergence_X

  SUBROUTINE FinalizeIncrement_Divergence_X

    DEALLOCATE( PositionIndexZ_F, PositionIndexZ_K )
    DEALLOCATE( IndexTableZ_F, IndexTableZ_K )
    DEALLOCATE( J_L, H_u_1_L, H_u_2_L, H_u_3_L )
    DEALLOCATE( J_R, H_u_1_R, H_u_2_R, H_u_3_R )
    DEALLOCATE( V_u_1_F, V_u_2_F, V_u_3_F )
    DEALLOCATE( NumericalFlux, Flux_q , dU_Z )
    ! DEALLOCATE( nIterations_L, nIterations_R, nIterations_K )

    NULLIFY( Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K )
    NULLIFY( Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F )
    NULLIFY( uV1_K, uV2_K, uV3_K )
    NULLIFY( uV1_L, uV2_L, uV3_L )
    NULLIFY( uV1_R, uV2_R, uV3_R )
    NULLIFY( E_K, F_d_1_K, F_d_2_K, F_d_3_K )
    NULLIFY( E_L, F_d_1_L, F_d_2_L, F_d_3_L )
    NULLIFY( E_R, F_d_1_R, F_d_2_R, F_d_3_R )

  END SUBROUTINE FinalizeIncrement_Divergence_X

  SUBROUTINE InitializeIncrement_ObserverCorrections( iZ_B0, iZ_E0 )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iNodeE
    INTEGER :: iX_K, iZ_K, iNodeX, iNodeZ
    INTEGER :: iX_F, iZ_F, iNode_E

    ! --- Geometry Fields ---

    ALLOCATE( uGF_K(1:nDOFX, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  , &
                    1:nGF) )

    Gm_dd_11_K(1:nNodesX_K) => uGF_K(:,:,:,iZ_B0(4):iZ_E0(4),iGF_Gm_dd_11)
    Gm_dd_22_K(1:nNodesX_K) => uGF_K(:,:,:,iZ_B0(4):iZ_E0(4),iGF_Gm_dd_22)
    Gm_dd_33_K(1:nNodesX_K) => uGF_K(:,:,:,iZ_B0(4):iZ_E0(4),iGF_Gm_dd_33)
    SqrtGm_K  (1:nNodesX_K) => uGF_K(:,:,:,iZ_B0(4):iZ_E0(4),iGF_SqrtGm  )

    ! --- Primitive Fluid Fields ---

    ALLOCATE( uPF_K(1:nDOFX, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  , &
                    1:nPF) )

    uV1_K(1:nNodesX_K) => uPF_K(:,:,:,iZ_B0(4):iZ_E0(4),iPF_V1)
    uV2_K(1:nNodesX_K) => uPF_K(:,:,:,iZ_B0(4):iZ_E0(4),iPF_V2)
    uV3_K(1:nNodesX_K) => uPF_K(:,:,:,iZ_B0(4):iZ_E0(4),iPF_V3)

    ! --- Conserved Two-Moment Fields ---

    ALLOCATE( uCM_K(1:nDOFZ, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  , &
                    1:nSpecies, &
                    iZ_B0(1)-1:iZ_E0(1)+1, &
                    1:nCM) )
    ALLOCATE( uCM_L(1:nDOF_E, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  , &
                    1:nSpecies, &
                    iZ_B0(1)  :iZ_E0(1)+1, &
                    1:nCM) )
    ALLOCATE( uCM_R(1:nDOF_E, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  , &
                    1:nSpecies, &
                    iZ_B0(1)  :iZ_E0(1)+1, &
                    1:nCM) )

    E_K (1:nNodesZ_K) => uCM_K(:,:,:,:,:,iZ_B0(1):iZ_E0(1),iCM_E )
    F_d_1_K(1:nNodesZ_K) => uCM_K(:,:,:,:,:,iZ_B0(1):iZ_E0(1),iCM_F1)
    F_d_2_K(1:nNodesZ_K) => uCM_K(:,:,:,:,:,iZ_B0(1):iZ_E0(1),iCM_F2)
    F_d_3_K(1:nNodesZ_K) => uCM_K(:,:,:,:,:,iZ_B0(1):iZ_E0(1),iCM_F3)

    E_L (1:nNodesZ_E) => uCM_L(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCM_E )
    F_d_1_L(1:nNodesZ_E) => uCM_L(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCM_F1)
    F_d_2_L(1:nNodesZ_E) => uCM_L(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCM_F2)
    F_d_3_L(1:nNodesZ_E) => uCM_L(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCM_F3)

    E_R (1:nNodesZ_E) => uCM_R(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCM_E )
    F_d_1_R(1:nNodesZ_E) => uCM_R(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCM_F1)
    F_d_2_R(1:nNodesZ_E) => uCM_R(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCM_F2)
    F_d_3_R(1:nNodesZ_E) => uCM_R(:,:,:,:,:,iZ_B0(1):iZ_E0(1)+1,iCM_F3)

    ! --- Primitive Two-Moment Fields ---

    ALLOCATE( J_L (nNodesZ_E) )
    ALLOCATE( H_u_1_L(nNodesZ_E) )
    ALLOCATE( H_u_2_L(nNodesZ_E) )
    ALLOCATE( H_u_3_L(nNodesZ_E) )

    ALLOCATE( J_R (nNodesZ_E) )
    ALLOCATE( H_u_1_R(nNodesZ_E) )
    ALLOCATE( H_u_2_R(nNodesZ_E) )
    ALLOCATE( H_u_3_R(nNodesZ_E) )

    ! --- Fluxes ---

    ALLOCATE( NumericalFlux (nDOF_E,nCM, &
                             iZ_B0(2)  :iZ_E0(2)  , &
                             iZ_B0(3)  :iZ_E0(3)  , &
                             iZ_B0(4)  :iZ_E0(4)  , &
                             nSpecies, &
                             iZ_B0(1)  :iZ_E0(1)+1) )
    ALLOCATE( Flux_q        (nDOFZ,nCM, &
                             iZ_B0(2)  :iZ_E0(2)  , &
                             iZ_B0(3)  :iZ_E0(3)  , &
                             iZ_B0(4)  :iZ_E0(4)  , &
                             nSpecies, &
                             iZ_B0(1)  :iZ_E0(1)  ) )

    ! --- Increment ---

    ALLOCATE( dU_Z          (nDOFZ,nCM, &
                             iZ_B0(2)  :iZ_E0(2)  , &
                             iZ_B0(3)  :iZ_E0(3)  , &
                             iZ_B0(4)  :iZ_E0(4)  , &
                             nSpecies, &
                             iZ_B0(1)  :iZ_E0(1)  ) )

    ! --- Eigenvalues ---
                             
    ALLOCATE( Alpha(nDOFX, &
                    iZ_B0(2)  :iZ_E0(2)  , &
                    iZ_B0(3)  :iZ_E0(3)  , &
                    iZ_B0(4)  :iZ_E0(4)  ) )

    ! --- Derivatives ---

    ALLOCATE( dV_u_dX1   (nDOFX,0:3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )
    ALLOCATE( dV_u_dX2   (nDOFX,0:3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )
    ALLOCATE( dV_u_dX3   (nDOFX,0:3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )

    ALLOCATE( dV_d_dX1   (nDOFX,0:3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )
    ALLOCATE( dV_d_dX2   (nDOFX,0:3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )
    ALLOCATE( dV_d_dX3   (nDOFX,0:3, &
                          iZ_B0(2)  :iZ_E0(2)  , &
                          iZ_B0(3)  :iZ_E0(3)  , &
                          iZ_B0(4)  :iZ_E0(4)  ) )

    ! ---

    ALLOCATE( nIterations_L(nNodesZ_E) )
    ALLOCATE( nIterations_R(nNodesZ_E) )
    ALLOCATE( nIterations_K(nNodesZ_K) )

    ALLOCATE( PositionIndexZ_F(nNodesZ_E) )
    ALLOCATE( PositionIndexZ_K(nNodesZ_K) )

    ALLOCATE( IndexTableZ_F(7,nNodesZ_E) )
    ALLOCATE( IndexTableZ_K(7,nNodesZ_K) )

    DO iZ1 = iZ_B0(1), iZ_E0(1)+1
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iNode_E = 1, nDOF_E

      iZ_F = iNode_E &
             + ( iZ2 - iZ_B0(2) ) * nDOF_E &
             + ( iZ3 - iZ_B0(3) ) * nDOF_E * nZ_E(2) &
             + ( iZ4 - iZ_B0(4) ) * nDOF_E * nZ_E(2) * nZ_E(3) &
             + ( iS  - 1        ) * nDOF_E * nZ_E(2) * nZ_E(3) * nZ_E(4) &
             + ( iZ1 - iZ_B0(1) ) * nDOF_E * nZ_E(2) * nZ_E(3) * nZ_E(4) * nSpecies

      iX_F = iNode_E &
             + ( iZ2 - iZ_B0(2) ) * nDOF_E &
             + ( iZ3 - iZ_B0(3) ) * nDOF_E * nZ_E(2) &
             + ( iZ4 - iZ_B0(4) ) * nDOF_E * nZ_E(2) * nZ_E(3)

      PositionIndexZ_F(iZ_F) = iX_F

      IndexTableZ_F(:,iZ_F) = [ 1, iNode_E, iZ2, iZ3, iZ4, iS, iZ1 ]

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iNodeZ = iNodeE &
               + ( iNodeX - 1 ) * nDOFE

      iZ_K = iNodeZ &
             + ( iZ2 - iZ_B0(2) ) * nDOFZ &
             + ( iZ3 - iZ_B0(3) ) * nDOFZ * nZ(2) &
             + ( iZ4 - iZ_B0(4) ) * nDOFZ * nZ(2) * nZ(3) &
             + ( iS  - 1        ) * nDOFZ * nZ(2) * nZ(3) * nZ(4) &
             + ( iZ1 - iZ_B0(1) ) * nDOFZ * nZ(2) * nZ(3) * nZ(4) * nSpecies

      iX_K = iNodeX &
             + ( iZ2 - iZ_B0(2) ) * nDOFX &
             + ( iZ3 - iZ_B0(3) ) * nDOFX * nZ(2) &
             + ( iZ4 - iZ_B0(4) ) * nDOFX * nZ(2) * nZ(3)

      PositionIndexZ_K(iZ_K) = iX_K

      IndexTableZ_K(:,iZ_K) = [ iNodeE, iNodeX, iZ2, iZ3, iZ4, iS, iZ1 ]

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    RETURN
  END SUBROUTINE InitializeIncrement_ObserverCorrections

  SUBROUTINE FinalizeIncrement_ObserverCorrections

    DEALLOCATE( PositionIndexZ_F, PositionIndexZ_K )
    DEALLOCATE( IndexTableZ_F, IndexTableZ_K )
    DEALLOCATE( uGF_K, uPF_K, uCM_K, uCM_L, uCM_R )
    DEALLOCATE( J_L, H_u_1_L, H_u_2_L, H_u_3_L )
    DEALLOCATE( J_R, H_u_1_R, H_u_2_R, H_u_3_R )
    DEALLOCATE( NumericalFlux, Flux_q, dU_Z, Alpha)
    DEALLOCATE( dV_u_dX1, dV_u_dX2, dV_u_dX3 )
    DEALLOCATE( dV_d_dX1, dV_d_dX2, dV_d_dX3 )
    DEALLOCATE( nIterations_L, nIterations_R, nIterations_K )

    NULLIFY( Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K )
    NULLIFY( uV1_K, uV2_K, uV3_K )
    NULLIFY( E_K, F_d_1_K, F_d_2_K, F_d_3_K )
    NULLIFY( E_L, F_d_1_L, F_d_2_L, F_d_3_L )
    NULLIFY( E_R, F_d_1_R, F_d_2_R, F_d_3_R )

  END SUBROUTINE FinalizeIncrement_ObserverCorrections

END MODULE TwoMoment_DiscretizationModule_Streaming_FMC