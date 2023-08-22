MODULE TwoMoment_DiscretizationModule_Streaming_FMC

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, &
    SqrtTiny, Pi
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, nDOFX_X2, nDOFX_X3
  USE ReferenceElementModule, ONLY: &
    nDOF_E, &
    nDOF_X1, nDOF_X2, nDOF_X3, &
    Weights_q
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE, iGE_Ep1, iGE_Ep2, iGE_Ep3
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nPF, iPF_V1, iPF_V2, iPF_V3
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nSpecies, &
    nCM, iCM_E, iCM_F1, iCM_F2, iCM_F3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit

  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    uFD_K, uS1_K, uS2_K, uS3_K, &
    uFD_L, uS1_L, uS2_L, uS3_L, &
    uFD_R, uS1_R, uS2_R, uS3_R

  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    E_K, F_d_1_K, F_d_2_K, F_d_3_K, &
    E_L, F_d_1_L, F_d_2_L, F_d_3_L, &
    E_R, F_d_1_R, F_d_2_R, F_d_3_R

  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    J_K, H_u_1_K, H_u_2_K, H_u_3_K

  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    V_u_1_K, V_u_2_K, V_u_3_K

  INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
  INTEGER :: nZ   (4), nK_X , nK_Z , nNodesX_K, nNodesZ_K
  INTEGER :: nZ_E (4), nE_X , nE_Z , nNodesZ_E
  INTEGER :: nZ_X1(4), nX1_X, nX1_Z, nNodesX_X1, nNodesZ_X1
  INTEGER :: nZ_X2(4), nX2_X, nX2_Z, nNodesX_X2, nNodesZ_X2
  INTEGER :: nZ_X3(4), nX3_X, nX3_Z, nNodesX_X3, nNodesZ_X3

CONTAINS

  SUBROUTINE ComputeIncrement_TwoMoment_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_M, dU_M )

    ! --- Input/Output variables ---
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: GE(1:nDOFE, iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: GX(1:nDOFX, &
                               iZ_B1(2):iZ_E1(2), &
                               iZ_B1(3):iZ_E1(3), &
                               iZ_B1(4):iZ_E1(4), &
                               1:nGF)
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

    ! CALL ApplyBoundaryConditions_TwoMoment &
    !        ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R )

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
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_M, dU_M )

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
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_M, dU_M )

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

    ! --- Conserved Fluid Fields ---

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

    INTEGER  :: iZP_B0(4), iZP_E0(4)

    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN

    ! --- Permuted Phase Space Limits ---

    iZP_B0(1) = iZ_B0(1) ; iZP_E0(1) = iZ_E0(1)
    iZP_B0(2) = iZ_B0(3) ; iZP_E0(2) = iZ_E0(3)
    iZP_B0(3) = iZ_B0(4) ; iZP_E0(3) = iZ_E0(4)
    iZP_B0(4) = iZ_B0(2) ; iZP_E0(4) = iZ_E0(2)

    CALL InitializeIncrement_Divergence_X &
           ( iZP_B0, iZP_E0, nDOFX_X1, nDOF_X1, &
             uPF_K, uPF_L, uPF_R, uCM_K, uCM_L, uCM_R )

  END SUBROUTINE ComputeIncrement_Divergence_X1

  SUBROUTINE InitializeIncrement_Divergence_X &
    ( iZP_B0, iZP_E0, nDOFX_X, nDOFZ_X, &
      uPF_K, uPF_L, uPF_R, uCM_K, uCM_L, uCM_R )

    INTEGER, INTENT(in) :: iZP_B0(4), iZP_E0(4) ! Permuted limits
    INTEGER, INTENT(in) :: nDOFX_X ! nDOFX_X1, ...
    INTEGER, INTENT(in) :: nDOFZ_X ! nDOFZ_X1, ...

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

  END SUBROUTINE InitializeIncrement_Divergence_X

END MODULE TwoMoment_DiscretizationModule_Streaming_FMC