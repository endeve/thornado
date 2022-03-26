MODULE MHD_DiscretizationModule_Relativistic

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    SqrtTiny, &
    Third, &
    Half, &
    One, &
    Two, &
    Three, &
    Four
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    WeightsX_X1, &
    WeightsX_X2, &
    WeightsX_X3, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    dLXdX2_q, &
    dLXdX3_q, &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Psi, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Phi_N, &
    CoordinateSystem
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi, &
    nPM, &
    iPM_D, &
    iPM_V1, &
    iPM_V2, &
    iPM_V3, &
    iPM_E, &
    iPM_Ne, &
    iPM_B1, &
    iPM_B2, &
    iPM_B3, &
    iPM_Chi, &
    nDM, &
    iDM_Sh_X1, &
    iDM_Sh_X2, &
    iDM_Sh_X3
  USE MHD_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_MHD
  USE MHD_UtilitiesModule, ONLY: &
    ComputePrimitive_MHD, &
    Eigenvalues_MHD, &
    Flux_X1_MHD, &
    NumericalFlux_MHD_X1
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_MHD_DG_Explicit

  LOGICAL  :: UseDivergenceCleaning
  REAL(DP) :: Time, DampingParameter

  REAL(DP), POINTER, CONTIGUOUS :: &
    Gm_dd_11_K(:), Gm_dd_22_K(:), Gm_dd_33_K(:), SqrtGm_K(:), &
    Gm_dd_11_F(:), Gm_dd_22_F(:), Gm_dd_33_F(:), SqrtGm_F(:), &
    Beta_1_K  (:), Beta_2_K  (:), Beta_3_K  (:), Alpha_K (:), &
    Beta_1_F  (:), Beta_2_F  (:), Beta_3_F  (:), Alpha_F (:)

  REAL(DP), POINTER, CONTIGUOUS :: &
    uD_K(:), uS1_K(:), uS2_K(:), uS3_K(:), uE_K(:), uNe_K(:), &
    uB1_K(:), uB2_K(:), uB3_K(:), uChi_K(:), &
    uD_L(:), uS1_L(:), uS2_L(:), uS3_L(:), uE_L(:), uNe_L(:), &
    uB1_L(:), uB2_L(:), uB3_L(:), uChi_L(:), &
    uD_R(:), uS1_R(:), uS2_R(:), uS3_R(:), uE_R(:), uNe_R(:), &
    uB1_R(:), uB2_R(:), uB3_R(:), uChi_R(:)
 
  REAL(DP), ALLOCATABLE :: &
    pD_K(:), pV1_K(:), pV2_K(:), pV3_K(:), pE_K(:), pNe_K(:), &
    pB1_K(:), pB2_K(:), pB3_K(:), pChi_K(:), &
    pD_L(:), pV1_L(:), pV2_L(:), pV3_L(:), pE_L(:), pNe_L(:), &
    pB1_L(:), pB2_L(:), pB3_L(:), pChi_L(:), & 
    pD_R(:), pV1_R(:), pV2_R(:), pV3_R(:), pE_R(:), pNe_R(:), &
    pB1_R(:), pB2_R(:), pB3_R(:), pChi_R(:)
 
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: &
    IndexTableX_F(:,:), IndexTableX_V(:,:)

  INTEGER :: nX   (3), nX_K,  nNodesX_K
  INTEGER :: nX_X1(3), nX1_X, nNodesX_X1
  INTEGER :: nX_X2(3), nX2_X, nNodesX_X2
  INTEGER :: nX_X3(3), nX3_X, nNodesX_X3
 
CONTAINS


  SUBROUTINE ComputeIncrement_MHD_DG_Explicit &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, &
      SuppressBC_Option, UseDivergenceCleaning_Option, &
      DampingParameter_Option )
    
    INTEGER,  INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)            :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in),  OPTIONAL :: &
      SuppressBC_Option
    LOGICAL,  INTENT(in),  OPTIONAL :: &
      UseDivergenceCleaning_Option
    REAL(DP), INTENT(in),  OPTIONAL :: &
      DampingParameter_Option
    REAL(DP), INTENT(inout)         :: &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)         :: &
      D (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out)           :: &
      dU(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3, iCM
    LOGICAL  :: SuppressBC
    REAL(DP) :: tau(nDOFX,iX_B1(1):iX_E1(1), &
                          iX_B1(2):iX_E1(2), &
                          iX_B1(3):iX_E1(3))

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    UseDivergenceCleaning = .FALSE.
    IF( PRESENT( UseDivergenceCleaning_Option ) ) &
      UseDivergenceCleaning = UseDivergenceCleaning_Option

    DampingParameter = 0.0_DP
    IF( UseDivergenceCleaning .AND. PRESENT( DampingParameter_Option ) )THEN
      DampingParameter = DampingParameter_Option
    END IF

    IF( .NOT. SuppressBC )THEN

      CALL ApplyBoundaryConditions_MHD &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    END IF

    CALL InitializeIncrement_MHD &
           ( iX_B0, iX_E0, iX_B1, iX_E1 )

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      tau(iNX,iX1,iX2,iX3) = One

    END DO
    END DO
    END DO
    END DO

    DO iCM = 1, nCM
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      dU(iNX,iX1,iX2,iX3,iCM) = Zero

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL ComputeIncrement_MHD_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    !PRINT*
    !PRINT*, 'Increments after divergence.'
    !PRINT*, 'B1: ', dU(:,:,:,:,iCM_B1)
    !PRINT*, 'S3: ', dU(:,:,:,:,iCM_S3)
    !PRINT*, 'B3: ', dU(:,:,:,:,iCM_B3)

    ! --- Multiply Inverse Mass Matrix ---

    DO iCM = 1, nCM
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      dU(iNX,iX1,iX2,iX3,iCM) &
        = dU(iNX,iX1,iX2,iX3,iCM) &
            * tau(iNX,iX1,iX2,iX3) &
            / ( WeightsX_q(iNX) * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                  * dX1(iX1) * dX2(iX2) * dX3(iX3) )

    END DO
    END DO
    END DO
    END DO
    END DO

    !PRINT*, 'Increments after mass matrix multiplication.'
    !PRINT*, 'B1: ', dU(:,:,:,:,iCM_B1)
    !PRINT*, 'S3: ', dU(:,:,:,:,iCM_S3)
    !PRINT*, 'B3: ', dU(:,:,:,:,iCM_B3)

    CALL ComputeIncrement_Geometry &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, dU )

    !PRINT*, 'Increments after geometry.'
    !PRINT*, 'B1: ', dU(:,:,:,:,iCM_B1)
    !PRINT*, 'S3: ', dU(:,:,:,:,iCM_S3)
    !PRINT*, 'B3: ', dU(:,:,:,:,iCM_B3)

    CALL ComputeIncrement_Gravity &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    !PRINT*, 'Increments after gravity.'
    !PRINT*, 'B1: ', dU(:,:,:,:,iCM_B1)
    !PRINT*, 'S3: ', dU(:,:,:,:,iCM_S3)
    !PRINT*, 'B3: ', dU(:,:,:,:,iCM_B3)

    CALL FinalizeIncrement_MHD

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_MHD_DG_Explicit


  SUBROUTINE ComputeIncrement_MHD_Divergence_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM), &
      D (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDM)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)

    INTEGER  :: iNX, iNX_X, iNX_K, iX1, iX2, iX3, iCM, iGF
    INTEGER  :: iXP_B0(3), iXP_E0(3)

    REAL(DP) :: AlphaMns, AlphaPls
    REAL(DP) :: P_L, P_R, Cs_L, Cs_R, P_K

    REAL(DP) :: EigVals_L(2), EigVals_R(2)
    REAL(DP) :: Flux_L   (nCM), Flux_R   (nCM)
    REAL(DP) :: Flux_F   (nCM), Flux_K   (nCM)
    REAL(DP) :: uCM_L_nCM(nCM), uCM_R_nCM(nCM)

    ! --- Geometry Fields ---

    REAL(DP) :: &
      G_K(nDOFX, &
          iX_B0(2)  :iX_E0(2),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(1)-1:iX_E0(1)+1, &
          nGF)
    REAL(DP) :: &
      G_F(nDOFX_X1, &
          iX_B0(2)  :iX_E0(2),   &
          iX_B0(3)  :iX_E0(3),   &
          iX_B0(1)  :iX_E0(1)+1, &
          nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP) :: &
      uCM_K(nDOFX, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)-1:iX_E0(1)+1, &
            nCM)
    REAL(DP) :: &
      uCM_L(nDOFX_X1, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)  :iX_E0(1)+1, &
            nCM)
    REAL(DP) :: &
      uCM_R(nDOFX_X1, &
            iX_B0(2)  :iX_E0(2),   &
            iX_B0(3)  :iX_E0(3),   &
            iX_B0(1)  :iX_E0(1)+1, &
            nCM)

    ! --- Diagnostic Fields ---

    REAL(DP) :: &
      uDM_L(iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            iX_B0(1):iX_E0(1)+1,2)
    REAL(DP) :: &
      uDM_R(iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            iX_B0(1):iX_E0(1)+1,2)

    ! --- Fluxes ---

    REAL(DP) :: &
      NumericalFlux(nDOFX_X1,nCM, &
                    iX_B0(2):iX_E0(2), &
                    iX_B0(3):iX_E0(3), &
                    iX_B0(1):iX_E0(1)+1)
    REAL(DP) :: &
      Flux_q(nDOFX,nCM, &
             iX_B0(2):iX_E0(2), &
             iX_B0(3):iX_E0(3), &
             iX_B0(1):iX_E0(1))

    ! --- X1 Increment ---

    REAL(DP) :: &
      dU_X1(nDOFX,nCM, &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            iX_B0(1):iX_E0(1))

    IF( iX_E0(1) .EQ. iX_B0(1) ) RETURN

    ! --- Permuted Limits ---

    iXP_B0(1) = iX_B0(2) ; iXP_E0(1) = iX_E0(2)
    iXP_B0(2) = iX_B0(3) ; iXP_E0(2) = iX_E0(3)
    iXP_B0(3) = iX_B0(1) ; iXP_E0(3) = iX_E0(1)

    ASSOCIATE( dX2 => MeshX(2) % Width, dX3 => MeshX(3) % Width )

    CALL InitializeIncrement_Divergence &
           ( iXP_B0, iXP_E0, nDOFX_X1, &
             G_K, G_F, uCM_K, uCM_L, uCM_R )

    ! --- Permute data ---

    DO iGF = 1           , nGF
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

      G_K(iNX,iX2,iX3,iX1,iGF) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iCM = 1           , nCM
    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

      uCM_K(iNX,iX2,iX3,iX1,iCM) = U(iNX,iX1,iX2,iX3,iCM)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX1 = iX_B0(1), iX_E0(1) + 1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      uDM_L(iX2,iX3,iX1,1) = D(1,iX1-1,iX2,iX3,iDM_Sh_X2)
      uDM_L(iX2,iX3,iX1,2) = D(1,iX1-1,iX2,iX3,iDM_Sh_X3)
      uDM_R(iX2,iX3,iX1,1) = D(1,iX1  ,iX2,iX3,iDM_Sh_X2)
      uDM_R(iX2,iX3,iX1,2) = D(1,iX1  ,iX2,iX3,iDM_Sh_X3)

    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale factor (X1) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_1), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_1), nDOFX_X1 )

    ! --- Scale factor (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_2), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_2), nDOFX_X1 )

    ! --- Scale factor (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_h_3), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_h_3), nDOFX_X1 )

    ! --- Lapse function ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_Alpha), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Alpha), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Alpha), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Alpha), nDOFX_X1 )

    ! --- Shift vector (X1) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_Beta_1), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_1), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_1), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_1), nDOFX_X1 )

    ! --- Shift vector (X2) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_Beta_2), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_2), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_2), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_2), nDOFX_X1 )

    ! --- Shift vector (X3) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iGF_Beta_3), nDOFX, Zero, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_3), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_3), nDOFX, Half, &
             G_F(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iGF_Beta_3), nDOFX_X1 )

    ! --- Compute metric and metric determinant on faces ---

    DO iX1   = iX_B0(1), iX_E0(1) + 1
    DO iX3   = iX_B0(3), iX_E0(3)
    DO iX2   = iX_B0(2), iX_E0(2)
    DO iNX_X = 1       , nDOFX_X1

      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_1) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1), SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_2) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_2), SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_h_3) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_3), SqrtTiny )

      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_11) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1     )**2, SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_22) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_2     )**2, SqrtTiny )
      G_F             (iNX_X,iX2,iX3,iX1,iGF_Gm_dd_33) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_3     )**2, SqrtTiny )

      G_F             (iNX_X,iX2,iX3,iX1,iGF_SqrtGm) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_h_1   ) &
                 * G_F(iNX_X,iX2,iX3,iX1,iGF_h_2   ) &
                 * G_F(iNX_X,iX2,iX3,iX1,iGF_h_3   ), SqrtTiny )

      G_F             (iNX_X,iX2,iX3,iX1,iGF_Alpha) &
        = MAX( G_F    (iNX_X,iX2,iX3,iX1,iGF_Alpha), SqrtTiny )

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iCM = 1, nCM

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               uCM_K(1,iX_B0(2),iX_B0(3),iX_B0(1)-1,iCM), nDOFX, Zero, &
               uCM_L(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               uCM_K(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX, Zero, &
               uCM_R(1,iX_B0(2),iX_B0(3),iX_B0(1)  ,iCM), nDOFX_X1 )

    END DO

    ! --- Numerical Flux ---

    !PRINT*, 'Computing the primitive variables for the left state.'

    CALL ComputePrimitive_MHD &
           ( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L, &
             uB1_L, uB2_L, uB3_L, uChi_L, &
             pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, &
             pB1_L, pB2_L, pB3_L, pChi_L, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             Alpha_F, Beta_1_F, Beta_2_F, Beta_3_F )

    !PRINT*, 'Computing the primitive variables for the right state.'

    CALL ComputePrimitive_MHD &
           ( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R, &
             uB1_R, uB2_R, uB3_R, uChi_R, &
             pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
             pB1_R, pB2_R, pB3_R, pChi_R, &
             Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, &
             Alpha_F, Beta_1_F, Beta_2_F, Beta_3_F )

    DO iNX_X = 1, nNodesX_X1

      ! --- Left state ---

      CALL ComputePressureFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_L(iNX_X), pE_L(iNX_X), pNe_L(iNX_X), Cs_L )

      EigVals_L &
        = Eigenvalues_MHD &
            ( pV1_L     (iNX_X), &
              Cs_L             , &
              Gm_dd_11_F(iNX_X), &
              pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              pB1_L     (iNX_X), &
              pB2_L     (iNX_X), &
              pB3_L     (iNX_X), & 
              pChi_L    (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      Flux_L &
        = Flux_X1_MHD &
            ( pD_L      (iNX_X), &
              pV1_L     (iNX_X), &
              pV2_L     (iNX_X), &
              pV3_L     (iNX_X), &
              pE_L      (iNX_X), &
              pNe_L     (iNX_X), &
              pB1_L     (iNX_X), &
              pB2_L     (iNX_X), &
              pB3_L     (iNX_X), &
              pChi_L    (iNX_X), &
              P_L              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      IF( .FALSE. )THEN

      !PRINT*
      !PRINT*, 'Left states and flux for iNX_X = ', iNX_X
      !PRINT*, 'D_L:      ', pD_L(iNX_X)
      !PRINT*, 'V1_L:     ', pV1_L(iNX_X)
      !PRINT*, 'V2_L:     ', pV2_L(iNX_X)
      !PRINT*, 'V3_L:     ', pV3_L(iNX_X)
      !PRINT*, 'E_L:      ', pE_L(iNX_X)
      !PRINT*, 'Ne_L:     ', pNe_L(iNX_X)
      !PRINT*, 'B1_L:     ', pB1_L(iNX_X)
      !PRINT*, 'B2_L:     ', pB2_L(iNX_X)
      !PRINT*, 'B3_L:     ', pB3_L(iNX_X)
      !PRINT*, 'Chi_L:    ', pChi_L(iNX_X)
      !PRINT*, 'P_L:      ', P_L
      !PRINT*, 'Gm_dd_11: ', Gm_dd_11_F(iNX_X)
      !PRINT*, 'Flux_L:   ', Flux_L

      END IF

      ! --- Right state ---

      CALL ComputePressureFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( pD_R(iNX_X), pE_R(iNX_X), pNe_R(iNX_X), Cs_R )

      EigVals_R &
        = Eigenvalues_MHD &
            ( pV1_R     (iNX_X), &
              Cs_R             , &
              Gm_dd_11_F(iNX_X), &
              pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              pB1_R     (iNX_X), &
              pB2_R     (iNX_X), &
              pB3_R     (iNX_X), & 
              pChi_R    (iNX_X), &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      Flux_R &
        = Flux_X1_MHD &
            ( pD_R      (iNX_X), &
              pV1_R     (iNX_X), &
              pV2_R     (iNX_X), &
              pV3_R     (iNX_X), &
              pE_R      (iNX_X), &
              pNe_R     (iNX_X), &
              pB1_R     (iNX_X), &
              pB2_R     (iNX_X), &
              pB3_R     (iNX_X), &
              pChi_R    (iNX_X), &
              P_R              , &
              Gm_dd_11_F(iNX_X), &
              Gm_dd_22_F(iNX_X), &
              Gm_dd_33_F(iNX_X), &
              Alpha_F   (iNX_X), &
              Beta_1_F  (iNX_X), &
              Beta_2_F  (iNX_X), &
              Beta_3_F  (iNX_X), &
              UseDivergenceCleaning )

      IF( .FALSE. )THEN
     
      !PRINT*
      !PRINT*, 'Right states and flux for iNX_X = ', iNX_X
      !PRINT*, 'D_R:      ', pD_R(iNX_X)
      !PRINT*, 'V1_R:     ', pV1_R(iNX_X)
      !PRINT*, 'V2_R:     ', pV2_R(iNX_X)
      !PRINT*, 'V3_R:     ', pV3_R(iNX_X)
      !PRINT*, 'E_R:      ', pE_R(iNX_X)
      !PRINT*, 'Ne_R:     ', pNe_R(iNX_X)
      !PRINT*, 'B1_R:     ', pB1_R(iNX_X)
      !PRINT*, 'B2_R:     ', pB2_R(iNX_X)
      !PRINT*, 'B3_R:     ', pB3_R(iNX_X)
      !PRINT*, 'Chi_R:    ', pChi_R(iNX_X)
      !PRINT*, 'P_R:      ', P_R
      !PRINT*, 'Gm_dd_11: ', Gm_dd_11_F(iNX_X)
      !PRINT*, 'Flux_R:   ', Flux_R

      END IF

      ! --- Numerical flux ---

      AlphaMns &
        = MAX( Zero, MAXVAL( - EigVals_L ), MAXVAL( - EigVals_R ) )

      AlphaPls &
        = MAX( Zero, MAXVAL( + EigVals_L ), MAXVAL( + EigVals_R ) )

      iNX = IndexTableX_F(1,iNX_X)
      iX2 = IndexTableX_F(2,iNX_X)
      iX3 = IndexTableX_F(3,iNX_X)
      iX1 = IndexTableX_F(4,iNX_X)

      DO iCM = 1, nCM

        uCM_L_nCM(iCM) = uCM_L(iNX,iX2,iX3,iX1,iCM)
        uCM_R_nCM(iCM) = uCM_R(iNX,iX2,iX3,iX1,iCM)

      END DO

      Flux_F &
        = NumericalFlux_MHD_X1 &
            ( uCM_L_nCM           , &
              uCM_R_nCM           , &
              Flux_L              , &
              Flux_R              , &
              AlphaPls            , &
              AlphaMns            )

      DO iCM = 1, nCM

        NumericalFlux(iNX,iCM,iX2,iX3,iX1) &
          = Flux_F(iCM) &
              * Alpha_F(iNX_X) * SqrtGm_F(iNX_X) * dX2(iX2) * dX3(iX3) &
              * WeightsX_X1(iNX)

      END DO ! iCM

    END DO ! iNX_X

    ! --- Surface Contribution ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX_X1, + One, LX_X1_Dn, nDOFX_X1, &
             NumericalFlux(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, &
             Zero, dU_X1, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX_X1, - One, LX_X1_Up, nDOFX_X1, &
             NumericalFlux(1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, &
             One,  dU_X1, nDOFX )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    !PRINT*, 'Computing the primitive variables for the volume term.'

    CALL ComputePrimitive_MHD &
           ( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K, &
             uB1_K, uB2_K, uB3_K, uChi_K, &
             pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
             pB1_K, pB2_K, pB3_K, pChi_K, &
             Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, &
             Alpha_K, Beta_1_K, Beta_2_K, Beta_3_K )

    DO iNX_K = 1, nNodesX_K

      CALL ComputePressureFromPrimitive &
             ( pD_K(iNX_K), pE_K(iNX_K), pNe_K(iNX_K), P_K )

      Flux_K &
        = Flux_X1_MHD &
            ( pD_K      (iNX_K), &
              pV1_K     (iNX_K), &
              pV2_K     (iNX_K), &
              pV3_K     (iNX_K), &
              pE_K      (iNX_K), &
              pNe_K     (iNX_K), &
              pB1_K     (iNX_K), &
              pB2_K     (iNX_K), &
              pB3_K     (iNX_K), &
              pChi_K    (iNX_K), &
              P_K              , &
              Gm_dd_11_K(iNX_K), &
              Gm_dd_22_K(iNX_K), &
              Gm_dd_33_K(iNX_K), &
              Alpha_K   (iNX_K), &
              Beta_1_K  (iNX_K), &
              Beta_2_K  (iNX_K), &
              Beta_3_K  (iNX_K), &
              UseDivergenceCleaning )

      iNX = IndexTableX_V(1,iNX_K)
      iX2 = IndexTableX_V(2,iNX_K)
      iX3 = IndexTableX_V(3,iNX_K)
      iX1 = IndexTableX_V(4,iNX_K)

      DO iCM = 1, nCM

        Flux_q(iNX,iCM,iX2,iX3,iX1) &
          = Flux_K(iCM) &
              * Alpha_K(iNX_K) * SqrtGm_K(iNX_K) &
              * dX2(iX2) * dX3(iX3) * WeightsX_q(iNX)

      END DO

    END DO

    ! --- Contribution from Volume ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nX_K*nCM, nDOFX, One, dLXdX1_q, nDOFX, &
             Flux_q, nDOFX, One, dU_X1, nDOFX )

    DO iCM = 1       , nCM
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      dU    (iNX,iX1,iX2,iX3,iCM) &
        = dU(iNX,iX1,iX2,iX3,iCM) &
            + dU_X1(iNX,iCM,iX2,iX3,iX1)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeIncrement_Divergence

    END ASSOCIATE

  END SUBROUTINE ComputeIncrement_MHD_Divergence_X1


  SUBROUTINE ComputeIncrement_Geometry &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in)    :: &
      tau(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    CALL ComputeIncrement_Geometry_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, dU )

  END SUBROUTINE ComputeIncrement_Geometry


  SUBROUTINE ComputeIncrement_Geometry_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, tau, dU )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF), &
      U (1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)
    REAL(DP), INTENT(in)           :: &
      tau(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))
    REAL(DP), INTENT(inout)        :: &
      dU(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCM)

    INTEGER :: iX1, iX2, iX3, iNX, iCM, iGF, i, k, iDim
    INTEGER :: nGF_K

    REAL(DP) :: P(nPM)
    REAL(DP) :: Pressure
    REAL(DP) :: W, B0u, bSq
    REAL(DP) :: PressureTensor(3,3, nDOFX,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3))
    REAL(DP) :: PressureTensorTrace(nDOFX,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3))

    REAL(DP) :: G_K_X1 (nDOFX,   nGF+1,iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: G_Dn_X1(nDOFX_X1,nGF+1,iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(1)  :iX_E0(1))
    REAL(DP) :: G_Up_X1(nDOFX_X1,nGF+1,iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(1)  :iX_E0(1))
    REAL(DP) :: dGdX1  (nDOFX,   nGF+1,iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(1)  :iX_E0(1))

    REAL(DP) :: G_K_X2 (nDOFX,   nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: G_Dn_X2(nDOFX_X2,nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(2)  :iX_E0(2))
    REAL(DP) :: G_Up_X2(nDOFX_X2,nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(2)  :iX_E0(2))
    REAL(DP) :: dGdX2  (nDOFX,   nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(3)  :iX_E0(3), &
                                       iX_B0(2)  :iX_E0(2))

    REAL(DP) :: G_K_X3 (nDOFX,   nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: G_Dn_X3(nDOFX_X3,nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3))
    REAL(DP) :: G_Up_X3(nDOFX_X3,nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3))
    REAL(DP) :: dGdX3  (nDOFX,   nGF+1,iX_B0(1)  :iX_E0(1), &
                                       iX_B0(2)  :iX_E0(2), &
                                       iX_B0(3)  :iX_E0(3))

    REAL(DP) :: EnergyDensitySourceTerms(7,nDOFX,iX_B0(1):iX_E0(1), &
                                                 iX_B0(2):iX_E0(2), &
                                                 iX_B0(3):iX_E0(3))

    REAL(DP) :: PressureTensor_ud(3,3)
    REAL(DP) :: DivGridVolume

    REAL(DP) :: GradPsiF(nDOFX)

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    nGF_K = ( nGF + 1 ) * PRODUCT( nX )

    ! --- Permute data ---

    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iGF = 1           , nGF
    DO iNX = 1           , nDOFX

      G_K_X1(iNX,iGF,iX2,iX3,iX1) = G(iNX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX1 = iX_B0(1) - 1, iX_E0(1) + 1
    DO iX3 = iX_B0(3)    , iX_E0(3)
    DO iX2 = iX_B0(2)    , iX_E0(2)
    DO iNX = 1           , nDOFX

      G_K_X1(iNX,nGF+1,iX2,iX3,iX1) &
        =   G_K_X1(iNX,iGF_SqrtGm,iX2,iX3,iX1) &
          * G_K_X1(iNX,iGF_Beta_1,iX2,iX3,iX1)

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate to faces ---

    ! --- X1 ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_K, nDOFX, One , LX_X1_Up, nDOFX_X1, &
             G_K_X1 (1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Zero, &
             G_Up_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_K, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K_X1 (1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX, Half, &
             G_Up_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_K, nDOFX, One , LX_X1_Up, nDOFX_X1, &
             G_K_X1 (1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), nDOFX, Zero, &
             G_Dn_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )
    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nGF_K, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             G_K_X1 (1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Half, &
             G_Dn_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

      ! --- Compute metric on faces from scale factors ---

    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iNX = 1       , nDOFX_X1

      G_Up_X1         (iNX,iGF_Gm_dd_11,iX2,iX3,iX1) &
        = MAX( G_Up_X1(iNX,iGF_h_1     ,iX2,iX3,iX1)**2, SqrtTiny )
      G_Up_X1         (iNX,iGF_Gm_dd_22,iX2,iX3,iX1) &
        = MAX( G_Up_X1(iNX,iGF_h_2     ,iX2,iX3,iX1)**2, SqrtTiny )
      G_Up_X1         (iNX,iGF_Gm_dd_33,iX2,iX3,iX1) &
        = MAX( G_Up_X1(iNX,iGF_h_3     ,iX2,iX3,iX1)**2, SqrtTiny )

      G_Up_X1        (iNX,iGF_SqrtGm   ,iX2,iX3,iX1) &
        = G_Up_X1    (iNX,iGF_h_1      ,iX2,iX3,iX1) &
            * G_Up_X1(iNX,iGF_h_2      ,iX2,iX3,iX1) &
            * G_Up_X1(iNX,iGF_h_3      ,iX2,iX3,iX1)

      G_Up_X1        (iNX,nGF+1        ,iX2,iX3,iX1) &
        = G_Up_X1    (iNX,iGF_SqrtGm   ,iX2,iX3,iX1) &
            * G_Up_X1(iNX,iGF_Beta_1   ,iX2,iX3,iX1)

      G_Dn_X1         (iNX,iGF_Gm_dd_11,iX2,iX3,iX1) &
        = MAX( G_Dn_X1(iNX,iGF_h_1     ,iX2,iX3,iX1)**2, SqrtTiny )
      G_Dn_X1         (iNX,iGF_Gm_dd_22,iX2,iX3,iX1) &
        = MAX( G_Dn_X1(iNX,iGF_h_2     ,iX2,iX3,iX1)**2, SqrtTiny )
      G_Dn_X1         (iNX,iGF_Gm_dd_33,iX2,iX3,iX1) &
        = MAX( G_Dn_X1(iNX,iGF_h_3     ,iX2,iX3,iX1)**2, SqrtTiny )

      G_Dn_X1        (iNX,iGF_SqrtGm   ,iX2,iX3,iX1) &
        = G_Dn_X1    (iNX,iGF_h_1      ,iX2,iX3,iX1) &
            * G_Dn_X1(iNX,iGF_h_2      ,iX2,iX3,iX1) &
            * G_Dn_X1(iNX,iGF_h_3      ,iX2,iX3,iX1)

      G_Dn_X1        (iNX,nGF+1        ,iX2,iX3,iX1) &
        = G_Dn_X1    (iNX,iGF_SqrtGm   ,iX2,iX3,iX1) &
            * G_Dn_X1(iNX,iGF_Beta_1   ,iX2,iX3,iX1)

    END DO
    END DO
    END DO
    END DO

    ! --- Compute derivatives (X1) ---

    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1       , nGF + 1
    DO iNX = 1       , nDOFX_X1

      G_Dn_X1(iNX,iGF,iX2,iX3,iX1) &
        = G_Dn_X1(iNX,iGF,iX2,iX3,iX1) * WeightsX_X1(iNX)

      G_Up_X1(iNX,iGF,iX2,iX3,iX1) &
        = G_Up_X1(iNX,iGF,iX2,iX3,iX1) * WeightsX_X1(iNX)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1       , nGF + 1
    DO iNX = 1       , nDOFX

      G_K_X1(iNX,iGF,iX2,iX3,iX1) &
        = G_K_X1(iNX,iGF,iX2,iX3,iX1) * WeightsX_q(iNX)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, +One, LX_X1_Up, nDOFX_X1, &
             G_Up_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1, Zero, &
             dGdX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX_X1, -One, LX_X1_Dn, nDOFX_X1, &
             G_Dn_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX_X1, One,  &
             dGdX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, nGF_K, nDOFX   , -One, dLXdX1_q, nDOFX   , &
             G_K_X1 (1,1,iX_B0(2),iX_B0(3),iX_B0(1)), nDOFX   , One , &
             dGdX1, nDOFX )

    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iGF = 1       , nGF + 1
    DO iNX = 1       , nDOFX

      dGdX1(iNX,iGF,iX2,iX3,iX1) &
        = dGdX1(iNX,iGF,iX2,iX3,iX1) / ( WeightsX_q(iNX) * dX1(iX1) )

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iGF = 1       , nGF + 1
    DO iNX = 1       , nDOFX

      dGdX2(iNX,iGF,iX1,iX3,iX2) = Zero

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iGF = 1       , nGF + 1
    DO iNX = 1       , nDOFX

      dGdX3(iNX,iGF,iX1,iX2,iX3) = Zero

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Contributions from time-independent metric ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      !PRINT*, 'In cell: ', iX1, iX2, iX3
      !PRINT*, 'In node: ', iNX

      !PRINT*, 'Computing the primitive variables for the time-independent metric source term.'

      CALL ComputePrimitive_MHD  &
             ( U(   iNX,iX1,iX2,iX3,iCM_D ),  &
               U(   iNX,iX1,iX2,iX3,iCM_S1),  &
               U(   iNX,iX1,iX2,iX3,iCM_S2),  &
               U(   iNX,iX1,iX2,iX3,iCM_S3),  &
               U(   iNX,iX1,iX2,iX3,iCM_E ),  &
               U(   iNX,iX1,iX2,iX3,iCM_Ne),  &
               U(   iNX,iX1,iX2,iX3,iCM_B1),  &
               U(   iNX,iX1,iX2,iX3,iCM_B2),  &
               U(   iNX,iX1,iX2,iX3,iCM_B3),  &
               U(   iNX,iX1,iX2,iX3,iCM_Chi), &
               P(iPM_D ),  &
               P(iPM_V1),  &
               P(iPM_V2),  &
               P(iPM_V3),  &
               P(iPM_E ),  &
               P(iPM_Ne),  &
               P(iPM_B1),  &
               P(iPM_B2),  &
               P(iPM_B3),  &
               P(iPM_Chi), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(   iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               G(   iNX,iX1,iX2,iX3,iGF_Alpha   ), &
               G(   iNX,iX1,iX2,iX3,iGF_Beta_1  ), &
               G(   iNX,iX1,iX2,iX3,iGF_Beta_2  ), &
               G(   iNX,iX1,iX2,iX3,iGF_Beta_3  ) )

      CALL ComputePressureFromPrimitive &
             ( P(iPM_D), P(iPM_E), P(iPM_Ne), Pressure )

      W   = One &
            / SQRT( One & 
                    - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) * P(iPM_V1)**2 &
                    - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) * P(iPM_V2)**2 &
                    - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) * P(iPM_V3)**2 )

      B0u = ( W / G(iNX,iX1,iX2,iX3,iGF_Alpha) ) &
            * ( G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                * P(iPM_V1) * U(iNX,iX1,iX2,iX3,iCM_B1) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                  * P(iPM_V2) * U(iNX,iX1,iX2,iX3,iCM_B2) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                  * P(iPM_V3) * U(iNX,iX1,iX2,iX3,iCM_B3) )

      bSq = ( 1.0_DP / W**2 ) &
            * ( G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) * U(iNX,iX1,iX2,iX3,iCM_B1)**2 &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) * U(iNX,iX1,iX2,iX3,iCM_B2)**2 &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) * U(iNX,iX1,iX2,iX3,iCM_B3)**2 ) &
            + ( G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                * U(iNX,iX1,iX2,iX3,iPM_V1) &
                * U(iNX,iX1,iX2,iX3,iCM_B1) &
              + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                * U(iNX,iX1,iX2,iX3,iPM_V2) &
                * U(iNX,iX1,iX2,iX3,iCM_B2) &
              + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) & 
                * U(iNX,iX1,iX2,iX3,iPM_V3) &
                * U(iNX,iX1,iX2,iX3,iCM_B3) )

      Pressure = Pressure + bSq / 2.0_DP

      PressureTensor(1,1,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S1) * P(iPM_V1) + Pressure ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
          - ( P(iPM_B1)**2 &
              + Two * B0u * P(iPM_B1) * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_1)**2 )
              
      PressureTensor(2,1,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S2) * P(iPM_V1) ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
          - ( P(iPM_B2) * P(iPM_B1) &
              + B0u * P(iPM_B2) * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
              + B0u * P(iPM_B1) * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_2) )
 
      PressureTensor(3,1,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S3) * P(iPM_V1) ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
          - ( P(iPM_B3) * P(iPM_B1) &
              + B0u * P(iPM_B3) * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
              + B0u * P(iPM_B1) * G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 * G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_1) )

      PressureTensor(1,2,iNX,iX1,iX2,iX3) &
        = PressureTensor(2,1,iNX,iX1,iX2,iX3)

      PressureTensor(2,2,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S2) * P(iPM_V2) + Pressure ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
          - ( P(iPM_B2)**2 &
              + Two * B0u * P(iPM_B2) * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_2)**2 )

      PressureTensor(3,2,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S3) * P(iPM_V2) ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
          - ( P(iPM_B3) * P(iPM_B2) &
              + B0u * P(iPM_B3) * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
              + B0u * P(iPM_B2) * G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 * G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_2) )

      PressureTensor(1,3,iNX,iX1,iX2,iX3) &
        = PressureTensor(3,1,iNX,iX1,iX2,iX3)

      PressureTensor(2,3,iNX,iX1,iX2,iX3) &
        = PressureTensor(3,2,iNX,iX1,iX2,iX3)

      PressureTensor(3,3,iNX,iX1,iX2,iX3) &
        = ( U(iNX,iX1,iX2,iX3,iCM_S3) * P(iPM_V3) + Pressure ) &
              / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
          - ( P(iPM_B3)**2 &
              + Two * B0u * P(iPM_B3) * G(iNX,iX1,iX2,iX3,iGF_Beta_3) &
              + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 &
                * G(iNX,iX1,iX2,iX3,iGF_Beta_3)**2 )

      PressureTensorTrace(iNX,iX1,iX2,iX3) &
        =   U(iNX,iX1,iX2,iX3,iCM_S1) * P(iPM_V1) &
          + U(iNX,iX1,iX2,iX3,iCM_S2) * P(iPM_V2) &
          + U(iNX,iX1,iX2,iX3,iCM_S3) * P(iPM_V3) &
          - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
            * ( P(iPM_B1)**2 & 
                + Two * B0u * P(iPM_B1) * G(iNX,iX1,iX2,iX3,iGF_Beta_1)  &
                + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 &
                  * G(iNX,iX1,iX2,iX3,iGF_Beta_1)**2 ) &
          - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
            * ( P(iPM_B2)**2 & 
                + Two * B0u * P(iPM_B2) * G(iNX,iX1,iX2,iX3,iGF_Beta_2)  &
                + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 &
                  * G(iNX,iX1,iX2,iX3,iGF_Beta_2)**2 ) &
          - G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
            * ( P(iPM_B3)**2 & 
                + Two * B0u * P(iPM_B3) * G(iNX,iX1,iX2,iX3,iGF_Beta_3)  &
                + ( G(iNX,iX1,iX2,iX3,iGF_Alpha) * B0u )**2 &
                  * G(iNX,iX1,iX2,iX3,iGF_Beta_3)**2 ) &
          + Three * Pressure

      ! --- X1 increments ---

      ! --- Momentum increment ---

      dU(iNX,iX1,iX2,iX3,iCM_S1) &
        = dU(iNX,iX1,iX2,iX3,iCM_S1) &
            + tau(iNX,iX1,iX2,iX3) &
                * ( G(iNX,iX1,iX2,iX3,iGF_Alpha) &
                      * (   PressureTensor(1,1,iNX,iX1,iX2,iX3) &
                              * G    (iNX,iX1,iX2,iX3,iGF_h_1) &
                              * dGdX1(iNX,iGF_h_1,iX2,iX3,iX1) &
                          + PressureTensor(2,2,iNX,iX1,iX2,iX3) &
                              * G    (iNX,iX1,iX2,iX3,iGF_h_2) &
                              * dGdX1(iNX,iGF_h_2,iX2,iX3,iX1) &
                          + PressureTensor(3,3,iNX,iX1,iX2,iX3) &
                              * G    (iNX,iX1,iX2,iX3,iGF_h_3) &
                              * dGdX1(iNX,iGF_h_3,iX2,iX3,iX1) ) &
                  + U(iNX,iX1,iX2,iX3,iCM_S1) &
                      * dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) &
                  + U(iNX,iX1,iX2,iX3,iCM_S2) &
                      * dGdX1(iNX,iGF_Beta_2,iX2,iX3,iX1) &
                  + U(iNX,iX1,iX2,iX3,iCM_S3) &
                      * dGdX1(iNX,iGF_Beta_3,iX2,iX3,iX1) &
                  - ( U(iNX,iX1,iX2,iX3,iCM_D) + U(iNX,iX1,iX2,iX3,iCM_E) ) &
                      * dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1) )

      ! --- Energy increment ---

      EnergyDensitySourceTerms(1,iNX,iX1,iX2,iX3) &
        = -U(iNX,iX1,iX2,iX3,iCM_S1) &
             / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
             * dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1)

      dU(iNX,iX1,iX2,iX3,iCM_E) &
        = dU(iNX,iX1,iX2,iX3,iCM_E) &
            + tau(iNX,iX1,iX2,iX3) * EnergyDensitySourceTerms(1,iNX,iX1,iX2,iX3)

      IF( UseDivergenceCleaning )THEN

        ! --- Eulerian Magnetic Field increment ---
  
        dU(iNX,iX1,iX2,iX3,iCM_B1) &
          = dU(iNX,iX1,iX2,iX3,iCM_B1) &
              - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_B1) &
                    * dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) ) &
              + U(iNX,iX1,iX2,iX3,iCM_Chi) &
                  * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                  * ( dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1) &
                        / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11 ) ) &
              + U(iNX,iX1,iX2,iX3,iCM_Chi) &
                  * G(iNX,iX1,iX2,iX3,iGF_Alpha) &
                  * ( - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                            / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11)**2 ) &
                        * dGdX1(iNX,iGF_Gm_dd_11,iX2,iX3,iX1) &
                      + ( dGdX1(iNX,iGF_SqrtGm,iX2,iX3,iX1) & 
                            / G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) ) )

        dU(iNX,iX1,iX2,iX3,iCM_B2) &
          = dU(iNX,iX1,iX2,iX3,iCM_B2) &
              - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_B1) &
                    * dGdX1(iNX,iGF_Beta_2,iX2,iX3,iX1) )
 
        dU(iNX,iX1,iX2,iX3,iCM_B3) &
          = dU(iNX,iX1,iX2,iX3,iCM_B3) &
              - ( G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_B1) &
                    * dGdX1(iNX,iGF_Beta_3,iX2,iX3,iX1) )
         
        ! --- Divergence violation field increment ---
  
        dU(iNX,iX1,iX2,iX3,iCM_Chi) &
          = dU(iNX,iX1,iX2,iX3,iCM_Chi) &
              - DampingParameter * U(iNX,iX1,iX2,iX3,iCM_Chi) &
              + ( U(iNX,iX1,iX2,iX3,iCM_B1) / G(iNX,iX1,iX2,iX3,iGF_Alpha) ) &
                * dGdX1(iNX,iGF_Alpha,iX2,iX3,iX1)
  
      END IF

    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, dX2, dX3

  END SUBROUTINE ComputeIncrement_Geometry_Relativistic


  SUBROUTINE ComputeIncrement_Gravity &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

   INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    CALL ComputeIncrement_Gravity_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

  END SUBROUTINE ComputeIncrement_Gravity


  SUBROUTINE ComputeIncrement_Gravity_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
  END SUBROUTINE ComputeIncrement_Gravity_Relativistic


  SUBROUTINE InitializeIncrement_MHD & 
               ( iX_B0, iX_E0, iX_B1, iX_E1 )

    INTEGER,  INTENT(in)           :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    nX        = iX_E0 - iX_B0 + 1 ! Number of Elements per Dimension
    nX_K      = PRODUCT( nX )     ! Number of Elements in Position Space
    nNodesX_K = nDOFX * nX_K      ! Number of Nodes in Elements

    nX_X1 = nX + [1,0,0] ! Number of X1 Faces per Dimension
    nX_X2 = nX + [0,1,0] ! Number of X2 Faces per Dimension
    nX_X3 = nX + [0,0,1] ! Number of X3 Faces per Dimension

    nX1_X = PRODUCT( nX_X1 ) ! Number of X1 Faces
    nX2_X = PRODUCT( nX_X2 ) ! Number of X2 Faces
    nX3_X = PRODUCT( nX_X3 ) ! Number of X3 Faces

    nNodesX_X1 = nDOFX_X1 * nX1_X ! Number of Nodes on X1 Faces
    nNodesX_X2 = nDOFX_X2 * nX2_X ! Number of Nodes on X2 Faces
    nNodesX_X3 = nDOFX_X3 * nX3_X ! Number of Nodes on X3 Faces

    ALLOCATE( pD_K  (nNodesX_K) )
    ALLOCATE( pV1_K (nNodesX_K) )
    ALLOCATE( pV2_K (nNodesX_K) )
    ALLOCATE( pV3_K (nNodesX_K) )
    ALLOCATE( pE_K  (nNodesX_K) )
    ALLOCATE( pNe_K (nNodesX_K) )
    ALLOCATE( pB1_K (nNodesX_K) )
    ALLOCATE( pB2_K (nNodesX_K) )
    ALLOCATE( pB3_K (nNodesX_K) )
    ALLOCATE( pChi_K(nNodesX_K) )

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE InitializeIncrement_MHD


  SUBROUTINE FinalizeIncrement_MHD

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    END ASSOCIATE ! dX1, etc.

    DEALLOCATE( pD_K, pV1_K, pV2_K, pV3_K, pE_K, pNe_K, &
                pB1_K, pB2_K, pB3_K, pChi_K )

  END SUBROUTINE FinalizeIncrement_MHD


  SUBROUTINE InitializeIncrement_Divergence &
    ( iXP_B0, iXP_E0, nDOFX_X, G_K, G_F, uCM_K, uCM_L, uCM_R )

    INTEGER, INTENT(in) :: iXP_B0(3), iXP_E0(3) ! Permuted limits
    INTEGER, INTENT(in) :: nDOFX_X              ! nDOFX_X1, ...

    ! --- Geometry Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      G_K (1:nDOFX, &
           iXP_B0(1)  :iXP_E0(1)  , &
           iXP_B0(2)  :iXP_E0(2)  , &
           iXP_B0(3)-1:iXP_E0(3)+1, &
           1:nGF), &
      G_F (1:nDOFX_X, &
           iXP_B0(1)  :iXP_E0(1)  , &
           iXP_B0(2)  :iXP_E0(2)  , &
           iXP_B0(3)  :iXP_E0(3)+1, &
           1:nGF)

    ! --- Conserved Fluid Fields ---

    REAL(DP), TARGET, INTENT(in) :: &
      uCM_K(1:nDOFX, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)-1:iXP_E0(3)+1, &
            1:nCM), &
      uCM_L(1:nDOFX_X, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)  :iXP_E0(3)+1, &
            1:nCM), &
      uCM_R(1:nDOFX_X, &
            iXP_B0(1)  :iXP_E0(1)  , &
            iXP_B0(2)  :iXP_E0(2)  , &
            iXP_B0(3)  :iXP_E0(3)+1, &
            1:nCM)

    INTEGER :: iNX, iX1, iX2, iX3, iCM
    INTEGER :: nXP(3), nXP_X(3), nX_X, nNodesX_X, iX_F, iX_V

    nXP       = iXP_E0 - iXP_B0 + 1
    nXP_X     = nXP + [0,0,1]
    nX_X      = PRODUCT( nXP_X )
    nNodesX_X = nDOFX_X * nX_X

    Gm_dd_11_K(1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Gm_dd_11)
    Gm_dd_22_K(1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Gm_dd_22)
    Gm_dd_33_K(1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Gm_dd_33)
    SqrtGm_K  (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_SqrtGm  )
    Beta_1_K  (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Beta_1  )
    Beta_2_K  (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Beta_2  )
    Beta_3_K  (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Beta_3  )
    Alpha_K   (1:nNodesX_K) => G_K(:,:,:,iXP_B0(3):iXP_E0(3),iGF_Alpha   )

    Gm_dd_11_F(1:nNodesX_X) => G_F(:,:,:,:,iGF_Gm_dd_11)
    Gm_dd_22_F(1:nNodesX_X) => G_F(:,:,:,:,iGF_Gm_dd_22)
    Gm_dd_33_F(1:nNodesX_X) => G_F(:,:,:,:,iGF_Gm_dd_33)
    SqrtGm_F  (1:nNodesX_X) => G_F(:,:,:,:,iGF_SqrtGm  )
    Beta_1_F  (1:nNodesX_X) => G_F(:,:,:,:,iGF_Beta_1  )
    Beta_2_F  (1:nNodesX_X) => G_F(:,:,:,:,iGF_Beta_2  )
    Beta_3_F  (1:nNodesX_X) => G_F(:,:,:,:,iGF_Beta_3  )
    Alpha_F   (1:nNodesX_X) => G_F(:,:,:,:,iGF_Alpha   )

    uD_K  (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_D )
    uS1_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_S1)
    uS2_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_S2)
    uS3_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_S3)
    uE_K  (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_E )
    uNe_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_Ne)
    uB1_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_B1)
    uB2_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_B2)
    uB3_K (1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_B3)
    uChi_K(1:nNodesX_K) => uCM_K(:,:,:,iXP_B0(3):iXP_E0(3),iCM_Chi)

    uD_L  (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_D )
    uS1_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_S1)
    uS2_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_S2)
    uS3_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_S3)
    uE_L  (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_E )
    uNe_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_Ne)
    uB1_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_B1)
    uB2_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_B2)
    uB3_L (1:nNodesX_X) => uCM_L(:,:,:,:,iCM_B3)
    uChi_L(1:nNodesX_X) => uCM_L(:,:,:,:,iCM_Chi)

    uD_R  (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_D )
    uS1_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_S1)
    uS2_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_S2)
    uS3_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_S3)
    uE_R  (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_E )
    uNe_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_Ne)
    uB1_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_B1)
    uB2_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_B2)
    uB3_R (1:nNodesX_X) => uCM_R(:,:,:,:,iCM_B3)
    uChi_R(1:nNodesX_X) => uCM_R(:,:,:,:,iCM_Chi)

    ALLOCATE( pD_L  (nNodesX_X) )
    ALLOCATE( pV1_L (nNodesX_X) )
    ALLOCATE( pV2_L (nNodesX_X) )
    ALLOCATE( pV3_L (nNodesX_X) )
    ALLOCATE( pE_L  (nNodesX_X) )
    ALLOCATE( pNe_L (nNodesX_X) )
    ALLOCATE( pB1_L (nNodesX_X) )
    ALLOCATE( pB2_L (nNodesX_X) )
    ALLOCATE( pB3_L (nNodesX_X) )
    ALLOCATE( pChi_L(nNodesX_X) )

    ALLOCATE( pD_R  (nNodesX_X) )
    ALLOCATE( pV1_R (nNodesX_X) )
    ALLOCATE( pV2_R (nNodesX_X) )
    ALLOCATE( pV3_R (nNodesX_X) )
    ALLOCATE( pE_R  (nNodesX_X) )
    ALLOCATE( pNe_R (nNodesX_X) )
    ALLOCATE( pB1_R (nNodesX_X) )
    ALLOCATE( pB2_R (nNodesX_X) )
    ALLOCATE( pB3_R (nNodesX_X) )
    ALLOCATE( pChi_R(nNodesX_X) )

    ALLOCATE( IndexTableX_F(4,nNodesX_X) )
    ALLOCATE( IndexTableX_V(4,nNodesX_K) )

    DO iX3 = iXP_B0(3), iXP_E0(3) + 1
    DO iX2 = iXP_B0(2), iXP_E0(2)
    DO iX1 = iXP_B0(1), iXP_E0(1)
    DO iNX = 1, nDOFX_X

      iX_F = iNX &
               + ( iX1 - iXP_B0(1) ) * nDOFX_X &
               + ( iX2 - iXP_B0(2) ) * nDOFX_X * nXP_X(1) &
               + ( iX3 - iXP_B0(3) ) * nDOFX_X * nXP_X(1) * nXP_X(2)

      IndexTableX_F(:,iX_F) = [ iNX, iX1, iX2, iX3 ]

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iXP_B0(3), iXP_E0(3)
    DO iX2 = iXP_B0(2), iXP_E0(2)
    DO iX1 = iXP_B0(1), iXP_E0(1)
    DO iNX = 1, nDOFX

      iX_V = iNX &
               + ( iX1 - iXP_B0(1) ) * nDOFX &
               + ( iX2 - iXP_B0(2) ) * nDOFX * nXP(1) &
               + ( iX3 - iXP_B0(3) ) * nDOFX * nXP(1) * nXP(2)

      IndexTableX_V(:,iX_V) = [ iNX, iX1, iX2, iX3 ]

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeIncrement_Divergence


  SUBROUTINE FinalizeIncrement_Divergence

    DEALLOCATE( pD_L, pV1_L, pV2_L, pV3_L, pE_L, pNe_L, & 
                pB1_L, pB2_L, pB3_L, pChi_L )
    DEALLOCATE( pD_R, pV1_R, pV2_R, pV3_R, pE_R, pNe_R, &
                pB1_R, pB2_R, pB3_R, pChi_R )
    DEALLOCATE( IndexTableX_F, IndexTableX_V )

    NULLIFY( Gm_dd_11_K, Gm_dd_22_K, Gm_dd_33_K, SqrtGm_K )
    NULLIFY( Gm_dd_11_F, Gm_dd_22_F, Gm_dd_33_F, SqrtGm_F )
    NULLIFY(   Beta_1_K,   Beta_2_K,   Beta_3_K,  Alpha_K )
    NULLIFY(   Beta_1_F,   Beta_2_F,   Beta_3_F,  Alpha_F )
    NULLIFY( uD_K, uS1_K, uS2_K, uS3_K, uE_K, uNe_K, &
             uB1_K, uB2_K, uB3_K, uChi_K )
    NULLIFY( uD_L, uS1_L, uS2_L, uS3_L, uE_L, uNe_L, &
             uB1_L, uB2_L, uB3_L, uChi_L )
    NULLIFY( uD_R, uS1_R, uS2_R, uS3_R, uE_R, uNe_R, &
             uB1_R, uB2_R, uB3_R, uChi_R )

  END SUBROUTINE FinalizeIncrement_Divergence


END MODULE MHD_DiscretizationModule_Relativistic
