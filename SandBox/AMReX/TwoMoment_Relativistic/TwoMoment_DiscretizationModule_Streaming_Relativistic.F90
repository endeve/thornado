MODULE TwoMoment_DiscretizationModule_Streaming_Relativistic

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, nDOFX_X2, nDOFX_X3, &
    WeightsX_q, &
    WeightsX_X1, &
    WeightsX_X2, &
    WeightsX_X3
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, LX_X1_Dn, LX_X1_Up, &
    dLXdX2_q, LX_X2_Dn, LX_X2_Up, &
    dLXdX3_q, LX_X3_Dn, LX_X3_Up
  USE ReferenceElementModule, ONLY: &
    nDOF_E, &
    nDOF_X1, nDOF_X2, nDOF_X3, &
    Weights_q, &
    Weights_E, &
    Weights_X1, &
    Weights_X2, &
    Weights_X3
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_E_Dn,  L_E_Up, &
    dLdE_q, &
    L_X1_Dn, L_X1_Up, &
    dLdX1_q, &
    L_X2_Dn, L_X2_Up, &
    dLdX2_q, &
    L_X3_Dn, L_X3_Up, &
    dLdX3_q
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE, &
    iGE_Ep2, &
    iGE_Ep3
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
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE TwoMoment_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_TwoMoment, &
    Flux_X1, &
    ComputeEddingtonTensorComponents_ud, &
    NumericalFlux_LLF





  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit

CONTAINS

  SUBROUTINE ComputeIncrement_TwoMoment_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R, Verbose_Option )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                   iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                   iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R (1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS
    LOGICAL :: Verbose

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option    
   
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )


    CALL ApplyBoundaryConditions_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U_F )

    CALL ApplyBoundaryConditions_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R )

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO



    CALL ComputeIncrement_Divergence_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R, &
             Verbose_Option = Verbose )
    ! --- Multiply Inverse Mass Matrix ---

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
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

    END ASSOCIATE ! dZ1, etc.


  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit

  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R, Verbose_Option )

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
           1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeZ, iNodeE, iNodeX
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF, iPF
    INTEGER  :: nZ(4), nZ_X1(4), nV_X1, nV, nX_X1
    INTEGER  :: nIterations
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF)
    REAL(DP) :: uPR_L(nPR), Flux_L(nCR)
    REAL(DP) :: uPR_R(nPR), Flux_R(nCR)
    REAL(DP) :: uPR_K(nPR), Flux_K(nCR)
    REAL(DP) :: &
      GX_K(nDOFX,nGF, &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B1(2):iZ_E1(2))
    REAL(DP) :: &
      GX_F(nDOFX_X1,nGF, &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      uCF_K(nDOFX,nCF, &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B1(2):iZ_E1(2))
    REAL(DP) :: &
      uCF_L(nDOFX_X1,nCF, &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      uCF_R(nDOFX_X1,nCF, &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      uPF_K(nDOFX,nPF, &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2))
    REAL(DP) :: &
      V_u(3,nDOFX_X1, &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      uCR_K(nDOFZ,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B1(2):iZ_E1(2))
    REAL(DP) :: &
      uCR_L(nDOF_X1,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      uCR_R(nDOF_X1,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      NumericalFlux &
        (nDOF_X1,nCR, &
         iZ_B0(1):iZ_E0(1), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies, &
         iZ_B0(2):iZ_E1(2))
    REAL(DP) :: &
      Flux_q(nDOFZ,nCR, &
             iZ_B0(1):iZ_E0(1), &
             iZ_B0(3):iZ_E0(3), &
             iZ_B0(4):iZ_E0(4), &
             nSpecies, &
             iZ_B0(2):iZ_E0(2))
    REAL(DP) :: &
      dU_X1(nDOFZ,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B0(2):iZ_E0(2))
    LOGICAL :: Verbose

    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
      PRINT*, "      ComputeIncrement_Divergence_X1"
    END IF

    nZ    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nZ_X1 = nZ + [0,1,0,0]    ! Number of X3 Faces per Phase Space Dimension
    nV    = nCR * nSpecies * PRODUCT( nZ )
    nV_X1 = nCR * nSpecies * PRODUCT( nZ_X1 )
    nX_X1 = PRODUCT( nZ_X1(2:4) ) ! Number of X1 Faces in Position Space

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )
    ! --- Permute Geometry Fields ---

    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iGF    = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iZ3,iZ4,iZ2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    ! --- Face States (Average of Left and Right States) ---


    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX_X1*nGF, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             GX_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1), nDOFX, Zero, &
             GX_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX_X1*nGF, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             GX_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX, Half, &
             GX_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )




    ! --- Permute Fluid Fields ---

    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iCF = 1, nCF
      DO iNodeX = 1, nDOFX
        
        uCF_K(iNodeX,iCF,iZ3,iZ4,iZ2) = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)
      END DO
      END DO

    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX_X1*nCF, nDOFX, One, LX_X1_Up, nDOFX_X1, &
             uCF_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1), nDOFX, Zero, &
             uCF_L(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX_X1*nCF, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
             uCF_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX, Zero, &
             uCF_R(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )

    ! --- Compute Face Velocity Components --




    DO iZ2 = iZ_B0(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

    END DO
    END DO
    END DO

    DO iZ2  = iZ_B0(2), iZ_E1(2)
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

!this is where the issue is

        ! --- Left State ---
        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_L(iNodeX,iCF_D ,iZ3,iZ4,iZ2), &
                 uCF_L(iNodeX,iCF_S1,iZ3,iZ4,iZ2), &
                 uCF_L(iNodeX,iCF_S2,iZ3,iZ4,iZ2), &
                 uCF_L(iNodeX,iCF_S3,iZ3,iZ4,iZ2), &
                 uCF_L(iNodeX,iCF_E ,iZ3,iZ4,iZ2), &
                 uCF_L(iNodeX,iCF_Ne,iZ3,iZ4,iZ2), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )
        ! --- Right State ---
        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_R(iNodeX,iCF_D ,iZ3,iZ4,iZ2), &
                 uCF_R(iNodeX,iCF_S1,iZ3,iZ4,iZ2), &
                 uCF_R(iNodeX,iCF_S2,iZ3,iZ4,iZ2), &
                 uCF_R(iNodeX,iCF_S3,iZ3,iZ4,iZ2), &
                 uCF_R(iNodeX,iCF_E ,iZ3,iZ4,iZ2), &
                 uCF_R(iNodeX,iCF_Ne,iZ3,iZ4,iZ2), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

        V_u(1:3,iNodeX,iZ3,iZ4,iZ2) &
          = FaceVelocity_X1 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )
      END DO

    END DO
    END DO
    END DO
    ! --- Permute Radiation Fields ---

    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iCR = 1, nCR

      DO iNodeZ = 1, nDOFZ

        uCR_K(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Radiation Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X1, nV_X1, nDOFZ, One, L_X1_Up, nDOF_X1, &
             uCR_K(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)-1), nDOFZ, Zero, &
             uCR_L(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), nDOF_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X1, nV_X1, nDOFZ, One, L_X1_Dn, nDOF_X1, &
             uCR_K(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), nDOFZ, Zero, &
             uCR_R(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), nDOF_X1 )


    ! --- Numerical Flux ---

    DO iZ2 = iZ_B0(2), iZ_E1(2)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)
      DO iNodeX = 1, nDOFX_X1
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE
        ! --- Left State Primitive ---
        CALL ComputePrimitive_TwoMoment &
               ( uCR_L(iNodeZ,iCR_N ,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_L(iNodeZ,iCR_G1,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_L(iNodeZ,iCR_G2,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_L(iNodeZ,iCR_G3,iZ1,iZ3,iZ4,iS,iZ2), &
                 uPR_L(iPR_D ), uPR_L(iPR_I1), &
                 uPR_L(iPR_I2), uPR_L(iPR_I3), &
                 V_u(1,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(2,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(3,iNodeX,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 0.0_DP, 0.0_DP, 0.0_DP,                &
                 GX_F(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2), &
                 nIterations )
        ! --- Left State Flux ---

        Flux_L &
          = Flux_X1 &
              ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                uPR_L(iPR_I2), uPR_L(iPR_I3), &
                V_u(1,iNodeX,iZ3,iZ4,iZ2), &
                V_u(2,iNodeX,iZ3,iZ4,iZ2), &
                V_u(3,iNodeX,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), & 
                GX_F(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2) )
        ! --- Right State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_R(iNodeZ,iCR_N ,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_R(iNodeZ,iCR_G1,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_R(iNodeZ,iCR_G2,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_R(iNodeZ,iCR_G3,iZ1,iZ3,iZ4,iS,iZ2), &
                 uPR_R(iPR_D ), uPR_R(iPR_I1), &
                 uPR_R(iPR_I2), uPR_R(iPR_I3), &
                 V_u(1,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(2,iNodeX,iZ3,iZ4,iZ2), &
                 V_u(3,iNodeX,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 0.0_DP,0.0_DP,0.0_DP,                  &
                 GX_F(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                 GX_F(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2), &
                 nIterations )
        ! --- Right State Flux ---

        Flux_R &
          = Flux_X1 &
              ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                uPR_R(iPR_I2), uPR_R(iPR_I3), &
                V_u(1,iNodeX,iZ3,iZ4,iZ2), &
                V_u(2,iNodeX,iZ3,iZ4,iZ2), &
                V_u(3,iNodeX,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                GX_F(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2) )
        ! --- Numerical Flux ---

        DO iCR = 1, nCR

          NumericalFlux(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = NumericalFlux_LLF &
                ( uCR_L(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2), &
                  uCR_R(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2), &
                  Flux_L(iCR), Flux_R(iCR), One )

          NumericalFlux(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
                * Weights_X1(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
                * GX_F(iNodeX,iGF_SqrtGm,iZ3,iZ4,iZ2) &
                * NumericalFlux(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2)

        END DO
      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_X1, + One, L_X1_Dn, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), &
             nDOF_X1, Zero, dU_X1, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_X1, - One, L_X1_Up, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)+1), &
             nDOF_X1, One,  dU_X1, nDOFZ )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    ! --- Compute Primitive Fluid in Spatial Elements ---

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF_K(iNodeX,iCF_D       ,iZ3,iZ4,iZ2), &
                 uCF_K(iNodeX,iCF_S1      ,iZ3,iZ4,iZ2), &
                 uCF_K(iNodeX,iCF_S2      ,iZ3,iZ4,iZ2), &
                 uCF_K(iNodeX,iCF_S3      ,iZ3,iZ4,iZ2), &
                 uCF_K(iNodeX,iCF_E       ,iZ3,iZ4,iZ2), &
                 uCF_K(iNodeX,iCF_Ne      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_D       ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_V1      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_V2      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_V3      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_E       ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_Ne      ,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

      END DO

    END DO
    END DO
    END DO

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNodeZ,iCR_N ,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_K(iNodeZ,iCR_G1,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_K(iNodeZ,iCR_G2,iZ1,iZ3,iZ4,iS,iZ2), &
                 uCR_K(iNodeZ,iCR_G3,iZ1,iZ3,iZ4,iS,iZ2), &
                 uPR_K(iPR_D ), uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), uPR_K(iPR_I3), &
                 uPF_K(iNodeX,iPF_V1      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_V2      ,iZ3,iZ4,iZ2), &
                 uPF_K(iNodeX,iPF_V3      ,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 0.0_DP, 0.0_DP, 0.0_DP,                &
                 GX_K(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                 GX_K(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                 GX_K(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                 GX_K(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2), &
                 nIterations )

        Flux_K &
          = Flux_X1 &
              ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                uPR_K(iPR_I2), uPR_K(iPR_I3), &
                uPF_K(iNodeX,iPF_V1      ,iZ3,iZ4,iZ2), &
                uPF_K(iNodeX,iPF_V2      ,iZ3,iZ4,iZ2), &
                uPF_K(iNodeX,iPF_V3      ,iZ3,iZ4,iZ2), &
                GX_K (iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                GX_K (iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                GX_K (iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                GX_K(iNodeX,iGF_Alpha   ,iZ3,iZ4,iZ2), &
                GX_K(iNodeX,iGF_Beta_1  ,iZ3,iZ4,iZ2), &
                GX_K(iNodeX,iGF_Beta_2  ,iZ3,iZ4,iZ2), &
                GX_K(iNodeX,iGF_Beta_3  ,iZ3,iZ4,iZ2) )

        DO iCR = 1, nCR

          Flux_q(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) &
                * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
                * GX_K(iNodeX,iGF_SqrtGm,iZ3,iZ4,iZ2) &
                * Flux_K(iCR)

        END DO

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOFZ, One, dLdX1_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_X1, nDOFZ )


    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_X1(iNodeZ,iCR,iZ1,iZ3,iZ4,iS,iZ2)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X1




  FUNCTION FaceVelocity_X1 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R )

    REAL(DP), INTENT(in) :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in) :: V1_R, V2_R, V3_R
    REAL(DP)             :: FaceVelocity_X1(1:3)

    ! --- Average Left and Right States ---

    FaceVelocity_X1(1) = Half * ( V1_L + V1_R )
    FaceVelocity_X1(2) = Half * ( V2_L + V2_R )
    FaceVelocity_X1(3) = Half * ( V3_L + V3_R )

    RETURN
  END FUNCTION FaceVelocity_X1


END MODULE TwoMoment_DiscretizationModule_Streaming_Relativistic
