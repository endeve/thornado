MODULE TwoMoment_DiscretizationModule_Streaming_OrderV

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Explicit, &
    Timer_Ex_Permute, &
    Timer_Ex_Interpolate, &
    Timer_Ex_Flux, &
    Timer_Ex_Increment, &
    Timer_Ex_In, &
    Timer_Ex_Div, &
    Timer_Ex_Geometry
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
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputePrimitive_TwoMoment, &
    Flux_E, &
    Flux_X1, &
    Flux_X2, &
    Flux_X3, &
    ComputeEddingtonTensorComponents_ud, &
    NumericalFlux_LLF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit

CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

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

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iNodeE, iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS

    CALL TimersStart( Timer_Explicit )

    PRINT*, "      ComputeIncrement_TwoMoment_Explicit (Start)"

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

    CALL TimersStart( Timer_Ex_In ) ! --- Hijacked Timer

    CALL ComputeIncrement_Divergence_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    CALL TimersStop( Timer_Ex_In )

    CALL TimersStart( Timer_Ex_Div  ) ! --- Hijacked Timer

    CALL ComputeIncrement_Divergence_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    CALL TimersStop( Timer_Ex_Div )

    CALL TimersStart( Timer_Ex_Geometry  ) ! --- Hijacked Timer

    CALL ComputeIncrement_Divergence_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    CALL TimersStop( Timer_Ex_Geometry )

    CALL ComputeIncrement_ObserverCorrections &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

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

    PRINT*, "      ComputeIncrement_TwoMoment_Explicit (End)"

    CALL TimersStop( Timer_Explicit )

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

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

    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN

    PRINT*, "      ComputeIncrement_Divergence_X1"

    nZ    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nZ_X1 = nZ + [0,1,0,0]    ! Number of X3 Faces per Phase Space Dimension
    nV    = nCR * nSpecies * PRODUCT( nZ )
    nV_X1 = nCR * nSpecies * PRODUCT( nZ_X1 )
    nX_X1 = PRODUCT( nZ_X1(2:4) ) ! Number of X1 Faces in Position Space

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

    CALL TimersStart( Timer_Ex_Permute )

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

    CALL TimersStop( Timer_Ex_Permute )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX_X1*nGF, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             GX_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1), nDOFX, Zero, &
             GX_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX_X1*nGF, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             GX_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX, Half, &
             GX_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1 )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Recompute Geometry from Scale Factors ---

    DO iZ2  = iZ_B0(2), iZ_E1(2)
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2) &
          = MAX( GX_F(iNodeX,iGF_h_1,iZ3,iZ4,iZ2)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2) &
          = MAX( GX_F(iNodeX,iGF_h_2,iZ3,iZ4,iZ2)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) &
          = MAX( GX_F(iNodeX,iGF_h_3,iZ3,iZ4,iZ2)**2, SqrtTiny )
        GX_F(iNodeX,iGF_SqrtGm,iZ3,iZ4,iZ2) &
          = SQRT(   GX_F(iNodeX,iGF_Gm_dd_11,iZ3,iZ4,iZ2) &
                  * GX_F(iNodeX,iGF_Gm_dd_22,iZ3,iZ4,iZ2) &
                  * GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Permute )

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

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Interpolate )

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

    CALL TimersStop( Timer_Ex_Interpolate )

    ! --- Compute Face Velocity Components ---

    DO iZ2  = iZ_B0(2), iZ_E1(2)
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        ! --- Left State ---

        CALL ComputePrimitive_Euler_NonRelativistic &
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

        CALL ComputePrimitive_Euler_NonRelativistic &
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

    CALL TimersStart( Timer_Ex_Permute )

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

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Interpolate )

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

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Flux )

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
                GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

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
                GX_F(iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

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

    CALL TimersStop( Timer_Ex_Flux )

    CALL TimersStart( Timer_Ex_Interpolate )

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

    CALL TimersStop( Timer_Ex_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    ! --- Compute Primitive Fluid in Spatial Elements ---

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_NonRelativistic &
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

    CALL TimersStart( Timer_Ex_Flux )

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
                GX_K (iNodeX,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

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

    CALL TimersStop( Timer_Ex_Flux )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOFZ, One, dLdX1_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_X1, nDOFZ )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Increment )

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

    CALL TimersStop( Timer_Ex_Increment )

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X1


  SUBROUTINE ComputeIncrement_Divergence_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

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

    INTEGER  :: iNode, iNodeZ, iNodeE, iNodeX
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER  :: nZ(4), nZ_X2(4), nV_X2, nV, nX_X2
    INTEGER  :: nIterations
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF)
    REAL(DP) :: uPR_L(nPR), Flux_L(nCR)
    REAL(DP) :: uPR_R(nPR), Flux_R(nCR)
    REAL(DP) :: uPR_K(nPR), Flux_K(nCR)
    REAL(DP) :: &
      GX_K(nDOFX,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B1(3):iZ_E1(3))
    REAL(DP) :: &
      GX_F(nDOFX_X2,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      uCF_K(nDOFX,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B1(3):iZ_E1(3))
    REAL(DP) :: &
      uCF_L(nDOFX_X2,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      uCF_R(nDOFX_X2,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      uPF_K(nDOFX,nPF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E0(3))
    REAL(DP) :: &
      V_u(3,nDOFX_X2, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(4):iZ_E0(4), &
          iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      uCR_K(nDOFZ,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B1(3):iZ_E1(3))
    REAL(DP) :: &
      uCR_L(nDOF_X2,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      uCR_R(nDOF_X2,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      NumericalFlux &
        (nDOF_X2,nCR, &
         iZ_B0(1):iZ_E0(1), &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies, &
         iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      Flux_q(nDOFZ,nCR, &
             iZ_B0(1):iZ_E0(1), &
             iZ_B0(2):iZ_E0(2), &
             iZ_B0(4):iZ_E0(4), &
             nSpecies, &
             iZ_B0(3):iZ_E0(3))
    REAL(DP) :: &
      dU_X2(nDOFZ,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B0(3):iZ_E0(3))

    IF( iZ_E0(3) .EQ. iZ_B0(3) ) RETURN

    PRINT*, "      ComputeIncrement_Divergence_X2"

    nZ    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nZ_X2 = nZ + [0,0,1,0]    ! Number of X2 Faces per Phase Space Dimension
    nV    = nCR * nSpecies * PRODUCT( nZ )
    nV_X2 = nCR * nSpecies * PRODUCT( nZ_X2 )
    nX_X2 = PRODUCT( nZ_X2(2:4) ) ! Number of X2 Faces in Position Space

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ4 => MeshX(3) % Width )

    ! --- Permute Geometry Fields ---

    CALL TimersStart( Timer_Ex_Permute )

    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iGF    = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iZ2,iZ4,iZ3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX_X2*nGF, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             GX_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1), nDOFX, Zero, &
             GX_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX_X2*nGF, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             GX_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX, Half, &
             GX_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Recompute Geometry from Scale Factors ---

    DO iZ3  = iZ_B0(3), iZ_E1(3)
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X2

        GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3) &
          = MAX( GX_F(iNodeX,iGF_h_1,iZ2,iZ4,iZ3)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3) &
          = MAX( GX_F(iNodeX,iGF_h_2,iZ2,iZ4,iZ3)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) &
          = MAX( GX_F(iNodeX,iGF_h_3,iZ2,iZ4,iZ3)**2, SqrtTiny )
        GX_F(iNodeX,iGF_SqrtGm,iZ2,iZ4,iZ3) &
          = SQRT(   GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3) &
                  * GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3) &
                  * GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Permute Fluid Fields ---

    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iCF    = 1, nCF
      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iCF,iZ2,iZ4,iZ3) = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Fluid Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX_X2*nCF, nDOFX, One, LX_X2_Up, nDOFX_X2, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1), nDOFX, Zero, &
             uCF_L(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX_X2*nCF, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX, Zero, &
             uCF_R(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )

    CALL TimersStop( Timer_Ex_Interpolate )

    ! --- Compute Face Velocity Components ---

    DO iZ3 = iZ_B0(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNode = 1, nDOFX_X2

        ! --- Left State ---

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_L(iNode,iCF_D ,iZ2,iZ4,iZ3), &
                 uCF_L(iNode,iCF_S1,iZ2,iZ4,iZ3), &
                 uCF_L(iNode,iCF_S2,iZ2,iZ4,iZ3), &
                 uCF_L(iNode,iCF_S3,iZ2,iZ4,iZ3), &
                 uCF_L(iNode,iCF_E ,iZ2,iZ4,iZ3), &
                 uCF_L(iNode,iCF_Ne,iZ2,iZ4,iZ3), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_F(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_F(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

        ! --- Right State ---

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_R(iNode,iCF_D ,iZ2,iZ4,iZ3), &
                 uCF_R(iNode,iCF_S1,iZ2,iZ4,iZ3), &
                 uCF_R(iNode,iCF_S2,iZ2,iZ4,iZ3), &
                 uCF_R(iNode,iCF_S3,iZ2,iZ4,iZ3), &
                 uCF_R(iNode,iCF_E ,iZ2,iZ4,iZ3), &
                 uCF_R(iNode,iCF_Ne,iZ2,iZ4,iZ3), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F(iNode,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_F(iNode,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_F(iNode,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

        ! --- Face Velocity ---

        V_u(1:3,iNode,iZ2,iZ4,iZ3) &
          = FaceVelocity_X2 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Permute Radiation Fields ---

    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iCR = 1, nCR

      DO iNodeZ = 1, nDOFZ

        uCR_K(iNodeZ,iCR,iZ1,iZ2,iZ4,iS,iZ3) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Radiation Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X2, nV_X2, nDOFZ, One, L_X2_Up, nDOF_X2, &
             uCR_K(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)-1), nDOFZ, Zero, &
             uCR_L(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)  ), nDOF_X2 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X2, nV_X2, nDOFZ, One, L_X2_Dn, nDOF_X2, &
             uCR_K(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)  ), nDOFZ, Zero, &
             uCR_R(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)  ), nDOF_X2 )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Flux )

    ! --- Numerical Flux ---

    DO iZ3 = iZ_B0(3), iZ_E1(3)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX_X2
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        ! --- Left State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_L(iNodeZ,iCR_N ,iZ1,iZ2,iZ4,iS,iZ3), &
                 uCR_L(iNodeZ,iCR_G1,iZ1,iZ2,iZ4,iS,iZ3), &
                 uCR_L(iNodeZ,iCR_G2,iZ1,iZ2,iZ4,iS,iZ3), &
                 uCR_L(iNodeZ,iCR_G3,iZ1,iZ2,iZ4,iS,iZ3), &
                 uPR_L(iPR_D ), uPR_L(iPR_I1), &
                 uPR_L(iPR_I2), uPR_L(iPR_I3), &
                 V_u(1,iNodeX,iZ2,iZ4,iZ3), &
                 V_u(2,iNodeX,iZ2,iZ4,iZ3), &
                 V_u(3,iNodeX,iZ2,iZ4,iZ3), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3), &
                 nIterations )

        ! --- Left State Flux ---

        Flux_L &
          = Flux_X2 &
              ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                uPR_L(iPR_I2), uPR_L(iPR_I3), &
                V_u(1,iNodeX,iZ2,iZ4,iZ3), &
                V_u(2,iNodeX,iZ2,iZ4,iZ3), &
                V_u(3,iNodeX,iZ2,iZ4,iZ3), &
                GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

        ! --- Right State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_R(iNodeZ,iCR_N ,iZ1,iZ2,iZ4,iS,iZ3), &
                 uCR_R(iNodeZ,iCR_G1,iZ1,iZ2,iZ4,iS,iZ3), &
                 uCR_R(iNodeZ,iCR_G2,iZ1,iZ2,iZ4,iS,iZ3), &
                 uCR_R(iNodeZ,iCR_G3,iZ1,iZ2,iZ4,iS,iZ3), &
                 uPR_R(iPR_D ), uPR_R(iPR_I1), &
                 uPR_R(iPR_I2), uPR_R(iPR_I3), &
                 V_u(1,iNodeX,iZ2,iZ4,iZ3), &
                 V_u(2,iNodeX,iZ2,iZ4,iZ3), &
                 V_u(3,iNodeX,iZ2,iZ4,iZ3), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3), &
                 nIterations )

        ! --- Right State Flux ---

        Flux_R &
          = Flux_X2 &
              ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                uPR_R(iPR_I2), uPR_R(iPR_I3), &
                V_u(1,iNodeX,iZ2,iZ4,iZ3), &
                V_u(2,iNodeX,iZ2,iZ4,iZ3), &
                V_u(3,iNodeX,iZ2,iZ4,iZ3), &
                GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

        ! --- Numerical Flux ---

        DO iCR = 1, nCR

          NumericalFlux(iNodeZ,iCR,iZ1,iZ2,iZ4,iS,iZ3) &
            = NumericalFlux_LLF &
                ( uCR_L(iNodeZ,iCR,iZ1,iZ2,iZ4,iS,iZ3), &
                  uCR_R(iNodeZ,iCR,iZ1,iZ2,iZ4,iS,iZ3), &
                  Flux_L(iCR), Flux_R(iCR), One )

          NumericalFlux(iNodeZ,iCR,iZ1,iZ2,iZ4,iS,iZ3) &
            = dZ1(iZ1) * dZ2(iZ2) * dZ4(iZ4) &
                * Weights_X2(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
                * GX_F(iNodeX,iGF_SqrtGm,iZ2,iZ4,iZ3) &
                * NumericalFlux(iNodeZ,iCR,iZ1,iZ2,iZ4,iS,iZ3)

        END DO

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_X2, + One, L_X2_Dn, nDOF_X2, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)  ), &
             nDOF_X2, Zero, dU_X2, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_X2, - One, L_X2_Up, nDOF_X2, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(4),1,iZ_B0(3)+1), &
             nDOF_X2, One,  dU_X2, nDOFZ )

    CALL TimersStop( Timer_Ex_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    ! --- Compute Primitive Fluid in Spatial Elements ---

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_K(iNodeX,iCF_D       ,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_S1      ,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_S2      ,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_S3      ,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_E       ,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_Ne      ,iZ2,iZ4,iZ3), &
                 uPF_K(iNodeX,iPF_D       ,iZ2,iZ4,iZ3), &
                 uPF_K(iNodeX,iPF_V1      ,iZ2,iZ4,iZ3), &
                 uPF_K(iNodeX,iPF_V2      ,iZ2,iZ4,iZ3), &
                 uPF_K(iNodeX,iPF_V3      ,iZ2,iZ4,iZ3), &
                 uPF_K(iNodeX,iPF_E       ,iZ2,iZ4,iZ3), &
                 uPF_K(iNodeX,iPF_Ne      ,iZ2,iZ4,iZ3), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Ex_Flux )

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNodeZ,iCR_N ,iZ1,iZ2,iZ4,iS,iZ3), &
                 uCR_K(iNodeZ,iCR_G1,iZ1,iZ2,iZ4,iS,iZ3), &
                 uCR_K(iNodeZ,iCR_G2,iZ1,iZ2,iZ4,iS,iZ3), &
                 uCR_K(iNodeZ,iCR_G3,iZ1,iZ2,iZ4,iS,iZ3), &
                 uPR_K(iPR_D ), uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), uPR_K(iPR_I3), &
                 uPF_K(iNodeX,iPF_V1      ,iZ2,iZ4,iZ3), &
                 uPF_K(iNodeX,iPF_V2      ,iZ2,iZ4,iZ3), &
                 uPF_K(iNodeX,iPF_V3      ,iZ2,iZ4,iZ3), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3), &
                 nIterations )

        Flux_K &
          = Flux_X2 &
              ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                uPR_K(iPR_I2), uPR_K(iPR_I3), &
                uPF_K(iNodeX,iPF_V1      ,iZ2,iZ4,iZ3), &
                uPF_K(iNodeX,iPF_V2      ,iZ2,iZ4,iZ3), &
                uPF_K(iNodeX,iPF_V3      ,iZ2,iZ4,iZ3), &
                GX_K (iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                GX_K (iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                GX_K (iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

        DO iCR = 1, nCR

          Flux_q(iNodeZ,iCR,iZ1,iZ2,iZ4,iS,iZ3) &
            = dZ1(iZ1) * dZ2(iZ2) * dZ4(iZ4) &
                * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
                * GX_K(iNodeX,iGF_SqrtGm,iZ2,iZ4,iZ3) &
                * Flux_K(iCR)

        END DO

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOFZ, One, dLdX2_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_X2, nDOFZ )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Increment )

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_X2(iNodeZ,iCR,iZ1,iZ2,iZ4,iS,iZ3)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Increment )

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X2


  SUBROUTINE ComputeIncrement_Divergence_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

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

    INTEGER  :: iNodeZ, iNodeE, iNodeX
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER  :: nZ(4), nZ_X3(4), nV_X3, nV, nX_X3
    INTEGER  :: nIterations
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF)
    REAL(DP) :: uPR_L(nPR), Flux_L(nCR)
    REAL(DP) :: uPR_R(nPR), Flux_R(nCR)
    REAL(DP) :: uPR_K(nPR), Flux_K(nCR)
    REAL(DP) :: &
      GX_K(nDOFX,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B1(4):iZ_E1(4))
    REAL(DP) :: &
      GX_F(nDOFX_X3,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      uCF_K(nDOFX,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B1(4):iZ_E1(4))
    REAL(DP) :: &
      uCF_L(nDOFX_X3,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      uCF_R(nDOFX_X3,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      uPF_K(nDOFX,nPF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      V_u(3,nDOFX_X3, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      uCR_K(nDOFZ,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            nSpecies, &
            iZ_B1(4):iZ_E1(4))
    REAL(DP) :: &
      uCR_L(nDOF_X3,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            nSpecies, &
            iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      uCR_R(nDOF_X3,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            nSpecies, &
            iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      NumericalFlux &
        (nDOF_X3,nCR, &
         iZ_B0(1):iZ_E0(1), &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         nSpecies, &
         iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      Flux_q(nDOFZ,nCR, &
             iZ_B0(1):iZ_E0(1), &
             iZ_B0(2):iZ_E0(2), &
             iZ_B0(3):iZ_E0(3), &
             nSpecies, &
             iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dU_X3(nDOFZ,nCR, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            nSpecies, &
            iZ_B0(4):iZ_E0(4))

    IF( iZ_E0(4) .EQ. iZ_B0(4) ) RETURN

    PRINT*, "      ComputeIncrement_Divergence_X3"

    nZ    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nZ_X3 = nZ + [0,0,0,1]    ! Number of X3 Faces per Phase Space Dimension
    nV    = nCR * nSpecies * PRODUCT( nZ )
    nV_X3 = nCR * nSpecies * PRODUCT( nZ_X3 )
    nX_X3 = PRODUCT( nZ_X3(2:4) ) ! Number of X3 Faces in Position Space

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width )

    ! --- Permute Geometry Fields ---

    CALL TimersStart( Timer_Ex_Permute )

    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iGF    = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX_X3*nGF, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             GX_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1), nDOFX, Zero, &
             GX_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX_X3*nGF, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             GX_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX, Half, &
             GX_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3 )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Recompute Geometry from Scale Factors ---

    DO iZ4  = iZ_B0(4), iZ_E1(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X3

        GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4) &
          = MAX( GX_F(iNodeX,iGF_h_1,iZ2,iZ3,iZ4)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4) &
          = MAX( GX_F(iNodeX,iGF_h_2,iZ2,iZ3,iZ4)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) &
          = MAX( GX_F(iNodeX,iGF_h_3,iZ2,iZ3,iZ4)**2, SqrtTiny )
        GX_F(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
          = SQRT(   GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4) &
                  * GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4) &
                  * GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Permute Fluid Fields ---

    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iCF    = 1, nCF
      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iCF,iZ2,iZ3,iZ4) = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Fluid Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX_X3*nCF, nDOFX, One, LX_X3_Up, nDOFX_X3, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1), nDOFX, Zero, &
             uCF_L(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX_X3*nCF, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX, Zero, &
             uCF_R(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3 )

    CALL TimersStop( Timer_Ex_Interpolate )

    ! --- Compute Face Velocity Components ---

    DO iZ4 = iZ_B0(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X3

        ! --- Left State ---

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_L(iNodeX,iCF_D ,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_S1,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_S2,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_S3,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_E ,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_Ne,iZ2,iZ3,iZ4), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Right State ---

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_R(iNodeX,iCF_D ,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_S1,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_S2,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_S3,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_E ,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_Ne,iZ2,iZ3,iZ4), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Face Velocity ---

        V_u(1:3,iNodeX,iZ2,iZ3,iZ4) &
          = FaceVelocity_X3 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Permute Radiation Fields ---

    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iS  = 1, nSpecies
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iCR = 1, nCR

      DO iNodeZ = 1, nDOFZ

        uCR_K(iNodeZ,iCR,iZ1,iZ2,iZ3,iS,iZ4) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Radiation Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X3, nV_X3, nDOFZ, One, L_X3_Up, nDOF_X3, &
             uCR_K(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)-1), nDOFZ, Zero, &
             uCR_L(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)  ), nDOF_X3 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X3, nV_X3, nDOFZ, One, L_X3_Dn, nDOF_X3, &
             uCR_K(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)  ), nDOFZ, Zero, &
             uCR_R(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)  ), nDOF_X3 )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Flux )

    ! --- Numerical Flux ---

    DO iZ4 = iZ_B0(4), iZ_E1(4)
    DO iS  = 1, nSpecies
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX_X3
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        ! --- Left State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_L(iNodeZ,iCR_N ,iZ1,iZ2,iZ3,iS,iZ4), &
                 uCR_L(iNodeZ,iCR_G1,iZ1,iZ2,iZ3,iS,iZ4), &
                 uCR_L(iNodeZ,iCR_G2,iZ1,iZ2,iZ3,iS,iZ4), &
                 uCR_L(iNodeZ,iCR_G3,iZ1,iZ2,iZ3,iS,iZ4), &
                 uPR_L(iPR_D ), uPR_L(iPR_I1), &
                 uPR_L(iPR_I2), uPR_L(iPR_I3), &
                 V_u(1,iNodeX,iZ2,iZ3,iZ4), &
                 V_u(2,iNodeX,iZ2,iZ3,iZ4), &
                 V_u(3,iNodeX,iZ2,iZ3,iZ4), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                 nIterations )

        ! --- Left State Flux ---

        Flux_L &
          = Flux_X3 &
              ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                uPR_L(iPR_I2), uPR_L(iPR_I3), &
                V_u(1,iNodeX,iZ2,iZ3,iZ4), &
                V_u(2,iNodeX,iZ2,iZ3,iZ4), &
                V_u(3,iNodeX,iZ2,iZ3,iZ4), &
                GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Right State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_R(iNodeZ,iCR_N ,iZ1,iZ2,iZ3,iS,iZ4), &
                 uCR_R(iNodeZ,iCR_G1,iZ1,iZ2,iZ3,iS,iZ4), &
                 uCR_R(iNodeZ,iCR_G2,iZ1,iZ2,iZ3,iS,iZ4), &
                 uCR_R(iNodeZ,iCR_G3,iZ1,iZ2,iZ3,iS,iZ4), &
                 uPR_R(iPR_D ), uPR_R(iPR_I1), &
                 uPR_R(iPR_I2), uPR_R(iPR_I3), &
                 V_u(1,iNodeX,iZ2,iZ3,iZ4), &
                 V_u(2,iNodeX,iZ2,iZ3,iZ4), &
                 V_u(3,iNodeX,iZ2,iZ3,iZ4), &
                 GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                 nIterations )

        ! --- Right State Flux ---

        Flux_R &
          = Flux_X3 &
              ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                uPR_R(iPR_I2), uPR_R(iPR_I3), &
                V_u(1,iNodeX,iZ2,iZ3,iZ4), &
                V_u(2,iNodeX,iZ2,iZ3,iZ4), &
                V_u(3,iNodeX,iZ2,iZ3,iZ4), &
                GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Numerical Flux ---

        DO iCR = 1, nCR

          NumericalFlux(iNodeZ,iCR,iZ1,iZ2,iZ3,iS,iZ4) &
            = NumericalFlux_LLF &
                ( uCR_L(iNodeZ,iCR,iZ1,iZ2,iZ3,iS,iZ4), &
                  uCR_R(iNodeZ,iCR,iZ1,iZ2,iZ3,iS,iZ4), &
                  Flux_L(iCR), Flux_R(iCR), One )

          NumericalFlux(iNodeZ,iCR,iZ1,iZ2,iZ3,iS,iZ4) &
            = dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) &
                * Weights_X3(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
                * GX_F(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
                * NumericalFlux(iNodeZ,iCR,iZ1,iZ2,iZ3,iS,iZ4)

        END DO

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_X3, + One, L_X3_Dn, nDOF_X3, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)  ), &
             nDOF_X3, Zero, dU_X3, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_X3, - One, L_X3_Up, nDOF_X3, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(2),iZ_B0(3),1,iZ_B0(4)+1), &
             nDOF_X3, One,  dU_X3, nDOFZ )

    CALL TimersStop( Timer_Ex_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    ! --- Compute Primitive Fluid in Spatial Elements ---

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_K(iNodeX,iCF_D       ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S1      ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S2      ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S3      ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_E       ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_Ne      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_D       ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V1      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V2      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V3      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_E       ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_Ne      ,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Ex_Flux )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iS  = 1, nSpecies
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNodeZ,iCR_N ,iZ1,iZ2,iZ3,iS,iZ4), &
                 uCR_K(iNodeZ,iCR_G1,iZ1,iZ2,iZ3,iS,iZ4), &
                 uCR_K(iNodeZ,iCR_G2,iZ1,iZ2,iZ3,iS,iZ4), &
                 uCR_K(iNodeZ,iCR_G3,iZ1,iZ2,iZ3,iS,iZ4), &
                 uPR_K(iPR_D ), uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), uPR_K(iPR_I3), &
                 uPF_K(iNodeX,iPF_V1      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V2      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V3      ,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                 nIterations )

        Flux_K &
          = Flux_X3 &
              ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                uPR_K(iPR_I2), uPR_K(iPR_I3), &
                uPF_K(iNodeX,iPF_V1      ,iZ2,iZ3,iZ4), &
                uPF_K(iNodeX,iPF_V2      ,iZ2,iZ3,iZ4), &
                uPF_K(iNodeX,iPF_V3      ,iZ2,iZ3,iZ4), &
                GX_K (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                GX_K (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                GX_K (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        DO iCR = 1, nCR

          Flux_q(iNodeZ,iCR,iZ1,iZ2,iZ3,iS,iZ4) &
            = dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) &
                * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
                * GX_K(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
                * Flux_K(iCR)

        END DO

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOFZ, One, dLdX3_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_X3, nDOFZ )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Increment )

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_X3(iNodeZ,iCR,iZ1,iZ2,iZ3,iS,iZ4)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Increment )

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X3


  SUBROUTINE ComputeIncrement_ObserverCorrections &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

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

    INTEGER  :: iNode, iNodeZ, iNodeE, iNodeX
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER  :: nK(4), nK_Z1(4), nV, nV_Z1
    REAL(DP) :: EdgeEnergyCubed, Alpha
    REAL(DP) :: k_ud_11_L, k_ud_12_L, k_ud_13_L
    REAL(DP) ::            k_ud_22_L, k_ud_23_L
    REAL(DP) ::                       k_ud_33_L
    REAL(DP) :: k_ud_11_R, k_ud_12_R, k_ud_13_R
    REAL(DP) ::            k_ud_22_R, k_ud_23_R
    REAL(DP) ::                       k_ud_33_R
    REAL(DP) :: uPR_K(nPR), Flux_K(nCR)
    REAL(DP) :: uPR_L(nPR), Flux_L(nCR)
    REAL(DP) :: uPR_R(nPR), Flux_R(nCR)
    REAL(DP) :: &
      uGF_K(nDOFX,nGF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      uCF_K(nDOFX,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      uPF_K(nDOFX,nPF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dV_u_dX1 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dV_d_dX1 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dV_u_dX2 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dV_d_dX2 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dV_u_dX3 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dV_d_dX3 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      uCR_K(nDOFZ,nCR, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            nSpecies, &
            iZ_B1(1):iZ_E1(1))
    REAL(DP) :: uCR_L &
                  (nDOF_E,nCR, &
                   iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                   iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(1):iZ_E1(1))
    REAL(DP) :: uCR_R &
                  (nDOF_E,nCR, &
                   iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                   iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(1):iZ_E1(1))
    REAL(DP) :: NumericalFlux &
                  (nDOF_E,nCR, &
                   iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                   iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(1):iZ_E1(1))
    REAL(DP) :: dU_E &
                  (nDOFZ,nCR, &
                   iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                   iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(1):iZ_E0(1))
    REAL(DP) :: Flux_q &
                  (nDOFZ,nCR, &
                   iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                   iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(1):iZ_E0(1))

    IF( iZ_E0(1) .EQ. iZ_B0(1) ) RETURN

    PRINT*, "      ComputeIncrement_ObserverCorrections"

    CALL ComputeWeakDerivatives_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, dV_u_dX1, dV_d_dX1 )

    CALL ComputeWeakDerivatives_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, dV_u_dX2, dV_d_dX2 )

    CALL ComputeWeakDerivatives_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, dV_u_dX3, dV_d_dX3 )

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_Z1 = nK + [1,0,0,0]    ! Number of Z1 Faces per Phase Space Dimension
    nV    = nCR * nSpecies * PRODUCT( nK )
    nV_Z1 = nCR * nSpecies * PRODUCT( nK_Z1 )

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, &
        dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Permute Geometry Fields ---

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iGF = 1, nGF

      DO iNodeX = 1, nDOFX

        uGF_K(iNodeX,iGF,iZ2,iZ3,iZ4) &
          = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Permute Fluid Fields ---

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iCF = 1, nCF

      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iCF,iZ2,iZ3,iZ4) &
          = U_F(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Permute Radiation Fields ---

    DO iZ1 = iZ_B1(1), iZ_E1(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iCR = 1, nCR

      DO iNodeZ = 1, nDOFZ

        uCR_K(iNodeZ,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Compute Primitive Fluid ---

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_K(iNodeX,iCF_D       ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S1      ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S2      ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S3      ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_E       ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_Ne      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_D       ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V1      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V2      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V3      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_E       ,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_Ne      ,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Radiation Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_E, nV_Z1, nDOFZ, One, L_E_Up, nDOF_E, &
             uCR_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)-1), nDOFZ, Zero, &
             uCR_L(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ), nDOF_E )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_E, nV_Z1, nDOFZ, One, L_E_Dn, nDOF_E, &
             uCR_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ), nDOFZ, Zero, &
             uCR_R(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ), nDOF_E )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Flux )

    ! --- Numerical Flux ---

    DO iZ1 = iZ_B0(1), iZ_E1(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNode = 1, nDOF_E ! = nDOFX

        ! --- Left State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_L(iNode,iCR_N       ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_L(iNode,iCR_G1      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_L(iNode,iCR_G2      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_L(iNode,iCR_G3      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uPR_L(iPR_D ), uPR_L(iPR_I1), &
                 uPR_L(iPR_I2), uPR_L(iPR_I3), &
                 uPF_K(iNode,iPF_V1      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNode,iPF_V2      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNode,iPF_V3      ,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Left State Flux ---

        Flux_L &
          = Flux_E( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                    uPR_L(iPR_I2), uPR_L(iPR_I3), &
                    dV_u_dX1(iNode,1           ,iZ2,iZ3,iZ4), &
                    dV_u_dX1(iNode,2           ,iZ2,iZ3,iZ4), &
                    dV_u_dX1(iNode,3           ,iZ2,iZ3,iZ4), &
                    dV_u_dX2(iNode,1           ,iZ2,iZ3,iZ4), &
                    dV_u_dX2(iNode,2           ,iZ2,iZ3,iZ4), &
                    dV_u_dX2(iNode,3           ,iZ2,iZ3,iZ4), &
                    dV_u_dX3(iNode,1           ,iZ2,iZ3,iZ4), &
                    dV_u_dX3(iNode,2           ,iZ2,iZ3,iZ4), &
                    dV_u_dX3(iNode,3           ,iZ2,iZ3,iZ4), &
                    uGF_K   (iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                    uGF_K   (iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                    uGF_K   (iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Right State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_R(iNode,iCR_N       ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_R(iNode,iCR_G1      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_R(iNode,iCR_G2      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_R(iNode,iCR_G3      ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uPR_R(iPR_D ), uPR_R(iPR_I1), &
                 uPR_R(iPR_I2), uPR_R(iPR_I3), &
                 uPF_K(iNode,iPF_V1      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNode,iPF_V2      ,iZ2,iZ3,iZ4), &
                 uPF_K(iNode,iPF_V3      ,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Right State Flux ---

        Flux_R &
          = Flux_E( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                    uPR_R(iPR_I2), uPR_R(iPR_I3), &
                    dV_u_dX1(iNode,1           ,iZ2,iZ3,iZ4), &
                    dV_u_dX1(iNode,2           ,iZ2,iZ3,iZ4), &
                    dV_u_dX1(iNode,3           ,iZ2,iZ3,iZ4), &
                    dV_u_dX2(iNode,1           ,iZ2,iZ3,iZ4), &
                    dV_u_dX2(iNode,2           ,iZ2,iZ3,iZ4), &
                    dV_u_dX2(iNode,3           ,iZ2,iZ3,iZ4), &
                    dV_u_dX3(iNode,1           ,iZ2,iZ3,iZ4), &
                    dV_u_dX3(iNode,2           ,iZ2,iZ3,iZ4), &
                    dV_u_dX3(iNode,3           ,iZ2,iZ3,iZ4), &
                    uGF_K   (iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                    uGF_K   (iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                    uGF_K   (iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Numerical Flux (Upwind) ---

        EdgeEnergyCubed &
          = ( MeshE % Center(iZ1) - Half * MeshE % Width(iZ1) )**3

        CALL ComputeEddingtonTensorComponents_ud &
               ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                 uPR_L(iPR_I2), uPR_L(iPR_I3), &
                 uGF_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                 k_ud_11_L, k_ud_12_L, k_ud_13_L, &
                            k_ud_22_L, k_ud_23_L, &
                                       k_ud_33_L )

        CALL ComputeEddingtonTensorComponents_ud &
               ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                 uPR_R(iPR_I2), uPR_R(iPR_I3), &
                 uGF_K(iNode,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNode,iGF_Gm_dd_33,iZ2,iZ3,iZ4), &
                 k_ud_11_R, k_ud_12_R, k_ud_13_R, &
                            k_ud_22_R, k_ud_23_R, &
                                       k_ud_33_R )

        ! --- Alpha := k^{i}_{j} d_{i}v^{j} ---

        Alpha &
          =   ( k_ud_11_L + k_ud_11_R )               &
                * (   dV_u_dX1(iNode,1,iZ2,iZ3,iZ4) ) &
            + ( k_ud_22_L + k_ud_22_R )               &
                * (   dV_u_dX2(iNode,2,iZ2,iZ3,iZ4) ) &
            + ( k_ud_33_L + k_ud_33_R )               &
                * (   dV_u_dX3(iNode,3,iZ2,iZ3,iZ4) ) &
            + ( k_ud_12_L + k_ud_12_R )               &
                * (   dV_u_dX1(iNode,2,iZ2,iZ3,iZ4)   &
                    + dV_u_dX2(iNode,1,iZ2,iZ3,iZ4) ) &
            + ( k_ud_13_L + k_ud_13_R )               &
                * (   dV_u_dX1(iNode,3,iZ2,iZ3,iZ4)   &
                    + dV_u_dX3(iNode,1,iZ2,iZ3,iZ4) ) &
            + ( k_ud_23_L + k_ud_23_R )               &
                * (   dV_u_dX2(iNode,3,iZ2,iZ3,iZ4)   &
                    + dV_u_dX3(iNode,2,iZ2,iZ3,iZ4) )

        Alpha = SIGN( Half, Alpha ) + Half

        DO iCR = 1, nCR

          NumericalFlux(iNode,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
            = (One-Alpha) * Flux_L(iCR) + Alpha * Flux_R(iCR)

          NumericalFlux(iNode,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
            = dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                * EdgeEnergyCubed * Weights_E(iNode) &
                * uGF_K(iNode,iGF_SqrtGm,iZ2,iZ3,iZ4) &
                * NumericalFlux(iNode,iCR,iZ2,iZ3,iZ4,iS,iZ1)

        END DO

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Surface Contributions ---

    ! --- Contributions from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_E, + One, L_E_Dn, nDOF_E, &
             NumericalFlux(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)  ), &
             nDOF_E, Zero, dU_E, nDOFZ )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOF_E, - One, L_E_Up, nDOF_E, &
             NumericalFlux(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4),1,iZ_B0(1)+1), &
             nDOF_E, One,  dU_E, nDOFZ )

    CALL TimersStop( Timer_Ex_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart( Timer_Ex_Flux )

    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNodeZ,iCR_N ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1), &
                 uPR_K(iPR_D ), &
                 uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), &
                 uPR_K(iPR_I3), &
                 uPF_K(iNodeX,iPF_V1,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V2,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V3,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        Flux_K &
          = Flux_E &
              ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                uPR_K(iPR_I2), uPR_K(iPR_I3), &
                dV_u_dX1(iNodeX,1           ,iZ2,iZ3,iZ4), &
                dV_u_dX1(iNodeX,2           ,iZ2,iZ3,iZ4), &
                dV_u_dX1(iNodeX,3           ,iZ2,iZ3,iZ4), &
                dV_u_dX2(iNodeX,1           ,iZ2,iZ3,iZ4), &
                dV_u_dX2(iNodeX,2           ,iZ2,iZ3,iZ4), &
                dV_u_dX2(iNodeX,3           ,iZ2,iZ3,iZ4), &
                dV_u_dX3(iNodeX,1           ,iZ2,iZ3,iZ4), &
                dV_u_dX3(iNodeX,2           ,iZ2,iZ3,iZ4), &
                dV_u_dX3(iNodeX,3           ,iZ2,iZ3,iZ4), &
                uGF_K   (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                uGF_K   (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                uGF_K   (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        DO iCR = 1, nCR

          Flux_q(iNodeZ,iCR,iZ2,iZ3,iZ4,iS,iZ1) &
            = dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) * Weights_q(iNodeZ) &
                * GE(iNodeE,iZ1,iGE_Ep3) &
                * uGF_K(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
                * Flux_K(iCR)

        END DO

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Flux )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFZ, nV, nDOFZ, One, dLdE_q, nDOFZ, &
             Flux_q, nDOFZ, One, dU_E, nDOFZ )

    CALL TimersStop( Timer_Ex_Interpolate )

    !--------------------------------------------------
    !--- Volume Sources (Number Flux Equation Only) ---
    !--------------------------------------------------

    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNodeZ,iCR_N ,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1), &
                 uCR_K(iNodeZ,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1), &
                 uPR_K(iPR_D ), &
                 uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), &
                 uPR_K(iPR_I3), &
                 uPF_K(iNodeX,iPF_V1,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V2,iZ2,iZ3,iZ4), &
                 uPF_K(iNodeX,iPF_V3,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 uGF_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        Flux_K &
          = Flux_E &
              ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                uPR_K(iPR_I2), uPR_K(iPR_I3), &
                dV_u_dX1(iNodeX,1           ,iZ2,iZ3,iZ4), &
                dV_u_dX1(iNodeX,2           ,iZ2,iZ3,iZ4), &
                dV_u_dX1(iNodeX,3           ,iZ2,iZ3,iZ4), &
                dV_u_dX2(iNodeX,1           ,iZ2,iZ3,iZ4), &
                dV_u_dX2(iNodeX,2           ,iZ2,iZ3,iZ4), &
                dV_u_dX2(iNodeX,3           ,iZ2,iZ3,iZ4), &
                dV_u_dX3(iNodeX,1           ,iZ2,iZ3,iZ4), &
                dV_u_dX3(iNodeX,2           ,iZ2,iZ3,iZ4), &
                dV_u_dX3(iNodeX,3           ,iZ2,iZ3,iZ4), &
                uGF_K   (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                uGF_K   (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                uGF_K   (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- iCR_G1 ---

        dU_E(iNodeZ,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1) &
          = dU_E(iNodeZ,iCR_G1,iZ2,iZ3,iZ4,iS,iZ1) &
            - dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * uGF_K(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
              * ( Flux_K(iCR_G1) &
                  + uPR_K(iPR_I1) * dV_d_dX1(iNodeX,1,iZ2,iZ3,iZ4) &
                  + uPR_K(iPR_I2) * dV_d_dX2(iNodeX,1,iZ2,iZ3,iZ4) &
                  + uPR_K(iPR_I3) * dV_d_dX3(iNodeX,1,iZ2,iZ3,iZ4) )

        ! --- iCR_G2 ---

        dU_E(iNodeZ,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1) &
          = dU_E(iNodeZ,iCR_G2,iZ2,iZ3,iZ4,iS,iZ1) &
            - dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * uGF_K(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
              * ( Flux_K(iCR_G2) &
                  + uPR_K(iPR_I1) * dV_d_dX1(iNodeX,2,iZ2,iZ3,iZ4) &
                  + uPR_K(iPR_I2) * dV_d_dX2(iNodeX,2,iZ2,iZ3,iZ4) &
                  + uPR_K(iPR_I3) * dV_d_dX3(iNodeX,2,iZ2,iZ3,iZ4) )

        ! --- iCR_G3 ---

        dU_E(iNodeZ,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1) &
          = dU_E(iNodeZ,iCR_G3,iZ2,iZ3,iZ4,iS,iZ1) &
            - dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
              * Weights_q(iNodeZ) * GE(iNodeE,iZ1,iGE_Ep2) &
              * uGF_K(iNodeX,iGF_SqrtGm,iZ2,iZ3,iZ4) &
              * ( Flux_K(iCR_G3) &
                  + uPR_K(iPR_I1) * dV_d_dX1(iNodeX,3,iZ2,iZ3,iZ4) &
                  + uPR_K(iPR_I2) * dV_d_dX2(iNodeX,3,iZ2,iZ3,iZ4) &
                  + uPR_K(iPR_I3) * dV_d_dX3(iNodeX,3,iZ2,iZ3,iZ4) )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Ex_Increment )

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_E(iNodeZ,iCR,iZ2,iZ3,iZ4,iS,iZ1)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Increment )

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeIncrement_ObserverCorrections


  SUBROUTINE ComputeWeakDerivatives_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, dV_u_dX1_Out, dV_d_dX1_Out )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)  :: &
      uCF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(out) :: &
      dV_u_dX1_Out &
         (1:nDOFX,1:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
      dV_d_dX1_Out &
         (1:nDOFX,1:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    INTEGER  :: nK(4), nK_X1(4), nX, nX_X1
    INTEGER  :: iNodeX
    INTEGER  :: i, iZ2, iZ3, iZ4, iCF, iGF
    REAL(DP) :: &
      uPF_K(nPF), uPF_L(nPF), uPF_R(nPF)
    REAL(DP) :: &
      V_u_X1(nDOFX_X1,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
             iZ_B0(2):iZ_E0(2)+1), &
      V_d_X1(nDOFX_X1,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
             iZ_B0(2):iZ_E0(2)+1), &
      V_u_K(nDOFX,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)), &
      V_d_K(nDOFX,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)), &
      dV_u_dX1(nDOFX,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
               iZ_B0(2):iZ_E0(2)), &
      dV_d_dX1(nDOFX,3,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
               iZ_B0(2):iZ_E0(2))
    REAL(DP) :: &
      GX_K(nDOFX   ,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B1(2):iZ_E1(2)  ,nGF), &
      GX_F(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)+1,nGF)
    REAL(DP) :: &
      uCF_K(nDOFX   ,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B1(2):iZ_E1(2)  ,nCF), &
      uCF_L(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)+1,nCF), &
      uCF_R(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
            iZ_B0(2):iZ_E0(2)+1,nCF)

    IF( iZ_E0(2) .EQ. iZ_B0(2) )THEN
      dV_u_dX1_Out = Zero
      dV_d_dX1_Out = Zero
      RETURN
    END IF

    PRINT*, "      ComputeWeakDerivatives_X1"

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X1 = nK + [0,1,0,0]    ! Number of X1 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X1 = PRODUCT( nK_X1(2:4) ) ! Number of X1 Faces in Position Space

    ! --- Permute Geometry Fields ---

    DO iGF = 1, nGF
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iZ3,iZ4,iZ2,iGF) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Geometry Fields ---

    DO iGF = iGF_h_1, iGF_h_3

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF), nDOFX, Zero, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX, Half, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

    END DO

    ! --- Compute Metric Components ---

    DO i   = 1, 3
    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11+i-1) &
          = MAX( GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_h_1+i-1)**2, SqrtTiny )

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Permute Fluid Fields ---

    DO iCF = 1, nCF
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF) = uCF(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    DO iCF = 1, nCF

      ! --- Interpolate Left State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               uCF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iCF), nDOFX, Zero, &
               uCF_L(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX_X1 )

      ! --- Interpolate Right State ---

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               uCF_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX, Zero, &
               uCF_R(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iCF), nDOFX_X1 )

    END DO

    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1

        ! --- Left States ---

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_D ), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_S1), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_S2), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_S3), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_E ), &
                 uCF_L(iNodeX,iZ3,iZ4,iZ2,iCF_Ne), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

        ! --- Right States ---

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_D ), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_S1), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_S2), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_S3), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_E ), &
                 uCF_R(iNodeX,iZ3,iZ4,iZ2,iCF_Ne), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
                 GX_F (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

        V_u_X1(iNodeX,1:3,iZ3,iZ4,iZ2) &
          = FaceVelocity_X1 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )

        V_u_X1(iNodeX,1,iZ3,iZ4,iZ2) &
          = V_u_X1(iNodeX,1,iZ3,iZ4,iZ2) * WeightsX_X1(iNodeX)

        V_u_X1(iNodeX,2,iZ3,iZ4,iZ2) &
          = V_u_X1(iNodeX,2,iZ3,iZ4,iZ2) * WeightsX_X1(iNodeX)

        V_u_X1(iNodeX,3,iZ3,iZ4,iZ2) &
          = V_u_X1(iNodeX,3,iZ3,iZ4,iZ2) * WeightsX_X1(iNodeX)

        V_d_X1(iNodeX,1,iZ3,iZ4,iZ2) &
          = GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
              * V_u_X1(iNodeX,1,iZ3,iZ4,iZ2)

        V_d_X1(iNodeX,2,iZ3,iZ4,iZ2) &
          = GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
              * V_u_X1(iNodeX,2,iZ3,iZ4,iZ2)

        V_d_X1(iNodeX,3,iZ3,iZ4,iZ2) &
          = GX_F(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) &
              * V_u_X1(iNodeX,3,iZ3,iZ4,iZ2)

      END DO

    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             V_u_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, Zero, &
             dV_u_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             V_d_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, Zero, &
             dV_d_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             V_u_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, One,  &
             dV_u_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             V_d_X1(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, One,  &
             dV_d_dX1, nDOFX )

    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_D ), &
                 uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_S1), &
                 uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_S2), &
                 uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_S3), &
                 uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_E ), &
                 uCF_K(iNodeX,iZ3,iZ4,iZ2,iCF_Ne), &
                 uPF_K(iPF_D ), &
                 uPF_K(iPF_V1), &
                 uPF_K(iPF_V2), &
                 uPF_K(iPF_V3), &
                 uPF_K(iPF_E ), &
                 uPF_K(iPF_Ne), &
                 GX_K (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11), &
                 GX_K (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22), &
                 GX_K (iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) )

        V_u_K(iNodeX,1,iZ3,iZ4,iZ2) &
          = uPF_K(iPF_V1) * WeightsX_q(iNodeX)

        V_u_K(iNodeX,2,iZ3,iZ4,iZ2) &
          = uPF_K(iPF_V2) * WeightsX_q(iNodeX)

        V_u_K(iNodeX,3,iZ3,iZ4,iZ2) &
          = uPF_K(iPF_V3) * WeightsX_q(iNodeX)

        V_d_K(iNodeX,1,iZ3,iZ4,iZ2) &
          = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_11) &
              * V_u_K(iNodeX,1,iZ3,iZ4,iZ2)

        V_d_K(iNodeX,2,iZ3,iZ4,iZ2) &
          = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_22) &
              * V_u_K(iNodeX,2,iZ3,iZ4,iZ2)

        V_d_K(iNodeX,3,iZ3,iZ4,iZ2) &
          = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF_Gm_dd_33) &
              * V_u_K(iNodeX,3,iZ3,iZ4,iZ2)

      END DO

    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             V_u_K, nDOFX, One, dV_u_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             V_d_K, nDOFX, One, dV_d_dX1, nDOFX )

    ASSOCIATE( dZ2 => MeshX(1) % Width )

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO i = 1, 3
      DO iNodeX = 1, nDOFX

        dV_u_dX1(iNodeX,i,iZ3,iZ4,iZ2) &
         = dV_u_dX1(iNodeX,i,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dV_d_dX1(iNodeX,i,iZ3,iZ4,iZ2) &
         = dV_d_dX1(iNodeX,i,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

      END DO
      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO i = 1, 3
      DO iNodeX = 1, nDOFX

        dV_u_dX1_Out(iNodeX,i,iZ2,iZ3,iZ4) &
          = dV_u_dX1(iNodeX,i,iZ3,iZ4,iZ2)

        dV_d_dX1_Out(iNodeX,i,iZ2,iZ3,iZ4) &
          = dV_d_dX1(iNodeX,i,iZ3,iZ4,iZ2)

      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeWeakDerivatives_X1


  SUBROUTINE ComputeWeakDerivatives_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, dV_u_dX2_Out, dV_d_dX2_Out )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in)  :: &
      uCF(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCF)
    REAL(DP), INTENT(out) :: &
      dV_u_dX2_Out &
         (1:nDOFX,1:3, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4)), &
      dV_d_dX2_Out &
         (1:nDOFX,1:3, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4))

    INTEGER  :: iNodeX
    INTEGER  :: i, iZ2, iZ3, iZ4, iGF, iCF
    INTEGER  :: nK(4), nK_X2(4), nX, nX_X2
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF), uPF_K(nPF)
    REAL(DP) :: &
      GX_K(nDOFX,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B1(3):iZ_E1(3))
    REAL(DP) :: &
      GX_F(nDOFX_X2,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      uCF_K(nDOFX,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B1(3):iZ_E1(3))
    REAL(DP) :: &
      uCF_L(nDOFX_X2,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      uCF_R(nDOFX_X2,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      V_u_F(nDOFX_X2,3, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      V_d_F(nDOFX_X2,3, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E1(3))
    REAL(DP) :: &
      V_u_K(nDOFX,3, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E0(3))
    REAL(DP) :: &
      V_d_K(nDOFX,3, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(4):iZ_E0(4), &
            iZ_B0(3):iZ_E0(3))
    REAL(DP) :: &
      dV_u_dX2 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(4):iZ_E0(4), &
         iZ_B0(3):iZ_E0(3))
    REAL(DP) :: &
      dV_d_dX2 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(4):iZ_E0(4), &
         iZ_B0(3):iZ_E0(3))

    IF( iZ_E0(3) .EQ. iZ_B0(3) )THEN
      dV_u_dX2_Out = Zero
      dV_d_dX2_Out = Zero
      RETURN
    END IF

    PRINT*, "      ComputeWeakDerivatives_X2"

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X2 = nK + [0,0,1,0]    ! Number of X2 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X2 = PRODUCT( nK_X2(2:4) ) ! Number of X2 Faces in Position Space

    ! --- Permute Geometry Fields ---

    CALL TimersStart( Timer_Ex_Permute )

    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iGF    = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iZ2,iZ4,iZ3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX_X2*nGF, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             GX_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1), nDOFX, Zero, &
             GX_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX_X2*nGF, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             GX_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX, Half, &
             GX_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Compute Metric Components from Scale Factors ---

    DO iZ3  = iZ_B0(3), iZ_E1(3)
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X2

        GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3) &
          = MAX( GX_F(iNodeX,iGF_h_1,iZ2,iZ4,iZ3)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3) &
          = MAX( GX_F(iNodeX,iGF_h_2,iZ2,iZ4,iZ3)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) &
          = MAX( GX_F(iNodeX,iGF_h_3,iZ2,iZ4,iZ3)**2, SqrtTiny )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Permute Fluid Fields ---

    CALL TimersStart( Timer_Ex_Permute )

    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iCF    = 1, nCF
      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iCF,iZ2,iZ4,iZ3) = uCF(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Interpolate Fluid Fields ---

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX_X2*nCF, nDOFX, One, LX_X2_Up, nDOFX_X2, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1), nDOFX, Zero, &
             uCF_L(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX_X2*nCF, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX, Zero, &
             uCF_R(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2 )

    CALL TimersStop( Timer_Ex_Interpolate )

    DO iZ3 = iZ_B0(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X2

        ! --- Left States ---

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_L(iNodeX,iCF_D ,iZ2,iZ4,iZ3), &
                 uCF_L(iNodeX,iCF_S1,iZ2,iZ4,iZ3), &
                 uCF_L(iNodeX,iCF_S2,iZ2,iZ4,iZ3), &
                 uCF_L(iNodeX,iCF_S3,iZ2,iZ4,iZ3), &
                 uCF_L(iNodeX,iCF_E ,iZ2,iZ4,iZ3), &
                 uCF_L(iNodeX,iCF_Ne,iZ2,iZ4,iZ3), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F (iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_F (iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_F (iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

        ! --- Right States ---

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_R(iNodeX,iCF_D ,iZ2,iZ4,iZ3), &
                 uCF_R(iNodeX,iCF_S1,iZ2,iZ4,iZ3), &
                 uCF_R(iNodeX,iCF_S2,iZ2,iZ4,iZ3), &
                 uCF_R(iNodeX,iCF_S3,iZ2,iZ4,iZ3), &
                 uCF_R(iNodeX,iCF_E ,iZ2,iZ4,iZ3), &
                 uCF_R(iNodeX,iCF_Ne,iZ2,iZ4,iZ3), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F (iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_F (iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_F (iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

        V_u_F(iNodeX,1:3,iZ2,iZ4,iZ3) &
          = FaceVelocity_X2 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )

        V_u_F(iNodeX,1,iZ2,iZ4,iZ3) &
          = V_u_F(iNodeX,1,iZ2,iZ4,iZ3) * WeightsX_X2(iNodeX)

        V_u_F(iNodeX,2,iZ2,iZ4,iZ3) &
          = V_u_F(iNodeX,2,iZ2,iZ4,iZ3) * WeightsX_X2(iNodeX)

        V_u_F(iNodeX,3,iZ2,iZ4,iZ3) &
          = V_u_F(iNodeX,3,iZ2,iZ4,iZ3) * WeightsX_X2(iNodeX)

        V_d_F(iNodeX,1,iZ2,iZ4,iZ3) &
          = GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3) &
              * V_u_F(iNodeX,1,iZ2,iZ4,iZ3)

        V_d_F(iNodeX,2,iZ2,iZ4,iZ3) &
          = GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3) &
              * V_u_F(iNodeX,2,iZ2,iZ4,iZ3)

        V_d_F(iNodeX,3,iZ2,iZ4,iZ3) &
          = GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) &
              * V_u_F(iNodeX,3,iZ2,iZ4,iZ3)

      END DO

    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             V_u_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2, Zero, &
             dV_u_dX2, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             V_d_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2, Zero, &
             dV_d_dX2, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             V_u_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)+1), nDOFX_X2, One,  &
             dV_u_dX2, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             V_d_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)+1), nDOFX_X2, One,  &
             dV_d_dX2, nDOFX )

    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_K(iNodeX,iCF_D ,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_S1,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_S2,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_S3,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_E ,iZ2,iZ4,iZ3), &
                 uCF_K(iNodeX,iCF_Ne,iZ2,iZ4,iZ3), &
                 uPF_K(iPF_D ), &
                 uPF_K(iPF_V1), &
                 uPF_K(iPF_V2), &
                 uPF_K(iPF_V3), &
                 uPF_K(iPF_E ), &
                 uPF_K(iPF_Ne), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) )

        V_u_K(iNodeX,1,iZ2,iZ4,iZ3) &
          = uPF_K(iPF_V1) * WeightsX_q(iNodeX)

        V_u_K(iNodeX,2,iZ2,iZ4,iZ3) &
          = uPF_K(iPF_V2) * WeightsX_q(iNodeX)

        V_u_K(iNodeX,3,iZ2,iZ4,iZ3) &
          = uPF_K(iPF_V3) * WeightsX_q(iNodeX)

        V_d_K(iNodeX,1,iZ2,iZ4,iZ3) &
          = GX_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ4,iZ3) &
              * V_u_K(iNodeX,1,iZ2,iZ4,iZ3)

        V_d_K(iNodeX,2,iZ2,iZ4,iZ3) &
          = GX_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ4,iZ3) &
              * V_u_K(iNodeX,2,iZ2,iZ4,iZ3)

        V_d_K(iNodeX,3,iZ2,iZ4,iZ3) &
          = GX_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ4,iZ3) &
              * V_u_K(iNodeX,3,iZ2,iZ4,iZ3)

      END DO

    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX, - One, dLXdX2_q, nDOFX, &
             V_u_K, nDOFX, One, dV_u_dX2, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX, - One, dLXdX2_q, nDOFX, &
             V_d_K, nDOFX, One, dV_d_dX2, nDOFX )

    ASSOCIATE( dZ3 => MeshX(2) % Width )

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO i      = 1, 3
      DO iNodeX = 1, nDOFX

        dV_u_dX2(iNodeX,i,iZ2,iZ4,iZ3) &
         = dV_u_dX2(iNodeX,i,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dV_d_dX2(iNodeX,i,iZ2,iZ4,iZ3) &
         = dV_d_dX2(iNodeX,i,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

      END DO
      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO i      = 1, 3
      DO iNodeX = 1, nDOFX

        dV_u_dX2_Out(iNodeX,i,iZ2,iZ3,iZ4) &
          = dV_u_dX2(iNodeX,i,iZ2,iZ4,iZ3)

        dV_d_dX2_Out(iNodeX,i,iZ2,iZ3,iZ4) &
          = dV_d_dX2(iNodeX,i,iZ2,iZ4,iZ3)

      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeWeakDerivatives_X2


  SUBROUTINE ComputeWeakDerivatives_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, uCF, dV_u_dX3_Out, dV_d_dX3_Out )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in)  :: &
      uCF(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCF)
    REAL(DP), INTENT(out) :: &
      dV_u_dX3_Out &
         (1:nDOFX,1:3, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4)), &
      dV_d_dX3_Out &
         (1:nDOFX,1:3, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4))

    INTEGER  :: iNodeX
    INTEGER  :: i, iZ2, iZ3, iZ4, iGF, iCF
    INTEGER  :: nK(4), nK_X3(4), nX, nX_X3
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF), uPF_K(nPF)
    REAL(DP) :: &
      GX_K(nDOFX,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B1(4):iZ_E1(4))
    REAL(DP) :: &
      GX_F(nDOFX_X3,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      uCF_K(nDOFX,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B1(4):iZ_E1(4))
    REAL(DP) :: &
      uCF_L(nDOFX_X3,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      uCF_R(nDOFX_X3,nCF, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      V_u_F(nDOFX_X3,3, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      V_d_F(nDOFX_X3,3, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E1(4))
    REAL(DP) :: &
      V_u_K(nDOFX,3, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      V_d_K(nDOFX,3, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dV_u_dX3 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dV_d_dX3 &
        (nDOFX,3, &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))

    IF( iZ_E0(4) .EQ. iZ_B0(4) )THEN
      dV_u_dX3_Out = Zero
      dV_d_dX3_Out = Zero
      RETURN
    END IF

    PRINT*, "      ComputeWeakDerivatives_X3"

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X3 = nK + [0,0,0,1]    ! Number of X2 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X3 = PRODUCT( nK_X3(2:4) ) ! Number of X2 Faces in Position Space

    ! --- Permute Geometry Fields ---

    CALL TimersStart( Timer_Ex_Permute )

    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iGF    = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX_X3*nGF, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             GX_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1), nDOFX, Zero, &
             GX_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX_X3*nGF, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             GX_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX, Half, &
             GX_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3 )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Permute )

    ! --- Compute Metric Components from Scale Factors ---

    DO iZ4  = iZ_B0(4), iZ_E1(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X3

        GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4) &
          = MAX( GX_F(iNodeX,iGF_h_1,iZ2,iZ3,iZ4)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4) &
          = MAX( GX_F(iNodeX,iGF_h_2,iZ2,iZ3,iZ4)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) &
          = MAX( GX_F(iNodeX,iGF_h_3,iZ2,iZ3,iZ4)**2, SqrtTiny )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    ! --- Permute Fluid Fields ---

    CALL TimersStart( Timer_Ex_Permute )

    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iCF    = 1, nCF
      DO iNodeX = 1, nDOFX

        uCF_K(iNodeX,iCF,iZ2,iZ3,iZ4) = uCF(iNodeX,iZ2,iZ3,iZ4,iCF)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Fluid Fields ---

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX_X3*nCF, nDOFX, One, LX_X3_Up, nDOFX_X3, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1), nDOFX, Zero, &
             uCF_L(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX_X3*nCF, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
             uCF_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX, Zero, &
             uCF_R(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3 )

    CALL TimersStop( Timer_Ex_Interpolate )

    ! --- Compute Face Velocity Components ---

    DO iZ4 = iZ_B0(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X3

        ! --- Left States ---

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_L(iNodeX,iCF_D ,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_S1,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_S2,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_S3,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_E ,iZ2,iZ3,iZ4), &
                 uCF_L(iNodeX,iCF_Ne,iZ2,iZ3,iZ4), &
                 uPF_L(iPF_D ), &
                 uPF_L(iPF_V1), &
                 uPF_L(iPF_V2), &
                 uPF_L(iPF_V3), &
                 uPF_L(iPF_E ), &
                 uPF_L(iPF_Ne), &
                 GX_F (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_F (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_F (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Right States ---

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_R(iNodeX,iCF_D ,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_S1,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_S2,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_S3,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_E ,iZ2,iZ3,iZ4), &
                 uCF_R(iNodeX,iCF_Ne,iZ2,iZ3,iZ4), &
                 uPF_R(iPF_D ), &
                 uPF_R(iPF_V1), &
                 uPF_R(iPF_V2), &
                 uPF_R(iPF_V3), &
                 uPF_R(iPF_E ), &
                 uPF_R(iPF_Ne), &
                 GX_F (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_F (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_F (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Face Velocity ---

        V_u_F(iNodeX,1:3,iZ2,iZ3,iZ4) &
          = FaceVelocity_X3 &
              ( uPF_L(iPF_V1), uPF_L(iPF_V2), uPF_L(iPF_V3), &
                uPF_R(iPF_V1), uPF_R(iPF_V2), uPF_R(iPF_V3) )

        V_u_F(iNodeX,1,iZ2,iZ3,iZ4) &
          = V_u_F(iNodeX,1,iZ2,iZ3,iZ4) * WeightsX_X3(iNodeX)

        V_u_F(iNodeX,2,iZ2,iZ3,iZ4) &
          = V_u_F(iNodeX,2,iZ2,iZ3,iZ4) * WeightsX_X3(iNodeX)

        V_u_F(iNodeX,3,iZ2,iZ3,iZ4) &
          = V_u_F(iNodeX,3,iZ2,iZ3,iZ4) * WeightsX_X3(iNodeX)

        V_d_F(iNodeX,1,iZ2,iZ3,iZ4) &
          = GX_F(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4) &
              * V_u_F(iNodeX,1,iZ2,iZ3,iZ4)

        V_d_F(iNodeX,2,iZ2,iZ3,iZ4) &
          = GX_F(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4) &
              * V_u_F(iNodeX,2,iZ2,iZ3,iZ4)

        V_d_F(iNodeX,3,iZ2,iZ3,iZ4) &
          = GX_F(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) &
              * V_u_F(iNodeX,3,iZ2,iZ3,iZ4)

      END DO

    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             V_u_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3, Zero, &
             dV_u_dX3, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             V_d_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3, Zero, &
             dV_d_dX3, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             V_u_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)+1), nDOFX_X3, One,  &
             dV_u_dX3, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             V_d_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)+1), nDOFX_X3, One,  &
             dV_d_dX3, nDOFX )

    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( uCF_K(iNodeX,iCF_D ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S1,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S2,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_S3,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_E ,iZ2,iZ3,iZ4), &
                 uCF_K(iNodeX,iCF_Ne,iZ2,iZ3,iZ4), &
                 uPF_K(iPF_D ), &
                 uPF_K(iPF_V1), &
                 uPF_K(iPF_V2), &
                 uPF_K(iPF_V3), &
                 uPF_K(iPF_E ), &
                 uPF_K(iPF_Ne), &
                 GX_K (iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4), &
                 GX_K (iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        V_u_K(iNodeX,1,iZ2,iZ3,iZ4) &
          = uPF_K(iPF_V1) * WeightsX_q(iNodeX)

        V_u_K(iNodeX,2,iZ2,iZ3,iZ4) &
          = uPF_K(iPF_V2) * WeightsX_q(iNodeX)

        V_u_K(iNodeX,3,iZ2,iZ3,iZ4) &
          = uPF_K(iPF_V3) * WeightsX_q(iNodeX)

        V_d_K(iNodeX,1,iZ2,iZ3,iZ4) &
          = GX_K(iNodeX,iGF_Gm_dd_11,iZ2,iZ3,iZ4) &
              * V_u_K(iNodeX,1,iZ2,iZ3,iZ4)

        V_d_K(iNodeX,2,iZ2,iZ3,iZ4) &
          = GX_K(iNodeX,iGF_Gm_dd_22,iZ2,iZ3,iZ4) &
              * V_u_K(iNodeX,2,iZ2,iZ3,iZ4)

        V_d_K(iNodeX,3,iZ2,iZ3,iZ4) &
          = GX_K(iNodeX,iGF_Gm_dd_33,iZ2,iZ3,iZ4) &
              * V_u_K(iNodeX,3,iZ2,iZ3,iZ4)

      END DO

    END DO
    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX, - One, dLXdX3_q, nDOFX, &
             V_u_K, nDOFX, One, dV_u_dX3, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nX, nDOFX, - One, dLXdX3_q, nDOFX, &
             V_d_K, nDOFX, One, dV_d_dX3, nDOFX )

    ASSOCIATE( dZ4 => MeshX(3) % Width )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO i      = 1, 3
      DO iNodeX = 1, nDOFX

        dV_u_dX3(iNodeX,i,iZ2,iZ3,iZ4) &
         = dV_u_dX3(iNodeX,i,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dV_d_dX3(iNodeX,i,iZ2,iZ3,iZ4) &
         = dV_d_dX3(iNodeX,i,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

      END DO
      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO i      = 1, 3
      DO iNodeX = 1, nDOFX

        dV_u_dX3_Out(iNodeX,i,iZ2,iZ3,iZ4) &
          = dV_u_dX3(iNodeX,i,iZ2,iZ3,iZ4)

        dV_d_dX3_Out(iNodeX,i,iZ2,iZ3,iZ4) &
          = dV_d_dX3(iNodeX,i,iZ2,iZ3,iZ4)

      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeWeakDerivatives_X3


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


  FUNCTION FaceVelocity_X2 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R )

    REAL(DP), INTENT(in) :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in) :: V1_R, V2_R, V3_R
    REAL(DP)             :: FaceVelocity_X2(1:3)

    ! --- Average Left and Right States ---

    FaceVelocity_X2(1) = Half * ( V1_L + V1_R )
    FaceVelocity_X2(2) = Half * ( V2_L + V2_R )
    FaceVelocity_X2(3) = Half * ( V3_L + V3_R )

    RETURN
  END FUNCTION FaceVelocity_X2


  FUNCTION FaceVelocity_X3 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R )

    REAL(DP), INTENT(in) :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in) :: V1_R, V2_R, V3_R
    REAL(DP)             :: FaceVelocity_X3(1:3)

    ! --- Average Left and Right States ---

    FaceVelocity_X3(1) = Half * ( V1_L + V1_R )
    FaceVelocity_X3(2) = Half * ( V2_L + V2_R )
    FaceVelocity_X3(3) = Half * ( V3_L + V3_R )

    RETURN
  END FUNCTION FaceVelocity_X3


END MODULE TwoMoment_DiscretizationModule_Streaming_OrderV
