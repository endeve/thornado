MODULE TwoMoment_DiscretizationModule_Streaming_OrderV

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOF
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Explicit, &
    Timer_Ex_Permute, &
    Timer_Ex_Interpolate, &
    Timer_Ex_Flux, &
    Timer_Ex_Increment
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, LX_X1_Up
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, &
    Weights_q, &
    Weights_X1
  USE ReferenceElementModule_Lagrange, ONLY: &
    dLdX1_q, &
    L_X1_Dn, L_X1_Up
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE
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
    nCF
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputePrimitive_TwoMoment, &
    Flux_X1, &
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
      GX  (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER :: iNodeX, iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS

    CALL TimersStart( Timer_Explicit )

    PRINT*, "      ComputeIncrement_TwoMoment_Explicit (Start)"

    ASSOCIATE &
      ( dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOF

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL ComputeIncrement_Divergence_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    ! --- Multiply Inverse Mass Matrix ---

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOF

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
              / ( Weights_q(iNodeZ) * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm) &
                    * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dZ2, etc.

    PRINT*, "      ComputeIncrement_TwoMoment_Explicit (End)"

    PRINT*, "        dU_R = ", MINVAL( dU_R ), MAXVAL( dU_R )

    CALL TimersStop( Timer_Explicit )

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U_F (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: iNode, iNodeZ, iNodeE, iNodeX
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF
    INTEGER  :: nZ(4), nZ_X1(4), nK, nF, nF_X1
    INTEGER  :: nIterations
    REAL(DP) :: uPR_L(nPR), Flux_L(nCR)
    REAL(DP) :: uPR_R(nPR), Flux_R(nCR)
    REAL(DP) :: uPR_K(nPR), Flux_K(nCR)
    REAL(DP) :: GX_K(nDOFX,   iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                              iZ_B1(2):iZ_E1(2),nGF)
    REAL(DP) :: GX_F(nDOFX_X1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                           iZ_B0(2):iZ_E1(2),nGF)
    REAL(DP) :: G_K(nDOF,nGF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                    iZ_B0(2):iZ_E0(2)  )
    REAL(DP) :: G_F(nDOF_X1,nGF,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                                iZ_B0(2):iZ_E1(2))
    REAL(DP) :: uCR_K(nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                      iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B1(2):iZ_E1(2))
    REAL(DP) :: uCR_L(nDOF_X1,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                      iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B0(2):iZ_E1(2))
    REAL(DP) :: uCR_R(nDOF_X1,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                      iZ_B0(4):iZ_E0(4),nCR,nSpecies,iZ_B0(2):iZ_E1(2))
    REAL(DP) :: NumericalFlux(nDOF_X1,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                      iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(2):iZ_E1(2))
    REAL(DP) :: dU_X1(nDOF,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                      iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(2):iZ_E0(2)  )

    REAL(DP) :: Flux_q(nDOF,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(3):iZ_E0(3), &
                       iZ_B0(4):iZ_E0(4),nSpecies,iZ_B0(2):iZ_E0(2)  )

    PRINT*, "      ComputeIncrement_Divergence_X1"

    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN

    nZ    = iZ_E0 - iZ_B0 + 1
    nZ_X1 = nZ + [0,1,0,0]
    nK    = nSpecies * nCR * PRODUCT( nZ )
    nF    = nSpecies * nCR * PRODUCT( nZ_X1 )
    nF_X1 = PRODUCT( nZ_X1(2:4) )

    ASSOCIATE &
      ( dZ3 => MeshX(2) % Width, &
        dZ4 => MeshX(3) % Width )

    CALL TimersStart( Timer_Ex_Permute )

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

    CALL TimersStop( Timer_Ex_Permute )

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Ex_Interpolate )

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    ! --- Scale Factors ---

    DO iGF = iGF_h_1, iGF_h_3

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nF_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1,iGF), nDOFX, Zero, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nF_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               GX_K(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX, Half, &
               GX_F(1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ,iGF), nDOFX_X1 )

    END DO

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Permute )

    DO iGF = iGF_h_1, iGF_h_3
    DO iZ2 = iZ_B0(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1)*nDOFE + iNodeE

        G_F(iNodeZ,iGF,iZ3,iZ4,iZ2) &
          = MAX( GX_F(iNodeX,iZ3,iZ4,iZ2,iGF), SqrtTiny )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO

    DO iZ2 = iZ_B0(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNode = 1, nDOF_X1

        G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2) &
          = MAX( G_F(iNode,iGF_h_1,iZ3,iZ4,iZ2)**2, SqrtTiny )
        G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2) &
          = MAX( G_F(iNode,iGF_h_2,iZ3,iZ4,iZ2)**2, SqrtTiny )
        G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) &
          = MAX( G_F(iNode,iGF_h_3,iZ3,iZ4,iZ2)**2, SqrtTiny )
        G_F(iNode,iGF_SqrtGm,iZ3,iZ4,iZ2) &
          =   G_F(iNode,iGF_h_1,iZ3,iZ4,iZ2) &
            * G_F(iNode,iGF_h_2,iZ3,iZ4,iZ2) &
            * G_F(iNode,iGF_h_3,iZ3,iZ4,iZ2)

      END DO

    END DO
    END DO
    END DO

    ! --- Permute Radiation Fields ---

    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOF

        uCR_K(iNodeZ,iZ1,iZ3,iZ4,iCR,iS,iZ2) &
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
           ( 'N', 'N', nDOF_X1, nF, nDOF, One, L_X1_Up, nDOF_X1, &
             uCR_K(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)-1), nDOF, Zero, &
             uCR_L(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)  ), nDOF_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOF_X1, nF, nDOF, One, L_X1_Dn, nDOF_X1, &
             uCR_K(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)  ), nDOF, Zero, &
             uCR_R(1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)  ), nDOF_X1 )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Flux )

    ! --- Numerical Flux ---

    DO iZ2 = iZ_B0(2), iZ_E1(2)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNode = 1, nDOF_X1

        ! --- Left State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_L(iNode,iZ1,iZ3,iZ4,iCR_N, iS,iZ2), &
                 uCR_L(iNode,iZ1,iZ3,iZ4,iCR_G1,iS,iZ2), &
                 uCR_L(iNode,iZ1,iZ3,iZ4,iCR_G2,iS,iZ2), &
                 uCR_L(iNode,iZ1,iZ3,iZ4,iCR_G3,iS,iZ2), &
                 uPR_L(iPR_D ), uPR_L(iPR_I1), &
                 uPR_L(iPR_I2), uPR_L(iPR_I3), &
                 Zero, Zero, Zero, &
                 G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 nIterations )

        ! --- Left State Flux ---

        Flux_L = Flux_X1 &
                   ( uPR_L(iPR_D ), uPR_L(iPR_I1), &
                     uPR_L(iPR_I2), uPR_L(iPR_I3), &
                     Zero, Zero, Zero, &
                     G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                     G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                     G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

        ! --- Right State Primitive ---

        CALL ComputePrimitive_TwoMoment &
               ( uCR_R(iNode,iZ1,iZ3,iZ4,iCR_N, iS,iZ2), &
                 uCR_R(iNode,iZ1,iZ3,iZ4,iCR_G1,iS,iZ2), &
                 uCR_R(iNode,iZ1,iZ3,iZ4,iCR_G2,iS,iZ2), &
                 uCR_R(iNode,iZ1,iZ3,iZ4,iCR_G3,iS,iZ2), &
                 uPR_R(iPR_D ), uPR_R(iPR_I1), &
                 uPR_R(iPR_I2), uPR_R(iPR_I3), &
                 Zero, Zero, Zero, &
                 G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 nIterations )

        ! --- Right State Flux ---

        Flux_R = Flux_X1 &
                   ( uPR_R(iPR_D ), uPR_R(iPR_I1), &
                     uPR_R(iPR_I2), uPR_R(iPR_I3), &
                     Zero, Zero, Zero, &
                     G_F(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                     G_F(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                     G_F(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

        DO iCR = 1, nCR

          NumericalFlux(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = NumericalFlux_LLF &
                ( uCR_L(iNode,iZ1,iZ3,iZ4,iCR,iS,iZ2), &
                  uCR_R(iNode,iZ1,iZ3,iZ4,iCR,iS,iZ2), &
                  Flux_L(iCR), Flux_R(iCR), One )

          NumericalFlux(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = dZ3(iZ3) * dZ4(iZ4) * Weights_X1(iNode) &
                * G_F(iNode,iGF_SqrtGm,iZ3,iZ4,iZ2) &
                * NumericalFlux(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2)

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
           ( 'T', 'N', nDOF, nK, nDOF_X1, + One, L_X1_Dn, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)  ), &
             nDOF_X1, Zero, dU_X1, nDOF )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOF, nK, nDOF_X1, - One, L_X1_Up, nDOF_X1, &
             NumericalFlux(1,1,iZ_B0(1),iZ_B0(3),iZ_B0(4),1,iZ_B0(2)+1), &
             nDOF_X1, One,  dU_X1, nDOF )

    CALL TimersStop( Timer_Ex_Interpolate )

    !--------------------
    ! --- Volume Term ---
    !--------------------

    CALL TimersStart( Timer_Ex_Permute )

    DO iGF = 1, nGF
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        G_K(iNodeZ,iGF,iZ3,iZ4,iZ2) = GX_K(iNodeX,iZ3,iZ4,iZ2,iGF)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Permute )

    CALL TimersStart( Timer_Ex_Flux )

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNode = 1, nDOF

        CALL ComputePrimitive_TwoMoment &
               ( uCR_K(iNode,iZ1,iZ3,iZ4,iCR_N, iS,iZ2), &
                 uCR_K(iNode,iZ1,iZ3,iZ4,iCR_G1,iS,iZ2), &
                 uCR_K(iNode,iZ1,iZ3,iZ4,iCR_G2,iS,iZ2), &
                 uCR_K(iNode,iZ1,iZ3,iZ4,iCR_G3,iS,iZ2), &
                 uPR_K(iPR_D ), uPR_K(iPR_I1), &
                 uPR_K(iPR_I2), uPR_K(iPR_I3), &
                 Zero, Zero, Zero, &
                 G_K(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                 G_K(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                 G_K(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2), &
                 nIterations )

        Flux_K &
          = Flux_X1 &
              ( uPR_K(iPR_D ), uPR_K(iPR_I1), &
                uPR_K(iPR_I2), uPR_K(iPR_I3), &
                Zero, Zero, Zero, &
                G_K(iNode,iGF_Gm_dd_11,iZ3,iZ4,iZ2), &
                G_K(iNode,iGF_Gm_dd_22,iZ3,iZ4,iZ2), &
                G_K(iNode,iGF_Gm_dd_33,iZ3,iZ4,iZ2) )

        DO iCR = 1, nCR

          Flux_q(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2) &
            = dZ3(iZ3) * dZ4(iZ4) * Weights_q(iNode) &
                * G_K(iNode,iGF_SqrtGm,iZ3,iZ4,iZ2) &
                * Flux_K(iCR)

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
           ( 'T', 'N', nDOF, nK, nDOF, One, dLdX1_q, nDOF, &
             Flux_q, nDOF, One, dU_X1, nDOF )

    CALL TimersStop( Timer_Ex_Interpolate )

    CALL TimersStart( Timer_Ex_Increment )

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNode = 1, nDOF

        dU_R(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = dU_R(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + dU_X1(iNode,iCR,iZ1,iZ3,iZ4,iS,iZ2)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Ex_Increment )

    END ASSOCIATE ! dZ3, etc.

  END SUBROUTINE ComputeIncrement_Divergence_X1


END MODULE TwoMoment_DiscretizationModule_Streaming_OrderV
