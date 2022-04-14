MODULE GeometryComputationModule_XCFC

  USE KindModule, ONLY: &
    DP, Zero, Half, Two, Three
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    iGF_h_1,      iGF_h_2,      iGF_h_3,      &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_Beta_1,   iGF_Beta_2,   iGF_Beta_3,   &
    iGF_Alpha, &
    iGF_K_dd_11, iGF_K_dd_12, iGF_K_dd_13, &
    iGF_K_dd_22, iGF_K_dd_23, iGF_K_dd_33, &
    nGF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeChristoffelSymbols_3D_XCFC

CONTAINS


  SUBROUTINE ComputeChristoffelSymbols_3D_XCFC &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, dGdX1, dGdX2, dGdX3, &
      Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3 )

      INTEGER,  INTENT(in) :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
      REAL(DP), INTENT(in) :: G    (nDOFX      ,iX_B1(1):iX_E1(1), &
                                                iX_B1(2):iX_E1(2), &
                                                iX_B1(3):iX_E1(3),1:nGF)
      REAL(DP), INTENT(in) :: dGdX1(nDOFX,nGF  ,iX_B0(2):iX_E0(2), &
                                                iX_B0(3):iX_E0(3), &
                                                iX_B0(1):iX_E0(1))
      REAL(DP), INTENT(in) :: dGdX2(nDOFX,nGF  ,iX_B0(1):iX_E0(1), &
                                                iX_B0(3):iX_E0(3), &
                                                iX_B0(2):iX_E0(2))
      REAL(DP), INTENT(in) :: dGdX3(nDOFX,nGF  ,iX_B0(1):iX_E0(1), &
                                                iX_B0(2):iX_E0(2), &
                                                iX_B0(3):iX_E0(3))

      REAL(DP), INTENT(out) :: Christoffel3D_X1(3,3,nDOFX,iX_B0(1):iX_E0(1), &
                                                          iX_B0(2):iX_E0(2), &
                                                          iX_B0(3):iX_E0(3))
      REAL(DP), INTENT(out) :: Christoffel3D_X2(3,3,nDOFX,iX_B0(1):iX_E0(1), &
                                                          iX_B0(2):iX_E0(2), &
                                                          iX_B0(3):iX_E0(3))
      REAL(DP), INTENT(out) :: Christoffel3D_X3(3,3,nDOFX,iX_B0(1):iX_E0(1), &
                                                          iX_B0(2):iX_E0(2), &
                                                          iX_B0(3):iX_E0(3))

      INTEGER :: iNX, iX1, iX2, iX3

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        ! --- X1 ---

        Christoffel3D_X1(1,1,iNX,iX1,iX2,iX3) &
          = dGdX1(iNX,iGF_h_1,iX2,iX3,iX1) / G(iNX,iX1,iX2,iX3,iGF_h_1)

        Christoffel3D_X1(2,1,iNX,iX1,iX2,iX3) &
          = dGdX2(iNX,iGF_h_1,iX1,iX3,iX2) / G(iNX,iX1,iX2,iX3,iGF_h_1)

        Christoffel3D_X1(3,1,iNX,iX1,iX2,iX3) &
          = dGdX3(iNX,iGF_h_1,iX1,iX2,iX3) / G(iNX,iX1,iX2,iX3,iGF_h_1)

        Christoffel3D_X1(1,2,iNX,iX1,iX2,iX3) &
          = Christoffel3D_X1(2,1,iNX,iX1,iX2,iX3)

        Christoffel3D_X1(2,2,iNX,iX1,iX2,iX3) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_2) / G(iNX,iX1,iX2,iX3,iGF_h_1)**2 &
              * dGdX1(iNX,iGF_h_2,iX2,iX3,iX1)

        Christoffel3D_X1(3,2,iNX,iX1,iX2,iX3) &
          = Zero

        Christoffel3D_X1(1,3,iNX,iX1,iX2,iX3) &
          = Christoffel3D_X1(3,1,iNX,iX1,iX2,iX3)

        Christoffel3D_X1(2,3,iNX,iX1,iX2,iX3) &
          = Christoffel3D_X1(3,2,iNX,iX1,iX2,iX3)

        Christoffel3D_X1(3,3,iNX,iX1,iX2,iX3) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_3) / G(iNX,iX1,iX2,iX3,iGF_h_1)**2 &
              * dGdX1(iNX,iGF_h_3,iX2,iX3,iX1)

        ! --- X2 ---

        Christoffel3D_X2(1,1,iNX,iX1,iX2,iX3) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_1) / G(iNX,iX1,iX2,iX3,iGF_h_2)**2 &
              * dGdX2(iNX,iGF_h_1,iX1,iX3,iX2)

        Christoffel3D_X2(2,1,iNX,iX1,iX2,iX3) &
          = dGdX1(iNX,iGF_h_2,iX2,iX3,iX1) / G(iNX,iX1,iX2,iX3,iGF_h_2)

        Christoffel3D_X2(3,1,iNX,iX1,iX2,iX3) &
          = Zero

        Christoffel3D_X2(1,2,iNX,iX1,iX2,iX3) &
          = Christoffel3D_X2(2,1,iNX,iX1,iX2,iX3)

        Christoffel3D_X2(2,2,iNX,iX1,iX2,iX3) &
          = dGdX2(iNX,iGF_h_2,iX1,iX3,iX2) / G(iNX,iX1,iX2,iX3,iGF_h_2)

        Christoffel3D_X2(3,2,iNX,iX1,iX2,iX3) &
          = dGdX3(iNX,iGF_h_2,iX1,iX2,iX3) / G(iNX,iX1,iX2,iX3,iGF_h_2)

        Christoffel3D_X2(1,3,iNX,iX1,iX2,iX3) &
          = Christoffel3D_X2(3,1,iNX,iX1,iX2,iX3)

        Christoffel3D_X2(2,3,iNX,iX1,iX2,iX3) &
          = Christoffel3D_X2(3,2,iNX,iX1,iX2,iX3)

        Christoffel3D_X2(3,3,iNX,iX1,iX2,iX3) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_3) / G(iNX,iX1,iX2,iX3,iGF_h_2)**2 &
              * dGdX2(iNX,iGF_h_3,iX1,iX3,iX2)

        ! --- X3 ---

        Christoffel3D_X3(1,1,iNX,iX1,iX2,iX3) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_1) / G(iNX,iX1,iX2,iX3,iGF_h_3)**2 &
              * dGdX3(iNX,iGF_h_1,iX1,iX2,iX3)

        Christoffel3D_X3(2,1,iNX,iX1,iX2,iX3) &
          = Zero

        Christoffel3D_X3(3,1,iNX,iX1,iX2,iX3) &
          = dGdX1(iNX,iGF_h_3,iX2,iX3,iX1) / G(iNX,iX1,iX2,iX3,iGF_h_3)

        Christoffel3D_X3(1,2,iNX,iX1,iX2,iX3) &
          = Christoffel3D_X3(2,1,iNX,iX1,iX2,iX3)

        Christoffel3D_X3(2,2,iNX,iX1,iX2,iX3) &
          = -G(iNX,iX1,iX2,iX3,iGF_h_2) / G(iNX,iX1,iX2,iX3,iGF_h_3)**2 &
              * dGdX3(iNX,iGF_h_2,iX1,iX2,iX3)

        Christoffel3D_X3(3,2,iNX,iX1,iX2,iX3) &
          = dGdX2(iNX,iGF_h_3,iX1,iX3,iX2) / G(iNX,iX1,iX2,iX3,iGF_h_3)

        Christoffel3D_X3(1,3,iNX,iX1,iX2,iX3) &
          = Christoffel3D_X3(3,1,iNX,iX1,iX2,iX3)

        Christoffel3D_X3(2,3,iNX,iX1,iX2,iX3) &
          = Christoffel3D_X3(3,2,iNX,iX1,iX2,iX3)

        Christoffel3D_X3(3,3,iNX,iX1,iX2,iX3) &
          = dGdX3(iNX,iGF_h_3,iX1,iX2,iX3) / G(iNX,iX1,iX2,iX3,iGF_h_3)

      END DO
      END DO
      END DO
      END DO

  END SUBROUTINE ComputeChristoffelSymbols_3D_XCFC


  SUBROUTINE ComputeExtrinsicCurvature_XCFC &
    ( iX_B0, iX_E0, iX_B1, iX_E1, dGdX1, dGdX2, dGdX3, &
      Christoffel3D_X1, Christoffel3D_X2, Christoffel3D_X3, G )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      dGdX1(1:,1:,iX_B0(2):,iX_B0(3):,iX_B0(1):), &
      dGdX2(1:,1:,iX_B0(1):,iX_B0(3):,iX_B0(2):), &
      dGdX3(1:,1:,iX_B0(1):,iX_B0(1):,iX_B0(3):)
    REAL(DP), INTENT(in)    :: &
      Christoffel3D_X1(1:,1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):), &
      Christoffel3D_X2(1:,1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):), &
      Christoffel3D_X3(1:,1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):)
    REAL(DP), INTENT(inout) :: G(:,:,:,:,:)

    INTEGER  :: iNX, iX1, iX2, iX3
    REAL(DP) :: DivGridVolume

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      DivGridVolume &
        =   dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) &
          + dGdX2(iNX,iGF_Beta_2,iX1,iX3,iX2) &
          + dGdX3(iNX,iGF_Beta_3,iX1,iX2,iX3) &
          + (   Christoffel3D_X1(1,1,iNX,iX1,iX2,iX3) &
              + Christoffel3D_X2(2,1,iNX,iX1,iX3,iX2) &
              + Christoffel3D_X3(3,1,iNX,iX1,iX2,iX3) ) &
            * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
          + (   Christoffel3D_X1(1,2,iNX,iX1,iX2,iX3) &
              + Christoffel3D_X2(2,2,iNX,iX1,iX3,iX2) &
              + Christoffel3D_X3(3,2,iNX,iX1,iX2,iX3) ) &
            * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
          + (   Christoffel3D_X1(1,3,iNX,iX1,iX2,iX3) &
              + Christoffel3D_X2(2,3,iNX,iX1,iX3,iX2) &
              + Christoffel3D_X3(3,3,iNX,iX1,iX2,iX3) ) &
            * G(iNX,iX1,iX2,iX3,iGF_Beta_3)

      G(iNX,iX1,iX2,iX3,iGF_K_dd_11) &
        = Half / G(iNX,iX1,iX2,iX3,iGF_Alpha) &
            * (   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                    * ( dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) &
                          + Christoffel3D_X1(1,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X1(1,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X1(1,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                    * ( dGdX1(iNX,iGF_Beta_1,iX2,iX3,iX1) &
                          + Christoffel3D_X1(1,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X1(1,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X1(1,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
                - Two / Three &
                    * G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) * DivGridVolume )

      G(iNX,iX1,iX2,iX3,iGF_K_dd_12) &
        = Half / G(iNX,iX1,iX2,iX3,iGF_Alpha) &
            * (   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                    * ( dGdX1(iNX,iGF_Beta_2,iX2,iX3,iX1) &
                          + Christoffel3D_X2(1,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X2(1,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X2(1,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                    * ( dGdX2(iNX,iGF_Beta_1,iX1,iX3,iX2) &
                          + Christoffel3D_X1(2,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X1(2,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X1(2,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) )

      G(iNX,iX1,iX2,iX3,iGF_K_dd_13) &
        = Half / G(iNX,iX1,iX2,iX3,iGF_Alpha) &
            * (   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                    * ( dGdX1(iNX,iGF_Beta_3,iX2,iX3,iX1) &
                          + Christoffel3D_X3(1,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X3(1,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X3(1,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                    * ( dGdX3(iNX,iGF_Beta_1,iX1,iX2,iX3) &
                          + Christoffel3D_X1(3,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X1(3,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X1(3,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) )

      G(iNX,iX1,iX2,iX3,iGF_K_dd_22) &
        = Half / G(iNX,iX1,iX2,iX3,iGF_Alpha) &
            * (   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                    * ( dGdX2(iNX,iGF_Beta_2,iX1,iX3,iX2) &
                          + Christoffel3D_X2(2,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X2(2,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X2(2,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                    * ( dGdX2(iNX,iGF_Beta_2,iX1,iX3,iX2) &
                          + Christoffel3D_X2(2,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X2(2,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X2(2,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
                - Two / Three * G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                    * DivGridVolume )

      G(iNX,iX1,iX2,iX3,iGF_K_dd_23) &
        = Half / G(iNX,iX1,iX2,iX3,iGF_Alpha) &
            * (   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                    * ( dGdX2(iNX,iGF_Beta_3,iX1,iX3,iX2) &
                          + Christoffel3D_X3(2,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X3(2,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X3(2,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                    * ( dGdX3(iNX,iGF_Beta_2,iX1,iX2,iX3) &
                          + Christoffel3D_X2(3,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X2(3,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X2(3,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) )

      G(iNX,iX1,iX2,iX3,iGF_K_dd_33) &
        = Half / G(iNX,iX1,iX2,iX3,iGF_Alpha) &
            * (   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                    * ( dGdX3(iNX,iGF_Beta_3,iX1,iX2,iX3) &
                          + Christoffel3D_X3(3,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X3(3,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X3(3,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
                + G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                    * ( dGdX3(iNX,iGF_Beta_3,iX1,iX2,iX3) &
                          + Christoffel3D_X3(3,1,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_1) &
                          + Christoffel3D_X3(3,2,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_2) &
                          + Christoffel3D_X3(3,3,iNX,iX1,iX2,iX3) &
                              * G(iNX,iX1,iX2,iX3,iGF_Beta_3) ) &
                - Two / Three * G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                    * DivGridVolume )

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeExtrinsicCurvature_XCFC


END MODULE GeometryComputationModule_XCFC
