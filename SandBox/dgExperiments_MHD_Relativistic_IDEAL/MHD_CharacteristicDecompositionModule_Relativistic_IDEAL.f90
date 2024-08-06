MODULE MHD_CharacteristicDecompositionModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    iCM_D,  &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E,  &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi
  USE MHD_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_MHD_Relativistic
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeEnthalpyFromPrimitive, &
    ComputeMagneticEnthalpyFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeCharacteristicDecomposition_MHD_Relativistic_IDEAL

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_MHD_Relativistic_IDEAL &
    ( iDim, G, U, EvolveOnlyMagnetic, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCM)
    LOGICAL , INTENT(in)  :: EvolveOnlyMagnetic
    REAL(DP), INTENT(out) :: R(nCM,nCM)
    REAL(DP), INTENT(out) :: invR(nCM,nCM)

    REAL(DP) :: D, V1, V2, V3, E, Ne, B1, B2, B3, Chi
    REAL(DP) :: Gmdd11, Gmdd22, Gmdd33, &
                              Lapse, Shift1, Shift2, Shift3
    REAL(DP) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq, W,   &
                              bu1, bu2, bu3, bu0, bd1, bd2, bd3, bd0, &
                              h, hb, hStar

    REAL(DP) :: dFdU(nCM,nCM)

    CALL ComputePrimitive_MHD_Relativistic &
           ( U(iCM_D  ), U(iCM_S1), U(iCM_S2), &
             U(iCM_S3 ), U(iCM_E ), U(iCM_Ne), &
             U(iCM_B1 ), U(iCM_B2), U(iCM_B3), &
             U(iCM_Chi),                               &
             D, V1, V2, V3, E, Ne, &
             B1, B2, B3, Chi,      &
             G(iGF_Gm_dd_11), &
             G(iGF_Gm_dd_22), &
             G(iGF_Gm_dd_33), &
             G(iGF_Alpha   ), &
             G(iGF_Beta_1  ), &
             G(iGF_Beta_2  ), &
             G(iGF_Beta_3  ), &
             EvolveOnlyMagnetic )

    Gmdd11 = G(iGF_Gm_dd_11)
    Gmdd22 = G(iGF_Gm_dd_22)
    Gmdd33 = G(iGF_Gm_dd_33)
    Lapse  = G(iGF_Alpha   )
    Shift1 = G(iGF_Beta_1  )
    Shift2 = G(iGF_Beta_2  )
    Shift3 = G(iGF_Beta_3  )

    CALL ComputeEnthalpyFromPrimitive &
           ( D, E, Ne, h )

    CALL ComputeMagneticEnthalpyFromPrimitive &
           ( D, V1, V2, V3, B1, B2, B3, Gmdd11, Gmdd22, Gmdd33, &
             Lapse, Shift1, Shift2, Shift3, hb )

    hStar = h + hb

    Vu1 = V1
    Vu2 = V2
    Vu3 = V3

    Vd1 = Gmdd11 * Vu1
    Vd2 = Gmdd22 * Vu2
    Vd3 = Gmdd33 * Vu3

    VSq = Vd1 * Vu1 + Vd2 * Vd2 + Vd3 * Vd3

    W = One / SQRT( One - VSq )

    bu1 = B1
    bu2 = B2
    bu3 = B3
    bu0 = ( Gmdd11 * Vu1 * bu1 &
            + Gmdd22 * Vu2     &
            + Gmdd33 * Vu3 )   &
          / ( Lapse - Gmdd11 * Vu1 * Shift1 &
                    - Gmdd22 * Vu2 * Shift2 &
                    - Gmdd33 * Vu3 * Shift3 )

    bd1 = Gmdd11 * bu1
    bd2 = Gmdd22 * bu2
    bd3 = Gmdd33 * bu3
    bd0 = bu0 * ( - Lapse**2 + Gmdd11 * Shift1**2   &
                             + Gmdd22 * Shift2**2   &
                             + Gmdd33 * Shift3**2 ) &
          + ( Gmdd11 * Shift1 * bu1   &
              + Gmdd22 * Shift2 * bu2 &
              + Gmdd33 * Shift3 * bu3 )

    SELECT CASE( iDim )
      CASE(1)
        CALL ComputeFluxJacobian_X1 &
               ( D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq, &
                 bu0, bu1, bu2, bu3, bd1, bd2, bd3,    &
                 W, h, hStar, dFdU )
    END SELECT

    !CALL ComputeCharacteristicDecomposition_Numeric( R, invR, dFdU )

  END SUBROUTINE ComputeCharacteristicDecomposition_MHD_Relativistic_IDEAL


  SUBROUTINE ComputeFluxJacobian_X1 &
    ( D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq,   &
      bu0, bu1, bu2, bu3, bd1, bd2, bd3, &
      W, h, hStar, dFdU_X1 )

    REAL(DP), INTENT(in)  :: D
    REAL(DP), INTENT(in)  :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq
    REAL(DP), INTENT(in)  :: bu0, bu1, bu2, bu3, bd1, bd2, bd3
    REAL(DP), INTENT(in)  :: W, h, hStar
    REAL(DP), INTENT(out) :: dFdU_X1(nCM,nCM)

 END SUBROUTINE ComputeFluxJacobian_X1


END MODULE MHD_CharacteristicDecompositionModule_Relativistic_IDEAL
