MODULE MHD_CharacteristicDecompositionModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One,  &
    Two
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
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeCharacteristicDecomposition_MHD_Relativistic_IDEAL

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_MHD_Relativistic_IDEAL &
    ( iDim, G, U, EvolveOnlyMagnetic, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(8)
    REAL(DP), INTENT(in)  :: U(nCM)
    LOGICAL , INTENT(in)  :: EvolveOnlyMagnetic
    REAL(DP), INTENT(out) :: R(nCM,nCM)
    REAL(DP), INTENT(out) :: invR(nCM,nCM)

    REAL(DP) :: D, V1, V2, V3, E, Ne, B1, B2, B3, Chi
    REAL(DP) :: Gmdd11, Gmdd22, Gmdd33, DetGm, &
                              Lapse, Shift1, Shift2, Shift3
    REAL(DP) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq, W,   &
                              bu1, bu2, bu3, bu0, bd1, bd2, bd3, bd0, &
                              h, hb, hStar, Ye

    REAL(DP) :: dUdV(nCM,nCM)
    REAL(DP) :: dFdV(nCM,nCM)
    REAL(DP) :: dFdU(nCM,nCM)

    CALL ComputePrimitive_MHD_Relativistic &
           ( U(iCM_D  ), U(iCM_S1), U(iCM_S2), &
             U(iCM_S3 ), U(iCM_E ), U(iCM_Ne), &
             U(iCM_B1 ), U(iCM_B2), U(iCM_B3), &
             U(iCM_Chi),                               &
             D, V1, V2, V3, E, Ne, &
             B1, B2, B3, Chi,      &
             G(1), &
             G(2), &
             G(3), &
             G(5), &
             G(6), &
             G(7), &
             G(8), &
             EvolveOnlyMagnetic )

    Gmdd11 = G(1)
    Gmdd22 = G(2)
    Gmdd33 = G(3)
    DetGm  = G(4)
    Lapse  = G(5)
    Shift1 = G(6)
    Shift2 = G(7)
    Shift3 = G(8)

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

    Ye = Ne / D

    CALL ComputePrimitiveConservedJacobian &
           ( D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq, &
             bu0, bu1, bu2, bu3, bd1, bd2, bd3,    &
             Gmdd11, Gmdd22, Gmdd33,               &
             W, h, hStar, Ye, dUdV )

    SELECT CASE( iDim )
      CASE(1)
        CALL ComputePrimitiveFluxJacobian_X1 &
               ( D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq, &
                 bu0, bu1, bu2, bu3, bd1, bd2, bd3,    &
                 Gmdd11, Gmdd22, Gmdd33,               &
                 W, h, hStar, dFdV )

    END SELECT

    CALL ComputeConservedFluxJacobian_Numeric &
           ( dUdV, dFdV, dFdU )

    CALL ComputeCharacteristicDecomposition_Numeric( R, invR, dFdU )

  END SUBROUTINE ComputeCharacteristicDecomposition_MHD_Relativistic_IDEAL


  SUBROUTINE ComputeCharacteristicDecomposition_Numeric( R, invR, dFdV )

    REAL(DP), INTENT(in)  :: dFdV(nCM,nCM)
    REAL(DP), INTENT(out) :: R(nCM,nCM)
    REAL(DP), INTENT(out) :: invR(nCM,nCM)
    REAL(DP)              :: Lambda(nCM,nCM)

    INTEGER  :: i, INFO, LWORK
    INTEGER  :: IPIV(nCM)

    REAL(DP)              :: WR(nCM)
    REAL(DP)              :: WI(nCM)
    REAL(DP)              :: TEMP(1)
    REAL(DP)              :: dFdV_Copy(nCM,nCM)
    REAL(DP)              :: invR_Copy(nCM,nCM)
    REAL(DP), ALLOCATABLE :: WORK1(:), WORK2(:)

    ! --- Copy to avoid overwriting dFdV ---

    dFdV_Copy = dFdV

    ! --- Necessary workplace query to get LWORK. ----

    CALL DGEEV( 'V', 'N', nCM, dFdV_Copy, nCM, WR, &
                WI, invR, nCM, 0, nCM, TEMP, &
                -1, INFO )

    LWORK = TEMP(1)
    ALLOCATE( WORK1(LWORK) )

    Lambda(:,:) = Zero

    ! --- Actual computation of eigedecomposition. ---

    CALL DGEEV( 'V', 'N', nCM, dFdV_Copy, nCM, WR, &
                WI, invR, nCM, 0, nCM, WORK1,  &
                LWORK, INFO )

    invR = TRANSPOSE( invR )

    invR_Copy = invR

    CALL DGETRF( nCM, nCM, invR_Copy, nCM, IPIV, INFO )

    LWORK = -1
    CALL DGETRI( nCM, invR_Copy, nCM, IPIV, TEMP, LWORK, INFO )

    LWORK = TEMP(1)
    ALLOCATE( WORK2(LWORK) )

    CALL DGETRI( nCM, invR_Copy, nCM, IPIV, WORK2, LWORK, INFO )

    R = invR_Copy

    IF ( ( INFO .NE. 0 ) .OR. ( ANY( ABS( WI ) > 1d-15 ) ) ) THEN

      PRINT*, 'INFO: ', INFO
      PRINT*, 'WR: ', WR
      PRINT*, 'WI: ', WI

      DO i = 1, nCM

        PRINT*, 'Lambda(i,:) : ', Lambda(i,:)

      END DO

    END IF

  END SUBROUTINE ComputeCharacteristicDecomposition_Numeric


  SUBROUTINE ComputePrimitiveConservedJacobian &
    ( D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq,   &
      bu0, bu1, bu2, bu3, bd1, bd2, bd3,      &
      Gmdd11, Gmdd22, Gmdd33,                 &
      W, h, hStar, Ye, dUdV )

    REAL(DP), INTENT(in)  :: D
    REAL(DP), INTENT(in)  :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq
    REAL(DP), INTENT(in)  :: bu0, bu1, bu2, bu3, bd1, bd2, bd3
    REAL(DP), INTENT(in)  :: Gmdd11, Gmdd22, Gmdd33
    REAL(DP), INTENT(in)  :: W, h, hStar, Ye
    REAL(DP), INTENT(out) :: dUdV(nCM,nCM)

    ! --- Upper 1st Index, Lower 2nd Index ---

    REAL(DP) :: b0b1, b0b2, b0b3
    REAL(DP) :: b1b1, b1b2, b1b3
    REAL(DP) :: b2b1, b2b2, b2b3
    REAL(DP) :: b3b1, b3b2, b3b3

    ! --- Upper 1st Index, Lower 2nd Index ---

    REAL(DP) :: v1v1, v1v2, v1v3
    REAL(DP) :: v2v1, v2v2, v2v3
    REAL(DP) :: v3v1, v3v2, v3v3

    ! --- Upper 1st Index, Upper 2nd Index ---

    REAL(DP) :: v1b1, v1b2, v1b3
    REAL(DP) :: v2b1, v2b2, v2b3
    REAL(DP) :: v3b1, v3b2, v3b3

    ! --- Upper 1st Index, Lower 2nd Index

    REAL(DP) :: b0v1, b0v2, b0v3

    ! --- Upper 1st Index, Lower 2nd Index ---

    REAL(DP) :: del11, del12, del13
    REAL(DP) :: del21, del22, del23
    REAL(DP) :: del31, del32, del33

    ! --- Lower 1st Index, Lower 2nd Index ---

    REAL(DP) :: gm11, gm12, gm13
    REAL(DP) :: gm21, gm22, gm23
    REAL(DP) :: gm31, gm32, gm33

    REAL(DP) :: Z, G

    b0b1 = bu0 * bd1
    b0b2 = bu0 * bd2
    b0b3 = bu0 * bd3

    b1b1 = bu1 * bd1
    b1b2 = bu1 * bd2
    b1b3 = bu1 * bd3
    b2b1 = bu2 * bd1
    b2b2 = bu2 * bd2
    b2b3 = bu2 * bd3
    b3b1 = bu3 * bd1
    b3b2 = bu3 * bd2
    b3b3 = bu3 * bd3

    v1v1 = vu1 * vd1
    v1v2 = vu1 * vd2
    v1v3 = vu1 * vd3
    v2v1 = vu2 * vd1
    v2v2 = vu2 * vd2
    v2v3 = vu2 * vd3
    v3v1 = vu3 * vd1
    v3v2 = vu3 * vd2
    v3v3 = vu3 * vd3

    v1b1 = vu1 * bu1
    v1b2 = vu1 * bu2
    v1b3 = vu1 * bu3
    v2b1 = vu2 * bu1
    v2b2 = vu2 * bu2
    v2b3 = vu2 * bu3
    v3b1 = vu3 * bu1
    v3b2 = vu3 * bu2
    v3b3 = vu3 * bu3

    del11 = One
    del12 = Zero
    del13 = Zero
    del21 = Zero
    del22 = One
    del23 = Zero
    del31 = Zero
    del32 = Zero
    del33 = One

    gm11 = Gmdd11
    gm12 = Zero
    gm13 = Zero
    gm21 = Zero
    gm22 = Gmdd22
    gm23 = Zero
    gm31 = Zero
    gm32 = Zero
    gm33 = Gmdd33

    ! --- Helper Expressions --

    Z = D * hStar * W**2

    G = Gamma_IDEAL / ( Gamma_IDEAL - One )

    ! --- Temp ---

    dUdV(1,:) &
      = [ W,              &
          D * W**3 * Vd1, &
          D * W**3 * Vd2, &
          D * W**3 * Vd3, &
          Zero,           &
          Zero,           &
          Zero,           &
          Zero,           &
          Zero,           &
          Zero ]

    dUdV(2,:) &
      = [ W**2 * vd1,                                          &
          - Two * b0b1 * W**2 * vd1 &
          + Z * ( Two * W**2 * vd1 * vd1 + gm11 ) - bd1 * bd1, &
          - Two * b0b2 * W**2 * vd1 &
          + Z * ( Two * W**2 * vd2 * vd1 + gm12 ) - bd2 * bd1, &
          - Two * b0b3 * W**2 * vd1 &
          + Z * ( Two * W**2 * vd3 * vd1 + gm13 ) - bd3 * bd1, &
          G * W**2 * vd1,                                      &
          Zero,                                                &
          ( - b0v1 + bd1 ) * ( Two * W**2 * vd1 ) &
          - bd1 * vd1 - gm11 * bu0,                            &
          ( - b0v2 + bd2 ) * ( Two * W**2 * vd1 ) &
          - bd1 * vd2 - gm12 * bu0,                            &
          ( - b0v3 + bd3 ) * ( Two * W**2 * vd1 ) &
          - bd1 * vd3 - gm13 * bu0,                            &
          Zero ]

    dUdV(3,:) &
      = [ W**2 * vd2,                                          &
          - Two * b0b1 * W**2 * vd2 &
          + Z * ( Two * W**2 * vd1 * vd2 + gm21 ) - bd1 * bd2, &
          - Two * b0b2 * W**2 * vd2 &
          + Z * ( Two * W**2 * vd2 * vd2 + gm22 ) - bd2 * bd2, &
          - Two * b0b3 * W**2 * vd2 &
          + Z * ( Two * W**2 * vd3 * vd2 + gm23 ) - bd3 * bd2, &
          G * W**2 * vd2,                                      &
          Zero,                                                &
          ( - b0v1 + bd1 ) * ( Two * W**2 * vd2 ) &
          - bd2 * vd1 - gm21 * bu0,                            &
          ( - b0v2 + bd2 ) * ( Two * W**2 * vd2 ) &
          - bd2 * vd2 - gm22 * bu0,                            &
          ( - b0v3 + bd3 ) * ( Two * W**2 * vd2 ) &
          - bd2 * vd3 - gm23 * bu0,                            &
          Zero ]

    dUdV(4,:) &
      = [ W**2 * vd3,                                          &
          - Two * b0b1 * W**2 * vd3 &
          + Z * ( Two * W**2 * vd1 * vd3 + gm31 ) - bd1 * bd3, &
          - Two * b0b2 * W**2 * vd3 &
          + Z * ( Two * W**2 * vd2 * vd3 + gm32 ) - bd2 * bd3, &
          - Two * b0b3 * W**2 * vd3 &
          + Z * ( Two * W**2 * vd3 * vd3 + gm33 ) - bd3 * bd3, &
          G * W**2 * vd3,                                      &
          Zero,                                                &
          ( - b0v1 + bd1 ) * ( Two * W**2 * vd3 ) &
          - bd3 * vd1 - gm31 * bu0,                            &
          ( - b0v2 + bd2 ) * ( Two * W**2 * vd3 ) &
          - bd3 * vd2 - gm32 * bu0,                            &
          ( - b0v3 + bd3 ) * ( Two * W**2 * vd3 ) &
          - bd3 * vd3 - gm33 * bu0,                            &
          Zero ]

    dUdV(5,:) &
      = [ ( W**2 - W ),                                             &
          Z * ( Two * W**2 * vd1 ) - D * W**3 * vd1 &
          - b0b1 * ( One + Two * W**2 ),                            &
          Z * ( Two * W**2 * vd2 ) - D * W**3 * vd2 &
          - b0b2 * ( One + Two * W**2 ),                            &
          Z * ( Two * W**2 * vd3 ) - D * W**3 * vd3 &
          - b0b3 * ( One + Two * W**2 ),                            &
          G * W**2 - One,                                           &
          Zero,                                                     &
          ( Two * W**2 - One ) * bd1 - ( Two * W**2 + One ) * b0v1, &
          ( Two * W**2 - One ) * bd2 - ( Two * W**2 + One ) * b0v2, &
          ( Two * W**2 - One ) * bd3 - ( Two * W**2 + One ) * b0v3, &
          Zero ]

    dUdV(6,:) &
      = [ Ye * W,              &
          D * Ye * W**3 * vd1, &
          D * Ye * W**3 * vd2, &
          D * Ye * W**3 * vd3, &
          Zero,                &
          D * W,               &
          Zero,                &
          Zero,                &
          Zero,                &
          Zero ]

    dUdV(7,:) &
      = [ Zero,                                                               &
          W**3 * vd1 * ( bu1 - bu0 * vu1 ) - W * ( bd1 * vu1 + bu0 * del11 ), &
          W**3 * vd2 * ( bu1 - bu0 * vu1 ) - W * ( bd2 * vu1 + bu0 * del12 ), &
          W**3 * vd3 * ( bu1 - bu0 * vu1 ) - W * ( bd3 * vu1 + bu0 * del13 ), &
          Zero,                                                               &
          Zero,                                                               &
          W * ( del11 - v1v1 ),                                               &
          W * ( del12 - v1v2 ),                                               &
          W * ( del13 - v1v3 ),                                               &
          Zero ]

    dUdV(8,:) &
      = [ Zero,                                                               &
          W**3 * vd1 * ( bu2 - bu0 * vu2 ) - W * ( bd1 * vu2 + bu0 * del21 ), &
          W**3 * vd2 * ( bu2 - bu0 * vu2 ) - W * ( bd2 * vu2 + bu0 * del22 ), &
          W**3 * vd3 * ( bu2 - bu0 * vu2 ) - W * ( bd3 * vu2 + bu0 * del23 ), &
          Zero,                                                               &
          Zero,                                                               &
          W * ( del21 - v2v1 ),                                               &
          W * ( del22 - v2v2 ),                                               &
          W * ( del23 - v2v3 ),                                               &
          Zero ]

    dUdV(9,:) &
      = [ Zero,                                                               &
          W**3 * vd1 * ( bu3 - bu0 * vu3 ) - W * ( bd1 * vu3 + bu0 * del31 ), &
          W**3 * vd2 * ( bu3 - bu0 * vu3 ) - W * ( bd2 * vu3 + bu0 * del32 ), &
          W**3 * vd3 * ( bu3 - bu0 * vu3 ) - W * ( bd3 * vu3 + bu0 * del33 ), &
          Zero,                                                               &
          Zero,                                                               &
          W * ( del31 - v3v1 ),                                               &
          W * ( del32 - v3v2 ),                                               &
          W * ( del33 - v3v3 ),                                               &
          Zero ]

    dUdV(10,:) &
      = [ Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          One ]

    ! --- ---

  END SUBROUTINE ComputePrimitiveConservedJacobian


  SUBROUTINE ComputePrimitiveFluxJacobian_X1 &
    ( D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq,   &
      bu0, bu1, bu2, bu3, bd1, bd2, bd3,      &
      Gmdd11, Gmdd22, Gmdd33,                 &
      W, h, hStar, dFdV_X1 )

    REAL(DP), INTENT(in)  :: D
    REAL(DP), INTENT(in)  :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq
    REAL(DP), INTENT(in)  :: bu0, bu1, bu2, bu3, bd1, bd2, bd3
    REAL(DP), INTENT(in)  :: Gmdd11, Gmdd22, Gmdd33
    REAL(DP), INTENT(in)  :: W, h, hStar
    REAL(DP), INTENT(out) :: dFdV_X1(nCM,nCM)

    ! --- Upper 1st Index, Lower 2nd Index ---

    REAL(DP) :: b0b1, b0b2, b0b3
    REAL(DP) :: b1b1, b1b2, b1b3
    REAL(DP) :: b2b1, b2b2, b2b3
    REAL(DP) :: b3b1, b3b2, b3b3

    ! --- Upper 1st Index, Lower 2nd Index ---

    REAL(DP) :: v1v1, v1v2, v1v3
    REAL(DP) :: v2v1, v2v2, v2v3
    REAL(DP) :: v3v1, v3v2, v3v3

    ! --- Upper 1st Index, Upper 2nd Index ---

    REAL(DP) :: v1b1, v1b2, v1b3
    REAL(DP) :: v2b1, v2b2, v2b3
    REAL(DP) :: v3b1, v3b2, v3b3

    ! --- Upper 1st Index, Lower 2nd Index

    REAL(DP) :: b0v1, b0v2, b0v3

    ! --- Upper 1st Index, Lower 2nd Index ---

    REAL(DP) :: del11, del12, del13
    REAL(DP) :: del21, del22, del23
    REAL(DP) :: del31, del32, del33

    ! --- Lower 1st Index, Lower 2nd Index ---

    REAL(DP) :: gm11, gm12, gm13
    REAL(DP) :: gm21, gm22, gm23
    REAL(DP) :: gm31, gm32, gm33

    REAL(DP) :: Z, G

    b0b1 = bu0 * bd1
    b0b2 = bu0 * bd2
    b0b3 = bu0 * bd3

    b1b1 = bu1 * bd1
    b1b2 = bu1 * bd2
    b1b3 = bu1 * bd3
    b2b1 = bu2 * bd1
    b2b2 = bu2 * bd2
    b2b3 = bu2 * bd3
    b3b1 = bu3 * bd1
    b3b2 = bu3 * bd2
    b3b3 = bu3 * bd3

    v1v1 = vu1 * vd1
    v1v2 = vu1 * vd2
    v1v3 = vu1 * vd3
    v2v1 = vu2 * vd1
    v2v2 = vu2 * vd2
    v2v3 = vu2 * vd3
    v3v1 = vu3 * vd1
    v3v2 = vu3 * vd2
    v3v3 = vu3 * vd3

    v1b1 = vu1 * bu1
    v1b2 = vu1 * bu2
    v1b3 = vu1 * bu3
    v2b1 = vu2 * bu1
    v2b2 = vu2 * bu2
    v2b3 = vu2 * bu3
    v3b1 = vu3 * bu1
    v3b2 = vu3 * bu2
    v3b3 = vu3 * bu3

    del11 = One
    del12 = Zero
    del13 = Zero
    del21 = Zero
    del22 = One
    del23 = Zero
    del31 = Zero
    del32 = Zero
    del33 = One

    gm11 = Gmdd11
    gm12 = Zero
    gm13 = Zero
    gm21 = Zero
    gm22 = Gmdd22
    gm23 = Zero
    gm31 = Zero
    gm32 = Zero
    gm33 = Gmdd33

    ! --- Helper Expressions --

    Z = D * hStar * W**2

    G = Gamma_IDEAL / ( Gamma_IDEAL - One )

    ! --- Temp ---

    dFdV_X1(1,:) &
      = [ W * Vu1,                      &
          D * W**3 * Vd1 * Vu1 + D * W, &
          D * W**3 * Vd2 * Vu1,         &
          D * W**3 * Vd3 * Vu1,         &
          Zero,                         &
          Zero,                         &
          Zero,                         &
          Zero,                         &
          Zero,                         &
          Zero ]

    dFdV_X1(2,:) &
      = [ W**2 * v1v1,                                                  &
          -b0b1 * ( del11 + Two * W**2 * v1v1 ) &
          + Z * ( Two * W**2 * vd1 * v1v1 + del11 * vd1 + gm11 * vu1 ), &
          -b0b2 * ( del11 + Two * W**2 * v1v1 ) &
          + Z * ( Two * W**2 * vd2 * v1v1 + del12 * vd1 + gm12 * vu1 ), &
          -b0b3 * ( del11 + Two * W**2 * v1v1 ) &
          + Z * ( Two * W**2 * vd3 * v1v1 + del13 * vd1 + gm13 * vu1 ), &
          del11 + G * W**2 * v1v1,                                      &
          Zero,                                                         &
          ( -b0v1 + bd1 ) * ( del11 + Two * W**2 * v1v1 ) &
          - bd1 * del11 - gm11 * bu1,                                   &
          ( -b0v2 + bd2 ) * ( del11 + Two * W**2 * v1v1 ) &
          - bd1 * del12 - gm12 * bu1,                                   &
          ( -b0v3 + bd3 ) * ( del11 + Two * W**2 * v1v1 ) &
          - bd1 * del13 - gm13 * bu1,                                   &
          Zero ]

    dFdV_X1(3,:) &
      = [ W**2 * v1v2,                                                  &
          -b0b1 * ( del12 + Two * W**2 * v1v2 ) &
          + Z * ( Two * W**2 * vd1 * v1v2 + del11 * vd2 + gm21 * vu1 ), &
          -b0b2 * ( del12 + Two * W**2 * v1v2 ) &
          + Z * ( Two * W**2 * vd2 * v1v2 + del12 * vd2 + gm22 * vu1 ), &
          -b0b3 * ( del12 + Two * W**2 * v1v2 ) &
          + Z * ( Two * W**2 * vd3 * v1v2 + del13 * vd2 + gm23 * vu1 ), &
          del12 + G * W**2 * v1v2,                                      &
          Zero,                                                         &
          ( -b0v1 + bd1 ) * ( del12 + Two * W**2 * v1v2 ) &
          - bd2 * del11 - gm21 * bu1,                                   &
          ( -b0v2 + bd2 ) * ( del12 + Two * W**2 * v1v2 ) &
          - bd2 * del12 - gm22 * bu1,                                   &
          ( -b0v3 + bd3 ) * ( del12 + Two * W**2 * v1v2 ) &
          - bd2 * del13 - gm23 * bu1,                                   &
          Zero ]

    dFdV_X1(4,:) &
      = [ W**2 * v1v1,                                                  &
          -b0b1 * ( del13 + Two * W**2 * v1v3 ) &
          + Z * ( Two * W**2 * vd1 * v1v3 + del11 * vd3 + gm31 * vu1 ), &
          -b0b2 * ( del13 + Two * W**2 * v1v3 ) &
          + Z * ( Two * W**2 * vd2 * v1v3 + del12 * vd3 + gm32 * vu1 ), &
          -b0b3 * ( del13 + Two * W**2 * v1v3 ) &
          + Z * ( Two * W**2 * vd3 * v1v3 + del13 * vd3 + gm33 * vu1 ), &
          del13 + G * W**2 * v1v3,                                      &
          Zero,                                                         &
          ( -b0v1 + bd1 ) * ( del13 + Two * W**2 * v1v3 ) &
          - bd3 * del11 - gm31 * bu1,                                   &
          ( -b0v2 + bd2 ) * ( del13 + Two * W**2 * v1v3 ) &
          - bd3 * del12 - gm32 * bu1,                                   &
          ( -b0v3 + bd3 ) * ( del13 + Two * W**2 * v1v3 ) &
          - bd3 * del13 - gm33 * bu1,                                   &
          Zero ]

    dFdV_X1(5,:) &
      = [ ( W**2 - W ) * vu1,            &
          Z * ( Two * W**2 * v1v1 + del11 ) - Two * W**2 * b0b1 * vu1 - b1b1 &
          - D * W**3 * v1v1 - D * W * del11, &
          Z * ( Two * W**2 * v1v2 + del12 ) - Two * W**2 * b0b2 * vu1 - b1b2 &
          - D * W**3 * v1v2 - D * W * del12, &
          Z * ( Two * W**2 * v1v3 + del13 ) - Two * W**2 * b0b3 * vu1 - b1b3 &
          - D * W**3 * v1v3 - D * W * del13, &
          G * W**2 * vu1,            &
          Zero,            &
          Two * ( -b0v1 + bd1 ) * W**2 * vu1 - bu1 * vd1 - bu0 * del11, &
          Two * ( -b0v2 + bd2 ) * W**2 * vu1 - bu1 * vd2 - bu0 * del12, &
          Two * ( -b0v3 + bd3 ) * W**2 * vu1 - bu1 * vd3 - bu0 * del13, &
          Zero ]

    dFdV_X1(6,:) &
      = [ Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero ]

    dFdV_X1(7,:) &
      = [ Zero,                                                             &
          W**3 * vd1 * ( v1b1 - v1b1 ) + W * ( del11 * bu1 - del11 * bu1 ), &
          W**3 * vd2 * ( v1b1 - v1b1 ) + W * ( del12 * bu1 - del12 * bu1 ), &
          W**3 * vd3 * ( v1b1 - v1b1 ) + W * ( del13 * bu1 - del13 * bu1 ), &
          Zero,                                                             &
          Zero,                                                             &
          W * ( vu1 * del11 - vu1 * del11 ),                                &
          W * ( vu1 * del12 - vu1 * del12 ),                                &
          W * ( vu1 * del13 - vu1 * del13 ),                                &
          Zero ]

    dFdV_X1(8,:) &
      = [ Zero,            &
          W**3 * vd1 * ( v1b2 - v2b1 ) + W * ( del11 * bu2 - del21 * bu1 ), &
          W**3 * vd2 * ( v1b2 - v2b1 ) + W * ( del12 * bu2 - del22 * bu1 ), &
          W**3 * vd3 * ( v1b2 - v2b1 ) + W * ( del13 * bu2 - del23 * bu1 ), &
          Zero,            &
          Zero,            &
          W * ( vu1 * del21 - vu2 * del11 ),            &
          W * ( vu1 * del22 - vu2 * del12 ),            &
          W * ( vu1 * del23 - vu2 * del13 ),            &
          Zero ]

    dFdV_X1(9,:) &
      = [ Zero,            &
          W**3 * vd1 * ( v1b3 - v3b1 ) + W * ( del11 * bu3 - del31 * bu1 ), &
          W**3 * vd2 * ( v1b3 - v3b1 ) + W * ( del12 * bu3 - del32 * bu1 ), &
          W**3 * vd3 * ( v1b3 - v3b1 ) + W * ( del13 * bu3 - del33 * bu1 ), &
          Zero,            &
          Zero,            &
          W * ( vu1 * del31 - vu3 * del11 ),            &
          W * ( vu1 * del32 - vu3 * del12 ),            &
          W * ( vu1 * del33 - vu3 * del13 ),            &
          Zero ]

    dFdV_X1(10,:) &
      = [ Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero, &
          Zero ]

    ! --- ---

  END SUBROUTINE ComputePrimitiveFluxJacobian_X1


  SUBROUTINE ComputeConservedFluxJacobian_Numeric( dUdV, dFdV, dFdU )

    REAL(DP), INTENT(in)  :: dUdV(nCM,nCM)
    REAL(DP), INTENT(in)  :: dFdV(nCM,nCM)
    REAL(DP), INTENT(out) :: dFdU(nCM,nCM)

    INTEGER  :: INFO, LWORK
    INTEGER  :: IPIV(nCM)
    REAL(DP) :: TEMP(1)
    REAL(DP) :: dVdU(nCM,nCM)

    ! --- Compute inverse of dUdV ---

    CALL DGETRF( nCM, nCM, dVdU, nCM, IPIV, INFO )

    LWORK = -1

    CALL DGETRI( nCM, dVdU, nCM, IPIV, TEMP, LWORK, INFO )

    CALL MatrixMatrixMultiply( 'N', 'N', nCM, nCM, nCM, One, dVdU, nCM, dFdV, nCM, Zero, dFdU, nCM )

  END SUBROUTINE ComputeConservedFluxJacobian_Numeric


END MODULE MHD_CharacteristicDecompositionModule_Relativistic_IDEAL
