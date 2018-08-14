MODULE CharacteristicDecompositionModule_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, SqrtTiny
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputePrimitive_GR, &
    Eigenvalues_GR
  USE EquationOfStateModule, ONLY: &
    ComputeSoundSpeedFromPrimitive_GR

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeCharacteristicDecomposition_GR

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_GR( iDim, G, U, R, LT )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,0:nCF-1)
    REAL(DP), INTENT(out) :: LT(nCF,0:nCF-1)

    INTEGER :: i
    REAL(DP), DIMENSION(1)   :: D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, E, Ne, P
    REAL(DP), DIMENSION(1)   :: Gmdd11, Gmdd22, Gmdd33, SqrtGm
    REAL(DP), DIMENSION(1)   :: W, h
    REAL(DP), DIMENSION(nCF) :: EigVals
    REAL(DP), DIMENSION(1)   :: ScriptK, GammaXX, DELTA, xi, Cs
    REAL(DP), DIMENSION(1)   :: LAMBDA_m, LAMBDA_p
    REAL(DP), DIMENSION(1)   :: ScriptA_m, ScriptA_p
    REAL(DP), DIMENSION(1)   :: ScriptV_m, ScriptV_p
    REAL(DP), DIMENSION(1)   :: ScriptN_m, ScriptN_p
    REAL(DP), DIMENSION(1)   :: ScriptC_m, ScriptC_p
    REAL(DP), DIMENSION(1)   :: Gamma = 4.0d0 / 3.0e0
    REAL(DP)                 :: L(nCF,0:nCF-1)

    CALL ComputePrimitive_GR &
           ( [U(iCF_D )], [U(iCF_S1)], [U(iCF_S2)], &
             [U(iCF_S3)], [U(iCF_E )], [U(iCF_Ne)], &
             D, Vu1, Vu2, Vu3, E, Ne, P, &
             [G(iGF_Gm_dd_11)], &
             [G(iGF_Gm_dd_22)], &
             [G(iGF_Gm_dd_33)] )

    CALL ComputeSoundSpeedFromPrimitive_GR( D, E, Ne, Cs )
     
    Gmdd11 = G(iGF_Gm_dd_11)
    Gmdd22 = G(iGF_Gm_dd_22)
    Gmdd33 = G(iGF_Gm_dd_33)
    SqrtGm = G(iGF_SqrtGm)

    Vd1 = Gmdd11 * Vu1
    Vd2 = Gmdd22 * Vu2
    Vd3 = Gmdd33 * Vu3

    ! --- Compute auxiliary quantities for right eigenvectors a la Rezzolla ---
    W = One / SQRT( One - ( Vu1*Vd1 + Vu2*Vd2 + Vu3*Vd3 ) )
    h = One + ( E + P ) / D
    ScriptK = ( Gamma - One ) / ( Gamma - One - Cs**2 ) ! --- Assumes ideal EOS ---

    SELECT CASE( iDim )

      CASE( 1 )

        EigVals = Eigenvalues_GR &
                    ( Vu1(1), Vu2(1), Vu3(1), Vu1(1), Cs(1), &
                      Gmdd11(1), Gmdd22(1), Gmdd33(1), Gmdd11(1), &
                      G(iGF_Alpha), G(iGF_Beta_1) )

        LAMBDA_m = ( EigVals(1) + G(iGF_Beta_1) ) / G(iGF_Alpha)
        LAMBDA_p = ( EigVals(3) + G(iGF_Beta_1) ) / G(iGF_Alpha)

        ScriptV_m = ( Vu1 - LAMBDA_m ) / ( One / Gmdd11 - Vu1 * LAMBDA_m )
        ScriptV_p = ( Vu1 - LAMBDA_p ) / ( One / Gmdd11 - Vu1 * LAMBDA_p )

        ScriptA_m = ( One / Gmdd11 - Vu1 * Vu1 ) &
                      / ( One / Gmdd11 - Vu1 * LAMBDA_m )
        ScriptA_p = ( One / Gmdd11 - Vu1 * Vu1 ) &
                      / ( One / Gmdd11 - Vu1 * LAMBDA_p )

        R(:,0) = [ One, &
                   h * W * ( Vd1 - ScriptV_m ), &
                   h * W * Vd2, &
                   h * W * Vd3, &
                   h * W * ScriptA_m - One, &
                   Zero ]
        R(:,1) = [ ScriptK / ( MAX( h * W, SqrtTiny ) ), &
                   Vd1, &
                   Vd2, &
                   Vd3, &
                   One - ScriptK / ( MAX( h * W, SqrtTiny ) ), &
                   Zero ]
        R(:,2) = [ W * Vd2, &
                   h *            Two * W**2 * Vd1 * Vd2,   &
                   h * ( Gmdd22 + Two * W**2 * Vd2 * Vd2 ), &
                   h *            Two * W**2 * Vd3 * Vd2,   &
                   W * Vd2 * ( Two * h * W - One ), &
                   Zero ]
        R(:,3) = [ W * Vd3, &
                   h *            Two * W**2 * Vd1 * Vd3,   &
                   h *            Two * W**2 * Vd2 * Vd3,   &
                   h * ( Gmdd33 + Two * W**2 * Vd3 * Vd3 ), &
                   W * Vd3 * ( Two * h * W - One ), &
                   Zero ]
        R(:,4) = [ One, &
                   h * W * ( Vd1 - ScriptV_p ), &
                   h * W * Vd2, &
                   h * W * Vd3, &
                   h * W * ScriptA_p - One, &
                   Zero ]
        R(:,5) = [ Zero, Zero, Zero, Zero, Zero, One ]


        ! --- Compute auxiliary quantities for left eigenvectors a la Rezzolla ---
        GammaXX = Gmdd22 * Gmdd33
        xi = GammaXX - SqrtGm * Vu1**2
        ScriptC_m = Vd1 - ScriptV_m
        ScriptC_p = Vd1 - ScriptV_p
        DELTA = h**3 * W * ( ScriptK - One ) * ( ScriptC_p - ScriptC_m ) * xi
        ScriptN_m = ( One - ScriptK ) &
                      * ( -SqrtGm * Vu1 + ScriptV_p * ( W**2 * xi - GammaXX ) ) &
                      - ScriptK * W**2 * ScriptV_p * xi
        ScriptN_p = ( One - ScriptK ) &
                      * ( -SqrtGm * Vu1 + ScriptV_m * ( W**2 * xi - GammaXX ) ) &
                      - ScriptK * W**2 * ScriptV_m * xi

        L(:,0) = h(1)**2 / DELTA(1) &
                      * [ h * W * ScriptV_p * xi + ScriptN_m, &
                          GammaXX * ( One - ScriptK * ScriptA_p ) &
                            + ( Two * ScriptK - One ) * ScriptN_p &
                            * ( W**2 * Vu1 * xi - GammaXX * Vu1 ), &
                          ( Two * ScriptK - One ) * ScriptN_p * W**2 * Vu2 * xi, &
                          ( Two * ScriptK - One ) * ScriptN_p * W**2 * Vu3 * xi, &
                          ScriptN_m, &
                          Zero ]
        L(:,1) = W(1) / ( ScriptK(1) - One ) &
                      * [ h - W, W * Vu1, W * Vu2, W * Vu3, -W, Zero ]
        L(:,2) = One / ( h(1) * xi(1) ) &
                      * [ -Gmdd33 * Vd2, &
                           Vu1 * Gmdd33 * Vd2, &
                           Gmdd33 * ( One - Vd1 * Vu1 ), &
                           Zero, &
                           -Gmdd33 * Vd2, &
                           Zero ]
        L(:,3) = One / ( h(1) * xi(1) ) &
                      * [ -Gmdd22 * Vd3, &
                           Vu1 * Gmdd22 * Vd3, &
                           Zero, &
                           Gmdd22 * ( One - Vd1 * Vu1 ), &
                           -Gmdd22 * Vd3, &
                           Zero ]
        L(:,4) = -h(1)**2 / DELTA(1) &
                      * [ h * W * ScriptV_m * xi + ScriptN_p, &
                          GammaXX * ( One - ScriptK * ScriptA_m ) &
                            + ( Two * ScriptK - One ) * ScriptN_m &
                            * ( W**2 * Vu1 * xi - GammaXX * Vu1 ), &
                          ( Two * ScriptK - One ) * ScriptN_m * W**2 * Vu2 * xi, &
                          ( Two * ScriptK - One ) * ScriptN_m * W**2 * Vu3 * xi, &
                          ScriptN_p, &
                          Zero ]
        L(:,5) = [ Zero, Zero, Zero, Zero, Zero, One ]

        LT = TRANSPOSE( L )

      CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A19,I2.2)') &
        '', 'Invalid dimension: ', iDim
      STOP

    END SELECT

  END SUBROUTINE ComputeCharacteristicDecomposition_GR


END MODULE CharacteristicDecompositionModule_GR
