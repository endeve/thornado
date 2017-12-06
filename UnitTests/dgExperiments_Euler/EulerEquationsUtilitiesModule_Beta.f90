MODULE EulerEquationsUtilitiesModule_Beta

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive
  PUBLIC :: ComputeConserved
  PUBLIC :: Eigenvalues
  PUBLIC :: AlphaPlus
  PUBLIC :: AlphaMinus
  PUBLIC :: AlphaMiddle
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3
  PUBLIC :: StressTensor_Diagonal
  PUBLIC :: NumericalFlux_X1_HLLC
  PUBLIC :: NumericalFlux_X2_HLLC
  PUBLIC :: NumericalFlux_X3_HLLC

CONTAINS


  SUBROUTINE ComputePrimitive &
               ( N, S_1, S_2, S_3, G, Ne, D, V_1, V_2, V_3, E, De, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:), INTENT(in)  :: N, S_1, S_2, S_3, G, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: D, V_1, V_2, V_3, E, De
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Three-Velocity: Index Up   ---
    ! --- Three-Momentum: Index Down ---

    D   = N
    V_1 = S_1 / ( Gm_dd_11 * N )
    V_2 = S_2 / ( Gm_dd_22 * N )
    V_3 = S_3 / ( Gm_dd_33 * N )
    E   = G - Half * ( S_1**2 / Gm_dd_11 &
                       + S_2**2 / Gm_dd_22 &
                       + S_3**2 / Gm_dd_33 ) / N
    De  = Ne

  END SUBROUTINE ComputePrimitive


  SUBROUTINE ComputeConserved &
               ( D, V_1, V_2, V_3, E, De, N, S_1, S_2, S_3, G, Ne, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, V_1, V_2, V_3, E, De
    REAL(DP), DIMENSION(:), INTENT(out) :: N, S_1, S_2, S_3, G, Ne
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Three-Velocity: Index Up   ---
    ! --- Three-Momentum: Index Down ---

    N   = D
    S_1 = D * V_1
    S_2 = D * V_2
    S_3 = D * V_3
    G   = E + Half * D * ( Gm_dd_11 * V_1**2 &
                           + Gm_dd_22 * V_2**2 &
                           + Gm_dd_33 * V_3**2 )
    Ne  = De

  END SUBROUTINE ComputeConserved


  PURE FUNCTION Eigenvalues( V, Cs )

    REAL(DP), INTENT(in)     :: V, Cs
    REAL(DP), DIMENSION(nCF) :: Eigenvalues

    Eigenvalues = [ V - Cs, V, V, V, V, V + Cs ]

    RETURN
  END FUNCTION Eigenvalues


  PURE ELEMENTAL REAL(DP) FUNCTION AlphaPlus &
    ( V_L, Cs_L, V_R, Cs_R, Gm_dd_ii )

    REAL(DP), INTENT(in) :: V_L, Cs_L, V_R, Cs_R, Gm_dd_ii

    REAL(DP) :: SqrtGm_uu_ii

    SqrtGm_uu_ii = SQRT( One / Gm_dd_ii )

    AlphaPlus = MAX( Zero, V_L + SqrtGm_uu_ii * Cs_L, &
                           V_R + SqrtGm_uu_ii * Cs_R )

    RETURN
  END FUNCTION AlphaPlus


  PURE ELEMENTAL REAL(DP) FUNCTION AlphaMinus &
    ( V_L, Cs_L, V_R, Cs_R, Gm_dd_ii )

    REAL(DP), INTENT(in) :: V_L, Cs_L, V_R, Cs_R, Gm_dd_ii

    REAL(DP) :: SqrtGm_uu_ii

    SqrtGm_uu_ii = SQRT( One / Gm_dd_ii )

    AlphaMinus = MAX( Zero, SqrtGm_uu_ii * Cs_L - V_L, &
                            SqrtGm_uu_ii * Cs_R - V_R )

    RETURN
  END FUNCTION AlphaMinus


  PURE ELEMENTAL REAL(DP) FUNCTION AlphaMiddle &
    ( D_L, FD_L, S_L, FS_L, D_R, FD_R, S_R, FS_R, aP, aM, Gm_dd_ii )

    REAL(DP), INTENT(in) :: D_L, FD_L, S_L, FS_L
    REAL(DP), INTENT(in) :: D_R, FD_R, S_R, FS_R
    REAL(DP), INTENT(in) :: aP, aM, Gm_dd_ii

    ! --- Middle Wavespeed as Suggested by Batten et al. (1997) ---
    ! --- (SIAM J. Sci. Comput., Vol. 18, No. 6, pp. 1553-1570) ---

    AlphaMiddle & ! --- Index Up
      = ( aP * S_R + aM * S_L - ( FS_R - FS_L ) ) &
        / ( aP * D_R + aM * D_L - ( FD_R - FD_L ) ) / Gm_dd_ii

    RETURN
  END FUNCTION AlphaMiddle


  PURE FUNCTION Flux_X1 &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Flux_X1(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: VSq

    VSq = Gm_dd_11 * V_1**2 + Gm_dd_22 * V_2**2 + Gm_dd_33 * V_3**2

    Flux_X1(iCF_D ) = D * V_1
    Flux_X1(iCF_S1) = D * Gm_dd_11 * V_1 * V_1 + P
    Flux_X1(iCF_S2) = D * Gm_dd_22 * V_2 * V_1
    Flux_X1(iCF_S3) = D * Gm_dd_33 * V_3 * V_1
    Flux_X1(iCF_E ) = ( E + Half * D * VSq + P ) * V_1
    Flux_X1(iCF_Ne) = Ne * V_1

    RETURN
  END FUNCTION Flux_X1


  PURE FUNCTION Flux_X2 &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Flux_X2(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: VSq

    VSq = Gm_dd_11 * V_1**2 + Gm_dd_22 * V_2**2 + Gm_dd_33 * V_3**2

    Flux_X2(iCF_D ) = D * V_2
    Flux_X2(iCF_S1) = D * Gm_dd_11 * V_1 * V_2
    Flux_X2(iCF_S2) = D * Gm_dd_22 * V_2 * V_2 + P
    Flux_X2(iCF_S3) = D * Gm_dd_33 * V_3 * V_2
    Flux_X2(iCF_E ) = ( E + Half * D * VSq + P ) * V_2
    Flux_X2(iCF_Ne) = Ne * V_2

    RETURN
  END FUNCTION Flux_X2


  PURE FUNCTION Flux_X3 &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Flux_X3(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    Flux_X3 = 0.0_DP

    RETURN
  END FUNCTION Flux_X3


  PURE FUNCTION StressTensor_Diagonal &
    ( D, V_1, V_2, V_3, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: StressTensor_Diagonal(1:3)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    StressTensor_Diagonal(1) = D * Gm_dd_11 * V_1**2 + P
    StressTensor_Diagonal(2) = D * Gm_dd_22 * V_2**2 + P
    StressTensor_Diagonal(3) = D * Gm_dd_33 * V_3**2 + P

    RETURN
  END FUNCTION StressTensor_Diagonal


  PURE FUNCTION NumericalFlux_X1_HLLC &
    ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd_11 )

    REAL(DP),                 INTENT(in) :: aP, aM, aC, Gm_dd_11
    REAL(DP), DIMENSION(nCF), INTENT(in) :: uL, FluxL, uR, FluxR
    REAL(DP), DIMENSION(nCF)             :: NumericalFlux_X1_HLLC

    REAL(DP) :: D, V1, V2, V3, P, E, Ne
    REAL(DP), DIMENSION(nCF) :: TMP

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X1_HLLC = FluxL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X1_HLLC = FluxR

    ELSE

      IF( aC .GE. Zero )THEN

        TMP = FluxL + aM * uL

        D  = TMP(iCF_D) / ( aC + aM )
        V1 = aC                       ! --- Index Up
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S1) - Gm_dd_11 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC + aM )
        Ne = TMP(iCF_Ne) / ( aC + aM )

      ELSE

        TMP = FluxR - aP * uR

        D  = TMP(iCF_D) / ( aC - aP )
        V1 = aC                       ! --- Index Up
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S1) - Gm_dd_11 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC - aP )
        Ne = TMP(iCF_Ne) / ( aC - aP )

      END IF

      NumericalFlux_X1_HLLC(iCF_D) &
        = D * V1
      NumericalFlux_X1_HLLC(iCF_S1) &
        = D * Gm_dd_11 * V1 * V1 + P
      NumericalFlux_X1_HLLC(iCF_S2) &
        = D * V2 * V1
      NumericalFlux_X1_HLLC(iCF_S3) &
        = D * V3 * V1
      NumericalFlux_X1_HLLC(iCF_E) &
        = ( E + P ) * V1
      NumericalFlux_X1_HLLC(iCF_Ne) &
        = Ne * V1

    END IF

    RETURN
  END FUNCTION NumericalFlux_X1_HLLC


  PURE FUNCTION NumericalFlux_X2_HLLC &
    ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd_22 )

    REAL(DP),                 INTENT(in) :: aP, aM, aC, Gm_dd_22
    REAL(DP), DIMENSION(nCF), INTENT(in) :: uL, FluxL, uR, FluxR
    REAL(DP), DIMENSION(nCF)             :: NumericalFlux_X2_HLLC

    REAL(DP) :: D, V1, V2, V3, P, E, Ne
    REAL(DP), DIMENSION(nCF) :: TMP

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X2_HLLC = FluxL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X2_HLLC = FluxR

    ELSE

      IF( aC .GE. Zero )THEN

        TMP = FluxL + aM * uL

        D  = TMP(iCF_D) / ( aC + aM )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = aC                       ! --- Index Up
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S2) - Gm_dd_22 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC + aM )
        Ne = TMP(iCF_Ne) / ( aC + aM )

      ELSE

        TMP = FluxR - aP * uR

        D  = TMP(iCF_D) / ( aC - aP )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = aC                       ! --- Index Up
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S2) - Gm_dd_22 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC - aP )
        Ne = TMP(iCF_Ne) / ( aC - aP )

      END IF

      NumericalFlux_X2_HLLC(iCF_D) &
        = D * V2
      NumericalFlux_X2_HLLC(iCF_S2) &
        = D * V1 * V2
      NumericalFlux_X2_HLLC(iCF_S2) &
        = D * Gm_dd_22 * V2 * V2 + P
      NumericalFlux_X2_HLLC(iCF_S3) &
        = D * V3 * V2
      NumericalFlux_X2_HLLC(iCF_E) &
        = ( E + P ) * V2
      NumericalFlux_X2_HLLC(iCF_Ne) &
        = Ne * V2

    END IF

    RETURN
  END FUNCTION NumericalFlux_X2_HLLC


  PURE FUNCTION NumericalFlux_X3_HLLC &
    ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd_33 )

    REAL(DP),                 INTENT(in) :: aP, aM, aC, Gm_dd_33
    REAL(DP), DIMENSION(nCF), INTENT(in) :: uL, FluxL, uR, FluxR
    REAL(DP), DIMENSION(nCF)             :: NumericalFlux_X3_HLLC

    REAL(DP) :: D, V1, V2, V3, P, E, Ne
    REAL(DP), DIMENSION(nCF) :: TMP

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X3_HLLC = FluxL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X3_HLLC = FluxR

    ELSE

      IF( aC .GE. Zero )THEN

        TMP = FluxL + aM * uL

        D  = TMP(iCF_D) / ( aC + aM )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = aC                       ! --- Index Up
        P  = TMP(iCF_S3) - Gm_dd_33 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC + aM )
        Ne = TMP(iCF_Ne) / ( aC + aM )

      ELSE

        TMP = FluxR - aP * uR

        D  = TMP(iCF_D) / ( aC - aP )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = aC                       ! --- Index Up
        P  = TMP(iCF_S3) - Gm_dd_33 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC - aP )
        Ne = TMP(iCF_Ne) / ( aC - aP )

      END IF

      NumericalFlux_X3_HLLC(iCF_D) &
        = D * V3
      NumericalFlux_X3_HLLC(iCF_S2) &
        = D * V1 * V3
      NumericalFlux_X3_HLLC(iCF_S2) &
        = D * V2 * V3
      NumericalFlux_X3_HLLC(iCF_S3) &
        = D * Gm_dd_33 * V3 * V3 + P
      NumericalFlux_X3_HLLC(iCF_E) &
        = ( E + P ) * V3
      NumericalFlux_X3_HLLC(iCF_Ne) &
        = Ne * V3

    END IF

    RETURN
  END FUNCTION NumericalFlux_X3_HLLC


END MODULE EulerEquationsUtilitiesModule_Beta
