MODULE Euler_UtilitiesModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, SqrtTiny
  USE UnitsModule, ONLY: &
    AtomicMassUnit
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_h_1, iGF_h_2, iGF_h_3, &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_P, iAF_Cs, iAF_E, iAF_Ye
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive
  PUBLIC :: ComputeConserved
  PUBLIC :: ComputeFromConserved
  PUBLIC :: ComputeTimeStep
  PUBLIC :: Eigenvalues
  PUBLIC :: AlphaPlus
  PUBLIC :: AlphaMinus
  PUBLIC :: AlphaMiddle
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3
  PUBLIC :: StressTensor_Diagonal
  PUBLIC :: NumericalFlux_HLL
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
    S_1 = D * Gm_dd_11 * V_1
    S_2 = D * Gm_dd_22 * V_2
    S_3 = D * Gm_dd_33 * V_3
    G   = E + Half * D * ( Gm_dd_11 * V_1**2 &
                           + Gm_dd_22 * V_2**2 &
                           + Gm_dd_33 * V_3**2 )
    Ne  = De

  END SUBROUTINE ComputeConserved


  SUBROUTINE ComputeFromConserved( iX_B, iX_E, G, U, P, A )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:), &
      U(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(DP), INTENT(out) :: &
      P(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:), &
      A(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)

    INTEGER :: iX1, iX2, iX3

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      CALL ComputePrimitive &
             ( U(:,iX1,iX2,iX3,iCF_D ), U(:,iX1,iX2,iX3,iCF_S1), &
               U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
               U(:,iX1,iX2,iX3,iCF_E ), U(:,iX1,iX2,iX3,iCF_Ne), &
               P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_V1), &
               P(:,iX1,iX2,iX3,iPF_V2), P(:,iX1,iX2,iX3,iPF_V3), &
               P(:,iX1,iX2,iX3,iPF_E ), P(:,iX1,iX2,iX3,iPF_Ne), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      ! Temporary Hack
      A(:,iX1,iX2,iX3,iAF_E)  = P(:,iX1,iX2,iX3,iPF_E) / P(:,iX1,iX2,iX3,iPF_D)

      A(:,iX1,iX2,iX3,iAF_Ye) = AtomicMassUnit * U(:,iX1,iX2,iX3,iCF_Ne) &
                                / U(:,iX1,iX2,iX3,iCF_D)

      CALL ComputePressureFromPrimitive &
             ( P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_E), &
               P(:,iX1,iX2,iX3,iPF_Ne), A(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeSoundSpeedFromPrimitive &
             ( P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_E ), &
               P(:,iX1,iX2,iX3,iPF_Ne), A(:,iX1,iX2,iX3,iAF_Cs) )

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved


  SUBROUTINE ComputeTimeStep( iX_B, iX_E, G, U, CFL, TimeStep )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:), &
      U(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep

    INTEGER  :: iX1, iX2, iX3
    REAL(DP) :: dX(3), dt_X(nDOFX,3)
    REAL(DP) :: P(nDOFX,nPF)
    REAL(DP) :: A(nDOFX,nAF)

    TimeStep = HUGE( One )

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      dX(1) = MeshX(1) % Width(iX1)
      dX(2) = MeshX(2) % Width(iX2)
      dX(3) = MeshX(3) % Width(iX3)

      CALL ComputePrimitive &
             ( U(:,iX1,iX2,iX3,iCF_D ), U(:,iX1,iX2,iX3,iCF_S1), &
               U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
               U(:,iX1,iX2,iX3,iCF_E ), U(:,iX1,iX2,iX3,iCF_Ne), &
               P(:,iPF_D), P(:,iPF_V1), P(:,iPF_V2), P(:,iPF_V3), &
               P(:,iPF_E), P(:,iPF_Ne), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( P(:,iPF_D), P(:,iPF_E), P(:,iPF_Ne), A(:,iAF_P) )

      CALL ComputeSoundSpeedFromPrimitive &
             ( P(:,iPF_D), P(:,iPF_E), P(:,iPF_Ne), A(:,iAF_Cs) )

      dt_X(:,1) &
        = dX(1) * G(:,iX1,iX2,iX3,iGF_h_1) &
            / MAX( ABS( P(:,iPF_V1) ) + A(:,iAF_Cs), SqrtTiny )
      dt_X(:,2) &
        = dX(2) * G(:,iX1,iX2,iX3,iGF_h_2) &
            / MAX( ABS( P(:,iPF_V2) ) + A(:,iAF_Cs), SqrtTiny )
      dt_X(:,3) &
        = dX(3) * G(:,iX1,iX2,iX3,iGF_h_3) &
            / MAX( ABS( P(:,iPF_V3) ) + A(:,iAF_Cs), SqrtTiny )

      TimeStep = MIN( TimeStep, MINVAL( dt_X ) )

    END DO
    END DO
    END DO

    TimeStep = CFL * TimeStep

  END SUBROUTINE ComputeTimeStep


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


  PURE FUNCTION NumericalFlux_HLL( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd )

    REAL(DP),                 INTENT(in) :: aP, aM, aC, Gm_dd
    REAL(DP), DIMENSION(nCF), INTENT(in) :: uL, FluxL, uR, FluxR
    REAL(DP), DIMENSION(nCF)             :: NumericalFlux_HLL

    NumericalFlux_HLL &
      = ( aP * FluxL + aM * FluxR - aP * aM * ( uR - uL ) ) / ( aP + aM )

    RETURN
  END FUNCTION NumericalFlux_HLL


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


END MODULE Euler_UtilitiesModule
