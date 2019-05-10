MODULE Euler_UtilitiesModule_NonRelativistic

  USE KindModule, ONLY: &
    DP, Zero, Half, One, SqrtTiny
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
    nAF, iAF_P, iAF_Cs
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_ComputePrimitive_NonRelativistic
  PUBLIC :: Euler_ComputeConserved_NonRelativistic
  PUBLIC :: Euler_ComputeFromConserved_NonRelativistic
  PUBLIC :: Euler_ComputeTimeStep_NonRelativistic
  PUBLIC :: Euler_Eigenvalues_NonRelativistic
  PUBLIC :: Euler_AlphaMiddle_NonRelativistic
  PUBLIC :: Euler_Flux_X1_NonRelativistic
  PUBLIC :: Euler_Flux_X2_NonRelativistic
  PUBLIC :: Euler_Flux_X3_NonRelativistic
  PUBLIC :: Euler_StressTensor_Diagonal_NonRelativistic
  PUBLIC :: Euler_NumericalFlux_HLL_NonRelativistic
  PUBLIC :: Euler_NumericalFlux_X1_HLLC_NonRelativistic
  PUBLIC :: Euler_NumericalFlux_X2_HLLC_NonRelativistic
  PUBLIC :: Euler_NumericalFlux_X3_HLLC_NonRelativistic

CONTAINS


  SUBROUTINE Euler_ComputePrimitive_NonRelativistic &
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

  END SUBROUTINE Euler_ComputePrimitive_NonRelativistic


  SUBROUTINE Euler_ComputeConserved_NonRelativistic &
               ( D, V_1, V_2, V_3, E, De, N, S_1, S_2, S_3, G, Ne, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, P )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, V_1, V_2, V_3, E, De
    REAL(DP), DIMENSION(:), INTENT(out) :: N, S_1, S_2, S_3, G, Ne
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Only used for relativistic code ---
    REAL(DP), INTENT(in) :: P(:)

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

  END SUBROUTINE Euler_ComputeConserved_NonRelativistic


  SUBROUTINE Euler_ComputeFromConserved_NonRelativistic &
               ( iX_B, iX_E, G, U, P, A )

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

      CALL Euler_ComputePrimitive_NonRelativistic &
             ( U(:,iX1,iX2,iX3,iCF_D ), U(:,iX1,iX2,iX3,iCF_S1), &
               U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
               U(:,iX1,iX2,iX3,iCF_E ), U(:,iX1,iX2,iX3,iCF_Ne), &
               P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_V1), &
               P(:,iX1,iX2,iX3,iPF_V2), P(:,iX1,iX2,iX3,iPF_V3), &
               P(:,iX1,iX2,iX3,iPF_E ), P(:,iX1,iX2,iX3,iPF_Ne), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_E), &
               P(:,iX1,iX2,iX3,iPF_Ne), A(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeSoundSpeedFromPrimitive &
             ( P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_E ), &
               P(:,iX1,iX2,iX3,iPF_Ne), A(:,iX1,iX2,iX3,iAF_Cs) )

    END DO
    END DO
    END DO

  END SUBROUTINE Euler_ComputeFromConserved_NonRelativistic


  SUBROUTINE Euler_ComputeTimeStep_NonRelativistic &
               ( iX_B, iX_E, G, U, CFL, TimeStep )

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

      CALL Euler_ComputePrimitive_NonRelativistic &
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

  END SUBROUTINE Euler_ComputeTimeStep_NonRelativistic


  PURE FUNCTION Euler_Eigenvalues_NonRelativistic &
    ( V, Cs, V1, V2, V3, Gm, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in)     :: V, Cs
    REAL(DP), DIMENSION(nCF) :: Euler_Eigenvalues_NonRelativistic

    ! --- Only used for relativistic code ---
    REAL(DP), INTENT(in) :: V1, V2, V3, Gm, Gm11, Gm22, Gm33, Lapse, Shift

    Euler_Eigenvalues_NonRelativistic = [ V - Cs, V, V, V, V, V + Cs ]

    RETURN
  END FUNCTION Euler_Eigenvalues_NonRelativistic


  REAL(DP) FUNCTION Euler_AlphaMiddle_NonRelativistic &
    ( D_L, S_L, E_L, FD_L, FS_L, FE_L, D_R, S_R, E_R, FD_R, FS_R, FE_R, &
      Gm_dd_ii, Lapse, Shift, aP, aM )

    REAL(DP), INTENT(in) :: D_L, S_L, E_L, FD_L, FS_L, FE_L, &
                            D_R, S_R, E_R, FD_R, FS_R, FE_R, &
                            Gm_dd_ii, aP, aM

    ! --- Only used for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

    ! --- Middle Wavespeed as Suggested by Batten et al. (1997) ---
    ! --- (SIAM J. Sci. Comput., Vol. 18, No. 6, pp. 1553-1570) ---

    Euler_AlphaMiddle_NonRelativistic & ! --- Index Up
      = ( aP * S_R + aM * S_L - ( FS_R - FS_L ) ) &
        / ( aP * D_R + aM * D_L - ( FD_R - FD_L ) ) / Gm_dd_ii

    RETURN
  END FUNCTION Euler_AlphaMiddle_NonRelativistic


  PURE FUNCTION Euler_Flux_X1_NonRelativistic &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33, Lapse, Shift )

    REAL(DP)             :: Euler_Flux_X1_NonRelativistic(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Only used for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

    REAL(DP) :: VSq

    VSq = Gm_dd_11 * V_1**2 + Gm_dd_22 * V_2**2 + Gm_dd_33 * V_3**2

    Euler_Flux_X1_NonRelativistic(iCF_D ) = D * V_1
    Euler_Flux_X1_NonRelativistic(iCF_S1) = D * Gm_dd_11 * V_1 * V_1 + P
    Euler_Flux_X1_NonRelativistic(iCF_S2) = D * Gm_dd_22 * V_2 * V_1
    Euler_Flux_X1_NonRelativistic(iCF_S3) = D * Gm_dd_33 * V_3 * V_1
    Euler_Flux_X1_NonRelativistic(iCF_E ) = ( E + Half * D * VSq + P ) * V_1
    Euler_Flux_X1_NonRelativistic(iCF_Ne) = Ne * V_1

    RETURN
  END FUNCTION Euler_Flux_X1_NonRelativistic


  PURE FUNCTION Euler_Flux_X2_NonRelativistic &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33, Lapse, Shift )

    REAL(DP)             :: Euler_Flux_X2_NonRelativistic(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Only used for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

    REAL(DP) :: VSq

    VSq = Gm_dd_11 * V_1**2 + Gm_dd_22 * V_2**2 + Gm_dd_33 * V_3**2

    Euler_Flux_X2_NonRelativistic(iCF_D ) = D * V_2
    Euler_Flux_X2_NonRelativistic(iCF_S1) = D * Gm_dd_11 * V_1 * V_2
    Euler_Flux_X2_NonRelativistic(iCF_S2) = D * Gm_dd_22 * V_2 * V_2 + P
    Euler_Flux_X2_NonRelativistic(iCF_S3) = D * Gm_dd_33 * V_3 * V_2
    Euler_Flux_X2_NonRelativistic(iCF_E ) = ( E + Half * D * VSq + P ) * V_2
    Euler_Flux_X2_NonRelativistic(iCF_Ne) = Ne * V_2

    RETURN
  END FUNCTION Euler_Flux_X2_NonRelativistic


  PURE FUNCTION Euler_Flux_X3_NonRelativistic &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33, Lapse, Shift )

    REAL(DP)             :: Euler_Flux_X3_NonRelativistic(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Only used for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

    Euler_Flux_X3_NonRelativistic = 0.0_DP

    RETURN
  END FUNCTION Euler_Flux_X3_NonRelativistic


  PURE FUNCTION Euler_StressTensor_Diagonal_NonRelativistic &
    ( S1, S2, S3, V1, V2, V3, P )

    REAL(DP), INTENT(in) :: S1, S2, S3, V1, V2, V3, P

    REAL(DP) :: Euler_StressTensor_Diagonal_NonRelativistic(3)

    Euler_StressTensor_Diagonal_NonRelativistic(1) = S1 * V1 + P
    Euler_StressTensor_Diagonal_NonRelativistic(2) = S2 * V2 + P
    Euler_StressTensor_Diagonal_NonRelativistic(3) = S3 * V3 + P

    RETURN
  END FUNCTION Euler_StressTensor_Diagonal_NonRelativistic


  PURE FUNCTION Euler_NumericalFlux_HLL_NonRelativistic &
    ( uL, uR, fL, fR, Gm_dd, vL, vR, pL, pR, Lapse, Shift, aP, aM, aC )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF)
    REAL(DP), INTENT(in) :: Gm_dd, aP, aM, aC

    ! --- Only used in relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_HLL_NonRelativistic(nCF)

    Euler_NumericalFlux_HLL_NonRelativistic &
      = ( aP * fL + aM * fR - aP * aM * ( uR - uL ) ) / ( aP + aM )

    RETURN
  END FUNCTION Euler_NumericalFlux_HLL_NonRelativistic


  PURE FUNCTION Euler_NumericalFlux_X1_HLLC_NonRelativistic &
    ( uL, uR, fL, fR, Gm_dd_11, vL, vR, pL, pR, Lapse, Shift_1, aP, aM, aC )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            Gm_dd_11, aP, aM, aC

    ! --- Only used for relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift_1

    REAL(DP) :: Euler_NumericalFlux_X1_HLLC_NonRelativistic(nCF)

    REAL(DP) :: D, V1, V2, V3, P, E, Ne, TMP(nCF)

    IF( aM .EQ. Zero )THEN

      Euler_NumericalFlux_X1_HLLC_NonRelativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      Euler_NumericalFlux_X1_HLLC_NonRelativistic = fR

    ELSE

      IF( aC .GE. Zero )THEN

        TMP = fL + aM * uL

        D  = TMP(iCF_D) / ( aC + aM )
        V1 = aC                       ! --- Index Up
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S1) - Gm_dd_11 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC + aM )
        Ne = TMP(iCF_Ne) / ( aC + aM )

      ELSE

        TMP = fR - aP * uR

        D  = TMP(iCF_D) / ( aC - aP )
        V1 = aC                       ! --- Index Up
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S1) - Gm_dd_11 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC - aP )
        Ne = TMP(iCF_Ne) / ( aC - aP )

      END IF

      Euler_NumericalFlux_X1_HLLC_NonRelativistic(iCF_D) &
        = D * V1
      Euler_NumericalFlux_X1_HLLC_NonRelativistic(iCF_S1) &
        = D * Gm_dd_11 * V1 * V1 + P
      Euler_NumericalFlux_X1_HLLC_NonRelativistic(iCF_S2) &
        = D * V2 * V1
      Euler_NumericalFlux_X1_HLLC_NonRelativistic(iCF_S3) &
        = D * V3 * V1
      Euler_NumericalFlux_X1_HLLC_NonRelativistic(iCF_E) &
        = ( E + P ) * V1
      Euler_NumericalFlux_X1_HLLC_NonRelativistic(iCF_Ne) &
        = Ne * V1

    END IF

    RETURN
  END FUNCTION Euler_NumericalFlux_X1_HLLC_NonRelativistic


  PURE FUNCTION Euler_NumericalFlux_X2_HLLC_NonRelativistic &
    ( uL, uR, fL, fR, Gm_dd_22, vL, vR, pL, pR, Lapse, Shift_2, aP, aM, aC )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            Gm_dd_22, aP, aM, aC

    ! --- Only used for relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift_2

    REAL(DP) :: Euler_NumericalFlux_X2_HLLC_NonRelativistic(nCF)
    
    REAL(DP) :: D, V1, V2, V3, P, E, Ne, TMP(nCF)

    IF( aM .EQ. Zero )THEN

      Euler_NumericalFlux_X2_HLLC_NonRelativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      Euler_NumericalFlux_X2_HLLC_NonRelativistic = fR

    ELSE

      IF( aC .GE. Zero )THEN

        TMP = fL + aM * uL

        D  = TMP(iCF_D) / ( aC + aM )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = aC                       ! --- Index Up
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S2) - Gm_dd_22 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC + aM )
        Ne = TMP(iCF_Ne) / ( aC + aM )

      ELSE

        TMP = fR - aP * uR

        D  = TMP(iCF_D) / ( aC - aP )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = aC                       ! --- Index Up
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S2) - Gm_dd_22 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC - aP )
        Ne = TMP(iCF_Ne) / ( aC - aP )

      END IF

      Euler_NumericalFlux_X2_HLLC_NonRelativistic(iCF_D) &
        = D * V2
      Euler_NumericalFlux_X2_HLLC_NonRelativistic(iCF_S2) &
        = D * V1 * V2
      Euler_NumericalFlux_X2_HLLC_NonRelativistic(iCF_S2) &
        = D * Gm_dd_22 * V2 * V2 + P
      Euler_NumericalFlux_X2_HLLC_NonRelativistic(iCF_S3) &
        = D * V3 * V2
      Euler_NumericalFlux_X2_HLLC_NonRelativistic(iCF_E) &
        = ( E + P ) * V2
      Euler_NumericalFlux_X2_HLLC_NonRelativistic(iCF_Ne) &
        = Ne * V2

    END IF

    RETURN
  END FUNCTION Euler_NumericalFlux_X2_HLLC_NonRelativistic


  PURE FUNCTION Euler_NumericalFlux_X3_HLLC_NonRelativistic &
    ( uL, uR, fL, fR, Gm_dd_33, vL, vR, pL, pR, Lapse, Shift_3, aP, aM, aC )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            Gm_dd_33, aP, aM, aC

    ! --- Only used for relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift_3

    REAL(DP) :: Euler_NumericalFlux_X3_HLLC_NonRelativistic(nCF)

    REAL(DP) :: D, V1, V2, V3, P, E, Ne, TMP(nCF)

    IF( aM .EQ. Zero )THEN

      Euler_NumericalFlux_X3_HLLC_NonRelativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      Euler_NumericalFlux_X3_HLLC_NonRelativistic = fR

    ELSE

      IF( aC .GE. Zero )THEN

        TMP = fL + aM * uL

        D  = TMP(iCF_D) / ( aC + aM )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = aC                       ! --- Index Up
        P  = TMP(iCF_S3) - Gm_dd_33 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC + aM )
        Ne = TMP(iCF_Ne) / ( aC + aM )

      ELSE

        TMP = fR - aP * uR

        D  = TMP(iCF_D) / ( aC - aP )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = aC                       ! --- Index Up
        P  = TMP(iCF_S3) - Gm_dd_33 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC - aP )
        Ne = TMP(iCF_Ne) / ( aC - aP )

      END IF

      Euler_NumericalFlux_X3_HLLC_NonRelativistic(iCF_D) &
        = D * V3
      Euler_NumericalFlux_X3_HLLC_NonRelativistic(iCF_S2) &
        = D * V1 * V3
      Euler_NumericalFlux_X3_HLLC_NonRelativistic(iCF_S2) &
        = D * V2 * V3
      Euler_NumericalFlux_X3_HLLC_NonRelativistic(iCF_S3) &
        = D * Gm_dd_33 * V3 * V3 + P
      Euler_NumericalFlux_X3_HLLC_NonRelativistic(iCF_E) &
        = ( E + P ) * V3
      Euler_NumericalFlux_X3_HLLC_NonRelativistic(iCF_Ne) &
        = Ne * V3

    END IF

    RETURN
  END FUNCTION Euler_NumericalFlux_X3_HLLC_NonRelativistic


END MODULE Euler_UtilitiesModule_NonRelativistic
