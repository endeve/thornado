MODULE Euler_UtilitiesModule_NonRelativistic

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, Half, One
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_h_1, iGF_h_2, iGF_h_3, &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_Alpha, iGF_Beta_1, iGF_Beta_2, iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, iAF_Gm, iAF_Cs
  USE EquationOfStateModule, ONLY: &
    ComputeSoundSpeedFromPrimitive, &
    ComputeAuxiliary_Fluid
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_ComputeTimeStep

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_Euler_NonRelativistic
  PUBLIC :: ComputeConserved_Euler_NonRelativistic
  PUBLIC :: ComputeFromConserved_Euler_NonRelativistic
  PUBLIC :: ComputeTimeStep_Euler_NonRelativistic
  PUBLIC :: Eigenvalues_Euler_NonRelativistic
  PUBLIC :: AlphaMiddle_Euler_NonRelativistic
  PUBLIC :: Flux_X1_Euler_NonRelativistic
  PUBLIC :: Flux_X2_Euler_NonRelativistic
  PUBLIC :: Flux_X3_Euler_NonRelativistic
  PUBLIC :: StressTensor_Diagonal_Euler_NonRelativistic
  PUBLIC :: NumericalFlux_HLL_Euler_NonRelativistic
  PUBLIC :: NumericalFlux_X1_HLLC_Euler_NonRelativistic
  PUBLIC :: NumericalFlux_X2_HLLC_Euler_NonRelativistic
  PUBLIC :: NumericalFlux_X3_HLLC_Euler_NonRelativistic

  INTERFACE ComputePrimitive_Euler_NonRelativistic
    MODULE PROCEDURE ComputePrimitive_Scalar
    MODULE PROCEDURE ComputePrimitive_Vector
  END INTERFACE ComputePrimitive_Euler_NonRelativistic

CONTAINS


  SUBROUTINE ComputePrimitive_Scalar &
    ( N, S_1, S_2, S_3, G, Ne, D, V_1, V_2, V_3, E, De, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in ) :: N, S_1, S_2, S_3, G, Ne
    REAL(DP), INTENT(out) :: D, V_1, V_2, V_3, E, De
    REAL(DP), INTENT(in ) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

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

  END SUBROUTINE ComputePrimitive_Scalar


  SUBROUTINE ComputePrimitive_Vector &
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

  END SUBROUTINE ComputePrimitive_Vector


  SUBROUTINE ComputeConserved_Euler_NonRelativistic &
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

  END SUBROUTINE ComputeConserved_Euler_NonRelativistic


  SUBROUTINE ComputeFromConserved_Euler_NonRelativistic &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      P(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      CALL ComputePrimitive_Euler_NonRelativistic &
             ( U(1:nDOFX,iX1,iX2,iX3,iCF_D),         &
               U(1:nDOFX,iX1,iX2,iX3,iCF_S1),        &
               U(1:nDOFX,iX1,iX2,iX3,iCF_S2),        &
               U(1:nDOFX,iX1,iX2,iX3,iCF_S3),        &
               U(1:nDOFX,iX1,iX2,iX3,iCF_E),         &
               U(1:nDOFX,iX1,iX2,iX3,iCF_Ne),        &
               P(1:nDOFX,iX1,iX2,iX3,iPF_D),         &
               P(1:nDOFX,iX1,iX2,iX3,iPF_V1),        &
               P(1:nDOFX,iX1,iX2,iX3,iPF_V2),        &
               P(1:nDOFX,iX1,iX2,iX3,iPF_V3),        &
               P(1:nDOFX,iX1,iX2,iX3,iPF_E),         &
               P(1:nDOFX,iX1,iX2,iX3,iPF_Ne),        &
               G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_11),  &
               G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_22),  &
               G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputeAuxiliary_Fluid &
             ( P(1:nDOFX,iX1,iX2,iX3,iPF_D ), &
               P(1:nDOFX,iX1,iX2,iX3,iPF_E ), &
               P(1:nDOFX,iX1,iX2,iX3,iPF_Ne), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_P ), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_T ), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_Ye), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_S ), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_E ), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_Gm), &
               A(1:nDOFX,iX1,iX2,iX3,iAF_Cs) )

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved_Euler_NonRelativistic


  SUBROUTINE ComputeTimeStep_Euler_NonRelativistic &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep

    INTEGER  :: iX1, iX2, iX3, iNodeX
    REAL(DP) :: dX(3), dt(3)
    REAL(DP) :: P(nDOFX,nPF)
    REAL(DP) :: Cs(nDOFX)
    REAL(DP) :: EigVals_X1(nCF,nDOFX), alpha_X1, &
                EigVals_X2(nCF,nDOFX), alpha_X2, &
                EigVals_X3(nCF,nDOFX), alpha_X3

    CALL TimersStart_Euler( Timer_Euler_ComputeTimeStep )

    TimeStep = HUGE( One )
    dt       = HUGE( One )

    ! --- Maximum wave-speeds ---

    alpha_X1 = -HUGE( One )
    alpha_X2 = -HUGE( One )
    alpha_X3 = -HUGE( One )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX(1) = MeshX(1) % Width(iX1)
      dX(2) = MeshX(2) % Width(iX2)
      dX(3) = MeshX(3) % Width(iX3)

      CALL ComputePrimitive_Euler_NonRelativistic &
             ( U(:,iX1,iX2,iX3,iCF_D ), U(:,iX1,iX2,iX3,iCF_S1), &
               U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
               U(:,iX1,iX2,iX3,iCF_E ), U(:,iX1,iX2,iX3,iCF_Ne), &
               P(:,iPF_D), P(:,iPF_V1), P(:,iPF_V2), P(:,iPF_V3), &
               P(:,iPF_E), P(:,iPF_Ne), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputeSoundSpeedFromPrimitive &
             ( P(:,iPF_D), P(:,iPF_E), P(:,iPF_Ne), Cs(:) )

      DO iNodeX = 1, nDOFX

        EigVals_X1(:,iNodeX) &
          = Eigenvalues_Euler_NonRelativistic &
              ( P (iNodeX,iPF_V1), &
                Cs(iNodeX),        &
                G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) )

      END DO

      alpha_X1 = MAX( alpha_X1, MAXVAL( ABS( EigVals_X1 ) ) )

      dt(1) = dX(1) / alpha_X1

      IF( nDimsX .GT. 1 )THEN

        DO iNodeX = 1, nDOFX

          EigVals_X2(:,iNodeX) &
            = Eigenvalues_Euler_NonRelativistic &
                ( P (iNodeX,iPF_V2), &
                  Cs(iNodeX),        &
                  G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) )

        END DO

        alpha_X2 = MAX( alpha_X2, MAXVAL( ABS( EigVals_X2 ) ) )

        dt(2) = dX(2) / alpha_X2

      END IF

      IF( nDimsX .GT. 2 )THEN

        DO iNodeX = 1, nDOFX

          EigVals_X3(:,iNodeX) &
            = Eigenvalues_Euler_NonRelativistic &
                ( P (iNodeX,iPF_V3), &
                  Cs(iNodeX),        &
                  G (iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

        END DO

        alpha_X3 = MAX( alpha_X3, MAXVAL( ABS( EigVals_X3 ) ) )

        dt(3) = dX(3) / alpha_X3

      END IF

      TimeStep = MIN( TimeStep, MINVAL( dt ) )

    END DO
    END DO
    END DO

    TimeStep = MAX( CFL * TimeStep, SqrtTiny )

    CALL TimersStop_Euler( Timer_Euler_ComputeTimeStep )

  END SUBROUTINE ComputeTimeStep_Euler_NonRelativistic


  PURE FUNCTION Eigenvalues_Euler_NonRelativistic &
    ( Vi, Cs, Gmii )

    ! --- Vi is the ith contravariant component of the three-velocity
    !     Gmii is the ith covariant component of the spatial three-metric ---

    REAL(DP), INTENT(in)     :: Vi, Cs, Gmii
    REAL(DP), DIMENSION(nCF) :: Eigenvalues_Euler_NonRelativistic

    Eigenvalues_Euler_NonRelativistic &
      = [ Vi - Cs / SQRT( Gmii ), Vi, Vi, Vi, Vi, Vi + Cs / SQRT( Gmii ) ]

    RETURN
  END FUNCTION Eigenvalues_Euler_NonRelativistic


  REAL(DP) FUNCTION AlphaMiddle_Euler_NonRelativistic &
    ( D_L, S_L, E_L, FD_L, FS_L, FE_L, D_R, S_R, E_R, FD_R, FS_R, FE_R, &
      Gmii, aP, aM )

    REAL(DP), INTENT(in) :: D_L, S_L, E_L, FD_L, FS_L, FE_L, &
                            D_R, S_R, E_R, FD_R, FS_R, FE_R, &
                            Gmii, aP, aM

    ! --- Middle Wavespeed as Suggested by Batten et al. (1997) ---
    ! --- (SIAM J. Sci. Comput., Vol. 18, No. 6, pp. 1553-1570) ---

    AlphaMiddle_Euler_NonRelativistic & ! --- Index Up
      = ( aP * S_R + aM * S_L - ( FS_R - FS_L ) ) &
        / ( aP * D_R + aM * D_L - ( FD_R - FD_L ) ) / Gmii

    RETURN
  END FUNCTION AlphaMiddle_Euler_NonRelativistic


  PURE FUNCTION Flux_X1_Euler_NonRelativistic &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Flux_X1_Euler_NonRelativistic(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: VSq

    VSq = Gm_dd_11 * V_1**2 + Gm_dd_22 * V_2**2 + Gm_dd_33 * V_3**2

    Flux_X1_Euler_NonRelativistic(iCF_D ) = D * V_1
    Flux_X1_Euler_NonRelativistic(iCF_S1) = D * Gm_dd_11 * V_1 * V_1 + P
    Flux_X1_Euler_NonRelativistic(iCF_S2) = D * Gm_dd_22 * V_2 * V_1
    Flux_X1_Euler_NonRelativistic(iCF_S3) = D * Gm_dd_33 * V_3 * V_1
    Flux_X1_Euler_NonRelativistic(iCF_E ) = ( E + Half * D * VSq + P ) * V_1
    Flux_X1_Euler_NonRelativistic(iCF_Ne) = Ne * V_1

    RETURN
  END FUNCTION Flux_X1_Euler_NonRelativistic


  PURE FUNCTION Flux_X2_Euler_NonRelativistic &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Flux_X2_Euler_NonRelativistic(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: VSq

    VSq = Gm_dd_11 * V_1**2 + Gm_dd_22 * V_2**2 + Gm_dd_33 * V_3**2

    Flux_X2_Euler_NonRelativistic(iCF_D ) = D * V_2
    Flux_X2_Euler_NonRelativistic(iCF_S1) = D * Gm_dd_11 * V_1 * V_2
    Flux_X2_Euler_NonRelativistic(iCF_S2) = D * Gm_dd_22 * V_2 * V_2 + P
    Flux_X2_Euler_NonRelativistic(iCF_S3) = D * Gm_dd_33 * V_3 * V_2
    Flux_X2_Euler_NonRelativistic(iCF_E ) = ( E + Half * D * VSq + P ) * V_2
    Flux_X2_Euler_NonRelativistic(iCF_Ne) = Ne * V_2

    RETURN
  END FUNCTION Flux_X2_Euler_NonRelativistic


  PURE FUNCTION Flux_X3_Euler_NonRelativistic &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Flux_X3_Euler_NonRelativistic(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    Flux_X3_Euler_NonRelativistic = 0.0_DP

    RETURN
  END FUNCTION Flux_X3_Euler_NonRelativistic


  PURE FUNCTION StressTensor_Diagonal_Euler_NonRelativistic &
    ( S1, S2, S3, V1, V2, V3, P )

    REAL(DP), INTENT(in) :: S1, S2, S3, V1, V2, V3, P

    REAL(DP) :: StressTensor_Diagonal_Euler_NonRelativistic(3)

    StressTensor_Diagonal_Euler_NonRelativistic(1) = S1 * V1 + P
    StressTensor_Diagonal_Euler_NonRelativistic(2) = S2 * V2 + P
    StressTensor_Diagonal_Euler_NonRelativistic(3) = S3 * V3 + P

    RETURN
  END FUNCTION StressTensor_Diagonal_Euler_NonRelativistic


  PURE FUNCTION NumericalFlux_HLL_Euler_NonRelativistic &
    ( uL, uR, fL, fR, aP, aM )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), aP, aM

    REAL(DP) :: NumericalFlux_HLL_Euler_NonRelativistic(nCF)

    NumericalFlux_HLL_Euler_NonRelativistic &
      = ( aP * fL + aM * fR - aP * aM * ( uR - uL ) ) / ( aP + aM )

    RETURN
  END FUNCTION NumericalFlux_HLL_Euler_NonRelativistic


  PURE FUNCTION NumericalFlux_X1_HLLC_Euler_NonRelativistic &
    ( uL, uR, fL, fR, aP, aM, aC, Gm11 )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm11

    REAL(DP) :: NumericalFlux_X1_HLLC_Euler_NonRelativistic(nCF)

    REAL(DP) :: D, V1, V2, V3, P, E, Ne, TMP(nCF)

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X1_HLLC_Euler_NonRelativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X1_HLLC_Euler_NonRelativistic = fR

    ELSE

      IF( aC .GE. Zero )THEN

        TMP = fL + aM * uL

        D  = TMP(iCF_D) / ( aC + aM )
        V1 = aC                       ! --- Index Up
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S1) - Gm11 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC + aM )
        Ne = TMP(iCF_Ne) / ( aC + aM )

      ELSE

        TMP = fR - aP * uR

        D  = TMP(iCF_D) / ( aC - aP )
        V1 = aC                       ! --- Index Up
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S1) - Gm11 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC - aP )
        Ne = TMP(iCF_Ne) / ( aC - aP )

      END IF

      NumericalFlux_X1_HLLC_Euler_NonRelativistic(iCF_D) &
        = D * V1
      NumericalFlux_X1_HLLC_Euler_NonRelativistic(iCF_S1) &
        = D * Gm11 * V1 * V1 + P
      NumericalFlux_X1_HLLC_Euler_NonRelativistic(iCF_S2) &
        = D * V2 * V1
      NumericalFlux_X1_HLLC_Euler_NonRelativistic(iCF_S3) &
        = D * V3 * V1
      NumericalFlux_X1_HLLC_Euler_NonRelativistic(iCF_E) &
        = ( E + P ) * V1
      NumericalFlux_X1_HLLC_Euler_NonRelativistic(iCF_Ne) &
        = Ne * V1

    END IF

    RETURN
  END FUNCTION NumericalFlux_X1_HLLC_Euler_NonRelativistic


  PURE FUNCTION NumericalFlux_X2_HLLC_Euler_NonRelativistic &
    ( uL, uR, fL, fR, aP, aM, aC, Gm22 )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm22

    REAL(DP) :: NumericalFlux_X2_HLLC_Euler_NonRelativistic(nCF)

    REAL(DP) :: D, V1, V2, V3, P, E, Ne, TMP(nCF)

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X2_HLLC_Euler_NonRelativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X2_HLLC_Euler_NonRelativistic = fR

    ELSE

      IF( aC .GE. Zero )THEN

        TMP = fL + aM * uL

        D  = TMP(iCF_D) / ( aC + aM )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = aC                       ! --- Index Up
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S2) - Gm22 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC + aM )
        Ne = TMP(iCF_Ne) / ( aC + aM )

      ELSE

        TMP = fR - aP * uR

        D  = TMP(iCF_D) / ( aC - aP )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = aC                       ! --- Index Up
        V3 = TMP(iCF_S3) / TMP(iCF_D) ! --- Index Down
        P  = TMP(iCF_S2) - Gm22 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC - aP )
        Ne = TMP(iCF_Ne) / ( aC - aP )

      END IF

      NumericalFlux_X2_HLLC_Euler_NonRelativistic(iCF_D) &
        = D * V2
      NumericalFlux_X2_HLLC_Euler_NonRelativistic(iCF_S1) &
        = D * V1 * V2
      NumericalFlux_X2_HLLC_Euler_NonRelativistic(iCF_S2) &
        = D * Gm22 * V2 * V2 + P
      NumericalFlux_X2_HLLC_Euler_NonRelativistic(iCF_S3) &
        = D * V3 * V2
      NumericalFlux_X2_HLLC_Euler_NonRelativistic(iCF_E) &
        = ( E + P ) * V2
      NumericalFlux_X2_HLLC_Euler_NonRelativistic(iCF_Ne) &
        = Ne * V2

    END IF

    RETURN
  END FUNCTION NumericalFlux_X2_HLLC_Euler_NonRelativistic


  PURE FUNCTION NumericalFlux_X3_HLLC_Euler_NonRelativistic &
    ( uL, uR, fL, fR, aP, aM, aC, Gm33 )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm33

    REAL(DP) :: NumericalFlux_X3_HLLC_Euler_NonRelativistic(nCF)

    REAL(DP) :: D, V1, V2, V3, P, E, Ne, TMP(nCF)

    IF( aM .EQ. Zero )THEN

      NumericalFlux_X3_HLLC_Euler_NonRelativistic = fL

    ELSEIF( aP .EQ. Zero )THEN

      NumericalFlux_X3_HLLC_Euler_NonRelativistic = fR

    ELSE

      IF( aC .GE. Zero )THEN

        TMP = fL + aM * uL

        D  = TMP(iCF_D) / ( aC + aM )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = aC                       ! --- Index Up
        P  = TMP(iCF_S3) - Gm33 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC + aM )
        Ne = TMP(iCF_Ne) / ( aC + aM )

      ELSE

        TMP = fR - aP * uR

        D  = TMP(iCF_D) / ( aC - aP )
        V1 = TMP(iCF_S1) / TMP(iCF_D) ! --- Index Down
        V2 = TMP(iCF_S2) / TMP(iCF_D) ! --- Index Down
        V3 = aC                       ! --- Index Up
        P  = TMP(iCF_S3) - Gm33 * aC * TMP(iCF_D)
        E  = ( TMP(iCF_E) - aC * P ) / ( aC - aP )
        Ne = TMP(iCF_Ne) / ( aC - aP )

      END IF

      NumericalFlux_X3_HLLC_Euler_NonRelativistic(iCF_D) &
        = D * V3
      NumericalFlux_X3_HLLC_Euler_NonRelativistic(iCF_S1) &
        = D * V1 * V3
      NumericalFlux_X3_HLLC_Euler_NonRelativistic(iCF_S2) &
        = D * V2 * V3
      NumericalFlux_X3_HLLC_Euler_NonRelativistic(iCF_S3) &
        = D * Gm33 * V3 * V3 + P
      NumericalFlux_X3_HLLC_Euler_NonRelativistic(iCF_E) &
        = ( E + P ) * V3
      NumericalFlux_X3_HLLC_Euler_NonRelativistic(iCF_Ne) &
        = Ne * V3

    END IF

    RETURN
  END FUNCTION NumericalFlux_X3_HLLC_Euler_NonRelativistic

END MODULE Euler_UtilitiesModule_NonRelativistic
