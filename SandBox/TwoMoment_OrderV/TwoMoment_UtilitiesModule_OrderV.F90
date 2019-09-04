MODULE TwoMoment_UtilitiesModule_OrderV

  USE KindModule, ONLY: &
    DP, Half, One, Three
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_TwoMoment
  PUBLIC :: ComputeConserved_TwoMoment

CONTAINS


  SUBROUTINE ComputePrimitive_TwoMoment &
    ( N, G1, G2, G3, D, I1, I2, I3, V1, V2, V3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in)  :: N, G1, G2, G3 ! --- Index Down
    REAL(DP), INTENT(out) :: D, I1, I2, I3 ! --- Index Up
    REAL(DP), INTENT(in)  ::    V1, V2, V3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

  END SUBROUTINE ComputePrimitive_TwoMoment


  SUBROUTINE ComputeConserved_TwoMoment &
    ( D, I_u_1, I_u_2, I_u_3, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in)  :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(out) :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_d_1, h_d_2, h_d_3
    REAL(DP) :: K_dd_11, K_dd_12, K_dd_13, K_dd_22, K_dd_23, K_dd_33

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    h_d_1 = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_u_3 / ( FF * D )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    ! --- Diagonal Eddington Tensor Components ---

    K_dd_11 = ( a * Gm_dd_11 + b * h_d_1 * h_d_1 ) * D
    K_dd_22 = ( a * Gm_dd_22 + b * h_d_2 * h_d_2 ) * D
    K_dd_33 = ( a * Gm_dd_33 + b * h_d_3 * h_d_3 ) * D

    ! --- Off-Diagonal Eddington Tensor Components ---

    K_dd_12 = b * h_d_1 * h_d_2 * D
    K_dd_13 = b * h_d_1 * h_d_3 * D
    K_dd_23 = b * h_d_2 * h_d_3 * D

    ! --- Conserved Number Density ---

    N = D + Gm_dd_11 * V_u_1 * I_u_1 &
          + Gm_dd_22 * V_u_2 * I_u_2 &
          + Gm_dd_33 * V_u_3 * I_u_3

    ! --- Conserved Number Flux Density (1) ---

    G_d_1 = Gm_dd_11 * I_u_1 &
              + V_u_1 * K_dd_11 + V_u_2 * K_dd_12 + V_u_3 * K_dd_13

    ! --- Conserved Number Flux Density (2) ---

    G_d_2 = Gm_dd_22 * I_u_2 &
              + V_u_1 * K_dd_12 + V_u_2 * K_dd_22 + V_u_3 * K_dd_23

    ! --- Conserved Number Flux Density (3) ---

    G_d_3 = Gm_dd_33 * I_u_3 &
              + V_u_1 * K_dd_13 + V_u_2 * K_dd_23 + V_u_3 * K_dd_33

  END SUBROUTINE ComputeConserved_TwoMoment


END MODULE TwoMoment_UtilitiesModule_OrderV
