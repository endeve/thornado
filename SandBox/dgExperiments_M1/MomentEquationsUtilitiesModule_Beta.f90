MODULE MomentEquationsUtilitiesModule_Beta

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    SqrtTiny, &
    Half, &
    One, &
    Three

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive
  PUBLIC :: ComputeConserved
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3
  PUBLIC :: NumericalFlux_LLF

CONTAINS


  SUBROUTINE ComputePrimitive &
    ( N, G_1, G_2, G_3, D, I_1, I_2, I_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:), INTENT(in)  :: N, G_1, G_2, G_3
    REAL(DP), DIMENSION(:), INTENT(out) :: D, I_1, I_2, I_3
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    D   = N
    I_1 = G_1 / Gm_dd_11
    I_2 = G_2 / Gm_dd_22
    I_3 = G_3 / Gm_dd_33

  END SUBROUTINE ComputePrimitive


  SUBROUTINE ComputeConserved &
    ( D, I_1, I_2, I_3, N, G_1, G_2, G_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, I_1, I_2, I_3
    REAL(DP), DIMENSION(:), INTENT(out) :: N, G_1, G_2, G_3
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    N   = D
    G_1 = Gm_dd_11 * I_1
    G_2 = Gm_dd_22 * I_2
    G_3 = Gm_dd_33 * I_3

  END SUBROUTINE ComputeConserved


  PURE FUNCTION Flux_X1 &
    ( D, I_1, I_2, I_3, FF, EF, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Flux_X1(1:4)
    REAL(DP), INTENT(in) :: D, I_1, I_2, I_3, FF, EF
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: h_u_1
    REAL(DP) :: h_d_1, h_d_2, h_d_3

    h_u_1 = I_1 / ( FF * D )

    h_d_1 = Gm_dd_11 * I_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_3 / ( FF * D )

    Flux_X1(1) &
      = I_1

    Flux_X1(2) &
      = D * Half * ( (Three*EF - One) * h_u_1 * h_d_1 + (One - EF) )

    Flux_X1(3) &
      = D * Half * ( (Three*EF - One) * h_u_1 * h_d_2 )

    Flux_X1(4) &
      = D * Half * ( (Three*EF - One) * h_u_1 * h_d_3 )

    RETURN
  END FUNCTION Flux_X1


  PURE FUNCTION Flux_X2 &
    ( D, I_1, I_2, I_3, FF, EF, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Flux_X2(1:4)
    REAL(DP), INTENT(in) :: D, I_1, I_2, I_3, FF, EF
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: h_u_2
    REAL(DP) :: h_d_1, h_d_2, h_d_3

    h_u_2 = I_2 / ( FF * D )

    h_d_1 = Gm_dd_11 * I_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_3 / ( FF * D )

    Flux_X2(1) &
      = I_2

    Flux_X2(2) &
      = D * Half * ( (Three*EF - One) * h_u_2 * h_d_1 )

    Flux_X2(3) &
      = D * Half * ( (Three*EF - One) * h_u_2 * h_d_2 + (One - EF) )

    Flux_X2(4) &
      = D * Half * ( (Three*EF - One) * h_u_2 * h_d_3 )

    RETURN
  END FUNCTION Flux_X2


  PURE FUNCTION Flux_X3 &
    ( D, I_1, I_2, I_3, FF, EF, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Flux_X3(1:4)
    REAL(DP), INTENT(in) :: D, I_1, I_2, I_3, FF, EF
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: h_u_3
    REAL(DP) :: h_d_1, h_d_2, h_d_3

    h_u_3 = I_3 / ( FF * D )

    h_d_1 = Gm_dd_11 * I_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_3 / ( FF * D )

    Flux_X3(1) &
      = I_3

    Flux_X3(2) &
      = D * Half * ( (Three*EF - One) * h_u_3 * h_d_1 )

    Flux_X3(3) &
      = D * Half * ( (Three*EF - One) * h_u_3 * h_d_2 )

    Flux_X3(4) &
      = D * Half * ( (Three*EF - One) * h_u_3 * h_d_3 + (One - EF) )

    RETURN
  END FUNCTION Flux_X3


  REAL(DP) PURE ELEMENTAL FUNCTION NumericalFlux_LLF &
    ( u_L, u_R, Flux_L, Flux_R, alpha )

    ! --- Local Lax-Friedrichs Flux ---

    REAL(DP), INTENT(in) :: u_L, u_R, flux_L, flux_R, alpha

    NumericalFlux_LLF &
      = Half * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF


END MODULE MomentEquationsUtilitiesModule_Beta
