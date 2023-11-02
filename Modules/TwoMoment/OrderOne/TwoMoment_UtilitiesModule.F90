MODULE TwoMoment_UtilitiesModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Three, SqrtTiny, Third
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_h_1, iGF_h_2, iGF_h_3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_TwoMoment
  PUBLIC :: ComputeConserved_TwoMoment
  PUBLIC :: ComputeTimeStep_TwoMoment
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3
  PUBLIC :: StressTensor_Diagonal
  PUBLIC :: NumericalFlux_LLF
  PUBLIC :: ComputeEddingtonTensorComponents_dd

  INTERFACE NumericalFlux_LLF
    MODULE PROCEDURE NumericalFlux_LLF_Scalar
    MODULE PROCEDURE NumericalFlux_LLF_Vector
  END INTERFACE

  INTERFACE Flux_X1
    MODULE PROCEDURE Flux_X1_Scalar
    MODULE PROCEDURE Flux_X1_Vector
  END INTERFACE

CONTAINS


  SUBROUTINE ComputePrimitive_TwoMoment &
    ( N, G_1, G_2, G_3, D, I_1, I_2, I_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), DIMENSION(:), INTENT(in)  :: N, G_1, G_2, G_3
    REAL(DP), DIMENSION(:), INTENT(out) :: D, I_1, I_2, I_3
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    D   = N
    I_1 = G_1 / Gm_dd_11
    I_2 = G_2 / Gm_dd_22
    I_3 = G_3 / Gm_dd_33

  END SUBROUTINE ComputePrimitive_TwoMoment


  SUBROUTINE ComputeConserved_TwoMoment &
    ( D, I_1, I_2, I_3, N, G_1, G_2, G_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, I_1, I_2, I_3
    REAL(DP), INTENT(out) :: N, G_1, G_2, G_3
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    N   = D
    G_1 = Gm_dd_11 * I_1
    G_2 = Gm_dd_22 * I_2
    G_3 = Gm_dd_33 * I_3

  END SUBROUTINE ComputeConserved_TwoMoment


  SUBROUTINE ComputeTimeStep_TwoMoment &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GX, CFL, TimeStep )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      GX(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep

    INTEGER  :: iX1, iX2, iX3, iNodeX
    REAL(DP) :: dX(3), dt(3)

    TimeStep = HUGE( One )
    dt       = HUGE( One )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX(1) = MeshX(1) % Width(iX1)
      dX(2) = MeshX(2) % Width(iX2)
      dX(3) = MeshX(3) % Width(iX3)

      dt(1) = dX(1) * MINVAL( GX(:,iX1,iX2,iX3,iGF_h_1) )

      IF( nDimsX .GT. 1 ) &
        dt(2) = dX(2) * MINVAL( GX(:,iX1,iX2,iX3,iGF_h_2) )

      IF( nDimsX .GT. 2 ) &
        dt(3) = dX(3) * MINVAL( GX(:,iX1,iX2,iX3,iGF_h_3) )

      TimeStep = MIN( TimeStep, MINVAL( dt ) )

    END DO
    END DO
    END DO

    TimeStep = MAX( CFL * TimeStep, SqrtTiny )

  END SUBROUTINE ComputeTimeStep_TwoMoment


  FUNCTION Flux_X1_Scalar &
    ( D, I_1, I_2, I_3, FF, EF, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
    RESULT( Flux_X1 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D, I_1, I_2, I_3, FF, EF
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP)             :: Flux_X1(1:4)

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
  END FUNCTION Flux_X1_Scalar

  FUNCTION Flux_X1_Vector &
    ( D, I_1, I_2, I_3, FF, EF, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
    RESULT( Flux_X1 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D(:), I_1(:), I_2(:), I_3(:), FF(:), EF(:)
    REAL(DP), INTENT(in) :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)
    REAL(DP)             :: Flux_X1(SIZE(D),1:4)

    REAL(DP) :: h_u_1(SIZE(D))
    REAL(DP) :: h_d_1(SIZE(D)), h_d_2(SIZE(D)), h_d_3(SIZE(D))

    h_u_1 = I_1 / ( FF * D )

    h_d_1 = Gm_dd_11 * I_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_3 / ( FF * D )

    Flux_X1(:,1) &
      = I_1

    Flux_X1(:,2) &
      = D * Half * ( (Three*EF - One) * h_u_1 * h_d_1 + (One - EF) )

    Flux_X1(:,3) &
      = D * Half * ( (Three*EF - One) * h_u_1 * h_d_2 )

    Flux_X1(:,4) &
      = D * Half * ( (Three*EF - One) * h_u_1 * h_d_3 )

    RETURN
  END FUNCTION Flux_X1_Vector


  FUNCTION Flux_X2 &
    ( D, I_1, I_2, I_3, FF, EF, Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

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


  FUNCTION Flux_X3 &
    ( D, I_1, I_2, I_3, FF, EF, Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

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


  FUNCTION StressTensor_Diagonal &
    ( D, I_1, I_2, I_3, FF, EF, Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP)             :: StressTensor_Diagonal(1:3)
    REAL(DP), INTENT(in) :: D, I_1, I_2, I_3, FF, EF
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: h_u_1, h_u_2, h_u_3
    REAL(DP) :: h_d_1, h_d_2, h_d_3

    h_u_1 = I_1 / ( FF * D )
    h_u_2 = I_2 / ( FF * D )
    h_u_3 = I_3 / ( FF * D )

    h_d_1 = Gm_dd_11 * h_u_1
    h_d_2 = Gm_dd_22 * h_u_2
    h_d_3 = Gm_dd_33 * h_u_3

    StressTensor_Diagonal(1) &
      = D * Half * ( (Three*EF - One) * h_u_1 * h_d_1 + (One - EF) )

    StressTensor_Diagonal(2) &
      = D * Half * ( (Three*EF - One) * h_u_2 * h_d_2 + (One - EF) )

    StressTensor_Diagonal(3) &
      = D * Half * ( (Three*EF - One) * h_u_3 * h_d_3 + (One - EF) )

    RETURN
  END FUNCTION StressTensor_Diagonal


  FUNCTION NumericalFlux_LLF_Scalar &
    ( u_L, u_R, Flux_L, Flux_R, alpha ) &
    RESULT( NumericalFlux_LLF )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Local Lax-Friedrichs Flux ---

    REAL(DP), INTENT(in) :: u_L, u_R, flux_L, flux_R, alpha
    REAL(DP) :: NumericalFlux_LLF

    NumericalFlux_LLF &
      = Half * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF_Scalar

  FUNCTION NumericalFlux_LLF_Vector &
    ( u_L, u_R, Flux_L, Flux_R, alpha ) &
    RESULT( NumericalFlux_LLF )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Local Lax-Friedrichs Flux ---

    REAL(DP), INTENT(in) :: u_L(:), u_R(:), flux_L(:), flux_R(:), alpha(:)
    REAL(DP) :: NumericalFlux_LLF(SIZE(u_L))

    NumericalFlux_LLF &
      = Half * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF_Vector


  SUBROUTINE ComputeEddingtonTensorComponents_dd &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(out) :: &
      k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33

    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_d_1, h_d_2, h_d_3

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    IF ( FF <= SqrtTiny ) THEN
      EF = Third
      a = Third
      b = Zero
    ELSE
      EF = EddingtonFactor( D, FF )
      a = Half * ( One - EF )
      b = Half * ( Three * EF - One )
    END IF

    h_d_1 = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_u_3 / ( FF * D )

    ! --- Diagonal Eddington Tensor Components ---

    k_dd_11 = a * Gm_dd_11 + b * h_d_1 * h_d_1
    k_dd_22 = a * Gm_dd_22 + b * h_d_2 * h_d_2
    k_dd_33 = a * Gm_dd_33 + b * h_d_3 * h_d_3

    ! --- Off-Diagonal Eddington Tensor Components ---

    k_dd_12 = b * h_d_1 * h_d_2
    k_dd_13 = b * h_d_1 * h_d_3
    k_dd_23 = b * h_d_2 * h_d_3

  END SUBROUTINE ComputeEddingtonTensorComponents_dd


END MODULE TwoMoment_UtilitiesModule
