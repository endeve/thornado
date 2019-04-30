MODULE Euler_UtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE FluidFieldsModule, ONLY: &
    nCF

#ifdef HYDRO_NONRELATIVISTIC

  USE Euler_UtilitiesModule_NonRelativistic

#elif HYDRO_RELATIVISTIC

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_ComputePrimitive
  PUBLIC :: Euler_ComputeConserved
  PUBLIC :: Euler_ComputeFromConserved
  PUBLIC :: Euler_ComputeTimeStep
  PUBLIC :: Euler_Eigenvalues
  PUBLIC :: Euler_AlphaPlus
  PUBLIC :: Euler_AlphaMinus
  PUBLIC :: Euler_AlphaMiddle
  PUBLIC :: Euler_Flux_X1
  PUBLIC :: Euler_Flux_X2
  PUBLIC :: Euler_Flux_X3
  PUBLIC :: Euler_StressTensor_Diagonal
  PUBLIC :: Euler_NumericalFlux_HLL
  PUBLIC :: Euler_NumericalFlux_X1_HLLC
  PUBLIC :: Euler_NumericalFlux_X2_HLLC
  PUBLIC :: Euler_NumericalFlux_X3_HLLC


CONTAINS


  SUBROUTINE Euler_ComputePrimitive &
    ( N, S_1, S_2, S_3, G, Ne, D, V_1, V_2, V_3, E, De, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:), INTENT(in)  :: N, S_1, S_2, S_3, G, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: D, V_1, V_2, V_3, E, De
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ComputePrimitive_NonRelativistic &
           ( N, S_1, S_2, S_3, G, Ne, D, V_1, V_2, V_3, E, De, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif HYDRO_RELATIVISTIC

#endif

  END SUBROUTINE Euler_ComputePrimitive


  SUBROUTINE Euler_ComputeConserved &
    ( D, V_1, V_2, V_3, E, De, N, S_1, S_2, S_3, G, Ne, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, V_1, V_2, V_3, E, De
    REAL(DP), DIMENSION(:), INTENT(out) :: N, S_1, S_2, S_3, G, Ne
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeConserved_NonRelativistic &
           ( D, V_1, V_2, V_3, E, De, N, S_1, S_2, S_3, G, Ne, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif HYDRO_RELATIVISTIC

#endif

  END SUBROUTINE Euler_ComputeConserved


  SUBROUTINE Euler_ComputeFromConserved( iX_B, iX_E, G, U, P, A )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:), &
      U(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(DP), INTENT(out) :: &
      P(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:), &
      A(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeFromConserved_NonRelativistic( iX_B, iX_E, G, U, P, A )

#elif HYDRO_RELATIVISTIC

#endif

  END SUBROUTINE Euler_ComputeFromConserved


  SUBROUTINE Euler_ComputeTimeStep( iX_B, iX_E, G, U, CFL, TimeStep )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:), &
      U(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeTimeStep_NonRelativistic &
           ( iX_B, iX_E, G, U, CFL, TimeStep )

#elif HYDRO_RELATIVISTIC

#endif

  END SUBROUTINE Euler_ComputeTimeStep


  PURE FUNCTION Euler_Eigenvalues( V, Cs )

    REAL(DP), INTENT(in) :: V, Cs
    REAL(DP), DIMENSION(nCF) :: Euler_Eigenvalues

#ifdef HYDRO_NONRELATIVISTIC

    Euler_Eigenvalues = Euler_Eigenvalues_NonRelativistic( V, Cs )

#elif HYDRO_RELATIVISTIC

#endif

  END FUNCTION Euler_Eigenvalues


  PURE ELEMENTAL REAL(DP) FUNCTION Euler_AlphaPlus &
    ( V_L, Cs_L, V_R, Cs_R, Gm_dd_ii )

    REAL(DP), INTENT(in) :: V_L, Cs_L, V_R, Cs_R, Gm_dd_ii

#ifdef HYDRO_NONRELATIVISTIC

    Euler_AlphaPlus = Euler_AlphaPlus_NonRelativistic &
                        ( V_L, Cs_L, V_R, Cs_R, Gm_dd_ii )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_AlphaPlus


  PURE ELEMENTAL REAL(DP) FUNCTION Euler_AlphaMinus &
    ( V_L, Cs_L, V_R, Cs_R, Gm_dd_ii )

    REAL(DP), INTENT(in) :: V_L, Cs_L, V_R, Cs_R, Gm_dd_ii

#ifdef HYDRO_NONRELATIVISTIC

    Euler_AlphaMinus = Euler_AlphaMinus_NonRelativistic &
                         ( V_L, Cs_L, V_R, Cs_R, Gm_dd_ii )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_AlphaMinus


  PURE ELEMENTAL REAL(DP) FUNCTION Euler_AlphaMiddle &
    ( D_L, FD_L, S_L, FS_L, D_R, FD_R, S_R, FS_R, aP, aM, Gm_dd_ii )

    REAL(DP), INTENT(in) :: D_L, FD_L, S_L, FS_L
    REAL(DP), INTENT(in) :: D_R, FD_R, S_R, FS_R
    REAL(DP), INTENT(in) :: aP, aM, Gm_dd_ii

#ifdef HYDRO_NONRELATIVISTIC

    Euler_AlphaMiddle = Euler_AlphaMiddle_NonRelativistic &
                          ( D_L, FD_L, S_L, FS_L, D_R, FD_R, S_R, FS_R, &
                            aP, aM, Gm_dd_ii )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_AlphaMiddle


  PURE FUNCTION Euler_Flux_X1 &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Euler_Flux_X1(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

#ifdef HYDRO_NONRELATIVISTIC

    Euler_Flux_X1 = Euler_Flux_X1_NonRelativistic &
                      ( D, V_1, V_2, V_3, E, Ne, P, &
                        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_Flux_X1


  PURE FUNCTION Euler_Flux_X2 &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Euler_Flux_X2(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

#ifdef HYDRO_NONRELATIVISTIC

    Euler_Flux_X2 = Euler_Flux_X2_NonRelativistic &
                      ( D, V_1, V_2, V_3, E, Ne, P, &
                        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_Flux_X2


  PURE FUNCTION Euler_Flux_X3 &
    ( D, V_1, V_2, V_3, E, Ne, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Euler_Flux_X3(1:nCF)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

#ifdef HYDRO_NONRELATIVISTIC

    Euler_Flux_X3 = Euler_Flux_X3_NonRelativistic &
                      ( D, V_1, V_2, V_3, E, Ne, P, &
                        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_Flux_X3


  PURE FUNCTION Euler_StressTensor_Diagonal &
    ( D, V_1, V_2, V_3, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP)             :: Euler_StressTensor_Diagonal(1:3)
    REAL(DP), INTENT(in) :: D, V_1, V_2, V_3, P
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

#ifdef HYDRO_NONRELATIVISTIC

    Euler_StressTensor_Diagonal &
      = Euler_StressTensor_Diagonal_NonRelativistic &
          ( D, V_1, V_2, V_3, P, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_StressTensor_Diagonal


  PURE FUNCTION Euler_NumericalFlux_HLL &
    ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd )

    REAL(DP),                 INTENT(in) :: aP, aM, aC, Gm_dd
    REAL(DP), DIMENSION(nCF), INTENT(in) :: uL, FluxL, uR, FluxR
    REAL(DP), DIMENSION(nCF)             :: Euler_NumericalFlux_HLL

#ifdef HYDRO_NONRELATIVISTIC

    Euler_NumericalFlux_HLL = Euler_NumericalFlux_HLL_NonRelativistic &
                                ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_HLL


  PURE FUNCTION Euler_NumericalFlux_X1_HLLC &
    ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd )

    REAL(DP),                 INTENT(in) :: aP, aM, aC, Gm_dd
    REAL(DP), DIMENSION(nCF), INTENT(in) :: uL, FluxL, uR, FluxR
    REAL(DP), DIMENSION(nCF)             :: Euler_NumericalFlux_X1_HLLC

#ifdef HYDRO_NONRELATIVISTIC

    Euler_NumericalFlux_X1_HLLC = Euler_NumericalFlux_X1_HLLC_NonRelativistic &
                                    ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_X1_HLLC


  PURE FUNCTION Euler_NumericalFlux_X2_HLLC &
    ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd )

    REAL(DP),                 INTENT(in) :: aP, aM, aC, Gm_dd
    REAL(DP), DIMENSION(nCF), INTENT(in) :: uL, FluxL, uR, FluxR
    REAL(DP), DIMENSION(nCF)             :: Euler_NumericalFlux_X2_HLLC

#ifdef HYDRO_NONRELATIVISTIC

    Euler_NumericalFlux_X2_HLLC = Euler_NumericalFlux_X2_HLLC_NonRelativistic &
                                    ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_X2_HLLC


  PURE FUNCTION Euler_NumericalFlux_X3_HLLC &
    ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd )

    REAL(DP),                 INTENT(in) :: aP, aM, aC, Gm_dd
    REAL(DP), DIMENSION(nCF), INTENT(in) :: uL, FluxL, uR, FluxR
    REAL(DP), DIMENSION(nCF)             :: Euler_NumericalFlux_X3_HLLC

#ifdef HYDRO_NONRELATIVISTIC

    Euler_NumericalFlux_X3_HLLC = Euler_NumericalFlux_X3_HLLC_NonRelativistic &
                                    ( uL, FluxL, uR, FluxR, aP, aM, aC, Gm_dd )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_X3_HLLC


END MODULE Euler_UtilitiesModule
