MODULE Euler_UtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE FluidFieldsModule, ONLY: &
    nCF

#ifdef HYDRO_NONRELATIVISTIC

  USE Euler_UtilitiesModule_NonRelativistic

#elif HYDRO_RELATIVISTIC

  USE Euler_UtilitiesModule_Relativistic

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
    ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33, &
      AF_P_Option )

    REAL(DP), INTENT(in)  :: CF_D(:), CF_S1(:), CF_S2(:), CF_S3(:), &
                             CF_E(:), CF_Ne(:)
    REAL(DP), INTENT(out) :: PF_D(:), PF_V1(:), PF_V2(:), PF_V3(:), &
                             PF_E(:), PF_Ne(:)
    REAL(DP), INTENT(in)  :: GF_Gm_dd_11(:), GF_Gm_dd_22(:), GF_Gm_dd_33(:)

    REAL(DP), INTENT(out), OPTIONAL :: AF_P_Option(:)

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ComputePrimitive_NonRelativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#elif HYDRO_RELATIVISTIC

    CALL Euler_ComputePrimitive_Relativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33, AF_P_Option )

#endif

  END SUBROUTINE Euler_ComputePrimitive


  SUBROUTINE Euler_ComputeConserved &
    ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
      AF_P )

    REAL(DP), INTENT(in)  :: PF_D(:), PF_V1(:), PF_V2(:), PF_V3(:), &
                             PF_E(:), PF_Ne(:)
    REAL(DP), INTENT(out) :: CF_D(:), CF_S1(:), CF_S2(:), CF_S3(:), &
                             CF_E(:), CF_Ne(:)
    REAL(DP), INTENT(in)  :: GF_Gm_dd_11(:), GF_Gm_dd_22(:), GF_Gm_dd_33(:)

    REAL(DP), INTENT(in), OPTIONAL :: AF_P(:)

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeConserved_NonRelativistic &
           ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#elif HYDRO_RELATIVISTIC

    CALL Euler_ComputeConserved_Relativistic &
           ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33, &
             AF_P )

#endif

  END SUBROUTINE Euler_ComputeConserved


  SUBROUTINE Euler_ComputeFromConserved( iX_B, iX_E, G, U, P, A )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      G(:,iX_B(1):,iX_B(2):,iX_B(3):,:), &
      U(:,iX_B(1):,iX_B(2):,iX_B(3):,:)
    REAL(DP), INTENT(out) :: &
      P(:,iX_B(1):,iX_B(2):,iX_B(3):,:), &
      A(:,iX_B(1):,iX_B(2):,iX_B(3):,:)

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeFromConserved_NonRelativistic( iX_B, iX_E, G, U, P, A )

#elif HYDRO_RELATIVISTIC

    CALL Euler_ComputeFromConserved_Relativistic( iX_B, iX_E, G, U, P, A )

#endif

  END SUBROUTINE Euler_ComputeFromConserved


  SUBROUTINE Euler_ComputeTimeStep &
    ( iX_B, iX_E, G, U, CFL, TimeStep, UseSourceTerm_Option )

    INTEGER,  INTENT(in)          :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)          :: &
      G(:,iX_B(1):,iX_B(2):,iX_B(3):,:), &
      U(:,iX_B(1):,iX_B(2):,iX_B(3):,:)
    REAL(DP), INTENT(in)          :: &
      CFL
    REAL(DP), INTENT(out)         :: &
      TimeStep
    LOGICAL, INTENT(in), OPTIONAL :: &
      UseSourceTerm_Option

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeTimeStep_NonRelativistic &
           ( iX_B, iX_E, G, U, CFL, TimeStep )

#elif HYDRO_RELATIVISTIC

    CALL Euler_ComputeTimeStep_Relativistic &
           ( iX_B, iX_E, G, U, CFL, TimeStep, UseSourceTerm_Option )

#endif

  END SUBROUTINE Euler_ComputeTimeStep


  PURE FUNCTION Euler_Eigenvalues &
    ( V, Cs, V1, V2, V3, Gm, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in) :: V, Cs
    REAL(DP), INTENT(in), OPTIONAL :: V1, V2, V3, Gm, Gm11, Gm22, Gm33, &
                                      Lapse, Shift

    REAL(DP) :: Euler_Eigenvalues(nCF)

#ifdef HYDRO_NONRELATIVISTIC

    Euler_Eigenvalues = Euler_Eigenvalues_NonRelativistic( V, Cs )

#elif HYDRO_RELATIVISTIC

    Euler_Eigenvalues = Euler_Eigenvalues_Relativistic &
                          ( V, Cs, V1, V2, V3, Gm, Gm11, Gm22, Gm33, &
                            Lapse, Shift )

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


  REAL(DP) FUNCTION Euler_AlphaMiddle &
    ( DL, SL, EL, F_DL, F_SL, F_EL, DR, SR, ER, F_DR, F_SR, F_ER, &
      aP, aM, Gm, Lapse, Shift )

    ! --- Gm is the covariant ii-component of the spatial three-metric
    !     Shift is the ith contravariant component of the shift-vector ---

    REAL(DP), INTENT(in) :: DL, SL, EL, F_DL, F_SL, F_EL, &
                            DR, SR, ER, F_DR, F_SR, F_ER, &
                            aP, aM, Gm

    REAL(DP), INTENT(in), OPTIONAL :: Lapse, Shift

#ifdef HYDRO_NONRELATIVISTIC

    Euler_AlphaMiddle = Euler_AlphaMiddle_NonRelativistic &
                          ( DL, F_DL, SL, F_SL, EL, F_EL, &
                            DR, F_DR, SR, F_SR, ER, F_ER, &
                            aP, aM, Gm )

#elif HYDRO_RELATIVISTIC

    Euler_AlphaMiddle = Euler_AlphaMiddle_Relativistic &
                          ( DL, F_DL, SL, F_SL, EL, F_EL, &
                            DR, F_DR, SR, F_SR, ER, F_ER, &
                            aP, aM, Gm, Lapse, Shift )

#endif

    RETURN
  END FUNCTION Euler_AlphaMiddle


  PURE FUNCTION Euler_Flux_X1 &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    ! --- Shift is the first contravariant component of the shift-vector ---

    REAL(DP)             :: Euler_Flux_X1(nCF)
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    REAL(DP), INTENT(in), OPTIONAL :: Lapse, Shift

#ifdef HYDRO_NONRELATIVISTIC

    Euler_Flux_X1 = Euler_Flux_X1_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#elif HYDRO_RELATIVISTIC

    Euler_Flux_X1 = Euler_Flux_X1_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )

#endif

    RETURN
  END FUNCTION Euler_Flux_X1


  PURE FUNCTION Euler_Flux_X2 &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    ! --- Shift is the second contravariant component of the shift-vector ---

    REAL(DP)             :: Euler_Flux_X2(nCF)
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    REAL(DP), INTENT(in), OPTIONAL :: Lapse, Shift

#ifdef HYDRO_NONRELATIVISTIC

    Euler_Flux_X2 = Euler_Flux_X2_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#elif HYDRO_RELATIVISTIC

    Euler_Flux_X2 = Euler_Flux_X2_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )

#endif

    RETURN
  END FUNCTION Euler_Flux_X2


  PURE FUNCTION Euler_Flux_X3 &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    ! --- Shift is the third contravariant component of the shift-vector ---

    REAL(DP)             :: Euler_Flux_X3(nCF)
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    REAL(DP), INTENT(in), OPTIONAL :: Lapse, Shift

#ifdef HYDRO_NONRELATIVISTIC

    Euler_Flux_X3 = Euler_Flux_X3_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#elif HYDRO_RELATIVISTIC

    Euler_Flux_X3 = Euler_Flux_X3_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )

#endif

    RETURN
  END FUNCTION Euler_Flux_X3


  PURE FUNCTION Euler_StressTensor_Diagonal( S1, S2, S3, V1, V2, V3, P )

    REAL(DP), INTENT(in) :: S1, S2, S3, V1, V2, V3, P

    REAL(DP) :: Euler_StressTensor_Diagonal(3)

#ifdef HYDRO_NONRELATIVISTIC

    Euler_StressTensor_Diagonal &
      = Euler_StressTensor_Diagonal_NonRelativistic( S1, S2, S3, V1, V2, V3, P )

#elif HYDRO_RELATIVISTIC

    Euler_StressTensor_Diagonal &
      = Euler_StressTensor_Diagonal_Relativistic( S1, S2, S3, V1, V2, V3, P )

#endif

    RETURN
  END FUNCTION Euler_StressTensor_Diagonal


  PURE FUNCTION Euler_NumericalFlux_HLL &
    ( uL, uR, fL, fR, aP, aM, aC, Gm, vL, vR, pL, pR, Lapse, Shift )

    REAL(DP), INTENT(in)           :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                                      aP, aM, aC, Gm
    REAL(DP), INTENT(in), OPTIONAL :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_HLL(nCF)

#ifdef HYDRO_NONRELATIVISTIC

    Euler_NumericalFlux_HLL = Euler_NumericalFlux_HLL_NonRelativistic &
                                ( uL, uR, fL, fR, aP, aM, aC, Gm )

#elif HYDRO_RELATIVISTIC

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_HLL


  PURE FUNCTION Euler_NumericalFlux_X1_HLLC &
    ( uL, uR, fL, fR, aP, aM, aC, Gm, vL, vR, pL, pR, Lapse, Shift )

    REAL(DP), INTENT(in)           :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                                      aP, aM, aC, Gm
    REAL(DP), INTENT(in), OPTIONAL :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_X1_HLLC(nCF)

#ifdef HYDRO_NONRELATIVISTIC

    Euler_NumericalFlux_X1_HLLC = Euler_NumericalFlux_X1_HLLC_NonRelativistic &
                                    ( uL, uR, fL, fR, aP, aM, aC, Gm )

#elif HYDRO_RELATIVISTIC

    Euler_NumericalFlux_X1_HLLC = Euler_NumericalFlux_X1_HLLC_Relativistic &
                                    ( uL, uR, fL, fR, aP, aM, aC, Gm, &
                                      vL, vR, pL, pR, Lapse, Shift )

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_X1_HLLC


  PURE FUNCTION Euler_NumericalFlux_X2_HLLC &
    ( uL, uR, fL, fR, aP, aM, aC, Gm, vL, vR, pL, pR, Lapse, Shift )

    REAL(DP), INTENT(in)           :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                                      aP, aM, aC, Gm
    REAL(DP), INTENT(in), OPTIONAL :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_X2_HLLC(nCF)

#ifdef HYDRO_NONRELATIVISTIC

    Euler_NumericalFlux_X2_HLLC = Euler_NumericalFlux_X2_HLLC_NonRelativistic &
                                    ( uL, uR, fL, fR, aP, aM, aC, Gm )

#elif HYDRO_RELATIVISTIC


    Euler_NumericalFlux_X2_HLLC = Euler_NumericalFlux_X2_HLLC_Relativistic &
                                    ( uL, uR, fL, fR, aP, aM, aC, Gm, &
                                      vL, vR, pL, pR, Lapse, Shift )

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_X2_HLLC


  PURE FUNCTION Euler_NumericalFlux_X3_HLLC &
    ( uL, uR, fL, fR, aP, aM, aC, Gm, vL, vR, pL, pR, Lapse, Shift )

    REAL(DP), INTENT(in)           :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                                      aP, aM, aC, Gm
    REAL(DP), INTENT(in), OPTIONAL :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_X3_HLLC(nCF)

#ifdef HYDRO_NONRELATIVISTIC

    Euler_NumericalFlux_X3_HLLC = Euler_NumericalFlux_X3_HLLC_NonRelativistic &
                                    ( uL, uR, fL, fR, aP, aM, aC, Gm )

#elif HYDRO_RELATIVISTIC


    Euler_NumericalFlux_X3_HLLC = Euler_NumericalFlux_X3_HLLC_Relativistic &
                                    ( uL, uR, fL, fR, aP, aM, aC, Gm, &
                                      vL, vR, pL, pR, Lapse, Shift )

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_X3_HLLC


END MODULE Euler_UtilitiesModule
