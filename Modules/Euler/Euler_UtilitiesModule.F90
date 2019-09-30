MODULE Euler_UtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule, ONLY: &
    nCF, nPF, nAF
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_ComputeTimeStep

#if defined HYDRO_NONRELATIVISTIC

  USE Euler_UtilitiesModule_NonRelativistic

#elif defined HYDRO_RELATIVISTIC

  USE Euler_UtilitiesModule_Relativistic

#else

  USE Euler_UtilitiesModule_NonRelativistic

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_ComputePrimitive
  PUBLIC :: Euler_ComputeConserved
  PUBLIC :: Euler_ComputeFromConserved
  PUBLIC :: Euler_ComputeTimeStep
  PUBLIC :: Euler_Eigenvalues
  PUBLIC :: Euler_AlphaMiddle
  PUBLIC :: Euler_Flux_X1
  PUBLIC :: Euler_Flux_X2
  PUBLIC :: Euler_Flux_X3
  PUBLIC :: Euler_StressTensor_Diagonal
  PUBLIC :: Euler_NumericalFlux_X1
  PUBLIC :: Euler_NumericalFlux_X2
  PUBLIC :: Euler_NumericalFlux_X3


CONTAINS


  SUBROUTINE Euler_ComputePrimitive &
    ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

    REAL(DP), INTENT(in)  :: CF_D(:), CF_S1(:), CF_S2(:), CF_S3(:), &
                             CF_E(:), CF_Ne(:)
    REAL(DP), INTENT(out) :: PF_D(:), PF_V1(:), PF_V2(:), PF_V3(:), &
                             PF_E(:), PF_Ne(:)
    REAL(DP), INTENT(in)  :: GF_Gm_dd_11(:), GF_Gm_dd_22(:), GF_Gm_dd_33(:)

#if defined HYDRO_NONRELATIVISTIC

    CALL Euler_ComputePrimitive_NonRelativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_ComputePrimitive_Relativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#else

    CALL Euler_ComputePrimitive_NonRelativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

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

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: AF_P(:)

#if defined HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeConserved_NonRelativistic &
           ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_ComputeConserved_Relativistic &
           ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33, &
             AF_P )

#else

    CALL Euler_ComputeConserved_NonRelativistic &
           ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#endif

  END SUBROUTINE Euler_ComputeConserved


  SUBROUTINE Euler_ComputeFromConserved &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      P(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

#if defined HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeFromConserved_NonRelativistic &
           ( iX_B0(1:3), iX_E0(1:3), iX_B1(1:3), iX_E1(1:3), &
             G(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nGF), &
             U(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nCF), &
             P(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nPF), &
             A(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nAF) )

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_ComputeFromConserved_Relativistic &
           ( iX_B0(1:3), iX_E0(1:3), iX_B1(1:3), iX_E1(1:3), &
             G(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nGF), &
             U(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nCF), &
             P(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nPF), &
             A(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nAF) )

#else

    CALL Euler_ComputeFromConserved_NonRelativistic &
           ( iX_B0(1:3), iX_E0(1:3), iX_B1(1:3), iX_E1(1:3), &
             G(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nGF), &
             U(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nCF), &
             P(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nPF), &
             A(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nAF) )


#endif

  END SUBROUTINE Euler_ComputeFromConserved


  SUBROUTINE Euler_ComputeTimeStep &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep )

    INTEGER,  INTENT(in)          :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)          :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)          :: &
      CFL
    REAL(DP), INTENT(out)         :: &
      TimeStep

    CALL TimersStart_Euler( Timer_Euler_ComputeTimeStep )

#if defined HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeTimeStep_NonRelativistic &
           ( iX_B0(1:3), iX_E0(1:3), iX_B1(1:3), iX_E1(1:3), &
             G(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nGF), &
             U(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nCF), &
             CFL, TimeStep )

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_ComputeTimeStep_Relativistic &
           ( iX_B0(1:3), iX_E0(1:3), iX_B1(1:3), iX_E1(1:3), &
             G(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nGF), &
             U(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nCF), &
             CFL, TimeStep )

#else

    CALL Euler_ComputeTimeStep_NonRelativistic &
           ( iX_B0(1:3), iX_E0(1:3), iX_B1(1:3), iX_E1(1:3), &
             G(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nGF), &
             U(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nCF), &
             CFL, TimeStep )

#endif

    CALL TimersStop_Euler( Timer_Euler_ComputeTimeStep )

  END SUBROUTINE Euler_ComputeTimeStep


  PURE FUNCTION Euler_Eigenvalues &
    ( V, Cs, Gmii, V1, V2, V3, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in) :: V, Cs, Gmii

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: V1, V2, V3, Gm11, Gm22, Gm33, Lapse, Shift

    REAL(DP) :: Euler_Eigenvalues(nCF)

#if defined HYDRO_NONRELATIVISTIC

    Euler_Eigenvalues = Euler_Eigenvalues_NonRelativistic &
                          ( V, Cs, Gmii )

#elif defined HYDRO_RELATIVISTIC

    Euler_Eigenvalues = Euler_Eigenvalues_Relativistic &
                          ( V, Cs, Gmii, V1, V2, V3, Gm11, Gm22, Gm33, &
                            Lapse, Shift )

#else

    Euler_Eigenvalues = Euler_Eigenvalues_NonRelativistic &
                          ( V, Cs, Gmii )

#endif

  END FUNCTION Euler_Eigenvalues


  REAL(DP) FUNCTION Euler_AlphaMiddle &
    ( DL, SL, EL, F_DL, F_SL, F_EL, DR, SR, ER, F_DR, F_SR, F_ER, &
      Gmii, aP, aM, Lapse, Shift )

    ! --- Gm is the covariant ii-component of the spatial three-metric
    !     Shift is the ith contravariant component of the shift-vector ---

    REAL(DP), INTENT(in) :: DL, SL, EL, F_DL, F_SL, F_EL, &
                            DR, SR, ER, F_DR, F_SR, F_ER, &
                            Gmii, aP, aM

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

#if defined HYDRO_NONRELATIVISTIC

    Euler_AlphaMiddle = Euler_AlphaMiddle_NonRelativistic &
                          ( DL, SL, EL, F_DL, F_SL, F_EL, &
                            DR, SR, ER, F_DR, F_SR, F_ER, &
                            Gmii, aP, aM )

#elif defined HYDRO_RELATIVISTIC

    Euler_AlphaMiddle = Euler_AlphaMiddle_Relativistic &
                          ( DL, SL, EL, F_DL, F_SL, F_EL, &
                            DR, SR, ER, F_DR, F_SR, F_ER, &
                            Gmii, aP, aM, Lapse, Shift )

#else

    Euler_AlphaMiddle = Euler_AlphaMiddle_NonRelativistic &
                          ( DL, SL, EL, F_DL, F_SL, F_EL, &
                            DR, SR, ER, F_DR, F_SR, F_ER, &
                            Gmii, aP, aM )

#endif

    RETURN
  END FUNCTION Euler_AlphaMiddle


  PURE FUNCTION Euler_Flux_X1 &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    ! --- Shift is the first contravariant component of the shift-vector ---

    REAL(DP)             :: Euler_Flux_X1(nCF)
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

#if defined HYDRO_NONRELATIVISTIC

    Euler_Flux_X1 = Euler_Flux_X1_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#elif defined HYDRO_RELATIVISTIC

    Euler_Flux_X1 = Euler_Flux_X1_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )

#else

    Euler_Flux_X1 = Euler_Flux_X1_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#endif

    RETURN
  END FUNCTION Euler_Flux_X1


  PURE FUNCTION Euler_Flux_X2 &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    ! --- Shift is the second contravariant component of the shift-vector ---

    REAL(DP)             :: Euler_Flux_X2(nCF)
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

#if defined HYDRO_NONRELATIVISTIC

    Euler_Flux_X2 = Euler_Flux_X2_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#elif defined HYDRO_RELATIVISTIC

    Euler_Flux_X2 = Euler_Flux_X2_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )

#else

    Euler_Flux_X2 = Euler_Flux_X2_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#endif

    RETURN
  END FUNCTION Euler_Flux_X2


  PURE FUNCTION Euler_Flux_X3 &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    ! --- Shift is the third contravariant component of the shift-vector ---

    REAL(DP)             :: Euler_Flux_X3(nCF)
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

#if defined HYDRO_NONRELATIVISTIC

    Euler_Flux_X3 = Euler_Flux_X3_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#elif defined HYDRO_RELATIVISTIC

    Euler_Flux_X3 = Euler_Flux_X3_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )

#else

    Euler_Flux_X3 = Euler_Flux_X3_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#endif

    RETURN
  END FUNCTION Euler_Flux_X3


  PURE FUNCTION Euler_StressTensor_Diagonal( S1, S2, S3, V1, V2, V3, P )

    REAL(DP), INTENT(in) :: S1, S2, S3, V1, V2, V3, P

    REAL(DP) :: Euler_StressTensor_Diagonal(3)

#if defined HYDRO_NONRELATIVISTIC

    Euler_StressTensor_Diagonal &
      = Euler_StressTensor_Diagonal_NonRelativistic( S1, S2, S3, V1, V2, V3, P )

#elif defined HYDRO_RELATIVISTIC

    Euler_StressTensor_Diagonal &
      = Euler_StressTensor_Diagonal_Relativistic( S1, S2, S3, V1, V2, V3, P )

#else

    Euler_StressTensor_Diagonal &
      = Euler_StressTensor_Diagonal_NonRelativistic( S1, S2, S3, V1, V2, V3, P )

#endif

    RETURN
  END FUNCTION Euler_StressTensor_Diagonal


  PURE FUNCTION Euler_NumericalFlux_X1 &
    ( uL, uR, fL, fR, aP, aM, aC, Gm11, vL, vR, pL, pR, Lapse, Shift )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm11

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_X1(nCF)

#if defined HYDRO_NONRELATIVISTIC

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X1 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X1 = Euler_NumericalFlux_X1_HLLC_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm11 )

#else

    Euler_NumericalFlux_X1 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#elif defined HYDRO_RELATIVISTIC

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X1 = Euler_NumericalFlux_HLL_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X1 = Euler_NumericalFlux_X1_HLLC_Relativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm11, &
                                 vL, vR, pL, pR, Lapse, Shift )

#else

    Euler_NumericalFlux_X1 = Euler_NumericalFlux_HLL_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#else

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X1 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X1 = Euler_NumericalFlux_X1_HLLC_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm11 )

#else

    Euler_NumericalFlux_X1 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_X1


  PURE FUNCTION Euler_NumericalFlux_X2 &
    ( uL, uR, fL, fR, aP, aM, aC, Gm22, vL, vR, pL, pR, Lapse, Shift )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm22

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_X2(nCF)

#if defined HYDRO_NONRELATIVISTIC

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X2 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X2 = Euler_NumericalFlux_X2_HLLC_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm22 )

#else

    Euler_NumericalFlux_X2 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#elif defined HYDRO_RELATIVISTIC

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X2 = Euler_NumericalFlux_HLL_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X2 = Euler_NumericalFlux_X2_HLLC_Relativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm22, &
                                 vL, vR, pL, pR, Lapse, Shift )

#else

    Euler_NumericalFlux_X2 = Euler_NumericalFlux_HLL_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#else

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X2 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X2 = Euler_NumericalFlux_X2_HLLC_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm22 )

#else

    Euler_NumericalFlux_X2 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_X2


  PURE FUNCTION Euler_NumericalFlux_X3 &
    ( uL, uR, fL, fR, aP, aM, aC, Gm33, vL, vR, pL, pR, Lapse, Shift )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm33

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_X3(nCF)

#if defined HYDRO_NONRELATIVISTIC

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X3 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X3 = Euler_NumericalFlux_X3_HLLC_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm33 )

#else

    Euler_NumericalFlux_X3 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#elif defined HYDRO_RELATIVISTIC

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X3 = Euler_NumericalFlux_HLL_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X3 = Euler_NumericalFlux_X3_HLLC_Relativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm33, &
                                 vL, vR, pL, pR, Lapse, Shift )

#else

    Euler_NumericalFlux_X3 = Euler_NumericalFlux_HLL_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#else

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X3 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X3 = Euler_NumericalFlux_X3_HLLC_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm33 )

#else

    Euler_NumericalFlux_X3 = Euler_NumericalFlux_HLL_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_X3


END MODULE Euler_UtilitiesModule
