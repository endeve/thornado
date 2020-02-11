MODULE Euler_UtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule, ONLY: &
    nCF, nPF, nAF

#if defined HYDRO_NONRELATIVISTIC

  USE Euler_UtilitiesModule_NonRelativistic

#elif defined HYDRO_RELATIVISTIC

  USE Euler_UtilitiesModule_Relativistic

#else

  USE Euler_UtilitiesModule_NonRelativistic

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_Euler
  PUBLIC :: ComputeConserved_Euler
  PUBLIC :: ComputeFromConserved_Euler
  PUBLIC :: ComputeTimeStep_Euler
  PUBLIC :: Eigenvalues_Euler
  PUBLIC :: AlphaMiddle_Euler
  PUBLIC :: Flux_X1_Euler
  PUBLIC :: Flux_X2_Euler
  PUBLIC :: Flux_X3_Euler
  PUBLIC :: StressTensor_Diagonal_Euler
  PUBLIC :: Euler_NumericalFlux_X1
  PUBLIC :: Euler_NumericalFlux_X2
  PUBLIC :: Euler_NumericalFlux_X3

  INTERFACE ComputePrimitive_Euler
    MODULE PROCEDURE ComputePrimitive_Scalar
    MODULE PROCEDURE ComputePrimitive_Vector
  END INTERFACE ComputePrimitive_Euler

CONTAINS


  SUBROUTINE ComputePrimitive_Scalar &
    ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

    REAL(DP), INTENT(in)            :: &
      CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne
    REAL(DP), INTENT(out)           :: &
      PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne
    REAL(DP), INTENT(in)            :: &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33

#if defined HYDRO_NONRELATIVISTIC

    CALL ComputePrimitive_Euler_NonRelativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#elif defined HYDRO_RELATIVISTIC

    CALL ComputePrimitive_Euler_Relativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#else

    CALL ComputePrimitive_Euler_NonRelativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#endif

  END SUBROUTINE ComputePrimitive_Scalar


  SUBROUTINE ComputePrimitive_Vector &
    ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
      PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

    REAL(DP), INTENT(in)            :: &
      CF_D(:), CF_S1(:), CF_S2(:), CF_S3(:), CF_E(:), CF_Ne(:)
    REAL(DP), INTENT(out)           :: &
      PF_D(:), PF_V1(:), PF_V2(:), PF_V3(:), PF_E(:), PF_Ne(:)
    REAL(DP), INTENT(in)            :: &
      GF_Gm_dd_11(:), GF_Gm_dd_22(:), GF_Gm_dd_33(:)

#if defined HYDRO_NONRELATIVISTIC

    CALL ComputePrimitive_Euler_NonRelativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#elif defined HYDRO_RELATIVISTIC

    CALL ComputePrimitive_Euler_Relativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#else

    CALL ComputePrimitive_Euler_NonRelativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#endif

  END SUBROUTINE ComputePrimitive_Vector


  SUBROUTINE ComputeConserved_Euler &
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

    CALL ComputeConserved_Euler_NonRelativistic &
           ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#elif defined HYDRO_RELATIVISTIC

    CALL ComputeConserved_Euler_Relativistic &
           ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33, &
             AF_P )

#else

    CALL ComputeConserved_Euler_NonRelativistic &
           ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

#endif

  END SUBROUTINE ComputeConserved_Euler


  SUBROUTINE ComputeFromConserved_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

    INTEGER,  INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)            :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out)           :: &
      P(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

#if defined HYDRO_NONRELATIVISTIC

    CALL ComputeFromConserved_Euler_NonRelativistic &
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

    CALL ComputeFromConserved_Euler_Relativistic &
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

    CALL ComputeFromConserved_Euler_NonRelativistic &
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

  END SUBROUTINE ComputeFromConserved_Euler


  SUBROUTINE ComputeTimeStep_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep )

    INTEGER,  INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)            :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)            :: &
      CFL
    REAL(DP), INTENT(out)           :: &
      TimeStep

#if defined HYDRO_NONRELATIVISTIC

    CALL ComputeTimeStep_Euler_NonRelativistic &
           ( iX_B0(1:3), iX_E0(1:3), iX_B1(1:3), iX_E1(1:3), &
             G(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nGF), &
             U(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nCF), &
             CFL, TimeStep )

#elif defined HYDRO_RELATIVISTIC

    CALL ComputeTimeStep_Euler_Relativistic &
           ( iX_B0(1:3), iX_E0(1:3), iX_B1(1:3), iX_E1(1:3), &
             G(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nGF), &
             U(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nCF), &
             CFL, TimeStep )

#else

    CALL ComputeTimeStep_Euler_NonRelativistic &
           ( iX_B0(1:3), iX_E0(1:3), iX_B1(1:3), iX_E1(1:3), &
             G(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nGF), &
             U(1:nDOFX,iX_B1(1):iX_E1(1),&
                       iX_B1(2):iX_E1(2),&
                       iX_B1(3):iX_E1(3),1:nCF), &
             CFL, TimeStep )

#endif

  END SUBROUTINE ComputeTimeStep_Euler


  PURE FUNCTION Eigenvalues_Euler &
    ( V, Cs, Gmii, V1, V2, V3, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in) :: V, Cs, Gmii

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: V1, V2, V3, Gm11, Gm22, Gm33, Lapse, Shift

    REAL(DP) :: Eigenvalues_Euler(nCF)

#if defined HYDRO_NONRELATIVISTIC

    Eigenvalues_Euler = Eigenvalues_Euler_NonRelativistic &
                          ( V, Cs, Gmii )

#elif defined HYDRO_RELATIVISTIC

    Eigenvalues_Euler = Eigenvalues_Euler_Relativistic &
                          ( V, Cs, Gmii, V1, V2, V3, Gm11, Gm22, Gm33, &
                            Lapse, Shift )

#else

    Eigenvalues_Euler = Eigenvalues_Euler_NonRelativistic &
                          ( V, Cs, Gmii )

#endif

  END FUNCTION Eigenvalues_Euler


  REAL(DP) FUNCTION AlphaMiddle_Euler &
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

    AlphaMiddle_Euler = AlphaMiddle_Euler_NonRelativistic &
                          ( DL, SL, EL, F_DL, F_SL, F_EL, &
                            DR, SR, ER, F_DR, F_SR, F_ER, &
                            Gmii, aP, aM )

#elif defined HYDRO_RELATIVISTIC

    AlphaMiddle_Euler = AlphaMiddle_Euler_Relativistic &
                          ( DL, SL, EL, F_DL, F_SL, F_EL, &
                            DR, SR, ER, F_DR, F_SR, F_ER, &
                            Gmii, aP, aM, Lapse, Shift )

#else

    AlphaMiddle_Euler = AlphaMiddle_Euler_NonRelativistic &
                          ( DL, SL, EL, F_DL, F_SL, F_EL, &
                            DR, SR, ER, F_DR, F_SR, F_ER, &
                            Gmii, aP, aM )

#endif

    RETURN
  END FUNCTION AlphaMiddle_Euler


  PURE FUNCTION Flux_X1_Euler &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    ! --- Shift is the first contravariant component of the shift-vector ---

    REAL(DP)             :: Flux_X1_Euler(nCF)
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

#if defined HYDRO_NONRELATIVISTIC

    Flux_X1_Euler = Flux_X1_Euler_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#elif defined HYDRO_RELATIVISTIC

    Flux_X1_Euler = Flux_X1_Euler_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )

#else

    Flux_X1_Euler = Flux_X1_Euler_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#endif

    RETURN
  END FUNCTION Flux_X1_Euler


  PURE FUNCTION Flux_X2_Euler &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    ! --- Shift is the second contravariant component of the shift-vector ---

    REAL(DP)             :: Flux_X2_Euler(nCF)
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

#if defined HYDRO_NONRELATIVISTIC

    Flux_X2_Euler = Flux_X2_Euler_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#elif defined HYDRO_RELATIVISTIC

    Flux_X2_Euler = Flux_X2_Euler_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )

#else

    Flux_X2_Euler = Flux_X2_Euler_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#endif

    RETURN
  END FUNCTION Flux_X2_Euler


  PURE FUNCTION Flux_X3_Euler &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    ! --- Shift is the third contravariant component of the shift-vector ---

    REAL(DP)             :: Flux_X3_Euler(nCF)
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift

#if defined HYDRO_NONRELATIVISTIC

    Flux_X3_Euler = Flux_X3_Euler_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#elif defined HYDRO_RELATIVISTIC

    Flux_X3_Euler = Flux_X3_Euler_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )

#else

    Flux_X3_Euler = Flux_X3_Euler_NonRelativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33 )

#endif

    RETURN
  END FUNCTION Flux_X3_Euler


  PURE FUNCTION StressTensor_Diagonal_Euler( S1, S2, S3, V1, V2, V3, P )

    REAL(DP), INTENT(in) :: S1, S2, S3, V1, V2, V3, P

    REAL(DP) :: StressTensor_Diagonal_Euler(3)

#if defined HYDRO_NONRELATIVISTIC

    StressTensor_Diagonal_Euler &
      = StressTensor_Diagonal_Euler_NonRelativistic( S1, S2, S3, V1, V2, V3, P )

#elif defined HYDRO_RELATIVISTIC

    StressTensor_Diagonal_Euler &
      = StressTensor_Diagonal_Euler_Relativistic( S1, S2, S3, V1, V2, V3, P )

#else

    StressTensor_Diagonal_Euler &
      = StressTensor_Diagonal_Euler_NonRelativistic( S1, S2, S3, V1, V2, V3, P )

#endif

    RETURN
  END FUNCTION StressTensor_Diagonal_Euler


  PURE FUNCTION Euler_NumericalFlux_X1 &
    ( uL, uR, fL, fR, aP, aM, aC, Gm11, vL, vR, pL, pR, Lapse, Shift )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            aP, aM, aC, Gm11

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_X1(nCF)

#if defined HYDRO_NONRELATIVISTIC

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X1 = NumericalFlux_HLL_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X1 = NumericalFlux_X1_HLLC_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm11 )

#else

    Euler_NumericalFlux_X1 = NumericalFlux_HLL_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#elif defined HYDRO_RELATIVISTIC

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X1 = NumericalFlux_HLL_Euler_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X1 = NumericalFlux_X1_HLLC_Euler_Relativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm11, &
                                 vL, vR, pL, pR, Lapse, Shift )

#else

    Euler_NumericalFlux_X1 = NumericalFlux_HLL_Euler_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#else

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X1 = NumericalFlux_HLL_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X1 = NumericalFlux_X1_HLLC_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm11 )

#else

    Euler_NumericalFlux_X1 = NumericalFlux_HLL_Euler_NonRelativistic &
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

    Euler_NumericalFlux_X2 = NumericalFlux_HLL_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X2 = NumericalFlux_X2_HLLC_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm22 )

#else

    Euler_NumericalFlux_X2 = NumericalFlux_HLL_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#elif defined HYDRO_RELATIVISTIC

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X2 = NumericalFlux_HLL_Euler_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X2 = NumericalFlux_X2_HLLC_Euler_Relativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm22, &
                                 vL, vR, pL, pR, Lapse, Shift )

#else

    Euler_NumericalFlux_X2 = NumericalFlux_HLL_Euler_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#else

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X2 = NumericalFlux_HLL_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X2 = NumericalFlux_X2_HLLC_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm22 )

#else

    Euler_NumericalFlux_X2 = NumericalFlux_HLL_Euler_NonRelativistic &
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

    Euler_NumericalFlux_X3 = NumericalFlux_HLL_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X3 = NumericalFlux_X3_HLLC_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm33 )

#else

    Euler_NumericalFlux_X3 = NumericalFlux_HLL_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#elif defined HYDRO_RELATIVISTIC

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X3 = NumericalFlux_HLL_Euler_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X3 = NumericalFlux_X3_HLLC_Euler_Relativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm33, &
                                 vL, vR, pL, pR, Lapse, Shift )

#else

    Euler_NumericalFlux_X3 = NumericalFlux_HLL_Euler_Relativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#else

#if defined HYDRO_RIEMANN_SOLVER_HLL

    Euler_NumericalFlux_X3 = NumericalFlux_HLL_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#elif defined HYDRO_RIEMANN_SOLVER_HLLC

    Euler_NumericalFlux_X3 = NumericalFlux_X3_HLLC_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM, aC, Gm33 )

#else

    Euler_NumericalFlux_X3 = NumericalFlux_HLL_Euler_NonRelativistic &
                               ( uL, uR, fL, fR, aP, aM )

#endif

#endif

    RETURN
  END FUNCTION Euler_NumericalFlux_X3


END MODULE Euler_UtilitiesModule
