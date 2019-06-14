MODULE Euler_UtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE FluidFieldsModule, ONLY: &
    nCF







  USE Euler_UtilitiesModule_Relativistic



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










    CALL Euler_ComputePrimitive_Relativistic &
           ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )



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



    CALL Euler_ComputeConserved_Relativistic &
           ( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
             CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33, &
             AF_P )



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







    CALL Euler_ComputeFromConserved_Relativistic( iX_B, iX_E, G, U, P, A )



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

    LOGICAL :: UseSourceTerm

    UseSourceTerm = .FALSE.
    IF( PRESENT( UseSourceTerm_Option ) ) &
      UseSourceTerm = UseSourceTerm_Option








    CALL Euler_ComputeTimeStep_Relativistic &
           ( iX_B, iX_E, G, U, CFL, TimeStep, UseSourceTerm )



  END SUBROUTINE Euler_ComputeTimeStep


  PURE FUNCTION Euler_Eigenvalues &
    ( V, Cs, V1, V2, V3, Gm, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in) :: V, Cs

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: V1, V2, V3, Gm, Gm11, Gm22, Gm33, &
                            Lapse, Shift

    REAL(DP) :: Euler_Eigenvalues(nCF)









    Euler_Eigenvalues = Euler_Eigenvalues_Relativistic &
                          ( V, Cs, V1, V2, V3, Gm, Gm11, Gm22, Gm33, &
                            Lapse, Shift )



  END FUNCTION Euler_Eigenvalues


  REAL(DP) FUNCTION Euler_AlphaMiddle &
    ( DL, SL, EL, F_DL, F_SL, F_EL, DR, SR, ER, F_DR, F_SR, F_ER, &
      Gm, Lapse, Shift, aP, aM )

    ! --- Gm is the covariant ii-component of the spatial three-metric
    !     Shift is the ith contravariant component of the shift-vector ---

    REAL(DP), INTENT(in) :: DL, SL, EL, F_DL, F_SL, F_EL, &
                            DR, SR, ER, F_DR, F_SR, F_ER, &
                            Gm, aP, aM

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift










    Euler_AlphaMiddle = Euler_AlphaMiddle_Relativistic &
                          ( DL, SL, EL, F_DL, F_SL, F_EL, &
                            DR, SR, ER, F_DR, F_SR, F_ER, &
                            Gm, Lapse, Shift, aP, aM )



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









    Euler_Flux_X1 = Euler_Flux_X1_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )



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









    Euler_Flux_X2 = Euler_Flux_X2_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )



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









    Euler_Flux_X3 = Euler_Flux_X3_Relativistic &
                      ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift )



    RETURN
  END FUNCTION Euler_Flux_X3


  PURE FUNCTION Euler_StressTensor_Diagonal( S1, S2, S3, V1, V2, V3, P )

    REAL(DP), INTENT(in) :: S1, S2, S3, V1, V2, V3, P

    REAL(DP) :: Euler_StressTensor_Diagonal(3)








    Euler_StressTensor_Diagonal &
      = Euler_StressTensor_Diagonal_Relativistic( S1, S2, S3, V1, V2, V3, P )



    RETURN
  END FUNCTION Euler_StressTensor_Diagonal


  PURE FUNCTION Euler_NumericalFlux_X1 &
    ( uL, uR, fL, fR, Gm, vL, vR, pL, pR, Lapse, Shift, aP, aM, aC )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            Gm, aP, aM, aC

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_X1(nCF)





    Euler_NumericalFlux_X1 = Euler_NumericalFlux_HLL_Relativistic &
                               ( uL, uR, fL, fR, Gm, vL, vR, pL, pR, &
                                 Lapse, Shift, aP, aM, aC )













    RETURN
  END FUNCTION Euler_NumericalFlux_X1


  PURE FUNCTION Euler_NumericalFlux_X2 &
    ( uL, uR, fL, fR, Gm, vL, vR, pL, pR, Lapse, Shift, aP, aM, aC )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            Gm, aP, aM, aC

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_X2(nCF)





    Euler_NumericalFlux_X2 = Euler_NumericalFlux_HLL_Relativistic &
                               ( uL, uR, fL, fR, Gm, vL, vR, pL, pR, &
                                 Lapse, Shift, aP, aM, aC )













    RETURN
  END FUNCTION Euler_NumericalFlux_X2


  PURE FUNCTION Euler_NumericalFlux_X3 &
    ( uL, uR, fL, fR, Gm, vL, vR, pL, pR, Lapse, Shift, aP, aM, aC )

    REAL(DP), INTENT(in) :: uL(nCF), uR(nCF), fL(nCF), fR(nCF), &
                            Gm, aP, aM, aC

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: vL, vR, pL, pR, Lapse, Shift

    REAL(DP) :: Euler_NumericalFlux_X3(nCF)





    Euler_NumericalFlux_X3 = Euler_NumericalFlux_HLL_Relativistic &
                               ( uL, uR, fL, fR, Gm, vL, vR, pL, pR, &
                                 Lapse, Shift, aP, aM, aC )













    RETURN
  END FUNCTION Euler_NumericalFlux_X3


END MODULE Euler_UtilitiesModule

