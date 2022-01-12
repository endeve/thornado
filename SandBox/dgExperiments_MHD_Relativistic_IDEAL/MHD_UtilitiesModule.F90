MODULE MHD_UtilitiesModule

  USE KindModule, ONLY: &
    DP, &
    Half
  USE MagnetofluidFieldsModule, ONLY: &
    nCM
  USE MHD_UtilitiesModule_Relativistic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_MHD
  PUBLIC :: ComputeConserved_MHD
  PUBLIC :: ComputeFromConserved_MHD
  PUBLIC :: ComputeTimeStep_MHD
  PUBLIC :: Eigenvalues_MHD
  PUBLIC :: Flux_X1_MHD
  PUBLIC :: NumericalFlux_MHD_X1

  INTERFACE ComputePrimitive_MHD
    MODULE PROCEDURE ComputePrimitive_Scalar
    MODULE PROCEDURE ComputePrimitive_Vector
  END INTERFACE ComputePrimitive_MHD

  INTERFACE ComputeConserved_MHD
    MODULE PROCEDURE ComputeConserved_Scalar
    MODULE PROCEDURE ComputeConserved_Vector
  END INTERFACE ComputeConserved_MHD


CONTAINS


  SUBROUTINE ComputePrimitive_Scalar &
    ( CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi,            &
      PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi,            &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
      GF_Alpha, GF_Beta_1, GF_Beta_2, GF_Beta_3 )

    REAL(DP), INTENT(in)  :: &
      CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi
    REAL(DP), INTENT(out) :: &
      PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi
    REAL(DP), INTENT(in)  :: &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33, &
      GF_Alpha, GF_Beta_1, GF_Beta_2, GF_Beta_3

    CALL ComputePrimitive_MHD_Relativistic &
           ( CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
             CM_B1, CM_B2, CM_B3, CM_Chi,            &
             PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
             PM_B1, PM_B2, PM_B3, PM_Chi,            &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
             GF_Alpha, GF_Beta_1, GF_Beta_2, GF_Beta_3 )

  END SUBROUTINE ComputePrimitive_Scalar


  SUBROUTINE ComputePrimitive_Vector &
    ( CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi,            &
      PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi,            &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
      GF_Alpha, GF_Beta_1, GF_Beta_2, GF_Beta_3 )
 
    REAL(DP), INTENT(in)  :: &
      CM_D(:), CM_S1(:), CM_S2(:), CM_S3(:), CM_E(:), CM_Ne(:), &
      CM_B1(:), CM_B2(:), CM_B3(:), CM_Chi(:)
    REAL(DP), INTENT(out) :: &
      PM_D(:), PM_V1(:), PM_V2(:), PM_V3(:), PM_E(:), PM_Ne(:), &
      PM_B1(:), PM_B2(:), PM_B3(:), PM_Chi(:)
    REAL(DP), INTENT(in)  :: &
      GF_Gm_dd_11(:), GF_Gm_dd_22(:), GF_Gm_dd_33(:), &
      GF_Alpha(:), GF_Beta_1(:), GF_Beta_2(:), GF_Beta_3(:)

    CALL ComputePrimitive_MHD_Relativistic &
           ( CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
             CM_B1, CM_B2, CM_B3, CM_Chi,            &
             PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
             PM_B1, PM_B2, PM_B3, PM_Chi,            &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
             GF_Alpha, GF_Beta_1, GF_Beta_2, GF_Beta_3 )

  END SUBROUTINE ComputePrimitive_Vector


  SUBROUTINE ComputeConserved_Scalar &
    ( PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi,            &
      CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi,            &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
      GF_Alpha, GF_Beta_1, GF_Beta_2, GF_Beta_3, &
      AM_P )

    REAL(DP), INTENT(in)  :: &
      PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi
    REAL(DP), INTENT(out) :: &
      CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi
    REAL(DP), INTENT(in)  :: &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33, &
      GF_Alpha, GF_Beta_1, GF_Beta_2, GF_Beta_3

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: AM_P

    CALL ComputeConserved_MHD_Relativistic &
           ( PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
             PM_B1, PM_B2, PM_B3, PM_Chi,            &
             CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
             CM_B1, CM_B2, CM_B3, CM_Chi,            &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
             GF_Alpha, GF_Beta_1, GF_Beta_2, GF_Beta_3, &
             AM_P )

  END SUBROUTINE ComputeConserved_Scalar


  SUBROUTINE ComputeConserved_Vector &
    ( PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
      PM_B1, PM_B2, PM_B3, PM_Chi,            &
      CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
      CM_B1, CM_B2, CM_B3, CM_Chi,            &
      GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
      GF_Alpha, GF_Beta_1, GF_Beta_2, GF_Beta_3, &
      AM_P )

    REAL(DP), INTENT(in)  :: &
      PM_D(:), PM_V1(:), PM_V2(:), PM_V3(:), PM_E(:), PM_Ne(:), &
      PM_B1(:), PM_B2(:), PM_B3(:), PM_Chi(:)
    REAL(DP), INTENT(out) :: &
      CM_D(:), CM_S1(:), CM_S2(:), CM_S3(:), CM_E(:), CM_Ne(:), &
      CM_B1(:), CM_B2(:), CM_B3(:), CM_Chi(:)
    REAL(DP), INTENT(in)  :: &
      GF_Gm_dd_11(:), GF_Gm_dd_22(:), GF_Gm_dd_33(:), &
      GF_Alpha(:), GF_Beta_1(:), GF_Beta_2(:), GF_Beta_3(:)

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: AM_P(:)

    CALL ComputeConserved_MHD_Relativistic &
           ( PM_D, PM_V1, PM_V2, PM_V3, PM_E, PM_Ne, &
             PM_B1, PM_B2, PM_B3, PM_Chi,            &
             CM_D, CM_S1, CM_S2, CM_S3, CM_E, CM_Ne, &
             CM_B1, CM_B2, CM_B3, CM_Chi,            &
             GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
             GF_Alpha, GF_Beta_1, GF_Beta_2, GF_Beta_3, &
             AM_P )

  END SUBROUTINE ComputeConserved_Vector


  SUBROUTINE ComputeFromConserved_MHD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      P(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CALL ComputeFromConserved_MHD_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

  END SUBROUTINE ComputeFromConserved_MHD


  SUBROUTINE ComputeTimeStep_MHD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CML, TimeStep )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)  :: &
      CML
    REAL(DP), INTENT(out) :: &
      TimeStep

    CALL ComputeTimeStep_MHD_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CML, TimeStep )

  END SUBROUTINE ComputeTimeStep_MHD


  FUNCTION Eigenvalues_MHD &
    ( Vi, Cs, Gmii, D, V1, V2, V3, E, Ne, &
      B1, B2, B3, Chi, &
      Gm11, Gm22, Gm33, &
      Lapse, Shift_X1, Shift_X2, Shift_X3 )

    REAL(DP), INTENT(in) :: Vi, Cs, Gmii

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, &
                            B1, B2, B3, Chi, &
                            Gm11, Gm22, Gm33, Lapse, & 
                            Shift_X1, Shift_X2, Shift_X3

    REAL(DP) :: Eigenvalues_MHD(2)

    Eigenvalues_MHD = Eigenvalues_MHD_Relativistic &
                          ( Vi, Cs, Gmii, D, V1, V2, V3, E, Ne, &
                            B1, B2, B3, Chi, &
                            Gm11, Gm22, Gm33, &
                            Lapse, Shift_X1, Shift_X2, Shift_X3 )

    RETURN
  END FUNCTION Eigenvalues_MHD


  FUNCTION Flux_X1_MHD &
    ( D, V1, V2, V3, E, Ne, B1, B2, B3, Chi, P, Gm11, Gm22, Gm33, Lapse, & 
      Shift_X1, Shift_X2, Shift_X3 )

    ! --- Shift is the first contravariant component of the shift-vector ---

    REAL(DP)             :: Flux_X1_MHD(nCM)
    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, &
                            B1, B2, B3, Chi, P
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    ! --- Only needed for relativistic code ---
    REAL(DP), INTENT(in) :: Lapse, Shift_X1, Shift_X2, Shift_X3

    Flux_X1_MHD = Flux_X1_MHD_Relativistic &
                      ( D, V1, V2, V3, E, Ne, B1, B2, B3, Chi, P, Gm11, Gm22, Gm33, &
                        Lapse, Shift_X1, Shift_X2, Shift_X3 )

    RETURN
  END FUNCTION Flux_X1_MHD


  FUNCTION NumericalFlux_MHD_X1 &
    ( uL, uR, fL, fR, aP, aM )

    REAL(DP), INTENT(in) :: uL(nCM), uR(nCM), fL(nCM), fR(nCM), &
                            aP, aM

    REAL(DP) :: NumericalFlux_MHD_X1(nCM)

    NumericalFlux_MHD_X1 &
      = NumericalFlux_MHD_HLL &
          ( uL, uR, fL, fR, aP, aM )

    RETURN
  END FUNCTION NumericalFlux_MHD_X1


  FUNCTION NumericalFlux_MHD_HLL( uL, uR, fL, fR, aP, aM )

    REAL(DP), INTENT(in) :: uL(nCM), uR(nCM), fL(nCM), fR(nCM), aP, aM

    REAL(DP) :: NumericalFlux_MHD_HLL(nCM)

    NumericalFlux_MHD_HLL = NumericalFlux_HLL_MHD_Relativistic &
                                ( uL, uR, fL, fR, aP, aM )

    RETURN
  END FUNCTION NumericalFlux_MHD_HLL


END MODULE MHD_UtilitiesModule
