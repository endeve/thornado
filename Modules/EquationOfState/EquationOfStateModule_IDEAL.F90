MODULE EquationOfStateModule_IDEAL

  USE KindModule, ONLY: &
    DP, Zero, One, Two

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PUBLIC :: &
    Gamma_IDEAL

  PUBLIC :: InitializeEquationOfState_IDEAL
  PUBLIC :: FinalizeEquationOfState_IDEAL
  PUBLIC :: ComputeInternalEnergyDensityFromPressure_IDEAL
  PUBLIC :: ComputePressureFromPrimitive_IDEAL
  PUBLIC :: ComputeMagneticPressureFromPrimitive_IDEAL
  PUBLIC :: ComputePressureFromSpecificInternalEnergy_IDEAL
  PUBLIC :: ComputeEnthalpyFromPrimitive_IDEAL
  PUBLIC :: ComputeMagneticEnthalpyFromPrimitive_IDEAL
  PUBLIC :: ComputeSoundSpeedFromPrimitive_IDEAL
  PUBLIC :: ComputeAlfvenSpeedFromPrimitive_IDEAL
  PUBLIC :: ComputeAuxiliary_Fluid_IDEAL
  PUBLIC :: ComputeAuxiliary_Magnetofluid_IDEAL

  INTERFACE ComputePressureFromPrimitive_IDEAL
    MODULE PROCEDURE ComputePressureFromPrimitive_IDEAL_Scalar
    MODULE PROCEDURE ComputePressureFromPrimitive_IDEAL_Vector
  END INTERFACE ComputePressureFromPrimitive_IDEAL

  INTERFACE ComputeMagneticPressureFromPrimitive_IDEAL
    MODULE PROCEDURE ComputeMagneticPressureFromPrimitive_IDEAL_Scalar
    MODULE PROCEDURE ComputeMagneticPressureFromPrimitive_IDEAL_Vector
  END INTERFACE ComputeMagneticPressureFromPrimitive_IDEAL

  INTERFACE ComputePressureFromSpecificInternalEnergy_IDEAL
    MODULE PROCEDURE ComputePressureFromSpecificInternalEnergy_IDEAL_Scalar
    MODULE PROCEDURE ComputePressureFromSpecificInternalEnergy_IDEAL_Vector
  END INTERFACE ComputePressureFromSpecificInternalEnergy_IDEAL

  INTERFACE ComputeInternalEnergyDensityFromPressure_IDEAL
    MODULE PROCEDURE ComputeInternalEnergyDensityFromPressure_IDEAL_Scalar
    MODULE PROCEDURE ComputeInternalEnergyDensityFromPressure_IDEAL_Vector
  END INTERFACE ComputeInternalEnergyDensityFromPressure_IDEAL

  INTERFACE ComputeEnthalpyFromPrimitive_IDEAL
    MODULE PROCEDURE ComputeEnthalpyFromPrimitive_IDEAL_Scalar
    MODULE PROCEDURE ComputeEnthalpyFromPrimitive_IDEAL_Vector
  END INTERFACE ComputeEnthalpyFromPrimitive_IDEAL

  INTERFACE ComputeMagneticEnthalpyFromPrimitive_IDEAL
    MODULE PROCEDURE ComputeMagneticEnthalpyFromPrimitive_IDEAL_Scalar
    MODULE PROCEDURE ComputeMagneticEnthalpyFromPrimitive_IDEAL_Vector
  END INTERFACE ComputeMagneticEnthalpyFromPrimitive_IDEAL

  INTERFACE ComputeSoundSpeedFromPrimitive_IDEAL
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_IDEAL_Scalar
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_IDEAL_Vector
  END INTERFACE ComputeSoundSpeedFromPrimitive_IDEAL

  INTERFACE ComputeAlfvenSpeedFromPrimitive_IDEAL
    MODULE PROCEDURE ComputeAlfvenSpeedFromPrimitive_IDEAL_Scalar
    MODULE PROCEDURE ComputeAlfvenSpeedFromPrimitive_IDEAL_Vector
  END INTERFACE ComputeAlfvenSpeedFromPrimitive_IDEAL

  INTERFACE ComputeAuxiliary_Fluid_IDEAL
    MODULE PROCEDURE ComputeAuxiliary_Fluid_IDEAL_Scalar
    MODULE PROCEDURE ComputeAuxiliary_Fluid_IDEAL_Vector
  END INTERFACE ComputeAuxiliary_Fluid_IDEAL

  INTERFACE ComputeAuxiliary_Magnetofluid_IDEAL
    MODULE PROCEDURE ComputeAuxiliary_Magnetofluid_IDEAL_Scalar
    MODULE PROCEDURE ComputeAuxiliary_Magnetofluid_IDEAL_Vector
  END INTERFACE ComputeAuxiliary_Magnetofluid_IDEAL


#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET( Gamma_IDEAL )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE( Gamma_IDEAL )
#endif

CONTAINS


  SUBROUTINE InitializeEquationOfState_IDEAL &
    ( Gamma_IDEAL_Option, Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Gamma_IDEAL_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    IF( PRESENT( Gamma_IDEAL_Option ) )THEN
      Gamma_IDEAL = Gamma_IDEAL_Option
    ELSE
      Gamma_IDEAL = 5.0_DP / 3.0_DP
    END IF

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A7,A13,ES10.3E3)') &
           '', 'Gamma_IDEAL: ', Gamma_IDEAL
      WRITE(*,*)

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( Gamma_IDEAL )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( Gamma_IDEAL )
#endif

  END SUBROUTINE InitializeEquationOfState_IDEAL


  SUBROUTINE FinalizeEquationOfState_IDEAL

  END SUBROUTINE FinalizeEquationOfState_IDEAL


  SUBROUTINE ComputeInternalEnergyDensityFromPressure_IDEAL_Scalar &
    ( D, P, Y, Ev )

    REAL(DP), INTENT(in)  :: D, P, Y
    REAL(DP), INTENT(out) :: Ev

    Ev = P / ( Gamma_IDEAL - 1.0_DP )

  END SUBROUTINE ComputeInternalEnergyDensityFromPressure_IDEAL_Scalar


  SUBROUTINE ComputeInternalEnergyDensityFromPressure_IDEAL_Vector &
    ( D, P, Y, Ev )

    REAL(DP), INTENT(in)  :: D(:), P(:), Y(:)
    REAL(DP), INTENT(out) :: Ev(:)

    Ev = P / ( Gamma_IDEAL - 1.0_DP )

  END SUBROUTINE ComputeInternalEnergyDensityFromPressure_IDEAL_Vector

#ifdef HYDRO_RELATIVISTIC

  SUBROUTINE ComputeEnthalpyFromPrimitive_IDEAL_Scalar &
    ( D, Ev, Ne, h )

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: h

    REAL(DP) :: P

    CALL ComputePressureFromPrimitive_IDEAL( D, Ev, Ne, P)

    h = One + (Ev + P) / D

  END SUBROUTINE ComputeEnthalpyFromPrimitive_IDEAL_Scalar


  SUBROUTINE ComputeEnthalpyFromPrimitive_IDEAL_Vector &
    ( D, Ev, Ne, h )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: h(:)

    INTEGER :: i

    DO i = 1, SIZE( D )

      CALL ComputeEnthalpyFromPrimitive_IDEAL_Scalar &
             ( D(i), Ev(i), Ne(i), h(i) )

    END DO

  END SUBROUTINE ComputeEnthalpyFromPrimitive_IDEAL_Vector

#else

  SUBROUTINE ComputeEnthalpyFromPrimitive_IDEAL_Scalar &
    ( D, Ev, Ne, h )

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: h

    REAL(DP) :: P

    CALL ComputePressureFromPrimitive_IDEAL( D, Ev, Ne, P)

    h = (Ev + P) / D

  END SUBROUTINE ComputeEnthalpyFromPrimitive_IDEAL_Scalar


  SUBROUTINE ComputeEnthalpyFromPrimitive_IDEAL_Vector &
    ( D, Ev, Ne, h )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: h(:)

    INTEGER :: i

    DO i = 1, SIZE( D )

      CALL ComputeEnthalpyFromPrimitive_IDEAL_Scalar &
             ( D(i), Ev(i), Ne(i), h(i) )

    END DO

  END SUBROUTINE ComputeEnthalpyFromPrimitive_IDEAL_Vector

#endif

  SUBROUTINE ComputeMagneticEnthalpyFromPrimitive_IDEAL_Scalar &
    ( D, V1, V2, V3,    &
      B1, B2, B3,       &
      Gm11, Gm22, Gm33, &
      Lapse, Shift1, Shift2, Shift3, hb )

    REAL(DP), INTENT(in)  :: D, V1, V2, V3,    &
                             B1, B2, B3,       &
                             Gm11, Gm22, Gm33, &
                             Lapse, Shift1, Shift2, Shift3
    REAL(DP), INTENT(out) :: hb

    REAL(DP) :: vdotb, vdotBeta, bdotBeta, bdotb, BetaSq, b0u, b0d, bSq

    vdotb    = Gm11 * V1 * B1 + Gm22 * V2 * B2 + Gm33 * V3 * B3

    bdotb    = Gm11 * B1**2 + Gm22 * B2**2 + Gm33 * B3**2

    vdotBeta = Gm11 * V1 * Shift1 + Gm22 * V2 * Shift2 + Gm33 * V3 * Shift3

    bdotBeta = Gm11 * B1 * Shift1 + Gm22 * B2 * Shift2 + Gm33 * B3 * Shift3

    BetaSq   = Gm11 * Shift1**2 + Gm22 * Shift2**2 + Gm33 * Shift3**2

    b0u = vdotb / ( Lapse - vdotBeta )

    b0d = b0u * ( BetaSq - Lapse**2 ) + bdotBeta

    bSq = b0d * b0u + b0u * bdotBeta + bdotb

    hb  = bSq / D

  END SUBROUTINE ComputeMagneticEnthalpyFromPrimitive_IDEAL_Scalar


  SUBROUTINE ComputeMagneticEnthalpyFromPrimitive_IDEAL_Vector &
    ( D, V1, V2, V3,    &
      B1, B2, B3,       &
      Gm11, Gm22, Gm33, &
      Lapse, Shift1, Shift2, Shift3, hb )

    REAL(DP), INTENT(in)  :: D(:), V1(:), V2(:), V3(:), &
                             B1(:), B2(:), B3(:),       &
                             Gm11(:), Gm22(:), Gm33(:), &
                             Lapse(:), Shift1(:), Shift2(:), Shift3(:)
    REAL(DP), INTENT(out) :: hb(:)

    INTEGER :: i

    DO i = 1, SIZE( D )

      CALL ComputeMagneticEnthalpyFromPrimitive_IDEAL_Scalar &
             ( D(i), V1(i), V2(i), V3(i), B1(i), B2(i), B3(i), &
               Gm11(i), Gm22(i), Gm33(i),                      &
               Lapse(i), Shift1(i), Shift2(i), Shift3(i), hb(i) )

    END DO

  END SUBROUTINE ComputeMagneticEnthalpyFromPrimitive_IDEAL_Vector


  SUBROUTINE ComputePressureFromPrimitive_IDEAL_Scalar( D, Ev, Ne, P )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: P

    P = ( Gamma_IDEAL - 1.0_DP ) * Ev

  END SUBROUTINE ComputePressureFromPrimitive_IDEAL_Scalar


  SUBROUTINE ComputePressureFromPrimitive_IDEAL_Vector( D, Ev, Ne, P )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: P(:)

    P(:) = ( Gamma_IDEAL - 1.0_DP ) * Ev(:)

  END SUBROUTINE ComputePressureFromPrimitive_IDEAL_Vector


  SUBROUTINE ComputeMagneticPressureFromPrimitive_IDEAL_Scalar &
               ( V1, V2, V3,       &
                 B1, B2, B3,       &
                 Gm11, Gm22, Gm33, &
                 Lapse, Shift1, Shift2, Shift3, Pb )

    REAL(DP), INTENT(in)  :: V1, V2, V3,       &
                             B1, B2, B3,       &
                             Gm11, Gm22, Gm33, &
                             Lapse, Shift1, Shift2, Shift3
    REAL(DP), INTENT(out) :: Pb

    REAL(DP) :: vdotb, vdotBeta, bdotBeta, bdotb, BetaSq, b0u, b0d, bSq

    vdotb    = Gm11 * V1 * B1 + Gm22 * V2 * B2 + Gm33 * V3 * B3

    bdotb    = Gm11 * B1**2 + Gm22 * B2**2 + Gm33 * B3**2

    vdotBeta = Gm11 * V1 * Shift1 + Gm22 * V2 * Shift2 + Gm33 * V3 * Shift3

    bdotBeta = Gm11 * B1 * Shift1 + Gm22 * B2 * Shift2 + Gm33 * B3 * Shift3

    BetaSq   = Gm11 * Shift1**2 + Gm22 * Shift2**2 + Gm33 * Shift3**2

    b0u = vdotb / ( Lapse - vdotBeta )

    b0d = b0u * ( BetaSq - Lapse**2 ) + bdotBeta

    bSq = b0d * b0u + b0u * bdotBeta + bdotb

    Pb = bSq / Two

  END SUBROUTINE ComputeMagneticPressureFromPrimitive_IDEAL_Scalar


  SUBROUTINE ComputeMagneticPressureFromPrimitive_IDEAL_Vector &
               ( V1, V2, V3,       &
                 B1, B2, B3,       &
                 Gm11, Gm22, Gm33, &
                 Lapse, Shift1, Shift2, Shift3, Pb )

    REAL(DP), INTENT(in)  :: V1(:), V2(:), V3(:),       &
                             B1(:), B2(:), B3(:),       &
                             Gm11(:), Gm22(:), Gm33(:), &
                             Lapse(:), Shift1(:), Shift2(:), Shift3(:)
    REAL(DP), INTENT(out) :: Pb(:)

    INTEGER :: i

    DO i = 1, SIZE( V1 )

      CALL ComputeMagneticPressureFromPrimitive_IDEAL_Scalar &
             ( V1(i), V2(i), V3(i),       &
               B1(i), B2(i), B3(i),       &
               Gm11(i), Gm22(i), Gm33(i), &
               Lapse(i), Shift1(i), Shift2(i), Shift3(i), Pb(i) )

    END DO

  END SUBROUTINE ComputeMagneticPressureFromPrimitive_IDEAL_Vector


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_IDEAL_Scalar &
    ( D, Em, Y, P )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Em, Y
    REAL(DP), INTENT(out) :: P

    P = ( Gamma_IDEAL - 1.0_DP ) * D * Em

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_IDEAL_Scalar


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_IDEAL_Vector &
    ( D, Em, Y, P )

    REAL(DP), INTENT(in)  :: D(:), Em(:), Y(:)
    REAL(DP), INTENT(out) :: P(:)

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputePressureFromSpecificInternalEnergy_IDEAL &
             ( D(iP), Em(iP), Y(iP), P(iP) )

    END DO

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_IDEAL_Vector


#ifdef HYDRO_RELATIVISTIC

  SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL_Scalar( D, Ev, Ne, Cs )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: Cs

    Cs = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - One ) * Ev &
                 / ( D + Gamma_IDEAL * Ev ) )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL_Scalar


  SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL_Vector( D, Ev, Ne, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: Cs(:)

    Cs = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - One ) * Ev &
                 / ( D + Gamma_IDEAL * Ev ) )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL_Vector

#else

  SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL_Scalar( D, Ev, Ne, Cs )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: Cs

    Cs = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - One ) * Ev / D )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL_Scalar


  SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL_Vector( D, Ev, Ne, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: Cs(:)

    Cs(:) = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - One ) * Ev(:) / D(:) )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL_Vector

#endif

  SUBROUTINE ComputeAlfvenSpeedFromPrimitive_IDEAL_Scalar &
    ( D, V1, V2, V3, E, Ne, &
      B1, B2, B3,           &
      Gm11, Gm22, Gm33,     &
      Lapse, Shift1, Shift2, Shift3, Ca )

    REAL(DP), INTENT(in)  :: D, V1, V2, V3, E, Ne, &
                             B1, B2, B3,           &
                             Gm11, Gm22, Gm33,     &
                             Lapse, Shift1, Shift2, Shift3
    REAL(DP), INTENT(out) :: Ca

    REAL(DP) :: h, hb

    CALL ComputeEnthalpyFromPrimitive_IDEAL( D, E, Ne, h )

    CALL ComputeMagneticEnthalpyFromPrimitive_IDEAL &
           ( D, V1, V2, V3, B1, B2, B3, &
             Gm11, Gm22, Gm33,          &
             Lapse, Shift1, Shift2, Shift3, hb )

    Ca = SQRT( D * hb / ( D * ( h + hb ) ) )

  END SUBROUTINE ComputeAlfvenSpeedFromPrimitive_IDEAL_Scalar


  SUBROUTINE ComputeAlfvenSpeedFromPrimitive_IDEAL_Vector &
    ( D, V1, V2, V3, E, Ne, B1, B2, B3, &
      Gm11, Gm22, Gm33,                 &
      Lapse, Shift1, Shift2, Shift3, Ca )

    REAL(DP), INTENT(in)  :: D(:), V1(:), V2(:), V3(:), E(:), Ne(:), &
                             B1(:), B2(:), B3(:),                    &
                             Gm11(:), Gm22(:), Gm33(:),              &
                             Lapse(:), Shift1(:), Shift2(:), Shift3(:)
    REAL(DP), INTENT(out) :: Ca(:)

    INTEGER :: i

    DO i = 1, SIZE( D )

      CALL ComputeAlfvenSpeedFromPrimitive_IDEAL_Scalar &
             ( D(i), V1(i), V2(i), V3(i), E(i), Ne(i), &
               B1(i), B2(i), B3(i),                    &
               Gm11(i), Gm22(i), Gm33(i),              &
               Lapse(i), Shift1(i), Shift2(i), Shift3(i), Ca(i) )

    END DO

  END SUBROUTINE ComputeAlfvenSpeedFromPrimitive_IDEAL_Vector


  SUBROUTINE ComputeAuxiliary_Fluid_IDEAL_Scalar &
    ( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: P, T, Y, S, Em, Gm, Cs

    P  = ( Gamma_IDEAL - 1.0_DP ) * Ev
    Gm = Gamma_IDEAL
    Em = Ev / D
    CALL ComputeSoundSpeedFromPrimitive_IDEAL( D, Ev, Ne, Cs )

    T = 0.0_DP
    Y = 0.0_DP
    S = 0.0_DP

  END SUBROUTINE ComputeAuxiliary_Fluid_IDEAL_Scalar


  SUBROUTINE ComputeAuxiliary_Fluid_IDEAL_Vector &
    ( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: P(:), T (:), Y (:), S(:), Em(:), Gm(:), Cs(:)

    INTEGER :: i

    DO i = 1, SIZE( D )

      CALL ComputeAuxiliary_Fluid_IDEAL_Scalar &
             ( D(i), Ev(i), Ne(i), P(i), T(i), Y(i), S(i), Em(i), Gm(i), Cs(i) )

    END DO

  END SUBROUTINE ComputeAuxiliary_Fluid_IDEAL_Vector


  SUBROUTINE ComputeAuxiliary_Magnetofluid_IDEAL_Scalar &
    ( D, V1, V2, V3, Ev, Ne,         &
      B1, B2, B3,                    &
      Gm11, Gm22, Gm33,              &
      Lapse, Shift1, Shift2, Shift3, &
      P, Pb, T, Y, S, Em, h, hb, Gm, Cs, Ca )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, V1, V2, V3, Ev, Ne, &
                             B1, B2, B3,            &
                             Gm11, Gm22, Gm33,      &
                             Lapse, Shift1, Shift2, Shift3
    REAL(DP), INTENT(out) :: P, Pb, T, Y, S, Em, h, hb, Gm, Cs, Ca

    P  = ( Gamma_IDEAL - 1.0_DP ) * Ev
    Pb = Zero
    Gm = Gamma_IDEAL
    Em = Ev / D
    CALL ComputeSoundSpeedFromPrimitive_IDEAL( D, Ev, Ne, Cs )
    CALL ComputeAlfvenSpeedFromPrimitive_IDEAL &
           ( D, V1, V2, V3, Ev, Ne, &
             B1, B2, B3,           &
             Gm11, Gm22, Gm33,     &
             Lapse, Shift1, Shift2, Shift3, Ca)
    CALL ComputeEnthalpyFromPrimitive_IDEAL( D, Ev, Ne, h )
    CALL ComputeMagneticEnthalpyFromPrimitive_IDEAL &
           ( D, V1, V2, V3, B1, B2, V3, &
             Gm11, Gm22, Gm33,          &
             Lapse, Shift1, Shift2, Shift3, hb )

    T = 0.0_DP
    Y = 0.0_DP
    S = 0.0_DP

  END SUBROUTINE ComputeAuxiliary_Magnetofluid_IDEAL_Scalar


  SUBROUTINE ComputeAuxiliary_Magnetofluid_IDEAL_Vector &
    ( D, V1, V2, V3, Ev, Ne,         &
      B1, B2, B3,                    &
      Gm11, Gm22, Gm33,              &
      Lapse, Shift1, Shift2, Shift3, &
      P, Pb, T, Y, S, Em, h, hb, Gm, Cs, Ca )

    REAL(DP), INTENT(in)  :: D(:), V1(:), V2(:), V3(:), Ev(:), Ne(:), &
                             B1(:), B2(:), B3(:),                     &
                             Gm11(:), Gm22(:), Gm33(:),               &
                             Lapse(:), Shift1(:), Shift2(:), Shift3(:)
    REAL(DP), INTENT(out) :: P(:), Pb(:), T (:), Y (:), S(:), Em(:), &
                             h(:), hb(:), Gm(:), Cs(:), Ca(:)

    INTEGER :: i

    DO i = 1, SIZE( D )

      CALL ComputeAuxiliary_Magnetofluid_IDEAL_Scalar &
             ( D(i), V1(i), V2(i), V3(i), Ev(i), Ne(i),   &
               B1(i), B2(i), B3(i),                       &
               Gm11(i), Gm22(i), Gm33(i),                 &
               Lapse(i), Shift1(i), Shift2(i), Shift3(i), &
               P(i), Pb(i), T(i), Y(i), S(i), Em(i), h(i), hb(i), &
               Gm(i), Cs(i), Ca(i) )

    END DO

  END SUBROUTINE ComputeAuxiliary_Magnetofluid_IDEAL_Vector


END MODULE EquationOfStateModule_IDEAL
