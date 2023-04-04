MODULE EquationOfStateModule_IDEAL

  USE KindModule, ONLY: &
    DP, One

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PUBLIC :: &
    Gamma_IDEAL

  PUBLIC :: InitializeEquationOfState_IDEAL
  PUBLIC :: FinalizeEquationOfState_IDEAL
  PUBLIC :: ComputeInternalEnergyDensityFromPressure_IDEAL
  PUBLIC :: ComputePressureFromPrimitive_IDEAL
  PUBLIC :: ComputePressureFromSpecificInternalEnergy_IDEAL
  PUBLIC :: ComputeSoundSpeedFromPrimitive_IDEAL
  PUBLIC :: ComputeAuxiliary_Fluid_IDEAL

  INTERFACE ComputePressureFromPrimitive_IDEAL
    MODULE PROCEDURE ComputePressureFromPrimitive_IDEAL_Scalar
    MODULE PROCEDURE ComputePressureFromPrimitive_IDEAL_Vector
  END INTERFACE ComputePressureFromPrimitive_IDEAL

  INTERFACE ComputePressureFromSpecificInternalEnergy_IDEAL
    MODULE PROCEDURE ComputePressureFromSpecificInternalEnergy_IDEAL_Scalar
    MODULE PROCEDURE ComputePressureFromSpecificInternalEnergy_IDEAL_Vector
  END INTERFACE ComputePressureFromSpecificInternalEnergy_IDEAL

  INTERFACE ComputeInternalEnergyDensityFromPressure_IDEAL
    MODULE PROCEDURE ComputeInternalEnergyDensityFromPressure_IDEAL_Scalar
    MODULE PROCEDURE ComputeInternalEnergyDensityFromPressure_IDEAL_Vector
  END INTERFACE ComputeInternalEnergyDensityFromPressure_IDEAL

  INTERFACE ComputeSoundSpeedFromPrimitive_IDEAL
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_IDEAL_Scalar
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_IDEAL_Vector
  END INTERFACE ComputeSoundSpeedFromPrimitive_IDEAL

  INTERFACE ComputeAuxiliary_Fluid_IDEAL
    MODULE PROCEDURE ComputeAuxiliary_Fluid_IDEAL_Scalar
    MODULE PROCEDURE ComputeAuxiliary_Fluid_IDEAL_Vector
  END INTERFACE ComputeAuxiliary_Fluid_IDEAL

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


END MODULE EquationOfStateModule_IDEAL
