MODULE EquationOfStateModule_IDEAL

  USE KindModule, ONLY: &
    DP, One
  USE FluidFieldsModule, ONLY: &
    ! --- Primitive Fluid Fields:
    iPF_D, iPF_E, iPF_Ne, nPF, &
    ! --- Auxiliary Fluid Fields:
    iAF_P, iAF_T, iAF_Ye, iAF_E, iAF_Gm, iAF_Cs, nAF

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
  PUBLIC :: Auxiliary_Fluid_IDEAL

  INTERFACE ComputePressureFromPrimitive_IDEAL
    MODULE PROCEDURE ComputePressureFromPrimitive_IDEAL_Scalar
    MODULE PROCEDURE ComputePressureFromPrimitive_IDEAL_Vector
  END INTERFACE ComputePressureFromPrimitive_IDEAL


CONTAINS


  SUBROUTINE InitializeEquationOfState_IDEAL( Gamma_IDEAL_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Gamma_IDEAL_Option

    IF( PRESENT( Gamma_IDEAL_Option ) )THEN
      Gamma_IDEAL = Gamma_IDEAL_Option
    ELSE
      Gamma_IDEAL = 5.0_DP / 3.0_DP
    END IF

  END SUBROUTINE InitializeEquationOfState_IDEAL


  SUBROUTINE FinalizeEquationOfState_IDEAL

  END SUBROUTINE FinalizeEquationOfState_IDEAL


  SUBROUTINE ComputeInternalEnergyDensityFromPressure_IDEAL( D, P, Y, Ev )

    REAL(DP), INTENT(in)  :: D(:), P(:), Y(:)
    REAL(DP), INTENT(out) :: Ev(:)

    Ev(:) = P(:) / ( Gamma_IDEAL - 1.0_DP )

  END SUBROUTINE ComputeInternalEnergyDensityFromPressure_IDEAL


  SUBROUTINE ComputePressureFromPrimitive_IDEAL_Scalar( D, Ev, Ne, P )

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: P

    P = ( Gamma_IDEAL - 1.0_DP ) * Ev

  END SUBROUTINE ComputePressureFromPrimitive_IDEAL_Scalar


  SUBROUTINE ComputePressureFromPrimitive_IDEAL_Vector( D, Ev, Ne, P )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: P(:)

    P(:) = ( Gamma_IDEAL - 1.0_DP ) * Ev(:)

  END SUBROUTINE ComputePressureFromPrimitive_IDEAL_Vector


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_IDEAL( D, Em, Y, P )

    REAL(DP), INTENT(in)  :: D(:), Em(:), Y(:)
    REAL(DP), INTENT(out) :: P(:)

    P(:) = ( Gamma_IDEAL - 1.0_DP ) * D(:) * Em(:)

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_IDEAL


#if defined HYDRO_NONRELATIVISTIC

  SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL( D, Ev, Ne, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: Cs(:)

    Cs(:) = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - One ) * Ev(:) / D(:) )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL

#elif defined HYDRO_RELATIVISTIC

  SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL( D, Ev, Ne, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: Cs(:)

    Cs = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - One ) * Ev &
                 / ( D + Gamma_IDEAL * Ev ) )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL

#else

  SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL( D, Ev, Ne, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: Cs(:)

    Cs(:) = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - One ) * Ev(:) / D(:) )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL

#endif


  SUBROUTINE ComputeAuxiliary_Fluid_IDEAL( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: P(:), T(:), Y(:), S(:), Em(:), Gm(:), Cs(:)

    P (:) = ( Gamma_IDEAL - 1.0_DP ) * Ev(:)
    Gm(:) = Gamma_IDEAL
    Em(:) = Ev(:) / D(:)
    CALL ComputeSoundSpeedFromPrimitive_IDEAL( D(:), Ev(:), Ne(:), Cs(:) )

  END SUBROUTINE ComputeAuxiliary_Fluid_IDEAL


  FUNCTION Auxiliary_Fluid_IDEAL( PF )

    REAL(DP), INTENT(in) :: PF(nPF)
    REAL(DP)             :: Auxiliary_Fluid_IDEAL(nAF)

    Auxiliary_Fluid_IDEAL(iAF_P)  &
      = ( Gamma_IDEAL - 1.0_DP ) * PF(iPF_E)
    Auxiliary_Fluid_IDEAL(iAF_Gm) &
      = Gamma_IDEAL
    Auxiliary_Fluid_IDEAL(iAF_E)  &
      = PF(iPF_E) / PF(iPF_D)
    Auxiliary_Fluid_IDEAL(iAF_Cs) &
      = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - 1.0_DP ) &
                * PF(iPF_E) / PF(iPF_D) )

    RETURN
  END FUNCTION Auxiliary_Fluid_IDEAL


END MODULE EquationOfStateModule_IDEAL
