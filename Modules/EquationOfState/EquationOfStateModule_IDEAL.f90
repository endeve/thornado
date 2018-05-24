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
  PUBLIC :: ComputeSoundSpeedFromPrimitive_GR_IDEAL
  PUBLIC :: ComputeAuxiliary_Fluid_IDEAL
  PUBLIC :: Auxiliary_Fluid_IDEAL

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

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, P, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Ev

    Ev(:) = P(:) / ( Gamma_IDEAL - 1.0_DP )

  END SUBROUTINE ComputeInternalEnergyDensityFromPressure_IDEAL


  SUBROUTINE ComputePressureFromPrimitive_IDEAL( D, Ev, Ne, P )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: P

    P(:) = ( Gamma_IDEAL - 1.0_DP ) * Ev(:)

  END SUBROUTINE ComputePressureFromPrimitive_IDEAL


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_IDEAL( D, Em, Y, P )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Em, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: P

    P(:) = ( Gamma_IDEAL - 1.0_DP ) * D(:) * Em(:)

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_IDEAL


  SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL( D, Ev, Ne, Cs )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: Cs

    Cs(:) = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - One ) * Ev(:) / D(:) )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_IDEAL


  SUBROUTINE ComputeSoundSpeedFromPrimitive_GR_IDEAL( D, Ev, Ne, Cs )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: Cs

    Cs(:) = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - One ) * Ev(:) &
                    / ( D(:) + Gamma_IDEAL * Ev(:) ) )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_GR_IDEAL


  SUBROUTINE ComputeAuxiliary_Fluid_IDEAL( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: P, T, Y, S, Em, Gm, Cs

    P (:) = ( Gamma_IDEAL - 1.0_DP ) * Ev(:)
    Gm(:) = Gamma_IDEAL
    Em(:) = Ev(:) / D(:)
    Cs(:) = SQRT( Gamma_IDEAL * P(:) / D(:) )

  END SUBROUTINE ComputeAuxiliary_Fluid_IDEAL


  FUNCTION Auxiliary_Fluid_IDEAL( PF )

    REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
    REAL(DP), DIMENSION(nAF)             :: Auxiliary_Fluid_IDEAL

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
