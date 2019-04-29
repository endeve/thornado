PROGRAM ComputeLeftState_NuclearEOS

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Centimeter, &
    Kilometer, &
    Gram, &
    Second, &
    Erg, &
    Kelvin
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    ComputeTemperatureFromPressure_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE, &
    ComputeSpecificInternalEnergy_TABLE
  USE UtilitiesModule_NuclearEOS, ONLY: &
    ComputeFVEC, &
    ComputeFJAC

  IMPLICIT NONE

  LOGICAL                  :: Converged
  INTEGER                  :: Iter, INFO
  INTEGER,  DIMENSION(4)   :: IPIV
  REAL(DP), PARAMETER      :: Tol = 1.0d-10
  REAL(DP), DIMENSION(1)   :: D_L, V_L, T_L, Y_L
  REAL(DP), DIMENSION(1)   :: P_L, E_L, Ev_L, N_L, S_L, Gm_L, Cs_L
  REAL(DP), DIMENSION(1)   :: D_R, V_R, T_R, Y_R
  REAL(DP), DIMENSION(1)   :: P_R, E_R, Ev_R, N_R, S_R, Gm_R, Cs_R
  REAL(DP), DIMENSION(1)   :: Mach, V_Sh
  REAL(DP), DIMENSION(1)   :: dPdD, dPdT, dPdY
  REAL(DP), DIMENSION(1)   :: dEdD, dEdT, dEdY
  REAL(DP), DIMENSION(4)   :: FVEC, UVEC, dUVEC
  REAL(DP), DIMENSION(4,4) :: FJAC

  CALL InitializeEquationOfState_TABLE( 'EquationOfStateTable.h5' )

  ! --- Right State ---

  D_R = 1.0d12 * Gram / Centimeter**3
  V_R = 0.0d00 * Kilometer / Second
  P_R = 1.0d31 * Erg / Centimeter**3
  Y_R = 0.3_DP

  CALL ComputeTemperatureFromPressure_TABLE &
         ( D_R, P_R, Y_R, T_R )
  CALL ComputeThermodynamicStates_Primitive_TABLE &
         ( D_R, T_R, Y_R, Ev_R, E_R, N_R )
  CALL ComputeAuxiliary_Fluid_TABLE &
         ( D_R, Ev_R, N_R, P_R, T_R, Y_R, S_R, E_R, Gm_R, Cs_R )

  WRITE(*,*)
  WRITE(*,'(A4,A)') '', 'Right State:'
  WRITE(*,*)
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'D_R  =', D_R / ( Gram / Centimeter**3 )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'V_R  =', V_R / ( Kilometer / Second )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'T_R  =', T_R / Kelvin
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'Y_R  =', Y_R
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'E_R  =', E_R / ( Erg / Gram )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'P_R  =', P_R / ( Erg / Centimeter**3 )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'Gm_R =', Gm_R
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'Cs_R =', Cs_R / ( Kilometer / Second )

  Mach = 3.0_DP
  V_Sh = Mach * Cs_R

  WRITE(*,*)
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'Mach =', Mach
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'V_Sh =', V_Sh / ( Kilometer / Second )

  ! --- Initial Guess for Left State ---

  D_L = 4.0_DP * D_R
  V_L = 0.10_DP
  P_L = 10.0_DP * P_R
  Y_L = Y_R

  CALL ComputeTemperatureFromPressure_TABLE &
         ( D_L, P_L, Y_L, T_L )
  CALL ComputeThermodynamicStates_Primitive_TABLE &
         ( D_L, T_L, Y_L, Ev_L, E_L, N_L )
  CALL ComputeAuxiliary_Fluid_TABLE &
         ( D_L, Ev_L, N_L, P_L, T_L, Y_L, S_L, E_L, Gm_L, Cs_L )

  WRITE(*,*)
  WRITE(*,'(A4,A)') '', 'Left State (Initial Guess):'
  WRITE(*,*)
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'D_L  =', D_L / ( Gram / Centimeter**3 )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'V_L  =', V_L / ( Kilometer / Second )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'T_L  =', T_L / Kelvin
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'Y_L  =', Y_L
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'E_L  =', E_L / ( Erg / Gram )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'P_L  =', P_L / ( Erg / Centimeter**3 )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'Gm_L =', Gm_L
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'Cs_L =', Cs_L / ( Kilometer / Second )

  Converged = .FALSE.
  Iter      = 0
  DO WHILE( .NOT. Converged )

    Iter = Iter + 1

    CALL ComputePressure_TABLE &
           ( D_L, T_L, Y_L, P_L, dPdD, dPdT, dPdY )

    CALL ComputeSpecificInternalEnergy_TABLE &
           ( D_L, T_L, Y_L, E_L, dEdD, dEdT, dEdY )

    CALL ComputeFVEC &
           ( D_R(1), V_R(1), P_R(1), E_R(1), Y_R(1), &
             D_L(1), V_L(1), P_L(1), E_L(1), Y_L(1), &
             V_Sh(1), FVEC )

    CALL ComputeFJAC &
           ( D_R(1), V_R(1), P_R(1), E_R(1), Y_R(1), &
             D_L(1), V_L(1), P_L(1), E_L(1), Y_L(1), &
             dPdD(1), dPdT(1), dPdY(1), &
             dEdD(1), dEdT(1), dEdY(1), &
             V_Sh(1), FJAC )

    UVEC  = [ D_L(1), V_L(1), T_L(1), Y_L(1) ]
    dUVEC = FVEC

    CALL DGESV( 4, 1, FJAC, 4, IPIV, dUVEC, 4, INFO )

    WRITE(*,*)
    WRITE(*,'(A6,A12,I2.2)') '', 'Iteration: ', Iter
    WRITE(*,'(A6,A12,4ES20.10E3)') '', '|FVEC| = ', &
      ABS( FVEC )
    WRITE(*,'(A6,A12,4ES20.10E3)') '', '|dU/U| = ', &
      ABS(  dUVEC / ( UVEC + TINY(1.0_DP) ) )
    WRITE(*,*)

    D_L(1) = UVEC(1) - dUVEC(1)
    V_L(1) = UVEC(2) - dUVEC(2)
    T_L(1) = UVEC(3) - dUVEC(3)
    Y_L(1) = UVEC(4) - dUVEC(4)

    IF( ALL( ABS(  dUVEC / ( UVEC + TINY(1.0_DP) ) ) < Tol ) ) &
      Converged = .TRUE.

  END DO

  CALL ComputePressure_TABLE &
         ( D_L, T_L, Y_L, P_L )
  CALL ComputeSpecificInternalEnergy_TABLE &
         ( D_L, T_L, Y_L, E_L )
  CALL ComputeThermodynamicStates_Primitive_TABLE &
         ( D_L, T_L, Y_L, Ev_L, E_L, N_L )
  CALL ComputeAuxiliary_Fluid_TABLE &
         ( D_L, Ev_L, N_L, P_L, T_L, Y_L, S_L, E_L, Gm_L, Cs_L )

  WRITE(*,*)
  WRITE(*,'(A4,A)') '', 'Left State:'
  WRITE(*,*)
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'D_L  =', D_L / ( Gram / Centimeter**3 )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'V_L  =', V_L / ( Kilometer / Second )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'T_L  =', T_L / Kelvin
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'Y_L  =', Y_L
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'E_L  =', E_L / ( Erg / Gram )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'P_L  =', P_L / ( Erg / Centimeter**3 )
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'Gm_L =', Gm_L
  WRITE(*,'(A4,A8,ES20.10E3)') '', 'Cs_L =', Cs_L / ( Kilometer / Second )
  WRITE(*,*)

END PROGRAM ComputeLeftState_NuclearEOS
