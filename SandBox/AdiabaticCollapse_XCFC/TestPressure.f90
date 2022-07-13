PROGRAM TestPressure

  USE KindModule, ONLY: &
    DP
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState,   &
    ComputePressureFromSpecificInternalEnergy
  USE UnitsModule, ONLY: &
    Gram, Centimeter, Erg

  IMPLICIT NONE

  INTEGER  :: i
  INTEGER, PARAMETER :: N = 1000

  REAL(DP) :: rhoMin, rhoMax, drho, rho(N)
  REAL(DP) :: epsMin, epsMax, deps, eps(N)
  REAL(DP) :: YeMin , YeMax , dYe , Ye(N)
  REAL(DP) :: p(N)

  CALL InitializeEquationOfState &
         ( EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = 'wl-EOS-SFHo-25-50-100.h5' )

rho(1) = 9.314925e13_DP * ( Gram / Centimeter**3 )
eps(1) = 5.465183e19_DP * ( Erg / Gram )
Ye (1) = 0.442_DP

rho(2) = 9.266702e13_DP * ( Gram / Centimeter**3 )
eps(2) = 5.299234e19_DP * ( Erg / Gram )
Ye (2) = 0.442_DP

    CALL ComputePressureFromSpecificInternalEnergy &
           ( rho(1), eps(1), Ye(1), p(1) )
    PRINT '(A,ES14.6E3)', 'p(1) = ', p(1) / ( Erg / Centimeter**3 )

    CALL ComputePressureFromSpecificInternalEnergy &
           ( rho(2), eps(2), Ye(2), p(2) )
    PRINT '(A,ES14.6E3)', 'p(2) = ', p(2) / ( Erg / Centimeter**3 )
stop

  rhoMin = 9.0e13_DP * ( Gram / Centimeter**3 )
  rhoMax = 1.0e14_DP * ( Gram / Centimeter**3 )
  drho   = ( rhoMax - rhoMin ) / DBLE( N )

  epsMin = 5.0e19_DP * ( Erg / Gram )
  epsMax = 6.0e19_DP * ( Erg / Gram )
  deps   = ( epsMax - epsMin ) / DBLE( N )

  YeMin  = 0.4_DP
  YeMax  = 0.5_DP
  dYe    = ( YeMax - YeMin ) / DBLE( N )

  rho(1) = rhoMin
  eps(1) = epsMin
  Ye (1) = YeMin

  DO i = 2, N

    rho(i) = rho(i-1) + drho
    eps(i) = eps(i-1) + deps
    Ye (i) = Ye (i-1) + dYe

  END DO

  DO i = 1, N

    CALL ComputePressureFromSpecificInternalEnergy &
           ( rho(i), eps(i), Ye(i), p(i) )

    PRINT '(A,ES14.6E3)', 'p(i) = ', p(i) / ( Erg / Centimeter**3 )

  END DO

  CALL FinalizeEquationOfState

END PROGRAM TestPressure
