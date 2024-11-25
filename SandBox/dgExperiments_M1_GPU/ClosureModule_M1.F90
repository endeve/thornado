MODULE ClosureModule_M1

  USE KindModule, ONLY: &
    DP, SqrtTiny, &
    One, Two, Three, &
    Fifth, Third

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMomentClosure
  PUBLIC :: FluxFactor
  PUBLIC :: EddingtonFactor

CONTAINS


  SUBROUTINE InitializeMomentClosure

#ifdef MOMENT_CLOSURE_MINERBO

    ! --- Maximum Entropy (ME) Minerbo Closure ---

    WRITE(*,*)
    WRITE(*,'(A6,A)') &
      '', 'Two-Moment Closure: Maximum Entropy (Minerbo)'

#elif  MOMENT_CLOSURE_MAXIMUM_ENTROPY_CB

    ! --- Cernohorsky-Bludman ME Closure ---

    WRITE(*,*)
    WRITE(*,'(A6,A)') &
      '', 'Two-Moment Closure: Maximum Entropy (Cernohorsky & Bludman)'

#elif  MOMENT_CLOSURE_MAXIMUM_ENTROPY_BL

    ! --- Banach-Larecki ME Closure ---

    WRITE(*,*)
    WRITE(*,'(A6,A)') &
      '', 'Two-Moment Closure: Maximum Entropy (Banach & Larecki)'

#elif  MOMENT_CLOSURE_KERSHAW_BL

    ! --- Banach-Larecki Kershaw Closure ---

    WRITE(*,*)
    WRITE(*,'(A6,A)') &
      '', 'Two-Moment Closure: Kershaw (Banach & Larecki)'

#else

    WRITE(*,*)
    WRITE(*,'(A6,A)') &
      '', 'Two-Moment Closure: Unknown'
    WRITE(*,*)
    STOP

#endif

  END SUBROUTINE InitializeMomentClosure


  ! --- Flux Factor ---


  PURE REAL(DP) ELEMENTAL FUNCTION FluxFactor &
    ( D, I_1, I_2, I_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in) :: D, I_1, I_2, I_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    FluxFactor &
      = MIN( MAX( SQRT( Gm_dd_11 * I_1**2 &
                        + Gm_dd_22 * I_2**2 &
                        + Gm_dd_33 * I_3**2 ) &
                    / MAX( D, SqrtTiny ), &
                  SqrtTiny ), &
             One )

    RETURN
  END FUNCTION FluxFactor


  ! --- Eddington Factor ---


  PURE REAL(DP) ELEMENTAL FUNCTION EddingtonFactor( D, FF )

    REAL(DP), INTENT(in) :: D, FF

#ifdef MOMENT_CLOSURE_MINERBO

    ! --- Maximum Entropy (ME) Minerbo Closure ---

    EddingtonFactor &
      = Third + Two * Third * ClosurePolynomial_ME_CB( FF )

#elif  MOMENT_CLOSURE_MAXIMUM_ENTROPY_CB

    ! --- Cernohorsky-Bludman ME Closure ---

    EddingtonFactor &
      = Third + Two * Third * ( One - D ) * ( One - Two * D ) &
          * ClosurePolynomial_ME_CB( FF / MAX( One - D, SqrtTiny ) )

#elif  MOMENT_CLOSURE_MAXIMUM_ENTROPY_BL

    ! --- Banach-Larecki ME Closure ---

    EddingtonFactor &
      = Third + Two * Third * ( One - D ) * ( One - Two * D ) &
          * ClosurePolynomial_ME_BL( FF / MAX( One - D, SqrtTiny ) )

#elif  MOMENT_CLOSURE_KERSHAW_BL

    ! --- Banach-Larecki Kershaw Closure ---

    EddingtonFactor &
      = Third + Two * Third * ( One - D ) * ( One - Two * D ) &
          * ClosurePolynomial_KE_BL( FF / MAX( One - D, SqrtTiny ) )

#endif

    RETURN
  END FUNCTION EddingtonFactor


  ! --- Closure Polynomials ---


  PURE REAL(DP) ELEMENTAL FUNCTION ClosurePolynomial_ME_CB( X )

    ! --- Cernohorsky-Bludman Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X

    ClosurePolynomial_ME_CB &
      = Three * Fifth * X**2 * ( One - Third * X + X**2 )

  END FUNCTION ClosurePolynomial_ME_CB


  PURE REAL(DP) ELEMENTAL FUNCTION ClosurePolynomial_ME_BL( X )

    ! --- Banach-Larecki Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X

    ClosurePolynomial_ME_BL &
      = ( 9.0_DP * X**2 - 5.0_DP &
          + SQRT( 33.0_DP * X**4 - 42.0_DP * X**2 + 25.0_DP ) ) / 8.0_DP

  END FUNCTION ClosurePolynomial_ME_BL


  PURE REAL(DP) ELEMENTAL FUNCTION ClosurePolynomial_KE_BL( X )

    ! --- Banach-Larecki Kershaw Closure ---

    REAL(DP), INTENT(in) :: X

    ClosurePolynomial_KE_BL &
      = X**2

  END FUNCTION ClosurePolynomial_KE_BL


END MODULE ClosureModule_M1
