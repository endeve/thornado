MODULE ClosureModule

  USE KindModule, ONLY: &
    DP, SqrtTiny, Zero, &
    One, Two, Three, &
    Four, Fifth, Third

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeEddingtonFactorDerivatives

CONTAINS


  SUBROUTINE ComputeEddingtonFactorDerivatives( D, FF, dEFdD_FF, dEFdFF_D )

    REAL(DP), INTENT(in)  :: D, FF
    REAL(DP), INTENT(out) :: dEFdD_FF ! --- wrt D, constant FF
    REAL(DP), INTENT(out) :: dEFdFF_D ! --- wrt FF, constant D
    
    REAL(DP) :: XX

    XX = FF / MAX( One - D, SqrtTiny )

#ifdef MOMENT_CLOSURE_MINERBO

    ! --- Maximum Entropy (ME) Minerbo Closure ---
   
    dEFdD_FF = Zero

    dEFdFF_D = Two * Third * ClosurePolynomialDerivative_ME_CB( FF ) 

#elif MOMENT_CLOSURE_MAXIMUM_ENTROPY_CB

    ! --- Cernohorsky-Bludman ME Closure ---

    dEFdD_FF &
      = Two * Third &
          * ( ClosurePolynomial_ME_CB( XX ) * ( Four * D - Three ) &
              + ClosurePolynomialDerivative_ME_CB( XX ) * XX * ( One - Two * D ) )

    dEFdFF_D &
      = Two * Third * ( One - Two * D ) * ClosurePolynomialDerivative_ME_CB( XX )

#elif MOMENT_CLOSURE_MAXIMUM_ENTROPY_BL

    ! --- Banach-Larecki ME Closure ---

    dEFdD_FF &
      = Two * Third &
          * ( ClosurePolynomial_ME_BL( XX ) * ( Four * D - Three ) &
              + ClosurePolynomialDerivative_ME_BL( XX ) * XX * ( One - Two * D ) )

    dEFdFF_D &
      = Two * Third * ( One - Two * D ) * ClosurePolynomialDerivative_ME_BL( XX )

#elif MOMENT_CLOSURE_KERSHAW_BL

    ! --- Banach-Larecki Kershaw Closure ---
 
    dEFdD_FF &
      = Two * Third &
          * ( ClosurePolynomial_KE_BL( XX ) * ( Four * D - Three ) &
              + ClosurePolynomialDerivative_KE_BL( XX ) * XX * ( One - Two * D ) )

    dEFdFF_D &
      = Two * Third * ( One - Two * D ) * ClosurePolynomialDerivative_KE_BL( XX )

#endif

  END SUBROUTINE ComputeEddingtonFactorDerivatives


  ! --- Derivative of Closure Polynomials ---


  PURE REAL(DP) ELEMENTAL FUNCTION ClosurePolynomialDerivative_ME_CB( X )

    ! --- Cernohorsky-Bludman Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X

    ClosurePolynomialDerivative_ME_CB &
      = Three * Fifth * X * ( Two - X + Four * X * X )

    RETURN
  END FUNCTION ClosurePolynomialDerivative_ME_CB


  PURE REAL(DP) ELEMENTAL FUNCTION ClosurePolynomialDerivative_ME_BL( X )

    ! --- Banach-Larecki Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X

    REAL(DP) :: aa

    aa = SQRT( 33.0_DP * X**4 - 42.0_DP * X * X + 25.0_DP )

    ClosurePolynomialDerivative_ME_BL &
      = Three * x * ( 11.0_DP * X * X + Three * X - 7.0_DP) / ( Four * aa )

    RETURN
  END FUNCTION ClosurePolynomialDerivative_ME_BL


  PURE REAL(DP) ELEMENTAL FUNCTION ClosurePolynomialDerivative_KE_BL( X )

    ! --- Banach-Larecki Kershaw Closure ---

    REAL(DP), INTENT(in) :: X

    ClosurePolynomialDerivative_KE_BL &
      = Two * X

    RETURN
  END FUNCTION ClosurePolynomialDerivative_KE_BL
  

! ====== COPIED FROM TwoMoment_ClosureModule ===== !
 
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
  
END MODULE ClosureModule
