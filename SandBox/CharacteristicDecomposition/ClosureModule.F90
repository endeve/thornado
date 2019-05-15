MODULE ClosureModule

  USE KindModule, ONLY: &
    DP, SqrtTiny, Zero, &
    One, Two, Three, &
    Four, Fifth, Third

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: PartialDerivativeEddingtonFactor

CONTAINS

  SUBROUTINE PartialDerivativeEddingtonFactor &
    ( D, FF, pEFpD, pEFpFF, Verbose_Option )

    REAL(DP), INTENT(in)  :: D, FF  
    REAL(DP), INTENT(out) :: pEFpD  ! wt FF constant
    REAL(DP), INTENT(out) :: pEFpFF ! wt D constant
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option
    
    LOGICAL :: Verbose
    REAL(DP) :: XX

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    XX = FF / MAX( One - D, SqrtTiny )

#ifdef MOMENT_CLOSURE_MINERBO

    ! --- Maximum Entropy (ME) Minerbo Closure ---
   
    IF ( Verbose ) THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'Minerbo : '
    END IF

    pEFpD  = Zero

    pEFpFF = Two * Third &
              * DerivativeClosurePolynomial_ME_CB &
                ( FF ) 

#elif MOMENT_CLOSURE_MAXIMUM_ENTROPY_CB

    ! --- Cernohorsky-Bludman ME Closure ---

    IF ( Verbose ) THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'Cernohorsky-Bludman :'
    END IF

    pEFpD  = Two * Third &
             * ( ClosurePolynomial_ME_CB( XX ) &
                 * ( Four * D - Three ) &
               + DerivativeClosurePolynomial_ME_CB( XX ) &
                 * xx * ( One - Two * D ) ) 
    pEFpFF = Two * Third * ( One - Two * D) &
              * DerivativeClosurePolynomial_ME_CB( XX )

#elif MOMENT_CLOSURE_MAXIMUM_ENTROPY_BL

    ! --- Banach-Larecki ME Closure ---

    IF ( Verbose ) THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'Banach-Larecki :'
    END IF

    pEFpD  = Two * Third &
             * ( ClosurePolynomial_ME_BL( XX ) &
                 * ( Four * D - Three ) &
               + DerivativeClosurePolynomial_ME_BL( XX ) &
                 * xx * ( One - Two * D ) )
    pEFpFF = Two * Third * ( One - Two * D) &
              * DerivativeClosurePolynomial_ME_BL( XX )

#elif MOMENT_CLOSURE_KERSHAW_BL

    ! --- Banach-Larecki Kershaw Closure ---
 
    IF ( Verbose ) THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'Kershaw :'
    END IF

    pEFpD  = Two * Third &
             * ( ClosurePolynomial_KE_BL( XX ) &
                 * ( Four * D - Three ) &
               + DerivativeClosurePolynomial_KE_BL( XX ) &
                 * xx * ( One - Two * D ) )
    pEFpFF = Two * Third * ( One - Two * D) &
              * DerivativeClosurePolynomial_KE_BL( XX )

#else

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'PartialDerivativeEddingtonFactor : Two-Moment Closure Unknown'
    WRITE(*,*)
    STOP
  
#endif

    IF ( Verbose ) THEN
      WRITE(*,'(A10,2A15)') '', 'pEFpD','pEFpFF'
      WRITE(*,'(A10,2ES15.4)') '', pEFpD, pEFpFF
    END IF

  END SUBROUTINE PartialDerivativeEddingtonFactor

  ! --- Derivative of Closure Polynomials ---


  PURE REAL(DP) ELEMENTAL FUNCTION DerivativeClosurePolynomial_ME_CB( X )

    ! --- Cernohorsky-Bludman Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X

    DerivativeClosurePolynomial_ME_CB &
      = Three * Fifth * X * ( Two - X + Four * X * X )

    RETURN

  END FUNCTION DerivativeClosurePolynomial_ME_CB


  PURE REAL(DP) ELEMENTAL FUNCTION DerivativeClosurePolynomial_ME_BL( X )

    ! --- Banach-Larecki Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X

    REAL(DP) :: aa

    aa = SQRT( 33.0_DP * X**4 - 42.0_DP * X * X + 25.0_DP )

    DerivativeClosurePolynomial_ME_BL &
      = Three * x * ( 11.0_DP * X * X + Three * X - 7.0_DP) / ( Four * aa )

  END FUNCTION DerivativeClosurePolynomial_ME_BL


  PURE REAL(DP) ELEMENTAL FUNCTION DerivativeClosurePolynomial_KE_BL( X )

    ! --- Banach-Larecki Kershaw Closure ---

    REAL(DP), INTENT(in) :: X

    DerivativeClosurePolynomial_KE_BL &
      = 2 * X

  END FUNCTION DerivativeClosurePolynomial_KE_BL
  

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
