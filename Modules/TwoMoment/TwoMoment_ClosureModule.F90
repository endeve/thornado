MODULE TwoMoment_ClosureModule

  USE KindModule, ONLY: &
    DP, SqrtTiny, &
    Zero, One, Two, Three, Four, &
    Fifth, Third

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeClosure_TwoMoment
  PUBLIC :: FluxFactor
  PUBLIC :: EddingtonFactor
  PUBLIC :: ComputeEddingtonFactorDerivatives

  INTERFACE FluxFactor
    MODULE PROCEDURE FluxFactor_Scalar
    MODULE PROCEDURE FluxFactor_Vector
  END INTERFACE

  INTERFACE EddingtonFactor
    MODULE PROCEDURE EddingtonFactor_Scalar
    MODULE PROCEDURE EddingtonFactor_Vector
  END INTERFACE

  INTERFACE ClosurePolynomial_ME_CB
    MODULE PROCEDURE ClosurePolynomial_ME_CB_Scalar
    MODULE PROCEDURE ClosurePolynomial_ME_CB_Vector
  END INTERFACE

  INTERFACE ClosurePolynomial_ME_BL
    MODULE PROCEDURE ClosurePolynomial_ME_BL_Scalar
    MODULE PROCEDURE ClosurePolynomial_ME_BL_Vector
  END INTERFACE

  INTERFACE ClosurePolynomial_KE_BL
    MODULE PROCEDURE ClosurePolynomial_KE_BL_Scalar
    MODULE PROCEDURE ClosurePolynomial_KE_BL_Vector
  END INTERFACE

  INTERFACE ClosurePolynomialDerivative_ME_CB
    MODULE PROCEDURE ClosurePolynomialDerivative_ME_CB_Scalar
    MODULE PROCEDURE ClosurePolynomialDerivative_ME_CB_Vector
  END INTERFACE

  INTERFACE ClosurePolynomialDerivative_ME_BL
    MODULE PROCEDURE ClosurePolynomialDerivative_ME_BL_Scalar
    MODULE PROCEDURE ClosurePolynomialDerivative_ME_BL_Vector
  END INTERFACE

  INTERFACE ClosurePolynomialDerivative_KE_BL
    MODULE PROCEDURE ClosurePolynomialDerivative_KE_BL_Scalar
    MODULE PROCEDURE ClosurePolynomialDerivative_KE_BL_Vector
  END INTERFACE

CONTAINS


  SUBROUTINE InitializeClosure_TwoMoment( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

#ifdef MOMENT_CLOSURE_MINERBO

    ! --- Maximum Entropy (ME) Minerbo Closure ---

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A6,A)') &
        '', 'Two-Moment Closure: Maximum Entropy (Minerbo)'

    END IF

#elif  MOMENT_CLOSURE_MAXIMUM_ENTROPY_CB

    ! --- Cernohorsky-Bludman ME Closure ---

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A6,A)') &
        '', 'Two-Moment Closure: Maximum Entropy (Cernohorsky & Bludman)'

    END IF

#elif  MOMENT_CLOSURE_MAXIMUM_ENTROPY_BL

    ! --- Banach-Larecki ME Closure ---

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A6,A)') &
        '', 'Two-Moment Closure: Maximum Entropy (Banach & Larecki)'

    END IF

#elif  MOMENT_CLOSURE_KERSHAW_BL

    ! --- Banach-Larecki Kershaw Closure ---

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A6,A)') &
        '', 'Two-Moment Closure: Kershaw (Banach & Larecki)'

    END IF

#else

    WRITE(*,*)
    WRITE(*,'(A6,A)') &
      '', 'Two-Moment Closure: Unknown'
    WRITE(*,*)
    STOP

#endif

  END SUBROUTINE InitializeClosure_TwoMoment


  ! --- Flux Factor ---


  FUNCTION FluxFactor_Scalar &
    ( D, I_1, I_2, I_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
    RESULT( FluxFactor )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D, I_1, I_2, I_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: FluxFactor

    FluxFactor &
      = MIN( MAX( SQRT( Gm_dd_11 * I_1**2 &
                        + Gm_dd_22 * I_2**2 &
                        + Gm_dd_33 * I_3**2 ) &
                    / MAX( D, SqrtTiny ), &
                  SqrtTiny ), &
             One )

    RETURN
  END FUNCTION FluxFactor_Scalar

  FUNCTION FluxFactor_Vector &
    ( D, I_1, I_2, I_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
    RESULT( FluxFactor )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D(:), I_1(:), I_2(:), I_3(:)
    REAL(DP), INTENT(in) :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)
    REAL(DP) :: FluxFactor(SIZE(D))

    FluxFactor &
      = MIN( MAX( SQRT( Gm_dd_11 * I_1**2 &
                        + Gm_dd_22 * I_2**2 &
                        + Gm_dd_33 * I_3**2 ) &
                    / MAX( D, SqrtTiny ), &
                  SqrtTiny ), &
             One )

    RETURN
  END FUNCTION FluxFactor_Vector


  ! --- Eddington Factor ---


  FUNCTION EddingtonFactor_Scalar( D, FF ) &
      RESULT( EddingtonFactor )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D, FF
    REAL(DP) :: EddingtonFactor

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
  END FUNCTION EddingtonFactor_Scalar

  FUNCTION EddingtonFactor_Vector( D, FF ) &
      RESULT( EddingtonFactor )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D(:), FF(:)
    REAL(DP) :: EddingtonFactor(SIZE(D))

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
  END FUNCTION EddingtonFactor_Vector


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
                + ClosurePolynomialDerivative_ME_CB( XX ) &
                    * XX * ( One - Two * D ) )

    dEFdFF_D &
      = Two * Third * ( One - Two * D ) &
          * ClosurePolynomialDerivative_ME_CB( XX )

#elif MOMENT_CLOSURE_MAXIMUM_ENTROPY_BL

    ! --- Banach-Larecki ME Closure ---

    dEFdD_FF &
      = Two * Third &
          * ( ClosurePolynomial_ME_BL( XX ) * ( Four * D - Three ) &
              + ClosurePolynomialDerivative_ME_BL( XX ) &
                  * XX * ( One - Two * D ) )

    dEFdFF_D &
      = Two * Third * ( One - Two * D ) &
          * ClosurePolynomialDerivative_ME_BL( XX )

#elif MOMENT_CLOSURE_KERSHAW_BL

    ! --- Banach-Larecki Kershaw Closure ---
 
    dEFdD_FF &
      = Two * Third &
          * ( ClosurePolynomial_KE_BL( XX ) * ( Four * D - Three ) &
              + ClosurePolynomialDerivative_KE_BL( XX ) &
                  * XX * ( One - Two * D ) )

    dEFdFF_D &
      = Two * Third * ( One - Two * D ) &
          * ClosurePolynomialDerivative_KE_BL( XX )

#endif

  END SUBROUTINE ComputeEddingtonFactorDerivatives


  ! --- Closure Polynomials ---


  FUNCTION ClosurePolynomial_ME_CB_Scalar( X ) &
      RESULT( ClosurePolynomial_ME_CB )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Cernohorsky-Bludman Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X
    REAL(DP) :: ClosurePolynomial_ME_CB

    ClosurePolynomial_ME_CB &
      = Three * Fifth * X**2 * ( One - Third * X + X**2 )

  END FUNCTION ClosurePolynomial_ME_CB_Scalar

  FUNCTION ClosurePolynomial_ME_CB_Vector( X ) &
      RESULT( ClosurePolynomial_ME_CB )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Cernohorsky-Bludman Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X(:)
    REAL(DP) :: ClosurePolynomial_ME_CB(SIZE(X))

    ClosurePolynomial_ME_CB &
      = Three * Fifth * X**2 * ( One - Third * X + X**2 )

  END FUNCTION ClosurePolynomial_ME_CB_Vector


  FUNCTION ClosurePolynomial_ME_BL_Scalar( X ) &
      RESULT( ClosurePolynomial_ME_BL )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Banach-Larecki Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X
    REAL(DP) :: ClosurePolynomial_ME_BL

    ClosurePolynomial_ME_BL &
      = ( 9.0_DP * X**2 - 5.0_DP &
          + SQRT( 33.0_DP * X**4 - 42.0_DP * X**2 + 25.0_DP ) ) / 8.0_DP

  END FUNCTION ClosurePolynomial_ME_BL_Scalar

  FUNCTION ClosurePolynomial_ME_BL_Vector( X ) &
      RESULT( ClosurePolynomial_ME_BL )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Banach-Larecki Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X(:)
    REAL(DP) :: ClosurePolynomial_ME_BL(SIZE(X))

    ClosurePolynomial_ME_BL &
      = ( 9.0_DP * X**2 - 5.0_DP &
          + SQRT( 33.0_DP * X**4 - 42.0_DP * X**2 + 25.0_DP ) ) / 8.0_DP

  END FUNCTION ClosurePolynomial_ME_BL_Vector


  FUNCTION ClosurePolynomial_KE_BL_Scalar( X ) &
      RESULT( ClosurePolynomial_KE_BL )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Banach-Larecki Kershaw Closure ---

    REAL(DP), INTENT(in) :: X
    REAL(DP) :: ClosurePolynomial_KE_BL

    ClosurePolynomial_KE_BL &
      = X**2

  END FUNCTION ClosurePolynomial_KE_BL_Scalar

  FUNCTION ClosurePolynomial_KE_BL_Vector( X ) &
      RESULT( ClosurePolynomial_KE_BL )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Banach-Larecki Kershaw Closure ---

    REAL(DP), INTENT(in) :: X(:)
    REAL(DP) :: ClosurePolynomial_KE_BL(SIZE(X))

    ClosurePolynomial_KE_BL &
      = X**2

  END FUNCTION ClosurePolynomial_KE_BL_Vector


  ! --- Derivative of Closure Polynomials ---


  FUNCTION ClosurePolynomialDerivative_ME_CB_Scalar( X ) &
      RESULT( ClosurePolynomialDerivative_ME_CB )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Cernohorsky-Bludman Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X
    REAL(DP) :: ClosurePolynomialDerivative_ME_CB

    ClosurePolynomialDerivative_ME_CB &
      = Three * Fifth * X * ( Two - X + Four * X * X )

  END FUNCTION ClosurePolynomialDerivative_ME_CB_Scalar


  FUNCTION ClosurePolynomialDerivative_ME_CB_Vector( X ) &
      RESULT( ClosurePolynomialDerivative_ME_CB )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Cernohorsky-Bludman Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X(:)
    REAL(DP) :: ClosurePolynomialDerivative_ME_CB(SIZE(X))

    ClosurePolynomialDerivative_ME_CB &
      = Three * Fifth * X * ( Two - X + Four * X**2 )

  END FUNCTION ClosurePolynomialDerivative_ME_CB_Vector


  FUNCTION ClosurePolynomialDerivative_ME_BL_Scalar( X ) &
      RESULT( ClosurePolynomialDerivative_ME_BL )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Banach-Larecki Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X
    REAL(DP) :: aa
    REAL(DP) :: ClosurePolynomialDerivative_ME_BL

    aa = SQRT( 33.0_DP * X**4 - 42.0_DP * X * X + 25.0_DP )

    ClosurePolynomialDerivative_ME_BL &
      = Three * x * ( 11.0_DP * X * X + Three * X - 7.0_DP) / ( Four * aa )

  END FUNCTION ClosurePolynomialDerivative_ME_BL_Scalar


  FUNCTION ClosurePolynomialDerivative_ME_BL_Vector( X ) &
      RESULT( ClosurePolynomialDerivative_ME_BL )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Banach-Larecki Maximum Entropy Closure ---

    REAL(DP), INTENT(in) :: X(:)
    REAL(DP) :: aa(SIZE(X))
    REAL(DP) :: ClosurePolynomialDerivative_ME_BL(SIZE(X))

    aa = SQRT( 33.0_DP * X**4 - 42.0_DP * X * X + 25.0_DP )

    ClosurePolynomialDerivative_ME_BL &
      = Three * x * ( 11.0_DP * X * X + Three * X - 7.0_DP) / ( Four * aa )

  END FUNCTION ClosurePolynomialDerivative_ME_BL_Vector


  FUNCTION ClosurePolynomialDerivative_KE_BL_Scalar( X ) &
      RESULT( ClosurePolynomialDerivative_KE_BL )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Banach-Larecki Kershaw Closure ---

    REAL(DP), INTENT(in) :: X
    REAL(DP) :: ClosurePolynomialDerivative_KE_BL

    ClosurePolynomialDerivative_KE_BL &
      = Two * X

  END FUNCTION ClosurePolynomialDerivative_KE_BL_Scalar


  FUNCTION ClosurePolynomialDerivative_KE_BL_Vector( X ) &
      RESULT( ClosurePolynomialDerivative_KE_BL )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Banach-Larecki Kershaw Closure ---

    REAL(DP), INTENT(in) :: X(:)
    REAL(DP) :: ClosurePolynomialDerivative_KE_BL(SIZE(X))

    ClosurePolynomialDerivative_KE_BL &
      = Two * X

  END FUNCTION ClosurePolynomialDerivative_KE_BL_Vector


END MODULE TwoMoment_ClosureModule
