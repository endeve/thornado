MODULE TwoMoment_ClosureModule

  USE KindModule, ONLY: &
    DP, SqrtTiny, &
    Zero, One, Two, Three, Four, &
    Fifth, Third
  USE, INTRINSIC :: ieee_arithmetic, ONLY: IEEE_IS_NAN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeClosure_TwoMoment
  PUBLIC :: FluxFactor
  PUBLIC :: FluxFactor_Relativistic
  PUBLIC :: EddingtonFactor
  PUBLIC :: HeatFluxFactor
  PUBLIC :: ComputeEddingtonFactorDerivatives

  INTERFACE FluxFactor
    MODULE PROCEDURE FluxFactor_Scalar
    MODULE PROCEDURE FluxFactor_Vector
  END INTERFACE FluxFactor

  INTERFACE FluxFactor_Relativistic
    MODULE PROCEDURE FluxFactor_Relativistic_Scalar
    MODULE PROCEDURE FluxFactor_Relativistic_Vector
  END INTERFACE FluxFactor_Relativistic

  INTERFACE EddingtonFactor
    MODULE PROCEDURE EddingtonFactor_Scalar
    MODULE PROCEDURE EddingtonFactor_Vector
  END INTERFACE EddingtonFactor

  INTERFACE HeatFluxFactor
    MODULE PROCEDURE HeatFluxFactor_Scalar
    MODULE PROCEDURE HeatFluxFactor_Vector
  END INTERFACE HeatFluxFactor

  INTERFACE ClosurePolynomial_ME_CB
    MODULE PROCEDURE ClosurePolynomial_ME_CB_Scalar
    MODULE PROCEDURE ClosurePolynomial_ME_CB_Vector
  END INTERFACE ClosurePolynomial_ME_CB

  INTERFACE ClosurePolynomial_ME_BL
    MODULE PROCEDURE ClosurePolynomial_ME_BL_Scalar
    MODULE PROCEDURE ClosurePolynomial_ME_BL_Vector
  END INTERFACE ClosurePolynomial_ME_BL

  INTERFACE ClosurePolynomial_KE_BL
    MODULE PROCEDURE ClosurePolynomial_KE_BL_Scalar
    MODULE PROCEDURE ClosurePolynomial_KE_BL_Vector
  END INTERFACE ClosurePolynomial_KE_BL

  INTERFACE ClosurePolynomialDerivative_ME_CB
    MODULE PROCEDURE ClosurePolynomialDerivative_ME_CB_Scalar
    MODULE PROCEDURE ClosurePolynomialDerivative_ME_CB_Vector
  END INTERFACE ClosurePolynomialDerivative_ME_CB

  INTERFACE ClosurePolynomialDerivative_ME_BL
    MODULE PROCEDURE ClosurePolynomialDerivative_ME_BL_Scalar
    MODULE PROCEDURE ClosurePolynomialDerivative_ME_BL_Vector
  END INTERFACE ClosurePolynomialDerivative_ME_BL

  INTERFACE ClosurePolynomialDerivative_KE_BL
    MODULE PROCEDURE ClosurePolynomialDerivative_KE_BL_Scalar
    MODULE PROCEDURE ClosurePolynomialDerivative_KE_BL_Vector
  END INTERFACE ClosurePolynomialDerivative_KE_BL

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

#elif  MOMENT_CLOSURE_KERSHAW

    ! --- Classical Kershaw Closure ---

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A6,A)') &
        '', 'Two-Moment Closure: Kershaw (Classical)'

    END IF

#elif  MOMENT_CLOSURE_KERSHAW_BL

    ! --- Banach-Larecki Kershaw Closure ---

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A6,A)') &
        '', 'Two-Moment Closure: Kershaw (Banach & Larecki)'

    END IF

#elif  MOMENT_CLOSURE_LEVERMORE

     ! --- Levermore Closure ---

     IF( Verbose )THEN

       WRITE(*,*)
       WRITE(*,'(A6,A)') &
         '', 'Two-Moment Closure: Levermore'
      
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


  FUNCTION FluxFactor_Relativistic_Scalar &
    ( D, I_1, I_2, I_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, v_u_1, v_u_2, v_u_3 ) &
    RESULT( FluxFactor_Relativistic )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D, I_1, I_2, I_3, alp, B_u_1, B_u_2, B_u_3, v_u_1, v_u_2, v_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: FluxFactor_Relativistic 
    REAL(DP) :: v_d_0, v_d_1, v_d_2, v_d_3, I_u_0, B_d_1, B_d_2, B_d_3, Gm_dd_00, Isq

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3
    
    v_d_1 = Gm_dd_11 * v_u_1
    v_d_2 = Gm_dd_22 * v_u_2
    v_d_3 = Gm_dd_33 * v_u_3

    v_d_0 = B_d_1 * v_u_1 +  B_d_2 * v_u_2 + B_d_3 * v_u_3
    I_u_0 = ( v_d_1 * I_1 + v_d_2 * I_2 + v_d_3 * I_3 ) / alp
    I_u_0 = I_u_0 * ( 1.0_DP / (1.0_DP - v_d_0 / alp ) ) 

    Gm_dd_00 = -alp**2 + B_d_1 * B_u_1 +  B_d_2 * B_u_2 + B_d_3 * B_u_3

    Isq = Gm_dd_00 * I_u_0 * I_u_0
    Isq = Isq + 2 * ( B_d_1 * I_1 * I_u_0 + B_d_2 * I_2 * I_u_0 + B_d_3 * I_3 * I_u_0 )
    Isq = Isq + Gm_dd_11 * I_1**2 +  Gm_dd_22 * I_2**2 + Gm_dd_33 * I_3**2 

    FluxFactor_Relativistic &
      = MIN( MAX( SQRT( Isq ) &
                    / MAX( D, SqrtTiny ), &
                  SqrtTiny ), &
             One )

    RETURN
  END FUNCTION FluxFactor_Relativistic_Scalar


  FUNCTION FluxFactor_Relativistic_Vector &
    ( D, I_1, I_2, I_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, v_u_1, v_u_2, v_u_3 ) &
    RESULT( FluxFactor_Relativistic )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D(:), I_1(:), I_2(:), I_3(:) 
    REAL(DP), INTENT(in) :: alp(:), B_u_1(:), B_u_2(:), B_u_3(:), v_u_1(:), v_u_2(:), v_u_3(:)
    REAL(DP), INTENT(in) :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)
    REAL(DP) :: FluxFactor_Relativistic(SIZE(D))
    REAL(DP) :: v_d_0(SIZE(D)), v_d_1(SIZE(D)), v_d_2(SIZE(D)), v_d_3(SIZE(D)), I_u_0(SIZE(D))
    REAL(DP) :: B_d_1(SIZE(D)), B_d_2(SIZE(D)), B_d_3(SIZE(D)), Gm_dd_00(SIZE(D)), Isq(SIZE(D))

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3
    
    v_d_1 = Gm_dd_11 * v_u_1
    v_d_2 = Gm_dd_22 * v_u_2
    v_d_3 = Gm_dd_33 * v_u_3

    v_d_0 = B_d_1 * v_u_1 +  B_d_2 * v_u_2 + B_d_3 * v_u_3
    I_u_0 = ( v_d_1 * I_1 + v_d_2 * I_2 + v_d_3 * I_3 ) / alp
    I_u_0 = I_u_0 * ( 1.0_DP / (1.0_DP - v_d_0 / alp ) ) 

    Gm_dd_00 = -alp**2 + B_d_1 * B_u_1 +  B_d_2 * B_u_2 + B_d_3 * B_u_3

    Isq = Gm_dd_00 * I_u_0 * I_u_0
    Isq = Isq + 2 * ( B_d_1 * I_1 * I_u_0 + B_d_2 * I_2 * I_u_0 + B_d_3 * I_3 * I_u_0 )
    Isq = Isq + Gm_dd_11 * I_1**2 +  Gm_dd_22 * I_2**2 + Gm_dd_33 * I_3**2 

    FluxFactor_Relativistic &
      = MIN( MAX( SQRT( Isq ) &
                    / MAX( D, SqrtTiny ), &
                  SqrtTiny ), &
             One )

    RETURN
  END FUNCTION FluxFactor_Relativistic_Vector


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

#ifdef THORNADO_DEBUG
    IF( IEEE_IS_NAN(D) ) STOP 'NAN in D when call EddingtonFactor_Scalar'
#endif

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

#elif  MOMENT_CLOSURE_KERSHAW

    ! --- Classical Kershaw Closure ---

    EddingtonFactor = Third * ( One + Two * FF**2 )

#elif  MOMENT_CLOSURE_KERSHAW_BL

    ! --- Banach-Larecki Kershaw Closure ---

    EddingtonFactor &
      = Third + Two * Third * ( One - D ) * ( One - Two * D ) &
          * ClosurePolynomial_KE_BL( FF / MAX( One - D, SqrtTiny ) )

#elif  MOMENT_CLOSURE_LEVERMORE

    ! --- Levermore Closure ---
  
    EddingtonFactor &
      = Third * ( 5.0_dp - Two * SQRT( Four - Three * FF * FF ) )

#else

    EddingtonFactor = 0.0_DP

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

#elif  MOMENT_CLOSURE_KERSHAW

    ! --- Classical Kershaw Closure ---

    EddingtonFactor = Third * ( One + Two * FF**2 )

#elif  MOMENT_CLOSURE_KERSHAW_BL

    ! --- Banach-Larecki Kershaw Closure ---

    EddingtonFactor &
      = Third + Two * Third * ( One - D ) * ( One - Two * D ) &
          * ClosurePolynomial_KE_BL( FF / MAX( One - D, SqrtTiny ) )

#elif  MOMENT_CLOSURE_LEVERMORE

    ! --- Levermore Closure ---

    EddingtonFactor &
      = Third * ( 5.0_dp - Two * SQRT( Four - Three * FF * FF ) )

#else

    EddingtonFactor = 0.0_DP

#endif

    RETURN
  END FUNCTION EddingtonFactor_Vector


  ! --- Heat Flux Factor ---


  FUNCTION HeatFluxFactor_Scalar( D, FF ) RESULT( HeatFluxFactor )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D, FF
    REAL(DP)             :: HeatFluxFactor

    REAL(DP) :: x

#ifdef MOMENT_CLOSURE_MINERBO

    ! --- Maximum Entropy (ME) Minerbo Closure -------------------
    ! --- Expression from Just et al. (2015), MNRAS, 453, 3386 ---

    HeatFluxFactor &
      = ( 45.0_DP + 10.0_DP * FF - 12.0 * FF**2 - 12.0_DP * FF**3 &
          + 38.0_DP * FF**4 - 12.0_DP * FF**5 + 18.0_DP * FF**6 ) &
        * FF / 75.0_DP

#elif  MOMENT_CLOSURE_MAXIMUM_ENTROPY_CB

    ! --- Expression derived from Richers (2020), PRD, 102 -------
    x = FF / ( One - D )
    HeatFluxFactor &
      = ( One - D ) &
        * ( ( ( One - Two * D + Two * D**2) - Three / 5.0_DP ) * x**6 &
            + Three * x / 5.0_DP )

#elif  MOMENT_CLOSURE_MAXIMUM_ENTROPY_BL

    HeatFluxFactor = 0.0_DP

#elif  MOMENT_CLOSURE_KERSHAW

    ! --- Classical Kershaw Closure ---

    HeatFluxFactor = FF * ( 4.0_DP * FF**2 + 5.0_DP ) / 9.0_DP

#elif  MOMENT_CLOSURE_KERSHAW_BL
    
    HeatFluxFactor &
    =((D + 1.0_DP) / (24.0_DP*D))*(((D*FF+(1.0_DP-D)**2)/(1.0_DP-D))**4 &
    -((D*FF-(1.0_DP-D)**2)/(1.0_DP-D))**4)+((2.0_DP-D)/3)*(FF**3+D**2*FF)

#elif  MOMENT_CLOSURE_LEVERMORE

    HeatFluxFactor = 0.0_DP

#else

    HeatFluxFactor = 0.0_DP

#endif

    RETURN
  END FUNCTION HeatFluxFactor_Scalar


  FUNCTION HeatFluxFactor_Vector( D, FF ) RESULT( HeatFluxFactor )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D(:), FF(:)
    REAL(DP)             :: HeatFluxFactor(SIZE(D))
    
    INTEGER :: i

#ifdef MOMENT_CLOSURE_MINERBO

    ! --- Maximum Entropy (ME) Minerbo Closure -------------------
    ! --- Expression from Just et al. (2015), MNRAS, 453, 3386 ---
    
    HeatFluxFactor &
      = ( 45.0_DP + 10.0_DP * FF - 12.0 * FF**2 - 12.0_DP * FF**3 &
          + 38.0_DP * FF**4 - 12.0_DP * FF**5 + 18.0_DP * FF**6 ) &
        * FF / 75.0_DP

#elif  MOMENT_CLOSURE_MAXIMUM_ENTROPY_CB

    DO i = 1, SIZE( D )
        HeatFluxFactor(i) = HeatFluxFactor_Scalar( D(i), FF(i) )
    END DO
#elif  MOMENT_CLOSURE_MAXIMUM_ENTROPY_BL

    HeatFluxFactor = 0.0_DP

#elif  MOMENT_CLOSURE_KERSHAW

    ! --- Classical Kershaw Closure ---

    HeatFluxFactor = FF * ( 4.0_DP * FF**2 + 5.0_DP ) / 9.0_DP

#elif  MOMENT_CLOSURE_KERSHAW_BL
    
    DO i = 1, SIZE( D )
        HeatFluxFactor(i) = HeatFluxFactor_Scalar( D(i), FF(i) )
    END DO
#elif  MOMENT_CLOSURE_LEVERMORE

    HeatFluxFactor = 0.0_DP

#else

    HeatFluxFactor = 0.0_DP

#endif

    RETURN
  END FUNCTION HeatFluxFactor_Vector


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

#elif  MOMENT_CLOSURE_KERSHAW

    ! --- Classical Kershaw Closure ---

    dEFdD_FF = Zero

    dEFdFF_D = Four * Third * FF

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

#elif MOMENT_CLOSURE_LEVERMORE

    ! --- Levermore Closure ---

    dEFdD_FF &
      = ZERO

    dEFdFF_D &
      = Two * FF / SQRT( Four - Three * FF * FF )

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
