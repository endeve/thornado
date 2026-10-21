#if defined(__FAST_MATH__)
#error Retained arithmetic must be compiled without fast-math reassociation.
#endif
MODULE TwoMoment_RetainedArithmeticModule

  USE KindModule, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  ! Two binary64 components, not two independent nonlinear unknowns.
  TYPE, PUBLIC :: RetainedReal
    REAL(DP) :: High, Low
  END TYPE RetainedReal

  PUBLIC :: Retain, RetainedValue, RetainedAbs, RetainedSqrt
  PUBLIC :: RetainedExp, RetainedFermi, RetainedAccumulate
  PUBLIC :: OPERATOR(+), OPERATOR(-), OPERATOR(*), OPERATOR(/)
  PUBLIC :: ASSIGNMENT(=)

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE AssignReal
  END INTERFACE
  INTERFACE OPERATOR(+)
    MODULE PROCEDURE AddPair, AddReal, RealAdd
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE SubtractPair, SubtractReal, RealSubtract, NegatePair
  END INTERFACE
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE MultiplyPair, MultiplyReal, RealMultiply
  END INTERFACE
  INTERFACE OPERATOR(/)
    MODULE PROCEDURE DividePair, DivideReal, RealDivide
  END INTERFACE

CONTAINS

  ELEMENTAL SUBROUTINE AssignReal( A, B )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(out) :: A
    REAL(DP), INTENT(in) :: B
    A % High = B
    A % Low = 0.0_DP
  END SUBROUTINE AssignReal

  ELEMENTAL FUNCTION Retain( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A
    REAL(DP), INTENT(in), OPTIONAL :: B
    TYPE(RetainedReal) :: C
    REAL(DP) :: X, Y, Z
    Y = 0.0_DP
    IF( PRESENT(B) ) Y = B
    X = A + Y
    Z = X - A
    C % High = X
    C % Low = ( A - ( X - Z ) ) + ( Y - Z )
  END FUNCTION Retain

  ELEMENTAL FUNCTION RetainedValue( A ) RESULT(B)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    REAL(DP) :: B
    B = A % High + A % Low
  END FUNCTION RetainedValue

  ELEMENTAL FUNCTION AddPair( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C
    REAL(DP) :: S, V, E
    S = A % High + B % High
    V = S - A % High
    E = ( A % High - ( S - V ) ) + ( B % High - V )
    C = Retain( S, ( E + A % Low ) + B % Low )
  END FUNCTION AddPair

  ELEMENTAL FUNCTION AddReal( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    REAL(DP), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    C = A + Retain(B)
  END FUNCTION AddReal

  ELEMENTAL FUNCTION RealAdd( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A
    TYPE(RetainedReal), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    C = Retain(A) + B
  END FUNCTION RealAdd

  ELEMENTAL FUNCTION NegatePair( A ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    TYPE(RetainedReal) :: C
    C % High = -A % High
    C % Low = -A % Low
  END FUNCTION NegatePair

  ELEMENTAL FUNCTION SubtractPair( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C
    C = A + (-B)
  END FUNCTION SubtractPair

  ELEMENTAL FUNCTION SubtractReal( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    REAL(DP), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    C = A - Retain(B)
  END FUNCTION SubtractReal

  ELEMENTAL FUNCTION RealSubtract( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A
    TYPE(RetainedReal), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    C = Retain(A) - B
  END FUNCTION RealSubtract

  ELEMENTAL FUNCTION MultiplyPair( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C
    REAL(DP) :: P, AH, AL, BH, BL, CA, CB, E
    ! Match the checked prototype on CPUs. Explicit association is required;
    ! compile this arithmetic without reassociation/unsafe fast-math.
    P = A % High * B % High
    CA = 134217729.0_DP * A % High
    CB = 134217729.0_DP * B % High
    AH = CA - ( CA - A % High )
    BH = CB - ( CB - B % High )
    AL = A % High - AH
    BL = B % High - BH
    E = ( ( AH * BH - P ) + AH * BL + AL * BH ) + AL * BL
    C = Retain( P, E + A % High * B % Low &
                      + A % Low * B % High + A % Low * B % Low )
  END FUNCTION MultiplyPair

  ELEMENTAL FUNCTION MultiplyReal( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    REAL(DP), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    C = A * Retain(B)
  END FUNCTION MultiplyReal

  ELEMENTAL FUNCTION RealMultiply( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A
    TYPE(RetainedReal), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    C = Retain(A) * B
  END FUNCTION RealMultiply

  ELEMENTAL FUNCTION DividePair( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C, R
    REAL(DP) :: Q
    Q = A % High / B % High
    R = A - B * Q
    C = Retain( Q, RetainedValue(R) / B % High )
  END FUNCTION DividePair

  ELEMENTAL FUNCTION DivideReal( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    REAL(DP), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    C = A / Retain(B)
  END FUNCTION DivideReal

  ELEMENTAL FUNCTION RealDivide( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A
    TYPE(RetainedReal), INTENT(in) :: B
    TYPE(RetainedReal) :: C
    C = Retain(A) / B
  END FUNCTION RealDivide

  ELEMENTAL FUNCTION RetainedAbs( A ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    TYPE(RetainedReal) :: C
    C = A
    IF( A % High < 0.0_DP ) C = -A
  END FUNCTION RetainedAbs

  ELEMENTAL FUNCTION RetainedSqrt( A ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    TYPE(RetainedReal) :: C
    INTEGER :: I
    C = SQRT(MAX(0.0_DP,A % High))
    IF( C % High == 0.0_DP ) RETURN
    DO I = 1, 2
      C = C + ( A - C * C ) / ( 2.0_DP * C )
    END DO
  END FUNCTION RetainedSqrt

  ELEMENTAL FUNCTION RetainedExp( A ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    TYPE(RetainedReal) :: C, R, Term, LnTwo
    INTEGER :: I, N
    ! This function is used for nonpositive FD-tail arguments only.
    IF( A % High < -745.0_DP )THEN
      C = 0.0_DP
      RETURN
    END IF
    LnTwo = Retain(LOG(2.0_DP),2.3190468138462995584e-17_DP)
    N = NINT(A % High / LnTwo % High)
    R = A - LnTwo * REAL(N,DP)
    C = 1.0_DP
    Term = 1.0_DP
    DO I = 1, 26
      Term = Term * R / REAL(I,DP)
      C = C + Term
    END DO
    C % High = SCALE(C % High,N)
    C % Low = SCALE(C % Low,N)
  END FUNCTION RetainedExp

  ELEMENTAL FUNCTION RetainedFermi( A ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A
    TYPE(RetainedReal) :: C, Tail
    Tail = RetainedExp(-RetainedAbs(A))
    IF( A % High >= 0.0_DP )THEN
      C = Tail / ( 1.0_DP + Tail )
    ELSE
      C = 1.0_DP / ( 1.0_DP + Tail )
    END IF
  END FUNCTION RetainedFermi

  PURE SUBROUTINE RetainedAccumulate( Sum, Term )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(inout) :: Sum
    TYPE(RetainedReal), INTENT(in) :: Term
    Sum = Sum + Term
  END SUBROUTINE RetainedAccumulate

END MODULE TwoMoment_RetainedArithmeticModule
