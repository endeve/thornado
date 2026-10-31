#if defined(__FAST_MATH__)
#error Retained arithmetic must be compiled without fast-math reassociation.
#endif
MODULE TwoMoment_RetainedArithmeticModule

  USE KindModule, ONLY: DP
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_FMA

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

    IF( .NOT. PRESENT(B) )THEN
      C % High = A+0.0_DP
      C % Low = 0.0_DP
      RETURN
    END IF
    Y = B
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
    REAL(DP) :: S, V, E
    S = A % High+B
    V = S-A % High
    E = (A % High-(S-V))+(B-V)
    C = Retain(S,E+A % Low)
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
    REAL(DP) :: S, V, E
    S = A+B % High
    V = S-A
    E = (A-(S-V))+(B % High-V)
    C = Retain(S,E+B % Low)
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
    REAL(DP) :: S, V, E
    S = A % High-B % High
    V = S-A % High
    E = (A % High-(S-V))-(B % High+V)
    C = Retain(S,(E+A % Low)-B % Low)
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
    REAL(DP) :: S, V, E
    S = A % High-B
    V = S-A % High
    E = (A % High-(S-V))-(B+V)
    C = Retain(S,E+A % Low)
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
    REAL(DP) :: S, V, E
    S = A-B % High
    V = S-A
    E = (A-(S-V))-(B % High+V)
    C = Retain(S,E-B % Low)
  END FUNCTION RealSubtract

  ELEMENTAL FUNCTION MultiplyPair( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(RetainedReal), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C
    REAL(DP) :: P,E
    ! Match the checked prototype on CPUs. Explicit association is required;
    ! compile this arithmetic without reassociation/unsafe fast-math.
    P = A % High * B % High
    E = IEEE_FMA(A % High,B % High,-P)
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
    REAL(DP) :: P,E
    P = A % High*B
    E = IEEE_FMA(A % High,B,-P)
    C = Retain(P,E+A % Low*B)
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
    REAL(DP) :: P,E
    P = A*B % High
    E = IEEE_FMA(A,B % High,-P)
    C = Retain(P,E+A*B % Low)
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
    TYPE(RetainedReal) :: C, R, P
    REAL(DP) :: Q
    Q = A % High/B
    P = RetainedProduct(B,Q)
    R = A-P
    C = Retain(Q,RetainedValue(R)/B)
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
    TYPE(RetainedReal) :: C, R
    TYPE(RetainedReal), PARAMETER :: LnStep = RetainedReal(4.33216987849965804E-02_DP,1.44940425865393722E-18_DP)
    TYPE(RetainedReal), PARAMETER :: Coefficient(0:13) = [ &
      RetainedReal(1.00000000000000000E+00_DP,0.00000000000000000E+00_DP), &
      RetainedReal(1.00000000000000000E+00_DP,0.00000000000000000E+00_DP), &
      RetainedReal(5.00000000000000000E-01_DP,0.00000000000000000E+00_DP), &
      RetainedReal(1.66666666666666657E-01_DP,9.25185853854297066E-18_DP), &
      RetainedReal(4.16666666666666644E-02_DP,2.31296463463574266E-18_DP), &
      RetainedReal(8.33333333333333322E-03_DP,1.15648231731787138E-19_DP), &
      RetainedReal(1.38888888888888894E-03_DP,-5.30054395437357706E-20_DP), &
      RetainedReal(1.98412698412698413E-04_DP,1.72095582934207053E-22_DP), &
      RetainedReal(2.48015873015873016E-05_DP,2.15119478667758816E-23_DP), &
      RetainedReal(2.75573192239858925E-06_DP,-1.85839327404647208E-22_DP), &
      RetainedReal(2.75573192239858883E-07_DP,2.37677146222502973E-23_DP), &
      RetainedReal(2.50521083854417202E-08_DP,-1.44881407093591197E-24_DP), &
      RetainedReal(2.08767569878681002E-09_DP,-1.20734505911325997E-25_DP), &
      RetainedReal(1.60590438368216133E-10_DP,1.25852945887520981E-26_DP) ]
    TYPE(RetainedReal), PARAMETER :: FractionalPower(0:15) = [ &
      RetainedReal(1.00000000000000000E+00_DP,0.00000000000000000E+00_DP), &
      RetainedReal(1.04427378242741375E+00_DP,8.55188970553796489E-17_DP), &
      RetainedReal(1.09050773266525769E+00_DP,-3.04678207981247115E-17_DP), &
      RetainedReal(1.13878863475669156E+00_DP,8.91281267602540778E-17_DP), &
      RetainedReal(1.18920711500272103E+00_DP,3.98201523146564611E-17_DP), &
      RetainedReal(1.24185781207348400E+00_DP,4.65802759183693679E-17_DP), &
      RetainedReal(1.29683955465100964E+00_DP,2.53825027948883150E-17_DP), &
      RetainedReal(1.35425554693689265E+00_DP,7.70094837980298946E-17_DP), &
      RetainedReal(1.41421356237309515E+00_DP,-9.66729331345291345E-17_DP), &
      RetainedReal(1.47682614593949935E+00_DP,-3.48399455689279580E-17_DP), &
      RetainedReal(1.54221082540794074E+00_DP,7.94983480969762086E-17_DP), &
      RetainedReal(1.61049033194925428E+00_DP,2.47071925697978879E-17_DP), &
      RetainedReal(1.68179283050742900E+00_DP,8.19901002058149652E-17_DP), &
      RetainedReal(1.75625216037329945E+00_DP,2.96014069544887331E-17_DP), &
      RetainedReal(1.83400808640934243E+00_DP,3.28310722424562720E-17_DP), &
      RetainedReal(1.91520656139714740E+00_DP,-1.06199460561959626E-16_DP) ]
    INTEGER :: I, N, J, Power

    IF( A % High < -745.0_DP )THEN
      C = 0.0_DP
      RETURN
    END IF
    ! exp(A) = 2**Power * 2**(J/16) * exp(R), |R| <= log(2)/32.
    ! The degree-13 Taylor remainder is below 5.87e-35.
    N = NINT(A % High/LnStep % High)
    J = MODULO(N,16)
    Power = (N-J)/16
    R = A-LnStep*REAL(N,DP)
    C = Coefficient(13)
    DO I = 12, 0, -1
      C = Coefficient(I)+R*C
    END DO
    C = C*FractionalPower(J)
    C % High = SCALE(C % High,Power)
    C % Low = SCALE(C % Low,Power)
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

  ELEMENTAL FUNCTION RetainedProduct( A, B ) RESULT(C)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: A, B
    TYPE(RetainedReal) :: C
    C % High = A*B
    C % Low = IEEE_FMA(A,B,-C % High)
  END FUNCTION RetainedProduct

END MODULE TwoMoment_RetainedArithmeticModule
