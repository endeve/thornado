MODULE UtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nNodesE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Locate
  PUBLIC :: NodeNumber
  PUBLIC :: NodeNumberX
  PUBLIC :: NodeNumberX_X1
  PUBLIC :: InitializeWeights
  PUBLIC :: Interpolate1D_Linear
  PUBLIC :: Interpolate1D_Log
  PUBLIC :: MapTo1D
  PUBLIC :: MapFrom1D
  PUBLIC :: GetRoots_Quadratic
  PUBLIC :: MinModB
  PUBLIC :: WriteVector
  PUBLIC :: WriteMatrix

  INTERFACE InitializeWeights
    MODULE PROCEDURE InitializeWeights3D
    MODULE PROCEDURE InitializeWeights4D
  END INTERFACE InitializeWeights

  INTERFACE MapTo1D
    MODULE PROCEDURE Map4DTo1D
  END INTERFACE MapTo1D

  INTERFACE MapFrom1D
    MODULE PROCEDURE Map1DTo4D
  END INTERFACE MapFrom1D

CONTAINS


  PURE INTEGER FUNCTION Locate( x, xx, n )

    REAL(DP), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    INTEGER :: il, im, iu

    il = 0
    iu = n+1
    DO WHILE ( iu - il > 1 )
      im = (iu+il)/2
      IF ((xx(n).ge.xx(1)).eqv.(x.ge.xx(im))) THEN
        il = im
      ELSE
        iu = im
      END IF
    END DO

    IF (x.eq.xx(1)) THEN
      Locate = 1
    ELSEIF (x.eq.xx(n)) THEN
      Locate = n-1
    ELSE
      Locate = il
    END IF

    RETURN
  END FUNCTION Locate


  PURE INTEGER FUNCTION NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

    INTEGER, INTENT(in) :: iNodeE, iNodeX1, iNodeX2, iNodeX3

    NodeNumber &
      = iNodeE &
        + ( (iNodeX1-1)  &
            + ( (iNodeX2-1)  &
                + (iNodeX3-1) * nNodesX(2) ) * nNodesX(1) ) * nNodesE

    RETURN
  END FUNCTION NodeNumber


  PURE INTEGER FUNCTION NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

    INTEGER, INTENT(in) :: iNodeX1, iNodeX2, iNodeX3

    NodeNumberX &
      = iNodeX1 &
        + nNodesX(1) * ( (iNodeX2-1) + (iNodeX3-1) * nNodesX(2) )

    RETURN
  END FUNCTION NodeNumberX


  PURE INTEGER FUNCTION NodeNumberX_X1( iNodeX2, iNodeX3 )

    INTEGER, INTENT(in) :: iNodeX2, iNodeX3

    NodeNumberX_X1 &
      = iNodeX2 + (iNodeX3-1) * nNodesX(2)

    RETURN
  END FUNCTION NodeNumberX_X1


  SUBROUTINE InitializeWeights3D( w1, w2, w3, w3D, w2D_1, w2D_2, w2D_3 )

    REAL(DP), DIMENSION(:), INTENT(in)  :: w1, w2, w3
    REAL(DP), DIMENSION(:), INTENT(out) :: w3D
    REAL(DP), DIMENSION(:), INTENT(out) :: w2D_1
    REAL(DP), DIMENSION(:), INTENT(out) :: w2D_2
    REAL(DP), DIMENSION(:), INTENT(out) :: w2D_3

    INTEGER :: i, i1, i2, i3

    DO i3 = 1, SIZE( w3 )
      DO i2 = 1, SIZE( w2 )
        DO i1 = 1, SIZE( w1 )

          i = NodeNumberX( i1, i2, i3 )

          w3D(i) = w1(i1) * w2(i2) * w3(i3)

        END DO
      END DO
    END DO

    DO i3 = 1, SIZE( w3 )
      DO i2 = 1, SIZE( w2 )

        i = NodeNumberX_X1( i2, i3 )

        w2D_1(i) = w2(i2) * w3(i3)

      END DO
    END DO

  END SUBROUTINE InitializeWeights3D


  SUBROUTINE InitializeWeights4D( w1, w2, w3, w4, w4D )

    REAL(DP), DIMENSION(:), INTENT(in)  :: w1, w2, w3, w4
    REAL(DP), DIMENSION(:), INTENT(out) :: w4D

    INTEGER :: i, i1, i2, i3, i4

    DO i4 = 1, SIZE( w4 )
      DO i3 = 1, SIZE( w3 )
        DO i2 = 1, SIZE( w2 )
          DO i1 = 1, SIZE( w1 )

            i = NodeNumber( i1, i2, i3, i4 )

            w4D(i) = w1(i1) * w2(i2) * w3(i3) * w4(i4)

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeWeights4D


  PURE REAL(DP) FUNCTION Interpolate1D_Linear( x, xL, xR, uL, uR )

    REAL(DP), INTENT(in) :: x, xL, xR, uL, uR

    Interpolate1D_Linear &
      = ( ( xR - x ) * uL + ( x - xL ) * uR ) / ( xR - xL )

    RETURN
  END FUNCTION Interpolate1D_Linear


  PURE REAL(DP) FUNCTION Interpolate1D_Log( x, xL, xR, uL, uR )

    REAL(DP), INTENT(in) :: x, xL, xR, uL, uR

    Interpolate1D_Log &
      = 10.0_DP**( ( LOG10( xR / x ) * LOG10( uL ) &
                     + LOG10( x / xL ) * LOG10( uR ) ) / LOG10( xR / xL ) )

    RETURN
  END FUNCTION Interpolate1D_Log


  SUBROUTINE Map4DTo1D( u4D, u1D )

    REAL(DP), DIMENSION(:,:,:,:), INTENT(in)  :: u4D
    REAL(DP), DIMENSION(:),       INTENT(out) :: u1D

    INTEGER :: i, i1, i2, i3, i4, n(4)

    n = SHAPE( u4D )

    i = 1
    DO i4 = 1, n(4)
      DO i3 = 1, n(3)
        DO i2 = 1, n(2)
          DO i1 = 1, n(1)

            u1D(i) = u4D(i1,i2,i3,i4)

            i = i + 1

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE Map4DTo1D


  SUBROUTINE Map1DTo4D( u4D, u1D )

    REAL(DP), DIMENSION(:,:,:,:), INTENT(out) :: u4D
    REAL(DP), DIMENSION(:),       INTENT(in)  :: u1D

    INTEGER :: i, i1, i2, i3, i4, n(4)

    n = SHAPE( u4D )

    i = 1
    DO i4 = 1, n(4)
      DO i3 = 1, n(3)
        DO i2 = 1, n(2)
          DO i1 = 1, n(1)

            u4D(i1,i2,i3,i4) = u1D(i)

            i = i + 1

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE Map1DTo4D


  SUBROUTINE GetRoots_Quadratic( a, b, c, r1, r2, Tol_Option )

    REAL(DP), INTENT(in)  :: a, b, c
    REAL(DP), INTENT(out) :: r1, r2
    REAL(DP), INTENT(in), OPTIONAL :: Tol_Option

    REAL(DP) :: delta, rq
    REAL(DP) :: Tol

    Tol = 1.0d-13
    IF( PRESENT( Tol_Option ) ) &
      Tol = Tol_Option

    IF( ABS( a ) <= Tol )THEN
      r1 = - c / b
      r2 = r1
    ELSE
      delta = b**2 - 4.0_DP * a * c
      IF( delta < 0.0_DP )THEN
        r1 = - SQRT( HUGE( 1.0_DP ) )
        r2 = - SQRT( HUGE( 1.0_DP ) )
      ELSE
        rq = - 0.5_DP * ( b + SIGN( 1.0_DP, b ) * SQRT( delta ) )
        r1 = rq / a
        r2 = c / rq
      END IF
    END IF

  END SUBROUTINE GetRoots_Quadratic


  REAL(DP) PURE ELEMENTAL FUNCTION MinMod2( a, b )

    REAL(DP), INTENT(in) :: a, b

    IF( a * b > 0.0_DP )THEN
      IF( ABS( a ) < ABS( b ) )THEN
        MinMod2 = a
      ELSE
        MinMod2 = b
      END IF
    ELSE
      MinMod2 = 0.0_DP
    END IF

    RETURN
  END FUNCTION MinMod2


  REAL(DP) PURE ELEMENTAL FUNCTION MinMod( a, b, c )

    REAL(DP), INTENT(in) :: a, b, c

    MinMod = MinMod2( a, MinMod2( b, c ) )

    RETURN
  END FUNCTION MinMod


  REAL(DP) PURE ELEMENTAL FUNCTION MinModB( a, b, c, dx, M )

    REAL(DP), INTENT(in) :: a, b, c, dx, M

    IF( ABS( a ) < M * dx**2 )THEN

      MinModB = a

    ELSE

      MinModB = MinMod( a, b, c )

    END IF

  END FUNCTION MinModB


  SUBROUTINE WriteVector( N, Vec, FileName )

    INTEGER,                INTENT(in) :: N
    REAL(DP), DIMENSION(N), INTENT(in) :: Vec
    CHARACTER(LEN=*),       INTENT(in) :: FileName

    INTEGER :: FUNIT

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) N, Vec(1:N)

    CLOSE( FUNIT )

  END SUBROUTINE WriteVector


  SUBROUTINE WriteMatrix( M, N, Mat, FileName )

    INTEGER,                  INTENT(in) :: M, N
    REAL(DP), DIMENSION(M,N), INTENT(in) :: Mat
    CHARACTER(LEN=*),         INTENT(in) :: FileName

    INTEGER :: FUNIT

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) M, N, Mat(1:M,1:N)

    CLOSE( FUNIT )

  END SUBROUTINE WriteMatrix


END MODULE UtilitiesModule
