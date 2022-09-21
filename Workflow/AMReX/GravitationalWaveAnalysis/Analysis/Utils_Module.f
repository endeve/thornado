c nodal points and weight for gauss-legendre
c
      SUBROUTINE GAULEG(X1,X2,X,W,N)
      IMPLICIT double precision (A-H,O-Z)
      integer n,I,J,M
      double precision x1,x2,x(N),w(N),xm,xl,z,p1,p2,p3,pp,z1,pi,eps
      PARAMETER (EPS=3.0D-14)

      pi = 4.0d0*datan(1.0d0)
c
      M=(N+1)/2
      XM=5.0D-1*(X2+X1)
      XL=5.0D-1*(X2-X1)
      DO 12 I=1,M
        Z=COS(pi*(I-2.5D-1)/(N+5.0D-1))
1       CONTINUE
          P1=1.0D0
          P2=0.0D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.0D0*J-1.0D0)*Z*P2-(J-1.0D0)*P3)/dble(J)
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.0D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.0D0*XL/((1.0D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END SUBROUTINE GAULEG

************************************************************************
************************************************************************

c      FUNCTION PLGNDR(L,M,X)
c     m=0 legendre
c     normalization is such that integration of (P_l)**2*sin(th) from
c     0 to pi/2 is unity

      double precision FUNCTION PLGNDR(L,th)
      implicit double precision (a-h,o-z)
	integer l,m
	double precision x,th
c	real*8 plgndr,x,th
	integer i,ll
	double precision fact,pll,pmm,pmmp1,somx2

      m = 0
      x = dcos(th)

      IF (M.LT.0.OR.M.GT.L.OR.DABS(X).GT.1.0d0) THEN
        WRITE(*,*) 'bad arguments'
        STOP
      ENDIF
      PMM=1.0d0
      IF(M.GT.0) THEN
        SOMX2=SQRT((1.0d0-X)*(1.0d0+X))
        FACT=1.0d0
        DO 11 I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.0d0
11      CONTINUE
      ENDIF
      IF(L.EQ.M) THEN
        PLGNDR=PMM
      ELSE
        PMMP1=X*dble(2*M+1)*PMM
        IF(L.EQ.M+1) THEN
          PLGNDR=PMMP1
        ELSE
          DO 12 LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          PLGNDR=PLL
        ENDIF
      ENDIF

      plgndr = dsqrt(5.0d-1*(2.0d0*L+1.0d0))*plgndr ! normalization

      RETURN
      END FUNCTION PLGNDR

************************************************************************
************************************************************************

      SUBROUTINE SPLC(X, N, Y, DF, IOPT, C, NC, IER)
************************************************************************
*  COMPUTE THE COEFFICIENTS OF THE CUBIC SPLINE.                       *
*  PARAMETERS                                                          *
*    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              *
*    (2) N: NUMBER OF KNOWN POINTS                                     *
*    (3) Y: 1-DIM. ARRAY FOR FUNCTION'S VALUES ON KNOWN POINTS         *
*    (4) DF: 1-DIM. ARRAY FOR DIFFERENTIALS AT END POINTS              *
*    (5) IOPT: 1-DIM. ARRAY SPECIFYING THE CONTENT OF DF               *
*    (6) C: 2-DIM. WORKING ARRAY                                       *
*    (7) NC: ROW SIZE OF THE ARRAY (C)                                 *
*    (8) IER: ERROR CODE                                               *
*  COPYRIGHT   T. OGUNI   JUNE 30 1989    VERSION 1.0                  *
************************************************************************
       
       INTEGER                        :: N, NC, IER
       double precision, dimension(N) :: X, Y
       double precision, dimension(2) :: DF, D
       INTEGER, dimension(2)          :: IOPT
       double precision, dimension(NC,3) :: C
       double precision, dimension(4) :: EC

       INTEGER :: I, I1, I2, IB, IDER, J, K, II, KS, KE
       double precision :: A1, A2, DY1, DY2, H, H1, H2, HH, HY
       double precision :: PIV, X1, X2, X3, Y1, Y2
C
      IF (N.LT.2 .OR. NC.LT.N-1 .OR. IOPT(1).LT.1 .OR. IOPT(1).GT.3
     *    .OR. IOPT(2).LT.1 .OR. IOPT(2).GT.3) THEN
       IER = 2
       WRITE(*,*) '(SUBR. SPLC) INVALID ARGUMENT.',N,NC,IOPT(1),IOPT(2)
       RETURN
      ENDIF
      DO 5 I=1,N-1
       IF (X(I) .GE. X(I+1)) THEN
        IER = 1
        WRITE(*,*) '(SUBR. SPLC) X SHOULD SATISFY UPWARD ORDER.'
        STOP
       ENDIF
    5 CONTINUE
      IER = 0
C  SET THE END CONDITIONS.
      II = 2
      KS = 1
      KE = MIN0(4,N)
      IDER = 1
      DO 70 I=1,2
       I1 = 2 * I - 1
       I2 = 2 * I
       IB = IOPT(I)
       GO TO (10, 20, 30), IB
   10  EC(I1) = 0.0D0
       EC(I2) = 2.0D0 * DF(I)
       GO TO 70
   20  D(I) = DF(I)
   25  IF (I .EQ. 2) II = N
       H = X(II) - X(II-1)
       EC(I1) = 1.0D0
       HY = Y(II) - Y(II-1)
       EC(I2) = 6.0D0 * (HY / H - D(I)) / H
       IF (I .EQ. 2) EC(I2) = - EC(I2)
       GO TO 70
   30  IF (I .NE. 1) THEN
        KS = MAX0(1,N-3)
        KE = N
        IDER = N
       ENDIF
       A2 = 0.0D0
       D(I) = 0.0D0
       DO 60 K=KS,KE
        IF (IDER .NE. K) THEN
         A1 = 1.0D0
         DO 50 J=KS,KE
          IF (J .NE. IDER .AND. J .NE. K) THEN
           X1 = X(IDER) - X(J)
           X2 = X(K) - X(J)
           A1 = A1 * X1 / X2
          ENDIF
   50    CONTINUE
         X3 = X(K) - X(IDER)
         D(I) = D(I) + A1 * Y(K) / X3
         A2 = A2 - 1.0D0 / X3
        ENDIF
   60  CONTINUE
       D(I) = D(I) + Y(IDER) * A2
       GO TO 25
   70 CONTINUE
C  SET THE ELEMENTS FOR THE SYMMETRIC TRIDIAGONAL EQUATION.
      IF (N .NE. 2) THEN
       H1 = X(2) - X(1)
       Y1 = Y(2) - Y(1)
       DO I=2,N-1
        H2 = X(I+1) - X(I)
        Y2 = Y(I+1) - Y(I)
        HH = H1 + H2
        C(I,1) = H2 / HH
        C(I,2) = 1.0D0 - C(I,1)
        C(I,3) = 6.0D0 * (Y2 / H2 - Y1 / H1) / HH
        H1 = H2
        Y1 = Y2
       END DO
      ENDIF
C  SOLVE THE EQUATION
      C(1,1) = - EC(1) * 0.5D0
      C(1,2) =   EC(2) * 0.5D0
      IF (N .NE. 2) THEN
       DO K=2,N-1
        PIV = 2.0D0 + C(K,2) * C(K-1,1)
        C(K,1) = - C(K,1) / PIV
        C(K,2) = (C(K,3) - C(K,2) * C(K-1,2)) / PIV
       END DO
      ENDIF
      DY1 = (EC(4) - EC(3) * C(N-1,2)) / (2.0D0 + EC(3) * C(N-1,1))
      DO I=1,N-1
       K = N - I
       DY2 = C(K,1) * DY1 + C(K,2)
       H = X(K+1) - X(K)
       C(K,3) = (DY1 - DY2) / (6.0D0 * H)
       C(K,2) = 0.5D0 * DY2
       C(K,1) = (Y(K+1) - Y(K)) / H - (C(K,2) + C(K,3) * H) * H
       DY1 = DY2
      END DO
C
      RETURN
      END SUBROUTINE SPLC

************************************************************************
************************************************************************

      SUBROUTINE SPLD(X, N, C, NC, V, M, D1, D2, IER)
************************************************************************
*  DIFFERENTIATION BY THE CUBIC SPLINE.                                *
*  PARAMETERS                                                          *
*    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              *
*    (2) N: NUMBER OF KNOWN POINTS                                     *
*    (3) C: 2-DIM. WORKING ARRAY                                       *
*    (4) NC: ROW SIZE OF THE ARRAY (C)                                 *
*    (5) V: 1-DIM. ARRAY FOR POINTS WHICH INTERPOLATION MUST BE MADE   *
*    (6) M: NUMBER OF POINTS FOR WHICH INTERPOLATION MUST BE MADE      *
*    (7) D1: 1-DIM. ARRAY FOR FIRST ORDER DIFFERENTIALS                *
*    (8) D2: 1-DIM. ARRAY FOR SECOND ORDER DIFFERENTIALS               *
*    (9) IER: ERROR CODE                                               *
*  COPYRIGHT   T. OGUNI   JUNE 30 1989    VERSION 1.0                  *
************************************************************************

      INTEGER                        :: N, NC, M, IER
      double precision, dimension(N) :: X
      double precision, dimension(M) :: V, D1, D2
      double precision, dimension(NC,3) :: C

      INTEGER                        :: I, K
      double precision               :: T, V1, V2
C
      IF (N .LT. 2 .OR. NC .LT. N-1 .OR. M .LT. 1) THEN
       IER = 2
       WRITE(*,*) '(SUBR. SPLD) INVALID ARGUMENT. ', N, NC, M
       RETURN
      ENDIF
      IER = 0
      I = 1
      DO 90 K=1,M
       V1 = V(K) - X(I)
       IF (V1) 10, 30, 40
   10  IF (I .GT. 1) GO TO 20
       IER = 1
       GO TO 80
   20  I = I - 1
       V1 = V(K) - X(I)
       IF (V1) 10, 30 ,80
   30  D1(K) = C(I,1)
       D2(K) = C(I,2) + C(I,2)
       GO TO 90
   40  IF (I .LT. N) GO TO 50
       IER = 1
       I = N - 1
       GO TO 80
   50  V2 = V(K) - X(I+1)
       IF (V2) 80, 60, 70
   60  IF (I .GE. N-1) GO TO 80
       I = I + 1
       GO TO 30
   70  I = I + 1
       V1 = V2
       GO TO 40
   80  T = C(I,2) + 3.0D0 * C(I,3) * V1
       D1(K) = C(I,1) + (T + C(I,2)) * V1
       D2(K) = T + T
   90 CONTINUE
      RETURN
      END SUBROUTINE SPLD

************************************************************************
************************************************************************

      SUBROUTINE SPLF(X, N, Y, C, NC, V, M, F, IER)
************************************************************************
*  INTERPOLATION BY THE CUBIC SPLINE.                                  *
*  PARAMETERS                                                          *
*    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              *
*    (2) N: NUMBER OF KNOWN POINTS                                     *
*    (3) Y: 1-DIM. ARRAY FOR FUNCTION'S VALUES ON KNOWN POINTS         *
*    (4) C: 2-DIM. WORKING ARRAY                                       *
*    (5) NC: ROW SIZE OF THE ARRAY (C)                                 *
*    (6) V: 1-DIM. ARRAY FOR POINTS WHICH INTERPOLATION MUST BE MADE   *
*    (7) M: NUMBER OF POINTS FOR WHICH INTERPOLATION MUST BE MADE      *
*    (8) F: 1-DIM. WORKING ARRAY                                       *
*    (9) IER: ERROR CODE                                               *
*  COPYRIGHT   T. OGUNI   JUNE 30 1989   VERSION 1.0                   *
************************************************************************

      INTEGER                        :: N, NC, M, IER
      double precision, dimension(N) :: X, Y
      double precision, dimension(M) :: V, F
      double precision, dimension(NC,3) :: C

      INTEGER                        :: I, K
      double precision               :: V1, V2

C
      IF (N .LT. 2 .OR. M .LT. 1 .OR. NC .LT. N-1) THEN
       IER = 2
       WRITE(*,*) '(SUBR. SPLF) INVALID ARGUMENT. ', N, NC, M
       RETURN
      ENDIF
      IER = 0
      I = 1
      DO 90 K=1,M
       V1 = V(K) - X(I)
       IF (V1) 10, 30, 40
   10  IF (I .GT. 1) GO TO 20
       IER = 1
       GO TO 80
   20  I = I - 1
       V1 = V(K) - X(I)
       IF (V1) 10, 30, 80
   30  F(K) = Y(I)
       GO TO 90
   40  IF (I .LT. N) GO TO 50
       IER = 1
       I = N - 1
       GO TO 80
   50  V2 = V(K) - X(I+1)
       IF (V2) 80, 60, 70
   60  I = I + 1
       GO TO 30
   70  I = I + 1
       V1 = V2
       GO TO 40
   80  F(K) = Y(I) + V1 * (C(I,1) + V1 * (C(I,2) + V1 * C(I,3)))
   90 CONTINUE
C
      RETURN
      END SUBROUTINE SPLF

************************************************************************
************************************************************************

      SUBROUTINE SPLQ(X, N, Y, C, NC, A, B, S, IER)
************************************************************************
*  INTERPOLATION BY THE CUBIC SPLINE.                                  *
*  PARAMETERS                                                          *
*    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              *
*    (2) N: NUMBER OF KNOWN POINTS                                     *
*    (3) Y: 1-DIM. ARRAY FOR FUNCTION'S VALUES ON KNOWN POINTS         *
*    (4) C: 2-DIM. WORKING ARRAY                                       *
*    (5) NC: ROW SIZE OF THE ARRAY (C)                                 *
*    (6) A: LOWER LIMIT OF THE INTERPOLATION INTERVAL                  *
*    (7) B: UPPER LIMIT OF THE INTERPOLATION INTERVAL                  *
*    (8) S: INTEGRAL OF THE FUNCTION ON THE INTERVAL                   *
*    (9) IER: ERROR CODE                                               *
*  COPYRIGHT   T. OGUNI    JUNE 30 1989    VERSION 1.0                 *
************************************************************************
 
      INTEGER                        :: N, NC, IER
      double precision, dimension(N) :: X, Y
      double precision, dimension(NC,3) :: C
      double precision               :: A, B, S

      INTEGER                        :: I, IND, IA, IFLA, IFLB, IB, IP
      double precision               :: XX, X1, X2, H, SA, SAB, SB, XA
      double precision               :: XB

C
      IF (N .LT. 2 .OR. NC .LT. N-1) THEN
       IER = 2
       WRITE(*,*) '(SUBR. SPLQ) INVALID ARGUMENT. ',N, NC
       RETURN
      ENDIF
      IER = 0
      IND = 1
      IA = 1
      IFLA = 0
      IFLB = 0
      XX = A
      IF (A .GT. B) XX = B
   10 X1 = XX - X(IA)
      DO 20 I=IA,N-1
       IP = I
       X2 = XX - X(I+1)
       IF (X2 .LT. 0.0) GO TO 30
       IF (I .LT. N-1) X1 = X2
   20 CONTINUE
      IP = N - 1
      IF (X2 .GT. 0.0) IER = 1
   30 IF (X1 .LT. 0.0) THEN
       IER = 1
       X1 = - X1
       IF (IND .EQ. 1) IFLA = 1
       IF (IND .EQ. 2) IFLB = 1
      ENDIF
      IF (IND .NE. 2) THEN
       IND = 2
       IA = IP
       XA = X1
       XX = B
       IF (A .GT. B) XX = A
       GO TO 10
      ENDIF
      IB = IP
      XB = X1
C  COMPUTE THE INTEGRAL FROM A TO B.
      SA = Y(IA)+XA*(C(IA,1)/2.0D0+XA*(C(IA,2)/3.0D0+XA*C(IA,3)/4.0D0))
      SA = SA * XA
      IF (IFLA .EQ. 1) SA = - SA
      SAB = 0.0D0
      IF (IB-1 .GE. IA) THEN
       DO I=IA,IB-1
        H = X(I+1) - X(I)
        SAB = SAB+H*(Y(I+1)+Y(I)-(C(I+1,2)+C(I,2))*H*H/6.0D0)/2.0D0
       END DO
      ENDIF
      SB = Y(IB)+XB*(C(IB,1)/2.0D0+XB*(C(IB,2)/3.0D0+XB*C(IB,3)/4.0D0))
      SB = SB * XB
      IF (IFLB .EQ. 1) SB = - SB
      S = SB + SAB - SA
      IF (A .GT. B) S =  - S
C
      RETURN
      END SUBROUTINE SPLQ

************************************************************************
************************************************************************

      SUBROUTINE basisf(theta,phi,f2m)

      IMPLICIT NONE

      double precision, INTENT(in)                 :: theta, phi
      COMPLEX(KIND=8), INTENT(out), DIMENSION(5,2) :: f2m

      !  f2m(m,k)      k=1: W2m part,  k=2: X2m part
      !                m = 2,1,0,-1,-2 corresponds to 1,2,3,4,5

      COMPLEX(KIND=8), PARAMETER :: ui = (0.0d0,1.0d0)
      double precision   :: alpha, root5, root15, root2pi, root4pi

      root15 = sqrt(1.5d1)
      root5 = sqrt(5.0d0)
      root2pi = sqrt(8.0d0*datan(1.0d0))
      root4pi = sqrt(1.6d1*datan(1.0d0))
      alpha = 1.0d0/(4.0d0*sqrt(3.0d0))

      !++++ W2m ++++++
      f2m(1,1) = alpha*5.0d-1*root15/root2pi*(cos(theta)**2+1.0d0)
     &           *exp(2.0d0*ui*phi)  ! m=2
      f2m(2,1) = alpha*5.0d-1*root15/root2pi*sin(2.0d0*theta)
     &           *exp(ui*phi) ! m=1
      f2m(3,1) = alpha*root5/root4pi*3.0d0*sin(theta)**2 ! m=0
      f2m(4,1) = -alpha*5.0d-1*root15/root2pi*sin(2.0d0*theta)
     &           *exp(-ui*phi) ! m=-1
      f2m(5,1) = alpha*5.0d-1*root15/root2pi*(cos(theta)**2+1.0d0)
     &           *exp(-2.0d0*ui*phi) ! m=-2

      !++++ X2m ++++++
      f2m(1,2) = ui*alpha*root15/root2pi*cos(theta)*exp(ui*2.0d0*phi)
      f2m(2,2) = ui*alpha*root15/root2pi*sin(theta)*exp(ui*phi) ! m=1
      f2m(3,2) = (0.0d0,0.0d0)
      f2m(4,2) = ui*alpha*root15/root2pi*sin(theta)*exp(-ui*phi) ! m=-1
      f2m(5,2) = -ui*alpha*root15/root2pi*cos(theta)*exp(-ui*2.0d0*phi)

      RETURN
      END SUBROUTINE basisf

