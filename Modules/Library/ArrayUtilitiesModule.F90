MODULE ArrayUtilitiesModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CreatePackIndex
  PUBLIC :: ArrayPack
  PUBLIC :: ArrayUnpack
  PUBLIC :: ArrayCopy

  INTERFACE ArrayPack
    MODULE PROCEDURE ArrayPack1D_1_i
    MODULE PROCEDURE ArrayPack1D_1
    MODULE PROCEDURE ArrayPack1D_2
    MODULE PROCEDURE ArrayPack1D_3
    MODULE PROCEDURE ArrayPack1D_4
    MODULE PROCEDURE ArrayPack1D_8
    MODULE PROCEDURE ArrayPack2D_1
    MODULE PROCEDURE ArrayPack2D_2
    MODULE PROCEDURE ArrayPack2D_3
    MODULE PROCEDURE ArrayPack2D_4
    MODULE PROCEDURE ArrayPack2D_8
    MODULE PROCEDURE ArrayPack3D_1
    MODULE PROCEDURE ArrayPack3D_2
    MODULE PROCEDURE ArrayPack3D_3
    MODULE PROCEDURE ArrayPack3D_4
    MODULE PROCEDURE ArrayPack3D_8
  END INTERFACE ArrayPack

  INTERFACE ArrayUnpack
    MODULE PROCEDURE ArrayUnpack1D_1_i
    MODULE PROCEDURE ArrayUnpack1D_1
    MODULE PROCEDURE ArrayUnpack1D_2
    MODULE PROCEDURE ArrayUnpack1D_3
    MODULE PROCEDURE ArrayUnpack1D_4
    MODULE PROCEDURE ArrayUnpack1D_8
    MODULE PROCEDURE ArrayUnpack2D_1
    MODULE PROCEDURE ArrayUnpack2D_2
    MODULE PROCEDURE ArrayUnpack2D_3
    MODULE PROCEDURE ArrayUnpack2D_4
    MODULE PROCEDURE ArrayUnpack2D_8
    MODULE PROCEDURE ArrayUnpack3D_1
    MODULE PROCEDURE ArrayUnpack3D_2
    MODULE PROCEDURE ArrayUnpack3D_3
    MODULE PROCEDURE ArrayUnpack3D_4
    MODULE PROCEDURE ArrayUnpack3D_8
  END INTERFACE ArrayUnpack

  INTERFACE ArrayCopy
    MODULE PROCEDURE ArrayCopy1D_1_i
    MODULE PROCEDURE ArrayCopy1D_1
    MODULE PROCEDURE ArrayCopy1D_2
    MODULE PROCEDURE ArrayCopy1D_3
    MODULE PROCEDURE ArrayCopy1D_4
    MODULE PROCEDURE ArrayCopy1D_5
    MODULE PROCEDURE ArrayCopy1D_8
    MODULE PROCEDURE ArrayCopy2D_1
    MODULE PROCEDURE ArrayCopy2D_2
    MODULE PROCEDURE ArrayCopy2D_3
    MODULE PROCEDURE ArrayCopy2D_4
    MODULE PROCEDURE ArrayCopy2D_5
    MODULE PROCEDURE ArrayCopy2D_8
    MODULE PROCEDURE ArrayCopy3D_1
    MODULE PROCEDURE ArrayCopy3D_2
    MODULE PROCEDURE ArrayCopy3D_3
    MODULE PROCEDURE ArrayCopy3D_4
    MODULE PROCEDURE ArrayCopy3D_5
    MODULE PROCEDURE ArrayCopy3D_8
  END INTERFACE ArrayCopy

CONTAINS


  SUBROUTINE CreatePackIndex &
    ( MASK, nP, PackIndex, UnpackIndex )

    LOGICAL, DIMENSION(1:), INTENT(in)  :: MASK
    INTEGER,                INTENT(out) :: nP
    INTEGER, DIMENSION(1:), INTENT(out) :: PackIndex, UnpackIndex

    INTEGER :: i, iPack, iUnpack

    ! --- Build Lookup Tables ---

    nP      = COUNT( MASK )
    iPack   = 0
    iUnpack = nP
    DO i = 1, SIZE( MASK )
      IF ( MASK(i) ) THEN
        iPack = iPack + 1
        PackIndex(i) = iPack
        UnpackIndex(iPack) = i
      ELSE
        iUnpack = iUnpack + 1
        PackIndex(i) = iUnpack
        UnpackIndex(iUnpack) = i
      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( PackIndex, UnpackIndex )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( PackIndex, UnpackIndex )
#endif

  END SUBROUTINE CreatePackIndex


  SUBROUTINE ArrayPack1D_1_i &
    ( nP, UnpackIndex, X1, X1_P )

    INTEGER,                 INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:), INTENT(in)    :: UnpackIndex
    INTEGER,  DIMENSION(1:), INTENT(in)    :: X1
    INTEGER,  DIMENSION(1:), INTENT(inout) :: X1_P

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

!#if defined(THORNADO_OMP_OL)
!      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
!      !$OMP PRIVATE( i )
!#elif defined(THORNADO_OACC)
!      !$ACC PARALLEL LOOP GANG VECTOR &
!      !$ACC PRIVATE( i ) &
!      !$ACC PRESENT( UnpackIndex, X1, X1_P )
!#elif defined(THORNADO_OMP)
!      !$OMP PARALLEL DO &
!      !$OMP PRIVATE( i )
!#endif
      DO iPack = 1, nP
        i = UnpackIndex(iPack)
        X1_P(iPack) = X1(i)
      END DO

    ELSE

      CALL ArrayCopy( X1, X1_P )

    END IF

  END SUBROUTINE ArrayPack1D_1_i


  SUBROUTINE ArrayPack1D_1 &
    ( nP, UnpackIndex, X1, X1_P )

    INTEGER,                 INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:), INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:), INTENT(in)    :: X1
    REAL(DP), DIMENSION(1:), INTENT(inout) :: X1_P

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X1_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
        i = UnpackIndex(iPack)
        X1_P(iPack) = X1(i)
      END DO

    ELSE

      CALL ArrayCopy( X1, X1_P )

    END IF

  END SUBROUTINE ArrayPack1D_1


  SUBROUTINE ArrayPack1D_2 &
    ( nP, UnpackIndex, X1, X2, X1_P, X2_P )

    INTEGER,                 INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:), INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:), INTENT(in)    :: X1, X2
    REAL(DP), DIMENSION(1:), INTENT(inout) :: X1_P, X2_P

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X2, X1_P, X2_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
        i = UnpackIndex(iPack)
        X1_P(iPack) = X1(i)
        X2_P(iPack) = X2(i)
      END DO

    ELSE

      CALL ArrayCopy( X1, X2, X1_P, X2_P )

    END IF

  END SUBROUTINE ArrayPack1D_2


  SUBROUTINE ArrayPack1D_3 &
    ( nP, UnpackIndex, X1, X2, X3, X1_P, X2_P, X3_P )

    INTEGER,                 INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:), INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:), INTENT(in)    :: X1, X2, X3
    REAL(DP), DIMENSION(1:), INTENT(inout) :: X1_P, X2_P, X3_P

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X2, X3, X1_P, X2_P, X3_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
        i = UnpackIndex(iPack)
        X1_P(iPack) = X1(i)
        X2_P(iPack) = X2(i)
        X3_P(iPack) = X3(i)
      END DO

    ELSE

      CALL ArrayCopy( X1, X2, X3, X1_P, X2_P, X3_P )

    END IF

  END SUBROUTINE ArrayPack1D_3


  SUBROUTINE ArrayPack1D_4 &
    ( nP, UnpackIndex, X1, X2, X3, X4, X1_P, X2_P, X3_P, X4_P )

    INTEGER,                 INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:), INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:), INTENT(in)    :: X1, X2, X3, X4
    REAL(DP), DIMENSION(1:), INTENT(inout) :: X1_P, X2_P, X3_P, X4_P

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X2, X3, X4, X1_P, X2_P, X3_P, X4_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
        i = UnpackIndex(iPack)
        X1_P(iPack) = X1(i)
        X2_P(iPack) = X2(i)
        X3_P(iPack) = X3(i)
        X4_P(iPack) = X4(i)
      END DO

    ELSE

      CALL ArrayCopy( X1, X2, X3, X4, X1_P, X2_P, X3_P, X4_P )

    END IF

  END SUBROUTINE ArrayPack1D_4


  SUBROUTINE ArrayPack1D_8 &
    ( nP, UnpackIndex, &
      X1, X2, X3, X4, X5, X6, X7, X8, &
      X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P )

    INTEGER,                 INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:), INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:), INTENT(in)    :: X1, X2, X3, X4, X5, X6, X7, X8
    REAL(DP), DIMENSION(1:), INTENT(inout) :: X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, &
      !$ACC          X1, X2, X3, X4, X5, X6, X7, X8, &
      !$ACC          X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
        i = UnpackIndex(iPack)
        X1_P(iPack) = X1(i)
        X2_P(iPack) = X2(i)
        X3_P(iPack) = X3(i)
        X4_P(iPack) = X4(i)
        X5_P(iPack) = X5(i)
        X6_P(iPack) = X6(i)
        X7_P(iPack) = X7(i)
        X8_P(iPack) = X8(i)
      END DO

    ELSE

      CALL ArrayCopy &
             ( X1, X2, X3, X4, X5, X6, X7, X8, &
               X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P )

    END IF

  END SUBROUTINE ArrayPack1D_8


  SUBROUTINE ArrayPack2D_1 &
    ( nP, UnpackIndex, X1, X1_P )

    INTEGER,                    INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:),    INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:,1:), INTENT(in)    :: X1
    REAL(DP), DIMENSION(1:,1:), INTENT(inout) :: X1_P

    INTEGER  :: i, iPack, j

    IF ( nP < SIZE(X1,2) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X1_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
      DO j = 1, SIZE(X1,1)
        i = UnpackIndex(iPack)
        X1_P(j,iPack) = X1(j,i)
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1, X1_P )

    END IF

  END SUBROUTINE ArrayPack2D_1


  SUBROUTINE ArrayPack2D_2 &
    ( nP, UnpackIndex, X1, X2, X1_P, X2_P )

    INTEGER,                    INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:),    INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:,1:), INTENT(in)    :: X1, X2
    REAL(DP), DIMENSION(1:,1:), INTENT(inout) :: X1_P, X2_P

    INTEGER  :: i, iPack, j

    IF ( nP < SIZE(X1,2) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X2, X1_P, X2_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
      DO j = 1, SIZE(X1,1)
        i = UnpackIndex(iPack)
        X1_P(j,iPack) = X1(j,i)
        X2_P(j,iPack) = X2(j,i)
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1, X2, X1_P, X2_P )

    END IF

  END SUBROUTINE ArrayPack2D_2


  SUBROUTINE ArrayPack2D_3 &
    ( nP, UnpackIndex, X1, X2, X3, X1_P, X2_P, X3_P )

    INTEGER,                    INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:),    INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:,1:), INTENT(in)    :: X1, X2, X3
    REAL(DP), DIMENSION(1:,1:), INTENT(inout) :: X1_P, X2_P, X3_P

    INTEGER  :: i, iPack, j

    IF ( nP < SIZE(X1,2) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X2, X3, X1_P, X2_P, X3_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
      DO j = 1, SIZE(X1,1)
        i = UnpackIndex(iPack)
        X1_P(j,iPack) = X1(j,i)
        X2_P(j,iPack) = X2(j,i)
        X3_P(j,iPack) = X3(j,i)
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1, X2, X3, X1_P, X2_P, X3_P )

    END IF

  END SUBROUTINE ArrayPack2D_3


  SUBROUTINE ArrayPack2D_4 &
    ( nP, UnpackIndex, X1, X2, X3, X4, X1_P, X2_P, X3_P, X4_P )

    INTEGER,                    INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:),    INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:,1:), INTENT(in)    :: X1, X2, X3, X4
    REAL(DP), DIMENSION(1:,1:), INTENT(inout) :: X1_P, X2_P, X3_P, X4_P

    INTEGER  :: i, iPack, j

    IF ( nP < SIZE(X1,2) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X2, X3, X4, X1_P, X2_P, X3_P, X4_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
      DO j = 1, SIZE(X1,1)
        i = UnpackIndex(iPack)
        X1_P(j,iPack) = X1(j,i)
        X2_P(j,iPack) = X2(j,i)
        X3_P(j,iPack) = X3(j,i)
        X4_P(j,iPack) = X4(j,i)
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1, X2, X3, X4, X1_P, X2_P, X3_P, X4_P )

    END IF

  END SUBROUTINE ArrayPack2D_4


  SUBROUTINE ArrayPack2D_8 &
    ( nP, UnpackIndex, &
      X1, X2, X3, X4, X5, X6, X7, X8, &
      X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P )

    INTEGER,                    INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:),    INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:,1:), INTENT(in)    :: X1, X2, X3, X4, X5, X6, X7, X8
    REAL(DP), DIMENSION(1:,1:), INTENT(inout) :: X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P

    INTEGER  :: i, iPack, j

    IF ( nP < SIZE(X1,2) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, &
      !$ACC          X1, X2, X3, X4, X5, X6, X7, X8, &
      !$ACC          X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
      DO j = 1, SIZE(X1,1)
        i = UnpackIndex(iPack)
        X1_P(j,iPack) = X1(j,i)
        X2_P(j,iPack) = X2(j,i)
        X3_P(j,iPack) = X3(j,i)
        X4_P(j,iPack) = X4(j,i)
        X5_P(j,iPack) = X5(j,i)
        X6_P(j,iPack) = X6(j,i)
        X7_P(j,iPack) = X7(j,i)
        X8_P(j,iPack) = X8(j,i)
      END DO
      END DO

    ELSE

      CALL ArrayCopy &
             ( X1, X2, X3, X4, X5, X6, X7, X8, &
               X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P )

    END IF

  END SUBROUTINE ArrayPack2D_8


  SUBROUTINE ArrayPack3D_1 &
    ( nP, UnpackIndex, X1, X1_P )

    INTEGER,                       INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:),       INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)    :: X1
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(inout) :: X1_P

    INTEGER  :: i, iPack, j, k

    IF ( nP < SIZE(X1,3) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X1_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
      DO j = 1, SIZE(X1,2)
      DO k = 1, SIZE(X1,1)
        i = UnpackIndex(iPack)
        X1_P(k,j,iPack) = X1(k,j,i)
      END DO
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1, X1_P )

    END IF

  END SUBROUTINE ArrayPack3D_1


  SUBROUTINE ArrayPack3D_2 &
    ( nP, UnpackIndex, X1, X2, X1_P, X2_P )

    INTEGER,                       INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:),       INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)    :: X1, X2
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(inout) :: X1_P, X2_P

    INTEGER  :: i, iPack, j, k

    IF ( nP < SIZE(X1,3) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X2, X1_P, X2_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
      DO j = 1, SIZE(X1,2)
      DO k = 1, SIZE(X1,1)
        i = UnpackIndex(iPack)
        X1_P(k,j,iPack) = X1(k,j,i)
        X2_P(k,j,iPack) = X2(k,j,i)
      END DO
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1, X2, X1_P, X2_P )

    END IF

  END SUBROUTINE ArrayPack3D_2


  SUBROUTINE ArrayPack3D_3 &
    ( nP, UnpackIndex, X1, X2, X3, X1_P, X2_P, X3_P )

    INTEGER,                       INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:),       INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)    :: X1, X2, X3
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(inout) :: X1_P, X2_P, X3_P

    INTEGER  :: i, iPack, j, k

    IF ( nP < SIZE(X1,3) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X2, X3, X1_P, X2_P, X3_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
      DO j = 1, SIZE(X1,2)
      DO k = 1, SIZE(X1,1)
        i = UnpackIndex(iPack)
        X1_P(k,j,iPack) = X1(k,j,i)
        X2_P(k,j,iPack) = X2(k,j,i)
        X3_P(k,j,iPack) = X3(k,j,i)
      END DO
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1, X2, X3, X1_P, X2_P, X3_P )

    END IF

  END SUBROUTINE ArrayPack3D_3


  SUBROUTINE ArrayPack3D_4 &
    ( nP, UnpackIndex, X1, X2, X3, X4, X1_P, X2_P, X3_P, X4_P )

    INTEGER,                       INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:),       INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)    :: X1, X2, X3, X4
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(inout) :: X1_P, X2_P, X3_P, X4_P

    INTEGER  :: i, iPack, j, k

    IF ( nP < SIZE(X1,3) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, X1, X2, X3, X4, X1_P, X2_P, X3_P, X4_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
      DO j = 1, SIZE(X1,2)
      DO k = 1, SIZE(X1,1)
        i = UnpackIndex(iPack)
        X1_P(k,j,iPack) = X1(k,j,i)
        X2_P(k,j,iPack) = X2(k,j,i)
        X3_P(k,j,iPack) = X3(k,j,i)
        X4_P(k,j,iPack) = X4(k,j,i)
      END DO
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1, X2, X3, X4, X1_P, X2_P, X3_P, X4_P )

    END IF

  END SUBROUTINE ArrayPack3D_4


  SUBROUTINE ArrayPack3D_8 &
    ( nP, UnpackIndex, &
      X1, X2, X3, X4, X5, X6, X7, X8, &
      X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P )

    INTEGER,                       INTENT(in)    :: nP
    INTEGER,  DIMENSION(1:),       INTENT(in)    :: UnpackIndex
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)    :: X1, X2, X3, X4, X5, X6, X7, X8
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(inout) :: X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P

    INTEGER  :: i, iPack, j, k

    IF ( nP < SIZE(X1,3) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( i )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( i ) &
      !$ACC PRESENT( UnpackIndex, &
      !$ACC          X1, X2, X3, X4, X5, X6, X7, X8, &
      !$ACC          X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( i )
#endif
      DO iPack = 1, nP
      DO j = 1, SIZE(X1,2)
      DO k = 1, SIZE(X1,1)
        i = UnpackIndex(iPack)
        X1_P(k,j,iPack) = X1(k,j,i)
        X2_P(k,j,iPack) = X2(k,j,i)
        X3_P(k,j,iPack) = X3(k,j,i)
        X4_P(k,j,iPack) = X4(k,j,i)
        X5_P(k,j,iPack) = X5(k,j,i)
        X6_P(k,j,iPack) = X6(k,j,i)
        X7_P(k,j,iPack) = X7(k,j,i)
        X8_P(k,j,iPack) = X8(k,j,i)
      END DO
      END DO
      END DO

    ELSE

      CALL ArrayCopy &
             ( X1, X2, X3, X4, X5, X6, X7, X8, &
               X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P )

    END IF

  END SUBROUTINE ArrayPack3D_8


  SUBROUTINE ArrayUnpack1D_1_i &
    ( nP, MASK, PackIndex, X1_P, X1 )

    INTEGER,                 INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:), INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:), INTENT(in)    :: PackIndex
    INTEGER,  DIMENSION(1:), INTENT(in)    :: X1_P
    INTEGER,  DIMENSION(1:), INTENT(inout) :: X1

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

!#if defined(THORNADO_OMP_OL)
!      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
!      !$OMP PRIVATE( iPack )
!#elif defined(THORNADO_OACC)
!      !$ACC PARALLEL LOOP GANG VECTOR &
!      !$ACC PRIVATE( iPack ) &
!      !$ACC PRESENT( PackIndex, X1_P, X1 )
!#elif defined(THORNADO_OMP)
!      !$OMP PARALLEL DO &
!      !$OMP PRIVATE( iPack )
!#endif
      DO i = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(i) = X1_P(iPack)
      END IF
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X1 )

    END IF

  END SUBROUTINE ArrayUnpack1D_1_i


  SUBROUTINE ArrayUnpack1D_1 &
    ( nP, MASK, PackIndex, X1_P, X1 )

    INTEGER,                 INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:), INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:), INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:), INTENT(in)    :: X1_P
    REAL(DP), DIMENSION(1:), INTENT(inout) :: X1

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X1 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(i) = X1_P(iPack)
      END IF
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X1 )

    END IF

  END SUBROUTINE ArrayUnpack1D_1


  SUBROUTINE ArrayUnpack1D_2 &
    ( nP, MASK, PackIndex, X1_P, X2_P, X1, X2 )

    INTEGER,                 INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:), INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:), INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:), INTENT(in)    :: X1_P, X2_P
    REAL(DP), DIMENSION(1:), INTENT(inout) :: X1, X2

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X2_P, X1, X2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(i) = X1_P(iPack)
        X2(i) = X2_P(iPack)
      END IF
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X2_P, X1, X2 )

    END IF

  END SUBROUTINE ArrayUnpack1D_2


  SUBROUTINE ArrayUnpack1D_3 &
    ( nP, MASK, PackIndex, X1_P, X2_P, X3_P, X1, X2, X3 )

    INTEGER,                 INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:), INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:), INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:), INTENT(in)    :: X1_P, X2_P, X3_P
    REAL(DP), DIMENSION(1:), INTENT(inout) :: X1, X2, X3

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X2_P, X3_P, X1, X2, X3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(i) = X1_P(iPack)
        X2(i) = X2_P(iPack)
        X3(i) = X3_P(iPack)
      END IF
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X2_P, X3_P, X1, X2, X3 )

    END IF

  END SUBROUTINE ArrayUnpack1D_3


  SUBROUTINE ArrayUnpack1D_4 &
    ( nP, MASK, PackIndex, X1_P, X2_P, X3_P, X4_P, X1, X2, X3, X4 )

    INTEGER,                 INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:), INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:), INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:), INTENT(in)    :: X1_P, X2_P, X3_P, X4_P
    REAL(DP), DIMENSION(1:), INTENT(inout) :: X1, X2, X3, X4

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X2_P, X3_P, X4_P, X1, X2, X3, X4 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(i) = X1_P(iPack)
        X2(i) = X2_P(iPack)
        X3(i) = X3_P(iPack)
        X4(i) = X4_P(iPack)
      END IF
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X2_P, X3_P, X4_P, X1, X2, X3, X4 )

    END IF

  END SUBROUTINE ArrayUnpack1D_4


  SUBROUTINE ArrayUnpack1D_8 &
    ( nP, MASK, PackIndex, &
      X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P, &
      X1, X2, X3, X4, X5, X6, X7, X8 )

    INTEGER,                 INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:), INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:), INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:), INTENT(in)    :: X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P
    REAL(DP), DIMENSION(1:), INTENT(inout) :: X1, X2, X3, X4, X5, X6, X7, X8

    INTEGER  :: i, iPack

    IF ( nP < SIZE(X1,1) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, &
      !$ACC          X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P, &
      !$ACC          X1, X2, X3, X4, X5, X6, X7, X8 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(i) = X1_P(iPack)
        X2(i) = X2_P(iPack)
        X3(i) = X3_P(iPack)
        X4(i) = X4_P(iPack)
        X5(i) = X5_P(iPack)
        X6(i) = X6_P(iPack)
        X7(i) = X7_P(iPack)
        X8(i) = X8_P(iPack)
      END IF
      END DO

    ELSE

      CALL ArrayCopy &
             ( X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P, &
               X1, X2, X3, X4, X5, X6, X7, X8 )

    END IF

  END SUBROUTINE ArrayUnpack1D_8


  SUBROUTINE ArrayUnpack2D_1 &
    ( nP, MASK, PackIndex, X1_P, X1 )

    INTEGER,                    INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:),    INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:),    INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:,1:), INTENT(in)    :: X1_P
    REAL(DP), DIMENSION(1:,1:), INTENT(inout) :: X1

    INTEGER  :: i, iPack, j

    IF ( nP < SIZE(X1,2) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X1 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,2)
      DO j = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(j,i) = X1_P(j,iPack)
      END IF
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X1 )

    END IF

  END SUBROUTINE ArrayUnpack2D_1


  SUBROUTINE ArrayUnpack2D_2 &
    ( nP, MASK, PackIndex, X1_P, X2_P, X1, X2 )

    INTEGER,                    INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:),    INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:),    INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:,1:), INTENT(in)    :: X1_P, X2_P
    REAL(DP), DIMENSION(1:,1:), INTENT(inout) :: X1, X2

    INTEGER  :: i, iPack, j

    IF ( nP < SIZE(X1,2) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X2_P, X1, X2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,2)
      DO j = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(j,i) = X1_P(j,iPack)
        X2(j,i) = X2_P(j,iPack)
      END IF
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X2_P, X1, X2 )

    END IF

  END SUBROUTINE ArrayUnpack2D_2


  SUBROUTINE ArrayUnpack2D_3 &
    ( nP, MASK, PackIndex, X1_P, X2_P, X3_P, X1, X2, X3 )

    INTEGER,                    INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:),    INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:),    INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:,1:), INTENT(in)    :: X1_P, X2_P, X3_P
    REAL(DP), DIMENSION(1:,1:), INTENT(inout) :: X1, X2, X3

    INTEGER  :: i, iPack, j

    IF ( nP < SIZE(X1,2) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X2_P, X3_P, X1, X2, X3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,2)
      DO j = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(j,i) = X1_P(j,iPack)
        X2(j,i) = X2_P(j,iPack)
        X3(j,i) = X3_P(j,iPack)
      END IF
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X2_P, X3_P, X1, X2, X3 )

    END IF

  END SUBROUTINE ArrayUnpack2D_3


  SUBROUTINE ArrayUnpack2D_4 &
    ( nP, MASK, PackIndex, X1_P, X2_P, X3_P, X4_P, X1, X2, X3, X4 )

    INTEGER,                    INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:),    INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:),    INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:,1:), INTENT(in)    :: X1_P, X2_P, X3_P, X4_P
    REAL(DP), DIMENSION(1:,1:), INTENT(inout) :: X1, X2, X3, X4

    INTEGER  :: i, iPack, j

    IF ( nP < SIZE(X1,2) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X2_P, X3_P, X4_P, X1, X2, X3, X4 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,2)
      DO j = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(j,i) = X1_P(j,iPack)
        X2(j,i) = X2_P(j,iPack)
        X3(j,i) = X3_P(j,iPack)
        X4(j,i) = X4_P(j,iPack)
      END IF
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X2_P, X3_P, X4_P, X1, X2, X3, X4 )

    END IF

  END SUBROUTINE ArrayUnpack2D_4


  SUBROUTINE ArrayUnpack2D_8 &
    ( nP, MASK, PackIndex, &
      X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P, &
      X1, X2, X3, X4, X5, X6, X7, X8 )

    INTEGER,                    INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:),    INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:),    INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:,1:), INTENT(in)    :: X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P
    REAL(DP), DIMENSION(1:,1:), INTENT(inout) :: X1, X2, X3, X4, X5, X6, X7, X8

    INTEGER  :: i, iPack, j

    IF ( nP < SIZE(X1,2) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, &
      !$ACC          X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P, &
      !$ACC          X1, X2, X3, X4, X5, X6, X7, X8 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,2)
      DO j = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(j,i) = X1_P(j,iPack)
        X2(j,i) = X2_P(j,iPack)
        X3(j,i) = X3_P(j,iPack)
        X4(j,i) = X4_P(j,iPack)
        X5(j,i) = X5_P(j,iPack)
        X6(j,i) = X6_P(j,iPack)
        X7(j,i) = X7_P(j,iPack)
        X8(j,i) = X8_P(j,iPack)
      END IF
      END DO
      END DO

    ELSE

      CALL ArrayCopy &
             ( X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P, &
               X1, X2, X3, X4, X5, X6, X7, X8 )

    END IF

  END SUBROUTINE ArrayUnpack2D_8


  SUBROUTINE ArrayUnpack3D_1 &
    ( nP, MASK, PackIndex, X1_P, X1 )

    INTEGER,                       INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:),       INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:),       INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)    :: X1_P
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(inout) :: X1

    INTEGER  :: i, iPack, j, k

    IF ( nP < SIZE(X1,3) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X1 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,3)
      DO j = 1, SIZE(X1,2)
      DO k = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(k,j,i) = X1_P(k,j,iPack)
      END IF
      END DO
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X1 )

    END IF

  END SUBROUTINE ArrayUnpack3D_1


  SUBROUTINE ArrayUnpack3D_2 &
    ( nP, MASK, PackIndex, X1_P, X2_P, X1, X2 )

    INTEGER,                       INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:),       INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:),       INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)    :: X1_P, X2_P
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(inout) :: X1, X2

    INTEGER  :: i, iPack, j, k

    IF ( nP < SIZE(X1,3) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X2_P, X1, X2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,3)
      DO j = 1, SIZE(X1,2)
      DO k = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(k,j,i) = X1_P(k,j,iPack)
        X2(k,j,i) = X2_P(k,j,iPack)
      END IF
      END DO
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X2_P, X1, X2 )

    END IF

  END SUBROUTINE ArrayUnpack3D_2


  SUBROUTINE ArrayUnpack3D_3 &
    ( nP, MASK, PackIndex, X1_P, X2_P, X3_P, X1, X2, X3 )

    INTEGER,                       INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:),       INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:),       INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)    :: X1_P, X2_P, X3_P
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(inout) :: X1, X2, X3

    INTEGER  :: i, iPack, j, k

    IF ( nP < SIZE(X1,3) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X2_P, X3_P, X1, X2, X3 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,3)
      DO j = 1, SIZE(X1,2)
      DO k = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(k,j,i) = X1_P(k,j,iPack)
        X2(k,j,i) = X2_P(k,j,iPack)
        X3(k,j,i) = X3_P(k,j,iPack)
      END IF
      END DO
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X2_P, X3_P, X1, X2, X3 )

    END IF

  END SUBROUTINE ArrayUnpack3D_3


  SUBROUTINE ArrayUnpack3D_4 &
    ( nP, MASK, PackIndex, X1_P, X2_P, X3_P, X4_P, X1, X2, X3, X4 )

    INTEGER,                       INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:),       INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:),       INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)    :: X1_P, X2_P, X3_P, X4_P
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(inout) :: X1, X2, X3, X4

    INTEGER  :: i, iPack, j, k

    IF ( nP < SIZE(X1,3) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, X1_P, X2_P, X3_P, X4_P, X1, X2, X3, X4 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,3)
      DO j = 1, SIZE(X1,2)
      DO k = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(k,j,i) = X1_P(k,j,iPack)
        X2(k,j,i) = X2_P(k,j,iPack)
        X3(k,j,i) = X3_P(k,j,iPack)
        X4(k,j,i) = X4_P(k,j,iPack)
      END IF
      END DO
      END DO
      END DO

    ELSE

      CALL ArrayCopy( X1_P, X2_P, X3_P, X4_P, X1, X2, X3, X4 )

    END IF

  END SUBROUTINE ArrayUnpack3D_4


  SUBROUTINE ArrayUnpack3D_8 &
    ( nP, MASK, PackIndex, &
      X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P, &
      X1, X2, X3, X4, X5, X6, X7, X8 )

    INTEGER,                       INTENT(in)    :: nP
    LOGICAL,  DIMENSION(1:),       INTENT(in)    :: MASK
    INTEGER,  DIMENSION(1:),       INTENT(in)    :: PackIndex
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)    :: X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(inout) :: X1, X2, X3, X4, X5, X6, X7, X8

    INTEGER  :: i, iPack, j, k

    IF ( nP < SIZE(X1,3) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( iPack )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( iPack ) &
      !$ACC PRESENT( PackIndex, &
      !$ACC          X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P, &
      !$ACC          X1, X2, X3, X4, X5, X6, X7, X8 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( iPack )
#endif
      DO i = 1, SIZE(X1,3)
      DO j = 1, SIZE(X1,2)
      DO k = 1, SIZE(X1,1)
      IF ( MASK(i) ) THEN
        iPack = PackIndex(i)
        X1(k,j,i) = X1_P(k,j,iPack)
        X2(k,j,i) = X2_P(k,j,iPack)
        X3(k,j,i) = X3_P(k,j,iPack)
        X4(k,j,i) = X4_P(k,j,iPack)
        X5(k,j,i) = X5_P(k,j,iPack)
        X6(k,j,i) = X6_P(k,j,iPack)
        X7(k,j,i) = X7_P(k,j,iPack)
        X8(k,j,i) = X8_P(k,j,iPack)
      END IF
      END DO
      END DO
      END DO

    ELSE

      CALL ArrayCopy &
             ( X1_P, X2_P, X3_P, X4_P, X5_P, X6_P, X7_P, X8_P, &
               X1, X2, X3, X4, X5, X6, X7, X8 )

    END IF

  END SUBROUTINE ArrayUnpack3D_8


  SUBROUTINE ArrayCopy1D_1_i &
    ( X1, Y1 )

    INTEGER,  DIMENSION(1:), INTENT(in)  :: X1
    INTEGER,  DIMENSION(1:), INTENT(out) :: Y1

    INTEGER  :: i

!#if defined(THORNADO_OMP_OL)
!    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
!#elif defined(THORNADO_OACC)
!    !$ACC PARALLEL LOOP GANG VECTOR &
!    !$ACC PRESENT( X1, Y1 )
!#elif defined(THORNADO_OMP)
!    !$OMP PARALLEL DO
!#endif
    DO i = 1, SIZE(X1,1)
      Y1(i) = X1(i)
    END DO

  END SUBROUTINE ArrayCopy1D_1_i


  SUBROUTINE ArrayCopy1D_1 &
    ( X1, Y1 )

    REAL(DP), DIMENSION(1:), INTENT(in)  :: X1
    REAL(DP), DIMENSION(1:), INTENT(out) :: Y1

    INTEGER  :: i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( X1, Y1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO i = 1, SIZE(X1,1)
      Y1(i) = X1(i)
    END DO

  END SUBROUTINE ArrayCopy1D_1


  SUBROUTINE ArrayCopy1D_2 &
    ( X1, X2, Y1, Y2 )

    REAL(DP), DIMENSION(1:), INTENT(in)  :: X1, X2
    REAL(DP), DIMENSION(1:), INTENT(out) :: Y1, Y2

    INTEGER  :: i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( X1, X2, Y1, Y2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO i = 1, SIZE(X1,1)
      Y1(i) = X1(i)
      Y2(i) = X2(i)
    END DO

  END SUBROUTINE ArrayCopy1D_2


  SUBROUTINE ArrayCopy1D_3 &
    ( X1, X2, X3, Y1, Y2, Y3 )

    REAL(DP), DIMENSION(1:), INTENT(in)  :: X1, X2, X3
    REAL(DP), DIMENSION(1:), INTENT(out) :: Y1, Y2, Y3

    INTEGER  :: i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( X1, X2, X3, Y1, Y2, Y3 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO i = 1, SIZE(X1,1)
      Y1(i) = X1(i)
      Y2(i) = X2(i)
      Y3(i) = X3(i)
    END DO

  END SUBROUTINE ArrayCopy1D_3


  SUBROUTINE ArrayCopy1D_4 &
    ( X1, X2, X3, X4, Y1, Y2, Y3, Y4 )

    REAL(DP), DIMENSION(1:), INTENT(in)  :: X1, X2, X3, X4
    REAL(DP), DIMENSION(1:), INTENT(out) :: Y1, Y2, Y3, Y4

    INTEGER  :: i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( X1, X2, X3, X4, Y1, Y2, Y3, Y4 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO i = 1, SIZE(X1,1)
      Y1(i) = X1(i)
      Y2(i) = X2(i)
      Y3(i) = X3(i)
      Y4(i) = X4(i)
    END DO

  END SUBROUTINE ArrayCopy1D_4


  SUBROUTINE ArrayCopy1D_5 &
    ( X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5 )

    REAL(DP), DIMENSION(1:), INTENT(in)  :: X1, X2, X3, X4, X5
    REAL(DP), DIMENSION(1:), INTENT(out) :: Y1, Y2, Y3, Y4, Y5

    INTEGER  :: i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO i = 1, SIZE(X1,1)
      Y1(i) = X1(i)
      Y2(i) = X2(i)
      Y3(i) = X3(i)
      Y4(i) = X4(i)
      Y5(i) = X5(i)
    END DO

  END SUBROUTINE ArrayCopy1D_5


  SUBROUTINE ArrayCopy1D_8 &
    ( X1, X2, X3, X4, X5, X6, X7, X8, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8 )

    REAL(DP), DIMENSION(1:), INTENT(in)  :: X1, X2, X3, X4, X5, X6, X7, X8
    REAL(DP), DIMENSION(1:), INTENT(out) :: Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8

    INTEGER  :: i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( X1, X2, X3, X4, X5, X6, X7, X8, &
    !$ACC          Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO i = 1, SIZE(X1,1)
      Y1(i) = X1(i)
      Y2(i) = X2(i)
      Y3(i) = X3(i)
      Y4(i) = X4(i)
      Y5(i) = X5(i)
      Y6(i) = X6(i)
      Y7(i) = X7(i)
      Y8(i) = X8(i)
    END DO

  END SUBROUTINE ArrayCopy1D_8


  SUBROUTINE ArrayCopy2D_1 &
    ( X1, Y1 )

    REAL(DP), DIMENSION(1:,1:), INTENT(in)  :: X1
    REAL(DP), DIMENSION(1:,1:), INTENT(out) :: Y1

    INTEGER  :: i, j

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( X1, Y1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO i = 1, SIZE(X1,2)
    DO j = 1, SIZE(X1,1)
      Y1(j,i) = X1(j,i)
    END DO
    END DO

  END SUBROUTINE ArrayCopy2D_1


  SUBROUTINE ArrayCopy2D_2 &
    ( X1, X2, Y1, Y2 )

    REAL(DP), DIMENSION(1:,1:), INTENT(in)  :: X1, X2
    REAL(DP), DIMENSION(1:,1:), INTENT(out) :: Y1, Y2

    INTEGER  :: i, j

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( X1, X2, Y1, Y2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO i = 1, SIZE(X1,2)
    DO j = 1, SIZE(X1,1)
      Y1(j,i) = X1(j,i)
      Y2(j,i) = X2(j,i)
    END DO
    END DO

  END SUBROUTINE ArrayCopy2D_2


  SUBROUTINE ArrayCopy2D_3 &
    ( X1, X2, X3, Y1, Y2, Y3 )

    REAL(DP), DIMENSION(1:,1:), INTENT(in)  :: X1, X2, X3
    REAL(DP), DIMENSION(1:,1:), INTENT(out) :: Y1, Y2, Y3

    INTEGER  :: i, j

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( X1, X2, X3, Y1, Y2, Y3 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO i = 1, SIZE(X1,2)
    DO j = 1, SIZE(X1,1)
      Y1(j,i) = X1(j,i)
      Y2(j,i) = X2(j,i)
      Y3(j,i) = X3(j,i)
    END DO
    END DO

  END SUBROUTINE ArrayCopy2D_3


  SUBROUTINE ArrayCopy2D_4 &
    ( X1, X2, X3, X4, Y1, Y2, Y3, Y4 )

    REAL(DP), DIMENSION(1:,1:), INTENT(in)  :: X1, X2, X3, X4
    REAL(DP), DIMENSION(1:,1:), INTENT(out) :: Y1, Y2, Y3, Y4

    INTEGER  :: i, j

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( X1, X2, X3, X4, Y1, Y2, Y3, Y4 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO i = 1, SIZE(X1,2)
    DO j = 1, SIZE(X1,1)
      Y1(j,i) = X1(j,i)
      Y2(j,i) = X2(j,i)
      Y3(j,i) = X3(j,i)
      Y4(j,i) = X4(j,i)
    END DO
    END DO

  END SUBROUTINE ArrayCopy2D_4


  SUBROUTINE ArrayCopy2D_5 &
    ( X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5 )

    REAL(DP), DIMENSION(1:,1:), INTENT(in)  :: X1, X2, X3, X4, X5
    REAL(DP), DIMENSION(1:,1:), INTENT(out) :: Y1, Y2, Y3, Y4, Y5

    INTEGER  :: i, j

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO i = 1, SIZE(X1,2)
    DO j = 1, SIZE(X1,1)
      Y1(j,i) = X1(j,i)
      Y2(j,i) = X2(j,i)
      Y3(j,i) = X3(j,i)
      Y4(j,i) = X4(j,i)
      Y5(j,i) = X5(j,i)
    END DO
    END DO

  END SUBROUTINE ArrayCopy2D_5


  SUBROUTINE ArrayCopy2D_8 &
    ( X1, X2, X3, X4, X5, X6, X7, X8, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8 )

    REAL(DP), DIMENSION(1:,1:), INTENT(in)  :: X1, X2, X3, X4, X5, X6, X7, X8
    REAL(DP), DIMENSION(1:,1:), INTENT(out) :: Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8

    INTEGER  :: i, j

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( X1, X2, X3, X4, X5, X6, X7, X8, &
    !$ACC          Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO i = 1, SIZE(X1,2)
    DO j = 1, SIZE(X1,1)
      Y1(j,i) = X1(j,i)
      Y2(j,i) = X2(j,i)
      Y3(j,i) = X3(j,i)
      Y4(j,i) = X4(j,i)
      Y5(j,i) = X5(j,i)
      Y6(j,i) = X6(j,i)
      Y7(j,i) = X7(j,i)
      Y8(j,i) = X8(j,i)
    END DO
    END DO

  END SUBROUTINE ArrayCopy2D_8


  SUBROUTINE ArrayCopy3D_1 &
    ( X1, Y1 )

    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)  :: X1
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(out) :: Y1

    INTEGER  :: i, j, k

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( X1, Y1 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO i = 1, SIZE(X1,3)
    DO j = 1, SIZE(X1,2)
    DO k = 1, SIZE(X1,1)
      Y1(k,j,i) = X1(k,j,i)
    END DO
    END DO
    END DO

  END SUBROUTINE ArrayCopy3D_1


  SUBROUTINE ArrayCopy3D_2 &
    ( X1, X2, Y1, Y2 )

    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)  :: X1, X2
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(out) :: Y1, Y2

    INTEGER  :: i, j, k

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( X1, X2, Y1, Y2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO i = 1, SIZE(X1,3)
    DO j = 1, SIZE(X1,2)
    DO k = 1, SIZE(X1,1)
      Y1(k,j,i) = X1(k,j,i)
      Y2(k,j,i) = X2(k,j,i)
    END DO
    END DO
    END DO

  END SUBROUTINE ArrayCopy3D_2


  SUBROUTINE ArrayCopy3D_3 &
    ( X1, X2, X3, Y1, Y2, Y3 )

    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)  :: X1, X2, X3
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(out) :: Y1, Y2, Y3

    INTEGER  :: i, j, k

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( X1, X2, X3, Y1, Y2, Y3 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO i = 1, SIZE(X1,3)
    DO j = 1, SIZE(X1,2)
    DO k = 1, SIZE(X1,1)
      Y1(k,j,i) = X1(k,j,i)
      Y2(k,j,i) = X2(k,j,i)
      Y3(k,j,i) = X3(k,j,i)
    END DO
    END DO
    END DO

  END SUBROUTINE ArrayCopy3D_3


  SUBROUTINE ArrayCopy3D_4 &
    ( X1, X2, X3, X4, Y1, Y2, Y3, Y4 )

    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)  :: X1, X2, X3, X4
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(out) :: Y1, Y2, Y3, Y4

    INTEGER  :: i, j, k

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( X1, X2, X3, X4, Y1, Y2, Y3, Y4 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO i = 1, SIZE(X1,3)
    DO j = 1, SIZE(X1,2)
    DO k = 1, SIZE(X1,1)
      Y1(k,j,i) = X1(k,j,i)
      Y2(k,j,i) = X2(k,j,i)
      Y3(k,j,i) = X3(k,j,i)
      Y4(k,j,i) = X4(k,j,i)
    END DO
    END DO
    END DO

  END SUBROUTINE ArrayCopy3D_4


  SUBROUTINE ArrayCopy3D_5 &
    ( X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5 )

    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)  :: X1, X2, X3, X4, X5
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(out) :: Y1, Y2, Y3, Y4, Y5

    INTEGER  :: i, j, k

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO i = 1, SIZE(X1,3)
    DO j = 1, SIZE(X1,2)
    DO k = 1, SIZE(X1,1)
      Y1(k,j,i) = X1(k,j,i)
      Y2(k,j,i) = X2(k,j,i)
      Y3(k,j,i) = X3(k,j,i)
      Y4(k,j,i) = X4(k,j,i)
      Y5(k,j,i) = X5(k,j,i)
    END DO
    END DO
    END DO

  END SUBROUTINE ArrayCopy3D_5


  SUBROUTINE ArrayCopy3D_8 &
    ( X1, X2, X3, X4, X5, X6, X7, X8, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8 )

    REAL(DP), DIMENSION(1:,1:,1:), INTENT(in)  :: X1, X2, X3, X4, X5, X6, X7, X8
    REAL(DP), DIMENSION(1:,1:,1:), INTENT(out) :: Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8

    INTEGER  :: i, j, k

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( X1, X2, X3, X4, X5, X6, X7, X8, &
    !$ACC          Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO i = 1, SIZE(X1,3)
    DO j = 1, SIZE(X1,2)
    DO k = 1, SIZE(X1,1)
      Y1(k,j,i) = X1(k,j,i)
      Y2(k,j,i) = X2(k,j,i)
      Y3(k,j,i) = X3(k,j,i)
      Y4(k,j,i) = X4(k,j,i)
      Y5(k,j,i) = X5(k,j,i)
      Y6(k,j,i) = X6(k,j,i)
      Y7(k,j,i) = X7(k,j,i)
      Y8(k,j,i) = X8(k,j,i)
    END DO
    END DO
    END DO

  END SUBROUTINE ArrayCopy3D_8


END MODULE ArrayUtilitiesModule
