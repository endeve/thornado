MODULE TwoMoment_BoundaryConditionsModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    bcZ, swZ, &
    nNodesZ, nDOF
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable4D
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_G1, iCR_G2, iCR_G3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_TwoMoment

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: U, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( U, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#endif

    CALL ApplyBC_TwoMoment_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    CALL ApplyBC_TwoMoment_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    CALL ApplyBC_TwoMoment_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: U ) &
    !$OMP MAP( release: iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( U ) &
    !$ACC DELETE( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#endif

  END SUBROUTINE ApplyBoundaryConditions_TwoMoment


  SUBROUTINE ApplyBC_TwoMoment_X1( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER :: iNode, iS, iCR, iZ1, iZ2, iZ3, iZ4
    INTEGER :: iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4
    INTEGER :: jNodeZ2, iNodeZ, jNodeZ

    SELECT CASE ( bcZ(2) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(7)
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = 1, swZ(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)
                  DO iNode = 1, nDOF

                    ! --- Inner Boundary ---

                    U(iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ_E0(2)-(iZ2-1),iZ3,iZ4,iCR,iS)

                    ! --- Outer Boundary ---

                    U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ_B0(2)+(iZ2-1),iZ3,iZ4,iCR,iS)

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 2 ) ! Homogeneous

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(7)
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = 1, swZ(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)
                  DO iNode = 1, nDOF

                    ! --- Inner Boundary ---

                    U(iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ_B0(2),iZ3,iZ4,iCR,iS)

                    ! --- Outer Boundary ---

                    U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR,iS)

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 3 ) ! Reflecting

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ2, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
      !$ACC PRIVATE( jNodeZ2, iNodeZ, jNodeZ ) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ2, iNodeZ, jNodeZ )
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = 1, swZ(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)

                  IF( iCR == iCR_G1 )THEN

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ2 = (nNodesZ(2)-iNodeZ2) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, jNodeZ2, iNodeZ3, iNodeZ4 )

                      ! --- Inner Boundary ---

                      U(iNodeZ,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                        = - U(jNodeZ,iZ1,iZ_B0(2),iZ3,iZ4,iCR,iS)

                      ! --- Outer Boundary ---

                      U(iNodeZ,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                        = - U(jNodeZ,iZ1,iZ_E0(2),iZ3,iZ4,iCR,iS)

                    END DO
                    END DO
                    END DO
                    END DO

                  ELSE

                    DO iNode = 1, nDOF

                      ! --- Inner Boundary ---

                      U(iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                        = U(iNode,iZ1,iZ_B0(2),iZ3,iZ4,iCR,iS)

                      ! --- Outer Boundary ---

                      U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                        = U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR,iS)

                    END DO

                  END IF

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 32 ) ! Reflecting (Inner), Homogeneous (Outer)

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ2, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
      !$ACC PRIVATE( jNodeZ2, iNodeZ, jNodeZ ) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ2, iNodeZ, jNodeZ )
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = 1, swZ(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)

                  ! --- Inner Boundary (Reflecting) ---

                  IF( iCR == iCR_G1 )THEN

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ2 = (nNodesZ(2)-iNodeZ2) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, jNodeZ2, iNodeZ3, iNodeZ4 )


                      U(iNodeZ,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                        = - U(jNodeZ,iZ1,iZ_B0(2),iZ3,iZ4,iCR,iS)

                    END DO
                    END DO
                    END DO
                    END DO

                  ELSE

                    DO iNode = 1, nDOF

                      U(iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                        = U(iNode,iZ1,iZ_B0(2),iZ3,iZ4,iCR,iS)

                    END DO

                  END IF

                  ! --- Outer Boundary (Homogeneous) ---

                  DO iNode = 1, nDOF

                    U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR,iS)

                  END DO

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A45,I2.2)') &
        '', 'Invalid Boundary Condition for TwoMoment X1: ', bcZ(2)
      STOP

    END SELECT

  END SUBROUTINE ApplyBC_TwoMoment_X1


  SUBROUTINE ApplyBC_TwoMoment_X2( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER :: iNode, iS, iCR, iZ1, iZ2, iZ3, iZ4
    INTEGER :: iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4
    INTEGER :: jNodeZ3, iNodeZ, jNodeZ

    SELECT CASE ( bcZ(3) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(7)
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = 1, swZ(3)
              DO iZ2 = iZ_B0(2), iZ_E0(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)
                  DO iNode = 1, nDOF

                    ! --- Inner Boundary ---

                    U(iNode,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ2,iZ_E0(3)-(iZ3-1),iZ4,iCR,iS)

                    ! --- Outer Boundary ---

                    U(iNode,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ2,iZ_B0(3)+(iZ3-1),iZ4,iCR,iS)

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 2 ) ! Homogeneous

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(7)
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = 1, swZ(3)
              DO iZ2 = iZ_B0(2), iZ_E0(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)
                  DO iNode = 1, nDOF

                    ! --- Inner Boundary ---

                    U(iNode,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ2,iZ_B0(3),iZ4,iCR,iS)

                    ! --- Outer Boundary ---

                    U(iNode,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ2,iZ_E0(3),iZ4,iCR,iS)

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 3 ) ! Reflecting

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ3, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
      !$ACC PRIVATE( jNodeZ3, iNodeZ, jNodeZ ) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ3, iNodeZ, jNodeZ )
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = 1, swZ(3)
              DO iZ2 = iZ_B0(2), iZ_E0(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)

                  IF( iCR == iCR_G2 )THEN

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ3 = ( nNodesZ(3) - iNodeZ3 ) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, jNodeZ3, iNodeZ4 )

                      ! --- Inner Boundary ---

                      U(iNodeZ,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                        = - U(jNodeZ,iZ1,iZ2,iZ_B0(3),iZ4,iCR,iS)

                      ! --- Outer Boundary ---

                      U(iNodeZ,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                        = - U(jNodeZ,iZ1,iZ2,iZ_E0(3),iZ4,iCR,iS)

                    END DO
                    END DO
                    END DO
                    END DO

                  ELSE

                    DO iNode = 1, nDOF

                      ! --- Inner Boundary ---

                      U(iNode,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                        = U(iNode,iZ1,iZ2,iZ_B0(3),iZ4,iCR,iS)

                      ! --- Outer Boundary ---

                      U(iNode,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                        = U(iNode,iZ1,iZ2,iZ_E0(3),iZ4,iCR,iS)

                    END DO

                  END IF

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 32 ) ! Reflecting (Inner), Homogeneous (Outer)

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ3, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
      !$ACC PRIVATE( jNodeZ3, iNodeZ, jNodeZ ) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ3, iNodeZ, jNodeZ )
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = 1, swZ(3)
              DO iZ2 = iZ_B0(2), iZ_E0(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)

                  ! --- Inner Boundary (Reflecting) ---

                  IF( iCR == iCR_G2 )THEN

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ3 = ( nNodesZ(3) - iNodeZ3 ) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, jNodeZ3, iNodeZ4 )

                      U(iNodeZ,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                        = - U(jNodeZ,iZ1,iZ2,iZ_B0(3),iZ4,iCR,iS)

                    END DO
                    END DO
                    END DO
                    END DO

                  ELSE

                    DO iNode = 1, nDOF

                      U(iNode,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                        = U(iNode,iZ1,iZ2,iZ_B0(3),iZ4,iCR,iS)

                    END DO

                  END IF

                  ! --- Outer Boundary (Homogeneous) ---

                  DO iNode = 1, nDOF

                    U(iNode,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ2,iZ_E0(3),iZ4,iCR,iS)

                  END DO

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A45,I2.2)') &
        '', 'Invalid Boundary Condition for TwoMoment X2: ', bcZ(3)
      STOP

    END SELECT

  END SUBROUTINE ApplyBC_TwoMoment_X2


  SUBROUTINE ApplyBC_TwoMoment_X3( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER :: iNode, iS, iCR, iZ1, iZ2, iZ3, iZ4
    INTEGER :: iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4
    INTEGER :: jNodeZ4, iNodeZ, jNodeZ

    SELECT CASE ( bcZ(4) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(7)
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = 1, swZ(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = iZ_B0(2), iZ_E0(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)
                  DO iNode = 1, nDOF

                    ! --- Inner Boundary ---

                    U(iNode,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ2,iZ3,iZ_E0(4)-(iZ4-1),iCR,iS)

                    ! --- Outer Boundary ---

                    U(iNode,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ2,iZ3,iZ_B0(4)+(iZ4-1),iCR,iS)

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 2 ) ! Homogeneous

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(7)
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = 1, swZ(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = iZ_B0(2), iZ_E0(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)
                  DO iNode = 1, nDOF

                    ! --- Inner Boundary ---

                    U(iNode,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ2,iZ3,iZ_B0(4),iCR,iS)

                    ! --- Outer Boundary ---

                    U(iNode,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ2,iZ3,iZ_E0(4),iCR,iS)

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 3 ) ! Reflecting

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ4, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
      !$ACC PRIVATE( jNodeZ4, iNodeZ, jNodeZ ) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ4, iNodeZ, jNodeZ )
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = 1, swZ(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = iZ_B0(2), iZ_E0(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)

                  IF( iCR == iCR_G3 )THEN

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ4 = (nNodesZ(4)-iNodeZ4) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, jNodeZ4 )

                      ! --- Inner Boundary ---

                      U(iNodeZ,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                        = - U(jNodeZ,iZ1,iZ2,iZ3,iZ_B0(4),iCR,iS)

                      ! --- Outer Boundary ---

                      U(iNodeZ,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                        = - U(jNodeZ,iZ1,iZ2,iZ3,iZ_E0(4),iCR,iS)

                    END DO
                    END DO
                    END DO
                    END DO

                  ELSE

                    DO iNode = 1, nDOF

                      ! --- Inner Boundary ---

                      U(iNode,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                        = U(iNode,iZ1,iZ2,iZ3,iZ_B0(4),iCR,iS)

                      ! --- Outer Boundary ---

                      U(iNode,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                        = U(iNode,iZ1,iZ2,iZ3,iZ_E0(4),iCR,iS)

                    END DO

                  END IF

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 32 ) ! Reflecting (Inner), Homogeneous (Outer)

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ4, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
      !$ACC PRIVATE( jNodeZ4, iNodeZ, jNodeZ ) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(6) &
      !$OMP PRIVATE( jNodeZ4, iNodeZ, jNodeZ )
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = 1, swZ(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = iZ_B0(2), iZ_E0(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)

                  ! --- Inner Boundary (Reflecting) ---

                  IF( iCR == iCR_G3 )THEN

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ4 = ( nNodesZ(4) - iNodeZ4 ) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, jNodeZ4 )


                      U(iNodeZ,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                        = - U(jNodeZ,iZ1,iZ2,iZ3,iZ_B0(4),iCR,iS)

                    END DO
                    END DO
                    END DO
                    END DO

                  ELSE

                    DO iNode = 1, nDOF

                      U(iNode,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                        = U(iNode,iZ1,iZ2,iZ3,iZ_B0(4),iCR,iS)

                    END DO

                  END IF

                  ! --- Outer Boundary (Homogeneous) ---

                  DO iNode = 1, nDOF

                    U(iNode,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                      = U(iNode,iZ1,iZ2,iZ3,iZ_E0(4),iCR,iS)

                  END DO

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A45,I2.2)') &
        '', 'Invalid Boundary Condition for TwoMoment X3: ', bcZ(4)
      STOP

    END SELECT

  END SUBROUTINE ApplyBC_TwoMoment_X3


END MODULE TwoMoment_BoundaryConditionsModule
