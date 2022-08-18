MODULE ScalarWave_BoundaryConditionsModule

  USE KindModule, ONLY: &
    DP
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ProgramHeaderModule, ONLY: &
    bcX, &
    swX, &
    nDOFX, &
    nNodesX
  USE ScalarFieldsModule, ONLY: &
    nSF, &
    iSF_u, &
    iSF_v
  USE ScalarWave_ErrorModule, ONLY: &
    DescribeError_ScalarWave

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_ScalarWave
  PUBLIC :: ApplyInnerBC_ScalarWave
  PUBLIC :: ApplyOuterBC_ScalarWave

  INTEGER, PARAMETER, PUBLIC :: iApplyBC_ScalarWave_Both  = 0
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_ScalarWave_Inner = 1
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_ScalarWave_Outer = 2
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_ScalarWave_None  = 3

  REAL(DP), PUBLIC :: ExpD
  REAL(DP), PUBLIC :: ExpE

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP ( ExpD, ExpE )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( ExpD, ExpE )
#endif


CONTAINS


  LOGICAL FUNCTION ApplyInnerBC_ScalarWave( iApplyBC )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iApplyBC

    ApplyInnerBC_ScalarWave = .FALSE.
    IF( iApplyBC .EQ. iApplyBC_ScalarWave_Inner .OR. &
        iApplyBC .EQ. iApplyBC_ScalarWave_Both ) &
    ApplyInnerBC_ScalarWave = .TRUE.

  END FUNCTION ApplyInnerBC_ScalarWave


  LOGICAL FUNCTION ApplyOuterBC_ScalarWave( iApplyBC )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iApplyBC

    ApplyOuterBC_ScalarWave = .FALSE.
    IF( iApplyBC .EQ. iApplyBC_ScalarWave_Outer .OR. &
        iApplyBC .EQ. iApplyBC_ScalarWave_Both ) &
    ApplyOuterBC_ScalarWave = .TRUE.

  END FUNCTION ApplyOuterBC_ScalarWave


  SUBROUTINE ApplyBoundaryConditions_ScalarWave &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    INTEGER :: iApplyBC(3)

    iApplyBC = iApplyBC_ScalarWave_Both

    IF( PRESENT( iApplyBC_Option ) ) &
      iApplyBC = iApplyBC_Option

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: U, iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  U, iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#endif

    CALL ApplyBC_ScalarWave_X1( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(1) )

    CALL ApplyBC_ScalarWave_X2( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(2) )

    CALL ApplyBC_ScalarWave_X3( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(3) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#endif

  END SUBROUTINE ApplyBoundaryConditions_ScalarWave


  SUBROUTINE ApplyBC_ScalarWave_X1( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iSF, iX1, iX2, iX3
    INTEGER :: iNode, iNodeX, iNodeX_0
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, jNodeX, jNodeX1
    REAL(DP) :: D_0, E_0, R_0, R_q

    SELECT CASE ( bcX(1) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#ifndef THORNADO_USE_AMREX

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNode = 1, nDOFX

          U(iNode,iX_B0(1)-iX1,iX2,iX3,iSF) &
            = U(iNode,iX_E0(1)-(iX1-1),iX2,iX3,iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNode = 1, nDOFX

          U(iNode,iX_E0(1)+iX1,iX2,iX3,iSF) &
            = U(iNode,iX_B0(1)+(iX1-1),iX2,iX3,iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

#endif

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNode = 1, nDOFX

          U(iNode,iX_B0(1)-iX1,iX2,iX3,iSF) &
            = U(iNode,iX_B0(1),iX2,iX3,iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNode = 1, nDOFX

          U(iNode,iX_E0(1)+iX1,iX2,iX3,iSF) &
            = U(iNode,iX_E0(1),iX2,iX3,iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

     ! --- Outer Boundary ---

     IF( ApplyOuterBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNode = 1, nDOFX

            U(iNode,iX_E0(1)+iX1,iX2,iX3,iSF) &
              = U(iNode,iX_E0(1),iX2,iX3,iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

      CALL DescribeError_ScalarWave( 01 )

    END SELECT

  END SUBROUTINE ApplyBC_ScalarWave_X1


  SUBROUTINE ApplyBC_ScalarWave_X2( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iSF, iX1, iX2, iX3
    INTEGER :: iNode, iNodeX, iNodeX1, iNodeX2, iNodeX3, jNodeX, jNodeX2

    SELECT CASE ( bcX(2) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#ifndef THORNADO_USE_AMREX

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNode = 1, nDOFX

          U(iNode,iX1,iX_B0(2)-iX2,iX3,iSF) &
            = U(iNode,iX1,iX_E0(2)-(iX2-1),iX3,iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNode = 1, nDOFX

          U(iNode,iX1,iX_E0(2)+iX2,iX3,iSF) &
            = U(iNode,iX1,iX_B0(2)+(iX2-1),iX3,iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

#endif

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNode = 1, nDOFX

          U(iNode,iX1,iX_B0(2)-iX2,iX3,iSF) &
            = U(iNode,iX1,iX_B0(2),iX3,iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNode = 1, nDOFX

          U(iNode,iX1,iX_E0(2)+iX2,iX3,iSF) &
            = U(iNode,iX1,iX_E0(2),iX3,iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
        !$OMP PRIVATE( iNodeX, jNodeX, jNodeX2 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRIVATE( iNodeX, jNodeX, jNodeX2 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5) &
        !$OMP PRIVATE( iNodeX, jNodeX, jNodeX2 )
#endif
        DO iSF = 1, nSF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNode = 1, nDOFX

          U(iNode,iX1,iX_E0(2)+iX2,iX3,iSF) &
            = U(iNode,iX1,iX_E0(2),iX3,iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

      CALL DescribeError_ScalarWave( 02 )

    END SELECT

  END SUBROUTINE ApplyBC_ScalarWave_X2


  SUBROUTINE ApplyBC_ScalarWave_X3( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iSF, iX1, iX2, iX3
    INTEGER :: iNode, iNodeX, iNodeX1, iNodeX2, iNodeX3, jNodeX, jNodeX3

    SELECT CASE ( bcX(3) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#ifndef THORNADO_USE_AMREX

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNode = 1, nDOFX

          U(iNode,iX1,iX2,iX_B0(3)-iX3,iSF) &
            = U(iNode,iX1,iX2,iX_E0(3)-(iX3-1),iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNode = 1, nDOFX

            U(iNode,iX1,iX2,iX_E0(3)+iX3,iSF) &
              = U(iNode,iX1,iX2,iX_B0(3)+(iX3-1),iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

#endif

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNode = 1, nDOFX

          U(iNode,iX1,iX2,iX_B0(3)-iX3,iSF) &
            = U(iNode,iX1,iX2,iX_B0(3),iSF)


        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNode = 1, nDOFX

          U(iNode,iX1,iX2,iX_E0(3)+iX3,iSF) &
            = U(iNode,iX1,iX2,iX_E0(3),iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_ScalarWave( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iSF = 1, nSF
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNode = 1, nDOFX

          U(iNode,iX1,iX2,iX_E0(3)+iX3,iSF) &
            = U(iNode,iX1,iX2,iX_E0(3),iSF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

      CALL DescribeError_ScalarWave( 03 )

    END SELECT

  END SUBROUTINE ApplyBC_ScalarWave_X3


END MODULE ScalarWave_BoundaryConditionsModule
