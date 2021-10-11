MODULE Euler_BoundaryConditionsModule

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
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_BoundaryConditions, &
    Timer_Euler_BC_ApplyBC, &
    Timer_Euler_BC_CopyIn, &
    Timer_Euler_BC_CopyOut
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_Euler
  PUBLIC :: ApplyInnerBC_Euler
  PUBLIC :: ApplyOuterBC_Euler

  INTEGER, PARAMETER, PUBLIC :: iApplyBC_Euler_Both  = 0
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_Euler_Inner = 1
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_Euler_Outer = 2
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_Euler_None  = 3

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


  LOGICAL FUNCTION ApplyInnerBC_Euler( iApplyBC )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iApplyBC

    ApplyInnerBC_Euler = .FALSE.
    IF( iApplyBC .EQ. iApplyBC_Euler_Inner .OR. &
        iApplyBC .EQ. iApplyBC_Euler_Both ) &
    ApplyInnerBC_Euler = .TRUE.

  END FUNCTION ApplyInnerBC_Euler


  LOGICAL FUNCTION ApplyOuterBC_Euler( iApplyBC )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iApplyBC

    ApplyOuterBC_Euler = .FALSE.
    IF( iApplyBC .EQ. iApplyBC_Euler_Outer .OR. &
        iApplyBC .EQ. iApplyBC_Euler_Both ) &
    ApplyOuterBC_Euler = .TRUE.

  END FUNCTION ApplyOuterBC_Euler


  SUBROUTINE ApplyBoundaryConditions_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    INTEGER :: iApplyBC(3)

    CALL TimersStart_Euler( Timer_Euler_BoundaryConditions )

    iApplyBC = iApplyBC_Euler_Both

    IF( PRESENT( iApplyBC_Option ) ) &
      iApplyBC = iApplyBC_Option

    CALL TimersStart_Euler( Timer_Euler_BC_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: U, iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  U, iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#endif

    CALL TimersStop_Euler( Timer_Euler_BC_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_BC_ApplyBC )

    CALL ApplyBC_Euler_X1( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(1) )

    CALL ApplyBC_Euler_X2( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(2) )

    CALL ApplyBC_Euler_X3( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(3) )

    CALL TimersStop_Euler( Timer_Euler_BC_ApplyBC )

    CALL TimersStart_Euler( Timer_Euler_BC_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, iApplyBC )
#endif

    CALL TimersStop_Euler( Timer_Euler_BC_CopyOut )

    CALL TimersStop_Euler( Timer_Euler_BoundaryConditions )

  END SUBROUTINE ApplyBoundaryConditions_Euler


  SUBROUTINE ApplyBC_Euler_X1( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iCF, iX1, iX2, iX3
    INTEGER  :: iNX, iNX_0
    INTEGER  :: iNX1, iNX2, iNX3, jNX, jNX1
    REAL(DP) :: D_0, E_0, R_0, R_q

    SELECT CASE ( bcX(1) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#ifndef USE_AMREX_TRUE

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

          U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF) &
            = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

          U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF) &
            = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

#endif

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

          U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF) &
            = U(iNX,iX_B0(1),iX2,iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

          U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF) &
            = U(iNX,iX_E0(1),iX2,iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 3 ) ! Reflecting

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX1 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX1 = ( nNodesX(1) - iNX1 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( jNX1, iNX2, iNX3 )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_D) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCF_D)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S1) &
              = - U(jNX,iX_B0(1),iX2,iX3,iCF_S1)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S2) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCF_S2)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S3) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCF_S3)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_E) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCF_E)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_Ne) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX1 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX1 = ( nNodesX(1) - iNX1 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( jNX1, iNX2, iNX3 )

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_D) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCF_D)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_S1) &
              = - U(jNX,iX_E0(1),iX2,iX3,iCF_S1)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_S2) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCF_S2)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_S3) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCF_S3)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_E) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCF_E)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_Ne) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 30 ) ! Reflecting (Inner), Zero (Outer)

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX1 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

            DO iNX3 = 1, nNodesX(3)
            DO iNX2 = 1, nNodesX(2)
            DO iNX1 = 1, nNodesX(1)

              jNX1 = ( nNodesX(1) - iNX1 ) + 1

              iNX = NodeNumberX( iNX1, iNX2, iNX3 )
              jNX = NodeNumberX( jNX1, iNX2, iNX3 )

              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_D) &
                = + U(jNX,iX_B0(1),iX2,iX3,iCF_D)
              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S1) &
                = - U(jNX,iX_B0(1),iX2,iX3,iCF_S1)
              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S2) &
                = + U(jNX,iX_B0(1),iX2,iX3,iCF_S2)
              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S3) &
                = + U(jNX,iX_B0(1),iX2,iX3,iCF_S3)
              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_E) &
                = + U(jNX,iX_B0(1),iX2,iX3,iCF_E)
              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_Ne) &
                = + U(jNX,iX_B0(1),iX2,iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 31 ) ! Reflecting (Inner), Homogeneous (Outer)

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX1 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX1 = ( nNodesX(1) - iNX1 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( jNX1, iNX2, iNX3 )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_D) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCF_D)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S1) &
              = - U(jNX,iX_B0(1),iX2,iX3,iCF_S1)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S2) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCF_S2)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S3) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCF_S3)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_E) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCF_E)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_Ne) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

     ! --- Outer Boundary ---

     IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF) &
              = U(iNX,iX_E0(1),iX2,iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 11 ) ! Custom BCs for Accretion Problem

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

        ASSOCIATE( X1_C  => MeshX(1) % Center, &
                   dX1   => MeshX(1) % Width,  &
                   eta_q => MeshX(1) % Nodes )

        R_0 = X1_C(1) + dX1(1) * eta_q(1)

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6), &
        !$OMP PRIVATE( iNX, iNX_0, D_0, E_0, R_q ) &
        !$OMP FIRSTPRIVATE( R_0 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC COPYIN( X1_C, dX1, eta_q ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX ) &
        !$ACC PRIVATE( iNX, iNX_0, D_0, E_0, R_q ) &
        !$ACC FIRSTPRIVATE( R_0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, iNX_0, D_0, E_0, R_q ) &
        !$OMP FIRSTPRIVATE( R_0 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            iNX   = NodeNumberX( iNX1, iNX2, iNX3 )
            iNX_0 = NodeNumberX( 1,       iNX2, iNX3 )

            D_0 = U(iNX_0,1,iX2,iX3,iCF_D)
            E_0 = U(iNX_0,1,iX2,iX3,iCF_E)

            R_q = NodeCoordinate &
                    ( X1_C( iX_B0(1) - iX1 ), dX1( iX_B0(1) - iX1 ), &
                      eta_q( iNX1 ) )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_D) &
              = D_0 * ( R_0 / R_q )**3

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_E) &
              = E_0 * ( R_0 / R_q )**4

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

        END ASSOCIATE

      END IF

    CASE ( 100 ) ! Custom BCs for Relativistic SAS Problem

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

        ASSOCIATE( X1_C  => MeshX(1) % Center, &
                   dX1   => MeshX(1) % Width,  &
                   eta_q => MeshX(1) % Nodes )

        R_0 = X1_C(1) + dX1(1) * eta_q(1)

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6), &
        !$OMP PRIVATE( iNX, iNX_0, D_0, E_0, R_q ) &
        !$OMP FIRSTPRIVATE( R_0 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC COPYIN( X1_C, dX1, eta_q ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX, ExpD, ExpE ) &
        !$ACC PRIVATE( iNX, iNX_0, D_0, E_0, R_q ) &
        !$ACC FIRSTPRIVATE( R_0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, iNX_0, D_0, E_0, R_q ) &
        !$OMP FIRSTPRIVATE( R_0 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            iNX   = NodeNumberX( iNX1, iNX2, iNX3 )
            iNX_0 = NodeNumberX( 1,       iNX2, iNX3 )

            D_0 = U(iNX_0,1,iX2,iX3,iCF_D)
            E_0 = U(iNX_0,1,iX2,iX3,iCF_E)

            R_q = NodeCoordinate &
                    ( X1_C( iX_B0(1) - iX1 ), dX1( iX_B0(1) - iX1 ), &
                      eta_q( iNX1 ) )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_D) &
              = D_0 * ( R_0 / R_q )**( ExpD )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_E) &
              = E_0 * ( R_0 / R_q )**( ExpE )

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

        END ASSOCIATE

      END IF

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

     ! --- Outer Boundary ---

     IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF) &
              = U(iNX,iX_E0(1),iX2,iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

      CALL DescribeError_Euler( 05 )

    END SELECT

  END SUBROUTINE ApplyBC_Euler_X1


  SUBROUTINE ApplyBC_Euler_X2( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iCF, iX1, iX2, iX3
    INTEGER :: iNX, iNX1, iNX2, iNX3, jNX, jNX2

    SELECT CASE ( bcX(2) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#ifndef USE_AMREX_TRUE

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF) &
            = U(iNX,iX1,iX_E0(2)-(iX2-1),iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF) &
            = U(iNX,iX1,iX_B0(2)+(iX2-1),iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

#endif

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF) &
            = U(iNX,iX1,iX_B0(2),iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF) &
            = U(iNX,iX1,iX_E0(2),iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 3 ) ! Reflecting

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX2 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_D) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_D)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S1) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_S1)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S2) &
              = - U(jNX,iX1,iX_B0(2),iX3,iCF_S2)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S3) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_S3)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_E) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_E)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_Ne) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX2 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_D) &
              = + U(jNX,iX1,iX_E0(2),iX3,iCF_D)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_S1) &
              = + U(jNX,iX1,iX_E0(2),iX3,iCF_S1)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_S2) &
              = - U(jNX,iX1,iX_E0(2),iX3,iCF_S2)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_S3) &
              = + U(jNX,iX1,iX_E0(2),iX3,iCF_S3)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_E) &
              = + U(jNX,iX1,iX_E0(2),iX3,iCF_E)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_Ne) &
              = + U(jNX,iX1,iX_E0(2),iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 30 ) ! Reflecting (Inner), Zero (Outer)

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX2 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_D) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_D)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S1) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_S1)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S2) &
              = - U(jNX,iX1,iX_B0(2),iX3,iCF_S2)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S3) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_S3)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_E) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_E)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_Ne) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 31 ) ! Reflecting (Inner), Homogeneous (Outer)

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX2 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_D) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_D)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S1) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_S1)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S2) &
              = - U(jNX,iX1,iX_B0(2),iX3,iCF_S2)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S3) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_S3)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_E) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_E)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_Ne) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRIVATE( iNX, jNX, jNX2 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF) &
            = U(iNX,iX1,iX_E0(2),iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRIVATE( iNX, jNX, jNX2 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF) &
            = U(iNX,iX1,iX_E0(2),iX3,iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

      CALL DescribeError_Euler( 06 )

    END SELECT

  END SUBROUTINE ApplyBC_Euler_X2


  SUBROUTINE ApplyBC_Euler_X3( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iCF, iX1, iX2, iX3
    INTEGER :: iNX

    SELECT CASE ( bcX(3) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#ifndef USE_AMREX_TRUE

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF) &
            = U(iNX,iX1,iX2,iX_E0(3)-(iX3-1),iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF) &
              = U(iNX,iX1,iX2,iX_B0(3)+(iX3-1),iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

#endif

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF) &
            = U(iNX,iX1,iX2,iX_B0(3),iCF)


        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF) &
            = U(iNX,iX1,iX2,iX_E0(3),iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_Euler( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0, swX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF) &
            = U(iNX,iX1,iX2,iX_E0(3),iCF)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

      CALL DescribeError_Euler( 07 )

    END SELECT

  END SUBROUTINE ApplyBC_Euler_X3


END MODULE Euler_BoundaryConditionsModule
