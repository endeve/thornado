MODULE TwoMoment_BoundaryConditionsModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    bcZ, swZ, nE, &
    nNodesZ, nDOF, nDOFE, nDOFX
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable4D
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR,iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE MeshModule,              ONLY: &
    MeshType,    &
    CreateMesh,  &
    DestroyMesh, &
    MeshX,          &
    MeshE,          &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE GeometryFieldsModuleE, ONLY: &
    nGE, iGE_Ep2
  USE ReferenceElementModule, ONLY: &
    nDOF_E, &
    nDOF_X1, &
    nDOF_X2, &
    nDOF_X3, &
    Weights_Q


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
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: U, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( U, iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#endif

    CALL ApplyBC_TwoMoment_E &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

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


  SUBROUTINE ApplyBC_TwoMoment_E( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

    INTEGER, PARAMETER :: nE_LeastSquares = 5
    INTEGER  :: i, iZ1, iZ2, iZ3, iZ4, iCR, iS, iNodeX, iNodeE, iNodeZ
    REAL(DP) :: E_K(2), E_R
    REAL(DP) :: &
      N_K(2, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          1:nSpecies), &
      a_K(iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          1:nSpecies), &
      b_K(iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          1:nSpecies)
    REAL(DP) :: L_lnN( nE_LeastSquares, nDOFE, nDOFX)
    REAL(DP) :: L_N( nE_LeastSquares, nDOFE, nDOFX)
    REAL(DP) :: N(nE_LeastSquares, nDOFE), lnN(nE_LeastSquares, nDOFE)
    REAL(DP) :: E(nE_LeastSquares, nDOFE), A_T(nDOFX), B_T(nDOFX)
    
    SELECT CASE ( bcZ(1) )

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
      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iZ1 = 1, swZ(1)

        DO iNodeZ = 1, nDOF

          ! --- Inner Boundary ---

          U(iNodeZ,iZ_B0(1)-iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            = U(iNodeZ,iZ_E0(1)-(iZ1-1),iZ2,iZ3,iZ4,iCR,iS)

          ! --- Outer Boundary ---

          U(iNodeZ,iZ_E0(1)+iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            = U(iNodeZ,iZ_B0(1)+(iZ1-1),iZ2,iZ3,iZ4,iCR,iS)

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
      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iZ1 = 1, swZ(1)

        DO iNodeZ = 1, nDOF

          ! --- Inner Boundary ---

          U(iNodeZ,iZ_B0(1)-iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            = U(iNodeZ,iZ_B0(1),iZ2,iZ3,iZ4,iCR,iS)

          ! --- Outer Boundary ---

          U(iNodeZ,iZ_E0(1)+iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            = U(iNodeZ,iZ_E0(1),iZ2,iZ3,iZ4,iCR,iS)

        END DO

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

    CASE ( 10 ) ! Custom

      ! --- Compute Cell Averaged Density in Last Two Elements ---

      ASSOCIATE &
        ( CenterE => MeshE % Center, &
          WidthE  => MeshE % Width )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET ENTER DATA &
      !$OMP MAP( alloc: N_K, a_K, b_K ) &
      !$OMP MAP( to: CenterE, WidthE )
#elif defined(THORNADO_OACC)
      !$ACC ENTER DATA &
      !$ACC CREATE( N_K, a_K, b_K ) &
      !$ACC COPYIN( CenterE, WidthE )
#endif

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, N_K, Weights_q )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iS  = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iZ1 = 1, 2

        N_K(iZ1,iZ2,iZ3,iZ4,iS) = Zero

        DO iNodeZ = 1, nDOF

          N_K(iZ1,iZ2,iZ3,iZ4,iS) &
            = N_K(iZ1,iZ2,iZ3,iZ4,iS) &
                + Weights_q(iNodeZ) &
                    * U(iNodeZ,iZ_E0(1)-(2-iZ1),iZ2,iZ3,iZ4,iCR_N,iS)

        END DO

      END DO
      END DO
      END DO
      END DO
      END DO

      ! --- Compute Exponential Extrapolation Coefficients ---

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( iZ_B0, iZ_E0, N_K, a_K, b_K, CenterE )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iS  = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)

        a_K(iZ2,iZ3,iZ4,iS) &
          = ( CenterE(iZ_E0(1)-1) * LOG( N_K(2,iZ2,iZ3,iZ4,iS) ) &
                - CenterE(iZ_E0(1)) * LOG( N_K(1,iZ2,iZ3,iZ4,iS) ) ) &
            / ( CenterE(iZ_E0(1)) - CenterE(iZ_E0(1)-1) )

        b_K(iZ2,iZ3,iZ4,iS) &
          = ( LOG( N_K(1,iZ2,iZ3,iZ4,iS) ) &
                - LOG( N_K(2,iZ2,iZ3,iZ4,iS) ) ) &
            / ( CenterE(iZ_E0(1)) - CenterE(iZ_E0(1)-1) )

      END DO
      END DO
      END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
      !$OMP PRIVATE( E_R )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ ) &
      !$ACC PRIVATE( E_R )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(7) &
      !$OMP PRIVATE( E_R )
#endif
      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iZ1 = 1, swZ(1)

        DO iNodeZ = 1, nDOF

          ! --- Inner Boundary ---

          U(iNodeZ,iZ_B0(1)-iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            = U(iNodeZ,iZ_B0(1),iZ2,iZ3,iZ4,iCR,iS)

          ! --- Outer Boundary ---

          IF( iCR == iCR_N )THEN

            E_R = CenterE(iZ_E0(1)) + Half * WidthE(iZ_E0(1))
            U(iNodeZ,iZ_E0(1)+iZ1,iZ2,iZ3,iZ4,iCR,iS) &
              = One / ( EXP(   a_K(iZ2,iZ3,iZ4,iS) &
                             + b_K(iZ2,iZ3,iZ4,iS) * E_R ) + One )

          ELSE

            U(iNodeZ,iZ_E0(1)+iZ1,iZ2,iZ3,iZ4,iCR,iS) &
              = Zero

          END IF

        END DO

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: N_K, a_K, b_K, CenterE, WidthE )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA &
      !$ACC DELETE( N_K, a_K, b_K, CenterE, WidthE )
#endif

      END ASSOCIATE

    CASE ( 11 ) ! Custom

      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iZ1 = 1, swZ(1)

        DO iNodeZ = 1, nDOF

          ! --- Inner Boundary ---

          U(iNodeZ,iZ_B0(1)-iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            = U(iNodeZ,iZ_B0(1),iZ2,iZ3,iZ4,iCR,iS)

          ! --- Outer Boundary ---

          IF( iCR == iCR_N )THEN

            U(iNodeZ,iZ_E0(1)+iZ1,iZ2,iZ3,iZ4,iCR,iS) = SqrtTiny

          ELSE

            U(iNodeZ,iZ_E0(1)+iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero

          END IF

        END DO

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

    CASE ( 22 ) ! Custom
      i = 1


      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iZ1 = iZ_B0(1), iZ_E0(1)

        IF ( iZ1 .GE. (nE - nE_LeastSquares + 1) ) THEN
         
          DO iNodeZ = 1, nDOF

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1


            L_N( i, iNodeE, iNodeX ) =  U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)
            
            IF(L_N(i,iNodeE,iNodeX) .EQ. 0.0_DP) THEN

              L_lnN( i, iNodeE, iNodeX ) = 0.0_DP

            ELSE

              L_lnN( i, iNodeE, iNodeX ) = LOG( U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) )

            END IF

            E(i,iNodeE) = NodeCoordinate( MeshE, iZ1, iNodeE )

          END DO

          i = i + 1
        END IF

      END DO 

      i = 1

      DO iNodeX = 1, nDOFX

        N(:,:) = L_N(:,:, iNodeX)
        
        lnN(:,:) = L_lnN(:,:,iNodeX)
        IF( ALL(N .EQ. 0.0_DP) )THEN


          A_T(iNodeX) = 0.0_DP

          B_T(iNodeX) = 0.0_DP

        ELSE

          B_T(iNodeX) &
            = ( SUM(N)*SUM(E*N*lnN)-SUM(E*N)*SUM(N*lnN) ) / (SUM(N)*SUM(E**2*N)-(SUM(E*N))**2)

          A_T(iNodeX) &
            = ( SUM(E**2*N)*SUM(N*lnN)-SUM(E*N)*SUM(E*N*lnN) ) / (SUM(N)*SUM(E**2*N)-(SUM(E*N))**2)
        
          A_T(iNodeX) = EXP(A_T(iNodeX) )

        END IF 


      END DO
      DO iZ1 = 1, swZ(1)


        DO iNodeZ = 1, nDOF

          iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1


          iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

       !   ! --- Inner Boundary ---
          U(iNodeZ,iZ_B0(1)-iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            = U(iNodeZ,iZ_B0(1),iZ2,iZ3,iZ4,iCR,iS)

          ! --- Outer Boundary ---

          E_R = MeshE % Center(iZ_E0(1)+iZ1) - 0.5_DP * MeshE % Width(iZ_E0(1)+iZ1)

          U(iNodeZ,iZ_E0(1)+iZ1,iZ2,iZ3,iZ4,iCR,iS) = A_T(iNodeX) * EXP( B_T(iNodeX) * E_R )

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
        '', 'Invalid Boundary Condition for TwoMoment E: ', bcZ(1)
      STOP
    END SELECT

  END SUBROUTINE ApplyBC_TwoMoment_E


  SUBROUTINE ApplyBC_TwoMoment_X1( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

    INTEGER :: iNode, iS, iCR, iZ1, iZ2, iZ3, iZ4
    INTEGER :: iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4
    INTEGER :: jNodeZ2, iNodeZ, jNodeZ, iNodeE
    REAL(DP):: E, Mu_0

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

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(7)
#endif
      DO iS = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = 1, swZ(2)
      DO iZ1 = iZ_B0(1), iZ_E0(1)

        ! --- Inner: No Boundary Condition ---

        ! --- Outer: Homogeneous ---

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

    CASE ( 32,31,30 ) ! Reflecting (Inner), Homogeneous (Outer)

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
                      = MAX( U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR,iS), 0.0_DP )

                  END DO

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 22 ) ! Custom Boundary Conditions for radiating inner and outflow outer
      Mu_0 = 0.9_DP
      DO iS = 1, nSpecies
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = 1, swZ(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)
                  DO iNode = 1, nDOF

                    iNodeE = MOD( (iNode-1)        , nDOFE ) + 1


                   E = NodeCoordinate( MeshE, iZ1, iNodeE )
                    ! --- Inner Boundary ---
                    U(iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR_N,iS)  &
                      != 0.5_DP * ( 1.0_DP - Mu_0 ) / ( EXP( E / 3.0_DP - 3.0_DP ) + 1.0_DP )
                      =   1.0_DP / ( EXP( E / 3.0_DP - 3.0_DP ) + 1.0_DP )
                      ! =   1.0_DP / ( EXP( E )- 1.0_DP )
                    U( iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR_G1,iS)  &
                      != 0.5_DP * ( 1.0_DP + Mu_0 ) * U(iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR_N,iS) 
                      = 0.99_DP * U(iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR_N,iS)
                    U( iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR_G2,iS) &
                      = 0.0_DP

                    U( iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR_G3,iS) &
                      = 0.0_DP

                    !print*,iNode, iZ1, iZ_B0(2)-iZ2, U(iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR_N,iS)


                    ! --- Outer Boundary ---

                    U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR_N,iS) &
                      = U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR_N,iS)

                    U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR_G1,iS) &
                      = U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR_G1,iS)

                    U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR_G2,iS) &
                      = U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR_G2,iS)

                    U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR_G3,iS) &
                      = U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR_G3,iS)
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
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

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

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(7)
#endif
      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = 1, swZ(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iZ1 = iZ_B0(1), iZ_E0(1)

        ! --- Inner: No Boundary Condition ---

        ! --- Outer: Homogeneous ---

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
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

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

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( U, iZ_B0, iZ_E0, swZ )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(7)
#endif
      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ4 = 1, swZ(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iZ1 = iZ_B0(1), iZ_E0(1)

        ! --- Inner: No Boundary Condition ---

        ! --- Outer: Homogeneous ---

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

  SUBROUTINE CheckRealizability(N, G1, V1,iZ2)

    REAL(DP), INTENT(in) ::  N, G1, V1
    INTEGER, INTENT(in) :: iZ2

    REAL(DP) :: gamma1, W

    W = 1.0_DP/SQRT(1.0_DP - V1**2)
    gamma1 = N - G1 / W
   
    IF(gamma1 .LT. 10d-26)THEN
    print*, iZ2, V1, N, G1, gamma1
    END IF
  END SUBROUTINE CheckRealizability



END MODULE TwoMoment_BoundaryConditionsModule
