MODULE SlopeLimiterModule_Euler_GR

  USE KindModule, ONLY: &
    DP, Zero, One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX
  USE UtilitiesModule, ONLY: &
    MinModB
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, L_X2, L_X3
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Fluid, &
    MapModalToNodal_Fluid
  USE ReferenceElementModuleX, ONLY: &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3, &
    WeightsX_q, &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_E, &
    Shock

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializeSlopeLimiter
  PUBLIC :: FinalizeSlopeLimiter
  PUBLIC :: ApplySlopeLimiter_Euler_GR

  REAL(DP) :: wTime

  LOGICAL  :: UseTroubledCellIndicator, UseSlopeLimiter
  REAL(DP) :: BetaTVD, BetaTVB
  REAL(DP) :: LimiterThreshold
  REAL(DP) :: SlopeTolerance
  REAL(DP), ALLOCATABLE :: WeightsX_X1_P(:), WeightsX_X1_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X2_P(:), WeightsX_X2_N(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X3_P(:), WeightsX_X3_N(:)

CONTAINS


  SUBROUTINE InitializeSlopeLimiter &
    ( BetaTVD_Option, BetaTVB_Option, LimiterThreshold_Option, &
      SlopeTolerance_Option, UseSlopeLimiter_Option, &
      UseTroubledCellIndicator_Option )

    REAL(DP), INTENT(in), OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option
    REAL(DP), INTENT(in), OPTIONAL :: &
      LimiterThreshold_Option
    REAL(DP), INTENT(in), OPTIONAL :: &
      SlopeTolerance_Option
    LOGICAL, INTENT(in), OPTIONAL :: &
      UseTroubledCellIndicator_Option, UseSlopeLimiter_Option

    INTEGER  :: iNodeX1, iNodex2, iNodeX3, iNode
    INTEGER  :: jNodeX1, jNodex2, jNodeX3, jNode
    REAL(DP) :: WeightX

    BetaTVD = One
    IF( PRESENT( BetaTVD_Option ) ) &
      BetaTVD = BetaTVD_Option

    BetaTVB = Zero
    IF( PRESENT( BetaTVB_Option ) ) &
      BetaTVB = BetaTVB_Option

    LimiterThreshold = 0.05_DP
    IF( PRESENT( LimiterThreshold_Option ) ) &
      LimiterThreshold = LimiterThreshold_Option

    SlopeTolerance = 1.0d-3
    IF( PRESENT( SlopeTolerance_Option ) ) &
      SlopeTolerance = SlopeTolerance_Option

    UseSlopeLimiter = .TRUE.
    IF( PRESENT( UseSlopeLimiter_Option ) ) &
         UseSlopeLimiter = UseSlopeLimiter_Option

    UseTroubledCellIndicator = .TRUE.
    IF( PRESENT( UseTroubledCellIndicator_Option ) ) &
         UseTroubledCellIndicator = UseTroubledCellIndicator_Option

    WRITE(*,*)
    WRITE(*,'(A31)') '  INFO: InitializeSlopeLimiter:'
    WRITE(*,'(A31)') '  -----------------------------'
    WRITE(*,*)
    WRITE(*,'(A30,L1)'   ) '    UseSlopeLimiter: ' , &
      UseSlopeLimiter
    WRITE(*,'(A30,L1)'   ) '    UseTroubledCellIndicator: ' , &
      UseTroubledCellIndicator
    WRITE(*,*)
    WRITE(*,'(A24,F5.3)' ) '      BetaTVD:          ' , &
      BetaTVD
    WRITE(*,'(A24,F5.3)' ) '      BetaTVB:          ' , &
      BetaTVB
    WRITE(*,'(A24,F5.3)' ) '      LimiterThreshold: ' , &
      LimiterThreshold
    WRITE(*,'(A24,F5.3)' ) '      SlopeTolerance:   ' , &
      SlopeTolerance
    WRITE(*,*)

    ALLOCATE( WeightsX_X1_P(nDOFX), WeightsX_X1_N(nDOFX) )
    ALLOCATE( WeightsX_X2_P(nDOFX), WeightsX_X2_N(nDOFX) )
    ALLOCATE( WeightsX_X3_P(nDOFX), WeightsX_X3_N(nDOFX) )

    DO jNode = 1, nDOFX

      jNodeX1 = NodeNumberTableX(1,jNode)
      jNodeX2 = NodeNumberTableX(2,jNode)
      jNodeX3 = NodeNumberTableX(3,jNode)

      WeightsX_X1_P(jNode) = Zero
      WeightsX_X1_N(jNode) = Zero
      WeightsX_X2_P(jNode) = Zero
      WeightsX_X2_N(jNode) = Zero
      WeightsX_X3_P(jNode) = Zero
      WeightsX_X3_N(jNode) = Zero

      DO iNode = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNode)
        iNodeX2 = NodeNumberTableX(2,iNode)
        iNodeX3 = NodeNumberTableX(3,iNode)

        WeightX = WeightsX1  (iNodeX1) &
                  * WeightsX2(iNodeX2) &
                  * WeightsX3(iNodeX3)

        WeightsX_X1_P(jNode) &
          = WeightsX_X1_P(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) + One ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X1_N(jNode) &
          = WeightsX_X1_N(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) - One ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X2_P(jNode) &
          = WeightsX_X2_P(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) + One ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X2_N(jNode) &
          = WeightsX_X2_N(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) - One ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) ) )

        WeightsX_X3_P(jNode) &
          = WeightsX_X3_P(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) + One ) )

        WeightsX_X3_N(jNode) &
          = WeightsX_X3_N(jNode) &
              + WeightX &
                * ( L_X1  (jNodeX1) % P( NodesX1(iNodeX1) ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) - One ) )

      END DO

    END DO

  END SUBROUTINE InitializeSlopeLimiter


  SUBROUTINE FinalizeSlopeLimiter

    DEALLOCATE( WeightsX_X1_P, WeightsX_X1_N )
    DEALLOCATE( WeightsX_X2_P, WeightsX_X2_N )
    DEALLOCATE( WeightsX_X3_P, WeightsX_X3_N )

  END SUBROUTINE FinalizeSlopeLimiter


  SUBROUTINE ApplySlopeLimiter_Euler_GR &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: SlopeDifference
    REAL(DP) :: dU (nCF,nDimsX)
    REAL(DP) :: U_M(nCF,0:2*nDimsX,nDOFX)

    IF( .NOT. UseSlopeLimiter ) RETURN
    
    IF( nDOFX == 1 ) RETURN

    wTime = MPI_WTIME( )

    CALL DetectTroubledCells &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    DO iX3 = iX_B0(3), iX_E0(3)
      dX3 = MeshX(3) % Width(iX3)
      DO iX2 = iX_B0(2), iX_E0(2)
        dX2 = MeshX(2) % Width(iX2)
        DO iX1 = iX_B0(1), iX_E0(1)
          dX1 = MeshX(1) % Width(iX1)

          IF( Shock(iX1,iX2,iX3) < LimiterThreshold ) CYCLE

          WRITE(*,'(A,1x,I3,1x,I3,1x,I3,1x,E13.6)') &
            "iX1,iX2,iX3,slope = ", iX1, iX2, iX3, Shock(iX1,iX2,iX3)

          DO iCF = 1, nCF

            CALL MapNodalToModal_Fluid &
                   ( U(:,iX1,  iX2,iX3,iCF), U_M(iCF,0,:) )
            CALL MapNodalToModal_Fluid &
                   ( U(:,iX1-1,iX2,iX3,iCF), U_M(iCF,1,:) )
            CALL MapNodalToModal_Fluid &
                   ( U(:,iX1+1,iX2,iX3,iCF), U_M(iCF,2,:) )

            IF( nDimsX > 1 )THEN

              CALL MapNodalToModal_Fluid &
                     ( U(:,iX1,iX2-1,iX3,iCF), U_M(iCF,3,:) )
              CALL MapNodalToModal_Fluid &
                     ( U(:,iX1,iX2+1,iX3,iCF), U_M(iCF,4,:) )

            END IF

            IF( nDimsX > 2 )THEN

              CALL MapNodalToModal_Fluid &
                     ( U(:,iX1,iX2,iX3-1,iCF), U_M(iCF,5,:) )
              CALL MapNodalToModal_Fluid &
                     ( U(:,iX1,iX2,iX3+1,iCF), U_M(iCF,6,:) )

            END IF

          END DO

          dU(:,1) = MinModB( U_M(:,0,2), &
                             BetaTVD * ( U_M(:,0,1) - U_M(:,1,1) ), &
                             BetaTVD * ( U_M(:,2,1) - U_M(:,0,1) ), &
                             dX1, BetaTVB )

          IF( nDimsX > 1 )THEN

            dU(:,2) = MinModB( U_M(:,0,3), &
                               BetaTVD * ( U_M(:,0,1) - U_M(:,3,1) ), &
                               BetaTVD * ( U_M(:,4,1) - U_M(:,0,1) ), &
                               dX2, BetaTVB )

          END IF

          IF( nDimsX > 2 )THEN

            dU(:,3) = MinModB( U_M(:,0,4), &
                               BetaTVD * ( U_M(:,0,1) - U_M(:,5,1) ), &
                               BetaTVD * ( U_M(:,6,1) - U_M(:,0,1) ), &
                               dX3, BetaTVB )

          END IF

!          WRITE(*,*) "dU(:,1)    = ", dU(:,1)
!          WRITE(*,*), "U_M(:,0,2) = ", U_M(:,0,2)

          DO iCF = 1, nCF

            SlopeDifference &
              = ABS( U_M(iCF,0,2) - dU(iCF,1) ) &
                / MAXVAL( [ ABS( U_M(iCF,0,2) ), &
                            ABS( dU (iCF,1) ), SqrtTiny ] )

            IF( nDimsX > 1 )THEN

              SlopeDifference &
                = MAX( SlopeDifference, &
                       ABS( U_M(iCF,0,3) - dU(iCF,2) ) &
                       / MAXVAL( [ ABS( U_M(iCF,0,3) ), &
                                   ABS( dU (iCF,2)), SqrtTiny ] ) )

            END IF

            IF( nDimsX > 2 )THEN

              SlopeDifference &
                = MAX( SlopeDifference, &
                       ABS( U_M(iCF,0,4) - dU(iCF,3) ) &
                       / MAXVAL( [ ABS( U_M(iCF,0,4) ), &
                                   ABS( dU (iCF,3)), SqrtTiny ] ) )

            END IF

            WRITE(*,*) "iCF, SlopeDifference = ", iCF, SlopeDifference

            IF( SlopeDifference > SlopeTolerance )THEN

              U_M(iCF,0,2:) = Zero

              U_M(iCF,0,2) = dU(iCF,1)

              IF( nDimsX > 1 ) &
                U_M(iCF,0,3) = dU(iCF,2)

              IF( nDimsX > 2 ) &
                U_M(iCF,0,4) = dU(iCF,3)

            END IF

            CALL MapModalToNodal_Fluid &
                   ( U(:,iX1,iX2,iX3,iCF), U_M(iCF,0,:) )

          END DO

        END DO
      END DO
   END DO

    wTime = MPI_WTIME( ) - wTime

  END SUBROUTINE ApplySlopeLimiter_Euler_GR


  SUBROUTINE DetectTroubledCells &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &         
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: V_K (0:2*nDimsX)
    REAL(DP) :: U_K (0:2*nDimsX,nCF)
    REAL(DP) :: U_K0(0:2*nDimsX,nCF)

    Shock = Zero
    
    IF ( .NOT. UseTroubledCellIndicator ) RETURN

    ! --- Using troubled-cell indicator from Fu, Shu (2017) ---
    
    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          V_K(0) = DOT_PRODUCT &
                     ( WeightsX_q(:), G(:,iX1,  iX2,iX3,iGF_SqrtGm) )

          V_K(1) = DOT_PRODUCT &
!                     ( WeightsX_q(:), G(:,iX1-1,iX2,iX3,iGF_SqrtGm) )
                     ( WeightsX_X1_N(:), G(:,iX1-1,iX2,iX3,iGF_SqrtGm) )

          V_K(2) = DOT_PRODUCT &
!                     ( WeightsX_q(:), G(:,iX1+1,iX2,iX3,iGF_SqrtGm) )
                     ( WeightsX_X1_P(:), G(:,iX1+1,iX2,iX3,iGF_SqrtGm) )

          DO iCF = 1, nCF

            U_K(0,iCF) &
              = DOT_PRODUCT &
                  ( WeightsX_q(:), &
                    G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                      * U(:,iX1,iX2,iX3,iCF) ) / V_K(0)

            U_K(1,iCF) &
              = DOT_PRODUCT &
!                  ( WeightsX_q(:), &
                  ( WeightsX_X1_N(:), &
                    G(:,iX1-1,iX2,iX3,iGF_SqrtGm) &
                      * U(:,iX1-1,iX2,iX3,iCF) ) / V_K(1)

            U_K0(1,iCF) &
              = DOT_PRODUCT &
!                  ( WeightsX_X1_P(:), &
                  ( WeightsX_q(:), &
                    G(:,iX1-1,iX2,iX3,iGF_SqrtGm) &
                      * U(:,iX1-1,iX2,iX3,iCF) ) / V_K(0)

            U_K(2,iCF) &
              = DOT_PRODUCT &
!                  ( WeightsX_q(:), &
                  ( WeightsX_X1_P(:), &
                    G(:,iX1+1,iX2,iX3,iGF_SqrtGm) &
                      * U(:,iX1+1,iX2,iX3,iCF) ) / V_K(2)

            U_K0(2,iCF) &
              = DOT_PRODUCT &
!                  ( WeightsX_X1_N(:), &
                  ( WeightsX_q(:), &
                    G(:,iX1+1,iX2,iX3,iGF_SqrtGm) &
                      * U(:,iX1+1,iX2,iX3,iCF) ) / V_K(0)

          END DO

          IF( nDimsX > 1 )THEN

            V_K(3) = DOT_PRODUCT &
                       ( WeightsX_q(:), G(:,iX1,iX2-1,iX3,iGF_SqrtGm) )

            V_K(4) = DOT_PRODUCT &
                       ( WeightsX_q(:), G(:,iX1,iX2+1,iX3,iGF_SqrtGm) )

            DO iCF = 1, nCF

              U_K(3,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_q(:), &
                      G(:,iX1,iX2-1,iX3,iGF_SqrtGm) &
                        * U(:,iX1,iX2-1,iX3,iCF) ) / V_K(3)

              U_K0(3,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_X2_P(:), &
                      G(:,iX1,iX2-1,iX3,iGF_SqrtGm) &
                        * U(:,iX1,iX2-1,iX3,iCF) ) / V_K(0)

              U_K(4,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_q(:), &
                      G(:,iX1,iX2+1,iX3,iGF_SqrtGm) &
                        * U(:,iX1,iX2+1,iX3,iCF) ) / V_K(4)

              U_K0(4,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_X2_N(:), &
                      G(:,iX1,iX2+1,iX3,iGF_SqrtGm) &
                        * U(:,iX1,iX2+1,iX3,iCF) ) / V_K(0)

            END DO

          END IF

          IF( nDimsX > 2 )THEN

            V_K(5) = DOT_PRODUCT &
                       ( WeightsX_q(:), G(:,iX1,iX2,iX3-1,iGF_SqrtGm) )

            V_K(6) = DOT_PRODUCT &
                       ( WeightsX_q(:), G(:,iX1,iX2,iX3+1,iGF_SqrtGm) )

            DO iCF = 1, nCF

              U_K(5,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_q(:), &
                      G(:,iX1,iX2,iX3-1,iGF_SqrtGm) &
                        * U(:,iX1,iX2,iX3-1,iCF) ) / V_K(5)

              U_K0(5,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_X3_P(:), &
                      G(:,iX1,iX2,iX3-1,iGF_SqrtGm) &
                        * U(:,iX1,iX2,iX3-1,iCF) ) / V_K(0)

              U_K(6,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_q(:), &
                      G(:,iX1,iX2,iX3+1,iGF_SqrtGm) &
                        * U(:,iX1,iX2,iX3+1,iCF) ) / V_K(6)

              U_K0(6,iCF) &
                = DOT_PRODUCT &
                    ( WeightsX_X3_N(:), &
                      G(:,iX1,iX2,iX3+1,iGF_SqrtGm) &
                        * U(:,iX1,iX2,iX3+1,iCF) ) / V_K(0)

            END DO

          END IF

          Shock(iX1,iX2,iX3) &
            = SUM( ABS( U_K(0,iCF_D) - U_K0(1:2*nDimsX,iCF_D) ) ) &
                / MAXVAL( ABS( U_K(0:2*nDimsX,iCF_D) ) )

          Shock(iX1,iX2,iX3) &
            = MAX( Shock(iX1,iX2,iX3), &
                   SUM( ABS( U_K(0,iCF_E) - U_K0(1:2*nDimsX,iCF_E) ) ) &
                   / MAXVAL( ABS( U_K(0:2*nDimsX,iCF_E) ) ) )

        END DO
      END DO
    END DO

  END SUBROUTINE DetectTroubledCells


END MODULE SlopeLimiterModule_Euler_GR
