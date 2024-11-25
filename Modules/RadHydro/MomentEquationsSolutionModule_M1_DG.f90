MODULE MomentEquationsSolutionModule_M1_DG

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nDOF, &
    nX, swX, nNodesX, &
    nE, nNodesE
  USE ReferenceElementModule, ONLY: &
    NodesX1, WeightsX1
  USE UtilitiesModule, ONLY: &
    NodeNumber
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, dL_X1
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    iCR_N, iCR_G1, iCR_G2, iCR_G3, nCR
  USE RiemannSolverModule, ONLY: &
    NumericalFlux_Radiation
  USE MomentEquationsUtilitiesModule, ONLY: &
    Eigenvalues, &
    AlphaP, &
    AlphaM, &
    Flux_X1, &
    GeometrySources, &
    FluxFactor, &
    EddingtonFactor

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  INTEGER,  DIMENSION(:,:,:,:), ALLOCATABLE :: NodeNumberTable
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: L_X1_Dn, L_X1_Up
  REAL(DP), DIMENSION(:,:),     ALLOCATABLE :: dL_X1_q

  PUBLIC :: ComputeRHS_M1_DG

  LOGICAL, PARAMETER :: DisplayTimers = .FALSE.
  REAL(DP) :: Timer_DIV
  REAL(DP) :: Timer_SRC

CONTAINS


  SUBROUTINE ComputeRHS_M1_DG( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    INTEGER, INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      U (1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
    REAL(DP), INTENT(out) :: &
      dU(1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,1:)

    CALL InitializeRHS

    CALL ComputeRHS_M1_DG_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    CALL ComputeRHS_M1_DG_GeometrySources &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    CALL FinalizeRHS

  END SUBROUTINE ComputeRHS_M1_DG


  SUBROUTINE ComputeRHS_M1_DG_X1( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      U (1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
    REAL(DP), INTENT(out) :: &
      dU(1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,1:)

    INTEGER :: iE, iX1, iX2, iX3, iS
    INTEGER :: iNodeE, iNodeX1, jNodeX1, iNodeX2, iNodeX3
    INTEGER :: iNode, jNode, iCR
    REAL(DP)                                  :: FF_q, FF_L, FF_R
    REAL(DP)                                  :: EF_q, EF_L, EF_R
    REAL(DP)                                  :: X1C, dX1
    REAL(DP)                                  :: Alpha, AlphaPls, AlphaMns
    REAL(DP), DIMENSION(1:nCR)                :: VolumeTerm 
    REAL(DP), DIMENSION(1:nCR)                :: Flux_L, Flux_R, Flux
    REAL(DP), DIMENSION(1:nCR)                :: uCR_L, uCR_R
    REAL(DP), DIMENSION(1:nCR)                :: Lambda_L, Lambda_R
    REAL(DP), DIMENSION(1:nX(1))              :: a_X1_L, a_X1_R
    REAL(DP), DIMENSION(1:nX(1))              :: b_X1_L, b_X1_R
    REAL(DP), DIMENSION(1:nNodesX(1))         :: x_q, w_q
    REAL(DP), DIMENSION(1:nNodesX(1),1:nX(1)) :: a_X1_q
    REAL(DP), DIMENSION(1:nNodesX(1),1:nX(1)) :: b_X1_q
    REAL(DP), DIMENSION(1:nDOF,1:nCR)         :: uCR_P, uCR_K, uCR_N
    REAL(DP), DIMENSION(1:nDOF,1:nCR)         :: Flux_X1_q

    CALL Timer_Start( Timer_DIV )

    x_q = NodesX1
    w_q = WeightsX1

    ! -- Precompute Metric Functions --

    DO iX1 = iX_B0(1), iX_E0(1)
      X1C = MeshX(1) % Center(iX1)
      dX1 = MeshX(1) % Width (iX1)
      a_X1_L(iX1) &
        = 1.0_DP!a( [ X1C - Half * dX1, Zero, Zero ] )
      a_X1_R(iX1) &
        = 1.0_DP!a( [ X1C + Half * dX1, Zero, Zero ] )
      b_X1_L(iX1) &
        = 1.0_DP!b( [ X1C - Half * dX1, Zero, Zero ] )
      b_X1_R(iX1) &
        = 1.0_DP!b( [ X1C + Half * dX1, Zero, Zero ] )
      DO iNodeX1 = 1, nNodesX(1)
        a_X1_q(iNodeX1,iX1) &
          = 1.0_DP!a( [ X1C + dX1 * x_q(iNodeX1), Zero, Zero ] )
        b_X1_q(iNodeX1,iX1) &
          = 1.0_DP!b( [ X1C + dX1 * x_q(iNodeX1), Zero, Zero ] )
      END DO
    END DO

    ! -- Compute Right-Hand Side for Moment Equations --

    DO iS = 1, nSpecies

      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          !$OMP PARALLEL DO PRIVATE &
          !$OMP&              ( iX1, iE, iNodeX1, iNodeX2, iNodeX3, iNodeE,  &
          !$OMP&                iNode, jNodeX1, jNode, uCR_P, uCR_K, uCR_N,  &
          !$OMP&                FF_q, EF_q, Flux_X1_q, VolumeTerm, dX1, iCR, &
          !$OMP&                uCR_L, uCR_R, AlphaPls, AlphaMns, Alpha,     &
          !$OMP&                FF_L, EF_L, Flux_L, Lambda_L, Flux,          &
          !$OMP&                FF_R, EF_R, Flux_R, Lambda_R )
          DO iX1 = iX_B0(1), iX_E0(1)
            dX1 = MeshX(1) % Width (iX1)
            DO iE = 1, nE

              uCR_P = U(1:nDOF,iE,iX1-1,iX2,iX3,1:nCR,iS) ! Previous Element
              uCR_K = U(1:nDOF,iE,iX1,  iX2,iX3,1:nCR,iS) ! This     Element
              uCR_N = U(1:nDOF,iE,iX1+1,iX2,iX3,1:nCR,iS) ! Next     Element

              DO iNode = 1, nDOF

                FF_q = FluxFactor &
                         ( uCR_K(iNode,iCR_N ), uCR_K(iNode,iCR_G1), &
                           uCR_K(iNode,iCR_G2), uCR_K(iNode,iCR_G3) )

                EF_q = EddingtonFactor( FF_q )

                Flux_X1_q(iNode,1:nCR) &
                  = Flux_X1 &
                      ( uCR_K(iNode,iCR_N ), uCR_K(iNode,iCR_G1), &
                        uCR_K(iNode,iCR_G2), uCR_K(iNode,iCR_G3), &
                        FF_q, EF_q )

              END DO

              dU(1:nDOF,iE,iX1,iX2,iX3,1:nCR,iS) = Zero

              DO iNodeX3 = 1, nNodesX(3)
                DO iNodeX2 = 1, nNodesX(2)
                  DO iNodeX1 = 1, nNodesX(1)
                    DO iNodeE = 1, nNodesE

                      iNode &
                        = NodeNumberTable &
                            ( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                      ! -----------------
                      ! -- Volume Term --
                      ! -----------------

                      VolumeTerm = Zero
                      DO jNodeX1 = 1, nNodesX(1)

                        jNode &
                          = NodeNumberTable &
                              ( iNodeE, jNodeX1, iNodeX2, iNodeX3 )

                        VolumeTerm(1:nCR) &
                          = VolumeTerm(1:nCR) &
                              + w_q(jNodeX1) &
                                  * a_X1_q(jNodeX1,iX1) &
                                  * b_X1_q(jNodeX1,iX1) &
                                  * Flux_X1_q(jNode,1:nCR) &
                                  * dL_X1_q(jNodeX1,iNodeX1)

                      END DO

                      dU(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                        = dU(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                            + VolumeTerm(1:nCR) &
                                / ( w_q(iNodeX1) * a_X1_q(iNodeX1,iX1) &
                                      * b_X1_q(iNodeX1,iX1) * dX1 )

                      ! -------------------
                      ! -- Surface Terms --
                      ! -------------------

                      ! -- Left Face --

                      ! -- Interpolate Left and Right State --

                      DO iCR = 1, nCR
                        uCR_L(iCR) = Zero
                        uCR_R(iCR) = Zero
                        DO jNodeX1 = 1, nNodesX(1)

                          jNode &
                            = NodeNumberTable &
                                ( iNodeE, jNodeX1, iNodeX2, iNodeX3 )

                          uCR_L(iCR) &
                            = uCR_L(iCR) &
                                + L_X1_Up(jNodeX1) * uCR_P(jNode,iCR)

                          uCR_R(iCR) &
                            = uCR_R(iCR) &
                                + L_X1_Dn(jNodeX1) * uCR_K(jNode,iCR)

                        END DO
                      END DO

                      ! -- Numerical Flux --

                      FF_L = FluxFactor &
                               ( uCR_L(iCR_N ), uCR_L(iCR_G1), &
                                 uCR_L(iCR_G2), uCR_L(iCR_G3) )

                      EF_L = EddingtonFactor( FF_L )

                      Flux_L &
                        = Flux_X1 &
                            ( uCR_L(iCR_N),  uCR_L(iCR_G1), &
                              uCR_L(iCR_G2), uCR_L(iCR_G3), FF_L, EF_L )

                      Lambda_L &
                        = Eigenvalues &
                            ( uCR_L(iCR_N ), uCR_L(iCR_G1), &
                              uCR_L(iCR_G2), uCR_L(iCR_G3), FF_L, EF_L )

                      FF_R = FluxFactor &
                               ( uCR_R(iCR_N ), uCR_R(iCR_G1), &
                                 uCR_R(iCR_G2), uCR_R(iCR_G3) )

                      EF_R = EddingtonFactor( FF_R )

                      Flux_R &
                        = Flux_X1 &
                            ( uCR_R(iCR_N),  uCR_R(iCR_G1), &
                              uCR_R(iCR_G2), uCR_R(iCR_G3), FF_R, EF_R )

                      Lambda_R &
                        = Eigenvalues &
                            ( uCR_R(iCR_N ), uCR_R(iCR_G1), &
                              uCR_R(iCR_G2), uCR_R(iCR_G3), FF_R, EF_R )

                      AlphaPls = AlphaP( Lambda_L, Lambda_R )
                      AlphaMns = AlphaM( Lambda_L, Lambda_R )
                      Alpha    = MAX( AlphaPls, AlphaMns )

                      Flux = NumericalFlux_Radiation &
                               ( uCR_L, uCR_R, Flux_L, Flux_R, &
                                 Alpha, AlphaPls, AlphaMns, One, nCR )

                      ! -- Contribution to Right-Hand Side --

                      dU(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                        = dU(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                            + a_X1_L(iX1) * b_X1_L(iX1) &
                                * Flux(1:nCR) * L_X1_Dn(iNodeX1) &
                                    / ( w_q(iNodeX1) * a_X1_q(iNodeX1,iX1) &
                                          * b_X1_q(iNodeX1,iX1) * dX1 )

                      ! -- Right Face --

                      ! -- Interpolate Left and Right State --

                      DO iCR = 1, nCR
                        uCR_L(iCR) = 0.0_DP
                        uCR_R(iCR) = 0.0_DP
                        DO jNodeX1 = 1, nNodesX(1)

                          jNode &
                            = NodeNumberTable &
                                ( iNodeE, jNodeX1, iNodeX2, iNodeX3 )

                          uCR_L(iCR) &
                            = uCR_L(iCR) &
                                + L_X1_Up(jNodeX1) * uCR_K(jNode,iCR)

                          uCR_R(iCR) &
                            = uCR_R(iCR) &
                                + L_X1_Dn(jNodeX1) * uCR_N(jNode,iCR)

                        END DO
                      END DO

                      ! -- Numerical Flux --

                      FF_L = FluxFactor &
                               ( uCR_L(iCR_N ), uCR_L(iCR_G1), &
                                 uCR_L(iCR_G2), uCR_L(iCR_G3) )

                      EF_L = EddingtonFactor( FF_L )

                      Flux_L &
                        = Flux_X1 &
                            ( uCR_L(iCR_N),  uCR_L(iCR_G1), &
                              uCR_L(iCR_G2), uCR_L(iCR_G3), FF_L, EF_L )

                      Lambda_L &
                        = Eigenvalues &
                            ( uCR_L(iCR_N ), uCR_L(iCR_G1), &
                              uCR_L(iCR_G2), uCR_L(iCR_G3), FF_L, EF_L )

                      FF_R = FluxFactor &
                               ( uCR_R(iCR_N ), uCR_R(iCR_G1), &
                                 uCR_R(iCR_G2), uCR_R(iCR_G3) )

                      EF_R = EddingtonFactor( FF_R )

                      Flux_R &
                        = Flux_X1 &
                            ( uCR_R(iCR_N),  uCR_R(iCR_G1), &
                              uCR_R(iCR_G2), uCR_R(iCR_G3), FF_R, EF_R )

                      Lambda_R &
                        = Eigenvalues &
                            ( uCR_R(iCR_N ), uCR_R(iCR_G1), &
                              uCR_R(iCR_G2), uCR_R(iCR_G3), FF_R, EF_R )

                      AlphaPls = AlphaP( Lambda_L, Lambda_R )
                      AlphaMns = AlphaM( Lambda_L, Lambda_R )
                      Alpha    = MAX( AlphaPls, AlphaMns )

                      Flux = NumericalFlux_Radiation &
                               ( uCR_L, uCR_R, Flux_L, Flux_R, &
                                 Alpha, AlphaPls, AlphaMns, One, nCR )

                      ! -- Contribution to Right-Hand Side --

                      dU(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                        = dU(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                            -  a_X1_R(iX1) * b_X1_R(iX1) &
                                 * Flux(1:nCR) * L_X1_Up(iNodeX1) &
                                     / ( w_q(iNodeX1) * a_X1_q(iNodeX1,iX1) &
                                           * b_X1_q(iNodeX1,iX1) * dX1 )

                    END DO ! iNodeE
                  END DO ! iNodeX1
                END DO ! iNodeX2
              END DO ! iNodeX3

            END DO ! iE
          END DO ! iX1
          !$OMP END PARALLEL DO
        END DO ! iX2
      END DO ! iX3

    END DO ! iS

    CALL Timer_Stop( Timer_DIV )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'ComputeDIV: ', Timer_DIV
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES16.6E2)') &
        '', 'SUM = ', SUM( dU )
      WRITE(*,*)

    END IF

  END SUBROUTINE ComputeRHS_M1_DG_X1


  SUBROUTINE ComputeRHS_M1_DG_GeometrySources &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      U (1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
    REAL(DP), INTENT(out) :: &
      dU(1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,1:)

    INTEGER                           :: iE, iX1, iX2, iX3, iS
    INTEGER                           :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER                           :: iNode
    REAL(DP)                          :: X1, X2, X3
    REAL(DP), DIMENSION(1:nDOF,1:nCR) :: uCR_K

    IF( TRIM( CoordinateSystem ) == 'CARTSEIAN' ) RETURN

    CALL Timer_Start( Timer_SRC )

    DO iS = 1, nSpecies

      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          !$OMP PARALLEL DO PRIVATE &
          !$OMP&              ( iX1, iE, iNodeX1, iNodeX2, iNodeX3, iNodeE, &
          !$OMP&                iNode, uCR_K, X1, X2, X3 )
          DO iX1 = iX_B0(1), iX_E0(1)
            DO iE = 1, nE

              uCR_K = U(1:nDOF,iE,iX1,iX2,iX3,1:nCR,iS) ! This Element

              DO iNodeX3 = 1, nNodesX(3)

                X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

                DO iNodeX2 = 1, nNodesX(2)

                  X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

                  DO iNodeX1 = 1, nNodesX(1)

                    X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                    DO iNodeE = 1, nNodesE

                      iNode = NodeNumberTable &
                                ( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                      dU(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                        = dU(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                            + GeometrySources &
                                ( uCR_K(iNode,iCR_N),  uCR_K(iNode,iCR_G1), &
                                  uCR_K(iNode,iCR_G2), uCR_K(iNode,iCR_G3), &
                                  [ X1, X2, X3 ] )

                    END DO
                  END DO
                END DO
              END DO

            END DO
          END DO
          !$OMP END PARALLEL DO
        END DO
      END DO

    END DO

    CALL Timer_Stop( Timer_SRC )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'ComputeSRC: ', Timer_SRC
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES16.6E2)') &
        '', 'SUM = ', SUM( dU )
      WRITE(*,*)

    END IF

  END SUBROUTINE ComputeRHS_M1_DG_GeometrySources


  SUBROUTINE InitializeRHS

    INTEGER :: iNode, iNodeE, iNodeX1, jNodeX1, iNodeX2, iNodeX3
    REAL(DP), DIMENSION(nNodesX(1)) :: x_q

    ALLOCATE( NodeNumberTable(nNodesE,nNodesX(1),nNodesX(2),nNodesX(3)) )

    iNode = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

             NodeNumberTable(iNodeE,iNodeX1,iNodeX2,iNodeX3) &
               = iNode

             iNode = iNode + 1

          END DO
        END DO
      END DO
    END DO

    ALLOCATE( L_X1_Dn(nNodesX(1)) )
    ALLOCATE( L_X1_Up(nNodesX(1)) )
    ALLOCATE( dL_X1_q(nNodesX(1),nNodesX(1)) )

    x_q = MeshX(1) % Nodes
    DO jNodeX1 = 1, nNodesX(1)
      L_X1_Dn(jNodeX1) &
        = L_X1(jNodeX1) % P( - Half )
      L_X1_Up(jNodeX1) &
        = L_X1(jNodeX1) % P( + Half )
      DO iNodeX1 = 1, nNodesX(1)
        dL_X1_q(iNodeX1,jNodeX1) &
          = dL_X1(jNodeX1) % P( x_q(iNodeX1) )
      END DO
    END DO

  END SUBROUTINE InitializeRHS


  SUBROUTINE FinalizeRHS

    DEALLOCATE( NodeNumberTable )

    DEALLOCATE( L_X1_Dn )
    DEALLOCATE( L_X1_Up )
    DEALLOCATE( dL_X1_q )

  END SUBROUTINE FinalizeRHS


  SUBROUTINE Timer_Start( Timer )

    REAL(DP) :: Timer

    IF( .NOT. DisplayTimers ) RETURN

    Timer = MPI_WTIME( )

  END SUBROUTINE Timer_Start


  SUBROUTINE Timer_Stop( Timer )

    REAL(DP) :: Timer

    IF( .NOT. DisplayTimers ) RETURN

    Timer = MPI_WTIME( ) - Timer

  END SUBROUTINE Timer_Stop


  SUBROUTINE Timer_Add( Timer, dT )

    REAL(DP) :: Timer, dT

    IF( .NOT. DisplayTimers ) RETURN    

    Timer = Timer + dT

  END SUBROUTINE Timer_Add


END MODULE MomentEquationsSolutionModule_M1_DG
