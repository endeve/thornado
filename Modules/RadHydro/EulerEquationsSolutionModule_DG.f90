MODULE EulerEquationsSolutionModule_DG

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE ReferenceElementModuleX, ONLY: &
    NodesX1, WeightsX1
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, dL_X1
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    uGF, nGF
  USE FluidFieldsModule, ONLY: &
    rhsCF, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, nCF, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, nPF, &
    uAF, iAF_P, iAF_T,  iAF_Ye, iAF_S,  iAF_E, iAF_Gm, iAF_Cs, nAF
  USE EquationOfStateModule, ONLY: &
    ComputeAuxiliary_Fluid
  USE RiemannSolverModule, ONLY: &
    NumericalFlux_Fluid
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputePrimitive, &
    Primitive, &
    Eigenvalues, &
    AlphaP, &
    AlphaM, &
    AlphaC, &
    Flux_X1, &
    GeometrySources, &
    ComputeGeometrySources_Gravity

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  INTEGER,  DIMENSION(:,:,:), ALLOCATABLE :: NodeNumberTableX
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: L_X1_Dn, L_X1_Up
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: dL_X1_q

  PUBLIC :: ComputeRHS_Euler_DG

  LOGICAL, PARAMETER :: DisplayTimers = .FALSE.
  REAL(DP) :: Timer_RHS
  REAL(DP) :: Timer_GEO
  REAL(DP) :: Timer_AUX
  REAL(DP) :: Timer_INIT
  REAL(DP) :: Timer_FIN

CONTAINS


  SUBROUTINE ComputeRHS_Euler_DG( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iX1, iX2, iX3

    CALL Timer_Start( Timer_AUX )

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        !$OMP PARALLEL DO PRIVATE ( iX1 )
        DO iX1 = iX_Begin(1), iX_End(1)

          CALL ComputePrimitive &
                 ( uCF(:,iX1,iX2,iX3,1:nCF), uPF(:,iX1,iX2,iX3,1:nPF) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uPF(:,iX1,iX2,iX3,iPF_E),  &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P),  &
                   uAF(:,iX1,iX2,iX3,iAF_T),  uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_S),  uAF(:,iX1,iX2,iX3,iAF_E),  &
                   uAF(:,iX1,iX2,iX3,iAF_Gm), uAF(:,iX1,iX2,iX3,iAF_Cs) )

        END DO
        !$OMP END PARALLEL DO
      END DO
    END DO

    CALL Timer_Stop( Timer_AUX )

    CALL Timer_Start( Timer_INIT )

    CALL InitializeRHS

    CALL Timer_Stop( Timer_INIT )

    CALL Timer_Start( Timer_RHS )

    CALL ComputeRHS_Euler_DG_X1( iX_Begin, iX_End )

    CALL Timer_Stop( Timer_RHS )

    CALL Timer_Start( Timer_GEO )

    CALL ComputeRHS_Euler_DG_GeometrySources( iX_Begin, iX_End )

    CALL Timer_Stop( Timer_GEO )

    CALL Timer_Start( Timer_FIN )

    CALL FinalizeRHS

    CALL Timer_Stop( Timer_FIN )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'ComputeRHS: ', Timer_RHS
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'ComputeGEO: ', Timer_GEO
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'ComputeAUX: ', Timer_AUX
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'Initialize: ', Timer_INIT
      WRITE(*,'(A4,A16,ES10.4E2)') &
        '', 'Finalize: ', Timer_FIN
      WRITE(*,*)
      WRITE(*,'(A4,A16,ES16.6E2)') &
        '', 'SUM = ', SUM( rhsCF )
      WRITE(*,*)

      STOP

    END IF

  END SUBROUTINE ComputeRHS_Euler_DG


  SUBROUTINE ComputeRHS_Euler_DG_X1( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iNodeX1, jNodeX1, iNodeX2, iNodeX3
    INTEGER :: iNodeX, jNodeX, iCF
    REAL(DP) :: X1C, dX1
    REAL(DP) :: Alpha, AlphaPls, AlphaMns, AlphaMdl
    REAL(DP), DIMENSION(1:nCF) :: VolumeTerm, Flux_L, Flux_R, Flux
    REAL(DP), DIMENSION(1:nCF) :: uCF_L, uCF_R
    REAL(DP), DIMENSION(1:nCF) :: Lambda_L, Lambda_R
    REAL(DP), DIMENSION(1:nPF) :: uPF_L, uPF_R
    REAL(DP), DIMENSION(1:nAF) :: uAF_L, uAF_R
    REAL(DP), DIMENSION(nX(1)) :: a_X1_L, a_X1_R
    REAL(DP), DIMENSION(nX(1)) :: b_X1_L, b_X1_R
    REAL(DP), DIMENSION(1:nNodesX(1)) :: x_q, w_q
    REAL(DP), DIMENSION(1:nNodesX(1),1:nX(1)) :: a_X1_q
    REAL(DP), DIMENSION(1:nNodesX(1),1:nX(1)) :: b_X1_q
    REAL(DP), DIMENSION(1:nDOFX,1:nCF) :: uCF_P, uCF_K, uCF_N
    REAL(DP), DIMENSION(1:nDOFX,1:nCF) :: Flux_X1_q
    REAL(DP), DIMENSION(1:nDOFX,1:nPF) :: uPF_K
    REAL(DP), DIMENSION(1:nDOFX,1:nAF) :: uAF_K

    x_q = NodesX1
    w_q = WeightsX1

    ! -- Precompute Metric Functions --

    DO iX1 = iX_Begin(1), iX_End(1)
      X1C = MeshX(1) % Center(iX1)
      dX1 = MeshX(1) % Width (iX1)
      a_X1_L(iX1) &
        = One!a( [ X1C - 0.5_DP * dX1, 0.0_DP, 0.0_DP ] )
      a_X1_R(iX1) &
        = One!a( [ X1C + 0.5_DP * dX1, 0.0_DP, 0.0_DP ] )
      b_X1_L(iX1) &
        = One!b( [ X1C - 0.5_DP * dX1, 0.0_DP, 0.0_DP ] )
      b_X1_R(iX1) &
        = One!b( [ X1C + 0.5_DP * dX1, 0.0_DP, 0.0_DP ] )
      DO iNodeX1 = 1, nNodesX(1)
        a_X1_q(iNodeX1,iX1) &
          = One!a( [ X1C + dX1 * x_q(iNodeX1), 0.0_DP, 0.0_DP ] )
        b_X1_q(iNodeX1,iX1) &
          = One!b( [ X1C + dX1 * x_q(iNodeX1), 0.0_DP, 0.0_DP ] )
      END DO
    END DO

    ! -- Compute Right-Hand Side for Euler Equations --

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        !$OMP PARALLEL DO PRIVATE &
        !$OMP&              ( iX1, iNodeX1, iNodeX2, iNodeX3, iNodeX, &
        !$OMP&                jNodeX1, jNodeX, uCF_P, uCF_K, uCF_N,   &
        !$OMP&                uPF_K, uAF_K, Flux_X1_q, VolumeTerm,    &
        !$OMP&                dX1, iCF, uCF_L, uCF_R, uPF_L, uPF_R,   &
        !$OMP&                uAF_L, uAF_R, Lambda_L, Lambda_R,       &
        !$OMP&                AlphaPls, AlphaMns, Alpha, AlphaMdl,    &
        !$OMP&                Flux_L, Flux_R, Flux )
        DO iX1 = iX_Begin(1), iX_End(1)

          dX1 = MeshX(1) % Width (iX1)

          uCF_P = uCF(:,iX1-1,iX2,iX3,:) ! Previous Element
          uCF_K = uCF(:,iX1  ,iX2,iX3,:) ! This     Element
          uCF_N = uCF(:,iX1+1,iX2,iX3,:) ! Next     Element

          uPF_K = uPF(:,iX1,iX2,iX3,:) ! Primitive This Element
          uAF_K = uAF(:,iX1,iX2,iX3,:) ! Auxiliary This Element

          DO iNodeX = 1, nDOFX

            Flux_X1_q(iNodeX,1:nCF) &
              = Flux_X1 &
                  ( uPF_K(iNodeX,iPF_D ), uPF_K(iNodeX,iPF_V1), &
                    uPF_K(iNodeX,iPF_V2), uPF_K(iNodeX,iPF_V3), &
                    uPF_K(iNodeX,iPF_E ), uAF_K(iNodeX,iAF_P ), &
                    uPF_K(iNodeX,iPF_Ne) )

          END DO

          rhsCF(1:nDOFX,iX1,iX2,iX3,1:nCF) = Zero

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX &
                  = NodeNumberTableX &
                      ( iNodeX1, iNodeX2, iNodeX3 )

                ! -----------------
                ! -- Volume Term --
                ! -----------------

                VolumeTerm = Zero
                DO jNodeX1 = 1, nNodesX(1)

                  jNodeX &
                    = NodeNumberTableX &
                        ( jNodeX1, iNodeX2, iNodeX3 )

                  VolumeTerm(1:nCF) &
                    = VolumeTerm(1:nCF) &
                        + w_q(jNodeX1) &
                            * a_X1_q(jNodeX1,iX1) &
                            * b_X1_q(jNodeX1,iX1) &
                            * Flux_X1_q(jNodeX,1:nCF) &
                            * dL_X1_q(jNodeX1,iNodeX1)

                END DO

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                  = rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                      + VolumeTerm(1:nCF) &
                          / ( w_q(iNodeX1) * a_X1_q(iNodeX1,iX1) &
                                * b_X1_q(iNodeX1,iX1) * dX1 )

                ! -------------------
                ! -- Surface Terms --
                ! -------------------

                ! -- Left Face -- 

                ! -- Interpolate Left and Right State -- 

                DO iCF = 1, nCF

                  uCF_L(iCF) = Zero
                  uCF_R(iCF) = Zero
                  DO jNodeX1 = 1, nNodesX(1)

                    jNodeX = NodeNumberTableX( jNodeX1, iNodeX2, iNodeX3 )

                    uCF_L(iCF) &
                      = uCF_L(iCF) &
                          + L_X1_Up(jNodeX1) * uCF_P(jNodeX,iCF)

                    uCF_R(iCF) &
                      = uCF_R(iCF) &
                          + L_X1_Dn(jNodeX1) * uCF_K(jNodeX,iCF)

                  END DO

                END DO

                ! -- Numerical Flux --

                uPF_L = Primitive( uCF_L )

                CALL ComputeAuxiliary_Fluid &
                       ( uPF_L(iPF_D ), uPF_L(iPF_E), uPF_L(iPF_Ne), &
                         uAF_L(iAF_P ), uAF_L(iAF_T), uAF_L(iAF_Ye), &
                         uAF_L(iAF_S ), uAF_L(iAF_E), uAF_L(iAF_Gm), &
                         uAF_L(iAF_Cs) )

                Flux_L &
                  = Flux_X1 &
                      ( uPF_L(iPF_D), uPF_L(iPF_V1), uPF_L(iPF_V2), &
                        uPF_L(iPF_V3), uPF_L(iPF_E), uAF_L(iAF_P),  &
                        uPF_L(iPF_Ne) )

                Lambda_L &
                  = Eigenvalues &
                      ( uPF_L(iPF_V1), uAF_L(iAF_Cs) )

                uPF_R = Primitive( uCF_R )

                CALL ComputeAuxiliary_Fluid &
                       ( uPF_R(iPF_D ), uPF_R(iPF_E), uPF_R(iPF_Ne), &
                         uAF_R(iAF_P ), uAF_R(iAF_T), uAF_R(iAF_Ye), &
                         uAF_R(iAF_S ), uAF_R(iAF_E), uAF_R(iAF_Gm), &
                         uAF_R(iAF_Cs) )

                Flux_R &
                  = Flux_X1 &
                      ( uPF_R(iPF_D), uPF_R(iPF_V1), uPF_R(iPF_V2), &
                        uPF_R(iPF_V3), uPF_R(iPF_E), uAF_R(iAF_P),  &
                        uPF_R(iPF_Ne) )

                Lambda_R &
                  = Eigenvalues &
                      ( uPF_R(iPF_V1), uAF_R(iAF_Cs) )

                AlphaPls = AlphaP( Lambda_L, Lambda_R )
                AlphaMns = AlphaM( Lambda_L, Lambda_R )
                Alpha    = MAX( AlphaPls, AlphaMns )

                AlphaMdl &
                  = AlphaC &
                      ( [ uCF_L (iCF_D), uCF_L (iCF_S1) ], &
                        [ uCF_R (iCF_D), uCF_R (iCF_S1) ], &
                        [ Flux_L(iCF_D), Flux_L(iCF_S1) ], &
                        [ Flux_R(iCF_D), Flux_R(iCF_S1) ], &
                        AlphaPls, AlphaMns )

                Flux &
                  = NumericalFlux_Fluid &
                      ( uCF_L, uCF_R, Flux_L, Flux_R, &
                        Alpha, AlphaPls, AlphaMns, AlphaMdl, nCF )

                ! -- Contribution to Right-Hand Side --

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                  = rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                      + a_X1_L(iX1) * b_X1_L(iX1) &
                          * Flux(1:nCF) * L_X1_Dn(iNodeX1) &
                              / ( w_q(iNodeX1) * a_X1_q(iNodeX1,iX1) &
                                    * b_X1_q(iNodeX1,iX1) * dX1 )

                ! -- Right Face --

                ! -- Interpolate Left and Right State --

                DO iCF = 1, nCF
                  uCF_L(iCF) = Zero
                  uCF_R(iCF) = Zero
                  DO jNodeX1 = 1, nNodesX(1)

                    jNodeX &
                      = NodeNumberTableX &
                          ( jNodeX1, iNodeX2, iNodeX3 )

                    uCF_L(iCF) &
                      = uCF_L(iCF) &
                          + L_X1_Up(jNodeX1) * uCF_K(jNodeX,iCF)

                    uCF_R(iCF) &
                      = uCF_R(iCF) &
                          + L_X1_Dn(jNodeX1) * uCF_N(jNodeX,iCF)

                  END DO
                END DO

                uPF_L = Primitive( uCF_L )

                CALL ComputeAuxiliary_Fluid &
                       ( uPF_L(iPF_D ), uPF_L(iPF_E), uPF_L(iPF_Ne), &
                         uAF_L(iAF_P ), uAF_L(iAF_T), uAF_L(iAF_Ye), &
                         uAF_L(iAF_S ), uAF_L(iAF_E), uAF_L(iAF_Gm), &
                         uAF_L(iAF_Cs) )

                Flux_L &
                  = Flux_X1 &
                      ( uPF_L(iPF_D), uPF_L(iPF_V1), uPF_L(iPF_V2), &
                        uPF_L(iPF_V3), uPF_L(iPF_E), uAF_L(iAF_P),  &
                        uPF_L(iPF_Ne) )

                Lambda_L &
                  = Eigenvalues &
                      ( uPF_L(iPF_V1), uAF_L(iAF_Cs) )

                uPF_R = Primitive( uCF_R )

                CALL ComputeAuxiliary_Fluid &
                       ( uPF_R(iPF_D ), uPF_R(iPF_E), uPF_R(iPF_Ne), &
                         uAF_R(iAF_P ), uAF_R(iAF_T), uAF_R(iAF_Ye), &
                         uAF_R(iAF_S ), uAF_R(iAF_E), uAF_R(iAF_Gm), &
                         uAF_R(iAF_Cs) )

                Flux_R &
                  = Flux_X1 &
                      ( uPF_R(iPF_D), uPF_R(iPF_V1), uPF_R(iPF_V2), &
                        uPF_R(iPF_V3), uPF_R(iPF_E), uAF_R(iAF_P),  &
                        uPF_R(iPF_Ne) )

                Lambda_R &
                  = Eigenvalues &
                      ( uPF_R(iPF_V1), uAF_R(iAF_Cs) )

                ! -- Numerical Flux --

                AlphaPls = AlphaP( Lambda_L, Lambda_R )
                AlphaMns = AlphaM( Lambda_L, Lambda_R )
                Alpha    = MAX( AlphaPls, AlphaMns )

                AlphaMdl &
                  = AlphaC &
                      ( [ uCF_L (iCF_D), uCF_L (iCF_S1) ], &
                        [ uCF_R (iCF_D), uCF_R (iCF_S1) ], &
                        [ Flux_L(iCF_D), Flux_L(iCF_S1) ], &
                        [ Flux_R(iCF_D), Flux_R(iCF_S1) ], &
                        AlphaPls, AlphaMns )

                Flux &
                  = NumericalFlux_Fluid &
                      ( uCF_L, uCF_R, Flux_L, Flux_R, &
                        Alpha, AlphaPls, AlphaMns, AlphaMdl, nCF )

                ! -- Contribution to Right-Hand Side --

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                  = rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                      - a_X1_R(iX1) * b_X1_R(iX1) &
                          * Flux(1:nCF) * L_X1_Up(iNodeX1) &
                              / ( w_q(iNodeX1) * a_X1_q(iNodeX1,iX1) &
                                    * b_X1_q(iNodeX1,iX1) * dX1 )

              END DO
            END DO
          END DO

        END DO
        !$OMP END PARALLEL DO
      END DO
    END DO

  END SUBROUTINE ComputeRHS_Euler_DG_X1


  SUBROUTINE ComputeRHS_Euler_DG_GeometrySources( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    INTEGER  :: iCF, iNodeG
    REAL(DP) :: X1, X2, X3
    REAL(DP), DIMENSION(nDOFX,nCF) :: GeometrySources_Gravity
    REAL(DP), DIMENSION(nNodesX(1)*nX(1)) :: X1_G, Fg_G
    REAL(DP), DIMENSION(nDOFX,nX(1)) :: X_1

    IF( TRIM( CoordinateSystem ) == 'CARTESIAN' ) RETURN

    ! --- Fictitious Forces ---

    iNodeG = 1
    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        DO iX1 = iX_Begin(1), iX_End(1)

          ASSOCIATE &
            ( uCF_K => uCF(:,iX1,iX2,iX3,:), & ! This Element Conserved
              uAF_K => uAF(:,iX1,iX2,iX3,:) )  ! This Element Auxiliary

          DO iNodeX3 = 1, nNodesX(3)

            X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

            DO iNodeX2 = 1, nNodesX(2)

              X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                  = rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                      + GeometrySources &
                          ( uCF_K(iNodeX,iCF_D),  uCF_K(iNodeX,iCF_S1), &
                            uCF_K(iNodeX,iCF_S2), uCF_K(iNodeX,iCF_S3), &
                            uAF_K(iNodeX,iAF_P), [ X1, X2, X3 ] )

                X1_G(iNodeG) = X1; iNodeG = iNodeG+1
                X_1(iNodeX,iX1) = X1

              END DO
            END DO
          END DO

          END ASSOCIATE ! uCF_K, etc.

        END DO
      END DO 
    END DO

    IF( TRIM( CoordinateSystem ) == 'CYLINDRICAL' ) RETURN

    ! --- Gravitational Sources ---

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        DO iX1 = iX_Begin(1), iX_End(1)

          CALL ComputeGeometrySources_Gravity  &
                 ( [ MeshX(1) % Width(iX1),    &
                     MeshX(2) % Width(iX2),    &
                     MeshX(3) % Width(iX3) ],  &
                   uCF(:,iX1,  iX2,iX3,1:nCF), &
                   uGF(:,iX1,  iX2,iX3,1:nGF), &
                   uGF(:,iX1-1,iX2,iX3,1:nGF), &
                   uGF(:,iX1+1,iX2,iX3,1:nGF), &
                   GeometrySources_Gravity(1:nDOFX,1:nCF) )

          DO iCF = 1, nCF

            rhsCF(:,iX1,iX2,iX3,iCF) &
              = rhsCF(:,iX1,iX2,iX3,iCF) &
                  + GeometrySources_Gravity(:,iCF)

          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE ComputeRHS_Euler_DG_GeometrySources


  SUBROUTINE InitializeRHS

    INTEGER :: iNodeX, iNodeX1, jNodeX1, iNodeX2, iNodeX3
    REAL(DP), DIMENSION(nNodesX(1)) :: x_q

    ALLOCATE( NodeNumberTableX(nNodesX(1),nNodesX(2),nNodesX(3)) )

    iNodeX = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          NodeNumberTableX(iNodeX1,iNodeX2,iNodeX3) &
            = iNodeX

          iNodeX = iNodeX + 1

        END DO
      END DO
    END DO

    ALLOCATE( L_X1_Dn(nNodesX(1)) )
    ALLOCATE( L_X1_Up(nNodesX(1)) )
    ALLOCATE( dL_X1_q(nNodesX(1),nNodesX(1)) )

    x_q = MeshX(1) % Nodes
    DO jNodeX1 = 1, nNodesX(1)
      L_X1_Dn(jNodeX1) &
        = L_X1(jNodeX1) % P( - 0.5_DP )
      L_X1_Up(jNodeX1) &
        = L_X1(jNodeX1) % P( + 0.5_DP )
      DO iNodeX1 = 1, nNodesX(1)
        dL_X1_q(iNodeX1,jNodeX1) &
          = dL_X1(jNodeX1) % P( x_q(iNodeX1) )
      END DO
    END DO

  END SUBROUTINE InitializeRHS


  SUBROUTINE FinalizeRHS

    DEALLOCATE( NodeNumberTableX )

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


END MODULE EulerEquationsSolutionModule_DG
