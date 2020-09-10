MODULE AccretionShockDiagnosticsModule

  USE KindModule,                  ONLY: &
    DP,    &
    Zero,  &
    One,   &
    Two,   &
    Three, &
    Five,  &
    TwoPi
  USE ProgramHeaderModule,         ONLY: &
    nDOFX,  &
    nNodes
  USE MeshModule,                  ONLY: &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule,           ONLY: &
    iPF_D, &
    iAF_P
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL
  USE UnitsModule,                 ONLY: &
    Erg,        &
    Gram,       &
    Centimeter, &
    UnitsDisplay
  USE QuadratureModule,            ONLY: &
    GetQuadrature

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeAccretionShockDiagnostics


CONTAINS


  SUBROUTINE ComputeAccretionShockDiagnostics &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, Power_Legendre )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: uPF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: uAF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: Power_Legendre(0:)

    CALL ComputePowerInLegendreModes &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, Power_Legendre )

  END SUBROUTINE ComputeAccretionShockDiagnostics


  SUBROUTINE ComputePowerInLegendreModes &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, Power )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: uPF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: uAF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: Power(0:)

    INTEGER  :: iX1, iX2, iX3, iNodeX, iNodeX1, iNodeX2, NodesX2(nNodes)
    REAL(DP) :: X1(nNodes), X2(nNodes), dX1, dX2
    REAL(DP) :: xQ(nNodes), wQ(nNodes)
    REAL(DP) :: P0(nNodes), P1(nNodes), P2(nNodes)
    REAL(DP) :: G(0:2,nNodes,iX_B0(1):iX_E0(1))
    REAL(DP) :: Field(nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))

! Debugging
logical::debug=.false.
character(len=16)::fmt='(ES24.16E3)'
real(dp)::LegendreProduct(6)

    CALL GetQuadrature( nNodes, xQ, wQ )

    ! --- Define Field used for Computing Power (e.g., entropy) ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        ! --- Entropy ---

        Field(iNodeX,iX1,iX2,iX3) &
          = LOG10( &
              ( uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Erg / Centimeter**3 ) ) &
                  / ( uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                  / ( Gram / Centimeter**3 ) )**( Gamma_IDEAL ) )

      END DO

    END DO
    END DO
    END DO

    Power = Zero
    G     = Zero

    ! --- Get nodes for integration over X2 ---

    IF     ( nNodes .EQ. 1 )THEN

      NodesX2 = [ 1 ]

    ELSE IF( nNodes .EQ. 2 )THEN

      NodesX2 = [ 1, 3 ]

    ELSE IF( nNodes .EQ. 3 )THEN

      NodesX2 = [ 1, 4, 7  ]

    END IF

    ! --- Loop over radii ---

    DO iX1 = iX_B0(1), iX_E0(1)

      dX1 = MeshX(1) % Width(iX1) / Centimeter

      DO iNodeX1 = 1, nNodes

        X1(iNodeX1) = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) / Centimeter

        ! --- For each radius, compute moments (G functions) ---

if(debug) LegendreProduct = Zero

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)

          dX2 = MeshX(2) % Width(iX2)

          DO iNodeX2 = 1, nNodes

            X2(iNodeX2) = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

          END DO

          CALL Legendre( X2, P0, P1, P2 )

if(debug)then
LegendreProduct(1) = LegendreProduct(1) + dX2 * SUM( wQ * P0 * P0 * SIN( X2 ) )
LegendreProduct(2) = LegendreProduct(2) + dX2 * SUM( wQ * P1 * P1 * SIN( X2 ) )
LegendreProduct(3) = LegendreProduct(3) + dX2 * SUM( wQ * P2 * P2 * SIN( X2 ) )
LegendreProduct(4) = LegendreProduct(4) + dX2 * SUM( wQ * P0 * P1 * SIN( X2 ) )
LegendreProduct(5) = LegendreProduct(5) + dX2 * SUM( wQ * P0 * P2 * SIN( X2 ) )
LegendreProduct(6) = LegendreProduct(6) + dX2 * SUM( wQ * P1 * P2 * SIN( X2 ) )
endif

          G(0,iNodeX1,iX1) &
            = G(0,iNodeX1,iX1) &
                + dX2 * SUM( wQ * Field(NodesX2,iX1,iX2,iX3) &
                               * P0 * SIN( X2 ) )

          G(1,iNodeX1,iX1) &
            = G(1,iNodeX1,iX1) &
                + dX2 * SUM( wQ * Field(NodesX2,iX1,iX2,iX3) &
                               * P1 * SIN( X2 ) )

          G(2,iNodeX1,iX1) &
            = G(2,iNodeX1,iX1) &
                + dX2 * SUM( wQ * Field(NodesX2,iX1,iX2,iX3) &
                               * P2 * SIN( X2 ) )

        END DO ! iX2
        END DO ! iX3

if(debug)then
do inodex=1,6
print fmt,LegendreProduct(inodex)
enddo
stop
print*
endif

      END DO ! iNodeX1

      Power(0) &
        = Power(0) + TwoPi * dX1 * SUM( wQ * G(0,:,iX1)**2 * X1**2 )

      Power(1) &
        = Power(1) + TwoPi * dX1 * SUM( wQ * G(1,:,iX1)**2 * X1**2 )

      Power(2) &
        = Power(2) + TwoPi * dX1 * SUM( wQ * G(2,:,iX1)**2 * X1**2 )

    END DO ! X1

if(debug)then
write(*,'(A,ES24.16E3)') 'Power(0) = ', Power(0)
write(*,'(A,ES24.16E3)') 'Power(1) = ', Power(1)
write(*,'(A,ES24.16E3)') 'Power(2) = ', Power(2)
endif

  END SUBROUTINE ComputePowerInLegendreModes


  SUBROUTINE Legendre( X2, P0, P1, P2 )

    REAL(DP), INTENT(in)  :: X2(nNodes)
    REAL(DP), INTENT(out) :: P0(nNodes), P1(nNodes), P2(nNodes)

    REAL(DP) :: X(nNodes)

    X = COS( X2 )

    P0 = SQRT( One   / Two )
    P1 = SQRT( Three / Two ) * X
    P2 = SQRT( Five  / Two ) * ( Three * X**2 - One ) / Two

  END SUBROUTINE Legendre


END MODULE AccretionShockDiagnosticsModule
