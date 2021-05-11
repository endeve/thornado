MODULE AccretionShockUtilitiesModule

  USE KindModule,                  ONLY: &
    DP,    &
    SqrtTiny, &
    Zero,  &
    One,   &
    Two,   &
    Three, &
    Five,  &
    TwoPi, &
    FourPi
  USE ProgramHeaderModule,         ONLY: &
    nDOFX,  &
    nNodesX, &
    nNodes, &
    swX, &
    nX
  USE ReferenceElementModuleX,     ONLY: &
    NodeNumberTableX3D, &
    NodeNumberTableX, &
    nDOFX_X1, &
    WeightsX_X1, &
    NodesX2
  USE MeshModule,                  ONLY: &
    MeshType, &
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
    Kilometer,  &
    UnitsDisplay
  USE QuadratureModule,            ONLY: &
    GetQuadrature

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeAccretionShockDiagnostics


CONTAINS


  SUBROUTINE ComputeAccretionShockDiagnostics &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, &
      Power_Legendre, MeshX, AngleAveragedShockRadius )

    INTEGER ,       INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP),       INTENT(in)    :: uPF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP),       INTENT(in)    :: uAF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    TYPE(MeshType), INTENT(in)    :: MeshX(3)
    REAL(DP),       INTENT(inout) :: Power_Legendre(0:)
    REAL(DP),       INTENT(inout) :: AngleAveragedShockRadius

    CALL ComputePowerInLegendreModes &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, Power_Legendre )

    CALL ComputeAngleAveragedShockRadius &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, &
             MeshX, AngleAveragedShockRadius )

  END SUBROUTINE ComputeAccretionShockDiagnostics


  SUBROUTINE ComputePowerInLegendreModes &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, Power )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: uPF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: uAF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: Power(0:)

    INTEGER  :: iX1, iX2, iX3, iNX, iNX1, iNX2, X2nodes(nNodes)
    REAL(DP) :: X1(nNodes), X2(nNodes), dX1, dX2
    REAL(DP) :: xQ(nNodes), wQ(nNodes)
    REAL(DP) :: P0(nNodes), P1(nNodes), P2(nNodes)
    REAL(DP) :: G(0:2,nNodes,iX_B0(1):iX_E0(1))
    REAL(DP) :: Field(nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))
    REAL(DP) :: EntropyUnits

! Debugging
logical::debug=.false.
character(len=16)::fmt='(ES24.16E3)'
real(dp)::LegendreProduct(6)

    EntropyUnits = Erg / Centimeter**3 &
                     / ( Gram / Centimeter**3 )**( Gamma_IDEAL )

    CALL GetQuadrature( nNodes, xQ, wQ )

    ! --- Define Field used for Computing Power (e.g., entropy) ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNX = 1, nDOFX

        ! --- Entropy ---

        Field(iNX,iX1,iX2,iX3) &
          = LOG10( uAF(iNX,iX1,iX2,iX3,iAF_P) &
                     / uPF(iNX,iX1,iX2,iX3,iPF_D)**( Gamma_IDEAL ) &
                     / EntropyUnits )

      END DO

    END DO
    END DO
    END DO

!    Power = Zero
    G     = Zero

    ! --- Get nodes for integration over X2 ---

    IF     ( nNodes .EQ. 1 )THEN

      X2nodes = [ 1 ]

    ELSE IF( nNodes .EQ. 2 )THEN

      X2nodes = [ 1, 3 ]

    ELSE IF( nNodes .EQ. 3 )THEN

      X2nodes = [ 1, 4, 7  ]

    END IF

    ! --- Loop over radii ---

    DO iX1 = iX_B0(1), iX_E0(1)

      dX1 = MeshX(1) % Width(iX1) / Centimeter

      DO iNX1 = 1, nNodes

        X1(iNX1) = NodeCoordinate( MeshX(1), iX1, iNX1 ) / Centimeter

        ! --- For each radius, compute moments (G functions) ---

if(debug) LegendreProduct = Zero

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)

          dX2 = MeshX(2) % Width(iX2)

          DO iNX2 = 1, nNodes

            X2(iNX2) = NodeCoordinate( MeshX(2), iX2, iNX2 )

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

          G(0,iNX1,iX1) &
            = G(0,iNX1,iX1) &
                + dX2 * SUM( wQ * Field(X2nodes,iX1,iX2,iX3) &
                               * P0 * SIN( X2 ) )

          G(1,iNX1,iX1) &
            = G(1,iNX1,iX1) &
                + dX2 * SUM( wQ * Field(X2nodes,iX1,iX2,iX3) &
                               * P1 * SIN( X2 ) )

          G(2,iNX1,iX1) &
            = G(2,iNX1,iX1) &
                + dX2 * SUM( wQ * Field(X2nodes,iX1,iX2,iX3) &
                               * P2 * SIN( X2 ) )

        END DO ! iX2
        END DO ! iX3

if(debug)then
do inx=1,6
print fmt,LegendreProduct(inx)
enddo
stop
print*
endif

      END DO ! iNX1

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


  SUBROUTINE LocateShockRadius &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, MeshX, ShockRadius )

    INTEGER ,       INTENT(in)  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP),       INTENT(in)  :: uPF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP),       INTENT(in)  :: uAF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    TYPE(MeshType), INTENT(in)  :: MeshX(3)
    REAL(DP),       INTENT(out) :: ShockRadius(1:,iX_B0(2):,iX_B0(3):)

    INTEGER  :: iNX, iNX1, iNX2, iNX3, iX1, iX2, iX3, iNX_X1
    REAL(DP) :: X1

    REAL(DP) :: Entropy(nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B0(2):iX_E0(2), &
                              iX_B0(3):iX_E0(3))

    REAL(DP) :: EntropyUnits
    REAL(DP) :: EntropyThreshold

    ! --- Compute entropy in elements ---

    EntropyUnits = Erg / Centimeter**3 &
                     / ( Gram / Centimeter**3 )**( Gamma_IDEAL )

    EntropyThreshold = LOG10( 1.0e15_DP * EntropyUnits )

    Entropy = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      Entropy(iNX,iX1,iX2,iX3) &
        = LOG10( uAF(iNX,iX1,iX2,iX3,iAF_P) &
                   / uPF(iNX,iX1,iX2,iX3,iPF_D)**( Gamma_IDEAL ) )

    END DO
    END DO
    END DO
    END DO

    ! --- Find shock radius as function of theta and phi ---

    ShockRadius = Zero

    ! --- If the grid doesn't contain the shock, return ---

    IF( ALL( Entropy .GT. EntropyThreshold ) &
          .OR. ALL( Entropy .LT. EntropyThreshold ) ) RETURN

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B1(1), iX_E0(1) ! include ghost cell w/o double counting

      DO iNX3 = 1, nNodesX(3)
      DO iNX2 = 1, nNodesX(2)
      DO iNX1 = 1, nNodesX(1)

        iNX_X1 = iNX2 + iNX3 - 1

        IF( ShockRadius(iNX_X1,iX2,iX3) .GT. SqrtTiny ) CYCLE

        iNX = NodeNumberTableX3D(iNX1,iNX2,iNX3)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

        IF( Entropy(iNX,iX1,iX2,iX3) .LT. EntropyThreshold ) &
          ShockRadius(iNX_X1,iX2,iX3) = X1

      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE LocateShockRadius


  SUBROUTINE ComputeAngleAveragedShockRadius &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, &
      MeshX, AngleAveragedShockRadius )

    INTEGER ,       INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP),       INTENT(in)    :: uPF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP),       INTENT(in)    :: uAF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    TYPE(MeshType), INTENT(in)    :: MeshX(3)
    REAL(DP),       INTENT(inout) :: AngleAveragedShockRadius

    REAL(DP) :: ShockRadius(nDOFX_X1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    INTEGER  :: iX2, iX3

    CALL LocateShockRadius &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uPF, uAF, MeshX, ShockRadius )

    ! --- Compute angle averaged shock radius ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      AngleAveragedShockRadius &
        = AngleAveragedShockRadius &
            + One / FourPi * MeshX(2) % Width(iX2) * MeshX(3) % Width(iX3) &
            * SUM( WeightsX_X1 * ShockRadius(:,iX2,iX3) &
                   * SIN( MeshX(2) % Center(iX2) &
                            + NodesX2 * MeshX(2) % Width(iX2) ) ) / Kilometer

    END DO
    END DO

  END SUBROUTINE ComputeAngleAveragedShockRadius


END MODULE AccretionShockUtilitiesModule
