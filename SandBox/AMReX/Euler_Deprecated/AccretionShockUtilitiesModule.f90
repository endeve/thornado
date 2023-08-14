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
    nDimsX, &
    nDOFX,  &
    nNodesX, &
    swX, &
    nX
  USE ReferenceElementModuleX,     ONLY: &
    NodeNumberTableX3D, &
    NodeNumberTableX, &
    nDOFX_X1, &
    WeightsX1, &
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
    UnitsDisplay

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeAccretionShockDiagnostics
  PUBLIC :: ComputePowerInLegendreModes
  PUBLIC :: ComputeAngleAveragedShockRadius


CONTAINS


  SUBROUTINE ComputeAccretionShockDiagnostics &
    ( iX_B0, iX_E0, uPF, uAF, &
      MeshX, PowerIntegrand, ShockRadius )

    INTEGER ,       INTENT(in)    :: iX_B0(3), iX_E0(3)
    REAL(DP),       INTENT(in)    :: uPF(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP),       INTENT(in)    :: uAF(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    TYPE(MeshType), INTENT(in)    :: MeshX(3)
    REAL(DP),       INTENT(inout) :: PowerIntegrand(0:,1:,iX_B0(1):)
    REAL(DP),       INTENT(inout) :: ShockRadius(1:,iX_B0(2):,iX_B0(3):)

    CALL ComputePowerIntegrand &
           ( iX_B0, iX_E0, uPF, uAF, MeshX, PowerIntegrand )

    CALL LocateShockRadius &
           ( iX_B0, iX_E0, uPF, uAF, MeshX, ShockRadius )

  END SUBROUTINE ComputeAccretionShockDiagnostics


  SUBROUTINE ComputePowerIntegrand &
    ( iX_B0, iX_E0, uPF, uAF, MeshX, PowerIntegrand )

    INTEGER ,       INTENT(in)    :: iX_B0(3), iX_E0(3)
    REAL(DP),       INTENT(in)    :: uPF(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP),       INTENT(in)    :: uAF(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    TYPE(MeshType), INTENT(in)    :: MeshX(3)
    REAL(DP),       INTENT(inout) :: PowerIntegrand(0:,1:,iX_B0(1):)

    INTEGER  :: iX1, iX2, iX3, iNX, iNX1, iNX2
    REAL(DP) :: X2(nNodesX(2)), dX2
    REAL(DP) :: P0(nDOFX_X1), P1(nDOFX_X1), P2(nDOFX_X1)
    REAL(DP) :: Field(nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))
    REAL(DP) :: EntropyUnits
    INTEGER, ALLOCATABLE :: nNodesX_X1(:)

    IF( nDimsX .EQ. 1 ) RETURN

    EntropyUnits = Erg / Centimeter**3 &
                     / ( Gram / Centimeter**3 )**( Gamma_IDEAL )

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

    ! --- Get nodes for integration over X1-Faces ---

    ALLOCATE( nNodesX_X1(nDOFX_X1) )

    nNodesX_X1(1) = 1

    DO iNX1 = 2, nDOFX_X1

      nNodesX_X1(iNX1) = nNodesX_X1(iNX1-1) + nNodesX(1)

    END DO

    ! --- Loop over radii ---

    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNX1 = 1, nNodesX(1)

        ! --- For each radius, compute moments (G functions) ---

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)

          dX2 = MeshX(2) % Width(iX2)

          DO iNX2 = 1, nNodesX(2)

            X2(iNX2) = NodeCoordinate( MeshX(2), iX2, iNX2 )

          END DO

          CALL Legendre( X2, P0, P1, P2 )

          PowerIntegrand(0,iNX1,iX1) &
            = PowerIntegrand(0,iNX1,iX1) &
                + dX2 &
                    * SUM( WeightsX_X1 * Field(nNodesX_X1,iX1,iX2,iX3) &
                             * P0 * SIN( X2 ) )

          PowerIntegrand(1,iNX1,iX1) &
            = PowerIntegrand(1,iNX1,iX1) &
                + dX2 &
                    * SUM( WeightsX_X1 * Field(nNodesX_X1,iX1,iX2,iX3) &
                             * P1 * SIN( X2 ) )

          PowerIntegrand(2,iNX1,iX1) &
            = PowerIntegrand(2,iNX1,iX1) &
                + dX2 &
                    * SUM( WeightsX_X1 * Field(nNodesX_X1,iX1,iX2,iX3) &
                             * P2 * SIN( X2 ) )

        END DO ! iX2
        END DO ! iX3

      END DO ! iNX1

    END DO ! iX1

    DEALLOCATE( nNodesX_X1 )

  END SUBROUTINE ComputePowerIntegrand


  SUBROUTINE ComputePowerInLegendreModes &
    ( iX1_B, iX1_E, MeshX1, PowerIntegrand, Power )

    INTEGER,        INTENT(in)  :: iX1_B, iX1_E
    TYPE(MeshType), INTENT(in)  :: MeshX1
    REAL(DP),       INTENT(in)  :: PowerIntegrand(0:,1:,iX1_B:)
    REAL(DP),       INTENT(out) :: Power(0:)

    INTEGER  :: iX1, iNX1
    REAL(DP) :: X1(nNodesX(1)), dX1

    Power = Zero

    DO iX1 = iX1_B, iX1_E

      DO iNX1 = 1, nNodesX(1)

        X1(iNX1) = NodeCoordinate( MeshX1, iX1, iNX1 ) / Centimeter

      END DO

      dX1 = MeshX1 % Width(iX1) / Centimeter

      Power(0) &
        = Power(0) &
            + TwoPi * dX1 &
                * SUM( WeightsX1 * PowerIntegrand(0,:,iX1)**2 * X1**2 )

      Power(1) &
        = Power(1) &
            + TwoPi * dX1 &
                * SUM( WeightsX1 * PowerIntegrand(1,:,iX1)**2 * X1**2 )

      Power(2) &
        = Power(2) &
            + TwoPi * dX1 &
                * SUM( WeightsX1 * PowerIntegrand(2,:,iX1)**2 * X1**2 )

    END DO

  END SUBROUTINE ComputePowerInLegendreModes


  SUBROUTINE Legendre( X2, P0, P1, P2 )

    REAL(DP), INTENT(in)  :: X2(nNodesX(2))
    REAL(DP), INTENT(out) :: P0(nDOFX_X1), P1(nDOFX_X1), P2(nDOFX_X1)

    REAL(DP) :: X(nNodesX(2))

    X = COS( X2 )

    P0 = SQRT( One   / Two )
    P1 = SQRT( Three / Two ) * X
    P2 = SQRT( Five  / Two ) * ( Three * X**2 - One ) / Two

  END SUBROUTINE Legendre


  SUBROUTINE LocateShockRadius &
    ( iX_B0, iX_E0, uPF, uAF, MeshX, ShockRadius )

    INTEGER ,       INTENT(in)    :: iX_B0(3), iX_E0(3)
    REAL(DP),       INTENT(in)    :: uPF(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP),       INTENT(in)    :: uAF(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    TYPE(MeshType), INTENT(in)    :: MeshX(3)
    REAL(DP),       INTENT(inout) :: ShockRadius(1:,iX_B0(2):,iX_B0(3):)

    INTEGER  :: iNX, iNX1, iNX2, iNX3, iX1, iX2, iX3, iNX_X1
    REAL(DP) :: X1

    REAL(DP) :: Entropy(nDOFX,iX_B0(1):iX_E0(1), &
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
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      Entropy(iNX,iX1,iX2,iX3) &
        = LOG10( uAF(iNX,iX1,iX2,iX3,iAF_P) &
                   / uPF(iNX,iX1,iX2,iX3,iPF_D)**( Gamma_IDEAL ) )

    END DO
    END DO
    END DO
    END DO

    ! --- Find shock radius as function of theta and phi ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNX3 = 1, nNodesX(3)
      DO iNX2 = 1, nNodesX(2)
      DO iNX1 = 1, nNodesX(1)

        iNX_X1 = iNX2 + iNX3 - 1

        iNX = NodeNumberTableX3D(iNX1,iNX2,iNX3)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

        IF( Entropy(iNX,iX1,iX2,iX3) .LT. EntropyThreshold &
              .AND. X1 .LT. ShockRadius(iNX_X1,iX2,iX3) ) &
            ShockRadius(iNX_X1,iX2,iX3) = X1

      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE LocateShockRadius


  SUBROUTINE ComputeAngleAveragedShockRadius &
    ( iX_B0, iX_E0, ShockRadius, AngleAveragedShockRadius )

    INTEGER , INTENT(in)  :: iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: ShockRadius(1:,iX_B0(1):,iX_B0(2):)
    REAL(DP), INTENT(out) :: AngleAveragedShockRadius

    INTEGER  :: iX2, iX3

    AngleAveragedShockRadius = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      AngleAveragedShockRadius &
        = AngleAveragedShockRadius &
            + One / FourPi * MeshX(2) % Width(iX2) * MeshX(3) % Width(iX3) &
            * SUM( WeightsX_X1 * ShockRadius(:,iX2,iX3) &
                   * SIN( MeshX(2) % Center(iX2) &
                            + NodesX2 * MeshX(2) % Width(iX2) ) )

    END DO
    END DO

  END SUBROUTINE ComputeAngleAveragedShockRadius


END MODULE AccretionShockUtilitiesModule
