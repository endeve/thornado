MODULE AccretionShockUtilitiesModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    Three, &
    Five, &
    TwoPi
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodes, &
    nNodesX, &
    nDimsX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    iPF_D, &
    iAF_P
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE UnitsModule, ONLY: &
    Erg, &
    Gram, &
    Centimeter
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WriteInitialConditionsToFile
  PUBLIC :: ComputeAccretionShockDiagnostics


CONTAINS


  SUBROUTINE WriteInitialConditionsToFile &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, FileName )

    INTEGER ,           INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP),           INTENT(in)    :: &
      uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP),           INTENT(inout) :: &
      uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    CHARACTER(LEN=128),   INTENT(in)  :: &
      FileName

    INTEGER  :: iNX, iX1
    REAL(DP) :: D(nNodesX(1)), V1(nNodesX(1)), V2(nNodesX(1)), V3(nNodesX(1)), &
                E(nNodesX(1)), Ne(nNodesX(1)), P (nNodesX(1))
    CHARACTER(LEN=16) :: FMT

    IF( nDimsX .GT. 1 )THEN

      WRITE(*,*) 'WriteInitialConditionsToFile not implemented for nDimsX > 1.'
      RETURN

    END IF

    WRITE(FMT,'(A3,I3.3,A10)') '(SP', nDOFX, 'ES25.16E3)'

    OPEN( 101, FILE = TRIM( FileName ) // '_D.dat' )
    OPEN( 102, FILE = TRIM( FileName ) // '_V.dat' )
    OPEN( 103, FILE = TRIM( FileName ) // '_P.dat' )

    WRITE(101,'(A16)') TRIM( FMT )
    WRITE(102,'(A16)') TRIM( FMT )
    WRITE(103,'(A16)') TRIM( FMT )

    CALL ApplyBoundaryConditions_Euler( iX_B0, iX_E0, iX_B1, iX_E1, uCF )

    DO iX1 = iX_B1(1), iX_E1(1)

      CALL ComputePrimitive_Euler_Relativistic &
             ( uCF(:,iX1,1,1,iCF_D ), &
               uCF(:,iX1,1,1,iCF_S1), &
               uCF(:,iX1,1,1,iCF_S2), &
               uCF(:,iX1,1,1,iCF_S3), &
               uCF(:,iX1,1,1,iCF_E ), &
               uCF(:,iX1,1,1,iCF_Ne), &
               D, V1, V2, V3, E, Ne, &
               uGF(:,iX1,1,1,iGF_Gm_dd_11), &
               uGF(:,iX1,1,1,iGF_Gm_dd_22), &
               uGF(:,iX1,1,1,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( D, E, Ne, P )

      WRITE(101,FMT) D
      WRITE(102,FMT) V1
      WRITE(103,FMT) P

    END DO

    CLOSE( 103 )
    CLOSE( 102 )
    CLOSE( 101 )

  END SUBROUTINE WriteInitialConditionsToFile


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

    ! --- Compute power in Legendre modes a la
    !     Blondin & Mezzacappa (2006), ApJ, 642, 401 ---

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: uPF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: uAF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: Power(0:)

    INTEGER  :: iX1, iX2, iX3, iNX, iNX1, iNX2, nX(3), nDOF(3), &
                NodesX2(nNodes), iX, iNX2_L, iNX2_U, iNX1_L, iNX1_U
    REAL(DP) :: dX1, dX2
    REAL(DP) :: xQ(nNodes), wQ(nNodes)
    REAL(DP) :: G(0:2,nNodes,iX_B0(1):iX_E0(1))
    REAL(DP) :: Field(nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))
    REAL(DP) :: FieldUnits
    REAL(DP), ALLOCATABLE :: X1(:), X2(:), P0(:), P1(:), P2(:)

    nX   = iX_E0 - iX_B0 + 1
    nDOF = nNodesX * nX

    FieldUnits = ( Erg / Centimeter**3 ) &
                   / ( Gram / Centimeter**3)**( Gamma_IDEAL )

    CALL GetQuadrature( nNodes, xQ, wQ )

    ALLOCATE( X1(iX_B0(1):iX_B0(1)+nDOF(1)-1) )
    ALLOCATE( X2(iX_B0(2):iX_B0(2)+nDOF(2)-1) )
    ALLOCATE( P0(iX_B0(2):iX_B0(2)+nDOF(2)-1) )
    ALLOCATE( P1(iX_B0(2):iX_B0(2)+nDOF(2)-1) )
    ALLOCATE( P2(iX_B0(2):iX_B0(2)+nDOF(2)-1) )

    ! --- Populate X1 array ---

    iX = 0

    DO iX1 = iX_B1(2), iX_E1(2)

      DO iNX1 = 1, nNodesX(1)

        iX = iX + 1

        X1(iX) = NodeCoordinate( MeshX(1), iX1, iNX1 ) / Centimeter

      END DO

    END DO

    ! --- Populate X2 array ---

    iX = 0

    DO iX2 = iX_B0(2), iX_E0(2)

      DO iNX2 = 1, nNodesX(2)

        iX = iX + 1

        X2(iX) = NodeCoordinate( MeshX(2), iX2, iNX2 )

      END DO

    END DO

    CALL ComputeLegendrePolynomials( X2, P0, P1, P2 )

    ! --- Define field used for computing power (e.g., entropy) ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      ! --- Entropy ---

      Field(iNX,iX1,iX2,iX3) &
        = LOG10( ( uAF(iNX,iX1,iX2,iX3,iAF_P) &
                     / uPF(iNX,iX1,iX2,iX3,iPF_D)**( Gamma_IDEAL ) ) &
                     / FieldUnits )

    END DO
    END DO
    END DO
    END DO

    Power = Zero
    G     = Zero

    ! --- Loop over radii ---

    DO iX1 = iX_B0(1), iX_E0(1)

      dX1 = MeshX(1) % Width(iX1) / Centimeter

      DO iNX1 = 1, nNodes

        CALL ComputeIntegrationNodesX2( iNX1, NodesX2 )

        ! --- For each radius, compute moments (G functions, Eq. (7)) ---

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)

          dX2 = MeshX(2) % Width(iX2)

          iNX2_L = iX_B0(2) + nNodesX(2) * ( iX2 - iX_B0(2) )
          iNX2_U = iNX2_L + nNodesX(2) - 1

          G(0,iNX1,iX1) &
            = G(0,iNX1,iX1) &
                + dX2 * SUM( wQ * Field(NodesX2,iX1,iX2,iX3) &
                               * P0(iNX2_L:iNX2_U) * SIN( X2(iNX2_L:iNX2_U) ) )

          G(1,iNX1,iX1) &
            = G(1,iNX1,iX1) &
                + dX2 * SUM( wQ * Field(NodesX2,iX1,iX2,iX3) &
                               * P1(iNX2_L:iNX2_U) * SIN( X2(iNX2_L:iNX2_U) ) )

          G(2,iNX1,iX1) &
            = G(2,iNX1,iX1) &
                + dX2 * SUM( wQ * Field(NodesX2,iX1,iX2,iX3) &
                               * P2(iNX2_L:iNX2_U) * SIN( X2(iNX2_L:iNX2_U) ) )

        END DO ! iX2
        END DO ! iX3

      END DO ! iNX1

      ! --- Compute powers (Eq. (8)) ---

      iNX1_L = iX_B0(1) + nNodesX(1) * ( iX1 - iX_B0(1) )
      iNX1_U = iNX1_L + nNodesX(1) - 1

      Power(0) &
        = Power(0) &
            + TwoPi * dX1 * SUM( wQ * G(0,:,iX1)**2 * X1(iNX1_L:iNX1_U)**2 )

      Power(1) &
        = Power(1) &
            + TwoPi * dX1 * SUM( wQ * G(1,:,iX1)**2 * X1(iNX1_L:iNX1_U)**2 )

      Power(2) &
        = Power(2) &
            + TwoPi * dX1 * SUM( wQ * G(2,:,iX1)**2 * X1(iNX1_L:iNX1_U)**2 )

    END DO ! X1

    DEALLOCATE( P2 )
    DEALLOCATE( P1 )
    DEALLOCATE( P0 )
    DEALLOCATE( X2 )
    DEALLOCATE( X1 )

  END SUBROUTINE ComputePowerInLegendreModes


  SUBROUTINE ComputeLegendrePolynomials( X2, P0, P1, P2 )

    REAL(DP), INTENT(in)  :: X2(:)
    REAL(DP), INTENT(out) :: P0(:), P1(:), P2(:)

    REAL(DP) :: X(SIZE(X2))

    X = COS( X2 )

    P0 = SQRT( One   / Two )
    P1 = SQRT( Three / Two ) * X
    P2 = SQRT( Five  / Two ) * Half * ( Three * X**2 - One )

  END SUBROUTINE ComputeLegendrePolynomials


  SUBROUTINE ComputeIntegrationNodesX2( iNX1, NodesX2 )

    INTEGER, INTENT(in)  :: iNX1
    INTEGER, INTENT(out) :: NodesX2(:)

    IF     ( nNodes .EQ. 1 )THEN

      NodesX2 = [ 1 ]

    ELSE IF( nNodes .EQ. 2 )THEN

      IF( iNX1 .EQ. 1 )THEN

        NodesX2 = [ 1, 3 ]

      ELSE

        NodesX2 = [ 2, 4 ]

      END IF

    ELSE IF( nNodes .EQ. 3 )THEN

      IF( iNX1 .EQ. 1 )THEN

        NodesX2 = [ 1, 4, 7  ]

      ELSE IF( iNX1 .EQ. 2 )THEN

        NodesX2 = [ 2, 5, 9 ]

      ELSE

        NodesX2 = [ 3, 6, 9 ]

      END IF

    ELSE

      PRINT*, 'Not implemented for nNodes > 3. Stopping...'
      STOP

    END IF

  END SUBROUTINE ComputeIntegrationNodesX2


END MODULE AccretionShockUtilitiesModule
