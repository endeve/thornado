MODULE XCFC_UtilitiesModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Two, &
    Pi
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOF, &
    nDimsX, &
    swE
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Psi, &
    iGF_K_dd_11, &
    iGF_K_dd_12, &
    iGF_K_dd_13, &
    iGF_K_dd_22, &
    iGF_K_dd_23, &
    iGF_K_dd_33
  USE GeometryBoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting, &
    ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MultiplyWithPsi6
  PUBLIC :: ApplyBoundaryConditions_Geometry_XCFC
  PUBLIC :: ComputeGravitationalMass
  PUBLIC :: UpdateConformalFactorAndMetric_XCFC
  PUBLIC :: UpdateLapseShiftCurvature_XCFC

  INTEGER, PUBLIC :: swX_GS(3)

  ! --- MF: Metric Fields ---

  INTEGER, PARAMETER, PUBLIC :: iMF_Psi     = 1
  INTEGER, PARAMETER, PUBLIC :: iMF_Alpha   = 2
  INTEGER, PARAMETER, PUBLIC :: iMF_Beta_1  = 3
  INTEGER, PARAMETER, PUBLIC :: iMF_Beta_2  = 4
  INTEGER, PARAMETER, PUBLIC :: iMF_Beta_3  = 5
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_11 = 6
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_12 = 7
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_13 = 8
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_22 = 9
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_23 = 10
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_33 = 11
  INTEGER, PARAMETER, PUBLIC :: nMF         = 11

  ! --- GS: Gravity/Geometry Sources ---

  INTEGER, PARAMETER, PUBLIC :: iGS_E  = 1
  INTEGER, PARAMETER, PUBLIC :: iGS_S1 = 2
  INTEGER, PARAMETER, PUBLIC :: iGS_S2 = 3
  INTEGER, PARAMETER, PUBLIC :: iGS_S3 = 4
  INTEGER, PARAMETER, PUBLIC :: iGS_S  = 5
  INTEGER, PARAMETER, PUBLIC :: iGS_Mg = 6
  INTEGER, PARAMETER, PUBLIC :: nGS    = 6

  INTERFACE MultiplyWithPsi6
    MODULE PROCEDURE MultiplyWithPsi6_X
    MODULE PROCEDURE MultiplyWithPsi6_Z
  END INTERFACE MultiplyWithPsi6

  ! https://amrex-codes.github.io/amrex/docs_html/Basics.html#fine-mask
  INTEGER, PARAMETER :: iLeaf    = 0
  INTEGER, PARAMETER :: iNotLeaf = 1

CONTAINS


  SUBROUTINE MultiplyWithPsi6_X &
    ( iX_B, iX_E, iX_B1, iX_E1, G, U, Power, Mask_Option )

    INTEGER , INTENT(in)    :: iX_B(3), iX_E(3), iX_B1(3), iX_E1(3), Power
    REAL(DP), INTENT(in)    :: G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    INTEGER , INTENT(in), OPTIONAL :: &
      Mask_Option(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP) :: Psi6
    INTEGER  :: iNX, iX1, iX2, iX3, iCF, nCF

    INTEGER :: Mask(iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1)

    IF( PRESENT( Mask_Option ) )THEN

      Mask = Mask_Option

    ELSE

      Mask = iLeaf

    END IF

    nCF = SIZE( U, DIM = 5 )

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

      DO iCF = 1, nCF
      DO iNX = 1, nDOFX

        Psi6 = G(iNX,iX1,iX2,iX3,iGF_Psi)**6

        U(iNX,iX1,iX2,iX3,iCF) &
          = U(iNX,iX1,iX2,iX3,iCF) * Psi6**( Power )

      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE MultiplyWithPsi6_X


  SUBROUTINE MultiplyWithPsi6_Z &
    ( iE_B1, iE_E1, iX_B, iX_E, iX_B1, iX_E1, G, U, Power, Mask_Option )

    INTEGER , INTENT(in)    :: &
      iE_B1, iE_E1, iX_B(3), iX_E(3), iX_B1(3), iX_E1(3), Power
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iE_B1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
    INTEGER , INTENT(in), OPTIONAL :: &
      Mask_Option(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP) :: Psi6
    INTEGER  :: iE_B0, iE_E0, nE
    INTEGER  :: iNZ, iE, iX1, iX2, iX3, iCR, iS, iNX, nCR, nS

    INTEGER :: Mask(iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1)

    IF( PRESENT( Mask_Option ) )THEN

      Mask = Mask_Option

    ELSE

      Mask = iLeaf

    END IF

    nCR = SIZE( U, DIM = 6 )
    nS  = SIZE( U, DIM = 7 )

    iE_B0 = iE_B1 + swE
    iE_E0 = iE_E1 - swE

    nE = iE_E0 - iE_B0 + 1

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

      DO iS  = 1    , nS
      DO iCR = 1    , nCR
      DO iE  = iE_B0, iE_E0
      DO iNZ = 1    , nDOF

        iNX = MOD( ( iNZ - 1 ) / nDOFE, nDOFX ) + 1

        Psi6 = G(iNX,iX1,iX2,iX3,iGF_Psi)**6

        U(iNZ,iE,iX1,iX2,iX3,iCR,iS) &
          = U(iNZ,iE,iX1,iX2,iX3,iCR,iS) * Psi6**( Power )

      END DO
      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE MultiplyWithPsi6_Z


  SUBROUTINE ApplyBoundaryConditions_Geometry_XCFC &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CALL ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G )

    CALL ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G )

  END SUBROUTINE ApplyBoundaryConditions_Geometry_XCFC


  SUBROUTINE ComputeGravitationalMass &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GS, GravitationalMass, Mask_Option )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      GS(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(inout) :: &
      GravitationalMass
    INTEGER , INTENT(in), OPTIONAL :: &
      Mask_Option(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3
    REAL(DP) :: d3X

    INTEGER :: Mask(iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1)

    IF( PRESENT( Mask_Option ) )THEN

      Mask = Mask_Option

    ELSE

      Mask = iLeaf

    END IF

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    ! --- Mg has been pre-multiplied with \alpha * \sqrt{ \gamma }

    ! --- Assuming 1D spherical symmetry ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

      IF( nDimsX .EQ. 1 )THEN

        d3X = Two / Pi * dX1(iX1) * dX2(iX2) * dX3(iX3)

      ELSE

        d3X = dX1(iX1) * dX2(iX2) * dX3(iX3)

      END IF

      GravitationalMass &
        = GravitationalMass + d3X                  &
            * SUM( WeightsX_q * GS(:,iX1,iX2,iX3,iGS_Mg) )

    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeGravitationalMass


  SUBROUTINE UpdateConformalFactorAndMetric_XCFC &
    ( iX_B0, iX_E0, iX_B1, iX_E1, M, G, Mask_Option )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: M(1:,iX_B0(1)-swX_GS(1):, &
                                    iX_B0(2)-swX_GS(2):, &
                                    iX_B0(3)-swX_GS(3):,1:)
    REAL(DP), INTENT(inout) :: G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    INTEGER , INTENT(in), OPTIONAL :: &
      Mask_Option(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iNX, iNX1, iNX2, iX_B(3), iX_E(3)
    REAL(DP) :: X1, X2, Psi

    INTEGER :: Mask(iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1)

    IF( PRESENT( Mask_Option ) )THEN

      Mask = Mask_Option

    ELSE

      Mask = iLeaf

    END IF

    iX_B = iX_B0 - swX_GS
    iX_E = iX_E0 + swX_GS

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

      DO iNX = 1, nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

        Psi = M(iNX,iX1,iX2,iX3,iMF_Psi)

        G(iNX,iX1,iX2,iX3,iGF_Psi) = Psi
        G(iNX,iX1,iX2,iX3,iGF_h_1) = Psi**2
        G(iNX,iX1,iX2,iX3,iGF_h_2) = Psi**2 * X1
        G(iNX,iX1,iX2,iX3,iGF_h_3) = Psi**2 * X1 * SIN( X2 )

      END DO ! iNX

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

    END DO ! iX1
    END DO ! iX2
    END DO ! iX3

  END SUBROUTINE UpdateConformalFactorAndMetric_XCFC


  SUBROUTINE UpdateLapseShiftCurvature_XCFC &
    ( iX_B0, iX_E0, iX_B1, iX_E1, M, G, Mask_Option )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: M(1:,iX_B0(1)-swX_GS(1):, &
                                    iX_B0(2)-swX_GS(2):, &
                                    iX_B0(3)-swX_GS(3):,1:)
    REAL(DP), INTENT(inout) :: G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    INTEGER , INTENT(in), OPTIONAL :: &
      Mask_Option(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX2, iX3, iX_B(3), iX_E(3)

    INTEGER :: Mask(iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1)

    IF( PRESENT( Mask_Option ) )THEN

      Mask = Mask_Option

    ELSE

      Mask = iLeaf

    END IF

    iX_B = iX_B0 - swX_GS
    iX_E = iX_E0 + swX_GS

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

      DO iNX = 1, nDOFX

        G(iNX,iX1,iX2,iX3,iGF_Alpha) = M(iNX,iX1,iX2,iX3,iMF_Alpha)

        G(iNX,iX1,iX2,iX3,iGF_Beta_1) = M(iNX,iX1,iX2,iX3,iMF_Beta_1)
        G(iNX,iX1,iX2,iX3,iGF_Beta_2) = M(iNX,iX1,iX2,iX3,iMF_Beta_2)
        G(iNX,iX1,iX2,iX3,iGF_Beta_3) = M(iNX,iX1,iX2,iX3,iMF_Beta_3)

        G(iNX,iX1,iX2,iX3,iGF_K_dd_11) = M(iNX,iX1,iX2,iX3,iMF_K_dd_11)
        G(iNX,iX1,iX2,iX3,iGF_K_dd_12) = M(iNX,iX1,iX2,iX3,iMF_K_dd_12)
        G(iNX,iX1,iX2,iX3,iGF_K_dd_13) = M(iNX,iX1,iX2,iX3,iMF_K_dd_13)
        G(iNX,iX1,iX2,iX3,iGF_K_dd_22) = M(iNX,iX1,iX2,iX3,iMF_K_dd_22)
        G(iNX,iX1,iX2,iX3,iGF_K_dd_23) = M(iNX,iX1,iX2,iX3,iMF_K_dd_23)
        G(iNX,iX1,iX2,iX3,iGF_K_dd_33) = M(iNX,iX1,iX2,iX3,iMF_K_dd_33)

      END DO ! iNX

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

    END DO ! iX1
    END DO ! iX2
    END DO ! iX3

  END SUBROUTINE UpdateLapseShiftCurvature_XCFC


  ! --- PRIVATE SUBROUTINES ---


  LOGICAL FUNCTION IsNotLeafElement( Element )

    INTEGER, INTENT(in) :: Element

    IF( Element .EQ. iNotLeaf )THEN
      IsNotLeafElement = .TRUE.
    ELSE
      IsNotLeafElement = .FALSE.
    END IF

    RETURN
  END FUNCTION IsNotLeafElement


END MODULE XCFC_UtilitiesModule
