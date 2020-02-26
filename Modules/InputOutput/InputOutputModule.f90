MODULE InputOutputModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nE, nNodesE, &
    nDOF
  USE UtilitiesModule, ONLY: &
    NodeNumber
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE PolynomialBasisModule_Lagrange, ONLY: &
    evalL, &
    evalLX
  USE MeshModule, ONLY: &
    MeshType, &
    MeshX, MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Phi_N, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, iAF_Me, iAF_Mp, &
    iAF_Mn, iAF_Xp, iAF_Xn, iAF_Xa, iAF_Xh, iAF_Gm, iAF_Cs, &
    Shock
  USE RadiationFieldsModule, ONLY: &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    Discontinuity

  IMPLICIT NONE
  PRIVATE

  CHARACTER(9),  PARAMETER :: &
    OutputDirectory = '../Output'
  CHARACTER(14), PARAMETER :: &
    GeometrySuffix = 'GeometryFields'
  CHARACTER(11), PARAMETER :: &
    FluidSuffix = 'FluidFields'
  CHARACTER(15), PARAMETER :: &
    RadiationSuffix = 'RadiationFields'
  INTEGER :: FileNumber = 0
  INTEGER :: RestartFileNumber = 0

  PUBLIC :: WriteFields1D
  PUBLIC :: WriteFieldsRestart1D
  PUBLIC :: ReadFluidFieldsRestart1D
  PUBLIC :: WriteRadiationFields

CONTAINS


  SUBROUTINE WriteFields1D &
               ( Time, WriteGeometryFields_Option, WriteFluidFields_Option, &
                 WriteRadiationFields_Option )

    REAL(DP), INTENT(in)           :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: WriteGeometryFields_Option
    LOGICAL,  INTENT(in), OPTIONAL :: WriteFluidFields_Option
    LOGICAL,  INTENT(in), OPTIONAL :: WriteRadiationFields_Option

    LOGICAL :: WriteGeometryFields
    LOGICAL :: WriteFluidFields
    LOGICAL :: WriteRadiationFields

    WriteGeometryFields = .FALSE.
    IF( PRESENT( WriteGeometryFields_Option ) ) &
      WriteGeometryFields = WriteGeometryFields_Option

    WriteFluidFields = .FALSE.
    IF( PRESENT( WriteFluidFields_Option ) ) &
      WriteFluidFields = WriteFluidFields_Option

    WriteRadiationFields = .FALSE.
    IF( PRESENT( WriteRadiationFields_Option ) ) &
      WriteRadiationFields = WriteRadiationFields_Option

    IF( WriteGeometryFields ) &
      CALL WriteGeometryFields1D( Time )

    IF( WriteFluidFields ) &
      CALL WriteFluidFields1D_N( Time )

    IF( WriteRadiationFields ) &
      CALL WriteRadiationFields1D_N( Time )

    FileNumber = FileNumber + 1

  END SUBROUTINE WriteFields1D


  SUBROUTINE WriteGeometryFields1D( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    INTEGER        :: FUNIT

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName = OutputDirectory // '/' // TRIM( ProgramName ) // '_' // &
                 GeometrySuffix // '_' // FileNumberString // '.dat'

    ASSOCIATE( X1 => MeshX(1) % Center(1:nX(1)), &
               U  => UnitsDisplay )

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) &
      Time / U % TimeUnit, nX(1), X1 / U % LengthX1Unit, &
      GeometryField1D( uGF(:,1:nX(1),1,1,iGF_Phi_N    ), nX(1) ) &
        / ( U % EnergyDensityUnit / U % MassDensityUnit ), &
      GeometryField1D( uGF(:,1:nX(1),1,1,iGF_Gm_dd_11 ), nX(1) ), &
      GeometryField1D( uGF(:,1:nX(1),1,1,iGF_Gm_dd_22 ), nX(1) ), &
      GeometryField1D( uGF(:,1:nX(1),1,1,iGF_Gm_dd_33 ), nX(1) )

    CLOSE( FUNIT )

    END ASSOCIATE ! X1, etc.

  END SUBROUTINE WriteGeometryFields1D


  FUNCTION GeometryField1D( u, nX1 )

    REAL(DP), DIMENSION(nX1)             :: GeometryField1D
    REAL(DP), DIMENSION(:,:), INTENT(in) :: u
    INTEGER,                  INTENT(in) :: nX1

    INTEGER :: iX1

    DO iX1 = 1, nX1
      GeometryField1D(iX1) &
        = evalLX( u(:,iX1), 0.0_DP, 0.0_DP, 0.0_DP )
    END DO

    RETURN
  END FUNCTION GeometryField1D


  SUBROUTINE WriteFluidFields1D( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    INTEGER        :: FUNIT

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName = OutputDirectory // '/' // TRIM( ProgramName ) // '_' // &
                 FluidSuffix // '_' // FileNumberString // '.dat'

    ASSOCIATE( X1 => MeshX(1) % Center(1:nX(1)), &
               U  => UnitsDisplay )

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) &
      Time / U % TimeUnit, nX(1), X1 / U % LengthX1Unit, &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_D ), nX(1) ) &
        / U % MassDensityUnit, &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_S1), nX(1) ) &
        / U % MomentumDensityUnit, &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_S2), nX(1) ) &
        / U % MomentumDensityUnit, &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_S3), nX(1) ) &
        / U % MomentumDensityUnit, &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_E ), nX(1) ) &
        / U % EnergyDensityUnit, &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_Ne), nX(1) ) &
        / U % ParticleDensityUnit, &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_D ), nX(1) ) &
        / U % MassDensityUnit, &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_V1), nX(1) ) &
        / U % VelocityX1Unit, &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_V2), nX(1) ) &
        / U % VelocityX2Unit, &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_V3), nX(1) ) &
        / U % VelocityX3Unit, &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_E ), nX(1) ) &
        / U % EnergyDensityUnit, &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_Ne), nX(1) ) &
        / U % ParticleDensityUnit, &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_P ), nX(1) ) &
        / U % PressureUnit, &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_T ), nX(1) ) &
        / U % TemperatureUnit, &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Ye), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_S ), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_E ), nX(1) ) &
        / ( U % EnergyDensityUnit / U % MassDensityUnit ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Me), nX(1) ) &
        / U % EnergyUnit, &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Mp), nX(1) ) &
        / U % EnergyUnit, &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Mn), nX(1) ) &
        / U % EnergyUnit, &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Xp), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Xn), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Xa), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Xh), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Gm), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Cs), nX(1) ) &
        / U % VelocityX1Unit, &
      Shock(1:nX(1),1,1)

    CLOSE( FUNIT )

    END ASSOCIATE ! X1, etc.

  END SUBROUTINE WriteFluidFields1D


  FUNCTION FluidField1D( u, nX1 )

    REAL(DP), DIMENSION(nX1)             :: FluidField1D
    REAL(DP), DIMENSION(:,:), INTENT(in) :: u
    INTEGER,                  INTENT(in) :: nX1

    INTEGER :: iX1

    DO iX1 = 1, nX1
      FluidField1D(iX1) &
        = evalLX( u(:,iX1), 0.0_DP, 0.0_DP, 0.0_DP )
    END DO

    RETURN
  END FUNCTION FluidField1D


  SUBROUTINE WriteFluidFields1D_N( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    INTEGER        :: FUNIT

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName = OutputDirectory // '/' // TRIM( ProgramName ) // '_' // &
                 FluidSuffix // '_' // FileNumberString // '.dat'

    ASSOCIATE( U  => UnitsDisplay )

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) &
      Time / U % TimeUnit, &
      nX(1), nNodesX(1), &
      NodeCoordinatesX1( nX(1), nNodesX(1) ) / U % LengthX1Unit, &
      FluidField1D_N( uCF(:,1:nX(1),1,1,iCF_D ), nX(1), nNodesX(1) ) &
        / U % MassDensityUnit, &
      FluidField1D_N( uCF(:,1:nX(1),1,1,iCF_S1), nX(1), nNodesX(1) ) &
        / U % MomentumDensityUnit, &
      FluidField1D_N( uCF(:,1:nX(1),1,1,iCF_S2), nX(1), nNodesX(1) ) &
        / U % MomentumDensityUnit, &
      FluidField1D_N( uCF(:,1:nX(1),1,1,iCF_S3), nX(1), nNodesX(1) ) &
        / U % MomentumDensityUnit, &
      FluidField1D_N( uCF(:,1:nX(1),1,1,iCF_E ), nX(1), nNodesX(1) ) &
        / U % EnergyDensityUnit, &
      FluidField1D_N( uCF(:,1:nX(1),1,1,iCF_Ne), nX(1), nNodesX(1) ) &
        / U % ParticleDensityUnit, &
      FluidField1D_N( uPF(:,1:nX(1),1,1,iPF_D ), nX(1), nNodesX(1) ) &
        / U % MassDensityUnit, &
      FluidField1D_N( uPF(:,1:nX(1),1,1,iPF_V1), nX(1), nNodesX(1) ) &
        / U % VelocityX1Unit, &
      FluidField1D_N( uPF(:,1:nX(1),1,1,iPF_V2), nX(1), nNodesX(1) ) &
        / U % VelocityX2Unit, &
      FluidField1D_N( uPF(:,1:nX(1),1,1,iPF_V3), nX(1), nNodesX(1) ) &
        / U % VelocityX3Unit, &
      FluidField1D_N( uPF(:,1:nX(1),1,1,iPF_E ), nX(1), nNodesX(1) ) &
        / U % EnergyDensityUnit, &
      FluidField1D_N( uPF(:,1:nX(1),1,1,iPF_Ne), nX(1), nNodesX(1) ) &
        / U % ParticleDensityUnit, &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_P ), nX(1), nNodesX(1) ) &
        / U % PressureUnit, &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_T ), nX(1), nNodesX(1) ) &
        / U % TemperatureUnit, &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_Ye), nX(1), nNodesX(1) ), &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_S ), nX(1), nNodesX(1) ), &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_E ), nX(1), nNodesX(1) ) &
        / ( U % EnergyDensityUnit / U % MassDensityUnit ), &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_Me), nX(1), nNodesX(1) ) &
        / U % EnergyUnit, &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_Mp), nX(1), nNodesX(1) ) &
        / U % EnergyUnit, &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_Mn), nX(1), nNodesX(1) ) &
        / U % EnergyUnit, &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_Xp), nX(1), nNodesX(1) ), &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_Xn), nX(1), nNodesX(1) ), &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_Xa), nX(1), nNodesX(1) ), &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_Xh), nX(1), nNodesX(1) ), &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_Gm), nX(1), nNodesX(1) ), &
      FluidField1D_N( uAF(:,1:nX(1),1,1,iAF_Cs), nX(1), nNodesX(1) ) &
        / U % VelocityX1Unit, &
      Shock(1:nX(1),1,1)

    CLOSE( FUNIT )

    END ASSOCIATE ! U

  END SUBROUTINE WriteFluidFields1D_N


  FUNCTION FluidField1D_N( u, nX1, nNX1 )

    REAL(DP), DIMENSION(nX1*nNX1)        :: FluidField1D_N
    REAL(DP), DIMENSION(:,:), INTENT(in) :: u
    INTEGER,                  INTENT(in) :: nX1, nNX1

    INTEGER :: i, iX1, iNX1

    i = 1
    DO iX1 = 1, nX1
      DO iNX1 = 1, nNX1
        FluidField1D_N(i) = u(iNX1,iX1)
        i = i + 1
      END DO
    END DO

    RETURN
  END FUNCTION FluidField1D_N


  SUBROUTINE WriteRadiationFields1D( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    INTEGER        :: FUNIT

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName = OutputDirectory // '/' // TRIM( ProgramName ) // '_' // &
                 RadiationSuffix // '_' // FileNumberString // '.dat'

    ASSOCIATE &
      ( X1 => MeshX(1) % Center(1:nX(1)), &
        E  => MeshE    % Center(1:nE), &
        U  => UnitsDisplay )

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) &
      Time / U % TimeUnit, &
      nE, E / U % EnergyUnit, nX(1), X1 / U % LengthX1Unit, &
      RadiationField1D( uCR(:,1:nE,1:nX(1),1,1,iCR_N, 1), nE, nX(1) ), &
      RadiationField1D( uCR(:,1:nE,1:nX(1),1,1,iCR_G1,1), nE, nX(1) ), &
      RadiationField1D( uCR(:,1:nE,1:nX(1),1,1,iCR_G2,1), nE, nX(1) ), &
      RadiationField1D( uCR(:,1:nE,1:nX(1),1,1,iCR_G3,1), nE, nX(1) ), &
      RadiationField1D( uPR(:,1:nE,1:nX(1),1,1,iPR_D, 1), nE, nX(1) ), &
      RadiationField1D( uCR(:,1:nE,1:nX(1),1,1,iPR_I1,1), nE, nX(1) ), &
      RadiationField1D( uCR(:,1:nE,1:nX(1),1,1,iPR_I2,1), nE, nX(1) ), &
      RadiationField1D( uCR(:,1:nE,1:nX(1),1,1,iPR_I3,1), nE, nX(1) ), &
      Discontinuity(1:nE,1:nX(1),1,1)

    CLOSE( FUNIT )

    END ASSOCIATE ! X1, E

  END SUBROUTINE WriteRadiationFields1D


  SUBROUTINE WriteRadiationFields( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    INTEGER        :: FUNIT

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        RadiationSuffix // '_' // &
        FileNumberString // '.dat'

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) &
      Time, &
      nE * nNodesE, &
      nX(1) * nNodesX(1), &
      nX(2) * nNodesX(2), &
      nX(3) * nNodesX(3), &
      NodeCoordinates( MeshE, nE, nNodesE ), &
      NodeCoordinates( MeshX(1), nX(1), nNodesX(1) ), &
      NodeCoordinates( MeshX(2), nX(2), nNodesX(2) ), &
      NodeCoordinates( MeshX(3), nX(3), nNodesX(3) )
!!$    WRITE( FUNIT, * ) &
!!$      RadiationField( uCR(:,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_N,1) )

    CLOSE( FUNIT )

    PRINT*, "Wrote ", TRIM( FileName )

  END SUBROUTINE WriteRadiationFields


  FUNCTION RadiationField( u )

    REAL(DP) :: &
      RadiationField &
        (1:nNodesE   *nE,   1:nNodesX(1)*nX(1), &
         1:nNodesX(2)*nX(2),1:nNodesX(3)*nX(3))
    REAL(DP), INTENT(in) :: &
      u(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3))

    INTEGER :: iE, iX1, iX2, iX3, iNode
    INTEGER :: iNodeE, iNodeX1, iNodeX2, iNodeX3

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE
            DO iNode = 1, nDOF

              iNodeE  = NodeNumberTable(1,iNode)
              iNodeX1 = NodeNumberTable(2,iNode)
              iNodeX2 = NodeNumberTable(3,iNode)
              iNodeX3 = NodeNumberTable(4,iNode)

              RadiationField &
                ((iE-1) *nNodesE   +iNodeE, (iX1-1)*nNodesX(1)+iNodeX1, &
                 (iX2-1)*nNodesX(2)+iNodeX2,(iX3-1)*nNodesX(3)+iNodeX3) &
                = u(iNode,iE,iX1,iX2,iX3)

            END DO
          END DO
        END DO
      END DO
    END DO

    RETURN
  END FUNCTION RadiationField


  FUNCTION RadiationField1D( u, nE, nX1 )

    REAL(DP), DIMENSION(nE,nX1)            :: RadiationField1D
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: u
    INTEGER,                    INTENT(in) :: nE, nX1

    INTEGER :: iE, iX1

    DO iX1 = 1, nX1
      DO iE = 1, nE
        RadiationField1D(iE,iX1) &
          = evalL( u(:,iE,iX1), 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP )
      END DO
    END DO

    RETURN
  END FUNCTION RadiationField1D


  SUBROUTINE WriteRadiationFields1D_N( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    INTEGER        :: FUNIT

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName = OutputDirectory // '/' // TRIM( ProgramName ) // '_' // &
                 RadiationSuffix // '_' // FileNumberString // '.dat'

    ASSOCIATE( U => UnitsDisplay )

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) &
      Time / U % TimeUnit, &
      nE, nNodesE, &
      NodeCoordinatesE( nE, nNodesE ) / U % EnergyUnit, &
      nX(1), nNodesX(1), &
      NodeCoordinatesX1( nX(1), nNodesX(1) ) / U % LengthX1Unit, &
      RadiationField1D_N &
        ( uCR(:,1:nE,1:nX(1),1,1,iCR_N, 1), nE, nNodesE, nX(1), nNodesX(1) ), &
      RadiationField1D_N &
        ( uCR(:,1:nE,1:nX(1),1,1,iCR_G1,1), nE, nNodesE, nX(1), nNodesX(1) ), &
      RadiationField1D_N &
        ( uCR(:,1:nE,1:nX(1),1,1,iCR_G2,1), nE, nNodesE, nX(1), nNodesX(1) ), &
      RadiationField1D_N &
        ( uCR(:,1:nE,1:nX(1),1,1,iCR_G3,1), nE, nNodesE, nX(1), nNodesX(1) ), &
      RadiationField1D_N &
        ( uPR(:,1:nE,1:nX(1),1,1,iPR_D, 1), nE, nNodesE, nX(1), nNodesX(1) ), &
      RadiationField1D_N &
        ( uPR(:,1:nE,1:nX(1),1,1,iPR_I1,1), nE, nNodesE, nX(1), nNodesX(1) ), &
      RadiationField1D_N &
        ( uPR(:,1:nE,1:nX(1),1,1,iPR_I2,1), nE, nNodesE, nX(1), nNodesX(1) ), &
      RadiationField1D_N &
        ( uPR(:,1:nE,1:nX(1),1,1,iPR_I3,1), nE, nNodesE, nX(1), nNodesX(1) ), &
      Discontinuity(1:nE,1:nX(1),1,1)

    CLOSE( FUNIT )

    END ASSOCIATE ! U

  END SUBROUTINE WriteRadiationFields1D_N


  FUNCTION RadiationField1D_N( u, nE, nNE, nX1, nNX1 )

    REAL(DP), DIMENSION(nE*nNE,nX1*nNX1)   :: RadiationField1D_N
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: u
    INTEGER,                    INTENT(in) :: nE,  nNE
    INTEGER,                    INTENT(in) :: nX1, nNX1

    INTEGER :: i, j, iNode
    INTEGER :: iE, iNE, iX1, iNX1

    j = 1
    DO iX1 = 1, nX1
      DO iNX1 = 1, nNX1
        i = 1
        DO iE = 1, nE
          DO iNE = 1, nNE
            iNode = NodeNumber(iNE,iNX1,1,1)
            RadiationField1D_N(i,j) &
              = u(iNode,iE,iX1)
            i = i + 1
          END DO
        END DO
        j = j + 1
      END DO
    END DO

    RETURN
  END FUNCTION RadiationField1D_N


  SUBROUTINE WriteFieldsRestart1D &
               ( Time, WriteGeometryFields_Option, WriteFluidFields_Option, &
                 WriteRadiationFields_Option )

    REAL(DP), INTENT(in)           :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: WriteGeometryFields_Option
    LOGICAL,  INTENT(in), OPTIONAL :: WriteFluidFields_Option
    LOGICAL,  INTENT(in), OPTIONAL :: WriteRadiationFields_Option

    LOGICAL :: WriteGeometryFields
    LOGICAL :: WriteFluidFields
    LOGICAL :: WriteRadiationFields

    WriteGeometryFields = .FALSE.
    IF( PRESENT( WriteGeometryFields_Option ) ) &
      WriteGeometryFields = WriteGeometryFields_Option

    WriteFluidFields = .FALSE.
    IF( PRESENT( WriteFluidFields_Option ) ) &
      WriteFluidFields = WriteFluidFields_Option

    WriteRadiationFields = .FALSE.
    IF( PRESENT( WriteRadiationFields_Option ) ) &
      WriteRadiationFields = WriteRadiationFields_Option

    IF( WriteFluidFields ) &
      CALL WriteFluidFieldsRestart1D( Time )

    RestartFileNumber = RestartFileNumber + 1

  END SUBROUTINE WriteFieldsRestart1D


  SUBROUTINE WriteFluidFieldsRestart1D( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    INTEGER        :: FUNIT

    WRITE( FileNumberString, FMT='(i6.6)') RestartFileNumber

    FileName = OutputDirectory // '/' // TRIM( ProgramName ) // '_' // &
                 FluidSuffix // '_Restart_' // FileNumberString // '.dat'

    ASSOCIATE( U => UnitsDisplay )

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) RestartFileNumber
    WRITE( FUNIT, * ) FileNumber
    WRITE( FUNIT, * ) Time / U % TimeUnit
    WRITE( FUNIT, * ) nX(1)
    WRITE( FUNIT, * ) nNodesX(1)

    WRITE( FUNIT, * ) &
      NodeCoordinatesX1( nX(1), nNodesX(1) ) / U % LengthX1Unit

    WRITE( FUNIT, * ) & ! uPF_D
      FluidFieldRestart1D_Out &
        ( uPF(:,1:nX(1),1,1,iPF_D ), nX(1), nNodesX(1) ) &
          / U % MassDensityUnit

    WRITE( FUNIT, * ) & ! uPF_V1
      FluidFieldRestart1D_Out &
        ( uPF(:,1:nX(1),1,1,iPF_V1), nX(1), nNodesX(1) ) &
          / U % VelocityX1Unit

    WRITE( FUNIT, * ) & ! uPF_V2
      FluidFieldRestart1D_Out &
        ( uPF(:,1:nX(1),1,1,iPF_V2), nX(1), nNodesX(1) ) &
          / U % VelocityX2Unit

    WRITE( FUNIT, * ) & ! uPF_V3
      FluidFieldRestart1D_Out &
        ( uPF(:,1:nX(1),1,1,iPF_V3), nX(1), nNodesX(1) ) &
          / U % VelocityX3Unit

    WRITE( FUNIT, * ) & ! uPF_E
      FluidFieldRestart1D_Out &
        ( uPF(:,1:nX(1),1,1,iPF_E ), nX(1), nNodesX(1) ) &
          / U % EnergyDensityUnit

    WRITE( FUNIT, * ) & ! uPF_Ne
      FluidFieldRestart1D_Out &
        ( uPF(:,1:nX(1),1,1,iPF_Ne), nX(1), nNodesX(1) ) &
          / U % ParticleDensityUnit

    WRITE( FUNIT, * ) & ! uAF_P
      FluidFieldRestart1D_Out &
        ( uAF(:,1:nX(1),1,1,iAF_P ), nX(1), nNodesX(1) ) &
          / U % PressureUnit

    WRITE( FUNIT, * ) & ! uAF_T
      FluidFieldRestart1D_Out &
        ( uAF(:,1:nX(1),1,1,iAF_T ), nX(1), nNodesX(1) ) &
          / U % TemperatureUnit

    WRITE( FUNIT, * ) & ! uAF_Ye
      FluidFieldRestart1D_Out &
        ( uAF(:,1:nX(1),1,1,iAF_Ye), nX(1), nNodesX(1) )

    CLOSE( FUNIT )

    END ASSOCIATE ! U

    WRITE(*,*)
    WRITE(*,'(A6,A20,A)') &
      '', 'Wrote Restart File: ', TRIM( FileName )
    WRITE(*,*)

  END SUBROUTINE WriteFluidFieldsRestart1D


  FUNCTION FluidFieldRestart1D_Out( u, nX1, nNodesX1 )

    REAL(DP), DIMENSION(nX1*nNodesX1)    :: FluidFieldRestart1D_Out
    REAL(DP), DIMENSION(:,:), INTENT(in) :: u
    INTEGER,                  INTENT(in) :: nX1
    INTEGER,                  INTENT(in) :: nNodesX1

    INTEGER :: iX1, iNodeX1, iNode

    iNode = 0
    DO iX1 = 1, nX1
      DO iNodeX1 = 1, nNodesX1
        iNode = iNode + 1
        FluidFieldRestart1D_Out(iNode) = u(iNodeX1,iX1)
      END DO
    END DO

    RETURN
  END FUNCTION FluidFieldRestart1D_Out


  SUBROUTINE ReadFluidFieldsRestart1D( RestartNumber, Time )

    INTEGER,  INTENT(in)  :: RestartNumber
    REAL(DP), INTENT(out) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    INTEGER        :: FUNIT
    INTEGER        :: nElements, nNodes
    REAL(DP), DIMENSION(:), ALLOCATABLE :: RealBuffer1D

    WRITE( FileNumberString, FMT='(i6.6)') RestartNumber

    FileName = OutputDirectory // '/' // TRIM( ProgramName ) // '_' // &
                 FluidSuffix // '_Restart_' // FileNumberString // '.dat'

    WRITE(*,*)
    WRITE(*,'(A4,A22,A)') &
      '', 'Reading Restart File: ', TRIM( FileName )
    WRITE(*,*)

    ASSOCIATE( U => UnitsDisplay )

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    READ( FUNIT, * ) RestartFileNumber
    RestartFileNumber &
      = RestartFileNumber + 1

    READ( FUNIT, * ) FileNumber
    FileNumber = FileNumber

    READ( FUNIT, * ) Time
    Time = Time * U % TimeUnit

    READ( FUNIT, * ) nElements
    READ( FUNIT, * ) nNodes

    IF( nElements /= nX(1) .OR. nNodes /= nNodesX(1) )THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'Error in ReadFluidFieldsRestart1D'
      WRITE(*,'(A6,A)') '', 'Incompatible Resolutions'
      WRITE(*,'(A6,A11,I4.4,A5,I2.2)') &
        '', 'nElements: ', nX(1),  ' vs. ', nElements
      WRITE(*,'(A6,A11,I4.4,A5,I2.2)') &
        '', '   nNodes: ', nNodes, ' vs. ', nNodesX(1)
      WRITE(*,*)
      STOP
    END IF

    ALLOCATE( RealBuffer1D(nElements*nNodes) )

    READ( FUNIT, * ) RealBuffer1D ! Coordinates (Not Needed)

    READ( FUNIT, * ) RealBuffer1D ! uPF_D
    uPF(:,1:nX(1),1,1,iPF_D) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes ) &
          * U % MassDensityUnit
    READ( FUNIT, * ) RealBuffer1D ! uPF_V1
    uPF(:,1:nX(1),1,1,iPF_V1) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes ) &
          * U % VelocityX1Unit
    READ( FUNIT, * ) RealBuffer1D ! uPF_V2
    uPF(:,1:nX(1),1,1,iPF_V2) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes ) &
          * U % VelocityX2Unit
    READ( FUNIT, * ) RealBuffer1D ! uPF_V3
    uPF(:,1:nX(1),1,1,iPF_V3) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes ) &
          * U % VelocityX3Unit
    READ( FUNIT, * ) RealBuffer1D ! uPF_E
    uPF(:,1:nX(1),1,1,iPF_E) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes ) &
          * U % EnergyDensityUnit
    READ( FUNIT, * ) RealBuffer1D ! uPF_Ne
    uPF(:,1:nX(1),1,1,iPF_Ne) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes ) &
          * U % ParticleDensityUnit

    READ( FUNIT, * ) RealBuffer1D ! uAF_P
    uAF(:,1:nX(1),1,1,iAF_P) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes ) &
          * U % PressureUnit

    READ( FUNIT, * ) RealBuffer1D ! uAF_T
    uAF(:,1:nX(1),1,1,iAF_T) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes ) &
          * U % TemperatureUnit

    READ( FUNIT, * ) RealBuffer1D ! uAF_Ye
    uAF(:,1:nX(1),1,1,iAF_Ye) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes )

    DEALLOCATE( RealBuffer1D )

    CLOSE( FUNIT )

    END ASSOCIATE ! U

  END SUBROUTINE ReadFluidFieldsRestart1D


  FUNCTION FluidFieldRestart1D_In( u, nX1, nNodesX1 )

    REAL(DP), DIMENSION(nNodesX1,nX1)  :: FluidFieldRestart1D_In
    REAL(DP), DIMENSION(:), INTENT(in) :: u
    INTEGER,                INTENT(in) :: nX1
    INTEGER,                INTENT(in) :: nNodesX1

    INTEGER :: iX1, iNodeX1, iNode

    iNode = 0
    DO iX1 = 1, nX1
      DO iNodeX1 = 1, nNodesX1
        iNode = iNode + 1
        FluidFieldRestart1D_In(iNodeX1,iX1) = u(iNode)
      END DO
    END DO

    RETURN
  END FUNCTION FluidFieldRestart1D_In


  FUNCTION NodeCoordinates( Mesh, nElements, nNodes )

    REAL(DP) :: NodeCoordinates(nElements*nNodes)
    TYPE(MeshType), INTENT(in) :: Mesh
    INTEGER,        INTENT(in) :: nElements
    INTEGER,        INTENT(in) :: nNodes

    INTEGER :: i, j, iNode

    iNode = 0
    DO j = 1, nElements
      DO i = 1, nNodes
        iNode = iNode + 1
        NodeCoordinates(iNode) &
          = NodeCoordinate( Mesh, j, i )
      END DO
    END DO

    RETURN
  END FUNCTION NodeCoordinates


  FUNCTION NodeCoordinatesE( nE, nNodesE )

    REAL(DP), DIMENSION(nE*nNodesE) :: NodeCoordinatesE
    INTEGER,             INTENT(in) :: nE
    INTEGER,             INTENT(in) :: nNodesE

    INTEGER :: iE, iNodeE, iNode

    iNode = 0
    DO iE = 1, nE
      DO iNodeE = 1, nNodesE
        iNode = iNode + 1
        NodeCoordinatesE(iNode) &
          = NodeCoordinate( MeshE, iE, iNodeE )
      END DO
    END DO

    RETURN
  END FUNCTION NodeCoordinatesE


  FUNCTION NodeCoordinatesX1( nX1, nNodesX1 )

    REAL(DP), DIMENSION(nX1*nNodesX1) :: NodeCoordinatesX1
    INTEGER,               INTENT(in) :: nX1
    INTEGER,               INTENT(in) :: nNodesX1

    INTEGER :: iX1, iNodeX1, iNode

    iNode = 0
    DO iX1 = 1, nX1
      DO iNodeX1 = 1, nNodesX1
        iNode = iNode + 1
        NodeCoordinatesX1(iNode) &
          = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
      END DO
    END DO

    RETURN
  END FUNCTION NodeCoordinatesX1


END MODULE InputOutputModule
