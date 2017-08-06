MODULE InputOutputModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, nE
  USE PolynomialBasisModule_Lagrange, ONLY: &
    evalL, &
    evalLX
  USE MeshModule, ONLY: &
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
      CALL WriteFluidFields1D( Time )

    IF( WriteRadiationFields ) &
      CALL WriteRadiationFields1D( Time )

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
      Time / U % TimeUnit, nX(1), X1 / U % LengthUnit, &
      GeometryField1D( uGF(:,1:nX(1),1,1,iGF_Phi_N    ), nX(1) ) &
        / ( U % EnergyDensityUnit / U % MassDensityUnit ), &
      GeometryField1D( uGF(:,1:nX(1),1,1,iGF_Gm_dd_11 ), nX(1) ), &
      GeometryField1D( uGF(:,1:nX(1),1,1,iGF_Gm_dd_22 ), nX(1) ), &
      GeometryField1D( uGF(:,1:nX(1),1,1,iGF_Gm_dd_33 ), nX(1) )

    CLOSE( FUNIT )

    END ASSOCIATE ! X1, etc.

  END SUBROUTINE WriteGeometryFields1D


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
      Time / U % TimeUnit, nX(1), X1 / U % LengthUnit, &
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
        / U % VelocityUnit, &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_V2), nX(1) ) &
        / U % VelocityUnit, &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_V3), nX(1) ) &
        / U % VelocityUnit, &
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
        / U % VelocityUnit, &
      Shock(1:nX(1),1,1)

    CLOSE( FUNIT )

    END ASSOCIATE ! X1, etc.

  END SUBROUTINE WriteFluidFields1D


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
      nE, E / U % EnergyUnit, nX(1), X1 / U % LengthUnit, &
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
      NodeCoordinatesX1( nX(1), nNodesX(1) ) / U % LengthUnit

    WRITE( FUNIT, * ) & ! uPF_D
      FluidFieldRestart1D_Out &
        ( uPF(:,1:nX(1),1,1,iPF_D ), nX(1), nNodesX(1) ) &
          / U % MassDensityUnit

    WRITE( FUNIT, * ) & ! uPF_V1
      FluidFieldRestart1D_Out &
        ( uPF(:,1:nX(1),1,1,iPF_V1), nX(1), nNodesX(1) ) &
          / U % VelocityUnit

    WRITE( FUNIT, * ) & ! uPF_V2
      FluidFieldRestart1D_Out &
        ( uPF(:,1:nX(1),1,1,iPF_V2), nX(1), nNodesX(1) ) &
          / U % VelocityUnit

    WRITE( FUNIT, * ) & ! uPF_V3
      FluidFieldRestart1D_Out &
        ( uPF(:,1:nX(1),1,1,iPF_V3), nX(1), nNodesX(1) ) &
          / U % VelocityUnit

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
          * U % VelocityUnit
    READ( FUNIT, * ) RealBuffer1D ! uPF_V2
    uPF(:,1:nX(1),1,1,iPF_V2) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes ) &
          * U % VelocityUnit
    READ( FUNIT, * ) RealBuffer1D ! uPF_V3
    uPF(:,1:nX(1),1,1,iPF_V3) &
      = FluidFieldRestart1D_In( RealBuffer1D, nElements, nNodes ) &
          * U % VelocityUnit
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


END MODULE InputOutputModule
