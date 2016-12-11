MODULE InputOutputModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nE
  USE PolynomialBasisModule_Lagrange, ONLY: &
    evalL, &
    evalLX
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
         iAF_Me, iAF_Mp, iAF_Mn, iAF_Gm, iAF_Cs
  USE RadiationFieldsModule, ONLY: &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3

  IMPLICIT NONE
  PRIVATE

  CHARACTER(9),  PARAMETER :: &
    OutputDirectory = '../Output'
  CHARACTER(11), PARAMETER :: &
    FluidSuffix = 'FluidFields'
  CHARACTER(15), PARAMETER :: &
    RadiationSuffix = 'RadiationFields'
  INTEGER :: FileNumber = 0

  PUBLIC :: WriteFields1D

CONTAINS


  SUBROUTINE WriteFields1D &
               ( Time, WriteFluidFields_Option, WriteRadiationFields_Option )

    REAL(DP), INTENT(in)           :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: WriteFluidFields_Option
    LOGICAL,  INTENT(in), OPTIONAL :: WriteRadiationFields_Option

    LOGICAL :: WriteFluidFields
    LOGICAL :: WriteRadiationFields

    WriteFluidFields = .TRUE.
    IF( PRESENT( WriteFluidFields_Option ) ) &
      WriteFluidFields = WriteFluidFields_Option

    WriteRadiationFields = .TRUE.
    IF( PRESENT( WriteRadiationFields_Option ) ) &
      WriteRadiationFields = WriteRadiationFields_Option

    IF( WriteFluidFields ) &
      CALL WriteFluidFields1D( Time )

    IF( WriteRadiationFields ) &
      CALL WriteRadiationFields1D( Time )

    FileNumber = FileNumber + 1

  END SUBROUTINE WriteFields1D


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
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Gm), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Cs), nX(1) ) &
        / U % VelocityUnit

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
      RadiationField1D( uCR(:,1:nE,1:nX(1),1,1,iPR_I3,1), nE, nX(1) )

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


END MODULE InputOutputModule
