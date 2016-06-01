MODULE InputOutputModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX
  USE PolynomialBasisModule_Lagrange, ONLY: &
    evalLX
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
         iAF_Me, iAF_Mp, iAF_Mn, iAF_Gm, iAF_Cs

  IMPLICIT NONE
  PRIVATE

  CHARACTER(9),  PARAMETER :: &
    OutputDirectory = '../Output'
  CHARACTER(11), PARAMETER :: &
    FluidSuffix = 'FluidFields'
  INTEGER :: FileNumber = 0

  PUBLIC :: WriteFields1D

CONTAINS


  SUBROUTINE WriteFields1D

    CALL WriteFluidFields1D

    FileNumber = FileNumber + 1

  END SUBROUTINE WriteFields1D


  SUBROUTINE WriteFluidFields1D

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    INTEGER        :: FUNIT

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName = OutputDirectory // '/' // TRIM( ProgramName ) // '_' // &
                 FluidSuffix // '_' // FileNumberString // '.dat'

    ASSOCIATE( X1 => MeshX(1) % Center(1:nX(1)) )

    OPEN( UNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) nX(1), X1, &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_D ), nX(1) ), &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_S1), nX(1) ), &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_S2), nX(1) ), &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_S3), nX(1) ), &
      FluidField1D( uCF(:,1:nX(1),1,1,iCF_E ), nX(1) ), &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_D ), nX(1) ), &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_V1), nX(1) ), &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_V2), nX(1) ), &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_V3), nX(1) ), &
      FluidField1D( uPF(:,1:nX(1),1,1,iPF_E ), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_P ), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_T ), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Ye), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_S ), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_E ), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Me), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Mp), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Mn), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Gm), nX(1) ), &
      FluidField1D( uAF(:,1:nX(1),1,1,iAF_Cs), nX(1) )

    CLOSE( FUNIT )

    END ASSOCIATE

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


END MODULE InputOutputModule
