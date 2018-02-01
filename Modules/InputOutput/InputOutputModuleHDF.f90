MODULE InputOutputModuleHDF

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nE, nNodesE, &
    nX, nNodesX, &
    nDOF
  USE ReferenceElementModule_Beta, ONLY: &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshE, MeshX
  USE InputOutputUtilitiesModule, ONLY: &
    NodeCoordinates, &
    Field4D
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, nCR, namesCR, &
    uPR, nPR, namesPR

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  CHARACTER(9),  PARAMETER :: &
    OutputDirectory = '../Output'
  CHARACTER(15), PARAMETER :: &
    RadiationSuffix = 'RadiationFields'
  INTEGER :: FileNumber = 0

  INTEGER :: HDFERR

  PUBLIC :: WriteFieldsHDF

CONTAINS


  SUBROUTINE WriteFieldsHDF &
              ( Time, WriteGF_Option, WriteFF_Option, WriteRF_Option )

    REAL(DP), INTENT(in) :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: WriteGF_Option
    LOGICAL,  INTENT(in), OPTIONAL :: WriteFF_Option
    LOGICAL,  INTENT(in), OPTIONAL :: WriteRF_Option

    LOGICAL :: WriteGF
    LOGICAL :: WriteFF
    LOGICAL :: WriteRF

    WriteGF = .FALSE.
    IF( PRESENT( WriteGF_Option ) ) &
      WriteGF = WriteGF_Option

    WriteFF = .FALSE.
    IF( PRESENT( WriteFF_Option ) ) &
      WriteFF = WriteFF_Option

    WriteRF = .FALSE.
    IF( PRESENT( WriteRF_Option ) ) &
      WriteRF = WriteRF_Option

    IF( WriteRF )THEN

      CALL WriteRadiationFieldsHDF( Time )

    END IF

    FileNumber = FileNumber + 1

  END SUBROUTINE WriteFieldsHDF


  SUBROUTINE WriteRadiationFieldsHDF( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(2)   :: String2
    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: GroupName
    CHARACTER(256) :: GroupNameSpecies
    CHARACTER(256) :: DatasetName
    INTEGER        :: iS, iRF
    INTEGER(HID_T) :: FILE_ID

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        RadiationSuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5OPEN_F( HDFERR )

    CALL H5FCREATE_F( TRIM( FileName ), H5F_ACC_TRUNC_F, FILE_ID, HDFERR )

    ! --- Write Time ---

    DatasetName = '/Time'

    CALL WriteDataset1DHDF &
           ( [ Time ], DatasetName, FILE_ID )

    ! --- Write Spatial Grid ---

    GroupName = 'Spatial Grid'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ) , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X1'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(1),nX(1),nNodesX(1)), &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X2'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(2),nX(2),nNodesX(2)), &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X3'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(3),nX(3),nNodesX(3)), &
             DatasetName, FILE_ID )

    ! --- Write Energy Grid ---

    GroupName = 'Energy Grid'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

    DatasetName = TRIM( GroupName ) // '/E'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshE,nE,nNodesE), DatasetName, FILE_ID )

    ! --- Write Radiation Variables ---

    GroupName = 'Radiation Fields'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

    DO iS = 1, nSpecies

      WRITE( String2, FMT='(i2.2)') iS

      GroupNameSpecies = 'Radiation Fields/' // 'Species_' // String2

      CALL CreateGroupHDF( FileName, TRIM( GroupNameSpecies ), FILE_ID )

      ! --- Conserved ---

      GroupName = TRIM( GroupNameSpecies ) // '/Conserved'

      CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

      DO iRF = 1, nCR

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesCR(iRF) )

        CALL WriteDataset4DHDF &
               ( Field4D &
                   ( uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iRF,iS), &
                     [ nE, nX(1), nX(2), nX(3) ],                     &
                     [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
                     nDOF, NodeNumberTable ), DatasetName, FILE_ID )

      END DO

      ! --- Primitive ---

      GroupName = TRIM( GroupNameSpecies  ) // '/Primitive'

      CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

      DO iRF = 1, nPR

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesPR(iRF) )

        CALL WriteDataset4DHDF &
               ( Field4D &
                   ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iRF,iS), &
                     [ nE, nX(1), nX(2), nX(3) ],                     &
                     [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
                     nDOF, NodeNumberTable ), DatasetName, FILE_ID )

      END DO

    END DO

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

  END SUBROUTINE WriteRadiationFieldsHDF


  SUBROUTINE CreateGroupHDF( FileName, GroupName, FILE_ID )

    CHARACTER(len=*), INTENT(in) :: FileName
    CHARACTER(len=*), INTENT(in) :: GroupName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HID_T) :: GROUP_ID

    CALL H5GCREATE_F( FILE_ID, TRIM( GroupName ), GROUP_ID, HDFERR )

    CALL H5GCLOSE_F( GROUP_ID, HDFERR )

  END SUBROUTINE CreateGroupHDF


  SUBROUTINE WriteDataset1DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(in) :: Dataset(:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HSIZE_T) :: DATASIZE(1)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SIZE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 1, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    call H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    call H5SCLOSE_F( DATASPACE_ID, HDFERR )

    call H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset1DHDF


  SUBROUTINE WriteDataset4DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(in) :: Dataset(:,:,:,:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HSIZE_T) :: DATASIZE(4)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SHAPE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 4, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    call H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    call H5SCLOSE_F( DATASPACE_ID, HDFERR )

    call H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset4DHDF


END MODULE InputOutputModuleHDF
