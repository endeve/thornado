MODULE ScalarWave_InputOutputModuleHDF

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nE, nNodesE, nDOFE, &
    nX, nNodesX, nDOFX, &
    nDOF
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshE, MeshX
  USE InputOutputUtilitiesModule, ONLY: &
    NodeCoordinates, &
    Field3D, &
    FromField3D, &
    Field4D
  USE ScalarFieldsModule, ONLY: &
    uSF, nSF, &
    unitsSF, &
    namesSF

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  CHARACTER (9), PARAMETER :: &
    OUTPUTDIRECTORY = '../Output'
  CHARACTER(10), PARAMETER :: &
    ScalarWaveSuffix  = 'ScalarWave'
  INTEGER :: FileNumber = 0

  INTEGER :: HDFERR

  PUBLIC :: WriteFieldsHDF
  PUBLIC :: ReadFieldsHDF

CONTAINS

  SUBROUTINE WriteFieldsHDF &
    ( Time, WriteSF_Option )

    REAL(DP), INTENT(in) :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: WriteSF_Option

    LOGICAL :: WriteSF

    IF( PRESENT( WriteSF_Option ) )THEN
      WriteSF = WriteSF_Option
    ELSE
      WriteSF = .FALSE.
    END IF

    IF( WriteSF )THEN

      CALL WriteScalarFieldsHDF( Time )

    END IF

    FileNumber = FileNumber + 1
  
  END SUBROUTINE WriteFieldsHDF


  SUBROUTINE ReadFieldsHDF &
    ( ReadFileNumber, Time, ReadSF_Option )

    INTEGER,  INTENT(in)  :: &
      ReadFileNumber
    REAL(DP), INTENT(out) :: &
      Time
    LOGICAL,  INTENT(in), OPTIONAL :: &
      ReadSF_Option

    LOGICAL :: ReadSF

    FileNumber = ReadFileNumber

    ReadSF = .FALSE.
    IF( PRESENT( ReadSF_Option ) ) ReadSF = ReadSF_Option

    IF( ReadSF ) CALL ReadScalarFieldsHDF( Time )

    FileNumber = FileNumber + 1

  END SUBROUTINE ReadFieldsHDF


  SUBROUTINE WriteScalarFieldsHDF( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: GroupName
    CHARACTER(256) :: DatasetName
    INTEGER        :: iSF
    INTEGER(HID_T) :: FILE_ID
    
    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        ScalarWaveSuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5OPEN_F( HDFERR )

    CALL H5FCREATE_F( TRIM( FileName ), H5F_ACC_TRUNC_F, FILE_ID, HDFERR )

    ASSOCIATE( U => UnitsDisplay )

    ! --- Write Time ---

    DatasetName = '/Time'

    CALL WriteDataset1DHDF &
           ( [ Time ] / U % TimeUnit, DatasetName, FILE_ID )

    ! --- Write Spatial Grid ---

    GroupName = 'Spatial Grid'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ) , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X1'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(1),nX(1),nNodesX(1)) &
               / U % LengthX1Unit, DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X2'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(2),nX(2),nNodesX(2)) &
               / U % LengthX2Unit, DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X3'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(3),nX(3),nNodesX(3)) &
               / U % LengthX3Unit, DatasetName, FILE_ID )

    END ASSOCIATE ! U

    ! --- Write Scalar Field Variables ---

    GroupName = 'Scalar Fields'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

    DO iSF = 1, nSF

      DatasetName = TRIM( GroupName ) // '/' // TRIM( namesSF(iSF) )

      CALL WriteDataset3DHDF &
             ( Field3D &
                 ( uSF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iSF), nX, nNodesX, &
                   nDOFX, NodeNumberTableX ) / unitsSF(iSF), &
               DatasetName, FILE_ID )

    END DO

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

  END SUBROUTINE WriteScalarFieldsHDF


  SUBROUTINE ReadScalarFieldsHDF( Time )

    REAL(DP), INTENT(out) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: DatasetName
    CHARACTER(256) :: GroupName
    INTEGER(HID_T) :: FILE_ID
    INTEGER        :: iSF
    REAL(DP)       :: Dataset1D(1)
    REAL(DP)       :: Dataset3D(nX(1)*nNodesX(1), &
                                nX(2)*nNodesX(2), &
                                nX(3)*nNodesX(3))

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        ScalarWaveSuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5OPEN_F( HDFERR )

    CALL H5FOPEN_F( TRIM( FileName ), H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

    ASSOCIATE( U => UnitsDisplay )

    ! --- Read Time ---

    DatasetName = '/Time'

    CALL ReadDataset1DHDF( Dataset1D, DatasetName, FILE_ID )

    Time = Dataset1D(1) * U % TimeUnit

    END ASSOCIATE

    ! --- Read Scalar Field Variables ---

    GroupName = 'Scalar Fields'

    DO iSF = 1, nSF

      DatasetName = TRIM( GroupName ) // '/' // TRIM( namesSF(iSF) )

      CALL ReadDataset3DHDF( Dataset3D, DatasetName, FILE_ID )

      uSF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iSF) &
        = FromField3D( Dataset3D, nX, nNodesX, nDOFX, NodeNumberTableX )

    END DO

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

  END SUBROUTINE ReadScalarFieldsHDF

  
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


  SUBROUTINE ReadDataset1DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(out) :: Dataset(:)
    CHARACTER(LEN=*), INTENT(in)  :: DatasetName
    INTEGER(HID_T),   INTENT(in)  :: FILE_ID

    INTEGER(HID_T) :: DATASET_ID
    INTEGER(HID_T) :: DATASIZE(1)

    DATASIZE = SHAPE( Dataset )

    CALL H5DOPEN_F( FILE_ID, TRIM( DatasetName ), DATASET_ID, HDFERR )

    CALL H5DREAD_F( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE ReadDataset1DHDF


  SUBROUTINE WriteDataset3DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(in) :: Dataset(:,:,:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HSIZE_T) :: DATASIZE(3)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SHAPE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 3, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    CALL H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5SCLOSE_F( DATASPACE_ID, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset3DHDF


  SUBROUTINE ReadDataset3DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(out) :: Dataset(:,:,:)
    CHARACTER(LEN=*), INTENT(in)  :: DatasetName
    INTEGER(HID_T),   INTENT(in)  :: FILE_ID

    INTEGER(HID_T) :: DATASET_ID
    INTEGER(HID_T) :: DATASIZE(3)

    DATASIZE = SHAPE( Dataset )

    CALL H5DOPEN_F( FILE_ID, TRIM( DatasetName ), DATASET_ID, HDFERR )

    CALL H5DREAD_F( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE ReadDataset3DHDF


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


END MODULE ScalarWave_InputOutputModuleHDF
