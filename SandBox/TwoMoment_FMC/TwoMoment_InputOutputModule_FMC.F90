MODULE TwoMoment_InputOutputModule_FMC

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay, &
    Centimeter,   &
    Second
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
    Field4D, &
    FromField4D
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nSpecies, &
    uCM, nCM, namesCM, unitsCM, &
    uPM, nPM, namesPM, unitsPM, &
    uAM, nAM, namesAM, unitsAM, &
    uGM, nGM, namesGM, unitsGM

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  CHARACTER(9),  PARAMETER :: &
    OutputDirectory = '../Output'
  CHARACTER(15), PARAMETER :: &
    TwoMomentSuffix = 'TwoMomentFields'
  INTEGER :: FileNumber = 0

  INTEGER :: HDFERR

  PUBLIC :: WriteTwoMomentFieldsHDF
  PUBLIC :: ReadTwoMomentFieldsHDF

CONTAINS


  SUBROUTINE WriteTwoMomentFieldsHDF( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(2)   :: String2
    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: GroupName
    CHARACTER(256) :: GroupNameSpecies
    CHARACTER(256) :: DatasetName
    INTEGER        :: iS, iMF
    INTEGER(HID_T) :: FILE_ID

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( uCM, uPM, uAM, uGM )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( uCM, uPM, uAM, uGM )
#endif

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        TwoMomentSuffix // '_' // &
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

    ! --- Write Cell Center Coordinates ---

    DatasetName = TRIM( GroupName ) // '/X1_C'

    CALL WriteDataset1DHDF &
           ( MeshX(1) % Center(1:nX(1)) / U % LengthX1Unit, &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X2_C'

    CALL WriteDataset1DHDF &
           ( MeshX(2) % Center(1:nX(2)) / U % LengthX2Unit, &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X3_C'

    CALL WriteDataset1DHDF &
           ( MeshX(3) % Center(1:nX(3)) / U % LengthX3Unit, &
             DatasetName, FILE_ID )

    ! --- Write Cell Widths ---

    DatasetName = TRIM( GroupName ) // '/dX1'

    CALL WriteDataset1DHDF &
           ( MeshX(1) % Width(1:nX(1)) / U % LengthX1Unit, &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/dX2'

    CALL WriteDataset1DHDF &
           ( MeshX(2) % Width(1:nX(2)) / U % LengthX2Unit, &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/dX3'

    CALL WriteDataset1DHDF &
           ( MeshX(3) % Width(1:nX(3)) / U % LengthX3Unit, &
             DatasetName, FILE_ID )

    ! --- Write Energy Grid ---

    GroupName = 'Energy Grid'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

    DatasetName = TRIM( GroupName ) // '/E'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshE,nE,nNodesE) &
               / U % EnergyUnit, DatasetName, FILE_ID )

    ! --- Write Cell Center Coordinates ---

    DatasetName = TRIM( GroupName ) // '/E_C'

    CALL WriteDataset1DHDF &
           ( MeshE % Center(1:nE) / U % EnergyUnit, &
             DatasetName, FILE_ID )

    ! --- Write Cell Widths ---

    DatasetName = TRIM( GroupName ) // '/dE'

    CALL WriteDataset1DHDF &
           ( MeshE % Width(1:nE) / U % EnergyUnit, &
             DatasetName, FILE_ID )

    END ASSOCIATE ! U

    ! --- Write Two-Moment Variables ---

    GroupName = 'Two-Moment Fields'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

    DO iS = 1, nSpecies

      WRITE( String2, FMT='(i2.2)') iS

      GroupNameSpecies = 'Two-Moment Fields/' // 'Species_' // String2

      CALL CreateGroupHDF( FileName, TRIM( GroupNameSpecies ), FILE_ID )

      ! --- Conserved ---

      GroupName = TRIM( GroupNameSpecies ) // '/Conserved'

      CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

      DO iMF = 1, nCM

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesCM(iMF) )

        CALL WriteDataset4DHDF &
               ( Field4D &
                   ( uCM(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iMF,iS), &
                     [ nE, nX(1), nX(2), nX(3) ],                     &
                     [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
                     nDOF, NodeNumberTable ), DatasetName, FILE_ID )

      END DO

      ! --- Primitive ---

      GroupName = TRIM( GroupNameSpecies  ) // '/Primitive'

      CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

      DO iMF = 1, nPM

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesPM(iMF) )

        CALL WriteDataset4DHDF &
               ( Field4D &
                   ( uPM(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iMF,iS), &
                     [ nE, nX(1), nX(2), nX(3) ],                     &
                     [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
                     nDOF, NodeNumberTable ), DatasetName, FILE_ID )

      END DO

      ! --- Auxiliary ---

      GroupName = TRIM( GroupNameSpecies  ) // '/Auxiliary'

      CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

      DO iMF = 1, nAM

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesAM(iMF) )

        CALL WriteDataset4DHDF &
               ( Field4D &
                   ( uAM(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iMF,iS), &
                     [ nE, nX(1), nX(2), nX(3) ],                     &
                     [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
                     nDOF, NodeNumberTable ), DatasetName, FILE_ID )

      END DO

      ! --- Gray ---

      GroupName = TRIM( GroupNameSpecies ) // '/Gray'

      CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

      DO iMF = 1, nGM

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesGM(iMF) )

        CALL WriteDataset3DHDF &
             ( Field3D &
                 ( uGM(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iMF,iS), &
                   nX, nNodesX, nDOFX, NodeNumberTableX ) / unitsGM(iMF), &
               DatasetName, FILE_ID )

      END DO

    END DO

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

    FileNumber = FileNumber + 1

  END SUBROUTINE WriteTwoMomentFieldsHDF


  SUBROUTINE ReadTwoMomentFieldsHDF( ReadFileNumber, Time )

    INTEGER,  INTENT(in)  :: ReadFileNumber
    REAL(DP), INTENT(out) :: Time

    CHARACTER(2)   :: String2
    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: DatasetName
    CHARACTER(256) :: GroupName
    CHARACTER(256) :: GroupNameSpecies
    INTEGER(HID_T) :: FILE_ID
    INTEGER        :: iS, iMF
    REAL(DP)       :: Dataset1D(1)
    REAL(DP)       :: Dataset4D(nE   *nNodesE,    &
                                nX(1)*nNodesX(1), &
                                nX(2)*nNodesX(2), &
                                nX(3)*nNodesX(3))

    FileNumber = ReadFileNumber

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        TwoMomentSuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5OPEN_F( HDFERR )

    CALL H5FOPEN_F( TRIM( FileName ), H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

    ASSOCIATE( U => UnitsDisplay )

    ! --- Read Time ---

    DatasetName = '/Time'

    CALL ReadDataset1DHDF( Dataset1D, DatasetName, FILE_ID )

    Time = Dataset1D(1) * U % TimeUnit

    END ASSOCIATE

    ! --- Read Two-Moment Variables ---

    DO iS = 1, nSpecies

      WRITE( String2, FMT='(i2.2)') iS

      GroupNameSpecies = 'Two-Moment Fields/' // 'Species_' // String2

      ! --- Conserved ---

      GroupName = TRIM( GroupNameSpecies ) // '/Conserved'

      DO iMF = 1, nCM

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesCM(iMF) )

        CALL ReadDataset4DHDF( Dataset4D, DatasetName, FILE_ID )

        uCM(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iMF,iS) &
          = FromField4D &
              ( Dataset4D, [ nE, nX(1), nX(2), nX(3) ], &
                [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
                nDOF, NodeNumberTable )

      END DO

      ! --- Primitive ---

      GroupName = TRIM( GroupNameSpecies ) // '/Primitive'

      DO iMF = 1, nPM

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesPM(iMF) )

        CALL ReadDataset4DHDF( Dataset4D, DatasetName, FILE_ID )

        uPM(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iMF,iS) &
          = FromField4D &
              ( Dataset4D, [ nE, nX(1), nX(2), nX(3) ], &
                [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
                nDOF, NodeNumberTable )

      END DO

    END DO

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( uCM, uPM )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( uCM, uPM )
#endif

    FileNumber = FileNumber + 1

  END SUBROUTINE ReadTwoMomentFieldsHDF


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


  SUBROUTINE ReadDataset4DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(out) :: Dataset(:,:,:,:)
    CHARACTER(LEN=*), INTENT(in)  :: DatasetName
    INTEGER(HID_T),   INTENT(in)  :: FILE_ID

    INTEGER(HID_T) :: DATASET_ID
    INTEGER(HID_T) :: DATASIZE(4)

    DATASIZE = SHAPE( Dataset )

    CALL H5DOPEN_F( FILE_ID, TRIM( DatasetName ), DATASET_ID, HDFERR )

    CALL H5DREAD_F( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE ReadDataset4DHDF


END MODULE TwoMoment_InputOutputModule_FMC
