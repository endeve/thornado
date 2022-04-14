MODULE ProgenitorModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Centimeter, &
    Gram, &
    Second, &
    Kelvin
  USE UtilitiesModule, ONLY: &
    WriteVector

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  INTEGER :: hdferr
  TYPE, PUBLIC :: ProgenitorType1D
    INTEGER :: nPoints
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Radius
    REAL(DP), DIMENSION(:), ALLOCATABLE :: MassDensity
    REAL(DP), DIMENSION(:), ALLOCATABLE :: RadialVelocity
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Temperature
    REAL(DP), DIMENSION(:), ALLOCATABLE :: ElectronFraction
  END type ProgenitorType1D

  PUBLIC :: ReadProgenitor1D

CONTAINS


  SUBROUTINE ReadProgenitor1D( FileName, P1D, Verbose_Option )

    CHARACTER(LEN=*),       INTENT(in)  :: FileName
    TYPE(ProgenitorType1D), INTENT(out) :: P1D
    LOGICAL,                INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER(HID_T) :: FileID
    INTEGER, DIMENSION(3) :: Dims
    LOGICAL :: Verbose

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(6x,A,A)') &
        'Reading Progenitor File: ', TRIM( FileName )
      WRITE(*,*)

    END IF

    CALL H5open_f( hdferr )

    CALL H5Fopen_f &
           ( TRIM( '../Progenitors/' // FileName ), &
             H5F_ACC_RDONLY_F, FileID, hdferr )

    CALL ReadInteger1D( FileID, './mesh/array_dimensions', Dims )

    CALL CreateProgenitor1D( Dims(1), P1D )

    CALL ReadFloat1D &
           ( FileID, './mesh/x_cf', P1D % Radius )
    CALL ReadFloat3Dto1D &
           ( FileID, './fluid/rho_c', Dims, P1D % MassDensity )
    CALL ReadFloat3Dto1D &
           ( FileID, './fluid/u_c',   Dims, P1D % RadialVelocity )
    CALL ReadFloat3Dto1D &
           ( FileID, './fluid/t_c',   Dims, P1D % Temperature )
    CALL ReadFloat3Dto1D &
           ( FileID, './fluid/ye_c',  Dims, P1D % ElectronFraction )

    CALL H5Fclose_f( FileID, hdferr )

    CALL H5close_f( hdferr )

    ! --- Convert to Code Units ---

    P1D % Radius &
      = P1D % Radius * Centimeter
    P1D % MassDensity &
      = P1D % MassDensity * Gram / Centimeter**3
    P1D % RadialVelocity &
      = P1D % RadialVelocity * Centimeter / Second
    P1D % Temperature &
      = P1D % Temperature * Kelvin
    P1D % ElectronFraction &
      = P1D % ElectronFraction

  END SUBROUTINE ReadProgenitor1D


  SUBROUTINE CreateProgenitor1D( nPoints, Progenitor1D )

    INTEGER,                INTENT(in)  :: nPoints
    TYPE(ProgenitorType1D), INTENT(out) :: Progenitor1D

    Progenitor1D % nPoints = nPoints
    ALLOCATE( Progenitor1D % Radius          (nPoints) )
    ALLOCATE( Progenitor1D % MassDensity     (nPoints) )
    ALLOCATE( Progenitor1D % RadialVelocity  (nPoints) )
    ALLOCATE( Progenitor1D % Temperature     (nPoints) )
    ALLOCATE( Progenitor1D % ElectronFraction(nPoints) )

  END SUBROUTINE CreateProgenitor1D


  SUBROUTINE ReadInteger1D( FileID, Name, Values )

    INTEGER(HID_T),        INTENT(in)  :: FileID
    CHARACTER(LEN=*),      INTENT(in)  :: Name
    INTEGER, DIMENSION(:), INTENT(out) :: values

    INTEGER(HID_T) :: &
      DatasetID
    INTEGER(HSIZE_T), DIMENSION(1) :: &
      Size

    Size = SHAPE( Values )

    CALL H5Dopen_f &
           ( FileID, Name, DatasetID, hdferr )
    CALL H5Dread_f &
           ( DatasetID, H5T_NATIVE_INTEGER, Values, Size, hdferr )
    CALL H5Dclose_f &
           ( DatasetID, hdferr )

  END SUBROUTINE ReadInteger1D


  SUBROUTINE ReadFloat1D( FileID, Name, Values )

    INTEGER(HID_T),         INTENT(in)  :: FileID
    CHARACTER(LEN=*),       INTENT(in)  :: Name
    REAL(DP), DIMENSION(:), INTENT(out) :: values

    INTEGER(HID_T) :: &
      DatasetID
    INTEGER(HSIZE_T), DIMENSION(1) :: &
      Size

    Size = SHAPE( Values )

    CALL H5Dopen_f &
           ( FileID, Name, DatasetID, hdferr )
    call H5Dread_f &
           ( DatasetID, H5T_NATIVE_DOUBLE, Values, Size, hdferr )
    CALL H5Dclose_f &
           ( DatasetID, hdferr )

  END SUBROUTINE ReadFloat1D


  SUBROUTINE ReadFloat3Dto1D( FileID, Name, Dims, Values )

    INTEGER(HID_T),         INTENT(in)  :: FileID
    CHARACTER(LEN=*),       INTENT(in)  :: Name
    INTEGER,  DIMENSION(3), INTENT(in)  :: Dims
    REAL(DP), DIMENSION(:), INTENT(out) :: values

    INTEGER(HID_T) :: &
      DatasetID
    INTEGER(HSIZE_T), DIMENSION(3) :: &
      Size
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
      Dataset3D

    Size = Dims

    ALLOCATE( Dataset3D(1:Dims(1),1:Dims(2),1:Dims(3)) )

    CALL H5Dopen_f &
           ( FileID, Name, DatasetID, hdferr )
    call H5Dread_f &
           ( DatasetID, H5T_NATIVE_DOUBLE, Dataset3D, Size, hdferr )
    CALL H5Dclose_f &
           ( DatasetID, hdferr )

    Values(1:Dims(1)) &
      = Dataset3D(1:Dims(1),1,1)

    DEALLOCATE( Dataset3D )

  END SUBROUTINE ReadFloat3Dto1D


END MODULE ProgenitorModule
