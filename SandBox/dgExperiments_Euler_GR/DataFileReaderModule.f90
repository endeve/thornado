MODULE DataFileReaderModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Meter, Kilogram, Second, Joule
  USE UtilitiesModule, ONLY: &
    Locate
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ReadParameters
  PUBLIC :: ReadData

CONTAINS

  SUBROUTINE ReadParameters( FILEIN, FluidFieldParameters )

    CHARACTER( LEN = * ), INTENT(in)   :: FILEIN
    REAL(DP), INTENT(out), ALLOCATABLE :: FluidFieldParameters(:)
    INTEGER                            :: i, nParams

    ! --- Get number of parameters ---
    nParams = 0
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    DO
      READ( 100, *, END = 10 )
      nParams = nParams + 1
    END DO
    10 CLOSE( 100 )

    ! --- Allocate and read in parameters ---
    ALLOCATE( FluidFieldParameters(nParams) )
    
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    DO i = 1, nParams
       READ( 100, '(ES23.16E2)' ) FluidFieldParameters(i)
    END DO
    CLOSE( 100 )

    ! --- Convert from physical-units to code-units ---
    FluidFieldParameters(1) = FluidFieldParameters(1) * Kilogram
    FluidFieldParameters(2) = FluidFieldParameters(2)
    FluidFieldParameters(3) = FluidFieldParameters(3) * Meter
    FluidFieldParameters(4) = FluidFieldParameters(4) * Meter
    FluidFieldParameters(5) = FluidFieldParameters(5) * Meter
    FluidFieldParameters(6) = FluidFieldParameters(6) * Meter
    FluidFieldParameters(7) = FluidFieldParameters(7) * Kilogram / Second

  END SUBROUTINE ReadParameters

  
  SUBROUTINE ReadData( FILEIN, nLines, FluidFieldData )

    CHARACTER( LEN = * ), INTENT(in)   :: FILEIN
    INTEGER,  INTENT(inout)            :: nLines
    REAL(DP), INTENT(out), ALLOCATABLE :: FluidFieldData(:,:)
    INTEGER                            :: i

    ! --- Get number of lines in data file ---
    nLines = 0
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    DO
      READ( 100, *, END = 10 )
      nLines = nLines + 1
    END DO
    10 CLOSE( 100 )

    ! --- Allocate and read in data ---
    ALLOCATE( FluidFieldData( 1:nLines, 4 ) )

    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    DO i = 1, nLines
       READ( 100, '(4ES23.16E2)' ) FluidFieldData(i,:)
    END DO
    CLOSE( 100 )

    !WRITE(*,*) nLines ! 20,000
    !STOP

    !WRITE(*,*) SHAPE(FluidFieldData)
    !WRITE(*,*) SIZE(FluidFieldData(:,1))
    !STOP

    ! --- Convert from physical-units to code-units ---
    FluidFieldData(:,1) = FluidFieldData(:,1) * Meter
    FluidFieldData(:,2) = FluidFieldData(:,2) * Kilogram / Meter**3
    FluidFieldData(:,3) = FluidFieldData(:,3) * Meter / Second
    FluidFieldData(:,4) = FluidFieldData(:,4) * Joule / Meter**3

  END SUBROUTINE ReadData

  
END MODULE DataFileReaderModule
