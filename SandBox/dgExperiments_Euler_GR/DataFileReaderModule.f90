MODULE DataFileReader

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ReadParameters
  PUBLIC :: ReadData

CONTAINS

  SUBROUTINE ReadParameters( FILEIN , FluidFieldParameters )

    CHARACTER( LEN = * ), INTENT(in)   :: FILEIN
    REAL(DP), INTENT(out), ALLOCATABLE :: FluidFieldParameters(:)
    INTEGER                            :: i, nParams

    ! --- Get number of parameters ---
    nParams = 0
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * )
    DO
      READ( 100, *, END = 10 )
      nParams = nParams + 1
    END DO
    10 CLOSE( 100 )

    ! --- Allocate and read in parameters ---
    ALLOCATE( FluidFieldParameters(nParams) )
    
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * )
    DO i = 1, nParams
       READ( 100, '(E17.10)' ) FluidFieldParameters(i)
    END DO
    CLOSE( 100 )

  END SUBROUTINE ReadParameters

  
  SUBROUTINE ReadData( FILEIN , FluidFieldData )

    CHARACTER( LEN = * ), INTENT(in)   :: FILEIN
    REAL(DP), INTENT(out), ALLOCATABLE :: FluidFieldData(:,:)
    INTEGER                            :: i, nLines

    ! --- Get number of lines in data file ---
    nLines = 0
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * )
    DO
      READ( 100, *, END = 10 )
      nLines = nLines + 1
    END DO
    10 CLOSE( 100 )

    ! --- Allocate and read in data ---
    ALLOCATE( FluidFieldData( nLines, 4 ) )

    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * )
    DO i = 1 , nLines
      READ( 100, '(4E17.10)' ) FluidFieldData( i , : )
    END DO
    CLOSE( 100 )

    WRITE(*,*) FluidFieldData( 1 , : )
  END SUBROUTINE ReadData

  
END MODULE DataFileReader
