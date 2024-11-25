MODULE ProgenitorModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Centimeter, &
    Gram, &
    Second, &
    Kelvin

  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: ProgenitorType1D
    INTEGER :: nPoints
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Radius
    REAL(DP), DIMENSION(:), ALLOCATABLE :: MassDensity
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Temperature
    REAL(DP), DIMENSION(:), ALLOCATABLE :: ElectronFraction
  END type ProgenitorType1D
 
  PUBLIC :: ReadProgenitor1D 

CONTAINS

  SUBROUTINE ReadProgenitor1D( FileName, P1D )

    CHARACTER(LEN=*),       INTENT(in)  :: FileName
    TYPE(ProgenitorType1D), INTENT(out) :: P1D 

    CHARACTER(LEN=8)      :: Format1 = "(A2,I10)"
    CHARACTER(LEN=6)      :: Format2 = "(4A13)"
    CHARACTER(LEN=9)      :: Format3 = "(4ES13.4)"

    CHARACTER(LEN=10)                   :: a
    INTEGER                             :: ii, datasize
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: database

    WRITE(*,*)
    WRITE(*,'(A4,A24,A1,A)') &
      '', 'Reading Progenitor File:', '', TRIM( FileName )
    WRITE(*,*)

    OPEN(1, FILE = FileName, FORM = "formatted", &
            ACTION = 'read')
    READ( 1, Format1 ) a, datasize
    ALLOCATE( database( datasize, 4) )
    READ( 1, Format2 )
    DO ii = 1, datasize
      READ( 1, Format3 ) database(ii,1:4)
    END DO
    CLOSE( 1, STATUS = 'keep')

    P1D % Radius           = database(:,1) * Centimeter
    P1D % MassDensity      = database(:,2) * Gram / Centimeter**3
    P1D % Temperature      = database(:,3) * Kelvin
    P1D % ElectronFraction = database(:,4)

  END SUBROUTINE ReadProgenitor1D

END MODULE ProgenitorModule
