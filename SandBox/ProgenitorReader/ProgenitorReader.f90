PROGRAM ProgenitorReader

  USE KindModule, ONLY: &
    DP
  USE UtilitiesModule, ONLY: &
    WriteVector

  IMPLICIT NONE

  CHARACTER(LEN=12) :: LabelR = '   R [km]'
  CHARACTER(LEN=12) :: LabelD = '   D [g/ccm]'
  CHARACTER(LEN=12) :: LabelT = '   T [mev]'
  CHARACTER(LEN=12) :: LabelY = '   Y'

  CALL ReadProgenitor &
         ( '../Progenitors/profile_data_000100ms.d', &
           Suffix_Option = '100ms' )

  CALL ReadProgenitor_VX &
         ( '../Progenitors/VX_G15/datafileVX-G15h+003.txt', &
           Suffix_Option = '003ms' )

  CALL ReadProgenitor_VX &
         ( '../Progenitors/VX_G15/datafileVX-G15h+010.txt', &
           Suffix_Option = '010ms' )

  CALL ReadProgenitor_VX &
         ( '../Progenitors/VX_G15/datafileVX-G15h+050.txt', &
           Suffix_Option = '050ms' )

  CALL ReadProgenitor_VX &
         ( '../Progenitors/VX_G15/datafileVX-G15h+100.txt', &
           Suffix_Option = '100ms' )

  CALL ReadProgenitor_VX &
         ( '../Progenitors/VX_G15/datafileVX-G15h+150.txt', &
           Suffix_Option = '150ms' )

  CALL ReadProgenitor_VX &
         ( '../Progenitors/VX_G15/datafileVX-G15h+250.txt', &
           Suffix_Option = '250ms' )

CONTAINS


  SUBROUTINE ReadProgenitor( FileName, Suffix_Option )

    CHARACTER(*), INTENT(in)           :: FileName
    CHARACTER(*), INTENT(in), OPTIONAL :: Suffix_Option

    CHARACTER(LEN=06)     :: Suffix
    CHARACTER(LEN=04)     :: Format0 = "(I5)"
    CHARACTER(LEN=27)     :: Format1 = "(I5,2ES12.3,ES14.6,6ES12.3)"
    CHARACTER(LEN=06)     :: Format2 = '(4A12)'
    CHARACTER(LEN=09)     :: Format3 = '(4ES12.3)'
    INTEGER               :: nPoints, iPoint
    INTEGER               :: LineNumber, Status
    REAL(DP), ALLOCATABLE :: Data(:,:)

    IF( PRESENT( Suffix_Option ) )THEN
      Suffix = TRIM( Suffix_Option )
    ELSE
      Suffix = ''
    END IF

    ! --- Count Lines ---

    OPEN( 1, FILE = TRIM( FileName ), FORM = "formatted", ACTION = 'read' )
    READ( 1, * )
    DO
      READ( 1, Format0, IOSTAT = Status ) nPoints
      IF( Status .NE. 0 ) EXIT
    END DO
    CLOSE( 1, STATUS = 'keep' )

    ! --- Read Data ---

    ALLOCATE( Data(nPoints,9) )

    OPEN( 1, FILE = TRIM( FileName ), FORM = "formatted", ACTION = 'read' )
    READ( 1, * )
    DO iPoint = 1, nPoints
      READ( 1, Format1, IOSTAT=Status) LineNumber, Data(iPoint,:)
    END DO
    CLOSE( 1, STATUS = 'keep' )

    CALL WriteVector( nPoints, Data(:,2), 'R.dat' )
    CALL WriteVector( nPoints, Data(:,5), 'D.dat' )
    CALL WriteVector( nPoints, Data(:,6), 'T.dat' )
    CALL WriteVector( nPoints, Data(:,9), 'Y.dat' )

    OPEN( 1, FILE = 'input_thornado_Chimera_' // TRIM( Suffix ) // '.dat', &
          FORM = 'formatted', ACTION = 'write' )

    WRITE( 1, Format2 ) LabelR, LabelD, LabelT, LabelY
    DO iPoint = 1, nPoints
      WRITE( 1, Format3 ) &
        Data(iPoint,2), Data(iPoint,5), Data(iPoint,6), Data(iPoint,9)
    END DO

    CLOSE( 1, STATUS = 'keep' )

    DEALLOCATE( Data )

  END SUBROUTINE ReadProgenitor


  SUBROUTINE ReadProgenitor_VX( FileName, Suffix_Option )

    CHARACTER(*), INTENT(in)           :: FileName
    CHARACTER(*), INTENT(in), OPTIONAL :: Suffix_Option

    CHARACTER(LEN=06)     :: Suffix
    CHARACTER(LEN=04)     :: Format0 = '(I3)'
    CHARACTER(LEN=12)     :: Format1 = '(I3,18E12.3)'
    CHARACTER(LEN=06)     :: Format2 = '(4A12)'
    CHARACTER(LEN=09)     :: Format3 = '(4ES12.3)'
    INTEGER               :: i, Status, LineNumber
    INTEGER               :: nPoints, iPoint
    REAL(DP), ALLOCATABLE :: Data(:,:)

    IF( PRESENT( Suffix_Option ) )THEN
      Suffix = TRIM( Suffix_Option )
    ELSE
      Suffix = ''
    END IF

    OPEN( 1, FILE = TRIM( FileName ), FORM = 'formatted', ACTION = 'read' )

    ! --- Skip 36 Lines

    DO i = 1, 36
      READ( 1, * )
    END DO

    ! --- Count Lines ---

    DO
      READ( 1, Format0, IOSTAT = Status ) nPoints
      IF( Status .NE. 0 ) EXIT
    END DO

    CLOSE( 1, STATUS = 'keep' )

    nPoints = nPoints + 1 ! Add one since counter starts at 0

    ! --- Read Data ---

    ALLOCATE( Data(nPoints,18) )

    OPEN( 1, FILE = TRIM( FileName ), FORM = 'formatted', ACTION = 'read' )

    ! --- Skip 36 Lines ---

    DO i = 1, 36
      READ( 1, * )
    END DO

    DO iPoint = 1, nPoints
      READ( 1, Format1, IOSTAT=Status) LineNumber, Data(iPoint,:)
    END DO

    CLOSE( 1, STATUS = 'keep' )

    CALL WriteVector( nPoints, Data(:,03), 'R_' // TRIM( Suffix ) // '.dat' )
    CALL WriteVector( nPoints, Data(:,05), 'D_' // TRIM( Suffix ) // '.dat' )
    CALL WriteVector( nPoints, Data(:,07), 'T_' // TRIM( Suffix ) // '.dat' )
    CALL WriteVector( nPoints, Data(:,10), 'Y_' // TRIM( Suffix ) // '.dat' )

    OPEN( 1, FILE = 'input_thornado_VX_' // TRIM( Suffix ) // '.dat', &
          FORM = 'formatted', ACTION = 'write' )

    WRITE( 1, Format2 ) LabelR, LabelD, LabelT, LabelY
    DO iPoint = 1, nPoints
      WRITE( 1, Format3 ) &
        Data(iPoint,03)/1.d5, Data(iPoint,05), &
        Data(iPoint,07)     , Data(iPoint,10)
    END DO

    CLOSE( 1, STATUS = 'keep' )

    DEALLOCATE( Data )

  END SUBROUTINE ReadProgenitor_VX


END PROGRAM ProgenitorReader
