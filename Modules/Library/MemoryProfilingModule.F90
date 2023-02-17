MODULE MemoryProfilingModule

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WriteMemoryUsage

CONTAINS


  SUBROUTINE WriteMemoryUsage( ounit, Label, iCycle )

  !-----------------------------------------------------------------------
  !
  !    Author:       O.E.B. Messer, ORNL
  !                  W.R. Hix, ORNL
  !
  !    Date:         8/12/11
  !
  !    Purpose:
  !      To determine and output the current memory usage
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !         Input Variables
  !-----------------------------------------------------------------------

  INTEGER         , INTENT(IN) :: ounit  ! Logical unit number for output
  CHARACTER(LEN=*), INTENT(IN) :: Label  ! Location label
  INTEGER         , INTENT(IN) :: iCycle ! Cycle number of simulation

#ifdef MEMORY_DISPLAY

  !-----------------------------------------------------------------------
  !         Local Variables
  !-----------------------------------------------------------------------

  CHARACTER(LEN=80) :: Line     ! Line from file
  INTEGER(KIND=8)   :: hwm, rss ! memory sizes
  INTEGER           :: istat    ! open file flag
  INTEGER           :: iunit    ! Logical unit number for input

  !-----------------------------------------------------------------------
  !         Find current memory size and high water mark
  !-----------------------------------------------------------------------

  hwm = 0
  rss = 0

  OPEN( NEWUNIT = iunit, FILE = '/proc/self/status', &
        STATUS = 'old', IOSTAT = istat )

  DO WHILE( .TRUE. )

    READ( iunit, '(A)', IOSTAT = istat ) Line

    IF( istat .LT. 0 ) EXIT

    IF( Line(1:6) .EQ. 'VmHWM:' ) READ( Line(8:80), * ) hwm

    IF( Line(1:6) .EQ. 'VmRSS:' ) READ( Line(8:80), * ) rss

  ENDDO

  CLOSE( iunit )

  !-----------------------------------------------------------------------
  !         Ouput current memory size and high water mark
  !-----------------------------------------------------------------------

  WRITE( ounit, '(3x,3A,I8,2A,I8,A,I8)' ) &
    'Current Memory usage ', TRIM( Label ), ' =', rss, ' kb; ', &
    'Maximum Memory Usage =', hwm, ' kB; iCycle ', iCycle

#else

    RETURN

#endif

  END SUBROUTINE WriteMemoryUsage

END MODULE MemoryProfilingModule
