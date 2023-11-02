MODULE MemoryProfilingModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WriteMemoryUsage

CONTAINS


  SUBROUTINE WriteMemoryUsage( ounit, Label, iCycle, Time, WriteMemory_Option )

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
    REAL(DP)        , INTENT(IN) :: Time   ! Current time of simulation
    LOGICAL         , INTENT(IN), OPTIONAL :: WriteMemory_Option

    !-----------------------------------------------------------------------
    !         Local Variables
    !-----------------------------------------------------------------------

    CHARACTER(LEN=80) :: Line     ! Line from file
    INTEGER(KIND=8)   :: hwm, rss ! memory sizes
    INTEGER           :: istat    ! open file flag
    INTEGER           :: iunit    ! Logical unit number for input
    LOGICAL           :: WriteMemory

    WriteMemory = .FALSE.
    IF( PRESENT( WriteMemory_Option ) ) &
      WriteMemory = WriteMemory_Option

    IF( WriteMemory )THEN

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

      END DO

      CLOSE( iunit )

      !-----------------------------------------------------------------------
      !         Ouput current memory size and high water mark
      !-----------------------------------------------------------------------

      WRITE( ounit, '(I8.8,1x,ES23.16E3,1x,I8.8,1x,I8.8)' ) &
        iCycle, Time, rss, hwm

    END IF ! WriteMemory

  END SUBROUTINE WriteMemoryUsage

END MODULE MemoryProfilingModule
