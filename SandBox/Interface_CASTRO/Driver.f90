PROGRAM Driver

  USE KindModule, ONLY: &
    DP
  USE ThornadoInitializationModule, ONLY: &
    InitThornado, &       ! --- To be called once
    InitThornado_Patch, & ! --- To be called once per patch
    FreeThornado_Patch    ! --- To be called once per parch

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: i
  INTEGER  :: mpierr
  REAL(DP) :: wTime

  CALL MPI_INIT( mpierr )

  wTime = MPI_WTIME( )

  CALL InitThornado( nDimsX = 3, nDimsE = 1 )

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A,A32,ES8.2E2)') '', 'InitThornado: ', wTime
  WRITE(*,*)

  wTime = MPI_WTIME( )

  DO i = 1, 2

    CALL InitThornado_Patch &
           ( nX  = [ 16, 16, 16 ], &
             swX = [ 2, 2, 2 ], &
             xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
             xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
             nE  = 12, &
             swE = 0, &
             eL  = 0.0_DP, &
             eR  = 1.0_DP, &
             nSpecies = 1 )

    CALL FreeThornado_Patch

  END DO

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A,A32,ES8.2E2)') '', 'InitThornado_Patch: ', wTime
  WRITE(*,*)

  CALL MPI_FINALIZE( mpierr )

END PROGRAM Driver
