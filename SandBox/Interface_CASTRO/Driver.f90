PROGRAM Driver

  USE KindModule, ONLY: &
    DP
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE FluidFieldsModule, ONLY: &
    uCF
  USE RadiationFieldsModule, ONLY: &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE ThornadoInitializationModule, ONLY: &
    InitThornado, &       ! --- To be called once
    InitThornado_Patch, & ! --- To be called once per patch
    FreeThornado_Patch    ! --- To be called once per parch
  USE TimeSteppingModule_Castro, ONLY: &
    Update_IMEX_PDARS

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: i
  INTEGER  :: mpierr
  REAL(DP) :: wTime
  REAL(DP) :: dt

  CALL MPI_INIT( mpierr )

  wTime = MPI_WTIME( )

  CALL InitThornado( nDimsX = 3, nE = 16, nSpeciesIn = 1 )

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
             swE = 0, &
             eL  = 0.0_DP, &
             eR  = 1.0_DP )

    dt = 1.0d-4

    uCR(:,:,:,:,:,iCR_N, :) = 1.0_DP
    uCR(:,:,:,:,:,iCR_G1,:) = 0.0_DP
    uCR(:,:,:,:,:,iCR_G2,:) = 0.0_DP
    uCR(:,:,:,:,:,iCR_G3,:) = 0.0_DP
!!$
!!$    CALL Update_IMEX_PDARS( dt, uCF, uCR )
!!$
    CALL FreeThornado_Patch

  END DO

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A,A32,ES8.2E2)') '', 'InitThornado_Patch: ', wTime
  WRITE(*,*)

  CALL MPI_FINALIZE( mpierr )

END PROGRAM Driver
