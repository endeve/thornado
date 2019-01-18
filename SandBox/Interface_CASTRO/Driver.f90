PROGRAM Driver

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer, &
    Millisecond, &
    MeV, &
    Erg
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
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

  CALL InitThornado &
         ( nDimsX = 3, nE = 10, swE = 0, eL_in = 0.0d0, eR_in = 1.0d2, &
           zoomE = 1.0_DP, nSpecies_in = 1 )

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A,A32,ES8.2E2)') '', 'InitThornado: ', wTime
  WRITE(*,*)

  wTime = MPI_WTIME( )

  DO i = 1, 2

    CALL InitThornado_Patch &
           ( nX    = [ 12, 12, 12 ], &
             swX   = [ 2, 2, 2 ], &
             xL    = [ 00.0_DP, 00.0_DP, 00.0_DP ] * Kilometer, &
             xR    = [ 16.0_DP, 16.0_DP, 16.0_DP ] * Kilometer )

    dt = 1.0d-4 * Millisecond

!    uCR(:,:,:,:,:,iCR_N, :) = 0.9_DP
!    uCR(:,:,:,:,:,iCR_G1,:) = 0.0_DP
!    uCR(:,:,:,:,:,iCR_G2,:) = 0.0_DP
!    uCR(:,:,:,:,:,iCR_G3,:) = 0.0_DP

!    uCF(:,:,:,:,iCF_D)  = 1.0d14 * Gram / Centimeter**3
!    uCF(:,:,:,:,iCF_S1) = 0.0_DP
!    uCF(:,:,:,:,iCF_S2) = 0.0_DP
!    uCF(:,:,:,:,iCF_S3) = 0.0_DP
!    uCF(:,:,:,:,iCF_E)  = 4.5d33 * Erg / Centimeter**3
!    uCF(:,:,:,:,iCF_Ne) = 2.0d37 / Centimeter**3

!    CALL Update_IMEX_PDARS( dt, uCF, uCR )

    CALL FreeThornado_Patch

  END DO

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A,A32,ES8.2E2)') '', 'InitThornado_Patch: ', wTime
  WRITE(*,*)

  CALL MPI_FINALIZE( mpierr )

END PROGRAM Driver
