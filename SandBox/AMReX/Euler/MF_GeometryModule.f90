MODULE MF_GeometryModule

  ! --- AMReX Modules ---

  USE amrex_base_module

  ! --- thornado Modules ---

  USE KindModule,                ONLY: &
    DP
  USE ProgramHeaderModule,       ONLY: &
    nDOFX
  USE GeometryFieldsModule,      ONLY: &
    nGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeGeometryX

CONTAINS


  SUBROUTINE MF_ComputeGeometryX( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF

    INTEGER :: iX1, iX2, iX3
    INTEGER :: lo(4), hi(4)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)

    IF( amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'MF_ComputeGeometryX'
    END IF

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = .TRUE. )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )

      BX = MFI % tilebox()
      lo = LBOUND( uGF )
      hi = UBOUND( uGF )

      ALLOCATE( G(1:nDOFX,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nGF) )

      CALL ComputeGeometryX &
             ( BX % lo(1:3), BX % hi(1:3), lo(1:3), hi(1:3), G )

      DO iX3 = lo(3), hi(3)
      DO iX2 = lo(2), hi(2)
      DO iX1 = lo(1), hi(1)

        uGF(iX1,iX2,iX3,:) &
          = RESHAPE( G(1:nDOFX,iX1,iX2,iX3,1:nGF), [hi(4)-lo(4)+1] )

      END DO
      END DO
      END DO

      DEALLOCATE( G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE MF_ComputeGeometryX


END MODULE MF_GeometryModule
