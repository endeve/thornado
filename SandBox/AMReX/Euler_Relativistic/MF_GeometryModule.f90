MODULE MF_GeometryModule

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---
  USE ProgramHeaderModule,                       ONLY: &
    nDOFX, swX
  USE GeometryFieldsModule,                      ONLY: &
    nGF, iGF_Phi_N
  USE GeometryComputationModule,                 ONLY: &
    ComputeGeometryX
  USE GravitySolutionModule_Newtonian_PointMass, ONLY: &
    ComputeGravitationalPotential

  ! --- Local Modules ---
  USE MyAmrModule,               ONLY: &
    nLevels

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeGeometryX
  PUBLIC :: MF_ComputeGravitationalPotential


CONTAINS


  SUBROUTINE MF_ComputeGeometryX( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels)

    INTEGER :: iX1, iX2, iX3, iLevel
    INTEGER :: lo(4), hi(4)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: G(:,:,:,:,:)

!!$    IF( amrex_parallel_ioprocessor() )THEN
!!$      WRITE(*,*)
!!$      WRITE(*,'(A4,A)') '', 'MF_ComputeGeometryX'
!!$    END IF

    DO iLevel = 0, nLevels

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

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

    END DO

  END SUBROUTINE MF_ComputeGeometryX


  SUBROUTINE MF_ComputeGravitationalPotential( MF_uGF, Mass )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels)
    REAL(amrex_real),     INTENT(in)    :: Mass

    INTEGER :: iX1, iX2, iX3, iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), ALLOCATABLE         :: G(:,:,:,:,:)

    DO iLevel = 0, nLevels

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()
        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        CALL ComputeGravitationalPotential &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                 G(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nGF), Mass )

        DO iX3 = iX_B1(3), iX_E1(3)
        DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)

          uGF(iX1,iX2,iX3,(iGF_Phi_N-1)*nDOFX+1:iGF_Phi_N*nDOFX) &
            = RESHAPE( G(1:nDOFX,iX1,iX2,iX3,iGF_Phi_N), [nDOFX] )

        END DO
        END DO
        END DO

        DEALLOCATE( G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ComputeGravitationalPotential


END MODULE MF_GeometryModule
