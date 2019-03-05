MODULE  MF_Euler_dgDiscretizationModule

  ! --- AMReX Modules ---
  USE amrex_base_module, ONLY: &
    amrex_multifab, &
    amrex_box,      &
    amrex_geometry, &
    amrex_mfiter,   &
    amrex_mfiter_build
  USE amrex_fort_module, ONLY: &
    amrex_real

  ! --- thornado Modules ---
  USE ProgramHeaderModule,      ONLY: &
    swX, nDOFX
  USE FluidFieldsModule,        ONLY: &
    nCF
  USE GeometryFieldsModule,     ONLY: &
    nGF
  USE Euler_dgDiscretizationModule, ONLY: &
    Euler_ComputeIncrement_DG_Explicit

  ! --- Local Modules ---
  USE MF_UtilitiesModule, ONLY: &
    AMReX2thornado, &
    thornado2AMReX
  USE MyAmrModule,        ONLY: &
    nLevels

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_Euler_ComputeIncrement


CONTAINS


  SUBROUTINE MF_Euler_ComputeIncrement( GEOM, MF_uGF, MF_uCF, MF_duCF )
 
    TYPE(amrex_geometry), INTENT(in)    :: GEOM   (0:nLevels)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF (0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF (0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCF(0:nLevels)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: duCF(:,:,:,:)

    REAL(amrex_real), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U (:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: dU(:,:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    DO iLevel = 0, nLevels

      ! --- Apply boundary conditions to geometry and conserved fluid ---
      !     Do we need to apply boundary conditions to geometry?
      !     If not, when is it applied? ---
      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )
      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_duCF(iLevel) % setval( 0.0_amrex_real )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF (iLevel) % DataPtr( MFI )
        uCF  => MF_uCF (iLevel) % DataPtr( MFI )
        duCF => MF_duCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        ALLOCATE( G (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nGF) )
        ALLOCATE( U (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nCF) )
        ALLOCATE( dU(1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nCF) )

        CALL AMReX2thornado &
               ( nGF, iX_B1, iX_E1, &
                 uGF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nGF), &
                 G(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nGF) )

        CALL AMReX2thornado &
               ( nCF, iX_B1, iX_E1, &
                 uCF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nCF), &
                 U(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nCF) )


        CALL Euler_ComputeIncrement_DG_Explicit &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                 G (1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF), &
                 U (1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF), &
                 dU(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nCF) )


        CALL thornado2AMReX &
               ( nCF, iX_B0, iX_E0, &
                 duCF(      iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nDOFX*nCF), &
                 dU(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nCF) )

        DEALLOCATE( G  )
        DEALLOCATE( U  )
        DEALLOCATE( dU )

      END DO

    END DO

  END SUBROUTINE MF_Euler_ComputeIncrement


END MODULE MF_Euler_dgDiscretizationModule
