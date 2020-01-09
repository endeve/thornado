MODULE MF_TimeSteppingModule_IMEX

  ! --- AMReX Modules ---
  USE amrex_fort_module,      ONLY: &
    amrex_real
  USE amrex_box_module,       ONLY: &
    amrex_box
  USE amrex_geometry_module,  ONLY: &
    amrex_geometry
  USE amrex_multifab_module,  ONLY: &
    amrex_multifab, &
    amrex_multifab_build, amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, amrex_mfiter_destroy
  USE amrex_boxarray_module,  ONLY: &
    amrex_boxarray, &
    amrex_boxarray_build, amrex_boxarray_destroy
  USE amrex_distromap_module, ONLY: &
    amrex_distromap, &
    amrex_distromap_build, amrex_distromap_destroy

  ! --- thornado Modules ---
  USE ProgramHeaderModule,              ONLY: &
    DescribeProgramHeaderX, &
    nDOFX,                  &
    nNodesX,                &
    nDOFE,                  &
    nNodesE,                &
    nDOFZ,                  &
    iZ_B0,                  &
    iZ_E0,                  &
    iZ_B1,                  &
    iZ_E1
  USE RadiationFieldsModule,            ONLY: &
    nPR,                    &
    iPR_D,                  &
    iPR_I1
  USE MF_UtilitiesModule,                ONLY: &
    AMReX2thornado, &
    thornado2AMReX

  ! --- Local modules ---
  USE MyAmrModule,                      ONLY: &
    t_end,                     &
    t,                         &
    dt,                        &
    t_wrt,                     &
    dt_wrt,                    &
    t_chk,                     &
    dt_chk,                    &
    CFL,                       &
    nNodes,                    &
    nStages,                   &
    nX,                        &
    nE,                        &
    swX,                       &
    swE,                       &
    bcX,                       &
    xL,                        &
    xR,                        &
    ProgramName,               &
    CoordSys,                  &
    StepNo,                    &
    nLevels,                   &
    iRestart,                  &
    MaxGridSizeX,              &
    BA,                        &
    DM,                        &
    GEOM,                      &
    StepNo,                    &
    nSpecies,                  &
    MyAmrInit


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_UpdateField_IMEX


CONTAINS

  SUBROUTINE MF_UpdateField_IMEX &
    ( t, dt, MF_uPR )

    REAL(amrex_real),     INTENT(in)    :: t(0:nLevels-1), dt(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uPR(0:nLevels-1)
    REAL(amrex_real), ALLOCATABLE :: P(:,:,:,:,:,:,:) 
    REAL(amrex_real), CONTIGUOUS, POINTER :: uPR (:,:,:,:)

    TYPE(amrex_mfiter)                    :: MFI

    INTEGER :: iLevel, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ
    DO iLevel = 0, nLevels-1
      
      CALL amrex_mfiter_build( MFI, MF_uPR(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )
        
        uPR  => MF_uPR (iLevel) % DataPtr( MFI )

        ALLOCATE( P(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                            iZ_B0(2):iZ_E0(2), &
                            iZ_B0(3):iZ_E0(3), & 
                            iZ_B0(4):iZ_E0(4), &
                            1:nPR, 1:nSpecies) )
        !Issue here with using the entire domain rather than just the part that is for the single worker still works with 1 worker 
        CALL AMReX2thornado( nPR, nSpecies, nE, iZ_B0(1:4), iZ_E0(1:4), & 
                              uPR(iZ_B0(2):iZ_E0(2), &
                                   iZ_B0(3):iZ_E0(3), & 
                                   iZ_B0(4):iZ_E0(4), &
                                   1:nDOFZ*nPR*nSpecies*nE),&
                              P(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                                   iZ_B0(2):iZ_E0(2), &
                                   iZ_B0(3):iZ_E0(3), & 
                                   iZ_B0(4):iZ_E0(4), &
                                   1:nPR, 1:nSpecies) )


         DO iS  = 1, nSpecies
         DO iZ4 = iZ_B0(4), iZ_E0(4)
         DO iZ3 = iZ_B0(3), iZ_E0(3)
         DO iZ2 = iZ_B0(2), iZ_E0(2)
         DO iZ1 = iZ_B0(1), iZ_E0(1)
         
           DO iNodeZ = 1, nDOFZ
            
             P(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
               = P(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) + 0.1_amrex_real

             P(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1 ,iS) &
               = P(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) + 0.1_amrex_real
           
           END DO

         END DO
         END DO
         END DO
         END DO
         END DO


         CALL thornado2AMReX( nPR, nSpecies, nE, iZ_B0(1:4), iZ_E0(1:4), & 
                              uPR(iZ_B0(2):iZ_E0(2), &
                                   iZ_B0(3):iZ_E0(3), & 
                                   iZ_B0(4):iZ_E0(4), &
                                   1:nDOFZ*nPR*nSpecies*nE),&
                              P(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                                   iZ_B0(2):iZ_E0(2), &
                                   iZ_B0(3):iZ_E0(3), & 
                                   iZ_B0(4):iZ_E0(4), &
                                   1:nPR, 1:nSpecies) )

      END DO

        CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_UpdateField_IMEX

END MODULE MF_TimeSteppingModule_IMEX
