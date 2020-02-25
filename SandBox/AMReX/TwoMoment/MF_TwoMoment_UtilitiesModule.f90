MODULE MF_TwoMoment_UtilitiesModule

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
    amrex_parallel_reduce_min, &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---
  USE ProgramHeaderModule,      ONLY: &
    swX, nDOFX, nDOFZ, swE, nDOFE, iE_B0, iE_E0, iE_B1, iE_E1
  USE RadiationFieldsModule,            ONLY: &
    nCR
  USE GeometryFieldsModule,     ONLY: &
    nGF

  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels, nSpecies, nE
  USE MF_UtilitiesModule,                ONLY: &
    AMReX2thornado_Euler, &
    AMReX2thornado

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeTimeStep_Fancy
  PUBLIC :: MF_ComputeTimeStep

CONTAINS




  SUBROUTINE MF_ComputeTimeStep_Fancy( MF_uGF, MF_uCR, CFL, TimeStepMin )

    TYPE(amrex_multifab),  INTENT(in)  :: MF_uGF(0:nLevels-1), MF_uCR(0:nLevels-1)
    REAL(amrex_real),     INTENT(in)  :: CFL
    REAL(amrex_real),     INTENT(out) :: TimeStepMin(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCR(:,:,:,:)


    REAL(amrex_real), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U(:,:,:,:,:,:,:)

    REAL(amrex_real) :: TimeStep(0:nLevels-1)
    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
 
    TimeStepMin = HUGE( 1.0e0_amrex_real )    

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )


        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0(1:3) = BX % lo(1:3)
        iX_E0(1:3) = BX % hi(1:3)
        iX_B1(1:3) = BX % lo(1:3) - swX(1:3)
        iX_E1(1:3) = BX % hi(1:3) + swX(1:3)


        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U (1:nDOFZ,iE_B1:iE_E1,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nCR,1:nSpecies) )

        CALL AMReX2thornado_Euler &
               ( nGF, iX_B1(1:3), iX_E1(1:3), &
                 uGF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nGF), &
                 G(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nGF) )

        CALL AMReX2thornado &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iX_B1, iX_E1,                    &
                 uCR(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFZ*nCR*nSpecies*nE), &
                 U(1:nDOFZ,iE_B0:iE_E0, &
                           iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nCR,1:nSpecies) )

!put compute timestep here

        TimeStep( iLevel ) = 0.01_amrex_real
 
        TimeStepMin( iLevel ) = MIN( TimeStepMin( iLevel ), TimeStep( iLevel ) )

        DEALLOCATE( U )
        DEALLOCATE( G )



      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

    END DO ! --- Loop over levels ---



  END SUBROUTINE MF_ComputeTimeStep_Fancy

  SUBROUTINE MF_ComputeTimeStep( nX, xR, xL, nNodes, CFL, TimeStepMin )

    INTEGER,              INTENT(in)  :: nX(:), nNodes
    REAL(amrex_real),     INTENT(in)  :: xR(:), xL(:), CFL
    REAL(amrex_real),     INTENT(out) :: TimeStepMin(0:nLevels-1)


    INTEGER          :: iLevel
    REAL(amrex_real) :: TimeStep(0:nLevels-1)


    TimeStepMin = HUGE( 1.0e0_amrex_real )    

    DO iLevel = 0, nLevels-1

      TimeStep( iLevel ) = CFL * MINVAL( (xR-xL) / DBLE(nX) ) &
                           / ( 2.0_amrex_real * DBLE(nNodes-1) + 1.0_amrex_real )

      TimeStepMin( iLevel ) = MIN( TimeStepMin( iLevel ), TimeStep( iLevel ) )

    END DO






  END SUBROUTINE MF_ComputeTimeStep


END MODULE MF_TwoMoment_UtilitiesModule
