MODULE  MF_TwoMoment_DiscretizationModule_Streaming_Relativistic

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter,   &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---
  USE ProgramHeaderModule,      ONLY: &
    swX, nDOFX, nDOFZ, swE, nDOFE, iE_B0, iE_E0, iE_B1, iE_E1
  USE GeometryFieldsModule,     ONLY: &
    nGF
  USE GeometryFieldsModuleE,     ONLY: &
    nGE, uGE
  USE RadiationFieldsModule,            ONLY: &
    nCR
  USE FluidFieldsModule,            ONLY: &
    nCF
  USE TwoMoment_DiscretizationModule_Streaming_Relativistic, ONLY: &
    ComputeIncrement_TwoMoment_Explicit

  ! --- Local Modules ---
  USE MF_UtilitiesModule,                ONLY: &
    AMReX2thornado_Euler, &
    thornado2AMReX_Euler, &
    AMReX2thornado, &
    thornado2AMReX
  USE MyAmrModule,                       ONLY: &
    nLevels, &
    nSpecies, &
    nE
  USE MF_TwoMoment_BoundaryConditionsModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap, &
    MF_ApplyBoundaryConditions_TwoMoment

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_TwoMoment_ComputeIncrement_Explicit


CONTAINS


  SUBROUTINE MF_TwoMoment_ComputeIncrement_Explicit( GEOM, MF_uGF, MF_uCF, MF_uCR, MF_duCR, Verbose_Option )

    TYPE(amrex_geometry), INTENT(in)    :: GEOM   (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCR(0:nLevels-1)
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCR (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: duCR(:,:,:,:)

    REAL(amrex_real), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: C (:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U (:,:,:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: dU(:,:,:,:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), i, iZ2, iZ1


    LOGICAL :: Verbose

    TYPE(EdgeMap) :: Edge_Map

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option



    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---
      CALL MF_uCR(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_duCR(iLevel) % setval( 0.0_amrex_real )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF (iLevel) % DataPtr( MFI )
        uCF  => MF_uCF (iLevel) % DataPtr( MFI )
        uCR  => MF_uCR (iLevel) % DataPtr( MFI )
        duCR => MF_duCR(iLevel) % DataPtr( MFI )
        
        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX
 
        i=1          

        DO WHILE (i<=4)
          
          IF (i==1) THEN
          
            iZ_B0(i)=iE_B0 
            iZ_E0(i)=iE_E0
            iZ_B1(i)=iE_B1
            iZ_E1(i)=iE_E1

          ELSE 

            iZ_B0(i)=iX_B0(i-1) 
            iZ_E0(i)=iX_E0(i-1)
            iZ_B1(i)=iX_B1(i-1)
            iZ_E1(i)=iX_E1(i-1)

          END IF
          i = i + 1 
        END DO


        ALLOCATE( G (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( C (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nCF) )

        ALLOCATE( U (1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                             iZ_B1(3):iZ_E1(3), &
                             iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )

        ALLOCATE( dU (1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                             iZ_B1(3):iZ_E1(3), &
                             iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )


        CALL AMReX2thornado_Euler &
               ( nGF, iX_B1, iX_E1, &
                 uGF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nGF), &
                 G(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nGF) )

        CALL AMReX2thornado_Euler &
               ( nCF, iX_B1, iX_E1, &
                 uCF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nCF), &
                 C(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nCF) )



        CALL AMReX2thornado &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iX_B1, iX_E1,                    &
                 uCR(      iZ_B1(2):iZ_E1(2), &
                           iZ_B1(3):iZ_E1(3), &
                           iZ_B1(4):iZ_E1(4),1:nDOFZ*nCR*nSpecies*nE), &
                 U(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                           iZ_B1(2):iZ_E1(2), &
                           iZ_B1(3):iZ_E1(3), &
                           iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )
        

        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )        

        CALL MF_ApplyBoundaryConditions_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, Edge_Map )

        CALL ComputeIncrement_TwoMoment_Explicit &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, G, C, U, dU, Verbose_Option = Verbose, &
               SuppressBC_Option = .TRUE.  )

        CALL thornado2AMReX &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, iX_B0, iX_E0, &
                 duCR(     iZ_B0(2):iZ_E0(2), &
                           iZ_B0(3):iZ_E0(3), &
                           iZ_B0(4):iZ_E0(4),1:nDOFZ*nCR*nSpecies*nE), &
                 dU(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                            iZ_B0(2):iZ_E0(2), &
                            iZ_B0(3):iZ_E0(3), &
                            iZ_B0(4):iZ_E0(4),1:nCR,1:nSpecies) )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO



  END SUBROUTINE MF_TwoMoment_ComputeIncrement_Explicit


END MODULE  MF_TwoMoment_DiscretizationModule_Streaming_Relativistic
