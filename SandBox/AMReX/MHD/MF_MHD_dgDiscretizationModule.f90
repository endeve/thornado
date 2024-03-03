MODULE  MF_MHD_dgDiscretizationModule

  ! --- AMReX Modules ---

  USE amrex_box_module,                   ONLY: &
    amrex_box
  USE amrex_geometry_module,              ONLY: &
    amrex_geometry
  USE amrex_multifab_module,              ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,              ONLY: &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE ProgramHeaderModule,                ONLY: &
    swX, &
    nDOFX
  USE MagnetofluidFieldsModule,                  ONLY: &
    nCM, &
    nDM
  USE GeometryFieldsModule,               ONLY: &
    nGF
  USE MHD_DiscretizationModule_Relativistic,       ONLY: &
    ComputeIncrement_MHD_DG_Explicit

  ! --- Local Modules ---

  USE MF_KindModule,                      ONLY: &
    DP, &
    Zero
  USE MF_UtilitiesModule,                 ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    WriteNodalDataToFile
  USE InputParsingModule,                 ONLY: &
    nLevels, &
    UseTiling, &
    DEBUG, &
    EvolveOnlyMagnetic, &
    UseDivergenceCleaning, &
    DampingParameter, &
    UsePowellSource, &
    WriteNodalData, &
    StepNo
  USE MF_MHD_BoundaryConditionsModule,  ONLY: &
    EdgeMap,          &
    ConstructEdgeMap, &
    MF_ApplyBoundaryConditions_MHD

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeIncrement_MHD


CONTAINS


  SUBROUTINE MF_ComputeIncrement_MHD( GEOM, MF_uGF, MF_uCM, MF_uDM, MF_duCM )

    TYPE(amrex_geometry), INTENT(in)    :: GEOM   (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCM (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uDM (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCM(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDM (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: duCM(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU(:,:,:,:,:)

    INTEGER :: iLevel, iCM
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    TYPE(EdgeMap) :: Edge_Map

    CHARACTER(LEN=19) :: NodalFileBaseName

!    !DO iLevel = 0, nLevels-1
!
!      ! --- Apply boundary conditions to interior domains ---
!
!      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )
!
!      CALL MF_uCM(iLevel) % Fill_Boundary( GEOM(iLevel) )
!
!      CALL MF_uDM(iLevel) % Fill_Boundary( GEOM(iLevel) )
!
!      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )
!
!      DO WHILE( MFI % next() )
!
!        uGF  => MF_uGF (iLevel) % DataPtr( MFI )
!        uCM  => MF_uCM (iLevel) % DataPtr( MFI )
!        uDM  => MF_uDM (iLevel) % DataPtr( MFI )
!
!        iLo_MF = LBOUND( uGF )
!
!        BX = MFI % tilebox()
!
!        iX_B0 = BX % lo
!        iX_E0 = BX % hi
!        iX_B1 = BX % lo - swX
!        iX_E1 = BX % hi + swX
!
!        ALLOCATE( G (1:nDOFX,iX_B1(1):iX_E1(1), &
!                             iX_B1(2):iX_E1(2), &
!                             iX_B1(3):iX_E1(3),1:nGF) )
!
!        ALLOCATE( U (1:nDOFX,iX_B1(1):iX_E1(1), &
!                             iX_B1(2):iX_E1(2), &
!                             iX_B1(3):iX_E1(3),1:nCM) )
!
!        ALLOCATE( D (1:nDOFX,iX_B1(1):iX_E1(1), &
!                             iX_B1(2):iX_E1(2), &
!                             iX_B1(3):iX_E1(3),1:nDM) )
!
!        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )
!
!        CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCM, U )
!
!        CALL amrex2thornado_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )
!
!        ! --- Apply boundary conditions to physical boundaries ---
!
!        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )
!
!        IF( DEBUG ) WRITE(*,'(A)') '    CALL MF_ApplyBoundaryConditions_MHD'
!
!        CALL MF_ApplyBoundaryConditions_MHD &
!               ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )
!
!        CALL thornado2amrex_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )
!
!        DEALLOCATE( D  )
!
!        DEALLOCATE( U  )
!
!        DEALLOCATE( G  )
!
!      END DO
!
!      CALL amrex_mfiter_destroy( MFI )
!
!    END DO

    DO iLevel = 0, nLevels-1

      ! --- Maybe don't need to apply boudnary conditions since
      !     they're applied in the shock detector ---

      ! --- Apply boundary conditions to interior domains ---

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCM(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uDM(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_duCM(iLevel) % setval( Zero )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF (iLevel) % DataPtr( MFI )
        uCM  => MF_uCM (iLevel) % DataPtr( MFI )
        uDM  => MF_uDM (iLevel) % DataPtr( MFI )
        duCM => MF_duCM(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

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
                             iX_B1(3):iX_E1(3),1:nCM) )

        ALLOCATE( D (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nDM) )

        ALLOCATE( dU(1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nCM) )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCM, U )

        CALL amrex2thornado_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )

        ! --- Apply boundary conditions to physical boundaries ---

        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )

        IF( DEBUG ) WRITE(*,'(A)') '    CALL MF_ApplyBoundaryConditions_MHD'

        CALL MF_ApplyBoundaryConditions_MHD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

        IF( DEBUG ) WRITE(*,'(A)') '    CALL ComputeIncrement_MHD_DG_Explicit'

        CALL ComputeIncrement_MHD_DG_Explicit &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, &
                 SuppressBC_Option = .TRUE., &
                 EvolveOnlyMagnetic_Option = EvolveOnlyMagnetic, &
                 UseDivergenceCleaning_Option = UseDivergenceCleaning, &
                 DampingParameter_Option = DampingParameter, &
                 UsePowellSource_Option = UsePowellSource )

        CALL thornado2amrex_X &
               ( nCM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, duCM, dU )

        DEALLOCATE( dU )

        DEALLOCATE( D  )

        DEALLOCATE( U  )

        DEALLOCATE( G  )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ComputeIncrement_MHD


END MODULE MF_MHD_dgDiscretizationModule
