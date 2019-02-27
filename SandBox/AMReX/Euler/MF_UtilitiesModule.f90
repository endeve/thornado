MODULE MF_UtilitiesModule

  ! --- AMReX Modules ---
  USE amrex_base_module, ONLY: &
    amrex_init, &
    amrex_finalize, &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_box, &
    amrex_boxarray, &
    amrex_boxarray_build, &
    amrex_boxarray_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_distromap, &
    amrex_distromap_build, &
    amrex_distromap_destroy
  USE amrex_geometry_module, ONLY: &
    amrex_geometry, &
    amrex_geometry_build, &
    amrex_geometry_destroy
  USE amrex_amr_module, ONLY: &
    amrex_amrcore_finalize
  USE amrex_fort_module, ONLY: &
    amrex_real

  ! --- thornado Modules ---
  USE ProgramHeaderModule,     ONLY: &
    InitializeProgramHeader
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE FluidFieldsModule,       ONLY: &
    nCF, nPF, nAF
  USE GeometryFieldsModule,    ONLY: &
    nGF

  USE MyAmrModule
  USE MyAmrDataModule
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_Plotfile, &
    ReadCheckpointFile


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: AMReX2thornado
  PUBLIC :: thornado2AMReX
  PUBLIC :: LinComb
  PUBLIC :: ShowVariableFromMultiFab
  PUBLIC :: MakeMF_Diff


CONTAINS


  SUBROUTINE AMReX2thornado( nVars, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,          INTENT(in)  :: nVars
    INTEGER,          INTENT(in)  :: iX_B(3), iX_E(3)
    REAL(amrex_real), INTENT(in)  :: &
      Data_amrex   (   iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(amrex_real), INTENT(out) :: &
      Data_thornado(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)

    INTEGER :: iX1, iX2, iX3, iVar

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      DO iVar = 1, nVars
        Data_thornado(1:nDOFX,iX1,iX2,iX3,iVar) &
          = Data_amrex(iX1,iX2,iX3,nDOFX*(iVar-1)+1:nDOFX*iVar)
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE AMReX2thornado


  SUBROUTINE thornado2AMReX( nVars, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,          INTENT(in)  :: nVars
    INTEGER,          INTENT(in)  :: iX_B(3), iX_E(3)
    REAL(amrex_real), INTENT(out) :: &
      Data_amrex   (   iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(amrex_real), INTENT(in)  :: &
      Data_thornado(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)

    INTEGER :: iX1, iX2, iX3, iVar

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      DO iVar = 1, nVars
        Data_amrex(iX1,iX2,iX3,nDOFX*(iVar-1)+1:nDOFX*iVar) &
          = Data_thornado(1:nDOFX,iX1,iX2,iX3,iVar)
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE thornado2AMReX


  SUBROUTINE LinComb( alpha, MF_U, beta, MF_D )

    TYPE(amrex_multifab), INTENT(inout) :: MF_U(0:nLevels)
    TYPE(amrex_multifab), INTENT(in)    :: MF_D(0:nLevels)
    REAL(amrex_real),     INTENT(in)    :: alpha, beta(0:nLevels)

    INTEGER                               :: iX1, iX2, iX3, iLevel
    INTEGER                               :: lo(4), hi(4)
    TYPE(amrex_box)                       :: BX
    TYPE(amrex_mfiter)                    :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: U(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: D(:,:,:,:)

    INTEGER :: nComp, iComp

    DO iLevel = 0, nLevels
      CALL amrex_mfiter_build( MFI, MF_U(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        nComp = MF_U(iLevel) % nComp()

        U => MF_U(iLevel) % DataPtr( MFI )
        D => MF_D(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo = LBOUND( U ); hi = UBOUND( U )

        DO iComp = 1, nComp

          DO iX3 = BX % lo(3), BX % hi(3)
          DO iX2 = BX % lo(2), BX % hi(2)
          DO iX1 = BX % lo(1), BX % hi(1)

            U(iX1,iX2,iX3,iComp) = alpha * U(iX1,iX2,iX3,iComp) &
                                   + beta(iLevel) &
                                   * D(iX1,iX2,iX3,iComp)

          END DO
          END DO
          END DO

        END DO

      END DO

    END DO

  END SUBROUTINE LinComb


  SUBROUTINE ShowVariableFromMultiFab( MF, swX, iComp )

    INTEGER,              INTENT(in) :: swX(3)
    TYPE(amrex_multifab), INTENT(in) :: MF(0:nLevels)
    INTEGER,              INTENT(in) :: iComp

    INTEGER                               :: iX1, iX2, iX3, iLevel
    INTEGER                               :: lo(4), hi(4)
    TYPE(amrex_box)                       :: BX
    TYPE(amrex_mfiter)                    :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: U(:,:,:,:)

    DO iLevel = 0, nLevels
      CALL amrex_mfiter_build( MFI, MF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        U => MF(iLevel) % DataPtr( MFI )
        BX = MFI % tilebox()

        lo = LBOUND( U ); hi = UBOUND( U )

        DO iX3 = BX % lo(3) - swX(3), BX % hi(3) + swX(3)
        DO iX2 = BX % lo(2) - swX(2), BX % hi(2) + swX(2)
        DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)

          WRITE(*,'(A,3I4.3,ES10.1E3)') &
            'iX1, iX2, iX3, Data: ',iX1, iX2, iX3, U(iX1,iX2,iX3,iComp)

        END DO
        END DO
        END DO

      END DO

    END DO
    WRITE(*,*)

  END SUBROUTINE ShowVariableFromMultiFab


  SUBROUTINE MakeMF_Diff( chk1, chk2 )

    INTEGER, INTENT(in) :: chk1, chk2

    TYPE(amrex_box)                    :: BX
    TYPE(amrex_boxarray),  ALLOCATABLE :: BA(:)
    TYPE(amrex_distromap), ALLOCATABLE :: DM(:)
    TYPE(amrex_geometry),  ALLOCATABLE :: GEOM(:)
    TYPE(amrex_multifab),  ALLOCATABLE :: MF_uGF_TEMP(:)
    TYPE(amrex_multifab),  ALLOCATABLE :: MF_uCF_TEMP(:)
    TYPE(amrex_multifab),  ALLOCATABLE :: MF_uPF_TEMP(:)
    TYPE(amrex_multifab),  ALLOCATABLE :: MF_uAF_TEMP(:)
    INTEGER                            :: iLevel, iComp

    ! --- Initialize AMReX ---
    CALL amrex_init()
    CALL amrex_amrcore_init()

    ! --- Parse parameter file ---
    CALL MyAmrInit

    CALL InitializeProgramHeader &
           ( ProgramName_Option = TRIM( ProgramName ), &
             nNodes_Option = nNodes, nX_Option = nX, swX_Option = swX, &
             xL_Option = xL, xR_Option = xR, bcX_Option = bcX, &
             Verbose_Option = .FALSE. )

    CALL InitializeReferenceElementX

    ALLOCATE( BA         (0:nLevels) )
    ALLOCATE( DM         (0:nLevels) )
    ALLOCATE( GEOM       (0:nLevels) )
    ALLOCATE( MF_uGF_TEMP(0:nLevels) )
    ALLOCATE( MF_uCF_TEMP(0:nLevels) )
    ALLOCATE( MF_uPF_TEMP(0:nLevels) )
    ALLOCATE( MF_uAF_TEMP(0:nLevels) )

    BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

    DO iLevel = 0, nLevels
      CALL amrex_boxarray_build( BA(iLevel), BX )
      CALL BA(iLevel) % maxSize( MaxGridSize )
      CALL amrex_geometry_build ( GEOM(iLevel), BX )
      CALL amrex_distromap_build( DM  (iLevel), BA(iLevel) )

      CALL amrex_multifab_build &
             ( MF_uGF_TEMP(iLevel), BA(iLevel), DM(iLevel), &
               nDOFX * nGF, swX(1) )
      CALL MF_uGF_TEMP(iLevel) % SetVal( 0.0d0 )
      CALL amrex_multifab_build &
             ( MF_uCF_TEMP(iLevel), BA(iLevel), DM(iLevel), &
               nDOFX * nCF, swX(1) )
      CALL MF_uCF_TEMP(iLevel) % SetVal( 0.0d0 )
      CALL amrex_multifab_build &
             ( MF_uPF_TEMP(iLevel), BA(iLevel), DM(iLevel), &
               nDOFX * nPF, swX(1) )
      CALL MF_uPF_TEMP(iLevel) % SetVal( 0.0d0 )
      CALL amrex_multifab_build &
             ( MF_uAF_TEMP(iLevel), BA(iLevel), DM(iLevel), &
               nDOFX * nAF, swX(1) )
      CALL MF_uAF_TEMP(iLevel) % SetVal( 0.0d0 )
    END DO

    CALL MyAmrFinalize
    CALL ReadCheckpointFile( chk1 )

    DO iLevel = 0, nLevels
      CALL MF_uGF_TEMP(iLevel) % PARALLEL_COPY( MF_uGF(iLevel), GEOM(iLevel) )
      CALL MF_uCF_TEMP(iLevel) % PARALLEL_COPY( MF_uCF(iLevel), GEOM(iLevel) )
      CALL MF_uPF_TEMP(iLevel) % PARALLEL_COPY( MF_uPF(iLevel), GEOM(iLevel) )
      CALL MF_uAF_TEMP(iLevel) % PARALLEL_COPY( MF_uAF(iLevel), GEOM(iLevel) )
    END DO

    CALL MyAmrFinalize
    CALL ReadCheckpointFile( chk2 )

    DO iLevel = 0, nLevels
      CALL MF_uGF_TEMP(iLevel) &
             % SUBTRACT( MF_uGF(iLevel), 1, 1, &
                         MF_uGF(iLevel) % nComp(), swX(1) )
      CALL MF_uCF_TEMP(iLevel) &
             % SUBTRACT( MF_uCF(iLevel), 1, 1, &
                         MF_uCF(iLevel) % nComp(), swX(1) )
      CALL MF_uPF_TEMP(iLevel) &
             % SUBTRACT( MF_uPF(iLevel), 1, 1, &
                         MF_uPF(iLevel) % nComp(), swX(1) )
      CALL MF_uAF_TEMP(iLevel) &
             % SUBTRACT( MF_uAF(iLevel), 1, 1, &
                         MF_uAF(iLevel) % nComp(), swX(1) )
    END DO

!!$    ! --- Write min/max of each component of each field ---
!!$    DO iLevel = 0, nLevels
!!$      DO iComp = 1, MF_uGF_TEMP(iLevel) % nComp()
!!$        WRITE(*,'(A,ES24.16E3)') 'MinG:   ', MF_uGF_TEMP(iLevel) % MIN(iComp)
!!$        WRITE(*,'(A,ES24.16E3)') 'MaxG:   ', MF_uGF_TEMP(iLevel) % MAX(iComp)
!!$      END DO
!!$      WRITE(*,*)
!!$      DO iComp = 1, MF_uCF_TEMP(iLevel) % nComp()
!!$        WRITE(*,'(A,ES24.16E3)') 'MinC:   ', MF_uCF_TEMP(iLevel) % MIN(iComp)
!!$        WRITE(*,'(A,ES24.16E3)') 'MaxC:   ', MF_uCF_TEMP(iLevel) % MAX(iComp)
!!$      END DO
!!$      WRITE(*,*)
!!$      DO iComp = 1, MF_uPF_TEMP(iLevel) % nComp()
!!$        WRITE(*,'(A,ES24.16E3)') 'MinP:   ', MF_uPF_TEMP(iLevel) % MIN(iComp)
!!$        WRITE(*,'(A,ES24.16E3)') 'MaxP:   ', MF_uPF_TEMP(iLevel) % MAX(iComp)
!!$      END DO
!!$      WRITE(*,*)
!!$      DO iComp = 1, MF_uAF_TEMP(iLevel) % nComp()
!!$        WRITE(*,'(A,ES24.16E3)') 'MinA:   ', MF_uAF_TEMP(iLevel) % MIN(iComp)
!!$        WRITE(*,'(A,ES24.16E3)') 'MaxA:   ', MF_uAF_TEMP(iLevel) % MAX(iComp)
!!$      END DO
!!$      WRITE(*,*)
!!$    END DO

    CALL WriteFieldsAMReX_Plotfile &
           ( t(0), GEOM, StepNo, &
             MF_uGF_Option = MF_uGF_TEMP, &
             MF_uCF_Option = MF_uCF_TEMP, &
             MF_uPF_Option = MF_uPF_TEMP, &
             MF_uAF_Option = MF_uAF_TEMP )

    CALL MyAmrFinalize

    CALL FinalizeReferenceElementX

    DO iLevel = 0, nLevels
      CALL amrex_geometry_destroy ( GEOM(iLevel) )
      CALL amrex_distromap_destroy( DM  (iLevel) )
      CALL amrex_boxarray_destroy ( BA  (iLevel) )
    END DO
    DEALLOCATE( GEOM )
    DEALLOCATE( DM   )
    DEALLOCATE( BA   )

    DO iLevel = 0, nLevels
      CALL amrex_multifab_destroy( MF_uAF_TEMP(iLevel) )
      CALL amrex_multifab_destroy( MF_uPF_TEMP(iLevel) )
      CALL amrex_multifab_destroy( MF_uCF_TEMP(iLevel) )
      CALL amrex_multifab_destroy( MF_uGF_TEMP(iLevel) )
    END DO
    DEALLOCATE( MF_uAF_TEMP )
    DEALLOCATE( MF_uPF_TEMP )
    DEALLOCATE( MF_uCF_TEMP )
    DEALLOCATE( MF_uGF_TEMP )

    CALL amrex_amrcore_finalize()
    CALL amrex_finalize()

    STOP 'MF_UtilitiesModule.f90'

  END SUBROUTINE MakeMF_Diff


END MODULE MF_UtilitiesModule
