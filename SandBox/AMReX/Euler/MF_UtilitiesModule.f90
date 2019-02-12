MODULE MF_UtilitiesModule

  ! --- AMReX Modules ---
  USE amrex_base_module, ONLY: &
    amrex_finalize, &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_box, &
    amrex_boxarray, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_distromap
  USE amrex_amr_module, ONLY: &
    amrex_amrcore_finalize
  USE amrex_fort_module, ONLY: &
    amrex_real

  ! --- thornado Modules ---
  USE ProgramHeaderModule, ONLY: &
    nDOFX, swX
  USE FluidFieldsModule, ONLY: &
    nCF, nPF, nAF
  USE GeometryFieldsModule, ONLY: &
    nGF

  USE MyAmrModule,            ONLY: &
    MyAmrInit, MyAmrFinalize
  USE MyAmrDataModule
  USE InputOutputModuleAMReX, ONLY: &
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


  SUBROUTINE LinComb( nLevels, alpha, MF_U, beta, MF_D )

    INTEGER,              INTENT(in)    :: nLevels
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


  SUBROUTINE ShowVariableFromMultiFab( nLevels, MF, iComp )

    INTEGER,              INTENT(in) :: nLevels
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


  SUBROUTINE MakeMF_Diff( nLevels, BA, DM, chk1, chk2 )

    INTEGER,               INTENT(in) :: nLevels, chk1, chk2
    TYPE(amrex_boxarray),  INTENT(in) :: BA(0:nLevels)
    TYPE(amrex_distromap), INTENT(in) :: DM(0:nLevels)

    TYPE(amrex_multifab) :: MF_uGF_TEMP(0:nLevels)
    TYPE(amrex_multifab) :: MF_uCF_TEMP(0:nLevels)
    TYPE(amrex_multifab) :: MF_uPF_TEMP(0:nLevels)
    TYPE(amrex_multifab) :: MF_uAF_TEMP(0:nLevels)
    INTEGER              :: iLevel

    DO iLevel = 0, nLevels
      CALL amrex_multifab_build &
             ( MF_uGF_TEMP(iLevel), BA(iLevel), DM(iLevel), &
               nDOFX * nGF, swX(1) )
      CALL amrex_multifab_build &
             ( MF_uCF_TEMP(iLevel), BA(iLevel), DM(iLevel), &
               nDOFX * nCF, swX(1) )
      CALL amrex_multifab_build &
             ( MF_uPF_TEMP(iLevel), BA(iLevel), DM(iLevel), &
               nDOFX * nPF, swX(1) )
      CALL amrex_multifab_build &
             ( MF_uAF_TEMP(iLevel), BA(iLevel), DM(iLevel), &
               nDOFX * nAF, swX(1) )
    END DO

    CALL MyAmrFinalize
    CALL ReadCheckpointFile( chk1 )

    DO iLevel = 0, nLevels
      CALL MF_uCF_TEMP(iLevel) &
             % COPY( MF_uCF(iLevel), 1, 1, MF_uCF(iLevel) % nComp(), swX(1) )
    END DO

    CALL MyAmrFinalize
    CALL ReadCheckpointFile( chk2 )

    DO iLevel = 0, nLevels
      CALL MF_uCF_TEMP(iLevel) &
             % SUBTRACT( MF_uCF(iLevel), 1, 1, &
                         MF_uCF(iLevel) % nComp(), swX(1) )
    END DO

    CALL ShowVariableFromMultifab( nLevels, MF_uCF_TEMP, 2 )

    CALL MyAmrFinalize

    DO iLevel = 0, nLevels
      CALL amrex_multifab_destroy( MF_uGF_TEMP(iLevel) )
      CALL amrex_multifab_destroy( MF_uCF_TEMP(iLevel) )
      CALL amrex_multifab_destroy( MF_uPF_TEMP(iLevel) )
      CALL amrex_multifab_destroy( MF_uAF_TEMP(iLevel) )
    END DO

    CALL amrex_amrcore_finalize()
    CALL amrex_finalize()

    STOP 'MF_UtilitiesModule.f90'

  END SUBROUTINE MakeMF_Diff


END MODULE MF_UtilitiesModule
