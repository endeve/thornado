MODULE MF_UtilitiesModule

  ! --- AMReX Modules ---
  USE amrex_base_module, ONLY: &
    amrex_multifab, &
    amrex_box, &
    amrex_mfiter, &
    amrex_mfiter_build
  USE amrex_fort_module, ONLY: &
    amrex_real

  ! --- thornado Modules ---
  USE ProgramHeaderModule, ONLY: &
    nDOFX, swX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: AMReX2thornado
  PUBLIC :: thornado2AMReX
  PUBLIC :: LinComb
  PUBLIC :: ShowVariableFromMultiFab


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


  SUBROUTINE LinComb( nLevels, alpha, MF_U, beta, MF_D, iS )

    INTEGER,              INTENT(in)    :: nLevels, iS
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
                                   * D(iX1,iX2,iX3,(iS-1)*nComp+iComp)

          END DO
          END DO
          END DO

        END DO

      END DO
    END DO

  END SUBROUTINE LinComb


  SUBROUTINE ShowVariableFromMultiFab( nLevels, MF, iComp )

    INTEGER,              INTENT(in)    :: nLevels
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

          WRITE(*,'(A,3I3.2,ES10.1E3)') &
            'iX1, iX2, iX3, Data: ',iX1, iX2, iX3, U(iX1,iX2,iX3,iComp)

        END DO
        END DO
        END DO

      END DO

    END DO
    WRITE(*,*)

  END SUBROUTINE ShowVariableFromMultiFab


END MODULE MF_UtilitiesModule
