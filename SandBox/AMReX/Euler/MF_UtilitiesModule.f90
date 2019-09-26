!> Module for operations on MultiFabs
!> @todo Modify data transfer subroutine to allow specification of looping over
!>       ghost cells
MODULE MF_UtilitiesModule

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build

  ! --- thornado Modules ---
  USE ProgramHeaderModule, ONLY: &
    nDOFX

  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: AMReX2thornado
  PUBLIC :: thornado2AMReX
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


END MODULE MF_UtilitiesModule
