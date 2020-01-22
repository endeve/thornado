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
    nDOFZ,                 &
    nDOFX

  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: AMReX2thornado
  PUBLIC :: thornado2AMReX
  PUBLIC :: AMReX2thornado_Euler
  PUBLIC :: thornado2AMReX_Euler

  CONTAINS

  SUBROUTINE AMReX2thornado( nVars, nS, nE, iZ_B, iZ_E, Data_amrex, Data_thornado )

    INTEGER,          INTENT(in)  :: nVars, nS, nE
    INTEGER,          INTENT(in)  :: iZ_B(4), iZ_E(4)
    REAL(amrex_real), INTENT(in)  :: &
      Data_amrex   (   iZ_B(2):,iZ_B(3):,iZ_B(4):,1:)
    REAL(amrex_real), INTENT(out) :: &
      Data_thornado(1:,iZ_B(1):,iZ_B(2):,iZ_B(3):,iZ_B(4):,1:,1:)
    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iVar

    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)
      
      DO iS = 1, nS
      DO iZ1 = iZ_B(1), iZ_E(1)
      DO iVar = 1, nVars
        Data_thornado(1:nDOFZ,iZ1,iZ2,iZ3,iZ4,iVar,iS) &
          = Data_amrex(iZ2,iZ3,iZ4,&
            nDOFZ*(iVar-1)+1+(iZ1-1)*nVars*nDOFZ+(iS-1)*nVars*nDOFZ*nE: &
            nDOFZ*iVar+(iZ1-1)*nVars*nDOFZ+(iS-1)*nVars*nDOFZ*nE)
      END DO
      END DO
      END DO
      
    END DO
    END DO
    END DO

  END SUBROUTINE AMReX2thornado


  SUBROUTINE thornado2AMReX( nVars, nS, nE, iZ_B, iZ_E, Data_amrex, Data_thornado )


    INTEGER,          INTENT(in)  :: nVars, nS, nE
    INTEGER,          INTENT(in)  :: iZ_B(4), iZ_E(4)
    REAL(amrex_real), INTENT(out)  :: &
      Data_amrex   (   iZ_B(2):,iZ_B(3):,iZ_B(4):,1:)
    REAL(amrex_real), INTENT(in) :: &
      Data_thornado(1:,iZ_B(1):,iZ_B(2):,iZ_B(3):,iZ_B(4):,1:,1:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iVar
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)
      
      DO iS = 1, nS
      DO iZ1 = iZ_B(1), iZ_E(1)
      DO iVar = 1, nVars
        Data_amrex(iZ2,iZ3,iZ4,&
        nDOFZ*(iVar-1)+1+(iZ1-1)*nVars*nDOFZ+(iS-1)*nVars*nDOFZ*nE: &
        nDOFZ*iVar+(iZ1-1)*nVars*nDOFZ+(iS-1)*nVars*nDOFZ*nE) & 
          =  Data_thornado(1:nDOFZ,iZ1,iZ2,iZ3,iZ4,iVar,iS) 
      END DO
      END DO
      END DO
      
    END DO
    END DO
    END DO

  END SUBROUTINE thornado2AMReX

  SUBROUTINE AMReX2thornado_Euler( nVars, iX_B, iX_E, Data_amrex, Data_thornado )

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

  END SUBROUTINE AMReX2thornado_Euler


  SUBROUTINE thornado2AMReX_Euler( nVars, iX_B, iX_E, Data_amrex, Data_thornado )

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

  END SUBROUTINE thornado2AMReX_Euler

END MODULE MF_UtilitiesModule
