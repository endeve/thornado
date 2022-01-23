!> Module for operations on MultiFabs
MODULE MF_UtilitiesModule

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_DataTransfer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: amrex2thornado_X
  PUBLIC :: thornado2amrex_X
  PUBLIC :: thornado2amrex_X_F
  PUBLIC :: amrex2thornado_X_F


CONTAINS


  SUBROUTINE amrex2thornado_X &
    ( nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      Data_amrex   (iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(out) :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_thornado(1:nDOFX,iX1,iX2,iX3,iFd) &
        = Data_amrex(iX1,iX2,iX3,nDOFX*(iFd-1)+1:nDOFX*iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

  END SUBROUTINE amrex2thornado_X


  SUBROUTINE thornado2amrex_X &
    ( nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(DP), INTENT(out) :: &
      Data_amrex   (iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(in)  :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_amrex(iX1,iX2,iX3,nDOFX*(iFd-1)+1:nDOFX*iFd) &
        = Data_thornado(1:nDOFX,iX1,iX2,iX3,iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

  END SUBROUTINE thornado2amrex_X


  SUBROUTINE thornado2amrex_X_F &
    ( nDOFX_X, nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, &
      Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nDOFX_X, nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(DP), INTENT(out) :: &
      Data_amrex   (iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(in)  :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_amrex(iX1,iX2,iX3,nDOFX_X*(iFd-1)+1:nDOFX_X*iFd) &
        = Data_thornado(1:nDOFX_X,iX1,iX2,iX3,iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

  END SUBROUTINE thornado2amrex_X_F


  SUBROUTINE amrex2thornado_X_F &
    ( nDOFX_X, nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, &
      Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nDOFX_X, nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      Data_amrex   (iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(out) :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_thornado(1:nDOFX_X,iX1,iX2,iX3,iFd) &
        = Data_amrex(iX1,iX2,iX3,nDOFX_X*(iFd-1)+1:nDOFX_X*iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

  END SUBROUTINE amrex2thornado_X_F


END MODULE MF_UtilitiesModule
