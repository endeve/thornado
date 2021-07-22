MODULE ProbModule

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshType, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne
  USE Euler_UtilitiesModule, ONLY: &
    ComputeConserved_Euler

  ! --- Local Modules ---

  USE AMReX_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    TwoPi
  USE MF_UtilitiesModule, ONLY: &
    thornado2amrex_X
  USE InputParsingModule, ONLY: &
    nX, &
    swX, &
    Gamma_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitProbData

CONTAINS


  SUBROUTINE InitProbData &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uCF, iLo, iHi, MeshX )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
                               iLo(4), iHi(4)
    REAL(DP), INTENT(inout) :: uCF(iLo(1):iHi(1),iLo(2):iHi(2), &
                                   iLo(3):iHi(3),iLo(4):iHi(4))

    TYPE(MeshType), INTENT(in) :: MeshX(3)

    INTEGER        :: iNX, iX1, iX2, iX3, iNX1, iNX2
    REAL(DP)       :: X1, X2
    REAL(DP)       :: uPF(1:nDOFX,1:nPF)
    REAL(DP)       :: U(1:nDOFX,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3), &
                        1:nCF)

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      iNX1 = NodeNumberTableX(1,iNX)
      X1   = NodeCoordinate( MeshX(1), iX1+1, iNX1 )

      iNX2 = NodeNumberTableX(2,iNX)
      X2   = NodeCoordinate( MeshX(2), iX2+1, iNX2 )

      uPF(iNX,iPF_D ) = 1.0_DP + 0.1_DP * SIN( TwoPi * X1 )
      uPF(iNX,iPF_V1) = 0.1_DP
      uPF(iNX,iPF_V2) = 0.0_DP
      uPF(iNX,iPF_V3) = 0.0_DP
      uPF(iNX,iPF_E ) = One / ( Gamma_IDEAL - One )
      uPF(iNX,iPF_Ne) = Zero

      CALL ComputeConserved_Euler &
             ( uPF(iNX,iPF_D ), &
               uPF(iNX,iPF_V1), &
               uPF(iNX,iPF_V2), &
               uPF(iNX,iPF_V3), &
               uPF(iNX,iPF_E ), &
               uPF(iNX,iPF_Ne), &
               U(iNX,iX1,iX2,iX3,iCF_D ), &
               U(iNX,iX1,iX2,iX3,iCF_S1), &
               U(iNX,iX1,iX2,iX3,iCF_S2), &
               U(iNX,iX1,iX2,iX3,iCF_S3), &
               U(iNX,iX1,iX2,iX3,iCF_E ), &
               U(iNX,iX1,iX2,iX3,iCF_Ne), &
               One, One, One, &
               One )

    END DO
    END DO
    END DO
    END DO

    CALL thornado2amrex_X &
      ( nCF, iX_B0, iX_E0, iLo, iX_B0, iX_E0, uCF, U )

  END SUBROUTINE InitProbData


END MODULE ProbModule
