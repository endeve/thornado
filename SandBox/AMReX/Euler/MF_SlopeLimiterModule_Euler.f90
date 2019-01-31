MODULE MF_SlopeLimiterModule_Euler

  USE amrex_base_module
  USE amrex_fort_module

  USE ProgramHeaderModule, ONLY: &
    swX, nDOFX
  USE FluidFieldsModule, ONLY: &
    nCF
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE SlopeLimiterModule_Euler, ONLY: &
    ApplySlopeLimiter_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ApplySlopeLimiter_Euler


CONTAINS


  SUBROUTINE MF_ApplySlopeLimiter_Euler( nLevels, MF_uGF, MF_uCF )

    INTEGER,              INTENT(in)    :: nLevels
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(amrex_real), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: G(:,:,:,:,:)

    INTEGER :: iLevel

    DO iLevel = 0, nLevels

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        ALLOCATE( U(1:nDOFX,BX%lo(1)-swX(1):BX%hi(1)+swX(1), &
                            BX%lo(2)-swX(2):BX%hi(2)+swX(2), &
                            BX%lo(3)-swX(3):BX%hi(3)+swX(3),1:nCF) )
        ALLOCATE( G(1:nDOFX,BX%lo(1)-swX(1):BX%hi(1)+swX(1), &
                            BX%lo(2)-swX(2):BX%hi(2)+swX(2), &
                            BX%lo(3)-swX(3):BX%hi(3)+swX(3),1:nGF) )

        CALL AMReX2thornado( nCF, BX, uCF, U )
        CALL AMReX2thornado( nGF, BX, uGF, G )

        CALL ApplySlopeLimiter_Euler &
               ( BX % lo, BX % hi, ( BX % lo ) - swX, ( BX % hi ) + swX, G, U )

        CALL thornado2AMReX( nCF, BX, uCF, U )
        CALL thornado2AMReX( nGF, BX, uGF, G )

        DEALLOCATE( G )
        DEALLOCATE( U )

      END DO

    END DO


  END SUBROUTINE MF_ApplySlopeLimiter_Euler


  SUBROUTINE AMReX2thornado( nVars, BX, Data_amrex, Data_thornado )

    INTEGER,          INTENT(in)    :: nVars
    TYPE(amrex_box),  INTENT(in)    :: BX
    REAL(amrex_real), INTENT(in)    :: &
      Data_amrex(BX%lo(1):,BX%lo(2):,BX%lo(3):,1:)
    REAL(amrex_real), INTENT(inout) :: &
      Data_thornado(1:nDOFX,BX%lo(1)-swX(1):BX%hi(1)+swX(1), &
                            BX%lo(2)-swX(2):BX%hi(2)+swX(2), &
                            BX%lo(3)-swX(3):BX%hi(3)+swX(3),1:nVars)
    INTEGER :: iX1, iX2, iX3, iY1, iY2, iY3, iVar

    DO iX3 = BX % lo(3) - swX(3), BX % hi(3) + swX(3)
      iY3 = iX3 + swX(3)
    DO iX2 = BX % lo(2) - swX(2), BX % hi(2) + swX(2)
      iY2 = iX2 + swX(2)
    DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)
      iY1 = iX1 + swX(1)

      DO iVar = 1, nVars
        Data_thornado(1:nDOFX,iX1,iX2,iX3,iVar) &
          = Data_amrex(iY1,iY2,iY3,nDOFX*(iVar-1)+1:nDOFX*iVar)
      END DO

    END DO
    END DO
    END DO


  END SUBROUTINE AMReX2thornado


  SUBROUTINE thornado2AMReX( nVars, BX, Data_amrex, Data_thornado )

    INTEGER,          INTENT(in)    :: nVars
    TYPE(amrex_box),  INTENT(in)    :: BX
    REAL(amrex_real), INTENT(inout) :: &
      Data_amrex(BX%lo(1):,BX%lo(2):,BX%lo(3):,1:)
    REAL(amrex_real), INTENT(in)    :: &
      Data_thornado(1:nDOFX,BX%lo(1)-swX(1):BX%hi(1)+swX(1), &
                            BX%lo(2)-swX(2):BX%hi(2)+swX(2), &
                            BX%lo(3)-swX(3):BX%hi(3)+swX(3),1:nVars)
    INTEGER :: iX1, iX2, iX3, iY1, iY2, iY3, iVar

    DO iX3 = BX % lo(3) - swX(3), BX % hi(3) + swX(3)
      iY3 = iX3 + swX(3)
    DO iX2 = BX % lo(2) - swX(2), BX % hi(2) + swX(2)
      iY2 = iX2 + swX(2)
    DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)
      iY1 = iX1 + swX(1)

      DO iVar = 1, nVars
        Data_amrex(iY1,iY2,iY3,nDOFX*(iVar-1)+1:nDOFX*iVar) &
          = Data_thornado(1:nDOFX,iX1,iX2,iX3,iVar)
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE thornado2AMReX


END MODULE MF_SlopeLimiterModule_Euler
