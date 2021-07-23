MODULE MF_GeometryModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MF_UtilitiesModule, ONLY: &
    thornado2amrex_X
  USE InputParsingModule, ONLY: &
    swX, &
    iOS_CPP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeGeometryX_MF


CONTAINS


  SUBROUTINE ComputeGeometryX_MF( MF_uGF )

    TYPE(amrex_multifab), INTENT(in) :: MF_uGF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER :: iNX, iX1, iX2, iX3, iGF
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: Gt(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iLo_G(3), iHi_G(3)

    REAL(DP) :: Mass

    Mass = Zero
    CALL amrex_parmparse_build( PP, 'thornado' )

      CALL PP % query( 'Mass', Mass )

    CALL amrex_parmparse_destroy( PP )

    CALL amrex_mfiter_build( MFI, MF_uGF )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      iLo_G = iX_B1 + iOS_CPP
      iHi_G = iX_E1 + iOS_CPP

      ALLOCATE( G (1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3), &
                   1:nGF) )

      ALLOCATE( Gt(1:nDOFX,iLo_G(1):iHi_G(1), &
                           iLo_G(2):iHi_G(2), &
                           iLo_G(3):iHi_G(3), &
                   1:nGF) )

#if defined HYDRO_RELATIVISTIC

      CALL ComputeGeometryX &
             ( iX_B0, iX_E0, iLo_G, iHi_G, Gt, Mass_Option = Mass )

#else

      CALL ComputeGeometryX &
             ( iX_B0, iX_E0, iLo_G, iHi_G, Gt )

#endif

      DO iGF = 1, nGF
      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)
      DO iNX = 1, nDOFX

        G(iNX,iX1,iX2,iX3,iGF) &
          = Gt(iNX,iX1+iOS_CPP(1),iX2+iOS_CPP(2),iX3+iOS_CPP(3),iGF)

      END DO
      END DO
      END DO
      END DO
      END DO

      CALL thornado2amrex_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

      DEALLOCATE( Gt )
      DEALLOCATE( G  )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ComputeGeometryX_MF


END MODULE MF_GeometryModule
