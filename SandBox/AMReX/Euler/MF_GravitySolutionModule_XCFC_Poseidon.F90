MODULE MF_GravitySolutionModule_XCFC_Poseidon

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE GravitySolutionModule_CFA_Poseidon, ONLY: &
    UpdateConformalFactorAndMetric

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    swX
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X
  USE TimersModule_AMReX_Euler, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler,  &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: UpdateConformalFactorAndMetric_MF

CONTAINS


  SUBROUTINE UpdateConformalFactorAndMetric_MF( MF_Psi, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_Psi(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iLo_G(4), iLo_Psi(4)

    REAL(DP), CONTIGUOUS, POINTER :: uPsi(:,:,:,:)
    REAL(DP), ALLOCATABLE         :: Psi (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G   (:,:,:,:,:)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uPsi => MF_Psi(iLevel) % DataPtr( MFI )
        uGF  => MF_uGF(iLevel) % DataPtr( MFI )

        iLo_Psi = LBOUND( uPsi )
        iLo_G   = LBOUND( uGF  )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        ALLOCATE( Psi(1:nDOFX,iX_B0(1):iX_E0(1), &
                              iX_B0(2):iX_E0(2), &
                              iX_B0(3):iX_E0(3)) )

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_X &
               ( iX_B0, iX_E0, iLo_Psi, iX_B0, iX_E0, uPsi, Psi )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_G, iX_B1, iX_E1, uGF, G )

        CALL UpdateConformalFactorAndMetric &
               ( iX_B0, iX_E0, iX_B1, iX_E1, Psi, G )

        CALL thornado2amrex_X &
               ( nGF, iX_B1, iX_E1, iLo_G, iX_B1, iX_E1, uGF, G )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( G   )
        DEALLOCATE( Psi )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE UpdateConformalFactorAndMetric_MF


END MODULE MF_GravitySolutionModule_XCFC_Poseidon
