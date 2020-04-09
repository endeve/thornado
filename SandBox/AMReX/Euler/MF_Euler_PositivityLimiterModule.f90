MODULE MF_Euler_PositivityLimiterModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,     ONLY: &
    AR => amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule,           ONLY: &
    swX, nDOFX
  USE FluidFieldsModule,             ONLY: &
    nCF
  USE GeometryFieldsModule,          ONLY: &
    nGF
  USE Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler

  ! --- Local Modules ---

  USE MF_UtilitiesModule,       ONLY: &
    AMReX2thornado, &
    thornado2AMReX
  USE MyAmrModule,              ONLY: &
    nLevels,              &
    UsePositivityLimiter, &
    DEBUG
  USE TimersModule_AMReX_Euler, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler,  &
    Timer_AMReX_Euler_DataTransfer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ApplyPositivityLimiter_Euler


CONTAINS


  SUBROUTINE MF_ApplyPositivityLimiter_Euler( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(AR), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: G(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF ) )

        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF ) )

        CALL AMReX2thornado( nGF, iX_B1, iX_E1, uGF, G )

        CALL AMReX2thornado( nCF, iX_B1, iX_E1, uCF, U )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

        IF( DEBUG ) WRITE(*,'(A)') '    CALL ApplyPositivityLimiter_Euler'

        CALL ApplyPositivityLimiter_Euler( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

        CALL thornado2AMReX( nCF, iX_B1, iX_E1, uCF, U )

        DEALLOCATE( U )

        DEALLOCATE( G )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ApplyPositivityLimiter_Euler


END MODULE MF_Euler_PositivityLimiterModule
