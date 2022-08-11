MODULE MF_Euler_PositivityLimiterModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator, &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDOFX
  USE FluidFieldsModule, ONLY: &
    nCF, &
    nDF
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X
  USE InputParsingModule, ONLY: &
    nLevels, &
    UsePositivityLimiter, &
    UseTiling, &
    DEBUG
!!$  USE AverageDownModule, ONLY: &
!!$    AverageDown
!!$  USE FillPatchModule, ONLY: &
!!$    FillPatch
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyPositivityLimiter_Euler_MF

  INTERFACE ApplyPositivityLimiter_Euler_MF
    MODULE PROCEDURE ApplyPositivityLimiter_Euler_MF_MultipleLevels
    MODULE PROCEDURE ApplyPositivityLimiter_Euler_MF_SingleLevel
  END INTERFACE ApplyPositivityLimiter_Euler_MF

CONTAINS


  SUBROUTINE ApplyPositivityLimiter_Euler_MF_MultipleLevels &
    ( MF_uGF, MF_uCF, MF_uDF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:)

    INTEGER, PARAMETER :: nCycles = 1

    INTEGER :: iCycle, iLevel, iErr

    DO iCycle = 1, nCycles

      DO iLevel = 0, nLevels-1

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,'(2x,A,I3.3)') &
              'CALL ApplyPositivityLimiter_Euler_MF_SingleLevel, iLevel: ', &
              iLevel

          END IF

        END IF ! DEBUG

        CALL ApplyPositivityLimiter_Euler_MF_SingleLevel &
               ( iLevel, MF_uGF, MF_uCF, MF_uDF )

!!$        CALL FillPatch( iLevel, 0.0_DP, MF_uGF, MF_uCF )

      END DO ! iLevel

      ! --- Ensure underlying coarse cells are consistent with
      !     cells on refined level ---

!!$      CALL AverageDown( MF_uGF, MF_uCF )

    END DO ! iCycle

  END SUBROUTINE ApplyPositivityLimiter_Euler_MF_MultipleLevels


  SUBROUTINE ApplyPositivityLimiter_Euler_MF_SingleLevel &
    ( iLevel, MF_uGF, MF_uCF, MF_uDF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D(:,:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF(iLevel) % DataPtr( MFI )
      uCF => MF_uCF(iLevel) % DataPtr( MFI )
      uDF => MF_uDF(iLevel) % DataPtr( MFI )

      iLo_MF = LBOUND( uGF )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = BX % lo - swX
      iX_E1 = BX % hi + swX

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                          iX_B1(2):iX_E1(2), &
                          iX_B1(3):iX_E1(3),1:nGF) )

      ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                          iX_B1(2):iX_E1(2), &
                          iX_B1(3):iX_E1(3),1:nCF) )

      ALLOCATE( D(1:nDOFX,iX_B1(1):iX_E1(1), &
                          iX_B1(2):iX_E1(2), &
                          iX_B1(3):iX_E1(3),1:nDF) )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

      CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

      CALL amrex2thornado_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )

      CALL ApplyPositivityLimiter_Euler( iX_B1, iX_E1, iX_B1, iX_E1, G, U, D )

      CALL thornado2amrex_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

      CALL thornado2amrex_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      DEALLOCATE( D )

      DEALLOCATE( U )

      DEALLOCATE( G )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ApplyPositivityLimiter_Euler_MF_SingleLevel


END MODULE MF_Euler_PositivityLimiterModule
