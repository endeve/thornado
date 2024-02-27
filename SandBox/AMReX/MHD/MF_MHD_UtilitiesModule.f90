MODULE MF_MHD_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_box_module,         ONLY: &
    amrex_box
  USE amrex_multifab_module,    ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,    ONLY: &
    amrex_parallel_reduce_min, &
    amrex_parallel_ioprocessor
  USE amrex_amr_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule,      ONLY: &
    nDOFX, &
    swX
  USE GeometryFieldsModule,     ONLY: &
    nGF
  USE MagnetofluidFieldsModule,        ONLY: &
    nCM, &
    nPM, &
    nAM, &
    nDM
  USE MHD_UtilitiesModule,    ONLY: &
    ComputeTimeStep_MHD, &
    ComputeFromConserved_MHD
  USE MHD_UtilitiesModule_Relativistic,    ONLY: &
    ComputeMagneticDivergence_MHD_Relativistic
  USE EquationOfStateModule,    ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  ! --- Local Modules ---

  USE MF_KindModule,            ONLY: &
    DP, &
    One
  USE InputParsingModule,       ONLY: &
    nLevels, &
    UseTiling, &
    EvolveOnlyMagnetic, &
    UseDivergenceCleaning
  USE MF_UtilitiesModule,       ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeFromConserved
  PUBLIC :: MF_ComputeTimeStep
  PUBLIC :: MF_ComputeMagneticDivergence

CONTAINS


  SUBROUTINE MF_ComputeFromConserved( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

    TYPE(amrex_multifab), INTENT(in)    :: &
      MF_uGF(0:nLevels-1), MF_uCM(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: &
      MF_uPM(0:nLevels-1), MF_uAM(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAM(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: A(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )
        uPM => MF_uPM(iLevel) % DataPtr( MFI )
        uAM => MF_uAM(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCM) )

        ALLOCATE( P(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nPM) )

        ALLOCATE( A(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nAM) )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCM, U )

        CALL amrex2thornado_X( nPM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uPM, P )

        CALL amrex2thornado_X( nAM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uAM, A )

        CALL ComputeFromConserved_MHD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A, EvolveOnlyMagnetic )

        CALL thornado2amrex_X( nPM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uPM, P )

        CALL thornado2amrex_X( nAM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uAM, A )

        DEALLOCATE( A )

        DEALLOCATE( P )

        DEALLOCATE( U )

        DEALLOCATE( G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ComputeFromConserved


  SUBROUTINE MF_ComputeTimeStep( MF_uGF, MF_uCM, CFL, TimeStepMin )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:nLevels-1), &
                                         MF_uCM(0:nLevels-1)

    REAL(DP),             INTENT(in)  :: CFL
    REAL(DP),             INTENT(out) :: TimeStepMin(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    REAL(DP) :: TimeStep(0:nLevels-1)

    TimeStepMin = HUGE( One )

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCM) )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCM, U )

        CALL ComputeTimeStep_MHD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep(iLevel), &
                 UseDivergenceCleaning, EvolveOnlyMagnetic )

        TimeStepMin(iLevel) = MIN( TimeStepMin(iLevel), TimeStep(iLevel) )

        DEALLOCATE( U )

        DEALLOCATE( G )

      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

    END DO ! --- Loop over levels ---

    CALL amrex_parallel_reduce_min( TimeStepMin, nLevels )

  END SUBROUTINE MF_ComputeTimeStep


  SUBROUTINE MF_ComputeMagneticDivergence( MF_uGF, MF_uCM, MF_uDM )

    TYPE(amrex_multifab), INTENT(in)    :: &
      MF_uGF(0:nLevels-1), MF_uCM(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: &
      MF_uDM(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDM(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    DO iLevel = 0, nLevels-1

      CALL MF_uGF(iLevel) % Fill_Boundary( amrex_geom(iLevel) )
      CALL MF_uCM(iLevel) % Fill_Boundary( amrex_geom(iLevel) )
      CALL MF_uDM(iLevel) % Fill_Boundary( amrex_geom(iLevel) )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )
        uDM => MF_uDM(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCM) )

        ALLOCATE( D(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nDM) )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCM, U )

        CALL amrex2thornado_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )

        CALL ComputeMagneticDivergence_MHD_Relativistic &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

        CALL thornado2amrex_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )

        DEALLOCATE( D )

        DEALLOCATE( U )

        DEALLOCATE( G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ComputeMagneticDivergence


END MODULE MF_MHD_UtilitiesModule
