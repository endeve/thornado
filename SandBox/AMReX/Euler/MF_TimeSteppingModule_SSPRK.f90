MODULE MF_TimeSteppingModule_SSPRK

  ! --- AMReX Modules ---

  USE amrex_fort_module,      ONLY: &
    AR => amrex_real
  USE amrex_box_module,       ONLY: &
    amrex_box
  USE amrex_geometry_module,  ONLY: &
    amrex_geometry
  USE amrex_multifab_module,  ONLY: &
    amrex_multifab,         &
    amrex_multifab_build,   &
    amrex_multifab_destroy, &
    amrex_mfiter,           &
    amrex_mfiter_build,     &
    amrex_mfiter_destroy
  USE amrex_boxarray_module,  ONLY: &
    amrex_boxarray,       &
    amrex_boxarray_build, &
    amrex_boxarray_destroy
  USE amrex_distromap_module, ONLY: &
    amrex_distromap,       &
    amrex_distromap_build, &
    amrex_distromap_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule,  ONLY: &
    swX,   &
    nDOFX, &
    nX
  USE FluidFieldsModule,    ONLY: &
    nCF
  USE GeometryFieldsModule, ONLY: &
    nGF

  ! --- Local Modules ---

  USE MF_Euler_SlopeLimiterModule,      ONLY: &
    MF_ApplySlopeLimiter_Euler
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    MF_ApplyPositivityLimiter_Euler
  USE MyAmrModule,                      ONLY: &
    nLevels, &
    DEBUG
  USE TimersModule_AMReX_Euler,         ONLY: &
    TimersStart_AMReX_Euler,       &
    TimersStop_AMReX_Euler,        &
    Timer_AMReX_Euler_UpdateFluid, &
    Timer_AMReX_Euler_InteriorBC,  &
    Timer_AMReX_Euler_CopyMultiFab

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFluid_SSPRK
  PUBLIC :: MF_UpdateFluid_SSPRK
  PUBLIC :: MF_FinalizeFluid_SSPRK

  INTEGER :: nStages_SSPRK
  REAL(AR), DIMENSION(:),   ALLOCATABLE :: c_SSPRK
  REAL(AR), DIMENSION(:),   ALLOCATABLE :: w_SSPRK
  REAL(AR), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  TYPE(amrex_multifab), DIMENSION(:),   ALLOCATABLE :: MF_U
  TYPE(amrex_multifab), DIMENSION(:,:), ALLOCATABLE :: MF_D

  LOGICAL :: Verbose

  REAL(AR), PARAMETER :: Zero = 0.0_AR
  REAL(AR), PARAMETER :: One  = 1.0_AR

  INTERFACE
    SUBROUTINE MF_Euler_Increment &
      ( GEOM, MF_uGF, MF_uCF, MF_uDF, MF_duCF )
      USE amrex_base_module, ONLY: &
        amrex_geometry, &
        amrex_multifab
      USE MyAmrModule, ONLY: &
        nLevels
      TYPE(amrex_geometry), INTENT(in   ) :: GEOM   (0:nLevels-1)
      TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF (0:nLevels-1)
      TYPE(amrex_multifab), INTENT(in   ) :: MF_uCF (0:nLevels-1)
      TYPE(amrex_multifab), INTENT(in   ) :: MF_uDF (0:nLevels-1)
      TYPE(amrex_multifab), INTENT(inout) :: MF_duCF(0:nLevels-1)
    END SUBROUTINE MF_Euler_Increment
  END INTERFACE


CONTAINS


  SUBROUTINE MF_InitializeFluid_SSPRK( nStages, BA, DM, Verbose_Option )

    INTEGER,               INTENT(in)           :: nStages
    TYPE(amrex_boxarray),  INTENT(in)           :: BA(0:nLevels-1)
    TYPE(amrex_distromap), INTENT(in)           :: DM(0:nLevels-1)
    LOGICAL,               INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER         :: iS, iLevel
    TYPE(amrex_box) :: BX

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    nStages_SSPRK = nStages

    CALL InitializeSSPRK( nStages )

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A,I1)') '', 'SSP RK Scheme: ', nStages

      WRITE(*,*)
      WRITE(*,'(A5,A)') '', 'Butcher Table:'
      WRITE(*,'(A5,A)') '', '--------------'
      DO iS = 1, nStages
        WRITE(*,'(A5,4ES14.4E3)') '', c_SSPRK(iS), a_SSPRK(iS,1:nStages)
      END DO
      WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_SSPRK(1:nStages)
      WRITE(*,*)

    END IF

    ALLOCATE( MF_U(0:nLevels-1) )
    ALLOCATE( MF_D(0:nLevels-1,1:nStages) )

    BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
        ( MF_U(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX(1) )

      DO iS = 1, nStages

        CALL amrex_multifab_build &
               ( MF_D(iLevel,iS), BA(iLevel), DM(iLevel), nDOFX * nCF, 0 )

      END DO

    END DO

  END SUBROUTINE MF_InitializeFluid_SSPRK


  SUBROUTINE MF_FinalizeFluid_SSPRK

    INTEGER :: iLevel, iS

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_U(iLevel) )

      DO iS = 1, nStages_SSPRK

        CALL amrex_multifab_destroy( MF_D(iLevel,iS) )

      END DO

    END DO

    DEALLOCATE( MF_U )
    DEALLOCATE( MF_D )

  END SUBROUTINE MF_FinalizeFluid_SSPRK


  SUBROUTINE InitializeSSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    INTEGER :: iS

    CALL AllocateButcherTables_SSPRK( nStages )

    SELECT CASE ( nStages )

      CASE ( 1 )

        a_SSPRK(1,1) = 0.0_AR
        w_SSPRK(1)   = 1.0_AR

      CASE ( 2 )

        a_SSPRK(1,1:2) = [ 0.0_AR, 0.0_AR ]
        a_SSPRK(2,1:2) = [ 1.0_AR, 0.0_AR ]
        w_SSPRK(1:2)   = [ 0.5_AR, 0.5_AR ]

      CASE ( 3 )

        a_SSPRK(1,1:3) = [ 0.00_AR, 0.00_AR, 0.00_AR ]
        a_SSPRK(2,1:3) = [ 1.00_AR, 0.00_AR, 0.00_AR ]
        a_SSPRK(3,1:3) = [ 0.25_AR, 0.25_AR, 0.00_AR ]
        w_SSPRK(1:3)   = [ 1.0_AR / 6.0_AR, &
                           1.0_AR / 6.0_AR, &
                           2.0_AR / 3.0_AR ]

    END SELECT

    DO iS = 1, nStages

      c_SSPRK(iS) = SUM( a_SSPRK(iS,1:iS-1) )

    END DO

  END SUBROUTINE InitializeSSPRK


  SUBROUTINE AllocateButcherTables_SSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    ALLOCATE( a_SSPRK(nStages,nStages) )
    ALLOCATE( c_SSPRK(nStages) )
    ALLOCATE( w_SSPRK(nStages) )

    a_SSPRK = Zero
    c_SSPRK = Zero
    w_SSPRK = Zero

  END SUBROUTINE AllocateButcherTables_SSPRK


  SUBROUTINE MF_UpdateFluid_SSPRK &
    ( t, dt, MF_uGF, MF_uCF, MF_uDF, GEOM, MF_Euler_ComputeIncrement )

    REAL(AR),             INTENT(in   ) :: t(0:nLevels-1), dt(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in   ) :: GEOM  (0:nLevels-1)
    PROCEDURE(MF_Euler_Increment)       :: MF_Euler_ComputeIncrement

    INTEGER :: iS, jS

    INTEGER                       :: iLevel
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:), U(:,:,:,:)

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_UpdateFluid )

    ! --- Set temporary MultiFabs U and dU to zero ---

    DO iLevel = 0, nLevels-1

      CALL MF_U(iLevel) % setval( Zero )

      DO iS = 1, nStages_SSPRK

        CALL MF_D(iLevel,iS) % setval( Zero )

      END DO

    END DO

    DO iS = 1, nStages_SSPRK

      ! --- Copy data from input MultiFab to temporary MultiFab ---

      DO iLevel = 0, nLevels-1

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_CopyMultiFab )

        CALL MF_U(iLevel) &
               % COPY( MF_uCF(iLevel), 1, 1, &
                       MF_uCF(iLevel) % nComp(), swX(1) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_CopyMultiFab )

        ! --- Apply boundary conditions to interior domains ---

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

        CALL MF_U(iLevel) % Fill_Boundary( GEOM(iLevel) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

        ! --- Copy ghost data from physical boundaries ---

        CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

        DO WHILE( MFI % next() )

          uCF => MF_uCF(iLevel) % DataPtr( MFI )
          U   => MF_U  (iLevel) % DataPtr( MFI )
          U   =  uCF

        END DO

        CALL amrex_mfiter_destroy( MFI )

      END DO

      DO iLevel = 0, nLevels-1

        DO jS = 1, iS-1

          IF( a_SSPRK(iS,jS) .NE. Zero ) &
            CALL MF_U(iLevel) &
                   % LinComb( One, MF_U(iLevel), 1, &
                              dt(iLevel) * a_SSPRK(iS,jS), MF_D(iLevel,jS), 1, &
                              1, MF_U(iLevel) % nComp(), 0 )

        END DO

      END DO

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        IF( DEBUG ) WRITE(*,'(A)') '  CALL MF_ApplySlopeLimiter_Euler (1)'

        CALL MF_ApplySlopeLimiter_Euler( MF_uGF, MF_U, MF_uDF, GEOM )

        IF( DEBUG ) WRITE(*,'(A)') '  CALL MF_ApplyPositivityLimiter_Euler (1)'

        CALL MF_ApplyPositivityLimiter_Euler( MF_uGF, MF_U )

        IF( DEBUG ) WRITE(*,'(A)') '  CALL MF_Euler_ComputeIncrement'

        CALL MF_Euler_ComputeIncrement( GEOM, MF_uGF, MF_U, MF_uDF, MF_D(:,iS) )

      END IF

    END DO

    DO iLevel = 0, nLevels-1

      DO iS = 1, nStages_SSPRK

        IF( w_SSPRK(iS) .NE. Zero ) &
          CALL MF_uCF(iLevel) &
                 % LinComb( One, MF_uCF(iLevel), 1, &
                            dt(iLevel) * w_SSPRK(iS), MF_D(iLevel,iS), 1, &
                            1, MF_uCF(iLevel) % nComp(), 0 )

      END DO

    END DO

    IF( DEBUG ) WRITE(*,'(A)') '  CALL MF_ApplySlopeLimiter_Euler (2)'

    CALL MF_ApplySlopeLimiter_Euler( MF_uGF, MF_uCF, MF_uDF, GEOM )

    IF( DEBUG ) WRITE(*,'(A)') '  CALL MF_ApplyPositivityLimiter_Euler (2)'

    CALL MF_ApplyPositivityLimiter_Euler( MF_uGF, MF_uCF )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_UpdateFluid )

  END SUBROUTINE MF_UpdateFluid_SSPRK


END MODULE MF_TimeSteppingModule_SSPRK
