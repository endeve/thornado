MODULE MF_TimeSteppingModule_SSPRK

  ! --- AMReX Modules ---

  USE amrex_box_module,                 ONLY: &
    amrex_box
  USE amrex_geometry_module,            ONLY: &
    amrex_geometry
  USE amrex_multifab_module,            ONLY: &
    amrex_multifab,         &
    amrex_multifab_build,   &
    amrex_multifab_destroy, &
    amrex_mfiter,           &
    amrex_mfiter_build,     &
    amrex_mfiter_destroy
  USE amrex_boxarray_module,            ONLY: &
    amrex_boxarray,       &
    amrex_boxarray_build, &
    amrex_boxarray_destroy
  USE amrex_distromap_module,           ONLY: &
    amrex_distromap,       &
    amrex_distromap_build, &
    amrex_distromap_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule,              ONLY: &
    swX,   &
    nDOFX, &
    nX
  USE MagnetofluidfieldsModule,                ONLY: &
    nCM
  USE GeometryFieldsModule,             ONLY: &
    nGF

  ! --- Local Modules ---

  USE MF_KindModule,                    ONLY: &
    DP, &
    Zero, &
    One
  USE InputParsingModule,               ONLY: &
    nLevels, &
    UseTiling, &
    DEBUG, &
    EvolveOnlyMagnetic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeMagnetofluid_SSPRK
  PUBLIC :: MF_UpdateMagnetofluid_SSPRK
  PUBLIC :: MF_FinalizeMagnetofluid_SSPRK

  INTEGER :: nStages_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: c_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: w_SSPRK
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  TYPE(amrex_multifab), DIMENSION(:),   ALLOCATABLE :: MF_U
  TYPE(amrex_multifab), DIMENSION(:,:), ALLOCATABLE :: MF_D

  LOGICAL :: Verbose

  INTERFACE
    SUBROUTINE MF_MHD_Increment &
      ( GEOM, MF_uGF, MF_uCM, MF_uDM, MF_duCM )
      USE amrex_geometry_module, ONLY: &
        amrex_geometry
      USE amrex_multifab_module, ONLY: &
        amrex_multifab
      USE InputParsingModule,    ONLY: &
        nLevels
      TYPE(amrex_geometry), INTENT(in)    :: GEOM   (0:nLevels-1)
      TYPE(amrex_multifab), INTENT(in)    :: MF_uGF (0:nLevels-1)
      TYPE(amrex_multifab), INTENT(in)    :: MF_uCM (0:nLevels-1)
      TYPE(amrex_multifab), INTENT(in)    :: MF_uDM (0:nLevels-1)
      TYPE(amrex_multifab), INTENT(inout) :: MF_duCM(0:nLevels-1)
    END SUBROUTINE MF_MHD_Increment
  END INTERFACE


CONTAINS


  SUBROUTINE MF_InitializeMagnetofluid_SSPRK( nStages, BA, DM, Verbose_Option )

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
        ( MF_U(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCM, swX )

      DO iS = 1, nStages

        CALL amrex_multifab_build &
               ( MF_D(iLevel,iS), BA(iLevel), DM(iLevel), nDOFX * nCM, swX )

      END DO

    END DO

  END SUBROUTINE MF_InitializeMagnetofluid_SSPRK


  SUBROUTINE MF_FinalizeMagnetofluid_SSPRK

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

  END SUBROUTINE MF_FinalizeMagnetofluid_SSPRK


  SUBROUTINE InitializeSSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    INTEGER :: iS

    CALL AllocateButcherTables_SSPRK( nStages )

    SELECT CASE ( nStages )

      CASE ( 1 )

        a_SSPRK(1,1) = 0.0_DP
        w_SSPRK(1)   = 1.0_DP

      CASE ( 2 )

        a_SSPRK(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_SSPRK(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_SSPRK(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 3 )

        a_SSPRK(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_SSPRK(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_SSPRK(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_SSPRK(1:3)   = [ 1.0_DP / 6.0_DP, &
                           1.0_DP / 6.0_DP, &
                           2.0_DP / 3.0_DP ]

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


  SUBROUTINE MF_UpdateMagnetofluid_SSPRK &
    ( t, dt, MF_uGF, MF_uCM, MF_uDM, GEOM, MF_ComputeIncrement_MHD )

    REAL(DP),             INTENT(in)    :: t(0:nLevels-1), dt(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDM(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)
    PROCEDURE(MF_MHD_Increment)       :: MF_ComputeIncrement_MHD

    INTEGER :: iS, jS

    INTEGER                       :: iLevel
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:), U(:,:,:,:)

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

        CALL MF_U(iLevel) &
               % COPY( MF_uCM(iLevel), 1, 1, &
                       MF_uCM(iLevel) % nComp(), swX )

        ! --- Apply boundary conditions to interior domains ---

        CALL MF_U(iLevel) % Fill_Boundary( GEOM(iLevel) )

        ! --- Copy ghost data from physical boundaries ---

        CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

        DO WHILE( MFI % next() )

          uCM => MF_uCM(iLevel) % DataPtr( MFI )
          U   => MF_U  (iLevel) % DataPtr( MFI )
          U   =  uCM

        END DO

        CALL amrex_mfiter_destroy( MFI )

      END DO

      DO iLevel = 0, nLevels-1

        DO jS = 1, iS-1

          IF( a_SSPRK(iS,jS) .NE. Zero ) &
            CALL MF_U(iLevel) &
                   % LinComb( One, MF_U(iLevel), 1, &
                              dt(iLevel) * a_SSPRK(iS,jS), MF_D(iLevel,jS), 1, &
                              1, MF_U(iLevel) % nComp(), swX )

        END DO

      END DO

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        IF( DEBUG ) WRITE(*,'(A)') '  CALL MF_ComputeIncrement_MHD'

        CALL MF_ComputeIncrement_MHD( GEOM, MF_uGF, MF_U, MF_uDM, MF_D(:,iS) )

      END IF

    END DO

    DO iLevel = 0, nLevels-1

      DO iS = 1, nStages_SSPRK

        IF( w_SSPRK(iS) .NE. Zero ) &
          CALL MF_uCM(iLevel) &
                 % LinComb( One, MF_uCM(iLevel), 1, &
                            dt(iLevel) * w_SSPRK(iS), MF_D(iLevel,iS), 1, &
                            1, MF_uCM(iLevel) % nComp(), swX )

      END DO

    END DO

  END SUBROUTINE MF_UpdateMagnetofluid_SSPRK


END MODULE MF_TimeSteppingModule_SSPRK
