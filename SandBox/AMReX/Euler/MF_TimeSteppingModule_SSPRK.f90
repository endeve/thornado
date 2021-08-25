MODULE MF_TimeSteppingModule_SSPRK

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_amrcore_module, ONLY: &
    amrex_max_level, &
    amrex_geom
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE FluidFieldsModule, ONLY: &
    nCF
  USE GeometryFieldsModule, ONLY: &
    nGF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two
  USE MF_Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_MF
!  USE MF_Euler_SlopeLimiterModule, ONLY: &
!    MF_ApplySlopeLimiter_Euler
!  USE MF_Euler_PositivityLimiterModule, ONLY: &
!    MF_ApplyPositivityLimiter_Euler
  USE InputParsingModule, ONLY: &
    UseTiling, &
    swX, &
    nX, &
    CFL, &
    nNodes, &
    nStages
  USE FillPatchModule, ONLY: &
    FillPatch_uCF
  USE MF_FieldsModule, ONLY: &
    MF_uGF_new, &
    MF_OffGridFlux_Euler
!  USE MF_Euler_TallyModule, ONLY: &
!    MF_IncrementOffGridTally_Euler
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_UpdateFluid, &
    Timer_AMReX_Euler_InteriorBC, &
    Timer_AMReX_Euler_CopyMultiFab

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFluid_SSPRK_MF
  PUBLIC :: UpdateFluid_SSPRK_MF
  PUBLIC :: FinalizeFluid_SSPRK_MF

  REAL(DP), DIMENSION(:),   ALLOCATABLE :: c_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: w_SSPRK
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  TYPE(amrex_multifab), ALLOCATABLE :: MF_U(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_D(:,:)

  LOGICAL :: Verbose


CONTAINS


  SUBROUTINE InitializeFluid_SSPRK_MF( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iS, iLevel

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL InitializeSSPRK

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A,I1)') '', 'SSP RK Scheme: ', nStages
      WRITE(*,'(A5,A,ES10.3E3)') '', 'CFL:           ', &
        CFL * ( DBLE( amrex_spacedim ) * ( Two * DBLE( nNodes ) - One ) )

      WRITE(*,*)
      WRITE(*,'(A5,A)') '', 'Butcher Table:'
      WRITE(*,'(A5,A)') '', '--------------'
      DO iS = 1, nStages
        WRITE(*,'(A5,4ES14.4E3)') '', c_SSPRK(iS), a_SSPRK(iS,1:nStages)
      END DO
      WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_SSPRK(1:nStages)
      WRITE(*,*)

    END IF

    ALLOCATE( MF_U(0:amrex_max_level) )
    ALLOCATE( MF_D(0:amrex_max_level,1:nStages) )

    DO iLevel = 0, amrex_max_level

      CALL amrex_multifab_build &
             ( MF_U(iLevel), MF_uGF_new(iLevel) % BA, &
                             MF_uGF_new(iLevel) % DM, nDOFX * nCF, swX )

      DO iS = 1, nStages

        CALL amrex_multifab_build &
               ( MF_D(iLevel,iS), MF_uGF_new(iLevel) % BA, &
                                  MF_uGF_new(iLevel) % DM, nDOFX * nCF, swX )

      END DO

    END DO

  END SUBROUTINE InitializeFluid_SSPRK_MF


  SUBROUTINE FinalizeFluid_SSPRK_MF

    INTEGER :: iLevel, iS

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

    DO iLevel = 0, amrex_max_level

      CALL amrex_multifab_destroy( MF_U(iLevel) )

      DO iS = 1, nStages

        CALL amrex_multifab_destroy( MF_D(iLevel,iS) )

      END DO

    END DO

    DEALLOCATE( MF_U )
    DEALLOCATE( MF_D )

  END SUBROUTINE FinalizeFluid_SSPRK_MF


  SUBROUTINE UpdateFluid_SSPRK_MF &
    ( t, dt, MF_uGF, MF_uCF, MF_uDF )

    REAL(DP),             INTENT(in)    :: t (0:amrex_max_level), &
                                           dt(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:amrex_max_level)

    INTEGER :: iS, jS

    INTEGER                       :: iLevel
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:), U(:,:,:,:)

    REAL(DP) :: dM_OffGrid_Euler(0:amrex_max_level,nCF)

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_UpdateFluid )

    dM_OffGrid_Euler = Zero

    ! --- Set temporary MultiFabs U and dU to zero ---

    DO iLevel = 0, amrex_max_level

      CALL MF_U(iLevel) % SetVal( Zero )

      DO iS = 1, nStages

        CALL MF_D(iLevel,iS) % SetVal( Zero )

      END DO

    END DO

    DO iS = 1, nStages

      ! --- Copy data from input MultiFab to temporary MultiFab ---

      DO iLevel = 0, amrex_max_level

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_CopyMultiFab )

        CALL MF_U(iLevel) &
               % COPY( MF_uCF(iLevel), 1, 1, &
                       MF_uCF(iLevel) % nComp(), swX )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_CopyMultiFab )

        ! --- Apply boundary conditions to interior domains ---

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

        CALL MF_U(iLevel) % Fill_Boundary( amrex_geom(iLevel) )

        CALL FillPatch_uCF( iLevel, t(iLevel), MF_U(iLevel) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

        ! --- Copy ghost data from physical boundaries ---

        CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

        DO WHILE( MFI % next() )

          uCF => MF_uCF(iLevel) % DataPtr( MFI )
          U   => MF_U  (iLevel) % DataPtr( MFI )
          U   =  uCF

        END DO

        CALL amrex_mfiter_destroy( MFI )

      END DO

      DO iLevel = 0, amrex_max_level

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

!        CALL MF_ApplySlopeLimiter_Euler( MF_uGF, MF_U, MF_uDF, GEOM )

!        CALL MF_ApplyPositivityLimiter_Euler( MF_uGF, MF_U, MF_uDF )

        CALL ComputeIncrement_Euler_MF( t, MF_uGF, MF_U, MF_uDF, MF_D(:,iS) )

!        DO iLevel = 0, amrex_max_level
!
!          dM_OffGrid_Euler(iLevel,:) &
!            = dM_OffGrid_Euler(iLevel,:) &
!                + dt(iLevel) * w_SSPRK(iS) * MF_OffGridFlux_Euler(iLevel,:)
!
!        END DO

      END IF

    END DO

    DO iLevel = 0, amrex_max_level

      DO iS = 1, nStages

        IF( w_SSPRK(iS) .NE. Zero ) &
          CALL MF_uCF(iLevel) &
                 % LinComb( One, MF_uCF(iLevel), 1, &
                            dt(iLevel) * w_SSPRK(iS), MF_D(iLevel,iS), 1, &
                            1, MF_uCF(iLevel) % nComp(), swX )

      END DO

    END DO

!    CALL MF_ApplySlopeLimiter_Euler( MF_uGF, MF_uCF, MF_uDF, GEOM )

!    CALL MF_ApplyPositivityLimiter_Euler( MF_uGF, MF_uCF, MF_uDF )

!    CALL MF_IncrementOffGridTally_Euler( dM_OffGrid_Euler )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_UpdateFluid )

  END SUBROUTINE UpdateFluid_SSPRK_MF


  ! --- PRIVATE Subroutines ---


  SUBROUTINE InitializeSSPRK

    INTEGER :: iS

    CALL AllocateButcherTables_SSPRK

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


  SUBROUTINE AllocateButcherTables_SSPRK

    ALLOCATE( a_SSPRK(nStages,nStages) )
    ALLOCATE( c_SSPRK(nStages) )
    ALLOCATE( w_SSPRK(nStages) )

    a_SSPRK = Zero
    c_SSPRK = Zero
    w_SSPRK = Zero

  END SUBROUTINE AllocateButcherTables_SSPRK


END MODULE MF_TimeSteppingModule_SSPRK
