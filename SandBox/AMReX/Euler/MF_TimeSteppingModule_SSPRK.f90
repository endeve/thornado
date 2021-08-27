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
    nStages, &
    StepNo
  USE FillPatchModule, ONLY: &
    FillPatch_uCF
  USE MF_FieldsModule, ONLY: &
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

  LOGICAL :: Verbose

  INTEGER, ALLOCATABLE :: nSubSteps(:)

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

    ALLOCATE( nSubSteps(0:amrex_max_level) )

    DO iLevel = 0, amrex_max_level

      nSubSteps(iLevel) = 2**iLevel

    END DO

  END SUBROUTINE InitializeFluid_SSPRK_MF


  SUBROUTINE FinalizeFluid_SSPRK_MF

    INTEGER :: iLevel, iS

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

    DEALLOCATE( nSubSteps )

  END SUBROUTINE FinalizeFluid_SSPRK_MF


  SUBROUTINE UpdateFluid_SSPRK_MF &
    ( t, dt, MF_uGF, MF_uCF, MF_uDF )

    REAL(DP),             INTENT(in)    :: t     (0:amrex_max_level)
    REAL(DP),             INTENT(in)    :: dt    (0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:amrex_max_level)

    TYPE(amrex_multifab) :: MF_U
    TYPE(amrex_multifab) :: MF_D(1:nStages)

    INTEGER :: iS, jS
    INTEGER :: iLevel, iSubStep

    REAL(DP) :: dM_OffGrid_Euler(0:amrex_max_level,nCF)

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_UpdateFluid )

    DO iLevel = 0, amrex_max_level

      DO iSubStep = 1, nSubSteps(iLevel)

        StepNo(iLevel) = StepNo(iLevel) + 1

        CALL amrex_multifab_build &
               ( MF_U, MF_uCF(iLevel) % BA, &
                       MF_uCF(iLevel) % DM, nDOFX * nCF, swX )

        CALL MF_U &
               % COPY( MF_uCF(iLevel), 1, 1, &
                       MF_uCF(iLevel) % nComp(), swX )

        DO iS = 1, nStages

          CALL amrex_multifab_build &
                 ( MF_D(iS), MF_uCF(iLevel) % BA, &
                             MF_uCF(iLevel) % DM, nDOFX * nCF, swX )

        END DO

        dM_OffGrid_Euler(iLevel,:) = Zero

        ! --- Initialize increment MultiFab to zero ---

        DO iS = 1, nStages

          CALL MF_D(iS) % SetVal( Zero )

        END DO

        DO iS = 1, nStages

          ! --- Copy data from input MultiFab to temporary MultiFab ---

          CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_CopyMultiFab )

          CALL MF_uCF(iLevel) &
                 % COPY( MF_U, 1, 1, MF_U % nComp(), swX )

          CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_CopyMultiFab )

          ! --- Apply boundary conditions to interior domains ---

          CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

          CALL FillPatch_uCF( iLevel, t(iLevel), MF_uCF(iLevel) )
          !CALL MF_U % Fill_Boundary( amrex_geom(iLevel) )

          CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

          DO jS = 1, iS-1

            IF( a_SSPRK(iS,jS) .NE. Zero ) &
              CALL MF_uCF(iLevel) &
                     % LinComb( One, MF_uCF(iLevel), 1, &
                                dt(iLevel) * a_SSPRK(iS,jS), &
                                MF_D(jS), 1, 1, &
                                MF_uCF(iLevel) % nComp(), swX )

          END DO

          IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
              .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

!!$            CALL MF_ApplySlopeLimiter_Euler &
!!$                   ( MF_uGF(iLevel), MF_U, MF_uDF(iLevel) )
!!$
!!$            CALL MF_ApplyPositivityLimiter_Euler &
!!$                   ( MF_uGF(iLevel), MF_U, MF_uDF(iLevel) )

            CALL ComputeIncrement_Euler_MF &
                   ( iLevel, t(iLevel), MF_uGF(iLevel), MF_uCF(iLevel), &
                     MF_uDF(iLevel), MF_D(iS) )

            dM_OffGrid_Euler(iLevel,:) &
              = dM_OffGrid_Euler(iLevel,:) &
                  + dt(iLevel) * w_SSPRK(iS) * MF_OffGridFlux_Euler(iLevel,:)

          END IF

        END DO ! iS

        DO iS = 1, nStages

          IF( w_SSPRK(iS) .NE. Zero ) &
            CALL MF_U &
                   % LinComb( One, MF_U, 1, &
                              dt(iLevel) * w_SSPRK(iS), MF_D(iS), 1, 1, &
                              MF_U % nComp(), swX )

        END DO

        CALL MF_uCF(iLevel) % COPY( MF_U, 1, 1, MF_U % nComp(), swX )

        DO iS = 1, nStages

          CALL amrex_multifab_destroy( MF_D(iS) )

        END DO

        CALL amrex_multifab_destroy( MF_U )

      END DO ! iSubStep

    END DO ! iLevel

!!$    CALL MF_ApplySlopeLimiter_Euler( MF_uGF, MF_uCF, MF_uDF )
!!$
!!$    CALL MF_ApplyPositivityLimiter_Euler( MF_uGF, MF_uCF, MF_uDF )
!!$
!!$    CALL MF_IncrementOffGridTally_Euler( dM_OffGrid_Euler )

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
