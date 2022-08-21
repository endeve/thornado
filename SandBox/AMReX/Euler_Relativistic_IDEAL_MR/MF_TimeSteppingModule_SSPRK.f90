MODULE MF_TimeSteppingModule_SSPRK

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator, &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE FluidFieldsModule, ONLY: &
    nCF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two
  USE MF_Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_MF
  USE MF_Euler_TallyModule, ONLY: &
    IncrementOffGridTally_Euler_MF
  USE MF_Euler_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_Euler_MF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler_MF
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uDF, &
    OffGridFlux_Euler_MF
  USE AverageDownModule, ONLY: &
    AverageDown
  USE InputParsingModule, ONLY: &
    nLevels, &
    nMaxLevels, &
    swX, &
    CFL, &
    nNodes, &
    nStages, &
    UseFluxCorrection_Euler, &
    t_new, &
    dt, &
    DEBUG
  USE FluxCorrectionModule_Euler, ONLY: &
    ApplyFluxCorrection_Euler_MF
  USE MF_GravitySolutionModule_XCFC_Poseidon, ONLY: &
    nGS, &
    MultiplyWithPsi6_MF, &
    ComputeConformalFactorSourcesAndMg_XCFC_MF, &
    ComputeConformalFactor_Poseidon_MF, &
    ComputePressureTensorTrace_XCFC_MF, &
    ComputeGeometry_Poseidon_MF
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

CONTAINS


  SUBROUTINE InitializeFluid_SSPRK_MF( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iS

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL InitializeSSPRK

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
        '    INFO: Timestepper'
      WRITE(*,'(A)') &
        '    -----------------'
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

  END SUBROUTINE InitializeFluid_SSPRK_MF


  SUBROUTINE FinalizeFluid_SSPRK_MF

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

  END SUBROUTINE FinalizeFluid_SSPRK_MF


  SUBROUTINE UpdateFluid_SSPRK_MF

    TYPE(amrex_multifab) :: MF_U(1:nStages,0:nMaxLevels-1)
    TYPE(amrex_multifab) :: MF_D(1:nStages,0:nMaxLevels-1)
    TYPE(amrex_multifab) :: MF_uGS(        0:nMaxLevels-1)

    INTEGER :: iS, jS, nCompCF
    INTEGER :: iLevel, iErr

    REAL(DP) :: dM_OffGrid_Euler(1:nCF,0:nMaxLevels-1)

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_UpdateFluid )

    dM_OffGrid_Euler = Zero

    nCompCF = nDOFX * nCF

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nGS, 0 )
      CALL MF_uGS(iLevel) % SetVal( Zero ) ! remove this after debugging

    END DO

    CALL MultiplyWithPsi6_MF( MF_uGF, +1, MF_uCF )

    DO iS = 1, nStages

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,*)
          WRITE(*,'(2x,A,I2.2)') 'iS: ', iS
          WRITE(*,*)

        END IF

      END IF ! DEBUG

      DO iLevel = 0, nLevels-1

        CALL amrex_multifab_build &
               ( MF_U(iS,iLevel), MF_uCF(iLevel) % BA, &
                 MF_uCF(iLevel) % DM, nCompCF, swX )

        CALL amrex_multifab_build &
               ( MF_D(iS,iLevel), MF_uCF(iLevel) % BA, &
                 MF_uCF(iLevel) % DM, nCompCF, swX )

        CALL MF_U(iS,iLevel) % COPY( MF_uCF(iLevel), 1, 1, nCompCF, swX )

      END DO ! iLevel

      DO jS = 1, iS-1

        DO iLevel = 0, nLevels-1

          IF( a_SSPRK(iS,jS) .NE. Zero )THEN

            CALL MF_U(iS,iLevel) &
                   % LinComb( One, MF_U(iS,iLevel), 1, &
                              dt(iLevel) * a_SSPRK(iS,jS), MF_D(jS,iLevel), 1, &
                              1, nCompCF, 0 )

          END IF

        END DO ! iLevel

!!$        IF( a_SSPRK(iS,jS) .NE. Zero ) &
!!$          CALL AverageDown( MF_uGF, MF_U(iS,:) )

      END DO ! jS

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        IF( iS .NE. 1 )THEN

          CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
                 ( MF_uGF, MF_U(iS,:), MF_uGS )

          CALL ComputeConformalFactor_Poseidon_MF( MF_uGS, MF_uGF )

          CALL MultiplyWithPsi6_MF( MF_uGF, -1, MF_U(iS,:) )

          CALL ApplySlopeLimiter_Euler_MF &
                 ( t_new, MF_uGF, MF_U(iS,:), MF_uDF )

          CALL ApplyPositivityLimiter_Euler_MF &
                 ( MF_uGF, MF_U(iS,:), MF_uDF )

          CALL MultiplyWithPsi6_MF( MF_uGF, +1, MF_U(iS,:) )

          CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
                 ( MF_uGF, MF_U(iS,:), MF_uGS )

          CALL ComputeConformalFactor_Poseidon_MF( MF_uGS, MF_uGF )

          CALL ComputePressureTensorTrace_XCFC_MF( MF_uGF, MF_U(iS,:), MF_uGS )

          CALL ComputeGeometry_Poseidon_MF( MF_uGS, MF_uGF )

        END IF

        CALL MultiplyWithPsi6_MF( MF_uGF, -1, MF_U(iS,:) )

        CALL ComputeIncrement_Euler_MF &
               ( t_new, MF_uGF, MF_U(iS,:), MF_uDF, MF_D(iS,:) )

        DO iLevel = 0, nLevels-1

          dM_OffGrid_Euler(:,iLevel) &
            = dM_OffGrid_Euler(:,iLevel) &
                + dt(iLevel) * w_SSPRK(iS) * OffGridFlux_Euler_MF(:,iLevel)

        END DO

        IF( nLevels .GT. 1 .AND. UseFluxCorrection_Euler ) &
          CALL ApplyFluxCorrection_Euler_MF( MF_uGF, MF_D(iS,:) )

      END IF ! a(:,iS) .NE. Zero .OR. w(iS) .NE. Zero

    END DO ! iS

    DO iS = 1, nStages

      DO iLevel = 0, nLevels-1

        IF( w_SSPRK(iS) .NE. Zero )THEN

          CALL MF_uCF(iLevel) &
                 % LinComb( One, MF_uCF(iLevel), 1, &
                            dt(iLevel) * w_SSPRK(iS), MF_D(iS,iLevel), 1, 1, &
                            nCompCF, 0 )

        END IF

      END DO ! iLevel

!!$      IF( w_SSPRK(iS) .NE. Zero ) &
!!$        CALL AverageDown( MF_uGF, MF_uCF )

    END DO ! iS

    DO iLevel = 0, nLevels-1

      DO iS = 1, nStages

        CALL amrex_multifab_destroy( MF_U(iS,iLevel) )
        CALL amrex_multifab_destroy( MF_D(iS,iLevel) )

      END DO ! iS

    END DO ! iLevel

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uGS )

    CALL ComputeConformalFactor_Poseidon_MF( MF_uGS, MF_uGF )

    CALL MultiplyWithPsi6_MF( MF_uGF, -1, MF_uCF )

    CALL ApplySlopeLimiter_Euler_MF &
           ( t_new, MF_uGF, MF_uCF, MF_uDF )

    CALL ApplyPositivityLimiter_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uDF )

    CALL MultiplyWithPsi6_MF( MF_uGF, +1, MF_uCF )

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uGS )

    CALL ComputeConformalFactor_Poseidon_MF( MF_uGS, MF_uGF )

    CALL ComputePressureTensorTrace_XCFC_MF( MF_uGF, MF_uCF, MF_uGS )

    CALL ComputeGeometry_Poseidon_MF( MF_uGS, MF_uGF )

    CALL MultiplyWithPsi6_MF( MF_uGF, -1, MF_uCF )

    CALL IncrementOffGridTally_Euler_MF( dM_OffGrid_Euler )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_uGS(iLevel) )

    END DO

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
