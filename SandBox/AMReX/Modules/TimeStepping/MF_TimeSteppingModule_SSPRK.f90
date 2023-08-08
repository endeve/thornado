MODULE MF_TimeSteppingModule_SSPRK

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
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
  USE XCFC_UtilitiesModule, ONLY: &
    nGS, &
    swXX, &
    MultiplyWithPsi6_MF
  USE MF_GravitySolutionModule_XCFC, ONLY: &
    ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF, &
    ComputeConformalFactor_MF, &
    ComputePressureTensorTrace_XCFC_Euler_MF, &
    ComputeGeometry_MF
  USE MF_TimersModule, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_UpdateFluid, &
    Timer_AMReX_GravitySolve

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFluid_SSPRK_MF
  PUBLIC :: UpdateFluid_SSPRK_MF
  PUBLIC :: FinalizeFluid_SSPRK_MF

  REAL(DP), DIMENSION(:)  , ALLOCATABLE :: c_SSPRK
  REAL(DP), DIMENSION(:)  , ALLOCATABLE :: w_SSPRK
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  LOGICAL :: Verbose
  LOGICAL :: EvolveGravity

CONTAINS


  SUBROUTINE InitializeFluid_SSPRK_MF( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iS

    TYPE(amrex_parmparse) :: PP

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    EvolveGravity = .FALSE.
    CALL amrex_parmparse_build( PP, 'TS' )
      CALL PP % query( 'EvolveGravity', EvolveGravity )
    CALL amrex_parmparse_destroy( PP )

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
      WRITE(*,'(A5,A,L)') '', 'EvolveGravity: ', EvolveGravity
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

    CALL TimersStart_AMReX( Timer_AMReX_UpdateFluid )

    dM_OffGrid_Euler = Zero

    nCompCF = nDOFX * nCF

    IF( EvolveGravity )THEN

      DO iLevel = 0, nLevels-1

        CALL amrex_multifab_build &
               ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
                 nDOFX * nGS, swXX )
        CALL MF_uGS(iLevel) % SetVal( Zero ) ! remove this after debugging

      END DO

    END IF ! EvolveGravity

    CALL MultiplyWithPsi6_MF( MF_uGF, +1, 1, 1, 1, 1, MF_uCF )

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

        CALL MF_U(iS,iLevel) % COPY( MF_uCF(iLevel), 1, 1, nCompCF, swX )

        CALL amrex_multifab_build &
               ( MF_D(iS,iLevel), MF_uCF(iLevel) % BA, &
                 MF_uCF(iLevel) % DM, nCompCF, swX )

      END DO ! iLevel = 0, nLevels-1

      DO jS = 1, iS-1

        DO iLevel = 0, nLevels-1

          IF( a_SSPRK(iS,jS) .NE. Zero ) &
            CALL MF_U(iS,iLevel) &
                   % LinComb( One, MF_U(iS,iLevel), 1, &
                              dt(iLevel) * a_SSPRK(iS,jS), MF_D(jS,iLevel), 1, &
                              1, nCompCF, 0 )

        END DO ! iLevel = 0, nLevels-1

!!$        IF( a_SSPRK(iS,jS) .NE. Zero ) &
!!$          CALL AverageDown &
!!$                 ( MF_uGF, MF_U(iS,:), &
!!$                   MF_uDF, ApplyPositivityLimiter_Option = .TRUE. )

      END DO ! jS = 1, iS-1

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        IF( iS .NE. 1 )THEN

          IF( EvolveGravity )THEN

            CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

            CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF &
                   ( MF_uGF, MF_U(iS,:), MF_uGS )

            CALL ComputeConformalFactor_MF( MF_uGS, MF_uGF )

            CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

          END IF ! EvolveGravity

          CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_U(iS,:) )

          CALL ApplySlopeLimiter_Euler_MF &
                 ( MF_uGF, MF_U(iS,:), MF_uDF )

          CALL ApplyPositivityLimiter_Euler_MF &
                 ( MF_uGF, MF_U(iS,:), MF_uDF )

          CALL MultiplyWithPsi6_MF( MF_uGF, +1, 1, 1, 1, 1, MF_U(iS,:) )

          IF( EvolveGravity )THEN

            CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

            CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF &
                   ( MF_uGF, MF_U(iS,:), MF_uGS )

            CALL ComputeConformalFactor_MF &
                   ( MF_uGS, MF_uGF )

            CALL ComputePressureTensorTrace_XCFC_Euler_MF &
                   ( MF_uGF, MF_U(iS,:), MF_uGS )

            CALL ComputeGeometry_MF &
                   ( MF_uGS, MF_uGF )

            CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

          END IF ! EvolveGravity

        END IF ! iS .NE. 1

        CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_U(iS,:) )

        ! Come in with U, leave with \psi^6 * dU
        CALL ComputeIncrement_Euler_MF &
               ( MF_uGF, MF_U(iS,:), MF_uDF, MF_D(iS,:) )

        DO iLevel = 0, nLevels-1

          dM_OffGrid_Euler(:,iLevel) &
            = dM_OffGrid_Euler(:,iLevel) &
                + dt(iLevel) * w_SSPRK(iS) * OffGridFlux_Euler_MF(:,iLevel)

        END DO

        IF( nLevels .GT. 1 .AND. UseFluxCorrection_Euler )THEN

          CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_D(iS,:) )

          CALL ApplyFluxCorrection_Euler_MF( MF_uGF, MF_D(iS,:) )

          CALL MultiplyWithPsi6_MF( MF_uGF, +1, 1, 1, 1, 1, MF_D(iS,:) )

        END IF

      END IF ! a(:,iS) .NE. Zero .OR. w(iS) .NE. Zero

    END DO ! iS = 1, nStages

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
!!$        CALL AverageDown &
!!$               ( MF_uGF, MF_uCF, &
!!$                 MF_uDF, ApplyPositivityLimiter_Option = .TRUE. )

    END DO ! iS

    DO iLevel = 0, nLevels-1

      DO iS = 1, nStages

        CALL amrex_multifab_destroy( MF_U(iS,iLevel) )
        CALL amrex_multifab_destroy( MF_D(iS,iLevel) )

      END DO

    END DO

    IF( EvolveGravity )THEN

      CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

      CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uGS )

      CALL ComputeConformalFactor_MF( MF_uGS, MF_uGF )

      CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

    END IF ! EvolveGravity

    CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_uCF )

    CALL ApplySlopeLimiter_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uDF )

    CALL ApplyPositivityLimiter_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uDF )

    CALL MultiplyWithPsi6_MF( MF_uGF, +1, 1, 1, 1, 1, MF_uCF )

    IF( EvolveGravity )THEN

      CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

      CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uGS )

      CALL ComputeConformalFactor_MF( MF_uGS, MF_uGF )

      CALL ComputePressureTensorTrace_XCFC_Euler_MF( MF_uGF, MF_uCF, MF_uGS )

      CALL ComputeGeometry_MF( MF_uGS, MF_uGF )

      DO iLevel = 0, nLevels-1

        CALL amrex_multifab_destroy( MF_uGS(iLevel) )

      END DO

      CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

    END IF ! EvolveGravity

    CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_uCF )

    CALL IncrementOffGridTally_Euler_MF( dM_OffGrid_Euler )

    CALL TimersStop_AMReX( Timer_AMReX_UpdateFluid )

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
