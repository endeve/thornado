MODULE MF_TimeSteppingModule_SSPRK_MHD

  ! --- AMReX Modules ---

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
    nDOFX, &
    nNodes, &
    nDimsX, &
    swX
  USE MagnetofluidFieldsModule, ONLY: &
    nCM
  USE XCFC_UtilitiesModule, ONLY: &
    nGS, &
    nMF, &
    swX_GS

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two
  USE InputParsingModule, ONLY: &
    nLevels, &
    nMaxLevels, &
    UseFluxCorrection_MHD, &
    dt, &
    DEBUG, &
    t_new
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_MHD, ONLY: &
    MF_uCM, &
    MF_uDM, &
    OffGridFlux_MHD_MF
  USE MF_MHD_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_MHD_MF
  USE MF_MHD_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_MHD_MF
  USE MF_MHD_dgDiscretizationModule, ONLY: &
    ComputeIncrement_MHD_MF
  USE MF_MHD_TallyModule, ONLY: &
    IncrementOffGridTally_MHD_MF
  USE MF_XCFC_UtilitiesModule, ONLY: &
    MultiplyWithPsi6_MF, &
    UpdateConformalFactorAndMetric_XCFC_MF, &
    UpdateLapseShiftCurvature_XCFC_MF, &
    ComputeConformalFactorSourcesAndMg_XCFC_MF, &
    ComputePressureTensorTrace_XCFC_MF
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF
  USE MF_GravitySolutionModule, ONLY: &
    EvolveGravity
  USE MF_GravitySolutionModule_XCFC, ONLY: &
    ComputeConformalFactor_XCFC_MF, &
    ComputeLapseShiftCurvature_XCFC_MF
  USE AverageDownModule_MHD, ONLY: &
    AverageDown
  USE FluxCorrectionModule_MHD, ONLY: &
    ApplyFluxCorrection_MHD_MF
  USE MF_TimersModule_MHD, ONLY: &
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

  INTEGER  :: nStages
  REAL(DP), PUBLIC :: CFL

CONTAINS


  SUBROUTINE InitializeFluid_SSPRK_MF( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iS

    TYPE(amrex_parmparse) :: PP

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL amrex_parmparse_build( PP, 'TS' )
      CALL PP % get  ( 'nStages', nStages )
      CALL PP % get  ( 'CFL'    , CFL )
    CALL amrex_parmparse_destroy( PP )

    CFL = CFL / ( DBLE( nDimsX ) * ( Two * DBLE( nNodes ) - One ) )

    CALL InitializeSSPRK

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
        '    INFO: Timestepper'
      WRITE(*,'(A)') &
        '    -----------------'
      WRITE(*,'(A5,A,I1)') '', 'SSP RK Scheme: ', nStages
      WRITE(*,'(A5,A,ES10.3E3)') '', 'CFL:           ', &
        CFL * ( DBLE( nDimsX ) * ( Two * DBLE( nNodes ) - One ) )
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
    TYPE(amrex_multifab) :: MF_uMF(        0:nMaxLevels-1)

    INTEGER :: iS, jS, nCompCF
    INTEGER :: iLevel, iErr

    REAL(DP) :: dM_OffGrid_MHD(1:nCM,0:nMaxLevels-1)

    CALL TimersStart_AMReX( Timer_AMReX_UpdateFluid )

    dM_OffGrid_MHD = Zero

    nCompCF = nDOFX * nCM

    IF( EvolveGravity )THEN

      DO iLevel = 0, nLevels-1

        CALL amrex_multifab_build &
               ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
                 nDOFX * nGS, 0 )
        CALL MF_uGS(iLevel) % SetVal( Zero ) ! remove this after debugging

        CALL amrex_multifab_build &
               ( MF_uMF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
                 nDOFX * nMF, swX_GS )
        CALL MF_uMF(iLevel) % SetVal( Zero ) ! remove this after debugging

      END DO

    END IF ! EvolveGravity

    CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCM, +1 )

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
               ( MF_U(iS,iLevel), MF_uCM(iLevel) % BA, &
                 MF_uCM(iLevel) % DM, nCompCF, swX )

        CALL MF_U(iS,iLevel) % COPY( MF_uCM(iLevel), 1, 1, nCompCF, swX )

        CALL amrex_multifab_build &
               ( MF_D(iS,iLevel), MF_uCM(iLevel) % BA, &
                 MF_uCM(iLevel) % DM, nCompCF, swX )

      END DO ! iLevel = 0, nLevels-1

      DO jS = 1, iS-1

        IF( a_SSPRK(iS,jS) .NE. Zero )THEN

          DO iLevel = 0, nLevels-1

              CALL MF_U(iS,iLevel) &
                     % LinComb( One, MF_U(iS,iLevel), 1, &
                                dt(iLevel) * a_SSPRK(iS,jS), &
                                MF_D(jS,iLevel), 1, &
                                1, nCompCF, 0 )

          END DO

        END IF ! a_SSPRK(iS,jS) .NE. Zero

      END DO ! jS = 1, iS-1

      IF( iS .GT. 1 )THEN

        CALL MultiplyWithPsi6_MF( MF_uGF, MF_U(iS,:), -1 )

        CALL AverageDown( MF_uGF, MF_U(iS,:) )
        CALL ApplyPositivityLimiter_MHD_MF &
               ( t_new, MF_uGF, MF_U(iS,:), MF_uDM )

        CALL MultiplyWithPsi6_MF( MF_uGF, MF_U(iS,:), +1 )

      END IF

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        IF( iS .NE. 1 )THEN

          IF( EvolveGravity )THEN

            CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

            CALL ComputeConformalFactor( MF_uGF, MF_U(iS,:), MF_uGS, MF_uMF )

            CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

          END IF ! EvolveGravity

          CALL MultiplyWithPsi6_MF( MF_uGF, MF_U(iS,:), -1 )

          CALL ApplySlopeLimiter_MHD_MF &
                 ( t_new, MF_uGF, MF_U(iS,:), MF_uDM )

          CALL ApplyPositivityLimiter_MHD_MF &
                 ( t_new, MF_uGF, MF_U(iS,:), MF_uDM )

          CALL MultiplyWithPsi6_MF( MF_uGF, MF_U(iS,:), +1 )

          IF( EvolveGravity )THEN

            CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

            CALL ComputeConformalFactor &
                   ( MF_uGF, MF_U(iS,:), MF_uGS, MF_uMF )

            CALL ComputeLapseShiftCurvature &
                   ( MF_uGF, MF_U(iS,:), MF_uGS, MF_uMF )

            CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

          END IF ! EvolveGravity

        END IF ! iS .NE. 1

        CALL MultiplyWithPsi6_MF( MF_uGF, MF_U(iS,:), -1 )

        ! Come in with U, leave with \psi^6 * dU
        CALL ComputeIncrement_MHD_MF &
               ( MF_uGF, MF_U(iS,:), MF_uDM, MF_D(iS,:) )

        DO iLevel = 0, nLevels-1

          dM_OffGrid_MHD(:,iLevel) &
            = dM_OffGrid_MHD(:,iLevel) &
                + dt(iLevel) * w_SSPRK(iS) * OffGridFlux_MHD_MF(:,iLevel)

        END DO

        IF( nLevels .GT. 1 .AND. UseFluxCorrection_MHD )THEN

          CALL MultiplyWithPsi6_MF( MF_uGF, MF_D(iS,:), -1 )

          CALL ApplyFluxCorrection_MHD_MF( MF_uGF, MF_D(iS,:) )

          CALL MultiplyWithPsi6_MF( MF_uGF, MF_D(iS,:), +1 )

        END IF

      END IF ! a(:,iS) .NE. Zero .OR. w(iS) .NE. Zero

    END DO ! iS = 1, nStages

    DO iS = 1, nStages

      IF( w_SSPRK(iS) .NE. Zero )THEN

        DO iLevel = 0, nLevels-1

          CALL MF_uCM(iLevel) &
                 % LinComb( One, MF_uCM(iLevel), 1, &
                            dt(iLevel) * w_SSPRK(iS), MF_D(iS,iLevel), 1, 1, &
                            nCompCF, 0 )

        END DO

      END IF ! w_SSPRK(iS) .NE. Zero

    END DO ! iS = 1, nStages

    CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCM, -1 )

    CALL AverageDown( MF_uGF, MF_uCM )
    CALL ApplyPositivityLimiter_MHD_MF &
           ( t_new, MF_uGF, MF_uCM, MF_uDM )

    CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCM, +1 )

    DO iLevel = 0, nLevels-1

      DO iS = 1, nStages

        CALL amrex_multifab_destroy( MF_U(iS,iLevel) )
        CALL amrex_multifab_destroy( MF_D(iS,iLevel) )

      END DO

    END DO

    IF( EvolveGravity )THEN

      CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

      CALL ComputeConformalFactor( MF_uGF, MF_uCM, MF_uGS, MF_uMF )

      CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

    END IF ! EvolveGravity

    CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCM, -1 )

    CALL ApplySlopeLimiter_MHD_MF &
           ( t_new, MF_uGF, MF_uCM, MF_uDM )

    CALL ApplyPositivityLimiter_MHD_MF &
           ( t_new, MF_uGF, MF_uCM, MF_uDM )

    CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCM, +1 )

    IF( EvolveGravity )THEN

      CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

      CALL ComputeConformalFactor( MF_uGF, MF_uCM, MF_uGS, MF_uMF )

      CALL ComputeLapseShiftCurvature( MF_uGF, MF_uCM, MF_uGS, MF_uMF )

      DO iLevel = 0, nLevels-1

        CALL amrex_multifab_destroy( MF_uMF(iLevel) )
        CALL amrex_multifab_destroy( MF_uGS(iLevel) )

      END DO

      CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

    END IF ! EvolveGravity

    CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCM, -1 )

    CALL IncrementOffGridTally_MHD_MF( dM_OffGrid_MHD )

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


  SUBROUTINE ComputeConformalFactor( MF_uGF, MF_uCM, MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCM, MF_uGS )

    CALL ComputeConformalFactor_XCFC_MF &
           ( MF_uGS, MF_uMF )

    CALL UpdateConformalFactorAndMetric_XCFC_MF &
           ( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF, UpdateSpatialMetric_Option = .TRUE. )

  END SUBROUTINE ComputeConformalFactor


  SUBROUTINE ComputeLapseShiftCurvature( MF_uGF, MF_uCM, MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

    CALL ComputePressureTensorTrace_XCFC_MF &
           ( MF_uGF, MF_uCM, MF_uGS )

    CALL ComputeLapseShiftCurvature_XCFC_MF &
           ( MF_uGS, MF_uMF )

    CALL UpdateLapseShiftCurvature_XCFC_MF &
           ( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF, UpdateSpatialMetric_Option = .TRUE. )

    CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )

  END SUBROUTINE ComputeLapseShiftCurvature


END MODULE MF_TimeSteppingModule_SSPRK_MHD
