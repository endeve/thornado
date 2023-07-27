MODULE MF_TimeSteppingModule_IMEX

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_amrcore_module, ONLY: &
    GEOM => amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ, &
    iE_B0, &
    iE_E0, &
    swX
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nCR
  USE MF_TwoMoment_DiscretizationModule_Streaming_Relativistic, ONLY: &
    ComputeIncrement_TwoMoment_Explicit_MF
  USE MF_TwoMoment_DiscretizationModule_Collisions_Neutrinos_GR, ONLY: &
    ComputeIncrement_TwoMoment_Implicit_Neutrinos_MF
  USE MF_FieldsModule_TwoMoment, ONLY: &
    OffGridFlux_TwoMoment_MF


  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uDF
  USE MF_FieldsModule_TwoMoment, ONLY: &
    MF_uCR
  USE InputParsingModule, ONLY: &
    t_new, &
    dt, &
    nLevels, &
    nMaxLevels, &
!    DEBUG, &
    nE, &
    nSpecies
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler_MF
  USE MF_Euler_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_Euler_MF
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment_MF
  USE MF_TwoMoment_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_TwoMoment_MF
  USE MF_Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_MF
  USE XCFC_UtilitiesModule, ONLY: &
    nGS, &
    MultiplyWithPsi6_MF
  USE MF_GravitySolutionModule_XCFC, ONLY: &
    ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF, &
    ComputeConformalFactor_MF, &
    ComputePressureTensorTrace_XCFC_TwoMoment_MF, &
    ComputeGeometry_MF
  USE MF_TimersModule, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_UpdateFluid, &
    Timer_AMReX_GravitySolve
  USE MF_TwoMoment_TallyModule, ONLY: &
    IncrementOffGridTally_TwoMoment_MF


  IMPLICIT NONE
  PRIVATE

  INTEGER               :: nStages
  REAL(DP), ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(DP), ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  LOGICAL               :: EvolveEuler
  LOGICAL               :: EvolveTwoMoment
  LOGICAL               :: Verbose

  PUBLIC :: Initialize_IMEX_RK_MF
  PUBLIC :: Finalize_IMEX_RK_MF
  PUBLIC :: Update_IMEX_RK_MF

  LOGICAL, PARAMETER :: DEBUG = .FALSE.

CONTAINS


  SUBROUTINE Update_IMEX_RK_MF

    TYPE(amrex_multifab) :: MF_R    (0:nMaxLevels-1          )
    TYPE(amrex_multifab) :: MF_F    (0:nMaxLevels-1          )
    TYPE(amrex_multifab) :: MF_DR_Im(0:nMaxLevels-1,1:nStages)
    TYPE(amrex_multifab) :: MF_DR_Ex(0:nMaxLevels-1,1:nStages)
    TYPE(amrex_multifab) :: MF_DF_Im(0:nMaxLevels-1,1:nStages)
    TYPE(amrex_multifab) :: MF_DF_Ex(0:nMaxLevels-1,1:nStages)
    TYPE(amrex_multifab) :: MF_uGS  (0:nMaxLevels-1          )

    TYPE(amrex_multifab) :: MF_R0   (0:nMaxLevels-1          )
    TYPE(amrex_multifab) :: MF_F0   (0:nMaxLevels-1          )

    INTEGER :: iS, jS, nCompCF, nCompCR, iLevel

    REAL(DP) :: dM_OffGrid_Euler    (1:nCF  ,0:nMaxLevels-1)
    REAL(DP) :: dM_OffGrid_TwoMoment(1:2*nCR,0:nMaxLevels-1)

    CALL TimersStart_AMReX( Timer_AMReX_UpdateFluid )


    iLevel = 0 ! temporary hack
    dM_OffGrid_Euler     = Zero
    dM_OffGrid_TwoMoment = Zero

    nCompCF = nDOFX * nCF
    nCompCR = nDOFZ * nCR * nE * nSpecies

    IF( DEBUG ) WRITE(*,'(A)') 'Entering Update_IMEX_RK_MF'

    CALL amrex_multifab_build &
           ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
             nDOFX * nGS, 0 )

    IF( DEBUG ) &
      CALL MF_uGS(iLevel) % SetVal( Zero )

    CALL MultiplyWithPsi6_MF( MF_uGF, +1, 1, 1, 1, 1, MF_uCF )
    CALL MultiplyWithPsi6_MF &
           ( MF_uGF, +1, nDOFE, iE_B0, iE_E0, nSpecies, MF_uCR )

    CALL amrex_multifab_build &
           ( MF_F0(iLevel), MF_uCF(iLevel) % BA, &
             MF_uCF(iLevel) % DM, nCompCF, swX )
    CALL MF_F0(iLevel) % COPY( MF_uCF(iLevel), 1, 1, nCompCF, swX )

    CALL amrex_multifab_build &
           ( MF_R0(iLevel), MF_uCR(iLevel) % BA, &
             MF_uCR(iLevel) % DM, nCompCR, swX )
    CALL MF_R0(iLevel) % COPY( MF_uCR(iLevel), 1, 1, nCompCR, swX )

    CALL amrex_multifab_build &
           ( MF_F(iLevel), MF_uCF(iLevel) % BA, &
             MF_uCF(iLevel) % DM, nCompCF, swX )

    CALL amrex_multifab_build &
           ( MF_R(iLevel), MF_uCR(iLevel) % BA, &
             MF_uCR(iLevel) % DM, nCompCR, swX )

    DO iS = 1, nStages

      IF( DEBUG ) WRITE(*,'(2x,A,I2.2)') 'iS: ', iS

      CALL MF_F(iLevel) % COPY( MF_F0(iLevel), 1, 1, nCompCF, swX )
      CALL MF_R(iLevel) % COPY( MF_R0(iLevel), 1, 1, nCompCR, swX )

      CALL amrex_multifab_build &
             ( MF_DF_Ex(iLevel,iS), MF_uCF(iLevel) % BA, &
               MF_uCF(iLevel) % DM, nCompCF, swX )

      CALL amrex_multifab_build &
             ( MF_DF_Im(iLevel,iS), MF_uCF(iLevel) % BA, &
               MF_uCF(iLevel) % DM, nCompCF, swX )

      CALL amrex_multifab_build &
             ( MF_DR_Ex(iLevel,iS), MF_uCR(iLevel) % BA, &
               MF_uCR(iLevel) % DM, nCompCR, swX )

      CALL amrex_multifab_build &
             ( MF_DR_Im(iLevel,iS), MF_uCR(iLevel) % BA, &
               MF_uCR(iLevel) % DM, nCompCR, swX )

      CALL MF_DF_Ex(iLevel,iS) % SetVal( Zero )
      CALL MF_DF_Im(iLevel,iS) % SetVal( Zero )
      CALL MF_DR_Ex(iLevel,iS) % SetVal( Zero )
      CALL MF_DR_Im(iLevel,iS) % SetVal( Zero )

      DO jS = 1, iS - 1

        IF( a_IM(iS,jS) .NE. Zero )THEN

          CALL MF_R(iLevel) &
                 % LinComb &
                     ( One                     , MF_R    (iLevel   ), 1, &
                       dt(iLevel) * a_IM(iS,jS), MF_DR_IM(iLevel,jS), 1, &
                       1, nCompCR, swX )

          CALL MF_F(iLevel) &
                 % LinComb &
                     ( One                     , MF_F    (iLevel   ), 1, &
                       dt(iLevel) * a_IM(iS,jS), MF_DF_IM(iLevel,jS), 1, &
                       1, nCompCF, swX )

        END IF ! a_IM(iS,jS) .NE. Zero

        IF( a_EX(iS,jS) .NE. Zero )THEN

          CALL MF_R(iLevel) &
                 % LinComb &
                     ( One                     , MF_R    (iLevel   ), 1, &
                       dt(iLevel) * a_EX(iS,jS), MF_DR_Ex(iLevel,jS), 1, &
                       1, nCompCR, swX )

          CALL MF_F(iLevel) &
                 % LinComb &
                     ( One                     , MF_F    (iLevel   ), 1, &
                       dt(iLevel) * a_EX(iS,jS), MF_DF_Ex(iLevel,jS), 1, &
                       1, nCompCF, swX )

        END IF ! a_EX(iS,jS) .NE. Zero

        IF( jS .EQ. iS - 1 )THEN

          IF( EvolveEuler )THEN

            CALL ApplySlopeLimiter_Euler_MF &
                   ( MF_uGF, MF_F, MF_uDF )

            CALL ApplyPositivityLimiter_Euler_MF &
                   ( MF_uGF, MF_F, MF_uDF )

          END IF ! EvolveEuler

          IF( EvolveTwoMoment )THEN
            CALL ApplySlopeLimiter_TwoMoment_MF &
                   ( GEOM, MF_uGF, MF_F, MF_R, &
                     Verbose_Option = .FALSE.  )

            CALL ApplyPositivityLimiter_TwoMoment_MF &
                   ( GEOM, MF_uGF, MF_F, MF_R, &
                     Verbose_Option = .FALSE.  )

          END IF ! EvolveTwoMoment

        END IF ! jS .EQ. iS - 1

      END DO ! jS = 1, iS - 1

      IF( DEBUG ) WRITE(*,'(4x,A,I2.2)') 'jS: ', jS

      IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

        IF( DEBUG )THEN

          WRITE(*,'(6x,A)') 'Computing implicit increment'
          WRITE(*,'(6x,A)') 'Adding implicit increment to stage data'

        END IF

        CALL ComputeIncrement_TwoMoment_Implicit_Neutrinos_MF &
             ( GEOM, MF_uGF, MF_F, MF_DF_Im(:,iS), MF_R, MF_DR_Im(:,iS), &
               dt(iLevel) * a_IM(iS,iS), Verbose_Option = .FALSE. )

        CALL MF_R(iLevel) &
               % LinComb &
                   ( One                     , MF_R    (iLevel   ), 1, &
                     dt(iLevel) * a_IM(iS,iS), MF_DR_Im(iLevel,iS), 1, &
                     1, nCompCR, swX )

        CALL MF_F(iLevel) &
               % LinComb &
                   ( One                     , MF_F    (iLevel   ), 1, &
                     dt(iLevel) * a_IM(iS,iS), MF_DF_Im(iLevel,iS), 1, &
                     1, nCompCF, swX )

        IF( EvolveEuler )THEN

          CALL ApplyPositivityLimiter_Euler_MF &
                 ( MF_uGF, MF_F, MF_uDF )

        END IF ! EvolveEuler

        IF( EvolveTwoMoment )THEN

          CALL ApplyPositivityLimiter_TwoMoment_MF &
                 ( GEOM, MF_uGF, MF_F, MF_R, &
                   Verbose_Option = .FALSE.  )

        END IF ! EvolveTwoMoment

      END IF ! ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero )

      IF( ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero ) )THEN

        IF( DEBUG ) &
          WRITE(*,'(6x,A)') 'Computing explicit increment'

        IF( EvolveEuler )THEN

          IF( iS .NE. 1 )THEN

            CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

            CALL ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF &
                   ( MF_uGF, MF_F, MF_R, MF_uGS )

            CALL ComputeConformalFactor_MF( MF_uGS, MF_uGF )

            CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_F )

            CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

            CALL ApplySlopeLimiter_Euler_MF &
                   ( MF_uGF, MF_F, MF_uDF )

            CALL ApplyPositivityLimiter_Euler_MF &
                   ( MF_uGF, MF_F, MF_uDF )

            CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

            CALL MultiplyWithPsi6_MF( MF_uGF, +1, 1, 1, 1, 1, MF_F )

            CALL ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF &
                   ( MF_uGF, MF_F, MF_R, MF_uGS )

            CALL ComputeConformalFactor_MF( MF_uGS, MF_uGF )

            CALL ComputePressureTensorTrace_XCFC_TwoMoment_MF &
                   ( MF_uGF, MF_F, MF_R, MF_uGS )

            CALL ComputeGeometry_MF( MF_uGS, MF_uGF )

            CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

          END IF ! iS .NE. 1

          CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_F )

          CALL ComputeIncrement_Euler_MF &
                 ( MF_uGF, MF_F, MF_uDF, MF_DF_Ex(:,iS) )

        END IF ! EvolveEuler

        IF( EvolveTwoMoment )THEN

          CALL ComputeIncrement_TwoMoment_Explicit_MF &
                 ( t_new, GEOM, MF_uGF, MF_F, MF_R, MF_DR_Ex(:,iS), &
                   Verbose_Option = .FALSE. )

          dM_OffGrid_TwoMoment(:,iLevel) &
            = dM_OffGrid_TwoMoment(:,iLevel) &
            + dt(iLevel) * w_EX(iS) * OffGridFlux_TwoMoment_MF(:,iLevel)



        END IF ! EvolveTwoMoment

      END IF ! ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero )

    END DO ! iS = 1, nStages

    ! --- Assembly Step ---

    IF( ANY( a_IM(nStages,:) .NE. w_IM(:) ) .OR. &
        ANY( a_EX(nStages,:) .NE. w_EX(:) ) )THEN

      IF( DEBUG ) WRITE(*,*) 'Assembly Step'

      CALL MF_F(iLevel) % COPY( MF_F0(iLevel), 1, 1, nCompCF, swX )
      CALL MF_R(iLevel) % COPY( MF_R0(iLevel), 1, 1, nCompCR, swX )

      DO iS = 1, nStages

        IF( w_IM(iS) .NE. Zero )THEN

          IF( DEBUG ) &
            WRITE(*,'(6x,A)') 'Adding implicit increment to original data'

          CALL MF_R(iLevel) &
                 % LinComb &
                     ( One                  , MF_R    (iLevel   ), 1, &
                       dt(iLevel) * w_IM(iS), MF_DR_Im(iLevel,iS), 1, &
                       1, nCompCR, swX )

          CALL MF_F(iLevel) &
                 % LinComb &
                     ( One                  , MF_F    (iLevel   ), 1, &
                       dt(iLevel) * w_IM(iS), MF_DF_Im(iLevel,iS), 1, &
                       1, nCompCF, swX )

        END IF ! w_IM(iS) .NE. Zero

        IF( w_EX(iS) .NE. Zero )THEN

          IF( DEBUG ) &
            WRITE(*,'(6x,A)') 'Adding explicit increment to original data'

          CALL MF_R(iLevel) &
                 % LinComb &
                     ( One                  , MF_R    (iLevel   ), 1, &
                       dt(iLevel) * w_EX(iS), MF_DR_Ex(iLevel,iS), 1, &
                       1, nCompCR, swX )

          CALL MF_F(iLevel) &
                 % LinComb &
                     ( One                  , MF_F    (iLevel   ), 1, &
                       dt(iLevel) * w_EX(iS), MF_DF_Ex(iLevel,iS), 1, &
                       1, nCompCF, swX )

        END IF ! w_EX(iS) .NE. Zero

      END DO ! iS = 1, nStages

      IF( EvolveEuler )THEN

        CALL ApplySlopeLimiter_Euler_MF &
               ( MF_uGF, MF_F, MF_uDF )

        CALL ApplyPositivityLimiter_Euler_MF &
               ( MF_uGF, MF_F, MF_uDF )

     END IF ! EvolveEuler

     IF( EvolveTwoMoment )THEN

        CALL ApplySlopeLimiter_TwoMoment_MF &
               ( GEOM, MF_uGF, MF_F, MF_R, &
                 Verbose_Option = .FALSE. )

        CALL ApplyPositivityLimiter_TwoMoment_MF &
               ( GEOM, MF_uGF, MF_F, MF_R, &
                 Verbose_Option = .FALSE. )

      END IF ! EvolveTwoMoment

    END IF ! ANY( a_IM(nStages,:) .NE. w_IM(:) ) .OR. &
           ! ANY( a_EX(nStages,:) .NE. w_EX(:) )

    CALL MF_uCF(iLevel) % COPY( MF_F(iLevel), 1, 1, nCompCF, swX )
    CALL MF_uCR(iLevel) % COPY( MF_R(iLevel), 1, 1, nCompCR, swX )

    CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

    CALL ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    CALL ComputeConformalFactor_MF( MF_uGS, MF_uGF )

    CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_uCF )

    CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

    CALL ApplySlopeLimiter_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uDF )

    CALL ApplyPositivityLimiter_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uDF )

    CALL TimersStart_AMReX( Timer_AMReX_GravitySolve )

    CALL MultiplyWithPsi6_MF( MF_uGF, +1, 1, 1, 1, 1, MF_uCF )

    CALL ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    CALL ComputeConformalFactor_MF( MF_uGS, MF_uGF )

    CALL ComputePressureTensorTrace_XCFC_TwoMoment_MF( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    CALL ComputeGeometry_MF( MF_uGS, MF_uGF )

    CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_uCF )
    CALL MultiplyWithPsi6_MF &
           ( MF_uGF, -1, nDOFE, iE_B0, iE_E0, nSpecies, MF_uCR )

    CALL TimersStop_AMReX( Timer_AMReX_GravitySolve )

    DO iS = 1, nStages

      CALL amrex_multifab_destroy( MF_DR_Im(iLevel,iS) )
      CALL amrex_multifab_destroy( MF_DR_Ex(iLevel,iS) )
      CALL amrex_multifab_destroy( MF_DF_Im(iLevel,iS) )
      CALL amrex_multifab_destroy( MF_DF_Ex(iLevel,iS) )

    END DO ! iS = 0, nStages


    CALL IncrementOffGridTally_TwoMoment_MF( dM_OffGrid_TwoMoment )

    CALL amrex_multifab_destroy( MF_R  (iLevel) )
    CALL amrex_multifab_destroy( MF_F  (iLevel) )
    CALL amrex_multifab_destroy( MF_R0 (iLevel) )
    CALL amrex_multifab_destroy( MF_F0 (iLevel) )
    CALL amrex_multifab_destroy( MF_uGS(iLevel) )

    IF( DEBUG ) WRITE(*,'(A)') 'Leaving Update_IMEX_RK_MF'

    CALL TimersStop_AMReX( Timer_AMReX_UpdateFluid )

  END SUBROUTINE Update_IMEX_RK_MF


  SUBROUTINE Initialize_IMEX_RK_MF &
    ( Scheme, EvolveEuler_Option, EvolveTwoMoment_Option, Verbose_Option )

    CHARACTER(LEN=*), INTENT(in) :: Scheme
    LOGICAL         , INTENT(in), OPTIONAL :: EvolveEuler_Option
    LOGICAL         , INTENT(in), OPTIONAL :: EvolveTwoMoment_Option
    LOGICAL         , INTENT(in), OPTIONAL :: Verbose_Option

    EvolveEuler = .FALSE.
    IF( PRESENT( EvolveEuler_Option ) ) &
      EvolveEuler = EvolveEuler_Option

    EvolveTwoMoment = .FALSE.
    IF( PRESENT( EvolveTwoMoment_Option ) ) &
      EvolveTwoMoment = EvolveTwoMoment_Option

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL Initialize_IMEX_RK( Scheme )

  END SUBROUTINE Initialize_IMEX_RK_MF


  SUBROUTINE Finalize_IMEX_RK_MF

    CALL Finalize_IMEX_RK

  END SUBROUTINE Finalize_IMEX_RK_MF


  SUBROUTINE Initialize_IMEX_RK( Scheme )

    CHARACTER(LEN=*), INTENT(in) :: Scheme

    INTEGER :: i

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A)')   'INFO: Time-stepper'
      WRITE(*,'(4x,A)')   '------------------'
      WRITE(*,*)
      WRITE(*,'(6x,A,A)') 'IMEX-RK Scheme: ', TRIM( Scheme )

    END IF

    SELECT CASE ( TRIM( Scheme ) )

      CASE ( 'BackwardEuler' )

        nStages = 1

        CALL AllocateButcherTables

        a_IM(1,1) = 1.0_DP
        w_IM(1)   = 1.0_DP

      CASE ( 'SSPRK1' )

        nStages = 1

        CALL AllocateButcherTables

        a_EX(1,1) = Zero
        w_EX(1)   = 1.0_DP

      CASE ( 'SSPRK2' )

        nStages = 2

        CALL AllocateButcherTables

        a_EX(1,1:2) = [ Zero, Zero ]
        a_EX(2,1:2) = [ 1.0_DP, Zero ]
        w_EX(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 'SSPRK3' )

        nStages = 3

        CALL AllocateButcherTables

        a_EX(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_EX(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_EX(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_EX(1:3)   = [ 1.00_DP, 1.00_DP, 4.00_DP ] / 6.0_DP

      CASE ( 'IMEX_ARS_111' )

        nStages = 2

        CALL AllocateButcherTables

        ! --- Coefficients from Ascher et al. (1997) ---

        a_EX(2,1) = 1.0_DP
        w_EX(1)   = a_EX(2,1)

        a_IM(2,2) = 1.0_DP
        w_IM(2)   = a_IM(2,2)

      CASE ( 'IMEX_PDARS' )

        nStages = 3

        CALL AllocateButcherTables

        a_EX(2,1) = 1.0_DP
        a_EX(3,1) = 0.5_DP
        a_EX(3,2) = 0.5_DP

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        a_IM(2,2) = 1.0_DP
        a_IM(3,2) = 0.5_DP
        a_IM(3,3) = 0.5_DP

        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

      CASE DEFAULT

        CALL DescribeError_MF( 105, Char_Option = [ TRIM( Scheme ) ] )

    END SELECT

    DO i = 1, nStages

      c_IM(i) = SUM( a_IM(i,1:i) )
      c_EX(i) = SUM( a_EX(i,1:i-1) )

    END DO

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'Implicit Butcher Table:'
      WRITE(*,'(A6,A)') '', '-----------------------'
      DO i = 1, nStages
        WRITE(*,'(A6,5ES14.4E3)') '', c_IM(i), a_IM(i,1:nStages)
      END DO
      WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_IM(1:nStages)

      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'Explicit Butcher Table:'
      WRITE(*,'(A6,A)') '', '-----------------------'
      DO i = 1, nStages
        WRITE(*,'(A6,5ES14.4E3)') '', c_EX(i), a_EX(i,1:nStages)
      END DO
      WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_EX(1:nStages)

    END IF

  END SUBROUTINE Initialize_IMEX_RK


  SUBROUTINE Finalize_IMEX_RK

    DEALLOCATE( c_IM, w_IM, a_IM )
    DEALLOCATE( c_EX, w_EX, a_EX )

  END SUBROUTINE Finalize_IMEX_RK


  SUBROUTINE AllocateButcherTables

    ! --- Implicit Coefficients ---

    ALLOCATE( c_IM(nStages) )
    ALLOCATE( w_IM(nStages) )
    ALLOCATE( a_IM(nStages,nStages) )

    c_IM = Zero
    w_IM = Zero
    a_IM = Zero

    ! --- Explicit Coefficients ---

    ALLOCATE( c_EX(nStages) )
    ALLOCATE( w_EX(nStages) )
    ALLOCATE( a_EX(nStages,nStages) )

    c_EX = Zero
    w_EX = Zero
    a_EX = Zero

  END SUBROUTINE AllocateButcherTables


END MODULE MF_TimeSteppingModule_IMEX
