MODULE MF_TimeSteppingModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray, &
    amrex_boxarray_build, &
    amrex_boxarray_destroy
  USE amrex_distromap_module, ONLY: &
    amrex_distromap, &
    amrex_distromap_build, &
    amrex_distromap_destroy
  USE amrex_amrcore_module, ONLY: &
    GEOM => amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFZ, &
    nDOFX, &
    nDOFE, &
    iZ_B0, &
    iZ_B1, &
    iZ_E0, &
    iZ_E1, &
    swX, &
    nX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nCR
  USE MF_TwoMoment_DiscretizationModule_Streaming_Relativistic, ONLY: &
    ComputeIncrement_TwoMoment_Explicit_MF
  USE MF_TwoMoment_DiscretizationModule_Collisions_Neutrinos_GR, ONLY: &
    ComputeIncrement_TwoMoment_Implicit_Neutrinos_MF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uDF
  USE MF_FieldsModule_TwoMoment, ONLY: &
    MF_uCR, &
    MF_uPR
  USE InputParsingModule, ONLY: &
    t_new, &
    dt, &
    nLevels, &
    nMaxLevels, &
!    DEBUG, &
    UseTiling, &
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
  USE MF_GravitySolutionModule_XCFC_Poseidon, ONLY: &
    nGS, &
    MultiplyWithPsi6_MF, &
    ComputeConformalFactorSourcesAndMg_XCFC_MF, &
    ComputeConformalFactor_Poseidon_MF, &
    ComputePressureTensorTrace_XCFC_MF, &
    ComputeGeometry_Poseidon_MF

  IMPLICIT NONE
  PRIVATE

  TYPE :: StageDataType
    REAL(DP), ALLOCATABLE :: dR_IM(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dR_EX(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dF_EX(:,:,:,:,:)
  END TYPE StageDataType

  INTEGER                           :: nStages
  REAL(DP)            , ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(DP)            , ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  REAL(DP)            , ALLOCATABLE :: U0(:,:,:,:,:,:,:)
  REAL(DP)            , ALLOCATABLE :: Ui(:,:,:,:,:,:,:)
  TYPE(StageDataType) , ALLOCATABLE :: StageData(:)

  PUBLIC :: Initialize_IMEX_RK_MF
  PUBLIC :: Finalize_IMEX_RK_MF
  PUBLIC :: Update_IMEX_RK_MF

  LOGICAL, PARAMETER :: DEBUG = .TRUE.

CONTAINS


  SUBROUTINE Update_IMEX_RK_MF

    TYPE(amrex_multifab) :: MF_R    (0:nMaxLevels-1,1:nStages)
    TYPE(amrex_multifab) :: MF_F    (0:nMaxLevels-1,1:nStages)
    TYPE(amrex_multifab) :: MF_DR_Im(0:nMaxLevels-1,1:nStages)
    TYPE(amrex_multifab) :: MF_DR_Ex(0:nMaxLevels-1,1:nStages)
    TYPE(amrex_multifab) :: MF_DF_Im(0:nMaxLevels-1,1:nStages)
    TYPE(amrex_multifab) :: MF_DF_Ex(0:nMaxLevels-1,1:nStages)
    TYPE(amrex_multifab) :: MF_uGS  (0:nMaxLevels-1          )

    INTEGER :: iS, jS, nCompCF, nCompCR, iLevel

    REAL(DP) :: dM_OffGrid_Euler    (1:nCF,0:nMaxLevels-1)
    REAL(DP) :: dM_OffGrid_TwoMoment(1:nCR,0:nMaxLevels-1)

    LOGICAL, PARAMETER :: Verbose = .FALSE.

    dM_OffGrid_Euler     = Zero
    dM_OffGrid_TwoMoment = Zero

    nCompCF = nDOFX * nCF
    nCompCR = nDOFZ * nCR * nE * nSpecies

    IF( DEBUG ) WRITE(*,'(A)') 'Entering Update_IMEX_RK_MF'

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nGS, 0 )

      IF( DEBUG ) &
        CALL MF_uGS(iLevel) % SetVal( Zero )

    END DO

    CALL MultiplyWithPsi6_MF( MF_uGF, +1, MF_uCF )
    !!$CALL MultiplyWithPsi6_MF( MF_uGF, +1, MF_uCR ) (needs mods for radiation)

    DO iS = 1, nStages

      IF( DEBUG ) WRITE(*,'(2x,A,I2.2)') 'iS: ', iS

      DO iLevel = 0, nLevels-1

        CALL amrex_multifab_build &
               ( MF_F(iLevel,iS), MF_uCF(iLevel) % BA, &
                 MF_uCF(iLevel) % DM, nCompCF, swX )

        CALL amrex_multifab_build &
               ( MF_DF_Ex(iLevel,iS), MF_uCF(iLevel) % BA, &
                 MF_uCF(iLevel) % DM, nCompCF, swX )

        CALL amrex_multifab_build &
               ( MF_DF_Im(iLevel,iS), MF_uCF(iLevel) % BA, &
                 MF_uCF(iLevel) % DM, nCompCF, swX )

        CALL amrex_multifab_build &
               ( MF_R(iLevel,iS), MF_uCR(iLevel) % BA, &
                 MF_uCR(iLevel) % DM, nCompCR, swX )

        CALL amrex_multifab_build &
               ( MF_DR_Ex(iLevel,iS), MF_uCR(iLevel) % BA, &
                 MF_uCR(iLevel) % DM, nCompCR, swX )

        CALL amrex_multifab_build &
               ( MF_DR_Im(iLevel,iS), MF_uCR(iLevel) % BA, &
                 MF_uCR(iLevel) % DM, nCompCR, swX )

        IF( DEBUG )THEN

          CALL MF_F    (iLevel,iS) % SetVal( Zero )
          CALL MF_DF_Ex(iLevel,iS) % SetVal( Zero )
          CALL MF_DF_Im(iLevel,iS) % SetVal( Zero )
          CALL MF_R    (iLevel,iS) % SetVal( Zero )
          CALL MF_DR_Ex(iLevel,iS) % SetVal( Zero )
          CALL MF_DR_Im(iLevel,iS) % SetVal( Zero )

        END IF ! DEBUG

        CALL MF_F(iLevel,iS) % COPY( MF_uCF(iLevel), 1, 1, nCompCF, swX )
        CALL MF_R(iLevel,iS) % COPY( MF_uCR(iLevel), 1, 1, nCompCR, swX )

      END DO ! iLevel = 0, nLevels-1

      DO jS = 1, iS-1

        IF( DEBUG ) WRITE(*,'(4x,A,I2.2)') 'jS: ', jS

        DO iLevel = 0, nLevels-1

          IF( a_EX(iS,jS) .NE. Zero )THEN

            IF( DEBUG ) &
              WRITE(*,'(6x,A)') 'Adding explicit increment to stage data'

            CALL MF_R(iLevel,iS) &
                   % LinComb &
                       ( One                     , MF_R    (iLevel,iS), 1, &
                         dt(iLevel) * a_EX(iS,jS), MF_DR_Ex(iLevel,jS), 1, &
                         1, nCompCR, swX )

            CALL MF_F(iLevel,iS) &
                   % LinComb &
                       ( One                     , MF_F    (iLevel,iS), 1, &
                         dt(iLevel) * a_EX(iS,jS), MF_DF_Ex(iLevel,jS), 1, &
                         1, nCompCF, swX )

          END IF

          IF( a_IM(iS,jS) .NE. Zero )THEN

            IF( DEBUG ) &
              WRITE(*,'(6x,A)') 'Adding implicit increment to stage data'

            CALL MF_R(iLevel,iS) &
                 % LinComb &
                     ( One                     , MF_R    (iLevel,iS), 1, &
                       dt(iLevel) * a_IM(iS,iS), MF_DR_Im(iLevel,jS), 1, &
                       1, nCompCR, swX )

            CALL MF_F(iLevel,iS) &
                 % LinComb &
                     ( One                     , MF_F    (iLevel,iS), 1, &
                       dt(iLevel) * a_IM(iS,jS), MF_DF_Im(iLevel,jS), 1, &
                       1, nCompCF, swX )

          END IF

        END DO ! iLevel = 0, nLevels-1

!!$        IF( a_SSPRK(iS,jS) .NE. Zero ) &
!!$          CALL AverageDown( MF_uGF, MF_F(iS,:) )

        IF( jS .EQ. iS - 1 )THEN

          CALL ApplySlopeLimiter_TwoMoment_MF &
                 ( GEOM, MF_uGF, MF_F(:,jS), MF_R(:,jS), &
                   Verbose_Option = Verbose  )

          CALL ApplyPositivityLimiter_TwoMoment_MF &
                 ( GEOM, MF_uGF, MF_F(:,jS), MF_R(:,jS), &
                   Verbose_Option = Verbose  )

          CALL ApplySlopeLimiter_Euler_MF &
                 ( t_new, MF_uGF, MF_F(:,jS), MF_uDF )

          CALL ApplyPositivityLimiter_Euler_MF &
                 ( MF_uGF, MF_F(:,jS), MF_uDF )

        END IF

      END DO ! jS = 1, iS-1

      IF( DEBUG ) WRITE(*,'(4x,A,I2.2)') 'jS: ', jS

      IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

        IF( DEBUG )THEN

          WRITE(*,'(6x,A)') 'Computing implicit increment'
          WRITE(*,'(6x,A)') 'Adding implicit increment to stage data'

        END IF

        DO iLevel = 0, nLevels-1

          CALL ComputeIncrement_TwoMoment_Implicit_Neutrinos_MF &
               ( GEOM, MF_uGF, MF_F(:,iS), MF_DF_Im(:,iS), &
                 MF_R(:,iS), MF_DR_Im(:,iS), &
                 dt(iLevel) * a_IM(iS,iS), Verbose_Option = Verbose )

          CALL MF_R(iLevel,iS) &
                 % LinComb &
                     ( One                     , MF_R    (iLevel,iS), 1, &
                       dt(iLevel) * a_IM(iS,iS), MF_DR_Im(iLevel,iS), 1, &
                       1, nCompCR, swX )

          CALL MF_F(iLevel,iS) &
                 % LinComb &
                     ( One                     , MF_F    (iLevel,iS), 1, &
                       dt(iLevel) * a_IM(iS,iS), MF_DF_Im(iLevel,iS), 1, &
                       1, nCompCF, swX )

        END DO ! iLevel = 0, nLevels-1

      END IF ! ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero )

      IF( ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero ) )THEN

        IF( DEBUG ) &
          WRITE(*,'(6x,A)') 'Computing explicit increment'

        CALL ComputeIncrement_TwoMoment_Explicit_MF &
               ( GEOM, MF_uGF, MF_F(:,iS), MF_R(:,iS), MF_DR_Ex(:,iS), &
                 Verbose_Option = Verbose )

        IF( iS .NE. 1 )THEN

          CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
                 ( MF_uGF, MF_F(:,iS), MF_uGS )

          CALL ComputeConformalFactor_Poseidon_MF( MF_uGS, MF_uGF )

          CALL MultiplyWithPsi6_MF( MF_uGF, -1, MF_F(:,iS) )

          CALL ApplySlopeLimiter_Euler_MF &
                 ( t_new, MF_uGF, MF_F(:,iS), MF_uDF )

          CALL ApplyPositivityLimiter_Euler_MF &
                 ( MF_uGF, MF_F(:,iS), MF_uDF )

          CALL MultiplyWithPsi6_MF( MF_uGF, +1, MF_F(:,iS) )

          CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
                 ( MF_uGF, MF_F(:,iS), MF_uGS )

          CALL ComputeConformalFactor_Poseidon_MF( MF_uGS, MF_uGF )

          CALL ComputePressureTensorTrace_XCFC_MF( MF_uGF, MF_F(:,iS), MF_uGS )

          CALL ComputeGeometry_Poseidon_MF( MF_uGS, MF_uGF )

        END IF ! iS .NE. 1

        CALL MultiplyWithPsi6_MF( MF_uGF, -1, MF_F(:,iS) )

        CALL ComputeIncrement_Euler_MF &
               ( t_new, MF_uGF, MF_F(:,iS), MF_uDF, MF_DF_Ex(:,iS) )

      END IF ! ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero )

    END DO ! iS = 1, nStages

    ! --- Assembly Step ---

    IF( DEBUG ) WRITE(*,*) 'Assembly Step'

    DO iLevel = 0, nLevels-1

      IF( ANY( a_IM(nStages,:) .NE. w_IM(:) ) .OR. &
          ANY( a_EX(nStages,:) .NE. w_EX(:) ) )THEN

!!$        U = uCR
!!$        F = uCF

        DO iS = 1, nStages

          IF( w_IM(iS) .NE. Zero )THEN

            IF( DEBUG ) &
              WRITE(*,'(6x,A)') 'Adding implicit increment to original data'

            CALL MF_uCR(iLevel) &
                   % LinComb &
                       ( One                  , MF_uCR  (iLevel   ), 1, &
                         dt(iLevel) * w_IM(iS), MF_DR_Im(iLevel,iS), 1, &
                         1, nCompCR, swX )

            CALL MF_uCF(iLevel) &
                   % LinComb &
                       ( One                  , MF_uCF  (iLevel   ), 1, &
                         dt(iLevel) * w_IM(iS), MF_DF_Im(iLevel,iS), 1, &
                         1, nCompCF, swX )

          END IF

          IF( w_EX(iS) .NE. Zero )THEN

            IF( DEBUG ) &
              WRITE(*,'(6x,A)') 'Adding explicit increment to original data'

            CALL MF_uCR(iLevel) &
                   % LinComb &
                       ( One                  , MF_uCR  (iLevel   ), 1, &
                         dt(iLevel) * w_EX(iS), MF_DR_Ex(iLevel,iS), 1, &
                         1, nCompCR, swX )

            CALL MF_uCF(iLevel) &
                   % LinComb &
                       ( One                  , MF_uCF  (iLevel   ), 1, &
                         dt(iLevel) * w_EX(iS), MF_DF_Ex(iLevel,iS), 1, &
                         1, nCompCF, swX )

          END IF

        END DO ! iS = 1, nStages

      END IF ! ANY( a_IM(nStages,:) .NE. w_IM(:) ) .OR. &
             ! ANY( a_EX(nStages,:) .NE. w_EX(:) )

    END DO ! iLevel = 0, nLevels-1

!!$    CALL ApplySlopeLimiter_TwoMoment_MF &
!!$           ( GEOM, MF_uGF, MF_F, MF_U, Verbose_Option = Verbose )
!!$
!!$    CALL ApplyPositivityLimiter_TwoMoment_MF &
!!$           ( GEOM, MF_uGF, MF_F, MF_U, Verbose_Option = Verbose )

!!$    uCR = U
!!$    uCF = F

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

    IF( DEBUG ) WRITE(*,'(A)') 'Leaving Update_IMEX_RK_MF'

  END SUBROUTINE Update_IMEX_RK_MF


  SUBROUTINE Initialize_IMEX_RK_MF( Scheme, Verbose_Option )

    CHARACTER(LEN=*), INTENT(in) :: Scheme
    LOGICAL         , INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL Initialize_IMEX_RK( Scheme, Verbose )

  END SUBROUTINE Initialize_IMEX_RK_MF


  SUBROUTINE Finalize_IMEX_RK_MF

    CALL Finalize_IMEX_RK

  END SUBROUTINE Finalize_IMEX_RK_MF


  SUBROUTINE Initialize_IMEX_RK &
    ( Scheme, Verbose, EvolveEuler_Option, EvolveTwoMoment_Option )

    CHARACTER(LEN=*), INTENT(in) :: Scheme
    LOGICAL         , INTENT(in) :: Verbose
    LOGICAL         , INTENT(in), OPTIONAL :: EvolveEuler_Option
    LOGICAL         , INTENT(in), OPTIONAL :: EvolveTwoMoment_Option

    INTEGER :: i
    LOGICAL :: EvolveEuler
    LOGICAL :: EvolveTwoMoment

    EvolveEuler = .FALSE.
    IF( PRESENT( EvolveEuler_Option ) ) &
      EvolveEuler = EvolveEuler_Option

    EvolveTwoMoment = .FALSE.
    IF( PRESENT( EvolveTwoMoment_Option ) ) &
      EvolveTwoMoment = EvolveTwoMoment_Option

    IF (Verbose) THEN
      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', 'IMEX-RK Scheme: ', TRIM( Scheme )
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

        WRITE(*,*)
        WRITE(*,'(A6,A,A)') &
          '', 'Unknown Time Stepping Scheme: ', TRIM( Scheme )
        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'Available Options:'
        WRITE(*,*)
        WRITE(*,'(A6,A)') '', 'BackwardEuler'
        WRITE(*,'(A6,A)') '', 'SSPRK1'
        WRITE(*,'(A6,A)') '', 'SSPRK2'
        WRITE(*,'(A6,A)') '', 'SSPRK3'
        WRITE(*,'(A6,A)') '', 'IMEX_ARS_111'
        WRITE(*,'(A6,A)') '', 'IMEX_PDARS'
        WRITE(*,*)
        STOP

    END SELECT

    DO i = 1, nStages
      c_IM(i) = SUM( a_IM(i,1:i) )
      c_EX(i) = SUM( a_EX(i,1:i-1) )
    END DO

    IF (Verbose) THEN
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

    CALL AllocateArray7D( U0 )
    CALL AllocateArray7D( Ui )

    ALLOCATE( StageData(nStages) )

    DO i = 1, nStages

      CALL AllocateArray7D( StageData(i) % dR_IM )
      CALL AllocateArray7D( StageData(i) % dR_EX )
      CALL AllocateArray5D( StageData(i) % dF_EX )

    END DO

  END SUBROUTINE Initialize_IMEX_RK


  SUBROUTINE Finalize_IMEX_RK

    INTEGER :: i

    DEALLOCATE( c_IM, w_IM, a_IM )
    DEALLOCATE( c_EX, w_EX, a_EX )

    CALL DeallocateArray7D( U0 )
    CALL DeallocateArray7D( Ui )

    DO i = 1, nStages

      CALL DeallocateArray7D( StageData(i) % dR_IM )
      CALL DeallocateArray7D( StageData(i) % dR_EX )
      CALL DeallocateArray5D( StageData(i) % dF_EX )

    END DO

    DEALLOCATE( StageData )

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


  SUBROUTINE AllocateArray7D( Array7D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    ALLOCATE &
      ( Array7D(nDOFZ, &
                iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
                nCR,nSpecies) )

    Array7D = Zero

  END SUBROUTINE AllocateArray7D


  SUBROUTINE AllocateArray5D( Array5D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array5D(:,:,:,:,:)

    ALLOCATE &
      ( Array5D(nDOFZ, &
                iZ_B1(2):iZ_E1(2), &
                iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
                nCF) )

    Array5D = Zero

  END SUBROUTINE AllocateArray5D


  SUBROUTINE DeallocateArray7D( Array7D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    DEALLOCATE( Array7D )

  END SUBROUTINE DeallocateArray7D


  SUBROUTINE DeallocateArray5D( Array5D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array5D(:,:,:,:,:)

    DEALLOCATE( Array5D )

  END SUBROUTINE DeallocateArray5D

END MODULE MF_TimeSteppingModule
