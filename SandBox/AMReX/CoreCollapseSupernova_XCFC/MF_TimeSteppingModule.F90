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
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    nSpecies
  USE MF_TwoMoment_DiscretizationModule_Streaming_Relativistic, ONLY: &
    MF_TwoMoment_ComputeIncrement_Explicit
!!$  USE MF_TwoMoment_DiscretizationModule_Collisions_Neutrinos_GR, ONLY: &
!!$    MF_TwoMoment_ComputeIncrement_Implicit_Neutrinos
  USE MF_TwoMoment_PositivityLimiter, ONLY: &
    MF_TwoMoment_ApplyPositivityLimiter
  USE MF_TwoMoment_SlopeLimiter, ONLY: &
    MF_TwoMoment_ApplySlopeLimiter

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE InputParsingModule, ONLY: &
    nLevels, &
    DEBUG, &
    UseTiling
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
    REAL(DP), ALLOCATABLE :: dU_IM(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_EX(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dF_EX(:,:,:,:,:)
  END TYPE StageDataType

  INTEGER                           :: nStages
  REAL(DP)            , ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(DP)            , ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_U(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_F(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_DU_Im(:,:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_DU_Ex(:,:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_DF_Im(:,:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_DF_Ex(:,:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uGS  (:,:)
  REAL(DP)            , ALLOCATABLE :: U0(:,:,:,:,:,:,:)
  REAL(DP)            , ALLOCATABLE :: Ui(:,:,:,:,:,:,:)
  TYPE(StageDataType) , ALLOCATABLE :: StageData(:)

  PUBLIC :: Initialize_IMEX_RK
  PUBLIC :: Finalize_IMEX_RK
  PUBLIC :: MF_InitializeField_IMEX_RK
  PUBLIC :: MF_FinalizeField_IMEX_RK
  PUBLIC :: MF_Update_IMEX_RK

CONTAINS


  SUBROUTINE MF_Update_IMEX_RK &
    ( t, dt, uGE, MF_uGF, MF_uCF, MF_uCR, GEOM, Verbose_Option )

    REAL(DP)            , INTENT(in)    :: t(0:nLevels-1), dt(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)
    REAL(DP)            , INTENT(in)    :: uGE(:,:,:)
    LOGICAL             , INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iLevel, iS, jS

    ! --- For physical boundary conditions ---
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:), U(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:), F(:,:,:,:)

    LOGICAL :: Verbose
    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    ! --- Set temporary MultiFabs U and dU to zero --
    DO iLevel = 0, nLevels-1

      CALL MF_U(iLevel) % setval( Zero )
      CALL MF_F(iLevel) % setval( Zero )

      DO iS = 1, nStages

        CALL MF_DU_Ex(iLevel,iS) % setval( Zero )
        CALL MF_DU_Im(iLevel,iS) % setval( Zero )
        CALL MF_DF_Ex(iLevel,iS) % setval( Zero )
        CALL MF_DF_Im(iLevel,iS) % setval( Zero )

      END DO

    END DO

    CALL MultiplyWithPsi6_MF( MF_uGF, +1, MF_uCF )
    CALL MultiplyWithPsi6_MF( MF_uGF, +1, MF_uCR )

    DO iS = 1, nStages

      ! --- Copy data from input MultiFab to temporary MultiFab ---

      DO iLevel = 0, nLevels-1

        CALL MF_U(iLevel) &
               % COPY( MF_uCR(iLevel), 1, 1, &
                       MF_uCR(iLevel) % nComp(), swX )

        CALL MF_F(iLevel) &
               % COPY( MF_uCF(iLevel), 1, 1, &
                       MF_uCF(iLevel) % nComp(), swX )

        ! --- Apply boundary conditions to interior domains ---

        CALL MF_U(iLevel) % Fill_Boundary( GEOM(iLevel) )
        CALL MF_F(iLevel) % Fill_Boundary( GEOM(iLevel) )

        ! --- Copy ghost data from physical boundaries ---

        CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

        DO WHILE( MFI % next() )

          uCR => MF_uCR(iLevel) % DataPtr( MFI )
          U   => MF_U  (iLevel) % DataPtr( MFI )
          U   =  uCR

          uCF => MF_uCF(iLevel) % DataPtr( MFI )
          F   => MF_F  (iLevel) % DataPtr( MFI )
          F   =  uCF

        END DO

        CALL amrex_mfiter_destroy( MFI )

      END DO

      DO iLevel = 0, nLevels-1

        DO jS = 1, iS - 1

          IF( a_EX(iS,jS) .NE. Zero )THEN

            CALL MF_U(iLevel) &
                   % LinComb( 1.0_DP,              MF_U(iLevel),    1, &
                              dt(iLevel) * a_EX(iS,jS), MF_DU_Ex(iLevel,jS), 1, &
                              1, MF_U(iLevel) % nComp(), swX )

            CALL MF_F(iLevel) &
                   % LinComb( 1.0_DP,              MF_F(iLevel),    1, &
                              dt(iLevel) * a_EX(iS,jS), MF_DF_Ex(iLevel,jS), 1, &
                              1, MF_F(iLevel) % nComp(), swX )

          END IF

          IF( a_IM(iS,jS) .NE. Zero )THEN

            CALL MF_U(iLevel) &
                 % LinComb( 1.0_DP,              MF_U(iLevel),    1, &
                            dt(iLevel) * a_IM(iS,iS), MF_DU_Im(iLevel,iS), 1, &
                            1, MF_U(iLevel) % nComp(), swX )

            CALL MF_F(iLevel) &
                 % LinComb( 1.0_DP,              MF_F(iLevel),    1, &
                            dt(iLevel) * a_IM(iS,jS), MF_DF_Im(iLevel,jS), 1, &
                            1, MF_F(iLevel) % nComp(), swX )

          END IF

          IF( jS == iS - 1 )THEN

            ! --- Apply Limiters ---

            CALL MF_TwoMoment_ApplySlopeLimiter &
                   ( GEOM, MF_uGF, MF_F, MF_U, Verbose_Option = Verbose  )

            CALL MF_TwoMoment_ApplyPositivityLimiter &
                   ( GEOM, MF_uGF, MF_F, MF_U, Verbose_Option = Verbose  )

          END IF

        END DO ! jS = 1, iS - 1

        IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

          IF (Verbose) THEN
            PRINT*, "    IMPLICIT: ", iS
          END IF

          CALL MF_TwoMoment_ComputeIncrement_Implicit_Neutrinos &
               ( GEOM, MF_uGF, MF_uCF, MF_DF_Im(:,iS), MF_U, MF_DU_Im(:,iS), &
                 dt(iLevel) * a_IM(iS,iS), Verbose_Option = Verbose )

          CALL MF_U(iLevel) &
                 % LinComb( 1.0_DP,              MF_U(iLevel),    1, &
                            dt(iLevel) * a_IM(iS,iS), MF_DU_Im(iLevel,iS), 1, &
                            1, MF_U(iLevel) % nComp(), swX )

          CALL MF_F(iLevel) &
                 % LinComb( 1.0_DP,              MF_F(iLevel),    1, &
                            dt(iLevel) * a_IM(iS,iS), MF_DF_Im(iLevel,iS), 1, &
                            1, MF_F(iLevel) % nComp(), swX )

        END IF

        IF( ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero ) )THEN

          ! --- Explicit Solve ---
          IF (Verbose) THEN
            PRINT*, "    EXPLICIT: ", iS
          END IF

          CALL MF_TwoMoment_ComputeIncrement_Explicit &
                 ( GEOM, MF_uGF, MF_F, MF_U, MF_DU_Ex(:,iS), &
                   Verbose_Option = Verbose )

          IF( iS .NE. 1 )THEN

            CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
                   ( MF_uGF, MF_F(iS,:), MF_uGS )

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

          END IF

          CALL MultiplyWithPsi6_MF( MF_uGF, -1, MF_F(:,iS) )

          CALL ComputeIncrement_Euler_MF &
                 ( t_new, MF_uGF, MF_F, MF_uDF, MF_DF_Ex(:,iS) )

        END IF

      END DO ! iLevel = 0, nLevels-1

    END DO ! iS = 1, nStages

    ! --- Assembly Step ---

    DO iLevel = 0, nLevels-1

      IF( ANY( a_IM(nStages,:) .NE. w_IM(:) ) .OR. &
          ANY( a_EX(nStages,:) .NE. w_EX(:) ) )THEN

        U = uCR
        F = uCF

        IF (Verbose) THEN
          PRINT*, "    ASSEMBLY:"
        END IF

        DO iS = 1, nStages

          IF( w_IM(iS) .NE. Zero )THEN

            CALL MF_U(iLevel) &
                   % LinComb( 1.0_DP,              MF_U(iLevel),    1, &
                              dt(iLevel) * w_IM(iS), MF_DU_Im(iLevel,iS), 1, &
                              1, MF_U(iLevel) % nComp(), swX )

            CALL MF_F(iLevel) &
                   % LinComb( 1.0_DP,              MF_F(iLevel),    1, &
                              dt(iLevel) * w_IM(iS), MF_DF_Im(iLevel,iS), 1, &
                              1, MF_F(iLevel) % nComp(), swX )

          END IF

          IF( w_EX(iS) .NE. Zero )THEN

            CALL MF_U(iLevel) &
                   % LinComb( 1.0_DP,              MF_U(iLevel),    1, &
                              dt(iLevel) * w_EX(iS), MF_DU_Ex(iLevel,iS), 1, &
                              1, MF_U(iLevel) % nComp(), swX )

            CALL MF_F(iLevel) &
                   % LinComb( 1.0_DP,              MF_F(iLevel),    1, &
                              dt(iLevel) * w_EX(iS), MF_DF_Ex(iLevel,iS), 1, &
                              1, MF_F(iLevel) % nComp(), swX )

          END IF

        END DO

        CALL MF_TwoMoment_ApplySlopeLimiter &
               ( GEOM, MF_uGF, MF_F, MF_U, Verbose_Option = Verbose )

        CALL MF_TwoMoment_ApplyPositivityLimiter &
               ( GEOM, MF_uGF, MF_F, MF_U, Verbose_Option = Verbose )

      END IF

    END DO

    uCR = U
    uCF = F

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

  END SUBROUTINE MF_Update_IMEX_RK


  SUBROUTINE MF_InitializeField_IMEX_RK &
    ( Scheme, BA, DM, Verbose_Option )

    CHARACTER(LEN=*)     , INTENT(in) :: Scheme
    TYPE(amrex_boxarray) , INTENT(in) :: BA(0:nLevels-1)
    TYPE(amrex_distromap), INTENT(in) :: DM(0:nLevels-1)
    LOGICAL              , INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER         :: iS, iLevel
    TYPE(amrex_box) :: BX

    CALL Initialize_IMEX_RK( Scheme, Verbose_Option )

    ALLOCATE( MF_U  (0:nLevels-1) )
    ALLOCATE( MF_F  (0:nLevels-1) )
    ALLOCATE( MF_uGS(0:nLevels-1) )
    ALLOCATE( MF_DU_Ex(0:nLevels-1,1:nStages) )
    ALLOCATE( MF_DU_Im(0:nLevels-1,1:nStages) )
    ALLOCATE( MF_DF_Ex(0:nLevels-1,1:nStages) )
    ALLOCATE( MF_DF_Im(0:nLevels-1,1:nStages) )

    BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_U(iLevel), BA(iLevel), DM(iLevel), &
               nDOFZ * nCR * ( iZ_E0( 1 ) - iZ_B0( 1 ) + 1 ) * nSpecies, swX )

      CALL amrex_multifab_build &
             ( MF_F(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX )

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), BA(iLevel), DM(iLevel), &
               nDOFX * nGS, 0 )
      CALL MF_uGS(iLevel) % SetVal( Zero ) ! remove this after debugging

      DO iS = 1, nStages

        CALL amrex_multifab_build &
               ( MF_DU_Ex(iLevel,iS), BA(iLevel), DM(iLevel), &
                 nDOFZ * nCR * ( iZ_E0( 1 ) - iZ_B0( 1 ) + 1 ) * nSpecies, swX )

        CALL amrex_multifab_build &
               ( MF_DU_Im(iLevel,iS), BA(iLevel), DM(iLevel), &
                 nDOFZ * nCR * ( iZ_E0( 1 ) - iZ_B0( 1 ) + 1 ) * nSpecies, swX )

        CALL amrex_multifab_build &
               ( MF_DF_Ex(iLevel,iS), BA(iLevel), DM(iLevel), &
                 nDOFX * nCF, swX )

        CALL amrex_multifab_build &
               ( MF_DF_Im(iLevel,iS), BA(iLevel), DM(iLevel), nDOFX * nCF, swX )

      END DO

    END DO

  END SUBROUTINE MF_InitializeField_IMEX_RK


  SUBROUTINE MF_FinalizeField_IMEX_RK

    INTEGER :: iLevel, iS

    CALL Finalize_IMEX_RK

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_U  (iLevel) )
      CALL amrex_multifab_destroy( MF_F  (iLevel) )
      CALL amrex_multifab_destroy( MF_uGS(iLevel) )

      DO iS = 1, nStages

        CALL amrex_multifab_destroy( MF_DU_Ex(iLevel,iS) )
        CALL amrex_multifab_destroy( MF_DU_Im(iLevel,iS) )
        CALL amrex_multifab_destroy( MF_DF_Ex(iLevel,iS) )
        CALL amrex_multifab_destroy( MF_DF_Im(iLevel,iS) )

      END DO

    END DO

    DEALLOCATE( MF_uGS )
    DEALLOCATE( MF_F )
    DEALLOCATE( MF_U )
    DEALLOCATE( MF_DU_Ex)
    DEALLOCATE( MF_DU_Im)
    DEALLOCATE( MF_DF_Ex)
    DEALLOCATE( MF_DF_Im)

  END SUBROUTINE MF_FinalizeField_IMEX_RK


  SUBROUTINE Initialize_IMEX_RK( Scheme, Verbose_Option )

    CHARACTER(LEN=*), INTENT(in) :: Scheme
    LOGICAL,               INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: i
    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

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

      CALL AllocateArray7D( StageData(i) % dU_IM )
      CALL AllocateArray7D( StageData(i) % dU_EX )
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

      CALL DeallocateArray7D( StageData(i) % dU_IM )
      CALL DeallocateArray7D( StageData(i) % dU_EX )
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
