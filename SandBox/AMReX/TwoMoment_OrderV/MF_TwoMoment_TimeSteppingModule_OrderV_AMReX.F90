MODULE MF_TwoMoment_TimeSteppingModule_OrderV_AMReX

  ! --- AMReX Modules ---

  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_fort_module,      ONLY: &
    AR => amrex_real
  USE amrex_box_module,       ONLY: &
    amrex_box
  USE amrex_geometry_module,  ONLY: &
    amrex_geometry
  USE amrex_multifab_module,  ONLY: &
    amrex_multifab, amrex_imultifab, &
    amrex_multifab_build, amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, amrex_mfiter_destroy
  USE amrex_boxarray_module,  ONLY: &
    amrex_boxarray, &
    amrex_boxarray_build, amrex_boxarray_destroy
  USE amrex_distromap_module, ONLY: &
    amrex_distromap, &
    amrex_distromap_build, amrex_distromap_destroy
  USE amrex_amrcore_module, ONLY: &
    GEOM => amrex_geom
  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE, nE, &
    iE_B0, iE_E0, &
    iE_B1, iE_E1, &
    iZ_B0, iZ_B1, iZ_E0, iZ_E1, swX, nX, nDimsX, nNodes
  USE GeometryFieldsModuleE, ONLY: &
    nGE, uGE
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nCR
  USE TwoMoment_DiscretizationModule_Collisions, ONLY: &
    ComputeIncrement_TwoMoment_Implicit
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
  USE thornado_amrex_timeintegrator_module, ONLY: &
    amrex_integrate_twomoment
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE FillPatchModule, ONLY: &
    FillPatch
  USE MeshModule, ONLY: &
    MeshX, MeshE


  ! --- Local Modules ---

  USE InputParsingModule,                      ONLY: &
    nLevels, DEBUG, UseTiling, nSpecies, UseFluxCorrection_TwoMoment
  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uDF, &
    OffGridFlux_Euler_MF
  USE MF_TwoMoment_DiscretizationModule_Streaming_OrderV, ONLY: &
    ComputeIncrement_TwoMoment_Explicit_MF
  USE MF_TwoMoment_DiscretizationModule_Collisions_OrderV, ONLY: &
    ComputeIncrement_TwoMoment_Implicit_MF
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment_MF
  USE MF_TwoMoment_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_TwoMoment_MF
  USE MF_FieldsModule_TwoMoment, ONLY: &
    OffGridFlux_TwoMoment_MF, &
    MF_uCR
  USE InputParsingModule, ONLY: &
    t_new, &
    dt, &
    nMaxLevels
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler_MF
  USE MF_Euler_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_Euler_MF
  USE MF_GravitySolutionModule, ONLY: &
    EvolveGravity
  USE AverageDownModule, ONLY: &
    AverageDown
  USE MF_TimersModule, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_UpdateFluid
  USE MF_Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_MF
  USE MF_UtilitiesModule
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MaskModule, ONLY: &
    CreateFineMask, DestroyFineMask, IsNotLeafElement

  USE FluxCorrectionModule_Euler, ONLY: &
    ApplyFluxCorrection_Euler_MF
  USE FluxCorrectionModule_TwoMoment, ONLY: &
    ApplyFluxCorrection_TwoMoment_MF
  USE MF_TwoMoment_TallyModule

  USE ISO_C_BINDING

  IMPLICIT NONE
  PRIVATE


  TYPE :: StageDataType
    REAL(AR), ALLOCATABLE :: dU_IM(:,:,:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: dU_EX(:,:,:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: dF_EX(:,:,:,:,:)
  END TYPE StageDataType

  INTEGER               :: nStages
  REAL(DP), ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(DP), ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  LOGICAL               :: EvolveEuler
  LOGICAL               :: EvolveTwoMoment
  LOGICAL               :: Verbose
  REAL(AR),            ALLOCATABLE :: U0(:,:,:,:,:,:,:)
  REAL(AR),            ALLOCATABLE :: Ui(:,:,:,:,:,:,:)
  TYPE(StageDataType), ALLOCATABLE :: StageData(:)

  REAL(DP)                 , PUBLIC :: CFL
  CHARACTER(:), ALLOCATABLE, PUBLIC :: Scheme

  ! FIXME
  PUBLIC :: Initialize_IMEX_RK_MF
  PUBLIC :: Finalize_IMEX_RK_MF
  PUBLIC :: Update_AMReX_TimeIntegrator

CONTAINS



  SUBROUTINE Initialize_IMEX_RK_MF &
    (BA, DM , Verbose_Option )

    LOGICAL         , INTENT(in), OPTIONAL :: Verbose_Option
    TYPE(amrex_boxarray),  INTENT(in)           :: BA(0:nLevels-1)
    TYPE(amrex_distromap), INTENT(in)           :: DM(0:nLevels-1)
    TYPE(amrex_parmparse) :: PP

    INTEGER         :: iS, iLevel

    EvolveEuler     = .TRUE.
    EvolveTwoMoment = .TRUE.
    CALL amrex_parmparse_build( PP, 'TS' )
      CALL PP % get  ( 'Scheme'         , Scheme )
      CALL PP % get  ( 'CFL'            , CFL )
      CALL PP % query( 'EvolveEuler'    , EvolveEuler )
      CALL PP % query( 'EvolveTwoMoment', EvolveTwoMoment )
    CALL amrex_parmparse_destroy( PP )

    CFL = CFL / ( DBLE( nDimsX ) * ( Two * DBLE( nNodes ) - One ) )

    !PRINT *, 'CFL: ', CFL
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
      WRITE(*,'(6x,A18,A)') 'IMEX-RK Scheme: ', TRIM( Scheme )
      WRITE(*,'(A5,A,ES10.3E3)') '', 'CFL:           ', &
        CFL * ( DBLE( nDimsX ) * ( Two * DBLE( nNodes ) - One ) )
      WRITE(*,'(6x,A18,L)') 'EvolveEuler: ', EvolveEuler
      WRITE(*,'(6x,A18,L)') 'EvolveTwoMoment: ', EvolveTwoMoment

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

  SUBROUTINE Update_AMReX_TimeIntegrator

    !! ! nCompCF, nCompCR: number of components (fluid and radiation),
    !! ! iLevel: index for refinement level
    INTEGER :: nCompCF, nCompCR, iLevel

    ! What is the purpose of these?
    REAL(DP) :: dM_OffGrid_Euler    (1:nCF  ,0:nMaxLevels-1)
    REAL(DP) :: dM_OffGrid_TwoMoment(1:2*nCR,0:nMaxLevels-1)

    TYPE(amrex_imultifab) :: iMF_FineMask

    CALL TimersStart_AMReX( Timer_AMReX_UpdateFluid )
    dM_OffGrid_Euler     = Zero
    dM_OffGrid_TwoMoment = Zero

    ! nCompGF = nDOFX * nGF
    ! nCompCF = nDOFX * nCF
    ! nCompCR = nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies!nDOFZ * nCR * nE * nSpecies

    !PRINT *, DEBUG

    IF( DEBUG ) WRITE(*,'(A)') 'Entering Update_AMReX_TimeIntegrator'

    !TODO: Call AMReX time integrator
    !TODO: Add options for time integrator, e.g. multirate or adaptive stepsize


    IF( EvolveEuler )THEN
      ! TODO: Implement
      ! CALL amrex_integrate_euler &
      !      ( GEOM, MF_uGF, MF_uCF, MF_uCR, t_new, dt, &
      !        Verbose_Option = Verbose )
    END IF ! EvolveEuler

    IF( EvolveTwoMoment )THEN
      DO iLevel = 0, nLevels-1

        ! --- Apply boundary conditions to interior domains ---

        CALL FillPatch &
             ( iLevel, MF_uGF,MF_uCR, ApplyBoundaryConditions_TwoMoment_Option = .TRUE.)

        CALL FillPatch( iLevel, MF_uGF, MF_uCF, ApplyBoundaryConditions_Euler_Option = .TRUE.)


        CALL FillPatch &
             ( iLevel, MF_uGF, &
               ApplyBoundaryConditions_Geometry_Option = .TRUE. )

        CALL ApplyPositivityLimiter_TwoMoment_MF &
               ( iLevel, MF_uGF(iLevel), MF_uCF(iLevel), MF_uCR(iLevel) )

        ! CALL MF_duCR(iLevel) % setval( 0.0_amrex_real )

        CALL CreateMesh_MF( iLevel, MeshX )

        CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

        CALL amrex_integrate_twomoment &
             ( GEOM(iLevel), uGE, MF_uGF(iLevel), MF_uCF(iLevel), MF_uCR(iLevel), &
               t_new(iLevel), dt(iLevel), &
               C_FUNLOC(ComputeIncrement_TwoMoment_Implicit), &
               C_FUNLOC(ComputeIncrement_TwoMoment_Explicit), &
               C_FUNLOC(ApplyBoundaryConditions_TwoMoment), &
               C_FUNLOC(ApplyBoundaryConditions_Euler), &
               swX, nGF, nCF, nCR, nDOFX, nDOFE, nSpecies, &
               iE_E0, iE_B0, iE_E1, iE_B1)
      END DO
    END IF ! EvolveTwoMoment

    CALL DestroyMesh_MF( MeshX )
    CALL DestroyFineMask( iMF_FineMask )

    !! IF( EvolveEuler )THEN
    !!   CALL ApplyPositivityLimiter_Euler_MF(MF_uGF, MF_uCF, MF_uDF)
    !! END IF ! EvolveEuler

    !! IF( EvolveTwoMoment )THEN
    !!   CALL ApplyPositivityLimiter_TwoMoment_MF(MF_uGF, MF_uCF, MF_uCR)
    !! END IF ! EvolveTwoMoment

    !! ! --- Explicit step ---

    !! IF( EvolveEuler )THEN
    !!   CALL ComputeIncrement_Euler_MF(MF_uGF, MF_uCF, MF_uDF, MF_DF_EX(:,iS))
    !! END IF ! EvolveEuler

    !! IF( EvolveTwoMoment )THEN
    !!   CALL ComputeIncrement_TwoMoment_Explicit_MF &
    !!          ( t_new, GEOM, MF_uGF, MF_uCF, MF_uCR, MF_DR_EX(:,iS), &
    !!            Verbose_Option = .FALSE. )

    !!   DO iLevel = 0, nLevels-1
    !!      ! What does this do??
    !!      dM_OffGrid_TwoMoment(:,iLevel) &
    !!         = dM_OffGrid_TwoMoment(:,iLevel) &
    !!         + dt(iLevel) * w_EX(iS) * OffGridFlux_TwoMoment_MF(:,iLevel)
    !!   END DO !DO iLevel = 0, nLevels-1

    !!   IF( nLevels .GT. 1 .AND. UseFluxCorrection_TwoMoment )THEN
    !!     CALL ApplyFluxCorrection_TwoMoment_MF( MF_uGF, MF_DR_EX(:,iS) )
    !!   END IF

    !! END IF ! EvolveTwoMoment

    !! CALL IncrementOffGridTally_TwoMoment_MF( dM_OffGrid_TwoMoment )

    !! IF( nLevels .GT. 1 ) &
    !!   CALL AverageDown( MF_uGF, MF_uCR, Radiation_Option = .TRUE. )

    !! IF( DEBUG ) WRITE(*,'(A)') 'Leaving Update_AMReX_TimeIntegrator'

    !! CALL TimersStop_AMReX( Timer_AMReX_UpdateFluid )

  END SUBROUTINE Update_AMReX_TimeIntegrator

END MODULE MF_TwoMoment_TimeSteppingModule_OrderV_AMReX
