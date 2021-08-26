MODULE InitializationModule

  ! --- AMReX Modules ---
  USE amrex_fort_module, ONLY: &
    AR => amrex_real, &
    amrex_spacedim
  USE amrex_amr_module, ONLY: &
    amrex_init, &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_init
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray,         &
    amrex_boxarray_build,   &
    amrex_boxarray_destroy, &
    amrex_print
  USE amrex_distromap_module, ONLY: &
    amrex_distromap,       &
    amrex_distromap_build, &
    amrex_distromap_destroy
  USE amrex_geometry_module, ONLY: &
    amrex_geometry, &
    amrex_geometry_build
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse,       &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---
  USE ProgramHeaderModule,              ONLY: &
    DescribeProgramHeaderX, &
    nDOFX,                  &
    nNodesX,                &
    nDOFE,                  &
    nNodesE,                &
    nDOFZ,                  &
    iZ_B0,                  &
    iZ_E0,                  &
    iZ_B1,                  &
    iZ_E1,                  &
    iE_B0,                  &
    iE_E0,                  &
    iE_B1,                  &
    iE_E1
  USE PolynomialBasisModuleX_Lagrange,  ONLY: &
    InitializePolynomialBasisX_Lagrange
  USE PolynomialBasisModuleX_Legendre,  ONLY: &
    InitializePolynomialBasisX_Legendre
  USE ReferenceElementModuleX,          ONLY: &
    InitializeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange
  USE ReferenceElementModule,           ONLY: &
    InitializeReferenceElement
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement
  USE EquationOfStateModule,            ONLY: &
    InitializeEquationOfState
  USE MeshModule,                       ONLY: &
    MeshX,      &
    MeshE,      &
    CreateMesh, &
    DestroyMesh
  USE RadiationFieldsModule,            ONLY: &
    nCR,      &
    nPR,      &
    CreateRadiationFields
  USE GeometryFieldsModule,             ONLY: &
    nGF,                     &
    CoordinateSystem,        &
    CreateGeometryFields
  USE GeometryFieldsModuleE, ONLY: &
    CreateGeometryFieldsE, &
    DestroyGeometryFieldsE, &
    uGE
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE FluidFieldsModule,                ONLY: &
    nCF,                     &
    nAF,                     &
    CreateFluidFields
  USE TwoMoment_OpacityModule_Relativistic,  ONLY: &
    CreateOpacities,         &
    SetOpacities
  USE OpacityModule_Table, ONLY:   &
    InitializeOpacities_TABLE
  USE PolynomialBasisMappingModule,     ONLY: &
    InitializePolynomialBasisMapping
  USE PolynomialBasisModule_Lagrange,   ONLY: &
    InitializePolynomialBasis_Lagrange
  USE PolynomialBasisModule_Legendre,   ONLY: &
    InitializePolynomialBasis_Legendre
  USE InputOutput,           ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint, &
    ReadCheckpointFile
  USE TwoMoment_PositivityLimiterModule_Relativistic, ONLY: &
    InitializePositivityLimiter_TwoMoment
  USE TwoMoment_SlopeLimiterModule_Relativistic, ONLY: &
    InitializeSlopeLimiter_TwoMoment
  ! --- Local modules ---
  USE MyAmrDataModule,                  ONLY: &
    MF_uGF, &
    MF_uCF, &
    MF_uAF, &
    MF_uPR, &
    MF_uCR
  USE MyAmrModule,                      ONLY: &
    t_end,                     &
    t,                         &
    dt,                        &
    t_wrt,                     &
    dt_wrt,                    &
    t_chk,                     &
    dt_chk,                    &
    CFL,                       &
    nNodes,                    &
    nStages,                   &
    nX,                        &
    nE,                        &
    swX,                       &
    swE,                       &
    eR,                        &
    eL,                        &
    zoomE,                     &
    bcX,                       &
    xL,                        &
    xR,                        &
    ProgramName,               &
    Scheme,                    &
    CoordSys,                  &
    StepNo,                    &
    nLevels,                   &
    iRestart,                  &
    MaxGridSizeX,              &
    BA,                        &
    DM,                        &
    GEOM,                      &
    StepNo,                    &
    nSpecies,                  &
    V_0,                       &
    Gamma_IDEAL,               &
    D_0,                       &
    Chi,                       &
    Sigma,                     &
    kT,                        &
    E0,                        &
    mu0,                       &
    R0,                        &
    EquationOfState,           &
    EosTableName,              &
    OpacityTableName_AbEm,     &
    OpacityTableName_Iso,     &
    OpacityTableName_NES,     &
    OpacityTableName_Pair,     &
    Min_1,                     &
    Min_2,                     &
    UsePositivityLimiter,      &
    UseSlopeLimiter,      &
    BetaTVD,      &
    MyAmrInit
  USE MF_InitializationModule,          ONLY: &
    MF_InitializeFields
  USE MF_GeometryModule,                ONLY: &
    MF_ComputeGeometryX
  USE MF_TwoMoment_TimeSteppingModule_Relativistic,  ONLY: &
    MF_InitializeField_IMEX_RK
  USE TwoMoment_ClosureModule,                       ONLY: &
    InitializeClosure_TwoMoment
  USE MF_TwoMoment_UtilitiesModule,     ONLY: & 
    MF_ComputeFromConserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

  LOGICAL, PUBLIC :: wrt, chk

  REAL(AR), PARAMETER :: Zero = 0.0_AR
  REAL(AR), PARAMETER :: One  = 1.0_AR
  REAL(AR), PARAMETER :: Two  = 2.0_AR

CONTAINS


  SUBROUTINE InitializeProgram

    INTEGER               :: iLevel, iDim
    TYPE(amrex_parmparse) :: PP
    TYPE(amrex_box)       :: BX
    REAL(AR)              :: Mass, W, Vad

    ! --- Initialize AMReX ---
    CALL amrex_init()

    CALL amrex_amrcore_init()

    ! --- Parse parameter file ---
    CALL MyAmrInit

    IF( iRestart .LT. 0 )THEN

      BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

      ALLOCATE( BA(0:nLevels-1) )

      DO iLevel = 0, nLevels-1

        CALL amrex_boxarray_build( BA(iLevel), BX )

      END DO

      DO iLevel = 0, nLevels-1

        CALL BA(iLevel) % maxSize( MaxGridSizeX )

      END DO

      ALLOCATE( GEOM(0:nLevels-1) )

      ALLOCATE( DM  (0:nLevels-1) )

      DO iLevel = 0, nLevels-1

        CALL amrex_geometry_build ( GEOM(iLevel), BX )
        CALL amrex_distromap_build( DM  (iLevel), BA(iLevel) )

      END DO


      DO iLevel = 0, nLevels-1
        CALL amrex_multifab_build &
               ( MF_uGF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nGF, swX )
        CALL MF_uGF(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uCF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX )
        CALL MF_uCF(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uAF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nAF, swX )
        CALL MF_uAF(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uPR(iLevel), BA(iLevel), DM(iLevel), &
                 nDOFZ * nPR * ( iZ_E0( 1 ) - iZ_B0( 1 ) + 1 ) * nSpecies, swX )

        CALL MF_uPR(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uCR(iLevel), BA(iLevel), DM(iLevel), &
                 nDOFZ * nCR * ( iZ_E0( 1 ) - iZ_B0( 1 ) + 1 ) * nSpecies, swX )

        CALL MF_uCR(iLevel) % SetVal( Zero )

      END DO

      t     = Zero
      dt    = Zero
      t_wrt = dt_wrt
      t_chk = dt_chk

    ELSE

      CALL ReadCheckpointFile( iRestart )

    END IF

    wrt = .FALSE.
    chk = .FALSE.


     DO iDim = 1, 3
      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
               xL(iDim), xR(iDim) )
    END DO

    CALL InitializeReferenceElementX

    CALL InitializeReferenceElementE

    CALL InitializeReferenceElement

    CALL InitializePolynomialBasisX_Lagrange
    CALL InitializePolynomialBasisX_Legendre

    CALL InitializePolynomialBasis_Lagrange
    CALL InitializePolynomialBasis_Legendre

    CALL InitializeReferenceElementX_Lagrange
    CALL InitializeReferenceElement_Lagrange

    CALL CreateRadiationFields( nX, swX, nE, swE, nSpecies_Option = nSpecies, &
                                Verbose_Option = amrex_parallel_ioprocessor()  )

!    CALL CreateFluidFields( nX, swX, CoordinateSystem_Option = 'CARTESIAN', &
!                              Verbose_Option = amrex_parallel_ioprocessor()  )
!
!    CALL CreateGeometryFields( nX, swX, CoordinateSystem_Option = 'CARTESIAN', &
!                               Verbose_Option = amrex_parallel_ioprocessor()  )
    
    


    CALL CreateFluidFields( nX, swX, CoordinateSystem_Option = 'SPHERICAL', &
                              Verbose_Option = amrex_parallel_ioprocessor()  )

    CALL CreateGeometryFields( nX, swX, CoordinateSystem_Option = 'SPHERICAL', &
                               Verbose_Option = amrex_parallel_ioprocessor()  )

    CALL CreateMesh &
           ( MeshE, nE, nNodesE, swE, eL, eR, zoomOption = zoomE )

    CALL CreateGeometryFieldsE &
           ( nE, swE, Verbose_Option = .FALSE. )

    CALL ComputeGeometryE &
           ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

    CALL MF_ComputeGeometryX( MF_uGF, 0.0_AR )


    #ifdef MICROPHYSICS_WEAKLIB

      CALL InitializeEquationOfState &
             ( EquationOfState_Option = EquationOfState, &
               EquationOfStateTableName_Option &
                 = EosTableName )

      CALL InitializeOpacities_TABLE &
        ( OpacityTableName_EmAb_Option = OpacityTableName_AbEm, &
          OpacityTableName_Iso_Option  = OpacityTableName_Iso,  &
          OpacityTableName_NES_Option  = OpacityTableName_NES,  &
          OpacityTableName_Pair_Option = OpacityTableName_Pair, &
          EquationOfStateTableName_Option = EosTableName, &
          Verbose_Option =  amrex_parallel_ioprocessor())

    #else


      CALL InitializeEquationOfState &
             ( EquationOfState_Option = EquationOfState, &
               Gamma_IDEAL_Option = Gamma_IDEAL, &
               Verbose_Option = amrex_parallel_ioprocessor()  )

      CALL CreateOpacities &
         ( nX, [ 1, 1, 1 ], nE, 1, Verbose_Option = amrex_parallel_ioprocessor() )

      CALL SetOpacities( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D_0, Chi, Sigma, kT, E0, mu0, R0, & 
                       Verbose_Option = amrex_parallel_ioprocessor()  )

    #endif

    CALL MF_InitializeFields( TRIM( ProgramName ), MF_uGF, MF_uCR, MF_uCF, V_0, &
                              Verbose_Option = amrex_parallel_ioprocessor() )
    DO iLevel = 0, nLevels-1

      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

    END DO

    CALL MF_InitializeField_IMEX_RK( Scheme, BA, DM, &
                                     Verbose_Option = amrex_parallel_ioprocessor() )


    CALL InitializeClosure_TwoMoment

    CALL InitializePositivityLimiter_TwoMoment &
         ( Min_1_Option = Min_1, &
           Min_2_Option = Min_2, &
           UsePositivityLimiter_Option &
             = UsePositivityLimiter, &
           Verbose_Option = amrex_parallel_ioprocessor() )

    CALL InitializeSlopeLimiter_TwoMoment &
           ( BetaTVD_Option &
               = BetaTVD, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter, &
             Verbose_Option &
               = amrex_parallel_ioprocessor()  )

    CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

    CALL WriteFieldsAMReX_PlotFile &
           ( t(0), StepNo, &
             MF_uCR_Option = MF_uCR, &
             MF_uPR_Option = MF_uPR )
!  IF (ProgramName == "SineWaveStreaming") THEN
!    W = 1.0_AR - DOT_PRODUCT(V_0,V_0)
!    W = 1.0_AR / SQRT(W)
!    Vad = (1.0_AR + V_0(1) * W ) / ( W + V_0(1) )
!    t_end = t_end / Vad
!  END IF
! 

  END SUBROUTINE InitializeProgram


END MODULE InitializationModule
