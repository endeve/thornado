  MODULE MF_TwoMoment_PositivityLimiterModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDOFX, &
    nDOFZ, &
    iE_B0, &
    iE_E0, &
    iE_B1, &
    iE_E1
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    nSpecies
  USE FluidFieldsModule, ONLY: &
    nCF
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment, &
    ApplyPositivityLimiter_TwoMoment
  USE MeshModule, ONLY: &
    MeshX
  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    amrex2thornado_Z, &
    thornado2amrex_Z, &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    nE, &
    DEBUG
  USE MF_EdgeMapModule, ONLY: &
    ConstructEdgeMap, &
    EdgeMap
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler_MF
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_TwoMoment_MF
  PUBLIC :: FinalizePositivityLimiter_TwoMoment_MF
  PUBLIC :: ApplyPositivityLimiter_TwoMoment_MF

  INTERFACE ApplyPositivityLimiter_TwoMoment_MF
    MODULE PROCEDURE ApplyPositivityLimiter_TwoMoment_MF_MultipleLevels
    MODULE PROCEDURE ApplyPositivityLimiter_TwoMoment_MF_SingleLevel
  END INTERFACE ApplyPositivityLimiter_TwoMoment_MF

  LOGICAL :: UsePositivityLimiter, UseEnergyLimiter

CONTAINS


  SUBROUTINE InitializePositivityLimiter_TwoMoment_MF

    TYPE(amrex_parmparse) :: PP

    REAL(DP) :: Min_1, Min_2

    UsePositivityLimiter = .TRUE.
    UseEnergyLimiter     = .FALSE.
    Min_1                = 1.0e-12_DP
    Min_2                = 1.0e-12_DP
    CALL amrex_parmparse_build( PP, 'PL' )
      CALL PP % query( 'UsePositivityLimiter_TwoMoment', &
                        UsePositivityLimiter )
      CALL PP % query( 'UseEnergyLimiter_TwoMoment', &
                        UseEnergyLimiter )
      CALL PP % query( 'Min_1_TwoMoment', &
                        Min_1 )
      CALL PP % query( 'Min_2_TwoMoment', &
                        Min_2 )
    CALL amrex_parmparse_destroy( PP )

    CALL InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option = Min_1, &
             Min_2_Option = Min_2, &
             UsePositivityLimiter_Option &
               = UsePositivityLimiter, &
             UseEnergyLimiter_Option &
               = UseEnergyLimiter, &
             Verbose_Option = amrex_parallel_ioprocessor() )

  END SUBROUTINE InitializePositivityLimiter_TwoMoment_MF


  SUBROUTINE FinalizePositivityLimiter_TwoMoment_MF

    CALL FinalizePositivityLimiter_TwoMoment

  END SUBROUTINE FinalizePositivityLimiter_TwoMoment_MF


  SUBROUTINE ApplyPositivityLimiter_TwoMoment_MF_SingleLevel &
    ( iLevel, MF_uGF, MF_uCF, MF_uCR, swX_Option )

    INTEGER             , INTENT(in)           :: iLevel
    TYPE(amrex_multifab), INTENT(in)           :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout)        :: MF_uCF
    TYPE(amrex_multifab), INTENT(inout)        :: MF_uCR
    INTEGER             , INTENT(in), OPTIONAL :: swX_Option(3)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR (:,:,:,:)

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: C (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U (:,:,:,:,:,:,:)

    INTEGER       :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX_B(3), iX_E(3), &
                     iLo_MF(4), swXX(3)
    INTEGER       :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), iZ_B(4), iZ_E(4)
    TYPE(EdgeMap) :: Edge_Map

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    swXX = 0
    IF( PRESENT( swX_Option ) ) &
      swXX = swX_Option

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

      ! --- Apply boundary conditions to interior domains ---

      DO WHILE( MFI % next() )

        uGF  => MF_uGF % DataPtr( MFI )
        uCF  => MF_uCF% DataPtr( MFI )
        uCR  => MF_uCR% DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX
        iX_B  = BX % lo - swXX
        iX_E  = BX % hi + swXX

        iZ_B0(1) = iE_B0
        iZ_B1(1) = iE_B1
        iZ_E0(1) = iE_E0
        iZ_E1(1) = iE_E1

        iZ_B0(2:4) = iX_B0
        iZ_B1(2:4) = iX_B1
        iZ_E0(2:4) = iX_E0
        iZ_E1(2:4) = iX_E1

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 C )

        CALL AllocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nCR     , &
                   nSpecies ], &
                 U )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uCF, C )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uCR, U )

      ! --- Apply boundary conditions to physical boundaries
      !     (needed for AMR) ---

        CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

        !CALL ApplyBoundaryConditions_Euler_MF &
             !( iX_B0, iX_E0, iX_B1, iX_E1, C, Edge_Map )

        CALL ApplyPositivityLimiter_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, G, C, U )

        CALL thornado2amrex_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uCR, U )

        CALL DeallocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nCR     , &
                   nSpecies ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 C )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! DO WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ApplyPositivityLimiter_TwoMoment_MF_SingleLevel


  SUBROUTINE ApplyPositivityLimiter_TwoMoment_MF_MultipleLevels &
    ( MF_uGF, MF_uCF, MF_uCR, swX_Option )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
    INTEGER             , INTENT(in), OPTIONAL :: swX_Option(3)

    INTEGER :: iLevel, iErr, swXX(3)

    swXX = 0
    IF( PRESENT( swX_Option ) ) &
      swXX = swX_Option

    DO iLevel = 0, nLevels-1

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(2x,A,I3.3)') &
            'CALL ApplyPositivityLimiter_TwoMoment_MF_SingleLevel, iLevel: ', &
            iLevel

        END IF

      END IF ! DEBUG

      CALL ApplyPositivityLimiter_TwoMoment_MF_SingleLevel &
             ( iLevel, MF_uGF(iLevel), MF_uCF(iLevel), MF_uCR(iLevel), &
               swX_Option = swXX )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ApplyPositivityLimiter_TwoMoment_MF_MultipleLevels


END MODULE MF_TwoMoment_PositivityLimiterModule
