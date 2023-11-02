MODULE MF_Euler_PerturbationModule

  ! ---AMReX Modules ---

  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse,       &
    amrex_parmparse_build, &
    amrex_parmparse_destroy


  ! --- thornado Modules ---

  USE ProgramHeaderModule,      ONLY: &
    swX, &
    nDOFX
  USE FluidFieldsModule,        ONLY: &
    nCF, &
    nDF
  USE GeometryFieldsModule,     ONLY: &
    nGF
  USE Euler_PerturbationModule, ONLY: &
    InitializeShellPerturbations, &
    ApplyShellPerturbations

  ! --- Local Modules ---

  USE MF_KindModule,                     ONLY: &
    DP
  USE MF_UtilitiesModule,                ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X
  USE InputParsingModule,                ONLY: &
    nLevels,          &
    DEBUG,            &
    UsePhysicalUnits, &
    UseTiling,        &
    GEOM

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializePerturbations_Euler, MF_ApplyPerturbations_Euler


CONTAINS


  SUBROUTINE MF_InitializePerturbations_Euler


    CHARACTER(LEN=:), ALLOCATABLE :: PerturbType

    INTEGER  :: NumShells
    INTEGER  :: BD_Order
    INTEGER  :: CO_V_Nodes
    INTEGER  :: MJ_D_Nodes
    INTEGER  :: MJ_V_RadNodes
    INTEGER  :: MJ_V_AngNodes
    REAL(DP) :: ShellWidth
    REAL(DP) :: BD_Amplitude
    REAL(DP) :: CO_V_Amplitude
    REAL(DP) :: MJ_D_Amplitude
    REAL(DP) :: MJ_V_Amplitude

    TYPE(amrex_parmparse) :: PP

    PerturbType = 'None'

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'PerturbType', PerturbType )
    CALL amrex_parmparse_destroy( PP )

    SELECT CASE( TRIM( PerturbType ) )

      CASE( 'None' )

        CALL InitializeShellPerturbations &
               ( PerturbType_Option = 'None' )

      CASE( 'Shell_BasicDensity' )

        CALL amrex_parmparse_build( PP, 'BD' )
          CALL PP % get     ( 'NumShells',   NumShells )
          CALL PP % get     ( 'ShellWidth',  ShellWidth )
          CALL PP % get     ( 'Order',       BD_Order )
          CALL PP % get     ( 'Amplitude',   BD_Amplitude )
        CALL amrex_parmparse_destroy( PP )

        CALL InitializeShellPerturbations &
               ( PerturbType_Option  = 'BasicDensity',    &
                 NumShells_Option    = NumShells,         &
                 ShellWidth_Option   = ShellWidth,        &
                 BD_Order_Option     = BD_Order,          &
                 BD_Amplitude_Option = BD_Amplitude,      &
                 UseUnits_Option     = UsePhysicalUnits )

      CASE( 'Shell_CouchOtt_Velocity' )

        CALL amrex_parmparse_build( PP, 'CO_V' )
          CALL PP % get   ( 'NumShells',  NumShells )
          CALL PP % get   ( 'ShellWidth', ShellWidth )
          CALL PP % get   ( 'Nodes',      CO_V_Nodes )
          CALL PP % get   ( 'Amplitude',  CO_V_Amplitude )
        CALL amrex_parmparse_destroy( PP )

        CALL InitializeShellPerturbations &
               ( PerturbType_Option    = 'CouchOtt_Velocity', &
                 NumShells_Option      = NumShells,           &
                 ShellWidth_Option     = ShellWidth,          &
                 CO_V_Nodes_Option     = CO_V_Nodes,          &
                 CO_V_Amplitude_Option = CO_V_Amplitude,      &
                 UseUnits_Option       = UsePhysicalUnits )

      CASE( 'Shell_MullerJanka_Density' )

        CALL amrex_parmparse_build( PP, 'MJ_D' )
          CALL PP % get     ( 'NumShells',   NumShells )
          CALL PP % get     ( 'ShellWidth',  ShellWidth )
          CALL PP % get     ( 'Nodes',       MJ_D_Nodes )
          CALL PP % get     ( 'Amplitude',   MJ_D_Amplitude )
        CALL amrex_parmparse_destroy( PP )

        CALL InitializeShellPerturbations &
               ( PerturbType_Option    = 'MullerJanka_Density', &
                 NumShells_Option      = NumShells,             &
                 ShellWidth_Option     = ShellWidth,            &
                 MJ_D_Nodes_Option     = MJ_D_Nodes,            &
                 MJ_D_Amplitude_Option = MJ_D_Amplitude,        &
                 UseUnits_Option       = UsePhysicalUnits )

      CASE( 'Shell_MullerJanka_Velocity' )

        CALL amrex_parmparse_build( PP, 'MJ_V' )
          CALL PP % get     ( 'NumShells',   NumShells )
          CALL PP % get     ( 'ShellWidth',  ShellWidth )
          CALL PP % get     ( 'RadNodes',    MJ_V_RadNodes )
          CALL PP % get     ( 'AngNodes',    MJ_V_AngNodes )
          CALL PP % get     ( 'Amplitude',   MJ_V_Amplitude )
        CALL amrex_parmparse_destroy( PP )

        CALL InitializeShellPerturbations &
               ( PerturbType_Option    = 'MullerJanka_Velocity', &
                 NumShells_Option      = NumShells,              &
                 ShellWidth_Option     = ShellWidth,             &
                 MJ_V_RadNodes_Option  = MJ_V_RadNodes,          &
                 MJ_V_AngNodes_Option  = MJ_V_AngNodes,          &
                 MJ_V_Amplitude_Option = MJ_V_Amplitude,         &
                 UseUnits_Option       = UsePhysicalUnits )

      CASE( 'MuellerJanka_RandVelocity' )

        CALL amrex_parmparse_build( PP )
        CALL amrex_parmparse_destroy( PP )

      CASE DEFAULT

        CALL InitializeShellPerturbations &
               ( PerturbType_Option = 'None' )

     END SELECT

  END SUBROUTINE MF_InitializePerturbations_Euler


  SUBROUTINE MF_ApplyPerturbations_Euler( Time, MF_uGF, MF_uCF, GEOM )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in   ) :: GEOM  (0:nLevels-1)

    REAL(DP)            , INTENT(in   ) :: Time(0:nLevels-1)

    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    INTEGER  :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF) )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

        IF( DEBUG ) WRITE(*,'(A)') '    CALL ApplyPerturbations_Euler'

        IF( iX_E1(1) .GT. GEOM(iLevel) % DOMAIN % hi(1) )THEN

          CALL ApplyShellPerturbations &
                 ( Time(iLevel), iX_B0, iX_E0, iX_B1, iX_E1, G, U )

        END IF

        CALL thornado2amrex_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

        DEALLOCATE( U )

        DEALLOCATE( G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ApplyPerturbations_Euler


END MODULE MF_Euler_PerturbationModule
