MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE UnitsModule, ONLY: &
    Kilometer, &
    Gram, &
    Centimeter, &
    Erg
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    NodeCoordinate, &
    MeshX
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    swX
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nCF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nPF
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    nGF
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL, &
    ComputePressureFromPrimitive_IDEAL
  USE Euler_UtilitiesModule, ONLY: &
    ComputeConserved_Euler

  ! --- Local Modules ---

  USE MF_UtilitiesModule, ONLY: &
    AllocateArray_X, &
    DeallocateArray_X, &
    amrex2thornado_X, &
    thornado2amrex_X
  USE MF_KindModule, ONLY: &
    DP, &
    One, &
    Zero, &
    Pi
  USE MF_EdgeMapModule, ONLY: &
    EdgeMap, &
    ConstructEdgeMap
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler_MF
  USE InputParsingModule, ONLY: &
    UseTiling, &
    xR

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF, MF_uPF, MF_uAF

    CALL InitializeFields_Polytrope_Newtonian_test( iLevel, MF_uGF, MF_uCF )

  END SUBROUTINE InitializeFields_MF


SUBROUTINE InitializeFields_Polytrope_Newtonian_test( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX_B(3), iX_E(3)
    REAL(DP) :: uPF(nDOFX,nPF)
    REAL(DP) :: Pressure

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    TYPE(EdgeMap) :: Edge_Map
    TYPE(amrex_parmparse) :: PP


    ! --- Problem-specific Parameters ---

    INTEGER  :: iNX1
    REAL(DP) :: X1
    REAL(DP) :: rho_c
    REAL(DP) :: P_c
    REAL(DP) :: K_polytrope, R_neutron_star

    CALL amrex_parmparse_build( PP, 'inputs' )
      CALL PP % get( 'rho_c'      , rho_c )
      CALL PP % get('R_neutron_star', R_neutron_star)
    CALL amrex_parmparse_destroy( PP )

    rho_c       = rho_c * ( Gram / Centimeter**3 )
    R_neutron_star = R_neutron_star * Kilometer
    K_polytrope = ( 2.0_DP / Pi ) * R_neutron_star**2
    P_c = K_polytrope * rho_c**Gamma_IDEAL

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(6x,A)')   'Poseidon Newtonian Polytrope - test'
      WRITE(*,'(6x,A,A)') '-----------------'
      WRITE(*,*)
      WRITE(*,'(8x,A,ES15.6E3)') 'rho_c in g/cc: ', rho_c / ( Gram / Centimeter**3 )

      WRITE(*,*)

    END IF

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      CALL amrex2thornado_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B0, iX_E0, uGF, G )

      iX_B = BX % lo
      iX_E = BX % hi

      IF( BX % hi(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) ) &
          iX_E(1) = iX_E(1) + swX(1)

      DO iX3 = iX_B(3), iX_E(3)
      DO iX2 = iX_B(2), iX_E(2)
      DO iX1 = iX_B(1), iX_E(1)
      DO iNX = 1      , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

  ! --- Density profile for n=1 polytrope, lane Emden Solution --- !

	      uPF(iNX,iPF_D) = rho_c * ( ( R_neutron_star / ( Pi * X1 ) ) * ( SIN( Pi * X1 / R_neutron_star ) ) )
        PRINT "(A,ES25.16E3)", 'uPF(iNX,iPF_D)', uPF(iNX,iPF_D)/ ( Gram / Centimeter**3 )
	      uPF(iNX,iPF_V1) = Zero
	      uPF(iNX,iPF_V2) = Zero
        uPF(iNX,iPF_V3) = Zero
        uPF(iNX,iPF_E ) = ( K_polytrope * uPF(iNX,iPF_D) ** Gamma_IDEAL ) / ( Gamma_IDEAL - One )
        PRINT *,'uPF(iNX,iPF_E)', uPF(iNX,iPF_E)
        uPF(iNX,iPF_Ne) = Zero

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPF(iNX,iPF_D ), uPF(iNX,iPF_E), &
                 uPF(iNX,iPF_Ne), Pressure )

        CALL ComputeConserved_Euler &
               ( uPF(iNX,iPF_D ), &
                 uPF(iNX,iPF_V1), &
                 uPF(iNX,iPF_V2), &
                 uPF(iNX,iPF_V3), &
                 uPF(iNX,iPF_E ), &
                 uPF(iNX,iPF_Ne), &
                 U  (iNX,iX1,iX2,iX3,iCF_D ), &
                 U  (iNX,iX1,iX2,iX3,iCF_S1), &
                 U  (iNX,iX1,iX2,iX3,iCF_S2), &
                 U  (iNX,iX1,iX2,iX3,iCF_S3), &
                 U  (iNX,iX1,iX2,iX3,iCF_E ), &
                 U  (iNX,iX1,iX2,iX3,iCF_Ne), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G  (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 Pressure )

      END DO
      END DO
      END DO
      END DO

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X &
             ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_Polytrope_Newtonian_test


END MODULE MF_InitializationModule
