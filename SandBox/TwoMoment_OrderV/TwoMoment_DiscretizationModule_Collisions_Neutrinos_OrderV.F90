MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_T, iAF_E , iAF_Ye
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Auxiliary_TABLE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit
  PUBLIC :: Initialize_TwoMoment_Collisions_Neutrinos
  PUBLIC :: Finalize_TwoMoment_Collisions_Neutrinos

  INTEGER               :: nE_G, nX_G
  INTEGER               :: nZ(4), nX(3), nE
  INTEGER               :: iE_B0, iE_E0, iX_B0(3), iX_E0(3)
  INTEGER               :: iE_B1, iE_E1, iX_B1(3), iX_E1(3)
  REAL(DP), ALLOCATABLE :: GE_N(:,:)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)
  REAL(DP), ALLOCATABLE :: CF_N(:,:)
  REAL(DP), ALLOCATABLE :: PF_N(:,:)
  REAL(DP), ALLOCATABLE :: AF_N(:,:)
  REAL(DP), ALLOCATABLE :: CR_N(:,:,:,:)

CONTAINS


  ! --- Public Subroutines ---


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(inout) :: &
      dU_F(1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)

    INTEGER :: iN_X

    CALL InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    CALL MapDataForCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    ! --- Compute Primitive Fluid ---

    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_NonRelativistic &
             ( CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), &
               CF_N(iN_X,iCF_S2), &
               CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), &
               CF_N(iN_X,iCF_Ne), &
               PF_N(iN_X,iPF_D ), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), &
               PF_N(iN_X,iPF_Ne), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

    END DO

    ! --- EOS Table Lookup ---

    CALL ComputeThermodynamicStates_Auxiliary_TABLE &
           ( PF_N(:,iPF_D), PF_N(:,iPF_E), PF_N(:,iPF_Ne), &
             AF_N(:,iAF_T), AF_N(:,iAF_E), AF_N(:,iAF_Ye) )

    CALL FinalizeCollisions

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit


  SUBROUTINE Initialize_TwoMoment_Collisions_Neutrinos

  END SUBROUTINE Initialize_TwoMoment_Collisions_Neutrinos


  SUBROUTINE Finalize_TwoMoment_Collisions_Neutrinos

  END SUBROUTINE Finalize_TwoMoment_Collisions_Neutrinos


  ! --- Private Subroutines ---


  SUBROUTINE InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)
    INTEGER, INTENT(in) :: iZ_B1(4), iZ_E1(4)

    iE_B0 = iZ_B0(1)  ; iE_E0 = iZ_E0(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iE_B1 = iZ_B1(1)  ; iE_E1 = iZ_E1(1)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    nZ = iZ_E0 - iZ_B0 + 1
    nX = iX_E0 - iX_B0 + 1
    nE = iE_E0 - iE_E0 + 1

    nX_G = nDOFX * PRODUCT( nX )
    nE_G = nDOFE * nE

    ALLOCATE( GE_N(nE_G,nGE) )
    ALLOCATE( GX_N(nX_G,nGF) )

    ALLOCATE( CF_N(nX_G,nCF) )
    ALLOCATE( PF_N(nX_G,nPF) )
    ALLOCATE( AF_N(nX_G,nAF) )

    ALLOCATE( CR_N(nE_G,nX_G,nSpecies,nCR) )

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

    DEALLOCATE( GE_N, GX_N )
    DEALLOCATE( CF_N, PF_N, AF_N )
    DEALLOCATE( CR_N )

  END SUBROUTINE FinalizeCollisions


  SUBROUTINE MapDataForCollisions &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in)  :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, &
           1:nSpecies)

    INTEGER :: iGE, iGF, iCF, iCR, iS
    INTEGER :: iE, iX1, iX2, iX3
    INTEGER :: iNodeE, iNodeX, iNodeZ
    INTEGER :: iN_E, iN_X

    ! --- Momentum Space Geometry ---

#if   defined(THORNADO_OMP_OL)
#elif defined(THORNADO_OACC)
#elif defined(THORNADO_OMP)
#endif
    DO iGE  = 1, nGE
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nE    ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      GE_N(iN_E,iGE) = GE(iNodeE,iE,iGE)

    END DO
    END DO

    ! --- Position Space Geometry ---

#if   defined(THORNADO_OMP_OL)
#elif defined(THORNADO_OACC)
#elif defined(THORNADO_OMP)
#endif
    DO iGF  = 1, nGF
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      GX_N(iN_X,iGF) = GX(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO

    ! --- Fluid Fields ---

#if   defined(THORNADO_OMP_OL)
#elif defined(THORNADO_OACC)
#elif defined(THORNADO_OMP)
#endif
    DO iCF  = 1, nCF
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      CF_N(iN_X,iCF) = U_F(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO

    ! --- Radiation Fields ---

#if   defined(THORNADO_OMP_OL)
#elif defined(THORNADO_OACC)
#elif defined(THORNADO_OMP)
#endif
    DO iS   = 1, nSpecies
    DO iCR  = 1, nCR
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

      CR_N(iN_E,iN_X,iS,iCR) = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapDataForCollisions


END MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV
