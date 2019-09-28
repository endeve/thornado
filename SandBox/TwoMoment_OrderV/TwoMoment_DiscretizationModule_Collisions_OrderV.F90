MODULE TwoMoment_DiscretizationModule_Collisions_OrderV

  USE KindModule, ONLY: &
    DP, Zero
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Im_MapForward, &
    Timer_Im_MapBackward
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_NonRelativistic_IDEAL, ONLY: &
    Euler_ComputePrimitive_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit

  INTEGER :: iE_B0,    iE_E0
  INTEGER :: iX_B0(3), iX_E0(3)
  INTEGER :: nZ(4), nE, nX(3), nE_G, nX_G

  REAL(DP)              :: PF_N(nPF)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)
  REAL(DP), ALLOCATABLE :: CF_N(:,:)
  REAL(DP), ALLOCATABLE :: CR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: dCR_N(:,:,:,:)

CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      dt
    REAL(DP), INTENT(in)  :: &
      GE  (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)  :: &
      GX  (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                   iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)  :: &
      U_F (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                   iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in)  :: &
      U_R (1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF
    INTEGER :: iX1, iX2, iX3, iE
    INTEGER :: iNodeZ, iNodeX, iNodeE, iN_X, iN_E

    PRINT*, "      ComputeIncrement_TwoMoment_Implicit (In)"

    CALL InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Im_MapForward )

    ! --- Rearrange Geometry Fields ---

    DO iN_X = 1, nX_G
    DO iGF  = 1, nGF

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      GX_N(iGF,iN_X) = GX(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO

    ! --- Rearrange Fluid Fields ---

    DO iN_X = 1, nX_G
    DO iCF  = 1, nCF

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      CF_N(iCF,iN_X) = U_F(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO

    ! --- Rearrange Radiation Fields ---

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

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      CR_N(iCR,iS,iN_E,iN_X) = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Im_MapForward )

    DO iN_X = 1, nX_G

      CALL Euler_ComputePrimitive_NonRelativistic &
             ( CF_N(iCF_D ,iN_X), &
               CF_N(iCF_S1,iN_X), &
               CF_N(iCF_S2,iN_X), &
               CF_N(iCF_S3,iN_X), &
               CF_N(iCF_E ,iN_X), &
               CF_N(iCF_Ne,iN_X), &
               PF_N(iPF_D ), &
               PF_N(iPF_V1), &
               PF_N(iPF_V2), &
               PF_N(iPF_V3), &
               PF_N(iPF_E ), &
               PF_N(iPF_Ne), &
               GX_N(iGF_Gm_dd_11,iN_X), &
               GX_N(iGF_Gm_dd_22,iN_X), &
               GX_N(iGF_Gm_dd_33,iN_X) )

      PRINT*, "V1 = ", PF_N(iPF_V1)
      PRINT*, "V2 = ", PF_N(iPF_V2)
      PRINT*, "V3 = ", PF_N(iPF_V3)

      DO iN_E = 1, nE_G
      DO iS   = 1, nSpecies
      DO iCR  = 1, nCR

        dCR_N(iCR,iS,iN_E,iN_X) = 0.0_DP

      END DO
      END DO
      END DO

    END DO

    CALL TimersStart( Timer_Im_MapBackward )

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0, iE_E0

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

        iN_X = iNodeX &
               + ( iX1 - iX_B0(1) ) * nDOFX &
               + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
               + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)
        iN_E = iNodeE &
             + ( iE  - iE_B0    ) * nDOFE

        dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) = dCR_N(iCR,iS,iN_E,iN_X)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Im_MapBackward )

    CALL FinalizeCollisions

    PRINT*, "      ComputeIncrement_TwoMoment_Implicit (Out)"

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit


  SUBROUTINE InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)

    nZ = iZ_E0 - iZ_B0 + 1
    nE = iE_E0 - iE_B0 + 1
    nX = iX_E0 - iX_B0 + 1

    nE_G = nDOFE * nE
    nX_G = nDOFX * PRODUCT( nX )

    ALLOCATE( GX_N(nGF,nX_G) )
    ALLOCATE( CF_N(nCF,nX_G) )
    ALLOCATE( CR_N (nCR,nSpecies,nE_G,nX_G) )
    ALLOCATE( dCR_N(nCR,nSpecies,nE_G,nX_G) )

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

    DEALLOCATE( GX_N, CF_N, CR_N, dCR_N )

  END SUBROUTINE FinalizeCollisions


END MODULE TwoMoment_DiscretizationModule_Collisions_OrderV
