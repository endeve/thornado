MODULE TwoMoment_DiscretizationModule_Collisions_OrderV

  USE KindModule, ONLY: &
    DP, Zero, One, Half, Three
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
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE TwoMoment_TimersModule_OrderV, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Collisions, &
    Timer_Collisions_PrimitiveFluid, &
    Timer_Collisions_Solve
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeEddingtonTensorComponents_dd
  USE TwoMoment_OpacityModule_OrderV, ONLY: &
    uOP, iOP_D0, iOP_Chi, iOP_Sigma, nOP
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor, &
    HeatFluxFactor
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit

  INTEGER :: iE_B0,    iE_E0
  INTEGER :: iX_B0(3), iX_E0(3)
  INTEGER :: nZ(4), nE, nX(3), nE_G, nX_G

  REAL(DP), ALLOCATABLE :: GX_N(:,:)
  REAL(DP), ALLOCATABLE :: PF_N(:,:)
  REAL(DP), ALLOCATABLE :: CF_N(:,:)
  REAL(DP), ALLOCATABLE :: CR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: dCR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: OP_N(:,:,:,:)

CONTAINS


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
    REAL(DP), INTENT(in) :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR,1:nSpecies)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF, iOP
    INTEGER :: iX1, iX2, iX3, iE
    INTEGER :: iNodeZ, iNodeX, iNodeE, iN_X, iN_E

    CALL TimersStart( Timer_Collisions )

    CALL InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )


#if defined(THORNADO_OMP_OL)
    !!$OMP TARGET ENTER DATA &
    !!$OMP MAP( to: GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !!$OMP MAP( alloc: dU_F, dU_R )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( GX, U_F, U_R, iZ_B0, iZ_E0, iZ_B1, iZ_E1 ) &
    !$ACC CREATE( dU_F, dU_R )
#endif



!PRINT*
!PRINT*, " BEFORE "
!PRINT*, "  uOP = ", uOP(1,1,2,1,2,1,1)
!PRINT*, "  OP_N = ", OP_N(:,1,1,1)
!PRINT*, "  dU_F = ", dU_F(:,2,1,1,1)


#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B1, iZ_E1 ) &
    !$OMP MAP( alloc: dU_F, dU_R )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( iZ_B1, iZ_E1 ) &
    !$ACC CREATE( dU_F, dU_R )
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( dU_F, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)

      DO iNodeX = 1, nDOFX

        dU_F(iNodeX,iZ2,iZ3,iZ4,iCF) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO



#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU_R, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(7)
#endif
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




#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: dU_F, dU_R ) &
    !$OMP MAP( release: iZ_B1, iZ_E1 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( dU_F, dU_R ) &
    !$ACC DELETE( iZ_B1, iZ_E1 )
#endif


#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: GX, U_F, U_R, uOP, iZ_B1, iZ_E1, iX_B0, nX, nZ ) &
    !$OMP MAP( alloc: GX_N, CF_N, CR_N, OP_N)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( GX, U_F, U_R, uOP, iZ_B1, iZ_E1, iX_B0, nX, nZ ) &
    !$ACC CREATE( GX_N, CF_N, CR_N, OP_N )
#endif


    ! --- Arrange Geometry Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, GX_N, GX )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iX1, iX2, iX3, iNodeX )
#endif
    DO iN_X = 1, nX_G
    DO iGF  = 1, nGF

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      GX_N(iGF,iN_X) = GX(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO

    ! --- Arrange Fluid Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, CF_N, U_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iX1, iX2, iX3, iNodeX )
#endif
    DO iN_X = 1, nX_G
    DO iCF  = 1, nCF

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      CF_N(iCF,iN_X) = U_F(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO

    ! --- Arrange Radiation Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iNodeE, iNodeZ, iNodeX, iX1, iX2, iX3, iE )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iNodeZ, iNodeX, iNodeE, iX1, iX2, iX3, iE ) &
    !$ACC PRESENT( nZ, nX, iX_B0, CR_N, U_R )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iE, iX1, iX2, iX3, iNodeE, iNodeX, iNodeZ )
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

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      CR_N(iCR,iS,iN_E,iN_X) = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO

    ! --- Arrange Opacities ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iNodeE, iNodeZ, iNodeX, iX1, iX2, iX3, iE )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iNodeZ, iNodeX, iNodeE, iX1, iX2, iX3, iE ) &
    !$ACC PRESENT( nZ, nX, iX_B0, OP_N, uOP )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iE, iX1, iX2, iX3, iNodeE, iNodeX, iNodeZ )
#endif
    DO iS   = 1, nSpecies
    DO iOP  = 1, nOP
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      OP_N(iOP,iS,iN_E,iN_X) = uOP(iNodeZ,iE,iX1,iX2,iX3,iOP,iS)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: GX_N, CF_N, CR_N, OP_N ) &
    !$OMP MAP( release: GX, U_F, U_R, uOP, iZ_B1, iZ_E1, iX_B0, nX )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( GX_N, CF_N, CR_N, OP_N ) &
    !$ACC DELETE( GX, U_F, U_R, uOP, iZ_B1, iZ_E1, iX_B0, nX )
#endif



!PRINT*
!PRINT*, " AFTER "
!PRINT*, "  uOP = ", uOP(1,1,2,1,2,1,1)
!PRINT*, "  OP_N = ", OP_N(:,1,1,1)
!PRINT*, "  dU_F = ", dU_F(:,2,1,1,1)


#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: CF_N, GX_N ) &
    !$OMP MAP( alloc: PF_N)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( CF_N, GX_N ) &
    !$ACC CREATE( PF_N )
#endif


    CALL TimersStart( Timer_Collisions_PrimitiveFluid )



#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( CF_N, PF_N, GX_N )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_NonRelativistic &
             ( CF_N(iCF_D       ,iN_X), &
               CF_N(iCF_S1      ,iN_X), &
               CF_N(iCF_S2      ,iN_X), &
               CF_N(iCF_S3      ,iN_X), &
               CF_N(iCF_E       ,iN_X), &
               CF_N(iCF_Ne      ,iN_X), &
               PF_N(iPF_D       ,iN_X), &
               PF_N(iPF_V1      ,iN_X), &
               PF_N(iPF_V2      ,iN_X), &
               PF_N(iPF_V3      ,iN_X), &
               PF_N(iPF_E       ,iN_X), &
               PF_N(iPF_Ne      ,iN_X), &
               GX_N(iGF_Gm_dd_11,iN_X), &
               GX_N(iGF_Gm_dd_22,iN_X), &
               GX_N(iGF_Gm_dd_33,iN_X) )

    END DO

    CALL TimersStop( Timer_Collisions_PrimitiveFluid )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: PF_N ) &
    !$OMP MAP( release: CF_N, GX_N )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( PF_N ) &
    !$ACC DELETE( CF_N, GX_N )
#endif




#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: CR_N, PF_N, GX_N, OP_N ) &
    !$OMP MAP( alloc: dCR_N )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( CR_N, PF_N, GX_N, OP_N ) &
    !$ACC CREATE( dCR_N )
#endif


    CALL TimersStart( Timer_Collisions_Solve )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( CR_N, PF_N, GX_N, OP_N, dCR_N )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G
    DO iS   = 1, nSpecies

      CALL ComputeIncrement_FixedPoint_Richardson &
             ( dt, &
               CR_N (iCR_N       ,iS,iN_E,iN_X), &
               CR_N (iCR_G1      ,iS,iN_E,iN_X), &
               CR_N (iCR_G2      ,iS,iN_E,iN_X), &
               CR_N (iCR_G3      ,iS,iN_E,iN_X), &
               PF_N (iPF_V1              ,iN_X), &
               PF_N (iPF_V2              ,iN_X), &
               PF_N (iPF_V3              ,iN_X), &
               GX_N (iGF_Gm_dd_11        ,iN_X), &
               GX_N (iGF_Gm_dd_22        ,iN_X), &
               GX_N (iGF_Gm_dd_33        ,iN_X), &
               OP_N (iOP_D0      ,iS,iN_E,iN_X), &
               OP_N (iOP_Chi     ,iS,iN_E,iN_X), &
               OP_N (iOP_Sigma   ,iS,iN_E,iN_X), &
               dCR_N(iCR_N       ,iS,iN_E,iN_X), &
               dCR_N(iCR_G1      ,iS,iN_E,iN_X), &
               dCR_N(iCR_G2      ,iS,iN_E,iN_X), &
               dCR_N(iCR_G3      ,iS,iN_E,iN_X) )

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Collisions_Solve )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: dCR_N ) &
    !$OMP MAP( release: CR_N, PF_N, GX_N, OP_N )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( dCR_N ) &
    !$ACC DELETE( CR_N, PF_N, GX_N, OP_N )
#endif


#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dCR_N, dU_R, nX, iX_B0, iX_E0 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( dCR_N, dU_R, nX, iX_B0, iX_E0)
#endif


    ! --- Revert Radiation Increment ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeE, iNodeX, iN_X, iN_E )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRIVATE( iNodeX, iNodeE, iN_X, iN_E ) &
    !$ACC PRESENT( nX, iX_B0, iX_E0, dU_R, dCR_N )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeE, iNodeX, iN_E, iN_X )
#endif
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
                 + (iX1-iX_B0(1)) * nDOFX &
                 + (iX2-iX_B0(2)) * nDOFX * nX(1) &
                 + (iX3-iX_B0(3)) * nDOFX * nX(1) * nX(2)

        iN_E = iNodeE &
                 + (iE-iE_B0) * nDOFE

        dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) = dCR_N(iCR,iS,iN_E,iN_X)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: dU_R ) &
    !$OMP MAP( release: dCR_N, nX, iX_B0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( dU_R ) &
    !$ACC DELETE( dCR_N, nX, iX_B0 )
#endif



    CALL FinalizeCollisions

    CALL TimersStop( Timer_Collisions )

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit


  SUBROUTINE ComputeIncrement_FixedPoint &
    ( dt, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, D_0, Chi, Sigma, &
      dN, dG_d_1, dG_d_2, dG_d_3 )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: dt
    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: D_0, Chi, Sigma
    REAL(DP), INTENT(out) :: dN, dG_d_1, dG_d_2, dG_d_3

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: LWORK = 2 * M
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, j, k, mk, INFO
    REAL(DP) :: D, I_d_1, I_d_2, I_d_3, Kappa
    REAL(DP) ::    I_u_1, I_u_2, I_u_3
    REAL(DP) :: A_d_1, A_d_2, A_d_3, k_dd(3,3)
    REAL(DP) :: D_00, D_ii
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: LMAT(4,4), DET, Alpha(M)
    REAL(DP) :: BVEC(4), AMAT(4,M), WORK(LWORK)

    Kappa = Chi + Sigma

    D_00 = One + dt * Chi
    D_ii = One + dt * Kappa

    ! --- Constant Vector ---

    CVEC = [ N + dt * Chi * D_0, G_d_1, G_d_2, G_d_3 ]

    ! --- Initial Guess ---

    D     = N
    I_u_1 = Zero
    I_u_2 = Zero
    I_u_3 = Zero

    I_d_1 = Gm_dd_11 * I_u_1
    I_d_2 = Gm_dd_22 * I_u_2
    I_d_3 = Gm_dd_33 * I_u_3

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIterations )

      k  = k + 1
      mk = MIN( M, k )

      UVEC = [ D, I_d_1, I_d_2, I_d_3 ]

      k_dd = EddingtonTensorComponents_dd &
               ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      A_d_1 = V_u_1 * k_dd(1,1) + V_u_2 * k_dd(1,2) + V_u_3 * k_dd(1,3)
      A_d_2 = V_u_1 * k_dd(1,2) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(2,3)
      A_d_3 = V_u_1 * k_dd(1,3) + V_u_2 * k_dd(2,3) + V_u_3 * k_dd(3,3)

      DET = ( D_00 * D_ii &
              - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 + V_u_3 * A_d_3 ) ) * D_ii**2

      LMAT(1,1) = D_ii**3
      LMAT(2,1) = - A_d_1 * D_ii**2
      LMAT(3,1) = - A_d_2 * D_ii**2
      LMAT(4,1) = - A_d_3 * D_ii**2

      LMAT(1,2) = - V_u_1 * D_ii**2
      LMAT(2,2) = D_00 * D_ii**2 &
                    - ( V_u_2 * A_d_2 + V_u_3 * A_d_3 ) * D_ii
      LMAT(3,2) = V_u_1 * A_d_2 * D_ii
      LMAT(4,2) = V_u_1 * A_d_3 * D_ii

      LMAT(1,3) = - V_u_2 * D_ii**2
      LMAT(2,3) = V_u_2 * A_d_1 * D_ii
      LMAT(3,3) = D_00 * D_ii**2 &
                    - ( V_u_1 * A_d_1 + V_u_3 * A_d_3 ) * D_ii
      LMAT(4,3) = V_u_2 * A_d_3 * D_ii

      LMAT(1,4) = - V_u_3 * D_ii**2
      LMAT(2,4) = V_u_3 * A_d_1 * D_ii
      LMAT(3,4) = V_u_3 * A_d_2 * D_ii
      LMAT(4,4) = D_00 * D_ii**2 &
                    - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 ) * D_ii

      LMAT = LMAT / DET

      GVEC(:,mk) = Zero

      DO j = 1, 4
      DO i = 1, 4

        GVEC(i,mk) = GVEC(i,mk) + LMAT(i,j) * CVEC(j)

      END DO
      END DO

      FVEC(:,mk) = GVEC(:,mk) - UVEC

      IF( mk == 1 )THEN

        ! --- Picard Iteration ---

        GVECm = GVEC(:,mk)

      ELSE

        ! --- Anderson Accelerated Fixed-Point ---

        Alpha = Alpha_LS( M, mk, FVEC )

        GVECm = Zero
        DO i = 1, mk

          GVECm = GVECm + Alpha(i) * GVEC(:,i)

        END DO

      END IF

      FVECm = GVECm - UVEC

      IF( ALL( ABS( FVECm ) <= Rtol * ABS( CVEC ) ) )THEN

        CONVERGED = .TRUE.

      END IF

      UVEC = GVECm

      IF( mk == M .AND. .NOT. CONVERGED )THEN

        FVEC = ShiftVec( M, mk, FVEC )
        GVEC = ShiftVec( M, mk, GVEC )

      END IF

      D     = UVEC(1)
      I_d_1 = UVEC(2); I_u_1 = I_d_1 / Gm_dd_11
      I_d_2 = UVEC(3); I_u_2 = I_d_2 / Gm_dd_22
      I_d_3 = UVEC(4); I_u_3 = I_d_3 / Gm_dd_33

    END DO

    dN     = Chi * ( D_0 - D )
    dG_d_1 = - Kappa * I_d_1
    dG_d_2 = - Kappa * I_d_2
    dG_d_3 = - Kappa * I_d_3

    ! IF( k == MaxIterations )THEN

    !   PRINT*
    !   PRINT*, "ComputeIncrement_FixedPoint"
    !   PRINT*
    !   PRINT*, "  N     = ", N
    !   PRINT*, "  G_d_1 = ", G_d_1
    !   PRINT*, "  G_d_2 = ", G_d_2
    !   PRINT*, "  G_d_3 = ", G_d_3
    !   PRINT*
    !   PRINT*, "  V_u_1 = ", V_u_1
    !   PRINT*, "  V_u_2 = ", V_u_2
    !   PRINT*, "  V_u_3 = ", V_u_3
    !   PRINT*

    !   PRINT*, "  Converged with k = ", k

    !   PRINT*
    !   PRINT*, "  FVECm = ", FVECm
    !   PRINT*

    !   PRINT*
    !   PRINT*, "  D     = ", D
    !   PRINT*, "  I_u_1 = ", I_u_1
    !   PRINT*, "  I_u_2 = ", I_u_2
    !   PRINT*, "  I_u_3 = ", I_u_3
    !   PRINT*

    ! END IF

  END SUBROUTINE ComputeIncrement_FixedPoint


  SUBROUTINE ComputeIncrement_FixedPoint_Richardson &
    ( dt, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, D_0, Chi, Sigma, &
      dN, dG_d_1, dG_d_2, dG_d_3 )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: dt
    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: D_0, Chi, Sigma
    REAL(DP), INTENT(out) :: dN, dG_d_1, dG_d_2, dG_d_3

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 1
    INTEGER,  PARAMETER :: LWORK = 2 * M
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, j, k, mk, INFO
    REAL(DP) :: D, I_d_1, I_d_2, I_d_3, Kappa
    REAL(DP) ::    I_u_1, I_u_2, I_u_3
    REAL(DP) :: k_dd(3,3)
    ! REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: D_00, D_ii
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: vMag, Omega, vI, vK, Alpha(M)
    REAL(DP) :: BVEC(4), AMAT(4,M), WORK(LWORK)

    Kappa = Chi + Sigma

    D_00 = One + dt * Chi
    D_ii = One + dt * Kappa

    ! --- Constant Vector ---

    CVEC = [ N, G_d_1, G_d_2, G_d_3 ]

    ! --- Initial Guess ---

    D     = N
    I_u_1 = Zero
    I_u_2 = Zero
    I_u_3 = Zero

    I_d_1 = Gm_dd_11 * I_u_1
    I_d_2 = Gm_dd_22 * I_u_2
    I_d_3 = Gm_dd_33 * I_u_3

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIterations )

      k  = k + 1
      mk = MIN( M, k )

      UVEC = [ D, I_d_1, I_d_2, I_d_3 ]

      ! CALL ComputeEddingtonTensorComponents_dd &
      !        ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      !          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )

      k_dd = EddingtonTensorComponents_dd &
               ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )


      vMag = SQRT(    V_u_1 * Gm_dd_11 * V_u_1 &
                    + V_u_2 * Gm_dd_22 * V_u_2 &
                    + V_u_3 * Gm_dd_33 * V_u_3 )

      Omega = One / ( One + vMag )

      vI = V_u_1 * UVEC(2) + V_u_2 *  UVEC(3) + V_u_3 *  UVEC(4)
      GVEC(1,mk) = (One - Omega) * UVEC(1) + &
                   Omega / D_00 * (CVEC(1) + dt * Chi * D_0 - vI)

      DO j = 1, 3
        vK = V_u_1 * k_dd(j,1) + V_u_2 * k_dd(j,2) + V_u_3 * k_dd(j,3)
        GVEC(j+1,mk) = (One - Omega) * UVEC(j+1) + &
                     Omega / D_ii * (CVEC(j+1) - vK * UVEC(1))
      END DO


      FVEC(:,mk) = GVEC(:,mk) - UVEC

      IF( mk == 1 )THEN

        ! --- Picard Iteration ---

        GVECm = GVEC(:,mk)

      ELSE

        Alpha = Alpha_LS( M, mk, FVEC )

        GVECm = Zero
        DO i = 1, mk

          GVECm = GVECm + Alpha(i) * GVEC(:,i)

        END DO

      END IF

      FVECm = GVECm - UVEC

      IF( ALL( ABS( FVECm ) <= Rtol * ABS( CVEC ) ) )THEN

        CONVERGED = .TRUE.

      END IF

      UVEC = GVECm

      IF( mk == M .AND. .NOT. CONVERGED )THEN

        FVEC = ShiftVec( M, mk, FVEC )
        GVEC = ShiftVec( M, mk, GVEC )

      END IF

      D     = UVEC(1)
      I_d_1 = UVEC(2); I_u_1 = I_d_1 / Gm_dd_11
      I_d_2 = UVEC(3); I_u_2 = I_d_2 / Gm_dd_22
      I_d_3 = UVEC(4); I_u_3 = I_d_3 / Gm_dd_33

    END DO

    dN     = Chi * ( D_0 - D )
    dG_d_1 = - Kappa * I_d_1
    dG_d_2 = - Kappa * I_d_2
    dG_d_3 = - Kappa * I_d_3

    ! IF( k == MaxIterations )THEN

    !   PRINT*
    !   PRINT*, "ComputeIncrement_FixedPoint"
    !   PRINT*
    !   PRINT*, "  N     = ", N
    !   PRINT*, "  G_d_1 = ", G_d_1
    !   PRINT*, "  G_d_2 = ", G_d_2
    !   PRINT*, "  G_d_3 = ", G_d_3
    !   PRINT*
    !   PRINT*, "  V_u_1 = ", V_u_1
    !   PRINT*, "  V_u_2 = ", V_u_2
    !   PRINT*, "  V_u_3 = ", V_u_3
    !   PRINT*

    !   PRINT*, "  Converged with k = ", k

    !   PRINT*
    !   PRINT*, "  FVECm = ", FVECm
    !   PRINT*

    !   PRINT*
    !   PRINT*, "  D     = ", D
    !   PRINT*, "  I_u_1 = ", I_u_1
    !   PRINT*, "  I_u_2 = ", I_u_2
    !   PRINT*, "  I_u_3 = ", I_u_3
    !   PRINT*

    ! END IF

  END SUBROUTINE ComputeIncrement_FixedPoint_Richardson


  FUNCTION Alpha_LS( M, mk, FVEC )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: M, mk
    REAL(DP), INTENT(in) :: FVEC(4,M)
    REAL(DP)             :: Alpha_LS(M)

    INTEGER  :: i
    REAL(DP) :: BVEC(4), AMAT(4,M)
    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA, SUM1

    BVEC = - FVEC(:,mk)

    DO i = 1, mk - 1

      AMAT(:,i) = FVEC(:,i) - FVEC(:,mk)

    END DO

    IF( mk == 2 )THEN

      AA11 = Zero
      AB1  = Zero

      DO i = 1, 4

        AA11 = AA11 + AMAT(i,1) * AMAT(i,1)
        AB1  = AB1  + AMAT(i,1) * BVEC(i)

      END DO

      BVEC(1) = AB1 / AA11

    ELSEIF( mk == 3 )THEN

      AA11 = Zero
      AA12 = Zero
      AA22 = Zero
      AB1  = Zero
      AB2  = Zero

      DO i = 1, 4

        AA11 = AA11 + AMAT(i,1) * AMAT(i,1)
        AA12 = AA12 + AMAT(i,1) * AMAT(i,2)
        AA22 = AA22 + AMAT(i,2) * AMAT(i,2)
        AB1  = AB1  + AMAT(i,1) * BVEC(i)
        AB2  = AB2  + AMAT(i,2) * BVEC(i)

      END DO

      DET_AA = AA11 * AA22 - AA12 * AA12

      BVEC(1) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
      BVEC(2) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

    ELSEIF( mk > 3 )THEN

      ! --- Not Implemented ---

    END IF

    SUM1 = Zero
    DO i = 1, mk - 1

      Alpha_LS(i) = BVEC(i)

      SUM1 = SUM1 + BVEC(i)

    END DO

    Alpha_LS(mk) = One - SUM1

    RETURN
  END FUNCTION Alpha_LS


  FUNCTION ShiftVec( M, mk, Vec )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: M, mk
    REAL(DP), INTENT(in) :: Vec(4,M)
    REAL(DP)             :: ShiftVec(4,M)

    INTEGER  :: i, j
    REAL(DP) :: VecTMP(4,M)

    DO j = 1, mk - 1
    DO i = 1, 4

      VecTMP(i,j) = Vec(i,j+1)

    END DO
    END DO

    DO j = 1, mk - 1
    DO i = 1, 4

      ShiftVec(i,j) = VecTMP(i,j)

    END DO
    END DO

    RETURN
  END FUNCTION ShiftVec




  SUBROUTINE SolveAlpha_LS( M, mk, FVEC, Alpha )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)    :: M, mk
    REAL(DP), INTENT(inout) :: FVEC(4,M), Alpha(M)

    INTEGER  :: i
    REAL(DP) :: BVEC(4), AMAT(4,M)
    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA, SUM1

    BVEC = - FVEC(:,mk)

    DO i = 1, mk - 1

      AMAT(:,i) = FVEC(:,i) - FVEC(:,mk)

    END DO

    IF( mk == 2 )THEN

      AA11 = Zero
      AB1  = Zero

      DO i = 1, 4

        AA11 = AA11 + AMAT(i,1) * AMAT(i,1)
        AB1  = AB1  + AMAT(i,1) * BVEC(i)

      END DO

      BVEC(1) = AB1 / AA11

    ELSEIF( mk == 3 )THEN

      AA11 = Zero
      AA12 = Zero
      AA22 = Zero
      AB1  = Zero
      AB2  = Zero

      DO i = 1, 4

        AA11 = AA11 + AMAT(i,1) * AMAT(i,1)
        AA12 = AA12 + AMAT(i,1) * AMAT(i,2)
        AA22 = AA22 + AMAT(i,2) * AMAT(i,2)
        AB1  = AB1  + AMAT(i,1) * BVEC(i)
        AB2  = AB2  + AMAT(i,2) * BVEC(i)

      END DO

      DET_AA = AA11 * AA22 - AA12 * AA12

      BVEC(1) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
      BVEC(2) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

    ELSEIF( mk > 3 )THEN

!      STOP

    END IF

    SUM1 = Zero
    DO i = 1, mk - 1

      Alpha(i) = BVEC(i)

      SUM1 = SUM1 + BVEC(i)

    END DO

    Alpha(mk) = One - SUM1

  END SUBROUTINE SolveAlpha_LS


  SUBROUTINE ShiftVectors( M, mk, FVEC, GVEC )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)    :: M, mk
    REAL(DP), INTENT(inout) :: FVEC(4,M), GVEC(4,M)

    INTEGER  :: i, j
    REAL(DP) :: FTMP(4,M), GTMP(4,M)

    DO j = 1, mk - 1
    DO i = 1, 4

      FTMP(i,j) = FVEC(i,j+1)
      GTMP(i,j) = GVEC(i,j+1)

    END DO
    END DO

    DO j = 1, mk - 1
    DO i = 1, 4

      FVEC(i,j) = FTMP(i,j)
      GVEC(i,j) = GTMP(i,j)

    END DO
    END DO

  END SUBROUTINE ShiftVectors

  FUNCTION EddingtonTensorComponents_dd &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP)              :: &
      EddingtonTensorComponents_dd(3,3)

    INTEGER  :: i, j
    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_d(3), Gm_dd(3,3)

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    h_d(1) = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d(2) = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d(3) = Gm_dd_33 * I_u_3 / ( FF * D )

    Gm_dd = Zero
    Gm_dd(1,1) = Gm_dd_11
    Gm_dd(2,2) = Gm_dd_22
    Gm_dd(3,3) = Gm_dd_33

    DO j = 1, 3
    DO i = 1, 3

      EddingtonTensorComponents_dd(i,j) &
        = a * Gm_dd(i,j) + b * h_d(i) * h_d(j)

    END DO
    END DO

    RETURN
  END FUNCTION EddingtonTensorComponents_dd

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
    ALLOCATE( PF_N(nPF,nX_G) )
    ALLOCATE( CF_N(nCF,nX_G) )
    ALLOCATE( CR_N (nCR,nSpecies,nE_G,nX_G) )
    ALLOCATE( dCR_N(nCR,nSpecies,nE_G,nX_G) )
    ALLOCATE( OP_N (nOP,nSpecies,nE_G,nX_G) )

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

    DEALLOCATE( GX_N, PF_N, CF_N, CR_N, dCR_N, OP_N )

  END SUBROUTINE FinalizeCollisions


END MODULE TwoMoment_DiscretizationModule_Collisions_OrderV
