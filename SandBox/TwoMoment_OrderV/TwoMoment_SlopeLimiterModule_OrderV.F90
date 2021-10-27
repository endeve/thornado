MODULE TwoMoment_SlopeLimiterModule_OrderV

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFE, nDOFX, nDimsX
  USE TwoMoment_TimersModule_OrderV, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_SL, &
    Timer_SL_Permute, &
    Timer_SL_LinearAlgebra, &
    Timer_SL_MinMod, &
    Timer_SL_ReplaceSlopes, &
    Timer_SL_Correction
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply
  USE UtilitiesModule, ONLY: &
    MinModB
  USE ReferenceElementModuleX, ONLY: &
    NodesX_q, WeightsX_q
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModalX, &
    MapModalToNodalX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
    DetectTroubledCells_TwoMoment

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_TwoMoment
  PUBLIC :: FinalizeSlopeLimiter_TwoMoment
  PUBLIC :: ApplySlopeLimiter_TwoMoment

  CHARACTER(4) :: SlopeLimiterMethod
  LOGICAL      :: UseSlopeLimiter
  INTEGER      :: nE, nE_G
  REAL(DP)     :: BetaTVD, BetaTVB
  REAL(DP)     :: SlopeTolerance
  REAL(DP), ALLOCATABLE :: N2M_Vec_0(:)
  REAL(DP), ALLOCATABLE :: N2M_Vec_1(:)
  REAL(DP), ALLOCATABLE :: N2M_Vec_2(:)
  REAL(DP), ALLOCATABLE :: N2M_Vec_3(:)
  REAL(DP), ALLOCATABLE :: M2N_Vec_0(:)
  REAL(DP), ALLOCATABLE :: M2N_Vec_1(:)
  REAL(DP), ALLOCATABLE :: M2N_Vec_2(:)
  REAL(DP), ALLOCATABLE :: M2N_Vec_3(:)

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE &
  !$OMP TARGET( N2M_Vec_0, N2M_Vec_1, N2M_Vec_2, N2M_Vec_3, &
  !$OMP         M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3, &
  !$OMP         BetaTVD, BetaTVB, SlopeTolerance )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE &
  !$ACC CREATE( N2M_Vec_0, N2M_Vec_1, N2M_Vec_2, N2M_Vec_3, &
  !$ACC         M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3, &
  !$ACC         BetaTVD, BetaTVB, SlopeTolerance )
#endif

CONTAINS


  SUBROUTINE InitializeSlopeLimiter_TwoMoment &
    ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
      UseSlopeLimiter_Option, SlopeLimiterMethod_Option, Verbose_Option )

    REAL(DP),     INTENT(in), OPTIONAL :: BetaTVD_Option
    REAL(DP),     INTENT(in), OPTIONAL :: BetaTVB_Option
    REAL(DP),     INTENT(in), OPTIONAL :: SlopeTolerance_Option
    LOGICAL,      INTENT(in), OPTIONAL :: UseSlopeLimiter_Option
    CHARACTER(*), INTENT(in), OPTIONAL :: SlopeLimiterMethod_Option
    LOGICAL,      INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose
    INTEGER :: iNodeX

    IF( PRESENT( BetaTVD_Option ) )THEN
      BetaTVD = BetaTVD_Option
    ELSE
      BetaTVD = One
    END IF

    IF( PRESENT( BetaTVB_Option ) )THEN
      BetaTVB = BetaTVB_Option
    ELSE
      BetaTVB = Zero
    END IF

    IF( PRESENT( SlopeTolerance_Option ) )THEN
      SlopeTolerance = SlopeTolerance_Option
    ELSE
      SlopeTolerance = 1.0d-6
    END IF

    IF( PRESENT( UseSlopeLimiter_Option ) )THEN
      UseSlopeLimiter = UseSlopeLimiter_Option
    ELSE
      UseSlopeLimiter = .TRUE.
    END IF

    IF( PRESENT( SlopeLimiterMethod_Option ) )THEN
      SlopeLimiterMethod = TRIM( SlopeLimiterMethod_Option )
    ELSE
      SlopeLimiterMethod = 'TVD'
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
        '  INFO: InitializeSlopeLimiter_TwoMoment:'
      WRITE(*,'(A)') &
        '  ---------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A32,L1)'       ) '', 'Use Slope Limiter: ' , &
        UseSlopeLimiter
      WRITE(*,*)
      WRITE(*,'(A4,A32,A)')         '', 'Method: ', &
        TRIM( SlopeLimiterMethod )
      WRITE(*,'(A4,A32,ES11.3E3)' ) '', 'BetaTVD: ' , &
        BetaTVD
      WRITE(*,'(A4,A32,ES11.3E3)' ) '', 'BetaTVB: ' , &
        BetaTVB
      WRITE(*,'(A4,A32,ES11.3E3)' ) '', 'SlopeTolerance: ' , &
        SlopeTolerance

    END IF

    ! --- For Computing Modal Coefficients from Nodal Values ---

    ALLOCATE( N2M_Vec_0(nDOFX) )
    ALLOCATE( N2M_Vec_1(nDOFX) )
    ALLOCATE( N2M_Vec_2(nDOFX) )
    ALLOCATE( N2M_Vec_3(nDOFX) )

    DO iNodeX = 1, nDOFX

      N2M_Vec_0(iNodeX) =           WeightsX_q(iNodeX)
      N2M_Vec_1(iNodeX) = 12.0_DP * WeightsX_q(iNodeX) * NodesX_q(1,iNodeX)
      N2M_Vec_2(iNodeX) = 12.0_DP * WeightsX_q(iNodeX) * NodesX_q(2,iNodeX)
      N2M_Vec_3(iNodeX) = 12.0_DP * WeightsX_q(iNodeX) * NodesX_q(3,iNodeX)

    END DO

    ! --- For Computing Nodal Values from Modal Coefficients ---

    ALLOCATE( M2N_Vec_0(nDOFX) )
    ALLOCATE( M2N_Vec_1(nDOFX) )
    ALLOCATE( M2N_Vec_2(nDOFX) )
    ALLOCATE( M2N_Vec_3(nDOFX) )

    DO iNodeX = 1, nDOFX

      M2N_Vec_0(iNodeX) = One
      M2N_Vec_1(iNodeX) = NodesX_q(1,iNodeX)
      M2N_Vec_2(iNodeX) = NodesX_q(2,iNodeX)
      M2N_Vec_3(iNodeX) = NodesX_q(3,iNodeX)

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE &
    !$OMP TO( N2M_Vec_0, N2M_Vec_1, N2M_Vec_2, N2M_Vec_3, &
    !$OMP     M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3, &
    !$OMP     BetaTVD, BetaTVB, SlopeTolerance )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE &
    !$ACC DEVICE( N2M_Vec_0, N2M_Vec_1, N2M_Vec_2, N2M_Vec_3, &
    !$ACC         M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3, &
    !$ACC         BetaTVD, BetaTVB, SlopeTolerance )
#endif

  END SUBROUTINE InitializeSlopeLimiter_TwoMoment


  SUBROUTINE FinalizeSlopeLimiter_TwoMoment

    DEALLOCATE( N2M_Vec_0, N2M_Vec_1, N2M_Vec_2, N2M_Vec_3 )
    DEALLOCATE( M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3 )

  END SUBROUTINE FinalizeSlopeLimiter_TwoMoment


  SUBROUTINE ApplySlopeLimiter_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, SuppressBC_Option )

    INTEGER, INTENT(in)            :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)           :: &
      GE (1:nDOFE, &
          iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)           :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)           :: &
      U_F(1:nDOFX, &
          iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout)        :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    LOGICAL :: SuppressBC

    IF( .NOT. UseSlopeLimiter .OR. nDOFX == 1 ) RETURN

    IF( PRESENT( SuppressBC_Option ) )THEN
      SuppressBC = SuppressBC_Option
    ELSE
      SuppressBC = .FALSE.
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B0, iZ_E0, GE, GX, U_F, U_R )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC COPYIN( iZ_B0, iZ_E0, GE, GX, U_F, U_R )
#endif

    IF( .NOT. SuppressBC )THEN

      CALL ApplyBoundaryConditions_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R )

    END IF

    SELECT CASE ( TRIM( SlopeLimiterMethod ) )

      CASE( 'TVD' )

        CALL ApplySlopeLimiter_TVD &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

      CASE( 'WENO' )

        CALL ApplySlopeLimiter_WENO &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

      CASE DEFAULT

        CALL ApplySlopeLimiter_TVD &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    END SELECT

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: iZ_B0, iZ_E0, GE, GX, U_F, U_R )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC DELETE( iZ_B0, iZ_E0, GE, GX, U_F, U_R )

    !$ACC WAIT
#endif

  END SUBROUTINE ApplySlopeLimiter_TwoMoment


  SUBROUTINE ApplySlopeLimiter_TVD &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE (1:nDOFE, &
          iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX, &
          iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: &
      iZ1, iZ2, iZ3, iZ4, iCR, iS, iE_G, &
      iNodeZ, iNodeE, iNodeX, nV_KX
    REAL(DP) :: &
      dSlope, Alpha
    LOGICAL  :: &
      TroubledCell &
           (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    LOGICAL  :: &
      Limited &
           (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nCR,1:nSpecies)
    REAL(DP) :: &
      C_0  (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nCR,1:nSpecies)
    REAL(DP) :: &
      C_X1 (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nCR,1:nSpecies)
    REAL(DP) :: &
      C_X2 (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nCR,1:nSpecies)
    REAL(DP) :: &
      C_X3 (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nCR,1:nSpecies)
    REAL(DP) :: &
      CL_X1(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nCR,1:nSpecies)
    REAL(DP) :: &
      CL_X2(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nCR,1:nSpecies)
    REAL(DP) :: &
      CL_X3(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nCR,1:nSpecies)
    REAL(DP) :: &
      uCR_K(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nCR,1:nSpecies)
    REAL(DP) :: &
      wSqrtGm &
           (1:nDOFX, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      uCR  (1:nDOFX, &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nCR,1:nSpecies)

    nE   = iZ_E0(1) - iZ_B0(1) + 1
    nE_G = nE * nDOFE

    nV_KX = PRODUCT( SHAPE( C_0 ) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B0, iZ_E0, GE, GX, U_F, U_R ) &
    !$OMP MAP( alloc: TroubledCell, Limited, &
    !$OMP             CL_X1, CL_X2, CL_X3, &
    !$OMP             C_0, C_X1, C_X2, C_X3, &
    !$OMP             uCR_K, wSqrtGm, uCR )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC COPYIN( iZ_B0, iZ_E0, GE, GX, U_F, U_R ) &
    !$ACC CREATE( TroubledCell, Limited, &
    !$ACC         CL_X1, CL_X2, CL_X3, &
    !$ACC         C_0, C_X1, C_X2, C_X3, &
    !$ACC         uCR_K, wSqrtGm, uCR )
#endif

    CALL DetectTroubledCells_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, TroubledCell )

#if   defined( THORNADO_OMP_OL )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE ASYNC DEVICE( TroubledCell )
#endif

    CALL TimersStart( Timer_SL )

    CALL ComputeLimitedSlopes_X1( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, CL_X1 )

    CALL ComputeLimitedSlopes_X2( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, CL_X2 )

    CALL ComputeLimitedSlopes_X3( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, CL_X3 )

    ! --- Permute Radiation Fields ---

    CALL TimersStart( Timer_SL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iZ1, iNodeE, iNodeZ )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
    !$ACC PRIVATE( iZ1, iNodeE, iNodeZ ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, U_R, uCR )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iZ1, iNodeE, iNodeZ )
#endif
    DO iS     = 1       , nSpecies
    DO iCR    = 1       , nCR
    DO iE_G   = 1       , nE_G
    DO iZ4    = iZ_B0(4), iZ_E0(4)
    DO iZ3    = iZ_B0(3), iZ_E0(3)
    DO iZ2    = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1       , nDOFX

      iZ1    = MOD( (iE_G-1) / nDOFE, nE    ) + iZ_B0(1)
      iNodeE = MOD( (iE_G-1)        , nDOFE ) + 1

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      uCR(iNodeX,iZ2,iZ3,iZ4,iE_G,iCR,iS) &
        = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, GX, WeightsX_q, wSqrtGm )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ4    = iZ_B0(4), iZ_E0(4)
    DO iZ3    = iZ_B0(3), iZ_E0(3)
    DO iZ2    = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1       , nDOFX

      wSqrtGm(iNodeX,iZ2,iZ3,iZ4) &
        = WeightsX_q(iNodeX) * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_SL_Permute )

    ! --- Compute Cell Average ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( Alpha )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC PRIVATE( Alpha ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, uCR_K, wSqrtGm, uCR )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( Alpha )
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      Alpha = Zero
      DO iNodeX = 1, nDOFX
        Alpha = Alpha + wSqrtGm(iNodeX,iZ2,iZ3,iZ4) * uCR(iNodeX,iZ2,iZ3,iZ4,iE_G,iCR,iS)
      END DO
      uCR_K(iZ2,iZ3,iZ4,iE_G,iCR,iS) = Alpha

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Compute Legendre Coefficients ---

    CALL TimersStart( Timer_SL_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_0, 1, Zero, C_0 , 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_1, 1, Zero, C_X1, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_2, 1, Zero, C_X2, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_3, 1, Zero, C_X3, 1 )

    CALL TimersStop( Timer_SL_LinearAlgebra )

    CALL TimersStart( Timer_SL_ReplaceSlopes )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( dSlope )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC PRIVATE( dSlope ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, TroubledCell, Limited, &
    !$ACC          CL_X1, CL_X2, CL_X3, C_0, C_X1, C_X2, C_X3, &
    !$ACC          uCR, M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( dSlope )
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      Limited(iZ2,iZ3,iZ4,iE_G,iCR,iS) = .FALSE.

      IF( TroubledCell(iZ2,iZ3,iZ4,iE_G,iS) )THEN

        dSlope &
          = MAX( ABS( CL_X1(iZ2,iZ3,iZ4,iE_G,iCR,iS) &
                      -C_X1(iZ2,iZ3,iZ4,iE_G,iCR,iS) ), &
                 ABS( CL_X2(iZ2,iZ3,iZ4,iE_G,iCR,iS) &
                      -C_X2(iZ2,iZ3,iZ4,iE_G,iCR,iS) ), &
                 ABS( CL_X3(iZ2,iZ3,iZ4,iE_G,iCR,iS) &
                      -C_X3(iZ2,iZ3,iZ4,iE_G,iCR,iS) ) )

        IF( dSlope > SlopeTolerance * ABS( C_0(iZ2,iZ3,iZ4,iE_G,iCR,iS) ) )THEN

          DO iNodeX = 1, nDOFX
            uCR(iNodeX,iZ2,iZ3,iZ4,iE_G,iCR,iS) &
              =   M2N_Vec_0(iNodeX) * C_0  (iZ2,iZ3,iZ4,iE_G,iCR,iS) &
                + M2N_Vec_1(iNodeX) * CL_X1(iZ2,iZ3,iZ4,iE_G,iCR,iS) &
                + M2N_Vec_2(iNodeX) * CL_X2(iZ2,iZ3,iZ4,iE_G,iCR,iS) &
                + M2N_Vec_3(iNodeX) * CL_X3(iZ2,iZ3,iZ4,iE_G,iCR,iS)
          END DO

          Limited(iZ2,iZ3,iZ4,iE_G,iCR,iS) = .TRUE.

        END IF

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_SL_ReplaceSlopes )

    ! --- Conservative Correction ---

    CALL TimersStart( Timer_SL_Correction )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( Alpha )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC PRIVATE( Alpha ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, Limited, uCR_K, wSqrtGm, uCR )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( Alpha )
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      IF( Limited(iZ2,iZ3,iZ4,iE_G,iCR,iS) )THEN

        Alpha = Zero
        DO iNodeX = 1, nDOFX
          Alpha = Alpha + wSqrtGm(iNodeX,iZ2,iZ3,iZ4) * uCR(iNodeX,iZ2,iZ3,iZ4,iE_G,iCR,iS)
        END DO

        IF( ABS( Alpha ) > Zero )THEN
          DO iNodeX = 1, nDOFX
            uCR(iNodeX,iZ2,iZ3,iZ4,iE_G,iCR,iS) &
              = uCR(iNodeX,iZ2,iZ3,iZ4,iE_G,iCR,iS) * uCR_K(iZ2,iZ3,iZ4,iE_G,iCR,iS) / Alpha
          END DO
        END IF

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_SL_Correction )

    CALL TimersStart( Timer_SL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, iNodeE, iE_G )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
    !$ACC PRIVATE( iNodeX, iNodeE, iE_G ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, U_R, uCR )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeX, iNodeE, iE_G )
#endif
    DO iS     = 1, nSpecies
    DO iCR    = 1, nCR
    DO iZ4    = iZ_B0(4), iZ_E0(4)
    DO iZ3    = iZ_B0(3), iZ_E0(3)
    DO iZ2    = iZ_B0(2), iZ_E0(2)
    DO iZ1    = iZ_B0(1), iZ_E0(1)
    DO iNodeZ = 1, nDOFZ

      iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
      iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

      iE_G   = iNodeE + ( iZ1 - iZ_B0(1) ) * nDOFE

      U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
        = uCR(iNodeX,iZ2,iZ3,iZ4,iE_G,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_SL_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: U_R ) &
    !$OMP MAP( release: iZ_B0, iZ_E0, GE, GX, U_F, U_R, &
    !$OMP               Limited, TroubledCell, &
    !$OMP               C_0, C_X1, C_X2, C_X3, &
    !$OMP               CL_X1, CL_X2, CL_X3,  &
    !$OMP               uCR_K, wSqrtGm, uCR )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC COPYOUT( U_R ) &
    !$ACC DELETE( iZ_B0, iZ_E0, GE, GX, U_F, &
    !$ACC         Limited, TroubledCell, &
    !$ACC         C_0, C_X1, C_X2, C_X3, &
    !$ACC         CL_X1, CL_X2, CL_X3,  &
    !$ACC         uCR_K, wSqrtGm, uCR )
#endif

    CALL TimersStop( Timer_SL )

  END SUBROUTINE ApplySlopeLimiter_TVD


  SUBROUTINE ComputeLimitedSlopes_X1( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, CL_X1 )

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      CL_X1(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nCR,1:nSpecies)

    INTEGER  :: &
      iE_G, iZ1, iZ2, iZ3, iZ4, iCR, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      C0(iZ_B0(2)-1:iZ_E0(2)+1, &
         iZ_B0(3)  :iZ_E0(3)  , &
         iZ_B0(4)  :iZ_E0(4)  , &
         1:nE_G,1:nCR,1:nSpecies), &
      C1(iZ_B0(2)-1:iZ_E0(2)+1, &
         iZ_B0(3)  :iZ_E0(3)  , &
         iZ_B0(4)  :iZ_E0(4)  , &
         1:nE_G,1:nCR,1:nSpecies)
    REAL(DP) :: &
      uCR(1:nDOFX, &
          iZ_B0(2)-1:iZ_E0(2)+1, &
          iZ_B0(3)  :iZ_E0(3)  , &
          iZ_B0(4)  :iZ_E0(4)  , &
          1:nE_G,1:nCR,1:nSpecies)

    nV_KX = PRODUCT( SHAPE( C0 ) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B0, iZ_E0, U_R ) &
    !$OMP MAP( alloc: CL_X1, C0, C1, uCR )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC COPYIN( iZ_B0, iZ_E0, U_R ) &
    !$ACC CREATE( CL_X1, C0, C1, uCR )
#endif

    ! --- Permute Radiation Fields ---

    CALL TimersStart( Timer_SL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iZ1, iNodeE, iNodeZ )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
    !$ACC PRIVATE( iZ1, iNodeE, iNodeZ ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, U_R, uCR )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iZ1, iNodeE, iNodeZ )
#endif
    DO iS     = 1         , nSpecies
    DO iCR    = 1         , nCR
    DO iE_G   = 1         , nE_G
    DO iZ4    = iZ_B0(4)  , iZ_E0(4)
    DO iZ3    = iZ_B0(3)  , iZ_E0(3)
    DO iZ2    = iZ_B0(2)-1, iZ_E0(2)+1
    DO iNodeX = 1         , nDOFX

      iZ1    = MOD( (iE_G-1) / nDOFE, nE    ) + iZ_B0(1)
      iNodeE = MOD( (iE_G-1)        , nDOFE ) + 1

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      uCR(iNodeX,iZ2,iZ3,iZ4,iE_G,iCR,iS) &
        = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_SL_Permute )

    ! --- Legendre Coefficients C0 and C1 ---

    CALL TimersStart( Timer_SL_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_0, 1, Zero, C0, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_1, 1, Zero, C1, 1 )

    CALL TimersStop( Timer_SL_LinearAlgebra )

    ! --- Limited Legendre Coefficient CL_X1 ---

    CALL TimersStart( Timer_SL_MinMod )

    ASSOCIATE( dX1 => MeshX(1) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP MAP( to: dX1 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC COPYIN( dX1 ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, CL_X1, C0, C1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      CL_X1(iZ2,iZ3,iZ4,iE_G,iCR,iS) &
        = MinModB &
            ( C1(iZ2,iZ3,iZ4,iE_G,iCR,iS), &
              BetaTVD * ( C0  (iZ2  ,iZ3,iZ4,iE_G,iCR,iS)    &
                          - C0(iZ2-1,iZ3,iZ4,iE_G,iCR,iS) ), &
              BetaTVD * ( C0  (iZ2+1,iZ3,iZ4,iE_G,iCR,iS)    &
                          - C0(iZ2  ,iZ3,iZ4,iE_G,iCR,iS) ), &
              dX1(iZ2), BetaTVB )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1

    CALL TimersStop( Timer_SL_MinMod )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: CL_X1 ) &
    !$OMP MAP( release: iZ_B0, iZ_E0, U_R, C0, C1, uCR )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC COPYOUT( CL_X1 ) &
    !$ACC DELETE( iZ_B0, iZ_E0, U_R, C0, C1, uCR )
#endif

  END SUBROUTINE ComputeLimitedSlopes_X1


  SUBROUTINE ComputeLimitedSlopes_X2( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, CL_X2 )

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      CL_X2(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nCR,1:nSpecies)

    INTEGER  :: &
      iE_G, iZ1, iZ2, iZ3, iZ4, iCR, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      C0 (iZ_B0(3)-1:iZ_E0(3)+1, &
          iZ_B0(2)  :iZ_E0(2)  , &
          iZ_B0(4)  :iZ_E0(4)  , &
          1:nE_G,1:nCR,1:nSpecies), &
      C2 (iZ_B0(3)-1:iZ_E0(3)+1, &
          iZ_B0(2)  :iZ_E0(2)  , &
          iZ_B0(4)  :iZ_E0(4)  , &
          1:nE_G,1:nCR,1:nSpecies), &
      CL2(iZ_B0(3)  :iZ_E0(3)  , &
          iZ_B0(2)  :iZ_E0(2)  , &
          iZ_B0(4)  :iZ_E0(4)  , &
          1:nE_G,1:nCR,1:nSpecies)
    REAL(DP) :: &
      uCR(1:nDOFX, &
          iZ_B0(3)-1:iZ_E0(3)+1, &
          iZ_B0(2)  :iZ_E0(2)  , &
          iZ_B0(4)  :iZ_E0(4)  , &
          1:nE_G,1:nCR,1:nSpecies)

    IF( iZ_E0(3) .EQ. iZ_B0(3) )THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
      !$ACC PRESENT( iZ_B0, iZ_E0, CL_X2 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(6)
#endif
      DO iS   = 1       , nSpecies
      DO iCR  = 1       , nCR
      DO iE_G = 1       , nE_G
      DO iZ4  = iZ_B0(4), iZ_E0(4)
      DO iZ3  = iZ_B0(3), iZ_E0(3)
      DO iZ2  = iZ_B0(2), iZ_E0(2)

        CL_X2(iZ2,iZ3,iZ4,iE_G,iCR,iS) = Zero

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

      RETURN

    END IF

    nV_KX = PRODUCT( SHAPE( C0 ) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B0, iZ_E0, U_R ) &
    !$OMP MAP( alloc: CL_X2, C0, C2, CL2, uCR )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC COPYIN( iZ_B0, iZ_E0, U_R ) &
    !$ACC CREATE( CL_X2, C0, C2, CL2, uCR )
#endif

    ! --- Permute Radiation Fields ---

    CALL TimersStart( Timer_SL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iZ1, iNodeE, iNodeZ )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
    !$ACC PRIVATE( iZ1, iNodeE, iNodeZ ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, U_R, uCR )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iZ1, iNodeE, iNodeZ )
#endif
    DO iS     = 1         , nSpecies
    DO iCR    = 1         , nCR
    DO iE_G   = 1         , nE_G
    DO iZ4    = iZ_B0(4)  , iZ_E0(4)
    DO iZ2    = iZ_B0(2)  , iZ_E0(2)
    DO iZ3    = iZ_B0(3)-1, iZ_E0(3)+1
    DO iNodeX = 1         , nDOFX

      iZ1    = MOD( (iE_G-1) / nDOFE, nE    ) + iZ_B0(1)
      iNodeE = MOD( (iE_G-1)        , nDOFE ) + 1

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      uCR(iNodeX,iZ3,iZ2,iZ4,iE_G,iCR,iS) &
        = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_SL_Permute )

    ! --- Legendre Coefficients C0 and C2 ---

    CALL TimersStart( Timer_SL_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_0, 1, Zero, C0, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_2, 1, Zero, C2, 1 )

    CALL TimersStop( Timer_SL_LinearAlgebra )

    ! --- Limited Legendre Coefficient CL_X2 ---

    CALL TimersStart( Timer_SL_MinMod )

    ASSOCIATE( dX2 => MeshX(2) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP MAP( to: dX2 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC COPYIN( dX2 ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, CL2, C0, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ2  = iZ_B0(2), iZ_E0(2)
    DO iZ3  = iZ_B0(3), iZ_E0(3)

      CL2(iZ3,iZ2,iZ4,iE_G,iCR,iS) &
        = MinModB &
            ( C2(iZ3,iZ2,iZ4,iE_G,iCR,iS), &
              BetaTVD * ( C0  (iZ3  ,iZ2,iZ4,iE_G,iCR,iS)    &
                          - C0(iZ3-1,iZ2,iZ4,iE_G,iCR,iS) ), &
              BetaTVD * ( C0  (iZ3+1,iZ2,iZ4,iE_G,iCR,iS)    &
                          - C0(iZ3  ,iZ2,iZ4,iE_G,iCR,iS) ), &
              dX2(iZ3), BetaTVB )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dX2

    CALL TimersStop( Timer_SL_MinMod )

    ! --- Permute Legendre Coefficient CL_X2 ---

    CALL TimersStart( Timer_SL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, CL_X2, CL2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      CL_X2(iZ2,iZ3,iZ4,iE_G,iCR,iS) &
        = CL2(iZ3,iZ2,iZ4,iE_G,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_SL_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: iZ_B0, iZ_E0, U_R, CL_X2, C0, C2, CL2, uCR )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC DELETE( iZ_B0, iZ_E0, U_R, CL_X2, C0, C2, CL2, uCR )
#endif

  END SUBROUTINE ComputeLimitedSlopes_X2


  SUBROUTINE ComputeLimitedSlopes_X3( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, CL_X3 )

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      CL_X3(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nCR,1:nSpecies)

    INTEGER  :: &
      iE_G, iZ1, iZ2, iZ3, iZ4, iCR, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      C0 (iZ_B0(4)-1:iZ_E0(4)+1, &
          iZ_B0(3)  :iZ_E0(3)  , &
          iZ_B0(2)  :iZ_E0(2)  , &
          1:nE_G,1:nCR,1:nSpecies), &
      C3 (iZ_B0(4)-1:iZ_E0(4)+1, &
          iZ_B0(3)  :iZ_E0(3)  , &
          iZ_B0(2)  :iZ_E0(2)  , &
          1:nE_G,1:nCR,1:nSpecies), &
      CL3(iZ_B0(4)  :iZ_E0(4)  , &
          iZ_B0(3)  :iZ_E0(3)  , &
          iZ_B0(2)  :iZ_E0(2)  , &
          1:nE_G,1:nCR,1:nSpecies)
    REAL(DP) :: &
      uCR(1:nDOFX, &
          iZ_B0(4)-1:iZ_E0(4)+1, &
          iZ_B0(3)  :iZ_E0(3)  , &
          iZ_B0(2)  :iZ_E0(2)  , &
          1:nE_G,1:nCR,1:nSpecies)

    IF( iZ_E0(4) .EQ. iZ_B0(4) )THEN

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, CL_X3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
      DO iS   = 1       , nSpecies
      DO iCR  = 1       , nCR
      DO iE_G = 1       , nE_G
      DO iZ4  = iZ_B0(4), iZ_E0(4)
      DO iZ3  = iZ_B0(3), iZ_E0(3)
      DO iZ2  = iZ_B0(2), iZ_E0(2)

        CL_X3(iZ2,iZ3,iZ4,iE_G,iCR,iS) = Zero

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

      RETURN

    END IF

    nV_KX = PRODUCT( SHAPE( C0 ) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B0, iZ_E0, U_R ) &
    !$OMP MAP( alloc: CL_X3, C0, C3, CL3, uCR )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC COPYIN( iZ_B0, iZ_E0, U_R ) &
    !$ACC CREATE( CL_X3, C0, C3, CL3, uCR )
#endif

    ! --- Permute Radiation Fields ---

    CALL TimersStart( Timer_SL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iZ1, iNodeE, iNodeZ )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
    !$ACC PRIVATE( iZ1, iNodeE, iNodeZ ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, U_R, uCR )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iZ1, iNodeE, iNodeZ )
#endif
    DO iS     = 1         , nSpecies
    DO iCR    = 1         , nCR
    DO iE_G   = 1         , nE_G
    DO iZ2    = iZ_B0(2)  , iZ_E0(2)
    DO iZ3    = iZ_B0(3)  , iZ_E0(3)
    DO iZ4    = iZ_B0(4)-1, iZ_E0(4)+1
    DO iNodeX = 1         , nDOFX

      iZ1    = MOD( (iE_G-1) / nDOFE, nE    ) + iZ_B0(1)
      iNodeE = MOD( (iE_G-1)        , nDOFE ) + 1

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      uCR(iNodeX,iZ4,iZ3,iZ2,iE_G,iCR,iS) &
        = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_SL_Permute )

    ! --- Legendre Coefficients C0 and C3 ---

    CALL TimersStart( Timer_SL_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_0, 1, Zero, C0, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_3, 1, Zero, C3, 1 )

    CALL TimersStop( Timer_SL_LinearAlgebra )

    ! --- Limited Legendre Coefficient CL_X3 ---

    CALL TimersStart( Timer_SL_MinMod )

    ASSOCIATE( dX3 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP MAP( to: dX3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC COPYIN( dX3 ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, CL3, C0, C3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iZ2  = iZ_B0(2), iZ_E0(2)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ4  = iZ_B0(4), iZ_E0(4)

      CL3(iZ4,iZ3,iZ2,iE_G,iCR,iS) &
        = MinModB &
            ( C3(iZ4,iZ3,iZ2,iE_G,iCR,iS), &
              BetaTVD * ( C0  (iZ4  ,iZ3,iZ2,iE_G,iCR,iS)    &
                          - C0(iZ4-1,iZ3,iZ2,iE_G,iCR,iS) ), &
              BetaTVD * ( C0  (iZ4+1,iZ3,iZ2,iE_G,iCR,iS)    &
                          - C0(iZ4  ,iZ3,iZ2,iE_G,iCR,iS) ), &
              dX3(iZ4), BetaTVB )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dX3

    CALL TimersStop( Timer_SL_MinMod )

    ! --- Permute Legendre Coefficient CL_X3 ---

    CALL TimersStart( Timer_SL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, CL_X3, CL3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      CL_X3(iZ2,iZ3,iZ4,iE_G,iCR,iS) &
        = CL3(iZ4,iZ3,iZ2,iE_G,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_SL_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: iZ_B0, iZ_E0, U_R, CL_X3, C0, C3, CL3, uCR )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC DELETE( iZ_B0, iZ_E0, U_R, CL_X3, C0, C3, CL3, uCR )
#endif

  END SUBROUTINE ComputeLimitedSlopes_X3


  SUBROUTINE ApplySLopeLimiter_WENO &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE (1:nDOFE, &
          iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX, &
          iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    PRINT*, "      ApplySlopeLimiter_WENO"

    STOP

  END SUBROUTINE ApplySLopeLimiter_WENO


END MODULE TwoMoment_SlopeLimiterModule_OrderV
