MODULE TwoMoment_SlopeLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFE, nDOFX, nDimsX
  USE TwoMoment_TimersModule, ONLY: &
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
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( always, to: &
    !$OMP      N2M_Vec_0, N2M_Vec_1, N2M_Vec_2, N2M_Vec_3, &
    !$OMP      M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3, &
    !$OMP      BetaTVD, BetaTVB, SlopeTolerance )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE &
    !$ACC DEVICE( N2M_Vec_0, N2M_Vec_1, N2M_Vec_2, N2M_Vec_3, &
    !$ACC         M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3, &
    !$ACC         BetaTVD, BetaTVB, SlopeTolerance )
#endif

  END SUBROUTINE InitializeSlopeLimiter_TwoMoment


  SUBROUTINE FinalizeSlopeLimiter_TwoMoment

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: &
    !$OMP      N2M_Vec_0, N2M_Vec_1, N2M_Vec_2, N2M_Vec_3, &
    !$OMP      M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3, &
    !$OMP      BetaTVD, BetaTVB, SlopeTolerance )
#endif

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

    CALL TimersStart( Timer_SL )

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

    IF ( TRIM( SlopeLimiterMethod ) == 'WENO' ) THEN

      CALL ApplySlopeLimiter_WENO &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    ELSE

      CALL ApplySlopeLimiter_TVD &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: iZ_B0, iZ_E0, GE, GX, U_F, U_R )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC DELETE( iZ_B0, iZ_E0, GE, GX, U_F, U_R )

    !$ACC WAIT
#endif

    CALL TimersStop( Timer_SL )

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
      iNodeZ, iNodeE, iNodeX
    REAL(DP) :: &
      dSlope, wSqrtGm, Alpha, uCR_K, &
      C_0, C_X1, C_X2, C_X3, &
      C0_L, C0_R, CL_X1, CL_X2, CL_X3
    REAL(DP) :: &
      uCR(1:nDOFX)
    LOGICAL  :: &
      TroubledCell &
           (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)

    nE   = iZ_E0(1) - iZ_B0(1) + 1
    nE_G = nE * nDOFE

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B0, iZ_E0, GE, GX, U_F, U_R ) &
    !$OMP MAP( alloc: TroubledCell )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC COPYIN( iZ_B0, iZ_E0, GE, GX, U_F, U_R ) &
    !$ACC CREATE( TroubledCell )
#endif

    CALL DetectTroubledCells_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, TroubledCell )

    CALL TimersStart( Timer_SL_ReplaceSlopes )

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP MAP( to: dX1, dX2, dX3 ) &
    !$OMP PRIVATE( iNodeZ, iNodeE, iZ1, dSlope, wSqrtGm, &
    !$OMP          C_0, C_X1, C_X2, C_X3, uCR, uCR_K, Alpha, &
    !$OMP          CL_X1, CL_X2, CL_X3, C0_L, C0_R )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC COPYIN( dX1, dX2, dX3 ) &
    !$ACC PRIVATE( iNodeZ, iNodeE, iZ1, dSlope, wSqrtGm, &
    !$ACC          C_0, C_X1, C_X2, C_X3, uCR, uCR_K, Alpha, &
    !$ACC          CL_X1, CL_X2, CL_X3, C0_L, C0_R ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, TroubledCell, WeightsX_q, GX, U_R, &
    !$ACC          M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ, iNodeE, iZ1, dSlope, wSqrtGm, &
    !$OMP          C_0, C_X1, C_X2, C_X3, uCR, uCR_K, Alpha, &
    !$OMP          CL_X1, CL_X2, CL_X3, C0_L, C0_R )
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      IF( TroubledCell(iZ2,iZ3,iZ4,iE_G,iS) )THEN

        iZ1    = MOD( (iE_G-1) / nDOFE, nE    ) + iZ_B0(1)
        iNodeE = MOD( (iE_G-1)        , nDOFE ) + 1

        DO iNodeX = 1, nDOFX
          iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE
          uCR(iNodeX) = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)
        END DO

        ! --- Compute Legendre Coefficients ---

        C_0  = Zero
        C_X1 = Zero
        C_X2 = Zero
        C_X3 = Zero
        DO iNodeX = 1, nDOFX
          C_0  = C_0  + N2M_Vec_0(iNodeX) * uCR(iNodeX)
          C_X1 = C_X1 + N2M_Vec_1(iNodeX) * uCR(iNodeX)
          C_X2 = C_X2 + N2M_Vec_2(iNodeX) * uCR(iNodeX)
          C_X3 = C_X3 + N2M_Vec_3(iNodeX) * uCR(iNodeX)
        END DO

        ! --- Limited Legendre Coefficients ---

        C0_L = Zero
        DO iNodeX = 1, nDOFX
          iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE
          C0_L = C0_L + N2M_Vec_0(iNodeX) * U_R(iNodeZ,iZ1,iZ2-1,iZ3,iZ4,iCR,iS)
        END DO
        C0_R = Zero
        DO iNodeX = 1, nDOFX
          iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE
          C0_R = C0_R + N2M_Vec_0(iNodeX) * U_R(iNodeZ,iZ1,iZ2+1,iZ3,iZ4,iCR,iS)
        END DO
        CL_X1 &
          = MinModB( C_X1, &
                     BetaTVD * ( C_0 - C0_L ), &
                     BetaTVD * ( C0_R - C_0 ), &
                     dX1(iZ2), BetaTVB )

        IF ( iZ_B0(3) /= iZ_E0(3) ) THEN
          C0_L = Zero
          DO iNodeX = 1, nDOFX
            iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE
            C0_L = C0_L + N2M_Vec_0(iNodeX) * U_R(iNodeZ,iZ1,iZ2,iZ3-1,iZ4,iCR,iS)
          END DO
          C0_R = Zero
          DO iNodeX = 1, nDOFX
            iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE
            C0_R = C0_R + N2M_Vec_0(iNodeX) * U_R(iNodeZ,iZ1,iZ2,iZ3+1,iZ4,iCR,iS)
          END DO
          CL_X2 &
            = MinModB( C_X2, &
                       BetaTVD * ( C_0 - C0_L ), &
                       BetaTVD * ( C0_R - C_0 ), &
                       dX2(iZ3), BetaTVB )
        ELSE
          CL_X2 = C_X2
        END IF

        IF ( iZ_B0(4) /= iZ_E0(4) ) THEN
          C0_L = Zero
          DO iNodeX = 1, nDOFX
            iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE
            C0_L = C0_L + N2M_Vec_0(iNodeX) * U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4-1,iCR,iS)
          END DO
          C0_R = Zero
          DO iNodeX = 1, nDOFX
            iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE
            C0_R = C0_R + N2M_Vec_0(iNodeX) * U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4+1,iCR,iS)
          END DO
          CL_X3 &
            = MinModB( C_X3, &
                       BetaTVD * ( C_0 - C0_L ), &
                       BetaTVD * ( C0_R - C_0 ), &
                       dX3(iZ4), BetaTVB )
        ELSE
          CL_X3 = C_X3
        END IF

        dSlope &
          = MAX( ABS( CL_X1 - C_X1 ), &
                 ABS( CL_X2 - C_X2 ), &
                 ABS( CL_X3 - C_X3 ) )

        IF( dSlope > SlopeTolerance * ABS( C_0 ) )THEN

          ! --- Conservative Correction ---

          uCR_K = Zero
          Alpha = Zero

          DO iNodeX = 1, nDOFX

            wSqrtGm = WeightsX_q(iNodeX) * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

            uCR_K = uCR_K + wSqrtGm * uCR(iNodeX)

            uCR(iNodeX) &
              =   M2N_Vec_0(iNodeX) * C_0 &
                + M2N_Vec_1(iNodeX) * CL_X1 &
                + M2N_Vec_2(iNodeX) * CL_X2 &
                + M2N_Vec_3(iNodeX) * CL_X3

            Alpha = Alpha + wSqrtGm * uCR(iNodeX)

          END DO

          IF( ABS( Alpha ) > Zero )THEN
            DO iNodeX = 1, nDOFX
              iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE
              U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) = uCR(iNodeX) * uCR_K / Alpha
            END DO
          END IF

        END IF

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

    CALL TimersStop( Timer_SL_ReplaceSlopes )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: iZ_B0, iZ_E0, GE, GX, U_F, U_R, &
    !$OMP               TroubledCell )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC DELETE( iZ_B0, iZ_E0, GE, GX, U_F, U_R, &
    !$ACC         TroubledCell )
#endif

  END SUBROUTINE ApplySlopeLimiter_TVD


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


END MODULE TwoMoment_SlopeLimiterModule
