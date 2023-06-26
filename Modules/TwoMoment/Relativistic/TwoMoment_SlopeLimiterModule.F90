MODULE TwoMoment_SlopeLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFE, nDOFX, nDimsX
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

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_TwoMoment
  PUBLIC :: FinalizeSlopeLimiter_TwoMoment
  PUBLIC :: ApplySlopeLimiter_TwoMoment

  CHARACTER(4) :: SlopeLimiterMethod
  LOGICAL      :: UseSlopeLimiter
  INTEGER      :: iX_B0(3), iX_E0(3)
  INTEGER      :: iE_B0, iE_E0, nE_G
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

  END SUBROUTINE InitializeSlopeLimiter_TwoMoment


  SUBROUTINE FinalizeSlopeLimiter_TwoMoment

    DEALLOCATE( N2M_Vec_0, N2M_Vec_1, N2M_Vec_2, N2M_Vec_3 )
    DEALLOCATE( M2N_Vec_0, M2N_Vec_1, M2N_Vec_2, M2N_Vec_3 )

  END SUBROUTINE FinalizeSlopeLimiter_TwoMoment


  SUBROUTINE ApplySlopeLimiter_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, SuppressBC_Option, Verbose_Option )

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
    LOGICAL,  INTENT(in), OPTIONAL :: &
      Verbose_Option

    LOGICAL :: SuppressBC
    LOGICAL :: Verbose

    IF( .NOT. UseSlopeLimiter .OR. nDOFX == 1 ) RETURN


    IF( PRESENT( SuppressBC_Option ) )THEN
      SuppressBC = SuppressBC_Option
    ELSE
      SuppressBC = .FALSE.
    END IF


    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    END IF



    IF (Verbose) THEN

      PRINT*, "      ApplySlopeLimiter_TwoMoment"

    END IF


!   IF( .NOT. SuppressBC )THEN
!
!      CALL ApplyBoundaryConditions_TwoMoment &
!             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R )
!
!    END IF

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
      iE, iX1, iX2, iX3, iCR, iS, iE_G, &
      iNodeZ, iNodeE, iNodeX, nV_KX
    REAL(DP) :: &
      dSlope, Alpha
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
      uCR_K(iZ_B1(2):iZ_E1(2), &
            iZ_B1(3):iZ_E1(3), &
            iZ_B1(4):iZ_E1(4), &
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


    iE_B0 = iZ_B0(1); iX_B0 = iZ_B0(2:4)
    iE_E0 = iZ_E0(1); iX_E0 = iZ_E0(2:4)

    nE_G = ( iZ_E0(1) - iZ_B0(1) + 1 ) * nDOFE
    CALL ComputeLimitedSlopes_X1( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, CL_X1 )

    CALL ComputeLimitedSlopes_X2( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, CL_X2 )

    CALL ComputeLimitedSlopes_X3( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, CL_X3 )

    ! --- Permute Radiation Fields ---


#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(8) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#endif
    DO iS     = 1       , nSpecies
    DO iCR    = 1       , nCR
    DO iE     = iE_B0   , iE_E0
    DO iNodeE = 1       , nDOFE
    DO iX3    = iX_B0(3), iX_E0(3)
    DO iX2    = iX_B0(2), iX_E0(2)
    DO iX1    = iX_B0(1), iX_E0(1)
    DO iNodeX = 1       , nDOFX

      iNodeZ = (iNodeX-1) * nDOFE + iNodeE 
      iE_G   = (iE    -1) * nDOFE + iNodeE

      uCR(iNodeX,iX1,iX2,iX3,iE_G,iCR,iS) &
        = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3    = iX_B0(3), iX_E0(3)
    DO iX2    = iX_B0(2), iX_E0(2)
    DO iX1    = iX_B0(1), iX_E0(1)
    DO iNodeX = 1       , nDOFX

      wSqrtGm(iNodeX,iX1,iX2,iX3) &
        = WeightsX_q(iNodeX) * GX(iNodeX,iX1,iX2,iX3,iGF_SqrtGm)

    END DO
    END DO
    END DO
    END DO


    ! --- Compute Cell Average ---

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      uCR_K(iX1,iX2,iX3,iE_G,iCR,iS) &
        = SUM( wSqrtGm(:,iX1,iX2,iX3) * uCR(:,iX1,iX2,iX3,iE_G,iCR,iS) )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Compute Legendre Coefficients ---


    nV_KX = PRODUCT( SHAPE( C_0 ) )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_0, 1, Zero, C_0 , 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_1, 1, Zero, C_X1, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_2, 1, Zero, C_X2, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_3, 1, Zero, C_X3, 1 )



#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( dSlope )
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      Limited(iX1,iX2,iX3,iE_G,iCR,iS) = .FALSE.


        dSlope &
          = MAX( ABS( CL_X1(iX1,iX2,iX3,iE_G,iCR,iS) &
                      -C_X1(iX1,iX2,iX3,iE_G,iCR,iS) ), &
                 ABS( CL_X2(iX1,iX2,iX3,iE_G,iCR,iS) &
                      -C_X2(iX1,iX2,iX3,iE_G,iCR,iS) ), &
                 ABS( CL_X3(iX1,iX2,iX3,iE_G,iCR,iS) &
                      -C_X3(iX1,iX2,iX3,iE_G,iCR,iS) ) )

        IF( dSlope > SlopeTolerance * ABS( C_0(iX1,iX2,iX3,iE_G,iCR,iS) ) )THEN

          uCR(:,iX1,iX2,iX3,iE_G,iCR,iS) &
            =   M2N_Vec_0(:) * C_0  (iX1,iX2,iX3,iE_G,iCR,iS) &
              + M2N_Vec_1(:) * CL_X1(iX1,iX2,iX3,iE_G,iCR,iS) &
              + M2N_Vec_2(:) * CL_X2(iX1,iX2,iX3,iE_G,iCR,iS) &
              + M2N_Vec_3(:) * CL_X3(iX1,iX2,iX3,iE_G,iCR,iS)

          Limited(iX1,iX2,iX3,iE_G,iCR,iS) = .TRUE.

        END IF


    END DO
    END DO
    END DO
    END DO
    END DO
    END DO


    ! --- Conservative Correction ---


    IF( ANY( Limited ) )THEN

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( Alpha )
#endif
      DO iS   = 1       , nSpecies
      DO iCR  = 1       , nCR
      DO iE_G = 1       , nE_G
      DO iX3  = iX_B0(3), iX_E0(3)
      DO iX2  = iX_B0(2), iX_E0(2)
      DO iX1  = iX_B0(1), iX_E0(1)

        IF( Limited(iX1,iX2,iX3,iE_G,iCR,iS) )THEN

          Alpha &
            = SUM( wSqrtGm(:,iX1,iX2,iX3) * uCR(:,iX1,iX2,iX3,iE_G,iCR,iS) )

          IF( ABS( Alpha ) > Zero )THEN

            Alpha = uCR_K(iX1,iX2,iX3,iE_G,iCR,iS) / Alpha

            uCR(:,iX1,iX2,iX3,iE_G,iCR,iS) &
              = Alpha * uCR(:,iX1,iX2,iX3,iE_G,iCR,iS)

          END IF

        END IF

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

    END IF



#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(8) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#endif
    DO iS     = 1, nSpecies
    DO iCR    = 1, nCR
    DO iX3    = iX_B0(3), iX_E0(3)
    DO iX2    = iX_B0(2), iX_E0(2)
    DO iX1    = iX_B0(1), iX_E0(1)
    DO iE     = iE_B0   , iE_E0
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iNodeZ = (iNodeX-1) * nDOFE + iNodeE
      iE_G   = (iE    -1) * nDOFE + iNodeE

      U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) &
        = uCR(iNodeX,iX1,iX2,iX3,iE_G,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO


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
      CL_X1(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nCR,1:nSpecies)

    INTEGER  :: &
      iE, iE_G, iX1, iX2, iX3, iCR, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      C0(iX_B0(1)-1:iX_E0(1)+1, &
         iX_B0(2)  :iX_E0(2)  , &
         iX_B0(3)  :iX_E0(3)  , &
         1:nE_G,1:nCR,1:nSpecies), &
      C1(iX_B0(1)-1:iX_E0(1)+1, &
         iX_B0(2)  :iX_E0(2)  , &
         iX_B0(3)  :iX_E0(3)  , &
         1:nE_G,1:nCR,1:nSpecies)
    REAL(DP) :: &
      uCR(1:nDOFX, &
          iX_B0(1)-1:iX_E0(1)+1, &
          iX_B0(2)  :iX_E0(2)  , &
          iX_B0(3)  :iX_E0(3)  , &
          1:nE_G,1:nCR,1:nSpecies)

    ! --- Permute Radiation Fields ---


#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(8) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#endif




    IF( iZ_E0(2) .EQ. iZ_B0(2) )THEN
      CL_X1 = Zero
 
      RETURN
    END IF



    DO iS     = 1         , nSpecies
    DO iCR    = 1         , nCR
    DO iE     = iE_B0     , iE_E0
    DO iNodeE = 1         , nDOFE
    DO iX3    = iX_B0(3)  , iX_E0(3)
    DO iX2    = iX_B0(2)  , iX_E0(2)
    DO iX1    = iX_B0(1)-1, iX_E0(1)+1
    DO iNodeX = 1         , nDOFX

      iNodeZ = (iNodeX-1) * nDOFE + iNodeE 
      iE_G   = (iE    -1) * nDOFE + iNodeE

      uCR(iNodeX,iX1,iX2,iX3,iE_G,iCR,iS) &
        = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO


    ! --- Legendre Coefficients C0 and C1 ---


    nV_KX = PRODUCT( SHAPE( C0 ) )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_0, 1, Zero, C0, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_1, 1, Zero, C1, 1 )


    ! --- Limited Legendre Coefficient CL_X1 ---


    ASSOCIATE( dX1 => MeshX(1) % Width )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      CL_X1(iX1,iX2,iX3,iE_G,iCR,iS) &
        = MinModB &
            ( C1(iX1,iX2,iX3,iE_G,iCR,iS), &
              BetaTVD * ( C0  (iX1  ,iX2,iX3,iE_G,iCR,iS)    &
                          - C0(iX1-1,iX2,iX3,iE_G,iCR,iS) ), &
              BetaTVD * ( C0  (iX1+1,iX2,iX3,iE_G,iCR,iS)    &
                          - C0(iX1  ,iX2,iX3,iE_G,iCR,iS) ), &
              dX1(iX1), BetaTVB )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1


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
      CL_X2(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nCR,1:nSpecies)

    INTEGER  :: &
      iE, iE_G, iX1, iX2, iX3, iCR, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      C0 (iX_B0(2)-1:iX_E0(2)+1, &
          iX_B0(1)  :iX_E0(1)  , &
          iX_B0(3)  :iX_E0(3)  , &
          1:nE_G,1:nCR,1:nSpecies), &
      C2 (iX_B0(2)-1:iX_E0(2)+1, &
          iX_B0(1)  :iX_E0(1)  , &
          iX_B0(3)  :iX_E0(3)  , &
          1:nE_G,1:nCR,1:nSpecies), &
      CL2(iX_B0(2)  :iX_E0(2)  , &
          iX_B0(1)  :iX_E0(1)  , &
          iX_B0(3)  :iX_E0(3)  , &
          1:nE_G,1:nCR,1:nSpecies)
    REAL(DP) :: &
      uCR(1:nDOFX, &
          iX_B0(2)-1:iX_E0(2)+1, &
          iX_B0(1)  :iX_E0(1)  , &
          iX_B0(3)  :iX_E0(3)  , &
          1:nE_G,1:nCR,1:nSpecies)

    ! --- Permute Radiation Fields ---


#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(8) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#endif
    IF( iZ_E0(3) .EQ. iZ_B0(3) )THEN
      CL_X2 = Zero
 
      RETURN
    END IF

    DO iS     = 1         , nSpecies
    DO iCR    = 1         , nCR
    DO iE     = iE_B0     , iE_E0
    DO iNodeE = 1         , nDOFE
    DO iX3    = iX_B0(3)  , iX_E0(3)
    DO iX1    = iX_B0(1)  , iX_E0(1)
    DO iX2    = iX_B0(2)-1, iX_E0(2)+1
    DO iNodeX = 1         , nDOFX

      iNodeZ = (iNodeX-1) * nDOFE + iNodeE 
      iE_G   = (iE    -1) * nDOFE + iNodeE

      uCR(iNodeX,iX2,iX1,iX3,iE_G,iCR,iS) &
        = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO



    nV_KX = PRODUCT( SHAPE( C0 ) )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_0, 1, Zero, C0, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_2, 1, Zero, C2, 1 )


    ! --- Limited Legendre Coefficient CL_X2 ---


    ASSOCIATE( dX2 => MeshX(2) % Width )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX1  = iX_B0(1), iX_E0(1)
    DO iX2  = iX_B0(2), iX_E0(2)

      CL2(iX2,iX1,iX3,iE_G,iCR,iS) &
        = MinModB &
            ( C2(iX2,iX1,iX3,iE_G,iCR,iS), &
              BetaTVD * ( C0  (iX2  ,iX1,iX3,iE_G,iCR,iS)    &
                          - C0(iX2-1,iX1,iX3,iE_G,iCR,iS) ), &
              BetaTVD * ( C0  (iX2+1,iX1,iX3,iE_G,iCR,iS)    &
                          - C0(iX2  ,iX1,iX3,iE_G,iCR,iS) ), &
              dX2(iX2), BetaTVB )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dX2


    ! --- Permute Legendre Coefficient CL_X2 ---


#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      CL_X2(iX1,iX2,iX3,iE_G,iCR,iS) &
        = CL2(iX2,iX1,iX3,iE_G,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO


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
      CL_X3(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nCR,1:nSpecies)

    INTEGER  :: &
      iE, iE_G, iX1, iX2, iX3, iCR, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      C0 (iX_B0(3)-1:iX_E0(3)+1, &
          iX_B0(2)  :iX_E0(2)  , &
          iX_B0(1)  :iX_E0(1)  , &
          1:nE_G,1:nCR,1:nSpecies), &
      C3 (iX_B0(3)-1:iX_E0(3)+1, &
          iX_B0(2)  :iX_E0(2)  , &
          iX_B0(1)  :iX_E0(1)  , &
          1:nE_G,1:nCR,1:nSpecies), &
      CL3(iX_B0(3)  :iX_E0(3)  , &
          iX_B0(2)  :iX_E0(2)  , &
          iX_B0(1)  :iX_E0(1)  , &
          1:nE_G,1:nCR,1:nSpecies)
    REAL(DP) :: &
      uCR(1:nDOFX, &
          iX_B0(3)-1:iX_E0(3)+1, &
          iX_B0(2)  :iX_E0(2)  , &
          iX_B0(1)  :iX_E0(1)  , &
          1:nE_G,1:nCR,1:nSpecies)

    ! --- Permute Radiation Fields ---


#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(8) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#endif
    IF( iZ_E0(4) .EQ. iZ_B0(4) )THEN
      CL_X3 = Zero
 
      RETURN
    END IF

    DO iS     = 1         , nSpecies
    DO iCR    = 1         , nCR
    DO iE     = iE_B0     , iE_E0
    DO iNodeE = 1         , nDOFE
    DO iX1    = iX_B0(1)  , iX_E0(1)
    DO iX2    = iX_B0(2)  , iX_E0(2)
    DO iX3    = iX_B0(3)-1, iX_E0(3)+1
    DO iNodeX = 1         , nDOFX

      iNodeZ = (iNodeX-1) * nDOFE + iNodeE 
      iE_G   = (iE    -1) * nDOFE + iNodeE

      uCR(iNodeX,iX3,iX2,iX1,iE_G,iCR,iS) &
        = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO


    ! --- Legendre Coefficients C0 and C3 ---


    nV_KX = PRODUCT( SHAPE( C0 ) )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_0, 1, Zero, C0, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, N2M_Vec_3, 1, Zero, C3, 1 )


    ! --- Limited Legendre Coefficient CL_X3 ---


    ASSOCIATE( dX3 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iX1  = iX_B0(1), iX_E0(1)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX3  = iX_B0(3), iX_E0(3)

      CL3(iX3,iX2,iX1,iE_G,iCR,iS) &
        = MinModB &
            ( C3(iX3,iX2,iX1,iE_G,iCR,iS), &
              BetaTVD * ( C0  (iX3  ,iX2,iX1,iE_G,iCR,iS)    &
                          - C0(iX3-1,iX2,iX1,iE_G,iCR,iS) ), &
              BetaTVD * ( C0  (iX3+1,iX2,iX1,iE_G,iCR,iS)    &
                          - C0(iX3  ,iX2,iX1,iE_G,iCR,iS) ), &
              dX3(iX3), BetaTVB )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dX3


    ! --- Permute Legendre Coefficient CL_X3 ---


#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS   = 1       , nSpecies
    DO iCR  = 1       , nCR
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      CL_X3(iX1,iX2,iX3,iE_G,iCR,iS) &
        = CL3(iX3,iX2,iX1,iE_G,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO


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


END MODULE TwoMoment_SlopeLimiterModule
