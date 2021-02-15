MODULE TwoMoment_SlopeLimiterModule_OrderV

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFE, nDOFX, nDimsX
  USE TwoMoment_TimersModule_OrderV, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_SlopeLimiter
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply
  USE UtilitiesModule, ONLY: &
    MinModB
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModalX, &
    MapModalToNodalX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF
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
  REAL(DP)     :: BetaTVD, BetaTVB
  REAL(DP)     :: SlopeTolerance

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

  END SUBROUTINE InitializeSlopeLimiter_TwoMoment


  SUBROUTINE FinalizeSlopeLimiter_TwoMoment

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

    CALL TimersStart( Timer_SlopeLimiter )

    IF( PRESENT( SuppressBC_Option ) )THEN
      SuppressBC = SuppressBC_Option
    ELSE
      SuppressBC = .FALSE.
    END IF

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

    CALL TimersStop( Timer_SlopeLimiter )

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
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER  :: &
      iX1, iX2, iX3, iZ1, iCR, iS, iE_G, nE_G, &
      iNodeZ, iNodeE, iNodeX, nV_KX, iDim
    REAL(DP) :: &
      dX1, dX2, dX3
    REAL(DP) :: &
      SlopeDifference(1:nCR)
    REAL(DP) :: &
      uCR_M(1:nCR,1:nDOFX), dCR(1:nCR,1:nDimsX)
    REAL(DP) :: &
      uCR  (1:nDOFX, &
            1:nCR,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: &
      uCR_K(1:nCR,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)

    iX_B0 = iZ_B0(2:4)
    iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4)
    iX_E1 = iZ_E1(2:4)

    nE_G = ( iZ_E0(1) - iZ_B0(1) + 1 ) * nDOFE

    ! --- Reorder Data ---

    DO iS     = 1, nSpecies
    DO iZ1    = iZ_B0(1), iZ_E0(1)
    DO iNodeE = 1, nDOFE
    DO iX3    = iX_B1(3), iX_E1(3)
    DO iX2    = iX_B1(2), iX_E1(2)
    DO iX1    = iX_B1(1), iX_E1(1)
    DO iCR    = 1, nCR
    DO iNodeX = 1, nDOFX

      iNodeZ = (iNodeX-1) * nDOFE + iNodeE 
      iE_G   = (iZ1   -1) * nDOFE + iNodeE

      uCR(iNodeX,iCR,iX1,iX2,iX3,iE_G,iS) &
        = U_R(iNodeZ,iZ1,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Compute Cell Averages ---

    ! --- Variables Per Spatial Element ---

    nV_KX = nCR * PRODUCT(iX_E1-iX_B1+1) * nE_G * nSpecies

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nV_KX, One, uCR, nDOFX, WeightsX_q, 1, Zero, uCR_K, 1 )

    DO iS   = 1, nSpecies
    DO iE_G = 1, nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      dX1 = MeshX(1) % Width(iX1)
      dX2 = MeshX(2) % Width(iX2)
      dX3 = MeshX(3) % Width(iX3)

      ! --- Map to Modal Representation ---

      DO iCR = 1, nCR

        CALL MapNodalToModalX &
               ( uCR(1:nDOFX,iCR,iX1,iX2,iX3,iE_G,iS), uCR_M(iCR,1:nDOFX) )

      END DO

      ! --- Compute Limited Slopes ---

      dCR(:,1) &
        = MinModB( uCR_M(:,2), &
                   BetaTVD * ( uCR_K  (:,iX1  ,iX2,iX3,iE_G,iS)    &
                               - uCR_K(:,iX1-1,iX2,iX3,iE_G,iS) ), &
                   BetaTVD * ( uCR_K  (:,iX1+1,iX2,iX3,iE_G,iS)    &
                               - uCR_K(:,iX1  ,iX2,iX3,iE_G,iS) ), &
                   dX1, BetaTVB )

      IF( nDimsX > 1 )THEN

        dCR(:,2) &
          = MinModB( uCR_M(:,3), &
                     BetaTVD * ( uCR_K  (:,iX1,iX2  ,iX3,iE_G,iS)    &
                                 - uCR_K(:,iX1,iX2-1,iX3,iE_G,iS) ), &
                     BetaTVD * ( uCR_K  (:,iX1,iX2+1,iX3,iE_G,iS)    &
                                 - uCR_K(:,iX1,iX2  ,iX3,iE_G,iS) ), &
                     dX2, BetaTVB )

      END IF

      IF( nDimsX > 2 )THEN

        dCR(:,3) &
          = MinModB( uCR_M(:,4), &
                     BetaTVD * ( uCR_K  (:,iX1,iX2,iX3  ,iE_G,iS)    &
                                 - uCR_K(:,iX1,iX2,iX3-1,iE_G,iS) ), &
                     BetaTVD * ( uCR_K  (:,iX1,iX2,iX3+1,iE_G,iS)    &
                                 - uCR_K(:,iX1,iX2,iX3  ,iE_G,iS) ), &
                     dX3, BetaTVB )

      END IF

      ! --- Compare Limited Slopes to Original Slopes ---

      SlopeDifference = Zero

      DO iDim = 1, nDimsX
      DO iCR  = 1, nCR

        SlopeDifference(iCR) &
          = MAX( SlopeDifference(iCR), &
                 ABS( dCR(iCR,iDim) - uCR_M(iCR,1+iDim) ) )

      END DO
      END DO

      ! --- Replace Slopes and Discard High-Order Components ---
      ! --- if Limited Slopes Deviate too Much from Original ---

      DO iCR = 1, nCR

        IF( SlopeDifference(iCR) &
              > SlopeTolerance * ABS( uCR_K(iCR,iX1,iX2,iX3,iE_G,iS) ) )THEN

          uCR_M(iCR,2:nDOFX) = Zero

          DO iDim = 1, nDimsX

            uCR_M(iCR,1+iDim) = dCR(iCR,iDim)

          END DO

          CALL MapModalToNodalX &
                 ( uCR(:,iCR,iX1,iX2,iX3,iE_G,iS), uCR_M(iCR,:) )

        END IF

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Reorder Data ---

    DO iS     = 1, nSpecies
    DO iCR    = 1, nCR
    DO iX3    = iX_B0(3), iX_E0(3)
    DO iX2    = iX_B0(2), iX_E0(2)
    DO iX1    = iX_B0(1), iX_E0(1)
    DO iZ1    = iZ_B0(1), iZ_E0(1)
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iNodeZ = (iNodeX-1) * nDOFE + iNodeE
      iE_G   = (iZ1   -1) * nDOFE + iNodeE

      U_R(iNodeZ,iZ1,iX1,iX2,iX3,iCR,iS) &
        = uCR(iNodeX,iCR,iX1,iX2,iX3,iE_G,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

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


END MODULE TwoMoment_SlopeLimiterModule_OrderV
