MODULE TwoMoment_TroubledCellIndicatorModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFE, nDOFX, nDimsX
  USE TwoMoment_TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_TCI, &
    Timer_TCI_Permute, &
    Timer_TCI_LinearAlgebra, &
    Timer_TCI_Compute
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3, &
    WeightsX_q
  USE PolynomialBasisModuleX_Lagrange, ONLY: &
    L_X1, L_X2, L_X3
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTroubledCellIndicator_TwoMoment
  PUBLIC :: FinalizeTroubledCellIndicator_TwoMoment
  PUBLIC :: DetectTroubledCells_TwoMoment

  LOGICAL               :: UseTroubledCellIndicator
  INTEGER               :: nE, nE_G
  REAL(DP)              :: C_TCI
  REAL(DP), ALLOCATABLE :: WeightsX_X1_Up(:), WeightsX_X1_Dn(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X2_Up(:), WeightsX_X2_Dn(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X3_Up(:), WeightsX_X3_Dn(:)

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE &
  !$OMP TARGET( WeightsX_X1_Up, WeightsX_X1_Dn, &
  !$OMP         WeightsX_X2_Up, WeightsX_X2_Dn, &
  !$OMP         WeightsX_X3_Up, WeightsX_X3_Dn, &
  !$OMP         C_TCI )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE &
  !$ACC CREATE( WeightsX_X1_Up, WeightsX_X1_Dn, &
  !$ACC         WeightsX_X2_Up, WeightsX_X2_Dn, &
  !$ACC         WeightsX_X3_Up, WeightsX_X3_Dn, &
  !$ACC         C_TCI )
#endif

CONTAINS


  SUBROUTINE InitializeTroubledCellIndicator_TwoMoment &
    ( UseTroubledCellIndicator_Option, C_TCI_Option, Verbose_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UseTroubledCellIndicator_Option
    REAL(DP), INTENT(in), OPTIONAL :: C_TCI_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL  :: Verbose
    INTEGER  :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: jNodeX, jNodeX1, jNodeX2, jNodeX3
    REAL(DP) :: WeightX

    IF( PRESENT( UseTroubledCellIndicator_Option ) )THEN
      UseTroubledCellIndicator = UseTroubledCellIndicator_Option
    ELSE
      UseTroubledCellIndicator = .TRUE.
    END IF

    IF( PRESENT( C_TCI_Option ) )THEN
      C_TCI = C_TCI_Option
    ELSE
      C_TCI = 1.0d-2
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') '  INFO: InitializeTroubledCellIndicator_TwoMoment:'
      WRITE(*,'(A)') '  ------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A32,L1)') &
        '', 'Use Troubled Cell Indicator: ', UseTroubledCellIndicator
      WRITE(*,*)
      WRITE(*,'(A4,A32,ES11.3E3)') &
        '', 'C_TCI: ', C_TCI

    END IF

    ! --- Weights for Cell Average of Neighbors into Target Cell ---

    ALLOCATE( WeightsX_X1_Up(nDOFX), WeightsX_X1_Dn(nDOFX) )
    ALLOCATE( WeightsX_X2_Up(nDOFX), WeightsX_X2_Dn(nDOFX) )
    ALLOCATE( WeightsX_X3_Up(nDOFX), WeightsX_X3_Dn(nDOFX) )

    WeightsX_X1_Up = Zero; WeightsX_X1_Dn = Zero
    WeightsX_X2_Up = Zero; WeightsX_X2_Dn = Zero
    WeightsX_X3_Up = Zero; WeightsX_X3_Dn = Zero

    DO jNodeX = 1, nDOFX

      jNodeX1 = NodeNumberTableX(1,jNodeX)
      jNodeX2 = NodeNumberTableX(2,jNodeX)
      jNodeX3 = NodeNumberTableX(3,jNodeX)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        WeightX =   WeightsX1(iNodeX1) &
                  * WeightsX2(iNodeX2) &
                  * WeightsX3(iNodeX3)

        WeightsX_X1_Up(jNodeX) &
          = WeightsX_X1_Up(jNodeX) &
              + WeightX &
                * (   L_X1(jNodeX1) % P( NodesX1(iNodeX1) + One ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2)       ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3)       ) )

        WeightsX_X1_Dn(jNodeX) &
          = WeightsX_X1_Dn(jNodeX) &
              + WeightX &
                * (   L_X1(jNodeX1) % P( NodesX1(iNodeX1) - One ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2)       ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3)       ) )

        WeightsX_X2_Up(jNodeX) &
          = WeightsX_X2_Up(jNodeX) &
              + WeightX &
                * (   L_X1(jNodeX1) % P( NodesX1(iNodeX1)       ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) + One ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3)       ) )

        WeightsX_X2_Dn(jNodeX) &
          = WeightsX_X2_Dn(jNodeX) &
              + WeightX &
                * (   L_X1(jNodeX1) % P( NodesX1(iNodeX1)       ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2) - One ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3)       ) )

        WeightsX_X3_Up(jNodeX) &
          = WeightsX_X3_Up(jNodeX) &
              + WeightX &
                * (   L_X1(jNodeX1) % P( NodesX1(iNodeX1)       ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2)       ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) + One ) )

        WeightsX_X3_Dn(jNodeX) &
          = WeightsX_X3_Dn(jNodeX) &
              + WeightX &
                * (   L_X1(jNodeX1) % P( NodesX1(iNodeX1)       ) &
                    * L_X2(jNodeX2) % P( NodesX2(iNodeX2)       ) &
                    * L_X3(jNodeX3) % P( NodesX3(iNodeX3) - One ) )

      END DO

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( always, to: &
    !$OMP      WeightsX_X1_Up, WeightsX_X1_Dn, &
    !$OMP      WeightsX_X2_Up, WeightsX_X2_Dn, &
    !$OMP      WeightsX_X3_Up, WeightsX_X3_Dn, &
    !$OMP      C_TCI )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE &
    !$ACC DEVICE( WeightsX_X1_Up, WeightsX_X1_Dn, &
    !$ACC         WeightsX_X2_Up, WeightsX_X2_Dn, &
    !$ACC         WeightsX_X3_Up, WeightsX_X3_Dn, &
    !$ACC         C_TCI )
#endif

  END SUBROUTINE InitializeTroubledCellIndicator_TwoMoment


  SUBROUTINE FinalizeTroubledCellIndicator_TwoMoment

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: &
    !$OMP      WeightsX_X1_Up, WeightsX_X1_Dn, &
    !$OMP      WeightsX_X2_Up, WeightsX_X2_Dn, &
    !$OMP      WeightsX_X3_Up, WeightsX_X3_Dn, &
    !$OMP      C_TCI )
#endif

    DEALLOCATE( WeightsX_X1_Up, WeightsX_X1_Dn )
    DEALLOCATE( WeightsX_X2_Up, WeightsX_X2_Dn )
    DEALLOCATE( WeightsX_X3_Up, WeightsX_X3_Dn )

  END SUBROUTINE FinalizeTroubledCellIndicator_TwoMoment


  SUBROUTINE DetectTroubledCells_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, TroubledCell, SuppressBC_Option )

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR,1:nSpecies)
    LOGICAL, INTENT(out)    :: &
      TroubledCell &
        (iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    LOGICAL, INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    LOGICAL  :: SuppressBC
    INTEGER  :: iZ2, iZ3, iZ4, iE_G, iS
    REAL(DP) :: TCI
    REAL(DP) :: & ! --- Cell Averaged Density in Target Cell ---
      N_K0 (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: & ! --- Cell Averaged Densities from West   Cell Data ---
      N_KW (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies), &
      N_KW0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: & ! --- Cell Averaged Densities from East   Cell Data ---
      N_KE (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies), &
      N_KE0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: & ! --- Cell Averaged Densities from South  Cell Data ---
      N_KS (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies), &
      N_KS0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: & ! --- Cell Averaged Densities from North  Cell Data ---
      N_KN (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies), &
      N_KN0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: & ! --- Cell Averaged Densities from Bottom Cell Data ---
      N_KB (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies), &
      N_KB0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: & ! --- Cell Averaged Densities from Top    Cell Data ---
      N_KT (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies), &
      N_KT0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)

    IF( .NOT. UseTroubledCellIndicator )THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
      !$ACC PRESENT( iZ_B0, iZ_E0, TroubledCell )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iS   = 1       , nSpecies
      DO iE_G = 1       , nE_G
      DO iZ4  = iZ_B0(4), iZ_E0(4)
      DO iZ3  = iZ_B0(3), iZ_E0(3)
      DO iZ2  = iZ_B0(2), iZ_E0(2)

        TroubledCell(iZ2,iZ3,iZ4,iE_G,iS) = .TRUE.

      END DO
      END DO
      END DO
      END DO
      END DO

      RETURN

    END IF

    CALL TimersStart( Timer_TCI )

    IF( PRESENT( SuppressBC_Option ) )THEN
      SuppressBC = SuppressBC_Option
    ELSE
      SuppressBC = .FALSE.
    END IF

    nE   = iZ_E0(1) - iZ_B0(1) + 1
    nE_G = nE * nDOFE ! --- Global Energy Points

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B0, iZ_E0, U_R ) &
    !$OMP MAP( alloc: TroubledCell, N_K0, &
    !$OMP             N_KW, N_KW0, N_KE, N_KE0, &
    !$OMP             N_KS, N_KS0, N_KN, N_KN0, &
    !$OMP             N_KB, N_KB0, N_KT, N_KT0 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC COPYIN( iZ_B0, iZ_E0, U_R ) &
    !$ACC CREATE( TroubledCell, N_K0, &
    !$ACC         N_KW, N_KW0, N_KE, N_KE0, &
    !$ACC         N_KS, N_KS0, N_KN, N_KN0, &
    !$ACC         N_KB, N_KB0, N_KT, N_KT0 )
#endif

    IF( .NOT. SuppressBC )THEN

      CALL ApplyBoundaryConditions_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R )

    END IF

    CALL ComputeCellAverages_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, N_K0, N_KW, N_KW0, N_KE, N_KE0 )

    CALL ComputeCellAverages_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, N_K0, N_KS, N_KS0, N_KN, N_KN0 )

    CALL ComputeCellAverages_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, N_K0, N_KB, N_KB0, N_KT, N_KT0 )

    CALL TimersStart( Timer_TCI_Compute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( TCI )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRIVATE( TCI ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, TroubledCell, N_K0, &
    !$ACC          N_KW, N_KW0, N_KE, N_KE0, &
    !$ACC          N_KS, N_KS0, N_KN, N_KN0, &
    !$ACC          N_KB, N_KB0, N_KT, N_KT0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( TCI )
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      TCI &
        = (   ABS( N_K0(iZ2,iZ3,iZ4,iE_G,iS) - N_KW0(iZ2,iZ3,iZ4,iE_G,iS) ) &
            + ABS( N_K0(iZ2,iZ3,iZ4,iE_G,iS) - N_KE0(iZ2,iZ3,iZ4,iE_G,iS) ) &
            + ABS( N_K0(iZ2,iZ3,iZ4,iE_G,iS) - N_KS0(iZ2,iZ3,iZ4,iE_G,iS) ) &
            + ABS( N_K0(iZ2,iZ3,iZ4,iE_G,iS) - N_KN0(iZ2,iZ3,iZ4,iE_G,iS) ) &
            + ABS( N_K0(iZ2,iZ3,iZ4,iE_G,iS) - N_KB0(iZ2,iZ3,iZ4,iE_G,iS) ) &
            + ABS( N_K0(iZ2,iZ3,iZ4,iE_G,iS) - N_KT0(iZ2,iZ3,iZ4,iE_G,iS) ) ) &
          / MAX( ABS( N_K0(iZ2,iZ3,iZ4,iE_G,iS) ), &
                 ABS( N_KW(iZ2,iZ3,iZ4,iE_G,iS) ), &
                 ABS( N_KE(iZ2,iZ3,iZ4,iE_G,iS) ), &
                 ABS( N_KS(iZ2,iZ3,iZ4,iE_G,iS) ), &
                 ABS( N_KN(iZ2,iZ3,iZ4,iE_G,iS) ), &
                 ABS( N_KB(iZ2,iZ3,iZ4,iE_G,iS) ), &
                 ABS( N_KT(iZ2,iZ3,iZ4,iE_G,iS) ) )

      IF( TCI > C_TCI )THEN

        TroubledCell(iZ2,iZ3,iZ4,iE_G,iS) = .TRUE.

      ELSE

        TroubledCell(iZ2,iZ3,iZ4,iE_G,iS) = .FALSE.

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Compute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: iZ_B0, iZ_E0, U_R, TroubledCell, N_K0, &
    !$OMP               N_KW, N_KW0, N_KE, N_KE0, &
    !$OMP               N_KS, N_KS0, N_KN, N_KN0, &
    !$OMP               N_KB, N_KB0, N_KT, N_KT0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC DELETE( iZ_B0, iZ_E0, U_R, TroubledCell, N_K0, &
    !$ACC         N_KW, N_KW0, N_KE, N_KE0, &
    !$ACC         N_KS, N_KS0, N_KN, N_KN0, &
    !$ACC         N_KB, N_KB0, N_KT, N_KT0 )
#endif

    CALL TimersStop( Timer_TCI )

  END SUBROUTINE DetectTroubledCells_TwoMoment


  SUBROUTINE ComputeCellAverages_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, N_K0, N_KW, N_KW0, N_KE, N_KE0 )

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_K0 (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KW (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KW0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KE (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KE0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)

    INTEGER  :: &
      iE, iE_G, iZ2, iZ3, iZ4, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      N  (1:nDOFX, &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          1:nE_G,1:nSpecies, &
          iZ_B0(2)-1:iZ_E0(2)+1)
    REAL(DP) :: &
      N_K(iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          1:nE_G,1:nSpecies, &
          iZ_B0(2):iZ_E0(2))

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: N, N_K )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC CREATE( N, N_K )
#endif

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
    !$ACC PRIVATE( iNodeZ, iE_G ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, N, U_R )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#endif
    DO iZ2    = iZ_B0(2)-1, iZ_E0(2)+1
    DO iS     = 1         , nSpecies
    DO iE     = iZ_B0(1)  , iZ_E0(1)
    DO iNodeE = 1         , nDOFE
    DO iZ4    = iZ_B0(4)  , iZ_E0(4)
    DO iZ3    = iZ_B0(3)  , iZ_E0(3)

      DO iNodeX = 1, nDOFX

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE
        iE_G   = (iE    -1) * nDOFE + iNodeE

        N(iNodeX,iZ3,iZ4,iE_G,iS,iZ2) &
          = U_R(iNodeZ,iE,iZ2,iZ3,iZ4,iCR_N,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    nV_KX = PRODUCT( SHAPE( N_K ) )

    ! --- Cell Average in Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)  ), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_K0, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_K0(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ3,iZ4,iE_G,iS,iZ2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in West Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)-1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KW, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KW(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ3,iZ4,iE_G,iS,iZ2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average of West Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)-1), &
             nDOFX, WeightsX_X1_Up, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KW0, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KW0(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ3,iZ4,iE_G,iS,iZ2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in East Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)+1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KE, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KE(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ3,iZ4,iE_G,iS,iZ2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in East Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(3),iZ_B0(4),1,1,iZ_B0(2)+1), &
             nDOFX, WeightsX_X1_Dn, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KE0, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KE0(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ3,iZ4,iE_G,iS,iZ2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: N, N_K )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC DELETE( N, N_K )
#endif

  END SUBROUTINE ComputeCellAverages_X1


  SUBROUTINE ComputeCellAverages_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, N_K0, N_KS, N_KS0, N_KN, N_KN0 )

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR,1:nSpecies)
    REAL(DP), INTENT(in)  :: &
      N_K0 (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KS (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KS0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KN (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KN0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)

    INTEGER  :: &
      iE, iE_G, iZ2, iZ3, iZ4, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      N  (1:nDOFX, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(4):iZ_E0(4), &
          1:nE_G,1:nSpecies, &
          iZ_B0(3)-1:iZ_E0(3)+1)
    REAL(DP) :: &
      N_K(iZ_B0(2):iZ_E0(2), &
          iZ_B0(4):iZ_E0(4), &
          1:nE_G,1:nSpecies, &
          iZ_B0(3):iZ_E0(3))

    IF( iZ_E0(3) .EQ. iZ_B0(3) )THEN
#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
      !$ACC PRESENT( iZ_B0, iZ_E0, N_K0, N_KS, N_KS0, N_KN, N_KN0 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
       DO iS   = 1       , nSpecies
       DO iE_G = 1       , nE_G
       DO iZ4  = iZ_B0(4), iZ_E0(4)
       DO iZ3  = iZ_B0(3), iZ_E0(3)
       DO iZ2  = iZ_B0(2), iZ_E0(2)

         N_KS (iZ2,iZ3,iZ4,iE_G,iS) = N_K0(iZ2,iZ3,iZ4,iE_G,iS)
         N_KS0(iZ2,iZ3,iZ4,iE_G,iS) = N_K0(iZ2,iZ3,iZ4,iE_G,iS)
         N_KN (iZ2,iZ3,iZ4,iE_G,iS) = N_K0(iZ2,iZ3,iZ4,iE_G,iS)
         N_KN0(iZ2,iZ3,iZ4,iE_G,iS) = N_K0(iZ2,iZ3,iZ4,iE_G,iS)

       END DO
       END DO
       END DO
       END DO
       END DO
       RETURN
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: N, N_K )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC CREATE( N, N_K )
#endif

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
    !$ACC PRIVATE( iNodeZ, iE_G ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, N, U_R )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#endif
    DO iZ3    = iZ_B0(3)-1, iZ_E0(3)+1
    DO iS     = 1         , nSpecies
    DO iE     = iZ_B0(1)  , iZ_E0(1)
    DO iNodeE = 1         , nDOFE
    DO iZ4    = iZ_B0(4)  , iZ_E0(4)
    DO iZ2    = iZ_B0(2)  , iZ_E0(2)

      DO iNodeX = 1, nDOFX

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE
        iE_G   = (iE    -1) * nDOFE + iNodeE

        N(iNodeX,iZ2,iZ4,iE_G,iS,iZ3) &
          = U_R(iNodeZ,iE,iZ2,iZ3,iZ4,iCR_N,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    nV_KX = PRODUCT( SHAPE( N_K ) )

    ! --- Cell Average in South Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(2),iZ_B0(4),1,1,iZ_B0(3)-1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KS, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KS(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ2,iZ4,iE_G,iS,iZ3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average of South Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(3),iZ_B0(4),1,1,iZ_B0(3)-1), &
             nDOFX, WeightsX_X2_Up, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KS0, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KS0(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ2,iZ4,iE_G,iS,iZ3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in North Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(2),iZ_B0(4),1,1,iZ_B0(3)+1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KN, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KN(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ2,iZ4,iE_G,iS,iZ3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in North Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(2),iZ_B0(4),1,1,iZ_B0(3)+1), &
             nDOFX, WeightsX_X2_Dn, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KN0, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KN0(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ2,iZ4,iE_G,iS,iZ3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: N, N_K )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC DELETE( N, N_K )
#endif

  END SUBROUTINE ComputeCellAverages_X2


  SUBROUTINE ComputeCellAverages_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, N_K0, N_KB, N_KB0, N_KT, N_KT0 )

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR,1:nSpecies)
    REAL(DP), INTENT(in)  :: &
      N_K0 (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KB (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KB0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KT (iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KT0(iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4), &
            1:nE_G,1:nSpecies)

    INTEGER  :: &
      iE, iE_G, iZ2, iZ3, iZ4, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      N  (1:nDOFX, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          1:nE_G,1:nSpecies, &
          iZ_B0(4)-1:iZ_E0(4)+1)
    REAL(DP) :: &
      N_K(iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          1:nE_G,1:nSpecies, &
          iZ_B0(4):iZ_E0(4))

    IF( iZ_E0(4) .EQ. iZ_B0(4) )THEN
#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
      !$ACC PRESENT( iZ_B0, iZ_E0, N_K0, N_KB, N_KB0, N_KT, N_KT0 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
       DO iS   = 1       , nSpecies
       DO iE_G = 1       , nE_G
       DO iZ4  = iZ_B0(4), iZ_E0(4)
       DO iZ3  = iZ_B0(3), iZ_E0(3)
       DO iZ2  = iZ_B0(2), iZ_E0(2)

         N_KB (iZ2,iZ3,iZ4,iE_G,iS) = N_K0(iZ2,iZ3,iZ4,iE_G,iS)
         N_KB0(iZ2,iZ3,iZ4,iE_G,iS) = N_K0(iZ2,iZ3,iZ4,iE_G,iS)
         N_KT (iZ2,iZ3,iZ4,iE_G,iS) = N_K0(iZ2,iZ3,iZ4,iE_G,iS)
         N_KT0(iZ2,iZ3,iZ4,iE_G,iS) = N_K0(iZ2,iZ3,iZ4,iE_G,iS)

       END DO
       END DO
       END DO
       END DO
       END DO
      RETURN
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: N, N_K )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC CREATE( N, N_K )
#endif

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
    !$ACC PRIVATE( iNodeZ, iE_G ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, N, U_R )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ, iE_G )
#endif
    DO iZ4    = iZ_B0(4)-1, iZ_E0(4)+1
    DO iS     = 1         , nSpecies
    DO iE     = iZ_B0(1)  , iZ_E0(1)
    DO iNodeE = 1         , nDOFE
    DO iZ3    = iZ_B0(3)  , iZ_E0(3)
    DO iZ2    = iZ_B0(2)  , iZ_E0(2)

      DO iNodeX = 1, nDOFX

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE
        iE_G   = (iE    -1) * nDOFE + iNodeE

        N(iNodeX,iZ2,iZ3,iE_G,iS,iZ4) &
          = U_R(iNodeZ,iE,iZ2,iZ3,iZ4,iCR_N,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    nV_KX = PRODUCT( SHAPE( N_K ) )

    ! --- Cell Average in Bottom Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(2),iZ_B0(3),1,1,iZ_B0(4)-1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KB, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KB(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ2,iZ3,iE_G,iS,iZ4)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average of Bottom Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(2),iZ_B0(3),1,1,iZ_B0(4)-1), &
             nDOFX, WeightsX_X3_Up, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KB0, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KB0(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ2,iZ3,iE_G,iS,iZ4)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in Top Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(2),iZ_B0(3),1,1,iZ_B0(4)+1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KT, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KT(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ2,iZ3,iE_G,iS,iZ4)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in Top Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iZ_B0(2),iZ_B0(3),1,1,iZ_B0(4)+1), &
             nDOFX, WeightsX_X3_Dn, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_KT0, N_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iZ4  = iZ_B0(4), iZ_E0(4)
    DO iZ3  = iZ_B0(3), iZ_E0(3)
    DO iZ2  = iZ_B0(2), iZ_E0(2)

      N_KT0(iZ2,iZ3,iZ4,iE_G,iS) = N_K(iZ2,iZ3,iE_G,iS,iZ4)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: N, N_K )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC DELETE( N, N_K )
#endif

  END SUBROUTINE ComputeCellAverages_X3


END MODULE TwoMoment_TroubledCellIndicatorModule
