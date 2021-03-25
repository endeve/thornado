MODULE TwoMoment_TroubledCellIndicatorModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFE, nDOFX, nDimsX
  USE TwoMoment_TimersModule_OrderV, ONLY: &
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
  INTEGER               :: iX_B0(3), iX_E0(3)
  INTEGER               :: iE_B0, iE_E0, nE_G
  REAL(DP)              :: C_TCI
  REAL(DP), ALLOCATABLE :: WeightsX_X1_Up(:), WeightsX_X1_Dn(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X2_Up(:), WeightsX_X2_Dn(:)
  REAL(DP), ALLOCATABLE :: WeightsX_X3_Up(:), WeightsX_X3_Dn(:)

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

  END SUBROUTINE InitializeTroubledCellIndicator_TwoMoment


  SUBROUTINE FinalizeTroubledCellIndicator_TwoMoment

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
    INTEGER  :: iX1, iX2, iX3, iE_G, iS
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
      TroubledCell = .FALSE.
      RETURN
    END IF

    CALL TimersStart( Timer_TCI )

    IF( PRESENT( SuppressBC_Option ) )THEN
      SuppressBC = SuppressBC_Option
    ELSE
      SuppressBC = .FALSE.
    END IF

    IF( .NOT. SuppressBC )THEN

      CALL ApplyBoundaryConditions_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R )

    END IF

    iE_B0 = iZ_B0(1); iX_B0 = iZ_B0(2:4)
    iE_E0 = iZ_B0(1); iX_E0 = iZ_E0(2:4)

    nE_G  = ( iE_E0 - iE_B0 + 1 ) * nDOFE ! --- Global Energy Points

    CALL ComputeCellAverages_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, N_K0, N_KW, N_KW0, N_KE, N_KE0 )

    CALL ComputeCellAverages_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R      , N_KS, N_KS0, N_KN, N_KN0 )

    CALL ComputeCellAverages_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R      , N_KB, N_KB0, N_KT, N_KT0 )

    CALL TimersStart( Timer_TCI_Compute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      TCI &
        = (   ABS( N_K0(iX1,iX2,iX3,iE_G,iS) - N_KW0(iX1,iX2,iX3,iE_G,iS) ) &
            + ABS( N_K0(iX1,iX2,iX3,iE_G,iS) - N_KE0(iX1,iX2,iX3,iE_G,iS) ) &
            + ABS( N_K0(iX1,iX2,iX3,iE_G,iS) - N_KS0(iX1,iX2,iX3,iE_G,iS) ) &
            + ABS( N_K0(iX1,iX2,iX3,iE_G,iS) - N_KN0(iX1,iX2,iX3,iE_G,iS) ) &
            + ABS( N_K0(iX1,iX2,iX3,iE_G,iS) - N_KB0(iX1,iX2,iX3,iE_G,iS) ) &
            + ABS( N_K0(iX1,iX2,iX3,iE_G,iS) - N_KT0(iX1,iX2,iX3,iE_G,iS) ) ) &
          / MAX( ABS( N_K0(iX1,iX2,iX3,iE_G,iS) ), &
                 ABS( N_KW(iX1,iX2,iX3,iE_G,iS) ), &
                 ABS( N_KE(iX1,iX2,iX3,iE_G,iS) ), &
                 ABS( N_KS(iX1,iX2,iX3,iE_G,iS) ), &
                 ABS( N_KN(iX1,iX2,iX3,iE_G,iS) ), &
                 ABS( N_KB(iX1,iX2,iX3,iE_G,iS) ), &
                 ABS( N_KT(iX1,iX2,iX3,iE_G,iS) ) )

      IF( TCI > C_TCI )THEN

        TroubledCell(iX1,iX2,iX3,iE_G,iS) = .TRUE.

      ELSE

        TroubledCell(iX1,iX2,iX3,iE_G,iS) = .FALSE.

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Compute )

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
      N_K0 (iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KW (iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KW0(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KE (iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KE0(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)

    INTEGER  :: &
      iE, iE_G, iX1, iX2, iX3, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      N  (1:nDOFX, &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3), &
          1:nE_G,1:nSpecies, &
          iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: &
      N_K(iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3), &
          1:nE_G,1:nSpecies, &
          iX_B0(1):iX_E0(1))

    CALL TimersStart( Timer_TCI_Permute )

    DO iX1    = iX_B0(1)-1, iX_E0(1)+1
    DO iS     = 1         , nSpecies
    DO iE     = 1         , nE_G
    DO iNodeE = iE_B0     , iE_E0
    DO iX3    = iX_B0(3)  , iX_E0(3)
    DO iX2    = iX_B0(2)  , iX_E0(2)

      DO iNodeX = 1, nDOFX

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE
        iE_G   = (iE    -1) * nDOFE + iNodeE

        N(iNodeX,iX2,iX3,iE_G,iS,iX1) &
          = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)

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
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(2),iX_B0(3),1,1,iX_B0(1)  ), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_K0(iX1,iX2,iX3,iE_G,iS) = N_K(iX2,iX3,iE_G,iS,iX1)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in West Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(2),iX_B0(3),1,1,iX_B0(1)-1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KW(iX1,iX2,iX3,iE_G,iS) = N_K(iX2,iX3,iE_G,iS,iX1)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average of West Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(2),iX_B0(3),1,1,iX_B0(1)-1), &
             nDOFX, WeightsX_X1_Up, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KW0(iX1,iX2,iX3,iE_G,iS) = N_K(iX2,iX3,iE_G,iS,iX1)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in East Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(2),iX_B0(3),1,1,iX_B0(1)+1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KE(iX1,iX2,iX3,iE_G,iS) = N_K(iX2,iX3,iE_G,iS,iX1)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in East Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(2),iX_B0(3),1,1,iX_B0(1)+1), &
             nDOFX, WeightsX_X1_Dn, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KE0(iX1,iX2,iX3,iE_G,iS) = N_K(iX2,iX3,iE_G,iS,iX1)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

  END SUBROUTINE ComputeCellAverages_X1


  SUBROUTINE ComputeCellAverages_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, N_KS, N_KS0, N_KN, N_KN0 )

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
      N_KS (iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KS0(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KN (iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KN0(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)

    INTEGER  :: &
      iE, iE_G, iX1, iX2, iX3, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      N  (1:nDOFX, &
          iX_B0(1):iX_E0(1), &
          iX_B0(3):iX_E0(3), &
          1:nE_G,1:nSpecies, &
          iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: &
      N_K(iX_B0(1):iX_E0(1), &
          iX_B0(3):iX_E0(3), &
          1:nE_G,1:nSpecies, &
          iX_B0(2):iX_E0(2))

    CALL TimersStart( Timer_TCI_Permute )

    DO iX2    = iX_B0(2)-1, iX_E0(2)+1
    DO iS     = 1         , nSpecies
    DO iE     = 1         , nE_G
    DO iNodeE = iE_B0     , iE_E0
    DO iX3    = iX_B0(3)  , iX_E0(3)
    DO iX1    = iX_B0(1)  , iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE
        iE_G   = (iE    -1) * nDOFE + iNodeE

        N(iNodeX,iX1,iX3,iE_G,iS,iX2) &
          = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)

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
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(1),iX_B0(3),1,1,iX_B0(2)-1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KS(iX1,iX2,iX3,iE_G,iS) = N_K(iX1,iX3,iE_G,iS,iX2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average of South Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(2),iX_B0(3),1,1,iX_B0(2)-1), &
             nDOFX, WeightsX_X2_Up, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KS0(iX1,iX2,iX3,iE_G,iS) = N_K(iX1,iX3,iE_G,iS,iX2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in North Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(1),iX_B0(3),1,1,iX_B0(2)+1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KN(iX1,iX2,iX3,iE_G,iS) = N_K(iX1,iX3,iE_G,iS,iX2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in North Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(1),iX_B0(3),1,1,iX_B0(2)+1), &
             nDOFX, WeightsX_X2_Dn, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KN0(iX1,iX2,iX3,iE_G,iS) = N_K(iX1,iX3,iE_G,iS,iX2)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

  END SUBROUTINE ComputeCellAverages_X2


  SUBROUTINE ComputeCellAverages_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, N_KB, N_KB0, N_KT, N_KT0 )

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
      N_KB (iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KB0(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KT (iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      N_KT0(iX_B0(1):iX_E0(1), &
            iX_B0(2):iX_E0(2), &
            iX_B0(3):iX_E0(3), &
            1:nE_G,1:nSpecies)

    INTEGER  :: &
      iE, iE_G, iX1, iX2, iX3, iS, &
      iNodeE, iNodeX, iNodeZ, nV_KX
    REAL(DP) :: &
      N  (1:nDOFX, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          1:nE_G,1:nSpecies, &
          iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: &
      N_K(iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          1:nE_G,1:nSpecies, &
          iX_B0(3):iX_E0(3))

    CALL TimersStart( Timer_TCI_Permute )

    DO iX3    = iX_B0(3)-1, iX_E0(3)+1
    DO iS     = 1         , nSpecies
    DO iE     = 1         , nE_G
    DO iNodeE = iE_B0     , iE_E0
    DO iX2    = iX_B0(2)  , iX_E0(2)
    DO iX1    = iX_B0(1)  , iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE
        iE_G   = (iE    -1) * nDOFE + iNodeE

        N(iNodeX,iX1,iX2,iE_G,iS,iX3) &
          = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)

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
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(1),iX_B0(2),1,1,iX_B0(3)-1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KB(iX1,iX2,iX3,iE_G,iS) = N_K(iX1,iX2,iE_G,iS,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average of Bottom Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(1),iX_B0(2),1,1,iX_B0(3)-1), &
             nDOFX, WeightsX_X3_Up, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KB0(iX1,iX2,iX3,iE_G,iS) = N_K(iX1,iX2,iE_G,iS,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in Top Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(1),iX_B0(2),1,1,iX_B0(3)+1), &
             nDOFX, WeightsX_q, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KT(iX1,iX2,iX3,iE_G,iS) = N_K(iX1,iX2,iE_G,iS,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

    ! --- Cell Average in Top Cell Extrapolated into Target Cell ---

    CALL TimersStart( Timer_TCI_LinearAlgebra )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFx, nV_KX, One, N(1,iX_B0(1),iX_B0(2),1,1,iX_B0(3)+1), &
             nDOFX, WeightsX_X3_Dn, 1, Zero, N_K, 1 )

    CALL TimersStop( Timer_TCI_LinearAlgebra )

    CALL TimersStart( Timer_TCI_Permute )

    DO iS   = 1       , nSpecies
    DO iE_G = 1       , nE_G
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      N_KT0(iX1,iX2,iX3,iE_G,iS) = N_K(iX1,iX2,iE_G,iS,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI_Permute )

  END SUBROUTINE ComputeCellAverages_X3


END MODULE TwoMoment_TroubledCellIndicatorModule
