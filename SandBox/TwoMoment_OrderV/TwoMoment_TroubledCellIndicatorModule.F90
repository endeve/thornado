MODULE TwoMoment_TroubledCellIndicatorModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFE, nDOFX, nDimsX
  USE TwoMoment_TimersModule_OrderV, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_TCI
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
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, SuppressBC_Option )

    INTEGER,  INTENT(in)          :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout)       :: &
      U(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)
    LOGICAL, INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    LOGICAL  :: SuppressBC
    INTEGER  :: iX1, iX2, iX3, iZ1, iE, iCR, iS
    INTEGER  :: iNodeZ, iNodeE, iNodeX, nE_G, nV_KX
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: N_K0 (1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: N_KW (1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies) ! West
    REAL(DP) :: N_KW0(1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: N_KE (1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies) ! East
    REAL(DP) :: N_KE0(1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: N_KS (1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies) ! South
    REAL(DP) :: N_KS0(1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: N_KN (1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies) ! North
    REAL(DP) :: N_KN0(1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: N_KB (1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies) ! Bottom
    REAL(DP) :: N_KB0(1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: N_KT (1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies) ! Top
    REAL(DP) :: N_KT0(1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies)
    REAL(DP) :: CR_N(1:nDOFX, &
                     1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies, &
                     iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4))
    REAL(DP) :: TCI (1:(iZ_E0(1)-iZ_B0(1)+1)*nDOFE,1:nSpecies, &
                     iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4))
    REAL(DP) :: uDR (1:nDOFZ, &
                     iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                     iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
                     1,1:nSpecies)

    IF( .NOT. UseTroubledCellIndicator ) RETURN

    CALL TimersStart( Timer_TCI )

    IF( PRESENT( SuppressBC_Option ) )THEN
      SuppressBC = SuppressBC_Option
    ELSE
      SuppressBC = .FALSE.
    END IF

    IF( .NOT. SuppressBC )THEN

      CALL ApplyBoundaryConditions_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    END IF

    iX_B0 = iZ_B0(2:4)
    iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4)
    iX_E1 = iZ_E1(2:4)

    nE_G  = ( iZ_E0(1) - iZ_B0(1) + 1 ) * nDOFE ! --- Global Energy Points

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iS     = 1, nSpecies
      DO iZ1    = iZ_B0(1), iZ_E0(1)
      DO iNodeE = 1, nDOFE
      DO iNodeX = 1, nDOFX

        iE     = (iZ1   -1) * nDOFE + iNodeE
        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        CR_N(iNodeX,iE,iS,iX1,iX2,iX3) &
          = U(iNodeZ,iZ1,iX1,iX2,iX3,iCR_N,iS)

      END DO
      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

    nV_KX = nE_G * nSpecies ! --- Variables Per Spatial Element

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      CALL MatrixVectorMultiply &
             ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1  ,iX2,iX3), &
               nDOFX, WeightsX_q    , 1, Zero, N_K0 , 1 )

      CALL MatrixVectorMultiply &
             ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1-1,iX2,iX3), &
               nDOFX, WeightsX_q    , 1, Zero, N_KW , 1 )

      CALL MatrixVectorMultiply &
             ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1-1,iX2,iX3), &
               nDOFX, WeightsX_X1_Up, 1, Zero, N_KW0, 1 )

      CALL MatrixVectorMultiply &
             ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1+1,iX2,iX3), &
               nDOFX, WeightsX_q    , 1, Zero, N_KE , 1 )

      CALL MatrixVectorMultiply &
             ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1+1,iX2,iX3), &
               nDOFX, WeightsX_X1_Dn, 1, Zero, N_KE0, 1 )

      IF( nDimsX > 1 )THEN

        CALL MatrixVectorMultiply &
               ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1,iX2-1,iX3), &
                 nDOFX, WeightsX_q    , 1, Zero, N_KS , 1 )

        CALL MatrixVectorMultiply &
               ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1,iX2-1,iX3), &
                 nDOFX, WeightsX_X2_Up, 1, Zero, N_KS0, 1 )

        CALL MatrixVectorMultiply &
               ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1,iX2+1,iX3), &
                 nDOFX, WeightsX_q    , 1, Zero, N_KN , 1 )

        CALL MatrixVectorMultiply &
               ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1,iX2+1,iX3), &
                 nDOFX, WeightsX_X2_Dn, 1, Zero, N_KN0, 1 )

      ELSE

        N_KS = N_K0; N_KS0 = N_K0
        N_KN = N_K0; N_KN0 = N_K0

      END IF

      IF( nDimsX > 2 )THEN

        CALL MatrixVectorMultiply &
               ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1,iX2,iX3-1), &
                 nDOFX, WeightsX_q    , 1, Zero, N_KB , 1 )

        CALL MatrixVectorMultiply &
               ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1,iX2,iX3-1), &
                 nDOFX, WeightsX_X3_Up, 1, Zero, N_KB0, 1 )

        CALL MatrixVectorMultiply &
               ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1,iX2,iX3+1), &
                 nDOFX, WeightsX_q    , 1, Zero, N_KT , 1 )

        CALL MatrixVectorMultiply &
               ( 'T', nDOFX, nV_KX, One, CR_N(1,1,1,iX1,iX2,iX3+1), &
                 nDOFX, WeightsX_X3_Dn, 1, Zero, N_KT0, 1 )

      ELSE

        N_KB = N_K0; N_KB0 = N_K0
        N_KT = N_K0; N_KT0 = N_K0

      END IF

      DO iS = 1, nSpecies
      DO iE = 1, nE_G

        TCI(iE,iS,iX1,iX2,iX3) &
          = (   ABS( N_K0(iE,iS) - N_KW0(iE,iS) ) &
              + ABS( N_K0(iE,iS) - N_KE0(iE,iS) ) &
              + ABS( N_K0(iE,iS) - N_KS0(iE,iS) ) &
              + ABS( N_K0(iE,iS) - N_KN0(iE,iS) ) &
              + ABS( N_K0(iE,iS) - N_KB0(iE,iS) ) &
              + ABS( N_K0(iE,iS) - N_KT0(iE,iS) ) ) &
            / MAX( ABS( N_K0(iE,iS) ), &
                   ABS( N_KW(iE,iS) ), &
                   ABS( N_KE(iE,iS) ), &
                   ABS( N_KS(iE,iS) ), &
                   ABS( N_KN(iE,iS) ), &
                   ABS( N_KB(iE,iS) ), &
                   ABS( N_KT(iE,iS) ) )

      END DO
      END DO

    END DO
    END DO
    END DO

    DO iS  = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE
        iE     = (iZ1   -1) * nDOFE + iNodeE

        uDR(iNodeZ,iZ1,iX1,iX2,iX3,1,iS) = TCI(iE,iS,iX1,iX2,iX3)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TCI )

  END SUBROUTINE DetectTroubledCells_TwoMoment


END MODULE TwoMoment_TroubledCellIndicatorModule
