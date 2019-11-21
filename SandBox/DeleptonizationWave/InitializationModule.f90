MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, TwoPi
  USE UnitsModule, ONLY: &
    Gram, Centimeter, &
    Kilometer, Kelvin, &
    MeV, BoltzmannConstant
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    iE_B0, iE_E0, &
    nDOF, nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModule, ONLY: &
    NodeNumbersX, &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshE, MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, &
    iAF_Xa, iAF_Xh, iAF_Gm
  USE RadiationFieldsModule, ONLY: &
    nSpecies, iNuE, iNuE_Bar, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule_TABLE, ONLY: &
    ApplyEquationOfState_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializeFields_DeleptonizationWave

CONTAINS


  SUBROUTINE InitializeFields_DeleptonizationWave &
    ( ProfileName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ProfileName_Option

    IF( PRESENT( ProfileName_Option ) )THEN

      IF( LEN_TRIM( ProfileName_Option ) .GT. 0 )THEN

        CALL InitializeFluidFields_DeleptonizationWave_Profile &
               ( ProfileName_Option )

      ELSE

        CALL InitializeFluidFields_DeleptonizationWave

      END IF

    ELSE

      CALL InitializeFluidFields_DeleptonizationWave

    END IF

    CALL InitializeRadiationFields_DeleptonizationWave

  END SUBROUTINE InitializeFields_DeleptonizationWave


  SUBROUTINE InitializeFluidFields_DeleptonizationWave

    ! --- Density Profile ---
    REAL(DP), PARAMETER :: MinD = 1.0d08 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: MaxD = 4.0d14 * Gram / Centimeter**3
!    REAL(DP), PARAMETER :: C_D  = 7.5_DP
    REAL(DP), PARAMETER :: R_D  = 2.0d01 * Kilometer
    REAL(DP), PARAMETER :: H_D  = 1.0d01 * Kilometer
!    REAL(DP), PARAMETER :: H_D  = 5.0d00 * Kilometer
    ! --- Temperature Profile ---
    REAL(DP), PARAMETER :: MinT = 5.0d09 * Kelvin
    REAL(DP), PARAMETER :: MaxT = 2.6d11 * Kelvin
!    REAL(DP), PARAMETER :: C_T  = 1.0_DP
    REAL(DP), PARAMETER :: R_T  = 2.5d01 * Kilometer
    REAL(DP), PARAMETER :: H_T  = 2.0d01 * Kilometer
!    REAL(DP), PARAMETER :: H_T  = 5.0d01 * Kilometer
    ! --- Electron Fraction Profile ---
    REAL(DP), PARAMETER :: MinY = 3.0d-1
    REAL(DP), PARAMETER :: MaxY = 4.6d-1
!    REAL(DP), PARAMETER :: C_Y  = 1.0_DP
    REAL(DP), PARAMETER :: R_Y  = 4.5d01 * Kilometer
    REAL(DP), PARAMETER :: H_Y  = 1.0d01 * Kilometer
!    REAL(DP), PARAMETER :: H_Y  = 5.0d01 * Kilometer

    INTEGER  :: iX1, iX2, iX3, iNodeX
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3, R

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = &!MeshX(1) % Center(iX1)
             NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = &!MeshX(2) % Center(iX2)
             NodeCoordinate( MeshX(2), iX2, iNodeX2 )
        X3 = &!MeshX(3) % Center(iX3)
             NodeCoordinate( MeshX(3), iX3, iNodeX3 )

        IF( TRIM( CoordinateSystem ) == 'SPHERICAL' )THEN

          R = X1

        ELSE

          R = SQRT( X1**2 + X2**2 + X3**2 )

        END IF

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = Half * ( MaxD * ( One - TANH( (R-R_D)/H_D ) ) &
                     + MinD * ( One - TANH( (R_D-R)/H_D ) ) )
!!$              = MaxD * C_D / ( C_D + ( R / H_D )**4 )

        uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
          = Half * ( MaxT * ( One - TANH( (R-R_T)/H_T ) ) &
                     + MinT * ( One - TANH( (R_T-R)/H_T ) ) )
!!$              = MaxT * C_T / ( C_T + ( R / H_T )**2 )

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
          = Half * ( MinY * ( One - TANH( (R-R_Y)/H_Y ) ) &
                     + MaxY * ( One - TANH( (R_Y-R)/H_Y ) ) )
!!$              = MinY * ( One + C_Y / ( C_Y + ( R / H_Y )**(-12) ) )

      END DO

      CALL ComputeThermodynamicStates_Primitive_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne) )

      uPF(:,iX1,iX2,iX3,iPF_V1) = Zero
      uPF(:,iX1,iX2,iX3,iPF_V2) = Zero
      uPF(:,iX1,iX2,iX3,iPF_V3) = Zero

      CALL ApplyEquationOfState_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
               uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
               uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
               uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
               uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )

      ! --- Conserved ---

      uCF(:,iX1,iX2,iX3,iCF_D ) = uPF(:,iX1,iX2,iX3,iPF_D)
      uCF(:,iX1,iX2,iX3,iCF_S1) = Zero
      uCF(:,iX1,iX2,iX3,iCF_S2) = Zero
      uCF(:,iX1,iX2,iX3,iCF_S3) = Zero
      uCF(:,iX1,iX2,iX3,iCF_E ) = uPF(:,iX1,iX2,iX3,iPF_E)
      uCF(:,iX1,iX2,iX3,iCF_Ne) = uPF(:,iX1,iX2,iX3,iPF_Ne)

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFluidFields_DeleptonizationWave


  SUBROUTINE InitializeFluidFields_DeleptonizationWave_Profile &
    ( ProfileName )

    CHARACTER(LEN=*), INTENT(in) :: ProfileName

    INTEGER  :: iX1, iX2, iX3, iNodeX, iR
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3, R
    REAL(DP), ALLOCATABLE :: R_P(:), D_P(:), T_P(:), Y_P(:)

    WRITE(*,*)
    WRITE(*,'(A6,A,A)') '', 'Initializing from Profile: ', TRIM( ProfileName )

    CALL ReadFluidProfile( ProfileName, R_P, D_P, T_P, Y_P )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

        IF( TRIM( CoordinateSystem ) == 'SPHERICAL' )THEN

          R = X1

        ELSE

          R = SQRT( X1**2 + X2**2 + X3**2 )

        END IF

        iR = MAX( Locate( R, R_P, SIZE( R_P ) ), 1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = Interpolate1D_Linear( R, R_P(iR), R_P(iR+1), D_P(iR), D_P(iR+1) )

        uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
          = Interpolate1D_Linear( R, R_P(iR), R_P(iR+1), T_P(iR), T_P(iR+1) )

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
          = Interpolate1D_Linear( R, R_P(iR), R_P(iR+1), Y_P(iR), Y_P(iR+1) )

      END DO

      CALL ComputeThermodynamicStates_Primitive_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne) )

      uPF(:,iX1,iX2,iX3,iPF_V1) = Zero
      uPF(:,iX1,iX2,iX3,iPF_V2) = Zero
      uPF(:,iX1,iX2,iX3,iPF_V3) = Zero

      CALL ApplyEquationOfState_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
               uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
               uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
               uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
               uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )

      ! --- Conserved ---

      uCF(:,iX1,iX2,iX3,iCF_D ) = uPF(:,iX1,iX2,iX3,iPF_D)
      uCF(:,iX1,iX2,iX3,iCF_S1) = Zero
      uCF(:,iX1,iX2,iX3,iCF_S2) = Zero
      uCF(:,iX1,iX2,iX3,iCF_S3) = Zero
      uCF(:,iX1,iX2,iX3,iCF_E ) = uPF(:,iX1,iX2,iX3,iPF_E)
      uCF(:,iX1,iX2,iX3,iCF_Ne) = uPF(:,iX1,iX2,iX3,iPF_Ne)

    END DO
    END DO
    END DO

    DEALLOCATE( R_P, D_P, T_P, Y_P )

  END SUBROUTINE InitializeFluidFields_DeleptonizationWave_Profile


  SUBROUTINE InitializeRadiationFields_DeleptonizationWave

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNode, iNodeE
    REAL(DP) :: kT(nDOF)
    REAL(DP) :: Mnu(nDOF), E
    REAL(DP) :: Gm_dd_11(nDOF)
    REAL(DP) :: Gm_dd_22(nDOF)
    REAL(DP) :: Gm_dd_33(nDOF)

    DO iS = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Gm_dd_11 &
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_11)

      Gm_dd_22 &
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_22)

      Gm_dd_33 &
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_33)

      kT = BoltzmannConstant &
             * uAF(NodeNumbersX,iX1,iX2,iX3,iAF_T)

      IF( iS .EQ. iNuE )THEN

        Mnu = + ( uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Me) &
                  + uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mp) &
                  - uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mn) )

      ELSEIF( iS .EQ. iNuE_Bar )THEN

        Mnu = - ( uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Me) &
                  + uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mp) &
                  - uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mn) )

      END IF

      DO iE = iE_B0, iE_E0
      DO iNode = 1, nDOF

        iNodeE = NodeNumberTable(1,iNode)

        E = NodeCoordinate( MeshE, iE, iNodeE )

        uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
          = MAX( One / ( EXP( (E-Mnu(iNode))/kT(iNode) ) + One ), 1.0d-99 )

        uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
          = Zero

        uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
          = Zero

        uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
          = Zero
        
      END DO

        CALL ComputeConserved_TwoMoment &
               ( uPR(:,iE,iX1,iX2,iX3,iPR_D, iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_N, iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      END DO

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeRadiationFields_DeleptonizationWave


  SUBROUTINE ReadFluidProfile( FileName, R, D, T, Y )

    CHARACTER(*),          INTENT(in)    :: FileName
    REAL(DP), ALLOCATABLE, INTENT(inout) :: R(:), D(:), T(:), Y(:)

    CHARACTER(LEN=9)      :: Format1 = '(4ES12.3)'
    INTEGER               :: nPoints, iPoint
    INTEGER               :: Status
    REAL(DP)              :: Buffer(4)
    REAL(DP), ALLOCATABLE :: Data(:,:)

    ! --- Count Lines ---

    nPoints = 0
    OPEN( 1, FILE = TRIM( FileName ), FORM = "formatted", ACTION = 'read' )
    READ( 1, * )
    DO
      READ( 1, Format1, IOSTAT=Status ) Buffer
      IF( Status .NE. 0 ) EXIT
      nPoints = nPoints + 1
    END DO
    CLOSE( 1, STATUS = 'keep' )

    ! --- Read Data ---

    ALLOCATE( Data(nPoints,4) )

    OPEN( 1, FILE = TRIM( FileName ), FORM = "formatted", ACTION = 'read' )
    READ( 1, * )
    DO iPoint = 1, nPoints
      READ( 1, Format1, IOSTAT=Status ) Data(iPoint,:)
    END DO
    CLOSE( 1, STATUS = 'keep' )

    ALLOCATE( R(nPoints), D(nPoints), T(nPoints), Y(nPoints) )

    R = Data(:,1) * Kilometer
    D = Data(:,2) * Gram / Centimeter**3
    T = Data(:,3) * MeV
    Y = Data(:,4)

    DEALLOCATE( Data )

  END SUBROUTINE ReadFluidProfile


END MODULE InitializationModule
