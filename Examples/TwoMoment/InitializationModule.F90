MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, One, Two, TwoPi
  USE UnitsModule, ONLY: &
    Centimeter, &
    Kilometer, &
    Gram, &
    MeV, &
    BoltzmannConstant, &
    SpeedOfLight
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    iE_B0, iE_E0, &
    nDOF, nDOFE, nDOFX
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, iAF_Me, iAF_Mp, iAF_Mn, &
    iAF_Xp, iAF_Xn, iAF_Xa, iAF_Xh, iAF_Gm
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment
  USE EquationOfStateModule_TABLE, ONLY: &
    ApplyEquationOfState_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields
  PUBLIC :: ComputeError

CONTAINS


  SUBROUTINE InitializeFields

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming

      CASE( 'Relaxation' )

        CALL InitializeFields_Relaxation

    END SELECT

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( uCF, uCR )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( uCF, uCR )
#endif

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_SineWaveStreaming

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNode, iNodeX1
    REAL(DP) :: X1, Ones(nDOF), L

    Ones = One
    L    = 1.0d2 * Kilometer

    DO iS  = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNode = 1, nDOF

        iNodeX1 = NodeNumberTable(2,iNode)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
          = 0.5_DP + 0.49_DP * SIN( TwoPi * X1 / L )

        uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
          = ( One - 1.0d-12 ) * uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS)

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
               Ones, Ones, Ones )

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_SineWaveStreaming


  SUBROUTINE InitializeFields_Relaxation

    REAL(DP), PARAMETER :: D_0 = 2.564d13 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: T_0 = 1.978d01 * MeV
    REAL(DP), PARAMETER :: Y_0 = 0.2876d0

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeE, iNodeX, iNodeZ
    REAL(DP) :: kT, E, Ones(nDOF)

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = D_0

        uAF(iNodeX,iX1,iX2,iX3,iAF_T ) = T_0

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = Y_0

        CALL ComputeThermodynamicStates_Primitive_TABLE &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

        CALL ApplyEquationOfState_TABLE &
             ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_S ), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_Me), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_Mp), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_Mn), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_Xp), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_Xn), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_Xa), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_Xh), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_Gm) )

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 One, One, One )

      END DO

    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

    Ones = One

    DO iS  = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

        kT = BoltzmannConstant * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

        E = NodeCoordinate( MeshE, iE, iNodeE )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D,iS) &
          = MAX( 0.99_DP * EXP( - ( E - Two * kT )**2 &
                                  / ( Two * (1.0d1*MeV)**2 ) ), 1.0d-99 )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) &
          = Zero

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) &
          = Zero

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) &
          = Zero

      END DO
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
               Ones, Ones, Ones )

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Relaxation


  SUBROUTINE ComputeError( t )

    REAL(DP), INTENT(in) :: t

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'SineWaveStreaming' )

        CALL ComputeError_SineWaveStreaming( t )

      CASE( 'Relaxation' )

        CALL ComputeError_Relaxation

    END SELECT

  END SUBROUTINE ComputeError


  SUBROUTINE ComputeError_SineWaveStreaming( t )

    REAL(DP), INTENT(in) :: t

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNode, iNodeX1
    REAL(DP) :: MaxError(nSpecies), L, X1, N_A

    L = 1.0d2 * Kilometer

    MaxError = Zero

    DO iS  = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNode = 1, nDOF

        iNodeX1 = NodeNumberTable(2,iNode)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        N_A = 0.5_DP + 0.49_DP * SIN( TwoPi * ( X1 - SpeedOfLight * t ) / L )

        MaxError(iS) &
          = MAX( ABS( N_A - uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) ), MaxError(iS) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    WRITE(*,*)
    WRITE(*,'(A2,A)') '', 'INFO: Error Check'
    WRITE(*,*)
    DO iS = 1, nSpecies
      WRITE(*,'(A4,A10,I2.2,A14,ES10.4E2)') &
      '', 'Species = ', iS, ', Inf Error = ', MaxError(iS)
    END DO
    WRITE(*,*)

  END SUBROUTINE ComputeError_SineWaveStreaming


  SUBROUTINE ComputeError_Relaxation

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeE, iNodeX, iNodeZ
    REAL(DP) :: MaxError(nSpecies), N0, kT, Mnu, E

    MaxError = Zero

    DO iS = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

        kT = BoltzmannConstant * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

        Mnu = ( - One )**(iS-1) &
              * (   uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                  + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                  - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn) )

        E = NodeCoordinate( MeshE, iE, iNodeE )

        N0 = One / ( EXP( ( E - Mnu ) / kT ) + One )

        MaxError(iS) &
          = MAX( ABS(N0-uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS))/N0, MaxError(iS) )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    WRITE(*,*)
    WRITE(*,'(A2,A)') '', 'INFO: Error Check'
    WRITE(*,*)
    DO iS = 1, nSpecies
      WRITE(*,'(A4,A10,I2.2,A14,ES10.4E2)') &
      '', 'Species = ', iS, ', Inf Error = ', MaxError(iS)
    END DO
    WRITE(*,*)

  END SUBROUTINE ComputeError_Relaxation


END MODULE InitializationModule
