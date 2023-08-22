MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, One, TwoPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFZ, nDOFX, nDOFE, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    xL, xR
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModuleZ, ONLY: &
    NodeNumberTableZ
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nSpecies, &
    uPM, iPM_J, iPM_H1, iPM_H2, iPM_H3, nPM, &
    uCM, iCM_E, iCM_F1, iCM_F2, iCM_F3, nCM
  USE TwoMoment_UtilitiesModule_FMC, ONLY: &
    ComputeConserved_TwoMoment_FMC

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming( V_0 )

    END SELECT

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( uCF, uPF, uCM, uPM )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( uCF, uPF, uCM, uPM )
#endif

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_SineWaveStreaming( V_0 )

    REAL(DP), INTENT(in) :: V_0(3)

    INTEGER :: iNodeX, iX1, iX2, iX3
    INTEGER :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER :: iNodeZ2, iNodeZ3, iNodeZ4
    Real(DP) :: X1, X2, X3, W

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_0(1)
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_0(2)
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_0(3)
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

      END DO

    END DO
    END DO
    END DO

    ! --- Two-Moment Fields ---

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE,nDOFX ) + 1

        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)
        iNodeZ3 = NodeNumberTableZ(3,iNodeZ)
        iNodeZ4 = NodeNumberTableZ(4,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeZ2)
        X2 = NodeCoordinate( MeshX(2), iZ3, iNodeZ3)
        X3 = NodeCoordinate( MeshX(3), iZ4, iNodeZ4)

        W = One / SQRT( One - uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V1)**2 )

        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS) &
          = 0.50_DP + 0.49_DP * SIN( TwoPi * X1)
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS) &
          = W * uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J,iS)
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS) &
          = 0.0_DP
        uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS) &
          = 0.0_DP

        CALL ComputeConserved_TwoMoment_FMC &
               ( uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 uPM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 uCM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V1), &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V2), &
                 uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V3), &
                 uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11),  &
                 uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22),  &
                 uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_SineWaveStreaming


END MODULE InitializationModule
