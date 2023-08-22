MODULE TwoMoment_DiscretizationModule_Collisions_FMC

  USE KindModule, ONLY: &
    DP, Zero, One, Half, Three
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nSpecies, &
    nCM, iCM_E, iCM_F1, iCM_F2, iCM_F3
  USE TwoMoment_UtilitiesModule_FMC, ONLY: &
    EddingtonTensorComponents_dd
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor, &
    HeatFluxFactor
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit

  CONTAINS

  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_M, dU_M)

    ! --- Input/Output variables ---
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: dt
    REAL(DP), INTENT(in) :: GE(1:nDOFE, iZ_B1(1):iZ_E1(1), 1:nGE)
    REAL(DP), INTENT(in) :: GX(1:nDOFX, &
                               iZ_B1(2):iZ_E1(2), &
                               iZ_B1(3):iZ_E1(3), &
                               iZ_B1(4):iZ_E1(4), &
                               1:nGF)
    REAL(DP), INTENT(in) :: U_M(1:nDOFZ, &
                                iZ_B1(1):iZ_E1(1), &
                                iZ_B1(2):iZ_E1(2), &
                                iZ_B1(3):iZ_E1(3), &
                                iZ_B1(4):iZ_E1(4), &
                                1:nCM, &
                                1:nSpecies)
    REAL(DP), INTENT(out) :: dU_M(1:nDOFZ, &
                                  iZ_B1(1):iZ_E1(1), &
                                  iZ_B1(2):iZ_E1(2), &
                                  iZ_B1(3):iZ_E1(3), &
                                  iZ_B1(4):iZ_E1(4), &
                                  1:nCM, &
                                  1:nSpecies)
    
    Write(*,*)
    print *,'ComputeIncrement_TwoMoment_Implicit'

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit

END MODULE TwoMoment_DiscretizationModule_Collisions_FMC