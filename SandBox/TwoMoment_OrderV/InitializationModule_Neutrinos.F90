MODULE InitializationModule_Neutrinos

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, Three, TwoPi, SqrtTiny
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer, &
    Millisecond, &
    MeV, &
    Kelvin, &
    Erg, &
    SpeedOfLight, &
    BoltzmannConstant
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    iE_B0, iE_E0, &
    nDOFX, nDOFE
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    WeightsX1
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshE, MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, &
    iAF_Xa, iAF_Xh, iAF_Gm
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, iNuE, iNuE_Bar, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeTemperatureFromPressure_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ApplyEquationOfState_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions_DG, &
    ComputeNeutrinoOpacities_EC, &
    SetOpacityMask, &
    nOp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields
  PUBLIC :: ComputeError
  PUBLIC :: OpaqueDiffusionDomainLength
  PUBLIC :: OpaqueDiffusionOpacity
  PUBLIC :: OpaqueDiffusionDecayRate
  PUBLIC :: OpaqueDiffusionAmplitude
  PUBLIC :: OpaqueDiffusionBackground
  PUBLIC :: VariableOpacityDiffusionMeanOpacity
  PUBLIC :: SetVariableOpacityDiffusionMeanOpticalDepth
  PUBLIC :: VariableOpacityDiffusionOpacity
  PUBLIC :: VariableOpacityDiffusionDensitySource
  PUBLIC :: ScatteringDiffusionDensitySource
  PUBLIC :: AbsorptionDiffusionAbsorptionOpacity
  PUBLIC :: AbsorptionDiffusionEquilibriumDensity
  PUBLIC :: CoupledDiffusionAbsorptionOpacity
  PUBLIC :: CoupledDiffusionEquilibriumSlope
  PUBLIC :: CoupledDiffusionDecayRate

  REAL(DP), PARAMETER :: OpaqueDiffusionDomainLength = One * Kilometer
  REAL(DP), PARAMETER :: OpaqueDiffusionOpticalDepth = 1.0d2
  REAL(DP), PARAMETER :: OpaqueDiffusionAmplitude    = 1.0d-4
  REAL(DP), PARAMETER :: OpaqueDiffusionBackground   = 2.5d-1
  REAL(DP), SAVE :: VariableOpacityDiffusionMeanOpticalDepth = 1.0d3
  REAL(DP), PARAMETER :: VariableOpacityDiffusionContrast = 9.0d-1
  REAL(DP), PARAMETER :: AbsorptionDiffusionAbsorptionFraction = 5.0d-1
  REAL(DP), PARAMETER :: CoupledDiffusionEquilibriumSlope = One
  REAL(DP), PARAMETER :: CoupledDiffusionTotalBackground = 5.0d-1
  REAL(DP), PARAMETER :: CoupledDiffusionTotalAmplitude = 5.0d-5
  REAL(DP), SAVE :: VariableOpacityDiffusionMeanOpacity &
    = 1.0d3 / OpaqueDiffusionDomainLength
  REAL(DP), SAVE :: ScatteringDiffusionSourceMean = Zero
  REAL(DP), PARAMETER :: &
    OpaqueDiffusionOpacity &
      = OpaqueDiffusionOpticalDepth / OpaqueDiffusionDomainLength

CONTAINS


  SUBROUTINE SetVariableOpacityDiffusionMeanOpticalDepth( OpticalDepth )

    REAL(DP), INTENT(in) :: OpticalDepth

    IF( OpticalDepth <= Zero ) &
      ERROR STOP 'Variable-opacity diffusion optical depth must be positive.'

    VariableOpacityDiffusionMeanOpticalDepth = OpticalDepth
    VariableOpacityDiffusionMeanOpacity &
      = VariableOpacityDiffusionMeanOpticalDepth &
          / OpaqueDiffusionDomainLength

  END SUBROUTINE SetVariableOpacityDiffusionMeanOpticalDepth


  SUBROUTINE InitializeFields( ProfileName, FileName, DnuMax_Option )

    CHARACTER(LEN=*), INTENT(in) :: ProfileName
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: FileName
    REAL(DP), INTENT(in), OPTIONAL :: DnuMax_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE( TRIM( ProgramName ) )

       CASE( 'HomogeneousEquilibrium' )

         IF( .NOT. PRESENT( DnuMax_Option ) ) &
           ERROR STOP 'HomogeneousEquilibrium requires DnuMax_Option.'

         CALL InitializeFields_LocalEquilibrium( '', DnuMax_Option )

       CASE( 'VaryingLocalEquilibrium' )

         IF( .NOT. PRESENT( DnuMax_Option ) ) &
           ERROR STOP 'VaryingLocalEquilibrium requires DnuMax_Option.'

         CALL InitializeFields_LocalEquilibrium &
                ( ProfileName, DnuMax_Option )

       CASE( 'ManufacturedOpaqueDiffusion' )

         IF( .NOT. PRESENT( DnuMax_Option ) ) &
           ERROR STOP 'ManufacturedOpaqueDiffusion requires DnuMax_Option.'

         CALL InitializeFields_ManufacturedOpaqueDiffusion &
                ( DnuMax_Option )

       CASE( 'ManufacturedVariableDiffusion' )

         IF( .NOT. PRESENT( DnuMax_Option ) ) &
           ERROR STOP &
             'ManufacturedVariableDiffusion requires DnuMax_Option.'

         CALL InitializeFields_ManufacturedVariableOpacityDiffusion &
                ( DnuMax_Option )

       CASE( 'ScatteringDiffusionLimit' )

         IF( .NOT. PRESENT( DnuMax_Option ) ) &
           ERROR STOP &
             'ScatteringDiffusionLimit requires DnuMax_Option.'

         CALL InitializeScatteringDiffusionSourceMean

         CALL InitializeFields_ManufacturedVariableOpacityDiffusion &
                ( DnuMax_Option )

       CASE( 'AbsorptionDiffusionBalance' )

         IF( .NOT. PRESENT( DnuMax_Option ) ) &
           ERROR STOP 'AbsorptionDiffusionBalance requires DnuMax_Option.'

         CALL InitializeFields_AbsorptionDiffusionBalance &
                ( DnuMax_Option )

       CASE( 'CoupledAbsorptionDiffusion' )

         IF( .NOT. PRESENT( DnuMax_Option ) ) &
           ERROR STOP 'CoupledAbsorptionDiffusion requires DnuMax_Option.'

         CALL InitializeFields_CoupledAbsorptionDiffusion &
                ( DnuMax_Option )

       CASE( 'Relaxation' )

         CALL InitializeFields_Relaxation(FileName)

       CASE( 'DeleptonizationWave1D' )

         CALL InitializeFields_DeleptonizationWave1D( ProfileName )

       CASE( 'EquilibriumAdvection' )

         CALL InitializeFields_EquilibriumAdvection

    END SELECT

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( uAF, uPF, uCF, uPR, uCR )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( uAF, uPF, uCF, uPR, uCR )
#endif

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_Relaxation(FileName)

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: FileName

    INTEGER  :: iE, iX1, iX2, iX3, iS, iNodeE, iNodeX, iNodeZ
    INTEGER  :: iCR, iN_E
    REAL(DP) :: kT, E, f_E

    REAL(DP) :: D_0, T_0, Y_0, E_0, V_u_1, V_u_2, V_u_3
    REAL(DP) :: Mu_0

    REAL(DP), DIMENSION(iE_E0*nDOFE) :: tmp_nu

    INTEGER  :: tmp_int, nE_G


    IF( PRESENT (FileName) ) THEN

      OPEN( UNIT=17,file=FileName,&
            ACCESS='sequential', FORM='formatted', ACTION='read')

      READ(17,*) D_0, T_0, Y_0, E_0, V_u_1, V_u_2, V_u_3
      READ(17,*) nE_G

      IF((iE_E0 - iE_B0 + 1) * nDOFE /= nE_G) THEN
        WRITE(*,*) 'The energy grids in the data file and simulation are different, aborting.'
        STOP
      ENDIF

      D_0   = D_0 * Gram / Centimeter**3
      T_0   = T_0 * MeV
      V_u_1 = V_u_1 * SpeedOfLight
      V_u_2 = V_u_2 * SpeedOfLight
      V_u_3 = V_u_3 * SpeedOfLight

      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
        READ(17,*) tmp_int, tmp_nu

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNodeX = 1, nDOFX
        DO iE  = iE_B0   , iE_E0

        DO iNodeE = 1, nDOFE

          iNodeZ = (iNodeX - 1 ) * nDOFE + iNodeE
          iN_E   = (iE     - 1 ) * nDOFE + iNodeE
           
          uPR(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) = tmp_nu(iN_E)
     
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO
      END DO
      END DO

      CLOSE(17)


      ! --- Fluid Fields Conserved Variables ---

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

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_u_1
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_u_2
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_u_3

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
                   uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

        END DO

      END DO
      END DO
      END DO

      ! --- Radiation Fields Conserved Variables ---

      DO iS  = 1       , nSpecies
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iE  = iE_B0   , iE_E0

        DO iNodeX = 1, nDOFX
        DO iNodeE = 1, nDOFE

          iNodeZ = (iNodeX-1) * nDOFE + iNodeE

          CALL ComputeConserved_TwoMoment &
                 ( uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS), &
                   uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS), &
                   uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS), &
                   uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS), &
                   uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
                   uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
                   uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
                   uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
                   uPF(iNodeX   ,iX1,iX2,iX3,iPF_V1),        &
                   uPF(iNodeX   ,iX1,iX2,iX3,iPF_V2),        &
                   uPF(iNodeX   ,iX1,iX2,iX3,iPF_V3),        &
                   uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_11),  &
                   uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_22),  &
                   uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_33) )

          END DO
          END DO

        END DO
        END DO
        END DO
        END DO
        END DO

    ELSE

      D_0   = 1.032d12 * Gram / Centimeter**3
      T_0   = 7.588d0 * MeV
      Y_0   = 0.1347_DP
      V_u_1 = 0.1_DP * SpeedOfLight
      V_u_2 = 0.0_DP * SpeedOfLight
      V_u_3 = 0.0_DP * SpeedOfLight
      Mu_0  = 0.0_DP ! \in [-1,1]

      ! --- Fluid Fields Conserved Variables ---

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

          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_u_1
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_u_2
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_u_3

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
                   uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

        END DO

      END DO
      END DO
      END DO

      ! --- Radiation Fields ---

      DO iS  = 1       , nSpecies
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iE  = iE_B0   , iE_E0

        DO iNodeX = 1, nDOFX
        DO iNodeE = 1, nDOFE

          iNodeZ = (iNodeX-1) * nDOFE + iNodeE

          kT = BoltzmannConstant * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

          E = NodeCoordinate( MeshE, iE, iNodeE )

          f_E = MAX( 0.99_DP * EXP( - ( E - Two*kT )**2 &
                                      / ( Two*(1.0d1*MeV)**2 ) ), 1.0d-99 )

          uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D,iS) &
            = f_E * 0.50_DP * ( One - Mu_0 )

          uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) &
            = f_E * 0.25_DP * ( One - Mu_0**2 )

          uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
          uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

          CALL ComputeConserved_TwoMoment &
                 ( uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS), &
                   uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS), &
                   uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS), &
                   uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS), &
                   uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
                   uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
                   uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
                   uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
                   uPF(iNodeX   ,iX1,iX2,iX3,iPF_V1),        &
                   uPF(iNodeX   ,iX1,iX2,iX3,iPF_V2),        &
                   uPF(iNodeX   ,iX1,iX2,iX3,iPF_V3),        &
                   uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_11),  &
                   uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_22),  &
                   uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_33) )

        END DO
        END DO

      END DO
      END DO
      END DO
      END DO
      END DO

    END IF
   
  END SUBROUTINE InitializeFields_Relaxation


  SUBROUTINE InitializeFields_DeleptonizationWave1D( ProfileName )

    CHARACTER(LEN=*), INTENT(in) :: ProfileName

    INTEGER  :: iX1, iX2, iX3, iS, iNodeE, iNodeX, iNodeX1, iNodeZ
    INTEGER  :: i, iR, iE, nR, nE
    REAL(DP) :: R, Tau
    REAL(DP), ALLOCATABLE :: R_P(:), D_P(:), T_P(:), Y_P(:)
    REAL(DP), ALLOCATABLE :: E_Nu(:), R_Nu(:,:), Chi(:,:,:), fEQ(:,:,:)
    REAL(DP), ALLOCATABLE :: D_Nu_P(:,:,:), I1_Nu_P(:,:,:)

    LOGICAL, ALLOCATABLE :: Opacity_Mask(:,:)

    WRITE(*,*)
    WRITE(*,'(A6,A,A)') '', &
      'Initializing from Profile: ', TRIM( ProfileName )

    CALL ReadFluidProfile( ProfileName, R_P, D_P, T_P, Y_P )

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        iR = MAX( Locate( R, R_P, SIZE( R_P ) ), 1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
          = Interpolate1D_Linear( R, R_P(iR), R_P(iR+1), D_P(iR), D_P(iR+1) )

        uAF(iNodeX,iX1,iX2,iX3,iAF_T ) &
          = Interpolate1D_Linear( R, R_P(iR), R_P(iR+1), T_P(iR), T_P(iR+1) )

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
          = Interpolate1D_Linear( R, R_P(iR), R_P(iR+1), Y_P(iR), Y_P(iR+1) )

        CALL ComputeThermodynamicStates_Primitive_TABLE &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) )

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

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

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
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

    nR = SIZE( R_P )
    nE = ( iE_E0 - iE_B0 + 1 ) * nDOFE

    ALLOCATE( E_Nu   (nE) )
    ALLOCATE( R_Nu   (nE,nSpecies) )
    ALLOCATE( Chi    (nE,nSpecies,nR) )
    ALLOCATE( fEQ    (nE,nSpecies,nR) )
    ALLOCATE( D_Nu_P (nE,nSpecies,nR) )
    ALLOCATE( I1_Nu_P(nE,nSpecies,nR) )

    ALLOCATE( Opacity_Mask(nOp,nR) )

    CALL SetOpacityMask &
           ( nR, 1, D_P, T_P, Y_P, Opacity_Mask )

    Chi = Zero

    ! --- Neutrino Energies ---

    DO iE = iE_B0, iE_E0
    DO iNodeE = 1, nDOFE

      E_Nu((iE-1)*nDOFE+iNodeE) = NodeCoordinate( MeshE, iE, iNodeE )

    END DO
    END DO

    ! --- Neutrino Absorption Opacities and Equilibrium Distributions ---

    CALL ComputeEquilibriumDistributions_DG &
           ( 1, nE, 1, nSpecies, 1, nR, E_Nu, D_P, T_P, Y_P, fEQ )

    CALL ComputeNeutrinoOpacities_EC &
           ( 1, nE, 1, nSpecies, 1, nR, E_Nu, D_P, T_P, Y_P, fEQ, Chi, Opacity_Mask )

    DO iR = 1, nR

      ! --- Prevent too large drop-off of the opacity -------
      ! --- This is mainly to prevent opacity for Nue_Bar ---
      ! --- be close to zero for low neutrino energies ------

      DO iS = 1, nSpecies
      DO iE = 1, nE-1

        IF( Chi(iE,iS,iR) < 1.d-16 * Chi(iE+1,iS,iR) )THEN

          Chi(1:iE,iS,iR) = Chi(iE+1,iS,iR) * ( E_Nu(1:iE) / E_Nu(iE+1) )**2

        END IF

      END DO
      END DO

    END DO

    ! --- Approximate Neutrino Sphere Radii ---

    DO iS = 1, nSpecies
    DO iE = 1, nE

      Tau = Zero
      DO iR = nR-1, 1, -1

        IF( Tau > 2.0_DP / 3.0_DP ) CYCLE

        Tau = Tau + Half * ( R_P(iR+1) - R_P(iR) ) &
                         * ( Chi(iE,iS,iR+1) + Chi(iE,iS,iR) )

        R_Nu(iE,iS) = MAX( R_P(iR), 1.0d1 * Kilometer )

      END DO

    END DO
    END DO

    ! --- Homogeneous Sphere Solution ---

    DO iR = 1, nR
    DO iS = 1, nSpecies
    DO iE = 1, nE

      ! --- Use Chi and fEQ at MIN( local radius, neutrino sphere radius ) ---

      i = MIN( iR, Locate( R_Nu(iE,iS), R_P, SIZE( R_P ) ) )

      CALL ComputeSphereSolution &
             ( R_Nu(iE,iS), Chi(iE,iS,i), fEQ(iE,iS,i), &
               R_P(iR), D_Nu_P(iE,iS,iR), I1_Nu_P(iE,iS,iR) )

      D_Nu_P(iE,iS,iR) = MAX( D_Nu_P(iE,iS,iR), SqrtTiny )

    END DO
    END DO
    END DO

    DO iS  = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        iNodeX1 = NodeNumberTable(2,iNodeZ)

        R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        iR = MAX( Locate( R, R_P, SIZE( R_P ) ), 1 )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) &
          = Interpolate1D_Linear &
              ( R, R_P(iR), R_P(iR+1), &
                D_Nu_P ((iE-1)*nDOFE+iNodeE,iS,iR  ), &
                D_Nu_P ((iE-1)*nDOFE+iNodeE,iS,iR+1) )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) &
          = Interpolate1D_Linear &
              ( R, R_P(iR), R_P(iR+1), &
                I1_Nu_P((iE-1)*nDOFE+iNodeE,iS,iR  ), &
                I1_Nu_P((iE-1)*nDOFE+iNodeE,iS,iR+1) )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V1),        &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V2),        &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V3),        &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_11),  &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_22),  &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    DEALLOCATE( R_P, D_P, T_P, Y_P )
    DEALLOCATE( E_Nu, R_Nu, Chi, fEQ )
    DEALLOCATE( D_Nu_P, I1_Nu_P )
    DEALLOCATE( Opacity_Mask )

  END SUBROUTINE InitializeFields_DeleptonizationWave1D


  SUBROUTINE InitializeFields_LocalEquilibrium( ProfileName, DnuMax )

    CHARACTER(LEN=*), INTENT(in) :: ProfileName
    REAL(DP), INTENT(in) :: DnuMax

    REAL(DP), PARAMETER :: D_0 = 3.0d14 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: T_0 = 1.2d11 * Kelvin
    REAL(DP), PARAMETER :: Y_0 = 0.3_DP

    INTEGER :: iE, iN_E, iN_X, iNodeE, iNodeX, iNodeX1, iNodeZ
    INTEGER :: iR, iS, iX1, iX2, iX3
    INTEGER :: nE_G, nX_G, nX_Local(3)
    REAL(DP) :: R
    REAL(DP), ALLOCATABLE :: E_N(:), D_N(:), T_N(:), Y_N(:)
    REAL(DP), ALLOCATABLE :: Dnu_EQ(:,:,:)
    REAL(DP), ALLOCATABLE :: R_P(:), D_P(:), T_P(:), Y_P(:)
    LOGICAL :: UseProfile

    UseProfile = LEN_TRIM( ProfileName ) > 0

    IF( UseProfile )THEN
      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', &
        'Initializing local equilibrium from Profile: ', TRIM( ProfileName )
      CALL ReadFluidProfile( ProfileName, R_P, D_P, T_P, Y_P )
    END IF

    nE_G = nDOFE * ( iE_E0 - iE_B0 + 1 )
    nX_Local = iX_E0 - iX_B0 + 1
    nX_G = nDOFX * PRODUCT( nX_Local )

    ALLOCATE( E_N(nE_G), D_N(nX_G), T_N(nX_G), Y_N(nX_G) )
    ALLOCATE( Dnu_EQ(nE_G,nSpecies,nX_G) )

    DO iE = iE_B0, iE_E0
    DO iNodeE = 1, nDOFE

      iN_E = ( iE - iE_B0 ) * nDOFE + iNodeE
      E_N(iN_E) = NodeCoordinate( MeshE, iE, iNodeE )

    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      iN_X = iNodeX &
             + ( iX1 - iX_B0(1) ) * nDOFX &
             + ( iX2 - iX_B0(2) ) * nDOFX * nX_Local(1) &
             + ( iX3 - iX_B0(3) ) * nDOFX &
                 * nX_Local(1) * nX_Local(2)

      IF( UseProfile )THEN

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        iR = MIN( MAX( Locate( R, R_P, SIZE(R_P) ), 1 ), SIZE(R_P)-1 )

        D_N(iN_X) = Interpolate1D_Linear &
          ( R, R_P(iR), R_P(iR+1), D_P(iR), D_P(iR+1) )
        T_N(iN_X) = Interpolate1D_Linear &
          ( R, R_P(iR), R_P(iR+1), T_P(iR), T_P(iR+1) )
        Y_N(iN_X) = Interpolate1D_Linear &
          ( R, R_P(iR), R_P(iR+1), Y_P(iR), Y_P(iR+1) )

      ELSE

        D_N(iN_X) = D_0
        T_N(iN_X) = T_0
        Y_N(iN_X) = Y_0

      END IF

      uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = D_N(iN_X)
      uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
      uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
      uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero
      uAF(iNodeX,iX1,iX2,iX3,iAF_T ) = T_N(iN_X)
      uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = Y_N(iN_X)

      CALL ComputeThermodynamicStates_Primitive_TABLE &
             ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
               uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) )

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
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO
    END DO

    CALL ComputeEquilibriumDistributions_DG &
           ( 1, nE_G, 1, nSpecies, 1, nX_G, &
             E_N, D_N, T_N, Y_N, Dnu_EQ, DnuMax )

    DO iS  = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0, iE_E0
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iN_E = ( iE - iE_B0 ) * nDOFE + iNodeE
      iN_X = iNodeX &
             + ( iX1 - iX_B0(1) ) * nDOFX &
             + ( iX2 - iX_B0(2) ) * nDOFX * nX_Local(1) &
             + ( iX3 - iX_B0(3) ) * nDOFX &
                 * nX_Local(1) * nX_Local(2)
      iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) = Dnu_EQ(iN_E,iS,iN_X)
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) = Zero
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

      CALL ComputeConserved_TwoMoment &
             ( uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS), &
               uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS), &
               uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS), &
               uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    DEALLOCATE( E_N, D_N, T_N, Y_N, Dnu_EQ )
    IF( UseProfile ) DEALLOCATE( R_P, D_P, T_P, Y_P )

  END SUBROUTINE InitializeFields_LocalEquilibrium


  SUBROUTINE InitializeFields_ManufacturedOpaqueDiffusion( DnuMax )

    REAL(DP), INTENT(in) :: DnuMax

    INTEGER :: iE, iS, iX1, iX2, iX3
    INTEGER :: iNodeE, iNodeX, iNodeX1, iNodeZ
    REAL(DP) :: X, WaveNumber, DecayRate, FluxAmplitude
    REAL(DP) :: Dnu, Inu_1

    CALL InitializeFields_LocalEquilibrium( '', DnuMax )

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength
    DecayRate = OpaqueDiffusionDecayRate()
    FluxAmplitude &
      = - DecayRate * OpaqueDiffusionBackground &
          * OpaqueDiffusionAmplitude / WaveNumber

    DO iS  = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0, iE_E0
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeZ  = ( iNodeX - 1 ) * nDOFE + iNodeE
      X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

      IF( iS == 1 )THEN
        Dnu = OpaqueDiffusionBackground &
          * ( One + OpaqueDiffusionAmplitude * COS( WaveNumber * X ) )
        Inu_1 = FluxAmplitude * SIN( WaveNumber * X )
      ELSE
        Dnu = Half * OpaqueDiffusionBackground
        Inu_1 = Zero
      END IF

      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) = Dnu
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) = Inu_1
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

      CALL ComputeConserved_TwoMoment &
             ( Dnu, Inu_1, Zero, Zero, &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    WRITE(*,* )
    WRITE(*,'(A)') 'Manufactured opaque-diffusion initial data:'
    WRITE(*,'(A,ES14.6E3)') '  domain optical depth         = ', &
      OpaqueDiffusionOpticalDepth
    WRITE(*,'(A,ES14.6E3)') '  mean free path [km]          = ', &
      ( One / OpaqueDiffusionOpacity ) / Kilometer
    WRITE(*,'(A,ES14.6E3)') '  perturbation amplitude       = ', &
      OpaqueDiffusionAmplitude
    WRITE(*,'(A,ES14.6E3)') '  slow-mode decay rate [1/ms]  = ', &
      OpaqueDiffusionDecayRate() * Millisecond
    WRITE(*,* )

  END SUBROUTINE InitializeFields_ManufacturedOpaqueDiffusion


  SUBROUTINE InitializeFields_ManufacturedVariableOpacityDiffusion &
    ( DnuMax )

    REAL(DP), INTENT(in) :: DnuMax

    INTEGER :: iE, iS, iX1, iX2, iX3
    INTEGER :: iNodeE, iNodeX, iNodeX1, iNodeZ
    REAL(DP) :: X, WaveNumber, Opacity
    REAL(DP) :: Dnu, Inu_1

    CALL InitializeFields_LocalEquilibrium( '', DnuMax )

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength

    DO iS  = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0, iE_E0
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeZ  = ( iNodeX - 1 ) * nDOFE + iNodeE
      X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

      IF( iS == 1 )THEN
        Opacity = VariableOpacityDiffusionOpacity( X )
        Dnu = OpaqueDiffusionBackground &
          * ( One + OpaqueDiffusionAmplitude * COS( WaveNumber * X ) )
        Inu_1 = OpaqueDiffusionBackground * OpaqueDiffusionAmplitude &
          * WaveNumber * SIN( WaveNumber * X ) / ( Three * Opacity )
      ELSE
        Dnu = Half * OpaqueDiffusionBackground
        Inu_1 = Zero
      END IF

      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) = Dnu
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) = Inu_1
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

      CALL ComputeConserved_TwoMoment &
             ( Dnu, Inu_1, Zero, Zero, &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    WRITE(*,* )
    IF( TRIM( ProgramName ) == 'ScatteringDiffusionLimit' )THEN
      WRITE(*,'(A)') &
        'Variable-opacity scattering diffusion-limit initial data:'
    ELSE
      WRITE(*,'(A)') 'Manufactured variable-opacity diffusion initial data:'
    END IF
    WRITE(*,'(A,ES14.6E3)') '  mean domain optical depth    = ', &
      VariableOpacityDiffusionMeanOpticalDepth
    WRITE(*,'(A,ES14.6E3)') '  opacity contrast             = ', &
      VariableOpacityDiffusionContrast
    WRITE(*,'(A,2(1X,ES14.6E3))') '  opacity min/max [1/km]       = ', &
      VariableOpacityDiffusionMeanOpacity &
        * ( One - VariableOpacityDiffusionContrast ) * Kilometer, &
      VariableOpacityDiffusionMeanOpacity &
        * ( One + VariableOpacityDiffusionContrast ) * Kilometer
    WRITE(*,* )

  END SUBROUTINE InitializeFields_ManufacturedVariableOpacityDiffusion


  SUBROUTINE InitializeFields_AbsorptionDiffusionBalance( DnuMax )

    REAL(DP), INTENT(in) :: DnuMax

    INTEGER :: iE, iS, iX1, iX2, iX3
    INTEGER :: iNodeE, iNodeX, iNodeX1, iNodeZ
    REAL(DP) :: X, WaveNumber, Opacity
    REAL(DP) :: CellWidthMin, CellWidthMax
    REAL(DP) :: Dnu, Inu_1

    CALL InitializeFields_LocalEquilibrium( '', DnuMax )

    ! The first species is initialized to the exact stationary balance
    !
    !   I_* = -(3*kappa_tr)^{-1} d_x D_*,
    !   d_x I_* = kappa_a ( D_0 - D_* ).

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength

    DO iS  = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0, iE_E0
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeZ  = ( iNodeX - 1 ) * nDOFE + iNodeE
      X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

      IF( iS == 1 )THEN
        Opacity = VariableOpacityDiffusionOpacity( X )
        Dnu = OpaqueDiffusionBackground &
          * ( One + OpaqueDiffusionAmplitude * COS( WaveNumber * X ) )
        Inu_1 = OpaqueDiffusionBackground * OpaqueDiffusionAmplitude &
          * WaveNumber * SIN( WaveNumber * X ) / ( Three * Opacity )
      ELSE
        Dnu = Half * OpaqueDiffusionBackground
        Inu_1 = Zero
      END IF

      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) = Dnu
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) = Inu_1
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

      CALL ComputeConserved_TwoMoment &
             ( Dnu, Inu_1, Zero, Zero, &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    WRITE(*,* )
    WRITE(*,'(A)') 'Absorption-diffusion balance initial data:'
    WRITE(*,'(A,ES14.6E3)') '  mean transport optical depth = ', &
      VariableOpacityDiffusionMeanOpticalDepth
    CellWidthMin = MINVAL &
      ( MeshX(1) % Width(iX_B0(1):iX_E0(1)) )
    CellWidthMax = MAXVAL &
      ( MeshX(1) % Width(iX_B0(1):iX_E0(1)) )
    WRITE(*,'(A,2(1X,ES14.6E3))') &
      '  cell optical-depth min/max = ', &
      VariableOpacityDiffusionMeanOpacity &
        * ( One - VariableOpacityDiffusionContrast ) * CellWidthMin, &
      VariableOpacityDiffusionMeanOpacity &
        * ( One + VariableOpacityDiffusionContrast ) * CellWidthMax
    WRITE(*,'(A,ES14.6E3)') '  opacity contrast             = ', &
      VariableOpacityDiffusionContrast
    WRITE(*,'(A,ES14.6E3)') '  absorption fraction          = ', &
      AbsorptionDiffusionAbsorptionFraction
    WRITE(*,* )

  END SUBROUTINE InitializeFields_AbsorptionDiffusionBalance


  SUBROUTINE InitializeFields_CoupledAbsorptionDiffusion( DnuMax )

    REAL(DP), INTENT(in) :: DnuMax

    INTEGER :: iE, iS, iX1, iX2, iX3
    INTEGER :: iNodeE, iNodeX, iNodeX1, iNodeZ
    REAL(DP) :: X, Angle, WaveNumber, DecayRate
    REAL(DP) :: RadiationFraction, EquilibriumRadiationFraction
    REAL(DP) :: Dnu, Matter, Inu_1

    IF( nSpecies /= 2 ) &
      ERROR STOP 'CoupledAbsorptionDiffusion requires exactly two species.'

    CALL InitializeFields_LocalEquilibrium( '', DnuMax )

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength
    DecayRate = CoupledDiffusionDecayRate()
    RadiationFraction = CoupledDiffusionRadiationFraction( DecayRate )
    EquilibriumRadiationFraction &
      = CoupledDiffusionEquilibriumSlope &
          / ( One + CoupledDiffusionEquilibriumSlope )

    DO iS  = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0, iE_E0
    DO iNodeX = 1, nDOFX
    DO iNodeE = 1, nDOFE

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeZ  = ( iNodeX - 1 ) * nDOFE + iNodeE
      X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
      Angle = WaveNumber * X

      IF( iS == 1 )THEN
        Dnu = EquilibriumRadiationFraction &
                * CoupledDiffusionTotalBackground &
              + RadiationFraction * CoupledDiffusionTotalAmplitude &
                  * COS( Angle )
        Inu_1 = - DecayRate * CoupledDiffusionTotalAmplitude &
          * SIN( Angle ) / WaveNumber
      ELSE
        Matter = ( One - EquilibriumRadiationFraction ) &
                   * CoupledDiffusionTotalBackground &
                 + ( One - RadiationFraction ) &
                     * CoupledDiffusionTotalAmplitude * COS( Angle )
        Dnu = Matter
        Inu_1 = Zero
      END IF

      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) = Dnu
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) = Inu_1
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
      uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

      CALL ComputeConserved_TwoMoment &
             ( Dnu, Inu_1, Zero, Zero, &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
               uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
               uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    WRITE(*,* )
    WRITE(*,'(A)') 'Coupled absorption-diffusion initial data:'
    WRITE(*,'(A,ES14.6E3)') '  transport optical depth      = ', &
      VariableOpacityDiffusionMeanOpticalDepth
    WRITE(*,'(A,ES14.6E3)') '  absorption fraction          = ', &
      AbsorptionDiffusionAbsorptionFraction
    WRITE(*,'(A,ES14.6E3)') '  equilibrium slope dD0/dM     = ', &
      CoupledDiffusionEquilibriumSlope
    WRITE(*,'(A,ES14.6E3)') '  exact slow decay rate [1/ms] = ', &
      DecayRate * Millisecond
    WRITE(*,* )

  END SUBROUTINE InitializeFields_CoupledAbsorptionDiffusion


  REAL(DP) FUNCTION VariableOpacityDiffusionOpacity( X ) RESULT( Opacity )

    REAL(DP), INTENT(in) :: X

    Opacity = VariableOpacityDiffusionMeanOpacity &
      * ( One + VariableOpacityDiffusionContrast &
          * COS( TwoPi * X / OpaqueDiffusionDomainLength ) )

  END FUNCTION VariableOpacityDiffusionOpacity


  REAL(DP) FUNCTION VariableOpacityDiffusionDensitySource( X ) &
    RESULT( Source )

    REAL(DP), INTENT(in) :: X

    REAL(DP) :: Angle, Denominator, WaveNumber

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength
    Angle = WaveNumber * X
    Denominator = One &
      + VariableOpacityDiffusionContrast * COS( Angle )

    Source = OpaqueDiffusionBackground * OpaqueDiffusionAmplitude &
      * WaveNumber**2 / ( Three * VariableOpacityDiffusionMeanOpacity ) &
      * (   COS( Angle ) / Denominator &
          + VariableOpacityDiffusionContrast * SIN( Angle )**2 &
              / Denominator**2 )

  END FUNCTION VariableOpacityDiffusionDensitySource


  REAL(DP) FUNCTION ScatteringDiffusionDensitySource( X ) &
    RESULT( Source )

    REAL(DP), INTENT(in) :: X

    Source = VariableOpacityDiffusionDensitySource( X ) &
      - ScatteringDiffusionSourceMean

  END FUNCTION ScatteringDiffusionDensitySource


  SUBROUTINE InitializeScatteringDiffusionSourceMean

    INTEGER :: iX1, iNodeX, iNodeX1
    REAL(DP) :: X, Weight, DomainVolume, SourceScale

    ScatteringDiffusionSourceMean = Zero
    DomainVolume = Zero

    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
      Weight = WeightsX1(iNodeX1) * MeshX(1) % Width(iX1)

      ScatteringDiffusionSourceMean &
        = ScatteringDiffusionSourceMean &
            + Weight * VariableOpacityDiffusionDensitySource( X )
      DomainVolume = DomainVolume + Weight

    END DO
    END DO

    ScatteringDiffusionSourceMean &
      = ScatteringDiffusionSourceMean / DomainVolume

    SourceScale = OpaqueDiffusionBackground * OpaqueDiffusionAmplitude &
      * ( TwoPi / OpaqueDiffusionDomainLength )**2 &
      / ( Three * VariableOpacityDiffusionMeanOpacity )

    WRITE(*,'(A,ES14.6E3)') &
      '  removed discrete source mean / source scale = ', &
      ScatteringDiffusionSourceMean / SourceScale
    WRITE(*,* )

  END SUBROUTINE InitializeScatteringDiffusionSourceMean


  REAL(DP) FUNCTION AbsorptionDiffusionAbsorptionOpacity( X ) &
    RESULT( Opacity )

    REAL(DP), INTENT(in) :: X

    Opacity = AbsorptionDiffusionAbsorptionFraction &
      * VariableOpacityDiffusionOpacity( X )

  END FUNCTION AbsorptionDiffusionAbsorptionOpacity


  REAL(DP) FUNCTION AbsorptionDiffusionEquilibriumDensity( X ) &
    RESULT( Dnu_0 )

    REAL(DP), INTENT(in) :: X

    REAL(DP) :: Dnu, WaveNumber

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength
    Dnu = OpaqueDiffusionBackground &
      * ( One + OpaqueDiffusionAmplitude * COS( WaveNumber * X ) )

    ! VariableOpacityDiffusionDensitySource is d_x I_*.
    Dnu_0 = Dnu + VariableOpacityDiffusionDensitySource( X ) &
      / AbsorptionDiffusionAbsorptionOpacity( X )

  END FUNCTION AbsorptionDiffusionEquilibriumDensity


  REAL(DP) FUNCTION CoupledDiffusionAbsorptionOpacity() RESULT( Opacity )

    Opacity = AbsorptionDiffusionAbsorptionFraction &
      * VariableOpacityDiffusionMeanOpacity

  END FUNCTION CoupledDiffusionAbsorptionOpacity


  REAL(DP) FUNCTION CoupledDiffusionDecayRate() RESULT( DecayRate )

    INTEGER, PARAMETER :: MaxIterations = 50
    INTEGER :: iIteration
    REAL(DP) :: WaveNumber, TransportOpacity, AbsorptionOpacity
    REAL(DP) :: EquilibriumFraction, A, F, dF, Correction

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength
    TransportOpacity = VariableOpacityDiffusionMeanOpacity
    AbsorptionOpacity = CoupledDiffusionAbsorptionOpacity()
    EquilibriumFraction = CoupledDiffusionEquilibriumSlope &
      / ( One + CoupledDiffusionEquilibriumSlope )
    A = AbsorptionOpacity &
      * ( One + CoupledDiffusionEquilibriumSlope )

    DecayRate = - EquilibriumFraction * WaveNumber**2 &
      / ( Three * TransportOpacity )

    DO iIteration = 1, MaxIterations

      F = Three * DecayRate * ( DecayRate + TransportOpacity ) &
            * ( DecayRate + A ) &
          + WaveNumber**2 &
              * ( DecayRate &
                  + AbsorptionOpacity * CoupledDiffusionEquilibriumSlope )
      dF = Three &
             * (   ( DecayRate + TransportOpacity ) * ( DecayRate + A ) &
                 + DecayRate * ( DecayRate + A ) &
                 + DecayRate * ( DecayRate + TransportOpacity ) ) &
           + WaveNumber**2
      Correction = F / dF
      DecayRate = DecayRate - Correction

      IF( ABS( Correction ) &
            <= 10.0_DP * EPSILON( One ) * MAX( ABS( DecayRate ), One ) ) &
        EXIT

    END DO

    IF( DecayRate >= Zero .OR. iIteration > MaxIterations ) &
      ERROR STOP 'Unable to determine coupled diffusion slow decay rate.'

  END FUNCTION CoupledDiffusionDecayRate


  REAL(DP) FUNCTION CoupledDiffusionRadiationFraction( DecayRate ) &
    RESULT( Fraction )

    REAL(DP), INTENT(in) :: DecayRate

    REAL(DP) :: AbsorptionOpacity

    AbsorptionOpacity = CoupledDiffusionAbsorptionOpacity()
    Fraction &
      = ( DecayRate &
          + AbsorptionOpacity * CoupledDiffusionEquilibriumSlope ) &
          / ( DecayRate &
              + AbsorptionOpacity &
                  * ( One + CoupledDiffusionEquilibriumSlope ) )

  END FUNCTION CoupledDiffusionRadiationFraction


  REAL(DP) FUNCTION OpaqueDiffusionDecayRate() RESULT( DecayRate )

    REAL(DP) :: WaveNumber, Q, Discriminant

    ! For d_t D + d_x I = 0, d_t I + d_x(D/3) = -kappa I,
    ! the slow periodic mode satisfies s^2 + kappa*s + k^2/3 = 0.
    ! This form avoids cancellation in the opaque limit.
    WaveNumber = TwoPi / OpaqueDiffusionDomainLength
    Q = WaveNumber**2 / Three
    Discriminant = OpaqueDiffusionOpacity**2 - 4.0_DP * Q

    IF( Discriminant <= Zero ) &
      ERROR STOP 'Manufactured opaque-diffusion mode is not overdamped.'

    DecayRate = - Two * Q &
      / ( OpaqueDiffusionOpacity + SQRT( Discriminant ) )

  END FUNCTION OpaqueDiffusionDecayRate


  SUBROUTINE InitializeFields_EquilibriumAdvection

    REAL(DP), PARAMETER :: D_0   = 3.0d14 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: T_0   = 1.2d11 * Kelvin
    REAL(DP), PARAMETER :: Y_0   = 0.3_DP
    REAL(DP), PARAMETER :: V_u_1 = 0.1_DP * SpeedOfLight
    REAL(DP), PARAMETER :: V_u_2 = 0.0_DP * SpeedOfLight
    REAL(DP), PARAMETER :: V_u_3 = 0.0_DP * SpeedOfLight

    INTEGER  :: iE, iX1, iX2, iX3, iS, iNodeE, iNodeX, iNodeZ
    REAL(DP) :: kT, Mnu, E, f_E

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

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_u_1
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_u_2
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_u_3

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
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

    DO iS  = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        kT = BoltzmannConstant * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

        Mnu = uAF  (iNodeX,iX1,iX2,iX3,iAF_Me) &
              + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
              - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn)

        IF( iS == iNuE_Bar )THEN

          Mnu = - Mnu

        END IF

        E = NodeCoordinate( MeshE, iE, iNodeE )

        f_E = MAX( One / ( EXP( ( E - Mnu ) / kT ) + One ), 1.0d-99 )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) = f_E
        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) = Zero
        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V1),        &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V2),        &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V3),        &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_11),  &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_22),  &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_EquilibriumAdvection


  SUBROUTINE ComputeSphereSolution( R0, Chi, f0, R, D, I )

    REAL(DP), INTENT(in)  :: R0, Chi, f0, R
    REAL(DP), INTENT(out) :: D, I

    INTEGER, PARAMETER :: nMu = 2048
    INTEGER            :: iMu
    REAL(DP)           :: Mu(nMu), Distribution(nMu)

    DO iMu = 1, nMu

      Mu(iMu) = - One + Two * DBLE(iMu-1)/DBLE(nMu-1)

      Distribution(iMu) = f_A( R, Mu(iMu), R0, f0, Chi )

    END DO

    D = Half * TRAPEZ( nMu, Mu, Distribution )
    I = Half * TRAPEZ( nMu, Mu, Distribution * Mu )

  END SUBROUTINE ComputeSphereSolution


  REAL(DP) FUNCTION f_A( R, Mu, R0, f0, Chi )

    REAL(DP), INTENT(in) :: R, Mu, R0, f0, Chi

    REAL(DP) :: s

    IF( R < R0 )THEN
      s = ( R * Mu + R0 * SQRT( One - ( R / R0 )**2 * ( One - Mu**2 ) ) )
    ELSE
      IF( Mu >= SQRT( One - ( R0 / R )**2 ) )THEN
        s = ( Two * R0 * SQRT( One - ( R / R0 )**2 * ( One - Mu**2 ) ) )
      ELSE
        s = Zero
      END IF
    END IF

    f_A = f0 * ( One - EXP( - Chi * s ) )

    RETURN
  END FUNCTION f_A


  REAL(DP) FUNCTION TRAPEZ( n, x, y )

    INTEGER,  INTENT(in) :: n
    REAL(dp), INTENT(in) :: x(n), y(n)

    INTEGER :: i

    TRAPEZ = 0.0_dp
    DO i = 1, n - 1
      TRAPEZ = TRAPEZ + 0.5_dp * ( x(i+1) - x(i) ) * ( y(i) + y(i+1) )
    END DO

    RETURN
  END FUNCTION TRAPEZ


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


  SUBROUTINE ComputeError( t )

    REAL(DP), INTENT(in) :: t

    SELECT CASE( TRIM( ProgramName ) )

       CASE( 'Relaxation' )

         CALL ComputeError_Relaxation

       CASE( 'DeleptonizationWave1D' )

         CALL ComputeError_DeleptonizationWave1D

       CASE( 'EquilibriumAdvection' )

         CALL ComputeError_EquilibriumAdvection

       CASE( 'ManufacturedOpaqueDiffusion' )

         CALL ComputeError_ManufacturedOpaqueDiffusion( t )

       CASE( 'ManufacturedVariableDiffusion' )

         CALL ComputeError_ManufacturedVariableOpacityDiffusion( t )

       CASE( 'ScatteringDiffusionLimit' )

         CALL ComputeError_ManufacturedVariableOpacityDiffusion( t )

       CASE( 'AbsorptionDiffusionBalance' )

         CALL ComputeError_AbsorptionDiffusionBalance( t )

       CASE( 'CoupledAbsorptionDiffusion' )

         CALL ComputeError_CoupledAbsorptionDiffusion( t )

    END SELECT

  END SUBROUTINE ComputeError


  SUBROUTINE ComputeError_Relaxation

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeE, iNodeX, iNodeZ
    INTEGER  :: nE, nX(3), nE_P, nX_P, iE_P, iX_P
    REAL(DP), ALLOCATABLE :: E_P(:), D_P(:), T_P(:), Y_P(:), f0_P(:,:,:)
    REAL(DP) :: MaxError(nSpecies), N0

    nE = iE_E0 - iE_B0 + 1
    nX = iX_E0 - iX_B0 + 1

    nE_P = nE * nDOFE
    nX_P = PRODUCT( nX ) * nDOFX

    ALLOCATE( E_P(nE_P) )

    DO iE = iE_B0, iE_E0
    DO iNodeE = 1, nDOFE
      iE_P = iNodeE + ( iE - iE_B0 ) * nDOFE
      E_P(iE_P) = NodeCoordinate( MeshE, iE, iNodeE )
    END DO
    END DO

    ALLOCATE( D_P(nX_P) )
    ALLOCATE( T_P(nX_P) )
    ALLOCATE( Y_P(nX_P) )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNodeX = 1, nDOFX
      iX_P = iNodeX &
             + ( iX1 - iX_B0(1) ) * nDOFX &
             + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
             + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)
      D_P(iX_P) = uPF(iNodeX,iX1,iX2,iX3,iPF_D )
      T_P(iX_P) = uAF(iNodeX,iX1,iX2,iX3,iAF_T )
      Y_P(iX_P) = uAF(iNodeX,iX1,iX2,iX3,iAF_Ye)
    END DO
    END DO
    END DO
    END DO

    ALLOCATE( f0_P(nE_P,nSpecies,nX_P) )

    CALL ComputeEquilibriumDistributions_DG &
           ( 1, nE_P, 1, nSpecies, 1, nX_P, E_P, D_P, T_P, Y_P, f0_P )

    MaxError = Zero

    DO iS = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iE_P = iNodeE + ( iE - iE_B0 ) * nDOFE
        iX_P = iNodeX &
               + ( iX1 - iX_B0(1) ) * nDOFX &
               + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
               + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)

        N0 = f0_P(iE_P,iS,iX_P)

        iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

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


  SUBROUTINE ComputeError_DeleptonizationWave1D

    INTEGER  :: iS
    REAL(DP) :: MaxError(nSpecies)

    MaxError = Zero

    WRITE(*,*)
    WRITE(*,'(A2,A)') '', 'INFO: Error Check'
    WRITE(*,*)
    DO iS = 1, nSpecies
      WRITE(*,'(A4,A10,I2.2,A14,ES10.4E2)') &
      '', 'Species = ', iS, ', Inf Error = ', MaxError(iS)
    END DO
    WRITE(*,*)

  END SUBROUTINE ComputeError_DeleptonizationWave1D


  SUBROUTINE ComputeError_EquilibriumAdvection

    INTEGER  :: iS
    REAL(DP) :: MaxError(nSpecies)

    MaxError = Zero

    WRITE(*,*)
    WRITE(*,'(A2,A)') '', 'INFO: Error Check'
    WRITE(*,*)
    DO iS = 1, nSpecies
      WRITE(*,'(A4,A10,I2.2,A14,ES10.4E2)') &
      '', 'Species = ', iS, ', Inf Error = ', MaxError(iS)
    END DO
    WRITE(*,*)

  END SUBROUTINE ComputeError_EquilibriumAdvection


  SUBROUTINE ComputeError_ManufacturedOpaqueDiffusion( t )

    REAL(DP), INTENT(in) :: t

    INTEGER :: iNodeX, iNodeX1, iNodeZ, iX1
    INTEGER :: OutputUnit
    REAL(DP) :: X, Weight, CellError
    REAL(DP) :: WaveNumber, DecayRate, DecayFactor
    REAL(DP) :: DiffusionDecayFactor
    REAL(DP) :: DensityAmplitude, FluxAmplitude
    REAL(DP) :: DiffusionFluxAmplitude
    REAL(DP) :: DnuNumerical, DnuExact, DnuDiffusion
    REAL(DP) :: InuNumerical, InuExact, InuDiffusion
    REAL(DP) :: DensityL1, DensityL2, DensityLInf
    REAL(DP) :: FluxL1, FluxL2, FluxLInf
    REAL(DP) :: DensityCosCoefficient, FluxSinCoefficient
    REAL(DP) :: MeanDensityError, CheckerboardAmplitude
    REAL(DP) :: DiffusionDecayRate, RateDifference
    REAL(DP) :: DensityModeErrorP1, DensityModeErrorDiffusion
    REAL(DP) :: FluxModeErrorP1, FluxModeErrorDiffusion

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength
    DecayRate = OpaqueDiffusionDecayRate()
    DecayFactor = EXP( DecayRate * t )
    DensityAmplitude &
      = OpaqueDiffusionBackground * OpaqueDiffusionAmplitude
    FluxAmplitude = - DecayRate * DensityAmplitude / WaveNumber
    DiffusionDecayRate = - WaveNumber**2 &
      / ( Three * OpaqueDiffusionOpacity )
    DiffusionDecayFactor = EXP( DiffusionDecayRate * t )
    DiffusionFluxAmplitude &
      = DensityAmplitude * WaveNumber / ( Three * OpaqueDiffusionOpacity )

    DensityL1 = Zero
    DensityL2 = Zero
    DensityLInf = Zero
    FluxL1 = Zero
    FluxL2 = Zero
    FluxLInf = Zero
    DensityCosCoefficient = Zero
    FluxSinCoefficient = Zero
    MeanDensityError = Zero
    CheckerboardAmplitude = Zero

    OPEN( NEWUNIT = OutputUnit, &
          FILE = 'manufactured_opaque_diffusion.dat', &
          STATUS = 'REPLACE', ACTION = 'WRITE' )
    WRITE(OutputUnit,'(A)') &
      '# x_km D_num D_P1 D_diff I1_num I1_P1 I1_diff D_error I1_error'

    DO iX1 = iX_B0(1), iX_E0(1)

      CellError = Zero

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeZ = ( iNodeX - 1 ) * nDOFE + 1
        X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        Weight = WeightsX1(iNodeX1) * MeshX(1) % Width(iX1)

        DnuNumerical = uPR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iPR_D,1)
        InuNumerical = uPR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iPR_I1,1)
        DnuExact = OpaqueDiffusionBackground &
          + DensityAmplitude * DecayFactor * COS( WaveNumber * X )
        InuExact = FluxAmplitude * DecayFactor * SIN( WaveNumber * X )
        DnuDiffusion = OpaqueDiffusionBackground &
          + DensityAmplitude * DiffusionDecayFactor * COS( WaveNumber * X )
        InuDiffusion &
          = DiffusionFluxAmplitude * DiffusionDecayFactor &
              * SIN( WaveNumber * X )

        DensityL1 = DensityL1 + Weight * ABS( DnuNumerical - DnuExact )
        DensityL2 = DensityL2 + Weight * ( DnuNumerical - DnuExact )**2
        DensityLInf = MAX( DensityLInf, ABS( DnuNumerical - DnuExact ) )
        FluxL1 = FluxL1 + Weight * ABS( InuNumerical - InuExact )
        FluxL2 = FluxL2 + Weight * ( InuNumerical - InuExact )**2
        FluxLInf = MAX( FluxLInf, ABS( InuNumerical - InuExact ) )

        DensityCosCoefficient = DensityCosCoefficient &
          + Weight * ( DnuNumerical - OpaqueDiffusionBackground ) &
              * COS( WaveNumber * X )
        FluxSinCoefficient = FluxSinCoefficient &
          + Weight * InuNumerical * SIN( WaveNumber * X )
        MeanDensityError = MeanDensityError &
          + Weight * ( DnuNumerical - OpaqueDiffusionBackground )

        CellError = CellError + Weight * ( DnuNumerical - DnuExact )

        WRITE(OutputUnit,'(9(ES24.16E3,1X))') &
          X / Kilometer, DnuNumerical, DnuExact, DnuDiffusion, &
          InuNumerical, InuExact, InuDiffusion, &
          DnuNumerical - DnuExact, InuNumerical - InuExact

      END DO

      CheckerboardAmplitude = CheckerboardAmplitude &
        + (-One)**(iX1-iX_B0(1)) * CellError

    END DO

    CLOSE( OutputUnit )

    DensityL1 = DensityL1 &
      / ( DensityAmplitude * OpaqueDiffusionDomainLength )
    DensityL2 = SQRT( DensityL2 / OpaqueDiffusionDomainLength ) &
      / DensityAmplitude
    DensityLInf = DensityLInf / DensityAmplitude
    FluxL1 = FluxL1 &
      / ( FluxAmplitude * OpaqueDiffusionDomainLength )
    FluxL2 = SQRT( FluxL2 / OpaqueDiffusionDomainLength ) &
      / FluxAmplitude
    FluxLInf = FluxLInf / FluxAmplitude

    DensityCosCoefficient &
      = Two * DensityCosCoefficient / OpaqueDiffusionDomainLength
    FluxSinCoefficient &
      = Two * FluxSinCoefficient / OpaqueDiffusionDomainLength
    MeanDensityError = MeanDensityError &
      / ( DensityAmplitude * OpaqueDiffusionDomainLength )
    CheckerboardAmplitude = ABS( CheckerboardAmplitude ) &
      / ( DensityAmplitude * OpaqueDiffusionDomainLength )

    RateDifference = ABS( ( DecayRate - DiffusionDecayRate ) / DecayRate )
    DensityModeErrorP1 = ABS &
      ( DensityCosCoefficient - DensityAmplitude * DecayFactor ) &
      / ABS( DensityAmplitude * DecayFactor )
    DensityModeErrorDiffusion = ABS &
      ( DensityCosCoefficient - DensityAmplitude * DiffusionDecayFactor ) &
      / ABS( DensityAmplitude * DiffusionDecayFactor )
    FluxModeErrorP1 = ABS &
      ( FluxSinCoefficient - FluxAmplitude * DecayFactor ) &
      / ABS( FluxAmplitude * DecayFactor )
    FluxModeErrorDiffusion = ABS &
      ( FluxSinCoefficient &
          - DiffusionFluxAmplitude * DiffusionDecayFactor ) &
      / ABS( DiffusionFluxAmplitude * DiffusionDecayFactor )

    WRITE(*,* )
    WRITE(*,'(A)') 'Manufactured opaque-diffusion error:'
    WRITE(*,'(A,ES14.6E3)') '  elapsed slow-mode e-foldings = ', &
      - DecayRate * t
    WRITE(*,'(A,ES14.6E3)') '  P1/diffusion rate difference = ', &
      RateDifference
    WRITE(*,'(A,3(1X,ES14.6E3))') '  density rel. L1/L2/Linf     = ', &
      DensityL1, DensityL2, DensityLInf
    WRITE(*,'(A,3(1X,ES14.6E3))') '  flux rel. L1/L2/Linf        = ', &
      FluxL1, FluxL2, FluxLInf
    WRITE(*,'(A,2(1X,ES14.6E3))') '  density cos coeff num/exact = ', &
      DensityCosCoefficient, DensityAmplitude * DecayFactor
    WRITE(*,'(A,2(1X,ES14.6E3))') '  density cos rel. err P1/diff= ', &
      DensityModeErrorP1, DensityModeErrorDiffusion
    WRITE(*,'(A,2(1X,ES14.6E3))') '  flux sin coeff num/exact    = ', &
      FluxSinCoefficient, FluxAmplitude * DecayFactor
    WRITE(*,'(A,2(1X,ES14.6E3))') '  flux sin rel. err P1/diff   = ', &
      FluxModeErrorP1, FluxModeErrorDiffusion
    WRITE(*,'(A,ES14.6E3)') '  mean density rel. error      = ', &
      MeanDensityError
    WRITE(*,'(A,ES14.6E3)') '  checkerboard error amplitude = ', &
      CheckerboardAmplitude
    WRITE(*,'(A)') &
      '  profile: manufactured_opaque_diffusion.dat'
    WRITE(*,* )

  END SUBROUTINE ComputeError_ManufacturedOpaqueDiffusion


  SUBROUTINE ComputeError_ManufacturedVariableOpacityDiffusion( t )

    REAL(DP), INTENT(in) :: t

    INTEGER :: iNodeX, iNodeX1, iNodeZ, iX1
    INTEGER :: OutputUnit
    REAL(DP) :: X, Weight, CellError
    REAL(DP) :: WaveNumber, Opacity
    REAL(DP) :: DensityAmplitude, FluxScale
    REAL(DP) :: DnuNumerical, DnuExact, InuNumerical, InuExact
    REAL(DP) :: DensityL1, DensityL2, DensityLInf
    REAL(DP) :: FluxL1, FluxL2, FluxLInf
    REAL(DP) :: MeanDensityError, CheckerboardAmplitude
    REAL(DP) :: FaceJumpL2, FaceJumpLInf
    LOGICAL :: ScatteringDiffusionLimitTest
    CHARACTER(64) :: ProfileFileName

    ScatteringDiffusionLimitTest &
      = TRIM( ProgramName ) == 'ScatteringDiffusionLimit'

    IF( ScatteringDiffusionLimitTest )THEN
      ProfileFileName = 'scattering_diffusion_limit.dat'
    ELSE
      ProfileFileName = 'manufactured_variable_opacity_diffusion.dat'
    END IF

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength
    DensityAmplitude &
      = OpaqueDiffusionBackground * OpaqueDiffusionAmplitude
    FluxScale = DensityAmplitude * WaveNumber &
      / ( Three * VariableOpacityDiffusionMeanOpacity )

    DensityL1 = Zero
    DensityL2 = Zero
    DensityLInf = Zero
    FluxL1 = Zero
    FluxL2 = Zero
    FluxLInf = Zero
    MeanDensityError = Zero
    CheckerboardAmplitude = Zero

    OPEN( NEWUNIT = OutputUnit, &
          FILE = TRIM( ProfileFileName ), &
          STATUS = 'REPLACE', ACTION = 'WRITE' )
    WRITE(OutputUnit,'(A)') &
      '# x_km kappa_per_km D_num D_exact I1_num I1_exact D_error I1_error'

    DO iX1 = iX_B0(1), iX_E0(1)

      CellError = Zero

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeZ = ( iNodeX - 1 ) * nDOFE + 1
        X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        Weight = WeightsX1(iNodeX1) * MeshX(1) % Width(iX1)
        Opacity = VariableOpacityDiffusionOpacity( X )

        DnuNumerical = uPR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iPR_D,1)
        InuNumerical = uPR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iPR_I1,1)
        DnuExact = OpaqueDiffusionBackground &
          + DensityAmplitude * COS( WaveNumber * X )
        InuExact = DensityAmplitude * WaveNumber * SIN( WaveNumber * X ) &
          / ( Three * Opacity )

        DensityL1 = DensityL1 + Weight * ABS( DnuNumerical - DnuExact )
        DensityL2 = DensityL2 + Weight * ( DnuNumerical - DnuExact )**2
        DensityLInf = MAX( DensityLInf, ABS( DnuNumerical - DnuExact ) )
        FluxL1 = FluxL1 + Weight * ABS( InuNumerical - InuExact )
        FluxL2 = FluxL2 + Weight * ( InuNumerical - InuExact )**2
        FluxLInf = MAX( FluxLInf, ABS( InuNumerical - InuExact ) )
        MeanDensityError = MeanDensityError &
          + Weight * ( DnuNumerical - DnuExact )
        CellError = CellError + Weight * ( DnuNumerical - DnuExact )

        WRITE(OutputUnit,'(8(ES24.16E3,1X))') &
          X / Kilometer, Opacity * Kilometer, &
          DnuNumerical, DnuExact, InuNumerical, InuExact, &
          DnuNumerical - DnuExact, InuNumerical - InuExact

      END DO

      CheckerboardAmplitude = CheckerboardAmplitude &
        + (-One)**(iX1-iX_B0(1)) * CellError

    END DO

    CLOSE( OutputUnit )

    DensityL1 = DensityL1 &
      / ( DensityAmplitude * OpaqueDiffusionDomainLength )
    DensityL2 = SQRT( DensityL2 / OpaqueDiffusionDomainLength ) &
      / DensityAmplitude
    DensityLInf = DensityLInf / DensityAmplitude
    FluxL1 = FluxL1 / ( FluxScale * OpaqueDiffusionDomainLength )
    FluxL2 = SQRT( FluxL2 / OpaqueDiffusionDomainLength ) / FluxScale
    FluxLInf = FluxLInf / FluxScale
    MeanDensityError = MeanDensityError &
      / ( DensityAmplitude * OpaqueDiffusionDomainLength )
    CheckerboardAmplitude = ABS( CheckerboardAmplitude ) &
      / ( DensityAmplitude * OpaqueDiffusionDomainLength )

    IF( ScatteringDiffusionLimitTest )THEN

      CALL ComputeScatteringDiffusionFaceJumps &
             ( DensityAmplitude, FaceJumpL2, FaceJumpLInf )

      WRITE(*,* )
      WRITE(*,'(A)') 'Variable-opacity scattering diffusion-limit error:'
      WRITE(*,'(A,ES14.6E3)') '  elapsed light-crossing times = ', &
        t / OpaqueDiffusionDomainLength
      WRITE(*,'(A,ES14.6E3)') '  nominal diffusion e-foldings= ', &
        t * ( TwoPi / OpaqueDiffusionDomainLength )**2 &
          / ( Three * VariableOpacityDiffusionMeanOpacity )

    ELSE

      FaceJumpL2   = Zero
      FaceJumpLInf = Zero

      WRITE(*,* )
      WRITE(*,'(A)') 'Manufactured variable-opacity diffusion error:'
      WRITE(*,'(A,ES14.6E3)') '  elapsed light-crossing times = ', &
        t / OpaqueDiffusionDomainLength

    END IF

    WRITE(*,'(A,3(1X,ES14.6E3))') '  density rel. L1/L2/Linf     = ', &
      DensityL1, DensityL2, DensityLInf
    WRITE(*,'(A,3(1X,ES14.6E3))') '  flux rel. L1/L2/Linf        = ', &
      FluxL1, FluxL2, FluxLInf
    WRITE(*,'(A,ES14.6E3)') '  mean density rel. error      = ', &
      MeanDensityError
    WRITE(*,'(A,ES14.6E3)') '  checkerboard error amplitude = ', &
      CheckerboardAmplitude
    IF( ScatteringDiffusionLimitTest )THEN
      WRITE(*,'(A,2(1X,ES14.6E3))') &
        '  density face-jump rel. L2/Linf = ', FaceJumpL2, FaceJumpLInf
      WRITE(*,'(A,2(1X,ES14.6E3))') &
        '  optical-depth-scaled jumps    = ', &
        VariableOpacityDiffusionMeanOpticalDepth * FaceJumpL2, &
        VariableOpacityDiffusionMeanOpticalDepth * FaceJumpLInf
      WRITE(*,'(A)') '  profile: scattering_diffusion_limit.dat'
      WRITE(*,'(A)') '  faces:   scattering_diffusion_limit_faces.dat'
    ELSE
      WRITE(*,'(A)') &
        '  profile: manufactured_variable_opacity_diffusion.dat'
    END IF
    WRITE(*,* )

  END SUBROUTINE ComputeError_ManufacturedVariableOpacityDiffusion


  SUBROUTINE ComputeScatteringDiffusionFaceJumps &
    ( DensityAmplitude, FaceJumpL2, FaceJumpLInf )

    USE ProgramHeaderModule, ONLY: &
      nDOFZ
    USE ReferenceElementModule_Lagrange, ONLY: &
      L_X1_Dn, L_X1_Up

    REAL(DP), INTENT(in)  :: DensityAmplitude
    REAL(DP), INTENT(out) :: FaceJumpL2, FaceJumpLInf

    INTEGER :: iX1, iX1_P, nFaces, OutputUnit
    REAL(DP) :: X_Face, Dnu_M, Dnu_P, Jump

    FaceJumpL2   = Zero
    FaceJumpLInf = Zero
    nFaces       = iX_E0(1) - iX_B0(1) + 1

    OPEN( NEWUNIT = OutputUnit, &
          FILE = 'scattering_diffusion_limit_faces.dat', &
          STATUS = 'REPLACE', ACTION = 'WRITE' )
    WRITE(OutputUnit,'(A)') '# x_face_km D_minus D_plus jump'

    DO iX1 = iX_B0(1), iX_E0(1)

      IF( iX1 < iX_E0(1) )THEN
        iX1_P = iX1 + 1
      ELSE
        iX1_P = iX_B0(1)
      END IF

      Dnu_M = DOT_PRODUCT &
        ( L_X1_Up(1,1:nDOFZ), &
          uCR(1:nDOFZ,iE_B0,iX1,  iX_B0(2),iX_B0(3),iCR_N,1) )
      Dnu_P = DOT_PRODUCT &
        ( L_X1_Dn(1,1:nDOFZ), &
          uCR(1:nDOFZ,iE_B0,iX1_P,iX_B0(2),iX_B0(3),iCR_N,1) )

      Jump = Dnu_P - Dnu_M
      X_Face = MeshX(1) % Center(iX1) &
        + Half * MeshX(1) % Width(iX1)

      FaceJumpL2   = FaceJumpL2 + Jump**2
      FaceJumpLInf = MAX( FaceJumpLInf, ABS( Jump ) )

      WRITE(OutputUnit,'(4(ES24.16E3,1X))') &
        X_Face / Kilometer, Dnu_M, Dnu_P, Jump

    END DO

    CLOSE( OutputUnit )

    FaceJumpL2 = SQRT( FaceJumpL2 / REAL( nFaces, DP ) ) &
      / DensityAmplitude
    FaceJumpLInf = FaceJumpLInf / DensityAmplitude

  END SUBROUTINE ComputeScatteringDiffusionFaceJumps


  SUBROUTINE ComputeError_AbsorptionDiffusionBalance( t )

    REAL(DP), INTENT(in) :: t

    INTEGER :: iNodeX, iNodeX1, iNodeZ, iX1
    INTEGER :: OutputUnit
    REAL(DP) :: X, Weight, CellSourceError
    REAL(DP) :: WaveNumber, TransportOpacity, AbsorptionOpacity
    REAL(DP) :: DensityAmplitude, FluxScale, SourceScale
    REAL(DP) :: DnuNumerical, DnuExact, Dnu_0
    REAL(DP) :: InuNumerical, InuExact
    REAL(DP) :: SourceNumerical, SourceExact
    REAL(DP) :: DensityL1, DensityL2, DensityLInf
    REAL(DP) :: FluxL1, FluxL2, FluxLInf
    REAL(DP) :: SourceL1, SourceL2, SourceLInf
    REAL(DP) :: MeanSourceError, CheckerboardAmplitude

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength
    DensityAmplitude &
      = OpaqueDiffusionBackground * OpaqueDiffusionAmplitude
    FluxScale = DensityAmplitude * WaveNumber &
      / ( Three * VariableOpacityDiffusionMeanOpacity )
    SourceScale = DensityAmplitude * WaveNumber**2 &
      / ( Three * VariableOpacityDiffusionMeanOpacity )

    DensityL1 = Zero
    DensityL2 = Zero
    DensityLInf = Zero
    FluxL1 = Zero
    FluxL2 = Zero
    FluxLInf = Zero
    SourceL1 = Zero
    SourceL2 = Zero
    SourceLInf = Zero
    MeanSourceError = Zero
    CheckerboardAmplitude = Zero

    OPEN( NEWUNIT = OutputUnit, &
          FILE = 'absorption_diffusion_balance.dat', &
          STATUS = 'REPLACE', ACTION = 'WRITE' )
    WRITE(OutputUnit,'(A)') &
      '# x_km kappa_tr_per_km kappa_a_per_km D_num D_exact I1_num ' // &
      'I1_exact source_num source_exact D_error I1_error source_error'

    DO iX1 = iX_B0(1), iX_E0(1)

      CellSourceError = Zero

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeZ = ( iNodeX - 1 ) * nDOFE + 1
        X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        Weight = WeightsX1(iNodeX1) * MeshX(1) % Width(iX1)
        TransportOpacity = VariableOpacityDiffusionOpacity( X )
        AbsorptionOpacity = AbsorptionDiffusionAbsorptionOpacity( X )

        DnuNumerical = uPR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iPR_D,1)
        InuNumerical = uPR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iPR_I1,1)
        DnuExact = OpaqueDiffusionBackground &
          + DensityAmplitude * COS( WaveNumber * X )
        InuExact = DensityAmplitude * WaveNumber * SIN( WaveNumber * X ) &
          / ( Three * TransportOpacity )
        Dnu_0 = AbsorptionDiffusionEquilibriumDensity( X )
        SourceNumerical = AbsorptionOpacity * ( Dnu_0 - DnuNumerical )
        SourceExact = VariableOpacityDiffusionDensitySource( X )

        DensityL1 = DensityL1 + Weight * ABS( DnuNumerical - DnuExact )
        DensityL2 = DensityL2 + Weight * ( DnuNumerical - DnuExact )**2
        DensityLInf = MAX( DensityLInf, ABS( DnuNumerical - DnuExact ) )
        FluxL1 = FluxL1 + Weight * ABS( InuNumerical - InuExact )
        FluxL2 = FluxL2 + Weight * ( InuNumerical - InuExact )**2
        FluxLInf = MAX( FluxLInf, ABS( InuNumerical - InuExact ) )
        SourceL1 = SourceL1 + Weight * ABS( SourceNumerical - SourceExact )
        SourceL2 = SourceL2 + Weight * ( SourceNumerical - SourceExact )**2
        SourceLInf = MAX( SourceLInf, ABS( SourceNumerical - SourceExact ) )
        MeanSourceError = MeanSourceError &
          + Weight * ( SourceNumerical - SourceExact )
        CellSourceError = CellSourceError &
          + Weight * ( SourceNumerical - SourceExact )

        WRITE(OutputUnit,'(12(ES24.16E3,1X))') &
          X / Kilometer, TransportOpacity * Kilometer, &
          AbsorptionOpacity * Kilometer, &
          DnuNumerical, DnuExact, InuNumerical, InuExact, &
          SourceNumerical, SourceExact, &
          DnuNumerical - DnuExact, InuNumerical - InuExact, &
          SourceNumerical - SourceExact

      END DO

      CheckerboardAmplitude = CheckerboardAmplitude &
        + (-One)**(iX1-iX_B0(1)) * CellSourceError

    END DO

    CLOSE( OutputUnit )

    DensityL1 = DensityL1 &
      / ( DensityAmplitude * OpaqueDiffusionDomainLength )
    DensityL2 = SQRT( DensityL2 / OpaqueDiffusionDomainLength ) &
      / DensityAmplitude
    DensityLInf = DensityLInf / DensityAmplitude
    FluxL1 = FluxL1 / ( FluxScale * OpaqueDiffusionDomainLength )
    FluxL2 = SQRT( FluxL2 / OpaqueDiffusionDomainLength ) / FluxScale
    FluxLInf = FluxLInf / FluxScale
    SourceL1 = SourceL1 / ( SourceScale * OpaqueDiffusionDomainLength )
    SourceL2 = SQRT( SourceL2 / OpaqueDiffusionDomainLength ) / SourceScale
    SourceLInf = SourceLInf / SourceScale
    MeanSourceError = MeanSourceError &
      / ( SourceScale * OpaqueDiffusionDomainLength )
    CheckerboardAmplitude = ABS( CheckerboardAmplitude ) &
      / ( SourceScale * OpaqueDiffusionDomainLength )

    WRITE(*,* )
    WRITE(*,'(A)') 'Absorption-diffusion one-step error:'
    WRITE(*,'(A,ES14.6E3)') '  elapsed light-crossing times = ', &
      t / OpaqueDiffusionDomainLength
    WRITE(*,'(A,3(1X,ES14.6E3))') '  density rel. L1/L2/Linf     = ', &
      DensityL1, DensityL2, DensityLInf
    WRITE(*,'(A,3(1X,ES14.6E3))') '  flux rel. L1/L2/Linf        = ', &
      FluxL1, FluxL2, FluxLInf
    WRITE(*,'(A,3(1X,ES14.6E3))') '  source rel. L1/L2/Linf      = ', &
      SourceL1, SourceL2, SourceLInf
    WRITE(*,'(A,ES14.6E3)') '  mean source rel. error       = ', &
      MeanSourceError
    WRITE(*,'(A,ES14.6E3)') '  source checkerboard amplitude= ', &
      CheckerboardAmplitude
    WRITE(*,'(A)') '  profile: absorption_diffusion_balance.dat'
    WRITE(*,* )

  END SUBROUTINE ComputeError_AbsorptionDiffusionBalance


  SUBROUTINE ComputeError_CoupledAbsorptionDiffusion( t )

    REAL(DP), INTENT(in) :: t

    INTEGER :: iNodeX, iNodeX1, iNodeZ, iX1
    INTEGER :: OutputUnit
    REAL(DP) :: X, Weight, CellExchangeError, CellTotalError
    REAL(DP) :: WaveNumber, DecayRate, DecayFactor
    REAL(DP) :: AbsorptionOpacity, RadiationFraction
    REAL(DP) :: EquilibriumRadiationFraction
    REAL(DP) :: RadiationNumerical, MatterNumerical, TotalNumerical
    REAL(DP) :: FluxNumerical, ExchangeNumerical
    REAL(DP) :: RadiationExact, MatterExact, TotalExact
    REAL(DP) :: FluxExact, ExchangeExact
    REAL(DP) :: TotalL1, TotalL2, TotalLInf
    REAL(DP) :: PartitionL1, PartitionL2, PartitionLInf
    REAL(DP) :: FluxL1, FluxL2, FluxLInf
    REAL(DP) :: ExchangeL1, ExchangeL2, ExchangeLInf
    REAL(DP) :: MeanTotalError, TotalCheckerboard, ExchangeCheckerboard
    REAL(DP) :: TotalAmplitude, FluxAmplitude, ExchangeAmplitude

    WaveNumber = TwoPi / OpaqueDiffusionDomainLength
    DecayRate = CoupledDiffusionDecayRate()
    DecayFactor = EXP( DecayRate * t )
    AbsorptionOpacity = CoupledDiffusionAbsorptionOpacity()
    RadiationFraction = CoupledDiffusionRadiationFraction( DecayRate )
    EquilibriumRadiationFraction &
      = CoupledDiffusionEquilibriumSlope &
          / ( One + CoupledDiffusionEquilibriumSlope )

    TotalAmplitude = CoupledDiffusionTotalAmplitude * DecayFactor
    FluxAmplitude = - DecayRate * TotalAmplitude / WaveNumber
    ExchangeAmplitude &
      = - DecayRate * ( One - RadiationFraction ) * TotalAmplitude

    TotalL1 = Zero
    TotalL2 = Zero
    TotalLInf = Zero
    PartitionL1 = Zero
    PartitionL2 = Zero
    PartitionLInf = Zero
    FluxL1 = Zero
    FluxL2 = Zero
    FluxLInf = Zero
    ExchangeL1 = Zero
    ExchangeL2 = Zero
    ExchangeLInf = Zero
    MeanTotalError = Zero
    TotalCheckerboard = Zero
    ExchangeCheckerboard = Zero

    OPEN( NEWUNIT = OutputUnit, &
          FILE = 'coupled_absorption_diffusion.dat', &
          STATUS = 'REPLACE', ACTION = 'WRITE' )
    WRITE(OutputUnit,'(A)') &
      '# x_km D_num D_exact M_num M_exact U_num U_exact I_num I_exact ' // &
      'Q_num Q_exact U_error D_error I_error Q_error'

    DO iX1 = iX_B0(1), iX_E0(1)

      CellTotalError = Zero
      CellExchangeError = Zero

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeZ = ( iNodeX - 1 ) * nDOFE + 1
        X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        Weight = WeightsX1(iNodeX1) * MeshX(1) % Width(iX1)

        RadiationNumerical &
          = uPR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iPR_D,1)
        MatterNumerical &
          = uPR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iPR_D,2)
        FluxNumerical &
          = uPR(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iPR_I1,1)
        TotalNumerical = RadiationNumerical + MatterNumerical
        ExchangeNumerical = AbsorptionOpacity &
          * ( CoupledDiffusionEquilibriumSlope * MatterNumerical &
              - RadiationNumerical )

        TotalExact = CoupledDiffusionTotalBackground &
          + TotalAmplitude * COS( WaveNumber * X )
        RadiationExact = EquilibriumRadiationFraction &
                * CoupledDiffusionTotalBackground &
              + RadiationFraction * TotalAmplitude * COS( WaveNumber * X )
        MatterExact = TotalExact - RadiationExact
        FluxExact = FluxAmplitude * SIN( WaveNumber * X )
        ExchangeExact = ExchangeAmplitude * COS( WaveNumber * X )

        TotalL1 = TotalL1 + Weight * ABS( TotalNumerical - TotalExact )
        TotalL2 = TotalL2 + Weight * ( TotalNumerical - TotalExact )**2
        TotalLInf = MAX( TotalLInf, ABS( TotalNumerical - TotalExact ) )
        PartitionL1 = PartitionL1 &
          + Weight * ABS( RadiationNumerical - RadiationExact )
        PartitionL2 = PartitionL2 &
          + Weight * ( RadiationNumerical - RadiationExact )**2
        PartitionLInf = MAX &
          ( PartitionLInf, ABS( RadiationNumerical - RadiationExact ) )
        FluxL1 = FluxL1 + Weight * ABS( FluxNumerical - FluxExact )
        FluxL2 = FluxL2 + Weight * ( FluxNumerical - FluxExact )**2
        FluxLInf = MAX( FluxLInf, ABS( FluxNumerical - FluxExact ) )
        ExchangeL1 = ExchangeL1 &
          + Weight * ABS( ExchangeNumerical - ExchangeExact )
        ExchangeL2 = ExchangeL2 &
          + Weight * ( ExchangeNumerical - ExchangeExact )**2
        ExchangeLInf = MAX &
          ( ExchangeLInf, ABS( ExchangeNumerical - ExchangeExact ) )
        MeanTotalError = MeanTotalError &
          + Weight * ( TotalNumerical - TotalExact )
        CellTotalError = CellTotalError &
          + Weight * ( TotalNumerical - TotalExact )
        CellExchangeError = CellExchangeError &
          + Weight * ( ExchangeNumerical - ExchangeExact )

        WRITE(OutputUnit,'(15(ES24.16E3,1X))') &
          X / Kilometer, RadiationNumerical, RadiationExact, &
          MatterNumerical, MatterExact, TotalNumerical, TotalExact, &
          FluxNumerical, FluxExact, ExchangeNumerical, ExchangeExact, &
          TotalNumerical - TotalExact, &
          RadiationNumerical - RadiationExact, &
          FluxNumerical - FluxExact, ExchangeNumerical - ExchangeExact

      END DO

      TotalCheckerboard = TotalCheckerboard &
        + (-One)**(iX1-iX_B0(1)) * CellTotalError
      ExchangeCheckerboard = ExchangeCheckerboard &
        + (-One)**(iX1-iX_B0(1)) * CellExchangeError

    END DO

    CLOSE( OutputUnit )

    TotalL1 = TotalL1 / ( TotalAmplitude * OpaqueDiffusionDomainLength )
    TotalL2 = SQRT( TotalL2 / OpaqueDiffusionDomainLength ) / TotalAmplitude
    TotalLInf = TotalLInf / TotalAmplitude
    PartitionL1 &
      = PartitionL1 / ( RadiationFraction * TotalAmplitude &
                          * OpaqueDiffusionDomainLength )
    PartitionL2 = SQRT( PartitionL2 / OpaqueDiffusionDomainLength ) &
      / ( RadiationFraction * TotalAmplitude )
    PartitionLInf = PartitionLInf / ( RadiationFraction * TotalAmplitude )
    FluxL1 = FluxL1 / ( FluxAmplitude * OpaqueDiffusionDomainLength )
    FluxL2 = SQRT( FluxL2 / OpaqueDiffusionDomainLength ) / FluxAmplitude
    FluxLInf = FluxLInf / FluxAmplitude
    ExchangeL1 &
      = ExchangeL1 / ( ExchangeAmplitude * OpaqueDiffusionDomainLength )
    ExchangeL2 = SQRT( ExchangeL2 / OpaqueDiffusionDomainLength ) &
      / ExchangeAmplitude
    ExchangeLInf = ExchangeLInf / ExchangeAmplitude
    MeanTotalError = MeanTotalError &
      / ( CoupledDiffusionTotalBackground * OpaqueDiffusionDomainLength )
    TotalCheckerboard = ABS( TotalCheckerboard ) &
      / ( TotalAmplitude * OpaqueDiffusionDomainLength )
    ExchangeCheckerboard = ABS( ExchangeCheckerboard ) &
      / ( ExchangeAmplitude * OpaqueDiffusionDomainLength )

    WRITE(*,* )
    WRITE(*,'(A)') 'Coupled absorption-diffusion error:'
    WRITE(*,'(A,ES14.6E3)') '  elapsed slow-mode e-foldings = ', &
      - DecayRate * t
    WRITE(*,'(A,3(1X,ES14.6E3))') '  total rel. L1/L2/Linf       = ', &
      TotalL1, TotalL2, TotalLInf
    WRITE(*,'(A,3(1X,ES14.6E3))') '  partition rel. L1/L2/Linf   = ', &
      PartitionL1, PartitionL2, PartitionLInf
    WRITE(*,'(A,3(1X,ES14.6E3))') '  flux rel. L1/L2/Linf        = ', &
      FluxL1, FluxL2, FluxLInf
    WRITE(*,'(A,3(1X,ES14.6E3))') '  exchange rel. L1/L2/Linf    = ', &
      ExchangeL1, ExchangeL2, ExchangeLInf
    WRITE(*,'(A,ES14.6E3)') '  mean total conservation err.= ', &
      MeanTotalError
    WRITE(*,'(A,ES14.6E3)') '  total checkerboard amplitude= ', &
      TotalCheckerboard
    WRITE(*,'(A,ES14.6E3)') '  exch. checkerboard amplitude= ', &
      ExchangeCheckerboard
    WRITE(*,'(A)') '  profile: coupled_absorption_diffusion.dat'
    WRITE(*,* )

  END SUBROUTINE ComputeError_CoupledAbsorptionDiffusion


END MODULE InitializationModule_Neutrinos
