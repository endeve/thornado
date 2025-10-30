MODULE InitializationModule_Neutrinos

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, SqrtTiny
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer, &
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
    NodeNumberTableX
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshE, MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, iCF_Nm, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, iPF_Nm, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_Ym, iAF_S, iAF_E, &
    iAF_Me, iAF_Mm, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, &
    iAF_Xa, iAF_Xh, iAF_Gm
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, iNuE, iNuE_Bar, iNuM, iNuM_Bar, &
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
    ComputeNeutrinoOpacities_EC

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields
  PUBLIC :: ComputeError

CONTAINS


  SUBROUTINE InitializeFields( ProfileName, FileName )

    CHARACTER(LEN=*), INTENT(in) :: ProfileName
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: FileName

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE( TRIM( ProgramName ) )

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

    REAL(DP) :: D_0, T_0, Ye_0, Ym_0, E_0, V_u_1, V_u_2, V_u_3
    REAL(DP) :: Mu_0

    REAL(DP), DIMENSION(iE_E0*nDOFE) :: tmp_nu

    INTEGER  :: tmp_int, nE_G


    IF( PRESENT (FileName) ) THEN

      OPEN( UNIT=17,file=FileName,&
            ACCESS='sequential', FORM='formatted', ACTION='read')

      READ(17,*) D_0, T_0, Ye_0, Ym_0, E_0, V_u_1, V_u_2, V_u_3
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
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = Ye_0
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ym) = Ym_0

          CALL ComputeThermodynamicStates_Primitive_TABLE &
                 ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Ym), &
                   uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                   uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                   uPF(iNodeX,iX1,iX2,iX3,iPF_Nm) )

          CALL ApplyEquationOfState_TABLE &
                 ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Ym), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_S ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Me), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Mm), &
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
                   uPF(iNodeX,iX1,iX2,iX3,iPF_Nm), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_Nm), &
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
      Ye_0  = 0.1347_DP
      Ym_0  = Zero
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
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = Ye_0
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ym) = Ym_0

          CALL ComputeThermodynamicStates_Primitive_TABLE &
                 ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Ym), &
                   uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                   uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                   uPF(iNodeX,iX1,iX2,iX3,iPF_Nm) )

          CALL ApplyEquationOfState_TABLE &
                 ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Ym), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_S ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Me), &
                   uAF(iNodeX,iX1,iX2,iX3,iAF_Mm), &
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
                   uPF(iNodeX,iX1,iX2,iX3,iPF_Nm), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                   uCF(iNodeX,iX1,iX2,iX3,iCF_Nm), &
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
    REAL(DP), ALLOCATABLE :: R_P(:), D_P(:), T_P(:), Ye_P(:), Ym_P(:)
    REAL(DP), ALLOCATABLE :: E_Nu(:), R_Nu(:,:), Chi(:,:,:), fEQ(:,:,:)
    REAL(DP), ALLOCATABLE :: D_Nu_P(:,:,:), I1_Nu_P(:,:,:)

    WRITE(*,*)
    WRITE(*,'(A6,A,A)') '', &
      'Initializing from Profile: ', TRIM( ProfileName )

    CALL ReadFluidProfile( ProfileName, R_P, D_P, T_P, Ye_P, Ym_P )

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
          = Interpolate1D_Linear( R, R_P(iR), R_P(iR+1), Ye_P(iR), Ye_P(iR+1) )

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ym) &
          = Interpolate1D_Linear( R, R_P(iR), R_P(iR+1), Ym_P(iR), Ym_P(iR+1) )

        CALL ComputeThermodynamicStates_Primitive_TABLE &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ym), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Nm) )

        CALL ApplyEquationOfState_TABLE &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ym), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_S ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Me), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Mm), &
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
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Nm), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Nm), &
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

    Chi = Zero

    ! --- Neutrino Energies ---

    DO iE = iE_B0, iE_E0
    DO iNodeE = 1, nDOFE

      E_Nu((iE-1)*nDOFE+iNodeE) = NodeCoordinate( MeshE, iE, iNodeE )

    END DO
    END DO

    ! --- Neutrino Absorption Opacities and Equilibrium Distributions ---

    CALL ComputeEquilibriumDistributions_DG &
           ( 1, nE, 1, nSpecies, 1, nR, E_Nu, D_P, T_P, Ye_P, Ym_P, fEQ )

    CALL ComputeNeutrinoOpacities_EC &
           ( 1, nE, 1, nSpecies, 1, nR, E_Nu, D_P, T_P, Ye_P, Ym_P, fEQ, Chi )

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

    DEALLOCATE( R_P, D_P, T_P, Ye_P, Ym_P )
    DEALLOCATE( E_Nu, R_Nu, Chi, fEQ )
    DEALLOCATE( D_Nu_P, I1_Nu_P )

  END SUBROUTINE InitializeFields_DeleptonizationWave1D


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
        uAF(iNodeX,iX1,iX2,iX3,iAF_Ym) = Zero

        CALL ComputeThermodynamicStates_Primitive_TABLE &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ym), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Nm) )

        CALL ApplyEquationOfState_TABLE &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ym), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_S ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Me), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Mn), &
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
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Nm), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Nm), &
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

        IF ( iS == iNuE ) THEN

          Mnu = + uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn)

        ELSE IF ( iS == iNuE_Bar ) THEN

          Mnu = - uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                - uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                + uAF(iNodeX,iX1,iX2,iX3,iAF_Mn)

        ELSE IF ( iS == iNuM ) THEN

          Mnu = + uAF(iNodeX,iX1,iX2,iX3,iAF_Mm) &
                + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn)

        ELSE IF ( iS == iNuM_Bar ) THEN

          Mnu = - uAF(iNodeX,iX1,iX2,iX3,iAF_Mm) &
                - uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                + uAF(iNodeX,iX1,iX2,iX3,iAF_Mn)

        ELSE

          Mnu = Zero

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


  SUBROUTINE ReadFluidProfile( FileName, R, D, T, Ye, Ym )

    CHARACTER(*),          INTENT(in)    :: FileName
    REAL(DP), ALLOCATABLE, INTENT(inout) :: R(:), D(:), T(:), Ye(:), Ym(:)

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

    ALLOCATE( R(nPoints), D(nPoints), T(nPoints), Ye(nPoints), Ym(nPoints) )

    R  = Data(:,1) * Kilometer
    D  = Data(:,2) * Gram / Centimeter**3
    T  = Data(:,3) * MeV
    Ye = Data(:,4)
    Ym = Zero

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

    END SELECT

  END SUBROUTINE ComputeError


  SUBROUTINE ComputeError_Relaxation

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeE, iNodeX, iNodeZ
    INTEGER  :: nE, nX(3), nE_P, nX_P, iE_P, iX_P
    REAL(DP), ALLOCATABLE :: E_P(:), D_P(:), T_P(:), Ye_P(:), Ym_P(:), f0_P(:,:,:)
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
    ALLOCATE( Ye_P(nX_P) )
    ALLOCATE( Ym_P(nX_P) )

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
      Ye_P(iX_P) = uAF(iNodeX,iX1,iX2,iX3,iAF_Ye)
      Ym_P(iX_P) = uAF(iNodeX,iX1,iX2,iX3,iAF_Ym)
    END DO
    END DO
    END DO
    END DO

    ALLOCATE( f0_P(nE_P,nSpecies,nX_P) )

    CALL ComputeEquilibriumDistributions_DG &
           ( 1, nE_P, 1, nSpecies, 1, nX_P, E_P, D_P, T_P, Ye_P, Ym_P, f0_P )

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


END MODULE InitializationModule_Neutrinos
