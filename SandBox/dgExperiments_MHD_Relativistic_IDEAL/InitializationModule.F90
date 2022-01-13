MODULE InitializationModule

  USE KindModule, ONLY: &
    DP,       &
    Zero,     &
    Half,     &
    One,      &
    Two,      &
    Three,    &
    Four,     &
    Pi,       &
    TwoPi,    &
    FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nNodesX,     &
    nDimsX,      &
    nDOFX,       &
    iX_B0,       &
    iX_B1,       &
    iX_E0,       &
    iX_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3
  USE MagnetofluidFieldsModule, ONLY: &
    nPM,     &
    uPM,     &
    iPM_D,   &
    iPM_V1,  &
    iPM_V2,  &
    iPM_V3,  &
    iPM_E,   &
    iPM_Ne,  &
    iPM_B1,  &
    iPM_B2,  &
    iPM_B3,  &
    iPM_Chi, &
    uCM,     &
    iCM_D,   &
    iCM_S1,  &
    iCM_S2,  &
    iCM_S3,  &
    iCM_E,   &
    iCM_Ne,  &
    iCM_B1,  &
    iCM_B2,  &
    iCM_B3,  &
    iCM_Chi, &
    uAM,     &
    iAM_P
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL, &
    ComputePressureFromPrimitive_IDEAL
  USE MHD_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_MHD_Relativistic
  USE UnitsModule, ONLY: &
    Gauss
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE PolynomialBasisModule_Lagrange, ONLY: &
    LagrangeP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Relativistic_MHD


CONTAINS


  SUBROUTINE InitializeFields_Relativistic_MHD &
               ( AdvectionProfile_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: AdvectionProfile_Option
    CHARACTER(LEN=64) :: AdvectionProfile = 'SineWave'

    uPM(:,:,:,:,iPM_Ne) = Zero

    IF( PRESENT( AdvectionProfile_Option ) ) &
      AdvectionProfile = TRIM( AdvectionProfile_Option )

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Advection' )

        CALL InitializeFields_Advection &
               ( TRIM( AdvectionProfile ) )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE InitializeFields_Relativistic_MHD


  SUBROUTINE InitializeFields_Advection( AdvectionProfile )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          ! Sine wave advection problem with comparable parameters
          ! to the Top Hat advection problem from Evans and Hawley (1988).

          CASE( 'SineWave' )

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.1_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0001_DP * SIN( TwoPi * X1 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

           ! Variant of Top Hat advection problem from Evans and 
           ! Hawley (1988).

           CASE( 'TopHat' )

            IF( X1 .GT. 0.05 .AND. X1 .LT. 0.55 )THEN

              uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.25_DP
 
            ELSE

              uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
 
            END IF

              uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = 1.0_DP 
              uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.1_DP
              uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
              uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
              uAM(iNodeX,iX1,iX2,iX3,iAM_P)  = 1.0_DP
              uPM(iNodeX,iX1,iX2,iX3,iPM_E)  &
                = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
              uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.25_DP
              uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP
              uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  SineWave'
            WRITE(*,'(A)') '  TopHat'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(:,iX1,iX2,iX3,iPM_D ), uPM(:,iX1,iX2,iX3,iPM_V1),  &
               uPM(:,iX1,iX2,iX3,iPM_V2), uPM(:,iX1,iX2,iX3,iPM_V3),  &
               uPM(:,iX1,iX2,iX3,iPM_E ), uPM(:,iX1,iX2,iX3,iPM_Ne),  &
               uPM(:,iX1,iX2,iX3,iPM_B1), uPM(:,iX1,iX2,iX3,iPM_B2),  &
               uPM(:,iX1,iX2,iX3,iPM_B3), uPM(:,iX1,iX2,iX3,iPM_Chi), &           
               uCM(:,iX1,iX2,iX3,iCM_D ), uCM(:,iX1,iX2,iX3,iCM_S1),  &
               uCM(:,iX1,iX2,iX3,iCM_S2), uCM(:,iX1,iX2,iX3,iCM_S3),  &
               uCM(:,iX1,iX2,iX3,iCM_E ), uCM(:,iX1,iX2,iX3,iCM_Ne),  &
               uCM(:,iX1,iX2,iX3,iCM_B1), uCM(:,iX1,iX2,iX3,iCM_B2),  &
               uCM(:,iX1,iX2,iX3,iCM_B3), uCM(:,iX1,iX2,iX3,iCM_Chi), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(:,iX1,iX2,iX3,iGF_Alpha   ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_1  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_2  ), &
               uGF(:,iX1,iX2,iX3,iGF_Beta_3  ), &
               uAM(:,iX1,iX2,iX3,iAM_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Advection


END MODULE InitializationModule
