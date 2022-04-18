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
  USE EquationOfStateModule, ONLY : &
      ComputePressureFromPrimitive
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
               ( AdvectionProfile_Option, SmoothProfile_Option, &
                 ConstantDensity_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: AdvectionProfile_Option
    LOGICAL,          INTENT(in), OPTIONAL :: SmoothProfile_Option
    LOGICAL,          INTENT(in), OPTIONAL :: ConstantDensity_Option

    CHARACTER(LEN=64) :: AdvectionProfile = 'MagneticSineWaveX1'
    LOGICAL           :: SmoothProfile
    LOGICAL           :: ConstantDensity = .TRUE.

    uPM(:,:,:,:,iPM_Ne) = Zero

    IF( PRESENT( AdvectionProfile_Option ) ) &
      AdvectionProfile = TRIM( AdvectionProfile_Option )

    IF( PRESENT( SmoothProfile_Option ) ) &
      SmoothProfile = SmoothProfile_Option

    IF( PRESENT( ConstantDensity_Option ) ) &
      ConstantDensity = ConstantDensity_Option

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Advection' )

        CALL InitializeFields_Advection &
               ( TRIM( AdvectionProfile ) )

      CASE( 'Advection2D' )

        CALL InitializeFields_Advection2D &
               ( TRIM( AdvectionProfile ) )

      CASE( 'Cleaning1D' )

        CALL InitializeFields_Cleaning1D &
               ( SmoothProfile )

      CASE( 'Cleaning2D' )

        CALL InitializeFields_Cleaning2D &
               ( ConstantDensity )

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
    REAL(DP) :: Eta, h, P, VA, W, k, VdotB, &
                V1_Transport, V2_Transport, V3_Transport

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

          ! Sine wave advection problem to test hydro portion of code.

          CASE( 'HydroSineWaveX1' )

            !PRINT*
            !PRINT*, 'In cell: ', iX1, iX2, iX3
            !PRINT*, 'In node: ', iNodeX
            !PRINT*, 'X1: ', X1

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One + 0.1_DP * SIN ( TwoPi * X1 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.1_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            CALL ComputePressureFromPrimitive &
                   ( uPM(iNodeX,iX1,iX2,iX3,iPM_D), &
                     uPM(iNodeX,iX1,iX2,iX3,iPM_E), &
                     0.0_DP, P )
              
            !PRINT*, 'PM_D: ',  uPM(iNodeX,iX1,iX2,iX3,iPM_D)
            !PRINT*, 'W: ', One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 )
            !PRINT*, 'h: ', One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                     * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) )
            !PRINT*, 'mu: ', One / ( ( One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 ) ) & 
            !                        * ( One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                                  * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) ) ) )
            !PRINT*

          CASE( 'MagneticSineWaveX1' )

            !PRINT*
            !PRINT*, 'In cell: ', iX1, iX2, iX3
            !PRINT*, 'In node: ', iNodeX
            !PRINT*, 'X1: ', X1

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.1_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0001_DP * SIN( TwoPi * X1 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0001_DP * COS( TwoPi * X1 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            !CALL ComputePressureFromPrimitive &
            !       ( uPM(iNodeX,iX1,iX2,iX3,iPM_D), &
            !         uPM(iNodeX,iX1,iX2,iX3,iPM_E), &
            !         0.0_DP, P )
    
            !PRINT*, 'PM_D: ',  uPM(iNodeX,iX1,iX2,iX3,iPM_D)
            !PRINT*, 'W: ', One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 )
            !PRINT*, 'h: ', One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                     * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) )
            !PRINT*, 'mu: ', One / ( ( One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 ) ) & 
            !                        * ( One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                                  * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) ) ) )
           !PRINT*

          ! Variant of Top Hat advection problem from Evans and 
          ! Hawley (1988).

          CASE( 'TopHat' )

          ! Circularly polarized Alfven wave with the exact solution
          ! from Del Zanna et al. (2007) and Mattia and Mignone (2022).

          CASE( 'CPAlfven' )

            !PRINT*
            !PRINT*, 'In cell: ', iX1, iX2, iX3
            !PRINT*, 'In node: ', iNodeX
            !PRINT*, 'X1: ', X1

            Eta = One
            k   = One
            h   = One + Gamma_IDEAL / ( Gamma_IDEAL - One )
            VA  = SQRT( ( Two / ( h + ( One + Eta**2 ) ) ) &
                        * ( One / ( One + SQRT( One - ( Two * Eta / ( h + ( One + Eta**2 ) ) )**2 ) ) ) )

            !PRINT*, 'h: ', h

            !PRINT*, 'VA: ', VA

            !PRINT*, 'One Period: ', Two * Pi / VA

            W = One / SQRT( One - VA**2 * Eta**2 )

            !PRINT*, 'Lorentz Factor: ', W

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = -VA * Eta * COS( k * X1 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = -VA * Eta * SIN( k * X1 )
            uAM(iNodeX,iX1,iX2,iX3,iAM_P)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )

            VdotB = uPM(iNodeX,iX1,iX2,iX3,iPM_V2) * Eta * COS( k * X1 ) &
                      + uPM(iNodeX,iX1,iX2,iX3,iPM_V3) * Eta * SIN( k * X1 ) 
          
            !PRINT*, 'VdotB: ', VdotB

            !PRINT*, 'Alpha:  ', uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha)
            !PRINT*, 'Beta 2: ', uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2)
            !PRINT*, 'Beta 3: ', uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3)

            V1_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V1) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_1) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) ) 
 
            V2_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V2) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) ) 

            V3_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V3) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) ) 

            !PRINT*, 'V1_Transport: ', V1_Transport
            !PRINT*, 'V2_Transport: ', V2_Transport
            !PRINT*, 'V3_Transport: ', V3_Transport

            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) &
              = W * VdotB * V1_Transport + ( One / W )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) &
              = W * VdotB * V2_Transport + Eta * COS( k * X1 ) / W 
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) &
              = W * VdotB * V3_Transport + Eta * SIN( k * X1 ) / W 
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            !PRINT*

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  HydroSineWaveX1'
            WRITE(*,'(A)') '  MagneticSineWaveX1'
            WRITE(*,'(A)') '  TopHat'
            WRITE(*,'(A)') '  CPAlfven'
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

    PRINT*, 'Finished initialization.'

  END SUBROUTINE InitializeFields_Advection


  SUBROUTINE InitializeFields_Advection2D( AdvectionProfile )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP) :: Eta, h, P, VA, W, k, VdotB, &
                V1_Transport, V2_Transport, V3_Transport

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          ! Sine wave advection problems to test hydro portion of code.

         CASE( 'HydroSineWaveX2' )

            !PRINT*
            !PRINT*, 'In cell: ', iX1, iX2, iX3
            !PRINT*, 'In node: ', iNodeX
            !PRINT*, 'X1: ', X1

            uPM(iNodeX,iX1,iX2,iX3,iPM_D) & 
              = One + 0.1_DP * SIN ( TwoPi * X2 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.1_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            CALL ComputePressureFromPrimitive &
                   ( uPM(iNodeX,iX1,iX2,iX3,iPM_D), &
                     uPM(iNodeX,iX1,iX2,iX3,iPM_E), &
                     0.0_DP, P )
              
            !PRINT*, 'PM_D: ',  uPM(iNodeX,iX1,iX2,iX3,iPM_D)
            !PRINT*, 'W: ', One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 )
            !PRINT*, 'h: ', One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                     * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) )
            !PRINT*, 'mu: ', One / ( ( One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 ) ) & 
            !                        * ( One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                                  * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) ) ) )
            !PRINT*

          !PRINT*

          CASE( 'HydroSineWaveX1X2' )

            !PRINT*
            !PRINT*, 'In cell: ', iX1, iX2, iX3
            !PRINT*, 'In node: ', iNodeX
            !PRINT*, 'X1: ', X1

            uPM(iNodeX,iX1,iX2,iX3,iPM_D) & 
              = One + 0.1_DP * SIN ( SQRT( Two ) * TwoPi * ( X1 + X2 ) )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.1_DP / SQRT( Two )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.1_DP / SQRT( Two )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            CALL ComputePressureFromPrimitive &
                   ( uPM(iNodeX,iX1,iX2,iX3,iPM_D), &
                     uPM(iNodeX,iX1,iX2,iX3,iPM_E), &
                     0.0_DP, P )
              
            !PRINT*, 'PM_D: ',  uPM(iNodeX,iX1,iX2,iX3,iPM_D)
            !PRINT*, 'W: ', One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 )
            !PRINT*, 'h: ', One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                     * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) )
            !PRINT*, 'mu: ', One / ( ( One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 ) ) & 
            !                        * ( One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                                  * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) ) ) )
            !PRINT*

          !PRINT*

         CASE( 'MagneticSineWaveX2' )

            !PRINT*
            !PRINT*, 'In cell: ', iX1, iX2, iX3
            !PRINT*, 'In node: ', iNodeX
            !PRINT*, 'X1: ', X1

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.1_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0001_DP * SIN( TwoPi * X2 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0001_DP * COS( TwoPi * X2 )
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            !CALL ComputePressureFromPrimitive &
            !       ( uPM(iNodeX,iX1,iX2,iX3,iPM_D), &
            !         uPM(iNodeX,iX1,iX2,iX3,iPM_E), &
            !         0.0_DP, P )
    
            !PRINT*, 'PM_D: ',  uPM(iNodeX,iX1,iX2,iX3,iPM_D)
            !PRINT*, 'W: ', One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 )
            !PRINT*, 'h: ', One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                     * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) )
            !PRINT*, 'mu: ', One / ( ( One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 ) ) & 
            !                        * ( One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                                  * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) ) ) )
           !PRINT*

          CASE( 'MagneticSineWaveX1X2' )

            !PRINT*
            !PRINT*, 'In cell: ', iX1, iX2, iX3
            !PRINT*, 'In node: ', iNodeX
            !PRINT*, 'X1: ', X1

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.1_DP / SQRT( Two )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.1_DP / SQRT( Two )
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0001_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0001_DP * SIN( SQRT( Two ) * TwoPi * ( X1 + X2 ) )
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            !CALL ComputePressureFromPrimitive &
            !       ( uPM(iNodeX,iX1,iX2,iX3,iPM_D), &
            !         uPM(iNodeX,iX1,iX2,iX3,iPM_E), &
            !         0.0_DP, P )
    
            !PRINT*, 'PM_D: ',  uPM(iNodeX,iX1,iX2,iX3,iPM_D)
            !PRINT*, 'W: ', One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 )
            !PRINT*, 'h: ', One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                     * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) )
            !PRINT*, 'mu: ', One / ( ( One / SQRT( One - uPM(iNodeX,iX1,iX2,iX3,iPM_V1)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V2)**2 &
            !                               - uPM(iNodeX,iX1,iX2,iX3,iPM_V3)**2 ) ) & 
            !                        * ( One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
            !                                  * ( P / uPM(iNodeX,iX1,iX2,iX3,iPM_D) ) ) )
           !PRINT*

         ! Loop advection problem from Mosta et al. (2013)

          CASE( 'LoopAdvection' )

            IF( SQRT( X1**2 + X2**2 ) .LE. 0.3 )THEN

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = One / Two
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = One / 24.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = Three
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )


            W = One / SQRT( One - ( ( One / Two )**2 + ( One / 24.0_DP )**2 ) )

            VdotB = uPM(iNodeX,iX1,iX2,iX3,iPM_V1) &
                      * ( -1.0d-3 * ( X2 / SQRT( X1**2 + X2**2 ) ) ) &
                      + uPM(iNodeX,iX1,iX2,iX3,iPM_V2) &
                          * ( 1.0d-3 * ( X1 / SQRT( X1**2 + X2**2 ) ) )

            V1_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V1) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_1) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) ) 
 
            V2_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V2) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) ) 

            V3_Transport = uPM(iNodeX,iX1,iX2,iX3,iPM_V3) &
                           - ( uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3) &
                               / uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha ) ) 

            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = W * VdotB * V1_Transport &
                                             + ( -1.0d-3 * ( X2 / SQRT( X1**2 + X2**2 ) ) ) / W
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = W * VdotB * V2_Transport & 
                                             + (  1.0d-3 * ( X1 / SQRT( X1**2 + X2**2 ) ) ) / W 
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = W * VdotB * V3_Transport
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            ELSE

            uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
            uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = One / Two
            uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = One / 24.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
            uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = Three
            uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
              = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP
            uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

            END IF

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  HydroSineWaveX2'
            WRITE(*,'(A)') '  HydroSineWaveX1X2'
            WRITE(*,'(A)') '  MagneticSineWaveX2'
            WRITE(*,'(A)') '  MagneticSineWaveX1X2'
            WRITE(*,'(A)') '  LoopAdvection'
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

    PRINT*, 'Finished initialization.'

  END SUBROUTINE InitializeFields_Advection2D


  SUBROUTINE InitializeFields_Cleaning1D( SmoothProfile )

    LOGICAL, INTENT(in) :: SmoothProfile

    ! 1D divergence cleaning test from Section 5.1 of Derigs et al. (2018)
    ! with option to use only the smooth part of the initial condition.

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        uPM(iNodeX,iX1,iX2,iX3,iPM_D)  = One
        uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.0_DP
        uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
        uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
        uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = One
        uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
          = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )

        IF( SmoothProfile )THEN

          IF( ( X1 >= -1.0_DP ) .AND. ( X1 <= -0.6_DP ) )THEN
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
          ELSE IF( ( X1 >= 0.6_DP ) .AND. ( X1 <= 1.0_DP ) )THEN
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
          ELSE
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = EXP( -( X1 / 0.11_DP )**2 / Two )
          END IF

        ELSE

          IF( ( X1 > -0.8_DP ) .AND. ( X1 <= -0.6_DP ) )THEN
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = -Two * ( X1 + 0.8_DP )
          ELSE IF( ( X1 > -0.6_DP ) .AND. ( X1 <= 0.6_DP ) )THEN
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = EXP( -( X1 / 0.11_DP )**2 / Two )
          ELSE IF( ( X1 > 0.6_DP ) .AND. ( X1 <= 1.0_DP ) )THEN
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.5_DP
          ELSE
            uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
          END IF

        END IF

        uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
        uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = 0.0_DP         
        uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

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

    PRINT*, 'Finished initialization.'

  END SUBROUTINE InitializeFields_Cleaning1D


  SUBROUTINE InitializeFields_Cleaning2D( ConstantDensity )

    ! 2D divergence cleaning test from Section 5.2 of
    ! Derigs et al. (2018) with option to use 
    ! constant density.

    LOGICAL :: ConstantDensity

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    
    REAL(DP) :: R, R0

    R0 = One / SQRT( 8.0_DP )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        R = SQRT( X1**2 + X2**2 )

        IF( ConstantDensity )THEN

          uPM(iNodeX,iX1,iX2,iX3,iPM_D) = 1.0_DP

        ELSE

          IF( X1 .LE. 0.5_DP )THEN
            uPM(iNodeX,iX1,iX2,iX3,iPM_D) = 1.0_DP
          ELSE
            uPM(iNodeX,iX1,iX2,iX3,iPM_D) = 2.0_DP
          END IF

        END IF

        uPM(iNodeX,iX1,iX2,iX3,iPM_V1) = 0.0_DP
        uPM(iNodeX,iX1,iX2,iX3,iPM_V2) = 0.0_DP
        uPM(iNodeX,iX1,iX2,iX3,iPM_V3) = 0.0_DP
        uAM(iNodeX,iX1,iX2,iX3,iAM_P ) = 6.0_DP
        uPM(iNodeX,iX1,iX2,iX3,iPM_E )  &
          = uAM(iNodeX,iX1,iX2,iX3,iAM_P) / ( Gamma_IDEAL - One )

        IF( R .LE. R0 )THEN

          uPM(iNodeX,iX1,iX2,iX3,iPM_B1) &
            = ( One / SQRT( Four * Pi ) ) &
              * ( ( R / R0 )**8 - Two * ( R / R0 )**4 + One ) 
    
        ELSE

          uPM(iNodeX,iX1,iX2,iX3,iPM_B1) = 0.0_DP
   
        END IF

        uPM(iNodeX,iX1,iX2,iX3,iPM_B2) = 0.0_DP
        uPM(iNodeX,iX1,iX2,iX3,iPM_B3) = One / SQRT( Four * Pi )         
        uPM(iNodeX,iX1,iX2,iX3,iPM_Chi) = 0.0_DP

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

    PRINT*, 'Finished initialization.'

  END SUBROUTINE InitializeFields_Cleaning2D


END MODULE InitializationModule
