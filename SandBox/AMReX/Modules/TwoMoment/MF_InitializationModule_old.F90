MODULE MF_InitializationModule
  ! --- AMReX Modules ---

  USE amrex_fort_module,      ONLY: &
    AR => amrex_real
  USE amrex_box_module,       ONLY: &
    amrex_box
  USE amrex_multifab_module,  ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,  ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse,       &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE EquationOfStateModule,   ONLY: &
    ComputeAuxiliary_Fluid, &
    ApplyEquationOfState, &
    ComputePressureFromPrimitive
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeTemperatureFromPressure_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ApplyEquationOfState_TABLE
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ProgramHeaderModule,              ONLY: &
    DescribeProgramHeaderX, &
    nDOFX,                  &
    nx,                     &
    swX,                    &
    nDimsX,                 &
    nNodesX,                &
    nDOFE,                  &
    nNodesE,                &
    nDOFZ,                  &
    iZ_B0,                  &
    iZ_E0,                  &
    iX_B0,                  &
    iX_B1,                  &
    iX_E0
  USE RadiationFieldsModule,       ONLY: &
    nCR,    &
    nPR,    &
    iPR_D,  &
    iPR_I1, &
    iPR_I2, &
    iPR_I3, &
    iCR_N,  &
    iCR_G1, &
    iCR_G2, &
    iCR_G3
  USE FluidFieldsModule,       ONLY: &
    nCF,    &
    iCF_D,  &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iPF_Ne, &
    nPF,    &
    iPF_D,  &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E,  &
    iCF_Ne, &
    nAF,    &
    iAF_P, &
    iAF_T, &
    iAF_Ye, &
    iAF_S, &
    iAF_E, &
    iAF_Gm, &
    iAF_Cs, &
    iAF_Me, &
    iAF_Mp, &
    iAF_Mn, &
    iAF_Xp, &
    iAF_Xn, &
    iAF_Xa, &
    iAF_Xh
  USE GeometryFieldsModule,    ONLY: &
    uGF,          &
    nGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3,   &
    iGF_h_1,   &
    iGF_h_2,   &
    iGF_h_3,   &
    iGF_Phi_N, &
    iGF_Psi,   &
    iGF_SqrtGm
  USE TwoMoment_OpacityModule_Relativistic,  ONLY: &
   uOP, iOP_D0, iOP_Chi, iOP_Sigma, nOP
  USE MeshModule,              ONLY: &
    MeshType,    &
    CreateMesh,  &
    DestroyMesh, &
    MeshX,          &
    MeshE,          &
    NodeCoordinate
  USE ReferenceElementModule, ONLY: &
    OuterProduct1D3D
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE TwoMoment_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_TwoMoment
  USE UnitsModule,             ONLY: &
    UnitsDisplay, &
    Kilometer,    &
    SolarMass,    &
    SpeedOfLight, &
    Erg,          &
    Gram,         &
    Centimeter,   &
    MeV,          &
    BoltzmannConstant, &
    GravitationalConstant
! --- Local Modules ---
  USE InputParsingModule, ONLY: &
    nLevels, &
    xL,      &
    xR,      &
    nE,      &
    nSpecies, &
    Mass,     &
    UseTiling, &
    R0


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields

  REAL(AR), PARAMETER :: Zero   = 0.0_AR
  REAL(AR), PARAMETER :: Half   = 0.5_AR
  REAL(AR), PARAMETER :: One    = 1.0_AR
  REAL(AR), PARAMETER :: Two    = 2.0_AR
  REAL(AR), PARAMETER :: Three  = 3.0_AR
  REAL(AR), PARAMETER :: Pi     = ACOS( -1.0_AR )
  REAL(AR), PARAMETER :: TwoPi  = 2.0_AR * Pi
  REAL(AR), PARAMETER :: FourPi = 4.0_AR * Pi

CONTAINS

  SUBROUTINE MF_InitializeFields &
    ( ProgramName, MF_uGF, MF_uCR, MF_uCF, V_0, Verbose_Option, Direction_Option )

    CHARACTER(LEN=*),     INTENT(in   ) :: ProgramName
    TYPE(amrex_multifab), INTENT(inout   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1) 
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    REAL(AR),             INTENT(in)    :: V_0(3)    
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option
    CHARACTER(1), INTENT(in), OPTIONAL :: Direction_Option

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option


    SELECT CASE ( TRIM( ProgramName ) )
      
      CASE ( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming( MF_uGF, MF_uCR, MF_uCF, V_0, Direction_Option ) 

      CASE ( 'SineWaveDiffusion' )

        CALL InitializeFields_SineWaveDiffusion( MF_uGF, MF_uCR, MF_uCF, V_0 ) 

      CASE( 'StreamingDopplerShift' )

        CALL InitializeFields_StreamingDopplerShift ( MF_uGF, MF_uCR, MF_uCF, V_0, Direction_Option )

      CASE( 'HomogeneousSphere1D' )

        CALL InitializeFields_HomogeneousSphere1D ( MF_uGF, MF_uCR, MF_uCF )
      
      CASE( 'HomogeneousSphereGR' )

        CALL InitializeFields_HomogeneousSphereGR ( MF_uGF, MF_uCR, MF_uCF )

      CASE( 'Relaxation' )

        CALL InitializeFields_Relaxation ( MF_uGF, MF_uCR, MF_uCF )

        IF (Verbose) THEN
          print*, ProgramName
        END IF

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(4x,A,A)') 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(4x,A)')   'Valid Options:'
          WRITE(*,'(6x,A)')     'SineWaveStreaming'
          WRITE(*,'(6x,A)')     'SineWaveDiffusion'
          WRITE(*,'(6x,A)')     'StreamingDopplerShift'
          WRITE(*,'(6x,A)')     'HomogeneousSphere1D'
          WRITE(*,'(6x,A)')     'HomogeneousSphereGR'
          STOP 'MF_InitializationModule_NonRelativistic_IDEAL.f90'
        END IF

    END SELECT

  END SUBROUTINE MF_InitializeFields

  SUBROUTINE InitializeFields_SineWaveStreaming( MF_uGF, MF_uCR, MF_uCF, V_0, Direction )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    REAL(AR),             INTENT(in)    :: V_0(3)
    CHARACTER(1),         INTENT(in)    :: Direction

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(AR)       :: X1, X2, X3
    REAL(AR)       :: uCR_K( nDOFZ, nE, nCR, nSpecies )
    REAL(AR)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(AR)       :: uGF_K( nDOFX, nGF )
    REAL(AR)       :: uPF_K( nDOFX, nPF )
    REAL(AR)       :: uCF_K( nDOFX, nCF )
    REAL(AR)       :: uAF_K( nDOFX, nAF )
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR)                      :: Ones(nDOFE), W

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    Ones=1.0_AR
   
    W = 1.0_AR - ( V_0(1)*V_0(1) + V_0(2)*V_0(2) + V_0(3)*V_0(3) )
    W = 1.0_AR / SQRT( W )
   
    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO
    
    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uCR(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )
        
        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_C = LBOUND( uCR )
        hi_C = UBOUND( uCR )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)         
          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            uPF_K(iNodeX,iPF_D ) = 1.0_AR
            uPF_K(iNodeX,iPF_V1) = V_0(1)
            uPF_K(iNodeX,iPF_V2) = V_0(2)
            uPF_K(iNodeX,iPF_V3) = V_0(3)
            uPF_K(iNodeX,iPF_E ) = 0.1_AR
            uPF_K(iNodeX,iPF_Ne) = 0.0_AR
          END DO
        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )   
      
          DO iNodeZ = 1, nDOFZ
                        
            DO iS = 1, nSpecies
            DO iZ1 = 1, nE          

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)
              iNodeZ3 = NodeNumberTable(3,iNodeZ)
              iNodeZ4 = NodeNumberTable(4,iNodeZ)

              X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )
              X2 = NodeCoordinate( MeshX(2), iX2, iNodeZ3 )
              X3 = NodeCoordinate( MeshX(3), iX3, iNodeZ4 )
print*, X1

              SELECT CASE( TRIM( Direction ) )
           
              CASE( 'X' )


                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_AR + 0.49_AR * SIN( TwoPi * X1 ) 
 
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = W * uPR_K( iNodeZ, iZ1, iPR_D, iS )
         
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_AR

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_AR

              CASE( 'Y' )


                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_AR + 0.49_AR * SIN( TwoPi * X2 )  
            
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_AR
         
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = W * uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_AR

              CASE( 'Z' )


                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_AR + 0.49_AR * SIN( TwoPi * X3 )  
            
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_AR
         
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_AR

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = W * uPR_K( iNodeZ, iZ1, iPR_D, iS )
             
          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A8,A)')    '', 'InitializeFields_SineWaveStreaming'
            WRITE(*,'(A8,A,A2)') '', 'Invalid Direction: ', TRIM( Direction )
            WRITE(*,*)
            STOP

          END SELECT

 
            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &                         
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_AR, 0.0_AR, 0.0_AR,     &
                     uGF_K(iNodeX,iGF_Alpha), &
                     uGF_K(iNodeX,iGF_Beta_1), &
                     uGF_K(iNodeX,iGF_Beta_2), &
                     uGF_K(iNodeX,iGF_Beta_3) )

              
            END DO
            END DO
          END DO 
           
            uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
              = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

 
        END DO
        END DO
        END DO

      END DO
      
      CALL amrex_mfiter_destroy( MFI )

    END DO
  
  END SUBROUTINE InitializeFields_SineWaveStreaming


  SUBROUTINE InitializeFields_SineWaveDiffusion( MF_uGF, MF_uCR, MF_uCF, V_0 )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(AR)       :: X1, X2, X3, V_0(3)
    REAL(AR)       :: uCR_K( nDOFZ, nE, nCR, nSpecies )
    REAL(AR)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(AR)       :: uGF_K( nDOFX, nGF )
    REAL(AR)       :: uPF_K( nDOFX, nPF )
    REAL(AR)       :: uCF_K( nDOFX, nCF )
    REAL(AR)       :: uAF_K( nDOFX, nAF )
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR)                      :: Ones(nDOFE), W, Third

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    Ones=1.0_AR

    Third = 1.0_AR / 3.0_AR
  
    W = 1.0_AR - ( V_0(1)*V_0(1) + V_0(2)*V_0(2) + V_0(3)*V_0(3) )
    W = 1.0_AR / SQRT( W )
   
    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO
    
    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uCR(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )
        
        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_C = LBOUND( uCR )
        hi_C = UBOUND( uCR )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)         
          
          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            uPF_K(iNodeX,iPF_D ) = 1.0_AR
            uPF_K(iNodeX,iPF_V1) = V_0(1)
            uPF_K(iNodeX,iPF_V2) = V_0(2)
            uPF_K(iNodeX,iPF_V3) = V_0(3)
            uPF_K(iNodeX,iPF_E ) = 0.1_AR
            uPF_K(iNodeX,iPF_Ne) = 0.0_AR

          END DO
        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )   
      
          DO iNodeZ = 1, nDOFZ
                        
            DO iS = 1, nSpecies
            DO iZ1 = 1, nE          

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)

              X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )
    
              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                =  0.49_AR * SIN( Third * Pi * X1 ) + 0.5_AR  
            
              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                = - W * ( 0.49_AR * Pi / ( 9.0_AR * uOP(iNodeZ,iZ1,iX1,iX2,iX3,iOP_Sigma,iS) ) ) &
                  * COS( Third * Pi * X1 ) 
              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_AR

              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_AR

              
            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &                         
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_AR, 0.0_AR, 0.0_AR,     &
                     
                     1.0_AR, 0.0_AR, 0.0_AR, 0.0_AR )

              
            END DO
            END DO
!Reshape here instead of up top look at Hydro example
          END DO 
           
            uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
              = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

 
        END DO
        END DO
        END DO

      END DO
      
      CALL amrex_mfiter_destroy( MFI )

    END DO
  
  END SUBROUTINE InitializeFields_SineWaveDiffusion

  SUBROUTINE InitializeFields_StreamingDopplerShift( MF_uGF, MF_uCR, MF_uCF, V_0, Direction )


    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    REAL(AR), INTENT(in) :: V_0(3)
    CHARACTER(1),         INTENT(in)    :: Direction

    REAL(AR), PARAMETER :: X_0 = 2.0_AR
    REAL(AR), PARAMETER :: X_1 = 3.5_AR
    REAL(AR), PARAMETER :: X_2 = 6.5_AR
    REAL(AR), PARAMETER :: X_3 = 8.0_AR
    REAL(AR), PARAMETER :: L_X = 6.0_AR



    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(AR)       :: X1, X2, X3
    REAL(AR)       :: uCR_K( nDOFZ, nE, nCR, nSpecies )
    REAL(AR)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(AR)       :: uGF_K( nDOFX, nGF )
    REAL(AR)       :: uPF_K( nDOFX, nPF )
    REAL(AR)       :: uCF_K( nDOFX, nCF )
    REAL(AR)       :: uAF_K( nDOFX, nAF )
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR)                      :: Ones(nDOFE), W, Third, E, Mu_0, S
    CHARACTER(len=40) :: name1
    CHARACTER(len=1):: nds
    CHARACTER(len=2)::nxn1
    CHARACTER(len=3)::nxn2

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    Ones=1.0_AR

    Third = 1.0_AR / 3.0_AR
  
    W = 1.0_AR - ( V_0(1)*V_0(1) + V_0(2)*V_0(2) + V_0(3)*V_0(3) )
    W = 1.0_AR / SQRT( W )
  
    Mu_0 = 0.8_AR

 
    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO
    
    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uCR(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )
        
        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_C = LBOUND( uCR )
        hi_C = UBOUND( uCR )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)         
          
          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX



            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)
            iNodeX3 = NodeNumberTableX(3,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
            X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

            uPF_K(iNodeX,iPF_D ) = 1.0_AR 

            SELECT CASE( TRIM( Direction ) )

            CASE( 'X' ) 

              IF( X1 .LT. X_0 )THEN
                uPF_K(iNodeX,iPF_V1) &
                  = 0.0_AR
              ELSEIF( X1 .GE. X_0 .AND. X1 .LT. X_1 )THEN
                uPF_K(iNodeX,iPF_V1) &
                  = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
              ELSEIF( X1 .GE. X_1 .AND. X1 .LT. X_2 )THEN
               uPF_K(iNodeX,iPF_V1) &
                  = V_0(1)
              ELSEIF( X1 .GE. X_2 .AND. X1 .LT. X_3 )THEN
                uPF_K(iNodeX,iPF_V1) &
                  = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
              ELSE
                uPF_K(iNodeX,iPF_V1) &
                 = 0.0_AR 
              END IF

              uPF_K(iNodeX,iPF_V2) = V_0(2)
              uPF_K(iNodeX,iPF_V3) = V_0(3)


            CASE( 'Y' )

              uPF_K(iNodeX,iPF_V1) = V_0(1)
              IF( X2 .LT. X_0 )THEN
                uPF_K(iNodeX,iPF_V2) &
                  = 0.0_AR
              ELSEIF( X2 .GE. X_0 .AND. X2 .LT. X_1 )THEN
                uPF_K(iNodeX,iPF_V2) &
                  = V_0(2) * SIN( TwoPi * ( X2 - X_0 ) / L_X )**2
              ELSEIF( X2 .GE. X_1 .AND. X2 .LT. X_2 )THEN
                uPF_K(iNodeX,iPF_V2) &
                  = V_0(2)
              ELSEIF( X2 .GE. X_2 .AND. X2 .LT. X_3 )THEN
                uPF_K(iNodeX,iPF_V2) &
                  = V_0(2) * SIN( TwoPi * ( X2 - X_0 ) / L_X )**2
              ELSE
                uPF_K(iNodeX,iPF_V2) &
                  = 0.0_AR
              END IF
              uPF_K(iNodeX,iPF_V3) = V_0(3)

            CASE( 'Z' )

              uPF_K(iNodeX,iPF_V1) = V_0(1)
              uPF_K(iNodeX,iPF_V2) = V_0(2)
              IF( X3 .LT. X_0 )THEN
                uPF_K(iNodeX,iPF_V3) &
                  = 0.0_AR
              ELSEIF( X3 .GE. X_0 .AND. X3 .LT. X_1 )THEN
                uPF_K(iNodeX,iPF_V3) &
                  = V_0(3) * SIN( TwoPi * ( X3 - X_0 ) / L_X )**2
              ELSEIF( X3 .GE. X_1 .AND. X3 .LT. X_2 )THEN
                uPF_K(iNodeX,iPF_V3) &
                  = V_0(3)
              ELSEIF( X3 .GE. X_2 .AND. X3 .LT. X_3 )THEN
                uPF_K(iNodeX,iPF_V3) &
                  = V_0(3) * SIN( TwoPi * ( X3 - X_0 ) / L_X )**2
              ELSE
                uPF_K(iNodeX,iPF_V3) &
                  = 0.0_AR
              END IF

            CASE DEFAULT

              WRITE(*,*)
              WRITE(*,'(A8,A)')    '', 'InitializeFields_StreamingDopplerShift'
              WRITE(*,'(A8,A,A2)') '', 'Invalid Direction: ', TRIM( Direction )
              WRITE(*,*)
              STOP

            END SELECT


            uPF_K(iNodeX,iPF_E ) = 0.1_AR
            uPF_K(iNodeX,iPF_Ne) = 0.0_AR

          END DO

        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )
          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )   
      
          DO iNodeZ = 1, nDOFZ
                        
            DO iS = 1, nSpecies
            DO iZ1 = 1, nE          

              iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)
           
 

              W = 1.0_AR - (uPF_K(iNodeX,iPF_V1)**2 +  uPF_K(iNodeX,iPF_V2)**2 + uPF_K(iNodeX,iPF_V3)**2 )
 
              W = 1.0_AR / SQRT( W )

              X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )
 
              E = NodeCoordinate( MeshE, iZ1, iNodeE )
              IF(     TRIM( Direction ) .EQ. 'X' )THEN

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 1.0_AR / ( EXP( E / 3.0_AR - 3.0_AR ) + 1.0_AR )
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.99_AR * W * uPR_K( iNodeZ, iZ1, iPR_D, iS )
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_AR
                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_AR
              ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 1.0_AR / ( EXP( E / 3.0_AR - 3.0_AR ) + 1.0_AR )
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_AR
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.99_AR * W * uPR_K( iNodeZ, iZ1, iPR_D, iS )
                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_AR
              ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN
              
                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 1.0_AR / ( EXP( E / 3.0_AR - 3.0_AR ) + 1.0_AR )
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_AR
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_AR
                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.99_AR * W * uPR_K( iNodeZ, iZ1, iPR_D, iS )

              END IF

              CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &                         
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_AR, 0.0_AR, 0.0_AR,     &                     
                     1.0_AR, 0.0_AR, 0.0_AR, 0.0_AR )
            END DO
            END DO
          END DO 
           
            uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
              = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

 
        END DO
        END DO
        END DO

      END DO
      
      CALL amrex_mfiter_destroy( MFI )

    END DO
  END SUBROUTINE InitializeFields_StreamingDopplerShift



  SUBROUTINE InitializeFields_HomogeneousSphere1D( MF_uGF, MF_uCR, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(AR)       :: X1, X2, X3, V_0(3)
    REAL(AR)       :: uCR_K( nDOFZ, nE, nCR, nSpecies )
    REAL(AR)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(AR)       :: uGF_K( nDOFX, nGF )
    REAL(AR)       :: uPF_K( nDOFX, nPF )
    REAL(AR)       :: uCF_K( nDOFX, nCF )
    REAL(AR)       :: uAF_K( nDOFX, nAF )
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR)                      :: Ones(nDOFE), W

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    Ones=1.0_AR
   
   
    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO
    
    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uCR(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )
        
        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_C = LBOUND( uCR )
        hi_C = UBOUND( uCR )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)         
          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            uPF_K(iNodeX,iPF_D ) = 1.0_AR
            uPF_K(iNodeX,iPF_V1) = 0.0_AR
            uPF_K(iNodeX,iPF_V2) = 0.0_AR
            uPF_K(iNodeX,iPF_V3) = 0.0_AR
            uPF_K(iNodeX,iPF_E ) = 0.1_AR
            uPF_K(iNodeX,iPF_Ne) = 0.0_AR

          END DO
        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )


          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )   
      
          DO iNodeZ = 1, nDOFZ
                        
            DO iS = 1, nSpecies
            DO iZ1 = 1, nE          

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)

              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                = 10d-8  
              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                 = 0.0_AR
              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_AR

              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_AR

              
            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &                         
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_AR, 0.0_AR, 0.0_AR,     &
                     uGF_K(iNodeX,iGF_Alpha), &
                     uGF_K(iNodeX,iGF_Beta_1), &
                     uGF_K(iNodeX,iGF_Beta_2), &
                     uGF_K(iNodeX,iGF_Beta_3) )

              
            END DO
            END DO
!Reshape here instead of up top look at Hydro example
          END DO 
           
            uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
              = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

 
        END DO
        END DO
        END DO

      END DO
      
      CALL amrex_mfiter_destroy( MFI )

    END DO
  END SUBROUTINE InitializeFields_HomogeneousSphere1D




  SUBROUTINE InitializeFields_HomogeneousSphereGR( MF_uGF, MF_uCR, MF_uCF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(AR)       :: X1, X2, X3, V_0(3)
    REAL(AR)       :: uCR_K( nDOFZ, nE, nCR, nSpecies )
    REAL(AR)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(AR)       :: uGF_K( nDOFX, nGF )
    REAL(AR)       :: uPF_K( nDOFX, nPF )
    REAL(AR)       :: uCF_K( nDOFX, nCF )
    REAL(AR)       :: uAF_K( nDOFX, nAF )
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR)                      :: Ones(nDOFE), W, R, Theta

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    Ones=1.0_AR
   
   
    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
               xL(iDim), xR(iDim) )

    END DO
    
    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uCR(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )
        
        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_C = LBOUND( uCR )
        hi_C = UBOUND( uCR )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)      
 
          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            iNodeX2 = NodeNumberTableX(2,iNodeX)

            R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
IF (iX1 .NE. 0 .AND. iX1 .NE. nX(1)+1) THEN
print*, R / UnitsDisplay % LengthX1Unit
END IF
            theta = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            uPF_K(iNodeX,iPF_D ) = 1.0d12 * Gram / Centimeter**3
            uPF_K(iNodeX,iPF_V1) = 0.0_AR
            uPF_K(iNodeX,iPF_V2) = 0.0_AR
            uPF_K(iNodeX,iPF_V3) = 0.0_AR
            uPF_K(iNodeX,iPF_E ) = 1.00d25 * Erg / Centimeter**3
            uPF_K(iNodeX,iPF_Ne) = 0.0_AR

            CALL ComputeAlphaPsi( Mass, R, R0, theta, uGF_K(iNodeX,:) )
          END DO
        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )


          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )   
      
          uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)) &
            = RESHAPE( uGF_K, [ hi_G(4) - lo_G(4) + 1 ] )   

          DO iNodeZ = 1, nDOFZ
                        
            DO iS = 1, nSpecies
            DO iZ1 = 1, nE          

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)

              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                = 10d-20  
              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                 = 0.0_AR
              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_AR

              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_AR

              
            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &                         
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_AR, 0.0_AR, 0.0_AR,     &
                     uGF_K(iNodeX,iGF_Alpha), &
                     uGF_K(iNodeX,iGF_Beta_1), &
                     uGF_K(iNodeX,iGF_Beta_2), &
                     uGF_K(iNodeX,iGF_Beta_3) )

              
            END DO
            END DO
!Reshape here instead of up top look at Hydro example
          END DO 
           
            uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
              = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

 
        END DO
        END DO
        END DO

      END DO
      
      CALL amrex_mfiter_destroy( MFI )

    END DO
  END SUBROUTINE InitializeFields_HomogeneousSphereGR


  SUBROUTINE InitializeFields_Relaxation( MF_uGF, MF_uCR, MF_uCF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(AR)       :: X1, X2, X3, V_0(3)
    REAL(AR)       :: uCR_K( nDOFZ, nE, nCR, nSpecies )
    REAL(AR)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(AR)       :: uGF_K( nDOFX, nGF )
    REAL(AR)       :: uPF_K( nDOFX, nPF )
    REAL(AR)       :: uCF_K( nDOFX, nCF )
    REAL(AR)       :: uAF_K( nDOFX, nAF )
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR)                      :: kT, f_E, E, Mu, E12



!    REAL(AR), PARAMETER :: D_0   = 1.032d11 * Gram / Centimeter**3
!    REAL(AR), PARAMETER :: T_0   = 4.484d0 * MeV
!    REAL(AR), PARAMETER :: Y_0   = 0.1689_AR

!------

    REAL(AR), PARAMETER :: D_0   = 1.032d12 * Gram / Centimeter**3
    REAL(AR), PARAMETER :: T_0   = 7.588d0 * MeV
    REAL(AR), PARAMETER :: Y_0   = 0.1347_AR

!------

!    REAL(AR), PARAMETER :: D_0   = 1.022d13 * Gram / Centimeter**3
!    REAL(AR), PARAMETER :: T_0   = 1.617d1 * MeV
!    REAL(AR), PARAMETER :: Y_0   = 0.1421_AR


    REAL(AR), PARAMETER :: V_u_1 = 0.0_AR * SpeedOfLight
    REAL(AR), PARAMETER :: V_u_2 = 0.0_AR * SpeedOfLight
    REAL(AR), PARAMETER :: V_u_3 = 0.0_AR * SpeedOfLight
    REAL(AR), PARAMETER :: Mu_0  = -1.0_AR ! \in [-1,1]
     

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

   
   
    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
               xL(iDim), xR(iDim) )

    END DO
    
    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uCR(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )
        
        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_C = LBOUND( uCR )
        hi_C = UBOUND( uCR )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)      
 
          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            iNodeX2 = NodeNumberTableX(2,iNodeX)


            uPF_K(iNodeX,iPF_D ) = D_0
            uAF_K(iNodeX,iAF_T)  = T_0
            uAF_K(iNodeX,iAF_Ye) = Y_0

        CALL ComputeThermodynamicStates_Primitive_TABLE &
               ( uPF_K(iNodeX,iPF_D ), &
                 uAF_K(iNodeX,iAF_T ), &
                 uAF_K(iNodeX,iAF_Ye), &
                 uPF_K(iNodeX,iPF_E ), &
                 uAF_K(iNodeX,iAF_E ), &
                 uPF_K(iNodeX,iPF_Ne) )

        CALL ApplyEquationOfState_TABLE &
               ( uPF_K(iNodeX,iPF_D ), &
                 uAF_K(iNodeX,iAF_T ), &
                 uAF_K(iNodeX,iAF_Ye), &
                 uAF_K(iNodeX,iAF_P ), &
                 uAF_K(iNodeX,iAF_S ), &
                 uAF_K(iNodeX,iAF_E ), &
                 uAF_K(iNodeX,iAF_Me), &
                 uAF_K(iNodeX,iAF_Mp), &
                 uAF_K(iNodeX,iAF_Mn), &
                 uAF_K(iNodeX,iAF_Xp), &
                 uAF_K(iNodeX,iAF_Xn), &
                 uAF_K(iNodeX,iAF_Xa), &
                 uAF_K(iNodeX,iAF_Xh), &
                 uAF_K(iNodeX,iAF_Gm) )

        uPF_K(iNodeX,iPF_V1) = V_u_1
        uPF_K(iNodeX,iPF_V2) = V_u_2
        uPF_K(iNodeX,iPF_V3) = V_u_3


          END DO
        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )


          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )   
      
          uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)) &
            = RESHAPE( uGF_K, [ hi_G(4) - lo_G(4) + 1 ] )   

          DO iNodeZ = 1, nDOFZ
                        
            DO iS = 1, nSpecies
            DO iZ1 = 1, nE          

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

              kT = BoltzmannConstant * uAF_K(iNodeX,iAF_T)

              !Mu = uAF_K(iNodeX,iAF_Me) + uAF_K(iNodeX,iAF_Mp) - uAF_K(iNodeX,iAF_Mn)
              E = NodeCoordinate( MeshE, iZ1, iNodeE )

              IF (iNodeE .EQ. 2) THEN
               print*, E / MeV
              END IF


              f_E = MAX( 0.99_AR * EXP( - ( E - Two*kT )**2 &
                                    / ( Two*(1.0d1*MeV)**2 ) ), 1.0d-99 )
              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                = f_E * 0.50_AR * ( One - Mu_0 )
              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                = f_E * 0.25_AR * ( One - Mu_0**2 )
              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_AR

              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_AR
              
            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &                         
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_AR, 0.0_AR, 0.0_AR,     &
                     uGF_K(iNodeX,iGF_Alpha), &
                     uGF_K(iNodeX,iGF_Beta_1), &
                     uGF_K(iNodeX,iGF_Beta_2), &
                     uGF_K(iNodeX,iGF_Beta_3) )

              
            END DO
            END DO
          END DO 
           
            uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
              = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

 
        END DO
        END DO
        END DO

      END DO
      CALL amrex_mfiter_destroy( MFI )

    END DO



  END SUBROUTINE InitializeFields_Relaxation


  SUBROUTINE ComputeAlphaPsi( M, R, R0, theta, G ) 

    REAL(AR), INTENT(in) :: M, R, R0, theta
    REAL(AR), INTENT(inout) :: G(nGF)


    REAL(AR) :: Phi, Psi, a, b

    a = GravitationalConstant * M  / ( 2.0_AR * ( R0 )**3 )
    b = GravitationalConstant * M 

      !Need to add constants here
    IF ( R .LE. R0 ) THEN
        
      G(iGF_Phi_N) =  a * ( ( R )**2 - 3.0_AR * ( R0 ) **2 ) 

    ELSE

      G(iGF_Phi_N) = - b / ( R )

    END IF
   
    G(iGF_Phi_N) = G(iGF_Phi_N) / ( SpeedOfLight )**2

    

    G(iGF_Alpha) = 1.0_AR + G(iGF_Phi_N)      

    G(iGF_Psi) = 1.0_AR - G(iGF_Phi_N) / 2.0_AR


    G(iGF_h_1) = G(iGF_Psi)**2
    G(iGF_h_2) = G(iGF_Psi)**2 * R 
    G(iGF_h_3) = G(iGF_Psi)**2 * R * sin(theta)


    G(iGF_Gm_dd_11) = G(iGF_h_1)**2
    G(iGF_Gm_dd_22) = G(iGF_h_2)**2 
    G(iGF_Gm_dd_33) = G(iGF_h_3)**2 

    G(iGF_SqrtGm) = SQRT( G(iGF_Gm_dd_11) * G(iGF_Gm_dd_22) * G(iGF_Gm_dd_33) )   


  END SUBROUTINE ComputeAlphaPsi

END MODULE MF_InitializationModule
