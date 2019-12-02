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
    iX_E0
  USE RadiationFieldsModule,       ONLY: &
    nCR,    &
    uCR,    &
    nPR,    &
    uPR,    &
    nSpecies, &
    iPR_D,  &
    iPR_I1, &
    iPR_I2, &
    iPR_I3, &
    iCR_N,  &
    iCR_G1, &
    iCR_G2, &
    iCR_G3
  USE MeshModule,              ONLY: &
    MeshType,    &
    CreateMesh,  &
    DestroyMesh, &
    NodeCoordinate, &
    MeshX,          &
    MeshE,          &
    NodeCoordinate
  USE ReferenceElementModuleZ, ONLY: &
    NodeNumberTableZ
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
!  USE FluidFieldsModule, ONLY: &
 !   uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
  !  uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
!  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
 !   ComputeConserved_Euler_NonRelativistic
  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels, &
    xL,      &
    xR

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
    ( ProgramName, MF_uPR, MF_uCR )

    CHARACTER(LEN=*),     INTENT(in   ) :: ProgramName
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uPR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)

    SELECT CASE ( TRIM( ProgramName ) )
      
      CASE ( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming( MF_uPR, MF_uCR ) 
        print*, ProgramName

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(4x,A,A)') 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(4x,A)')   'Valid Options:'
          WRITE(*,'(6x,A)')     'IsentropicVortex'
          WRITE(*,'(6x,A)')     'Sod'
          WRITE(*,'(6x,A)')     'SphericalSod'
          WRITE(*,'(6x,A)')     'TopHatAdvection'
          WRITE(*,'(6x,A)')     'Implosion'
          WRITE(*,'(6x,A)')     'StandingAccretionShock'
          STOP 'MF_InitializationModule_NonRelativistic_IDEAL.f90'
        END IF

    END SELECT

  END SUBROUTINE MF_InitializeFields

  SUBROUTINE InitializeFields_SineWaveStreaming( MF_uPR, MF_uCR )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uPR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    !print*, 'hi'

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2
    REAL(AR)       :: X1, X2, X3, V_0(3)
    REAL(AR)       :: uCR_K( nDOFZ, nNodesE, nCR, nSpecies )
    REAL(AR)       :: uPR_K( nDOFZ, nNodesE, nPR, nSpecies )
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_P(4), hi_P(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uPR(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCR(:,:,:,:)

    V_0(1) = 0.1
    V_0(2) = 0.0
    V_0(3) = 0.0

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX
!
!        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_AR
!        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_0(1)
!        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_0(2)
!        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_0(3)
!        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_AR
!        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_AR
!
      END DO

!      CALL ComputeConserved_Euler_NonRelativistic &
!             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
!               uPF(:,iX1,iX2,iX3,iPF_V1), &
!               uPF(:,iX1,iX2,iX3,iPF_V2), &
!               uPF(:,iX1,iX2,iX3,iPF_V3), &
!               uPF(:,iX1,iX2,iX3,iPF_E ), &
!               uPF(:,iX1,iX2,iX3,iPF_Ne), &
!               uCF(:,iX1,iX2,iX3,iCF_D ), &
!               uCF(:,iX1,iX2,iX3,iCF_S1), &
!               uCF(:,iX1,iX2,iX3,iCF_S2), &
!               uCF(:,iX1,iX2,iX3,iCF_S3), &
!               uCF(:,iX1,iX2,iX3,iCF_E ), &
!               uCF(:,iX1,iX2,iX3,iCF_Ne), &
!               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
!               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
!               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO


    uCR_K = Zero
    uPR_K = Zero

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uPR(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uPR => MF_uPR(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_P = LBOUND( uPR )
        hi_P = UBOUND( uPR )

        lo_C = LBOUND( uCR )
        hi_C = UBOUND( uCR )
        
        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)
                    
          !uPR_K &
            != RESHAPE( uPR(iX1,iX2,iX3,lo_P(4):hi_P(4)), [ nDOFZ, nNodesE, nPR, nSpecies ] )
          DO iNodeZ = 1, nDOFZ
            
            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeZ2 = NodeNumberTableZ(2,iNodeZ)

            X1 = NodeCoordinate( MeshX(1), iX2, iNodeZ2 )
            
            uPR_K( iNodeZ, :, iPR_D, : ) &
              = 0.50_AR + 0.49_AR * SIN( TwoPi * X1 )  
            
            uPR_K( iNodeZ, :, iPR_I1, : ) &
              = uPR_K( iNodeZ, :, iPR_D, : )
         
            uPR_K( iNodeZ, :, iPR_I2, : ) &
              = 0.0_AR

            uPR_K( iNodeZ, :, iPR_I3, : ) &
              = 0.0_AR
            

          END DO         


 
        END DO
        END DO
        END DO
   


      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO


  END SUBROUTINE InitializeFields_SineWaveStreaming

END MODULE MF_InitializationModule
