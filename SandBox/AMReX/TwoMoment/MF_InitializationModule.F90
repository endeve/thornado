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
  USE ReferenceElementModule, ONLY: &
    OuterProduct1D3D
  USE ReferenceElementModuleZ, ONLY: &
    NodeNumberTableZ
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment
  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels, &
    xL,      &
    xR,      &
    nE,      &
    nSpecies


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
    ( ProgramName, MF_uGF, MF_uCR, MF_uCF )

    CHARACTER(LEN=*),     INTENT(in   ) :: ProgramName
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1) 
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    SELECT CASE ( TRIM( ProgramName ) )
      
      CASE ( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming( MF_uGF, MF_uCR, MF_uCF ) 
        print*, ProgramName

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(4x,A,A)') 'Unknown Program: ', TRIM( ProgramName )
          WRITE(*,'(4x,A)')   'Valid Options:'
          WRITE(*,'(6x,A)')     'SineWaveStreaming'
          STOP 'MF_InitializationModule_NonRelativistic_IDEAL.f90'
        END IF

    END SELECT

  END SUBROUTINE MF_InitializeFields

  SUBROUTINE InitializeFields_SineWaveStreaming( MF_uGF, MF_uCR, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    !print*, 'hi'

    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(AR)       :: X1, X2, X3, V_0(3)
    REAL(AR)       :: uCR_K( nDOFZ, nE, nCR, nSpecies )
    REAL(AR)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    TYPE(MeshType) :: MeshX(3)

    ! --- AMReX ---
    INTEGER                       :: iLevel
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_P(4), hi_P(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uPR(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(AR)                      :: Ones(nDOFE)
    REAL(AR)     :: Gm_dd_11(nDOFZ)
    REAL(AR)     :: Gm_dd_22(nDOFZ)
    REAL(AR)     :: Gm_dd_33(nDOFZ)

    uCR_K = Zero
    uPR_K = Zero

    print*, nSpecies
   
    Ones=1.0_AR

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uCR => MF_uCR(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_P = LBOUND( uPR )
        hi_P = UBOUND( uPR )

        lo_C = LBOUND( uCR )
        hi_C = UBOUND( uCR )
        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)         

          Gm_dd_11 &
            = OuterProduct1D3D &
                ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), nDOFX )

          Gm_dd_22 &
            = OuterProduct1D3D &
                ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), nDOFX )

          Gm_dd_33 &
            = OuterProduct1D3D &
                ( Ones, nDOFE, uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), nDOFX )

          DO iNodeZ = 1, nDOFZ
            
            DO iS = 1, nSpecies
            DO iZ1 = 1, nE          

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTableZ(2,iNodeZ)

              X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )
    
              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                = iZ1  !0.50_AR + 0.49_AR * SIN( TwoPi * X1 )  
            
              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                = iS/10.0_AR !uPR_K( iNodeZ, iZ1, iPR_D, iS )
         
              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_AR

              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_AR


            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(:,iZ1,iPR_D,iS), &
                     uPR_K(:,iZ1,iPR_I1,iS), &
                     uPR_K(:,iZ1,iPR_I2,iS), &
                     uPR_K(:,iZ1,iPR_I3,iS), &
                     uCR_K(:,iZ1,iCR_N,iS), &
                     uCR_K(:,iZ1,iCR_G1,iS), &
                     uCR_K(:,iZ1,iCR_G2,iS), &
                     uCR_K(:,iZ1,iCR_G3,iS), &
                     Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:))

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
    
  END SUBROUTINE InitializeFields_SineWaveStreaming

END MODULE MF_InitializationModule
