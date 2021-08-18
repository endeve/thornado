MODULE MF_InitializationModule_Relativistic_IDEAL

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne
  USE Euler_UtilitiesModule, ONLY: &
    ComputeConserved_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    TwoPi
  USE MF_UtilitiesModule, ONLY: &
    thornado2amrex_X, &
    amrex2thornado_X
  USE InputParsingModule, ONLY: &
    nX, &
    swX, &
    Gamma_IDEAL, &
    ProgramName
  USE MF_Euler_ErrorModule, ONLY: &
    DescribeError_Euler_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Euler_MF


CONTAINS


  SUBROUTINE InitializeFields_Euler_MF( iLevel, MF_uGF, MF_uCF )

    INTEGER,              INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'Advection2D' )

        CALL InitializeFields_Advection2D( iLevel, MF_uGF, MF_uCF )

      CASE DEFAULT

        CALL DescribeError_Euler_MF &
               ( 01, Message_Option &
                       = 'Invalid ProgramName: ' // TRIM( ProgramName ) )

    END SELECT

  END SUBROUTINE InitializeFields_Euler_MF


  SUBROUTINE InitializeFields_Advection2D( iLevel, MF_uGF, MF_uCF )

    INTEGER,              INTENT(in) :: iLevel
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF, MF_uCF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER  :: iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: uPF(nDOFX,nPF)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-specific parameters ---

    INTEGER  :: iNX1, iNX2
    REAL(DP) :: X1, X2

    CHARACTER(:), ALLOCATABLE :: AdvectionProfile

    ! --- Gaussian ---

    REAL(DP), PARAMETER :: D_0    = 2.0_DP
    REAL(DP), PARAMETER :: X1_0   = 0.5_DP
    REAL(DP), PARAMETER :: X2_0   = 0.5_DP
    REAL(DP), PARAMETER :: Radius = 0.1_DP
    REAL(DP), PARAMETER :: sigma  = 0.1_DP

    AdvectionProfile = 'Gaussian'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'AdvectionProfile', AdvectionProfile )
    CALL amrex_parmparse_destroy( PP )

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(2x,A,A)') 'Advection Profile: ', TRIM( AdvectionProfile )

    END IF

    CALL amrex_mfiter_build( MFI, MF_uGF )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      ALLOCATE( G (1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3), &
                   1:nGF) )

      ALLOCATE( U (1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3), &
                   1:nCF) )

      CALL amrex2thornado_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        X1   = NodeCoordinate( MeshX(1), iX1, iNX1 )

        iNX2 = NodeNumberTableX(2,iNX)
        X2   = NodeCoordinate( MeshX(2), iX2, iNX2 )

        IF( TRIM( AdvectionProfile ) .EQ. 'Gaussian' )THEN

          IF( ( X1 - X1_0 )**2 + ( X2 - X2_0 )**2 .LT. Radius**2 )THEN

            uPF(iNX,iPF_D) &
              = D_0 * EXP( ( X1 - X1_0 )**2 / ( 2.0_DP * sigma**2 ) ) &
                    * EXP( ( X2 - X2_0 )**2 / ( 2.0_DP * sigma**2 ) )

          ELSE

            uPF(iNX,iPF_D) = 1.0_DP

          END IF

          uPF(iNX,iPF_V1) = 0.5_DP
          uPF(iNX,iPF_V2) = 0.0_DP
          uPF(iNX,iPF_V3) = 0.0_DP
          uPF(iNX,iPF_E ) = One / ( Gamma_IDEAL - One )
          uPF(iNX,iPF_Ne) = Zero

        ELSE IF( TRIM( AdvectionProfile ) .EQ. 'SineWaveX1' )THEN

          uPF(iNX,iPF_D ) = 1.0_DP + 0.1_DP * SIN( TwoPi * X1 )
          uPF(iNX,iPF_V1) = 0.1_DP
          uPF(iNX,iPF_V2) = 0.0_DP
          uPF(iNX,iPF_V3) = 0.0_DP
          uPF(iNX,iPF_E ) = One / ( Gamma_IDEAL - One )
          uPF(iNX,iPF_Ne) = Zero

        END IF

        CALL ComputeConserved_Euler &
               ( uPF(iNX,iPF_D ), &
                 uPF(iNX,iPF_V1), &
                 uPF(iNX,iPF_V2), &
                 uPF(iNX,iPF_V3), &
                 uPF(iNX,iPF_E ), &
                 uPF(iNX,iPF_Ne), &
                 U(iNX,iX1,iX2,iX3,iCF_D ), &
                 U(iNX,iX1,iX2,iX3,iCF_S1), &
                 U(iNX,iX1,iX2,iX3,iCF_S2), &
                 U(iNX,iX1,iX2,iX3,iCF_S3), &
                 U(iNX,iX1,iX2,iX3,iCF_E ), &
                 U(iNX,iX1,iX2,iX3,iCF_Ne), &
                 G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 One )

      END DO
      END DO
      END DO
      END DO

      CALL thornado2amrex_X &
             ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

      DEALLOCATE( U )
      DEALLOCATE( G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_Advection2D


END MODULE MF_InitializationModule_Relativistic_IDEAL
