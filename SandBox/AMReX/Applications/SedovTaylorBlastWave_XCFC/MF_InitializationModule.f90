MODULE MF_InitializationModule

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFX, &
    swX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    MeshType, &
    NodeCoordinate
  USE EquationOfStateModule_IDEAL, ONLY: &
    ComputePressureFromPrimitive_IDEAL, &
    Gamma_IDEAL
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    CoordinateSystem
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    iAF_P
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Two, &
    Three, &
    FourPi, &
    Pi
  USE InputParsingModule, ONLY: &
    UseTiling, &
    nMaxLevels, &
    t_end, &
    xR
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

  REAL(DP), PUBLIC :: IntE_Min_Euler_PL

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF, MF_uPF, MF_uAF

    LOGICAL :: Verbose

    TYPE(amrex_parmparse) :: PP

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAF(:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3), iX1, iX2, iX3, iNX

    CHARACTER(LEN=32) :: FMT

    ! --- Problem-Dependent Parameters ---

    INTEGER  :: iNX1, iNX2, iNX3, nDetCells
    REAL(DP) :: X1, X2, X3, Radius, dX1, dX2, dX3, &
                sigma, X_D, Edet, eCentral, rho0

    REAL(DP), PARAMETER :: RhoEratio = 1.0e2_DP
    REAL(DP), PARAMETER :: e1e2ratio = 1.0e-10_DP

    TYPE(MeshType) :: MeshXX(3)

    Verbose = .FALSE.
    IF( amrex_parallel_ioprocessor() .AND. iLevel .EQ. 0 ) Verbose = .TRUE.

    X_D = -1.0_DP
    CALL amrex_parmparse_build( PP, 'Sedov' )
      CALL PP % get  ( 'Edet'     , Edet )
      CALL PP % get  ( 'nDetCells', nDetCells )
      CALL PP % query( 'X_D'      , X_D )
    CALL amrex_parmparse_destroy( PP )

    CALL CreateMesh_MF( nMaxLevels-1, MeshXX )

    IF( X_D .LT. Zero )THEN

      IF( TRIM( CoordinateSystem ) .EQ. 'SPHERICAL' )THEN

        dX1 = MeshXX(1) % Width(0)

        X_D = SQRT( dX1**2 )

      ELSE IF( TRIM( CoordinateSystem ) .EQ. 'CYLINDRICAL' )THEN

        dX1 = MeshXX(1) % Width(0)
        dX2 = MeshXX(2) % Width(0)

        X_D = SQRT( dX1**2 + dX2**2 )

      ELSE

        dX1 = MeshXX(1) % Width(0)
        dX2 = MeshXX(2) % Width(0)
        dX3 = MeshXX(3) % Width(0)

        X_D = SQRT( dX1**2 + dX2**2 + dX3**2 )

      END IF

    END IF

    CALL DestroyMesh_MF( MeshXX )

    sigma = nDetCells * X_D

    eCentral &
      = Edet &
         / ( Pi * sigma**3 &
               * ( SQRT( Pi ) * ERF( xR(1) / sigma ) &
                     - Two * xR(1) / sigma * EXP( -xR(1)**2 / sigma**2 ) ) )

    ! --- Enforce non-relativistic specific enthalpy ---
    rho0 = RhoEratio * eCentral

    IF( Verbose )THEN

      WRITE(FMT,'(A)') '(A20,ES24.16E3)'

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'Initializing: ', TRIM( ProgramName )
      WRITE(*,'(4x,A)')   '-------------'
      WRITE(*,*)
      WRITE(*,'(A20,I4.4)') 'nDetCells: ', nDetCells
      WRITE(*,TRIM(FMT)) 'X_D: ', X_D
      WRITE(*,TRIM(FMT)) 'Edet: ', Edet
      WRITE(*,TRIM(FMT)) 'eCentral: ', eCentral
      WRITE(*,TRIM(FMT)) 't_end = ', t_end
      WRITE(*,TRIM(FMT)) 'rho0  = ', rho0
      WRITE(*,TRIM(FMT)) 'ener0 = ', e1e2ratio * eCentral
      WRITE(*,TRIM(FMT)) 'pres0 = ', e1e2ratio * eCentral * ( Gamma_IDEAL - 1.0_DP )
      WRITE(*,TRIM(FMT)) 'cs0   = ', &
        SQRT( Gamma_IDEAL * ( Gamma_IDEAL - 1.0e0_DP ) &
                * e1e2ratio / RhoEratio  )
      WRITE(*,*)

    END IF

    ! --- Map to computational domain ---

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )
      uPF => MF_uPF % DataPtr( MFI )
      uAF => MF_uAF % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1) + swX(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)
        iNX3 = NodeNumberTableX(3,iNX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

        IF( TRIM( CoordinateSystem ) .EQ. 'CARTESIAN' )THEN

          Radius = SQRT( X1**2 + X2**2 + X3**2 )

        ELSE IF( TRIM( CoordinateSystem ) .EQ. 'CYLINDRICAL' )THEN

          Radius = SQRT( X1**2 + X2**2 )

        ELSE

          Radius = SQRT( X1**2 )

        END IF

        uPF(iX1,iX2,iX3,nDOFX*(iPF_D -1)+iNX) = rho0
        uPF(iX1,iX2,iX3,nDOFX*(iPF_V1-1)+iNX) = Zero
        uPF(iX1,iX2,iX3,nDOFX*(iPF_V2-1)+iNX) = Zero
        uPF(iX1,iX2,iX3,nDOFX*(iPF_V3-1)+iNX) = Zero
        uPF(iX1,iX2,iX3,nDOFX*(iPF_Ne-1)+iNX) = Zero

        uPF(iX1,iX2,iX3,nDOFX*(iPF_E-1)+iNX) &
          = MAX( eCentral * EXP( -Radius**2 / sigma**2 ), &
                 1.0e2_DP * TINY( 1.0_DP ) ) ! avoid zero internal energy

        IntE_Min_Euler_PL &
          = MIN( IntE_Min_Euler_PL, &
                 1.0e-1_DP * uPF(iX1,iX2,iX3,nDOFX*(iPF_E-1)+iNX) )

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPF(iX1,iX2,iX3,nDOFX*(iPF_D -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_E -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_Ne-1)+iNX), &
                 uAF(iX1,iX2,iX3,nDOFX*(iAF_P -1)+iNX) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF(iX1,iX2,iX3,nDOFX*(iPF_D       -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_V1      -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_V2      -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_V3      -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_E       -1)+iNX), &
                 uPF(iX1,iX2,iX3,nDOFX*(iPF_Ne      -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_D       -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_S1      -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_S2      -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_S3      -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_E       -1)+iNX), &
                 uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne      -1)+iNX), &
                 uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                 uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                 uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                 uAF(iX1,iX2,iX3,nDOFX*(iAF_P       -1)+iNX) )

      END DO
      END DO
      END DO
      END DO

    END DO ! WHILE MFI % next()

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_MF


END MODULE MF_InitializationModule
