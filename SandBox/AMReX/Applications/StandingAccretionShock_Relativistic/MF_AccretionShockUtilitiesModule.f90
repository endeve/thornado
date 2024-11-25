MODULE MF_AccretionShockUtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator, &
    amrex_parallel_reduce_sum, &
    amrex_parallel_myproc
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_imultifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE UtilitiesModule, ONLY: &
    thornado_abort
  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDimsX, &
    nDOFX, &
    nNodesX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX1, &
    WeightsX2, &
    WeightsX_q, &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX2_q
  USE MeshModule, ONLY: &
    MeshType, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm, &
    iGF_Psi
  USE FluidFieldsModule, ONLY: &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iAF_P
  USE Euler_BoundaryConditionsModule, ONLY: &
    ExpD, &
    ExpE
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond, &
    Centimeter, &
    Second

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    FourPi
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_UtilitiesModule, ONLY: &
    ShowVariableFromMultiFab, &
    IndLo_X, &
    IndHi_X
  USE InputParsingModule, ONLY: &
    nLevels, &
    t_new, &
    StepNo, &
    PlotFileNameRoot
  USE FillPatchModule, ONLY: &
    FillPatch
  USE AverageDownModule, ONLY: &
    AverageDown
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler_MF
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WriteNodal1DICToFile_SAS
  PUBLIC :: ComputePowerInLegendreModes

  LOGICAL,          PUBLIC              :: WriteNodal1DIC_SAS
  CHARACTER(LEN=:), PUBLIC, ALLOCATABLE :: FileName_Nodal1DIC_SAS

  INTEGER, PARAMETER :: nLeg = 5 ! number of Legendre polynomials
  INTEGER, PARAMETER :: SL = 256

  LOGICAL :: ComputePowerInSitu

CONTAINS


  SUBROUTINE WriteNodal1DICToFile_SAS

    TYPE(amrex_parmparse) :: PP

    INTEGER :: iLevel, iErr

    CHARACTER(256) :: FileName, DirName

    IF( nDimsX .GT. 1 ) RETURN

    DO iLevel = 0, nLevels - 1

      CALL FillPatch( iLevel, MF_uGF )
      CALL FillPatch( iLevel, MF_uGF, MF_uCF )

    END DO

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF, &
             swXX_Option = swX )

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get( 'FileName_Nodal1DIC_SAS', &
                      FileName_Nodal1DIC_SAS )
    CALL amrex_parmparse_destroy( PP )

    WRITE(FileName,'(A)') TRIM( FileName_Nodal1DIC_SAS ) // '_PF_D'
    CALL ShowVariableFromMultiFab &
           ( MF_uPF, iPF_D, swXX_Option = swX, UseFineMask_Option = .FALSE., &
             WriteToFile_Option = .TRUE., &
             FileNameBase_Option = TRIM( FileName ) )

    WRITE(FileName,'(A)') TRIM( FileName_Nodal1DIC_SAS ) // '_PF_V1'
    CALL ShowVariableFromMultiFab &
           ( MF_uPF, iPF_V1, swXX_Option = swX, UseFineMask_Option = .FALSE., &
             WriteToFile_Option = .TRUE., &
             FileNameBase_Option = TRIM( FileName ) )

    WRITE(FileName,'(A)') TRIM( FileName_Nodal1DIC_SAS ) // '_AF_P'
    CALL ShowVariableFromMultiFab &
           ( MF_uAF, iAF_P, swXX_Option = swX, UseFineMask_Option = .FALSE., &
             WriteToFile_Option = .TRUE., &
             FileNameBase_Option = TRIM( FileName ) )

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(DirName,'(A,I8.8,A)') &
        TRIM( PlotFileNameRoot ), StepNo(0), '_nodal/'

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() ) &
        CALL SYSTEM( 'mkdir -p ' // TRIM( DirName ) // ' 2>/dev/null' )

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      OPEN( UNIT = 101, &
            FILE = TRIM( DirName ) &
                     // TRIM( FileName_Nodal1DIC_SAS ) // '_BC.dat' )

      WRITE(101,'(ES24.16E3)') ExpD
      WRITE(101,'(ES24.16E3)') ExpE

    END IF

  END SUBROUTINE WriteNodal1DICToFile_SAS


  SUBROUTINE ComputePowerInLegendreModes( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uPF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uAF(0:)

    TYPE(amrex_multifab) :: DivV2(0:nLevels-1)

    INTEGER :: iLevel, nX

    REAL(DP) :: Power(0:nLeg-1)
    REAL(DP), ALLOCATABLE :: G(:,:)
    REAL(DP), ALLOCATABLE :: Psi(:)

    TYPE(amrex_parmparse) :: PP

    ComputePowerInSitu = .FALSE.
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % query( 'ComputePower', ComputePowerInSitu )
    CALL amrex_parmparse_destroy( PP )

    IF( .NOT. ComputePowerInSitu ) RETURN

    CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )
    CALL ApplyBoundaryConditions_Euler_MF   ( MF_uCF )

    DO iLevel = 0, nLevels - 1

      CALL FillPatch( iLevel, MF_uGF )
      CALL FillPatch( iLevel, MF_uGF, MF_uCF )

    END DO

    DO iLevel = 0, nLevels-1

      CALL MF_uPF(iLevel) % SetVal( Zero )
      CALL MF_uAF(iLevel) % SetVal( Zero )

    END DO

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF, swXX_Option = [ 0, 1, 0 ] )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( DivV2(iLevel), &
               MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, nDOFX, 0 )
      CALL DivV2(iLevel) % SetVal( Zero )

    END DO

    CALL ComputeDivV2( MF_uPF, DivV2 )

    nX = amrex_geom(0) % domain % hi(1) - amrex_geom(0) % domain % lo(1) + 1

    ALLOCATE( G(0:nLeg-1,1:nX*nNodesX(1)) )
    ALLOCATE( Psi       (1:nX*nNodesX(1)) )

    CALL ComputePowerDensity( MF_uGF, DivV2, G, Psi )

    CALL ComputePower( G, Psi, Power )

    CALL WriteToFile( Power )

    DEALLOCATE( Psi )
    DEALLOCATE( G   )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( DivV2(iLevel) )

    END DO

  END SUBROUTINE ComputePowerInLegendreModes


  SUBROUTINE ComputeDivV2( MF_uPF, DivV2 )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uPF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: DivV2 (0:)

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_imultifab) :: iMF_FineMask

    INTEGER :: iX_B0(3), iX_E0(3), iX1, iX2, iX3, iNX, iNX2
    INTEGER :: iLevel, iLo, iHi

    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: dV2     (:,:,:,:)

    REAL(DP) :: X2, dV2dX2(nDOFX)

    TYPE(MeshType) :: MeshX(3)

    iLo = IndLo_X( nDOFX, iPF_V2 )
    iHi = IndHi_X( nDOFX, iPF_V2 )

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uPF % BA, MF_uPF % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uPF(iLevel) )

      DO WHILE( MFI % next() )

        FineMask => iMF_FineMask         % DataPtr( MFI )
        uPF      => MF_uPF      (iLevel) % DataPtr( MFI )
        dV2      => DivV2       (iLevel) % DataPtr( MFI )

        BX = MFI % TileBox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

          dV2dX2 &
            = One / MeshX(2) % Width(iX2) &
                * MATMUL( dLXdX2_q, uPF(iX1,iX2,iX3,iLo:iHi) )

          DO iNX = 1, nDOFX

            iNX2 = NodeNumberTableX(2,iNX)

            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

            dV2(iX1,iX2,iX3,iNX) &
              = dV2(iX1,iX2,iX3,iNX) &
                  + dV2dX2(iNX) * SIN( X2 ) &
                  + uPF(iX1,iX2,iX3,iLo+iNX-1) * COS( X2 )

          END DO

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel

    CALL AverageDown( DivV2 )

  END SUBROUTINE ComputeDivV2


  SUBROUTINE ComputePowerDensity( MF_uGF, DivV2, G, Psi )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)  :: DivV2 (0:)
    REAL(DP)            , INTENT(out) :: G(0:,1:)
    REAL(DP)            , INTENT(out) :: Psi (1:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    INTEGER :: iX_B0(3), iX_E0(3), iX1, iX2, iNX1, jNX(nNodesX(2))
    INTEGER :: ell, iDOFX1

    REAL(DP), CONTIGUOUS, POINTER :: dV2(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    REAL(DP) :: dX2, X2(nNodesX(2)), Pell(nNodesX(2))

    TYPE(MeshType) :: MeshX(3)

real(dp)::Pell2(0:nLeg-1,nNodesX(2)),PP(0:nLeg-1,0:nLeg-1)
integer::ell1,ell2
logical::debug
debug=.false.
if(debug)pp=zero

    Psi = Zero
    G   = Zero

    CALL CreateMesh_MF( 0, MeshX )

    CALL amrex_mfiter_build( MFI, DivV2(0) )

    DO WHILE( MFI % next() )

      dV2 => DivV2 (0) % DataPtr( MFI )
      uGF => MF_uGF(0) % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi

      DO ell = 0, nLeg-1

        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX1 = 1, nNodesX(1)

            iDOFX1 &
              = ( iX1 - amrex_geom(0) % domain % lo(1) ) * nNodesX(1) + iNX1

            ! --- Assume psi is spherically-symmetric ---
            IF( ell .EQ.0 ) &
              Psi(iDOFX1) = uGF(iX1,iX_B0(2),iX_B0(3),nDOFX*(iGF_Psi-1)+iNX1)

            CALL GetjNX( iNX1, jNX )

            DO iX2 = iX_B0(2), iX_E0(2)

              dX2 = MeshX(2) % Width (iX2)
              X2  = MeshX(2) % Center(iX2) + MeshX(2) % Nodes * dX2

              CALL ComputeLegendrePolynomial( ell, COS( X2 ), Pell )

              G(ell,iDOFX1) &
                = G(ell,iDOFX1) &
                    + dX2 * SUM( WeightsX2 &
                                 * dV2(iX1,iX2,iX_B0(3),jNX) * Pell )

if(debug)then
if(ix1.eq.amrex_geom(0)%domain%lo(1).and.inx1.eq.1.and.ell.eq.0)then
do ell1=0,nleg-1
do ell2=0,nleg-1
CALL ComputeLegendrePolynomial( ell1, COS( X2 ), Pell2(ell1,:) )
CALL ComputeLegendrePolynomial( ell2, COS( X2 ), Pell2(ell2,:) )
PP(ell1,ell2) &
  = PP(ell1,ell2) &
      + dX2 * SUM( WeightsX2 * Pell2(ell1,:) * Pell2(ell2,:) * SIN( X2 ) )
enddo
enddo
endif
endif

            END DO ! iX2

          END DO ! iNX1

        END DO ! iX1

      END DO ! ell

    END DO ! WHILE

    CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshX )

    DO iX1 = amrex_geom(0) % domain % lo(1), amrex_geom(0) % domain % hi(1)

      DO iNX1 = 1, nNodesX(1)

        iDOFX1 &
          = ( iX1 - amrex_geom(0) % domain % lo(1) ) * nNodesX(1) + iNX1

        CALL amrex_parallel_reduce_sum( Psi(iDOFX1) )

        DO ell = 0, nLeg-1

          CALL amrex_parallel_reduce_sum( G(ell,iDOFX1) )

        END DO

      END DO

    END DO

if(debug)then
do ell1 = 0, nleg-1
do ell2 = 0, nleg-1
call amrex_parallel_reduce_sum( pp(ell1,ell2) )
end do
end do
if(amrex_parallel_ioprocessor())then
do ell1 = 0, nleg-1
print'(5ES13.3E3)',PP(ell1,:)
enddo
call thornado_abort
endif
endif
  END SUBROUTINE ComputePowerDensity


  SUBROUTINE ComputePower( G, Psi, Power )

    REAL(DP) , INTENT(in)  :: G(0:,1:), Psi(1:)
    REAL(DP) , INTENT(out) :: Power(0:)

    INTEGER :: iX1, ell, iLo, iHi

    REAL(DP) :: X1(nNodesX(1)), dX1, ra, rb, ShockRadius

    TYPE(MeshType) :: MeshX(3)

    TYPE(amrex_parmparse) :: PP

    IF( .NOT. amrex_parallel_ioprocessor() ) RETURN

    ra = 0.8_DP
    rb = 0.9_DP
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % query( 'ShellBottom', ra )
      CALL PP % query( 'ShellTop'   , rb )
      CALL PP % get  ( 'ShockRadius', ShockRadius )
    CALL amrex_parmparse_destroy( PP )

    ra = ra * ShockRadius * Kilometer
    rb = rb * ShockRadius * Kilometer

    Power = Zero

    CALL CreateMesh_MF( 0, MeshX )

    DO iX1 = amrex_geom(0) % domain % lo(1), amrex_geom(0) % domain % hi(1)

      IF( .NOT. ( MeshX(1) % Center(iX1) .GT. ra &
            .AND. MeshX(1) % Center(iX1) .LT. rb ) ) CYCLE

      dX1 = MeshX(1) % Width (iX1)
      X1  = MeshX(1) % Center(iX1) + MeshX(1) % Nodes * dX1

      iLo = ( iX1 - amrex_geom(0) % domain % lo(1) ) * nNodesX(1) + 1
      iHi = iLo + nNodesX(1) - 1

      DO ell = 0, nLeg-1

        Power(ell) &
          = Power(ell) &
              + FourPi &
                  * SUM( WeightsX1 * G(ell,iLo:iHi)**2 * Psi(iLo:iHi)**6 &
                           * X1**2 ) * dX1

     END DO ! ell

    END DO ! iX1

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ComputePower


  SUBROUTINE WriteToFile( Power )

    REAL(DP), INTENT(in) :: Power(0:nLeg-1)

    TYPE(amrex_parmparse) :: PP

    INTEGER, PARAMETER :: FileNameLen = 256
    CHARACTER(256) :: FMT
    CHARACTER(:), ALLOCATABLE :: PowerFileName_T
    CHARACTER(FileNameLen) :: PowerFileName

    IF( .NOT. amrex_parallel_ioprocessor() ) RETURN

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get( 'PowerFileName', PowerFileName_T )
    CALL amrex_parmparse_destroy( PP )

    WRITE(FMT,'(A3,I2.2,A10)') '(SP', nLeg+1, 'ES25.16E3)'

    PowerFileName = TRIM( PowerFileName_T )

    CALL CheckFileExistenceAndAppend( FileNameLen, PowerFileName )

    OPEN( 100, FILE = TRIM( PowerFileName ), POSITION = 'APPEND' )

    WRITE(100,TRIM(FMT)) &
      t_new(0) / Millisecond, Power / ( Centimeter**3 / Second**2 )

    CLOSE( 100 )

  END SUBROUTINE WriteToFile


  SUBROUTINE ComputeLegendrePolynomial( ell, x, Pell )

    INTEGER , INTENT(in)   :: ell
    REAL(DP), INTENT(in)   :: x   (nNodesX(2))
    REAL(DP), INTENT(out)  :: Pell(nNodesX(2))

    SELECT CASE( ell )

    CASE( 0 )

      Pell = SQRT( 1.0_DP / 2.0_DP )

    CASE( 1 )

      Pell = SQRT( 3.0_DP / 2.0_DP ) * x

    CASE( 2 )

      Pell = SQRT( 5.0_DP / 2.0_DP ) * ( 3.0_DP * x**2 - 1.0_DP ) / 2.0_DP

    CASE( 3 )

      Pell = SQRT( 7.0_DP / 2.0_DP ) &
               * 1.0_DP / 2.0_DP * ( 5.0_DP * x**3 - 3.0_DP * x )

    CASE( 4 )

      Pell = SQRT( 9.0_DP / 2.0_DP ) &
               * 1.0_DP / 8.0_DP * ( 35.0_DP * x**4 - 30.0_DP * x**2 + 3.0_DP )

    END SELECT

  END SUBROUTINE ComputeLegendrePolynomial


  SUBROUTINE GetjNX( iNX1, jNX )

    INTEGER, INTENT(in)  :: iNX1
    INTEGER, INTENT(out) :: jNX(nNodesX(2))

    INTEGER :: i

    jNX(1) = iNX1

    DO i = 2, nNodesX(2)

      jNX(i) = jNX(i-1) + nNodesX(2)

    END DO

  END SUBROUTINE GetjNX


  RECURSIVE SUBROUTINE CheckFileExistenceAndAppend &
    ( FileNameLen, FileName, IntSuffix_Option )

    INTEGER                   , INTENT(in)    :: FileNameLen
    CHARACTER(LEN=FileNameLen), INTENT(inout) :: FileName
    INTEGER                   , INTENT(inout), OPTIONAL :: IntSuffix_Option

    LOGICAL :: IsFile
    INTEGER :: IntSuffix
    INTEGER :: SL_T

    IntSuffix = 1
    IF( PRESENT( IntSuffix_Option ) ) &
      IntSuffix = IntSuffix_Option

    SL_T = LEN( TRIM( FileName ) )

    INQUIRE( FILE = TRIM( FileName ), EXIST = IsFile )

    IF( IsFile )THEN

      IF( FileName(SL_T-3:SL_T) .EQ. '.dat' )THEN

        WRITE(FileName,'(A,A,I2.2)') TRIM( FileName ), '_', IntSuffix

      ELSE

        WRITE(FileName(SL_T-1:SL_T),'(I2.2)') IntSuffix

      END IF

      IntSuffix = IntSuffix + 1

      CALL CheckFileExistenceAndAppend &
             ( FileNameLen, FileName, IntSuffix_Option = IntSuffix )

    END IF

  END SUBROUTINE CheckFileExistenceAndAppend


END MODULE MF_AccretionShockUtilitiesModule
