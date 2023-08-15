!> Module for operations on MultiFabs
MODULE MF_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Psi, &
    iGF_SqrtGm, &
    unitsGF
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
    iPF_Ne, &
    unitsPF, &
    nDF, &
    unitsAF, &
    nAF, &
    iAF_P
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE UtilitiesModule, ONLY: &
    IsCornerCell

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE InputParsingModule, ONLY: &
    nLevels, &
    nX, &
    swX, &
    UseTiling
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    EdgeMap, &
    ConstructEdgeMap, &
    MF_ApplyBoundaryConditions_Euler
  USE TimersModule_AMReX_Euler, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_DataTransfer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: amrex2thornado_X
  PUBLIC :: thornado2amrex_X
  PUBLIC :: amrex2thornado_X_Global
  PUBLIC :: ShowVariableFromMultiFab
  PUBLIC :: WriteNodalDataToFile
  PUBLIC :: CombineGridData


CONTAINS


  SUBROUTINE amrex2thornado_X &
    ( nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      Data_amrex   (iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(out) :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_thornado(1:nDOFX,iX1,iX2,iX3,iFd) &
        = Data_amrex(iX1,iX2,iX3,nDOFX*(iFd-1)+1:nDOFX*iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

  END SUBROUTINE amrex2thornado_X


  SUBROUTINE thornado2amrex_X &
    ( nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(DP), INTENT(out) :: &
      Data_amrex   (iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(DP), INTENT(in)  :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_amrex(iX1,iX2,iX3,nDOFX*(iFd-1)+1:nDOFX*iFd) &
        = Data_thornado(1:nDOFX,iX1,iX2,iX3,iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

  END SUBROUTINE thornado2amrex_X


  SUBROUTINE amrex2thornado_X_Global( GEOM, MF, nF, U, ApplyBC_Option )

    TYPE(amrex_geometry), INTENT(in)  :: GEOM(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)  :: MF(0:nLevels-1)
    INTEGER,              INTENT(in)  :: nF
    REAL(DP),             INTENT(out) :: U(1:,1-swX(1):,1-swX(2):,1-swX(3):,1:)
    LOGICAL,              INTENT(in), OPTIONAL :: ApplyBC_Option

    INTEGER                       :: iNX, iX1, iX2, iX3, iFd, iLevel, &
                                     iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
                                     iX_B (3), iX_E (3)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    TYPE(EdgeMap)                 :: Edge_Map
    LOGICAL                       :: ApplyBC
    REAL(DP), CONTIGUOUS, POINTER :: uA(:,:,:,:)
    REAL(DP), ALLOCATABLE         :: uT(:,:,:,:,:)

    ApplyBC = .FALSE.
    IF( PRESENT( ApplyBC_Option ) ) &
      ApplyBC = ApplyBC_Option

    DO iLevel = 0, nLevels-1

      CALL MF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL amrex_mfiter_build( MFI, MF(iLevel), tiling = UseTiling )

      U = Zero

      DO WHILE( MFI % next() )

        uA => MF(iLevel) % DataPtr( MFI )
        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        IF( ApplyBC )THEN

          ALLOCATE( uT(nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3), &
                       nF) )

          CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )

          CALL amrex2thornado_X &
                 ( nF, iX_B1, iX_E1, LBOUND( uA ), iX_B1, iX_E1, uA, uT )

          CALL MF_ApplyBoundaryConditions_Euler &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, uT, Edge_Map )

          CALL thornado2amrex_X &
                 ( nF, iX_B1, iX_E1, LBOUND( uA ), iX_B1, iX_E1, uA, uT )

          DEALLOCATE( uT )

        END IF

        iX_B = iX_B0
        iX_E = iX_E0

        IF( BX % lo(1) .EQ. 1     ) iX_B(1) = 1     - swX(1)
        IF( BX % lo(2) .EQ. 1     ) iX_B(2) = 1     - swX(2)
        IF( BX % lo(3) .EQ. 1     ) iX_B(3) = 1     - swX(3)
        IF( BX % hi(1) .EQ. nX(1) ) iX_E(1) = nX(1) + swX(1)
        IF( BX % hi(2) .EQ. nX(2) ) iX_E(2) = nX(2) + swX(2)
        IF( BX % hi(3) .EQ. nX(3) ) iX_E(3) = nX(3) + swX(3)

        DO iFd = 1, nF
        DO iX3 = iX_B(3), iX_E(3)
        DO iX2 = iX_B(2), iX_E(2)
        DO iX1 = iX_B(1), iX_E(1)

          U(1:nDOFX,iX1,iX2,iX3,iFd) &
            = uA(iX1,iX2,iX3,nDOFX*(iFd-1)+1:nDOFX*iFd)

        END DO
        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

      DO iFd = 1, nF
      DO iX3 = 1-swX(3), nX(3)+swX(3)
      DO iX2 = 1-swX(2), nX(2)+swX(2)
      DO iX1 = 1-swX(1), nX(1)+swX(1)
      DO iNX = 1, nDOFX

        CALL amrex_parallel_reduce_sum( U(iNX,iX1,iX2,iX3,iFd) )

      END DO
      END DO
      END DO
      END DO
      END DO

    END DO

  END SUBROUTINE amrex2thornado_X_Global


  SUBROUTINE ShowVariableFromMultiFab( MF, swXX, iComp )

    TYPE(amrex_multifab), INTENT(in) :: MF(0:nLevels-1)
    INTEGER,              INTENT(in) :: swXX(3)
    INTEGER,              INTENT(in) :: iComp

    INTEGER                       :: iX1, iX2, iX3, iLevel
    INTEGER                       :: lo(4), hi(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: U(:,:,:,:)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        U => MF(iLevel) % DataPtr( MFI )
        BX = MFI % tilebox()

        lo = LBOUND( U ); hi = UBOUND( U )

        DO iX3 = BX % lo(3) - swXX(3), BX % hi(3) + swXX(3)
        DO iX2 = BX % lo(2) - swXX(2), BX % hi(2) + swXX(2)
        DO iX1 = BX % lo(1) - swXX(1), BX % hi(1) + swXX(1)

          WRITE(*,'(A,3I4.3,ES10.1E3)') &
            'iX1, iX2, iX3, Data: ',iX1, iX2, iX3, U(iX1,iX2,iX3,iComp)

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    WRITE(*,*)

  END SUBROUTINE ShowVariableFromMultiFab


  SUBROUTINE WriteNodalDataToFile( GEOM, MF_uGF, MF_uCF, MF_uDF, FileNameBase )

    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uDF(0:nLevels-1)
    CHARACTER(LEN=*)    , INTENT(in) :: FileNameBase

    INTEGER           :: iLo(3), iHi(3), iNX, iNX1, iX1, iX2, iX3
    CHARACTER(LEN=16) :: FMT

    REAL(DP) :: P(1:nDOFX,1:nPF)
    REAL(DP) :: A(1:nDOFX,1:nAF)
    REAL(DP) :: G(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                          1-swX(2):nX(2)+swX(2), &
                          1-swX(3):nX(3)+swX(3), &
                  1:nGF)
    REAL(DP) :: U(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                          1-swX(2):nX(2)+swX(2), &
                          1-swX(3):nX(3)+swX(3), &
                  1:nCF)
    REAL(DP) :: D(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                          1-swX(2):nX(2)+swX(2), &
                          1-swX(3):nX(3)+swX(3), &
                  1:nDF)

    CALL amrex2thornado_X_Global &
           ( GEOM, MF_uGF, nGF, G, ApplyBC_Option = .FALSE. )

    CALL amrex2thornado_X_Global &
           ( GEOM, MF_uCF, nCF, U, ApplyBC_Option = .TRUE. )

    CALL amrex2thornado_X_Global &
           ( GEOM, MF_uDF, nDF, D, ApplyBC_Option = .FALSE. )

    IF( amrex_parallel_ioprocessor() )THEN

      iLo = 1  - swX
      iHi = nX + swX

      OPEN( UNIT = 100, FILE = TRIM( FileNameBase ) // 'r.dat'      )
      OPEN( UNIT = 101, FILE = TRIM( FileNameBase ) // 'Alpha.dat'  )
      OPEN( UNIT = 102, FILE = TRIM( FileNameBase ) // 'Psi.dat'    )
      OPEN( UNIT = 103, FILE = TRIM( FileNameBase ) // 'SqrtGm.dat' )
      OPEN( UNIT = 104, FILE = TRIM( FileNameBase ) // 'D.dat'      )
      OPEN( UNIT = 105, FILE = TRIM( FileNameBase ) // 'V.dat'      )
      OPEN( UNIT = 106, FILE = TRIM( FileNameBase ) // 'E.dat'      )
      OPEN( UNIT = 107, FILE = TRIM( FileNameBase ) // 'P.dat'      )

      WRITE(FMT,'(A3,I3.3,A10)') '(SP', nDOFX, 'ES25.16E3)'

      DO iX3 = iLo(3), iHi(3)
      DO iX2 = iLo(2), iHi(2)
      DO iX1 = iLo(1), iHi(1)

        IF( IsCornerCell( iLo, iHi, iX1, iX2, iX3 ) ) CYCLE

        CALL ComputePrimitive_Euler &
               ( U   (:,iX1,iX2,iX3,iCF_D ), &
                 U   (:,iX1,iX2,iX3,iCF_S1), &
                 U   (:,iX1,iX2,iX3,iCF_S2), &
                 U   (:,iX1,iX2,iX3,iCF_S3), &
                 U   (:,iX1,iX2,iX3,iCF_E ), &
                 U   (:,iX1,iX2,iX3,iCF_Ne), &
                 P   (:,iPF_D ), &
                 P   (:,iPF_V1), &
                 P   (:,iPF_V2), &
                 P   (:,iPF_V3), &
                 P   (:,iPF_E ), &
                 P   (:,iPF_Ne), &
                 G   (:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G   (:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G   (:,iX1,iX2,iX3,iGF_Gm_dd_33) )

        CALL ComputePressureFromPrimitive &
               ( P(:,iPF_D ), P(:,iPF_E ), P(:,iPF_Ne), A(:,iAF_P) )

        DO iNX = 1, nDOFX

          iNX1 = NodeNumberTableX(1,iNX)

          WRITE(100,'(ES24.16E3,1x)',ADVANCE='NO') &
            NodeCoordinate( MeshX(1), iX1, iNX1 )

        END DO

        WRITE(100,*)

        WRITE(101,FMT) &
          G(:,iX1,iX2,iX3,iGF_Alpha) &
            / unitsGF(iGF_Alpha)

        WRITE(102,FMT) &
          G(:,iX1,iX2,iX3,iGF_Psi) &
            / unitsGF(iGF_Psi)

        WRITE(103,FMT) &
          G(:,iX1,iX2,iX3,iGF_SqrtGm) &
            / unitsGF(iGF_SqrtGm)

        WRITE(104,FMT) &
          P(:,iPF_D) &
            / unitsPF(iPF_D)

        WRITE(105,FMT) &
          P(:,iPF_V1) &
            / unitsPF(iPF_V1)

        WRITE(106,FMT) &
          P(:,iPF_E) &
            / unitsPF(iPF_E)

        WRITE(107,FMT) &
          A(:,iAF_P) &
            / unitsAF(iAF_P)

      END DO
      END DO
      END DO

      CLOSE( 107 )
      CLOSE( 106 )
      CLOSE( 105 )
      CLOSE( 104 )
      CLOSE( 103 )
      CLOSE( 102 )
      CLOSE( 101 )
      CLOSE( 100 )

    END IF

  END SUBROUTINE WriteNodalDataToFile


  SUBROUTINE CombineGridData( MF, nF, iField, U )

    TYPE(amrex_multifab), INTENT(in)  :: MF(0:nLevels-1)
    INTEGER,              INTENT(in)  :: nF, iField
    REAL(DP),             INTENT(out) :: U(1:,1-swX(1):,1-swX(2):,1-swX(3):)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: U_P(:,:,:,:)
    REAL(DP)                      :: U_K(nDOFX,nF)

    REAL(DP) :: U0(nDOFX,1-swX(1):nX(1)+swX(1), &
                         1-swX(2):nX(2)+swX(2), &
                         1-swX(3):nX(3)+swX(3),0:nLevels-1)

    INTEGER  :: iNodeX, iX1, iX2, iX3, iLevel, iX_B(3), iX_E(3), lo(4), hi(4)

    U0 = Zero

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        U_P => MF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo = LBOUND( U_P )
        hi = UBOUND( U_P )

        iX_B = BX % lo
        iX_E = BX % hi

        ! -- Get physical ghost cells right ---

        IF( BX % lo(1) .EQ. 1     ) iX_B(1) = iX_B(1) - swX(1)
        IF( BX % lo(2) .EQ. 1     ) iX_B(2) = iX_B(2) - swX(2)
        IF( BX % lo(3) .EQ. 1     ) iX_B(3) = iX_B(3) - swX(3)
        IF( BX % hi(1) .EQ. nX(1) ) iX_E(1) = iX_E(1) + swX(1)
        IF( BX % hi(2) .EQ. nX(2) ) iX_E(2) = iX_E(2) + swX(2)
        IF( BX % hi(3) .EQ. nX(3) ) iX_E(3) = iX_E(3) + swX(3)

        DO iX3 = iX_B(3), iX_E(3)
        DO iX2 = iX_B(2), iX_E(2)
        DO iX1 = iX_B(1), iX_E(1)

          U_K = RESHAPE( U_P(iX1,iX2,iX3,lo(4):hi(4)), [ nDOFX, nF ] )

          U0(:,iX1,iX2,iX3,iLevel) = U_K(:,iField)

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    ! --- Combine data from different grids ---

    DO iX3 = 1-swX(3), nX(3)+swX(3)
    DO iX2 = 1-swX(2), nX(2)+swX(2)
    DO iX1 = 1-swX(1), nX(1)+swX(1)

      DO iNodeX = 1, nDOFX

        CALL amrex_parallel_reduce_sum( U0(iNodeX,iX1,iX2,iX3,:), nLevels )

        U(iNodeX,iX1,iX2,iX3) = U0(iNodeX,iX1,iX2,iX3,0)

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE CombineGridData


END MODULE MF_UtilitiesModule
