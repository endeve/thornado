!> Module for operations on MultiFabs
MODULE MF_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_multifab, &
    amrex_imultifab

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MakeFineMaskModule, ONLY: &
    MakeFineMask, &
    iCoarse_MFM
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_DataTransfer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ShowVariableFromMultiFab
  PUBLIC :: MultiplyWithMetric
  PUBLIC :: amrex2thornado_X
  PUBLIC :: thornado2amrex_X
  PUBLIC :: thornado2amrex_X_F
  PUBLIC :: amrex2thornado_X_F

  INTERFACE ShowVariableFromMultiFab
    MODULE PROCEDURE ShowVariableFromMultiFab_Single
    MODULE PROCEDURE ShowVariableFromMultiFab_Vector
  END INTERFACE ShowVariableFromMultiFab

CONTAINS


  SUBROUTINE ShowVariableFromMultiFab_Single &
    ( iLevel, MF, iField, iMF_Mask, &
      swXX_Option, WriteToFile_Option, FileName_Option )

    INTEGER              , INTENT(in) :: iLevel, iField
    TYPE(amrex_multifab) , INTENT(in) :: MF
    TYPE(amrex_imultifab), INTENT(in) :: iMF_Mask
    INTEGER              , INTENT(in), OPTIONAL :: swXX_Option(3)
    LOGICAL              , INTENT(in), OPTIONAL :: WriteToFile_Option
    CHARACTER(*)         , INTENT(in), OPTIONAL :: FileName_Option

    INTEGER                       :: iX1, iX2, iX3, iNX
    INTEGER                       :: lo(4), hi(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: U(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: Mask(:,:,:,:)
    INTEGER                       :: swXX(3)
    LOGICAL                       :: WriteToFile
    CHARACTER(128)                :: FMT
    CHARACTER(128)                :: FileName

    REAL(DP) :: NodesX1(nNodesX(1))
    REAL(DP) :: NodesX2(nNodesX(2))
    REAL(DP) :: NodesX3(nNodesX(3))

    swXX = 0
    IF( PRESENT( swXX_Option ) ) swXX = swXX_Option

    WriteToFile = .FALSE.
    IF( PRESENT( WriteToFile_Option ) ) WriteToFile = WriteToFile_Option

    FileName = ''
    IF( PRESENT( FileName_Option ) ) FileName = TRIM( FileName_Option )

    WRITE(FMT,'(A,I2.2,A,I2.2,A,I2.2,A,I3.3,A)') &
      '(I2.2,3I4.3,3ES25.16E3,', &
      nNodesX(1),  'ES25.16E3,', &
      nNodesX(2),  'ES25.16E3,', &
      nNodesX(3),  'ES25.16E3,', &
      nDOFX     ,  'ES25.16E3)'

    IF( WriteToFile ) OPEN( 100, FILE = TRIM( FileName ), POSITION = 'APPEND' )

    CALL amrex_mfiter_build( MFI, MF, tiling = UseTiling )

    CALL CreateMesh_MF( iLevel, MeshX )

    DO WHILE( MFI % next() )

      IF( nLevels .GT. 1 .AND. iLevel .LT. nLevels-1 ) &
        Mask => iMF_Mask % DataPtr( MFI )

      U => MF % DataPtr( MFI )
      BX = MFI % tilebox()

      lo = LBOUND( U ); hi = UBOUND( U )

      DO iX3 = BX % lo(3) - swXX(3), BX % hi(3) + swXX(3)
      DO iX2 = BX % lo(2) - swXX(2), BX % hi(2) + swXX(2)
      DO iX1 = BX % lo(1) - swXX(1), BX % hi(1) + swXX(1)

        IF( nLevels .GT. 1 .AND. iLevel .LT. nLevels-1 )THEN

          IF( Mask(iX1,iX2,iX3,1) .NE. iCoarse_MFM ) CYCLE

        END IF

        DO iNX = 1, nNodesX(1)
          NodesX1(iNX) = NodeCoordinate( MeshX(1), iX1, iNX )
        END DO
        DO iNX = 1, nNodesX(2)
          NodesX2(iNX) = NodeCoordinate( MeshX(2), iX2, iNX )
        END DO
        DO iNX = 1, nNodesX(3)
          NodesX3(iNX) = NodeCoordinate( MeshX(3), iX3, iNX )
        END DO

        IF( WriteToFile )THEN

          WRITE(100,TRIM(FMT)) &
            iLevel, iX1, iX2, iX3, &
            MeshX(1) % Width(iX1), &
            MeshX(2) % Width(iX2), &
            MeshX(3) % Width(iX3), &
            NodesX1, NodesX2, NodesX3, &
            U(iX1,iX2,iX3,1+nDOFX*(iField-1):nDOFX*iField)

        ELSE

          WRITE(*,TRIM(FMT)) &
            iLevel, iX1, iX2, iX3, &
            MeshX(1) % Width(iX1), &
            MeshX(2) % Width(iX2), &
            MeshX(3) % Width(iX3), &
            NodesX1, NodesX2, NodesX3, &
            U(iX1,iX2,iX3,1+nDOFX*(iField-1):nDOFX*iField)

        END IF

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshX )

    IF( WriteToFile) CLOSE(100)

    WRITE(*,*)

  END SUBROUTINE ShowVariableFromMultiFab_Single


  SUBROUTINE ShowVariableFromMultiFab_Vector &
    ( MF, iField, swXX_Option, WriteToFile_Option, FileName_Option )

    INTEGER             , INTENT(in) :: iField
    TYPE(amrex_multifab), INTENT(in) :: MF(0:nLevels-1)
    INTEGER             , INTENT(in), OPTIONAL :: swXX_Option(3)
    LOGICAL             , INTENT(in), OPTIONAL :: WriteToFile_Option
    CHARACTER(*)        , INTENT(in), OPTIONAL :: FileName_Option

    INTEGER :: iLevel

    INTEGER        :: swXX(3)
    LOGICAL        :: WriteToFile
    CHARACTER(128) :: FileName

    TYPE(amrex_imultifab) :: iMF_Mask

    swXX = 0
    IF( PRESENT( swXX_Option ) ) swXX = swXX_Option

    WriteToFile = .FALSE.
    IF( PRESENT( WriteToFile_Option ) ) WriteToFile = WriteToFile_Option

    FileName = ''
    IF( PRESENT( FileName_Option ) ) FileName = TRIM( FileName_Option )

    DO iLevel = 0, nLevels-1

      IF( nLevels .GT. 1 .AND. iLevel .LT. nLevels-1 ) &
        CALL MakeFineMask &
               ( iMF_Mask, MF(iLevel) % BA, MF(iLevel) % DM, &
                 MF(iLevel+1) % BA )

      CALL ShowVariableFromMultiFab_Single &
             ( iLevel, MF(iLevel), iField, iMF_Mask, &
               swXX_Option = swXX, &
               WriteToFile_Option = WriteToFile, &
               FileName_Option = TRIM( FileName ) )

    END DO

  END SUBROUTINE ShowVariableFromMultiFab_Vector


  SUBROUTINE MultiplyWithMetric( MF_uGF, MF_uCF, Power )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF
    INTEGER             , INTENT(in)    :: Power

    INTEGER                       :: iX1, iX2, iX3, iNX, iCF
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_U(4), hi_U(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: G(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U(:,:,:,:)
    REAL(DP)                      :: G_K(nDOFX,nGF)
    REAL(DP)                      :: U_K(nDOFX,nCF)

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      G => MF_uGF % DataPtr( MFI )
      U => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( G ); hi_G = UBOUND( G )
      lo_U = LBOUND( U ); hi_U = UBOUND( U )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        G_K(1:nDOFX,1:nGF) &
          = RESHAPE( G(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        U_K(1:nDOFX,1:nCF) &
          = RESHAPE( U(iX1,iX2,iX3,lo_U(4):hi_U(4)), [ nDOFX, nCF ] )

        DO iCF = 1, nCF
        DO iNX = 1, nDOFX

          U_K(iNX,iCF) = U_K(iNX,iCF) * G_K(iNX,iGF_SqrtGm)**( Power )

        END DO
        END DO

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE MultiplyWithMetric


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


  SUBROUTINE thornado2amrex_X_F &
    ( nDOFX_X, nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, &
      Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nDOFX_X, nFields
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

      Data_amrex(iX1,iX2,iX3,nDOFX_X*(iFd-1)+1:nDOFX_X*iFd) &
        = Data_thornado(1:nDOFX_X,iX1,iX2,iX3,iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

  END SUBROUTINE thornado2amrex_X_F


  SUBROUTINE amrex2thornado_X_F &
    ( nDOFX_X, nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, &
      Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nDOFX_X, nFields
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

      Data_thornado(1:nDOFX_X,iX1,iX2,iX3,iFd) &
        = Data_amrex(iX1,iX2,iX3,nDOFX_X*(iFd-1)+1:nDOFX_X*iFd)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

  END SUBROUTINE amrex2thornado_X_F


END MODULE MF_UtilitiesModule
