MODULE MF_TwoMoment_BoundaryConditionsModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,     ONLY: &
    AR => amrex_real, &
    amrex_spacedim
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment, &
    iApplyBC_TwoMoment_Inner, &
    iApplyBC_TwoMoment_Both
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDOFE, nDOFX, &
    nE, swX, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    ProgramName
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF

  ! --- Local Modules ---

  USE InputParsingModule,              ONLY: &
    DEBUG, &
    UseTiling, &
    nLevels
  USE MF_EdgeMapModule,                ONLY: &
    EdgeMap, &
    ConstructEdgeMap
  USE MF_UtilitiesModule,              ONLY: &
    amrex2thornado_Z, &
    thornado2amrex_Z

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_TwoMoment_MF

  REAL(AR), ALLOCATABLE, PUBLIC :: uCR_iBC(:,:,:,:) ! (nDOF,iE_B0:iE_E0,nCR,nSpecies)

  INTERFACE ApplyBoundaryConditions_TwoMoment_MF
    MODULE PROCEDURE ApplyBoundaryConditions_TwoMoment_MF_MultiLevel
    MODULE PROCEDURE ApplyBoundaryConditions_TwoMoment_MF_SingleLevel
    MODULE PROCEDURE ApplyBoundaryConditions_TwoMoment_MF_SingleLevel_Box
  END INTERFACE ApplyBoundaryConditions_TwoMoment_MF


CONTAINS


  LOGICAL FUNCTION UseInnerBC_TwoMoment()

    UseInnerBC_TwoMoment &
      = ALLOCATED( uCR_iBC ) .AND. &
        ( TRIM( ProgramName ) .EQ. 'TransparentVortex'           .OR. &
          TRIM( ProgramName ) .EQ. 'TransparentVortex_Spherical' )

    RETURN
  END FUNCTION UseInnerBC_TwoMoment


  SUBROUTINE ApplyBoundaryConditions_TwoMoment_MF_MultiLevel( MF_uCR )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      CALL ApplyBoundaryConditions_TwoMoment_MF( iLevel, MF_uCR(iLevel) )

    END DO

  END SUBROUTINE ApplyBoundaryConditions_TwoMoment_MF_MultiLevel


  SUBROUTINE ApplyBoundaryConditions_TwoMoment_MF_SingleLevel( iLevel, MF_uCR )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX
    TYPE(EdgeMap)      :: Edge_Map

    REAL(AR), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(AR), ALLOCATABLE         :: U  (:,:,:,:,:,:,:)

    INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    INTEGER :: iZ_B (4), iZ_E (4), iLo_MF(4)
    INTEGER :: iZ2, iZ3, iZ4, iCR, iS

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = .FALSE. )

    DO WHILE( MFI % next() )

      uCR => MF_uCR % DataPtr( MFI )

      iLo_MF = LBOUND( uCR )

      BX = MFI % tilebox()

      iZ_B0 = [ iE_B0, BX % lo(1)       , BX % lo(2)       , BX % lo(3)        ]
      iZ_E0 = [ iE_E0, BX % hi(1)       , BX % hi(2)       , BX % hi(3)        ]
      iZ_B1 = [ iE_B1, BX % lo(1)-swX(1), BX % lo(2)-swX(2), BX % lo(3)-swX(3) ]
      iZ_E1 = [ iE_E1, BX % hi(1)+swX(1), BX % hi(2)+swX(2), BX % hi(3)+swX(3) ]

      ALLOCATE( U(1:nDOF, &
                  iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
                  iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
                  1:nCR, 1:nSpecies) )

      iZ_B = [ iE_B0, iZ_B1(2), iZ_B1(3), iZ_B1(4) ]
      iZ_E = [ iE_E0, iZ_E1(2), iZ_E1(3), iZ_E1(4) ]

      CALL amrex2thornado_Z &
             ( nCR, nSpecies, nE, iE_B0, iE_E0, &
               iZ_B1, iZ_E1, iLo_MF, iZ_B, iZ_E, uCR, U )

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_TwoMoment_MF &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, Edge_Map )

      IF( UseInnerBC_TwoMoment() .AND. &
          BX % lo(1) .EQ. amrex_geom(iLevel) % domain % lo(1) )THEN

        DO iS  = 1, nSpecies
        DO iCR = 1, nCR
        DO iZ4 = iZ_B1(4), iZ_E1(4)
        DO iZ3 = iZ_B1(3), iZ_E1(3)
        DO iZ2 = iZ_B1(2), iZ_B0(2) - 1

          U(:,iE_B0:iE_E0,iZ2,iZ3,iZ4,iCR,iS) = uCR_iBC(:,:,iCR,iS)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      CALL thornado2amrex_Z &
             ( nCR, nSpecies, nE, iE_B0, iE_E0, &
               iZ_B1, iZ_E1, iLo_MF, iZ_B, iZ_E, uCR, U )

      DEALLOCATE( U )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ApplyBoundaryConditions_TwoMoment_MF_SingleLevel


  SUBROUTINE ApplyBoundaryConditions_TwoMoment_MF_SingleLevel_Box &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, Edge_Map )

    INTEGER,       INTENT(in   ) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(AR),      INTENT(inout) :: &
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)
    TYPE(EdgeMap), INTENT(in   ) :: &
      Edge_Map

    INTEGER :: iApplyBC(3)

    INTEGER :: iZ2, iZ3, iZ4, iCR, iS

    CALL Edge_Map % GetBC( iApplyBC )

    IF( DEBUG ) WRITE(*,'(A)') '      CALL ApplyBoundaryConditions_TwoMoment'

    CALL ApplyBoundaryConditions_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, &
             iApplyBC_Option = iApplyBC )


    IF( UseInnerBC_TwoMoment() .AND. &
        ( iApplyBC(1) .EQ. iApplyBC_TwoMoment_Inner .OR. &
          iApplyBC(1) .EQ. iApplyBC_TwoMoment_Both ) )THEN

      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ4 = iZ_B1(4), iZ_E1(4)
      DO iZ3 = iZ_B1(3), iZ_E1(3)
      DO iZ2 = iZ_B1(2), iZ_B0(2) - 1

        U(:,iE_B0:iE_E0,iZ2,iZ3,iZ4,iCR,iS) = uCR_iBC(:,:,iCR,iS)

      END DO
      END DO
      END DO
      END DO
      END DO

    END IF

  END SUBROUTINE ApplyBoundaryConditions_TwoMoment_MF_SingleLevel_Box


END MODULE MF_TwoMoment_BoundaryConditionsModule