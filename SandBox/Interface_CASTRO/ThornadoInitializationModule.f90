MODULE ThornadoInitializationModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    InitializeProgramHeader, &
    nNodesX, nNodesE
  USE QuadratureModule, ONLY: &
    InitializeQuadratures
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement
  USE PolynomialBasisModule_Lagrange, ONLY: &
    InitializePolynomialBasis_Lagrange
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange
  USE MeshModule, ONLY: &
    MeshX, MeshE, &
    CreateMesh, &
    DestroyMesh
  USE GeometryFieldsModule, ONLY: &
    CreateGeometryFields, &
    DestroyGeometryFields
  USE FluidFieldsModule, ONLY: &
    CreateFluidFields, &
    DestroyFluidFields
  USE GeometryFieldsModuleE, ONLY: &
    CreateGeometryFieldsE, &
    DestroyGeometryFieldsE
  USE RadiationFieldsModule, ONLY: &
    CreateRadiationFields, &
    DestroyRadiationFields

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitThornado
  PUBLIC :: InitThornado_Patch
  PUBLIC :: FreeThornado_Patch

CONTAINS


  SUBROUTINE InitThornado( nDimsX, nDimsE )

    INTEGER, INTENT(in) :: nDimsX, nDimsE

    INTEGER :: nX(3), nE, i

    nX = 1
    DO i = 1, nDimsX
      nX(i) = nX(i) + 1
    END DO

    nE = 1 + nDimsE

    CALL InitializeProgramHeader &
           ( ProgramName_Option = '', nNodes_Option = 2, &
             nX_Option = nX, nE_Option = nE )

    CALL InitializeQuadratures

    CALL InitializeReferenceElementX

    CALL InitializeReferenceElementE

    CALL InitializeReferenceElement

    CALL InitializePolynomialBasis_Lagrange

    CALL InitializeReferenceElementX_Lagrange

    CALL InitializeReferenceElementE_Lagrange

    CALL InitializeReferenceElement_Lagrange

  END SUBROUTINE InitThornado


  SUBROUTINE InitThornado_Patch( nX, swX, xL, xR, nE, swE, eL, eR, nSpecies )

    INTEGER,  INTENT(in) :: nX(3), swX(3)
    REAL(DP), INTENT(in) :: xL(3), xR(3)
    INTEGER,  INTENT(in) :: nE, swE
    REAL(DP), INTENT(in) :: eL, eR
    INTEGER,  INTENT(in) :: nSpecies

    INTEGER :: iDim

    CALL InitializeProgramHeader &
           ( ProgramName_Option = '', nNodes_Option = 2, &
             nX_Option = nX, swX_Option = swX, &
             xL_Option = xL, xR_Option  = xR,  &
             nE_Option = nE, swE_Option = swE, &
             eL_Option = eL, eR_Option  = eR )

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), &
               swX(iDim), xL(iDim), xR(iDim) )

    END DO

    CALL CreateMesh( MeshE, nE, nNodesE, swE, eL, eR )

    CALL CreateGeometryFields &
           ( nX, swX, CoordinateSystem_Option = 'CARTESIAN', &
             Verbose_Option = .FALSE. )

    CALL CreateFluidFields &
           ( nX, swX, Verbose_Option = .FALSE. )

    CALL CreateGeometryFieldsE &
           ( nE, swE, Verbose_Option = .FALSE. )

    CALL CreateRadiationFields &
           ( nX, swX, nE, swE, nSpecies_Option = nSpecies, &
             Verbose_Option = .FALSE. )

  END SUBROUTINE InitThornado_Patch


  SUBROUTINE FreeThornado_Patch

    INTEGER :: iDim

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

    CALL DestroyMesh( MeshE )

    CALL DestroyGeometryFields

    CALL DestroyFluidFields

    CALL DestroyGeometryFieldsE

    CALL DestroyRadiationFields

  END SUBROUTINE FreeThornado_Patch


END MODULE ThornadoInitializationModule
