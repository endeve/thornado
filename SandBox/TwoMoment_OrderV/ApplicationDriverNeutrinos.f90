PROGRAM ApplicationDriverNeutrinos

  USE KindModule, ONLY: &
    DP, Zero, One
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE

  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: TimeSteppingScheme
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UsePositivityLimiter
  INTEGER       :: nSpecies
  INTEGER       :: nNodes
  INTEGER       :: nE, bcE, nX(3), bcX(3)
  REAL(DP)      :: xL(3), xR(3), ZoomX(3) = One
  REAL(DP)      :: eL, eR, ZoomE = One
  REAL(DP)      :: t_end

  ProgramName = 'Relaxation'

  CoordinateSystem = 'CARTESIAN'

  SELECT CASE( TRIM( ProgramName ) )

    CASE( 'Relaxation' )

      nSpecies = 2
      nNodes   = 2

      nX  = [ 1, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer
      bcX = [ 0, 0, 0 ]

      nE    = 16
      eL    = 0.0d0 * MeV
      eR    = 3.0d2 * MeV
      bcE   = 0
      ZoomE = 1.266038160710160_DP

      TimeSteppingScheme = 'BackwardEuler'

      t_end = 1.0d0 * Millisecond

      UseSlopeLimiter = .FALSE.

      UsePositivityLimiter = .FALSE.

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', 'Unknown Program Name: ', TRIM( ProgramName )
      WRITE(*,*)
      STOP

  END SELECT

  CALL InitializeDriver

  CALL FinalizeDriver

CONTAINS


  SUBROUTINE InitializeDriver

    USE TwoMoment_TimersModule_OrderV, ONLY: &
      InitializeTimers
    USE ProgramInitializationModule, ONLY: &
      InitializeProgram
    USE ReferenceElementModuleX, ONLY: &
      InitializeReferenceElementX
    USE ReferenceElementModuleX_Lagrange, ONLY: &
      InitializeReferenceElementX_Lagrange
    USE GeometryComputationModule, ONLY: &
      ComputeGeometryX
    USE ReferenceElementModuleE, ONLY: &
      InitializeReferenceElementE
    USE ReferenceElementModuleE_Lagrange, ONLY: &
      InitializeReferenceElementE_Lagrange
    USE GeometryComputationModuleE, ONLY: &
      ComputeGeometryE
    USE ReferenceElementModuleZ, ONLY: &
      InitializeReferenceElementZ
    USE ReferenceElementModule, ONLY: &
      InitializeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      InitializeReferenceElement_Lagrange

    CALL InitializeTimers

    CALL InitializeProgram &
           ( ProgramName_Option &
               = TRIM( ProgramName ), &
             nX_Option &
               = nX, &
             swX_Option &
               = [ 1, 1, 1 ], &
             bcX_Option &
               = bcX, &
             xL_Option &
               = xL, &
             xR_Option &
               = xR, &
             zoomX_Option &
               = zoomX, &
             nE_Option &
               = nE, &
             swE_Option &
               = 1, &
             bcE_Option &
               = bcE, &
             eL_Option &
               = eL, &
             eR_Option &
               = eR, &
             zoomE_Option &
               = zoomE, &
             nNodes_Option &
               = nNodes, &
             CoordinateSystem_Option &
               = TRIM( CoordinateSystem ), &
             ActivateUnits_Option &
               = .TRUE., &
             nSpecies_Option &
               = nSpecies, &
             BasicInitialization_Option &
               = .TRUE. )

    ! --- Position Space Reference Element and Geometry ---

    CALL InitializeReferenceElementX

    CALL InitializeReferenceElementX_Lagrange

    CALL ComputeGeometryX &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    ! --- Energy Space Reference Element and Geometry ---

    CALL InitializeReferenceElementE

    CALL InitializeReferenceElementE_Lagrange

    CALL ComputeGeometryE &
           ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

    ! --- Phase Space Reference Element ---

    CALL InitializeReferenceElementZ

    CALL InitializeReferenceElement

    CALL InitializeReferenceElement_Lagrange

  END SUBROUTINE InitializeDriver


  SUBROUTINE FinalizeDriver

    USE ReferenceElementModuleX, ONLY: &
      FinalizeReferenceElementX
    USE ReferenceElementModuleX_Lagrange, ONLY: &
      FinalizeReferenceElementX_Lagrange
    USE ReferenceElementModuleE, ONLY: &
      FinalizeReferenceElementE
    USE ReferenceElementModuleE_Lagrange, ONLY: &
      FinalizeReferenceElementE_Lagrange
    USE ReferenceElementModuleZ, ONLY: &
      FinalizeReferenceElementZ
    USE ReferenceElementModule, ONLY: &
      FinalizeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      FinalizeReferenceElement_Lagrange
    USE ProgramInitializationModule, ONLY: &
      FinalizeProgram
    USE TwoMoment_TimersModule_OrderV, ONLY: &
      FinalizeTimers

    CALL FinalizeReferenceElementX

    CALL FinalizeReferenceElementX_Lagrange

    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementE_Lagrange

    CALL FinalizeReferenceElementZ

    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElement_Lagrange

    CALL FinalizeProgram

    CALL FinalizeTimers

  END SUBROUTINE FinalizeDriver


END PROGRAM ApplicationDriverNeutrinos
