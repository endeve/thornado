MODULE GeometryFieldsModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX

  IMPLICIT NONE
  PRIVATE

  CHARACTER(16), PUBLIC :: &
    CoordinateSystem = 'CARTESIAN'

  LOGICAL :: Verbose

  ! --- Spatial Geometry Fields ---

  INTEGER, PUBLIC, PARAMETER :: iGF_Phi_N    = 1  ! Newtonian Potential
  INTEGER, PUBLIC, PARAMETER :: iGF_h_1      = 2  ! Spatial Scale Factor 1
  INTEGER, PUBLIC, PARAMETER :: iGF_h_2      = 3  ! Spatial Scale Factor 2
  INTEGER, PUBLIC, PARAMETER :: iGF_h_3      = 4  ! Spatial Scale Factor 3
  INTEGER, PUBLIC, PARAMETER :: iGF_Gm_dd_11 = 5  ! Spatial Metric Component 11
  INTEGER, PUBLIC, PARAMETER :: iGF_Gm_dd_22 = 6  ! Spatial Metric Component 22
  INTEGER, PUBLIC, PARAMETER :: iGF_Gm_dd_33 = 7  ! Spatial Metric Component 33
  INTEGER, PUBLIC, PARAMETER :: iGF_SqrtGm   = 8  ! Sqrt of Metric Determinant
  INTEGER, PUBLIC, PARAMETER :: iGF_Alpha    = 9  ! Lapse Function
  INTEGER, PUBLIC, PARAMETER :: iGF_Beta_1   = 10 ! Shift Vector 1
  INTEGER, PUBLIC, PARAMETER :: iGF_Beta_2   = 11 ! Shift Vector 2
  INTEGER, PUBLIC, PARAMETER :: iGF_Beta_3   = 12 ! Shift Vector 3
  INTEGER, PUBLIC, PARAMETER :: iGF_Psi      = 13 ! Conformal Factor
  INTEGER, PUBLIC, PARAMETER :: nGF          = 13 ! n Geometry Fields

  CHARACTER(32), DIMENSION(nGF), PUBLIC, PARAMETER :: &
    namesGF = [ 'Newtonian Potential                         ', &
                'Spatial Scale Factor (1)                    ', &
                'Spatial Scale Factor (2)                    ', &
                'Spatial Scale Factor (3)                    ', &
                'Spatial Metric Component (11)               ', &
                'Spatial Metric Component (22)               ', &
                'Spatial Metric Component (33)               ', &
                'Sqrt of Spatial Metric Determinant          ', &
                'Lapse Function                              ', &
                'Shift Vector (1)                            ', &
                'Shift Vector (2)                            ', &
                'Shift Vector (3)                            ', &
                'Conformal Factor                            ' ]

  REAL(DP), DIMENSION(nGF), PUBLIC :: unitsGF
  REAL(DP), ALLOCATABLE,    PUBLIC :: uGF(:,:,:,:,:)

  PUBLIC :: CreateGeometryFields
  PUBLIC :: DestroyGeometryFields

CONTAINS


  SUBROUTINE CreateGeometryFields &
    ( nX, swX, CoordinateSystem_Option, Verbose_Option )

    INTEGER,      INTENT(in)           :: nX(3), swX(3)
    CHARACTER(*), INTENT(in), OPTIONAL :: CoordinateSystem_Option
    LOGICAL,      INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iGF

    IF( PRESENT( CoordinateSystem_Option ) )THEN
      CoordinateSystem = CoordinateSystem_Option
    ELSE
      CoordinateSystem = 'CARTESIAN'
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A15)') '', 'Geometry Fields'
      WRITE(*,*)
      WRITE(*,'(A5,A,A)') &
        '', 'Coordinate System: ', TRIM( CoordinateSystem )
      WRITE(*,*)
      DO iGF = 1, nGF
        WRITE(*,'(A5,A32)') '', TRIM( namesGF(iGF) )
      END DO
    END IF

    ALLOCATE &
      ( uGF(1:nDOFX, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nGF) )

    CALL SetUnitsGeometryFields

    ! --- Initialize to Flat Spacetime (Cartesian) ---

    uGF(:,:,:,:,iGF_Phi_N)    = 0.0_DP
    uGF(:,:,:,:,iGF_h_1)      = 1.0_DP
    uGF(:,:,:,:,iGF_h_2)      = 1.0_DP
    uGF(:,:,:,:,iGF_h_3)      = 1.0_DP
    uGF(:,:,:,:,iGF_Gm_dd_11) = 1.0_DP
    uGF(:,:,:,:,iGF_Gm_dd_22) = 1.0_DP
    uGF(:,:,:,:,iGF_Gm_dd_33) = 1.0_DP
    uGF(:,:,:,:,iGF_SqrtGm)   = 1.0_DP
    uGF(:,:,:,:,iGF_Alpha)    = 1.0_DP
    uGF(:,:,:,:,iGF_Beta_1)   = 0.0_DP
    uGF(:,:,:,:,iGF_Beta_2)   = 0.0_DP
    uGF(:,:,:,:,iGF_Beta_3)   = 0.0_DP
    uGF(:,:,:,:,iGF_Psi)      = 1.0_DP

  END SUBROUTINE CreateGeometryFields


  SUBROUTINE DestroyGeometryFields

    DEALLOCATE( uGF )

  END SUBROUTINE DestroyGeometryFields


  SUBROUTINE SetUnitsGeometryFields

    USE UnitsModule, ONLY: &
      UnitsActive, &
      Erg, &
      Gram

    IF( UnitsActive )THEN

      unitsGF(iGF_Phi_N)    = Erg / Gram
      unitsGF(iGF_h_1)      = 1.0_DP
      unitsGF(iGF_h_2)      = 1.0_DP
      unitsGF(iGF_h_3)      = 1.0_DP
      unitsGF(iGF_Gm_dd_11) = 1.0_DP
      unitsGF(iGF_Gm_dd_22) = 1.0_DP
      unitsGF(iGF_Gm_dd_33) = 1.0_DP
      unitsGF(iGF_SqrtGm)   = 1.0_DP
      unitsGF(iGF_Alpha)    = 1.0_DP
      unitsGF(iGF_Beta_1)   = 1.0_DP
      unitsGF(iGF_Beta_2)   = 1.0_DP
      unitsGF(iGF_Beta_3)   = 1.0_DP
      unitsGF(iGF_Psi)      = 1.0_DP

    ELSE

      unitsGF = 1.0_DP

    END IF

  END SUBROUTINE SetUnitsGeometryFields


END MODULE GeometryFieldsModule
