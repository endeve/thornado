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
  INTEGER, PUBLIC, PARAMETER :: iGF_K_dd_11  = 14 ! Extrinsic Curvature 11
  INTEGER, PUBLIC, PARAMETER :: iGF_K_dd_12  = 15 ! Extrinsic Curvature 12
  INTEGER, PUBLIC, PARAMETER :: iGF_K_dd_13  = 16 ! Extrinsic Curvature 13
  INTEGER, PUBLIC, PARAMETER :: iGF_K_dd_22  = 17 ! Extrinsic Curvature 22
  INTEGER, PUBLIC, PARAMETER :: iGF_K_dd_23  = 18 ! Extrinsic Curvature 23
  INTEGER, PUBLIC, PARAMETER :: iGF_K_dd_33  = 19 ! Extrinsic Curvature 33
  INTEGER, PUBLIC, PARAMETER :: nGF          = 19 ! n Geometry Fields


  CHARACTER(32), DIMENSION(nGF), PUBLIC, PARAMETER :: &
    namesGF = [ 'Newtonian Potential             ', &
                'Spatial Scale Factor (1)        ', &
                'Spatial Scale Factor (2)        ', &
                'Spatial Scale Factor (3)        ', &
                'Spatial Metric Component (11)   ', &
                'Spatial Metric Component (22)   ', &
                'Spatial Metric Component (33)   ', &
                'Sqrt Spatial Metric Determinant ', &
                'Lapse Function                  ', &
                'Shift Vector (1)                ', &
                'Shift Vector (2)                ', &
                'Shift Vector (3)                ', &
                'Conformal Factor                ', &
                'Extrinsic Curvature Comp. (11)  ', &
                'Extrinsic Curvature Comp. (12)  ', &
                'Extrinsic Curvature Comp. (13)  ', &
                'Extrinsic Curvature Comp. (22)  ', &
                'Extrinsic Curvature Comp. (23)  ', &
                'Extrinsic Curvature Comp. (33)  ' ]

  CHARACTER(10), DIMENSION(nGF), PUBLIC, PARAMETER :: &
   ShortNamesGF = [ 'GF_Phi_N  ', &
                    'GF_h_1    ', &
                    'GF_h_2    ', &
                    'GF_h_3    ', &
                    'GF_Gm_11  ', &
                    'GF_Gm_22  ', &
                    'GF_Gm_33  ', &
                    'GF_SqrtGm ', &
                    'GF_Alpha  ', &
                    'GF_Beta_1 ', &
                    'GF_Beta_2 ', &
                    'GF_Beta_3 ', &
                    'GF_Psi    ', &
                    'GF_K_11   ', &
                    'GF_K_12   ', &
                    'GF_K_13   ', &
                    'GF_K_22   ', &
                    'GF_K_23   ', &
                    'GF_K_33   ' ]

  REAL(DP), DIMENSION(nGF), PUBLIC :: unitsGF
  REAL(DP), ALLOCATABLE,    PUBLIC :: uGF(:,:,:,:,:)

  PUBLIC :: CreateGeometryFields
  PUBLIC :: DestroyGeometryFields
  PUBLIC :: SetUnitsGeometryFields
  PUBLIC :: DescribeGeometryFields

CONTAINS


  SUBROUTINE CreateGeometryFields &
    ( nX, swX, CoordinateSystem_Option, Verbose_Option )

    INTEGER,      INTENT(in)           :: nX(3), swX(3)
    CHARACTER(*), INTENT(in), OPTIONAL :: CoordinateSystem_Option
    LOGICAL,      INTENT(in), OPTIONAL :: Verbose_Option

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

    CALL DescribeGeometryFields( CoordinateSystem, Verbose )

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
    uGF(:,:,:,:,iGF_K_dd_11)  = 0.0_DP
    uGF(:,:,:,:,iGF_K_dd_12)  = 0.0_DP
    uGF(:,:,:,:,iGF_K_dd_13)  = 0.0_DP
    uGF(:,:,:,:,iGF_K_dd_22)  = 0.0_DP
    uGF(:,:,:,:,iGF_K_dd_23)  = 0.0_DP
    uGF(:,:,:,:,iGF_K_dd_33)  = 0.0_DP

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uGF )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( uGF )
#endif

  END SUBROUTINE CreateGeometryFields


  SUBROUTINE DescribeGeometryFields( CoordinateSystem, Verbose )

    CHARACTER(*), INTENT(in) :: CoordinateSystem
    LOGICAL,      INTENT(in) :: Verbose

    INTEGER :: iGF

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

  END SUBROUTINE DescribeGeometryFields


  SUBROUTINE DestroyGeometryFields

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: uGF )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( uGF )
#endif

    DEALLOCATE( uGF )

  END SUBROUTINE DestroyGeometryFields


  SUBROUTINE SetUnitsGeometryFields

    USE UnitsModule, ONLY: &
      UnitsActive, &
      Erg, &
      Gram, &
      UnitsDisplay

    ASSOCIATE( U => UnitsDisplay )

    IF( UnitsActive )THEN

      unitsGF(iGF_Phi_N) = Erg / Gram

      IF     ( TRIM( CoordinateSystem ) .EQ. 'CARTESIAN' )THEN

        unitsGF(iGF_h_1) = 1.0_DP
        unitsGF(iGF_h_2) = 1.0_DP
        unitsGF(iGF_h_3) = 1.0_DP

        unitsGF(iGF_K_dd_11) = 1.0_DP / U % LengthX1Unit
        unitsGF(iGF_K_dd_12) = 1.0_DP / U % LengthX1Unit
        unitsGF(iGF_K_dd_13) = 1.0_DP / U % LengthX1Unit
        unitsGF(iGF_K_dd_22) = 1.0_DP / U % LengthX1Unit
        unitsGF(iGF_K_dd_23) = 1.0_DP / U % LengthX1Unit
        unitsGF(iGF_K_dd_33) = 1.0_DP / U % LengthX1Unit

      ELSE IF( TRIM( CoordinateSystem ) .EQ. 'CYLINDRICAL' )THEN

        unitsGF(iGF_h_1) = 1.0_DP
        unitsGF(iGF_h_2) = 1.0_DP
        unitsGF(iGF_h_3) = U % LengthX1Unit

        unitsGF(iGF_K_dd_11) = 1.0_DP / U % LengthX1Unit
        unitsGF(iGF_K_dd_12) = 1.0_DP / U % LengthX1Unit
        unitsGF(iGF_K_dd_13) = 1.0_DP
        unitsGF(iGF_K_dd_22) = 1.0_DP / U % LengthX1Unit
        unitsGF(iGF_K_dd_23) = 1.0_DP
        unitsGF(iGF_K_dd_33) = 1.0_DP

      ELSE IF( TRIM( CoordinateSystem ) .EQ. 'SPHERICAL' )THEN

        unitsGF(iGF_h_1) = 1.0_DP
        unitsGF(iGF_h_2) = U % LengthX1Unit
        unitsGF(iGF_h_3) = U % LengthX1Unit

        unitsGF(iGF_K_dd_11) = 1.0_DP / U % LengthX1Unit
        unitsGF(iGF_K_dd_12) = 1.0_DP
        unitsGF(iGF_K_dd_13) = 1.0_DP
        unitsGF(iGF_K_dd_22) = U % LengthX1Unit
        unitsGF(iGF_K_dd_23) = U % LengthX1Unit
        unitsGF(iGF_K_dd_33) = U % LengthX1Unit

      ELSE

        WRITE(*,*) 'Invalid choice of CoordinateSystem: ', &
                   TRIM( CoordinateSystem )

      END IF

      unitsGF(iGF_Gm_dd_11) = unitsGF(iGF_h_1)**2
      unitsGF(iGF_Gm_dd_22) = unitsGF(iGF_h_2)**2
      unitsGF(iGF_Gm_dd_33) = unitsGF(iGF_h_3)**2
      unitsGF(iGF_SqrtGm)   = unitsGF(iGF_h_1) &
                                * unitsGF(iGF_h_2) &
                                * unitsGF(iGF_h_3)
      unitsGF(iGF_Beta_1)   = U % VelocityX1Unit
      unitsGF(iGF_Beta_2)   = U % VelocityX2Unit
      unitsGF(iGF_Beta_3)   = U % VelocityX3Unit
      unitsGF(iGF_Alpha)    = 1.0_DP
      unitsGF(iGF_Psi)      = 1.0_DP

    ELSE

      unitsGF = 1.0_DP

    END IF

    END ASSOCIATE

  END SUBROUTINE SetUnitsGeometryFields


END MODULE GeometryFieldsModule
