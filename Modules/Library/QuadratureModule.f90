MODULE QuadratureModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PUBLIC, DIMENSION(1:1) :: xG1, wG1 ! 1-Point Gauss
  REAL(DP), PUBLIC, DIMENSION(1:2) :: xG2, wG2 ! 2-Point Gauss
  REAL(DP), PUBLIC, DIMENSION(1:3) :: xG3, wG3 ! 3-Point Gauss
  REAL(DP), PUBLIC, DIMENSION(1:4) :: xG4, wG4 ! 4-Point Gauss
  REAL(DP), PUBLIC, DIMENSION(1:5) :: xG5, wG5 ! 5-Point Gauss

  REAL(DP), PUBLIC, DIMENSION(1:2) :: xL2, wL2 ! 2-Point Lobatto
  REAL(DP), PUBLIC, DIMENSION(1:3) :: xL3, wL3 ! 3-Point Lobatto
  REAL(DP), PUBLIC, DIMENSION(1:4) :: xL4, wL4 ! 4-Point Lobatto
  REAL(DP), PUBLIC, DIMENSION(1:5) :: xL5, wL5 ! 5-Point Lobatto

  PUBLIC :: &
    InitializeQuadratures, &
    GetQuadrature

CONTAINS


  SUBROUTINE InitializeQuadratures

    ! Quadratures Defined on [-0.5, 0.5]

    !**********************************************
    ! 1-Point Gaussian Quadrature
    !**********************************************

    xG1(1) = 0.0_DP

    wG1(1) = 1.0_DP

    !**********************************************
    ! 2-Point Gaussian Quadrature
    !**********************************************

    xG2(1) = - SQRT( 1.0_DP / 12.0_DP )
    xG2(2) = + SQRT( 1.0_DP / 12.0_DP )

    wG2(1) = 0.5_DP
    wG2(2) = 0.5_DP

    !**********************************************
    ! 3-Point Gaussian Quadrature
    !**********************************************

    xG3(1) = - SQRT( 15.0_DP ) / 10.0_DP
    xG3(2) = 0.0_DP
    xG3(3) = + SQRT( 15.0_DP ) / 10.0_DP

    wG3(1) = 5.0_DP / 18.0_DP
    wG3(2) = 8.0_DP / 18.0_DP
    wG3(3) = 5.0_DP / 18.0_DP

    !**********************************************
    ! 4-Point Gaussian Quadrature
    !**********************************************

    xG4(1) = - SQRT( (3.0_DP+2.0_DP*SQRT(6.0_DP/5.0_DP)) / 7.0_DP ) / 2.0_DP
    xG4(2) = - SQRT( (3.0_DP-2.0_DP*SQRT(6.0_DP/5.0_DP)) / 7.0_DP ) / 2.0_DP
    xG4(3) = + SQRT( (3.0_DP-2.0_DP*SQRT(6.0_DP/5.0_DP)) / 7.0_DP ) / 2.0_DP
    xG4(4) = + SQRT( (3.0_DP+2.0_DP*SQRT(6.0_DP/5.0_DP)) / 7.0_DP ) / 2.0_DP

    wG4(1) = ( 18.0_DP - SQRT( 30.0_DP ) ) / 72.0_DP
    wG4(2) = ( 18.0_DP + SQRT( 30.0_DP ) ) / 72.0_DP
    wG4(3) = ( 18.0_DP + SQRT( 30.0_DP ) ) / 72.0_DP
    wG4(4) = ( 18.0_DP - SQRT( 30.0_DP ) ) / 72.0_DP

    !**********************************************
    ! 5-Point Gaussian Quadrature
    !**********************************************

    xG5(1) = - SQRT( 245.0_DP + 14.0_DP * SQRT( 70.0_DP ) ) / 21.0_DP / 2.0_DP
    xG5(2) = - SQRT( 245.0_DP - 14.0_DP * SQRT( 70.0_DP ) ) / 21.0_DP / 2.0_DP
    xG5(3) = 0.0_DP
    xG5(4) = - xG5(2)
    xG5(5) = - xG5(1)

    wG5(1) = ( 322.0_DP - 13.0_DP * SQRT( 70.0_DP ) ) / 900.0_DP / 2.0_DP
    wG5(2) = ( 322.0_DP + 13.0_DP * SQRT( 70.0_DP ) ) / 900.0_DP / 2.0_DP
    wG5(3) = 128.0_DP / 225.0_DP / 2.0_DP
    wG5(4) = wG5(2)
    wG5(5) = wG5(1)

    !*************************************************
    ! 2-Point Gauss-Lobatto-Quadrature
    !*************************************************

    xL2(1) = - 0.5_DP
    xL2(2) = + 0.5_DP

    wL2(1) = 0.5_DP
    wL2(2) = 0.5_DP

    !*************************************************
    ! 3-Point Gauss-Lobatto-Quadrature
    !*************************************************

    xL3(1) = - 0.5_DP
    xL3(2) =   0.0_DP
    xL3(3) =   0.5_DP

    wL3(1) = 1.0_DP / 6.0_DP
    wL3(2) = 2.0_DP / 3.0_DP
    wL3(3) = 1.0_DP / 6.0_DP

    !*************************************************
    ! 4-Point Gauss-Lobatto-Quadrature
    !*************************************************

    xL4(1) = - 0.5_DP
    xL4(2) = - SQRT( 5.0_DP ) / 10.0_DP
    xL4(3) =   SQRT( 5.0_DP ) / 10.0_DP
    xL4(4) =   0.5_DP

    wL4(1) = 1.0_DP / 12.0_DP
    wL4(2) = 5.0_DP / 12.0_DP
    wL4(3) = 5.0_DP / 12.0_DP
    wL4(4) = 1.0_DP / 12.0_DP

    !*************************************************
    ! 5-Point Gauss-Lobatto-Quadrature
    !*************************************************

    xL5(1) = - 0.5_DP
    xL5(2) = - SQRT( 21.0_DP ) / 14.0_DP
    xL5(3) =   0.0_DP
    xL5(4) =   SQRT( 21.0_DP ) / 14.0_DP
    xL5(5) =   0.5_DP

    wL5(1) =  1.0_DP /  20.0_DP
    wL5(2) = 49.0_DP / 180.0_DP
    wL5(3) = 32.0_DP /  90.0_DP
    wL5(4) = 49.0_DP / 180.0_DP
    wL5(5) =  1.0_DP /  20.0_DP

  END SUBROUTINE InitializeQuadratures


  SUBROUTINE GetQuadrature( nQ, xQ, wQ, QuadratureNameOption )

    INTEGER, INTENT(in)                    :: nQ
    REAL(DP), DIMENSION(nQ), INTENT(inout) :: xQ, wQ
    CHARACTER(32), INTENT(in), OPTIONAL    :: QuadratureNameOption

    CHARACTER(32) :: QuadratureName

    CALL InitializeQuadratures

    QuadratureName = 'Gaussian'
    IF( PRESENT( QuadratureNameOption ) ) &
      QuadratureName = TRIM( QuadratureNameOption )

    SELECT CASE ( TRIM( QuadratureName ) )
      CASE ( 'Gaussian' )
        SELECT CASE ( nQ )
          CASE ( 1 )
            xQ = xG1; wQ = wG1
          CASE ( 2 )
            xQ = xG2; wQ = wG2
          CASE ( 3 )
            xQ = xG3; wQ = wG3
          CASE ( 4 )
            xQ = xG4; wQ = wG4
          CASE ( 5 )
            xQ = xG5; wQ = wG5
          CASE DEFAULT
            WRITE(*,*)
            WRITE(*,'(A5,A45,I2.2)') &
              '', 'Gaussian Quadrature Not Implemented for nQ = ', nQ
            STOP
        END SELECT
      CASE ( 'Lobatto' )
        SELECT CASE ( nQ )
          CASE ( 2 )
            xQ = xL2; wQ = wL2
          CASE ( 3 )
            xQ = xL3; wQ = wL3
          CASE ( 4 )
            xQ = xL4; wQ = wL4
          CASE ( 5 )
            xQ = xL5; wQ = wL5
          CASE DEFAULT
            WRITE(*,*)
            WRITE(*,'(A5,A44,I2.2)') &
              '', 'Lobatto Quadrature Not Implemented for nQ = ', nQ
            STOP
         END SELECT
      CASE DEFAULT
        WRITE(*,*)
        WRITE(*,'(A5,A28,A)') &
          '', 'Quadrature Not Implemented: ', TRIM( QuadratureName )
    END SELECT

  END SUBROUTINE GetQuadrature


END MODULE QuadratureModule
