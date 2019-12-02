MODULE GravitySolutionModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX, &
    iX_B0, iX_E0, iX_B1, iX_E1
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumberX_X1
  USE FluidFieldsModule, ONLY: &
    uCF
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE PolynomialBasisModule_Lagrange, ONLY: &
     L_X1,  L_X2,  L_X3, &
    dL_X1, dL_X2, dL_X3
  USE MeshModule, ONLY: &
    MeshX
  USE GravitySolutionModule_Newtonian_PointMass_Old, ONLY: &
    InitializeGravitySolver_Newtonian_PointMass, &
    FinalizeGravitySolver_Newtonian_PointMass, &
    SolveGravity_Newtonian_PointMass
  USE GravitySolutionModule_Newtonian_Poseidon_Old, ONLY: &
    InitializeGravitySolver_Newtonian_Poseidon, &
    FinalizeGravitySolver_Newtonian_Poseidon, &
    SolveGravity_Newtonian_Poseidon

  IMPLICIT NONE
  PRIVATE

  CHARACTER(32) :: GravitySolver
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X1_L ! Basis Function
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X1_H ! Basis Function
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: V_X1_L ! Test  Function
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: V_X1_H ! Test  Function
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dVdX1  ! Test  Derivatives
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dVdX2  ! Test  Derivatives
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dVdX3  ! Test  Derivatives

  PROCEDURE (GS), POINTER, PUBLIC :: &
    SolveGravity => NULL()

  INTERFACE
    SUBROUTINE GS
    END SUBROUTINE GS
  END INTERFACE

  PUBLIC :: InitializeGravitySolver
  PUBLIC :: FinalizeGravitySolver

CONTAINS


  SUBROUTINE InitializeGravitySolver( GravitySolver_Option, PointMass_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: GravitySolver_Option
    REAL(DP),         INTENT(in), OPTIONAL :: PointMass_Option

    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    INTEGER :: jNodeX1, jNodeX2, jNodeX3, jNodeX

    GravitySolver = 'Dummy'
    IF( PRESENT( GravitySolver_Option ) ) &
      GravitySolver = GravitySolver_Option

    WRITE(*,*)
    WRITE(*,'(A5,A16,A)') &
      '', 'Gravity Solver: ', TRIM( GravitySolver )
    WRITE(*,'(A5,A16)') &
      '', '--------------- '

    SELECT CASE ( TRIM( GravitySolver ) )

      CASE ( 'Newtonian_PointMass' )

        CALL InitializeGravitySolver_Newtonian_PointMass &
               ( PointMass_Option = PointMass_Option )
        SolveGravity &
          => SolveGravity_Newtonian_PointMass

      CASE ( 'Newtonian_Poseidon' )

        CALL InitializeGravitySolver_Newtonian_Poseidon
        SolveGravity &
          => SolveGravity_Newtonian_Poseidon

      CASE DEFAULT

        SolveGravity &
          => SolveGravity_Dummy

    END SELECT

    ALLOCATE &
      ( dVdX1(nDOFX,nDOFX), &
        dVdX2(nDOFX,nDOFX), &
        dVdX3(nDOFX,nDOFX) )

    DO jNodeX3 = 1, nNodesX(3)
      DO jNodeX2 = 1, nNodesX(2)
        DO jNodeX1 = 1, nNodesX(1)

          jNodeX = NodeNumberX( jNodeX1, jNodeX2, jNodeX3 )

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                dVdX1(iNodeX,jNodeX) &
                  = dL_X1(jNodeX1) % P( MeshX(1) % Nodes(iNodeX1) ) &
                      * L_X2(jNodeX2) % P( MeshX(2) % Nodes(iNodeX2) ) &
                      * L_X3(jNodeX3) % P( MeshX(3) % Nodes(iNodeX3) )

                dVdX2(iNodeX,jNodeX) &
                  = dL_X2(jNodeX2) % P( MeshX(2) % Nodes(iNodeX2) ) &
                      * L_X1(jNodeX1) % P( MeshX(1) % Nodes(iNodeX1) ) &
                      * L_X3(jNodeX3) % P( MeshX(3) % Nodes(iNodeX3) )

                dVdX3(iNodeX,jNodeX) &
                  = dL_X3(jNodeX3) % P( MeshX(3) % Nodes(iNodeX3) ) &
                      * L_X1(jNodeX1) % P( MeshX(1) % Nodes(iNodeX1) ) &
                      * L_X2(jNodeX2) % P( MeshX(2) % Nodes(iNodeX2) )

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

    ALLOCATE &
      ( L_X1_L(nDOFX,nNodesX(2)*nNodesX(3)), &
        L_X1_H(nDOFX,nNodesX(2)*nNodesX(3)) )

    ALLOCATE &
      ( V_X1_L(nNodesX(2)*nNodesX(3),nDOFX), &
        V_X1_H(nNodesX(2)*nNodesX(3),nDOFX) )

    DO jNodeX3 = 1, nNodesX(3)
      DO jNodeX2 = 1, nNodesX(2)

        jNodeX1 = NodeNumberX_X1( jNodeX2, jNodeX3 )

        DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
            DO iNodeX1 = 1, nNodesX(1)

              iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

              L_X1_L(iNodeX,jNodeX1) &
                = L_X1(iNodeX1) % P( - 0.5_DP ) &
                    * L_X2(iNodeX2) % P( MeshX(2) % Nodes(jNodeX2) ) &
                    * L_X3(iNodeX3) % P( MeshX(3) % Nodes(jNodeX3) )

              V_X1_L(jNodeX1,iNodeX) &
                = L_X1_L(iNodeX,jNodeX1)

              L_X1_H(iNodeX,jNodeX1) &
                = L_X1(iNodeX1) % P( + 0.5_DP ) &
                    * L_X2(iNodeX2) % P( MeshX(2) % Nodes(jNodeX2) ) &
                    * L_X3(iNodeX3) % P( MeshX(3) % Nodes(jNodeX3) )

              V_X1_H(jNodeX1,iNodeX) &
                = L_X1_H(iNodeX,jNodeX1)

            END DO
          END DO
        END DO

      END DO
    END DO

  END SUBROUTINE InitializeGravitySolver


  SUBROUTINE FinalizeGravitySolver

    SELECT CASE ( TRIM( GravitySolver ) )

      CASE ( 'Newtonian_PointMass' )

        CALL FinalizeGravitySolver_Newtonian_PointMass

      CASE ( 'Newtonian_Poseidon' )

        CALL FinalizeGravitySolver_Newtonian_Poseidon

      CASE DEFAULT

    END SELECT

    NULLIFY( SolveGravity )

    DEALLOCATE( L_X1_L, L_X1_H )
    DEALLOCATE( V_X1_L, V_X1_H )
    DEALLOCATE( dVdX1, dVdX2, dVdX3 )

  END SUBROUTINE FinalizeGravitySolver


  SUBROUTINE SolveGravity_Dummy

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'GravitySolutionModule: SolveGravity_Dummy'
    WRITE(*,*)

  END SUBROUTINE SolveGravity_Dummy


END MODULE GravitySolutionModule
