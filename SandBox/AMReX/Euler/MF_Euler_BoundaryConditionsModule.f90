MODULE MF_Euler_BoundaryConditionsModule

  USE ISO_C_BINDING

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_Euler_ApplyBoundaryConditions

  INTERFACE
    SUBROUTINE MF_Euler_ApplyBoundaryConditions( pMF, pGEOM, BC ) BIND(c)
      IMPORT
      IMPLICIT NONE
      TYPE(c_ptr),    INTENT(inout) :: pMF
      TYPE(c_ptr),    INTENT(in)    :: pGEOM
      INTEGER(c_int), INTENT(in)    :: BC(*)
    END SUBROUTINE MF_Euler_ApplyBoundaryConditions
  END INTERFACE

END MODULE MF_Euler_BoundaryConditionsModule
