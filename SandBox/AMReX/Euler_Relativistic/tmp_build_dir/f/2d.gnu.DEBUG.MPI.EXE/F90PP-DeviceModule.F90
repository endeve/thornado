MODULE DeviceModule

  USE, INTRINSIC :: ISO_C_BINDING
  USE KindModule, ONLY: &
    DP



































  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  INTEGER, PUBLIC :: mydevice, ndevices

  INTERFACE QueryOnGPU
    MODULE PROCEDURE QueryOnGPU_3D_DP_1
    MODULE PROCEDURE QueryOnGPU_3D_DP_2
    MODULE PROCEDURE QueryOnGPU_3D_DP_3
    MODULE PROCEDURE QueryOnGPU_2D_DP_1
    MODULE PROCEDURE QueryOnGPU_DP_1
    MODULE PROCEDURE QueryOnGPU_DP_2
    MODULE PROCEDURE QueryOnGPU_DP_3
    MODULE PROCEDURE QueryOnGPU_DP_4
    MODULE PROCEDURE QueryOnGPU_DP_5
    MODULE PROCEDURE QueryOnGPU_DP_6
    MODULE PROCEDURE QueryOnGPU_DP_7
  END INTERFACE

  PUBLIC :: InitializeDevice
  PUBLIC :: FinalizeDevice
  PUBLIC :: device_is_present
  PUBLIC :: get_device_num
  PUBLIC :: QueryOnGpu

CONTAINS


  SUBROUTINE InitializeDevice

    INTEGER :: ierr, myrank, nranks


    mydevice = -1
    ndevices = 0

























    RETURN
  END SUBROUTINE InitializeDevice


  SUBROUTINE FinalizeDevice
    RETURN
  END SUBROUTINE FinalizeDevice


  LOGICAL FUNCTION device_is_present( hostptr, device, bytes )
    TYPE(C_PTR), INTENT(in) :: hostptr
    INTEGER, INTENT(in) :: device
    INTEGER(C_SIZE_T), INTENT(in) :: bytes





    device_is_present = .false.

    RETURN
  END FUNCTION device_is_present


  INTEGER FUNCTION get_device_num()





    get_device_num = -1

    RETURN
  END FUNCTION get_device_num


  FUNCTION QueryOnGPU_3D_DP_1( X1 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET :: X1
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 )

  END FUNCTION QueryOnGPU_3D_DP_1


  FUNCTION QueryOnGPU_3D_DP_2( X1, X2 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET :: X1, X2
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 )

  END FUNCTION QueryOnGPU_3D_DP_2


  FUNCTION QueryOnGPU_3D_DP_3( X1, X2, X3 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET :: X1, X2, X3
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 )

  END FUNCTION QueryOnGPU_3D_DP_3


  FUNCTION QueryOnGPU_2D_DP_1( X1 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:,:), INTENT(in), TARGET :: X1
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 )

  END FUNCTION QueryOnGPU_2D_DP_1


  FUNCTION QueryOnGPU_DP_1( X1 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 )

  END FUNCTION QueryOnGPU_DP_1


  FUNCTION QueryOnGPU_DP_2( X1, X2 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 )

  END FUNCTION QueryOnGPU_DP_2


  FUNCTION QueryOnGPU_DP_3( X1, X2, X3 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2, X3
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 )

  END FUNCTION QueryOnGPU_DP_3


  FUNCTION QueryOnGPU_DP_4( X1, X2, X3, X4 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2, X3, X4
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3
    INTEGER(C_SIZE_T) :: SizeOf_X4

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)
    SizeOf_X4 = SIZE(X4) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 ) &
           .AND. device_is_present( C_LOC( X4 ), mydevice, SizeOf_X4 )

  END FUNCTION QueryOnGPU_DP_4


  FUNCTION QueryOnGPU_DP_5( X1, X2, X3, X4, X5 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2, X3, X4, X5
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3
    INTEGER(C_SIZE_T) :: SizeOf_X4
    INTEGER(C_SIZE_T) :: SizeOf_X5

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)
    SizeOf_X4 = SIZE(X4) * C_SIZEOF(0.0_DP)
    SizeOf_X5 = SIZE(X5) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 ) &
           .AND. device_is_present( C_LOC( X4 ), mydevice, SizeOf_X4 ) &
           .AND. device_is_present( C_LOC( X5 ), mydevice, SizeOf_X5 )

  END FUNCTION QueryOnGPU_DP_5


  FUNCTION QueryOnGPU_DP_6( X1, X2, X3, X4, X5, X6 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2, X3, X4, X5, X6
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3
    INTEGER(C_SIZE_T) :: SizeOf_X4
    INTEGER(C_SIZE_T) :: SizeOf_X5
    INTEGER(C_SIZE_T) :: SizeOf_X6

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)
    SizeOf_X4 = SIZE(X4) * C_SIZEOF(0.0_DP)
    SizeOf_X5 = SIZE(X5) * C_SIZEOF(0.0_DP)
    SizeOf_X6 = SIZE(X6) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 ) &
           .AND. device_is_present( C_LOC( X4 ), mydevice, SizeOf_X4 ) &
           .AND. device_is_present( C_LOC( X5 ), mydevice, SizeOf_X5 ) &
           .AND. device_is_present( C_LOC( X6 ), mydevice, SizeOf_X6 )

  END FUNCTION QueryOnGPU_DP_6


  FUNCTION QueryOnGPU_DP_7( X1, X2, X3, X4, X5, X6, X7 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2, X3, X4, X5, X6, X7
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3
    INTEGER(C_SIZE_T) :: SizeOf_X4
    INTEGER(C_SIZE_T) :: SizeOf_X5
    INTEGER(C_SIZE_T) :: SizeOf_X6
    INTEGER(C_SIZE_T) :: SizeOf_X7

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)
    SizeOf_X4 = SIZE(X4) * C_SIZEOF(0.0_DP)
    SizeOf_X5 = SIZE(X5) * C_SIZEOF(0.0_DP)
    SizeOf_X6 = SIZE(X6) * C_SIZEOF(0.0_DP)
    SizeOf_X7 = SIZE(X7) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 ) &
           .AND. device_is_present( C_LOC( X4 ), mydevice, SizeOf_X4 ) &
           .AND. device_is_present( C_LOC( X5 ), mydevice, SizeOf_X5 ) &
           .AND. device_is_present( C_LOC( X6 ), mydevice, SizeOf_X6 ) &
           .AND. device_is_present( C_LOC( X7 ), mydevice, SizeOf_X7 )

  END FUNCTION QueryOnGPU_DP_7


END MODULE DeviceModule

