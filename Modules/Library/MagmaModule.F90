!
!   -- MAGMA (version 2.3.0) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      @date November 2017
!
! JAH: this file combines magma2.f90, magma2_common.f90, and magma2_dfortran.f90 from magma 2.3.0

module MagmaModule

use iso_c_binding

implicit none

integer(c_int) :: magma_device

integer :: mymagma_queue
type(c_ptr) :: magma_queue
!$omp threadprivate(magma_queue,mymagma_queue)

integer :: nmagma_queue
type(c_ptr) :: magma_queue_array_cptr
type(c_ptr), pointer :: magma_queue_array(:)

!! =============================================================================
!! Parameter constants from magma_types.h
integer(c_int), parameter ::   &
    MagmaFalse         = 0,    &
    MagmaTrue          = 1,    &
    
    MagmaRowMajor      = 101,  &
    MagmaColMajor      = 102,  &
    
    MagmaNoTrans       = 111,  &
    MagmaTrans         = 112,  &
    MagmaConjTrans     = 113,  &
    
    MagmaUpper         = 121,  &
    MagmaLower         = 122,  &
    MagmaGeneral       = 123,  &
    MagmaFull          = 123,  &  !! deprecated, use MagmaGeneral
    
    MagmaNonUnit       = 131,  &
    MagmaUnit          = 132,  &
    
    MagmaLeft          = 141,  &
    MagmaRight         = 142,  &
    MagmaBothSides     = 143
!! todo all the rest

!! =====================================================================
!! Parameter constants
real(c_float),             parameter :: sdummy = 0
real(c_double),            parameter :: ddummy = 0
complex(c_float_complex),  parameter :: cdummy = 0
complex(c_double_complex), parameter :: zdummy = 0
integer(c_int),            parameter :: idummy = 0
type(c_ptr),               parameter :: ptr_dummy = c_null_ptr

!! Intel ifort chokes on c_sizeof here, so use extension sizeof
!! see https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/495001
integer(c_size_t), parameter :: &
    sizeof_real      = sizeof(sdummy), &
    sizeof_double    = sizeof(ddummy), &
    sizeof_complex   = sizeof(cdummy), &
    sizeof_complex16 = sizeof(zdummy), &
    sizeof_int       = sizeof(idummy)!, &
    !sizeof_ptr       = sizeof(ptr_dummy)


!! =============================================================================
!! Fortran interfaces to C functions
interface

    !! -------------------------------------------------------------------------
    !! initialize
    subroutine magma_init() &
    bind(C, name="magma_init")
        use iso_c_binding
    end subroutine

    subroutine magma_finalize() &
    bind(C, name="magma_finalize")
        use iso_c_binding
    end subroutine

    !! -------------------------------------------------------------------------
    !! version
    subroutine magma_version( major, minor, micro ) &
    bind(C, name="magma_version")
        use iso_c_binding
        integer(c_int), target :: major, minor, micro
    end subroutine

    subroutine magma_print_environment() &
    bind(C, name="magma_print_environment")
        use iso_c_binding
    end subroutine

    !! -------------------------------------------------------------------------
    !! timing
    real(c_double) function magma_wtime() &
    bind(C, name="magma_wtime")
        use iso_c_binding
    end function

    real(c_double) function magma_sync_wtime( queue ) &
    bind(C, name="magma_sync_wtime")
        use iso_c_binding
        type(c_ptr), value :: queue
    end function

    !! -------------------------------------------------------------------------
    !! device support
    integer(c_int) function magma_num_gpus() &
    bind(C, name="magma_num_gpus")
        use iso_c_binding
    end function

    integer(c_int) function magma_getdevice_arch() &
    bind(C, name="magma_getdevice_arch")
        use iso_c_binding
    end function

    subroutine magma_getdevice( dev ) &
    bind(C, name="magma_getdevice")
        use iso_c_binding
        integer(c_int), target :: dev
    end subroutine

    subroutine magma_setdevice( dev ) &
    bind(C, name="magma_setdevice")
        use iso_c_binding
        integer(c_int), value :: dev
    end subroutine

    integer(c_size_t) function magma_mem_size( queue ) &
    bind(C, name="magma_mem_size")
        use iso_c_binding
        type(c_ptr), value :: queue
    end function

    !! -------------------------------------------------------------------------
    !! queue support
    subroutine magma_queue_create_internal( dev, queue_ptr, func, file, line ) &
    bind(C, name="magma_queue_create_internal")
        use iso_c_binding
        integer(c_int), value :: dev
        type(c_ptr), target :: queue_ptr  !! queue_t*
        character(c_char) :: func, file
        integer(c_int), value :: line
    end subroutine

    subroutine magma_queue_destroy_internal( queue, func, file, line ) &
    bind(C, name="magma_queue_destroy_internal")
        use iso_c_binding
        type(c_ptr), value :: queue  !! queue_t
        character(c_char) :: func, file
        integer(c_int), value :: line
    end subroutine

#if defined(THORNADO_CUDA)
    subroutine magma_queue_create_from_cuda_internal( dev, stream, cublas_handle, &
    cusparse_handle, queue_ptr, func, file, line ) &
    bind(C, name="magma_queue_create_from_cuda_internal")
      use iso_c_binding
      integer(c_int), value :: dev
      type(c_ptr), value :: stream, cublas_handle, cusparse_handle
      type(c_ptr), target :: queue_ptr  !! queue_t*
      character(c_char) :: func, file
      integer(c_int), value :: line
    end subroutine magma_queue_create_from_cuda_internal
#elif defined(THORNADO_HIP)
    subroutine magma_queue_create_from_hip_internal( dev, stream, hipblas_handle, &
    hipsparse_handle, queue_ptr, func, file, line ) &
    bind(C, name="magma_queue_create_from_hip_internal")
      use iso_c_binding
      integer(c_int), value :: dev
      type(c_ptr), value :: stream, hipblas_handle, hipsparse_handle
      type(c_ptr), target :: queue_ptr  !! queue_t*
      character(c_char) :: func, file
      integer(c_int), value :: line
    end subroutine magma_queue_create_from_hip_internal
#endif

    subroutine magma_queue_sync_internal( queue, func, file, line ) &
    bind(C, name="magma_queue_sync_internal")
        use iso_c_binding
        type(c_ptr), value :: queue  !! queue_t
        character(c_char) :: func, file
        integer(c_int), value :: line
    end subroutine

    integer(c_int) function magma_queue_get_device( queue ) &
    bind(C, name="magma_queue_get_device")
        use iso_c_binding
        type(c_ptr), value :: queue  !! queue_t
    end function

    !! -------------------------------------------------------------------------
    !! offsets pointers -- 1D vectors with inc
    !! see offset.c
    type(c_ptr) function magma_soffset_1d( ptr, inc, i ) &
    bind(C, name="magma_soffset_1d")
        use iso_c_binding
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function

    type(c_ptr) function magma_doffset_1d( ptr, inc, i ) &
    bind(C, name="magma_doffset_1d")
        use iso_c_binding
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function

    type(c_ptr) function magma_coffset_1d( ptr, inc, i ) &
    bind(C, name="magma_coffset_1d")
        use iso_c_binding
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function

    type(c_ptr) function magma_zoffset_1d( ptr, inc, i ) &
    bind(C, name="magma_zoffset_1d")
        use iso_c_binding
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function

    type(c_ptr) function magma_ioffset_1d( ptr, inc, i ) &
    bind(C, name="magma_ioffset_1d")
        use iso_c_binding
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function

    !! -------------------------------------------------------------------------
    !! offsets pointers -- 2D matrices with lda
    !! see offset.c
    type(c_ptr) function magma_soffset_2d( ptr, lda, i, j ) &
    bind(C, name="magma_soffset_2d")
        use iso_c_binding
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function

    type(c_ptr) function magma_doffset_2d( ptr, lda, i, j ) &
    bind(C, name="magma_doffset_2d")
        use iso_c_binding
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function

    type(c_ptr) function magma_coffset_2d( ptr, lda, i, j ) &
    bind(C, name="magma_coffset_2d")
        use iso_c_binding
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function

    type(c_ptr) function magma_zoffset_2d( ptr, lda, i, j ) &
    bind(C, name="magma_zoffset_2d")
        use iso_c_binding
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function

    type(c_ptr) function magma_ioffset_2d( ptr, lda, i, j ) &
    bind(C, name="magma_ioffset_2d")
        use iso_c_binding
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function

    !! -------------------------------------------------------------------------
    !! magma_malloc (GPU memory)
    integer(c_int) function magma_malloc( ptr, bytes ) &
    bind(C, name="magma_malloc")
        use iso_c_binding
        type(c_ptr), target :: ptr  !! void**
        integer(c_size_t), value :: bytes
    end function

    !! todo imalloc

    integer(c_int) function magma_free_internal( ptr, func, file, line ) &
    bind(C, name="magma_free_internal")
        use iso_c_binding
        type(c_ptr), value :: ptr  !! void*
        character(c_char) :: func, file
        integer(c_int), value :: line
    end function

    !! -------------------------------------------------------------------------
    !! magma_malloc_cpu (CPU main memory)
    !! these are aligned to 32-byte boundary
    integer(c_int) function magma_malloc_cpu( ptr, bytes ) &
    bind(C, name="magma_malloc_cpu")
        use iso_c_binding
        type(c_ptr), target :: ptr  !! void**
        integer(c_size_t), value :: bytes
    end function

    !! todo imalloc_cpu

    integer(c_int) function magma_free_cpu( ptr ) &
    bind(C, name="magma_free_cpu")
        use iso_c_binding
        type(c_ptr), value :: ptr  !! void*
    end function

    !! -------------------------------------------------------------------------
    !! magma_malloc_pinned (pinned CPU main memory)
    integer(c_int) function magma_malloc_pinned( ptr, bytes ) &
    bind(C, name="magma_malloc_pinned")
        use iso_c_binding
        type(c_ptr), target :: ptr  !! void**
        integer(c_size_t), value :: bytes
    end function

    !! todo imalloc_pinned

    integer(c_int) function magma_free_pinned_internal( ptr, func, file, line ) &
    bind(C, name="magma_free_pinned_internal")
        use iso_c_binding
        type(c_ptr), value :: ptr  !! void*
        character(c_char), value :: func, file
        integer(c_int), value :: line
    end function

    !! -------------------------------------------------------------------------
    !! set/get
    subroutine magma_setmatrix_internal( &
        m, n, elemsize, hA_src, lda, dB_dst, ldb, queue, func, file, line ) &
    bind(C, name="magma_setmatrix_internal")
        use iso_c_binding
        integer(c_int),    value  :: m, n, elemsize, lda, ldb
        type(c_ptr),       value  :: hA_src
        type(c_ptr),       value  :: dB_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine

    subroutine magma_getmatrix_internal( &
        m, n, elemsize, dA_src, lda, hB_dst, ldb, queue, func, file, line ) &
    bind(C, name="magma_getmatrix_internal")
        use iso_c_binding
        integer(c_int),    value  :: m, n, elemsize, lda, ldb
        type(c_ptr),       value  :: dA_src
        type(c_ptr),       value  :: hB_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine
    
    subroutine magma_setmatrix_async_internal( &
        m, n, elemsize, hA_src, lda, dB_dst, ldb, queue, func, file, line ) &
    bind(C, name="magma_setmatrix_async_internal")
        use iso_c_binding
        integer(c_int),    value  :: m, n, elemsize, lda, ldb
        type(c_ptr),       value  :: hA_src
        type(c_ptr),       value  :: dB_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine

    subroutine magma_getmatrix_async_internal( &
        m, n, elemsize, dA_src, lda, hB_dst, ldb, queue, func, file, line ) &
    bind(C, name="magma_getmatrix_async_internal")
        use iso_c_binding
        integer(c_int),    value  :: m, n, elemsize, lda, ldb
        type(c_ptr),       value  :: dA_src
        type(c_ptr),       value  :: hB_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine
    
    subroutine magma_setvector_internal( &
        n, elemsize, hx_src, incx, dy_dst, incy, queue, func, file, line ) &
    bind(C, name="magma_setvector_internal")
        use iso_c_binding
        integer(c_int),    value  :: n, elemsize, incx, incy
        type(c_ptr),       value  :: hx_src
        type(c_ptr),       value  :: dy_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine

    subroutine magma_getvector_internal( &
        n, elemsize, dx_src, incx, hy_dst, incy, queue, func, file, line ) &
    bind(C, name="magma_getvector_internal")
        use iso_c_binding
        integer(c_int),    value  :: n, elemsize, incx, incy
        type(c_ptr),       value  :: dx_src
        type(c_ptr),       value  :: hy_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine
    
    subroutine magma_setvector_async_internal( &
        n, elemsize, hx_src, incx, dy_dst, incy, queue, func, file, line ) &
    bind(C, name="magma_setvector_async_internal")
        use iso_c_binding
        integer(c_int),    value  :: n, elemsize, incx, incy
        type(c_ptr),       value  :: hx_src
        type(c_ptr),       value  :: dy_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine

    subroutine magma_getvector_async_internal( &
        n, elemsize, dx_src, incx, hy_dst, incy, queue, func, file, line ) &
    bind(C, name="magma_getvector_async_internal")
        use iso_c_binding
        integer(c_int),    value  :: n, elemsize, incx, incy
        type(c_ptr),       value  :: dx_src
        type(c_ptr),       value  :: hy_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine

    !! -------------------------------------------------------------------------
    !! CPU interfaces (matrix in CPU memory)
    subroutine magma_dgetrf( m, n, A, lda, ipiv, info ) &
    bind(C, name="magma_dgetrf")
        use iso_c_binding
        integer(c_int),            value  :: m, n, lda
        real(c_double), target :: A(lda,*)
        integer(c_int),            target :: ipiv(*)
        integer(c_int),            target :: info  !! int*
    end subroutine

    subroutine magma_dpotrf( uplo, n, A, lda, info ) &
    bind(C, name="magma_dpotrf")
        use iso_c_binding
        integer(c_int),            value  :: uplo
        integer(c_int),            value  :: n, lda
        real(c_double), target :: A(lda,*)
        integer(c_int),            target :: info  !! int*
    end subroutine

    !! -------------------------------------------------------------------------
    !! GPU interfaces (matrix in GPU memory)
    subroutine magma_dgetrf_gpu( m, n, dA, lda, ipiv, info ) &
    bind(C, name="magma_dgetrf_gpu")
        use iso_c_binding
        integer(c_int), value  :: m, n, lda
        type(c_ptr),    value  :: dA
        integer(c_int), target :: ipiv(*)
        integer(c_int), target :: info  !! int*
    end subroutine

    subroutine magma_dpotrf_gpu( uplo, n, dA, lda, info ) &
    bind(C, name="magma_dpotrf_gpu")
        use iso_c_binding
        integer(c_int), value  :: uplo, n, lda
        type(c_ptr),    value  :: dA
        integer(c_int), target :: info  !! int*
    end subroutine

    subroutine magma_dgels_gpu( trans, m, n, nrhs, dA, lda, dB, ldb, hwork, lwork, info ) &
    bind(C, name="magma_dgels_gpu")
        use iso_c_binding
        integer(c_int), value  :: trans, m, n, nrhs, lda, ldb
        type(c_ptr),    value  :: dA, dB
        type(c_ptr),    value  :: hwork
        integer(c_int), value  :: lwork
        integer(c_int), target :: info  !! int*
    end subroutine

    subroutine magmablas_dtranspose( m, n, dA, ldda, dAT, lddat, queue ) &
    bind(C, name="magmablas_dtranspose")
        use iso_c_binding
        integer(c_int), value  :: m, n, ldda, lddat
        type(c_ptr),    value  :: dA, dAT
        type(c_ptr),    value  :: queue  !! queue_t
    end subroutine

    subroutine magmablas_dlacpy( uplo, m, n, dA, ldda, dB, lddb, queue ) &
    bind(C, name="magmablas_dlacpy")
        use iso_c_binding
        integer(c_int), value  :: uplo, m, n, ldda, lddb
        type(c_ptr),    value  :: dA, dB
        type(c_ptr),    value  :: queue  !! queue_t
    end subroutine

    subroutine magmablas_dlascl2( uplo, m, n, dD, dA, ldda, queue, info ) &
    bind(C, name="magmablas_dlascl2")
        use iso_c_binding
        integer(c_int), value  :: uplo, m, n, ldda
        type(c_ptr),    value  :: dD, dA
        type(c_ptr),    value  :: queue  !! queue_t
        integer(c_int), target :: info  !! int*
    end subroutine

    subroutine magmablas_dgeadd2( &
        m, n, &
        alpha, dA, ldda, &
        beta,  dB, lddb, &
        queue ) &
    bind(C, name="magmablas_dgeadd2")
        use iso_c_binding
        integer(c_int), value  :: m, n, ldda, lddb
        real(c_double), value  :: alpha, beta
        type(c_ptr),    value  :: dA, dB
        type(c_ptr),    value  :: queue  !! queue_t
    end subroutine

    !! -------------------------------------------------------------------------
    !! batched GPU interfaces (all arrays in GPU memory)
    subroutine magma_dgetrf_batched( &
        m, n, dA_array, lda, dipiv_array, dinfo_array, batchcount, queue ) &
    bind(C, name="magma_dgetrf_batched")
        use iso_c_binding
        integer(c_int), value  :: m, n, lda, batchcount
        type(c_ptr),    value  :: dA_array    !! double_real**
        type(c_ptr),    value  :: dipiv_array !! int**
        type(c_ptr),    value  :: dinfo_array !! int*
        type(c_ptr),    value  :: queue
    end subroutine

    subroutine magma_dgetrs_batched( &
        trans, n, nrhs, dA_array, ldda, dipiv_array, dB_array, lddb, batchcount, queue ) &
    bind(C, name="magma_dgetrs_batched")
        use iso_c_binding
        integer(c_int), value  :: trans, n, nrhs, ldda, lddb, batchcount
        type(c_ptr),    value  :: dA_array    !! double_real**
        type(c_ptr),    value  :: dipiv_array !! int**
        type(c_ptr),    value  :: dB_array    !! double_real**
        type(c_ptr),    value  :: queue
    end subroutine

    subroutine magmablas_dgemm_batched_strided( &
        transA, transB, m, n, k, &
        alpha, dA, lda, strideA, &
               dB, ldb, strideB, &
        beta,  dC, ldc, strideC, &
        batchcount, queue ) &
    bind(C, name="magmablas_dgemm_batched_strided")
        use iso_c_binding
        integer(c_int), value :: transA, transB, m, n, k, lda, ldb, ldc, batchcount
        integer(c_int), value :: strideA, strideB, strideC
        real(c_double), value :: alpha, beta
        type(c_ptr),    value :: dA, dB, dC
        type(c_ptr),    value :: queue  !! queue_t
    end subroutine

    !! -------------------------------------------------------------------------
    !! BLAS (matrices in GPU memory)
    real(c_double) function magma_dnrm2( &
        n, &
        dx, incx, &
        queue ) &
    bind(C, name="magma_dnrm2")
        use iso_c_binding
        integer(c_int),         value :: n, incx
        type(c_ptr),            value :: dx
        type(c_ptr),            value :: queue  !! queue_t
    end function

    subroutine magma_daxpy( &
        n, &
        alpha, dx, incx, &
               dy, incy, &
        queue ) &
    bind(C, name="magma_daxpy")
        use iso_c_binding
        integer(c_int),         value :: n, incx, incy
        real(c_double), value :: alpha
        type(c_ptr),            value :: dx, dy
        type(c_ptr),            value :: queue  !! queue_t
    end subroutine

    subroutine magma_dgemv( &
        transA, m, n, &
        alpha, dA, lda, &
               dx, incx, &
        beta,  dy, incy, &
        queue ) &
    bind(C, name="magma_dgemv")
        use iso_c_binding
        integer(c_int),         value :: transA, m, n, lda, incx, incy
        real(c_double), value :: alpha, beta
        type(c_ptr),            value :: dA, dx, dy
        type(c_ptr),            value :: queue  !! queue_t
    end subroutine

    subroutine magmablas_dgemv( &
        transA, m, n, &
        alpha, dA, lda, &
               dx, incx, &
        beta,  dy, incy, &
        queue ) &
    bind(C, name="magmablas_dgemv")
        use iso_c_binding
        integer(c_int),         value :: transA, m, n, lda, incx, incy
        real(c_double), value :: alpha, beta
        type(c_ptr),            value :: dA, dx, dy
        type(c_ptr),            value :: queue  !! queue_t
    end subroutine

    subroutine magma_dgemm( &
        transA, transB, m, n, k, &
        alpha, dA, lda, &
               dB, ldb, &
        beta,  dC, ldc, &
        queue ) &
    bind(C, name="magma_dgemm")
        use iso_c_binding
        integer(c_int),         value :: transA, transB, m, n, k, lda, ldb, ldc
        real(c_double), value :: alpha, beta
        type(c_ptr),            value :: dA, dB, dC
        type(c_ptr),            value :: queue  !! queue_t
    end subroutine

end interface

!! =============================================================================
!! Fortran routines & functions
contains

    !! -------------------------------------------------------------------------
    !! queue support
    subroutine magma_queue_create( dev, queue_ptr )
        use iso_c_binding
        integer(c_int), value :: dev
        type(c_ptr), target :: queue_ptr  !! queue_t*
        
        call magma_queue_create_internal( &
                dev, queue_ptr, &
                "magma_queue_create" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_queue_destroy( queue )
        use iso_c_binding
        type(c_ptr), value :: queue  !! queue_t
        
        call magma_queue_destroy_internal( &
                queue, &
                "magma_queue_destroy" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

#if defined(THORNADO_CUDA)
    subroutine magma_queue_create_from_cuda( dev, stream, cublas_handle, cusparse_handle, queue_ptr )
        use iso_c_binding
        integer(c_int), value :: dev
        type(c_ptr), value :: stream, cublas_handle, cusparse_handle
        type(c_ptr), target :: queue_ptr  !! queue_t*

        call magma_queue_create_from_cuda_internal( &
                dev, stream, cublas_handle, cusparse_handle, queue_ptr, &
                "magma_queue_create_from_cuda" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine
#elif defined(THORNADO_HIP)
    subroutine magma_queue_create_from_hip( dev, stream, hipblas_handle, hipsparse_handle, queue_ptr )
        use iso_c_binding
        integer(c_int), value :: dev
        type(c_ptr), value :: stream, hipblas_handle, hipsparse_handle
        type(c_ptr), target :: queue_ptr  !! queue_t*

        call magma_queue_create_from_hip_internal( &
                dev, stream, hipblas_handle, hipsparse_handle, queue_ptr, &
                "magma_queue_create_from_hip" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine
#endif

    subroutine magma_queue_sync( queue )
        use iso_c_binding
        type(c_ptr), value :: queue  !! queue_t
        
        call magma_queue_sync_internal( &
                queue, &
                "magma_queue_sync" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    !! -------------------------------------------------------------------------
    !! malloc wrappers
    integer(c_int) function magma_imalloc( ptr, n )
        use iso_c_binding
        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_imalloc = magma_malloc( ptr, n*sizeof_int )
    end function

    integer(c_int) function magma_imalloc_cpu( ptr, n )
        use iso_c_binding
        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_imalloc_cpu = magma_malloc_cpu( ptr, n*sizeof_int )
    end function

    integer(c_int) function magma_imalloc_pinned( ptr, n )
        use iso_c_binding
        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_imalloc_pinned = magma_malloc_pinned( ptr, n*sizeof_int )
    end function

    !! -------------------------------------------------------------------------
    !! magma_free wrappers
    integer(c_int) function magma_free( ptr )
        type(c_ptr) :: ptr
        
        magma_free = magma_free_internal( &
                ptr, &
                "magma_free" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end function

    integer(c_int) function magma_free_pinned( ptr )
        type(c_ptr) :: ptr
        
        magma_free_pinned = magma_free_internal( &
                ptr, &
                "magma_free_pinned" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end function

    !! -------------------------------------------------------------------------
    !! set/get wrappers
    subroutine magma_setmatrix( &
        m, n, elemsize, hA_src, lda, dB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int),    value  :: m, n, elemsize, lda, ldb
        type(c_ptr),       value  :: hA_src
        type(c_ptr),       value  :: dB_dst
        type(c_ptr),       value  :: queue
        
        call magma_setmatrix_internal( &
                m, n, elemsize, hA_src, lda, dB_dst, ldb, queue, &
                "magma_setmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_getmatrix( &
        m, n, elemsize, dA_src, lda, hB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int),    value  :: m, n, elemsize, lda, ldb
        type(c_ptr),       value  :: dA_src
        type(c_ptr),       value  :: hB_dst
        type(c_ptr),       value  :: queue
        
        call magma_getmatrix_internal( &
                m, n, elemsize, dA_src, lda, hB_dst, ldb, queue, &
                "magma_getmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine
    
    subroutine magma_setvector( &
        n, elemsize, hx_src, incx, dy_dst, incy, queue )
        use iso_c_binding
        integer(c_int),    value  :: n, elemsize, incx, incy
        type(c_ptr),       value  :: hx_src
        type(c_ptr),       value  :: dy_dst
        type(c_ptr),       value  :: queue
        
        call magma_setvector_internal( &
                n, elemsize, hx_src, incx, dy_dst, incy, queue, &
                "magma_setvector" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_getvector( &
        n, elemsize, dx_src, incx, hy_dst, incy, queue )
        use iso_c_binding
        integer(c_int),    value  :: n, elemsize, incx, incy
        type(c_ptr),       value  :: dx_src
        type(c_ptr),       value  :: hy_dst
        type(c_ptr),       value  :: queue
        
        call magma_getvector_internal( &
                n, elemsize, dx_src, incx, hy_dst, incy, queue, &
                "magma_getvector" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    !! -------------------------------------------------------------------------
    !! set/get wrappers
    !! matrices & vectors of integers
    subroutine magma_isetmatrix( &
        m, n, hA_src, lda, dB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int), value  :: m, n, lda, ldb
        integer(c_int), target :: hA_src(lda,*)
        type(c_ptr),    value  :: dB_dst
        type(c_ptr),    value  :: queue
        
        call magma_setmatrix_internal( &
                m, n, int(sizeof_int), c_loc(hA_src), lda, dB_dst, ldb, queue, &
                "magma_isetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_igetmatrix( &
        m, n, dA_src, lda, hB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int), value  :: m, n, lda, ldb
        type(c_ptr),    value  :: dA_src
        integer(c_int), target :: hB_dst(ldb,*)
        type(c_ptr),    value  :: queue
        
        call magma_getmatrix_internal( &
                m, n, int(sizeof_int), dA_src, lda, c_loc(hB_dst), ldb, queue, &
                "magma_igetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine
    
    subroutine magma_isetvector( &
        n, hx_src, incx, dy_dst, incy, queue )
        use iso_c_binding
        integer(c_int), value  :: n, incx, incy
        integer(c_int), target :: hx_src(*)
        type(c_ptr),    value  :: dy_dst
        type(c_ptr),    value  :: queue
        
        call magma_setvector_internal( &
                n, int(sizeof_int), c_loc(hx_src), incx, dy_dst, incy, queue, &
                "magma_isetvector" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_igetvector( &
        n, dx_src, incx, hy_dst, incy, queue )
        use iso_c_binding
        integer(c_int), value  :: n, incx, incy
        type(c_ptr),    value  :: dx_src
        integer(c_int), target :: hy_dst(*)
        type(c_ptr),    value  :: queue
        
        call magma_getvector_internal( &
                n, int(sizeof_int), dx_src, incx, c_loc(hy_dst), incy, queue, &
                "magma_igetvector" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    !! -------------------------------------------------------------------------
    !! set/get wrappers
    !! matrices & vectors of c_ptr pointers
    !subroutine magma_psetmatrix( &
    !    m, n, hA_src, lda, dB_dst, ldb, queue )
    !    use iso_c_binding
    !    integer(c_int), value  :: m, n, lda, ldb
    !    type(c_ptr),    target :: hA_src(lda,*)
    !    type(c_ptr),    value  :: dB_dst
    !    type(c_ptr),    value  :: queue
    !    
    !    call magma_setmatrix_internal( &
    !            m, n, int(sizeof_ptr), c_loc(hA_src), lda, dB_dst, ldb, queue, &
    !            "magma_psetmatrix" // c_null_char, &
    !            __FILE__ // c_null_char, &
    !            __LINE__ )
    !end subroutine

    !subroutine magma_pgetmatrix( &
    !    m, n, dA_src, lda, hB_dst, ldb, queue )
    !    use iso_c_binding
    !    integer(c_int), value  :: m, n, lda, ldb
    !    type(c_ptr),    value  :: dA_src
    !    type(c_ptr),    target :: hB_dst(ldb,*)
    !    type(c_ptr),    value  :: queue
    !    
    !    call magma_getmatrix_internal( &
    !            m, n, int(sizeof_ptr), dA_src, lda, c_loc(hB_dst), ldb, queue, &
    !            "magma_pgetmatrix" // c_null_char, &
    !            __FILE__ // c_null_char, &
    !            __LINE__ )
    !end subroutine
    
    !subroutine magma_psetvector( &
    !    n, hx_src, incx, dy_dst, incy, queue )
    !    use iso_c_binding
    !    integer(c_int), value  :: n, incx, incy
    !    type(c_ptr),    target :: hx_src(*)
    !    type(c_ptr),    value  :: dy_dst
    !    type(c_ptr),    value  :: queue
    !    
    !    call magma_setvector_internal( &
    !            n, int(sizeof_ptr), c_loc(hx_src), incx, dy_dst, incy, queue, &
    !            "magma_psetvector" // c_null_char, &
    !            __FILE__ // c_null_char, &
    !            __LINE__ )
    !end subroutine

    !subroutine magma_pgetvector( &
    !    n, dx_src, incx, hy_dst, incy, queue )
    !    use iso_c_binding
    !    integer(c_int), value  :: n, incx, incy
    !    type(c_ptr),    value  :: dx_src
    !    type(c_ptr),    target :: hy_dst(*)
    !    type(c_ptr),    value  :: queue
    !    
    !    call magma_getvector_internal( &
    !            n, int(sizeof_ptr), dx_src, incx, c_loc(hy_dst), incy, queue, &
    !            "magma_pgetvector" // c_null_char, &
    !            __FILE__ // c_null_char, &
    !            __LINE__ )
    !end subroutine

    !! -------------------------------------------------------------------------
    !! malloc wrappers
    integer(c_int) function magma_dmalloc( ptr, n )
        use iso_c_binding
        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_dmalloc = magma_malloc( ptr, n*sizeof_double )
    end function

    integer(c_int) function magma_dmalloc_cpu( ptr, n )
        use iso_c_binding
        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_dmalloc_cpu = magma_malloc_cpu( ptr, n*sizeof_double )
    end function

    integer(c_int) function magma_dmalloc_pinned( ptr, n )
        use iso_c_binding
        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_dmalloc_pinned = magma_malloc_pinned( ptr, n*sizeof_double )
    end function

    !! -------------------------------------------------------------------------
    !! set/get wrappers
    subroutine magma_dsetmatrix( &
        m, n, hA_src, lda, dB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int),            value  :: m, n, lda, ldb
        real(c_double), target :: hA_src(lda,*)
        type(c_ptr),               value  :: dB_dst
        type(c_ptr),               value  :: queue
        
        call magma_setmatrix_internal( &
                m, n, int(sizeof_double), c_loc(hA_src), lda, dB_dst, ldb, queue, &
                "magma_dsetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_dgetmatrix( &
        m, n, dA_src, lda, hB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int),            value  :: m, n, lda, ldb
        type(c_ptr),               value  :: dA_src
        real(c_double), target :: hB_dst(ldb,*)
        type(c_ptr),               value  :: queue
        
        call magma_getmatrix_internal( &
                m, n, int(sizeof_double), dA_src, lda, c_loc(hB_dst), ldb, queue, &
                "magma_dgetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_dsetmatrix_async( &
        m, n, hA_src, lda, dB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int),            value  :: m, n, lda, ldb
        real(c_double), target :: hA_src(lda,*)
        type(c_ptr),               value  :: dB_dst
        type(c_ptr),               value  :: queue
        
        call magma_setmatrix_async_internal( &
                m, n, int(sizeof_double), c_loc(hA_src), lda, dB_dst, ldb, queue, &
                "magma_dsetmatrix_async" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_dgetmatrix_async( &
        m, n, dA_src, lda, hB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int),            value  :: m, n, lda, ldb
        type(c_ptr),               value  :: dA_src
        real(c_double), target :: hB_dst(ldb,*)
        type(c_ptr),               value  :: queue
        
        call magma_getmatrix_async_internal( &
                m, n, int(sizeof_double), dA_src, lda, c_loc(hB_dst), ldb, queue, &
                "magma_dgetmatrix_async" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine
    
    subroutine magma_dsetvector( &
        n, hx_src, incx, dy_dst, incy, queue )
        use iso_c_binding
        integer(c_int), value  :: n, incx, incy
        real(c_double), target :: hx_src(*)
        type(c_ptr),    value  :: dy_dst
        type(c_ptr),    value  :: queue
        
        call magma_setvector_internal( &
                n, int(sizeof_double), c_loc(hx_src), incx, dy_dst, incy, queue, &
                "magma_dsetvector" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_dgetvector( &
        n, dx_src, incx, hy_dst, incy, queue )
        use iso_c_binding
        integer(c_int), value  :: n, incx, incy
        type(c_ptr),    value  :: dx_src
        real(c_double), target :: hy_dst(*)
        type(c_ptr),    value  :: queue
        
        call magma_getvector_internal( &
                n, int(sizeof_double), dx_src, incx, c_loc(hy_dst), incy, queue, &
                "magma_dgetvector" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine
    
    subroutine magma_dsetvector_async( &
        n, hx_src, incx, dy_dst, incy, queue )
        use iso_c_binding
        integer(c_int), value  :: n, incx, incy
        real(c_double), target :: hx_src(*)
        type(c_ptr),    value  :: dy_dst
        type(c_ptr),    value  :: queue
        
        call magma_setvector_async_internal( &
                n, int(sizeof_double), c_loc(hx_src), incx, dy_dst, incy, queue, &
                "magma_dsetvector_async" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_dgetvector_async( &
        n, dx_src, incx, hy_dst, incy, queue )
        use iso_c_binding
        integer(c_int), value  :: n, incx, incy
        type(c_ptr),    value  :: dx_src
        real(c_double), target :: hy_dst(*)
        type(c_ptr),    value  :: queue
        
        call magma_getvector_async_internal( &
                n, int(sizeof_double), dx_src, incx, c_loc(hy_dst), incy, queue, &
                "magma_dgetvector_async" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

end module MagmaModule
