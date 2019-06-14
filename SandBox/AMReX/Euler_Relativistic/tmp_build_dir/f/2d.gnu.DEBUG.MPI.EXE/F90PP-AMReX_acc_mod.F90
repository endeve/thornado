
module amrex_acc_module

  implicit none

contains

  subroutine amrex_initialize_acc (id) bind(c,name='amrex_initialize_acc')



    integer, intent(in), value :: id




  end subroutine amrex_initialize_acc

  subroutine amrex_finalize_acc () bind(c,name='amrex_finalize_acc')




  end subroutine amrex_finalize_acc

end module amrex_acc_module

