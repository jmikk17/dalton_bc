module memory_parallel

#include "mpi_mod.h"
  implicit none
#include "mpi_header.h"

#ifdef VAR_MPI
  public memallocmpi
  public memfreempi

  private

  integer(kind=MPI_INTEGER_KIND) :: ierr

contains

  subroutine memallocmpi(nelement, ptr)
  integer(kind=mpi_address_kind), intent(inout) :: ptr
  integer(kind=mpi_address_kind), intent(in)    :: nelement
  call mpi_alloc_mem(nelement, mpi_info_null, ptr, ierr)
  end subroutine memallocmpi

  subroutine memfreempi(buf)
  real(8), intent(inout) :: buf(*)
  call mpi_free_mem(buf, ierr)
  end subroutine memfreempi
#endif

end module
