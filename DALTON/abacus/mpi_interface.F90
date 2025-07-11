! This module can not be used in gen1int, as it is compiled after gen1int.
! In that case, the include mpi_header and mpi_mod files are used.

module mpi_interface
#if defined(VAR_MPI)
#ifdef USE_MPI_MOD_F08 
  use mpi_f08
#elif defined(USE_MPI_MOD_F90)
  use mpi
#else
  #include "mpif.h"
#endif
#endif
end module mpi_interface
