!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2013 all rights reserved
!================================================================
      module omp_lib
c     Dummy module for stupid PGF90
      external omp_get_thread_num, omp_get_num_threads
      integer  omp_get_thread_num, omp_get_num_threads
      external omp_in_parallel
      logical  omp_in_parallel
      end module
