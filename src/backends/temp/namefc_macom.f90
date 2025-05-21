program nmefc_macom
   use mod_misc
   use mod_csp
   use mod_csp_init
   USE mod_mpi_interfaces
   USE mod_mpi_variables
   USE mod_mpi_csp_io
   USE mod_mpi_test
   !YY
   USE mitice
   USE mitice_parameters
   USE mitice_utility
   USE mitice_init
   USE mitice_ave
   implicit none

   CALL mpi_process_init

   IF (mpi_rank < mpi_comp_procs) THEN

      call misc_namelist_read

!YY: ifdef of seaice module TO DO LATER
 !#ifdef SEAICE
      if (mitice_on) then
         !open an runtime file for screen output and model parameter summary
         call mitice_run_info_open
         !read seaice parameter namelist /namelist.mitice/
         call mitice_read_params
      endif
!#endif

      ! send information to io processes from computational processes
      IF (mpi_rank == 0) THEN
         CALL mpi_send_info_comp_to_io
      END IF

      call misc_run_info_open

      call csp_init

      if (mitice_on)then
         call mitice_init_allocate   !allocate host and device variables
         call mitice_init_fixed   !prepare a few parameters
         call gpu_seaice_parameters_update   !upload parameters to GPU device
         call mitice_init_vars      !initialized model variables and read initial fields
      endif

      if (restart_in) then
         call mpi_csp_io_restart_read
         if (assim_in) then
            call csp_init_asm
         end if
      end if
      
      IF (mpi_rank == 0) THEN
         ! send addtional infomation from comp communicator to io communicator
         CALL mpi_send_2nd_info_comp_to_io
      END IF

      ! ---- OPENACC: UPDATE VARS TO GPU FOR TIME INTEGRATION
      call GPU_VARS_INIT
      ! ----

      if (.not. restart_in) then
         ! send output to IO process
         CALL mpi_csp_io_send_output
      end if

      do myIter = nIter0, nIterMax
         IF (mpi_rank == 0) write(*,*) "myIter:", myIter
                          
         call csp

      end do

      call misc_run_info_close

   ELSE

      CALL mpi_csp_io_init
      CALL mpi_csp_io_main

   END IF

   if(mitice_on) then
      call mitice_run_info_close   !close seaice runing message file
      ! release memory for time-averaged vars
      IF ( SEAICE_taveFreq.GT.0.0_wp ) THEN
         call mitice_ave_release
      ENDIF
   endif

   CALL mpi_final_operations

end program nmefc_macom
