module mod_csp
   use mod_csp_misc
   use mod_csp_force
   use mod_csp_dyn_ab
   use mod_csp_thdyn
   use mod_csp_vdiff_coef_gls
   use mod_mpi_csp_io
   use mod_mpi_csp_io_sbc
   use mod_mpi_csp_io_bdy
   !!YY
   use mitice
   use mitice_parameters
   implicit none

contains
!==============================================================================
   subroutine csp
!==============================================================================
      implicit none
      integer :: i, k
      integer :: clock_count1, clock_count2, clock_rate, clock_max
      real(dp) :: sshtemp, sshgloavg, rhsMax, rhsMin
      real(wp) :: maxtfld, mintfld

      call system_clock(clock_count1, clock_rate, clock_max)

      call csp_init_step

      !YY,ZY: consider the calling sequence of misc_rhoInSitu and seaice_main
      call csp_misc_rhoInSitu

      call mpi_csp_io_rdforce

      if (ln_bdy .and. mpi_check_bdy) then
         call mpi_csp_io_rdbdy
      end if

      call csp_force_core

      if (mitice_on) then
         call mitice_main        !main interface for seaice module
         call freeze_surface      !check SST, whether any point is below -1.9 degree
      end if

      if (ln_EmPmRevise) call csp_force_EmPmR_revise

      call csp_vdiff_coef_gls

      call csp_dyn_ab

      call csp_thdyn

      if (restart_out .and. (mod(myIter, nResFreq) .eq. 0)) then
         call mpi_csp_io_restart_send
         if (mpi_rank == 0) then
            write (*, "(a,i8)") " send restart file at step ", myIter
         end if
      end if

      ! send output to IO process
      CALL mpi_csp_io_send_output

      !YY seaice i/o, including frame dump and time-averaged frames dump
      !note: rank0 is responsible for seaice i/o. TO DO LATER FOR parallel I/O

      if(mitice_on) call seaice_output

      call system_clock(clock_count2, clock_rate, clock_max)

      IF (mpi_rank == 0) THEN
         write(run_out_unit, "(a,i8,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,f7.3,a)") &
            ' Integral for iter ', myIter, ', date ', nowyear, '-', nowmon, '-', nowday, ' ', &
            nowhour, ':', nowmin, ':', nowsec, &
            ', Elapsed time ', real(clock_count2 - clock_count1)/real(clock_rate), ' secs'
      END IF

   end subroutine csp
end module mod_csp
