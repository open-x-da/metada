module mod_mpi_csp_io
   use mod_misc
   use mod_io_netcdf
   use mod_mpi_variables
   use mod_mpi_interfaces
   use mod_mpi_csp_io_netcdf4
   use mod_mpi_csp_io_restart
   use mod_mpi_test
   implicit none

contains

!==============================================================================
   subroutine mpi_csp_io_send_output
!==============================================================================
      implicit none
      integer :: i, k, djulian
      integer(i8) :: julian_st
      logical :: flag_write_out, out_avg

      djulian = 17

      flag_write_out = .false.

      if (trim(fOutFreq) .eq. "daily") then
         call greg2jul(0, 0, 0, nowday, nowmon, nowyear, julian_st)
         djulian = int(nowjulian - julian_st, 4)
      else if (trim(fOutFreq) .eq. "monthly") then
         call greg2jul(0, 0, 0, 1, nowmon, nowyear, julian_st)
         djulian = int(nowjulian - julian_st, 4)
      else if (trim(fOutFreq) .eq. "yearly") then
         call greg2jul(0, 0, 0, 1, 1, nowyear, julian_st)
         djulian = int(nowjulian - julian_st, 4)
      else
         write (run_out_unit, *) ' ERROR ! : Model only support daily, monthly or yearly now !'
         stop
      end if

      if (mod(djulian, nOutFreq) .eq. 0) flag_write_out = .true.

      out_avg = .true.

      if (out_avg) then
         !$ACC kernels present(tFld_out,uFld_out,vFld_out,wFld_out,KappaRT_out,KappaRM_out, &
         !$ACC                tFld,uFld,vFld,wFld,KappaRT,KappaRM,  &
         !$ACC       ssh_out,etaH_out,utau_out,vtau_out,qns_out,qsr_out,EmPmR_out, &
         !$ACC       ssh,etaH,utau,vtau,qns,qsr,EmPmR)

         !$ACC loop collapse(2) independent
         do k = 1, nk
            do i = 1, loc_nlpb
               tFld_out(i, k, 1) = tFld_out(i, k, 1) + tFld(i, k, 1)
               tFld_out(i, k, 2) = tFld_out(i, k, 2) + tFld(i, k, 2)
               uFld_out(i, k) = uFld_out(i, k) + uFld(i, k)
               vFld_out(i, k) = vFld_out(i, k) + vFld(i, k)
               wFld_out(i, k + 1) = wFld_out(i, k + 1) + wFld(i, k + 1)
               !KappaRT_out(i, k + 1) = KappaRT_out(i, k + 1) + KappaRT(i, k + 1)
               !KappaRM_out(i, k + 1) = KappaRM_out(i, k + 1) + KappaRM(i, k + 1)
            end do
         end do

         !$ACC loop
         do i = 1, loc_nlpb
            ssh_out(i) = ssh_out(i) + ssh(i)
            etaH_out(i) = etaH_out(i) + etaH(i)
            utau_out(i) = utau_out(i) + utau(i)
            vtau_out(i) = vtau_out(i) + vtau(i)
            qns_out(i) = qns_out(i) + qns(i)
            qsr_out(i) = qsr_out(i) + qsr(i)
            EmPmR_out(i) = EmPmR_out(i) + EmPmR(i)
         end do
         !$ACC end kernels

         output_nmuber = output_nmuber + 1.0_wp
      end if

      if (flag_write_out) then

         if (out_avg) then
            !$ACC update device(output_nmuber)

# if defined (EXDIAG)
            !$ACC kernels present(tFld_out,uFld_out,vFld_out,wFld_out,KappaRT_out,KappaRM_out, &
            !$ACC       ssh_out,etaH_out,utau_out,vtau_out,qns_out,qsr_out,EmPmR_out,pbt_base, &
            !$ACC       tFld_adv, sFld_adv, tFld_vdiff, sFld_vdiff, tFld_hdiff, sFld_hdiff, tFld_force, sFld_force)
# else 
            !$ACC kernels present(tFld_out,uFld_out,vFld_out,wFld_out,KappaRT_out,KappaRM_out, &
            !$ACC       ssh_out,etaH_out,utau_out,vtau_out,qns_out,qsr_out,EmPmR_out,pbt_base)
# endif
            !$ACC loop collapse(2) independent
            do k = 1, nk
               do i = 1, loc_nlpb
                  tFld_out(i, k, 1) = tFld_out(i, k, 1)/ output_nmuber
                  tFld_out(i, k, 2) = tFld_out(i, k, 2)/ output_nmuber
                  uFld_out(i, k) = uFld_out(i, k)/ output_nmuber
                  vFld_out(i, k) = vFld_out(i, k)/ output_nmuber
                  wFld_out(i, k + 1) = wFld_out(i, k + 1)/ output_nmuber
                  !KappaRT_out(i, k + 1) = KappaRT_out(i, k + 1)/ output_nmuber
                  !KappaRM_out(i, k + 1) = KappaRM_out(i, k + 1)/ output_nmuber
               end do
            end do

            !$ACC loop
            do i = 1, loc_nlpb
               ssh_out(i) = ssh_out(i)/ output_nmuber
               etaH_out(i) = etaH_out(i)/ output_nmuber + pbt_base(i)
               utau_out(i) = utau_out(i)/ output_nmuber
               vtau_out(i) = vtau_out(i)/ output_nmuber
               qns_out(i) = qns_out(i)/ output_nmuber
               qsr_out(i) = qsr_out(i)/ output_nmuber
               EmPmR_out(i) = EmPmR_out(i)/ output_nmuber
            end do

# if defined (EXDIAG)
            !$ACC loop collapse(2) independent
            do k = 1, nk
               do i = 1, loc_nlpb
                  tFld_adv(i, k) = tFld_adv(i, k)/output_nmuber
                  sFld_adv(i, k) = sFld_adv(i, k)/output_nmuber
                  tFld_vdiff(i, k) = tFld_vdiff(i, k)/output_nmuber
                  sFld_vdiff(i, k) = sFld_vdiff(i, k)/output_nmuber
                  tFld_hdiff(i, k) = tFld_hdiff(i, k)/output_nmuber
                  sFld_hdiff(i, k) = sFld_hdiff(i, k)/output_nmuber
                  tFld_force(i, k) = tFld_force(i, k)/output_nmuber
                  sFld_force(i, k) = sFld_force(i, k)/output_nmuber
               end do
            end do
# endif

            !$ACC end kernels
         else
            !$ACC kernels present(tFld_out,uFld_out,vFld_out,wFld_out,KappaRT_out,KappaRM_out, &
            !$ACC       tFld,uFld,vFld,wFld,KappaRT,KappaRM,pbt_base,  &
            !$ACC       ssh_out,etaH_out,utau_out,vtau_out,qns_out,qsr_out,EmPmR_out, &
            !$ACC       ssh,etaH,utau,vtau,qns,qsr,EmPmR)

            !$ACC loop collapse(2) independent
            do k = 1, nk
               do i = 1, loc_nlpb
                  tFld_out(i, k, 1) = tFld(i, k, 1)
                  tFld_out(i, k, 2) = tFld(i, k, 2)
                  uFld_out(i, k) = uFld(i, k)
                  vFld_out(i, k) = vFld(i, k)
                  wFld_out(i, k + 1) = wFld(i, k + 1)
                  !KappaRT_out(i, k + 1) = KappaRT(i, k + 1)
                  !KappaRM_out(i, k + 1) = KappaRM(i, k + 1)
                  !if (k .eq. nk) write(*,*) tFld(i, k, 1)
               end do
            end do

            !$ACC loop
            do i = 1, loc_nlpb
               ssh_out(i) = ssh(i)
               etaH_out(i) = etaH(i) + pbt_base(i)
               utau_out(i) = utau(i)
               vtau_out(i) = vtau(i)
               qns_out(i) = qns(i)
               qsr_out(i) = qsr(i)
               EmPmR_out(i) = EmPmR(i)
            end do

            !$ACC end kernels
         end if

         ! IF (mpi_rank == 0)  write (run_out_unit, "(a,i8)") &
         !    " send output file at step ", myIter

         ! ---- OPENACC: OFFLOAD DEVICE VARS TO HOST FOR I/O
         !wFld_out,KappaRT_out,KappaRM_out,
# if defined (EXDIAG)
         !$ACC update self(tFld_out,uFld_out,vFld_out, &
         !$ACC       ssh_out,etaH_out,utau_out,vtau_out,qns_out,qsr_out,EmPmR_out, &
         !$ACC       tFld_adv, sFld_adv, tFld_vdiff, sFld_vdiff, tFld_hdiff, sFld_hdiff, tFld_force, sFld_force)
# else
         !$ACC update self(tFld_out,uFld_out,vFld_out, &
         !$ACC       ssh_out,etaH_out,utau_out,vtau_out,qns_out,qsr_out,EmPmR_out)
# endif
! ----
         CALL mpi_csp_io_send_output_1d(ssh_out, mpi_comp_buf_1d_ssh, mpi_req_ssh, &
                                        mpi_flag_ssh)
         CALL mpi_csp_io_send_output_2d(tFld_out(:, :, 1), mpi_comp_buf_2d_tFld, &
                                        mpi_req_tFld, mpi_flag_tFld)
         CALL mpi_csp_io_send_output_2d(tFld_out(:, :, 2), mpi_comp_buf_2d_sFld, &
                                        mpi_req_sFld, mpi_flag_sFld)
         CALL mpi_csp_io_send_output_2d(uFld_out, mpi_comp_buf_2d_uFld, &
                                        mpi_req_uFld, mpi_flag_uFld)
         CALL mpi_csp_io_send_output_2d(vFld_out, mpi_comp_buf_2d_vFld, &
                                        mpi_req_vFld, mpi_flag_vFld)
         CALL mpi_csp_io_send_output_2d(wFld_out(1:loc_nlpb, 2:nkp1), mpi_comp_buf_2d_wFld, &
                                        mpi_req_wFld, mpi_flag_wFld)
         !CALL mpi_csp_io_send_output_2d(KappaRT_out(1:loc_nlpb, 2:nkp1), mpi_comp_buf_2d_KappaRT, &
         !                               mpi_req_KappaRT, mpi_flag_KappaRT)
         !CALL mpi_csp_io_send_output_2d(KappaRM_out(1:loc_nlpb, 2:nkp1), mpi_comp_buf_2d_KappaRM, &
         !                               mpi_req_KappaRM, mpi_flag_KappaRM)
         CALL mpi_csp_io_send_output_1d(etaH_out, mpi_comp_buf_1d_etaH, mpi_req_etaH, &
                                        mpi_flag_etaH)
         CALL mpi_csp_io_send_output_1d(utau_out, mpi_comp_buf_1d_utau, mpi_req_utau, &
                                        mpi_flag_utau)
         CALL mpi_csp_io_send_output_1d(vtau_out, mpi_comp_buf_1d_vtau, mpi_req_vtau, &
                                        mpi_flag_vtau)
         CALL mpi_csp_io_send_output_1d(qns_out, mpi_comp_buf_1d_qns, mpi_req_qns, &
                                        mpi_flag_qns)
         CALL mpi_csp_io_send_output_1d(qsr_out, mpi_comp_buf_1d_qsr, mpi_req_qsr, &
                                        mpi_flag_qsr)
         CALL mpi_csp_io_send_output_1d(EmPmR_out, mpi_comp_buf_1d_EmPmR, mpi_req_EmPmR, &
                                        mpi_flag_EmPmR)

# if defined (EXDIAG)
         CALL mpi_csp_io_send_output_2d(tFld_adv, mpi_comp_buf_2d_tFld_adv, &
                                        mpi_req_tFld_adv, mpi_flag_tFld_adv)
         CALL mpi_csp_io_send_output_2d(sFld_adv, mpi_comp_buf_2d_sFld_adv, &
                                        mpi_req_sFld_adv, mpi_flag_sFld_adv)
         CALL mpi_csp_io_send_output_2d(tFld_vdiff, mpi_comp_buf_2d_tFld_vdiff, &
                                        mpi_req_tFld_vdiff, mpi_flag_tFld_vdiff)
         CALL mpi_csp_io_send_output_2d(sFld_vdiff, mpi_comp_buf_2d_sFld_vdiff, &
                                        mpi_req_sFld_vdiff, mpi_flag_sFld_vdiff)
         CALL mpi_csp_io_send_output_2d(tFld_hdiff, mpi_comp_buf_2d_tFld_hdiff, &
                                        mpi_req_tFld_hdiff, mpi_flag_tFld_hdiff)
         CALL mpi_csp_io_send_output_2d(sFld_hdiff, mpi_comp_buf_2d_sFld_hdiff, &
                                        mpi_req_sFld_hdiff, mpi_flag_sFld_hdiff)
         CALL mpi_csp_io_send_output_2d(tFld_force, mpi_comp_buf_2d_tFld_force, &
                                        mpi_req_tFld_force, mpi_flag_tFld_force)
         CALL mpi_csp_io_send_output_2d(sFld_force, mpi_comp_buf_2d_sFld_force, &
                                        mpi_req_sFld_force, mpi_flag_sFld_force)
# endif
         
         if (out_avg) then
            !wFld_out,KappaRT_out,KappaRM_out,
# if defined (EXDIAG)
            !$ACC kernels present(tFld_out,uFld_out,vFld_out, &
            !$ACC       ssh_out,etaH_out,utau_out,vtau_out,qns_out,qsr_out,EmPmR_out, &
            !$ACC       tFld_adv, sFld_adv, tFld_vdiff, sFld_vdiff, tFld_hdiff, sFld_hdiff, tFld_force, sFld_force)
# else
            !$ACC kernels present(tFld_out,uFld_out,vFld_out, &
            !$ACC       ssh_out,etaH_out,utau_out,vtau_out,qns_out,qsr_out,EmPmR_out)
# endif
            !$ACC loop collapse(2) independent
            do k = 1, nk
               do i = 1, loc_nlpb
                  tFld_out(i, k, 1) = 0.0_wp
                  tFld_out(i, k, 2) = 0.0_wp
                  uFld_out(i, k) = 0.0_wp
                  vFld_out(i, k) = 0.0_wp
                  wFld_out(i, k + 1) = 0.0_wp
                  !KappaRT_out(i, k + 1) = 0.0_wp
                  !KappaRM_out(i, k + 1) = 0.0_wp
               end do
            end do

            !$ACC loop
            do i = 1, loc_nlpb
               ssh_out(i) = 0.0_wp
               etaH_out(i) = 0.0_wp
               utau_out(i) = 0.0_wp
               vtau_out(i) = 0.0_wp
               qns_out(i) = 0.0_wp
               qsr_out(i) = 0.0_wp
               EmPmR_out(i) = 0.0_wp
            end do
            !$ACC end kernels

# if defined (EXDIAG)
            !$ACC loop collapse(2) independent
            do k = 1, nk
               do i = 1, loc_nlpb
                  tFld_adv(i, k) = 0.0_wp
                  sFld_adv(i, k) = 0.0_wp
                  tFld_vdiff(i, k) = 0.0_wp
                  sFld_vdiff(i, k) = 0.0_wp
                  tFld_hdiff(i, k) = 0.0_wp
                  sFld_hdiff(i, k) = 0.0_wp
                  tFld_force(i, k) = 0.0_wp
                  sFld_force(i, k) = 0.0_wp
               end do
            end do
# endif

            output_nmuber = 0.0_wp
         end if
      end if

   end subroutine mpi_csp_io_send_output

!==============================================================================
   SUBROUTINE mpi_csp_io_recv_output
!==============================================================================
      implicit none

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_ssh, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_ssh)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_tFld(:, 1), &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_tFld)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_tFld(:, 2), &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_sFld)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_uFld, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_uFld)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_vFld, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_vFld)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_wFld, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_wFld)

      !CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_KappaRT, &
      !                                  mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
      !                                  mpi_comp_io_comm, mpi_req_KappaRT)

      !CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_KappaRM, &
      !                                  mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
      !                                  mpi_comp_io_comm, mpi_req_KappaRM)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_etaH, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_etaH)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_utau, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_utau)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_vtau, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_vtau)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_qns, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_qns)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_qsr, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_qsr)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_1d_EmPmR, mpi_nlpb_counts_all, &
                                        mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm, mpi_req_EmPmR)

# if defined (EXDIAG)
      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_tFld_adv, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_tFld_adv) 
                                        
      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_sFld_adv, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_sFld_adv)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_tFld_vdiff, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_tFld_vdiff) 
                                        
      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_sFld_vdiff, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_sFld_vdiff)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_tFld_hdiff, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_tFld_hdiff) 
                                        
      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_sFld_hdiff, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_sFld_hdiff)

      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_tFld_force, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_tFld_force) 
                                        
      CALL mpi_igatherv_root_comp_to_io(mpi_io_buf_2d_sFld_force, &
                                        mpi_2d_nlpb_nk_counts_all, mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, &
                                        mpi_comp_io_comm, mpi_req_sFld_force)
# endif
   END SUBROUTINE mpi_csp_io_recv_output

!==============================================================================
   SUBROUTINE mpi_csp_io_root_scatterv_output(filename, ts)
!==============================================================================
      character(lc), INTENT(IN) :: filename
      integer, INTENT(IN) :: ts

      CALL MPI_WAIT(mpi_req_ssh, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_ssh, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_1d, mpi_io_buf_1d_send_counts, &
                        mpi_io_buf_1d_send_displs, mpi_real_wp, mpi_io_buf_seperate_1d, &
                        mpi_io_buf_seperate_1d_count, mpi_real_wp, mpi_root_io, &
                        mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'ssh', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_WAIT(mpi_req_tFld, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_tFld(:, 1), mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 't', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_sFld, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_tFld(:, 2), mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 's', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_uFld, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_uFld, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'u', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_vFld, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_vFld, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'v', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_wFld, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_wFld, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'w', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)
      
      !CALL MPI_WAIT(mpi_req_KappaRT, mpi_stat, mpi_err)
      !CALL mpi_2d_data_adjust(mpi_io_buf_2d_KappaRT, mpi_io_buf_adjust_2d, &
      !                        mpi_comp_procs, mpi_nlpb_block_1d_all, &
      !                        mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
      !                        mpi_nlpb_block_1d_displs_all)
      !CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
      !                  mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
      !                  mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      !CALL mpi_netcdf_write(filename, 'kappa_h', mpi_io_buf_seperate_2d, 'sp', &
      !                      mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
      !                      mpi_io_buf_seperate_2d_count_nk, ts)
      !
      !CALL MPI_WAIT(mpi_req_KappaRM, mpi_stat, mpi_err)
      !CALL mpi_2d_data_adjust(mpi_io_buf_2d_KappaRM, mpi_io_buf_adjust_2d, &
      !                        mpi_comp_procs, mpi_nlpb_block_1d_all, &
      !                        mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
      !                        mpi_nlpb_block_1d_displs_all)
      !CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
      !                  mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
      !                  mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      !CALL mpi_netcdf_write(filename, 'kappa_m', mpi_io_buf_seperate_2d, 'sp', &
      !                      mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
      !                      mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_etaH, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_etaH, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_1d, mpi_io_buf_1d_send_counts, &
                        mpi_io_buf_1d_send_displs, mpi_real_wp, mpi_io_buf_seperate_1d, &
                        mpi_io_buf_seperate_1d_count, mpi_real_wp, mpi_root_io, &
                        mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'pbt', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_WAIT(mpi_req_utau, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_utau, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_1d, mpi_io_buf_1d_send_counts, &
                        mpi_io_buf_1d_send_displs, mpi_real_wp, mpi_io_buf_seperate_1d, &
                        mpi_io_buf_seperate_1d_count, mpi_real_wp, mpi_root_io, &
                        mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'utau', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_WAIT(mpi_req_vtau, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_vtau, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_1d, mpi_io_buf_1d_send_counts, &
                        mpi_io_buf_1d_send_displs, mpi_real_wp, mpi_io_buf_seperate_1d, &
                        mpi_io_buf_seperate_1d_count, mpi_real_wp, mpi_root_io, &
                        mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'vtau', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_WAIT(mpi_req_qns, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_qns, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_1d, mpi_io_buf_1d_send_counts, &
                        mpi_io_buf_1d_send_displs, mpi_real_wp, mpi_io_buf_seperate_1d, &
                        mpi_io_buf_seperate_1d_count, mpi_real_wp, mpi_root_io, &
                        mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'qns', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_WAIT(mpi_req_qsr, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_qsr, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_1d, mpi_io_buf_1d_send_counts, &
                        mpi_io_buf_1d_send_displs, mpi_real_wp, mpi_io_buf_seperate_1d, &
                        mpi_io_buf_seperate_1d_count, mpi_real_wp, mpi_root_io, &
                        mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'qsr', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_WAIT(mpi_req_EmPmR, mpi_stat, mpi_err)
      CALL mpi_1d_data_adjust(mpi_io_buf_1d_EmPmR, mpi_io_buf_adjust_1d, mpi_comp_procs, &
                              mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_1d, mpi_io_buf_1d_send_counts, &
                        mpi_io_buf_1d_send_displs, mpi_real_wp, mpi_io_buf_seperate_1d, &
                        mpi_io_buf_seperate_1d_count, mpi_real_wp, mpi_root_io, &
                        mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'empmr', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

# if defined (EXDIAG)                            
      CALL MPI_WAIT(mpi_req_tFld_adv, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_tFld_adv, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 't_adv', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_sFld_adv, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_sFld_adv, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 's_adv', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_tFld_vdiff, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_tFld_vdiff, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 't_vdiff', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_sFld_vdiff, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_sFld_vdiff, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 's_vdiff', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_tFld_hdiff, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_tFld_hdiff, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 't_hdiff', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_sFld_hdiff, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_sFld_hdiff, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 's_hdiff', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_tFld_force, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_tFld_force, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, &
                        mpi_io_buf_2d_displs, mpi_real_wp, mpi_io_buf_seperate_2d, &
                        mpi_total_nlpb*mpi_io_nk, mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 't_force', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_WAIT(mpi_req_sFld_force, mpi_stat, mpi_err)
      CALL mpi_2d_data_adjust(mpi_io_buf_2d_sFld_force, mpi_io_buf_adjust_2d, &
                              mpi_comp_procs, mpi_nlpb_block_1d_all, &
                              mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, &
                              mpi_nlpb_block_1d_displs_all)
      CALL MPI_SCATTERV(mpi_io_buf_adjust_2d, mpi_io_buf_2d_send_counts, mpi_io_buf_2d_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 's_force', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)
# endif
   END SUBROUTINE mpi_csp_io_root_scatterv_output

!==============================================================================
   SUBROUTINE mpi_csp_io_scatterv_output(filename, ts)
!==============================================================================
      character(lc), INTENT(IN) :: filename
      integer, INTENT(IN) :: ts
      REAL(wp) :: mpi_fake_send_buf(1)
      INTEGER :: mpi_fake_send_count(1)
      INTEGER :: mpi_fake_displs(1)

      mpi_fake_send_buf(1) = 0.0_wp
      mpi_fake_send_count(1) = 0
      mpi_fake_displs(1) = 0

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_1d, mpi_io_buf_seperate_1d_count, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'ssh', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 't', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 's', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'u', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'v', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'w', mpi_io_buf_seperate_2d, 'sp', &
                            mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                            mpi_io_buf_seperate_2d_count_nk, ts)
      
      !CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
      !                  mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
      !                  mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      !CALL mpi_netcdf_write(filename, 'kappa_h', mpi_io_buf_seperate_2d, 'sp', &
      !                      mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
      !                      mpi_io_buf_seperate_2d_count_nk, ts)
      !
      !CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
      !                  mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
      !                  mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      !CALL mpi_netcdf_write(filename, 'kappa_m', mpi_io_buf_seperate_2d, 'sp', &
      !                      mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
      !                      mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_1d, mpi_io_buf_seperate_1d_count, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'pbt', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_1d, mpi_io_buf_seperate_1d_count, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'utau', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_1d, mpi_io_buf_seperate_1d_count, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'vtau', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_1d, mpi_io_buf_seperate_1d_count, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'qns', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_1d, mpi_io_buf_seperate_1d_count, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'qsr', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_1d, mpi_io_buf_seperate_1d_count, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 'empmr', mpi_io_buf_seperate_1d, 'sp', &
                            mpi_total_nlpb, mpi_io_buf_seperate_1d_count, mpi_io_buf_seperate_1d_start, &
                            mpi_io_buf_seperate_1d_count, ts)

# if defined (EXDIAG)  
      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 't_adv', mpi_io_buf_seperate_2d, 'sp', &
                           mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                           mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 's_adv', mpi_io_buf_seperate_2d, 'sp', &
                           mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                           mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 't_vdiff', mpi_io_buf_seperate_2d, 'sp', &
                           mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                           mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 's_vdiff', mpi_io_buf_seperate_2d, 'sp', &
                           mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                           mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 't_hdiff', mpi_io_buf_seperate_2d, 'sp', &
                           mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                           mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 's_hdiff', mpi_io_buf_seperate_2d, 'sp', &
                           mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                           mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 't_force', mpi_io_buf_seperate_2d, 'sp', &
                           mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                           mpi_io_buf_seperate_2d_count_nk, ts)

      CALL MPI_SCATTERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_fake_displs, &
                        mpi_real_wp, mpi_io_buf_seperate_2d, mpi_total_nlpb*mpi_io_nk, &
                        mpi_real_wp, mpi_root_io, mpi_io_comm, mpi_err)
      CALL mpi_netcdf_write(filename, 's_force', mpi_io_buf_seperate_2d, 'sp', &
                           mpi_total_nlpb, nk, mpi_total_nlpb, mpi_io_nk, mpi_io_buf_seperate_2d_start_nk, &
                           mpi_io_buf_seperate_2d_count_nk, ts)
# endif

   END SUBROUTINE mpi_csp_io_scatterv_output

!==============================================================================
   SUBROUTINE mpi_csp_io_init
!==============================================================================
      implicit none

   END SUBROUTINE mpi_csp_io_init

!==============================================================================
   SUBROUTINE mpi_csp_io_main
!==============================================================================
      implicit none
      character(lc) :: filename
      integer :: ts, djulian, djulianf
      integer(i8) :: julian_st, julianf, julian_stf
      logical :: flag_write_out
      type(nc_attr) :: nc_attr1

      !---vvv receve basic information form compute root process
      IF (mpi_rank == mpi_comp_procs) THEN
         CALL mpi_recv_info_io_from_comp
         CALL mpi_send_info_io
      ELSE
         CALL mpi_recv_info_io
      END IF
      !---^^^ end receve basic information

      nc_attr1%var_units = 'seconds since 1949-10-01 00:00:00'

      !---vvv write initial field, this will only execute once at model run first
      call greg2jul(nowsec, nowmin, nowhour, nowday, nowmon, nowyear, nowjulian)
      ! set for the filename and real record time
      julianf = nowjulian - nOutFreq/2
      call jul2greg(nowsecf, nowminf, nowhourf, nowdayf, nowmonf, nowyearf, julianf)

      djulian = 17

      flag_write_out = .false.

      if (trim(fOutFreq) .eq. "daily") then
         write (filename, "(a,i4.4,i2.2,i2.2,a,i2.2,a)") &
            trim(fout_dir)//"MaCOM_"//trim(fOutFreq)//"_", nowyearf, nowmonf, &
            nowdayf, "_initial_d", idom, ".nc"
         call greg2jul(0, 0, 0, nowday, nowmon, nowyear, julian_st)
         call greg2jul(0, 0, 0, nowdayf, nowmonf, nowyearf, julian_stf)
         djulian = int(nowjulian - julian_st, 4)
         djulianf = int(julianf - julian_stf, 4)
      else if (trim(fOutFreq) .eq. "monthly") then
         write (filename, "(a,i4.4,i2.2,a,i2.2,a)") &
            trim(fout_dir)//"MaCOM_"//trim(fOutFreq)//"_", nowyearf, nowmonf, &
            "_initial_d", idom, ".nc"
         call greg2jul(0, 0, 0, 1, nowmon, nowyear, julian_st)
         call greg2jul(0, 0, 0, 1, nowmonf, nowyearf, julian_stf)
         djulian = int(nowjulian - julian_st, 4)
         djulianf = int(julianf - julian_stf, 4)
      else if (trim(fOutFreq) .eq. "yearly") then
         write (filename, "(a,i4.4,a,i2.2,a)") trim(fout_dir)//"MaCOM_"// &
            trim(fOutFreq)//"_", nowyearf, "_initial_d", idom, ".nc"
         call greg2jul(0, 0, 0, 1, 1, nowyear, julian_st)
         call greg2jul(0, 0, 0, 1, 1, nowyearf, julian_stf)
         djulian = int(nowjulian - julian_st, 4)
         djulianf = int(julianf - julian_stf, 4)
      else
         write (*, *) ' ERROR! : Model only support daily, monthly or yearly now!, my rank is ', mpi_rank
         stop
      end if

      if ((mod(djulian, nOutFreq) .eq. 0) .and. (.not. restart_in)) then
         flag_write_out = .true.
      end if

      if (flag_write_out) ts = 1

      IF (mpi_rank == mpi_comp_procs) THEN
         CALL mpi_recv_2nd_info_io_from_comp
         CALL mpi_send_2nd_info_io
         CALL mpi_csp_io_root_init

         if (flag_write_out) then
            CALL mpi_csp_io_recv_output

            call netcdf_write(trim(filename), "time", dble(julianf), &
                              'dp', ts=ts, nc_attr1=nc_attr1)

            CALL mpi_csp_io_root_scatterv_output(filename, ts)
         end if
      ELSE
         CALL mpi_recv_2nd_info_io
         CALL mpi_csp_io_non_root_init

         if (flag_write_out) then
            CALL mpi_csp_io_scatterv_output(filename, ts)
         end if
      END IF
      !---^^^ write initial field finished

      call greg2jul(nowsec, nowmin, nowhour, nowday, nowmon, nowyear, nowjulian)

      !---vvv write regular output field
      do myIterm = nIter0, nIterMax

         nowjulian = nowjulian + int(dTtracer, 8)

         call jul2greg(nowsec, nowmin, nowhour, nowday, nowmon, nowyear, nowjulian)
         julianf = nowjulian - nOutFreq/2
         call jul2greg(nowsecf, nowminf, nowhourf, nowdayf, nowmonf, nowyearf, julianf)

         if (mpi_rank == mpi_comp_procs) then
            if (restart_out .and. (mod(myIterm, nResFreq) .eq. 0)) then
               ! write (*, "(a,i8)") " recieve restart file and write at step ", myIterm
               call mpi_csp_io_restart_main
            end if
         end if

         djulian = 17

         flag_write_out = .false.

         if (trim(fOutFreq) .eq. "daily") then
            write (filename, "(a,i4.4,i2.2,i2.2,a,i2.2,a)") &
               trim(fout_dir)//"MaCOM_"//trim(fOutFreq)//"_", nowyearf, nowmonf, &
               nowdayf, "_d", idom, ".nc"
            call greg2jul(0, 0, 0, nowday, nowmon, nowyear, julian_st)
            call greg2jul(0, 0, 0, nowdayf, nowmonf, nowyearf, julian_stf)
            djulian = int(nowjulian - julian_st, 4)
            djulianf = int(julianf - julian_stf, 4)
         else if (trim(fOutFreq) .eq. "monthly") then
            write (filename, "(a,i4.4,i2.2,a,i2.2,a)") &
               trim(fout_dir)//"MaCOM_"//trim(fOutFreq)//"_", nowyearf, nowmonf, &
               "_d", idom, ".nc"
            call greg2jul(0, 0, 0, 1, nowmon, nowyear, julian_st)
            call greg2jul(0, 0, 0, 1, nowmonf, nowyearf, julian_stf)
            djulian = int(nowjulian - julian_st, 4)
            djulianf = int(julianf - julian_stf, 4)
         else if (trim(fOutFreq) .eq. "yearly") then
            write (filename, "(a,i4.4,a,i2.2,a)") trim(fout_dir)//"MaCOM_"// &
               trim(fOutFreq)//"_", nowyearf, "_d", idom, ".nc"
            call greg2jul(0, 0, 0, 1, 1, nowyear, julian_st)
            call greg2jul(0, 0, 0, 1, 1, nowyearf, julian_stf)
            djulian = int(nowjulian - julian_st, 4)
            djulianf = int(julianf - julian_stf, 4)
         else
            write (*, *) ' ERROR! : Model only support daily, monthly or', &
               ' yearly now!, my rank is ', mpi_rank
            stop
         end if

         if (mod(djulian, nOutFreq) .eq. 0) flag_write_out = .true.

         if (flag_write_out) ts = int(djulianf/nOutFreq, 4) + 1

         if (flag_write_out) then

            IF (mpi_rank == mpi_comp_procs) THEN

               CALL mpi_csp_io_recv_output

               call netcdf_write(trim(filename), "time", dble(julianf), &
                                 'dp', ts=ts, nc_attr1=nc_attr1)

               CALL mpi_csp_io_root_scatterv_output(filename, ts)
            ELSE
               CALL mpi_csp_io_scatterv_output(filename, ts)
            END IF

         end if
      end do
      !---^^^ write regular output field finished

   END SUBROUTINE mpi_csp_io_main

end module mod_mpi_csp_io
