module mod_mpi_csp_io_misc
   use mod_misc
   use mod_mpi_variables
   use mod_mpi_interfaces
   implicit none

contains

!==============================================================================
   SUBROUTINE mpi_csp_io_comp_init
!==============================================================================
      mpi_total_2d_nlpb_nk = loc_nlpb * nk
      ALLOCATE (mpi_comp_buf_1d_ssh(loc_nlpb))
      ALLOCATE (mpi_comp_buf_1d_etaH(loc_nlpb))
      ALLOCATE (mpi_comp_buf_1d_qns(loc_nlpb))
      ALLOCATE (mpi_comp_buf_1d_qsr(loc_nlpb))
      ALLOCATE (mpi_comp_buf_1d_utau(loc_nlpb))
      ALLOCATE (mpi_comp_buf_1d_vtau(loc_nlpb))
      ALLOCATE (mpi_comp_buf_1d_EmPmR(loc_nlpb))
      ALLOCATE (mpi_comp_buf_2d_tFld(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_sFld(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_uFld(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_vFld(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_wFld(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_KappaRT(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_KappaRM(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_1d_restart(loc_nlpb))
      ALLOCATE (mpi_comp_buf_2d_restart(loc_nlpb, nk))

      ALLOCATE (mpi_comp_buf_2d_tFld_adv(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_sFld_adv(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_tFld_hdiff(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_sFld_hdiff(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_tFld_vdiff(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_sFld_vdiff(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_tFld_force(loc_nlpb, nk))
      ALLOCATE (mpi_comp_buf_2d_sFld_force(loc_nlpb, nk))

   END SUBROUTINE mpi_csp_io_comp_init

!==============================================================================
   SUBROUTINE mpi_csp_io_root_init
!==============================================================================
      INTEGER :: i
      mpi_total_2d_nlpb_nk = mpi_total_nlpb*nk
      ALLOCATE (mpi_io_buf_2d_tFld(mpi_total_2d_nlpb_nk, ntracer))
      ALLOCATE (mpi_io_buf_2d_uFld(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_vFld(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_wFld(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_KappaRT(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_KappaRM(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_1d_ssh(mpi_total_nlpb))
      ALLOCATE (mpi_io_buf_1d_etaH(mpi_total_nlpb))
      ALLOCATE (mpi_io_buf_1d_qns(mpi_total_nlpb))
      ALLOCATE (mpi_io_buf_1d_qsr(mpi_total_nlpb))
      ALLOCATE (mpi_io_buf_1d_EmPmR(mpi_total_nlpb))
      ALLOCATE (mpi_io_buf_1d_utau(mpi_total_nlpb))
      ALLOCATE (mpi_io_buf_1d_vtau(mpi_total_nlpb))
      ALLOCATE (mpi_io_buf_adjust_1d(mpi_total_nlpb))
      ALLOCATE (mpi_io_buf_adjust_2d(mpi_total_nlpb, nk))
      ALLOCATE (mpi_1d_nlpb_counts_displs_all(mpi_comp_procs + 1))
      ALLOCATE (mpi_2d_nlpb_nk_counts_displs_all(mpi_comp_procs + 1))
      ALLOCATE (mpi_2d_nlpb_nk_counts_all(mpi_comp_procs + 1))
      ALLOCATE (mpi_io_buf_1d_send_counts(mpi_io_procs))
      ALLOCATE (mpi_io_buf_1d_send_displs(mpi_io_procs))
      ALLOCATE (mpi_io_buf_2d_restart(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_1d_restart(mpi_total_nlpb))

      ALLOCATE (mpi_io_buf_2d_tFld_adv(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_sFld_adv(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_tFld_hdiff(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_sFld_hdiff(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_tFld_vdiff(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_sFld_vdiff(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_tFld_force(mpi_total_2d_nlpb_nk))
      ALLOCATE (mpi_io_buf_2d_sFld_force(mpi_total_2d_nlpb_nk))

      DO i = 0, mpi_io_procs - 1, 1
         IF (i < MOD(mpi_total_nlpb, mpi_io_procs)) THEN
            mpi_io_buf_1d_send_counts(i + 1) = mpi_total_nlpb/mpi_io_procs + 1
         ELSE
            mpi_io_buf_1d_send_counts(i + 1) = mpi_total_nlpb/mpi_io_procs
         END IF
      END DO

      mpi_2d_nlpb_nk_counts_all(1:mpi_comp_procs) = mpi_nlpb_counts_all(1:mpi_comp_procs)*nk
      ! the last element is io root process, not needed to receive data
      mpi_2d_nlpb_nk_counts_all(mpi_comp_procs + 1) = 0

      mpi_2d_nlpb_nk_counts_displs_all(1) = 0
      mpi_1d_nlpb_counts_displs_all(1) = 0
      DO i = 2, mpi_comp_procs, 1
         mpi_2d_nlpb_nk_counts_displs_all(i) = mpi_2d_nlpb_nk_counts_displs_all(i - 1) + mpi_nlpb_counts_all(i - 1)*nk
         mpi_1d_nlpb_counts_displs_all(i) = mpi_1d_nlpb_counts_displs_all(i - 1) + mpi_nlpb_counts_all(i - 1)
      END DO
      ! this element is fake displs
      mpi_2d_nlpb_nk_counts_displs_all(mpi_comp_procs + 1) = mpi_2d_nlpb_nk_counts_displs_all(mpi_comp_procs)
      mpi_1d_nlpb_counts_displs_all(mpi_comp_procs + 1) = mpi_1d_nlpb_counts_displs_all(mpi_comp_procs)

      IF (mpi_io_rank < MOD(nk, mpi_io_procs)) THEN
         mpi_io_nk = nk/mpi_io_procs + 1
      ELSE
         mpi_io_nk = nk/mpi_io_procs
      END IF

      !allocate data in 1 dimension
      mpi_io_buf_seperate_1d_count = mpi_io_buf_1d_send_counts(1)
      mpi_io_buf_seperate_1d_start = 1
      ALLOCATE (mpi_io_buf_seperate_1d(mpi_io_buf_seperate_1d_count))

      ALLOCATE (mpi_io_buf_seperate_2d(mpi_total_nlpb, mpi_io_nk))

      ALLOCATE (mpi_io_buf_2d_send_counts(mpi_io_procs))
      DO i = 0, mpi_io_procs - 1, 1
         IF (i < MOD(nk, mpi_io_procs)) THEN
            mpi_io_buf_2d_send_counts(i + 1) = (nk/mpi_io_procs + 1)*mpi_total_nlpb
         ELSE
            mpi_io_buf_2d_send_counts(i + 1) = (nk/mpi_io_procs)*mpi_total_nlpb
         END IF
      END DO

      ALLOCATE (mpi_io_buf_2d_displs(mpi_io_procs))
      mpi_io_buf_2d_displs(1) = 0
      mpi_io_buf_1d_send_displs(1) = 0
      IF (mpi_io_procs > 1) THEN
         DO i = 2, mpi_io_procs, 1
            mpi_io_buf_1d_send_displs(i) = mpi_io_buf_1d_send_displs(i - 1) &
                                           + mpi_io_buf_1d_send_counts(i - 1)
            mpi_io_buf_2d_displs(i) = mpi_io_buf_2d_displs(i - 1) &
                                      + mpi_io_buf_2d_send_counts(i - 1)
         END DO
      END IF

      mpi_io_buf_seperate_2d_start_nk(1) = 1
      IF (mpi_io_rank < MOD(nk, mpi_io_procs)) THEN
         mpi_io_buf_seperate_2d_start_nk(2) = mpi_io_rank*(nk/mpi_io_procs) &
                                              + mpi_io_rank + 1
      ELSE
         mpi_io_buf_seperate_2d_start_nk(2) = mpi_io_rank*(nk/mpi_io_procs) &
                                              + MOD(nk, mpi_io_procs) + 1
      END IF

      mpi_io_buf_seperate_2d_count_nk(1) = mpi_total_nlpb
      mpi_io_buf_seperate_2d_count_nk(2) = mpi_io_nk

   END SUBROUTINE mpi_csp_io_root_init

!==============================================================================
   SUBROUTINE mpi_csp_io_non_root_init
!==============================================================================

      INTEGER :: i
      IF (mpi_io_rank < MOD(nk, mpi_io_procs)) THEN
         mpi_io_nk = nk/mpi_io_procs + 1
      ELSE
         mpi_io_nk = nk/mpi_io_procs
      END IF

      IF (mpi_io_rank < MOD(mpi_total_nlpb, mpi_io_procs)) THEN
         mpi_io_buf_seperate_1d_count = mpi_total_nlpb/mpi_io_procs + 1
         mpi_io_buf_seperate_1d_start = mpi_io_rank*(mpi_total_nlpb/mpi_io_procs) &
                                        + mpi_io_rank + 1
      ELSE
         mpi_io_buf_seperate_1d_count = mpi_total_nlpb/mpi_io_procs
         mpi_io_buf_seperate_1d_start = mpi_io_rank*(mpi_total_nlpb/mpi_io_procs) &
                                        + MOD(mpi_total_nlpb, mpi_io_procs) + 1
      END IF

      ALLOCATE (mpi_io_buf_seperate_1d(mpi_io_buf_seperate_1d_count))
      ALLOCATE (mpi_io_buf_seperate_2d(mpi_total_nlpb, mpi_io_nk))

      ALLOCATE (mpi_io_buf_2d_send_counts(mpi_io_procs))
      DO i = 0, mpi_io_procs - 1, 1
         IF (i < MOD(nk, mpi_io_procs)) THEN
            mpi_io_buf_2d_send_counts(i + 1) = (nk/mpi_io_procs + 1)*mpi_total_nlpb
         ELSE
            mpi_io_buf_2d_send_counts(i + 1) = (nk/mpi_io_procs)*mpi_total_nlpb
         END IF
      END DO

      ALLOCATE (mpi_io_buf_2d_displs(mpi_io_procs))
      mpi_io_buf_2d_displs(1) = 0
      IF (mpi_io_procs > 1) THEN
         DO i = 2, mpi_io_procs, 1
            mpi_io_buf_2d_displs(i) = mpi_io_buf_2d_displs(i - 1) &
                                      + mpi_io_buf_2d_send_counts(i - 1)
         END DO
      END IF

      mpi_io_buf_seperate_2d_start_nk(1) = 1
      IF (mpi_io_rank < MOD(nk, mpi_io_procs)) THEN
         mpi_io_buf_seperate_2d_start_nk(2) = mpi_io_rank*(nk/mpi_io_procs) &
                                              + mpi_io_rank + 1
      ELSE
         mpi_io_buf_seperate_2d_start_nk(2) = mpi_io_rank*(nk/mpi_io_procs) &
                                              + MOD(nk, mpi_io_procs) + 1
      END IF

      mpi_io_buf_seperate_2d_count_nk(1) = mpi_total_nlpb
      mpi_io_buf_seperate_2d_count_nk(2) = mpi_io_nk

   END SUBROUTINE mpi_csp_io_non_root_init

!==============================================================================
   SUBROUTINE mpi_send_info_comp_to_io
!==============================================================================
      mpi_buf_size = max_domain*100*8
      ALLOCATE (mpi_buf(mpi_buf_size))
      mpi_position = 0
      CALL MPI_PACK(ndomain, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(idom, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nOutFreq, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nResFreq, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nRatio, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nIter0, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nIterMax, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(dTsurf, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(dTmom, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(dTtracer, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(restart_in, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(restart_out, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(startyear, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(startmon, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(startday, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(starthour, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(startmin, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(startsec, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(fOutFreq, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(fout_dir, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_SEND(mpi_position, 1, MPI_INTEGER, mpi_root_comp_io, mpi_comp_io_tag, mpi_comp_io_comm, mpi_err)
      CALL MPI_SEND(mpi_buf, mpi_position, MPI_PACKED, mpi_root_comp_io, mpi_comp_io_tag, mpi_comp_io_comm, mpi_err)
      DEALLOCATE (mpi_buf)

   END SUBROUTINE mpi_send_info_comp_to_io

!==============================================================================
   SUBROUTINE mpi_recv_info_io_from_comp
!==============================================================================
      INTEGER :: i = 0

      CALL MPI_RECV(mpi_buf_size, 1, MPI_INTEGER, mpi_sub_root_comp_io, &
                    mpi_comp_io_tag, mpi_comp_io_comm, mpi_comp_io_status, mpi_err)
      ALLOCATE (mpi_buf(mpi_buf_size))
      CALL MPI_RECV(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_sub_root_comp_io, &
                    mpi_comp_io_tag, mpi_comp_io_comm, mpi_comp_io_status, mpi_err)
      mpi_position = 0
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ndomain, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, idom, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nOutFreq, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nResFreq, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nRatio, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nIter0, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nIterMax, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dTsurf, 1, mpi_real_wp, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dTmom, 1, mpi_real_wp, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dTtracer, 1, mpi_real_wp, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, restart_in, 1, MPI_LOGICAL, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, restart_out, 1, MPI_LOGICAL, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowyear, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowmon, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowday, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowhour, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowmin, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowsec, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fOutFreq, lc, MPI_CHARACTER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fout_dir, lc, MPI_CHARACTER, mpi_comp_io_comm, mpi_err)
      DEALLOCATE (mpi_buf)

      ! the last element of mpi_nlpb_counts_all is not needed, but you have to create it for MPI_GATHER
      ! the follow MPI_GATHER and MPI_GATHERV is pair with last code in mpi_init_1d_index_buf
      ALLOCATE (mpi_nlpb_counts_all(mpi_comp_io_procs))
      CALL MPI_RECV(mpi_nlpb_counts_all, mpi_comp_procs, MPI_INTEGER, 0, mpi_tag1, MPI_COMM_WORLD, mpi_stat, mpi_err)
      mpi_nlpb_counts_all(mpi_comp_io_procs) = 0

      ! mpi_nlpb_block_1d_all is not used in mpi_comp_io_comm, so the size is equal to mpi_comp_procs not mpi_comp_procs + 1
      ALLOCATE (mpi_nlpb_block_1d_all(mpi_comp_procs))
      ALLOCATE (mpi_nlpb_block_1d_displs_all(mpi_comp_procs))
      ! i is not needed, but have to send with MPI_GATHER
      CALL MPI_RECV(mpi_nlpb_block_1d_all, mpi_comp_procs, MPI_INTEGER, 0, mpi_tag1, MPI_COMM_WORLD, mpi_stat, mpi_err)

      mpi_nlpb_block_1d_displs_all(1) = 0
      DO i = 1, mpi_comp_procs - 1, 1
         mpi_nlpb_block_1d_displs_all(i + 1) = mpi_nlpb_block_1d_displs_all(i) + mpi_nlpb_block_1d_all(i)
      END DO

      mpi_nlpb_block_1d_sum = mpi_nlpb_block_1d_displs_all(mpi_comp_procs) + mpi_nlpb_block_1d_all(mpi_comp_procs)

      ALLOCATE (mpi_nlpb_block_starts_1d_all(mpi_nlpb_block_1d_sum))
      ALLOCATE (mpi_nlpb_block_counts_1d_all(mpi_nlpb_block_1d_sum))

      CALL MPI_RECV(mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_1d_sum, MPI_INTEGER, 0, mpi_tag1, MPI_COMM_WORLD, &
                    mpi_stat, mpi_err)
      CALL MPI_RECV(mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_sum, MPI_INTEGER, 0, mpi_tag1, MPI_COMM_WORLD, &
                    mpi_stat, mpi_err)


   END SUBROUTINE mpi_recv_info_io_from_comp

!==============================================================================
   SUBROUTINE mpi_send_info_io
!==============================================================================
      mpi_buf_size = max_domain*100*8
      ALLOCATE (mpi_buf(mpi_buf_size))
      mpi_position = 0
      CALL MPI_PACK(ndomain, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(idom, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nOutFreq, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nRatio, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nIter0, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nIterMax, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(dTsurf, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(dTmom, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(dTtracer, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(restart_in, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowyear, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowmon, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowday, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowhour, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowmin, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowsec, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(fOutFreq, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(fout_dir, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_BCAST(mpi_buf_size, 1, MPI_INTEGER, mpi_root_io, mpi_io_comm, mpi_err)
      CALL MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_root_io, mpi_io_comm, mpi_err)
      DEALLOCATE (mpi_buf)

   END SUBROUTINE mpi_send_info_io

!==============================================================================
   SUBROUTINE mpi_recv_info_io
!==============================================================================
      CALL MPI_BCAST(mpi_buf_size, 1, MPI_INTEGER, mpi_root_io, mpi_io_comm, mpi_err)
      ALLOCATE (mpi_buf(mpi_buf_size))
      CALL MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_root_io, mpi_io_comm, mpi_err)
      mpi_position = 0
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ndomain, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, idom, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nOutFreq, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nRatio, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nIter0, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nIterMax, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dTsurf, 1, mpi_real_wp, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dTmom, 1, mpi_real_wp, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dTtracer, 1, mpi_real_wp, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, restart_in, 1, MPI_LOGICAL, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowyear, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowmon, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowday, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowhour, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowmin, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowsec, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fOutFreq, lc, MPI_CHARACTER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fout_dir, lc, MPI_CHARACTER, mpi_io_comm, mpi_err)
      DEALLOCATE (mpi_buf)

   END SUBROUTINE mpi_recv_info_io

!==============================================================================
   SUBROUTINE mpi_send_2nd_info_comp_to_io
!==============================================================================
      mpi_buf_size = 100*8
      ALLOCATE (mpi_buf(mpi_buf_size))
      mpi_position = 0
      CALL MPI_PACK(nk, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(ntracer, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nowyear, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nowmon, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nowday, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nowhour, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nowmin, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(nowsec, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_PACK(mpi_total_nlpb, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_io_comm, mpi_err)
      CALL MPI_SEND(mpi_position, 1, MPI_INTEGER, mpi_root_comp_io, mpi_comp_io_tag, mpi_comp_io_comm, mpi_err)
      CALL MPI_SEND(mpi_buf, mpi_position, MPI_PACKED, mpi_root_comp_io, mpi_comp_io_tag, mpi_comp_io_comm, mpi_err)
      DEALLOCATE (mpi_buf)

   END SUBROUTINE mpi_send_2nd_info_comp_to_io

!==============================================================================
   SUBROUTINE mpi_recv_2nd_info_io_from_comp
!==============================================================================
      CALL MPI_RECV(mpi_buf_size, 1, MPI_INTEGER, mpi_sub_root_comp_io, &
                    mpi_comp_io_tag, mpi_comp_io_comm, mpi_comp_io_status, mpi_err)
      ALLOCATE (mpi_buf(mpi_buf_size))
      CALL MPI_RECV(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_sub_root_comp_io, &
                    mpi_comp_io_tag, mpi_comp_io_comm, mpi_comp_io_status, mpi_err)
      mpi_position = 0
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nk, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ntracer, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowyear, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowmon, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowday, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowhour, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowmin, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowsec, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, mpi_total_nlpb, 1, MPI_INTEGER, mpi_comp_io_comm, mpi_err)
      DEALLOCATE (mpi_buf)

   END SUBROUTINE mpi_recv_2nd_info_io_from_comp

!==============================================================================
   SUBROUTINE mpi_send_2nd_info_io
!==============================================================================
      mpi_buf_size = 100*8
      ALLOCATE (mpi_buf(mpi_buf_size))
      mpi_position = 0
      CALL MPI_PACK(nk, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(ntracer, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowyear, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowmon, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowday, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowhour, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowmin, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(nowsec, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_PACK(mpi_total_nlpb, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_io_comm, mpi_err)
      CALL MPI_BCAST(mpi_buf_size, 1, MPI_INTEGER, mpi_root_io, mpi_io_comm, mpi_err)
      CALL MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_root_io, mpi_io_comm, mpi_err)
      DEALLOCATE (mpi_buf)
   END SUBROUTINE mpi_send_2nd_info_io

!==============================================================================
   SUBROUTINE mpi_recv_2nd_info_io
!==============================================================================
      CALL MPI_BCAST(mpi_buf_size, 1, MPI_INTEGER, mpi_root_io, mpi_io_comm, mpi_err)
      ALLOCATE (mpi_buf(mpi_buf_size))
      CALL MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_root_io, mpi_io_comm, mpi_err)
      mpi_position = 0
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nk, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ntracer, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowyear, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowmon, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowday, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowhour, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowmin, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nowsec, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, mpi_total_nlpb, 1, MPI_INTEGER, mpi_io_comm, mpi_err)
      DEALLOCATE (mpi_buf)

   END SUBROUTINE mpi_recv_2nd_info_io

!==============================================================================
   SUBROUTINE mpi_csp_io_send_output_1d(mpi_buf_source, mpi_io_buf_1d, mpi_req, &
                                        mpi_flag)
!==============================================================================
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_io_buf_1d
      REAL(wp), DIMENSION(:), INTENT(IN) :: mpi_buf_source
      INTEGER, INTENT(INOUT) :: mpi_req
      LOGICAL, INTENT(INOUT) :: mpi_flag

      IF (mpi_flag) THEN
         CALL MPI_WAIT(mpi_req, mpi_stat, mpi_err)
         mpi_flag = .FALSE.
      END IF

      mpi_io_buf_1d(:) = mpi_buf_source(1:loc_nlpb)

      CALL mpi_igatherv_comp_to_io_1d(mpi_io_buf_1d, loc_nlpb, mpi_root_comp_io, &
                                      mpi_req, mpi_comp_io_comm)

      ! after sending , set mpi_flag=.TRUE., so we need to check whether it is
      ! finished before send the second round of tFld
      mpi_flag = .TRUE.

   END SUBROUTINE mpi_csp_io_send_output_1d

!==============================================================================
   SUBROUTINE mpi_igatherv_comp_to_io_1d(mpi_send_buf, mpi_send_count, &
                                         mpi_igatherv_root, mpi_igatherv_req, mpi_igatherv_comm)
!==============================================================================
      REAL(wp), INTENT(IN), DIMENSION(:) :: mpi_send_buf
      INTEGER, INTENT(IN) :: mpi_send_count
      INTEGER, INTENT(IN) :: mpi_igatherv_root
      INTEGER, INTENT(IN) :: mpi_igatherv_comm
      INTEGER, INTENT(INOUT) :: mpi_igatherv_req
      REAL(wp) :: mpi_fake_recvbuf(1)
      INTEGER :: mpi_fake_recvcounts(1)
      INTEGER :: mpi_fake_displs(1)

      mpi_fake_recvbuf(1) = 0.0_wp
      mpi_fake_recvcounts(1) = 0
      mpi_fake_displs(1) = 0

      CALL MPI_IGATHERV(mpi_send_buf, mpi_send_count, mpi_real_wp, &
                        mpi_fake_recvbuf, mpi_fake_recvcounts, mpi_fake_displs, mpi_real_wp, &
                        mpi_igatherv_root, mpi_igatherv_comm, mpi_igatherv_req, mpi_err)
   END SUBROUTINE mpi_igatherv_comp_to_io_1d

!==============================================================================
   SUBROUTINE mpi_csp_io_send_output_2d(mpi_buf_source, mpi_io_buf_2d, mpi_req, mpi_flag)
!==============================================================================
      REAL(wp), INTENT(INOUT), DIMENSION(:, :) :: mpi_io_buf_2d
      REAL(wp), INTENT(IN), DIMENSION(:, :) :: mpi_buf_source
      INTEGER, INTENT(INOUT) :: mpi_req
      LOGICAL, INTENT(INOUT) :: mpi_flag

      IF (mpi_flag) THEN
         CALL MPI_WAIT(mpi_req, mpi_stat, mpi_err)
         mpi_flag = .FALSE.
      END IF

      mpi_io_buf_2d(:, :) = mpi_buf_source(1:loc_nlpb, 1:nk)
      CALL mpi_igatherv_comp_to_io(mpi_io_buf_2d, mpi_total_2d_nlpb_nk, &
                                   mpi_root_comp_io, mpi_req, mpi_comp_io_comm)

      ! after sending fFld, set mpi_flag_tFld=.TRUE., so we need to check whether
      ! it is finished before send the second round of tFld
      mpi_flag = .TRUE.

   END SUBROUTINE mpi_csp_io_send_output_2d

!==============================================================================
   SUBROUTINE mpi_igatherv_comp_to_io(mpi_send_buf, mpi_send_count, &
                                      mpi_igatherv_root, mpi_igatherv_req, mpi_igatherv_comm)
!==============================================================================
      implicit none

      REAL(wp), INTENT(IN), DIMENSION(:, :) :: mpi_send_buf
      INTEGER, INTENT(IN) :: mpi_send_count
      INTEGER, INTENT(IN) :: mpi_igatherv_root
      INTEGER, INTENT(IN) :: mpi_igatherv_comm
      INTEGER, INTENT(INOUT) :: mpi_igatherv_req

      REAL(wp) :: mpi_fake_recvbuf(1)
      INTEGER :: mpi_fake_recvcounts(1)
      INTEGER :: mpi_fake_displs(1)

      mpi_fake_recvbuf(1) = 0.0_wp
      mpi_fake_recvcounts(1) = 0
      mpi_fake_displs(1) = 0

      CALL MPI_IGATHERV(mpi_send_buf, mpi_send_count, mpi_real_wp, &
                        mpi_fake_recvbuf, mpi_fake_recvcounts, mpi_fake_displs, mpi_real_wp, &
                        mpi_igatherv_root, mpi_igatherv_comm, mpi_igatherv_req, mpi_err)
   END SUBROUTINE mpi_igatherv_comp_to_io

!==============================================================================
   SUBROUTINE mpi_igatherv_root_comp_to_io(mpi_recv_buf, mpi_recv_count, &
                                           mpi_recv_displs, mpi_igatherv_root, mpi_igatherv_comm, mpi_igatherv_req)
!==============================================================================
      REAL(wp), INTENT(IN), DIMENSION(:) :: mpi_recv_buf
      INTEGER, INTENT(IN), DIMENSION(:) :: mpi_recv_count
      INTEGER, INTENT(IN) :: mpi_igatherv_root
      INTEGER, INTENT(IN), DIMENSION(:) :: mpi_recv_displs
      INTEGER, INTENT(IN) :: mpi_igatherv_comm
      INTEGER, INTENT(INOUT) :: mpi_igatherv_req
      REAL(wp) :: mpi_fake_send_buf = 0.0
      INTEGER :: mpi_fake_send_count = 0
      CALL MPI_IGATHERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_real_wp, &
                        mpi_recv_buf, mpi_recv_count, mpi_recv_displs, mpi_real_wp, &
                        mpi_igatherv_root, mpi_igatherv_comm, mpi_igatherv_req, mpi_err)
   END SUBROUTINE mpi_igatherv_root_comp_to_io

!==============================================================================
   SUBROUTINE mpi_2d_data_adjust(buf_original, buf_adjust, buf_num, buf_blocks, &
                                 buf_starts, buf_counts, buf_displs)
!==============================================================================
      REAL(wp), INTENT(IN), DIMENSION(:) :: buf_original
      REAL(wp), INTENT(INOUT), DIMENSION(:, :) :: buf_adjust
      INTEGER, INTENT(IN) :: buf_num
      INTEGER, INTENT(IN), DIMENSION(:) :: buf_blocks
      INTEGER, INTENT(IN), DIMENSION(:) :: buf_starts
      INTEGER, INTENT(IN), DIMENSION(:) :: buf_counts
      INTEGER, INTENT(IN), DIMENSION(:) :: buf_displs
      INTEGER :: i, j, k, index_original_start, index_original_end, index_adjust_start, index_adjust_end

      index_original_start = 1
      DO j = 1, buf_num, 1
         DO i = 1, nk, 1
            DO k = 1, buf_blocks(j), 1
               index_adjust_start = buf_starts(buf_displs(j) + k)
               index_adjust_end = index_adjust_start + buf_counts(buf_displs(j) + k) - 1
               index_original_end = index_original_start + buf_counts(buf_displs(j) + k) - 1
               buf_adjust(index_adjust_start:index_adjust_end, i) = &
                  buf_original(index_original_start:index_original_end)
               index_original_start = index_original_end + 1
            END DO
         END DO
      END DO

   END SUBROUTINE mpi_2d_data_adjust

!==============================================================================
   SUBROUTINE mpi_1d_data_adjust(buf_original, buf_adjust, buf_num, buf_blocks, &
                                 buf_starts, buf_counts, buf_displs)
!==============================================================================
      REAL(wp), INTENT(IN), DIMENSION(:) :: buf_original
      REAL(wp), INTENT(INOUT), DIMENSION(:) :: buf_adjust
      INTEGER, INTENT(IN) :: buf_num
      INTEGER, INTENT(IN), DIMENSION(:) :: buf_blocks
      INTEGER, INTENT(IN), DIMENSION(:) :: buf_starts
      INTEGER, INTENT(IN), DIMENSION(:) :: buf_counts
      INTEGER, INTENT(IN), DIMENSION(:) :: buf_displs
      INTEGER ::  j, k, index_original_start, index_original_end, index_adjust_start, index_adjust_end

      index_original_start = 1
      DO j = 1, buf_num, 1
         DO k = 1, buf_blocks(j), 1
            index_adjust_start = buf_starts(buf_displs(j) + k)
            index_adjust_end = index_adjust_start + buf_counts(buf_displs(j) + k) - 1
            index_original_end = index_original_start + buf_counts(buf_displs(j) + k) - 1
            buf_adjust(index_adjust_start:index_adjust_end) = &
               buf_original(index_original_start:index_original_end)
            index_original_start = index_original_end + 1
         END DO
      END DO

   END SUBROUTINE mpi_1d_data_adjust

end module mod_mpi_csp_io_misc
