MODULE mod_mpi_test
   use mod_mpi_test_variables
   USE mod_mpi_variables
   use mod_mpi_interfaces
   use mod_mpi_csp_io_misc
   use mod_io_netcdf
   use mpi
   use mod_csp_basic, only: nowjulian
   IMPLICIT NONE


CONTAINS

   ! initializing operation for mpi debug and test
   SUBROUTINE mpi_test_debug_init
      ALLOCATE (mpi_variable_3d_nk_2(nlpbz, nk, 2))
      ALLOCATE (mpi_variable_2d_nk(nlpbz, nk))
      ALLOCATE (mpi_variable_1d(nlpbz))
      ALLOCATE (mpi_dp_1d(nlpb))
      ALLOCATE (mpi_dp_2d_nk(nlpb, nk))
   END SUBROUTINE

   ! debug netcdf reading and exchanging functions
   SUBROUTINE mpi_debug_netcdf_read_exchange
      CHARACTER(LEN = lc) :: filename
      filename ="uFld.nc"
      CALL mpi_netcdf_read_exchange(filename, 'uFld', uFld, 3)
      WRITE (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_maskW_1.txt"
!      CALL mpi_write_file_2d_dp(mpi_debug_filename,uFld,nlpb,nk)

      maskW(1:nlpb,1:nk)= INT(uFld(1:nlpb,1:nk))
      CALL mpi_write_file_2d_i2(mpi_debug_filename,maskW,nlpb,nk)
      maskW(loc_nlpb + 1:nlpb,1:nk) = - 1
      WRITE (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_maskW_2.txt"
      CALL mpi_write_file_2d_i2(mpi_debug_filename,maskW,nlpb,nk)
      CALL mpi_write_file_2d_dp(mpi_debug_filename,uFld,nlpb,nk)
    END SUBROUTINE

   ! time test
   SUBROUTINE mpi_test_time_start
! test in computing processes
   IF (mpi_rank == 0) THEN
      WRITE (mpi_test_filename, '(I3.3,A1,A,A)') mpi_procs, '_', TRIM(mpi_time), "_comp_test.txt"
      OPEN (UNIT = MPI_TEST_UNIT, FILE = TRIM(mpi_test_filename), STATUS ='REPLACE', ACTION ='WRITE')
      WRITE (MPI_TEST_UNIT, *) "Number of whole parallel processors is:", mpi_procs
      WRITE (MPI_TEST_UNIT, *) "Number of whole computing processors is:", mpi_comp_procs
      CALL SYSTEM_CLOCK(mpi_clock_start,mpi_clock_rate,mpi_clock_count_max)
      CALL SYSTEM_CLOCK(mpi_clock_start)
      CALL SYSTEM_CLOCK(mpi_clock_init_start)
   END IF
! test in io processes
!   IF (mpi_rank == mpi_comp_procs) THEN
!      CALL DATE_AND_TIME(TIME = mpi_time)
!      WRITE (mpi_test_filename, '(I3.3,A1,A,A)'), mpi_procs, '_', TRIM(mpi_time), "_io_test.txt"
!      OPEN (UNIT = MPI_TEST_IO_UNIT, FILE = TRIM(mpi_test_filename), STATUS ='REPLACE', ACTION ='WRITE')
!      WRITE (MPI_TEST_IO_UNIT, *) "Number of whole parallel processors is:", mpi_procs
!      WRITE (MPI_TEST_IO_UNIT, *) "Number of whole IO processors is:", mpi_io_procs
!      mpi_time_real_io_write = 0
!   END IF
   END SUBROUTINE

   ! time test
   SUBROUTINE mpi_test_time_end
   IF (mpi_rank == 0) THEN
     CALL SYSTEM_CLOCK(mpi_clock_end)
      WRITE (MPI_TEST_UNIT, *) "the whole execution time (seconds) for proc 0 is:"
      WRITE (MPI_TEST_UNIT, *) dble(mpi_clock_end - mpi_clock_start)/ dble(mpi_clock_rate)
      WRITE (MPI_TEST_UNIT, *) "the whole computing time (seconds) for proc 0 is:"
      WRITE (MPI_TEST_UNIT, *) dble(mpi_clock_compute_end - mpi_clock_compute_start)/ dble(mpi_clock_rate)
      WRITE (MPI_TEST_UNIT, *) "the whole initializing time (seconds) for proc 0 is:"
      WRITE (MPI_TEST_UNIT, *) dble(mpi_clock_init_end - mpi_clock_init_start)/ dble(mpi_clock_rate)
!      WRITE (MPI_TEST_UNIT, *) "the whole IO - reading time (seconds) for proc 0 is:"
!      WRITE (MPI_TEST_UNIT, *) mpi_time_io_read
!      WRITE (MPI_TEST_UNIT, *) "the whole IO - writing time (seconds) for proc 0 is:"
!      WRITE (MPI_TEST_UNIT, *) mpi_time_io_write
      CLOSE (MPI_TEST_UNIT)
   END IF
!   IF (mpi_rank == mpi_comp_procs) THEN
!      WRITE (MPI_TEST_IO_UNIT, *) "the real whole IO - writing time (seconds) for proc ", mpi_comp_procs, " is:", mpi_time_real_io_write
!      CLOSE (MPI_TEST_IO_UNIT)
!   END IF
   END SUBROUTINE

   ! record the begining time of computation
   SUBROUTINE mpi_test_compute_time_start
      IF (mpi_rank == 0) THEN
        CALL SYSTEM_CLOCK(mpi_clock_compute_start)
      END IF
   END SUBROUTINE
   ! test initializing time
   SUBROUTINE mpi_test_init_time_end
      IF (mpi_rank == 0) THEN
        CALL SYSTEM_CLOCK(mpi_clock_init_end)
      END IF
   END SUBROUTINE
   ! compute time of computation
   SUBROUTINE mpi_test_compute_time_end
      IF (mpi_rank == 0) THEN
        CALL SYSTEM_CLOCK(mpi_clock_compute_end)
      END IF
   END SUBROUTINE

   ! check difference between variables
   SUBROUTINE mpi_debug_diff_uv(uv)
      REAL(wp), DIMENSION(:, :, :) :: uv
      mpi_variable_3d_nk_2(1:nlpbz, 1:nk, 1:2) = uv(1:nlpbz, 1:nk, 1:2)
      CALL mpi_diff_wp_3d(uv, mpi_variable_3d_nk_2, nk, 2, TRIM("uv"))
      WRITE (*, *) "finished comparison for uv at rank:", mpi_rank
   END SUBROUTINE

   ! compare data difference of 3 dimension
   SUBROUTINE mpi_diff_wp_3d(mpi_array1, mpi_array2, dimen1, dimen2, mpi_varible_name)
      REAL(wp), DIMENSION(:, :, :), INTENT(IN) :: mpi_array1, mpi_array2
      INTEGER, INTENT(IN) :: dimen1, dimen2
      CHARACTER(len =*), INTENT(IN) :: mpi_varible_name
      INTEGER :: i, j, k
      DO i = 1, dimen2
         DO j = 1, dimen1
            DO k = 1, nlpbz
            if (mpi_array1(k, j, i) /= mpi_array2(k, j, i)) then
               WRITE (*, *) TRIM(mpi_varible_name), "(", k, ",", j, ",", i, ") is not equal, my rank id is", mpi_rank &
                  , ",the first value is", mpi_array1(k, j, i), ",the second value is", mpi_array2(k, j, i)
               STOP
            end if
            END DO
         END DO
      END DO
   END SUBROUTINE

   ! compare data difference of 2 dimension
   SUBROUTINE mpi_diff_wp_2d(mpi_array1, mpi_array2, dimen1, mpi_varible_name)
      REAL(wp), DIMENSION(:, :), INTENT(IN) :: mpi_array1, mpi_array2
      INTEGER, INTENT(IN) :: dimen1
      CHARACTER(len =*), INTENT(IN) :: mpi_varible_name
      INTEGER :: i, j
      DO i = 1, dimen1
         DO j = 1, nlpbz
            if (mpi_array1(j, i) /= mpi_array2(j, i)) then
               WRITE (*, *) TRIM(mpi_varible_name), "(", j, ",", i, ") is not equal, my rank id is", mpi_rank &
                  , ",the first value is", mpi_array1(j, i), ",the second value is", mpi_array2(j, i)
               STOP
            end if
         END DO
      END DO
   END SUBROUTINE

   ! compare data difference of 1 dimension
   SUBROUTINE mpi_diff_wp_1d(mpi_array1, mpi_array2, mpi_varible_name)
      REAL(wp), DIMENSION(:), INTENT(IN) :: mpi_array1, mpi_array2
      CHARACTER(len =*), INTENT(IN) :: mpi_varible_name
      INTEGER :: i
      DO i = 1, nlpbz
         if (mpi_array1(i) /= mpi_array2(i)) then
            WRITE (*, *) TRIM(mpi_varible_name), "(", i, ") is not equal, my rank id is", mpi_rank &
               , ",the first value is", mpi_array1(i), ",the second value is", mpi_array2(i)
            STOP
         end if
      END DO
   END SUBROUTINE

   ! output integer data of 1 dimension
   SUBROUTINE mpi_gather_output_integer(mpi_input, mpi_filename)
      INTEGER, INTENT(IN) :: mpi_input
      CHARACTER(lc), INTENT(IN) :: mpi_filename
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_int_array
      INTEGER :: i
      IF (mpi_rank == 0) THEN
         ALLOCATE (mpi_int_array(mpi_comp_procs))
      END IF

      CALL MPI_GATHER(mpi_input, 1, MPI_INTEGER, mpi_int_array, 1, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
      IF (mpi_rank == 0) THEN
         OPEN (UNIT = MPI_TEST_UNIT, FILE = TRIM(mpi_filename), STATUS ='REPLACE', ACTION ='WRITE')
         DO i = 1, mpi_comp_procs
            WRITE (MPI_TEST_UNIT, *) mpi_int_array(i)
         END DO
         CLOSE (MPI_TEST_UNIT)
      END IF
   END SUBROUTINE

   ! send dp data of 1 dimension from computing process to IO process
   subroutine mpi_send_output_1d_test(file_name, var, d1)
      character(lc), intent(in) :: file_name
      integer, intent(in) :: d1
      real(dp), intent(in) :: var(d1)
      real(dp) :: send_var(d1)

      IF (mpi_rank == 0) THEN
         CALL MPI_SEND(file_name, lc, MPI_CHARACTER, mpi_root_comp_io, mpi_comp_io_tag, mpi_comp_io_comm, mpi_err)
      END IF

      send_var(:) = var(1:d1)

      CALL mpi_gatherv_comp_to_io_1d(send_var, d1, mpi_root_comp_io, mpi_comp_io_comm)

   end subroutine mpi_send_output_1d_test

   ! send dp data of 2 dimension from computing process to IO process
   subroutine mpi_send_output_2d_test(file_name, var, d1, d2)
      character(lc), intent(in) :: file_name
      integer, intent(in) :: d1, d2
      real(dp), intent(in) :: var(d1, d2)
      real(dp) :: send_var(d1, d2)
      integer :: d1_d2

      IF (mpi_rank == 0) THEN
         CALL MPI_SEND(file_name, lc, MPI_CHARACTER, mpi_root_comp_io, mpi_comp_io_tag, mpi_comp_io_comm, mpi_err)
      END IF

      d1_d2 = d1 * d2
      send_var(:, :) = var(1:d1, 1:d2)

      CALL mpi_gatherv_comp_to_io(send_var, d1_d2, mpi_root_comp_io, mpi_comp_io_comm)

   end subroutine mpi_send_output_2d_test

   ! gather data from computing process to IO process
   SUBROUTINE mpi_gatherv_root_comp_to_io(mpi_recv_buf, mpi_recv_count, mpi_recv_displs, mpi_gatherv_root, mpi_gatherv_comm)
      REAL(wp), INTENT(IN), DIMENSION(:) :: mpi_recv_buf
      INTEGER, INTENT(IN), DIMENSION(:) :: mpi_recv_count
      INTEGER, INTENT(IN) :: mpi_gatherv_root
      INTEGER, INTENT(IN), DIMENSION(:) :: mpi_recv_displs
      INTEGER, INTENT(IN) :: mpi_gatherv_comm
      REAL(wp) :: mpi_fake_send_buf = 0.0
      INTEGER :: mpi_fake_send_count = 0
      CALL MPI_GATHERV(mpi_fake_send_buf, mpi_fake_send_count, mpi_real_wp, mpi_recv_buf, mpi_recv_count, mpi_recv_displs, &
                       mpi_real_wp, mpi_gatherv_root, mpi_gatherv_comm, mpi_err)
   END SUBROUTINE mpi_gatherv_root_comp_to_io

   ! gather data from computing process to IO process
   SUBROUTINE mpi_gatherv_comp_to_io(mpi_send_buf, mpi_send_count, mpi_gatherv_root, mpi_gatherv_comm)
      REAL(DP), INTENT(IN), DIMENSION(:, :) :: mpi_send_buf
      INTEGER, INTENT(IN) :: mpi_send_count
      INTEGER, INTENT(IN) :: mpi_gatherv_root
      INTEGER, INTENT(IN) :: mpi_gatherv_comm
      REAL(DP) :: mpi_fake_recvbuf(1) = 0.0
      INTEGER :: mpi_fake_recvcounts(1) = 0
      INTEGER :: mpi_fake_displs(1) = 0
      CALL MPI_GATHERV(mpi_send_buf, mpi_send_count, MPI_DOUBLE_PRECISION, mpi_fake_recvbuf, mpi_fake_recvcounts, mpi_fake_displs, &
                       MPI_DOUBLE_PRECISION, mpi_gatherv_root, mpi_gatherv_comm, mpi_err)
   END SUBROUTINE mpi_gatherv_comp_to_io

   ! gather data from computing process to IO process in 1 dimension
   SUBROUTINE mpi_gatherv_comp_to_io_1d(mpi_send_buf, mpi_send_count, mpi_gatherv_root, mpi_gatherv_comm)
      REAL(DP), INTENT(IN), DIMENSION(:) :: mpi_send_buf
      INTEGER, INTENT(IN) :: mpi_send_count
      INTEGER, INTENT(IN) :: mpi_gatherv_root
      INTEGER, INTENT(IN) :: mpi_gatherv_comm
      REAL(DP) :: mpi_fake_recvbuf(1) = 0.0
      INTEGER :: mpi_fake_recvcounts(1) = 0
      INTEGER :: mpi_fake_displs(1) = 0
      CALL MPI_GATHERV(mpi_send_buf, mpi_send_count, MPI_DOUBLE_PRECISION, mpi_fake_recvbuf, mpi_fake_recvcounts, mpi_fake_displs, &
                       MPI_DOUBLE_PRECISION, mpi_gatherv_root, mpi_gatherv_comm, mpi_err)
   END SUBROUTINE mpi_gatherv_comp_to_io_1d

   ! output variables
   SUBROUTINE mpi_debug_output_rhoInSitu_latC
!      mpi_var_filename = "second_rhoInSitu.txt"
!      CALL mpi_send_output_2d_test(mpi_var_filename, rhoInSitu(1:nlpb,1:nk), nlpb, nk)
      CALL mpi_data_2d_output_comp(TRIM("_rhoInSitu.txt"), rhoInSitu, nk)
!      mpi_var_filename = "second_latC.txt"
!      CALL mpi_send_output_1d_test(mpi_var_filename, latC, nlpb)
      CALL mpi_data_1d_output_comp(TRIM("_latC.txt"), latC)
   ! output data collecting data from computing processes
   END SUBROUTINE
   ! output variables
   SUBROUTINE mpi_debug_output_dPhiHydX_dPhiHydY(k)
     INTEGER,INTENT(IN)::k
     WRITE (mpi_debug_filename, '(I2.2,A)') k, "_dPhiHydX"
     CALL mpi_data_1d_output_comp(TRIM(mpi_debug_filename), dPhiHydX)
     WRITE (mpi_debug_filename, '(I2.2,A)') k, "_dphiSurfX"
     CALL mpi_data_1d_output_comp(TRIM(mpi_debug_filename), dphiSurfX)
   END SUBROUTINE

   SUBROUTINE mpi_data_1d_output_comp(file_name, send_var)
      character(*), intent(in) :: file_name
      real(wp), intent(in), dimension(:) :: send_var
      character(lc) :: real_name
      real(wp) :: recv_buf(mpi_total_nlpb)
      REAL(wp), ALLOCATABLE, DIMENSION(:) :: adjust_buf
      integer :: i

      CALL MPI_GATHERV(send_var, loc_nlpb, mpi_real_wp, recv_buf, mpi_nlpb_counts_all, &
                       mpi_nlpb_displs_all, mpi_real_wp, 0, mpi_comp_comm, mpi_err)

      IF (mpi_rank == 0) THEN
         ALLOCATE (adjust_buf(mpi_total_nlpb))
         CALL mpi_1d_data_adjust(recv_buf, adjust_buf, mpi_comp_procs, mpi_nlpb_block_1d_all, &
                                 mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         WRITE (real_name, '(I4.4,A)') myIter, TRIM(file_name)
!         OPEN (UNIT = 98, FILE = TRIM(real_name), STATUS ='REPLACE', ACTION ='WRITE')
!!     WRITE (98, *) adjust_buf
!         DO i = 1, mpi_total_nlpb
!            WRITE (98, *) adjust_buf(i)
!         END DO
!         CLOSE (98)
        CALL netcdf_write(TRIM(real_name),"var_name",adjust_buf,'dp',d1 = mpi_total_nlpb,ts = 1)
      END IF
   END SUBROUTINE

   SUBROUTINE mpi_data_1d_output_snapshot(prefix,varname, send_var)
       !YY: output single frame (1d)
      character(*), intent(in) :: prefix,varname
      real(wp), intent(in), dimension(:) :: send_var
      character(lc) :: real_name
      real(wp) :: recv_buf(mpi_total_nlpb)
      REAL(wp), ALLOCATABLE, DIMENSION(:) :: adjust_buf
      integer :: i

      CALL MPI_GATHERV(send_var, loc_nlpb, mpi_real_wp, recv_buf, mpi_nlpb_counts_all, &
                       mpi_nlpb_displs_all, mpi_real_wp, 0, mpi_comp_comm, mpi_err)

      IF (mpi_rank == 0) THEN
         ALLOCATE (adjust_buf(mpi_total_nlpb))
         CALL mpi_1d_data_adjust(recv_buf, adjust_buf, mpi_comp_procs, mpi_nlpb_block_1d_all, &
                                 mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         write(real_name,'(A,i4.4,i2.2,i2.2,i2.2,i2.2,i2.2,A)') trim(adjustl(prefix))//trim(adjustl(varname)),    &
             nowyearf,nowmonf,nowdayf,nowhourf,nowminf,nowsecf,'.nc'

        CALL netcdf_write(TRIM(real_name),varname,adjust_buf,'dp',d1 = mpi_total_nlpb,ts = 1)
      END IF
   END SUBROUTINE mpi_data_1d_output_snapshot

   SUBROUTINE mpi_data_1d_output_ice_ts(prefix,varname, send_var,ts)
      !YY: output frames in a NC file (1d)
      character(*), intent(in) :: prefix,varname
      real(wp), intent(in), dimension(:) :: send_var
      character(lc) :: real_name
      real(wp) :: recv_buf(mpi_total_nlpb)
      REAL(wp), ALLOCATABLE, DIMENSION(:) :: adjust_buf
      integer :: i, ts
      type(nc_attr) :: nc_attr1
      
      nc_attr1%var_units = 'seconds since 1949-10-01 00:00:00'

      CALL MPI_GATHERV(send_var, loc_nlpb, mpi_real_wp, recv_buf, mpi_nlpb_counts_all, &
                       mpi_nlpb_displs_all, mpi_real_wp, 0, mpi_comp_comm, mpi_err)

      IF (mpi_rank == 0) THEN
         ALLOCATE (adjust_buf(mpi_total_nlpb))
         CALL mpi_1d_data_adjust(recv_buf, adjust_buf, mpi_comp_procs, mpi_nlpb_block_1d_all, &
                                 mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         write(real_name,'(A)') trim(adjustl(prefix))//trim(adjustl(varname))//'.nc'

        call netcdf_write(trim(real_name), "time", dble(nowjulian), 'dp', ts=ts, nc_attr1=nc_attr1)
        CALL netcdf_write(TRIM(real_name),varname,adjust_buf,'dp',d1 = mpi_total_nlpb,ts = ts)
      END IF
   END SUBROUTINE mpi_data_1d_output_ice_ts

   ! output variables
   SUBROUTINE mpi_debug_output_gTracer_diffKhu_tracer_etc
      CALL mpi_data_2d_output_comp(TRIM("_gTracer.txt"), gTracer, nk)
      CALL mpi_data_2d_output_comp(TRIM("_diffKhu.txt"), diffKhu, nk)
      CALL mpi_data_2d_output_comp(TRIM("_tracer.txt"), tracer, nk)
      CALL mpi_data_2d_output_comp(TRIM("_tFld.txt"), tFld(:, :, 1), nk)
   END SUBROUTINE

   ! check difference between two variables
   SUBROUTINE mpi_debug_diff_uFld
      mpi_variable_2d_nk(1:nlpbz, 1:nk) = uFld(1:nlpbz, 1:nk)
      !CALL mpi_data_exchange(mpi_sendbuf_uFld, mpi_variable_2d_nk, nk)
      CALL mpi_diff_wp_2d(uFld, mpi_variable_2d_nk, nk, TRIM("uFld"))
      WRITE (*, *) "finished comparison for uFld at rank:", mpi_rank
   END SUBROUTINE

   ! check difference between two variables
   SUBROUTINE mpi_debug_diff_etaH
      mpi_variable_1d(1:nlpbz) = etaH(1:nlpbz)
      !CALL mpi_data_exchange(mpi_sendbuf_etaH, mpi_variable_1d)
      CALL mpi_diff_wp_1d(etaH, mpi_variable_1d, TRIM("etaH"))
      WRITE (*, *) "finished comparison for etaH at rank:", mpi_rank
   END SUBROUTINE

   ! output variables
   SUBROUTINE mpi_debug_output_EmPmR_tileEmP_rAc
      INTEGER :: i
      CALL mpi_data_1d_output_comp(TRIM("_EmPmR.txt"), EmPmR)
      mpi_dp_1d(1:nlpb) = dble(EmPmR(1:nlpb))
      CALL mpi_data_1d_output_comp(TRIM("_dble_EmPmR.txt"), mpi_dp_1d(1:nlpb))
      do i = 1, nlpb
         mpi_dp_1d(i) = dble(rAc(i))* dble(EmPmR(i))* maskC(i, nk)
      end do
      CALL mpi_data_1d_output_comp(TRIM("_tileEmP.txt"), mpi_dp_1d(1:nlpb))
      CALL mpi_data_1d_output_comp(TRIM("_rAc.txt"), rAc)
!      CALL mpi_data_2d_output_comp(TRIM("_maskC.txt"),maskC,nk)
   ! output variables
   END SUBROUTINE

   SUBROUTINE mpi_debug_output_guExt_phi0surf_dphiSurfX_etc
      CALL mpi_data_1d_output_comp(TRIM("_guExt.txt"), guExt)
      CALL mpi_data_1d_output_comp(TRIM("_phi0surf.txt"), phi0surf)
      CALL mpi_data_1d_output_comp(TRIM("_dphiSurfX.txt"), dphiSurfX)
      CALL mpi_data_1d_output_comp(TRIM("_uDragBottom.txt"), uDragBottom)
      CALL mpi_data_1d_output_comp(TRIM("_recip_hFacZ.txt"), recip_hFacZ)
      CALL mpi_data_1d_output_comp(TRIM("_dPhiHydY.txt"), dPhiHydY)
      CALL mpi_data_1d_output_comp(TRIM("_hDiv.txt"), hDiv)
      CALL mpi_data_1d_output_comp(TRIM("_vort3.txt"), vort3)
      CALL mpi_data_1d_output_comp(TRIM("_omega3.txt"), omega3)
      CALL mpi_data_1d_output_comp(TRIM("_del2u.txt"), del2u)
      CALL mpi_data_1d_output_comp(TRIM("_uCorTerm.txt"), uCorTerm)
      CALL mpi_data_1d_output_comp(TRIM("_dKEdx.txt"), dKEdx)
      CALL mpi_data_1d_output_comp(TRIM("_viscA4_D.txt"), viscA4_D)
      CALL mpi_data_1d_output_comp(TRIM("_uDissip.txt"), uDissip)
      CALL mpi_data_1d_output_comp(TRIM("_uDragTerms.txt"), uDragTerms)
      CALL mpi_data_1d_output_comp(TRIM("_uShearTerm.txt"), uShearTerm)
      CALL mpi_data_2d_output_comp(TRIM("_gUnow.txt"), gUnow, nk)
      CALL mpi_data_2d_output_comp(TRIM("_gTracer.txt"), gTracer, nk)
      CALL mpi_data_1d_output_comp(TRIM("_etaN_second.txt"), etaN)
      CALL mpi_data_2d_output_comp(TRIM("_uFld.txt"), uFld, nk)
      CALL mpi_data_2d_output_comp(TRIM("_vFld.txt"), vFld, nk)
      CALL mpi_data_2d_output_comp(TRIM("_wFld.txt"), wFld(:, 2:nkp1), nk)
      CALL mpi_data_2d_output_comp(TRIM("_hFacC.txt"), hFacC, nk)
    END SUBROUTINE

    ! output variables
    SUBROUTINE mpi_debug_output_diss_en_hmxl(diss_en)
      real(wp), INTENT(IN),dimension(:, :) :: diss_en
      CALL mpi_data_2d_output_comp(TRIM("_diss_en.txt"), diss_en, nk)
      CALL mpi_data_2d_output_comp(TRIM("_hmxl.txt"), hmxl, nk)
      CALL mpi_data_2d_output_comp(TRIM("_en.txt"), en, nk)
    END SUBROUTINE

    ! output variables
    SUBROUTINE mpi_debug_output_c_mu_work_c(work_c)
      real(wp), INTENT(IN),dimension(:, :) :: work_c
!      CALL mpi_data_2d_output_comp(TRIM("_c_mu.txt"), mpi_dp_2d_nk, nk)
      CALL mpi_data_2d_output_comp(TRIM("_work_c.txt"), work_c, nk)
    END SUBROUTINE

    ! output variables
    SUBROUTINE mpi_debug_output_y_c_d_u(y,c,d,u)
      real(wp), INTENT(IN),dimension(:, :) :: y,c,d,u
      CALL mpi_data_2d_output_comp(TRIM("_y.txt"), y, nk)
      CALL mpi_data_2d_output_comp(TRIM("_c.txt"), c, nk)
      CALL mpi_data_2d_output_comp(TRIM("_d.txt"), d, nk)
      CALL mpi_data_2d_output_comp(TRIM("_u.txt"), u, nk)
    END SUBROUTINE

   ! output data collecting data of 2 dimension from computing processes
   SUBROUTINE mpi_data_2d_output_comp(file_name, send_var, d1)
      character(*), intent(in) :: file_name
      integer, intent(in) :: d1
      real(wp), intent(in), dimension(:, :) :: send_var
      character(lc) :: real_name
      real(wp) :: send_buf(loc_nlpb, d1)
      real(wp) :: recv_buf(mpi_total_nlpb * d1)
      REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: adjust_buf
      integer :: i, j, send_count

      send_buf(1:loc_nlpb, 1:d1) = send_var(1:loc_nlpb, 1:d1)
      send_count = loc_nlpb * d1

      CALL MPI_GATHERV(send_buf, send_count, mpi_real_wp, recv_buf, mpi_2d_nlpb_nk_counts_all, &
                       mpi_2d_nlpb_nk_counts_displs_all, mpi_real_wp, 0, mpi_comp_comm, mpi_err)

      IF (mpi_rank == 0) THEN
         ALLOCATE (adjust_buf(mpi_total_nlpb, d1))
         CALL mpi_2d_data_adjust(recv_buf, adjust_buf, mpi_comp_procs, mpi_nlpb_block_1d_all, &
                                 mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
         WRITE (real_name, '(I4.4,A)') myIter, TRIM(file_name)

!         OPEN (UNIT = 98, FILE = TRIM(real_name), STATUS ='REPLACE', ACTION ='WRITE')
!!     WRITE (98, *) adjust_buf
!         DO j = 1, d1
!            DO i = 1, mpi_total_nlpb
!               WRITE (98, *) adjust_buf(i, j)
!            END DO
!         END DO
!         CLOSE (98)
        CALL netcdf_write(TRIM(real_name),"var_name",adjust_buf,'dp',d1 = mpi_total_nlpb, d2 = d1, ts = 1)
      END IF
   END SUBROUTINE

   ! output data in io process of 2 dimension from computing processes
   SUBROUTINE mpi_recv_write_output_2d_test(d1, d2)
      character(lc) :: file_name
      integer, intent(in) :: d1, d2
      real(dp) :: var(d1 * d2)
      real(dp) :: var_adjust(d1, d2)
      integer i, j
      ! first receive file name and variable name, then decide how to scatter them
      CALL MPI_RECV(file_name, lc, MPI_CHARACTER, 0, mpi_comp_io_tag, mpi_comp_io_comm, mpi_stat, mpi_err)
      CALL mpi_gatherv_root_comp_to_io(var, mpi_2d_nlpb_nk_counts_all, &
         mpi_2d_nlpb_nk_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm)

      CALL mpi_2d_data_adjust(var, var_adjust, mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      OPEN (UNIT = 98, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO j = 1, d2
         DO i = 1, d1
            WRITE (98, *) var_adjust(i, j)
         END DO
      END DO
      CLOSE (98)

   END SUBROUTINE mpi_recv_write_output_2d_test

   ! output data in io process of 1 dimension from computing processes
   SUBROUTINE mpi_recv_write_output_1d_test(d1)
      character(lc) :: file_name
      integer, intent(in) :: d1
      real(dp) :: var(d1)
      real(dp) :: var_adjust(d1)
      integer ::i
      ! first receive file name and variable name, then decide how to scatter them
      CALL MPI_RECV(file_name, lc, MPI_CHARACTER, 0, mpi_comp_io_tag, mpi_comp_io_comm, mpi_stat, mpi_err)
      CALL mpi_gatherv_root_comp_to_io(var, mpi_nlpb_counts_all, mpi_1d_nlpb_counts_displs_all, mpi_root_comp_io, mpi_comp_io_comm)

      CALL mpi_1d_data_adjust(var, var_adjust, mpi_comp_procs, mpi_nlpb_block_1d_all, mpi_nlpb_block_starts_1d_all, &
                              mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_displs_all)
      OPEN (UNIT = 98, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, d1
         WRITE (98, *) var_adjust(i)
      END DO
      CLOSE (98)

   END SUBROUTINE mpi_recv_write_output_1d_test

   ! write files
   SUBROUTINE mpi_debug_write_ssh
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_ssh.txt"
      call mpi_write_file_1d_dp(mpi_debug_filename, mpi_io_buf_adjust_1d, mpi_total_nlpb)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_t
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_t.txt"
      call mpi_write_file_2d_dp(mpi_debug_filename, mpi_io_buf_adjust_2d, mpi_total_nlpb, nk)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_s
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_s.txt"
      call mpi_write_file_2d_dp(mpi_debug_filename, mpi_io_buf_adjust_2d, mpi_total_nlpb, nk)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_u
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_u.txt"
      call mpi_write_file_2d_dp(mpi_debug_filename, mpi_io_buf_adjust_2d, mpi_total_nlpb, nk)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_v
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_v.txt"
      call mpi_write_file_2d_dp(mpi_debug_filename, mpi_io_buf_adjust_2d, mpi_total_nlpb, nk)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_w
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_w.txt"
      call mpi_write_file_2d_dp(mpi_debug_filename, mpi_io_buf_adjust_2d, mpi_total_nlpb, nk)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_kappa_h
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_kappa_h.txt"
      call mpi_write_file_2d_dp(mpi_debug_filename, mpi_io_buf_adjust_2d, mpi_total_nlpb, nk)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_kappa_m
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_kappa_m.txt"
      call mpi_write_file_2d_dp(mpi_debug_filename, mpi_io_buf_adjust_2d, mpi_total_nlpb, nk)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_pbt
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_pbt.txt"
      call mpi_write_file_1d_dp(mpi_debug_filename, mpi_io_buf_adjust_1d, mpi_total_nlpb)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_utau
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_utau.txt"
      call mpi_write_file_1d_dp(mpi_debug_filename, mpi_io_buf_adjust_1d, mpi_total_nlpb)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_vtau
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_vtau.txt"
      call mpi_write_file_1d_dp(mpi_debug_filename, mpi_io_buf_adjust_1d, mpi_total_nlpb)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_qns
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_qns.txt"
      call mpi_write_file_1d_dp(mpi_debug_filename, mpi_io_buf_adjust_1d, mpi_total_nlpb)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_qsr
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_qsr.txt"
      call mpi_write_file_1d_dp(mpi_debug_filename, mpi_io_buf_adjust_1d, mpi_total_nlpb)
   END SUBROUTINE

   ! write files
   SUBROUTINE mpi_debug_write_empmr
      write (mpi_debug_filename, '(I4.4,A1,I8.8,A)') mpi_rank, "_", myIterm, "_empmr.txt"
      call mpi_write_file_1d_dp(mpi_debug_filename, mpi_io_buf_adjust_1d, mpi_total_nlpb)
   END SUBROUTINE

   ! output mpi_nlpb_block information
   SUBROUTINE mpi_debug_output_block_info
      if (mpi_rank == mpi_comp_procs) then
         write (mpi_debug_filename, '(I4.4, A)') mpi_rank, "_mpi_nlpb_block_1d_all.txt"
         call mpi_write_file_1d_integer(mpi_debug_filename, mpi_nlpb_block_1d_all, mpi_comp_procs)

         write (mpi_debug_filename, '(I4.4, A)') mpi_rank, "_mpi_nlpb_block_1d_displs_all.txt"
         call mpi_write_file_1d_integer(mpi_debug_filename, mpi_nlpb_block_1d_displs_all, mpi_comp_procs)

         write (mpi_debug_filename, '(I4.4, A)') mpi_rank, "_mpi_nlpb_block_starts_1d_all.txt"
         call mpi_write_file_1d_integer(mpi_debug_filename, mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_1d_sum)

         write (mpi_debug_filename, '(I4.4, A)') mpi_rank, "_mpi_nlpb_block_counts_1d_all.txt"
         call mpi_write_file_1d_integer(mpi_debug_filename, mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_sum)
      end if
    END SUBROUTINE

   ! write data of 1 dimension to file
   SUBROUTINE mpi_write_file_1d_dp(file_name, var, var_size)
      character(lc), intent(in) :: file_name
      real(dp), intent(in), dimension(:) ::var
      integer, intent(in) :: var_size
      integer :: i
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, var_size
         WRITE (99, *) var(i)
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_write_file_1d_dp

   ! write data of 2 dimension to file
   SUBROUTINE mpi_write_file_2d_dp(file_name, var, var_size1, var_size2)
      character(lc), intent(in) :: file_name
      real(dp), intent(in), dimension(:, :) ::var
      integer, intent(in) :: var_size1, var_size2
      integer :: i, j
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, var_size2
         DO j = 1, var_size1
            WRITE (99, *) var(j, i)
         END DO
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_write_file_2d_dp

   ! write data of 2 dimension to file
   SUBROUTINE mpi_write_file_2d_i2(file_name, var, var_size1, var_size2)
      character(lc), intent(in) :: file_name
      integer(i2), intent(in), dimension(:, :) ::var
      integer, intent(in) :: var_size1, var_size2
      integer :: i, j
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, var_size2
         DO j = 1, var_size1
            WRITE (99, *) var(j, i)
         END DO
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_write_file_2d_i2

   ! write data of 3 dimension to file
   SUBROUTINE mpi_write_file_3d_dp(file_name, var, var_size1, var_size2, var_size3)
      character(lc), intent(in) :: file_name
      real(dp), intent(in), dimension(:, :,:) ::var
      integer, intent(in) :: var_size1, var_size2, var_size3
      integer :: i, j, k
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO k = 1, var_size3
        DO i = 1, var_size2
           DO j = 1, var_size1
              WRITE (99, *) var(j, i, k)
           END DO
        END DO
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_write_file_3d_dp

   ! write integer data of 1 dimension to file
   SUBROUTINE mpi_write_file_1d_integer(file_name, var, var_size)
      character(lc), intent(in) :: file_name
      integer, intent(in), dimension(:) ::var
      integer, intent(in) :: var_size
      integer :: i
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, var_size
         WRITE (99, *) var(i)
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_write_file_1d_integer

   ! test indexes in mpi initialization
   SUBROUTINE mpi_check_indexes(i,j)
     INTEGER,INTENT(INOUT) :: i,j
     DO i = 1,mpi_graph_indegree
       DO j = i + 1,mpi_graph_indegree
         IF(mpi_graph_sources(i)>= mpi_graph_sources(j)) THEN
           WRITE(*,*)"error in mpi_graph_source(",i,"),value is:",mpi_graph_sources(i),&
           ",bigger than mpi_graph_source(",j,"),value is:",mpi_graph_sources(j),&
           ",rankID:",mpi_rank
           STOP
         END IF
       END DO
     END DO

     DO i = 1,mpi_graph_outdegree
       DO j = i + 1,mpi_graph_outdegree
         IF(mpi_graph_dests(i)>= mpi_graph_dests(j)) THEN
           WRITE(*,*)"error in mpi_graph_dests(",i,"),value is:",mpi_graph_dests(i),&
           ",bigger than mpi_graph_dests(",j,"),value is:",mpi_graph_dests(j),&
           ",rankID:",mpi_rank
           STOP
         END IF
       END DO
     END DO

     IF (nlpb - loc_nlpb .ne. mpi_recv_indexes_1d_counts(mpi_graph_indegree) + &
        mpi_recv_indexes_1d_displs(mpi_graph_indegree)) THEN
        WRITE(*,*) "mpi_rank:", mpi_rank, ", nlpb - loc_nlpb:", nlpb - loc_nlpb, "not equal to received size:", &
           mpi_recv_indexes_1d_counts(mpi_graph_indegree) + mpi_recv_indexes_1d_displs(mpi_graph_indegree)
        STOP
     END IF
     DO i = 1,mpi_send_indexes_block_1d_sum
        IF (mpi_send_indexes_1d(i)< 1) THEN
           WRITE(*,*) "error in mpi_send_indexes_1d(", i, "),value:", mpi_send_indexes_1d(i), ",rankID:", mpi_rank
        !  STOP
        END IF
     END DO

     DO i = 1, mpi_send_indexes_block_1d_sum
        IF (mpi_send_indexes_block_1d_starts(i)< 1) THEN
           WRITE(*,*) "error in mpi_send_indexes_block_1d_starts(", i, "),value:", &
              mpi_send_indexes_block_1d_starts(i), ",rankID:", mpi_rank
           WRITE(*,*) "error in mpi_send_indexes_block_1d_starts(", i - 1, "),value:", &
              mpi_send_indexes_block_1d_starts(i - 1), ",rankID:", mpi_rank
        !  STOP
        END IF
        IF (mpi_send_indexes_block_1d_counts(i)< 1) THEN
           WRITE(*,*) "error in mpi_send_indexes_block_1d_counts(", i, "),value:", &
              mpi_send_indexes_block_1d_counts(i), ",rankID:", mpi_rank
           WRITE(*,*) "error in mpi_send_indexes_block_1d_counts(", i - 1, "),value:", &
              mpi_send_indexes_block_1d_counts(i - 1), ",rankID:", mpi_rank
        !  STOP
        END IF
     END DO
   END SUBROUTINE

   ! test indexes in mpi initialization
   SUBROUTINE mpi_check_send_indexes(current_p,i,k)
     INTEGER,INTENT(IN) :: i, k,current_p

     IF(current_p + mpi_send_indexes_block_1d_counts(mpi_send_indexes_block_1d_all_displs(i) + k) - 1 &
        > mpi_send_indexes_1d_counts_sum * nk) THEN
        WRITE(*,*) "mpirankID:",mpi_rank,"access over mpi_sendbuf boundary,the position:", &
           current_p + mpi_send_indexes_block_1d_counts(mpi_send_indexes_block_1d_all_displs(i) + k) - 1
        WRITE(*,*) "the maximum boundary of mpi_sendbuf:",mpi_send_indexes_1d_counts_sum * nk
     END IF
     IF(mpi_send_indexes_block_1d_starts(mpi_send_indexes_block_1d_all_displs(i) + k) + &
        mpi_send_indexes_block_1d_counts(mpi_send_indexes_block_1d_all_displs(i) + k) - 1 > nlpb) THEN
        WRITE(*,*) "mpirankID:",mpi_rank,"access over mpi_source_buf boundary,the position:", &
           mpi_send_indexes_block_1d_starts(mpi_send_indexes_block_1d_all_displs(i) + k) + &
           mpi_send_indexes_block_1d_counts(mpi_send_indexes_block_1d_all_displs(i) + k) - 1
        WRITE(*,*) "the maximum boundary of mpi_source_buf:",nlpb
     END IF
     IF(mpi_send_indexes_block_1d_starts(mpi_send_indexes_block_1d_all_displs(i) + k) < 1) THEN
       WRITE(*,*) "mpirankID:", mpi_rank, "mpi_send_indexes_block_1d_starts:", &
          mpi_send_indexes_block_1d_starts(mpi_send_indexes_block_1d_all_displs(i) + k), &
          ",mpi_send_indexes_block_1d_all_displs:", mpi_send_indexes_block_1d_all_displs(i), &
          ",i:", i, ",k:", k
     END IF
     IF(mpi_send_indexes_block_1d_starts(mpi_send_indexes_block_1d_all_displs(i) + k)+ &
        mpi_send_indexes_block_1d_counts(mpi_send_indexes_block_1d_all_displs(i) + k) - 1 < 1) THEN
     WRITE(*,*) "mpirankID:", mpi_rank, "mpi_send_indexes_block_1d_counts:", &
        mpi_send_indexes_block_1d_counts(mpi_send_indexes_block_1d_all_displs(i) + k), &
        ",mpi_send_indexes_block_1d_all_displs:", mpi_send_indexes_block_1d_all_displs(i), &
        ",i:", i, ",k:", k
     END IF
   END SUBROUTINE
END MODULE
