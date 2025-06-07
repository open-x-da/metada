MODULE mod_mpi_interfaces
   USE mod_csp_basic
   USE mod_mpi_reallocate
   USE mod_mpi_variables
   use mod_mpi_test_variables
   use mod_io_netcdf
!!YY: seaice module
   use mitice_parameters
   use mitice_vars
# if defined (OPENACC)
   use openacc
# endif
   IMPLICIT NONE

   ! exchanging variables in computing processes in blocking communication or non - blocking communication
   INTERFACE mpi_data_exchange !数据交换！！！！！！！！！！！！！！！！！！！
!      MODULE PROCEDURE mpi_i2_1d_block_exchange
      MODULE PROCEDURE mpi_i2_2d_block_exchange
      MODULE PROCEDURE mpi_wp_1d_block_exchange
      MODULE PROCEDURE mpi_wp_2d_block_exchange
      MODULE PROCEDURE mpi_wp_3d_block_exchange
      MODULE PROCEDURE mpi_wp_1d_nonblock_exchange
      MODULE PROCEDURE mpi_wp_2d_nonblock_exchange
      MODULE PROCEDURE mpi_wp_3d_nonblock_exchange
   END INTERFACE

   ! check and wait until non - blocking communication is finished
   INTERFACE mpi_data_check_and_wait
      MODULE PROCEDURE mpi_dp_1d_check_and_wait
      MODULE PROCEDURE mpi_dp_2d_check_and_wait
      MODULE PROCEDURE mpi_dp_3d_check_and_wait
   END INTERFACE

   ! used for netcdf input reading and exchanging
   INTERFACE mpi_netcdf_read_exchange
      MODULE PROCEDURE mpi_netcdf_read_exchange_1d_wp
      MODULE PROCEDURE mpi_netcdf_read_exchange_1d_integer
      MODULE PROCEDURE mpi_netcdf_read_exchange_1d_i2
      MODULE PROCEDURE mpi_netcdf_read_exchange_2d_integer
      MODULE PROCEDURE mpi_netcdf_read_exchange_2d_i2
      MODULE PROCEDURE mpi_netcdf_read_exchange_2d_wp
   END INTERFACE

   ! used for the root computing process to send data from input netcdf files
   INTERFACE mpi_netcdf_root_exchange
      MODULE PROCEDURE mpi_netcdf_root_exchange_1d_wp
      MODULE PROCEDURE mpi_netcdf_root_exchange_1d_integer
      MODULE PROCEDURE mpi_netcdf_root_exchange_1d_i2
      MODULE PROCEDURE mpi_netcdf_root_exchange_2d_integer
      MODULE PROCEDURE mpi_netcdf_root_exchange_2d_i2
      MODULE PROCEDURE mpi_netcdf_root_exchange_2d_wp
   END INTERFACE

   ! used for the non - root computing processes to receive data sent from root computing process
   INTERFACE mpi_netcdf_nonroot_exchange
      MODULE PROCEDURE mpi_netcdf_nonroot_exchange_1d_wp
      MODULE PROCEDURE mpi_netcdf_nonroot_exchange_1d_integer
      MODULE PROCEDURE mpi_netcdf_nonroot_exchange_1d_i2
      MODULE PROCEDURE mpi_netcdf_nonroot_exchange_2d_integer
      MODULE PROCEDURE mpi_netcdf_nonroot_exchange_2d_i2
      MODULE PROCEDURE mpi_netcdf_nonroot_exchange_2d_wp
   END INTERFACE


CONTAINS

  ! deallocate buffer after reading netcdf input files
   SUBROUTINE mpi_netcdf_read_finalize
      IF (mpi_rank == 0) THEN
         DEALLOCATE (mpi_buf_1d_nlpb_wp)
         DEALLOCATE (mpi_buf_1d_nlpb_nlpbz_wp)
         DEALLOCATE (mpi_buf_1d_nlpb_nlpbz_integer)
         DEALLOCATE (mpi_buf_1d_nlpb_integer)
         DEALLOCATE (mpi_buf_1d_nlpb_i2)
         DEALLOCATE (mpi_adjust_buf_1d_nlpb_shadow_wp)
         DEALLOCATE (mpi_adjust_buf_1d_nlpb_shadow_nlpbz_wp)
         DEALLOCATE (mpi_adjust_buf_1d_nlpb_shadow_nlpbz_integer)
         DEALLOCATE (mpi_adjust_buf_1d_nlpb_shadow_integer)
         DEALLOCATE (mpi_adjust_buf_1d_nlpb_shadow_i2)
         DEALLOCATE (mpi_nlpb_shadow_nlpbz_block_1d_all)
         DEALLOCATE (mpi_nlpb_shadow_nlpbz_block_starts_1d_all)
         DEALLOCATE (mpi_nlpb_shadow_nlpbz_block_counts_1d_all)
         DEALLOCATE (mpi_nlpb_shadow_nlpbz_counts_all)
         DEALLOCATE (mpi_nlpb_shadow_nlpbz_block_1d_displs_all)
         DEALLOCATE (mpi_nlpb_shadow_nlpbz_displs_all)

         DEALLOCATE (mpi_buf_2d_nlpb_ni_integer)
         DEALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_ni_integer)
         DEALLOCATE (mpi_buf_2d_nlpb_nlpbz_ni_integer)
         DEALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_nlpbz_ni_integer)
         DEALLOCATE (mpi_buf_2d_nlpb_nk_i2)
         DEALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_nk_integer)
         DEALLOCATE (mpi_buf_2d_nlpb_nk_wp)
         DEALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_nk_wp)
         DEALLOCATE (mpi_buf_2d_nlpb_nlpbz_nk_i2)
         DEALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_nlpbz_nk_integer)
         DEALLOCATE (mpi_adjust_2d_nlpb_shadow_ni_counts_all)
         DEALLOCATE (mpi_adjust_2d_nlpb_shadow_ni_displs_all)
         DEALLOCATE (mpi_adjust_2d_nlpb_shadow_nk_counts_all)
         DEALLOCATE (mpi_adjust_2d_nlpb_shadow_nk_displs_all)
         DEALLOCATE (mpi_adjust_2d_nlpb_shadow_nlpbz_ni_counts_all)
         DEALLOCATE (mpi_adjust_2d_nlpb_shadow_nlpbz_ni_displs_all)
         DEALLOCATE (mpi_adjust_2d_nlpb_shadow_nlpbz_nk_counts_all)
         DEALLOCATE (mpi_adjust_2d_nlpb_shadow_nlpbz_nk_displs_all)

         if (ln_bdy) then
            DEALLOCATE (mpi_original_idx_bdy_all)
            DEALLOCATE (mpi_bdy_dests_all)
            DEALLOCATE (mpi_bdy_counts_1d_all)
            DEALLOCATE (mpi_bdy_displs_1d_all)
            DEALLOCATE (mpi_original_to_ordered_bdy_all)
            DEALLOCATE (mpi_ordered_idx_bdy_all)

            IF (ALLOCATED(mpi_buf_1d_nlbdy_wp))  DEALLOCATE(mpi_buf_1d_nlbdy_wp)
            IF(ALLOCATED(mpi_adjust_buf_1d_nlbdy_wp)) DEALLOCATE(mpi_adjust_buf_1d_nlbdy_wp)
            IF (ALLOCATED(mpi_buf_1d_nlbdy_integer)) DEALLOCATE(mpi_buf_1d_nlbdy_integer)
            IF (ALLOCATED(mpi_adjust_buf_1d_nlbdy_integer)) DEALLOCATE(mpi_adjust_buf_1d_nlbdy_integer)
         end if

      END IF

      IF (ln_bdy .and. mpi_check_bdy) THEN
        DEALLOCATE(mpi_idx_bdy)
      END IF
   END SUBROUTINE mpi_netcdf_read_finalize

   ! root process broadcast information to non - root processes among computing processes
   SUBROUTINE mpi_netcdf_root_bcast_first
      CHARACTER, ALLOCATABLE, DIMENSION(:) :: mpi_buf_tmp
      INTEGER :: mpi_buf_size_tmp
      INTEGER :: mpi_position_tmp

      mpi_position_tmp = 0
      mpi_buf_size_tmp = sizeof(mpi_total_nlpb)* 7
      ALLOCATE (mpi_buf_tmp(mpi_buf_size_tmp))

      CALL MPI_PACK(mpi_total_nlpb, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(mpi_total_nl, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(mpi_total_nlpbz, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(mpi_total_nlz, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nk, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nkp1, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ni, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)

      CALL MPI_BCAST(mpi_buf_tmp, mpi_buf_size_tmp, MPI_PACKED, mpi_root_comp, mpi_comp_comm, mpi_err)
      DEALLOCATE (mpi_buf_tmp)
   END SUBROUTINE

   ! non - root processes receive information from the root processes
   SUBROUTINE mpi_netcdf_nonroot_bcast_first
      CHARACTER, ALLOCATABLE, DIMENSION(:) :: mpi_buf_tmp
      INTEGER :: mpi_buf_size_tmp
      INTEGER :: mpi_position_tmp

      mpi_position_tmp = 0
      mpi_buf_size_tmp = sizeof(mpi_total_nlpb)* 7
      ALLOCATE (mpi_buf_tmp(mpi_buf_size_tmp))

      CALL MPI_BCAST(mpi_buf_tmp, mpi_buf_size_tmp, MPI_PACKED, mpi_root_comp, mpi_comp_comm, mpi_err)

      CALL MPI_UNPACK(mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_total_nlpb, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_total_nl, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_total_nlpbz, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_total_nlz, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, nk, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, nkp1, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, ni, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      DEALLOCATE (mpi_buf_tmp)

   END SUBROUTINE

   ! root process broadcast information to non - root processes among computing processes at 2nd time
   SUBROUTINE mpi_netcdf_root_bcast_second
      CHARACTER, ALLOCATABLE, DIMENSION(:) :: mpi_buf_tmp
      INTEGER :: mpi_buf_size_tmp
      INTEGER :: mpi_position_tmp

      mpi_position_tmp = 0
      mpi_buf_size_tmp = sizeof(mpi_total_nlpb)* 7
      ALLOCATE (mpi_buf_tmp(mpi_buf_size_tmp))

      CALL MPI_PACK(mpi_total_nlpb, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nl, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(mpi_total_nlpbz, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nlz, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nk, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nkp1, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ni, 1, MPI_INTEGER, mpi_buf_tmp, mpi_buf_size_tmp, mpi_position_tmp, mpi_comp_comm, mpi_err)

      CALL MPI_BCAST(mpi_buf_tmp, mpi_buf_size_tmp, MPI_PACKED, mpi_root_comp, mpi_comp_comm, mpi_err)
      DEALLOCATE (mpi_buf_tmp)

   END SUBROUTINE

   ! initialize grid index and graph communication, and send index to I / O process
   SUBROUTINE mpi_index_graph_init_and_send
     ! initialize grid index in 1 dimension
      CALL mpi_init_1d_index_buf
      CALL mpi_send_index_to_io
      CALL mpi_init_2d_index_buf
      CALL mpi_get_graph_source
      CALL mpi_get_graph_destination
      if (ln_bdy) then
         CALL mpi_init_bdy
      end if
      CALL mpi_init_graph_index
      CALL mpi_init_arrays

      ! ---------------------------------------------------------------------------
      !  VARS IN THIS SUBROUTINE ARE ACC-DECLARED IN MOD_MPI_VARIABLES
      !  BEFORE UPDATE, THEY HAS BEEN ALLOCATED AND INIT-ED
      ! ---------------------------------------------------------------------------

      !$ACC UPDATE DEVICE(mpi_send_indexes_1d, mpi_send_indexes_1d_counts,    &
      !$ACC mpi_send_indexes_1d_displs, mpi_recv_indexes_1d_counts,           &
      !$ACC mpi_recv_indexes_1d_displs, mpi_sendbuf_1d, mpi_sendbuf_2d_nk_f,  &
      !$ACC mpi_sendbuf_2d_nkp1, mpi_sendbuf_2d_nk_int, mpi_global_grids)

   END SUBROUTINE

   ! initialize index and communicator for reading netcdf file and scattering to some processes with boundary
   SUBROUTINE mpi_init_bdy
      CHARACTER(lc) :: filename
      INTEGER, ALLOCATABLE,DIMENSION(:) :: bdy_counts_tmp
      INTEGER, ALLOCATABLE,DIMENSION(:) :: bdy_dests_tmp
!      INTEGER, ALLOCATABLE,DIMENSION(:) :: bdy_processes
      INTEGER, ALLOCATABLE,DIMENSION(:) :: mpi_buf
      INTEGER, DIMENSION(:,:),ALLOCATABLE :: mpi_nlbdy_shadow_all ! all number of nlbdy and shadow, only for process 0
      INTEGER :: i,pos
      IF (mpi_rank == 0) THEN
        ALLOCATE(mpi_original_idx_bdy_all(mpi_total_nlbdy))
        filename = "nmefc_macom_csp_grid_d01.nc"
        CALL netcdf_read(trim(filename), 'idx_bdy', mpi_original_idx_bdy_all,d1 = mpi_total_nlbdy)

        filename = "nmefc_macom_par_dec.nc"
        ALLOCATE(mpi_nlbdy_shadow_all(2,mpi_comp_procs))
        CALL netcdf_read(trim(filename), 'nlbdy_layer', mpi_nlbdy_shadow_all,d1 = 2,d2 = mpi_comp_procs)

        ALLOCATE(mpi_buf(mpi_comp_procs))
        mpi_buf(:)= mpi_nlbdy_shadow_all(1,:)
        ! scatter number local boundary on each process
        CALL MPI_SCATTER(mpi_buf,1,MPI_INTEGER,loc_nlbdy,1,MPI_INTEGER,0,mpi_comp_comm,mpi_err)

        mpi_buf(:)= mpi_nlbdy_shadow_all(2,:)
        ! scatter number local boundary and their shadow on each process
        CALL MPI_SCATTER(mpi_buf,1,MPI_INTEGER,nlbdy,1,MPI_INTEGER,0,mpi_comp_comm,mpi_err)
        DEALLOCATE(mpi_buf)

        mpi_total_nlbdy_shadow = 0
        mpi_bdy_outdegree = 0
        DO i = 1,mpi_comp_procs
          IF(mpi_nlbdy_shadow_all(2,i)> 0) THEN
            mpi_total_nlbdy_shadow = mpi_total_nlbdy_shadow + mpi_nlbdy_shadow_all(2,i)
            mpi_bdy_outdegree = mpi_bdy_outdegree + 1
          END IF
        END DO

        pos = 0
        ALLOCATE(mpi_bdy_dests_all(mpi_bdy_outdegree))
        ALLOCATE(mpi_bdy_counts_1d_all(mpi_bdy_outdegree))
        DO i = 1,mpi_comp_procs
          IF(mpi_nlbdy_shadow_all(2,i)> 0) THEN
            pos = pos + 1
            mpi_bdy_dests_all(pos)= i - 1
            mpi_bdy_counts_1d_all(pos)= mpi_nlbdy_shadow_all(2,i)
          END IF
        END DO

        DEALLOCATE(mpi_nlbdy_shadow_all)

        IF(mpi_bdy_dests_all(1) .ne. 0) then
          ! add process 0 into mpi_bdy_comm, process 0 only is in charge of reading data
          mpi_bdy_outdegree = mpi_bdy_outdegree + 1
          CALL reallocate(mpi_bdy_dests_all,mpi_bdy_outdegree)
          CALL reallocate(mpi_bdy_counts_1d_all,mpi_bdy_outdegree)
          mpi_bdy_dests_all(2:mpi_bdy_outdegree)= mpi_bdy_dests_all(1:mpi_bdy_outdegree - 1)
          mpi_bdy_dests_all(1)= 0
          mpi_bdy_counts_1d_all(2:mpi_bdy_outdegree)= mpi_bdy_counts_1d_all(1:mpi_bdy_outdegree - 1)
          mpi_bdy_counts_1d_all(1)= 0
        END IF

        ALLOCATE(mpi_bdy_displs_1d_all(mpi_bdy_outdegree))
        mpi_bdy_displs_1d_all(1)= 0
        DO i = 2,mpi_bdy_outdegree
          mpi_bdy_displs_1d_all(i)= mpi_bdy_displs_1d_all(i - 1)+ mpi_bdy_counts_1d_all(i - 1)
        END DO

        ALLOCATE(mpi_original_to_ordered_bdy_all(mpi_total_nlbdy_shadow))
        ALLOCATE(mpi_ordered_idx_bdy_all(mpi_total_nlbdy_shadow))
!        ALLOCATE(bdy_processes(mpi_total_nlbdy_shadow))
        CALL netcdf_read(trim(filename), 'idx_bdy_ordered', mpi_ordered_idx_bdy_all,d1 = mpi_total_nlbdy_shadow)
!        CALL netcdf_read(trim(filename), 'idx_bdy_process', bdy_processes,d1 = mpi_total_nlbdy_shadow)

        DO i = 1,mpi_total_nlbdy_shadow
          pos = mpi_findloc(mpi_original_idx_bdy_all,mpi_ordered_idx_bdy_all(i))
          IF(pos > 0) THEN
            mpi_original_to_ordered_bdy_all(i)= pos
          ELSE
            WRITE(*,*) "there is no match between idx_bdy in nmefc_macom_csp_grid_d01.nc &
               and idx_bdy_ordered in nmefc_macom_par_dec.nc"
            STOP
          END IF
        END DO

        CALL MPI_BCAST(mpi_bdy_outdegree,1,MPI_INTEGER,0,mpi_comp_comm,mpi_err)

        CALL MPI_BCAST(mpi_bdy_dests_all,mpi_bdy_outdegree,MPI_INTEGER,0,mpi_comp_comm,mpi_err)

      ELSE

        ! scatter number local boundary on each process
        CALL MPI_SCATTER(mpi_buf,1,MPI_INTEGER,loc_nlbdy,1,MPI_INTEGER,0,mpi_comp_comm,mpi_err)

        ! scatter number local boundary and their shadow on each process
        CALL MPI_SCATTER(mpi_buf,1,MPI_INTEGER,nlbdy,1,MPI_INTEGER,0,mpi_comp_comm,mpi_err)

        CALL MPI_BCAST(mpi_bdy_outdegree,1,MPI_INTEGER,0,mpi_comp_comm,mpi_err)

        ALLOCATE(mpi_bdy_dests_all(mpi_bdy_outdegree))
        CALL MPI_BCAST(mpi_bdy_dests_all,mpi_bdy_outdegree,MPI_INTEGER,0,mpi_comp_comm,mpi_err)


      END IF

      ! check whether rank in MPI_COMM_WORLD is same with mpi_comp_comm, it should be same,
      ! otherwise there will be error in creating mpi_bdy_comm
      IF(mpi_rank /= mpi_comp_rank) THEN
        WRITE(*,*) "mpi_rank =",mpi_rank," is not equal to mpi_comp_rank =",mpi_comp_rank
        STOP
      END IF

      ! create group and communicator for sending data in boundary
      CALL MPI_GROUP_INCL(mpi_group_all, mpi_bdy_outdegree, mpi_bdy_dests_all, mpi_group_bdy, mpi_err)
      CALL MPI_COMM_CREATE(mpi_comp_comm, mpi_group_bdy, mpi_bdy_comm, mpi_err)

      DO i = 1,mpi_bdy_outdegree
        IF(mpi_rank == mpi_bdy_dests_all(i)) THEN
          mpi_check_bdy =.TRUE.
          EXIT
        END IF
      END DO

      IF(mpi_check_bdy) THEN
        IF(mpi_rank == 0) THEN
          ! it is possible that there is no boundary for process 0, so nlbdy = 0
          ALLOCATE(mpi_idx_bdy(nlbdy))
!
          CALL MPI_SCATTERV(mpi_ordered_idx_bdy_all,mpi_bdy_counts_1d_all,mpi_bdy_displs_1d_all,MPI_INTEGER, &
                            mpi_idx_bdy,nlbdy,MPI_INTEGER,0,mpi_bdy_comm,mpi_err)

          IF(nlbdy > 0)THEN
             CALL mpi_adjust_sub_indexes_ori(mpi_idx_bdy, nlbdy)
          END IF
        ELSE
          ALLOCATE(mpi_idx_bdy(nlbdy))
          CALL MPI_SCATTERV(mpi_ordered_idx_bdy_all,mpi_bdy_counts_1d_all,mpi_bdy_displs_1d_all,MPI_INTEGER, &
                            mpi_idx_bdy,nlbdy,MPI_INTEGER,0,mpi_bdy_comm,mpi_err)
          CALL mpi_adjust_sub_indexes_ori(mpi_idx_bdy, nlbdy)
        END IF
        !$acc update device(mpi_idx_bdy)
      END IF

   END SUBROUTINE

   ! send data from root of mpi_comp_comm to the root of mpi_io_comm
   SUBROUTINE mpi_send_index_to_io
     IF (mpi_rank == 0) THEN
       CALL MPI_SEND(mpi_nlpb_counts_all, mpi_comp_procs, MPI_INTEGER, mpi_comp_procs, mpi_tag1, MPI_COMM_WORLD, mpi_err)
       CALL MPI_SEND(mpi_nlpb_block_1d_all, mpi_comp_procs, MPI_INTEGER, mpi_comp_procs, mpi_tag1, MPI_COMM_WORLD, mpi_err)
       CALL MPI_SEND(mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_1d_sum, MPI_INTEGER, mpi_comp_procs, mpi_tag1, &
                      MPI_COMM_WORLD, mpi_err)
       CALL MPI_SEND(mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_sum, MPI_INTEGER, mpi_comp_procs, mpi_tag1, &
                      MPI_COMM_WORLD, mpi_err)
     END IF

     END SUBROUTINE

   ! initialize data index for data exchange in 1 dimension
   SUBROUTINE mpi_init_1d_index_buf
      INTEGER :: mpi_i
      INTEGER :: mpi_j
      ! only for root computing process, buffer for all global local, shadow and z grid index
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mpi_global_grids_all
!      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_starts ! buffer for temporary data exchange
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mpi_process_index_all ! only for root computing process, save all mpi_process_index
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mpi_zone_index_all ! only for root computing process, save all mpi_zone_index
      ! only for root computing process,  real data exclude null value from mpi_global_grids_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_global_grids_all_sendbuf ! buffer for temporary data exchange
      ! only for root computing process,  real data exclude null value from mpi_process_index_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_process_index_all_sendbuf
      ! only for root computing process,  bufffer for sending data of mpi_zone_index_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_zone_index_all_sendbuf
      ! only for root computing process, array for size of nlpb, nlpb + shadow, nlpb + shadow + z grids
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mpi_nlpb_shadow_nlpbz_all ! only for root computing process
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_starts_tmp ! temporay buffer for mpi_nlpb_shadow_nlpbz_starts
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_counts_tmp ! temporay buffer for mpi_nlpb_shadow_nlpbz_counts_all
      ! for all computing processes, size of nlpb, nlpb + shadow, nlpb + shadow + z grids
      INTEGER :: mpi_nlpb_shadow_nlpbz(3)

      CHARACTER(lc) :: filename

      IF (mpi_rank == 0) THEN
        filename = "nmefc_macom_par_dec.nc"
        ALLOCATE(mpi_nlpb_shadow_nlpbz_all(3,mpi_comp_procs))

        CALL netcdf_read(trim(filename), 'nlloc', mpi_max_nlpb_shadow_nlpbz)
        CALL netcdf_read(trim(filename), 'loc_nlpb', mpi_nlpb_shadow_nlpbz_all, 3, mpi_comp_procs)

        ! scatter size of nlpb, nlpb + shadow, nlpb + shadow + z grids
        CALL MPI_SCATTER(mpi_nlpb_shadow_nlpbz_all,3,MPI_INTEGER,mpi_nlpb_shadow_nlpbz,3,MPI_INTEGER,0,mpi_comp_comm,mpi_err)

        ! index for scattering data read from input netcdf files
        ALLOCATE(mpi_nlpb_counts_all(mpi_comp_procs))
        ALLOCATE(mpi_nlpb_displs_all(mpi_comp_procs))
        ALLOCATE(mpi_nlpb_shadow_counts_all(mpi_comp_procs))
        ALLOCATE(mpi_nlpb_shadow_displs_all(mpi_comp_procs))
        ALLOCATE(mpi_nlpb_shadow_nlpbz_counts_all(mpi_comp_procs))
        ALLOCATE(mpi_nlpb_shadow_nlpbz_displs_all(mpi_comp_procs))
        mpi_nlpb_displs_all(1) = 0
        mpi_nlpb_shadow_displs_all(1) = 0
        mpi_nlpb_shadow_nlpbz_displs_all(1) = 0
        mpi_nlpb_counts_all(1) = mpi_nlpb_shadow_nlpbz_all(1,1)
        mpi_nlpb_shadow_counts_all(1) = mpi_nlpb_shadow_nlpbz_all(2,1)
        mpi_nlpb_shadow_nlpbz_counts_all(1) = mpi_nlpb_shadow_nlpbz_all(3,1)
        mpi_total_nlpb_shadow = mpi_nlpb_shadow_nlpbz_all(2,1)
        mpi_total_nlpb_shadow_nlpbz = mpi_nlpb_shadow_nlpbz_all(3,1)
        DO mpi_i = 2, mpi_comp_procs, 1
          mpi_nlpb_counts_all(mpi_i) = mpi_nlpb_shadow_nlpbz_all(1,mpi_i)
          mpi_nlpb_displs_all(mpi_i) = mpi_nlpb_displs_all(mpi_i - 1) + mpi_nlpb_counts_all(mpi_i - 1)
          mpi_nlpb_shadow_counts_all(mpi_i) = mpi_nlpb_shadow_nlpbz_all(2,mpi_i)
          mpi_nlpb_shadow_displs_all(mpi_i) = mpi_nlpb_shadow_displs_all(mpi_i - 1) + &
                                                    mpi_nlpb_shadow_counts_all(mpi_i - 1)
          mpi_nlpb_shadow_nlpbz_counts_all(mpi_i) = mpi_nlpb_shadow_nlpbz_all(3,mpi_i)
          mpi_nlpb_shadow_nlpbz_displs_all(mpi_i) = mpi_nlpb_shadow_nlpbz_displs_all(mpi_i - 1) + &
                                                    mpi_nlpb_shadow_nlpbz_counts_all(mpi_i - 1)
          mpi_total_nlpb_shadow = mpi_total_nlpb_shadow + mpi_nlpb_shadow_nlpbz_all(2,mpi_i)
          mpi_total_nlpb_shadow_nlpbz = mpi_total_nlpb_shadow_nlpbz + mpi_nlpb_shadow_nlpbz_all(3,mpi_i)
        END DO
        mpi_nlpb_counts_sum = mpi_nlpb_displs_all(mpi_comp_procs) + mpi_nlpb_counts_all(mpi_comp_procs)
        mpi_nlpb_shadow_counts_sum = mpi_nlpb_shadow_displs_all(mpi_comp_procs) + &
                                           mpi_nlpb_shadow_counts_all(mpi_comp_procs)
        mpi_nlpb_shadow_nlpbz_counts_sum = mpi_nlpb_shadow_nlpbz_displs_all(mpi_comp_procs) + &
                                           mpi_nlpb_shadow_nlpbz_counts_all(mpi_comp_procs)


        ! scatter global grid index to each process
        ALLOCATE(mpi_global_grids_all_sendbuf(mpi_nlpb_shadow_nlpbz_counts_sum))
        ALLOCATE(mpi_global_grids_all(mpi_max_nlpb_shadow_nlpbz,mpi_comp_procs))
        ALLOCATE(mpi_global_grids(mpi_nlpb_shadow_nlpbz(3)))
        CALL netcdf_read(trim(filename), 'loc_idx1', mpi_global_grids_all, mpi_max_nlpb_shadow_nlpbz, mpi_comp_procs)
        CALL mpi_integer_2d_prepare_sendbuf_for_par_dec(mpi_global_grids_all,mpi_global_grids_all_sendbuf,mpi_comp_procs, &
                                                        mpi_nlpb_shadow_nlpbz_counts_all)
        CALL MPI_SCATTERV(mpi_global_grids_all_sendbuf, mpi_nlpb_shadow_nlpbz_counts_all, &
                          mpi_nlpb_shadow_nlpbz_displs_all, MPI_INTEGER, mpi_global_grids, &
                          mpi_nlpb_shadow_nlpbz(3), MPI_INTEGER, 0, mpi_comp_comm,mpi_err)

        ! scatter process index to each process
        ALLOCATE(mpi_process_index_all(mpi_max_nlpb_shadow_nlpbz,mpi_comp_procs))
        ALLOCATE(mpi_process_index_all_sendbuf(mpi_nlpb_shadow_nlpbz_counts_sum))
        ALLOCATE(mpi_process_index(mpi_nlpb_shadow_nlpbz(3)))
        CALL netcdf_read(trim(filename), 'loc_idx2', mpi_process_index_all, mpi_max_nlpb_shadow_nlpbz, mpi_comp_procs)
        CALL mpi_integer_2d_prepare_sendbuf_for_par_dec(mpi_process_index_all,mpi_process_index_all_sendbuf,mpi_comp_procs, &
                                                        mpi_nlpb_shadow_nlpbz_counts_all)
        CALL MPI_SCATTERV(mpi_process_index_all_sendbuf, mpi_nlpb_shadow_nlpbz_counts_all, &
                          mpi_nlpb_shadow_nlpbz_displs_all, MPI_INTEGER, &
                          mpi_process_index, mpi_nlpb_shadow_nlpbz(3), &
                          MPI_INTEGER, 0, mpi_comp_comm,mpi_err)

        ! scatter zone index to each process
        ALLOCATE(mpi_zone_index_all(mpi_max_nlpb_shadow_nlpbz,mpi_comp_procs))
        ALLOCATE(mpi_zone_index_all_sendbuf(mpi_nlpb_shadow_nlpbz_counts_sum))
        ALLOCATE(mpi_zone_index(mpi_nlpb_shadow_nlpbz(3)))
        CALL netcdf_read(trim(filename), 'loc_idx4', mpi_zone_index_all, mpi_max_nlpb_shadow_nlpbz, mpi_comp_procs)
        CALL mpi_integer_2d_prepare_sendbuf_for_par_dec(mpi_zone_index_all,mpi_zone_index_all_sendbuf,mpi_comp_procs, &
                                                        mpi_nlpb_shadow_nlpbz_counts_all)
        CALL MPI_SCATTERV(mpi_zone_index_all_sendbuf, mpi_nlpb_shadow_nlpbz_counts_all, &
                          mpi_nlpb_shadow_nlpbz_displs_all, MPI_INTEGER, &
                          mpi_zone_index, mpi_nlpb_shadow_nlpbz(3), &
                          MPI_INTEGER, 0, mpi_comp_comm,mpi_err)

      ELSE
        ! receive size of nlpb, nlpb + shadow, nlpb + shadow + z grids
        CALL MPI_SCATTER(mpi_nlpb_shadow_nlpbz_all,3,MPI_INTEGER,mpi_nlpb_shadow_nlpbz,3,MPI_INTEGER,0,mpi_comp_comm,mpi_err)

        ! receive global grid index from root process
        ALLOCATE(mpi_global_grids(mpi_nlpb_shadow_nlpbz(3)))
        CALL MPI_SCATTERV(mpi_global_grids_all_sendbuf, mpi_nlpb_shadow_nlpbz_counts_all, &
                          mpi_nlpb_shadow_nlpbz_displs_all, MPI_INTEGER, &
                          mpi_global_grids, mpi_nlpb_shadow_nlpbz(3), &
                          MPI_INTEGER, 0, mpi_comp_comm,mpi_err)

        ! receive process index from root process
        ALLOCATE(mpi_process_index(mpi_nlpb_shadow_nlpbz(3)))
        CALL MPI_SCATTERV(mpi_process_index_all_sendbuf, mpi_nlpb_shadow_nlpbz_counts_all, &
                          mpi_nlpb_shadow_nlpbz_displs_all, MPI_INTEGER, &
                          mpi_process_index, mpi_nlpb_shadow_nlpbz(3), &
                          MPI_INTEGER, 0, mpi_comp_comm,mpi_err)

        ! receive zone index to each process
        ALLOCATE(mpi_zone_index(mpi_nlpb_shadow_nlpbz(3)))
        CALL MPI_SCATTERV(mpi_zone_index_all_sendbuf, mpi_nlpb_shadow_nlpbz_counts_all, &
                          mpi_nlpb_shadow_nlpbz_displs_all, MPI_INTEGER, &
                          mpi_zone_index, mpi_nlpb_shadow_nlpbz(3), &
                          MPI_INTEGER, 0, mpi_comp_comm,mpi_err)

      ENDIF

      loc_nlpb = mpi_nlpb_shadow_nlpbz(1)
      ! for nlpb = mpi_nlpb_shadow_nlpbz(2), is used to be consistent to old codes, so the name is nlpb_shadow is better
      nlpb = mpi_nlpb_shadow_nlpbz(2)
      nlpbz = mpi_nlpb_shadow_nlpbz(3)

      ! update nl to nlpb
      nl = nlpb
      ! update nlz to nlpbz
      nlz = nlpbz

      ! check process value in nmefc_macom_par_dec.nc
      DO mpi_i = 1,nlpb
        IF(mpi_process_index(mpi_i) < 0) THEN
          WRITE(*,*) "process ID in process:", mpi_rank, " is wrong, position is:", &
             mpi_i, ", value is:", mpi_process_index(mpi_i)
          STOP
        END IF
      END DO


      ! get starting global index, size of contiguous blocks in local grids(nlpb)
      ALLOCATE (mpi_nlpb_shadow_nlpbz_starts_tmp(mpi_nlpb_shadow_nlpbz(3)))
      ALLOCATE (mpi_nlpb_shadow_nlpbz_counts_tmp(mpi_nlpb_shadow_nlpbz(3)))
      mpi_i = 1
      mpi_j = 2
      mpi_nlpb_block_1d = 1
      mpi_nlpb_shadow_nlpbz_starts_tmp(1) = mpi_global_grids(1)
      DO WHILE (mpi_i < mpi_nlpb_shadow_nlpbz(1))
        IF (mpi_global_grids(mpi_i) + 1 /= mpi_global_grids(mpi_j)) THEN
           ! get the size of each contiguous block
           mpi_nlpb_shadow_nlpbz_counts_tmp(mpi_nlpb_block_1d) = mpi_global_grids(mpi_i) - &
                                                                 mpi_nlpb_shadow_nlpbz_starts_tmp(mpi_nlpb_block_1d) + 1
           ! get number of blocks
           mpi_nlpb_block_1d = mpi_nlpb_block_1d + 1
           ! get starting index number of each block
           mpi_nlpb_shadow_nlpbz_starts_tmp(mpi_nlpb_block_1d) = mpi_global_grids(mpi_j)
        END IF
        mpi_i = mpi_i + 1
        mpi_j = mpi_j + 1
      END DO
      mpi_nlpb_shadow_nlpbz_counts_tmp(mpi_nlpb_block_1d) = mpi_global_grids(mpi_i) - &
        mpi_nlpb_shadow_nlpbz_starts_tmp(mpi_nlpb_block_1d) + 1
      ! copy starting index number of each block, and size of each contiguous block
      ALLOCATE (mpi_nlpb_block_starts_1d(mpi_nlpb_block_1d))
      ALLOCATE (mpi_nlpb_block_counts_1d(mpi_nlpb_block_1d))
      mpi_nlpb_block_starts_1d(1:mpi_nlpb_block_1d) = mpi_nlpb_shadow_nlpbz_starts_tmp(1:mpi_nlpb_block_1d)
      mpi_nlpb_block_counts_1d(1:mpi_nlpb_block_1d) = mpi_nlpb_shadow_nlpbz_counts_tmp(1:mpi_nlpb_block_1d)

      ! add starting global index, size of contiguous blocks in shadow grids
      mpi_i = mpi_nlpb_shadow_nlpbz(1) + 1
      mpi_j = mpi_nlpb_shadow_nlpbz(1) + 2
      mpi_nlpb_shadow_block_1d = mpi_nlpb_block_1d + 1
      mpi_nlpb_shadow_nlpbz_starts_tmp(mpi_nlpb_shadow_block_1d) = mpi_global_grids(mpi_nlpb_shadow_nlpbz(1) + 1)
      DO WHILE (mpi_i < mpi_nlpb_shadow_nlpbz(2))
        IF (mpi_global_grids(mpi_i) + 1 /= mpi_global_grids(mpi_j)) THEN
           ! get the size of each contiguous block
           mpi_nlpb_shadow_nlpbz_counts_tmp(mpi_nlpb_shadow_block_1d) = mpi_global_grids(mpi_i) - &
                                                                 mpi_nlpb_shadow_nlpbz_starts_tmp(mpi_nlpb_shadow_block_1d) + 1
           ! get number of blocks
           mpi_nlpb_shadow_block_1d = mpi_nlpb_shadow_block_1d + 1
           ! get starting index number of each block
           mpi_nlpb_shadow_nlpbz_starts_tmp(mpi_nlpb_shadow_block_1d) = mpi_global_grids(mpi_j)
        END IF
        mpi_i = mpi_i + 1
        mpi_j = mpi_j + 1
      END DO
      ! the previous cicle doesn't count the last block size and number
      mpi_nlpb_shadow_nlpbz_counts_tmp(mpi_nlpb_shadow_block_1d) = mpi_global_grids(mpi_i) - &
        mpi_nlpb_shadow_nlpbz_starts_tmp(mpi_nlpb_shadow_block_1d) + 1
      ! copy starting index number of each block, and size of each contiguous block
      ALLOCATE (mpi_nlpb_shadow_block_starts_1d(mpi_nlpb_shadow_block_1d))
      ALLOCATE (mpi_nlpb_shadow_block_counts_1d(mpi_nlpb_shadow_block_1d))
      mpi_nlpb_shadow_block_starts_1d(1:mpi_nlpb_shadow_block_1d) = mpi_nlpb_shadow_nlpbz_starts_tmp(1:mpi_nlpb_shadow_block_1d)
      mpi_nlpb_shadow_block_counts_1d(1:mpi_nlpb_shadow_block_1d) = mpi_nlpb_shadow_nlpbz_counts_tmp(1:mpi_nlpb_shadow_block_1d)

      ! add extra two z grid indexes
      IF(mpi_nlpb_shadow_nlpbz(2) == mpi_nlpb_shadow_nlpbz(3)) THEN
        ! no z grids
        mpi_nlpb_shadow_nlpbz_block_1d = mpi_nlpb_shadow_block_1d
        ALLOCATE(mpi_nlpb_shadow_nlpbz_block_starts_1d(mpi_nlpb_shadow_nlpbz_block_1d))
        ALLOCATE(mpi_nlpb_shadow_nlpbz_block_counts_1d(mpi_nlpb_shadow_nlpbz_block_1d))
        mpi_nlpb_shadow_nlpbz_block_starts_1d(1:mpi_nlpb_shadow_block_1d) = &
          mpi_nlpb_shadow_block_starts_1d(1:mpi_nlpb_shadow_block_1d)
        mpi_nlpb_shadow_nlpbz_block_counts_1d(1:mpi_nlpb_shadow_block_1d) = &
          mpi_nlpb_shadow_block_counts_1d(1:mpi_nlpb_shadow_block_1d)
      ELSE
        IF(mpi_nlpb_shadow_nlpbz(2) + 1 == mpi_nlpb_shadow_nlpbz(3)) THEN
          ! only 1 z grid
          mpi_nlpb_shadow_nlpbz_block_1d = mpi_nlpb_shadow_block_1d + 1
          ALLOCATE(mpi_nlpb_shadow_nlpbz_block_starts_1d(mpi_nlpb_shadow_nlpbz_block_1d))
          ALLOCATE(mpi_nlpb_shadow_nlpbz_block_counts_1d(mpi_nlpb_shadow_nlpbz_block_1d))
          mpi_nlpb_shadow_nlpbz_block_starts_1d(1:mpi_nlpb_shadow_block_1d) = &
            mpi_nlpb_shadow_block_starts_1d(1:mpi_nlpb_shadow_block_1d)
          mpi_nlpb_shadow_nlpbz_block_counts_1d(1:mpi_nlpb_shadow_block_1d) = &
            mpi_nlpb_shadow_block_counts_1d(1:mpi_nlpb_shadow_block_1d)
          mpi_nlpb_shadow_nlpbz_block_starts_1d(mpi_nlpb_shadow_nlpbz_block_1d) = mpi_global_grids(mpi_nlpb_shadow_nlpbz(3))
          mpi_nlpb_shadow_nlpbz_block_counts_1d(mpi_nlpb_shadow_nlpbz_block_1d) = 1
        ELSE
          IF(mpi_nlpb_shadow_nlpbz(2) + 2 == mpi_nlpb_shadow_nlpbz(3)) THEN
            ! two z grids are continuous
            IF(mpi_global_grids(mpi_nlpb_shadow_nlpbz(2) + 1) + 1 == mpi_global_grids(mpi_nlpb_shadow_nlpbz(3))) THEN
              mpi_nlpb_shadow_nlpbz_block_1d = mpi_nlpb_shadow_block_1d + 1
              ALLOCATE(mpi_nlpb_shadow_nlpbz_block_starts_1d(mpi_nlpb_shadow_nlpbz_block_1d))
              ALLOCATE(mpi_nlpb_shadow_nlpbz_block_counts_1d(mpi_nlpb_shadow_nlpbz_block_1d))
              mpi_nlpb_shadow_nlpbz_block_starts_1d(1:mpi_nlpb_shadow_block_1d) = &
                 mpi_nlpb_shadow_block_starts_1d(1:mpi_nlpb_shadow_block_1d)
              mpi_nlpb_shadow_nlpbz_block_counts_1d(1:mpi_nlpb_shadow_block_1d) = &
                 mpi_nlpb_shadow_block_counts_1d(1:mpi_nlpb_shadow_block_1d)
              mpi_nlpb_shadow_nlpbz_block_starts_1d(mpi_nlpb_shadow_nlpbz_block_1d) = &
                 mpi_global_grids(mpi_nlpb_shadow_nlpbz(2) + 1)
              mpi_nlpb_shadow_nlpbz_block_counts_1d(mpi_nlpb_shadow_nlpbz_block_1d) = 2
            ELSE
              ! two z grids are not continuous
              mpi_nlpb_shadow_nlpbz_block_1d = mpi_nlpb_shadow_block_1d + 2
              ALLOCATE(mpi_nlpb_shadow_nlpbz_block_starts_1d(mpi_nlpb_shadow_nlpbz_block_1d))
              ALLOCATE(mpi_nlpb_shadow_nlpbz_block_counts_1d(mpi_nlpb_shadow_nlpbz_block_1d))
              mpi_nlpb_shadow_nlpbz_block_starts_1d(1:mpi_nlpb_shadow_block_1d) = &
                 mpi_nlpb_shadow_block_starts_1d(1:mpi_nlpb_shadow_block_1d)
              mpi_nlpb_shadow_nlpbz_block_counts_1d(1:mpi_nlpb_shadow_block_1d) = &
                 mpi_nlpb_shadow_block_counts_1d(1:mpi_nlpb_shadow_block_1d)
              mpi_nlpb_shadow_nlpbz_block_starts_1d(mpi_nlpb_shadow_block_1d + 1) = &
                 mpi_global_grids(mpi_nlpb_shadow_nlpbz(2) + 1)
              mpi_nlpb_shadow_nlpbz_block_starts_1d(mpi_nlpb_shadow_block_1d + 2) = &
                 mpi_global_grids(mpi_nlpb_shadow_nlpbz(2) + 2)
              mpi_nlpb_shadow_nlpbz_block_counts_1d(mpi_nlpb_shadow_block_1d + 1) = 1
              mpi_nlpb_shadow_nlpbz_block_counts_1d(mpi_nlpb_shadow_block_1d + 2) = 1
            END IF
          ELSE
            ! mpi_nlpb_shadow_nlpbz(3)- mpi_nlpb_shadow_nlpbz(2) should be 0 or 1 or 2
            WRITE(*,*)"size of loc_nlpb(3) in process of ",mpi_rank,"is wrong, please check nmefc_macom_par_dec.nc"
            STOP
          END IF
        END IF
      END IF

      ! prepare index for communication and I / O
      IF (mpi_rank == 0) THEN

        ! gather index in nlpb grids
         ALLOCATE (mpi_nlpb_block_1d_all(mpi_comp_procs))
         ALLOCATE (mpi_nlpb_block_1d_displs_all(mpi_comp_procs))

         CALL MPI_GATHER(mpi_nlpb_block_1d, 1, MPI_INTEGER, mpi_nlpb_block_1d_all, &
                         1, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         mpi_nlpb_block_1d_displs_all(1) = 0
         DO mpi_i = 2, mpi_comp_procs, 1
            mpi_nlpb_block_1d_displs_all(mpi_i) = mpi_nlpb_block_1d_displs_all(mpi_i-1) + mpi_nlpb_block_1d_all(mpi_i-1)
         END DO

         mpi_nlpb_block_1d_sum = mpi_nlpb_block_1d_displs_all(mpi_comp_procs) + mpi_nlpb_block_1d_all(mpi_comp_procs)

         ALLOCATE (mpi_nlpb_block_starts_1d_all(mpi_nlpb_block_1d_sum))
         ALLOCATE (mpi_nlpb_block_counts_1d_all(mpi_nlpb_block_1d_sum))

         CALL MPI_GATHERV(mpi_nlpb_block_starts_1d, mpi_nlpb_block_1d, MPI_INTEGER, &
                          mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_1d_all, &
                          mpi_nlpb_block_1d_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         CALL MPI_GATHERV(mpi_nlpb_block_counts_1d, mpi_nlpb_block_1d, MPI_INTEGER, &
                          mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_all, &
                          mpi_nlpb_block_1d_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         ! gather index in nlpb + shadow grids
         ALLOCATE(mpi_nlpb_shadow_block_1d_all(mpi_comp_procs))
         ALLOCATE(mpi_nlpb_shadow_block_1d_displs_all(mpi_comp_procs))
         CALL MPI_GATHER(mpi_nlpb_shadow_block_1d, 1, MPI_INTEGER, mpi_nlpb_shadow_block_1d_all, &
                         1, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         mpi_nlpb_shadow_block_1d_displs_all(1) = 0
         DO mpi_i = 2, mpi_comp_procs, 1
            mpi_nlpb_shadow_block_1d_displs_all(mpi_i) = mpi_nlpb_shadow_block_1d_displs_all(mpi_i - 1) + &
                                                 mpi_nlpb_shadow_block_1d_all(mpi_i - 1)
         END DO

         mpi_nlpb_shadow_block_1d_sum = mpi_nlpb_shadow_block_1d_displs_all(mpi_comp_procs) + &
                                             mpi_nlpb_shadow_block_1d_all(mpi_comp_procs)

         ALLOCATE (mpi_nlpb_shadow_block_starts_1d_all(mpi_nlpb_shadow_block_1d_sum))
         ALLOCATE (mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_sum))

         CALL MPI_GATHERV(mpi_nlpb_shadow_block_starts_1d, mpi_nlpb_shadow_block_1d, &
                          MPI_INTEGER, mpi_nlpb_shadow_block_starts_1d_all, &
                          mpi_nlpb_shadow_block_1d_all, mpi_nlpb_shadow_block_1d_displs_all, &
                          MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         CALL MPI_GATHERV(mpi_nlpb_shadow_block_counts_1d, mpi_nlpb_shadow_block_1d, &
                          MPI_INTEGER, mpi_nlpb_shadow_block_counts_1d_all, &
                          mpi_nlpb_shadow_block_1d_all, mpi_nlpb_shadow_block_1d_displs_all, &
                          MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         ! gather index in nlpb + shadow + z grids
         ALLOCATE(mpi_nlpb_shadow_nlpbz_block_1d_all(mpi_comp_procs))
         ALLOCATE(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(mpi_comp_procs))
         CALL MPI_GATHER(mpi_nlpb_shadow_nlpbz_block_1d, 1, MPI_INTEGER, &
                         mpi_nlpb_shadow_nlpbz_block_1d_all, 1, MPI_INTEGER, 0, &
                         mpi_comp_comm, mpi_err)

         mpi_nlpb_shadow_nlpbz_block_1d_displs_all(1) = 0
         DO mpi_i = 2, mpi_comp_procs, 1
            mpi_nlpb_shadow_nlpbz_block_1d_displs_all(mpi_i) = mpi_nlpb_shadow_nlpbz_block_1d_displs_all(mpi_i - 1) + &
                                                 mpi_nlpb_shadow_nlpbz_block_1d_all(mpi_i - 1)
         END DO

         mpi_nlpb_shadow_nlpbz_block_1d_sum = mpi_nlpb_shadow_nlpbz_block_1d_displs_all(mpi_comp_procs) + &
                                             mpi_nlpb_shadow_nlpbz_block_1d_all(mpi_comp_procs)

         ALLOCATE (mpi_nlpb_shadow_nlpbz_block_starts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_sum))
         ALLOCATE (mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_sum))

         CALL MPI_GATHERV(mpi_nlpb_shadow_nlpbz_block_starts_1d, mpi_nlpb_shadow_nlpbz_block_1d, MPI_INTEGER, &
                          mpi_nlpb_shadow_nlpbz_block_starts_1d_all, mpi_nlpb_shadow_nlpbz_block_1d_all, &
                          mpi_nlpb_shadow_nlpbz_block_1d_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         CALL MPI_GATHERV(mpi_nlpb_shadow_nlpbz_block_counts_1d, mpi_nlpb_shadow_nlpbz_block_1d, MPI_INTEGER, &
                          mpi_nlpb_shadow_nlpbz_block_counts_1d_all, mpi_nlpb_shadow_nlpbz_block_1d_all, &
                          mpi_nlpb_shadow_nlpbz_block_1d_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         DEALLOCATE(mpi_global_grids_all_sendbuf)
         DEALLOCATE(mpi_process_index_all)
         DEALLOCATE(mpi_process_index_all_sendbuf)
         DEALLOCATE(mpi_nlpb_shadow_nlpbz_all)
         DEALLOCATE(mpi_global_grids_all)

      ELSE

        ! send indexes in nlpb grids
         CALL MPI_GATHER(mpi_nlpb_block_1d, 1, MPI_INTEGER, mpi_nlpb_block_1d_all, &
                         1, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         CALL MPI_GATHERV(mpi_nlpb_block_starts_1d, mpi_nlpb_block_1d, MPI_INTEGER, &
                          mpi_nlpb_block_starts_1d_all, mpi_nlpb_block_1d_all, &
                          mpi_nlpb_block_1d_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         CALL MPI_GATHERV(mpi_nlpb_block_counts_1d, mpi_nlpb_block_1d, MPI_INTEGER, &
                          mpi_nlpb_block_counts_1d_all, mpi_nlpb_block_1d_all, &
                          mpi_nlpb_block_1d_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         ! send indexes in nlpb and shadow grids
         CALL MPI_GATHER(mpi_nlpb_shadow_block_1d, 1, MPI_INTEGER, mpi_nlpb_shadow_block_1d_all, 1, &
                         MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         CALL MPI_GATHERV(mpi_nlpb_shadow_block_starts_1d, mpi_nlpb_shadow_block_1d, &
                          MPI_INTEGER, mpi_nlpb_shadow_block_starts_1d_all, &
                          mpi_nlpb_shadow_block_1d_all, mpi_nlpb_shadow_block_1d_displs_all, &
                          MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         CALL MPI_GATHERV(mpi_nlpb_shadow_block_counts_1d, mpi_nlpb_shadow_block_1d, &
                          MPI_INTEGER, mpi_nlpb_shadow_block_counts_1d_all, &
                          mpi_nlpb_shadow_block_1d_all, mpi_nlpb_shadow_block_1d_displs_all, &
                          MPI_INTEGER, 0, mpi_comp_comm, mpi_err)


         ! send indexes in nlpb, shadow and z grids
         CALL MPI_GATHER(mpi_nlpb_shadow_nlpbz_block_1d, 1, MPI_INTEGER, mpi_nlpb_shadow_nlpbz_block_1d_all, 1, &
                         MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         CALL MPI_GATHERV(mpi_nlpb_shadow_nlpbz_block_starts_1d, mpi_nlpb_shadow_nlpbz_block_1d, MPI_INTEGER, &
                          mpi_nlpb_shadow_nlpbz_block_starts_1d_all, mpi_nlpb_shadow_nlpbz_block_1d_all, &
                          mpi_nlpb_shadow_nlpbz_block_1d_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         CALL MPI_GATHERV(mpi_nlpb_shadow_nlpbz_block_counts_1d, mpi_nlpb_shadow_nlpbz_block_1d, MPI_INTEGER, &
                          mpi_nlpb_shadow_nlpbz_block_counts_1d_all, mpi_nlpb_shadow_nlpbz_block_1d_all, &
                          mpi_nlpb_shadow_nlpbz_block_1d_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
      END IF

      DEALLOCATE(mpi_nlpb_shadow_nlpbz_starts_tmp)
      DEALLOCATE(mpi_nlpb_shadow_nlpbz_counts_tmp)
   END SUBROUTINE mpi_init_1d_index_buf

   ! global reduction for sum operation
   SUBROUTINE mpi_comp_dp_allsum(mpi_value, mpi_reduce_value)
      REAL(dp), INTENT(INOUT) :: mpi_value
      REAL(dp), OPTIONAL, INTENT(OUT) :: mpi_reduce_value
      REAL(dp) :: mpi_global_value
      CALL MPI_ALLREDUCE(mpi_value, mpi_global_value, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comp_comm, mpi_err)
      !? temporary change, you need to delete SNGL in the future
      IF (PRESENT(mpi_reduce_value)) THEN
         mpi_reduce_value = SNGL(mpi_global_value)
      ELSE
         mpi_value = SNGL(mpi_global_value)
      END IF
   END SUBROUTINE

   ! initialize index for data exchange in 2 dimensions
   SUBROUTINE mpi_init_2d_index_buf
      INTEGER :: mpi_i

      IF (mpi_rank == 0) THEN
         ALLOCATE (mpi_buf_2d_nlpb_ni_integer(mpi_total_nlpb, ni))
         ALLOCATE (mpi_buf_2d_nlpb_nlpbz_ni_integer(mpi_total_nlpbz, ni))
         ALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_ni_integer(mpi_total_nlpb_shadow * ni))
         ALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_nlpbz_ni_integer(mpi_total_nlpb_shadow_nlpbz * ni))
         ALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_nk_integer(mpi_total_nlpb_shadow * nk))
         ALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_nlpbz_nk_integer(mpi_total_nlpb_shadow_nlpbz * nk))
         ALLOCATE (mpi_buf_2d_nlpb_nk_i2(mpi_total_nlpb, nk))
         ALLOCATE (mpi_buf_2d_nlpb_nlpbz_nk_i2(mpi_total_nlpbz, nk))
         ALLOCATE (mpi_buf_2d_nlpb_nk_wp(mpi_total_nlpb, nk))
         ALLOCATE (mpi_buf_2d_nlpb_sbct_wp(mpi_total_nlpb, sbct))
         ALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_nk_wp(mpi_total_nlpb_shadow * nk))
         ALLOCATE (mpi_adjust_buf_2d_nlpb_shadow_sbct_wp(mpi_total_nlpb_shadow * sbct))

         ALLOCATE (mpi_adjust_2d_nlpb_shadow_ni_counts_all(mpi_comp_procs))
         ALLOCATE (mpi_adjust_2d_nlpb_shadow_ni_displs_all(mpi_comp_procs))
         ALLOCATE (mpi_adjust_2d_nlpb_shadow_nk_counts_all(mpi_comp_procs))
         ALLOCATE (mpi_adjust_2d_nlpb_shadow_nk_displs_all(mpi_comp_procs))
         ALLOCATE (mpi_adjust_2d_nlpb_shadow_nlpbz_ni_counts_all(mpi_comp_procs))
         ALLOCATE (mpi_adjust_2d_nlpb_shadow_nlpbz_ni_displs_all(mpi_comp_procs))
         ALLOCATE (mpi_adjust_2d_nlpb_shadow_nlpbz_nk_counts_all(mpi_comp_procs))
         ALLOCATE (mpi_adjust_2d_nlpb_shadow_nlpbz_nk_displs_all(mpi_comp_procs))
         ALLOCATE (mpi_adjust_2d_nlpb_shadow_sbct_counts_all(mpi_comp_procs))
         ALLOCATE (mpi_adjust_2d_nlpb_shadow_sbct_displs_all(mpi_comp_procs))

         ALLOCATE(mpi_2d_nlpb_nk_counts_all(mpi_comp_procs))
         ALLOCATE(mpi_2d_nlpb_nk_counts_displs_all(mpi_comp_procs))

         mpi_adjust_2d_nlpb_shadow_ni_counts_all(1:mpi_comp_procs) = mpi_nlpb_shadow_counts_all(1:mpi_comp_procs)* ni
         mpi_adjust_2d_nlpb_shadow_nlpbz_ni_counts_all(1:mpi_comp_procs) = mpi_nlpb_shadow_nlpbz_counts_all(1:mpi_comp_procs)* ni
         mpi_adjust_2d_nlpb_shadow_nk_counts_all(1:mpi_comp_procs) = mpi_nlpb_shadow_counts_all(1:mpi_comp_procs)* nk
         mpi_adjust_2d_nlpb_shadow_nlpbz_nk_counts_all(1:mpi_comp_procs) = mpi_nlpb_shadow_nlpbz_counts_all(1:mpi_comp_procs)* nk
         mpi_adjust_2d_nlpb_shadow_sbct_counts_all(1:mpi_comp_procs) = mpi_nlpb_shadow_counts_all(1:mpi_comp_procs)* sbct
         mpi_2d_nlpb_nk_counts_all(1:mpi_comp_procs)= mpi_nlpb_counts_all(1:mpi_comp_procs)* nk

         mpi_adjust_2d_nlpb_shadow_ni_displs_all(1) = 0
         mpi_adjust_2d_nlpb_shadow_nlpbz_ni_displs_all(1) = 0
         mpi_adjust_2d_nlpb_shadow_nk_displs_all(1) = 0
         mpi_adjust_2d_nlpb_shadow_nlpbz_nk_displs_all(1) = 0
         mpi_adjust_2d_nlpb_shadow_sbct_displs_all(1) = 0
         mpi_2d_nlpb_nk_counts_displs_all(1) = 0

         DO mpi_i = 2, mpi_comp_procs, 1
            mpi_adjust_2d_nlpb_shadow_ni_displs_all(mpi_i) = mpi_adjust_2d_nlpb_shadow_ni_displs_all(mpi_i - 1) + &
                                                              mpi_adjust_2d_nlpb_shadow_ni_counts_all(mpi_i - 1)
            mpi_adjust_2d_nlpb_shadow_nlpbz_ni_displs_all(mpi_i) = mpi_adjust_2d_nlpb_shadow_nlpbz_ni_displs_all(mpi_i - 1) + &
                                                      mpi_adjust_2d_nlpb_shadow_nlpbz_ni_counts_all(mpi_i - 1)
            mpi_adjust_2d_nlpb_shadow_nk_displs_all(mpi_i) = mpi_adjust_2d_nlpb_shadow_nk_displs_all(mpi_i - 1) + &
                                                              mpi_adjust_2d_nlpb_shadow_nk_counts_all(mpi_i - 1)
            mpi_adjust_2d_nlpb_shadow_nlpbz_nk_displs_all(mpi_i) = mpi_adjust_2d_nlpb_shadow_nlpbz_nk_displs_all(mpi_i - 1) + &
                                                                    mpi_adjust_2d_nlpb_shadow_nlpbz_nk_counts_all(mpi_i - 1)
            mpi_adjust_2d_nlpb_shadow_sbct_displs_all(mpi_i) = mpi_adjust_2d_nlpb_shadow_sbct_displs_all(mpi_i - 1) + &
                                                              mpi_adjust_2d_nlpb_shadow_sbct_counts_all(mpi_i - 1)
            mpi_2d_nlpb_nk_counts_displs_all(mpi_i) = mpi_2d_nlpb_nk_counts_displs_all(mpi_i - 1) + &
                                                        mpi_2d_nlpb_nk_counts_all(mpi_i - 1)
         END DO
      END IF
   END SUBROUTINE

   ! increase array size for data exchange
   SUBROUTINE mpi_init_arrays

      INTEGER :: dimen(2), dimen_3d(3),i


      ! allocate buffer for io read and adjust
      IF (mpi_rank == 0) THEN
        ALLOCATE (mpi_buf_1d_nlpb_wp(mpi_total_nlpb))
        ALLOCATE (mpi_buf_1d_nlpb_nlpbz_wp(mpi_total_nlpbz))
        ALLOCATE (mpi_buf_1d_nlpb_nlpbz_integer(mpi_total_nlpbz))
        ALLOCATE (mpi_buf_1d_nlpb_integer(mpi_total_nlpb))
        ALLOCATE (mpi_buf_1d_nlpb_i2(mpi_total_nlpb))
        ALLOCATE (mpi_adjust_buf_1d_nlpb_shadow_wp(mpi_total_nlpb_shadow))
        ALLOCATE (mpi_adjust_buf_1d_nlpb_shadow_nlpbz_wp(mpi_total_nlpb_shadow_nlpbz))
        ALLOCATE (mpi_adjust_buf_1d_nlpb_shadow_nlpbz_integer(mpi_total_nlpb_shadow_nlpbz))
        ALLOCATE (mpi_adjust_buf_1d_nlpb_shadow_integer(mpi_total_nlpb_shadow))
        ALLOCATE (mpi_adjust_buf_1d_nlpb_shadow_i2(mpi_total_nlpb_shadow))

        if (ln_bdy) then
           ALLOCATE (mpi_buf_1d_nlbdy_wp(mpi_total_nlbdy))
           ALLOCATE (mpi_adjust_buf_1d_nlbdy_wp(mpi_total_nlbdy_shadow))
           ALLOCATE (mpi_buf_1d_nlbdy_integer(mpi_total_nlbdy))
           ALLOCATE (mpi_adjust_buf_1d_nlbdy_integer(mpi_total_nlbdy_shadow))

           ALLOCATE (mpi_buf_2d_nlbdy_nk_wp(mpi_total_nlbdy,nk))
           ALLOCATE (mpi_adjust_buf_2d_nlbdy_nk_wp(mpi_total_nlbdy_shadow * nk))
           ALLOCATE (mpi_adjust_2d_nlbdy_nk_counts_all(mpi_bdy_outdegree))
           ALLOCATE (mpi_adjust_2d_nlbdy_nk_displs_all(mpi_bdy_outdegree))

           DO i = 1,mpi_bdy_outdegree
              mpi_adjust_2d_nlbdy_nk_counts_all(i) = mpi_bdy_counts_1d_all(i)* nk
           END DO

           mpi_adjust_2d_nlbdy_nk_displs_all(1)= 0
           DO i = 2,mpi_bdy_outdegree
              mpi_adjust_2d_nlbdy_nk_displs_all(i)= mpi_adjust_2d_nlbdy_nk_displs_all(i-1) &
                 + mpi_adjust_2d_nlbdy_nk_counts_all(i-1)
           END DO
        end if

      ELSE
        ! VALUE IS NOT USED, ONLY FOR AVOIDING MPI WARNINGS IN SOME CASES

        ALLOCATE(mpi_nlpb_shadow_counts_all(1))
        ALLOCATE(mpi_adjust_buf_1d_nlpb_shadow_wp(1))
        ALLOCATE(mpi_nlpb_shadow_displs_all(1))
      END IF

      ! allocate send buffer for each exchaning arrays
      ALLOCATE (mpi_sendbuf_1d(mpi_send_indexes_1d_counts_sum))
      ALLOCATE (mpi_sendbuf_2d_nk_f(mpi_send_indexes_1d_counts_sum*nk))
      ALLOCATE (mpi_sendbuf_2d_nk_int(mpi_send_indexes_1d_counts_sum*nk))
      ALLOCATE (mpi_sendbuf_2d_nkp1(mpi_send_indexes_1d_counts_sum*nkp1))
#ifdef SEAICE_ITD
      if(mitice_on) &
        ALLOCATE (mpi_sendbuf_2d_itd_f(mpi_send_indexes_1d_counts_sum*nITD))
#endif

   END SUBROUTINE mpi_init_arrays

   ! broadcast information in namelist from root process to other processes among computing processes
   SUBROUTINE mpi_namelist_pack_and_bcast
      mpi_buf_size = 8000

      ALLOCATE (mpi_buf(mpi_buf_size), STAT = mpi_mem_status)

      CALL mpi_namelist_pack

      CALL MPI_BCAST(mpi_position, 1, MPI_INTEGER, mpi_root_comp, mpi_comp_comm, mpi_err)
      CALL MPI_BCAST(mpi_buf, mpi_position, MPI_PACKED, mpi_root_comp, mpi_comp_comm, mpi_err)
      ! mpi_position is the real size of buffer in mpi_root_comp
      DEALLOCATE (mpi_buf)
   END SUBROUTINE mpi_namelist_pack_and_bcast

   ! broadcast information in namelist from root process to other processes among computing processes
   SUBROUTINE mitice_mpi_initfixed_pack_and_bcast
      mpi_buf_size = 128
      mpi_position=0
      ALLOCATE (mpi_buf(mpi_buf_size), STAT = mpi_mem_status)
      CALL MPI_PACK(SWFracB, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_mcPheePiston, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
#ifdef SEAICE_ITD
      CALL MPI_PACK(Hlimit, nITD+1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
#endif
      ! becasue size of mpi_buf is small, so we send the whole buffer at one
      ! time, delete one time of sending operation to send size of buffer first
      CALL MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_root_comp, mpi_comp_comm, mpi_err)
      DEALLOCATE (mpi_buf)
   END SUBROUTINE mitice_mpi_initfixed_pack_and_bcast

   ! receive information in namelist from root process for other processes among computing processes
   SUBROUTINE mitice_mpi_initfixed_unpack_and_bcast
      mpi_buf_size=128
      mpi_position=0
      ALLOCATE (mpi_buf(mpi_buf_size))
      CALL MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_root_comp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SWFracB, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_mcPheePiston, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
#ifdef SEAICE_ITD
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, Hlimit, nITD+1, mpi_real_wp, mpi_comp_comm, mpi_err)
#endif
      DEALLOCATE (mpi_buf)
   END SUBROUTINE mitice_mpi_initfixed_unpack_and_bcast

   ! broadcast information in namelist in mitice from root process to other processes among computing processes
   SUBROUTINE mitice_mpi_namelist_pack_and_bcast
      mpi_buf_size = 8000
      mpi_position=0
      ALLOCATE (mpi_buf(mpi_buf_size), STAT = mpi_mem_status)
      CALL MPI_PACK(SEAICEuseDYNAMICS, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEuseEVP, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEuseFREEDRIFT, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEuseTEM, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEuseTilt, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(useHB87stressCoupling, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(usePW79thermodynamics, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      !CALL MPI_PACK(useMaykutSatVapPoly, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEuseFlooding, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_growMeltByConv, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEadvHeff, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEadvArea, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEadvSnow, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEadvSalt, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEaddSnowMass, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEmomAdvection, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_clipVelocities, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_no_slip, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEetaZmethod, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(IMAX_TICE, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(postSolvTempIter, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEdiffKhHeff, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEdiffKhSnow, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEdiffKhArea, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEdiffKhSalt, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_deltaTtherm, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_deltaTdyn, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_deltaTevp, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_elasticParm, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_evpTauRelax, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_evpDampC, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEnEVPstarSteps, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_evpAlpha, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_evpBeta, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEaEVPcoeff, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEaEVPcStar, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEaEVPalphaMin, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_zetaMin, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_zetaMaxFac, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEuseLinRemapITD, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(useHibler79IceStrength, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEpartFunc, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEredistFunc, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEridgingIterMax, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEsimpleRidging, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEsnowFracRidge, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEgStar, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEhStar, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEaStar, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEshearParm, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEmuRidging, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEmaxRaft, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_cf, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEpresH0, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEpresPow0, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEpresPow1, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_initialHEFF, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_areaGainFormula, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_areaLossFormula, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_doOpenWaterGrowth, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_doOpenWaterMelt, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_rhoAir, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_rhoIce, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_rhoSnow, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ICE2WATR, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_cpAir, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEscaleSurfStress, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_drag, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_waterDrag, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEdWatMin, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_dryIceAlb, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_wetIceAlb, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_drySnowAlb, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_wetSnowAlb, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(HO, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_drag_south, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_waterDrag_south, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_dryIceAlb_south, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_wetIceAlb_south, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_drySnowAlb_south, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_wetSnowAlb_south, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(HO_south, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_cBasalStar, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEbasalDragU0, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEbasalDragK1, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEbasalDragK2, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_wetAlbTemp, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_strength, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_cStar, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_eccen, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEpressReplFac, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_tensilFac, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_tensilDepth, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_lhFusion, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_lhEvap, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_dalton, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_iceConduct, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_snowConduct, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_emissivity, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_ice_emiss, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_snow_emiss, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_snowThick, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_shortwave, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(OCEAN_drag, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_tempFrz0, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_dTempFrz_dS, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_salt0, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_saltFrac, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEstressFactor, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
!      CALL MPI_PACK(SEAICE_availHeatTaper, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
!      CALL MPI_PACK(SEAICE_mcPheePiston, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_frazilFrac, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_mcPheeTaper, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_mcPheeStepFunc, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
!      CALL MPI_PACK(SEAICE_gamma_t, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
!      CALL MPI_PACK(SEAICE_gamma_t_frz, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_availHeatFrac, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
!      CALL MPI_PACK(SEAICE_availHeatFracFrz, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_PDF, nITD, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(AreaFile, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(HeffFile, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(uIceFile, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(vIceFile, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(HsnowFile, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(HsaltFile, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEheatConsFix, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_multDim, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_useMultDimSnow, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_deltaMin, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_area_reg, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_hice_reg, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_area_floor, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_area_max, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_tauAreaObsRelax, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_airTurnAngle, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_waterTurnAngle, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(MIN_ATEMP, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(MIN_LWDOWN, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(MIN_TICE, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_EPS, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_EPS_SQ, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEwriteState, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEuseEVPpickup, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEuseEVPstar, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICEuseEVPrev, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_monFreq, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_dumpFreq, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(SEAICE_taveFreq, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(dumpFileIntv, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)

      CALL MPI_BCAST(mpi_position, 1, MPI_INTEGER, mpi_root_comp, mpi_comp_comm, mpi_err)
      CALL MPI_BCAST(mpi_buf, mpi_position, MPI_PACKED, mpi_root_comp, mpi_comp_comm, mpi_err)
      ! mpi_position is the real size of buffer in mpi_root_comp
      DEALLOCATE (mpi_buf)
   END SUBROUTINE

   ! receive information in namelist from root process for other processes among computing processes
   SUBROUTINE mpi_namelist_unpack_and_bcast
      CALL MPI_BCAST(mpi_buf_size, 1, MPI_INTEGER, mpi_root_comp, mpi_comp_comm, mpi_err)

      ALLOCATE (mpi_buf(mpi_buf_size))
      CALL MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_root_comp, mpi_comp_comm, mpi_err)

      CALL mpi_namelist_unpack
      DEALLOCATE (mpi_buf)
   END SUBROUTINE mpi_namelist_unpack_and_bcast

   ! receive information in namelist from root process for other processes among computing processes
   SUBROUTINE mitice_mpi_namelist_unpack_and_bcast
      CALL MPI_BCAST(mpi_buf_size, 1, MPI_INTEGER, mpi_root_comp, mpi_comp_comm, mpi_err)

      ALLOCATE (mpi_buf(mpi_buf_size))
      CALL MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, mpi_root_comp, mpi_comp_comm, mpi_err)

      mpi_position=0
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEuseDYNAMICS, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEuseEVP, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEuseFREEDRIFT, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEuseTEM, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEuseTilt, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, useHB87stressCoupling, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, usePW79thermodynamics, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      !CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, useMaykutSatVapPoly, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEuseFlooding, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_growMeltByConv, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEadvHeff, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEadvArea, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEadvSnow, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEadvSalt, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEaddSnowMass, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEmomAdvection, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_clipVelocities, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_no_slip, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEetaZmethod, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, IMAX_TICE, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, postSolvTempIter, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEdiffKhHeff, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEdiffKhSnow, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEdiffKhArea, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEdiffKhSalt, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_deltaTtherm, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_deltaTdyn, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_deltaTevp, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_elasticParm, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_evpTauRelax, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_evpDampC, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEnEVPstarSteps, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_evpAlpha, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_evpBeta, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEaEVPcoeff, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEaEVPcStar, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEaEVPalphaMin, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_zetaMin, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_zetaMaxFac, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEuseLinRemapITD, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, useHibler79IceStrength, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEpartFunc, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEredistFunc, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEridgingIterMax, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEsimpleRidging, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEsnowFracRidge, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEgStar, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEhStar, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEaStar, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEshearParm, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEmuRidging, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEmaxRaft, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_cf, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEpresH0, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEpresPow0, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEpresPow1, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_initialHEFF, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_areaGainFormula, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_areaLossFormula, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_doOpenWaterGrowth, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_doOpenWaterMelt, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_rhoAir, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_rhoIce, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_rhoSnow, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ICE2WATR, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_cpAir, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEscaleSurfStress, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_drag, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_waterDrag, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEdWatMin, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_dryIceAlb, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_wetIceAlb, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_drySnowAlb, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_wetSnowAlb, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, HO, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_drag_south, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_waterDrag_south, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_dryIceAlb_south, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_wetIceAlb_south, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_drySnowAlb_south, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_wetSnowAlb_south, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, HO_south, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_cBasalStar, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEbasalDragU0, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEbasalDragK1, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEbasalDragK2, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_wetAlbTemp, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_strength, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_cStar, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_eccen, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEpressReplFac, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_tensilFac, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_tensilDepth, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_lhFusion, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_lhEvap, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_dalton, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_iceConduct, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_snowConduct, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_emissivity, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_ice_emiss, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_snow_emiss, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_snowThick, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_shortwave, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, OCEAN_drag, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_tempFrz0, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_dTempFrz_dS, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_salt0, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_saltFrac, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEstressFactor, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
!      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_availHeatTaper, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
!      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_mcPheePiston, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_frazilFrac, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_mcPheeTaper, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_mcPheeStepFunc, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
!      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_gamma_t, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
!      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_gamma_t_frz, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_availHeatFrac, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
!      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_availHeatFracFrz, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_PDF, nITD, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, AreaFile, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, HeffFile, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, uIceFile, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, vIceFile, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, HsnowFile, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, HsaltFile, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEheatConsFix, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_multDim, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_useMultDimSnow, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_deltaMin, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_area_reg, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_hice_reg, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_area_floor, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_area_max, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_tauAreaObsRelax, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_airTurnAngle, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_waterTurnAngle, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, MIN_ATEMP, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, MIN_LWDOWN, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, MIN_TICE, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_EPS, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_EPS_SQ, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEwriteState, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEuseEVPpickup, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEuseEVPstar, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICEuseEVPrev, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_monFreq, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_dumpFreq, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, SEAICE_taveFreq, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dumpFileIntv, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      DEALLOCATE (mpi_buf)
   END SUBROUTINE 

   ! pack information in namelist in root process
   SUBROUTINE mpi_namelist_pack
      CALL MPI_PACK(ndomain, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(idom, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(parent, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nRatio, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nIter0, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nIterMax, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nOutFreq, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fOutFreq, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nResFreq, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(startyear, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(startmon, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(startday, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(starthour, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(startmin, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(startsec, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(dTtracer, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(dTmom, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(dTsurf, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(gammas, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(gammat, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ffsfrs, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ffsfrt, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nffsfrs, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nffsfrt, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fsfrsavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fsfrtavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(sbc_cli, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(sbc_cyc, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ffuv, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fftq, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ffrad, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ffprc, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ffslp, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ffrunoff, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nffuv, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nfftq, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nffrad, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nffprc, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nffslp, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nffrunoff, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fuvavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ftqavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fradavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fprcavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fslpavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(frunoffavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ln_qdew, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ln_rain, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ln_EmPmRevise, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(zqt, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(zuv, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(rhoref, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ln_bdy, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fbvec, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fbt, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fbs, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fbpbt, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nfbvec, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nfbt, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nfbs, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(nfbpbt, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(bvecavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(btavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(bsavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(bpbtavg, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fsbc_dir, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fbdy_dir, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(fout_dir, lc, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(restart_in, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(assim_in, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(restart_out, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(date_res, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ViscAhDCon, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ViscAhZCon, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ViscA4DCon, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ViscA4ZCon, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(harmonic, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(biharmonic, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(bottomDrag, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(BottomDragMax, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(KappaRMCon, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(abEps, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(rStarFacLow, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(rStarFacUp, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ln_pbt_base, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ln_bous, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ln_tide, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(ntracer, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(rn_abs, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(rn_si0, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(KappaRTCon, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(harmonicT, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(biharmonicT, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(diffKhCon, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      CALL MPI_PACK(diffK4Con, 1, mpi_real_wp, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)
      !YY: ice module on/off
      CALL MPI_PACK(mitice_on, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, mpi_comp_comm, mpi_err)

   END SUBROUTINE mpi_namelist_pack

   ! unpack information in namelist for non - root processes
   SUBROUTINE mpi_namelist_unpack
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ndomain, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, idom, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, parent, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nRatio, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nIter0, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nIterMax, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nOutFreq, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fOutFreq, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nResFreq, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, startyear, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, startmon, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, startday, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, starthour, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, startmin, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, startsec, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dTtracer, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dTmom, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dTsurf, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, gammas, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, gammat, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ffsfrs, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ffsfrt, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nffsfrs, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nffsfrt, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fsfrsavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fsfrtavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, sbc_cli, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, sbc_cyc, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ffuv, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fftq, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ffrad, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ffprc, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ffslp, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ffrunoff, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nffuv, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nfftq, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nffrad, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nffprc, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nffslp, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nffrunoff, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fuvavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ftqavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fradavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fprcavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fslpavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, frunoffavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ln_qdew, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ln_rain, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ln_EmPmRevise, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, zqt, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, zuv, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, rhoref, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ln_bdy, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fbvec, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fbt, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fbs, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fbpbt, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nfbvec, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nfbt, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nfbs, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nfbpbt, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, bvecavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, btavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, bsavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, bpbtavg, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fsbc_dir, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fbdy_dir, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, fout_dir, lc, MPI_CHARACTER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, restart_in, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, assim_in, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, restart_out, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, date_res, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ViscAhDCon, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ViscAhZCon, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ViscA4DCon, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ViscA4ZCon, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, harmonic, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, biharmonic, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, bottomDrag, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, BottomDragMax, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, KappaRMCon, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, abEps, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, rStarFacLow, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, rStarFacUp, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ln_pbt_base, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ln_bous, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ln_tide, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, ntracer, 1, MPI_INTEGER, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, rn_abs, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, rn_si0, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, KappaRTCon, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, harmonicT, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, biharmonicT, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, diffKhCon, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, diffK4Con, 1, mpi_real_wp, mpi_comp_comm, mpi_err)
      CALL MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, mitice_on, 1, MPI_LOGICAL, mpi_comp_comm, mpi_err)

   END SUBROUTINE mpi_namelist_unpack

   ! prepare sending buffer of dp for each process,for all arrays which have same position for collecting sending buffer
   SUBROUTINE mpi_wp_1d_prepare_sendbuf(mpi_source_buf, mpi_send_buf)
      REAL(wp), DIMENSION(:), INTENT(IN) :: mpi_source_buf
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_send_buf
      INTEGER :: i
      !$acc kernels present(mpi_send_buf,mpi_source_buf)
      !$acc loop independent
      DO i = 1, mpi_send_indexes_1d_counts_sum
         mpi_send_buf(i) = mpi_source_buf(mpi_send_indexes_1d(i))
      END DO
      !$acc end kernels
   END SUBROUTINE mpi_wp_1d_prepare_sendbuf

   ! prepare sending buffer of i2 for each process,for all arrays which have same position for collecting sending buffer
   SUBROUTINE mpi_i2_1d_prepare_sendbuf(mpi_source_buf, mpi_send_buf)
      INTEGER(i2), DIMENSION(:), INTENT(IN) :: mpi_source_buf
      INTEGER(i2), DIMENSION(:), INTENT(INOUT) :: mpi_send_buf
      INTEGER :: i
      !$acc kernels present(mpi_send_buf,mpi_source_buf)
      !$acc loop independent
      DO i = 1, mpi_send_indexes_1d_counts_sum
         mpi_send_buf(i) = mpi_source_buf(mpi_send_indexes_1d(i))
      END DO
      !$acc end kernels
   END SUBROUTINE mpi_i2_1d_prepare_sendbuf

   ! prepare sending buffer of dp with 3 dimensions for each process
   SUBROUTINE mpi_dp_3d_prepare_sendbuf(mpi_source_buf, mpi_send_buf, dimen1, dimen2, mpi_send_indexes_1d_counts_tmp, &
                                        mpi_send_indexes_1d_displs_tmp)
      REAL(wp), DIMENSION(:, :, :), INTENT(IN) :: mpi_source_buf
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_send_buf
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, INTENT(IN) :: dimen2
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_indexes_1d_displs_tmp
      INTEGER :: i, j, k, m, kcount, ij, startpos

      mpi_send_indexes_1d_counts_tmp(1) = mpi_send_indexes_1d_counts(1)* dimen1 * dimen2
      mpi_send_indexes_1d_displs_tmp(1) = 0

      DO i = 2, mpi_graph_outdegree
         mpi_send_indexes_1d_counts_tmp(i) = mpi_send_indexes_1d_counts(i)* dimen1 * dimen2
         mpi_send_indexes_1d_displs_tmp(i) = mpi_send_indexes_1d_displs_tmp(i - 1) + mpi_send_indexes_1d_counts_tmp(i - 1)
      END DO

      DO i = 1, mpi_graph_outdegree
         startpos = mpi_send_indexes_1d_displs_tmp(i)
         !$acc parallel present(mpi_send_indexes_1d_counts,mpi_send_buf,mpi_source_buf,mpi_send_indexes_1d, &
         !$acc                 mpi_send_indexes_1d_displs), &
         !$acc         copyin(i,startpos)
         !$acc loop collapse(3) private(kcount, ij) independent
         DO m = 1, dimen2
            DO j = 1, dimen1
               DO k = 1,mpi_send_indexes_1d_counts(i)
                 kcount = startpos + (j - 1)* mpi_send_indexes_1d_counts(i) +&
                     (m - 1)* dimen1 * mpi_send_indexes_1d_counts(i) + k
                 ij = mpi_send_indexes_1d(mpi_send_indexes_1d_displs(i) + k)
                 mpi_send_buf(kcount) = mpi_source_buf(ij,j,m)
               END DO
            END DO
         END DO
         !$acc end parallel
      END DO

   END SUBROUTINE mpi_dp_3d_prepare_sendbuf

   ! prepare sending buffer of i2 with 2 dimensions for each process
   SUBROUTINE mpi_i2_2d_prepare_sendbuf(mpi_source_buf, mpi_send_buf, dimen1, mpi_send_indexes_1d_counts_tmp, &
                                        mpi_send_indexes_1d_displs_tmp)
      INTEGER(i2), DIMENSION(:, :), INTENT(IN) :: mpi_source_buf
      INTEGER(i2), DIMENSION(:), INTENT(INOUT) :: mpi_send_buf
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_indexes_1d_displs_tmp
      INTEGER :: i, j, k, kcount, ij, startpos

      mpi_send_indexes_1d_counts_tmp(1) = mpi_send_indexes_1d_counts(1)* dimen1
      mpi_send_indexes_1d_displs_tmp(1) = 0
      DO i = 2, mpi_graph_outdegree
         mpi_send_indexes_1d_counts_tmp(i) = mpi_send_indexes_1d_counts(i)* dimen1
         mpi_send_indexes_1d_displs_tmp(i) = mpi_send_indexes_1d_displs_tmp(i - 1) + mpi_send_indexes_1d_counts_tmp(i - 1)
      END DO

      DO i = 1, mpi_graph_outdegree
         startpos = mpi_send_indexes_1d_displs_tmp(i)
         !$acc parallel present(mpi_send_indexes_1d_counts,mpi_send_buf,mpi_source_buf,mpi_send_indexes_1d, &
         !$acc                 mpi_send_indexes_1d_displs), &
         !$acc         copyin(i,startpos)
         !$acc loop collapse(2) private(kcount, ij) independent
         DO j = 1, dimen1
            DO k = 1, mpi_send_indexes_1d_counts(i)
               kcount = startpos + (j - 1)* mpi_send_indexes_1d_counts(i) + k
               ij = mpi_send_indexes_1d(mpi_send_indexes_1d_displs(i) + k)
               mpi_send_buf(kcount) = mpi_source_buf(ij,j)
            END DO
         END DO
         !$acc end parallel
      END DO

   END SUBROUTINE mpi_i2_2d_prepare_sendbuf

   ! prepare sending buffer of dp with 2 dimensions for each process
   SUBROUTINE mpi_wp_2d_prepare_sendbuf(mpi_source_buf, mpi_send_buf, dimen1, mpi_send_indexes_1d_counts_tmp, &
                                        mpi_send_indexes_1d_displs_tmp)
      REAL(wp), DIMENSION(:, :), INTENT(IN) :: mpi_source_buf
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_send_buf
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_indexes_1d_displs_tmp
      INTEGER :: i, j, k, kcount, ij, startpos

      mpi_send_indexes_1d_counts_tmp(1) = mpi_send_indexes_1d_counts(1)* dimen1
      mpi_send_indexes_1d_displs_tmp(1) = 0
      DO i = 2, mpi_graph_outdegree
         mpi_send_indexes_1d_counts_tmp(i) = mpi_send_indexes_1d_counts(i)* dimen1
         mpi_send_indexes_1d_displs_tmp(i) = mpi_send_indexes_1d_displs_tmp(i - 1) + mpi_send_indexes_1d_counts_tmp(i - 1)
      END DO

      DO i = 1, mpi_graph_outdegree
         startpos = mpi_send_indexes_1d_displs_tmp(i)
         !$acc parallel present(mpi_send_indexes_1d_counts,mpi_send_buf,mpi_source_buf,mpi_send_indexes_1d, &
         !$acc                 mpi_send_indexes_1d_displs), &
         !$acc         copyin(i,startpos)
         !$acc loop collapse(2) private(kcount, ij) independent
         DO j = 1, dimen1
            DO k = 1, mpi_send_indexes_1d_counts(i)
              mpi_send_buf(startpos + (j - 1)* mpi_send_indexes_1d_counts(i) + k) = &
                mpi_source_buf(mpi_send_indexes_1d(mpi_send_indexes_1d_displs(i) + k),j)
            END DO
         END DO
         !$acc end parallel
      END DO

   END SUBROUTINE mpi_wp_2d_prepare_sendbuf

   ! prepare sending buffer of wp with 2 dimensions for each process for reading nmefc_macom_par_dec.nc
   SUBROUTINE mpi_integer_2d_prepare_sendbuf_for_par_dec(mpi_source_buf, &
      mpi_send_buf, dimen1, mpi_send_indexes_1d_counts_tmp)

      INTEGER, DIMENSION(:, :), INTENT(IN) :: mpi_source_buf
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_buf
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, DIMENSION(:), INTENT(IN) :: mpi_send_indexes_1d_counts_tmp
      INTEGER :: i
      INTEGER :: current_p

      current_p = 1
      DO i = 1, dimen1
        mpi_send_buf(current_p:current_p + mpi_send_indexes_1d_counts_tmp(i) - 1) = &
          mpi_source_buf(1:mpi_send_indexes_1d_counts_tmp(i), i)
        current_p = current_p + mpi_send_indexes_1d_counts_tmp(i)
      END DO

   END SUBROUTINE

   ! prepare index counts and displacement for data exchange in 2 dimension
   SUBROUTINE mpi_2d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, dimen1)
      INTEGER, DIMENSION(:), INTENT(INOUT) ::mpi_recv_indexes_counts_tmp
      INTEGER, DIMENSION(:), INTENT(INOUT) ::mpi_recv_indexes_displs_tmp
      INTEGER, INTENT(IN) :: dimen1
      INTEGER :: i
      mpi_recv_indexes_counts_tmp(1) = mpi_recv_indexes_1d_counts(1)* dimen1
      mpi_recv_indexes_displs_tmp(1) = 0
      DO i = 2, mpi_graph_indegree
         mpi_recv_indexes_counts_tmp(i) = mpi_recv_indexes_1d_counts(i)* dimen1
         mpi_recv_indexes_displs_tmp(i) = mpi_recv_indexes_displs_tmp(i - 1) + mpi_recv_indexes_counts_tmp(i - 1)
      END DO

   END SUBROUTINE

   ! prepare index counts and displacement for data exchange in 3 dimension
   SUBROUTINE mpi_dp_3d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, dimen1, dimen2)
      INTEGER, DIMENSION(:), INTENT(INOUT) ::mpi_recv_indexes_counts_tmp
      INTEGER, DIMENSION(:), INTENT(INOUT) ::mpi_recv_indexes_displs_tmp
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, INTENT(IN) :: dimen2
      INTEGER :: i
      mpi_recv_indexes_counts_tmp(1) = mpi_recv_indexes_1d_counts(1)*dimen1*dimen2
      mpi_recv_indexes_displs_tmp(1) = 0
      DO i = 2, mpi_graph_indegree
         mpi_recv_indexes_counts_tmp(i) = mpi_recv_indexes_1d_counts(i)*dimen1*dimen2
         mpi_recv_indexes_displs_tmp(i) = mpi_recv_indexes_displs_tmp(i - 1) + mpi_recv_indexes_counts_tmp(i - 1)
      END DO

   END SUBROUTINE

   ! adjust receiving data of dp in the right order in 2 dimension
   SUBROUTINE mpi_wp_2d_prepare_recvbuf(mpi_recv_buf_tmp, mpi_recv_data, dimen1)
      REAL(wp), DIMENSION(:), INTENT(IN) :: mpi_recv_buf_tmp
      REAL(wp), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_data
      INTEGER, INTENT(IN) :: dimen1
      INTEGER :: i, j, k, current_p, ij, kcount

      !$acc data present(loc_nlpb,mpi_recv_buf_tmp,mpi_recv_data,mpi_recv_indexes_1d_displs,mpi_recv_indexes_1d_counts)
      current_p = 1
      DO i = 1, mpi_graph_indegree
         !$acc parallel copyin(current_p,i)
         !$acc loop collapse(2) private(kcount,ij) independent
         DO k = 1, dimen1
           DO j = 1, mpi_recv_indexes_1d_counts(i)
              ij = loc_nlpb + mpi_recv_indexes_1d_displs(i) + j
              kcount = current_p+(k-1)*mpi_recv_indexes_1d_counts(i)+(j-1)
              mpi_recv_data(ij, k) = mpi_recv_buf_tmp(kcount)
            END DO
         END DO
         !$acc end parallel
         current_p = current_p + mpi_recv_indexes_1d_counts(i)*dimen1
      END DO
      !$acc end data
   END SUBROUTINE

   ! adjust receiving data of i2 in the right order in 2 dimension
   SUBROUTINE mpi_i2_2d_prepare_recvbuf(mpi_recv_buf_tmp, mpi_recv_data, dimen1)
      INTEGER(i2), DIMENSION(:), INTENT(IN) :: mpi_recv_buf_tmp
      INTEGER(i2), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_data
      INTEGER, INTENT(IN) :: dimen1
      INTEGER :: i, j, k, current_p, ij, kcount

      !$acc data present(loc_nlpb,mpi_recv_buf_tmp,mpi_recv_data,mpi_recv_indexes_1d_displs,mpi_recv_indexes_1d_counts)
      current_p = 1
      DO i = 1, mpi_graph_indegree
         !$acc parallel copyin(current_p,i)
         !$acc loop collapse(2) private(kcount,ij) independent
         DO k = 1, dimen1
           DO j = 1, mpi_recv_indexes_1d_counts(i)
              ij = loc_nlpb + mpi_recv_indexes_1d_displs(i) + j
              kcount = current_p+(k-1)*mpi_recv_indexes_1d_counts(i)+(j-1)
              mpi_recv_data(ij, k) = mpi_recv_buf_tmp(kcount)
            END DO
         END DO
         !$acc end parallel
         current_p = current_p + mpi_recv_indexes_1d_counts(i)*dimen1
      END DO
      !$acc end data
   END SUBROUTINE

   ! adjust receiving data of dp in the right order in 3 dimension
   SUBROUTINE mpi_wp_3d_prepare_recvbuf(mpi_recv_buf_tmp, mpi_recv_data, dimen1, dimen2)
      REAL(wp), DIMENSION(:), INTENT(IN) :: mpi_recv_buf_tmp
      REAL(wp), DIMENSION(:, :, :), INTENT(INOUT) :: mpi_recv_data
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, INTENT(IN) :: dimen2
      INTEGER :: i, j, k, l, current_p, ij, kcount

      !$acc data present(loc_nlpb,mpi_recv_buf_tmp,mpi_recv_data,mpi_recv_indexes_1d_displs,mpi_recv_indexes_1d_counts)
      current_p = 1
      DO i = 1, mpi_graph_indegree
         !$acc parallel copyin(i,current_p)
         !$acc loop collapse(3) private(kcount,ij) independent
         DO j = 1, dimen2
            DO k = 1, dimen1
               DO l = 1, mpi_recv_indexes_1d_counts(i)
                  ij = loc_nlpb + mpi_recv_indexes_1d_displs(i) + l
                  kcount = current_p + (j-1)*dimen1*mpi_recv_indexes_1d_counts(i) &
                         + (k-1)*mpi_recv_indexes_1d_counts(i) + l-1
                  mpi_recv_data(ij, k, j) = mpi_recv_buf_tmp(kcount)
               END DO
            END DO
         END DO
         !$acc end parallel
         current_p = current_p + mpi_recv_indexes_1d_counts(i)*dimen1*dimen2
      END DO
      !$acc end data
   END SUBROUTINE

   ! blocking data exchange for wp of 1 dimension
   SUBROUTINE mpi_wp_1d_block_exchange(mpi_send_data, mpi_recv_data)
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_send_data
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_recv_data

      CALL mpi_wp_1d_prepare_sendbuf(mpi_recv_data, mpi_send_data)
      ! !$acc host_data use_device(mpi_send_data, mpi_recv_data)
      ! CALL MPI_NEIGHBOR_ALLTOALLV(mpi_send_data, mpi_send_indexes_1d_counts, mpi_send_indexes_1d_displs, mpi_real_wp, &
      !                             mpi_recv_data(loc_nlpb + 1), mpi_recv_indexes_1d_counts, mpi_recv_indexes_1d_displs, &
      !                             mpi_real_wp, mpi_graph_comm, mpi_err)
      ! !$acc end host_data

   END SUBROUTINE mpi_wp_1d_block_exchange

   ! blocking data exchange for i2 of 1 dimension
   SUBROUTINE mpi_i2_1d_block_exchange(mpi_send_data, mpi_recv_data)
      INTEGER(i2), DIMENSION(:), INTENT(INOUT) :: mpi_send_data
      INTEGER(i2), DIMENSION(:), INTENT(INOUT) :: mpi_recv_data
      CALL mpi_i2_1d_prepare_sendbuf(mpi_recv_data, mpi_send_data)
      !$acc host_data use_device(mpi_send_data,mpi_recv_data)
      ! CALL MPI_NEIGHBOR_ALLTOALLV(mpi_send_data, mpi_send_indexes_1d_counts, &
      !    mpi_send_indexes_1d_displs, MPI_I2, mpi_recv_data(loc_nlpb + 1), &
      !    mpi_recv_indexes_1d_counts, mpi_recv_indexes_1d_displs, MPI_I2, &
      !    mpi_graph_comm, mpi_err)
      ! !$acc end host_data

   END SUBROUTINE mpi_i2_1d_block_exchange

   ! blocking data exchange for wp of 3 dimension
   SUBROUTINE mpi_wp_3d_block_exchange(mpi_send_data, mpi_recv_data, dimen1, dimen2)
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_send_data
      REAL(wp), DIMENSION(:, :, :), INTENT(INOUT) :: mpi_recv_data
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, INTENT(IN) :: dimen2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_displs_tmp
      REAL(wp), ALLOCATABLE, DIMENSION(:) :: mpi_recv_buf_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_displs_tmp
      ALLOCATE (mpi_send_indexes_1d_counts_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_send_indexes_1d_displs_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_recv_buf_tmp(mpi_recv_indexes_1d_counts_sum*dimen1*dimen2))
      ALLOCATE (mpi_recv_indexes_counts_tmp(mpi_graph_indegree))
      ALLOCATE (mpi_recv_indexes_displs_tmp(mpi_graph_indegree))

      CALL mpi_dp_3d_prepare_sendbuf(mpi_recv_data, mpi_send_data, dimen1, dimen2, mpi_send_indexes_1d_counts_tmp, &
                                      mpi_send_indexes_1d_displs_tmp)
      CALL mpi_dp_3d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, dimen1, dimen2)
      !$acc data create(mpi_recv_buf_tmp),  &
      !$acc      present(mpi_send_data,mpi_recv_data)
      !$acc host_data use_device(mpi_send_data,mpi_recv_buf_tmp)
      ! CALL MPI_NEIGHBOR_ALLTOALLV(mpi_send_data, mpi_send_indexes_1d_counts_tmp, &
      !    mpi_send_indexes_1d_displs_tmp, mpi_real_wp, mpi_recv_buf_tmp, &
      !    mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, mpi_real_wp, &
      !    mpi_graph_comm, mpi_err)
      !$acc end host_data
      CALL mpi_wp_3d_prepare_recvbuf(mpi_recv_buf_tmp, mpi_recv_data, dimen1, dimen2)
      !$acc end data

      DEALLOCATE (mpi_send_indexes_1d_counts_tmp)
      DEALLOCATE (mpi_send_indexes_1d_displs_tmp)
      DEALLOCATE (mpi_recv_buf_tmp)
      DEALLOCATE (mpi_recv_indexes_counts_tmp)
      DEALLOCATE (mpi_recv_indexes_displs_tmp)

   END SUBROUTINE mpi_wp_3d_block_exchange

   ! non - blocking data exchange for wp of 3 dimension
   SUBROUTINE mpi_wp_3d_nonblock_exchange(mpi_source_data, mpi_send_data, mpi_recv_buf, dimen1, dimen2, mpi_req)
      REAL(wp), DIMENSION(:, :, :), INTENT(INOUT) :: mpi_source_data
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_send_data
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, INTENT(IN) :: dimen2
      INTEGER, INTENT(INOUT) :: mpi_req
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_displs_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_displs_tmp
      ALLOCATE (mpi_send_indexes_1d_counts_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_send_indexes_1d_displs_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_recv_indexes_counts_tmp(mpi_graph_indegree))
      ALLOCATE (mpi_recv_indexes_displs_tmp(mpi_graph_indegree))

      CALL mpi_dp_3d_prepare_sendbuf(mpi_source_data, mpi_send_data, dimen1, &
         dimen2, mpi_send_indexes_1d_counts_tmp, mpi_send_indexes_1d_displs_tmp)
      CALL mpi_dp_3d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, &
         mpi_recv_indexes_displs_tmp, dimen1, dimen2)
      ! !$acc data present(mpi_recv_buf,mpi_send_data)
      ! !$acc host_data use_device(mpi_send_data,mpi_recv_buf)
      ! CALL MPI_INEIGHBOR_ALLTOALLV(mpi_send_data, mpi_send_indexes_1d_counts_tmp, &
      !    mpi_send_indexes_1d_displs_tmp, mpi_real_wp, mpi_recv_buf, &
      !    mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, mpi_real_wp, &
      !    mpi_graph_comm, mpi_req, mpi_err)
      ! !$acc end host_data
      ! !$acc end data

      DEALLOCATE (mpi_send_indexes_1d_counts_tmp)
      DEALLOCATE (mpi_send_indexes_1d_displs_tmp)
      DEALLOCATE (mpi_recv_indexes_counts_tmp)
      DEALLOCATE (mpi_recv_indexes_displs_tmp)

   END SUBROUTINE mpi_wp_3d_nonblock_exchange

   ! check and wait for finishing non - blocking communication of 3 dimension, and adjust receiving data
   SUBROUTINE mpi_dp_3d_check_and_wait(mpi_recv_data, mpi_recv_buf, dimen1, dimen2, mpi_req)
      REAL(wp), DIMENSION(:, :, :), INTENT(INOUT) :: mpi_recv_data
      REAL(wp), DIMENSION(:), INTENT(IN) :: mpi_recv_buf
      INTEGER, INTENT(INOUT) :: mpi_req
      INTEGER, INTENT(IN) :: dimen1, dimen2
      CALL MPI_WAIT(mpi_req, mpi_stat, mpi_err)
      CALL mpi_wp_3d_prepare_recvbuf(mpi_recv_buf, mpi_recv_data, dimen1, dimen2)
   END SUBROUTINE mpi_dp_3d_check_and_wait

   ! blocking data exchange for wp of 2 dimension
   SUBROUTINE mpi_wp_2d_block_exchange(mpi_send_data, mpi_recv_data, dimen1)
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_send_data
      REAL(wp), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_data
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_displs_tmp
      REAL(wp), ALLOCATABLE, DIMENSION(:) :: mpi_recv_buf_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_displs_tmp
      ALLOCATE (mpi_send_indexes_1d_counts_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_send_indexes_1d_displs_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_recv_buf_tmp(mpi_recv_indexes_1d_counts_sum*dimen1))
      ALLOCATE (mpi_recv_indexes_counts_tmp(mpi_graph_indegree))
      ALLOCATE (mpi_recv_indexes_displs_tmp(mpi_graph_indegree))

      CALL mpi_wp_2d_prepare_sendbuf(mpi_recv_data, mpi_send_data, dimen1, &
         mpi_send_indexes_1d_counts_tmp, mpi_send_indexes_1d_displs_tmp)
      CALL mpi_2d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, &
         mpi_recv_indexes_displs_tmp, dimen1)
      !$acc data create(mpi_recv_buf_tmp),  &
      !$acc      present(mpi_send_data,mpi_recv_data)
      !$acc host_data use_device(mpi_send_data,mpi_recv_buf_tmp)
      ! CALL MPI_NEIGHBOR_ALLTOALLV(mpi_send_data, mpi_send_indexes_1d_counts_tmp, &
      !    mpi_send_indexes_1d_displs_tmp, mpi_real_wp, mpi_recv_buf_tmp, &
      !    mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, mpi_real_wp, &
      !    mpi_graph_comm, mpi_err)
      !$acc end host_data
      CALL mpi_wp_2d_prepare_recvbuf(mpi_recv_buf_tmp, mpi_recv_data, dimen1)
      !$acc end data

      DEALLOCATE (mpi_send_indexes_1d_counts_tmp)
      DEALLOCATE (mpi_send_indexes_1d_displs_tmp)
      DEALLOCATE (mpi_recv_buf_tmp)
      DEALLOCATE (mpi_recv_indexes_counts_tmp)
      DEALLOCATE (mpi_recv_indexes_displs_tmp)

   END SUBROUTINE mpi_wp_2d_block_exchange

   ! non - blocking data exchange for wp of 2 dimension
   SUBROUTINE mpi_wp_2d_nonblock_exchange(mpi_source_data, mpi_send_data, mpi_recv_buf, dimen1, mpi_req)
      REAL(wp), DIMENSION(:, :), INTENT(INOUT) :: mpi_source_data
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_send_data
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(INOUT) :: mpi_req
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_displs_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_displs_tmp
      ALLOCATE (mpi_send_indexes_1d_counts_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_send_indexes_1d_displs_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_recv_indexes_counts_tmp(mpi_graph_indegree))
      ALLOCATE (mpi_recv_indexes_displs_tmp(mpi_graph_indegree))

      CALL mpi_wp_2d_prepare_sendbuf(mpi_source_data, mpi_send_data, dimen1, &
         mpi_send_indexes_1d_counts_tmp, mpi_send_indexes_1d_displs_tmp)
      CALL mpi_2d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, dimen1)
      ! !$acc data present(mpi_send_data,mpi_recv_buf)
      ! !$acc host_data use_device(mpi_send_data,mpi_recv_buf)
      ! CALL MPI_INEIGHBOR_ALLTOALLV(mpi_send_data, mpi_send_indexes_1d_counts_tmp, &
      !    mpi_send_indexes_1d_displs_tmp, mpi_real_wp, mpi_recv_buf, &
      !    mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, mpi_real_wp, &
      !    mpi_graph_comm, mpi_req, mpi_err)
      ! !$acc end host_data
      ! !$acc end data

      DEALLOCATE (mpi_send_indexes_1d_counts_tmp)
      DEALLOCATE (mpi_send_indexes_1d_displs_tmp)
      DEALLOCATE (mpi_recv_indexes_counts_tmp)
      DEALLOCATE (mpi_recv_indexes_displs_tmp)

   END SUBROUTINE mpi_wp_2d_nonblock_exchange

   ! check and wait for finishing non - blocking communication of 2 dimension, and adjust receiving data
   SUBROUTINE mpi_dp_2d_check_and_wait(mpi_recv_data, mpi_recv_buf, dimen1, mpi_req)
      ! mpi_source_data is final position for writing data
      REAL(wp), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_data
      REAL(wp), DIMENSION(:), INTENT(IN) :: mpi_recv_buf
      INTEGER, INTENT(INOUT) :: mpi_req
      INTEGER, INTENT(IN) :: dimen1
      CALL MPI_WAIT(mpi_req, mpi_stat, mpi_err)
      CALL mpi_wp_2d_prepare_recvbuf(mpi_recv_buf, mpi_recv_data, dimen1)
   END SUBROUTINE mpi_dp_2d_check_and_wait

   ! blocking data exchange for i2 of 2 dimension
   SUBROUTINE mpi_i2_2d_block_exchange(mpi_send_data, mpi_recv_data, dimen1)
      INTEGER(i2), DIMENSION(:), INTENT(INOUT) :: mpi_send_data
      INTEGER(i2), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_data
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_displs_tmp
      INTEGER(i2), ALLOCATABLE, DIMENSION(:) :: mpi_recv_buf_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_displs_tmp
      ALLOCATE (mpi_send_indexes_1d_counts_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_send_indexes_1d_displs_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_recv_buf_tmp(mpi_recv_indexes_1d_counts_sum * dimen1))
      ALLOCATE (mpi_recv_indexes_counts_tmp(mpi_graph_indegree))
      ALLOCATE (mpi_recv_indexes_displs_tmp(mpi_graph_indegree))

      CALL mpi_i2_2d_prepare_sendbuf(mpi_recv_data, mpi_send_data, dimen1, &
         mpi_send_indexes_1d_counts_tmp, mpi_send_indexes_1d_displs_tmp)
      CALL mpi_2d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, &
         mpi_recv_indexes_displs_tmp, dimen1)
      !$acc data create(mpi_recv_buf_tmp),  &
      !$acc      present(mpi_send_data,mpi_recv_data)
      !$acc host_data use_device(mpi_send_data,mpi_recv_buf_tmp)
      ! CALL MPI_NEIGHBOR_ALLTOALLV(mpi_send_data, mpi_send_indexes_1d_counts_tmp, &
      !    mpi_send_indexes_1d_displs_tmp, MPI_I2, mpi_recv_buf_tmp, &
      !    mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, MPI_I2, mpi_graph_comm, mpi_err)
      !$acc end host_data
      CALL mpi_i2_2d_prepare_recvbuf(mpi_recv_buf_tmp, mpi_recv_data, dimen1)
      !$acc end data

      DEALLOCATE (mpi_send_indexes_1d_counts_tmp)
      DEALLOCATE (mpi_send_indexes_1d_displs_tmp)
      DEALLOCATE (mpi_recv_buf_tmp)
      DEALLOCATE (mpi_recv_indexes_counts_tmp)
      DEALLOCATE (mpi_recv_indexes_displs_tmp)

   END SUBROUTINE mpi_i2_2d_block_exchange

   ! non - blocking data exchange for wp of 1 dimension
   SUBROUTINE mpi_wp_1d_nonblock_exchange(mpi_source_data, mpi_send_data, mpi_recv_buf, mpi_req)

      REAL(wp), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: mpi_source_data
      REAL(wp), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: mpi_send_data
      REAL(wp), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(INOUT) :: mpi_req
      CALL mpi_wp_1d_prepare_sendbuf(mpi_source_data, mpi_send_data)
      ! !$acc data present(mpi_send_data,mpi_recv_buf)
      ! !$acc host_data use_device(mpi_send_data,mpi_recv_buf)
      ! CALL MPI_INEIGHBOR_ALLTOALLV(mpi_send_data, mpi_send_indexes_1d_counts, &
      !    mpi_send_indexes_1d_displs, mpi_real_wp, mpi_recv_buf, &
      !    mpi_recv_indexes_1d_counts, mpi_recv_indexes_1d_displs, mpi_real_wp, &
      !    mpi_graph_comm, mpi_req, mpi_err)
      ! !$acc end host_data
      ! !$acc end data

   END SUBROUTINE mpi_wp_1d_nonblock_exchange

   ! check and wait for finishing non - blocking communication of 1 dimension, and copy receiving data
   SUBROUTINE mpi_dp_1d_check_and_wait(mpi_recv_data, mpi_recv_buf, mpi_req)
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_recv_data
      REAL(wp), DIMENSION(:), INTENT(IN) :: mpi_recv_buf
      INTEGER, INTENT(INOUT) :: mpi_req
      INTEGER :: i

      CALL MPI_WAIT(mpi_req, mpi_stat, mpi_err)
      !$acc kernels present(loc_nlpb,mpi_recv_data,mpi_recv_buf)
      !$acc loop independent
      DO i = 1,mpi_recv_indexes_1d_counts_sum
         mpi_recv_data(loc_nlpb + i) = mpi_recv_buf(i)
      END DO
      !$acc end kernels
   END SUBROUTINE mpi_dp_1d_check_and_wait

   ! Inintialize mpi operations
   SUBROUTINE mpi_process_init
      CHARACTER(LEN = 16) :: arg
      CHARACTER(LEN = 16) :: env_value
      INTEGER :: i
      INTEGER :: ranks1(1), ranks2(1)
# if defined (OPENACC)
      !Multiple GPUs
      INTEGER :: InNodeComm, InNodeRank
      CHARACTER (len=MPI_MAX_PROCESSOR_NAME) :: hostname
      INTEGER :: namelength
# endif

      ! 首先尝试从环境变量获取I/O处理器数量
      CALL GET_ENVIRONMENT_VARIABLE("MACOM_IO_PROCS", env_value)
      IF (LEN_TRIM(env_value) > 0) THEN
         ! 环境变量存在，使用它的值
         READ(env_value, *) mpi_io_procs
      ELSE
         ! 环境变量不存在，尝试从命令行获取
         CALL GET_COMMAND_ARGUMENT(1, arg)

         IF (LEN_TRIM(arg) == 0) THEN
            WRITE (*, *) "please add the number (between 1 and 128) of I / O processors following the executable file"
            STOP
         END IF

         READ (arg, *) mpi_io_procs
      END IF

      IF (mpi_io_procs < 0 .OR. mpi_io_procs > 128) THEN
         WRITE (*, *) "please change the number of I/O processors to the &
            arrange between 1 and 128 following the executable file"
         STOP
      END IF

      CALL MPI_INIT(mpi_err)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_procs, mpi_err)

      IF (mpi_procs < mpi_io_procs + 1) THEN
         WRITE (*, *) "the number or computing processors should be more than 2,", &
            "please decrease the number of I/O processors or increase the number of all processors"
         STOP
      END IF

      mpi_comp_procs = mpi_procs - mpi_io_procs
      mpi_comp_io_procs = mpi_comp_procs + 1
      ALLOCATE (mpi_comp_all_ranks(mpi_comp_procs))
      ALLOCATE (mpi_io_all_ranks(mpi_io_procs))
      ALLOCATE (mpi_comp_io_mixed_ranks(mpi_comp_io_procs))

      mpi_comp_all_ranks = [(i, i = 0, mpi_comp_procs - 1, 1)]
      mpi_comp_io_mixed_ranks = [(i, i = 0, mpi_comp_io_procs - 1, 1)]
      mpi_io_all_ranks = [(i, i = mpi_comp_procs, mpi_procs - 1, 1)]

      CALL MPI_COMM_GROUP(MPI_COMM_WORLD, mpi_group_all, mpi_err)
      CALL MPI_GROUP_INCL(mpi_group_all, mpi_comp_procs, mpi_comp_all_ranks, mpi_group_comp, mpi_err)
      CALL MPI_GROUP_INCL(mpi_group_all, mpi_comp_io_procs, mpi_comp_io_mixed_ranks, mpi_group_comp_io, mpi_err)
      CALL MPI_GROUP_INCL(mpi_group_all, mpi_io_procs, mpi_io_all_ranks, mpi_group_io, mpi_err)

      CALL MPI_COMM_CREATE(MPI_COMM_WORLD, mpi_group_comp, mpi_comp_comm, mpi_err)
      CALL MPI_COMM_CREATE(MPI_COMM_WORLD, mpi_group_comp_io, mpi_comp_io_comm, mpi_err)
      CALL MPI_COMM_CREATE(MPI_COMM_WORLD, mpi_group_io, mpi_io_comm, mpi_err)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_err)

      ! initialize value, - 1 means null and not valid process rank
      mpi_comp_rank = - 1
      mpi_comp_io_rank = - 1
      mpi_io_rank = - 1

      ! you have to judge whether the process is in the communicator of mpi_comp_comm, otherwise there is an error of "Invalid communicator"
      IF (mpi_rank < mpi_comp_procs) THEN
         CALL MPI_COMM_RANK(mpi_comp_comm, mpi_comp_rank, mpi_err)
      END IF
      IF (mpi_rank < mpi_comp_io_procs) THEN
         CALL MPI_COMM_RANK(mpi_comp_io_comm, mpi_comp_io_rank, mpi_err)
      END IF
      IF (mpi_rank >= mpi_comp_procs .AND. mpi_rank < mpi_procs) THEN
         CALL MPI_COMM_RANK(mpi_io_comm, mpi_io_rank, mpi_err)
      END IF

      ! the process ID (= mpi_comp_procs) is root process in mpi_comp_comm, but we don't know root rank id in mpi_comp_comm
      ! through this method to get the root id in mpi_comp_comm
      ranks1(1) = 0
      mpi_root_comp = - 1
      ranks2(1) = mpi_root_comp
      CALL MPI_GROUP_TRANSLATE_RANKS(mpi_group_all, 1, ranks1, mpi_group_comp, ranks2, mpi_err)
      mpi_root_comp = ranks2(1)

      ! the process ID (= mpi_comp_procs) is root process in mpi_comp_io_comm from I / O communicator, but we don't know root rank id in mpi_comp_io_comm
      ! through this method to get the root id in mpi_comp_io_comm
      ranks1(1) = mpi_comp_procs
      mpi_root_comp_io = - 1
      ranks2(1) = mpi_root_comp_io
      CALL MPI_GROUP_TRANSLATE_RANKS(mpi_group_all, 1, ranks1, mpi_group_comp_io, ranks2, mpi_err)
      mpi_root_comp_io = ranks2(1)

      ! the process ID (= mpi_comp_procs) is root process in mpi_comp_io_comm from computational communicator, but we don't know root rank id in mpi_comp_io_comm
      ! through this method to get the root id in mpi_comp_io_comm
      ranks1(1) = 0
      mpi_sub_root_comp_io = - 1
      ranks2(1) = mpi_sub_root_comp_io
      CALL MPI_GROUP_TRANSLATE_RANKS(mpi_group_all, 1, ranks1, mpi_group_comp_io, ranks2, mpi_err)
      mpi_sub_root_comp_io = ranks2(1)

      ! the process ID (= mpi_comp_procs) is root process in mpi_comp_io_comm, but we don't know root rank id in mpi_io_comm
      ! through this method to get the root id in mpi_io_comm
      ranks1(1) = mpi_comp_procs
      mpi_root_io = - 1
      ranks2(1) = mpi_root_io
      CALL MPI_GROUP_TRANSLATE_RANKS(mpi_group_all, 1, ranks1, mpi_group_io, ranks2, mpi_err)
      mpi_root_io = ranks2(1)

      DEALLOCATE (mpi_comp_all_ranks)
      DEALLOCATE (mpi_io_all_ranks)
      DEALLOCATE (mpi_comp_io_mixed_ranks)

! -----------------------------------------------------------------------------------
# if defined (OPENACCGPU)
! Multiple GPUs
! FIND GPU DEVICES IN EACH NODE AND BIND GPUS WITH LOCAL RANKS (IN EACH NODE)
! THIS IS DONE BY USING MPI_COMM_SPLIT_TYPE WITH ATTRIBUTE MPI_COMM_TYPE_SHARED
! -----------
      ! ngpus per node
      ngpus = ACC_GET_NUM_DEVICES(ACC_DEVICE_NVIDIA)
      call MPI_GET_PROCESSOR_NAME(hostname, namelength, mpi_err)
      if (ngpus .le. 0) then
         if (mpi_rank .eq. 0) then
            write(*,'(a)'), '***GPU*** No NVIDIA GPUs available. STOP!'
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_err)
         end if
      else
         call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,   &
              MPI_INFO_NULL, InNodeComm, mpi_err)
         call MPI_COMM_RANK(InNodeComm, InNodeRank, mpi_err)
         if (mpi_rank .ge. mpi_comp_procs) then
            write(*, "(A25,A20)") "***GPU*** I/O processors at Node:", hostname(1:namelength)
         else
            if (InNodeRank .ge. ngpus) then
               write (*,*) "***GPU*** WARNING: CPU processors: ", InNodeRank, &
                  "more than device number: ", ngpus
            else
               call ACC_SET_DEVICE_NUM(InNodeRank, ACC_DEVICE_NVIDIA)
               write(*,"(A18,I4,A13,I2,A8,A20)") "***GPU*** RANK: ", mpi_rank, &
                  " USING GPU: ", InNodeRank, " AT: ", hostname(1:namelength)
            end if
         end if
      end if
      call MPI_COMM_FREE(InNodeComm, mpi_err)
# endif
! -----------------------------------------------------------------------------------

      IF (mpi_rank == 0) THEN
         WRITE (*, *) "Finish initializing mpi operations, the number of total processes is ", mpi_procs, &
            ", the number of computational processes is ", mpi_comp_procs, &
            ", the number of I / O processes is", mpi_io_procs
      END IF

      CALL MPI_TYPE_CREATE_F90_INTEGER(4, MPI_I2, mpi_err)

      if (wp == dp) mpi_real_wp = MPI_DOUBLE_PRECISION
      if (wp == sp) mpi_real_wp = MPI_REAL

   END SUBROUTINE

   ! reading and exchanging data of wp for 1 dimension from netcdf input data
   SUBROUTINE mpi_netcdf_read_exchange_1d_wp(filename, varname, var, options,optional_ts)
     CHARACTER(*), INTENT(IN) :: filename, varname
     REAL(wp), DIMENSION(:), INTENT(INOUT) :: var
     INTEGER, INTENT(IN) :: options
     INTEGER, INTENT(IN),OPTIONAL :: optional_ts
     ! 1 for nlpb,  2 for nlpbz
     SELECT CASE (options)
     CASE (1)

       IF(mpi_rank == 0) THEN
         IF(PRESENT(optional_ts)) THEN
           CALL netcdf_read(filename, varname, mpi_buf_1d_nlpb_wp, d1 = mpi_total_nlpb,ts = optional_ts)
         ELSE
           CALL netcdf_read(filename, varname, mpi_buf_1d_nlpb_wp, d1 = mpi_total_nlpb)
         END IF
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF
     CASE (2)
       IF(mpi_rank == 0) THEN
         CALL netcdf_read(filename, varname, mpi_buf_1d_nlpb_nlpbz_wp, d1 = mpi_total_nlpbz)
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF

     CASE (5)
       IF(mpi_rank == 0) THEN
         IF(PRESENT(optional_ts)) THEN
           CALL netcdf_read(filename, varname, mpi_buf_1d_nlbdy_wp, d1 = mpi_total_nlbdy,ts = optional_ts)
         ELSE
           CALL netcdf_read(filename, varname, mpi_buf_1d_nlbdy_wp, d1 = mpi_total_nlbdy)
         END IF
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF

     CASE DEFAULT
        WRITE (*, *) "wrong options for netcdf reading"
        STOP
     END SELECT
   END SUBROUTINE

   ! reading and exchanging data of integer for 1 dimension from netcdf input data
   SUBROUTINE mpi_netcdf_read_exchange_1d_integer(filename, varname, var, options,optional_ts)
     CHARACTER(*), INTENT(IN) :: filename, varname
     INTEGER, DIMENSION(:), INTENT(INOUT) :: var
     INTEGER, INTENT(IN) :: options
     INTEGER, INTENT(IN),OPTIONAL :: optional_ts
     ! 1 for nlpb,  2 for nlpbz,5 for nlbdy
     SELECT CASE (options)
     CASE (1)
       IF(mpi_rank == 0) THEN
         CALL netcdf_read(filename, varname, mpi_buf_1d_nlpb_integer, d1 = mpi_total_nlpb)
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF
     CASE (2)
       IF(mpi_rank == 0) THEN
         CALL netcdf_read(filename, varname, mpi_buf_1d_nlpb_nlpbz_integer, d1 = mpi_total_nlpbz)
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF
     CASE (5)
       IF(mpi_rank == 0) THEN
         IF(PRESENT(optional_ts)) THEN
           CALL netcdf_read(filename, varname, mpi_buf_1d_nlbdy_integer, d1 = mpi_total_nlbdy,ts = optional_ts)
         ELSE
           CALL netcdf_read(filename, varname, mpi_buf_1d_nlbdy_integer, d1 = mpi_total_nlbdy)
         END IF
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF
     CASE DEFAULT
        WRITE (*, *) "wrong options for netcdf reading"
        STOP
     END SELECT
   END SUBROUTINE

   ! reading and exchanging data of integer for 2 dimension from netcdf input data
   SUBROUTINE mpi_netcdf_read_exchange_2d_integer(filename, varname, var, options)
     CHARACTER(*), INTENT(IN) :: filename, varname
     INTEGER, DIMENSION(:,:), INTENT(INOUT) :: var
     INTEGER, INTENT(IN) :: options
     ! 1 for nlpb,  2 for nlpbz
     SELECT CASE (options)
     CASE (1)
       IF(mpi_rank == 0) THEN
         CALL netcdf_read(filename, varname, mpi_buf_2d_nlpb_ni_integer, d1 = mpi_total_nlpb, d2 = ni)
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF
     CASE (2)
       IF(mpi_rank == 0) THEN
         CALL netcdf_read(filename, varname, mpi_buf_2d_nlpb_nlpbz_ni_integer, d1 = mpi_total_nlpbz, d2 = ni)
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF
     CASE DEFAULT
        WRITE (*, *) "wrong options for netcdf reading"
        STOP
     END SELECT
   END SUBROUTINE

   ! reading and exchanging data of wp for 2 dimension from netcdf input data
   SUBROUTINE mpi_netcdf_read_exchange_2d_wp(filename, varname, var, options, optional_dim1, optional_ts)
     CHARACTER(*), INTENT(IN) :: filename, varname
     REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: var
     INTEGER, INTENT(IN) :: options
     INTEGER, INTENT(IN),OPTIONAL :: optional_dim1
     INTEGER, INTENT(IN),OPTIONAL :: optional_ts
     SELECT CASE (options)
     CASE (3)
       IF(mpi_rank == 0) THEN
         IF(PRESENT(optional_ts)) THEN
           CALL netcdf_read(filename, varname, mpi_buf_2d_nlpb_nk_wp, d1 = mpi_total_nlpb, d2 = nk, ts = optional_ts)
         ELSE
           CALL netcdf_read(filename, varname, mpi_buf_2d_nlpb_nk_wp, d1 = mpi_total_nlpb, d2 = nk)
         END IF
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF
!     CASE (4)
!       IF(mpi_rank == 0) THEN
!         CALL netcdf_read(filename, varname, mpi_buf_2d_nlpb_nlpbz_nk_wp, d1 = mpi_total_nlpbz, d2 = nk)
!         CALL mpi_netcdf_root_exchange(var, options)
!       ELSE
!         CALL mpi_netcdf_nonroot_exchange(var, options)
!       END IF
     CASE (5)
       IF(mpi_rank == 0) THEN
         IF(PRESENT(optional_ts)) THEN
           CALL netcdf_read(filename, varname, mpi_buf_2d_nlpb_sbct_wp, d1 = mpi_total_nlpb, d2 = sbct, ts = optional_ts)
         ELSE
           CALL netcdf_read(filename, varname, mpi_buf_2d_nlpb_sbct_wp, d1 = mpi_total_nlpb, d2 = sbct)
         END IF
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF

     CASE (6)
       IF(mpi_rank == 0) THEN
         IF(PRESENT(optional_ts)) THEN
           CALL netcdf_read(filename, varname, mpi_buf_2d_nlbdy_nk_wp, d1 = mpi_total_nlbdy, d2 = nk, ts = optional_ts)
         ELSE
           CALL netcdf_read(filename, varname, mpi_buf_2d_nlbdy_nk_wp, d1 = mpi_total_nlbdy, d2 = nk)
         END IF
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF

     CASE (8)
       IF(.NOT. PRESENT(optional_dim1)) THEN
         WRITE(*,*) "please specify the length of the second dimension to be read and exchange"
       END IF
       IF(mpi_rank == 0) THEN
         ALLOCATE(mpi_buf_2d_wp(mpi_total_nlbdy,optional_dim1))
         IF(PRESENT(optional_ts)) THEN
           CALL netcdf_read(filename, varname, mpi_buf_2d_wp, d1 = mpi_total_nlbdy, d2 = optional_dim1, ts = optional_ts)
         ELSE
           CALL netcdf_read(filename, varname, mpi_buf_2d_wp, d1 = mpi_total_nlbdy, d2 = optional_dim1)
         END IF
         CALL mpi_netcdf_root_exchange(var, options,optional_dim1)
         DEALLOCATE(mpi_buf_2d_wp)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options,optional_dim1)
       END IF

     CASE (9)
       IF(.NOT. PRESENT(optional_dim1)) THEN
         WRITE(*,*) "please specify the length of the second dimension to be read and exchange"
       END IF
       IF(mpi_rank == 0) THEN
         ALLOCATE(mpi_buf_2d_wp(mpi_total_nlpb,optional_dim1))
         IF(PRESENT(optional_ts)) THEN
           CALL netcdf_read(filename, varname, mpi_buf_2d_wp, d1 = mpi_total_nlpb, d2 = optional_dim1, ts = optional_ts)
         ELSE
           CALL netcdf_read(filename, varname, mpi_buf_2d_wp, d1 = mpi_total_nlpb, d2 = optional_dim1)
         END IF
         CALL mpi_netcdf_root_exchange(var, options,optional_dim1)
         DEALLOCATE(mpi_buf_2d_wp)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options,optional_dim1)
       END IF

     CASE DEFAULT
        WRITE (*, *) "wrong options for netcdf reading"
        STOP
     END SELECT
   END SUBROUTINE

   ! reading and exchanging data of integer for 2 dimension from netcdf input data
   SUBROUTINE mpi_netcdf_read_exchange_2d_i2(filename, varname, var, options)
     CHARACTER(*), INTENT(IN) :: filename, varname
     INTEGER(i2), DIMENSION(:,:), INTENT(INOUT) :: var
     INTEGER, INTENT(IN) :: options
     SELECT CASE (options)
     CASE (3)
       IF(mpi_rank == 0) THEN
         CALL netcdf_read(filename, varname, mpi_buf_2d_nlpb_nk_i2, d1 = mpi_total_nlpb, d2 = nk)
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF
     CASE (4)
       IF(mpi_rank == 0) THEN
         CALL netcdf_read(filename, varname, mpi_buf_2d_nlpb_nlpbz_nk_i2, d1 = mpi_total_nlpbz, d2 = nk)
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF
     CASE DEFAULT
        WRITE (*, *) "wrong options for netcdf reading"
        STOP
     END SELECT
   END SUBROUTINE

   ! reading and exchanging data of i2 for 1 dimension from netcdf input data
   SUBROUTINE mpi_netcdf_read_exchange_1d_i2(filename, varname, var, options)
     CHARACTER(*), INTENT(IN) :: filename, varname
     INTEGER(i2), DIMENSION(:), INTENT(INOUT) :: var
     INTEGER, INTENT(IN) :: options
     ! 1 for nlpb,  2 for nlpbz
     SELECT CASE (options)
     CASE (1)
       IF(mpi_rank == 0) THEN
         CALL netcdf_read(filename, varname, mpi_buf_1d_nlpb_i2, d1 = mpi_total_nlpb)
         CALL mpi_netcdf_root_exchange(var, options)
       ELSE
         CALL mpi_netcdf_nonroot_exchange(var, options)
       END IF
!     CASE (2)
!       IF(mpi_rank == 0) THEN
!         CALL netcdf_read(filename, varname, mpi_buf_1d_nlpb_nlpbz_integer, d1 = mpi_total_nlpbz)
!         CALL mpi_netcdf_root_exchange(var, options)
!       ELSE
!         CALL mpi_netcdf_nonroot_exchange(var, options)
!       END IF
     CASE DEFAULT
        WRITE (*, *) "wrong options for netcdf reading"
        STOP
     END SELECT
   END SUBROUTINE

   ! sending data of wp for 1 dimension from netcdf input data for root process
   SUBROUTINE mpi_netcdf_root_exchange_1d_wp(mpi_recv_buf, mpi_options)
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(IN) :: mpi_options
      SELECT CASE (mpi_options)
      CASE (1)
         CALL mpi_1d_io_read_adjust_wp(mpi_buf_1d_nlpb_wp, mpi_adjust_buf_1d_nlpb_shadow_wp, mpi_options)
         CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlpb_shadow_wp,mpi_nlpb_shadow_counts_all, mpi_nlpb_shadow_displs_all, &
                           mpi_real_wp, mpi_recv_buf, nlpb, mpi_real_wp, 0, mpi_comp_comm, mpi_err)
      CASE (2)
         CALL mpi_1d_io_read_adjust_wp(mpi_buf_1d_nlpb_nlpbz_wp, mpi_adjust_buf_1d_nlpb_shadow_nlpbz_wp, mpi_options)
         CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlpb_shadow_nlpbz_wp, mpi_nlpb_shadow_nlpbz_counts_all, &
            mpi_nlpb_shadow_nlpbz_displs_all, mpi_real_wp, mpi_recv_buf, nlpbz, &
            mpi_real_wp, 0, mpi_comp_comm, mpi_err)
      CASE (5)
         CALL mpi_1d_io_read_adjust_wp(mpi_buf_1d_nlbdy_wp, mpi_adjust_buf_1d_nlbdy_wp, mpi_options)
         CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlbdy_wp,mpi_bdy_counts_1d_all,mpi_bdy_displs_1d_all, mpi_real_wp, &
                            mpi_recv_buf, nlbdy, mpi_real_wp, 0, mpi_bdy_comm, mpi_err)

      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE

   ! sending data of integer for 1 dimension from netcdf input data for root process
   SUBROUTINE mpi_netcdf_root_exchange_1d_integer(mpi_recv_buf, mpi_options)
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(IN) :: mpi_options
      SELECT CASE (mpi_options)
      CASE (1)
         CALL mpi_1d_io_read_adjust_integer(mpi_buf_1d_nlpb_integer, mpi_adjust_buf_1d_nlpb_shadow_integer, mpi_options)
         CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlpb_shadow_integer, mpi_nlpb_shadow_counts_all, mpi_nlpb_shadow_displs_all, &
                          MPI_INTEGER, mpi_recv_buf, nlpb, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
      CASE (2)
         CALL mpi_1d_io_read_adjust_integer(mpi_buf_1d_nlpb_nlpbz_integer, &
            mpi_adjust_buf_1d_nlpb_shadow_nlpbz_integer, mpi_options)
         CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlpb_shadow_nlpbz_integer, &
            mpi_nlpb_shadow_nlpbz_counts_all, mpi_nlpb_shadow_nlpbz_displs_all, &
            MPI_INTEGER, mpi_recv_buf, nlpbz, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

      CASE (5)
         CALL mpi_1d_io_read_adjust_integer(mpi_buf_1d_nlbdy_integer, mpi_adjust_buf_1d_nlbdy_integer, mpi_options)
         CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlbdy_integer,mpi_bdy_counts_1d_all,mpi_bdy_displs_1d_all, MPI_INTEGER, &
                            mpi_recv_buf, nlbdy, MPI_INTEGER, 0, mpi_bdy_comm, mpi_err)
      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE

   ! sending data of integer for 2 dimension from netcdf input data for root process
   SUBROUTINE mpi_netcdf_root_exchange_2d_integer(mpi_recv_buf, mpi_options)
      INTEGER, DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(IN) :: mpi_options
      SELECT CASE (mpi_options)
      CASE (1)
         CALL mpi_2d_io_read_adjust_integer(mpi_buf_2d_nlpb_ni_integer, &
            mpi_adjust_buf_2d_nlpb_shadow_ni_integer, mpi_options)
         CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_ni_integer, mpi_adjust_2d_nlpb_shadow_ni_counts_all, &
                           mpi_adjust_2d_nlpb_shadow_ni_displs_all, MPI_INTEGER, mpi_recv_buf, nlpb * ni, &
                           MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
      CASE (2)
         CALL mpi_2d_io_read_adjust_integer(mpi_buf_2d_nlpb_nlpbz_ni_integer, &
            mpi_adjust_buf_2d_nlpb_shadow_nlpbz_ni_integer, &
                                            mpi_options)
         CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_nlpbz_ni_integer, &
            mpi_adjust_2d_nlpb_shadow_nlpbz_ni_counts_all, &
            mpi_adjust_2d_nlpb_shadow_nlpbz_ni_displs_all, MPI_INTEGER, &
            mpi_recv_buf, nlpbz * ni, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE

   ! sending data of wp for 2 dimension from netcdf input data for root process
   SUBROUTINE mpi_netcdf_root_exchange_2d_wp(mpi_recv_buf, mpi_options,optional_dim1)
      REAL(wp), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(IN) :: mpi_options
      INTEGER, INTENT(IN),OPTIONAL :: optional_dim1
      REAL(wp), DIMENSION(:, :), ALLOCATABLE:: mpi_adjust_buf_2d
      SELECT CASE (mpi_options)
      CASE (3)
         CALL mpi_2d_io_read_adjust_wp(mpi_buf_2d_nlpb_nk_wp, &
                                       mpi_adjust_buf_2d_nlpb_shadow_nk_wp, mpi_options)

         CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_nk_wp, mpi_adjust_2d_nlpb_shadow_nk_counts_all, &
                           mpi_adjust_2d_nlpb_shadow_nk_displs_all, mpi_real_wp, mpi_recv_buf, nlpb * nk, &
                           mpi_real_wp, 0, mpi_comp_comm, mpi_err)

      CASE (5)
         CALL mpi_2d_io_read_adjust_wp(mpi_buf_2d_nlpb_sbct_wp, &
                                       mpi_adjust_buf_2d_nlpb_shadow_sbct_wp, mpi_options)

         CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_sbct_wp, mpi_adjust_2d_nlpb_shadow_sbct_counts_all, &
                           mpi_adjust_2d_nlpb_shadow_sbct_displs_all, mpi_real_wp, mpi_recv_buf, nlpb * sbct, &
                           mpi_real_wp, 0, mpi_comp_comm, mpi_err)
      CASE (6)
         CALL mpi_2d_io_read_adjust_wp(mpi_buf_2d_nlbdy_nk_wp, &
                                       mpi_adjust_buf_2d_nlbdy_nk_wp, mpi_options)

         ! communicator is mpi_bdy_comm
         CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlbdy_nk_wp, mpi_adjust_2d_nlbdy_nk_counts_all, &
                           mpi_adjust_2d_nlbdy_nk_displs_all, mpi_real_wp, mpi_recv_buf, nlbdy * nk, &
                           mpi_real_wp, 0, mpi_bdy_comm, mpi_err)

      CASE (8)
         ALLOCATE(mpi_buf_adjust_2d_wp(mpi_total_nlbdy_shadow * optional_dim1))
         ALLOCATE(mpi_adjust_2d_user_counts_all(mpi_comp_procs))
         ALLOCATE(mpi_adjust_2d_user_displs_all(mpi_comp_procs))
         mpi_adjust_2d_user_counts_all(1:mpi_comp_procs) = mpi_bdy_counts_1d_all(1:mpi_comp_procs)* optional_dim1
         mpi_adjust_2d_user_displs_all(1:mpi_comp_procs)= mpi_bdy_displs_1d_all(1:mpi_comp_procs)* optional_dim1
         CALL mpi_2d_io_read_adjust_wp(mpi_buf_2d_wp,mpi_buf_adjust_2d_wp, mpi_options,optional_dim1)

         CALL MPI_SCATTERV(mpi_buf_adjust_2d_wp, mpi_adjust_2d_user_counts_all, &
                           mpi_adjust_2d_user_displs_all, mpi_real_wp, mpi_recv_buf, nlbdy * optional_dim1, &
                           mpi_real_wp, 0, mpi_comp_comm, mpi_err)
         DEALLOCATE(mpi_buf_adjust_2d_wp)
         DEALLOCATE(mpi_adjust_2d_user_counts_all)
         DEALLOCATE(mpi_adjust_2d_user_displs_all)

      CASE (9)
         ALLOCATE(mpi_buf_adjust_2d_wp(mpi_total_nlpb_shadow * optional_dim1))
         ALLOCATE(mpi_adjust_2d_user_counts_all(mpi_comp_procs))
         ALLOCATE(mpi_adjust_2d_user_displs_all(mpi_comp_procs))
         mpi_adjust_2d_user_counts_all(1:mpi_comp_procs) = mpi_nlpb_shadow_counts_all(1:mpi_comp_procs)* optional_dim1
         mpi_adjust_2d_user_displs_all(1:mpi_comp_procs)= mpi_nlpb_shadow_displs_all(1:mpi_comp_procs)* optional_dim1
         CALL mpi_2d_io_read_adjust_wp(mpi_buf_2d_wp,mpi_buf_adjust_2d_wp, mpi_options,optional_dim1)

         CALL MPI_SCATTERV(mpi_buf_adjust_2d_wp, mpi_adjust_2d_user_counts_all, &
                           mpi_adjust_2d_user_displs_all, mpi_real_wp, mpi_recv_buf, nlpb * optional_dim1, &
                           mpi_real_wp, 0, mpi_comp_comm, mpi_err)
         DEALLOCATE(mpi_buf_adjust_2d_wp)
         DEALLOCATE(mpi_adjust_2d_user_counts_all)
         DEALLOCATE(mpi_adjust_2d_user_displs_all)
      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE mpi_netcdf_root_exchange_2d_wp

   ! reading and sending data of wp for 2 dimension from netcdf input data for root process
   SUBROUTINE mpi_netcdf_root_exchange_2d_i2(mpi_recv_buf, mpi_options)
      INTEGER(i2), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: mpi_recv_buf_integer
      INTEGER, INTENT(IN) :: mpi_options
      SELECT CASE (mpi_options)
      CASE (3)
         ALLOCATE (mpi_recv_buf_integer(nlpb, nk))
         CALL mpi_2d_io_read_adjust_i2(mpi_buf_2d_nlpb_nk_i2, mpi_adjust_buf_2d_nlpb_shadow_nk_integer, mpi_options)
         CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_nk_integer, &
            mpi_adjust_2d_nlpb_shadow_nk_counts_all, mpi_adjust_2d_nlpb_shadow_nk_displs_all, &
            MPI_INTEGER, mpi_recv_buf_integer, nlpb * nk, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         mpi_recv_buf(1:nlpb, 1:nk) = INT(mpi_recv_buf_integer(1:nlpb, 1:nk), KIND(mpi_recv_buf))
         DEALLOCATE (mpi_recv_buf_integer)

      CASE (4)
         ALLOCATE (mpi_recv_buf_integer(nlpbz, nk))
         CALL mpi_2d_io_read_adjust_i2(mpi_buf_2d_nlpb_nlpbz_nk_i2, &
            mpi_adjust_buf_2d_nlpb_shadow_nlpbz_nk_integer, mpi_options)
         CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_nlpbz_nk_integer, &
            mpi_adjust_2d_nlpb_shadow_nlpbz_nk_counts_all, &
            mpi_adjust_2d_nlpb_shadow_nlpbz_nk_displs_all, MPI_INTEGER, &
            mpi_recv_buf_integer, nlpbz * nk, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         mpi_recv_buf(1:nlpbz, 1:nk) = INT(mpi_recv_buf_integer(1:nlpbz, 1:nk), KIND(mpi_recv_buf))
         DEALLOCATE (mpi_recv_buf_integer)
      CASE DEFAULT
         STOP
         WRITE (*, *) "wrong options for netcdf reading"
      END SELECT

   END SUBROUTINE

   ! sending data of i2 for 1 dimension from netcdf input data for root process
   SUBROUTINE mpi_netcdf_root_exchange_1d_i2(mpi_recv_buf_i2, mpi_options)
      INTEGER(i2), DIMENSION(:), INTENT(INOUT) :: mpi_recv_buf_i2
      INTEGER, DIMENSION(:), ALLOCATABLE :: mpi_recv_buf_integer
      INTEGER, INTENT(IN) :: mpi_options
      SELECT CASE (mpi_options)
      CASE (1)
         CALL mpi_1d_io_read_adjust_i2(mpi_buf_1d_nlpb_i2, mpi_adjust_buf_1d_nlpb_shadow_integer, mpi_options)
!         ! pay attention size of adjust buffer is not mpi_total_nlpb, it should be bigger as mpi_total_nlpb_shadow
!         mpi_adjust_buf_1d_nlpb_shadow_integer(1:mpi_total_nlpb_shadow) = INT(mpi_adjust_buf_1d_nlpb_shadow_i2(1:mpi_total_nlpb_shadow), &
!                                                                        KIND(mpi_adjust_buf_1d_nlpb_shadow_integer))
         ALLOCATE (mpi_recv_buf_integer(nlpb))
         CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlpb_shadow_integer, &
            mpi_nlpb_shadow_counts_all, mpi_nlpb_shadow_displs_all, MPI_INTEGER, &
            mpi_recv_buf_integer, nlpb, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         mpi_recv_buf_i2(1:nlpb) = INT(mpi_recv_buf_integer(1:nlpb), KIND(mpi_recv_buf_i2))
         DEALLOCATE (mpi_recv_buf_integer)

      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE

   ! receiving data of wp for 1 dimension from netcdf input data for non - root processes
   SUBROUTINE mpi_netcdf_nonroot_exchange_1d_wp(mpi_recv_buf, mpi_options)
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(IN) :: mpi_options
      SELECT CASE (mpi_options)
      CASE (1)
        CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlpb_shadow_wp, &
           mpi_nlpb_shadow_counts_all, mpi_nlpb_shadow_displs_all, mpi_real_wp, &
           mpi_recv_buf, nlpb, mpi_real_wp, 0, mpi_comp_comm, mpi_err)
      CASE (2)
        CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlpb_shadow_nlpbz_wp, &
           mpi_nlpb_shadow_nlpbz_counts_all, mpi_nlpb_shadow_nlpbz_displs_all, &
           mpi_real_wp, mpi_recv_buf, nlpbz, mpi_real_wp, 0, mpi_comp_comm, mpi_err)
      CASE (5)
         CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlbdy_wp, mpi_bdy_counts_1d_all, &
            mpi_bdy_displs_1d_all, mpi_real_wp, mpi_recv_buf, nlbdy, &
            mpi_real_wp, 0, mpi_bdy_comm, mpi_err)
      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE

   ! receiving data of integer for 1 dimension from netcdf input data for non - root processes
   SUBROUTINE mpi_netcdf_nonroot_exchange_1d_integer(mpi_recv_buf, mpi_options)
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(IN) :: mpi_options
      SELECT CASE (mpi_options)
      CASE (1)
        CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlpb_shadow_integer, &
           mpi_nlpb_shadow_counts_all, mpi_nlpb_shadow_displs_all, MPI_INTEGER, &
           mpi_recv_buf, nlpb, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
      CASE (2)
        CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlpb_shadow_nlpbz_integer, &
           mpi_nlpb_shadow_nlpbz_counts_all, mpi_nlpb_shadow_nlpbz_displs_all, &
           MPI_INTEGER, mpi_recv_buf, nlpbz, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
      CASE (5)
         CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlbdy_integer, &
           mpi_bdy_counts_1d_all, mpi_bdy_displs_1d_all, MPI_INTEGER, &
           mpi_recv_buf, nlbdy, MPI_INTEGER, 0, mpi_bdy_comm, mpi_err)

      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE

   ! receiving data of integer for 2 dimension from netcdf input data for non - root processes
   SUBROUTINE mpi_netcdf_nonroot_exchange_2d_integer(mpi_recv_buf, mpi_options)
      INTEGER, DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(IN) :: mpi_options
      SELECT CASE (mpi_options)
      CASE (1)
       CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_ni_integer, &
          mpi_adjust_2d_nlpb_shadow_ni_counts_all, mpi_adjust_2d_nlpb_shadow_ni_displs_all, MPI_INTEGER, &
          mpi_recv_buf, nlpb * ni, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
      CASE (2)
       CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_nlpbz_ni_integer, &
          mpi_adjust_2d_nlpb_shadow_nlpbz_ni_counts_all, mpi_adjust_2d_nlpb_shadow_nlpbz_ni_displs_all, &
          MPI_INTEGER, mpi_recv_buf, nlpbz * ni, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE

   ! receiving data of wp for 2 dimension from netcdf input data for non - root processes
   SUBROUTINE mpi_netcdf_nonroot_exchange_2d_wp(mpi_recv_buf, mpi_options,optional_dim1)
      REAL(wp), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(IN) :: mpi_options
      INTEGER, INTENT(IN),OPTIONAL :: optional_dim1
      SELECT CASE (mpi_options)

      CASE (3)
       CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_nk_wp, mpi_adjust_2d_nlpb_shadow_nk_counts_all, &
                          mpi_adjust_2d_nlpb_shadow_nk_displs_all, mpi_real_wp, mpi_recv_buf, &
                          nlpb * nk, mpi_real_wp, 0, mpi_comp_comm, mpi_err)

      CASE (5)
       CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_sbct_wp, mpi_adjust_2d_nlpb_shadow_sbct_counts_all, &
                          mpi_adjust_2d_nlpb_shadow_sbct_displs_all, mpi_real_wp, mpi_recv_buf, nlpb * sbct, &
                          mpi_real_wp, 0, mpi_comp_comm, mpi_err)

      CASE (6)
       ! communicator is mpi_bdy_comm
       IF(mpi_check_bdy) THEN
         CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlbdy_nk_wp, mpi_adjust_2d_nlbdy_nk_counts_all, &
                           mpi_adjust_2d_nlbdy_nk_displs_all, mpi_real_wp, mpi_recv_buf, nlbdy * nk, &
                           mpi_real_wp, 0, mpi_bdy_comm, mpi_err)
       END IF

      CASE (8)
       ! communicator is mpi_bdy_comm
       IF(mpi_check_bdy) THEN
         CALL MPI_SCATTERV(mpi_buf_adjust_2d_wp, mpi_adjust_2d_user_counts_all, &
                           mpi_adjust_2d_user_displs_all, mpi_real_wp, mpi_recv_buf, nlbdy * optional_dim1, &
                           mpi_real_wp, 0, mpi_comp_comm, mpi_err)
       END IF

      CASE (9)
        CALL MPI_SCATTERV(mpi_buf_adjust_2d_wp, mpi_adjust_2d_user_counts_all, &
                           mpi_adjust_2d_user_displs_all, mpi_real_wp, mpi_recv_buf, nlpb * optional_dim1, &
                           mpi_real_wp, 0, mpi_comp_comm, mpi_err)
      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE

   ! receiving data of i2 for 2 dimension from netcdf input data for non - root processes
   SUBROUTINE mpi_netcdf_nonroot_exchange_2d_i2(mpi_recv_buf, mpi_options)
      INTEGER(i2), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_buf
      INTEGER, INTENT(IN) :: mpi_options
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: mpi_recv_buf_integer
      SELECT CASE (mpi_options)
      CASE (1)
      CASE (2)

      CASE (3)
         ALLOCATE (mpi_recv_buf_integer(nlpb, nk))
         CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_ni_integer, &
            mpi_adjust_2d_nlpb_shadow_nk_counts_all, mpi_adjust_2d_nlpb_shadow_nk_displs_all, &
            MPI_INTEGER, mpi_recv_buf_integer, nlpb * nk, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         mpi_recv_buf(1:nlpb, 1:nk) = INT(mpi_recv_buf_integer(1:nlpb, 1:nk), KIND(mpi_recv_buf))
         DEALLOCATE (mpi_recv_buf_integer)

      CASE (4)

         ALLOCATE (mpi_recv_buf_integer(nlpbz, nk))
         CALL MPI_SCATTERV(mpi_adjust_buf_2d_nlpb_shadow_nlpbz_nk_integer,mpi_adjust_2d_nlpb_shadow_nlpbz_nk_counts_all, &
                           mpi_adjust_2d_nlpb_shadow_nlpbz_nk_displs_all, MPI_INTEGER, &
                           mpi_recv_buf_integer, nlpbz * nk, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         mpi_recv_buf(1:nlpbz, 1:nk) = INT(mpi_recv_buf_integer(1:nlpbz, 1:nk), KIND(mpi_recv_buf))
         DEALLOCATE (mpi_recv_buf_integer)
      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE

   ! receiving data of i2 for 1 dimension from netcdf input data for non - root processes
   SUBROUTINE mpi_netcdf_nonroot_exchange_1d_i2(mpi_recv_buf_i2, mpi_options)
      INTEGER(i2), DIMENSION(:), INTENT(INOUT) :: mpi_recv_buf_i2
      INTEGER, DIMENSION(:), ALLOCATABLE :: mpi_recv_buf_integer
      INTEGER, INTENT(IN) :: mpi_options
      SELECT CASE (mpi_options)
      CASE (1)
         ALLOCATE (mpi_recv_buf_integer(nlpb))
         CALL MPI_SCATTERV(mpi_adjust_buf_1d_nlpb_shadow_integer, &
            mpi_nlpb_shadow_counts_all, mpi_nlpb_shadow_displs_all, MPI_INTEGER, &
            mpi_recv_buf_integer, nlpb, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         mpi_recv_buf_i2(1:nlpb) = INT(mpi_recv_buf_integer(1:nlpb), KIND(mpi_recv_buf_i2))
         DEALLOCATE (mpi_recv_buf_integer)


      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE

   ! adjust sending buffer of wp of 1 dimension for data exchange
   SUBROUTINE mpi_1d_io_read_adjust_wp(mpi_source_buf, mpi_adjust_buf, mpi_options)
      REAL(wp), DIMENSION(:), INTENT(IN) :: mpi_source_buf
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_adjust_buf
      INTEGER, INTENT(IN) :: mpi_options
      INTEGER ::  j, k, index_original_start, index_original_end, index_adjust_start, index_adjust_end

      SELECT CASE (mpi_options)
      CASE (1)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO k = 1, mpi_nlpb_shadow_block_1d_all(j), 1
               index_original_start = &
                  mpi_nlpb_shadow_block_starts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k)
               index_original_end = index_original_start + &
                  mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
               index_adjust_end = index_adjust_start + &
                  mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
               mpi_adjust_buf(index_adjust_start:index_adjust_end) = &
                  mpi_source_buf(index_original_start:index_original_end)
               index_adjust_start = index_adjust_end + 1
            END DO
         END DO

      CASE (2)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO k = 1, mpi_nlpb_shadow_nlpbz_block_1d_all(j), 1
               index_original_start = &
                  mpi_nlpb_shadow_nlpbz_block_starts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k)
               index_original_end = index_original_start + &
                  mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k) - 1
               index_adjust_end = index_adjust_start + &
                  mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k) - 1
               mpi_adjust_buf(index_adjust_start:index_adjust_end) = &
                  mpi_source_buf(index_original_start:index_original_end)
               index_adjust_start = index_adjust_end + 1
            END DO
         END DO

       CASE (5)
         index_adjust_start = 1
         DO j = 1, mpi_bdy_outdegree, 1
           DO k = 1, mpi_bdy_counts_1d_all(j), 1
              index_original_start = mpi_original_to_ordered_bdy_all(mpi_bdy_displs_1d_all(j)+ k)
              mpi_adjust_buf(index_adjust_start) = mpi_source_buf(index_original_start)
              index_adjust_start = index_adjust_start + 1
           END DO
         END DO


      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE mpi_1d_io_read_adjust_wp

   ! adjust sending buffer of integer of 1 dimension for data exchange
   SUBROUTINE mpi_1d_io_read_adjust_integer(mpi_source_buf, mpi_adjust_buf, mpi_options)
      INTEGER, DIMENSION(:), INTENT(IN) :: mpi_source_buf
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_adjust_buf
      INTEGER, INTENT(IN) :: mpi_options
      INTEGER ::  j, k, index_original_start, index_original_end, index_adjust_start, index_adjust_end

      SELECT CASE (mpi_options)
      CASE (1)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO k = 1, mpi_nlpb_shadow_block_1d_all(j), 1
               index_original_start = &
                  mpi_nlpb_shadow_block_starts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k)
               index_original_end = index_original_start + &
                  mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
               index_adjust_end = index_adjust_start + &
                  mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
               mpi_adjust_buf(index_adjust_start:index_adjust_end) = &
                  mpi_source_buf(index_original_start:index_original_end)
               index_adjust_start = index_adjust_end + 1
            END DO
         END DO

      CASE (2)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO k = 1, mpi_nlpb_shadow_nlpbz_block_1d_all(j), 1
               index_original_start = &
                  mpi_nlpb_shadow_nlpbz_block_starts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k)
               index_original_end = index_original_start + &
                  mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k) - 1
               index_adjust_end = index_adjust_start + &
                  mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k) - 1
               mpi_adjust_buf(index_adjust_start:index_adjust_end) = &
                  mpi_source_buf(index_original_start:index_original_end)
               index_adjust_start = index_adjust_end + 1
            END DO
         END DO

      CASE (5)
        index_adjust_start = 1
        DO j = 1, mpi_bdy_outdegree, 1
          DO k = 1, mpi_bdy_counts_1d_all(j), 1
             index_original_start = mpi_original_to_ordered_bdy_all(mpi_bdy_displs_1d_all(j)+ k)
             mpi_adjust_buf(index_adjust_start) = mpi_source_buf(index_original_start)
             index_adjust_start = index_adjust_start + 1
          END DO
        END DO

      CASE DEFAULT
         WRITE (*, *) "wrong options for netcdf reading"
         STOP
      END SELECT

   END SUBROUTINE mpi_1d_io_read_adjust_integer

   ! adjust sending buffer of integer of 2 dimension for data exchange
   SUBROUTINE mpi_2d_io_read_adjust_integer(mpi_source_buf, mpi_adjust_buf, mpi_options)
      INTEGER, DIMENSION(:, :), INTENT(IN) :: mpi_source_buf
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_adjust_buf
      INTEGER, INTENT(IN) :: mpi_options
      INTEGER ::  i, j, k, index_original_start, index_original_end, index_adjust_start, index_adjust_end

      SELECT CASE (mpi_options)
      CASE (1)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO i = 1, ni, 1
               DO k = 1, mpi_nlpb_shadow_block_1d_all(j), 1
                  index_original_start = &
                     mpi_nlpb_shadow_block_starts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k)
                  index_original_end = index_original_start + &
                     mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
                  index_adjust_end = index_adjust_start + &
                     mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
                  mpi_adjust_buf(index_adjust_start:index_adjust_end) = &
                     mpi_source_buf(index_original_start:index_original_end, i)
                  index_adjust_start = index_adjust_end + 1
               END DO
            END DO
         END DO

      CASE (2)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO i = 1, ni, 1
               DO k = 1, mpi_nlpb_shadow_nlpbz_block_1d_all(j), 1
                  index_original_start = &
                     mpi_nlpb_shadow_nlpbz_block_starts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k)
                  index_original_end = index_original_start + &
                     mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k) - 1
                  index_adjust_end = index_adjust_start + &
                     mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k) - 1
                  mpi_adjust_buf(index_adjust_start:index_adjust_end) = &
                     mpi_source_buf(index_original_start:index_original_end, i)
                  index_adjust_start = index_adjust_end + 1
               END DO
            END DO
         END DO

      CASE (3)

      CASE (4)

      CASE DEFAULT
         STOP
         WRITE (*, *) "wrong options for netcdf reading"
      END SELECT

   END SUBROUTINE mpi_2d_io_read_adjust_integer

   ! adjust sending buffer of wp of 2 dimension for data exchange
   SUBROUTINE mpi_2d_io_read_adjust_wp(mpi_source_buf, mpi_adjust_buf, mpi_options,optional_dim1)
      REAL(wp), DIMENSION(:, :), INTENT(IN) :: mpi_source_buf
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: mpi_adjust_buf
      INTEGER, INTENT(IN) :: mpi_options
      INTEGER, INTENT(IN),OPTIONAL :: optional_dim1
      INTEGER ::  i, j, k, index_original_start, index_original_end, index_adjust_start, index_adjust_end

      SELECT CASE (mpi_options)

      CASE (3)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO i = 1, nk, 1
               DO k = 1, mpi_nlpb_shadow_block_1d_all(j), 1
                  index_original_start = mpi_nlpb_shadow_block_starts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k)
                  index_original_end = index_original_start + &
                                      mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
                  index_adjust_end = index_adjust_start + &
                                      mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
                  mpi_adjust_buf(index_adjust_start:index_adjust_end) = mpi_source_buf(index_original_start:index_original_end, i)
                  index_adjust_start = index_adjust_end + 1
               END DO
            END DO
         END DO

      CASE (5)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO i = 1, sbct, 1
               DO k = 1, mpi_nlpb_shadow_block_1d_all(j), 1
                  index_original_start = mpi_nlpb_shadow_block_starts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k)
                  index_original_end = index_original_start + &
                                      mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
                  index_adjust_end = index_adjust_start + &
                                      mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
                  mpi_adjust_buf(index_adjust_start:index_adjust_end) = mpi_source_buf(index_original_start:index_original_end, i)
                  index_adjust_start = index_adjust_end + 1
               END DO
            END DO
         END DO

      CASE (6)
         index_adjust_start = 1
         DO j = 1, mpi_bdy_outdegree, 1
            DO i = 1, nk, 1
               DO k = 1, mpi_bdy_counts_1d_all(j), 1
                  index_original_start = mpi_original_to_ordered_bdy_all(mpi_bdy_displs_1d_all(j)+ k)
                  mpi_adjust_buf(index_adjust_start) = mpi_source_buf(index_original_start, i)
                  index_adjust_start = index_adjust_start + 1
               END DO
            END DO
         END DO

      CASE (8)
         index_adjust_start = 1
         DO j = 1, mpi_bdy_outdegree, 1
            DO i = 1, optional_dim1, 1
               DO k = 1, mpi_bdy_counts_1d_all(j), 1
                  index_original_start = mpi_original_to_ordered_bdy_all(mpi_bdy_displs_1d_all(j)+ k)
                  mpi_adjust_buf(index_adjust_start) = mpi_source_buf(index_original_start, i)
                  index_adjust_start = index_adjust_start + 1
               END DO
            END DO
         END DO

      CASE (9)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO i = 1, optional_dim1, 1
               DO k = 1, mpi_nlpb_shadow_block_1d_all(j), 1
                  index_original_start = mpi_nlpb_shadow_block_starts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k)
                  index_original_end = index_original_start + &
                                      mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
                  index_adjust_end = index_adjust_start + &
                                      mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
                  mpi_adjust_buf(index_adjust_start:index_adjust_end) = mpi_source_buf(index_original_start:index_original_end, i)
                  index_adjust_start = index_adjust_end + 1
               END DO
            END DO
         END DO

      CASE DEFAULT
         STOP
         WRITE (*, *) "wrong options for netcdf reading"
      END SELECT

   END SUBROUTINE mpi_2d_io_read_adjust_wp

   ! adjust sending buffer of i2 of 2 dimension for data exchange
   SUBROUTINE mpi_2d_io_read_adjust_i2(mpi_source_buf, mpi_adjust_buf, mpi_options)
      INTEGER(i2), DIMENSION(:, :), INTENT(IN) :: mpi_source_buf
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_adjust_buf
      INTEGER, INTENT(IN) :: mpi_options
      INTEGER ::  i, j, k, index_original_start, index_original_end, index_adjust_start, index_adjust_end

      SELECT CASE (mpi_options)
      CASE (1)

      CASE (2)

      CASE (3)

         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO i = 1, nk, 1
               DO k = 1, mpi_nlpb_shadow_block_1d_all(j), 1
                  index_original_start = mpi_nlpb_shadow_block_starts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k)
                  index_original_end = index_original_start + &
                                       mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
                  index_adjust_end = index_adjust_start + &
                                      mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
                  mpi_adjust_buf(index_adjust_start:index_adjust_end) = &
                                                                  INT(mpi_source_buf(index_original_start:index_original_end,i), &
                                                                  KIND(mpi_adjust_buf))
                  index_adjust_start = index_adjust_end + 1
               END DO
            END DO
         END DO
      CASE (4)

         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO i = 1, nk, 1
               DO k = 1, mpi_nlpb_shadow_nlpbz_block_1d_all(j), 1
                  index_original_start = mpi_nlpb_shadow_nlpbz_block_starts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k)
                  index_original_end = index_original_start + &
                     mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k) - 1
                  index_adjust_end = index_adjust_start + &
                     mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k) - 1
                  mpi_adjust_buf(index_adjust_start:index_adjust_end) = &
                     INT(mpi_source_buf(index_original_start:index_original_end,i),KIND(mpi_adjust_buf))
                  index_adjust_start = index_adjust_end + 1
               END DO
            END DO
         END DO
      CASE DEFAULT
         STOP
         WRITE (*, *) "wrong options for netcdf reading"
      END SELECT

   END SUBROUTINE mpi_2d_io_read_adjust_i2

   ! adjust sending buffer of i2 of 1 dimension for data exchange
   SUBROUTINE mpi_1d_io_read_adjust_i2(mpi_source_buf, mpi_adjust_buf, mpi_options)
      INTEGER(i2), DIMENSION(:), INTENT(IN) :: mpi_source_buf
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_adjust_buf
      INTEGER, INTENT(IN) :: mpi_options
      INTEGER ::  j, k, index_original_start, index_original_end, index_adjust_start, index_adjust_end

      SELECT CASE (mpi_options)
      CASE (1)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO k = 1, mpi_nlpb_shadow_block_1d_all(j), 1
               index_original_start = mpi_nlpb_shadow_block_starts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k)
               index_original_end = index_original_start + &
                                    mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
               index_adjust_end = index_adjust_start + &
                                    mpi_nlpb_shadow_block_counts_1d_all(mpi_nlpb_shadow_block_1d_displs_all(j) + k) - 1
               mpi_adjust_buf(index_adjust_start:index_adjust_end) = &
                                    INT(mpi_source_buf(index_original_start:index_original_end),KIND(mpi_adjust_buf))
               index_adjust_start = index_adjust_end + 1
            END DO
         END DO

      CASE (2)
         index_adjust_start = 1
         DO j = 1, mpi_comp_procs, 1
            DO k = 1, mpi_nlpb_shadow_nlpbz_block_1d_all(j), 1
               index_original_start = mpi_nlpb_shadow_nlpbz_block_starts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k)
               index_original_end = index_original_start + &
                                    mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k) - 1
               index_adjust_end = index_adjust_start + &
                                  mpi_nlpb_shadow_nlpbz_block_counts_1d_all(mpi_nlpb_shadow_nlpbz_block_1d_displs_all(j) + k) - 1
               mpi_adjust_buf(index_adjust_start:index_adjust_end) = &
                                  INT(mpi_source_buf(index_original_start:index_original_end),KIND(mpi_adjust_buf))
               index_adjust_start = index_adjust_end + 1
            END DO
         END DO
      CASE (3)

      CASE (4)

      CASE DEFAULT
         STOP
         WRITE (*, *) "wrong options for netcdf reading"
      END SELECT

   END SUBROUTINE mpi_1d_io_read_adjust_i2

   ! mpi finalizing operations
   SUBROUTINE mpi_final_operations

      ! deallocate some arrays used in netcdf reading
      ! CALL mpi_netcdf_read_finalize
      call mpi_barrier(mpi_comm_world,mpi_err)  !!!caoyu add for mpich
      CALL MPI_FINALIZE(mpi_err)

   END SUBROUTINE mpi_final_operations


   ! searching sourcing process ID for all computing processes in graph neighbouring communication
   SUBROUTINE mpi_get_graph_source
     INTEGER,ALLOCATABLE,DIMENSION(:,:) :: mpi_local_grids_all ! only for root computing process
     INTEGER,ALLOCATABLE,DIMENSION(:) :: mpi_local_grids_all_sendbuf ! only for root computing process
     INTEGER,ALLOCATABLE,DIMENSION(:) :: mpi_recv_indexes_1d_counts_tmp
     INTEGER,ALLOCATABLE,DIMENSION(:) :: mpi_recv_indexes_1d_displs_tmp
     INTEGER,ALLOCATABLE,DIMENSION(:) :: mpi_graph_sources_tmp
     CHARACTER(lc) :: filename
     INTEGER :: i,j

     ALLOCATE (mpi_recv_indexes_1d_counts_tmp(mpi_comp_procs - 1))
     ALLOCATE (mpi_recv_indexes_1d_displs_tmp(mpi_comp_procs - 1))
     ALLOCATE (mpi_graph_sources_tmp(mpi_comp_procs - 1))
     ALLOCATE(mpi_local_grids(nlpbz))

     IF (mpi_rank == 0) THEN
       filename = "nmefc_macom_par_dec.nc"
       ! scatter local grid index to each process
       ALLOCATE(mpi_local_grids_all_sendbuf(mpi_nlpb_shadow_nlpbz_counts_sum))
       ALLOCATE(mpi_local_grids_all(mpi_max_nlpb_shadow_nlpbz,mpi_comp_procs))
       CALL netcdf_read(trim(filename), 'loc_idx3', mpi_local_grids_all, mpi_max_nlpb_shadow_nlpbz, mpi_comp_procs)
       CALL mpi_integer_2d_prepare_sendbuf_for_par_dec(mpi_local_grids_all,mpi_local_grids_all_sendbuf,mpi_comp_procs, &
                                                       mpi_nlpb_shadow_nlpbz_counts_all)
       CALL MPI_SCATTERV(mpi_local_grids_all_sendbuf, mpi_nlpb_shadow_nlpbz_counts_all, &
                         mpi_nlpb_shadow_nlpbz_displs_all,MPI_INTEGER, &
                         mpi_local_grids, nlpbz, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
       DEALLOCATE(mpi_local_grids_all)


     ELSE
       CALL MPI_SCATTERV(mpi_local_grids_all_sendbuf, mpi_nlpb_shadow_nlpbz_counts_all, &
                         mpi_nlpb_shadow_nlpbz_displs_all, MPI_INTEGER, &
                         mpi_local_grids, nlpbz, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

     END IF

     ! check value in nmefc_macom_par_dec.nc for loc_idx3
     DO i = 1,loc_nlpb
      IF (mpi_local_grids(i)< 1 .OR. mpi_local_grids(i)> loc_nlpb) THEN
        WRITE(*,*)"loc_idx3(",i,") = ",mpi_local_grids(i),"in nmefc_macom_par_dec.nc is wrong"
        STOP
      END IF
     END DO

     DO i = loc_nlpb + 1,nlpb
      IF (mpi_local_grids(i)< 1 .OR. mpi_local_grids(i)> mpi_total_nlpb) THEN
        WRITE(*,*)"loc_idx3(",i,") = ",mpi_local_grids(i),"in nmefc_macom_par_dec.nc is wrong"
        STOP
      END IF
     END DO

#ifdef MPI_DEBUG_NOW
    WRITE(filename,*)mpi_rank,"_mpi_process_index.txt"
    CALL mpi_output_file_1d_integer(filename,mpi_process_index,nlpbz)
#endif
     ! get number of processes to receive data, size of received data and starting position of received data in buffer
     ! mpi_recv_indexes_1d_displs_tmp is displacement from mpi_process_index(1)
     mpi_graph_indegree = 0
     i = loc_nlpb + 1
     mpi_recv_indexes_1d_displs_tmp(1) = loc_nlpb
     DO j = loc_nlpb + 2,nlpb
       IF(mpi_process_index(i) .ne. mpi_process_index(j)) THEN
         mpi_graph_indegree = mpi_graph_indegree + 1
         mpi_graph_sources_tmp(mpi_graph_indegree) = mpi_process_index(i)
         mpi_recv_indexes_1d_counts_tmp(mpi_graph_indegree) = i - mpi_recv_indexes_1d_displs_tmp(mpi_graph_indegree)
         mpi_recv_indexes_1d_displs_tmp(mpi_graph_indegree + 1) = i
       END IF
       i = i + 1
     END DO
     mpi_graph_indegree = mpi_graph_indegree + 1
     mpi_graph_sources_tmp(mpi_graph_indegree) = mpi_process_index(i)
     mpi_recv_indexes_1d_counts_tmp(mpi_graph_indegree)= i - mpi_recv_indexes_1d_displs_tmp(mpi_graph_indegree)

     ALLOCATE(mpi_graph_sources(mpi_graph_indegree))
     ALLOCATE(mpi_recv_indexes_1d_counts(mpi_graph_indegree))
     ALLOCATE(mpi_recv_indexes_1d_displs(mpi_graph_indegree))
     mpi_graph_sources(1:mpi_graph_indegree) = mpi_graph_sources_tmp(1:mpi_graph_indegree)
     mpi_recv_indexes_1d_counts(1:mpi_graph_indegree)= mpi_recv_indexes_1d_counts_tmp(1:mpi_graph_indegree)

     ! mpi_recv_indexes_1d is displacement from mpi_process_index(loc_nlpb + 1)
     mpi_recv_indexes_1d_displs(1) = 0
     mpi_recv_indexes_1d_counts_sum = mpi_recv_indexes_1d_counts(1)
     DO i = 2, mpi_graph_indegree
        mpi_recv_indexes_1d_counts_sum = mpi_recv_indexes_1d_counts_sum + mpi_recv_indexes_1d_counts(i)
        mpi_recv_indexes_1d_displs(i) = mpi_recv_indexes_1d_displs(i - 1) + mpi_recv_indexes_1d_counts(i - 1)
     END DO

     DEALLOCATE(mpi_graph_sources_tmp)
     DEALLOCATE(mpi_recv_indexes_1d_counts_tmp)
     DEALLOCATE(mpi_recv_indexes_1d_displs_tmp)
   END SUBROUTINE mpi_get_graph_source

   ! change globle index to local index for index arrays
   SUBROUTINE mpi_adjust_indexes
      CALL mpi_adjust_sub_indexes(tw, nlpb)
      CALL mpi_adjust_sub_indexes(ts, nlpb)
      CALL mpi_adjust_sub_indexes(te, nlpb)
      CALL mpi_adjust_sub_indexes(tn, nlpb)
      CALL mpi_adjust_sub_indexes(zw, nlpbz)
      CALL mpi_adjust_sub_indexes(zs, nlpbz)
      CALL mpi_adjust_sub_indexes(ze, nlpbz)
      CALL mpi_adjust_sub_indexes(zn, nlpbz)
      CALL mpi_adjust_sub_indexes(auz1, nlpb)
      CALL mpi_adjust_sub_indexes(auz2, nlpb)
      CALL mpi_adjust_sub_indexes(auz3, nlpb)
      CALL mpi_adjust_sub_indexes(auz4, nlpb)
      CALL mpi_adjust_sub_indexes(auz5, nlpb)
      CALL mpi_adjust_sub_indexes(auz6, nlpb)
      CALL mpi_adjust_sub_indexes(avz1, nlpb)
      CALL mpi_adjust_sub_indexes(avz2, nlpb)
      CALL mpi_adjust_sub_indexes(avz3, nlpb)
      CALL mpi_adjust_sub_indexes(avz4, nlpb)
      CALL mpi_adjust_sub_indexes(avz5, nlpb)
      CALL mpi_adjust_sub_indexes(avz6, nlpb)
      CALL mpi_adjust_sub_indexes(uw(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(us(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(ue(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(un(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(vw(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(vs(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(ve(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(vn(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(z1(:, 1), nlpbz)
      CALL mpi_adjust_sub_indexes(z2(:, 1), nlpbz)
      CALL mpi_adjust_sub_indexes(z3(:, 1), nlpbz)
      CALL mpi_adjust_sub_indexes(z4(:, 1), nlpbz)
      CALL mpi_adjust_sub_indexes(au1(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(au2(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(au3(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(au4(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(av1(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(av2(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(av3(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(av4(:, 1), nlpb)
      CALL mpi_adjust_sub_indexes(tw2, nlpb)
      CALL mpi_adjust_sub_indexes(ts2, nlpb)
      CALL mpi_adjust_sub_indexes(te2, nlpb)
      CALL mpi_adjust_sub_indexes(tn2, nlpb)
   END SUBROUTINE mpi_adjust_indexes

   ! find grid index
# if defined (OPENACC)
   SUBROUTINE mpi_adjust_sub_indexes(source_indexes, indexes_size)
      INTEGER, INTENT(INOUT), DIMENSION(:) :: source_indexes
      INTEGER, INTENT(IN) :: indexes_size
      INTEGER :: i, j
      LOGICAL :: find_index

      !$acc data present(source_indexes, mpi_global_grids)
      !$acc kernels
      !$acc loop private(find_index)
      DO i = 1, indexes_size
         find_index = .FALSE.
         ! search from position where the last index is to be faster
         !$acc loop seq
         DO j = 1, nlpbz
            IF (source_indexes(i) == mpi_global_grids(j)) THEN
               source_indexes(i) = j
               find_index = .TRUE.
               EXIT
            END IF
         END DO

         IF (.NOT. find_index) THEN
            IF (i <= loc_nlpb) THEN
               !WRITE(*,*) "local index:",i,"can not be found in the local and shadow grids"
               STOP
            ELSE
               source_indexes(i) = i
            END IF
         END IF

      END DO
      !$acc end kernels
      !$acc end data

   END SUBROUTINE mpi_adjust_sub_indexes
# else
   SUBROUTINE mpi_adjust_sub_indexes(source_indexes, indexes_size)
      INTEGER, INTENT(INOUT), DIMENSION(:) :: source_indexes
      INTEGER, INTENT(IN) :: indexes_size
      INTEGER :: i, j
      INTEGER :: start
      LOGICAL :: find_index

      start = 1

      DO i = 1, indexes_size
        find_index = .FALSE.
        ! search from position where the last index is to be faster
        DO j = start, nlpbz
           IF (source_indexes(i) == mpi_global_grids(j)) THEN
              source_indexes(i) = j
              start = j
              find_index = .TRUE.
              EXIT
           END IF
        END DO

        ! if not find, we search from 1 to the postion where the last index is
        IF ( .NOT. find_index) THEN
          DO j = 1,start
            IF (source_indexes(i) == mpi_global_grids(j)) THEN
               source_indexes(i) = j
               start = j
               find_index = .TRUE.
               EXIT
            END IF
          END DO
         END IF

         IF (.NOT. find_index) THEN
           IF (i <= loc_nlpb) THEN
             WRITE(*,*) mpi_rank, "local index: ", i, source_indexes(i), " can not be found in the local and shadow grids"
             STOP
           ELSE
             source_indexes(i) = i
!             WRITE(*,*) "shadow index:",i," in zone:",mpi_zone_index(i)," can not be found, it's index changed to itself'"
           END IF
         END IF

      END DO

   END SUBROUTINE mpi_adjust_sub_indexes
# endif

SUBROUTINE mpi_adjust_sub_indexes_ori(source_indexes, indexes_size)
   INTEGER, INTENT(INOUT), DIMENSION(:) :: source_indexes
   INTEGER, INTENT(IN) :: indexes_size
   INTEGER :: i, j
   INTEGER :: start
   LOGICAL :: find_index

   start = 1

   DO i = 1, indexes_size
     find_index = .FALSE.
     ! search from position where the last index is to be faster
     DO j = start, nlpbz
        IF (source_indexes(i) == mpi_global_grids(j)) THEN
           source_indexes(i) = j
           start = j
           find_index = .TRUE.
           EXIT
        END IF
     END DO

     ! if not find, we search from 1 to the postion where the last index is
     IF ( .NOT. find_index) THEN
       DO j = 1,start
         IF (source_indexes(i) == mpi_global_grids(j)) THEN
            source_indexes(i) = j
            start = j
            find_index = .TRUE.
            EXIT
         END IF
       END DO
      END IF

      IF (.NOT. find_index) THEN
        IF (i <= loc_nlpb) THEN
          WRITE(*,*) mpi_rank, "local index: ", i, source_indexes(i), " can not be found in the local and shadow grids"
          STOP
        ELSE
          source_indexes(i) = i
!             WRITE(*,*) "shadow index:",i," in zone:",mpi_zone_index(i)," can not be found, it's index changed to itself'"
        END IF
      END IF

   END DO

END SUBROUTINE mpi_adjust_sub_indexes_ori

   ! initializing graph neighbouring communication
   SUBROUTINE mpi_init_graph_index
      INTEGER :: i,j,mpi_start,mpi_end
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_block_1d_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_block_1d_starts_tmp
      ! create graph topology
      INTEGER, ALLOCATABLE :: weights(:)
      ALLOCATE(weights(mpi_graph_indegree))
      weights = 0
      CALL MPI_DIST_GRAPH_CREATE_ADJACENT(mpi_comp_comm, mpi_graph_indegree, mpi_graph_sources, weights, &
                                          mpi_graph_outdegree, mpi_graph_dests, weights, &
                                          MPI_INFO_NULL, .false., mpi_graph_comm, mpi_err)

      ! preparing data counts and displs for sending index that will be used in data exchange
      ALLOCATE (mpi_send_indexes_block_1d_counts_tmp(mpi_send_indexes_1d_counts_sum))
      ALLOCATE (mpi_send_indexes_block_1d_starts_tmp(mpi_send_indexes_1d_counts_sum))
      ALLOCATE (mpi_send_indexes_block_1d_all(mpi_graph_outdegree))
      mpi_send_indexes_block_1d_sum = 0
      mpi_send_indexes_block_1d_all = 0
      DO i = 1,mpi_graph_outdegree
        mpi_start = 1
        mpi_send_indexes_block_1d_sum = mpi_send_indexes_block_1d_sum + 1
        mpi_send_indexes_block_1d_starts_tmp(mpi_send_indexes_block_1d_sum) = &
          mpi_send_indexes_1d(mpi_send_indexes_1d_displs(i) + 1)
        DO mpi_end = 2,mpi_send_indexes_1d_counts(i)
          IF (mpi_send_indexes_1d(mpi_send_indexes_1d_displs(i)+ mpi_start) + 1 /= &
              mpi_send_indexes_1d(mpi_send_indexes_1d_displs(i)+ mpi_end)) THEN
           ! size is equal to the end index of continuous block minus the start index
            mpi_send_indexes_block_1d_counts_tmp(mpi_send_indexes_block_1d_sum) = &
              mpi_send_indexes_1d(mpi_send_indexes_1d_displs(i) + mpi_start) - &
              mpi_send_indexes_block_1d_starts_tmp(mpi_send_indexes_block_1d_sum) + 1
           ! set start index of next block and add next block to number of total blocks
            mpi_send_indexes_block_1d_sum = mpi_send_indexes_block_1d_sum + 1
            mpi_send_indexes_block_1d_starts_tmp(mpi_send_indexes_block_1d_sum) = &
              mpi_send_indexes_1d(mpi_send_indexes_1d_displs(i)+ mpi_end)
            mpi_send_indexes_block_1d_all(i)= mpi_send_indexes_block_1d_all(i) + 1
          END IF
          mpi_start = mpi_start + 1
        END DO
        ! add size of the last block in each outdegree
        mpi_send_indexes_block_1d_counts_tmp(mpi_send_indexes_block_1d_sum) = &
          mpi_send_indexes_1d(mpi_send_indexes_1d_displs(i) + mpi_start) - &
          mpi_send_indexes_block_1d_starts_tmp(mpi_send_indexes_block_1d_sum) + 1

        mpi_send_indexes_block_1d_all(i)= mpi_send_indexes_block_1d_all(i) + 1
      END DO

      ALLOCATE(mpi_send_indexes_block_1d_counts(mpi_send_indexes_block_1d_sum))
      ALLOCATE(mpi_send_indexes_block_1d_starts(mpi_send_indexes_block_1d_sum))
      mpi_send_indexes_block_1d_counts = mpi_send_indexes_block_1d_counts_tmp(1:mpi_send_indexes_block_1d_sum)
      mpi_send_indexes_block_1d_starts = mpi_send_indexes_block_1d_starts_tmp(1:mpi_send_indexes_block_1d_sum)

      ALLOCATE(mpi_send_indexes_block_1d_all_displs(mpi_graph_outdegree))
      mpi_send_indexes_block_1d_all_displs(1) = 0
      DO i = 2,mpi_graph_outdegree
        mpi_send_indexes_block_1d_all_displs(i)= mpi_send_indexes_block_1d_all_displs(i - 1)+ mpi_send_indexes_block_1d_all(i - 1)
      END DO

#ifdef MPI_DEBUG
      CALL mpi_check_indexes(i,j)
#endif

   END SUBROUTINE mpi_init_graph_index

   ! searching destination process ID for all computing processes
   SUBROUTINE mpi_get_graph_destination
      INTEGER :: i, k, dest_process,dest_position,start,dest_postion_end,source_start,source_end
      INTEGER :: mpi_graph_indegree_counts_sum
      INTEGER :: mpi_graph_outdegree_counts_sum
      INTEGER :: mpi_graph_total_recv_indexes_all_current_p
      INTEGER :: mpi_recv_indexes_1d_counts_sum_all_sum
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_graph_indegree_displs_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_1d_counts_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_graph_outdegree_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_graph_outdegree_displs_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_graph_dests_all_current_p
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_all_current_p
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_1d_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_1d_counts_sum_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_1d_counts_sum_displs_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_sum_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_sum_displs_all

      IF (mpi_rank == 0) THEN

         ! collecting graph information in process 0, then scatter them to other processes
         ALLOCATE (mpi_graph_indegree_all(mpi_comp_procs))
         ALLOCATE (mpi_graph_indegree_displs_all(mpi_comp_procs))
         CALL MPI_GATHER(mpi_graph_indegree, 1, MPI_INTEGER, mpi_graph_indegree_all, 1, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         ! collecting neighbouring process number for each process
         mpi_graph_indegree_counts_sum = mpi_graph_indegree_all(1)
         mpi_graph_indegree_displs_all(1) = 0
         DO i = 2, mpi_comp_procs
            mpi_graph_indegree_counts_sum = mpi_graph_indegree_counts_sum + mpi_graph_indegree_all(i)
            mpi_graph_indegree_displs_all(i) = mpi_graph_indegree_displs_all(i - 1) + mpi_graph_indegree_all(i - 1)
         END DO
         ALLOCATE (mpi_graph_sources_all(mpi_graph_indegree_counts_sum))
         CALL MPI_GATHERV(mpi_graph_sources, mpi_graph_indegree, MPI_INTEGER, mpi_graph_sources_all, mpi_graph_indegree_all, &
                          mpi_graph_indegree_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         ALLOCATE (mpi_recv_indexes_1d_counts_all(mpi_graph_indegree_counts_sum))
         CALL MPI_GATHERV(mpi_recv_indexes_1d_counts, mpi_graph_indegree, MPI_INTEGER, mpi_recv_indexes_1d_counts_all, &
                          mpi_graph_indegree_all, mpi_graph_indegree_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         ! get number of destitaton processes
         ALLOCATE (mpi_graph_outdegree_all(mpi_comp_procs))
         mpi_graph_outdegree_all = 0
         ! total indegrees are equal to total outdegrees
         mpi_graph_outdegree_counts_sum = mpi_graph_indegree_counts_sum
         DO i = 1, mpi_graph_outdegree_counts_sum
           mpi_graph_outdegree_all(mpi_graph_sources_all(i) + 1) = mpi_graph_outdegree_all(mpi_graph_sources_all(i)+ 1)+ 1
         END DO
         CALL MPI_SCATTER(mpi_graph_outdegree_all, 1, MPI_INTEGER, mpi_graph_outdegree, 1, MPI_INTEGER, 0, mpi_comp_comm,mpi_err)

         ! according to process order, put all destitaton process index and sending index counts together in the following codes
         ALLOCATE (mpi_graph_outdegree_displs_all(mpi_comp_procs))
         mpi_graph_outdegree_displs_all(1) = 0
         DO i = 2, mpi_comp_procs
            mpi_graph_outdegree_displs_all(i) = mpi_graph_outdegree_displs_all(i - 1) + mpi_graph_outdegree_all(i - 1)
         END DO
         ALLOCATE (mpi_graph_dests_all(mpi_graph_outdegree_counts_sum))
         ALLOCATE (mpi_graph_dests_all_current_p(mpi_comp_procs))

         ! according to received process ID and received number of grids, get sending process ID and sending number of grids
         ALLOCATE (mpi_send_indexes_1d_counts_all(mpi_graph_outdegree_counts_sum))
         ALLOCATE(mpi_recv_indexes_1d_counts_sum_all(mpi_comp_procs))
         mpi_recv_indexes_1d_counts_sum_all = 0
         mpi_graph_dests_all_current_p = 1
         start = 1
         DO i = 1,mpi_comp_procs
           DO k = 1,mpi_graph_indegree_all(i)
             dest_process = mpi_graph_sources_all(start) + 1
             dest_position = mpi_graph_outdegree_displs_all(dest_process) + mpi_graph_dests_all_current_p(dest_process)
             mpi_graph_dests_all(dest_position) = i - 1
             mpi_send_indexes_1d_counts_all(dest_position) = mpi_recv_indexes_1d_counts_all(start)
             mpi_recv_indexes_1d_counts_sum_all(i) = mpi_recv_indexes_1d_counts_sum_all(i) + mpi_recv_indexes_1d_counts_all(start)
             start = start + 1
             mpi_graph_dests_all_current_p(dest_process) = mpi_graph_dests_all_current_p(dest_process) + 1
           END DO
         END DO

         ALLOCATE (mpi_graph_dests(mpi_graph_outdegree))
         CALL MPI_SCATTERV(mpi_graph_dests_all, mpi_graph_outdegree_all, mpi_graph_outdegree_displs_all, MPI_INTEGER, &
                            mpi_graph_dests, mpi_graph_outdegree, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         ALLOCATE (mpi_send_indexes_1d_counts(mpi_graph_outdegree))
         CALL MPI_SCATTERV(mpi_send_indexes_1d_counts_all, mpi_graph_outdegree_all, mpi_graph_outdegree_displs_all, &
                            MPI_INTEGER, mpi_send_indexes_1d_counts, mpi_graph_outdegree, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

        ! gather receiving index ID
        ALLOCATE(mpi_recv_indexes_1d_counts_sum_displs_all(mpi_comp_procs))
        mpi_recv_indexes_1d_counts_sum_displs_all(1) = 0
        mpi_recv_indexes_1d_counts_sum_all_sum = mpi_recv_indexes_1d_counts_sum_all(1)
        DO i = 2,mpi_comp_procs
          mpi_recv_indexes_1d_counts_sum_displs_all(i)= mpi_recv_indexes_1d_counts_sum_displs_all(i - 1)+ &
            mpi_recv_indexes_1d_counts_sum_all(i - 1)
          mpi_recv_indexes_1d_counts_sum_all_sum = mpi_recv_indexes_1d_counts_sum_all_sum + mpi_recv_indexes_1d_counts_sum_all(i)
        END DO
        ALLOCATE(mpi_recv_indexes_1d_all(mpi_recv_indexes_1d_counts_sum_all_sum))
        CALL MPI_GATHERV(mpi_local_grids(loc_nlpb + 1), mpi_recv_indexes_1d_counts_sum, MPI_INTEGER, mpi_recv_indexes_1d_all, &
                          mpi_recv_indexes_1d_counts_sum_all, mpi_recv_indexes_1d_counts_sum_displs_all, MPI_INTEGER, 0, &
                          mpi_comp_comm, mpi_err)


        ! according to received process ID and index ID to get sending index ID
        ALLOCATE(mpi_send_indexes_1d_counts_sum_all(mpi_comp_procs))
        ALLOCATE(mpi_send_indexes_1d_counts_sum_displs_all(mpi_comp_procs))
        start = 1
        mpi_send_indexes_1d_counts_sum_all = 0
        DO k = 1,mpi_graph_outdegree_all(1)
          mpi_send_indexes_1d_counts_sum_all(1) = mpi_send_indexes_1d_counts_sum_all(1)+ mpi_send_indexes_1d_counts_all(start)
          start = start + 1
        END DO
        mpi_send_indexes_1d_counts_sum_displs_all(1) = 0
        DO i = 2,mpi_comp_procs
          DO k = 1,mpi_graph_outdegree_all(i)
            mpi_send_indexes_1d_counts_sum_all(i) = mpi_send_indexes_1d_counts_sum_all(i)+ mpi_send_indexes_1d_counts_all(start)
            start = start + 1
          END DO
          mpi_send_indexes_1d_counts_sum_displs_all(i)= mpi_send_indexes_1d_counts_sum_displs_all(i - 1) + &
            mpi_send_indexes_1d_counts_sum_all(i - 1)
        END DO
        ALLOCATE(mpi_send_indexes_1d_all(mpi_recv_indexes_1d_counts_sum_all_sum))
        ALLOCATE(mpi_send_indexes_1d_all_current_p(mpi_comp_procs))
        mpi_send_indexes_1d_all_current_p = 1
        start = 1
        source_start = 1
        DO i = 1,mpi_comp_procs
          DO k = 1,mpi_graph_indegree_all(i)
            dest_process = mpi_graph_sources_all(start) + 1
            dest_position = mpi_send_indexes_1d_counts_sum_displs_all(dest_process) + &
              mpi_send_indexes_1d_all_current_p(dest_process)
            dest_postion_end = dest_position + mpi_recv_indexes_1d_counts_all(start) - 1
            source_end = source_start + mpi_recv_indexes_1d_counts_all(start) - 1
            mpi_send_indexes_1d_all(dest_position:dest_postion_end)= mpi_recv_indexes_1d_all(source_start:source_end)
            mpi_send_indexes_1d_all_current_p(dest_process) = mpi_send_indexes_1d_all_current_p(dest_process) + &
              mpi_recv_indexes_1d_counts_all(start)
            start = start + 1
            source_start = source_end + 1
          END DO
        END DO
        ! initialize number and displacement of send indexes
        ALLOCATE(mpi_send_indexes_1d_displs(mpi_graph_outdegree))
        mpi_send_indexes_1d_displs(1) = 0
        mpi_send_indexes_1d_counts_sum = mpi_send_indexes_1d_counts(1)
        DO i = 2, mpi_graph_outdegree
           mpi_send_indexes_1d_displs(i) = mpi_send_indexes_1d_displs(i - 1) + mpi_send_indexes_1d_counts(i - 1)
           mpi_send_indexes_1d_counts_sum = mpi_send_indexes_1d_counts_sum + mpi_send_indexes_1d_counts(i)
        END DO
        ! scatter send index ID
        ALLOCATE(mpi_send_indexes_1d(mpi_send_indexes_1d_counts_sum))
        CALL MPI_SCATTERV(mpi_send_indexes_1d_all, mpi_send_indexes_1d_counts_sum_all, mpi_send_indexes_1d_counts_sum_displs_all, &
                           MPI_INTEGER, mpi_send_indexes_1d, mpi_send_indexes_1d_counts_sum, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         DEALLOCATE(mpi_graph_outdegree_all)
         DEALLOCATE(mpi_graph_outdegree_displs_all)
         DEALLOCATE(mpi_send_indexes_1d_counts_all)
         DEALLOCATE(mpi_graph_dests_all_current_p)
         DEALLOCATE(mpi_graph_indegree_displs_all)
         DEALLOCATE(mpi_recv_indexes_1d_counts_all)
         DEALLOCATE(mpi_recv_indexes_1d_counts_sum_all)
         DEALLOCATE(mpi_recv_indexes_1d_counts_sum_displs_all)
         DEALLOCATE(mpi_recv_indexes_1d_all)
         DEALLOCATE(mpi_send_indexes_1d_all_current_p)
         DEALLOCATE(mpi_send_indexes_1d_all)
         DEALLOCATE(mpi_send_indexes_1d_counts_sum_all)

      ELSE

         CALL MPI_GATHER(mpi_graph_indegree, 1, MPI_INTEGER, mpi_graph_indegree_all, 1, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         CALL MPI_GATHERV(mpi_graph_sources, mpi_graph_indegree, MPI_INTEGER, mpi_graph_sources_all, mpi_graph_indegree_all, &
                          mpi_graph_indegree_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         CALL MPI_GATHERV(mpi_recv_indexes_1d_counts, mpi_graph_indegree, MPI_INTEGER, mpi_recv_indexes_1d_counts_all, &
                          mpi_graph_indegree_all, mpi_graph_indegree_displs_all, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         CALL MPI_SCATTER(mpi_graph_outdegree_all, mpi_comp_procs, MPI_INTEGER, mpi_graph_outdegree, 1, &
                          MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         ALLOCATE (mpi_graph_dests(mpi_graph_outdegree))
         CALL MPI_SCATTERV(mpi_graph_dests_all, mpi_graph_outdegree_all, mpi_graph_outdegree_displs_all, MPI_INTEGER, &
                            mpi_graph_dests, mpi_graph_outdegree, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
         ALLOCATE (mpi_send_indexes_1d_counts(mpi_graph_outdegree))
         CALL MPI_SCATTERV(mpi_send_indexes_1d_counts_all, mpi_graph_outdegree_all, &
                           mpi_graph_outdegree_displs_all, MPI_INTEGER, mpi_send_indexes_1d_counts, &
                           mpi_graph_outdegree, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)

         CALL MPI_GATHERV(mpi_local_grids(loc_nlpb + 1), mpi_recv_indexes_1d_counts_sum, MPI_INTEGER, mpi_recv_indexes_1d_all, &
                          mpi_recv_indexes_1d_counts_sum_all, mpi_recv_indexes_1d_counts_sum_displs_all, MPI_INTEGER, 0, &
                          mpi_comp_comm, mpi_err)

         ! initialize number and displacement of send indexes
         ALLOCATE(mpi_send_indexes_1d_displs(mpi_graph_outdegree))
         mpi_send_indexes_1d_displs(1) = 0
         mpi_send_indexes_1d_counts_sum = mpi_send_indexes_1d_counts(1)
         DO i = 2, mpi_graph_outdegree
            mpi_send_indexes_1d_displs(i) = mpi_send_indexes_1d_displs(i - 1) + mpi_send_indexes_1d_counts(i - 1)
            mpi_send_indexes_1d_counts_sum = mpi_send_indexes_1d_counts_sum + mpi_send_indexes_1d_counts(i)
         END DO
        ! scatter send index ID
        ALLOCATE(mpi_send_indexes_1d(mpi_send_indexes_1d_counts_sum))
        CALL MPI_SCATTERV(mpi_send_indexes_1d_all, mpi_send_indexes_1d_counts_sum_all, mpi_send_indexes_1d_counts_sum_displs_all, &
                           MPI_INTEGER, mpi_send_indexes_1d, mpi_send_indexes_1d_counts_sum, MPI_INTEGER, 0, mpi_comp_comm, mpi_err)
      END IF


   END SUBROUTINE mpi_get_graph_destination

   ! write data of 1 dimension to file
   SUBROUTINE mpi_output_file_1d_dp(file_name, var, var_size)
      character(*), intent(in) :: file_name
      real(dp), intent(in), dimension(:) ::var
      integer, intent(in) :: var_size
      integer :: i
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, var_size
         WRITE (99, *) var(i)
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_output_file_1d_dp

   ! write data of 1 dimension to file
   SUBROUTINE mpi_output_file_1d_integer(file_name, var, var_size)
      character(*), intent(in) :: file_name
      integer, intent(in), dimension(:) ::var
      integer, intent(in) :: var_size
      integer :: i
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, var_size
         WRITE (99, *) var(i)
      END DO
      CLOSE (99)
   END SUBROUTINE

   ! find index position of first value which is equal to the searching value
   INTEGER FUNCTION mpi_findloc(a, val) RESULT(i)
      INTEGER, DIMENSION(:) :: a
      INTEGER :: val
      INTEGER :: test
      test = size(a)
      do i = 1, size(a)
         if (a(i) == val) exit
      enddo
      if (i > size(a)) i = 0
      return
   END FUNCTION mpi_findloc

END MODULE mod_mpi_interfaces
