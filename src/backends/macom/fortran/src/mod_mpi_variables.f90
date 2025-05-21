MODULE mod_mpi_variables
   use mpi
   use mod_misc_basic
   IMPLICIT NONE

# if defined (OPENACC)
   !For Multiple GPUs use
   INTEGER, PUBLIC :: ngpus    !ngpu:number of gpus per node. see subroutine mpi_init
# endif

   INTEGER, PUBLIC :: mpi_err ! mpi error parameter for all processes
   INTEGER, PUBLIC :: mpi_tag1 = 1 ! mpi tag parameter for all processes
   INTEGER, PUBLIC :: mpi_rank ! Process ID in all processes
   INTEGER, PUBLIC :: mpi_procs ! Total number of all processes
   INTEGER, PUBLIC :: mpi_group_all ! mpi group for all processes
   INTEGER, PUBLIC :: mpi_stat(MPI_STATUS_SIZE) ! non - blocking communication status
   INTEGER, PUBLIC :: mpi_real_wp ! data type for double or real
   INTEGER :: MPI_I2 ! equal to i2 type defined in the program
   INTEGER, PUBLIC, PARAMETER  :: DI = MPI_OFFSET_KIND       ! double integer

   INTEGER, PUBLIC :: mpi_comp_rank ! Process ID for computation
   INTEGER, ALLOCATABLE, DIMENSION(:), PUBLIC :: mpi_comp_all_ranks ! All process ID for computation
   INTEGER, PUBLIC :: mpi_comp_procs ! Total number of Processes for computation
   INTEGER, PUBLIC :: mpi_comp_comm ! mpi communicator for all computation
   INTEGER, PUBLIC :: mpi_graph_comm ! mpi communicator for all computation in graph for neighbouring communication
   INTEGER, PUBLIC :: mpi_group_comp ! mpi group for computation
   INTEGER, PUBLIC :: mpi_root_comp ! root process ID for computation
   INTEGER, PUBLIC :: mpi_comp_tag = 1 ! communication tag for computation

   INTEGER, PUBLIC :: mpi_io_rank ! Process ID in I / O
   INTEGER, ALLOCATABLE, DIMENSION(:), PUBLIC :: mpi_io_all_ranks ! All process ID in I / O
   INTEGER, PUBLIC :: mpi_io_procs ! Total number of Processes in I / O
   INTEGER, PUBLIC :: mpi_io_comm ! mpi communicator for all I / O
   INTEGER, PUBLIC :: mpi_group_io ! mpi group for I / O
   INTEGER, PUBLIC :: mpi_root_io ! root process ID for I / O
   INTEGER, PUBLIC :: mpi_io_tag = 2 ! communication tag for I / O

   INTEGER, PUBLIC :: mpi_group_bdy ! mpi group for boundary communication
   INTEGER, PUBLIC :: mpi_bdy_comm ! mpi communicator for boundary communication
   LOGICAL, PUBLIC :: mpi_check_bdy = .FALSE. ! judge whether the process contains boundary or the root process of 0

   INTEGER, PUBLIC :: mpi_comp_io_rank ! Process ID in computation and I / O
   INTEGER, ALLOCATABLE, DIMENSION(:), PUBLIC :: mpi_comp_io_mixed_ranks ! All processes ID in computation and one I / O root process
   INTEGER, PUBLIC :: mpi_comp_io_procs ! Total number of Processes in computation and I / O
   INTEGER, PUBLIC :: mpi_comp_io_comm ! mpi communicator for all computation and one I / O
   INTEGER, PUBLIC :: mpi_group_comp_io ! mpi group for computation and I / O
   INTEGER, PUBLIC :: mpi_root_comp_io ! root process ID from I / O communicator to exchange data in  computation and I / O
   INTEGER, PUBLIC :: mpi_sub_root_comp_io ! root process ID from computation communicator to exchange data in computation and I / O
   INTEGER, PUBLIC :: mpi_comp_io_tag = 3 ! communication tag for computation and I / O
   INTEGER, PUBLIC :: mpi_comp_io_status(MPI_STATUS_SIZE) ! mpi status for I / O
   INTEGER, PUBLIC :: mpi_req_comp_io ! request for non - block comunication in computation and I / O

   CHARACTER, ALLOCATABLE, DIMENSION(:) :: mpi_buf ! buffer for data sharing
   INTEGER :: mpi_position = 0 ! real size of data in mpi_pack
   INTEGER :: mpi_buf_size ! size of mpi_buf,it can be updated
   INTEGER :: mpi_mem_status ! status of memory allocated
   INTEGER, PUBLIC, PARAMETER :: MPI_LOG_UNIT = 90 ! unit for writing log files
   INTEGER, PUBLIC, PARAMETER :: MPI_INPUT_UNIT = 91 ! unit for reading log files
   INTEGER, PUBLIC :: mpi_io_stat ! status parameter in MPI functions
!   INTEGER, PUBLIC :: mpi_nlpbz  ! nlpbz + mpi_recv_indexes_1d_counts_sum
   CHARACTER(LEN = lc), PUBLIC :: mpi_io_msg ! message parameter in MPI functions

   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_global_grids ! save global grid index, including shadow and z grids
   !$ACC DECLARE CREATE(mpi_global_grids)
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_idx_bdy ! save global grid index for boundary
   !$acc declare create(mpi_idx_bdy)
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_ordered_idx_bdy_all ! save all ordered global grid index for boundary
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_original_idx_bdy_all ! save all global grid index for boundary
   ! the sequential index of ordered array in position of the original array
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_original_to_ordered_bdy_all
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_local_grids ! save local grid index, including shadow and z grids
   ! only for root computing process, maximium size of array for size of nlpb, nlpb + shadow, nlpb + shadow + z grids in a process
   INTEGER, PUBLIC :: mpi_max_nlpb_shadow_nlpbz
   INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_process_index ! process rank ID for receiving data
   INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_zone_index ! zone ID for computing grid layer or shadow layer

   ! number of uncontiguous blocks for 1 dimension for nlpb
   INTEGER, PUBLIC :: mpi_nlpb_block_1d
   ! starting postion of uncontiguous blocks for 1 dimension for nlpb
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_block_starts_1d
   ! size of each uncontiguous blocks for 1 dimension for nlpb
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_block_counts_1d
   ! only for root computing process, displacement of mpi_nlpb_shadow_nlpbz(1) of each process
   INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_displs_all
   ! only for root computing process and root I / O process, mpi_nlpb_shadow_nlpbz(1) of each process
   ! in root computing process, dimension size is mpi_comp_procs
   ! in root I / O process, dimension size is mpi_comp_procs + 1
   INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_counts_all
   INTEGER :: mpi_nlpb_counts_sum ! sum of mpi_nlpb_counts_all
   ! only for root computing process, number of uncontiguous blocks for 1 dimension for nlpb
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_block_1d_all
   ! only for root computing process, sum of mpi_nlpb_block_1d_all
   INTEGER, PUBLIC :: mpi_nlpb_block_1d_sum
   ! only for root computing process, displacement of mpi_nlpb_block_1d
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_block_1d_displs_all
   ! only for root computing process,starting postion of uncontiguous blocks for 1 dimension for nlpb
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_block_starts_1d_all
   ! only for root computing process,size of each uncontiguous blocks for 1 dimension for nlpb
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_block_counts_1d_all

   ! only for root computing process,size of sending data to each process for 1 dimension for bdy
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_bdy_counts_1d_all
   ! only for root computing process,displacement of sending data to each process for 1 dimension for bdy
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_bdy_displs_1d_all

   ! number of uncontiguous blocks for 1 dimension for nlpb and shadow
   INTEGER, PUBLIC :: mpi_nlpb_shadow_block_1d
   ! starting postion of uncontiguous blocks for 1 dimension for nlpb and shadow
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_block_starts_1d
   ! size of each uncontiguous blocks for 1 dimension for nlpb and shadow
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_block_counts_1d
   ! only for root computing process, displacement of mpi_nlpb_shadow_nlpbz(2) of each process
   INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_displs_all
   ! only for root computing process, mpi_nlpb_shadow_nlpbz(2) of each process
   INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_counts_all
   INTEGER :: mpi_nlpb_shadow_counts_sum ! sum of mpi_nlpb_shadow_counts_all
   ! only for root computing process, number of uncontiguous blocks for 1 dimension for nlpb and shadow
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_block_1d_all
   ! only for root computing process, sum of mpi_nlpb_shadow_block_1d_all
   INTEGER, PUBLIC :: mpi_nlpb_shadow_block_1d_sum
   ! only for root computing process, displacement of mpi_nlpb_shadow_block_1d
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_block_1d_displs_all
   ! only for root computing process,starting postion of uncontiguous blocks for 1 dimension for nlpb and shadow
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_block_starts_1d_all
   ! only for root computing process,size of each uncontiguous blocks for 1 dimension for nlpb and shadow
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_block_counts_1d_all

   ! number of uncontiguous blocks for 1 dimension for nlpb, shadow and z grids
   INTEGER, PUBLIC :: mpi_nlpb_shadow_nlpbz_block_1d
   ! starting postion of uncontiguous blocks for 1 dimension for nlpb, shadow and z grids
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_block_starts_1d
   ! size of each uncontiguous blocks for 1 dimension for nlpb, shadow and z grids
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_block_counts_1d
   ! only for root computing process, displacement of mpi_nlpb_shadow_nlpbz(3) of each process
   INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_displs_all
   ! only for root computing process, mpi_nlpb_shadow_nlpbz(3) of each process
   INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_counts_all
   INTEGER :: mpi_nlpb_shadow_nlpbz_counts_sum ! sum of mpi_nlpb_shadow_nlpbz_counts_all
   ! only for root computing process,number of uncontiguous blocks for 1 dimension for nlpb, shadow and z grids
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_block_1d_all
   ! only for root computing process, sum of mpi_nlpb_shadow_nlpbz_block_1d_all
   INTEGER, PUBLIC :: mpi_nlpb_shadow_nlpbz_block_1d_sum
   ! only for root computing process, displacement of mpi_nlpb_shadow_nlpbz_block_1d
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_block_1d_displs_all
   ! only for root computing process,starting postion of uncontiguous blocks for 1 dimension for nlpb, shadow and z grids
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_block_starts_1d_all
   ! only for root computing process, size of each uncontiguous blocks for 1 dimension for nlpb, shadow and z grids
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_nlpb_shadow_nlpbz_block_counts_1d_all

   INTEGER, PUBLIC :: mpi_total_nlpb ! number of all grids for nlpb
   INTEGER, PUBLIC :: mpi_total_nlpb_shadow ! number of all grids for nlpb and shadow in each process
   INTEGER, PUBLIC :: mpi_total_nlpb_shadow_nlpbz ! number of all grids for nlpb, shadow,and z grids in each process
   INTEGER, PUBLIC :: mpi_total_nlpbz ! number of all grids for nlpbz
   INTEGER, PUBLIC :: mpi_total_nl ! number of all grids for nl
   INTEGER, PUBLIC :: mpi_total_nlz ! number of all grids for nlz
   INTEGER, PUBLIC :: mpi_total_nlbdy ! number of all nlbdy, only for process 0
   INTEGER, PUBLIC :: mpi_total_nlbdy_shadow ! number of all nlbdy and shadow, only for process 0
   INTEGER, PUBLIC :: mpi_nlbdy ! number of nlbdy on each process
   INTEGER, PUBLIC :: mpi_nlbdy_shadow ! number of nlbdy and shadow on each process


   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlpb_shadow_ni_counts_all ! all nlpb and shadow in I / O reading
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlpb_shadow_ni_displs_all ! nlpb displs in I / O reading
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlpb_shadow_nk_counts_all ! all nlpb in I / O reading with nk
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlpb_shadow_nk_displs_all ! nlpb displs in I / O reading with nk
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlpb_shadow_sbct_counts_all ! all nlpb in I / O reading with sbct
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlpb_shadow_sbct_displs_all ! nlpb displs in I / O reading with sbct
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_user_counts_all ! all nlpb in I / O reading with sbct
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_user_displs_all ! nlpb displs in I / O reading with sbct
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlbdy_nk_counts_all ! number of mpi_bdy_counts_1d_all * nk in I / O reading for boundary
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlbdy_nk_displs_all ! displacement of mpi_bdy_counts_1d_all * nk in I / O reading for boundary

   ! save displacement of mpi_grid_2d_blocks_nlpbz for each process in data buffer in I / O reading
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlpb_shadow_nlpbz_ni_displs_all
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlpb_shadow_nlpbz_ni_counts_all ! all nlpbz in two dimension in I / O reading
   ! save displacement of nlpbz in data buffer in I / O reading with nlpbz and nk
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlpb_shadow_nlpbz_nk_displs_all
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_2d_nlpb_shadow_nlpbz_nk_counts_all ! nlpbz in I / O reading with nk

   ! used for buffer to read netcdf input of wp,only for computing root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_buf_1d_nlpb_wp
   ! used for buffer to read netcdf input of integer,only for computing root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_buf_1d_nlpb_integer
   ! used for buffer to read netcdf input of i2,only for computing root process
   INTEGER(i2), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_buf_1d_nlpb_i2
   ! used for adjusted buffer to read netcdf input of wp,only for computing root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_1d_nlpb_shadow_wp
   ! used for adjusted buffer to read netcdf input of integer,only computing root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_1d_nlpb_shadow_integer
   ! used for adjusted buffer to read netcdf input of i2,only for computing root process
   INTEGER(i2), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_1d_nlpb_shadow_i2
   ! used for buffer to read netcdf input of wp for nlpbz,only for computing root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_buf_1d_nlpb_nlpbz_wp
   ! used for buffer to read netcdf input of integer for nlpbz,only for computing root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_buf_1d_nlpb_nlpbz_integer
   ! used for adjusted buffer to read netcdf input of wp for nlpbz,only for computing root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_1d_nlpb_shadow_nlpbz_wp
   ! used for buffer to read netcdf input of integerfor nlpbz ,only for computing root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_1d_nlpb_shadow_nlpbz_integer
   ! used for buffer to read netcdf input of wp for boundary,only for computing root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_buf_1d_nlbdy_wp
   ! used for buffer to read netcdf input of wp for boundary and its halo, only for computing root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_1d_nlbdy_wp
   ! used for buffer to read netcdf input of wp for boundary,only for computing root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_buf_1d_nlbdy_integer
   ! used for buffer to read netcdf input of wp for boundary and its halo, only for computing root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_1d_nlbdy_integer

   ! used for buffer to read netcdf input of integer for 2 dimenstions,only for io root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_buf_2d_nlpb_ni_integer
   ! used for buffer to read netcdf input of i2 for 2 dimensions,only for io root process
   INTEGER(i2), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_buf_2d_nlpb_nk_i2
   ! used for buffer to read netcdf input of i2 for 2 dimensions for nlbpz,only for io root process
   INTEGER(i2), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_buf_2d_nlpb_nlpbz_nk_i2
   ! used for buffer to read netcdf input of wp for 2 dimensions,only for io root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_buf_2d_nlpb_nk_wp
   ! used for buffer to read netcdf input of wp for 2 dimensions,only for io root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_buf_2d_nlpb_sbct_wp
   ! used for buffer to read netcdf input of integer for 2 dimensions for nlpbz,only for io root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_buf_2d_nlpb_nlpbz_ni_integer
   ! used for buffer to read netcdf input of wp for 2 dimensions,only for io root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_buf_2d_nlbdy_nk_wp
   ! used for buffer to read netcdf input of wp for 2 dimensions,the second dimension is user - defined,
   ! so we need to allocate and deallocate according to requirements of user, only for io root process,
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_buf_2d_wp
   ! used for buffer to read netcdf input of wp for 2 dimensions,the second dimension is user - defined,
   ! so we need to allocate and deallocate according to requirements of user, only for io root process,
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_buf_adjust_2d_wp

   ! used for adjusted buffer to read netcdf input of integer for 2 dimensions with mpi_total_nlpb and ni,only for io root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_2d_nlpb_shadow_ni_integer
   ! used for adjusted buffer to read netcdf input of integer for 2 dimensions with mpi_total_nlpb and nk,only for io root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_2d_nlpb_shadow_nk_integer
   ! used for adjusted buffer to read netcdf input of integer for 2 dimensions with mpi_total_nlpbz and nk,only for io root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_2d_nlpb_shadow_nlpbz_nk_integer
   ! used for adjusted buffer to read netcdf input of integer for 2 dimensions with mpi_total_nlpbz and ni,only for io root process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_2d_nlpb_shadow_nlpbz_ni_integer
   ! used for adjusted buffer to read netcdf input of wp for 2 dimensions,only for io root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_2d_nlpb_shadow_nk_wp
   ! used for adjusted buffer to read netcdf input of wp for 2 dimensions,only for io root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_2d_nlpb_shadow_sbct_wp
   ! used for adjusted buffer to read netcdf input of wp for 2 dimensions,only for io root process
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_adjust_buf_2d_nlbdy_nk_wp

   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_1d_send_counts
   ! counts for 1 dimension in netcdf for writing for all processes
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_1d_send_displs
   ! displs for 1 dimension in netcdf for writing for all processes
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_send_counts
   ! counts for 2 dimensions in netcdf for writing for all processes
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_recv_count
   ! counts for 2 dimensions in netcdf for writing for only self - owned process
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_displs ! counts for 2 dimensions in netcdf for writing
   INTEGER, PUBLIC :: mpi_io_nk ! io writing is seperated by nk in each process


   ! following variables are used for non - overlapping communication for netcdf writing to receive or send data,
   ! only for io root process, tFld for all tracer, others for themself
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_io_buf_2d_tFld
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_uFld
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_vFld
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_wFld
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_KappaRT
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_KappaRM
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_1d_ssh
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_1d_etaH
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_1d_utau
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_1d_vtau
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_1d_qns
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_1d_qsr
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_1d_EmPmR
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_restart
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_1d_restart
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_tFld_adv
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_sFld_adv
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_tFld_hdiff
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_sFld_hdiff
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_tFld_vdiff
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_sFld_vdiff
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_tFld_force
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_2d_sFld_force
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_rhoLev
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_io_sbc
   !YY seaice buffer
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_io_nITD

   ! following variables are used to receive or send data for all computing processes
   ! the shape of variables are same with relating variables in computing processes
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_tFld
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_sFld
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_uFld
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_vFld
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_wFld
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_KappaRT
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_KappaRM
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_restart
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_comp_buf_1d_ssh
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_comp_buf_1d_etaH
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_comp_buf_1d_utau
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_comp_buf_1d_vtau
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_comp_buf_1d_qns
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_comp_buf_1d_qsr
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_comp_buf_1d_EmPmR
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_comp_buf_1d_restart
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_tFld_adv
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_sFld_adv
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_tFld_hdiff
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_sFld_hdiff
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_tFld_vdiff
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_sFld_vdiff
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_tFld_force
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_comp_buf_2d_sFld_force

   ! following variables are only used in root process, adjust data in right order
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_adjust_1d
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_io_buf_adjust_2d

   ! following variables are used for netcdf writing in each process, the order is same with final netcdf output files
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_io_buf_seperate_2d
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_io_buf_seperate_1d

   INTEGER, PUBLIC :: mpi_io_buf_seperate_1d_start ! real starts for netcdf writing in each process
   INTEGER, PUBLIC :: mpi_io_buf_seperate_1d_count ! real count for netcdf writing in each process
   INTEGER, PUBLIC :: mpi_io_buf_seperate_2d_start_nk(2) ! real starts for netcdf writing for 2 dimensions in each process
   INTEGER, PUBLIC :: mpi_io_buf_seperate_2d_count_nk(2) ! real count for netcdf writing for 2 dimensions in each process

   ! following variables are request in I / O for non - blocking coummication to check status of communication
   INTEGER, PUBLIC :: mpi_req_tFld
   INTEGER, PUBLIC :: mpi_req_sFld
   INTEGER, PUBLIC :: mpi_req_uFld
   INTEGER, PUBLIC :: mpi_req_vFld
   INTEGER, PUBLIC :: mpi_req_wFld
   INTEGER, PUBLIC :: mpi_req_KappaRT
   INTEGER, PUBLIC :: mpi_req_KappaRM
   INTEGER, PUBLIC :: mpi_req_ssh
   INTEGER, PUBLIC :: mpi_req_etaH
   INTEGER, PUBLIC :: mpi_req_qns
   INTEGER, PUBLIC :: mpi_req_qsr
   INTEGER, PUBLIC :: mpi_req_utau
   INTEGER, PUBLIC :: mpi_req_vtau
   INTEGER, PUBLIC :: mpi_req_EmPmR
   INTEGER, PUBLIC :: mpi_req_restart
   INTEGER, PUBLIC :: mpi_req_tFld_adv
   INTEGER, PUBLIC :: mpi_req_sFld_adv
   INTEGER, PUBLIC :: mpi_req_tFld_vdiff
   INTEGER, PUBLIC :: mpi_req_sFld_vdiff
   INTEGER, PUBLIC :: mpi_req_tFld_hdiff
   INTEGER, PUBLIC :: mpi_req_sFld_hdiff
   INTEGER, PUBLIC :: mpi_req_tFld_force
   INTEGER, PUBLIC :: mpi_req_sFld_force

   ! following variables are request in non - blocking neibouring communication in mpi_comp_comm
   INTEGER, PUBLIC :: mpi_comp_req_dxS
   INTEGER, PUBLIC :: mpi_comp_req_dyW
   INTEGER, PUBLIC :: mpi_comp_req_etaH
   INTEGER, PUBLIC :: mpi_comp_req_etaN
   INTEGER, PUBLIC :: mpi_comp_req_KE
   INTEGER, PUBLIC :: mpi_comp_req_hFacZ
   INTEGER, PUBLIC :: mpi_comp_req_recip_hFacZ
   INTEGER, PUBLIC :: mpi_comp_req_Vort3Loc
   INTEGER, PUBLIC :: mpi_comp_req_hDivLoc
   INTEGER, PUBLIC :: mpi_comp_req_del2T
   INTEGER, PUBLIC :: mpi_comp_req_uDragBottom
   INTEGER, PUBLIC :: mpi_comp_req_uFld
   INTEGER, PUBLIC :: mpi_comp_req_vFld
   INTEGER, PUBLIC :: mpi_comp_req_wFld
   INTEGER, PUBLIC :: mpi_comp_req_hFacC
   INTEGER, PUBLIC :: mpi_comp_req_hFacW
   INTEGER, PUBLIC :: mpi_comp_req_hFacS
   INTEGER, PUBLIC :: mpi_comp_req_rStarFacC
   INTEGER, PUBLIC :: mpi_comp_req_af
   INTEGER, PUBLIC :: mpi_comp_req_tl
   INTEGER, PUBLIC :: mpi_comp_req_gTracer
   INTEGER, PUBLIC :: mpi_comp_req_tracer
   INTEGER, PUBLIC :: mpi_comp_req_dTdl
   INTEGER, PUBLIC :: mpi_comp_req_dfl
   INTEGER, PUBLIC :: mpi_comp_req_uv
   INTEGER, PUBLIC :: mpi_comp_req_Bo_surf
   INTEGER, PUBLIC :: mpi_comp_req_hDiv
   INTEGER, PUBLIC :: mpi_comp_req_omega3

   ! following variables are used for flags of request in non - blocking communication
   LOGICAL, PUBLIC :: mpi_flag_tFld = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_sFld = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_uFld = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_vFld = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_wFld = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_KappaRT = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_KappaRM = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_ssh = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_etaH = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_qns = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_qsr = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_utau = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_vtau = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_EmPmR = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_restart = .FALSE.

   LOGICAL, PUBLIC :: mpi_flag_tFld_adv = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_sFld_adv = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_tFld_hdiff = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_sFld_hdiff = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_tFld_vdiff = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_sFld_vdiff = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_tFld_force = .FALSE.
   LOGICAL, PUBLIC :: mpi_flag_sFld_force = .FALSE.

   ! following variables are used to save name of output files
   CHARACTER(lc) :: mpi_file_var_name_ssh(2)
   CHARACTER(lc) :: mpi_file_var_name_tFld(2)

   ! it is in I / O processes, equal to myIter in computitional processes,
   ! there is no direct define of myIter in serial program, so we define it here
   INTEGER, PUBLIC :: mpi_io_myIter = 1

   INTEGER, PUBLIC :: mpi_total_2d_nlpb_nk ! whole size for 2D netcdf data
   ! displacement for each process on 2D netcdf data,only for root of comp_io_comm
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_1d_nlpb_counts_displs_all
   ! displacement for each process on 2D netcdf data,only for root of comp_io_comm
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_2d_nlpb_nk_counts_displs_all
   ! element counts for each process on 2D netcdf data,only for root of comp_io_comm
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_2d_nlpb_nk_counts_all

   INTEGER, PUBLIC :: mpi_send_indexes_1d_counts_sum ! all number of mpi_send_indexes_1d_counts
   INTEGER, PUBLIC :: mpi_recv_indexes_1d_counts_sum ! all number of mpi_recv_indexes_1d_counts
   INTEGER, PUBLIC :: mpi_send_indexes_block_1d_sum ! used for preparing sending buffer
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_block_1d_counts ! used for preparing size of sending buffer
   ! used for all computing processes, the number of blocks for each processes to be sent
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_block_1d_all
   ! used for preparing starting point of sending buffer
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_block_1d_starts
   ! used for preparing displacement of sending buffer
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_block_1d_all_displs
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d ! local 1 dimension indexes of sending
   !$ACC DECLARE CREATE(mpi_send_indexes_1d)
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_recv_total_indexes
   ! all global 1 dimension indexes of receiving,only for process 0
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts ! sending counts for each process
   !$ACC DECLARE CREATE(mpi_send_indexes_1d_counts)
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_1d_counts ! receiving counts for each process
   !$ACC DECLARE CREATE(mpi_recv_indexes_1d_counts)
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_displs ! sending displs for each process
   !$ACC DECLARE CREATE(mpi_send_indexes_1d_displs)
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_1d_displs ! receiving displs for each process
   !$ACC DECLARE CREATE(mpi_recv_indexes_1d_displs)

   ! receiving counts for each process,only for process 0
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_recv_total_indexes_counts

   ! source process ID of graph neighbour communication
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_sources
   ! total number of source process ID of graph neighbour communication,only for process 0
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_sources_all
   ! total destination process ID of graph neighbour communication,only for process 0
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_dests_all
   ! destination process ID of graph neighbour communication
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_dests
   ! number of mpi_graph_sources
   INTEGER, PUBLIC :: mpi_graph_indegree
   ! number of mpi_graph_dests_all
   INTEGER, PUBLIC :: mpi_graph_outdegree
   ! all indegree in mpi_rank order,only for process 0
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_indegree_all
   ! all indegree in mpi_rank order,only for process 0
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_outdegree_all


   ! all process ID for receiving data of boundary
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_bdy_dests_all

   ! number of processes for receiving data of boundary
   INTEGER, PUBLIC :: mpi_bdy_outdegree

   ! following variables are buffer for data exchange
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_sendbuf_1d
   !$ACC DECLARE CREATE(mpi_sendbuf_1d)
   !YY: sendbuf_2d_itd_f is used in SEAICE module for 2d vars (nlpb,itd)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_sendbuf_2d_nk_f, mpi_sendbuf_2d_itd_f
   !$ACC DECLARE CREATE(mpi_sendbuf_2d_nk_f,mpi_sendbuf_2d_itd_f)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_sendbuf_2d_nkp1
   !$ACC DECLARE CREATE(mpi_sendbuf_2d_nkp1)
   INTEGER(i2), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_sendbuf_2d_nk_int
   !$ACC DECLARE CREATE(mpi_sendbuf_2d_nk_int)

END MODULE mod_mpi_variables
