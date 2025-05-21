module mod_csp_tide
   use mod_csp_basic
   use mod_csp_misc
   use mod_misc
   implicit none

   character(len=4), public, save :: tide_name(tide_total)
   real(wp), public, save :: tide_freq(tide_total)

contains

!==============================================================================
  subroutine csp_tide_init
!==============================================================================
     implicit none
     integer :: i, j, tideDimIdx
     character(lc) :: filename, tideName
     integer :: tide_start(6), chkidx

     tide_name = (/'  Z0', '  SA', ' SSA', & !  1
           ' MSM', '  MM', ' MSF', '  MF', & !  2
           'ALP1', ' 2Q1', 'SIG1', '  Q1', & !  3
           'RHO1', '  O1', 'TAU1', 'BET1', & !  4
           ' NO1', 'CHI1', ' PI1', '  P1', & !  5
           '  S1', '  K1', 'PSI1', 'PHI1', & !  6
           'THE1', '  J1', '2PO1', ' SO1', & !  7
           ' OO1', 'UPS1', 'ST36', '2NS2', & !  8
           'ST37', ' ST1', ' OQ2', 'EPS2', & !  9
           ' ST2', ' ST3', '  O2', ' 2N2', & ! 10
           ' MU2', 'SNK2', '  N2', ' NU2', & ! 11
           ' ST4', ' OP2', 'GAM2', '  H1', & ! 12
           '  M2', '  H2', 'MKS2', ' ST5', & ! 13
           ' ST6', 'LDA2', '  L2', '2SK2', & ! 14
           '  T2', '  S2', '  R2', '  K2', & ! 15
           'MSN2', 'ETA2', ' ST7', '2SM2', & ! 16
           'ST38', 'SKM2', '2SN2', ' NO3', & ! 17
           ' MO3', '  M3', ' NK3', ' SO3', & ! 18
           ' MK3', ' SP3', ' SK3', ' ST8', & ! 19
           '  N4', '3MS4', 'ST39', ' MN4', & ! 20
           ' ST9', 'ST40', '  M4', 'ST10', & ! 21
           ' SN4', ' KN4', ' MS4', ' MK4', & ! 22
           ' SL4', '  S4', ' SK4', 'MNO5', & ! 23
           '2MO5', '3MP5', 'MNK5', '2MP5', & ! 24
           '2MK5', 'MSK5', '3KM5', '2SK5', & ! 25
           'ST11', '2NM6', 'ST12', '2MN6', & ! 26
           'ST13', 'ST41', '  M6', 'MSN6', & ! 27
           'MKN6', 'ST42', '2MS6', '2MK6', & ! 28
           'NSK6', '2SM6', 'MSK6', '  S6', & ! 29
           'ST14', 'ST15', '  M7', 'ST16', & ! 30
           '3MK7', 'ST17', 'ST18', '3MN8', & ! 31
           'ST19', '  M8', 'ST20', 'ST21', & ! 32
           '3MS8', '3MK8', 'ST22', 'ST23', & ! 33
           'ST24', 'ST25', 'ST26', '4MK9', & ! 34
           'ST27', 'ST28', ' M10', 'ST29', & ! 35
           'ST30', 'ST31', 'ST32', 'ST33', & ! 36
           ' M12', 'ST34', 'ST35', '  M1'/)  ! 37

      tide_freq = (/   0.0000000000_wp, 0.0001140741_wp, 0.0002281591_wp, & !  1
      0.0013097808_wp, 0.0015121518_wp, 0.0028219327_wp, 0.0030500918_wp, & !  2
      0.0343965699_wp, 0.0357063507_wp, 0.0359087218_wp, 0.0372185026_wp, & !  3
      0.0374208736_wp, 0.0387306544_wp, 0.0389588136_wp, 0.0400404353_wp, & !  4
      0.0402685944_wp, 0.0404709654_wp, 0.0414385130_wp, 0.0415525871_wp, & !  5
      0.0416666721_wp, 0.0417807462_wp, 0.0418948203_wp, 0.0420089053_wp, & !  6
      0.0430905270_wp, 0.0432928981_wp, 0.0443745198_wp, 0.0446026789_wp, & !  7
      0.0448308380_wp, 0.0463429898_wp, 0.0733553835_wp, 0.0746651643_wp, & !  8
      0.0748675353_wp, 0.0748933234_wp, 0.0759749451_wp, 0.0761773161_wp, & !  9
      0.0764054753_wp, 0.0772331498_wp, 0.0774613089_wp, 0.0774870970_wp, & ! 10
      0.0776894680_wp, 0.0787710897_wp, 0.0789992488_wp, 0.0792016198_wp, & ! 11
      0.0794555670_wp, 0.0802832416_wp, 0.0803090296_wp, 0.0803973266_wp, & ! 12
      0.0805114007_wp, 0.0806254748_wp, 0.0807395598_wp, 0.0809677189_wp, & ! 13
      0.0815930224_wp, 0.0818211815_wp, 0.0820235525_wp, 0.0831051742_wp, & ! 14
      0.0832192592_wp, 0.0833333333_wp, 0.0834474074_wp, 0.0835614924_wp, & ! 15
      0.0848454852_wp, 0.0850736443_wp, 0.0853018034_wp, 0.0861552660_wp, & ! 16
      0.0863576370_wp, 0.0863834251_wp, 0.0876674179_wp, 0.1177299033_wp, & ! 17
      0.1192420551_wp, 0.1207671010_wp, 0.1207799950_wp, 0.1220639878_wp, & ! 18
      0.1222921469_wp, 0.1248859204_wp, 0.1251140796_wp, 0.1566887168_wp, & ! 19
      0.1579984976_wp, 0.1582008687_wp, 0.1592824904_wp, 0.1595106495_wp, & ! 20
      0.1597388086_wp, 0.1607946422_wp, 0.1610228013_wp, 0.1612509604_wp, & ! 21
      0.1623325821_wp, 0.1625607413_wp, 0.1638447340_wp, 0.1640728931_wp, & ! 22
      0.1653568858_wp, 0.1666666667_wp, 0.1668948258_wp, 0.1982413039_wp, & ! 23
      0.1997534558_wp, 0.1999816149_wp, 0.2012913957_wp, 0.2025753884_wp, & ! 24
      0.2028035475_wp, 0.2056254802_wp, 0.2058536393_wp, 0.2084474129_wp, & ! 25
      0.2372259056_wp, 0.2385098983_wp, 0.2387380574_wp, 0.2400220501_wp, & ! 26
      0.2402502093_wp, 0.2413060429_wp, 0.2415342020_wp, 0.2428439828_wp, & ! 27
      0.2430721419_wp, 0.2441279756_wp, 0.2443561347_wp, 0.2445842938_wp, & ! 28
      0.2458940746_wp, 0.2471780673_wp, 0.2474062264_wp, 0.2500000000_wp, & ! 29
      0.2787527046_wp, 0.2802906445_wp, 0.2817899023_wp, 0.2830867891_wp, & ! 30
      0.2833149482_wp, 0.2861368809_wp, 0.3190212990_wp, 0.3205334508_wp, & ! 31
      0.3207616099_wp, 0.3220456027_wp, 0.3233553835_wp, 0.3235835426_wp, & ! 32
      0.3248675353_wp, 0.3250956944_wp, 0.3264054753_wp, 0.3276894680_wp, & ! 33
      0.3279176271_wp, 0.3608020452_wp, 0.3623141970_wp, 0.3638263489_wp, & ! 34
      0.3666482815_wp, 0.4010448515_wp, 0.4025570033_wp, 0.4038667841_wp, & ! 35
      0.4053789360_wp, 0.4069168759_wp, 0.4082008687_wp, 0.4471596822_wp, & ! 36
      0.4830684040_wp, 0.4858903367_wp, 0.4874282766_wp, 0.0402557003_wp/)  ! 37

      if (mpi_comp_rank == 0) then
         CALL MPI_BCAST(tide_name_use, tide_total*4, mpi_character, mpi_root_comp, mpi_comp_comm, mpi_err)
      else
         CALL MPI_BCAST(tide_name_use, tide_total*4, mpi_character, mpi_root_comp, mpi_comp_comm, mpi_err)
      end if

      ! first find how many tide we will use
      do i = 1, tide_total
         inner1: do j = 1, tide_total
            if ( adjustl(trim(tide_name(j))) ==  adjustl(trim(tide_name_use(i)))) then
               ntide = ntide + 1
               exit inner1
            end if
         end do inner1
      end do

      ! second find index of tide we use in tide parameter table
      allocate(tide_idx(ntide))
      ntide = 0
      do i = 1, tide_total
         inner2: do j = 1, tide_total
            if ( adjustl(trim(tide_name(j))) ==  adjustl(trim(tide_name_use(i)))) then
               ntide = ntide + 1
               tide_idx(ntide) = j
               exit inner2
            end if
         end do inner2
      end do

      ! pick up tide consituent we used for calculate speed
      deallocate(tide_name_use)
      allocate(tide_name_use(ntide), tide_freq_use(ntide))
      do i = 1, ntide
         tide_name_use(i) = adjustl(trim(tide_name(tide_idx(i))))
         tide_freq_use(i) = tide_freq(tide_idx(i))
      end do

      !write(*,*) "myid = ", mpi_comp_rank, tide_idx, tide_name_use, tide_freq_use

      if (tide2d) then
         tideDim = nlpb
         tideDimIdx = 9
      else
         tideDim = nlbdy
         tideDimIdx = 5
      end if

      allocate(tide_amp_ssh(tideDim, ntide), tide_pha_ssh(tideDim, ntide))
      allocate(tide_u_ssh(tideDim, ntide), tide_v_ssh(tideDim, ntide), tide_f_ssh(tideDim, ntide))
      allocate(tide_amp_u(tideDim, ntide), tide_pha_u(tideDim, ntide))
      allocate(tide_u_curu(tideDim, ntide), tide_v_curu(tideDim, ntide), tide_f_curu(tideDim, ntide))
      allocate(tide_amp_v(tideDim, ntide), tide_pha_v(tideDim, ntide))
      allocate(tide_u_curv(tideDim, ntide), tide_v_curv(tideDim, ntide), tide_f_curv(tideDim, ntide))

      allocate(tide_ssh(tideDim), tide_u(tideDim), tide_v(tideDim))

      tide_ssh = 0.0_wp
      tide_u = 0.0_wp
      tide_v = 0.0_wp

      filename = "nmefc_macom_tide_harmonic.nc"

      IF (mpi_rank == 0) THEN
         call netcdf_read(filename, 'dateNodal', tide_start, d1 = 6)
         CALL MPI_BCAST(tide_start, 6, mpi_integer, mpi_root_comp, mpi_comp_comm, mpi_err)
      else
         CALL MPI_BCAST(tide_start, 6, mpi_integer, mpi_root_comp, mpi_comp_comm, mpi_err)
      end if

      call greg2jul(tide_start(6), tide_start(5), tide_start(4), tide_start(3), &
           tide_start(2), tide_start(1), dateNodalJulian)

      CALL mpi_netcdf_read_exchange(filename, "tide_amp_ssh", tide_amp_ssh, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_pha_ssh", tide_pha_ssh, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_u_ssh", tide_u_ssh, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_v_ssh", tide_v_ssh, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_f_ssh", tide_f_ssh, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_amp_u", tide_amp_u, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_pha_u", tide_pha_u, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_u_curu", tide_u_curu, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_v_curu", tide_v_curu, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_f_curu", tide_f_curu, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_amp_v", tide_amp_v, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_pha_v", tide_pha_v, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_u_curv", tide_u_curv, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_v_curv", tide_v_curv, tideDimIdx, optional_dim1 = ntide)
      CALL mpi_netcdf_read_exchange(filename, "tide_f_curv", tide_f_curv, tideDimIdx, optional_dim1 = ntide)

   end subroutine csp_tide_init

!==============================================================================
  subroutine csp_tide_update
!==============================================================================
     implicit none
     integer :: i, j, k
     complex(wp) :: ap, am, tide_values
     real(wp) :: temp, rpifreq
     real(dp) :: tim

     !$ACC data present(tide_ssh, tide_u, tide_v, tide_freq_use, &
     !$ACC              tide_amp_ssh, tide_pha_ssh, tide_u_ssh, tide_v_ssh, tide_f_ssh, &
     !$ACC              tide_amp_u, tide_pha_u, tide_u_curu, tide_v_curu, tide_f_curu, &
     !$ACC              tide_amp_v, tide_pha_v, tide_u_curv, tide_v_curv, tide_f_curv)
     !$ACC kernels
     !$ACC loop
     do i = 1, nlpb
        tide_ssh(i) = 0.0_wp
        tide_u(i) = 0.0_wp
        tide_v(i) = 0.0_wp
     end do
     !$ACC end kernels

     tim = (dble(nowjulian) - dble(dateNodalJulian))/3600.0_wp

     !$ACC kernels copyin(tim)
     !$ACC loop independent private(rpifreq, ap, am, tide_values, temp)
     do i = 1, tideDim
        !$ACC loop seq
        do k = 1, ntide
           rpifreq = 2.0_wp*rpi*tide_freq_use(k)*tim

           ! update tide ssh
           ap = 0.5_wp*tide_amp_ssh(i,k)*exp(unitConj*tide_pha_ssh(i,k)*radPi)
           am = conjg(ap)
           temp = 2.0_wp*rpi*(tide_u_ssh(i,k)+tide_v_ssh(i,k))
           ap = ap*tide_f_ssh(i,k)*exp(unitC   *temp)
           am = am*tide_f_ssh(i,k)*exp(unitConj*temp)
           tide_values = exp(unitC*rpifreq)*ap + exp(unitConj*rpifreq)*am
           tide_ssh(i) = tide_ssh(i) + real(tide_values)

           ! update tide u
           ap = 0.5_wp*tide_amp_u(i,k)*exp(unitConj*tide_pha_u(i,k)*radPi)
           am = conjg(ap)
           temp = 2.0_wp*rpi*(tide_u_curu(i,k)+tide_v_curu(i,k))
           ap = ap*tide_f_curu(i,k)*exp(unitC   *temp)
           am = am*tide_f_curu(i,k)*exp(unitConj*temp)
           tide_values = exp(unitC*rpifreq)*ap + exp(unitConj*rpifreq)*am
           tide_u(i) = tide_u(i) + real(tide_values)

           ! update tide v
           ap = 0.5_wp*tide_amp_v(i,k)*exp(unitConj*tide_pha_v(i,k)*radPi)
           am = conjg(ap)
           temp = 2.0_wp*rpi*(tide_u_curv(i,k)+tide_v_curv(i,k))
           ap = ap*tide_f_curv(i,k)*exp(unitC   *temp)
           am = am*tide_f_curv(i,k)*exp(unitConj*temp)
           tide_values = exp(unitC*rpifreq)*ap + exp(unitConj*rpifreq)*am
           tide_v(i) = tide_v(i) + real(tide_values)
        end do
     end do
     !$ACC end kernels

     !$ACC end data

  end subroutine csp_tide_update

end module mod_csp_tide
