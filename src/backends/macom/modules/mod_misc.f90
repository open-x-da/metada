module mod_misc
   use mod_misc_basic
   use mod_csp_basic
   USE mod_mpi_interfaces
   USE mod_mpi_variables
   implicit none

contains

!==============================================================================
   subroutine misc_run_info_open
!==============================================================================
      implicit none

      integer :: istat
      character(lc) :: msg, filename, domid_str

      run_out_unit = idom*10
      sol_out_unit = idom*10 + 1

      write (domid_str, "(i2.2)") idom

      IF (mpi_rank == 0) THEN
         filename = "nmefc_macom_run_information_d"//trim(domid_str)//".output"
         open (file=trim(filename), unit=run_out_unit, action='write', iostat=istat, iomsg=msg)
         if (istat /= 0) then
            write (*, *) 'nmefc_macom_run_information.output open failed, error message = ', msg
            stop
         end if
         write (unit=run_out_unit, fmt=*) '         NMEFC MaCOM run information'

         filename = "nmefc_macom_ssh_solver_information_d"//trim(domid_str)//".output"
         open (file=trim(filename), unit=sol_out_unit, action='write', iostat=istat, iomsg=msg)
         if (istat /= 0) then
            write (*, *) 'nmefc_macom_ssh_solver_information.output open failed, error message = ', msg
            stop
         end if
      END IF

   end subroutine misc_run_info_open

!==============================================================================
   subroutine misc_run_info_close
!==============================================================================
      implicit none
      integer :: istat
      character(lc) :: msg, time, date, zone, timestamp

      IF (mpi_rank == 0) THEN
         call date_and_time(date=date, time=time, zone=zone)
         timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                           time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//zone
         write (run_out_unit,"(a)") ' All Integral Step Success Finished !'
         write (run_out_unit,"(a)") trim(timestamp)
         close (unit=run_out_unit, status='keep', iostat=istat, iomsg=msg)
         if (istat /= 0) then
            write (*, *) 'nmefc_macom_run_information.output close failed, error message = ', msg
         end if

         close (unit=sol_out_unit, status='keep', iostat=istat, iomsg=msg)
         if (istat /= 0) then
            write (*, *) 'nmefc_macom_ssh_solver_information.output close failed, error message = ', msg
         end if
      END IF

   end subroutine misc_run_info_close

!==============================================================================
   subroutine misc_namelist_read
!==============================================================================
      implicit none

      integer :: namelist_unit
      integer :: istat
      character(lc) :: msg

      namelist /model_run_general/ ndomain
      namelist /model_run_general/ idom
      namelist /model_run_general/ parent
      namelist /model_run_general/ nRatio
      namelist /model_run_general/ nIter0
      namelist /model_run_general/ nIterMax
      namelist /model_run_general/ nOutFreq
      namelist /model_run_general/ fOutFreq
      namelist /model_run_general/ nResFreq
      namelist /model_run_general/ startyear
      namelist /model_run_general/ startmon
      namelist /model_run_general/ startday
      namelist /model_run_general/ starthour
      namelist /model_run_general/ startmin
      namelist /model_run_general/ startsec
      namelist /model_run_general/ dTtracer
      namelist /model_run_general/ dTmom
      namelist /model_run_general/ dTsurf
      namelist /model_run_general/ gammas
      namelist /model_run_general/ gammat
      namelist /model_run_general/ ffsfrs
      namelist /model_run_general/ ffsfrt
      namelist /model_run_general/ nffsfrs
      namelist /model_run_general/ nffsfrt
      namelist /model_run_general/ fsfrsavg
      namelist /model_run_general/ fsfrtavg
      namelist /model_run_general/ sbc_cli
      namelist /model_run_general/ sbc_cyc
      namelist /model_run_general/ ffuv
      namelist /model_run_general/ fftq
      namelist /model_run_general/ ffrad
      namelist /model_run_general/ ffprc
      namelist /model_run_general/ ffslp
      namelist /model_run_general/ ffrunoff
      namelist /model_run_general/ nffuv
      namelist /model_run_general/ nfftq
      namelist /model_run_general/ nffrad
      namelist /model_run_general/ nffprc
      namelist /model_run_general/ nffslp
      namelist /model_run_general/ nffrunoff
      namelist /model_run_general/ fuvavg
      namelist /model_run_general/ ftqavg
      namelist /model_run_general/ fradavg
      namelist /model_run_general/ fprcavg
      namelist /model_run_general/ fslpavg
      namelist /model_run_general/ frunoffavg
      namelist /model_run_general/ ln_qdew
      namelist /model_run_general/ ln_rain
      namelist /model_run_general/ ln_EmPmRevise
      namelist /model_run_general/ zqt
      namelist /model_run_general/ zuv
      namelist /model_run_general/ rhoref
      namelist /model_run_general/ ln_bdy
      namelist /model_run_general/ fbvec
      namelist /model_run_general/ fbt
      namelist /model_run_general/ fbs
      namelist /model_run_general/ fbpbt
      namelist /model_run_general/ nfbvec
      namelist /model_run_general/ nfbt
      namelist /model_run_general/ nfbs
      namelist /model_run_general/ nfbpbt
      namelist /model_run_general/ bvecavg
      namelist /model_run_general/ btavg
      namelist /model_run_general/ bsavg
      namelist /model_run_general/ bpbtavg
      namelist /model_run_general/ fsbc_dir
      namelist /model_run_general/ fbdy_dir
      namelist /model_run_general/ fout_dir
      namelist /model_run_general/ restart_in
      namelist /model_run_general/ assim_in
      namelist /model_run_general/ restart_out
      namelist /model_run_general/ date_res
!YY
      namelist /model_run_general/ mitice_on

      namelist /model_dynamic_csp/ ViscAhDCon
      namelist /model_dynamic_csp/ ViscAhZCon
      namelist /model_dynamic_csp/ ViscA4DCon
      namelist /model_dynamic_csp/ ViscA4ZCon
      namelist /model_dynamic_csp/ harmonic
      namelist /model_dynamic_csp/ biharmonic
      namelist /model_dynamic_csp/ bottomDrag
      namelist /model_dynamic_csp/ BottomDragMax
      namelist /model_dynamic_csp/ KappaRMCon
      namelist /model_dynamic_csp/ abEps
      namelist /model_dynamic_csp/ rStarFacLow
      namelist /model_dynamic_csp/ rStarFacUp
      namelist /model_dynamic_csp/ ln_pbt_base
      namelist /model_dynamic_csp/ ln_bous
      namelist /model_dynamic_csp/ ln_tide
      namelist /model_dynamic_csp/ tide_name_use

      namelist /model_thermodynamic_csp/ ntracer
      namelist /model_thermodynamic_csp/ rn_abs
      namelist /model_thermodynamic_csp/ rn_si0
      namelist /model_thermodynamic_csp/ KappaRTCon
      namelist /model_thermodynamic_csp/ harmonicT
      namelist /model_thermodynamic_csp/ biharmonicT
      namelist /model_thermodynamic_csp/ diffKhCon
      namelist /model_thermodynamic_csp/ diffK4Con

      allocate(tide_name_use(tide_total))
      tide_name_use(:) = 'AAAA'

      IF (mpi_rank == 0) THEN
         namelist_unit = 15

         open (file='namelist.nmefc_macom', unit=namelist_unit, action='read', &
               iostat=istat, iomsg=msg)
         if (istat /= 0) then
            write (*, *) 'namelist.nmefc_macom open failed, error message = ', msg
            write (*, *) 'we use default parameters'
         else
            read (unit=namelist_unit, nml=model_run_general, iostat=istat, iomsg=msg)
            if (istat /= 0) then
               write (*, *) 'model_run_general read failed, error message = ', msg
               write (*, *) 'we use default parameters'
            else
               write (*, *) 'read parameters for model_run_general successful'
            end if

            read (unit=namelist_unit, nml=model_dynamic_csp, iostat=istat, iomsg=msg)
            if (istat /= 0) then
               write (*, *) 'model_dynamic_csp read failed, error message = ', msg
               write (*, *) 'we use default parameters'
            else
               write (*, *) 'read parameters for model_dynamic_csp successful'
            end if

            read (unit=namelist_unit, nml=model_thermodynamic_csp, iostat=istat, iomsg=msg)
            if (istat /= 0) then
               write (*, *) 'model_thermodynamic_csp read failed, error message = ', msg
               write (*, *) 'we use default parameters'
            else
               write (*, *) 'read parameters for model_thermodynamic_csp successful'
            end if
         end if

         close (unit=namelist_unit, status='keep', iostat=istat, iomsg=msg)
         if (istat /= 0) then
            write (*, *) 'namelist.nmefc_macom close failed, error message = ', msg
         end if

         CALL mpi_namelist_pack_and_bcast

      ELSE

         CALL mpi_namelist_unpack_and_bcast

      END IF

      DEALLOCATE (mpi_buf, STAT=mpi_mem_status)

   end subroutine misc_namelist_read

!==============================================================================
   subroutine jul2greg(ksec, kminut, khour, kday, kmonth, kyear, prelday, krefdate)
!==============================================================================
    !!-----------------------------------------------------------------------
    !!
    !!                     ***  ROUTINE jul2greg  ***
    !!
    !! ** Purpose : Take the relative time in days and re-express in terms of
    !!              seconds, minutes, hours, days, month, year.
    !!
    !! ** Method  : Reference date : 19650101, 19491001
    !!
    !! ** Action  :
    !!
    !! History
    !!      ! 06-04  (A. Vidard) Original
    !!      ! 06-05  (A. Vidard) Reformatted and refdate
    !!      ! 06-10  (A. Weaver) Cleanup
    !!      ! 2014-09 (D. Lea) Change to use FLOOR to deal with negative prelday
    !!      ! 2020-09 (Yu Zhang) Reference from NEMO
    !!-----------------------------------------------------------------------
      implicit none
      ! * Arguments
      integer, intent(in), optional :: krefdate
      integer, intent(out) :: ksec, kminut, khour, kday, kmonth, kyear
      integer(8), intent(in) :: prelday

    !! * Local declarations
      integer, parameter :: jpgreg = 2299161, jporef = 2433191, jparef = 2438762
      integer :: ijulian, ij1, ija, ijb, ijc, ijd, ije
      integer :: iref, zday, ksec2

      iref = jporef

      ! Main computation
      if (present(krefdate)) then
         select case (krefdate)

         case (0)
            iref = jpgreg

         case (19650101)
            iref = jparef

         case default
            write (*, '(a,i8.8)') 'jul2greg: Unknown krefdate:', krefdate
            stop
         end select
      end if

      zday = int(prelday/86400, 4)

      ksec = int(mod(prelday, 86400), 4)
      ksec2 = abs(ksec)

      if (ksec < 0) ksec = 86400 + ksec

      khour = ksec/3600
      kminut = (ksec - 3600*khour)/60
      ksec = mod(ksec, 60)

      ijulian = iref + zday
      if ((zday <= 0) .and. (prelday < 0)) then
         if (ksec2 == 0) then
            continue
         else
            ijulian = ijulian - 1
         end if
      end if

      ! If input date after 10/15/1582 :
      if (ijulian >= jpgreg) then
         ij1 = int((dble(ijulian - 1867216) - 0.25)/36524.25)
         ija = ijulian + 1 + ij1 - int((0.25*ij1))
      else
         ija = ijulian
      end if

      ijb = ija + 1524
      ijc = int(6680.+(dble(ijb - 2439870) - 122.1)/365.25)
      ijd = 365*ijc + int(0.25*ijc)
      ije = int((ijb - ijd)/30.6001)
      kday = ijb - ijd - int(30.6001*ije)
      kmonth = ije - 1
      if (kmonth > 12) kmonth = kmonth - 12
      kyear = ijc - 4715
      if (kmonth > 2) kyear = kyear - 1
      if (kyear <= 0) kyear = kyear - 1

   end subroutine jul2greg

!==============================================================================
   subroutine greg2jul(ksec, kmin, khour, kday, kmonth, kyear, pjulian, krefdate)
!==============================================================================
    !!-----------------------------------------------------------------------
    !!
    !!                     ***  ROUTINE greg2jul  ***
    !!
    !! ** Purpose : Produce the time relative to the current date and time.
    !!
    !! ** Method  : The units are days, so hours and minutes transform to
    !!              fractions of a day.
    !!
    !!              Reference date : 19650101, 19491001
    !! ** Action  :
    !!
    !! History :
    !!      ! 06-04  (A. Vidard) Original
    !!      ! 06-04  (A. Vidard) Reformatted
    !!      ! 06-10  (A. Weaver) Cleanup
    !!      ! 09-20  (Yu Zhang) Reference from NEMO
    !!-----------------------------------------------------------------------
      implicit none
      ! * Arguments
      integer, intent(in) :: ksec, kmin, khour, kday, kmonth, kyear
      integer(8), intent(out) :: pjulian
      integer, intent(in), optional :: krefdate

    !! * Local declarations
      integer, parameter :: &
         jpgreg = 15 + 31*(10 + 12*1582), &     ! Gregorian calendar introduction date
         jporef = 2433191, &     ! Julian reference date: 19491001
         jparef = 2438762, &     ! Julian reference date: 19650101
         jpgref = 2299161                             ! Julian reference date start of Gregorian calender
      integer :: ija, ijy, ijm, ijultmp, ijyear, iref

      iref = jporef

      if (present(krefdate)) then
         select case (krefdate)

         case (0)
            iref = jpgref

         case (19650101)
            iref = jparef

         case default
            write (*, '(a,i8.8)') 'greg2jul: Unknown krefdate:', krefdate
            stop
         end select
      end if

      ! Main computation
      ijyear = kyear

      if (ijyear < 0) ijyear = ijyear + 1
      if (kmonth > 2) then
         ijy = ijyear
         ijm = kmonth + 1
      else
         ijy = ijyear - 1
         ijm = kmonth + 13
      end if
      ijultmp = int(365.25*ijy) + int(30.6001*ijm) + kday + 1720995
      if (kday + 31*(kmonth + 12*ijyear) >= jpgreg) then
         ija = int(0.01*ijy)
         ijultmp = ijultmp + 2 - ija + int(0.25*ija)
      end if
      pjulian = (int(ijultmp, 8) - int(iref, 8))*86400 + int(((60*khour + kmin)*60 + ksec), 8)

   end subroutine greg2jul

end module mod_misc
