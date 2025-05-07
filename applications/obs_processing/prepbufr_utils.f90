module prepbufr_utils
  use prepbufr_constants
  use prepbufr_data
  implicit none

  private
  public :: readpb, virtmp, process_config_line

contains

  ! Process a configuration file line based on its type
  subroutine process_config_line(line_type, crec, inlength)
    character(len=*), intent(in) :: line_type
    character(len=*), intent(in) :: crec
    integer, intent(in) :: inlength
    
    integer :: j

    select case (line_type)
      case ("LATS")
        read(crec, *) id, lat1, lat2
        print *, "Latitude Values: ", lat1, lat2
      
      case ("LONS")
        read(crec, *) id, lon1, lon2
        print *, "Longitude Values: ", lon1, lon2
      
      case ("SAID")
        ns = inlength / 6
        read(crec, *) id, (said(j), j=1, ns)
        print *, "Stations: ", (trim(said(j))//" ", j=1, ns)
      
      case ("PARM")
        np = inlength / 2
        read(crec, *) id, (parm(j), j=1, np)
        print *, "Parameters: ", (parm(j)//" ", j=1, np)
      
      case ("TYPE")
        nt = inlength / 7
        read(crec, *) id, (type(j), j=1, nt)
        print *, "Report Types: ", (trim(type(j))//" ", j=1, nt)
      
      case ("PLAT")
        nplat = inlength / 4
        read(crec, *) id, (plat(j), j=1, nplat)
        print *, "Reporting platforms: ", (plat(j), j=1, nplat)
      
      case ("VTMP")
        read(crec, *) id, vtmp
        if (vtmp == "FALSE") then
          print *, "*** Skipping virtual temperature ***"
          skip_vtmp = .true.
        else
          print *, "*** Retaining virtual temperature ***"
          skip_vtmp = .false.
        end if
    end select
  end subroutine process_config_line

  ! Read and combine mass and wind subsets from prepbufr file
  subroutine readpb(lunit, subset_out, idate_out, iret)
    integer, intent(in) :: lunit
    character(len=*), intent(out) :: subset_out
    integer, intent(out) :: idate_out
    integer, intent(out) :: iret
    
    ! Local variables
    character(len=MXSTRL) :: head
    character(len=MXSTRL) :: ostr(MXR8VT)
    real(kind=8) :: hdr2(NHR8PM)
    real(kind=8) :: evns2(MXR8PM, MXR8LV, MXR8VN, MXR8VT)
    real(kind=8) :: r8sid, r8sid2, pob1, pob2
    character(len=8) :: csid, csid2, subst2
    logical :: match = .true.
    integer :: jret, ii, jj, kk, lv, lv2, nlev2
    integer :: idate2
    
    ! Initialize header and data strings
    head = 'SID XOB YOB ELV DHR TYP T29 ITP'
    ostr(1) = 'POB PQM PPC PRC PFC PAN POE CAT'
    ostr(2) = 'QOB QQM QPC QRC QFC QAN QOE CAT'
    ostr(3) = 'TOB TQM TPC TRC TFC TAN TOE CAT'
    ostr(4) = 'ZOB ZQM ZPC ZRC ZFC ZAN ZOE CAT'
    ostr(5) = 'UOB WQM WPC WRC UFC UAN WOE CAT'
    ostr(6) = 'VOB WQM WPC WRC VFC VAN WOE CAT'
    
    ! Initialize return code
    iret = 0
    
    ! If previous call did not yield matching subsets, READNS is already
    ! pointing at an unmatched subset. Otherwise, advance to next subset.
    if (match) then
      call READNS(lunit, subset, idate, jret)
      if (jret /= 0) then
        iret = -1
        return
      end if
    else
      subset = subst2
      idate = idate2
    end if
    
    ! Read HDR and EVNS data for current subset
    call UFBINT(lunit, hdr, NHR8PM, 1, jret, head)
    do ii = 1, MXR8VT
      call UFBEVN(lunit, evns(1, 1, 1, ii), MXR8PM, MXR8LV, MXR8VN, nlev, ostr(ii))
    end do
    
    ! Advance to next subset and read its HDR data
    call READNS(lunit, subst2, idate2, jret)
    if (jret /= 0) then
      iret = 1
      subset_out = subset
      idate_out = idate
      return
    end if
    
    call UFBINT(lunit, hdr2, NHR8PM, 1, jret, head)
    
    ! Check if subsets match (same station)
    match = .true.
    
    if (subset /= subst2) then
      match = .false.
      subset_out = subset
      idate_out = idate
      return
    end if
    
    ! Convert real SID to character and compare
    write(csid, '(F8.0)') hdr(1)
    write(csid2, '(F8.0)') hdr2(1)
    
    if (csid /= csid2) then
      match = .false.
      subset_out = subset
      idate_out = idate
      return
    end if
    
    do ii = 2, 5
      if (hdr(ii) /= hdr2(ii)) then
        match = .false.
        subset_out = subset
        idate_out = idate
        return
      end if
    end do
    
    ! Read EVNS data for second matching subset
    do ii = 1, MXR8VT
      call UFBEVN(lunit, evns2(1, 1, 1, ii), MXR8PM, MXR8LV, MXR8VN, nlev2, ostr(ii))
    end do
    
    ! Merge EVNS2 into EVNS
    merge_levels: do lv2 = 1, nlev2
      do lv = 1, nlev
        pob1 = evns(1, lv, 1, 1)
        pob2 = evns2(1, lv2, 1, 1)
        
        if (pob1 == pob2) then
          ! This level exists in both subsets, overwrite missing data in EVNS
          do kk = 1, MXR8VT
            do jj = 1, MXR8VN
              do ii = 1, MXR8PM
                if (evns(ii, lv, jj, kk) == R8BFMS) then
                  evns(ii, lv, jj, kk) = evns2(ii, lv2, jj, kk)
                end if
              end do
            end do
          end do
          cycle merge_levels
          
        else if ((pob2 > pob1) .or. (lv == nlev)) then
          ! Add this level from second subset to EVNS
          nlev = nlev + 1
          do kk = 1, MXR8VT
            do jj = 1, MXR8VN
              do ii = 1, MXR8PM
                evns(ii, nlev, jj, kk) = evns2(ii, lv2, jj, kk)
              end do
            end do
          end do
          cycle merge_levels
        end if
      end do
    end do merge_levels
    
    subset_out = subset
    idate_out = idate
  end subroutine readpb

  ! Check for virtual temperature observation
  subroutine virtmp(lev, k, idx, flag)
    integer, intent(in) :: lev, k
    integer, intent(out) :: idx, flag
    
    integer :: j
    
    idx = 0
    flag = 1
    
    do j = 1, MXR8VN
      if (evns(3, lev, j, k) == VIRTMP_PROG_CODE) then
        idx = j
        
        ! Skip if reason code = 3
        if (evns(4, lev, j, k) == VIRTMP_REASON_CODE) then
          flag = -1
          return
        end if
      end if
    end do
  end subroutine virtmp

end module prepbufr_utils 