program read_prepbufr
  use prepbufr_constants
  use prepbufr_data
  use prepbufr_utils
  implicit none
  
  integer :: i, ii, jj, kk, lv, jret, mm
  
  ! Get command line arguments
  call getarg(1, argv)
  inf = argv
  call getarg(2, argv)
  outf = argv
  call getarg(3, argv)
  config = argv
  
  ! Initialize variables
  np = 0
  nt = 0
  nplat = 0
  ns = 0
  lon1 = 0.0
  lon2 = 360.0
  lat1 = 90.0
  lat2 = -90.0
  
  ! Print input/output information
  print *, "infile = ", trim(inf)
  print *, "outfile = ", trim(outf)
  print *, "config file = ", trim(config)
  
  ! Read configuration file
  open(unit=10, file=trim(config), form='formatted')
  do i = 1, 10
    read(10, '(A100)', iostat=io) crec
    if (io < 0) exit
    inlength = len_trim(crec) - 4
    call process_config_line(crec(1:4), crec, inlength)
  end do
  close(10)
  
  ! Open output files based on report types
  if (nt > 0) then
    do kk = 1, nt
      do ii = 1, NFILO
        if (type(kk) == filo(ii)) then
          open(unit=iunso(ii), file=trim(outf)//'.'//filo(ii))
          write(unit=iunso(ii), fmt='("#", 148("-"))')
          write(unit=iunso(ii), fmt='("#",a4,a11,a7,a9,a8,a9,a8,a7,a5,a6,8a9)') &
            'SID','XOB','YOB','ELV','DHR','TYP','T29','ITP', &
            'lev','var','OB','QM','PC','RC','FC','AN','OE','CAT'
          write(unit=iunso(ii), fmt='("#", 148("-"))')
          exit
        end if
      end do
    end do
  else
    do ii = 1, NFILO
      open(unit=iunso(ii), file=trim(outf)//'.'//filo(ii))
      write(unit=iunso(ii), fmt='("#", 148("-"))')
      write(unit=iunso(ii), fmt='("#",a4,a11,a7,a9,a8,a9,a8,a7,a5,a6,8a9)') &
        'SID','XOB','YOB','ELV','DHR','TYP','T29','ITP', &
        'lev','var','OB','QM','PC','RC','FC','AN','OE','CAT'
      write(unit=iunso(ii), fmt='("#", 148("-"))')
    end do
  end if
  
  ! Open PREPBUFR input file
  open(unit=11, file=trim(inf), form='unformatted')
  call OPENBF(11, 'IN', 11)
  call DATELEN(10)
  
  ! Main processing loop
  main_loop: do
    ! Get next station report
    call readpb(11, subset, idate, ierrpb)
    if (ierrpb == -1) then
      write(6,*) 'All subsets read in and processed. Exiting.'
      exit main_loop
    end if
    
    ! PREPBUFR data type subsetting filter
    if (nt > 0) then
      found = .false.
      k = 1
      do while ((.not. found) .and. (k <= nt))
        if (subset(1:6) == type(k)) then
          found = .true.
        else
          k = k + 1
        end if
      end do
      
      if ((.not. found) .and. (ierrpb == 0)) then
        cycle main_loop
      end if
      
      if ((.not. found) .and. (ierrpb == 1)) then
        exit main_loop
      end if
    end if
    
    ! Reporting platform subsetting filter
    if (nplat > 0) then
      found = .false.
      k = 1
      do while ((.not. found) .and. (k <= nplat))
        if (int(hdr(7)) == plat(k)) then
          found = .true.
        else
          k = k + 1
        end if
      end do
      
      if ((.not. found) .and. (ierrpb == 0)) then
        cycle main_loop
      end if
      
      if ((.not. found) .and. (ierrpb == 1)) then
        exit main_loop
      end if
    end if
    
    ! Station ID subsetting filter
    if (ns > 0) then
      found = .false.
      k = 1
      write(sid, '(F8.0)') hdr(1)
      sid = adjustl(sid)
      
      do while ((.not. found) .and. (k <= ns))
        if (sid == said(k)) then
          found = .true.
        else
          k = k + 1
        end if
      end do
      
      if ((.not. found) .and. (ierrpb == 0)) then
        cycle main_loop
      end if
      
      if ((.not. found) .and. (ierrpb == 1)) then
        exit main_loop
      end if
    end if
    
    ! Longitude/latitude subsetting filter
    if (ns == 0) then
      found = .false.
      
      ! Case lon1 < lon2
      if (lon1 < lon2) then
        if ((hdr(2) >= lon1) .and. (hdr(2) <= lon2)) then
          if ((hdr(3) <= lat1) .and. (hdr(3) >= lat2)) then
            found = .true.
          end if
        end if
      else
        ! Case lon1 > lon2
        if ((hdr(2) >= lon1) .or. (hdr(2) <= lon2)) then
          if ((hdr(3) <= lat1) .and. (hdr(3) >= lat2)) then
            found = .true.
          end if
        end if
      end if
      
      if ((.not. found) .and. (ierrpb == 0)) then
        cycle main_loop
      end if
      
      if ((.not. found) .and. (ierrpb == 1)) then
        exit main_loop
      end if
    end if
    
    ! Set appropriate output file unit number
    found = .false.
    ii = 1
    do while ((.not. found) .and. (ii <= NFILO))
      if (subset(1:6) == filo(ii)) then
        found = .true.
        iuno = iunso(ii)
      else
        ii = ii + 1
      end if
    end do
    
    if ((.not. found) .and. (ierrpb == 0)) then
      cycle main_loop
    end if
    
    ! Loop through event data array EVNS
    level_loop: do lv = 1, nlev
      var_loop: do kk = 1, MXR8VT
        ! Parameter subsetting filter
        if (np > 0) then
          found = .false.
          p = 1
          do while ((.not. found) .and. (p <= np))
            if (var(kk) == parm(p)) then
              found = .true.
            else
              p = p + 1
            end if
          end do
          
          if (.not. found) then
            cycle var_loop
          end if
        end if
        
        ! Check for virtual temperature
        tvflag = 1
        if ((var(kk) == 'T') .and. (skip_vtmp)) then
          call virtmp(lv, kk, tv_ev_idx, tvflag)
          if (tvflag == -1) then
            cycle var_loop
          end if
        end if
        
        ! Write header and EVNS data to output file
        event_loop: do jj = 1, MXR8VN
          ! Skip virtual temperature at tv_ev_idx
          if ((var(kk) == 'T') .and. (jj <= tv_ev_idx)) then
            cycle event_loop
          end if
          
          write(outstg, '(A8,1X,2F7.2,1X,F8.1,1X,F7.3,1X,F8.1,1X,F7.1,1X,F6.1,1X,I4,1X,A5,8(1X,F8.1))') &
            subset(1:8), (hdr(ii), ii=1,8), lv, var(kk), (evns(ii,lv,jj,kk), ii=1,8)
          
          ! Replace * with spaces
          do mm = 1, 150
            if (outstg(mm:mm) == '*') then
              outstg(mm:mm) = ' '
            end if
          end do
          
          ! Write non-blank lines
          if (outstg(77:137) /= ' ') then
            write(unit=iuno, fmt='(A150)') outstg
          end if
        end do event_loop
      end do var_loop
    end do level_loop
    
    ! Continue if not last report
    if (ierrpb == 0) then
      cycle main_loop
    end if
  end do main_loop
  
  write(6,*) 'All subsets read in and processed. End of program.'
  
end program read_prepbufr