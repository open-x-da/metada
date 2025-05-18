      INCLUDE         'readpb.prm'

      PARAMETER (STRLN = 180)
      parameter (maxtype=15)
      parameter (maxplat=15)
      parameter	(maxparm=5)
      parameter	(maxsaid=15)
C*
      CHARACTER  outstg*(150), subset*8, 
     +           inf*300, outf*300, config*300, argv*300,
     +           crec*101, type(maxtype)*6,
     +  		 parm(maxparm),id*4,idatec*10,
     +           said(maxsaid)*5,sid*5,vtmp*5
     
      LOGICAL skip_vtmp
 
C*
      INTEGER    plat(maxplat)
      
      PARAMETER  ( NFILO = 15 ) 
      INTEGER    iunso ( NFILO )
     +		/   51,   52,   53,   54,   55,
     +		    56,   57,   58,   59,   60,
     +		    61,   62,   63,   64,   65  /

	  CHARACTER*6 filo ( NFILO )
     +		/ 'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR',
     +		  'VADWND', 'SATBOG', 'SATEMP', 'ADPSFC', 'SFCSHP',
     +		  'SFCBOG', 'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND'  /

      CHARACTER  var ( MXR8VT )
     +        /'P','Q','T','Z','U','V'/

      LOGICAL    found 

      INTEGER    io,stat,n,inlength,np,nt,nplat,
     +           count,k,flag,pflag,p,ns,s,platflag
	
   	  REAL*8	 lat1,lat2,lon1,lon2,platform
        
      INTEGER    tv_ev_idx, tvflag

C*
C-----------------------------------------------------------------------
C*      Input arguments:
C*      1. PREPBUFR input file (path + file name)
C*      2. output file (path + file name)
C*      3. configuration file (path + file name)

      call getarg(1, argv)
      inf=argv
      call getarg(2, argv)
      outf=argv
      call getarg(3, argv)
      config=argv
  
      np=0
      nt=0
      nplat=0
      ns=0
      lon1=0.0
      lon2=360.0
      lat1=90.0
      lat2=-90.0
      
c* Omit virtual temperature in output (default)
      skip_vtmp = .false.

      print*, "infile = ",inf
      print*, "outfile = ",outf
      print*, "config file = ",config

c*	open the configuration file
      open (unit=10, file=config, form='formatted')
      do i = 1,10
        read(10, '(A100)', iostat=io) crec
        if (io < 0) exit
        inlength=len_trim(crec)-4
        select case (crec(1:4))
          case ("LATS")
            read (crec,*) id,lat1,lat2
            print *,"Latitude Values: ",lat1,lat2
          case ("LONS")
            read (crec,*) id,lon1,lon2
            print *,"Longitude Values: ",lon1,lon2
          case ("SAID")
            ns=inlength/6
            read (crec,*) id,(said(j), j=1,ns)
            print *,"Stations: ",(said(j)," ", j=1,ns)
          case ("PARM")
            np=inlength/2
            read (crec,*) id,(parm(j), j=1,np)
            print *,"Parameters: ",(parm(j)," ", j=1,np)
          case ("TYPE")
            nt=inlength/7
            read (crec,*) id,(type(j), j=1,nt)
            print *,"Report Types: ",(type(j)," ", j=1,nt)
          case ("PLAT")
            nplat=inlength/4
            read (crec,*) id,(plat(j), j=1,nplat)
            print *,"Reporting platforms: ",(plat(j)," ", j=1,nplat)
          case ("VTMP")
            read (crec,*) id,vtmp
            if (vtmp .eq. "FALSE") then
               print *,"*** Skipping virtual temperature ***"
               skip_vtmp = .true.
            else
               print *,"*** Retaining virtual temperature ***"
               skip_vtmp = .false.
            end if
        end select
      end do

c-----7---------------------------------------------------------------72
C*    Open the output file(s)
c-----7---------------------------------------------------------------72
	  if(nt .gt. 0) then
	    do kk = 1, nt
	    do ii = 1, NFILO
	      if (type(kk) .eq. filo(ii)) then
	        open (unit=iunso(ii), file=trim(outf) // '.' // filo(ii))
            WRITE (UNIT=iunso(ii), FMT=15)
            WRITE (UNIT=iunso(ii), FMT=20)
     +      'SID','XOB','YOB','ELV','DHR','TYP','T29','ITP',
     +      'lev','var','OB','QM', 'PC', 'RC', 'FC','AN','OE','CAT'
            WRITE (UNIT=iunso(ii), FMT=15)
	        exit
	      end if
	    end do
	    end do
	  else
	    do ii = 1, NFILO
	      open (unit=iunso(ii), file=trim(outf) // '.' // filo(ii))
            WRITE (UNIT=iunso(ii), FMT=15)
            WRITE (UNIT=iunso(ii), FMT=20)
     +      'SID','XOB','YOB','ELV','DHR','TYP','T29','ITP',
     +      'lev','var','OB','QM', 'PC', 'RC', 'FC','AN','OE','CAT'
            WRITE (UNIT=iunso(ii), FMT=15)
	    end do	  
	  end if
  15  FORMAT ("#", 148("-"))
  20  FORMAT ("#",a4,a11,a7,a9,a8,a9,a8,a7,a5,a6,8a9)
	  
c-----7---------------------------------------------------------------72
C*    Open the PREPBUFR input file
c-----7---------------------------------------------------------------72
      OPEN ( UNIT = 11, FILE = inf, FORM = 'UNFORMATTED' )
      CALL OPENBF  ( 11, 'IN', 11 )
      CALL DATELEN  ( 10 )

c-----7---------------------------------------------------------------72
C*   Get the next station report from the input file.
c-----7---------------------------------------------------------------72
  10  CALL READPB  ( 11, subset, idate, ierrpb )
      IF ( ierrpb .eq. -1 )  THEN
        WRITE(6,*) 'All subsets read in and processed. Exiting.'
        STOP
      END IF

c-----7---------------------------------------------------------------72
C*	PREPBUFR data type subsetting filter
c-----7---------------------------------------------------------------72
      if (nt .gt. 0) then
        k = 1
        found = .false.
        do while ((.not. found) .and. (k .le. nt))
          if(subset(1:6) .eq. type(k)) then
            found = .true.
          else
            k = k + 1
          end if
        end do
	    if ((.not. found) .and. (ierrpb .eq. 0)) then
          go to 10
        end if
c       End program if no subset match found and we are processing
c       the last report in the prepbufr file (ierrpb = 1)
        if ((.not. found) .and. (ierrpb .eq. 1)) then
          go to 9998
        end if
      end if

c-----7---------------------------------------------------------------72
c    Reporting platform (input report type) subsetting filter
c-----7---------------------------------------------------------------72
      if (nplat .gt. 0) then
        k = 1
        found = .false.
        do while ((.not. found) .and. (k .le. nplat))
          if(hdr(7) .eq. plat(k)) then
            found = .true.
          else
            k = k + 1
          end if
        end do
	    if ((.not. found) .and. (ierrpb .eq. 0)) then
          go to 10
        end if
c       End program if no subset match found and we are processing
c       the last report in the prepbufr file (ierrpb = 1)
        if ((.not. found) .and. (ierrpb .eq. 1)) then
          go to 9998
        end if
      end if

c-----7---------------------------------------------------------------72
c     Station ID subsetting filter
c-----7---------------------------------------------------------------72
      if (ns .gt. 0) then
        k = 1
        found = .false.       
        write(unit=sid,fmt='(a5)') hdr(1)
        do while ((.not. found) .and. (k .le. ns))
          if(sid .eq. said(k)) then
            found = .true.
          else
            k = k + 1
          end if
        end do
	    if ((.not. found) .and. (ierrpb .eq. 0)) then
          go to 10
        end if
c       End program if no subset match found and we are processing
c       the last report in the prepbufr file (ierrpb = 1)
        if ((.not. found) .and. (ierrpb .eq. 1)) then
          go to 9998
        end if
      end if
 
c-----7---------------------------------------------------------------72
c     Longitude/latitude subsetting filters.  This is ignored if
c     station ID subsetting is selected.
c-----7---------------------------------------------------------------72
      if (ns .eq. 0) then
        found = .false.
C*      Case lon1 < lon2
        if(lon1 .lt. lon2) then 
          if ((hdr(2) .ge. lon1) .and. (hdr(2) .le. lon2)) then 
            if ((hdr(3) .le. lat1) .and. (hdr(3) .ge. lat2)) then
              found = .true.
            end if
          end if
        else
C*      Case lon1 > lon2
          if ((hdr(2) .ge. lon1) .or. (hdr(2) .le. lon2)) then
            if ((hdr(3) .le. lat1) .and. (hdr(3) .ge. lat2)) then
              found = .true.
            end if
          end if
        end if
	    if ((.not. found) .and. (ierrpb .eq. 0)) then
          go to 10
        end if
c       End program if no subset match found and we are processing
c       the last report in the prepbufr file (ierrpb = 1)
        if ((.not. found) .and. (ierrpb .eq. 1)) then
          go to 9998
        end if
      end if
      
c-----7---------------------------------------------------------------72
C*	Set the appropriate output file unit number.
c-----7---------------------------------------------------------------72
	  ii = 1
	  found = .false.
	  DO WHILE ((.not. found) .and. (ii .le. NFILO))
	    IF (subset (1:6) .eq. filo (ii))  THEN
		  found = .true.
		  iuno = iunso(ii)
	    ELSE 
		  ii = ii + 1
	    END IF
	  END DO
	  IF ((.not. found) .and. (ierrpb .eq. 0)) THEN
	    GO TO 10
	  END IF

c-----7---------------------------------------------------------------72
c     Loop through the event data array EVNS
c-----7---------------------------------------------------------------72
      DO lv = 1, nlev
      DO kk = 1, MXR8VT

c-----7---------------------------------------------------------------72
c     Parameter subsetting filter (P, Q, T, Z, U, V)
c-----7---------------------------------------------------------------72
        if(np .gt. 0) then
          p = 1
          found = .false.
          do while ((.not. found) .and. (p .le. np))
            if(var(kk) .eq. parm(p)) then
              found = .true.
            else
              p = p + 1
            end if
          end do
          if (.not. found) then
            cycle
          end if
        end if

c-----7---------------------------------------------------------------72
c   Check for virtual temperature
c-----7---------------------------------------------------------------72
        tvflag=1
        if ((var(kk) .eq. 'T') .and. (skip_vtmp)) then
          call virtmp(lv, kk, tv_ev_idx, tvflag)
          if (tvflag .eq. -1) then
            cycle
          end if
        end if
                
c-----7---------------------------------------------------------------72                
C*    Write the header and EVNS data for this station report to the
C*    output file.
c-----7---------------------------------------------------------------72

        DO jj = 1, MXR8VN
C*        Skip virtual temperature at tv_ev_idx
          if ((var(kk) .eq. 'T') .and. (jj .le. tv_ev_idx)) then
            cycle
          endif
          WRITE (UNIT=outstg, FMT=500)
     +           (hdr (ii), ii = 1, 8),
     +           lv, var(kk), 
     +           (evns(ii,lv,jj,kk),ii=1,8)

          DO mm = 1, 150
            IF (outstg (mm:mm) .eq. '*') THEN
              outstg (mm:mm) = ' '
            END IF
          END DO
          if (outstg(77:137) .ne. ' ') then
            WRITE (UNIT=iuno, FMT='(A150)')  outstg
          endif
        END DO  ! End jj = 1, MXR8VN loop
      END DO  ! End kk = 1, MXR8VT loop
      END DO  ! End lv = 1, nlev loop

c-----7---------------------------------------------------------------72
      IF  ( ierrpb .eq. 0 )  THEN
        GO TO 10
      END IF

 9998 CONTINUE
      WRITE(6,*) 'All subsets read in and processed. End of program.'
        
C*    Format specifier for outstg
  500 FORMAT (A8, 1X, 2F7.2, 1X, F8.1, 1X, F7.3, 1X, F8.1, 1X, F7.1, 
     + 1X, F6.1, 1X, I4, 1X, A5, 8(1X,F8.1))
C* 
      STOP
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        SUBROUTINE READPB  ( lunit, subset, idate, iret )
C
C*      This subroutine will read and combine the mass and wind subsets
C*      of the next station report in the prepbufr file.  It is styled
C*      after entry point READNS, and it only requires the prepbufr file
C*      to be opened for reading with OPENBF.  The combined station
C*      report is returned to the caller in COMMON /PREPBC/.
C*      This common area contains the number of levels in the report,
C*      a one dimensional array with the header information, and a four
C*      dimensional array containing all events from the variables POB,
C*      QOB, TOB, ZOB, UOB, and VOB for the report.
C*
C*      The header array contains the following list of mnemonics:
C*
C*         HDR(1)  Station identification (SID)
C*         HDR(2)  Longitude (XOB)
C*         HDR(3)  Latitude  (YOB)
C*         HDR(4)  Elevation (ELV)
C*         HDR(5)  Observation time minus cycle time (DHR)
C*         HDR(6)  PREPBUFR report type (TYP)
C*         HDR(7) Input report type (T29)
C*         HDR(8) Instrument type (ITP)
C*
C*      The 4-D array of data, EVNS ( ii, lv, jj, kk ), is indexed
C*      as follows:
C*
C*      "ii" indexes the event data types; these consist of:
C*          1) OBservation        (e.g., POB, ZOB, UOB, VOB, TOB, QOB, PWO)
C*          2) Quality Mark       (e.g., PQM, ZRM, WQM, TQM, QQM, PWQ)
C*          3) Program Code       (e.g., PPC, ZPC, WPC, TPC, QPC, PWP)
C*          4) Reason Code        (e.g., PRC, ZRC, WRC, TRC, QRC, PWR)
C*          5) ForeCast value     (e.g., PFC, ZFC, UFC, VFC, TFC, QFC, PWF)
C*          6) ANalysed value     (e.g., PAN, ZAN, UAN, VAN, TAN, QAN, PWA)
C*          7) Observation Error  (e.g., POE, ZOE, WOE, TOE, QOE, PWO)
C*          8) PREPBUFR data level category (CAT)
C*      "lv" indexes the levels of the report
C*          1) Lowest level
C*      "jj" indexes the event stacks
C*          1) N'th event
C*          2) (N-1)'th event (if present)
C*          3) (N-2)'th event (if present)
C*                ...
C*         10) (N-9)'th event (if present)
C*      "kk" indexes the variable types
C*          1) Pressure
C*          2) Specific humidity
C*          3) Temperature
C*          4) Height
C*          5) U-component wind
C*          6) V-component wind
C*
C*      Note that the structure of this array is identical to one
C*      returned from UFBEVN, with an additional (4th) dimension to
C*      include the six variable types into the same array.
C*
C*      The return codes are as follows:
C*      iret =  0 - normal return
C*           =  1 - the station report within COMMON /PREPBC/ contains the
C*                  last available subset from within the prepbufr file
C*           = -1 - there are no more subsets available from within the
C*                  prepbufr file       
C*
        INCLUDE         'readpb.prm'
C*
        CHARACTER*(*)   subset
C* 
        CHARACTER*(MXSTRL)      head
   
        CHARACTER*(MXSTRL)      ostr ( MXR8VT ) 
       
        DATA head  / 'SID XOB YOB ELV DHR TYP T29 ITP' /
C*
        DATA ostr / 'POB PQM PPC PRC PFC PAN POE CAT',
     +              'QOB QQM QPC QRC QFC QAN QOE CAT',
     +              'TOB TQM TPC TRC TFC TAN TOE CAT',
     +              'ZOB ZQM ZPC ZRC ZFC ZAN ZOE CAT',
     +              'UOB WQM WPC WRC UFC UAN WOE CAT',
     +              'VOB WQM WPC WRC VFC VAN WOE CAT'  /
C*
        REAL*8          hdr2 ( NHR8PM ),
     +                  evns2 ( MXR8PM, MXR8LV, MXR8VN, MXR8VT )
C*
        REAL*8          r8sid, r8sid2, pob1, pob2
C*
        CHARACTER*8     csid, csid2, subst2
C*
        LOGICAL         match 
        
        DATA match / .true. /
C*
        EQUIVALENCE     ( r8sid, csid ), ( r8sid2, csid2 )
C*
        SAVE            match, subst2, idate2
C-----------------------------------------------------------------------
        iret = 0
C*
C*      If the previous call to this subroutine did not yield matching
C*      mass and wind subsets, then READNS is already pointing at an
C*      unmatched subset.  Otherwise, call READNS to advance the subset
C*      pointer to the next subset.
C*
        IF  ( match )  THEN
            CALL READNS  ( lunit, subset, idate, jret )
            IF  ( jret .ne. 0 )  THEN
                iret = -1
                RETURN
            END IF
        ELSE
            subset = subst2
            idate = idate2
        END IF
C*
C*      Read the HDR and EVNS data for the subset that is currently
C*      being pointed to.
C*
        CALL UFBINT  ( lunit, hdr, NHR8PM, 1, jret, head )
        DO ii = 1, MXR8VT
            CALL UFBEVN  ( lunit, evns ( 1, 1, 1, ii ), MXR8PM, MXR8LV,
     +                     MXR8VN, nlev, ostr (ii) )
        END DO
C
C*      Now, advance the subset pointer to the following subset and
C*      read its HDR data.
C
        CALL READNS  ( lunit, subst2, idate2, jret )
        IF  ( jret .ne. 0 )  THEN
            iret = 1
            RETURN
        END IF
        CALL UFBINT  ( lunit, hdr2, NHR8PM, 1, jret, head )
C 
C*      Check whether these two subsets have identical SID, YOB, XOB,
C*      ELV, and DHR values.  If so, then they are matching mass and
C*      wind subsets for a single station report.
C
        match = .true.
C
        IF  ( subset .ne. subst2 )  THEN
            match = .false.
            RETURN
        END IF
C 
        r8sid = hdr (1)
        r8sid2 = hdr2 (1)
        IF  ( csid .ne. csid2 )  THEN
            match = .false.
            RETURN
        END IF
C 
        DO ii = 2, 5
            IF  ( hdr (ii) .ne. hdr2 (ii) )  THEN
                match = .false.
                RETURN
            END IF
        END DO
C
C*      Read the EVNS data for the second of the two matching subsets.
C 
        DO ii = 1, MXR8VT
            CALL UFBEVN  ( lunit, evns2 ( 1, 1, 1, ii ), MXR8PM, MXR8LV,
     +                     MXR8VN, nlev2, ostr (ii) )
        ENDDO
C
C*      Combine the EVNS data for the two matching subsets into a
C*      single 4-D array.  Do this by merging the EVNS2 array into
C*      the EVNS array.
C
        DO 10 lv2 = 1, nlev2
            DO lv = 1, nlev
                pob1 = evns ( 1, lv, 1, 1 )
                pob2 = evns2 ( 1, lv2, 1, 1 )
                IF  ( pob1 .eq. pob2 )  THEN
C
C*                This pressure level from the second subset also exists
C*                in the first subset, so overwrite any "missing" piece
C*                of data for this pressure level in the first subset
C*                with the corresponding piece of data from the second
C*                subset (since this results in no net loss of data!).
C
                  DO kk = 1, MXR8VT
                    DO jj = 1, MXR8VN
                      DO ii = 1, MXR8PM
                        IF  ( evns ( ii, lv, jj, kk ) .eq. R8BFMS ) THEN
                          evns ( ii, lv, jj, kk ) =
     +                          evns2 ( ii, lv2, jj, kk )
                        END IF
                      END DO
                    END DO
                  END DO
                  GO TO 10
                ELSE IF  (  ( pob2 .gt. pob1 )  .or.
     +                       ( lv .eq. nlev )  )  THEN
C
C*                Either all remaining pressure levels within the first
C*                subset are less than this pressure level from the
C*                second subset (since levels within each subset are
C*                guaranteed to be in descending order wrt pressure!)
C*                *OR* there are more total levels within the second
C*                subset than in the first subset.  In either case, we
C*                should now add this second subset level to the end of
C*                the EVNS array.
C
                  nlev = nlev + 1
                  DO kk = 1, MXR8VT
                    DO jj = 1, MXR8VN
                      DO ii = 1, MXR8PM
                        evns ( ii, nlev, jj, kk ) =
     +                        evns2 ( ii, lv2, jj, kk )
                      END DO
                    END DO
                  END DO
                  GOTO 10
                END IF
            END DO
   10   END DO
C* 
        RETURN
        END

c----7----------------------------------------------------------------72
      SUBROUTINE virtmp(lev,k,idx,flag)
c
c   // Do not write virtual temperature observations.
c   // PREPBUFR Table 14 describes the VIRTMP processing step:
c   //    http://www.emc.ncep.noaa.gov/mmb/data_processing/prepbufr.doc/table_14.htm
c   // 
c   // For VIRTMP program code 8 with reason code 3, do not use this
c   // this observation of virtual temperature.
c   // For VIRTMP program code 8 with any other reason code, step down the
c   // event stack index to find sensible temperature.
c   //

      INCLUDE 'readpb.prm'

      PARAMETER (virtmp_prog_code = 8.0)
      PARAMETER (virtmp_reason_code = 3.0)
      INTEGER lev,j,k
      INTEGER idx, flag

      idx = 0
      flag = 1
      
      do j = 1, MXR8VN
        if (evns(3,lev,j,k) .eq. virtmp_prog_code) then
          idx = j
                    
c Skip if reason code = 3
          if (evns(4,lev,j,k) .eq. virtmp_reason_code) then
            flag = -1
            return
          endif
                    
        endif
      enddo
      
      return
      end
      