C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        SUBROUTINE READPB  ( lunit, subset, idate, hdr_out, iret )
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
        CHARACTER*(8)   subset
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

        REAL*8          hdr_out ( NHR8PM )  ! Output parameter for header data
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

C*      Copy the header data to the output parameter
        DO ii = 1, NHR8PM
            hdr_out(ii) = hdr(ii)
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