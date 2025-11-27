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