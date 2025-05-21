module mitice_ave
!  *===========================================================*
!   used for time-averaging vars giving a period
!   note: only interior grids can be used for averaging (except halo)
!  *===========================================================*

use mitice_vars
use mitice_parameters
use mod_misc_basic
use mod_csp_basic
   implicit none
       real(wp)  ::  SEAICE_timeAve
       real(wp), allocatable, dimension(:)  ::  UTAUtave
       real(wp), allocatable, dimension(:)  ::  VTAUtave
       real(wp), allocatable, dimension(:)  ::  EmPmRtave
       real(wp), allocatable, dimension(:)  ::  QNETtave
       real(wp), allocatable, dimension(:)  ::  QSWtave
       real(wp), allocatable, dimension(:)  ::  UICEtave
       real(wp), allocatable, dimension(:)  ::  VICEtave
       real(wp), allocatable, dimension(:)  ::  HEFFtave
       real(wp), allocatable, dimension(:)  ::  AREAtave
       !$acc declare create(UTAUtave, VTAUtave,          &
       !$acc EmPmRtave, QNETtave, QSWtave,               &
       !$acc UICEtave, VICEtave, HEFFtave, AREAtave)
       
#ifdef SEAICE_ITD
       real(wp), allocatable, dimension(:,:) ::  HEFFITDtave
       real(wp), allocatable, dimension(:,:) ::  AREAITDtave
       !$acc declare create(HEFFITDtave,AREAITDtave)
#endif

contains
   subroutine mitice_ave_allocate
     implicit none 
     allocate(UTAUtave(loc_nlpb),VTAUtave(loc_nlpb))
     allocate(EmPmRtave(loc_nlpb),QNETtave(loc_nlpb),QSWtave(loc_nlpb))
     allocate(UICEtave(loc_nlpb),VICEtave(loc_nlpb))
     allocate(HEFFtave(loc_nlpb),AREAtave(loc_nlpb))
#ifdef SEAICE_ITD
     allocate(HEFFITDtave(loc_nlpb,nITD),AREAITDtave(loc_nlpb,nITD))
#endif


   end subroutine mitice_ave_allocate

   subroutine mitice_ave_release
     implicit none 
     deallocate(UTAUtave,VTAUtave)
     deallocate(EmPmRtave,QNETtave,QSWtave)
     deallocate(UICEtave,VICEtave)
     deallocate(HEFFtave,AREAtave)
#ifdef SEAICE_ITD
     deallocate(HEFFITDtave,AREAITDtave)
#endif
   end subroutine mitice_ave_release


   subroutine mitice_ave_reset
!  *===========================================================*
!    initiate the time-averaged vars to zero
!  *===========================================================*
     implicit none 
!$acc data present(UTAUtave,VTAUtave,EmPmRtave,QNETtave,       &
!$acc        QSWtave,UICEtave,VICEtave,HEFFtave,AREAtave)

!     Initialize averages to zero
     CALL timeave_reset_1d( UTAUtave    )
     CALL timeave_reset_1d( VTAUtave    )
     CALL timeave_reset_1d( EmPmRtave )
     CALL timeave_reset_1d( QNETtave  )
     CALL timeave_reset_1d( QSWtave   )
     CALL timeave_reset_1d( UICEtave  )
     CALL timeave_reset_1d( VICEtave  )
     CALL timeave_reset_1d( HEFFtave  )
     CALL timeave_reset_1d( AREAtave  )
!$acc end data

!not necessary, put timeave_reset to GPU
!!UPLOAD TO DEVICE
!!$ACC UPDATE DEVICE(UTAUtave,VTAUtave,EmPmRtave,QNETtave,QSWtave,   &
!!$ACC              UICEtave,VICEtave,HEFFtave,AREAtave)

#ifdef SEAICE_ITD
     CALL timeave_reset_2d( HEFFITDtave,nITD )
     CALL timeave_reset_2d( AREAITDtave,nITD )
#endif

#ifdef SEAICE_ITD
!$ACC UPDATE DEVICE(HEFFITDtave,AREAITDtave)
#endif

     SEAICE_timeAve = ZERO
   end subroutine mitice_ave_reset

   subroutine mitice_ave_cumulate
!  *===========================================================*
!    sum up seaice vars multiplied by dt
!  *===========================================================*
       implicit none

       call timeave_cumulate_1d(UTAUtave,utau(1:loc_nlpb),SEAICE_deltaTtherm)
       call timeave_cumulate_1d(VTAUtave,vtau(1:loc_nlpb),SEAICE_deltaTtherm)
       call timeave_cumulate_1d(EmPmRtave,EmPmR(1:loc_nlpb),SEAICE_deltaTtherm)
       call timeave_cumulate_1d(QNETtave,QNET(1:loc_nlpb),SEAICE_deltaTtherm)
       call timeave_cumulate_1d(QSWtave,QSW(1:loc_nlpb),SEAICE_deltaTtherm)
       call timeave_cumulate_1d(UICEtave,uIce(1:loc_nlpb),SEAICE_deltaTtherm)
       call timeave_cumulate_1d(VICEtave,vIce(1:loc_nlpb),SEAICE_deltaTtherm)
       call timeave_cumulate_1d(HEFFtave,HEFF(1:loc_nlpb),SEAICE_deltaTtherm)
       call timeave_cumulate_1d(AREAtave,AREA(1:loc_nlpb),SEAICE_deltaTtherm)
#ifdef SEAICE_ITD
       call timeave_cumulate_2d(HEFFITDtave,HEFFITD(1:loc_nlpb,:),nITD,SEAICE_deltaTtherm)
       call timeave_cumulate_2d(AREAITDtave,AREAITD(1:loc_nlpb,:),nITD,SEAICE_deltaTtherm)
#endif
       SEAICE_timeAve = SEAICE_timeAve+SEAICE_deltaTtherm
       
   end subroutine mitice_ave_cumulate

   subroutine mitice_ave_normalize(timeAve)
!  *===========================================================*
!   sum-up value / averaging-period = time-averaged value
!  *===========================================================*
       implicit none
       !interface var
       real(wp)  :: timeAve
       call timeave_normalize_1d( UTAUtave,timeAve    )
       call timeave_normalize_1d( VTAUtave,timeAve    )
       call timeave_normalize_1d( EmPmRtave,timeAve )
       call timeave_normalize_1d( QNETtave,timeAve  )
       call timeave_normalize_1d( QSWtave,timeAve   )
       call timeave_normalize_1d( UICEtave,timeAve  )
       call timeave_normalize_1d( VICEtave, timeAve )
       call timeave_normalize_1d( HEFFtave, timeAve )
       call timeave_normalize_1d( AREAtave, timeAve )
#ifdef SEAICE_ITD
       call timeave_normalize_2d( HEFFITDtave, timeAve, nITD )
       call timeave_normalize_2d( AREAITDtave, timeAve, nITD )
#endif
   end subroutine mitice_ave_normalize
   

   subroutine timeave_cumulate_2d(fldtave,fld, Ksize, deltaTloc)
!  *==========================================================*
!   o Sum over time for 2d vars in MaCOM                           
!  *==========================================================*
      implicit none

!     fldtave - time averaged Field
!     fld  - Input Field
!     Ksize - k dimension
      INTEGER   ::  Ksize
      real(wp)  ::  fld(:,:)
      real(wp)  ::  fldtave(:,:)
      real(wp)  ::  deltaTloc
      INTEGER i, k
 
!$acc kernels loop collapse(2) present(fldtave,fld) copyin(deltaTloc)
      DO k=1,Ksize
        do i = 1,loc_nlpb
          fldtave(i,k)=fldtave(i,k)+fld(i,k)*deltaTloc
        enddo
      ENDDO
!$acc end kernels
 
   end subroutine timeave_cumulate_2d

   subroutine timeave_cumulate_1d(fldtave,fld, deltaTloc)
!  *==========================================================*
!    sum over time for 1d vars in MaCOM
!  *==========================================================*
      implicit none

!     fldtave - time averaged Field
!     fld  - Input Field
      real(wp)  ::   fld(:)
      real(wp)  ::   fldtave(:)
      real(wp)  ::   deltaTloc
      INTEGER i
 
!$acc kernels loop present(fldtave,fld) copyin(deltaTloc)
      do i = 1,loc_nlpb
        fldtave(i)=fldtave(i)+fld(i)*deltaTloc
      enddo
!$acc end kernels
 
   end subroutine timeave_cumulate_1d

   subroutine timeave_normalize_2d(fldtave,timeave_cumul, Ksize)
!  *==========================================================*
!   divided by averaging period (must be a multiple of deltaT)
!  *==========================================================*

      implicit none
!     fldtave       :: time averaged Field
!     timeave_cumul :: cumulated time for average
!     Ksize - k dimension
      integer  ::  Ksize
      real(wp) ::  fldtave(:,:)
      real(wp) ::  timeave_cumul
      integer i, k

      IF ( timeave_cumul .NE. 0.0_wp ) THEN
!$acc kernels present(fldtave) copyin(timeave_cumul)
!$acc loop collapse(2)
         DO k=1,Ksize
           do i = 1,loc_nlpb
             fldtave(i,k) = fldtave(i,k) / timeave_cumul
           enddo
         ENDDO
!$acc end kernels
      ENDIF

   end subroutine timeave_normalize_2d

   subroutine timeave_normalize_1d(fldtave,timeave_cumul)
!  *==========================================================*
!   divided by averaging period (must be a multiple of deltaT)
!  *==========================================================*

      implicit none
      real(wp) ::  fldtave(:)
      real(wp) ::  timeave_cumul
      integer i

      IF ( timeave_cumul .NE. 0.0_wp ) THEN
!$acc kernels present(fldtave) copyin(timeave_cumul)
!$acc loop
         do i = 1,loc_nlpb
           fldtave(i) = fldtave(i) / timeave_cumul
         enddo
!$acc end kernels
      ENDIF

   end subroutine timeave_normalize_1d

   subroutine timeave_reset_1d(fldtave)
!  *==========================================================*
!    Initialize arrays for averaging  (Ksize=1 for 2dvar)
!  *==========================================================*
      implicit none

!     fldtave - time averaged Field
      real(wp) :: fldtave(:)
      integer i

!$acc kernels loop present(fldtave) 
      do i = 1,loc_nlpb
        fldtave(i) = 0.0_wp
      enddo
!$acc end kernels
   end subroutine timeave_reset_1d

   subroutine timeave_reset_2d(fldtave,Ksize)
!  *==========================================================*
!    Initialize arrays for averaging  (Ksize=1 for 2dvar)
!  *==========================================================*
      implicit none

      integer  ::  Ksize
!     fldtave - time averaged Field
      real(wp) :: fldtave(:,:)
      integer i,k
 
!$acc kernels loop collapse(2) present(fldtave) 
      do k = 1,Ksize
        do i = 1,loc_nlpb
          fldtave(i,k) = 0.0_wp
        enddo
      enddo
!$acc end kernels
   end subroutine timeave_reset_2d
 
end module mitice_ave
