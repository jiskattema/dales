!> \file modtimestat.f90
!!  Timestat calculates timeseries of several variables

!>
!! Timestat calculates timeseries of several variables
!>
!! Timeseries of the most relevant parameters. Written to tmser1.expnr and tmsurf.expnr
!! If netcdf is true, this module leads the tmser.expnr.nc output
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!


module modtimestat


  use modglobal, only : longint

implicit none
! private
! PUBLIC :: inittimestat, timestat
save
!NetCDF variables
  !integer,parameter :: nvar = 28
  integer :: nvar
  integer :: ncid,nrec = 0
  character(80) :: fname = 'tmser.xxx.nc'
  !character(80),dimension(nvar,4) :: ncname
  character(80), allocatable, dimension(:,:)    :: ncname
  character(40) :: name

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: ltimestat= .false. !<switch for timestatistics (on/off)
  real    :: zi,ziold=-1, we
  integer, parameter :: iblh_flux = 1, iblh_grad = 2, iblh_thres = 3
  integer, parameter :: iblh_thv = -1,iblh_thl = -2, iblh_qt = -3
  integer :: iblh_meth = iblh_grad, iblh_var = iblh_thv
  integer :: blh_nsamp = 4
  real    :: blh_thres=-1 ,blh_sign=1.0
  real   :: zbaseav, ztopav, ztopmax,zbasemin
  real   :: qlintav, qlintmax, tke_tot
  real   :: cc, wmax, qlmax
  real   :: qlint
  logical:: store_zi = .false.
  integer:: nvar0
  real ::   cl_c,ci_c,c_totc                                       & !#sb3 START variables for mixed-phase microphysics 
           ,zt_cl_av,zt_ci_av,zt_cl_max,zt_ci_max                  &
           ,zb_cl_av,zb_ci_av,zb_cl_min,zb_ci_min                  &
           ,qcl_intav,qci_intav,qhr_intav,qhs_intav,qhg_intav      &
           ,qcl_intmax,qci_intmax,qhr_intmax,qhs_intmax,qhg_intmax &
           ,sfc_prec_av, sfc_precw_av                               !#sb3 END
  !Variables for heterogeneity
  real, allocatable :: u0av_patch (:,:)     ! patch averaged um    at full level
  real, allocatable :: v0av_patch (:,:)     ! patch averaged vm    at full level
  real, allocatable :: w0av_patch (:,:)     ! patch averaged wm    at full level
  real,allocatable, dimension(:,:) :: zbase_field, ztop_field, cc_field, qlint_field, tke_tot_field
  real,allocatable, dimension(:,:) :: zbase_patch, ztop_patch, zbasemin_patch, zbasemin_patchl
  real,allocatable, dimension(:,:) :: cc_patch, qlint_patch, qlintmax_patch, qlintmax_patchl, tke_tot_patch
  real,allocatable, dimension(:,:) :: wmax_patch, wmax_patchl, qlmax_patch, qlmax_patchl, ztopmax_patch, ztopmax_patchl
  real,allocatable, dimension(:,:) :: ust_patch, qst_patch, tst_patch, wthls_patch, wqls_patch, wthvs_patch
  !In combination with isurf = 1
  real,allocatable, dimension(:,:) :: Qnet_patch, H_patch, LE_patch, G0_patch, tendskin_patch,rs_patch,ra_patch
  real,allocatable, dimension(:,:) :: cliq_patch, wl_patch, rsveg_patch, rssoil_patch, tskin_patch, obl_patch
  real,allocatable, dimension(:,:) :: zi_patch,ziold_patch,we_patch, zi_field


contains
!> Initializing Timestat. Read out the namelist, initializing the variables
  subroutine inittimestat
    use mpi
    use modmpi,    only : my_real,myid,comm3d,mpi_logical,mpierr,mpi_integer
    use modglobal, only : ifnamopt, fname_options,cexpnr,dtmax,ifoutput,dtav_glob,tres,&
                          ladaptive,k1,kmax,rd,rv,dt_lim,btime,i1,j1,lwarmstart,       &
                          checknamelisterror, eps1 ! #sb3
    use modfields, only : thlprof,qtprof,svprof
    use modsurfdata, only : isurf, lhetero, xpatches, ypatches
    use modstat_nc, only : lnetcdf, open_nc, define_nc, ncinfo, nctiminfo
    use modmicrodata, only : imicro,imicro_bulk3               ! #sb3
    use modmicrodata3, only: iq_cl, iq_ci, iq_hr, iq_hs, iq_hg  ! #sb3
    implicit none
    integer :: ierr,k,location = 1
    real :: gradient = 0.0
    real, allocatable,dimension(:) :: profile
    integer :: i,j
    character(len=1000) :: line

    namelist/NAMTIMESTAT/ & !< namelist
    dtav,ltimestat,blh_thres,iblh_meth,iblh_var,blh_nsamp !! namelist contents

    dtav=dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMTIMESTAT,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMTIMESTAT')
      write(6 ,NAMTIMESTAT)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ltimestat  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(blh_thres  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(iblh_meth  ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iblh_var   ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(blh_nsamp  ,1,MPI_INTEGER,0,comm3d,mpierr)
    idtav = dtav/tres

    tnext = idtav+btime

    if(.not.(ltimestat)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'TIMESTAT: dtav should be a integer multiple of dtmax'
    end if

    allocate(profile(1:k1))
    select case (iblh_var)
    case(iblh_qt)
      profile = qtprof
    case(iblh_thl)
      profile = thlprof
    case(iblh_thv)
      do k=1,k1
        profile(k) = thlprof(k)*(1+(rv/rd-1)*qtprof(k))
      end do
    case(1:)
      profile = svprof(:,iblh_var)
    end select
    blh_sign = sign(1.0,profile(kmax)-profile(1))

    select case(iblh_meth)
    case (iblh_flux)
    case (iblh_grad)
    case (iblh_thres)
      if (blh_thres<0) then
        do k=kmax,2,-1
          if (blh_sign*(profile(k+1) - profile(k-1)) > gradient) then
            location = k
            gradient = blh_sign*(profile(k+1) - profile(k-1))
          endif
        enddo
        blh_thres=profile(location)
        if (myid==0) write (*,*) 'TIMESTAT: blh_tres =',blh_thres
      end if
    case default
      stop 'TIMESTAT: Incorrect iblh_meth'
    end select
    deallocate(profile)

    if(myid==0) then
       if (.not. lwarmstart) then
          !tmser1
          open (ifoutput,file='tmser1.'//cexpnr,status='replace',position='append')
          write(ifoutput,'(2a)') &
               '#  time      cc     z_cbase    z_ctop_avg  z_ctop_max      zi         we', &
               '   <<ql>>  <<ql>>_max   w_max   tke     ql_max'
          close(ifoutput)
          !tmsurf
          open (ifoutput,file='tmsurf.'//cexpnr,status='replace',position='append')
          write(ifoutput,'(2a)') &
               '#  time        ust        tst        qst         obukh', &
               '      thls        z0        wthls      wthvs      wqls '
          close(ifoutput)
          if(isurf == 1) then
             open (ifoutput,file='tmlsm.'//cexpnr,status='replace',position='append')
             write(ifoutput,'(3a)') &
                  '#     time      Qnet        H          LE         G0  ', &
                  '   tendskin     rs         ra        tskin        cliq  ', &
                  '    Wl          rssoil     rsveg       Resp       wco2         An', &
                  '    gcco2'
             write(ifoutput,'(3a)') &
                  '#      [s]     [W/m2]     [W/m2]     [W/m2]     [W/m2]', &
                  '   [W/m2]      [s/m]       [s/m]     [K]          [-]   ', &
                  '   [m]          [s/m]      [s/m]   [mgCm2/s]               [mgCm2/s]',&
                  '   [m/s]  '
             close(ifoutput)
          end if
          
          if(lhetero) then
             do i=1,xpatches
                do j=1,ypatches
                   name = 'tmser1patchiiixjjj.'//cexpnr
                   write (name(12:14),'(i3.3)') i
                   write (name(16:18),'(i3.3)') j
                   open (ifoutput,file=name,status='replace',position='append')
                   write(ifoutput,'(2a)') &
                        '#  time      cc     z_cbase    z_ctop_avg  z_ctop_max      zi         we', &
                        '   <<ql>>  <<ql>>_max   w_max   tke     ql_max'
                   close(ifoutput)
                   
                   name = 'tmsurfpatchiiixjjj.'//cexpnr
                   write (name(12:14),'(i3.3)') i
                   write (name(16:18),'(i3.3)') j
                   open (ifoutput,file=name,status='replace',position='append')
                   write(ifoutput,'(2a)') &
                        '#  time        ust        tst        qst         obukh', &
                        '      thls        z0        wthls      wthvs      wqls '
                   close(ifoutput)
                   
                   if(isurf == 1) then
                      name = 'tmlsmpatchiiixjjj.'//cexpnr
                      write (name(11:13),'(i3.3)') i
                      write (name(15:17),'(i3.3)') j
                      open (ifoutput,file=name,status='replace',position='append')
                      write(ifoutput,'(3a)') &
                           '#     time      Qnet        H          LE         G0  ', &
                           '   tendskin     rs         ra        tskin        cliq  ', &
                           '    Wl          rssoil     rsveg'
                      write(ifoutput,'(3a)') &
                           '#      [s]     [W/m2]     [W/m2]     [W/m2]     [W/m2]', &
                           '   [W/m2]      [s/m]       [s/m]     [K]          [-]   ', &
                           '   [m]          [s/m]      [s/m]'
                      close(ifoutput)
                   endif
                enddo
             enddo
          endif
       endif

      if (lnetcdf) then
        if(isurf == 1) then
          nvar = 32
        else
          nvar = 21
        end if
        if (imicro.eq.imicro_bulk3)then ! #sb3 START
           nvar0 = nvar
           nvar = nvar0 +23  ! +21 number of statistics for 3-phase SB microphysics
        endif ! #sb3 END


        allocate(ncname(nvar,4))

        fname(7:9) = cexpnr
        call nctiminfo(ncname(1,:))
        call ncinfo(ncname( 2,:),'cfrac','Cloud fraction','-','time')
        call ncinfo(ncname( 3,:),'zb','Cloud-base height','m','time')
        call ncinfo(ncname( 4,:),'zc_av','Average Cloud-top height','m','time')
        call ncinfo(ncname( 5,:),'zc_max','Maximum Cloud-top height','m','time')
        call ncinfo(ncname( 6,:),'zi','Boundary layer height','m','time')
        call ncinfo(ncname( 7,:),'we','Entrainment velocity','m/s','time')
        call ncinfo(ncname( 8,:),'lwp_bar','Liquid-water path','kg/m^2','time')
        call ncinfo(ncname( 9,:),'lwp_max','Maximum Liquid-water path','kg/m^2','time')
        call ncinfo(ncname(10,:),'wmax','Maximum vertical velocity','m/s','time')
        call ncinfo(ncname(11,:),'vtke','Vertical integral of total TKE','kg/s','time')
        call ncinfo(ncname(12,:),'lmax','Maximum liquid water specific humidity','kg/kg','time')
        call ncinfo(ncname(13,:),'ustar','Surface friction velocity','m/s','time')
        call ncinfo(ncname(14,:),'tstr','Turbulent temperature scale','K','time')
        call ncinfo(ncname(15,:),'qtstr','Turbulent humidity scale','K','time')
        call ncinfo(ncname(16,:),'obukh','Obukhov Length','m','time')
        call ncinfo(ncname(17,:),'thlskin','Surface liquid water potential temperature','K','time')
        call ncinfo(ncname(18,:),'z0','Roughness height','m','time')
        call ncinfo(ncname(19,:),'wtheta','Surface kinematic temperature flux','K m/s','time')
        call ncinfo(ncname(20,:),'wthetav','Surface kinematic virtual temperature flux','K m/s','time')
        call ncinfo(ncname(21,:),'wq','Surface kinematic moisture flux','kg/kg m/s','time')

        if(isurf==1) then
          call ncinfo(ncname(22,:),'Qnet','Net radiation','W/m^2','time')
          call ncinfo(ncname(23,:),'H','Sensible heat flux','W/m^2','time')
          call ncinfo(ncname(24,:),'LE','Latent heat flux','W/m^2','time')
          call ncinfo(ncname(25,:),'G0','Ground heat flux','W/m^2','time')
          call ncinfo(ncname(26,:),'tendskin','Skin tendency','W/m^2','time')
          call ncinfo(ncname(27,:),'rs','Surface resistance','s/m','time')
          call ncinfo(ncname(28,:),'ra','Aerodynamic resistance','s/m','time')
          call ncinfo(ncname(29,:),'cliq','Fraction of vegetated surface covered with liquid water','-','time')
          call ncinfo(ncname(30,:),'Wl','Liquid water reservoir','m','time')
          call ncinfo(ncname(31,:),'rssoil','Soil evaporation resistance','s/m','time')
          call ncinfo(ncname(32,:),'rsveg','Vegitation resistance','s/m','time')
        end if
        if (imicro.eq.imicro_bulk3)then        ! #sb3 START
          call ncinfo(ncname( nvar0+1,:),'cl_frac','Liquid Cloud fraction','-','time')
          call ncinfo(ncname( nvar0+2,:),'ci_frac','Ice Cloud fraction','-','time')
          call ncinfo(ncname( nvar0+3,:),'ctot_frac','Total Cloud fraction','-','time')
          call ncinfo(ncname( nvar0+4,:),'zb_l_av','Average Liquid Cloud-base height','m','time')
          call ncinfo(ncname( nvar0+5,:),'zb_l_min','Minimal Liquid Cloud-base height','m','time')
          call ncinfo(ncname( nvar0+6,:),'zc_l_av','Average Liquid Cloud-top height','m','time')
          call ncinfo(ncname( nvar0+7,:),'zc_l_max','Maximum Liquid Cloud-top height','m','time')
          call ncinfo(ncname( nvar0+8,:),'zb_i_av','Average Ice Cloud-base height','m','time')
          call ncinfo(ncname( nvar0+9,:),'zb_i_min','Minimal Ice Cloud-base height','m','time')
          call ncinfo(ncname( nvar0+10,:),'zc_i_av','Average Ice Cloud-top height','m','time')
          call ncinfo(ncname( nvar0+11,:),'zc_i_max','Maximum Ice Cloud-top height','m','time')          
          call ncinfo(ncname( nvar0+12,:),'clwp_bar','Cloud Liquid-water path','kg/m^2','time')
          call ncinfo(ncname( nvar0+13,:),'clwp_max','Maximum Cloud Liquid-water path','kg/m^2','time')
          call ncinfo(ncname( nvar0+14,:),'rlwp_bar','Rain Liquid-water path','kg/m^2','time')
          call ncinfo(ncname( nvar0+15,:),'rlwp_max','Maximum Rain Liquid-water path','kg/m^2','time') 
          call ncinfo(ncname( nvar0+16,:),'icwp_bar','Cloud Ice-water path','kg/m^2','time')
          call ncinfo(ncname( nvar0+17,:),'icwp_max','Maximum Cloud Ice-water path','kg/m^2','time')
          call ncinfo(ncname( nvar0+18,:),'siwp_bar','Snow Ice-water path','kg/m^2','time')
          call ncinfo(ncname( nvar0+19,:),'siwp_max','Maximum Snow Ice-water path','kg/m^2','time')            
          call ncinfo(ncname( nvar0+20,:),'giwp_bar','Graupel Ice-water path','kg/m^2','time')
          call ncinfo(ncname( nvar0+21,:),'giwp_max','Graupel Snow Ice-water path','kg/m^2','time')  
          call ncinfo(ncname( nvar0+22,:),'sfc_precw_av','Average surface precipitation','kg/m^2 /s','time')
          call ncinfo(ncname( nvar0+23,:),'sfc_prec_av','Average surface precipitation','W/m^2','time')
           ! - to add more
        endif      ! #sb3 END     
        call open_nc(fname,  ncid,nrec)
        if(nrec==0) call define_nc( ncid, NVar, ncname)
      end if
    end if

    if (lhetero) then
      allocate(zbase_field  (2:i1,2:j1))
      allocate(ztop_field   (2:i1,2:j1))
      allocate(cc_field     (2:i1,2:j1))
      allocate(qlint_field  (2:i1,2:j1))
      allocate(tke_tot_field(2:i1,2:j1))

      allocate(zbase_patch    (xpatches,ypatches))
      allocate(ztop_patch     (xpatches,ypatches))
      allocate(zbasemin_patch (xpatches,ypatches))
      allocate(zbasemin_patchl(xpatches,ypatches))
      allocate(cc_patch       (xpatches,ypatches))
      allocate(qlint_patch    (xpatches,ypatches))
      allocate(qlintmax_patch (xpatches,ypatches))
      allocate(qlintmax_patchl(xpatches,ypatches))
      allocate(tke_tot_patch  (xpatches,ypatches))
      allocate(wmax_patch     (xpatches,ypatches))
      allocate(wmax_patchl    (xpatches,ypatches))
      allocate(qlmax_patch    (xpatches,ypatches))
      allocate(qlmax_patchl   (xpatches,ypatches))
      allocate(ztopmax_patch  (xpatches,ypatches))
      allocate(ztopmax_patchl (xpatches,ypatches))

      allocate(u0av_patch(xpatches,ypatches))
      allocate(v0av_patch(xpatches,ypatches))
      allocate(w0av_patch(xpatches,ypatches))

      allocate(ust_patch (xpatches,ypatches))
      allocate(qst_patch (xpatches,ypatches))
      allocate(tst_patch (xpatches,ypatches))
      allocate(wthls_patch (xpatches,ypatches))
      allocate(wqls_patch(xpatches,ypatches))
      allocate(wthvs_patch(xpatches,ypatches))

      allocate(Qnet_patch    (xpatches,ypatches))
      allocate(H_patch       (xpatches,ypatches))
      allocate(LE_patch      (xpatches,ypatches))
      allocate(G0_patch      (xpatches,ypatches))
      allocate(tendskin_patch(xpatches,ypatches))
      allocate(rs_patch      (xpatches,ypatches))
      allocate(ra_patch      (xpatches,ypatches))
      allocate(cliq_patch    (xpatches,ypatches))
      allocate(wl_patch      (xpatches,ypatches))
      allocate(rsveg_patch   (xpatches,ypatches))
      allocate(rssoil_patch  (xpatches,ypatches))
      allocate(tskin_patch   (xpatches,ypatches))
      allocate(obl_patch     (xpatches,ypatches))

      allocate(zi_patch      (xpatches,ypatches))
      allocate(zi_field      (2:i1,2:j1))
      allocate(ziold_patch   (xpatches,ypatches))
      ziold_patch = -1
      allocate(we_patch      (xpatches,ypatches))
    endif

  end subroutine inittimestat

!>Run timestat. Calculate and write the statistics
  subroutine timestat

    use modglobal,  only : i1,j1,kmax,zf,dzf,cu,cv,rv,rd,eps1,&
                          ijtot,timee,rtimee,dt_lim,rk3step,cexpnr,ifoutput &
                          ,rlv   ! #sb3
!
    use modfields,  only : um,vm,wm,e12m,ql0,u0av,v0av,rhof,u0,v0,w0, sv0 !#sb3
    use modsurfdata,only : wtsurf, wqsurf, isurf,ustar,thlflux,qtflux,z0,oblav,qts,thls,&
                           Qnet, H, LE, G0, rs, ra, tskin, tendskin, &
                           cliq,rsveg,rssoil,Wl, &
                           lhetero, xpatches, ypatches, qts_patch, wt_patch, wq_patch, thls_patch,obl,z0mav_patch, wco2av, Anav, Respav,gcco2av
    use modsurface, only : patchxnr,patchynr
    use mpi
    use modmpi,     only : my_real,mpi_sum,mpi_max,mpi_min,comm3d,mpierr,myid
    use modstat_nc,  only : lnetcdf, writestat_nc,nc_fillvalue
    use modmicrodata,  only : imicro,imicro_bulk3                          ! #sb3
    use modmicrodata3, only : iq_cl,iq_ci,iq_hr, iq_hs,iq_hg,rlvi       &  ! #sb3
                             ,precep_l, precep_i                        &  ! #sb3
                             ,q_cl_statmin,q_ci_statmin                    ! #sb3

    implicit none

    real   :: zbaseavl, ztopavl, ztopmaxl, ztop,zbaseminl
    real   :: qlintavl, qlintmaxl, tke_totl
    real   :: ccl, wmaxl, qlmaxl
    real   :: ust,tst,qst,ustl,tstl,qstl,thlfluxl,qtfluxl
    real   :: usttst, ustqst
    real   :: wts, wqls,wthvs
    real   :: c1,c2 !Used to calculate wthvs
    real,dimension(nvar) :: vars

    ! lsm variables
    real   :: Qnetavl, Havl, LEavl, G0avl, tendskinavl, rsavl, raavl, tskinavl,Wlavl,cliqavl,rsvegavl,rssoilavl
    real   :: Qnetav, Hav, LEav, G0av, tendskinav, rsav, raav, tskinav,Wlav,cliqav,rsvegav,rssoilav
    integer:: i, j, k

    ! heterogeneity variables
    integer:: patchx, patchy
    
    !  stats for SB mixed-phase microphysics #sb3 START 
    real :: clcl,cicl,ctcl                                     &
           ,zb_cl_avl,zb_cl_minl,zt_cl_avl,zt_cl_maxl          &
           ,zb_ci_avl,zb_ci_minl,zt_ci_avl,zt_ci_maxl          & 
           ,qcl_intavl,qci_intavl,qhr_intavl,qhs_intavl        &
           ,qhg_intavl                                         &
           ,qcl_intmaxl,qci_intmaxl,qhr_intmaxl,qhs_intmaxl    &
           ,qhg_intmaxl                                        &
           ,qcl_int,qci_int,qhr_int,qhs_int,qhg_int            &
           ,ztopcl,ztopci                                      &
           ,sfc_precw_avl, sfc_prec_avl
    !- add more if needed
    ! #sb3 end


    if (.not.(ltimestat)) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    if (lhetero) then
      zbase_field    = 0
      ztop_field     = 0
      cc_field       = 0
      qlint_field    = 0
      tke_tot_field  = 0

      zbase_patch    = 0
      ztop_patch     = 0
      zbasemin_patch = zf(kmax)
      zbasemin_patchl= zf(kmax)
      cc_patch       = 0
      qlint_patch    = 0
      qlintmax_patch = 0
      qlintmax_patchl= 0
      tke_tot_patch  = 0
      wmax_patch     = 0
      wmax_patchl    = 0
      qlmax_patch    = 0
      qlmax_patchl   = 0
      ztopmax_patch  = 0
      ztopmax_patchl = 0

      ust_patch      = 0
      qst_patch      = 0
      tst_patch      = 0
      wthls_patch      = 0
      wqls_patch     = 0
      wthvs_patch     = 0

      Qnet_patch     = 0
      H_patch        = 0
      LE_patch       = 0
      G0_patch       = 0
      tendskin_patch = 0
      rs_patch       = 0
      ra_patch       = 0
      cliq_patch     = 0
      wl_patch       = 0
      rsveg_patch    = 0
      rssoil_patch   = 0
      tskin_patch    = 0
      obl_patch      = 0

      zi_patch       = 0
      zi_field       = 0
      we_patch       = 0
    endif

    !      -----------------------------------------------------------
  !     1     EVALUATION OF CLOUD COVER, CLOUD BASE, ETC.
  !    -----------------------------------------------------------

  !     -----------------------------------------------------
  !     1.   Set A:  entrainment and time evolution
  !     -----------------------------------------------------

    zbaseavl = 0.0
    ztopavl = 0.0
    zbaseminl = zf(kmax)

    store_zi = .true.

    call calcblheight

    store_zi = .false.

  !     --------------------------------------------------------------
  !     9.2  liq. waterpath, cloudcover, cloudbase and cloudtop
  !     --------------------------------------------------------------

    ccl      = 0.0
    qlintavl = 0.0
    qlintmaxl= 0.0
    tke_totl = 0.0

    do j=2,j1
      if (lhetero) then
        patchy = patchynr(j)
      endif
      do i=2,i1
        if (lhetero) then
          patchx = patchxnr(i)
        endif

        qlint     = 0.0
        do k=1,kmax
          qlint = qlint + ql0(i,j,k)*rhof(k)*dzf(k)
        end do
        if (qlint>0.) then
          ccl      = ccl      + 1.0
          qlintavl = qlintavl + qlint
          qlintmaxl = max(qlint,qlintmaxl)
          if (lhetero) then
            cc_field(i,j)                  = 1.0
            qlint_field(i,j)               = qlint
            qlintmax_patchl(patchx,patchy) = max(qlintmax_patchl(patchx,patchy),qlint)
          endif
        end if

        do k=1,kmax
          if (ql0(i,j,k) > 0.) then
            zbaseavl = zbaseavl + zf(k)
            zbaseminl = min(zf(k),zbaseminl)
            if (lhetero) then
              zbase_field(i,j)               = zf(k)
              zbasemin_patchl(patchx,patchy) = min(zbasemin_patchl(patchx,patchy),zf(k))
            endif
            exit
          end if
        end do

      end do
    end do

    call MPI_ALLREDUCE(ccl   , cc   , 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qlintavl, qlintav, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qlintmaxl, qlintmax, 1,    MY_REAL, &
                          MPI_MAX, comm3d,mpierr)
    call MPI_ALLREDUCE(zbaseavl, zbaseav, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(zbaseminl, zbasemin, 1,    MY_REAL, &
                          MPI_MIN, comm3d,mpierr)

    if (lhetero) then
      cc_patch    = patchsum_1level(cc_field   )
      qlint_patch = patchsum_1level(qlint_field)
      zbase_patch = patchsum_1level(zbase_field)
      call MPI_ALLREDUCE(qlintmax_patchl, qlintmax_patch, xpatches*ypatches,  MY_REAL,  MPI_MAX, comm3d,mpierr)
      call MPI_ALLREDUCE(zbasemin_patchl, zbasemin_patch, xpatches*ypatches,  MY_REAL,  MPI_MIN, comm3d,mpierr)
    endif
    
   !     --------------------------------------------------------------
   !     9.2 B  ice. waterpath, cloudcover, cloudbase and cloudtop
   !     --------------------------------------------------------------
   ! #sb3 START
     if (imicro.eq.imicro_bulk3) then 
     
     zb_cl_avl = 0.0
     zb_cl_minl = zf(kmax)
     zt_cl_avl = 0.0
     zt_cl_maxl = 0.0
     
     zb_ci_avl = 0.0
     zb_ci_minl = zf(kmax)
     zt_ci_avl = 0.0
     zt_ci_maxl = 0.0   
     
     clcl      = 0.0 ! ccl      = 0.0
     cicl      = 0.0 
     ctcl      = 0.0
     qcl_intavl = 0.0 ! qlintavl = 0.0 
     qci_intavl = 0.0 
     qhr_intavl = 0.0 
     qhs_intavl = 0.0
     qhg_intavl = 0.0
     qcl_intmaxl= 0.0 ! qlintmaxl= 0.0
     qci_intmaxl = 0.0 
     qhr_intmaxl = 0.0 
     qhs_intmaxl = 0.0
     qhg_intmaxl = 0.0   
     !
     sfc_prec_avl = 0.0
     sfc_precw_avl = 0.0
     !tke_totl = 0.0 ! tke_totl = 0.0 
     
     do j=2,j1
       !if (lhetero) then
       !  patchy = patchynr(j)
       !endif
       do i=2,i1
        ! if (lhetero) then
        !   patchx = patchxnr(i)
        ! endif
 
         qcl_int     = 0.0 ! qlint     = 0.0
         qci_int     = 0.0
         qhr_int     = 0.0
         qhs_int     = 0.0
         qhg_int     = 0.0
         ztopcl      = 0.0
         ztopci      = 0.0
         
         do k=1,kmax
          qcl_int=qcl_int+sv0(i,j,k,iq_cl)*rhof(k)*dzf(k) ! qlint = qlint + ql0(i,j,k)*rhof(k)*dzf(k)
          qci_int=qci_int+sv0(i,j,k,iq_ci)*rhof(k)*dzf(k)
          qhr_int=qhr_int+sv0(i,j,k,iq_hr)*rhof(k)*dzf(k)
          qhs_int=qhs_int+sv0(i,j,k,iq_hs)*rhof(k)*dzf(k)
          qhg_int=qhg_int+sv0(i,j,k,iq_hg)*rhof(k)*dzf(k)
         end do
         if (qcl_int.gt.0.) then
           qcl_intavl = qcl_intavl + qcl_int
           qcl_intmaxl = max(qcl_int,qcl_intmaxl)
           !if (lhetero) then
           !  cc_field(i,j)                  = 1.0
           !  qlint_field(i,j)               = qlint
           !  qlintmax_patchl(patchx,patchy) = max(qlintmax_patchl(patchx,patchy),qlint)
           !endif
         end if
         if (qcl_int.gt.q_cl_statmin) then ! used to be 0.0
           clcl      = clcl      + 1.0
         endif
         if (qci_int.gt.q_ci_statmin) then ! used to be 0.0
           cicl      = cicl      + 1.0    
         endif
         if (qci_int.gt.0.) then
           qci_intavl = qci_intavl + qci_int
           qci_intmaxl = max(qci_int,qci_intmaxl)
           !if (lhetero) then
           !  cc_field(i,j)                  = 1.0
           !  qlint_field(i,j)               = qlint
           !  qlintmax_patchl(patchx,patchy) = max(qlintmax_patchl(patchx,patchy),qlint)
           !endif
         end if
         if ((qci_int.gt.q_ci_statmin).or.(qcl_int.gt.q_cl_statmin)) then
           ctcl = ctcl +1.0
         endif
         if (qhr_int>0.) then
           ! ci_c      = ci_c      + 1.0
           qhr_intavl = qhr_intavl + qhr_int
           qhr_intmaxl = max(qhr_int,qhr_intmaxl)
           !if (lhetero) then
           !  cc_field(i,j)                  = 1.0
           !  qlint_field(i,j)               = qlint
           !  qlintmax_patchl(patchx,patchy) = max(qlintmax_patchl(patchx,patchy),qlint)
           !endif
         end if        
         if (qhs_int>0.) then
           ! ci_c      = ci_c      + 1.0
           qhs_intavl = qhs_intavl + qhs_int
           qhs_intmaxl = max(qhs_int,qhs_intmaxl)
           !if (lhetero) then
           !  cc_field(i,j)                  = 1.0
           !  qlint_field(i,j)               = qlint
           !  qlintmax_patchl(patchx,patchy) = max(qlintmax_patchl(patchx,patchy),qlint)
           !endif
         end if        
         if (qhg_int>0.) then
           ! ci_c      = ci_c      + 1.0
           qhg_intavl = qhg_intavl + qhg_int
           qhg_intmaxl = max(qhg_int,qhg_intmaxl)
           !if (lhetero) then
           !  cc_field(i,j)                  = 1.0
           !  qlint_field(i,j)               = qlint
           !  qlintmax_patchl(patchx,patchy) = max(qlintmax_patchl(patchx,patchy),qlint)
           !endif
         end if
         ! loop for liquid cloud top
         do k=1,kmax 
           if (sv0(i,j,k,iq_cl) > q_cl_statmin) then ! (sv0(i,j,k,iq_cl) > 0.) then
             zb_cl_avl = zb_cl_avl + zf(k) ! zbaseavl = zbaseavl + zf(k)
             zb_cl_minl = min(zf(k),zb_cl_minl)! zbaseminl = min(zf(k),zbaseminl)
             ! if (lhetero) then
             !   zbase_field(i,j)               = zf(k)
             !   zbasemin_patchl(patchx,patchy) = min(zbasemin_patchl(patchx,patchy),zf(k))
             ! endif
             exit
           end if
         end do
         ! loop for ice cloud base
         do k=1,kmax 
           if (sv0(i,j,k,iq_ci) > q_ci_statmin) then ! (sv0(i,j,k,iq_ci) > 0.) then
             zb_ci_avl = zb_ci_avl + zf(k) ! zbaseavl = zbaseavl + zf(k)
             zb_ci_minl = min(zf(k),zb_ci_minl)! zbaseminl = min(zf(k),zbaseminl)
             ! if (lhetero) then
             !   zbase_field(i,j)               = zf(k)
             !   zbasemin_patchl(patchx,patchy) = min(zbasemin_patchl(patchx,patchy),zf(k))
             ! endif
             exit
           end if
         end do
         !
         ! cloud tops
          do  k=1,kmax
            if (sv0(i,j,k,iq_cl) >q_cl_statmin) then  ! (sv0(i,j,k,iq_cl) > 0) then    !if (ql0(i,j,k) > 0) then
              ztopcl = zf(k)               ! ztop = zf(k)
            endif                        !endif
            ! wmaxl = max(wm(i,j,k),wmaxl) !wmaxl = max(wm(i,j,k),wmaxl)
            ! qlmaxl = max(ql0(i,j,k),qlmaxl) !qlmaxl = max(ql0(i,j,k),qlmaxl)
            if (sv0(i,j,k,iq_ci) > q_ci_statmin ) then ! (sv0(i,j,k,iq_ci) > 0) then   !if (ql0(i,j,k) > 0) then
              ztopci = zf(k)               ! ztop = zf(k)
            endif                        !endif           
          end do
          !if (lhetero) then
          !do  k=1,kmax
          !   if (ql0(i,j,k) > 0) ztop_field(i,j) = zf(k)
          !   wmax_patchl(patchx,patchy)  = max(wmax_patchl (patchx,patchy),wm (i,j,k))
          !   qlmax_patchl(patchx,patchy) = max(qlmax_patchl(patchx,patchy),ql0(i,j,k))
          !enddo
          !endif
 
         zt_cl_avl = zt_cl_avl + ztopcl
         zt_ci_avl = zt_ci_avl + ztopci
         if (ztopcl > zt_cl_maxl) then
           zt_cl_maxl = ztopcl
         endif 
         if (ztopci > zt_ci_maxl) then
           zt_ci_maxl = ztopci
         endif 
         
         ! if (lhetero) then
         !  ztop_field = ztop
         !  if (ztop > ztopmax_patchl(patchx,patchy)) ztopmax_patchl(patchx,patchy) = ztop
         ! endif
       end do
     end do
     
     ! and th surface precipitation 
     ! send to other processors
      do j=2,j1
      do i=2,i1
       ! calculate surface precipitation
       sfc_precw_avl = sfc_precw_avl+1000.0*rhof(1)         &
                       *(precep_l(i,j)+precep_i(i,j))
       sfc_prec_avl = sfc_prec_avl                          &
                      +rlv*rhof(1)*precep_l(i,j)            &
                      +rlvi*rhof(1)*precep_i(i,j)
      enddo 
      enddo
      
     !call MPI_ALLREDUCE(ccl   , cc   , 1,    MY_REAL, &
     !                      MPI_SUM, comm3d,mpierr)
     call MPI_ALLREDUCE(clcl   , cl_c   , 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
     call MPI_ALLREDUCE(cicl   , ci_c   , 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr) 
      call MPI_ALLREDUCE(ctcl   , c_totc   , 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr)                     
     !call MPI_ALLREDUCE(qlintavl, qlintav, 1,    MY_REAL, &
     !                      MPI_SUM, comm3d,mpierr)
     call MPI_ALLREDUCE(qcl_intavl, qcl_intav, 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr) 
     call MPI_ALLREDUCE(qci_intavl, qci_intav, 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr) 
     call MPI_ALLREDUCE(qhr_intavl, qhr_intav, 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr) 
     call MPI_ALLREDUCE(qhs_intavl, qhs_intav, 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr) 
     call MPI_ALLREDUCE(qhg_intavl, qhg_intav, 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr)  
                           
     !call MPI_ALLREDUCE(qlintmaxl, qlintmax, 1,    MY_REAL, &
     !                      MPI_MAX, comm3d,mpierr)
     call MPI_ALLREDUCE(qcl_intmaxl, qcl_intmax, 1,    MY_REAL, &
                           MPI_MAX, comm3d,mpierr)
     call MPI_ALLREDUCE(qci_intmaxl, qci_intmax, 1,    MY_REAL, &
                           MPI_MAX, comm3d,mpierr)
     call MPI_ALLREDUCE(qhr_intmaxl, qhr_intmax, 1,    MY_REAL, &
                           MPI_MAX, comm3d,mpierr)
     call MPI_ALLREDUCE(qhs_intmaxl, qhs_intmax, 1,    MY_REAL, &
                           MPI_MAX, comm3d,mpierr)  
     call MPI_ALLREDUCE(qhg_intmaxl, qhg_intmax, 1,    MY_REAL, &
                           MPI_MAX, comm3d,mpierr)                          
                        
     !call MPI_ALLREDUCE(zbaseavl, zbaseav, 1,    MY_REAL, &
     !                      MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(zb_cl_avl, zb_cl_av, 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(zb_ci_avl, zb_ci_av, 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr)                          
                           
     !call MPI_ALLREDUCE(zbaseminl, zbasemin, 1,    MY_REAL, &
     !                      MPI_MIN, comm3d,mpierr)
      call MPI_ALLREDUCE(zb_cl_minl, zb_cl_min, 1,    MY_REAL, &
                           MPI_MIN, comm3d,mpierr)                         
      call MPI_ALLREDUCE(zb_ci_minl, zb_ci_min, 1,    MY_REAL, &
                           MPI_MIN, comm3d,mpierr)                         
 
     !if (lhetero) then
     !  cc_patch    = patchsum_1level(cc_field   )
     !  qlint_patch = patchsum_1level(qlint_field)
     !  zbase_patch = patchsum_1level(zbase_field)
     !  call MPI_ALLREDUCE(qlintmax_patchl, qlintmax_patch, xpatches*ypatches,  MY_REAL,  MPI_MAX, comm3d,mpierr)
     !  call MPI_ALLREDUCE(zbasemin_patchl, zbasemin_patch, xpatches*ypatches,  MY_REAL,  MPI_MIN, comm3d,mpierr)
     !endif
     !endif
     ! call MPI_ALLREDUCE(ztopavl, ztopav, 1,    MY_REAL, &
     !                      MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(zt_cl_avl, zt_cl_av, 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(zt_ci_avl, zt_ci_av, 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
                           
     ! call MPI_ALLREDUCE(ztopmaxl, ztopmax, 1,    MY_REAL, &
     !                      MPI_MAX, comm3d,mpierr)
      call MPI_ALLREDUCE(zt_cl_maxl, zt_cl_max, 1,    MY_REAL, &
                           MPI_MAX, comm3d,mpierr)
      call MPI_ALLREDUCE(zt_ci_maxl, zt_ci_max, 1,    MY_REAL, &
                           MPI_MAX, comm3d,mpierr)    
      
     ! sum precipitation  
      call MPI_ALLREDUCE(sfc_precw_avl,sfc_precw_av, 1, MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(sfc_prec_avl, sfc_prec_av, 1,  MY_REAL, &
                           MPI_SUM, comm3d,mpierr)               
                           
     endif ! imicro #sb3    

  !     ---------------------------------------
  !     9.3  determine maximum ql_max and w_max
  !     ---------------------------------------

    wmaxl  = 0.0
    qlmaxl = 0.0
    ztopavl = 0.0
    ztopmaxl = 0.0

    do  j=2,j1
      if (lhetero) then
        patchy = patchynr(j)
      endif
      do  i=2,i1
        if (lhetero) then
          patchx = patchxnr(i)
        endif
         ztop  = 0.0

         do  k=1,kmax
           if (ql0(i,j,k) > 0) ztop = zf(k)
           wmaxl = max(wm(i,j,k),wmaxl)
           qlmaxl = max(ql0(i,j,k),qlmaxl)
         end do

         if (lhetero) then
         do  k=1,kmax
            if (ql0(i,j,k) > 0) ztop_field(i,j) = zf(k)
            wmax_patchl(patchx,patchy)  = max(wmax_patchl (patchx,patchy),wm (i,j,k))
            qlmax_patchl(patchx,patchy) = max(qlmax_patchl(patchx,patchy),ql0(i,j,k))
         enddo
         endif

        ztopavl = ztopavl + ztop
        if (ztop > ztopmaxl) ztopmaxl = ztop

        if (lhetero) then
          ztop_field = ztop
          if (ztop > ztopmax_patchl(patchx,patchy)) ztopmax_patchl(patchx,patchy) = ztop
        endif
      end do
    end do

    call MPI_ALLREDUCE(wmaxl   , wmax   , 1,    MY_REAL, &
                          MPI_MAX, comm3d,mpierr)
    call MPI_ALLREDUCE(qlmaxl, qlmax, 1,    MY_REAL, &
                          MPI_MAX, comm3d,mpierr)
    call MPI_ALLREDUCE(ztopavl, ztopav, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(ztopmaxl, ztopmax, 1,    MY_REAL, &
                          MPI_MAX, comm3d,mpierr)

    if (lhetero) then
      call MPI_ALLREDUCE(wmax_patchl ,      wmax_patch, xpatches*ypatches,  MY_REAL,  MPI_MAX, comm3d,mpierr)
      call MPI_ALLREDUCE(qlmax_patchl,     qlmax_patch, xpatches*ypatches,  MY_REAL,  MPI_MAX, comm3d,mpierr)
      call MPI_ALLREDUCE(ztopmax_patchl, ztopmax_patch, xpatches*ypatches,  MY_REAL,  MPI_MAX, comm3d,mpierr)
      ztop_patch = patchsum_1level(ztop_field)
    endif
  !     -------------------------
  !     9.4  normalise the fields
  !     -------------------------

    if (cc > 0.0) then
      zbaseav = zbaseav/cc
      ztopav  = ztopav/cc
    else
      zbaseav = 0.0
      ztopav = 0.0
    end if

    cc      = cc/ijtot
    qlintav = qlintav / ijtot !domain averaged liquid water path

    if (lhetero) then
      do j=1,ypatches
         do i=1,xpatches
           if (cc_patch(i,j) > 0.0) then
             zbase_patch(i,j) = zbase_patch(i,j)/cc_patch(i,j)
             ztop_patch (i,j) = ztop_patch(i,j) /cc_patch(i,j)
           else
             zbase_patch(i,j) = 0.0
             ztop_patch (i,j) = 0.0
           endif
           cc_patch    = cc_patch    * (xpatches*ypatches/ijtot)
           qlint_patch = qlint_patch * (xpatches*ypatches/ijtot)
        enddo
      enddo
    endif
    
   !     -------------------------
   !     9.4 B normalise the fields
   !     -------------------------
    if (imicro.eq.imicro_bulk3)then ! #sb3 START
     !if (lhetero) then
     !  do j=1,ypatches
     !     do i=1,xpatches
     !       if (cc_patch(i,j) > 0.0) then
     !         zbase_patch(i,j) = zbase_patch(i,j)/cc_patch(i,j)
     !         ztop_patch (i,j) = ztop_patch(i,j) /cc_patch(i,j)
     !       else
     !         zbase_patch(i,j) = 0.0
     !         ztop_patch (i,j) = 0.0
     !       endif
     !       cc_patch    = cc_patch    * (xpatches*ypatches/ijtot)
     !       qlint_patch = qlint_patch * (xpatches*ypatches/ijtot)
     !    enddo
     !  enddo
     !  endif

        if (cl_c > 0.0) then        ! if (cc > 0.0) then
         zb_cl_av = zb_cl_av/cl_c    ! zbaseav = zbaseav/cc
         zt_cl_av = zt_cl_av/cl_c    ! ztopav  = ztopav/cc
       else
         zb_cl_av = 0.0            ! zbaseav = 0.0
         zt_cl_av = 0.0            ! ztopav = 0.0
       end if
       if (ci_c > 0.0) then        ! if (cc > 0.0) then
         zb_ci_av = zb_ci_av/ci_c    ! zbaseav = zbaseav/cc
         zt_ci_av = zt_ci_av/ci_c    ! ztopav  = ztopav/cc
       else
         zb_ci_av = 0.0            ! zbaseav = 0.0
         zt_ci_av = 0.0            ! ztopav = 0.0
       end if     
       
       cl_c      = cl_c/ijtot      ! cc      = cc/ijtot
       ci_c      = ci_c/ijtot      !
       c_totc    = c_totc/ijtot
       !domain averaged  water path
       qcl_intav = qcl_intav / ijtot   ! qlintav = qlintav / ijtot
       qci_intav = qci_intav / ijtot
       qhr_intav = qhr_intav / ijtot
       qhs_intav = qhs_intav / ijtot
       qhg_intav = qhg_intav / ijtot
       
       ! domain averaged precipitation       
       sfc_precw_av = sfc_precw_av / ijtot
       sfc_prec_av  = sfc_prec_av / ijtot
 
       !if (lhetero) then
       !   do j=1,ypatches
       !   do i=1,xpatches
       !     if (cc_patch(i,j) > 0.0) then
       !       zbase_patch(i,j) = zbase_patch(i,j)/cc_patch(i,j)
       !       ztop_patch (i,j) = ztop_patch(i,j) /cc_patch(i,j)
       !     else
       !       zbase_patch(i,j) = 0.0
       !       ztop_patch (i,j) = 0.0
       !     endif
       !     cc_patch    = cc_patch    * (xpatches*ypatches/ijtot)
       !     qlint_patch = qlint_patch * (xpatches*ypatches/ijtot)
       !  enddo
       !  enddo
       !endif      
    endif ! #sb3 END


  !     -------------------------
  !     9.5  Domain Averaged TKE
  !     -------------------------

    do  k=1,kmax
      if (lhetero) then
        u0av_patch = patchsum_1level(u0(2:i1,2:j1,k)) * (xpatches*ypatches/ijtot)
        v0av_patch = patchsum_1level(v0(2:i1,2:j1,k)) * (xpatches*ypatches/ijtot)
        w0av_patch = patchsum_1level(w0(2:i1,2:j1,k)) * (xpatches*ypatches/ijtot)
      endif
      do  j=2,j1
        if (lhetero) then
          patchy = patchynr(j)
        endif
      do  i=2,i1
        if (lhetero) then
          patchx = patchxnr(i)
        endif

        tke_totl = tke_totl + 0.5*( &
                            (0.5*(um(i,j,k)+um(i+1,j,k))+cu-u0av(k))**2 &
                            +(0.5*(vm(i,j,k)+vm(i,j+1,k))+cv-v0av(k))**2 &
                            +(0.5*(wm(i,j,k)+wm(i,j,k+1))           )**2 &
                                  ) + e12m(i,j,k)**2

        if (lhetero) then
          tke_tot_field(i,j) = tke_tot_field(i,j) + 0.5*( &
                                   (0.5*(um(i,j,k)+um(i+1,j,k))+cu-u0av_patch(patchx,patchy))**2 + &
                                   (0.5*(vm(i,j,k)+vm(i,j+1,k))+cv-v0av_patch(patchx,patchy))**2 + &
                                   (0.5*(wm(i,j,k)+wm(i,j,k+1))   -w0av_patch(patchx,patchy))**2 &
                                       ) + e12m(i,j,k)**2
        endif
      end do
      end do
    end do

    call MPI_ALLREDUCE(tke_totl, tke_tot, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)

    tke_tot = tke_tot/ijtot

    if (lhetero) then
      tke_tot_patch = patchsum_1level(tke_tot_field) * (xpatches*ypatches/ijtot)
    endif

!     -------------------------
!     9.6  Horizontally  Averaged ustar, tstar and obl
!     -------------------------
    ustl=sum(ustar(2:i1,2:j1))
    tstl=sum(- thlflux(2:i1,2:j1) / ustar(2:i1,2:j1))
    qstl=sum(- qtflux (2:i1,2:j1) / ustar(2:i1,2:j1))

    call MPI_ALLREDUCE(ustl, ust, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(tstl, tst, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qstl, qst, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)

    ust = ust / ijtot
    tst = tst / ijtot
    qst = qst / ijtot

    if (lhetero) then
      ust_patch = patchsum_1level(ustar(2:i1,2:j1)) * (xpatches*ypatches/ijtot)
      tst_patch = patchsum_1level(- thlflux(2:i1,2:j1) / ustar(2:i1,2:j1)) * (xpatches*ypatches/ijtot)
      qst_patch = patchsum_1level(-  qtflux(2:i1,2:j1) / ustar(2:i1,2:j1)) * (xpatches*ypatches/ijtot)
    endif

    if(isurf < 3) then
      thlfluxl = sum(thlflux(2:i1, 2:j1))
      qtfluxl  = sum(qtflux (2:i1, 2:j1))

      call MPI_ALLREDUCE(thlfluxl, usttst, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(qtfluxl,  ustqst, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)

      usttst = -usttst / ijtot
      ustqst = -ustqst / ijtot
    end if

    !Constants c1 and c2
    c1   = 1.+(rv/rd-1)*qts
    c2   = (rv/rd-1)

    if(isurf >= 3) then
      wts  = wtsurf
      wqls = wqsurf
      wthvs = c1*wts + c2*thls*wqls
    else
      wts  = -usttst
      wqls = -ustqst
      wthvs = c1*wts + c2*thls*wqls
    end if

    if (lhetero) then
      if(isurf < 3) then
        wthls_patch  = patchsum_1level(thlflux(2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        wqls_patch = patchsum_1level( qtflux(2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
      else
        wthls_patch  = wt_patch
        wqls_patch = wq_patch
      endif
      wthvs_patch = (1.+(rv/rd-1)*qts_patch) * wthls_patch + c2 * (thls_patch) * wq_patch
      obl_patch  = patchsum_1level(obl(2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
    endif

  !  9.7  Create statistics for the land surface scheme
    if(isurf == 1) then
      Qnetavl      = sum(Qnet(2:i1,2:j1))
      Havl         = sum(H(2:i1,2:j1))
      LEavl        = sum(LE(2:i1,2:j1))
      G0avl        = sum(G0(2:i1,2:j1))
      tendskinavl  = sum(tendskin(2:i1,2:j1))
      rsavl        = sum(rs(2:i1,2:j1))
      raavl        = sum(ra(2:i1,2:j1))
      cliqavl      = sum(cliq(2:i1,2:j1))
      Wlavl        = sum(wl(2:i1,2:j1))
      rsvegavl     = sum(rsveg(2:i1,2:j1))
      rssoilavl    = sum(rssoil(2:i1,2:j1))
      tskinavl     = sum(tskin(2:i1,2:j1))

      call MPI_ALLREDUCE(Qnetavl,     Qnetav,     1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(Havl,        Hav,        1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(LEavl,       LEav,       1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(G0avl,       G0av,       1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(tendskinavl, tendskinav, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(rsavl,       rsav,       1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(raavl,       raav,       1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(cliqavl,     cliqav,     1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(wlavl,       wlav,       1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(rsvegavl,    rsvegav,    1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(rssoilavl,   rssoilav,   1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(tskinavl,    tskinav,    1,  MY_REAL,MPI_SUM, comm3d,mpierr)

      Qnetav        = Qnetav      / ijtot
      Hav           = Hav         / ijtot
      LEav          = LEav        / ijtot
      G0av          = G0av        / ijtot
      tendskinav    = tendskinav  / ijtot
      rsav          = rsav        / ijtot
      raav          = raav        / ijtot
      cliqav        = cliqav      / ijtot
      wlav          = wlav        / ijtot
      rsvegav       = rsvegav     / ijtot
      rssoilav      = rssoilav    / ijtot
      tskinav       = tskinav     / ijtot

      if (lhetero) then
        Qnet_patch     = patchsum_1level(Qnet    (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        H_patch        = patchsum_1level(H       (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        LE_patch       = patchsum_1level(LE      (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        G0_patch       = patchsum_1level(G0      (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        tendskin_patch = patchsum_1level(tendskin(2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        rs_patch       = patchsum_1level(rs      (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        ra_patch       = patchsum_1level(ra      (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        cliq_patch     = patchsum_1level(cliq    (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        wl_patch       = patchsum_1level(wl      (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        rsveg_patch    = patchsum_1level(rsveg   (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        rssoil_patch   = patchsum_1level(rssoil  (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        tskin_patch    = patchsum_1level(tskin   (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
      endif
    end if

  !  9.8  write the results to output file
  !     ---------------------------------------

    if(myid==0)then
       !tmser1
      open (ifoutput,file='tmser1.'//cexpnr,position='append')
      write( ifoutput,'(f10.2,f6.3,4f12.3,f10.4,5f9.3)') &
          rtimee, &
          cc, &
          zbaseav, &
          ztopav, &
          ztopmax, &
          zi, &
          we, &
          qlintav*1000., &
          qlintmax*1000., &
          wmax, &
          tke_tot*dzf(1), &
          qlmax*1000.
      close(ifoutput)

      !tmsurf
      open (ifoutput,file='tmsurf.'//cexpnr,position='append')
      write( ifoutput,'(f10.2,4e11.3,f11.3,4e11.3)') &
          rtimee   ,&
          ust     ,&
          tst     ,&
          qst     ,&
          oblav   ,&
          thls    ,&
          z0      ,&
          wts     ,&
          wthvs    ,&
          wqls
      close(ifoutput)

      if (isurf == 1) then
        !tmlsm
        open (ifoutput,file='tmlsm.'//cexpnr,position='append')
        write(ifoutput,'(f10.2,9f11.3,e13.3, 5f11.3,e13.3)') &
            rtimee       ,&
            Qnetav      ,&
            Hav         ,&
            LEav        ,&
            G0av        ,&
            tendskinav  ,&
            rsav        ,&
            raav        ,&
            tskinav     ,&
            cliqav      ,&
            wlav        ,&
            rssoilav    ,&
            rsvegav     ,&
            Respav      ,&
            wco2av      ,&
            Anav        ,&
            gcco2av
        close(ifoutput)
      end if
      if (lnetcdf) then
        vars( 1) = rtimee
        vars( 2) = cc
        if (vars(2)<eps1) vars(2) = nc_fillvalue
        vars( 3) = zbaseav
        if (vars(3)<eps1) vars(3) = nc_fillvalue
        vars( 4) = ztopav
        if (vars(4)<eps1) vars(4) = nc_fillvalue
        vars( 5) = ztopmax
        if (vars(5)<eps1) vars(5) = nc_fillvalue
        vars( 6) = zi
        vars( 7) = we
        vars( 8) = qlintav
        vars( 9) = qlintmax
        vars(10) = wmax
        vars(11) = tke_tot*dzf(1)
        vars(12) = qlmax
        vars(13) = ust
        vars(14) = tst
        vars(15) = qst
        vars(16) = oblav
        vars(17) = thls
        vars(18) = z0
        vars(19) = wts
        vars(20) = wthvs
        vars(21) = wqls
        if (isurf == 1) then
          vars(22) = Qnetav
          vars(23) = Hav
          vars(24) = LEav
          vars(25) = G0av
          vars(26) = tendskinav
          vars(27) = rsav
          vars(28) = raav
          vars(29) = cliqav
          vars(30) = wlav
          vars(31) = rssoilav
          vars(32) = rsvegav
        end if
        if (imicro.eq.imicro_bulk3)then ! #sb3 START
           vars(nvar0+1)  = cl_c          ! ('cl_frac','Liquid Cloud fraction','-','time')
           vars(nvar0+2)  = ci_c          ! ('ci_frac','Ice Cloud fraction','-','time')
           vars(nvar0+3)  = c_totc        ! ('ctot_frac','Total Cloud fraction','-','time')
           vars(nvar0+4)  = zb_cl_av      ! ('zb_l_av','Average Liquid Cloud-base height','m','time')
           vars(nvar0+5)  = zb_cl_min     ! ('zb_l_min','Minimal Liquid Cloud-base height','m','time')
           vars(nvar0+6)  = zt_cl_av      ! ('zc_l_av','Average Liquid Cloud-top height','m','time')
           vars(nvar0+7)  = zt_cl_max     ! ('zc_l_max','Maximum Liquid Cloud-top height','m','time')
           vars(nvar0+8)  = zb_ci_av      ! ('zb_i_av','Average Ice Cloud-base height','m','time')
           vars(nvar0+9)  = zb_ci_min     ! ('zb_i_min','Minimal Ice Cloud-base height','m','time')
           vars(nvar0+10) = zt_ci_av      ! ('zc_i_av','Average Ice Cloud-top height','m','time')
           vars(nvar0+11) = zt_ci_max     ! ('zc_i_max','Maximum Ice Cloud-top height','m','time')          
           vars(nvar0+12) = qcl_intav     ! ('clwp_bar','Cloud Liquid-water path','kg/m^2','time')
           vars(nvar0+13) = qcl_intmax    ! ('clwp_max','Maximum Cloud Liquid-water path','kg/m^2','time')
           vars(nvar0+14) = qhr_intav     ! ('rlwp_bar','Rain Liquid-water path','kg/m^2','time')
           vars(nvar0+15) = qhr_intmax    ! ('rlwp_max','Maximum Rain Liquid-water path','kg/m^2','time') 
           vars(nvar0+16) = qci_intav     ! ('icwp_bar','Cloud Ice-water path','kg/m^2','time')
           vars(nvar0+17) = qci_intmax    ! ('icwp_max','Maximum Cloud Ice-water path','kg/m^2','time')
           vars(nvar0+18) = qhs_intav     ! ('siwp_bar','Snow Ice-water path','kg/m^2','time')
           vars(nvar0+19) = qhs_intmax    ! ('siwp_max','Maximum Snow Ice-water path','kg/m^2','time')            
           vars(nvar0+20) = qhg_intav     ! ('giwp_bar','Graupel Ice-water path','kg/m^2','time')
           vars(nvar0+21) = qhg_intmax    ! ('giwp_max','Graupel Snow Ice-water path','kg/m^2','time')
           vars(nvar0+22) = sfc_precw_av
           vars(nvar0+23) = sfc_prec_av
           ! vars(nvar0+22) = 
           ! - to add more    
        endif ! #sb3 END

        call writestat_nc(ncid,nvar,ncname,vars,nrec,.true.)
      end if

      if(lhetero) then
        do i=1,xpatches
          do j=1,ypatches
            name = 'tmser1patchiiixjjj.'//cexpnr
            write (name(12:14),'(i3.3)') i
            write (name(16:18),'(i3.3)') j
            open (ifoutput,file=name,position='append')
            write( ifoutput,'(f10.2,f6.3,4f12.3,f10.4,5f9.3)') &
              rtimee, &
              cc_patch(i,j), &
              zbase_patch(i,j), &
              ztop_patch(i,j), &
              ztopmax_patch(i,j), &
              zi_patch(i,j), &
              we_patch(i,j), &
              qlint_patch(i,j)*1000., &
              qlintmax_patch(i,j)*1000., &
              wmax_patch(i,j), &
              tke_tot_patch(i,j)*dzf(1), &
              qlmax_patch(i,j)*1000.
            close(ifoutput)

            name = 'tmsurfpatchiiixjjj.'//cexpnr
            write (name(12:14),'(i3.3)') i
            write (name(16:18),'(i3.3)') j
            open (ifoutput,file=name,position='append')
            write( ifoutput,'(f10.2,4e11.3,f11.3,4e11.3)') &
              rtimee      ,&
              ust_patch(i,j)   ,&
              tst_patch(i,j)   ,&
              qst_patch(i,j)   ,&
              obl_patch(i,j)   ,&
              thls_patch(i,j)  ,&
              z0mav_patch(i,j) ,&
              wthls_patch(i,j)   ,&
              wthvs_patch(i,j)  ,&
              wqls_patch(i,j)
            close(ifoutput)

            if(isurf == 1) then
              name = 'tmlsmpatchiiixjjj.'//cexpnr
              write (name(11:13),'(i3.3)') i
              write (name(15:17),'(i3.3)') j
              open (ifoutput,file=name,position='append')
              write(ifoutput,'(f10.2,9f11.3,e13.3, 2f11.3)') &
                rtimee           ,&
                Qnet_patch(i,j)       ,&
                H_patch(i,j)          ,&
                LE_patch(i,j)         ,&
                G0_patch(i,j)         ,&
                tendskin_patch(i,j)   ,&
                rs_patch(i,j)         ,&
                ra_patch(i,j)         ,&
                tskin_patch(i,j)      ,&
                cliq_patch(i,j)       ,&
                wl_patch(i,j)         ,&
                rssoil_patch(i,j)     ,&
                rsveg_patch(i,j)
              close(ifoutput)
            endif
          enddo
        enddo
      endif

    end if

  end subroutine timestat

!>Calculate the boundary layer height
!!
!! There are 3 available ways to calculate the boundary layer height:
!! - By determining the minimum flux in some scalar, e.g. buoyancy
!! - By determining the minimum local gradient of some scalar, averaged over a definable number of columns
!! - By monitoring a threshold value of some scalar, averaged over a definable number of columns
  subroutine calcblheight

    use modglobal,  only : ih,i1,jh,j1,kmax,k1,cp,rlv,imax,rd,zh,dzh,zf,dzf,rv,ijtot,iadv_sv,iadv_kappa
    use modfields,  only : w0,qt0,qt0h,ql0,thl0,thl0h,thv0h,sv0,exnf,whls
    use modsurfdata,only : svs, lhetero, xpatches, ypatches
    use modsurface, only : patchxnr,patchynr
    use mpi
    use modmpi,     only : mpierr, comm3d,mpi_sum,my_real
    implicit none
    real    :: zil, dhdt,locval,oldlocval
    integer :: location,i,j,k,nsamp,stride
    real, allocatable,dimension(:,:,:) :: blh_fld,  sv0h, blh_fld2
    real, allocatable, dimension(:) :: profile, gradient, dgrad
    allocate(blh_fld(2-ih:i1+ih,2-jh:j1+jh,k1),sv0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(profile(k1),gradient(k1),dgrad(k1))
    if (lhetero) then
      allocate(blh_fld2(2:i1,2:j1,k1))
    endif
    zil = 0.0
    gradient = 0.0
    dgrad = 0.0
    select case (iblh_meth)
    case (iblh_flux)
      select case (iblh_var)
      case(iblh_qt)
        blh_fld = w0*qt0h
      case(iblh_thl)
        blh_fld = w0*thl0h
      case(iblh_thv)
        blh_fld = w0*thv0h
      case(1:)
        if (iadv_sv(iblh_var)==iadv_kappa) then
          call halflev_kappa(sv0(:,:,:,iblh_var),sv0h)
        else
          do  k=2,k1
          do  j=2,j1
          do  i=2,i1
            sv0h(i,j,k) = (sv0(i,j,k,iblh_var)*dzf(k-1)+sv0(i,j,k-1,iblh_var)*dzf(k))/(2*dzh(k))
          enddo
          enddo
          enddo
          sv0h(2:i1,2:j1,1) = svs(iblh_var)
          blh_fld = w0*sv0h
        end if
      end select
    case (iblh_grad,iblh_thres)
      select case (iblh_var)
      case(iblh_qt)
        blh_fld = qt0
      case(iblh_thl)
        blh_fld = thl0
      case(iblh_thv)
        do k=1,k1
          blh_fld(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                      *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
        end do
      case(1:)
        blh_fld = sv0(2:i1,2:j1,1:k1,iblh_var)
      end select
    end select

    select case (iblh_meth)

    case (iblh_flux)
      stride = ceiling(real(imax)/real(blh_nsamp))
      do i=2,stride+1
        nsamp =  ceiling(real(i1-i+1)/real(stride))
        do j=2,j1
          zil = zil + nsamp*zh(minloc(sum(blh_fld(i:i1:stride,j,:),1),1))
        end do
      end do

    case (iblh_grad)
      stride = ceiling(real(imax)/real(blh_nsamp))
      do i=2,stride+1
        nsamp =  ceiling(real(i1-i+1)/real(stride))
        do j=2,j1
          profile  = sum(blh_fld(i:i1:stride,j,:),1)
          select case (iblh_var)
          case(iblh_qt) !Water vapour gradients near the inversion layer can be either positive or negative
            gradient(2:k1) = abs(profile(2:k1) - profile(1:kmax))/dzh(2:k1)
          case(iblh_thl,iblh_thv) !temperature jumps near the inversion layer are always positive
            gradient(2:k1) = (profile(2:k1) - profile(1:kmax))/dzh(2:k1)
          case default
            gradient(2:k1) = (profile(2:k1) - profile(1:kmax))/dzh(2:k1)
          end select
          dgrad(2:kmax)    = (gradient(3:k1) - gradient(2:kmax))/dzf(2:kmax)
          location = maxloc(gradient,1)
          zil  = zil + nsamp*(zh(location-1) - dzh(location)*dgrad(location-1)/(dgrad(location)-dgrad(location-1) + 1.e-8))
        enddo
      enddo
    case (iblh_thres)
      stride = ceiling(real(imax)/real(blh_nsamp))
      do i=2,stride+1
        nsamp =  ceiling(real(i1-i+1)/real(stride))
        do j=2,j1
          locval = 0.0
          do k=kmax,1,-1
            oldlocval = locval
            locval = blh_sign*sum(blh_fld(i:i1:stride,j,k))/nsamp
            if (locval < blh_sign*blh_thres) then
              zil = zil + nsamp *(zf(k) +  (blh_sign*blh_thres-locval) &
                        *dzh(k+1)/(oldlocval-locval))
              exit
            endif
          enddo
        enddo
      enddo

    end select

    if (lhetero) then !Not using stride, instead use adjacent grid points in x-direction (to prevent processor communication)

      select case (iblh_meth)

      case (iblh_flux)
        blh_fld2(2,2:j1,:)        = blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)        + blh_fld(3,2:j1,:)
        blh_fld2(3:(i1-1),2:j1,:) = blh_fld(2:(i1-2),2:j1,:) + blh_fld(3:(i1-1),2:j1,:) + blh_fld(4:i1,2:j1,:)
        blh_fld2(i1,2:j1,:)       = blh_fld(i1-1,2:j1,:)     + blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)

        do i=2,i1
          do j=2,j1
            zi_field(i,j) = zh(minloc(blh_fld2(i,j,:),1))
          end do
        end do

      case (iblh_grad)
        blh_fld2(2,2:j1,:)        = blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)        + blh_fld(3,2:j1,:)
        blh_fld2(3:(i1-1),2:j1,:) = blh_fld(2:(i1-2),2:j1,:) + blh_fld(3:(i1-1),2:j1,:) + blh_fld(4:i1,2:j1,:)
        blh_fld2(i1,2:j1,:)       = blh_fld(i1-1,2:j1,:)     + blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)

        do i=2,i1
          do j=2,j1
            profile  = blh_fld2(i,j,:)
            select case (iblh_var)
            case(iblh_qt) !Water vapour gradients near the inversion layer can be either positive or negative
              gradient(2:k1) = abs(profile(2:k1) - profile(1:kmax))/dzh(2:k1)
            case(iblh_thl,iblh_thv) !temperature jumps near the inversion layer are always positive
              gradient(2:k1) = (profile(2:k1) - profile(1:kmax))/dzh(2:k1)
            case default
              gradient(2:k1) = (profile(2:k1) - profile(1:kmax))/dzh(2:k1)
            end select
            dgrad(2:kmax)    = (gradient(3:k1) - gradient(2:kmax))/dzf(2:kmax)
            location = maxloc(gradient,1)
            zi_field(i,j) = (zh(location-1) - dzh(location)*dgrad(location-1)/(dgrad(location)-dgrad(location-1) + 1.e-8))
          enddo
        enddo

      case (iblh_thres)
        blh_fld2(2,2:j1,:)        = blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)        + blh_fld(3,2:j1,:)
        blh_fld2(3:(i1-1),2:j1,:) = blh_fld(2:(i1-2),2:j1,:) + blh_fld(3:(i1-1),2:j1,:) + blh_fld(4:i1,2:j1,:)
        blh_fld2(i1,2:j1,:)       = blh_fld(i1-1,2:j1,:)     + blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)

        do i=2,i1
          do j=2,j1
            locval = 0.0
            do k=kmax,1,-1
              oldlocval = locval
              locval = blh_sign*blh_fld2(i,j,k)/3
              if (locval < blh_sign*blh_thres) then
                zi_field(i,j) = (zf(k) +  (blh_sign*blh_thres-locval) * dzh(k+1)/(oldlocval-locval))
                exit
              endif
            enddo
          enddo
        enddo

      end select

    endif

    call MPI_ALLREDUCE(zil, zi, 1, MY_REAL, MPI_SUM, comm3d,mpierr)
    zi = zi / ijtot

    if (lhetero) then
      zi_patch = patchsum_1level(zi_field) * (xpatches*ypatches/ijtot)
    endif

    if (ziold< 0) ziold = zi
    dhdt = (zi-ziold)/dtav
    if(store_zi) ziold = zi

    k=2
    do while (zh(k)<zi .and. k < kmax)
      k=k+1
    end do
    we = dhdt - whls (k)   !include for large-scale vertical velocity

    if (lhetero) then
      do j=1,ypatches
      do i=1,xpatches
        if (ziold_patch(i,j)<0) ziold_patch(i,j) = zi_patch(i,j)
        k=2
        do while (zh(k)<zi_patch(i,j) .and. k < kmax)
          k=k+1
        end do
        we_patch(i,j) = ((zi_patch(i,j)-ziold_patch(i,j))/dtav) - whls(k)
      enddo
      enddo
      if(store_zi) ziold_patch = zi_patch
    endif

    deallocate(blh_fld,sv0h)
    if (lhetero) then
      deallocate(blh_fld2)
    endif
    deallocate(profile,gradient,dgrad)

  end subroutine calcblheight

!> Clean up when leaving the run
  subroutine exittimestat
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modsurfdata,only :lhetero
    implicit none

    if(ltimestat .and. lnetcdf .and. myid==0) call exitstat_nc(ncid)

    if (lhetero) then
      deallocate(zbase_field  )
      deallocate(ztop_field   )
      deallocate(cc_field     )
      deallocate(qlint_field  )
      deallocate(tke_tot_field)

      deallocate(zbase_patch    )
      deallocate(ztop_patch     )
      deallocate(zbasemin_patch )
      deallocate(zbasemin_patchl)
      deallocate(cc_patch       )
      deallocate(qlint_patch    )
      deallocate(qlintmax_patch )
      deallocate(qlintmax_patchl)
      deallocate(tke_tot_patch  )
      deallocate(wmax_patch     )
      deallocate(wmax_patchl    )
      deallocate(qlmax_patch    )
      deallocate(qlmax_patchl   )
      deallocate(ztopmax_patch  )
      deallocate(ztopmax_patchl )

      deallocate(u0av_patch)
      deallocate(v0av_patch)
      deallocate(w0av_patch)

      deallocate(ust_patch )
      deallocate(qst_patch )
      deallocate(tst_patch )
      deallocate(wthls_patch )
      deallocate(wqls_patch)
      deallocate(wthvs_patch)

      deallocate(Qnet_patch    )
      deallocate(H_patch       )
      deallocate(LE_patch      )
      deallocate(G0_patch      )
      deallocate(tendskin_patch)
      deallocate(rs_patch      )
      deallocate(ra_patch      )
      deallocate(cliq_patch    )
      deallocate(wl_patch      )
      deallocate(rsveg_patch   )
      deallocate(rssoil_patch  )
      deallocate(tskin_patch   )
      deallocate(obl_patch     )

      deallocate(zi_patch      )
      deallocate(zi_field      )
      deallocate(ziold_patch   )
      deallocate(we_patch      )
    endif
  end subroutine exittimestat

 function patchsum_1level(x)
   use modglobal,  only : imax,jmax
   use modsurface, only : patchxnr, patchynr
   use modsurfdata,only : xpatches,ypatches
   use mpi
   use modmpi,     only : mpierr,comm3d,my_real,mpi_sum
   implicit none
   real                :: patchsum_1level(xpatches,ypatches),xl(xpatches,ypatches)
   real, intent(in)    :: x(imax,jmax)
   integer             :: i,j,iind,jind

   patchsum_1level = 0
   xl              = 0

   do j=1,jmax
     jind  = patchynr(j)

     do i=1,imax
       iind  = patchxnr(i)

       xl(iind,jind) = xl(iind,jind) + x(i,j)
     enddo
   enddo

  call MPI_ALLREDUCE(xl,patchsum_1level, xpatches*ypatches,MY_REAL,MPI_SUM, comm3d,mpierr)

  end function

end module modtimestat
