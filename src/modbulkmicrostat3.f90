!> \file modbulkmicrostat3.f90
!!  Calculates profiles coming from the bulkmicrophysics
!!  in case of modbulkmicro3
!>
!! Profiles coming from the bulkmicrophysics. Written to precep.expnr for the
!! rain rates etc., and to qlptend.expnr, nptend.expnr and qtptend.expnr for the
!! tendencies is rain water content, droplet number, and total water content,
!! respectively.
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Olivier Geoffroy, KNMI
!!  \author Johan van de Dussen, TU Delft
!!  \author Jan Chylik, IGMK
!!  \author Jisk Attema, NLeSC
!!
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
!                 2020 NLeSC
!
module modbulkmicrostat3
  use modglobal, only : longint

implicit none
private
PUBLIC  :: initbulkmicrostat3, bulkmicrostat3, exitbulkmicrostat3     &
          ,bulkmicrotend3
save
!NetCDF variables
  integer,parameter :: nvar = 23
  integer :: nvar_mphys = 0 ! adjusted later
  integer :: ncid_mphys = 0
  integer :: nrec_mphys = 0
  character(80) :: fname_mphys = 'mphysprofiles.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  character(80),allocatable, dimension(:,:) :: ncmname
  character(80),dimension(1,4) :: tncmname
  real          :: dtav, timeav
  integer(kind=longint):: idtav, itimeav, tnext, tnextwrite
  integer          :: nsamples
  logical          :: lmicrostat = .false.
  integer, parameter      :: nrfields = 5 , &
                             iauto    = 2 , &
                             iaccr    = 3 , &
                             ievap    = 4 , &
                             ised     = 5
  integer          :: nmphyout,nbaspout,ncolout
  real, allocatable, dimension(:,:)  :: Npav   , &
                                        Npmn   , &
                                        qlpav  , &
                                        qlpmn  , &
                                        qtpav  , &
                                        qtpmn
  real, allocatable, dimension(:)    :: precmn       , &
                                        preccountmn  , &
                                        prec_prcmn   , &
                                        cloudcountavl, &
                                        cloudcountav , &
                                        cloudcountmn , &
                                        raincountavl , &
                                        raincountav  , &
                                        raincountmn  , &
                                        Nrrainavl    , &
                                        Nrrainav     , &
                                        Nrrainmn     , &
                                        qravl        , &
                                        qrav         , &
                                        qrmn         , &
                                        Dvravl       , &
                                        Dvrav        , &
                                        Dvrmn

    real, allocatable, dimension(:,:) :: tend  ! tend(ntends, k)
real, allocatable, dimension(:) :: cl_countav, hr_countav, ci_countav, hs_countav, hg_countav
real, allocatable, dimension(:) :: cl_countavl,hr_countavl,ci_countavl,hs_countavl,hg_countavl
real, allocatable, dimension(:) :: cl_countmn,hr_countmn,ci_countmn,hs_countmn,hg_countmn 
real, allocatable, dimension(:) :: prec_r_mn,prec_i_mn,prec_s_mn,prec_g_mn
real, allocatable, dimension(:) :: cfrac_tot,cfrac_l,cfrac_i

contains

!> Initialization routine, reads namelists and inits variables
subroutine initbulkmicrostat3
    use modmpi,    only  : myid, mpi_logical, my_real, comm3d, mpierr
    use modglobal, only  : ifnamopt, fname_options, cexpnr, ifoutput     &
              ,dtav_glob, timeav_glob, ladaptive, k1, dtmax,btime,tres   &
              ,lwarmstart,checknamelisterror, kmax
    use modstat_nc, only : lnetcdf,open_nc,define_nc,ncinfo,nctiminfo    &
                          ,writestat_dims_nc
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid
    use modmicrodata,  only : imicro, imicro_bulk, imicro_sice, imicro_bulk3
    implicit none
    integer      :: ierr

    ! read nambulkmicrostat instead of initbulkmicrostat in modbulkmicro.f90
    namelist/NAMBULKMICROSTAT/ &
    lmicrostat, dtav, timeav


    ! if ((imicro /=imicro_bulk) .and. (imicro /= imicro_sice)) return #sb3
    ! - since we want to generate statistics even for imicro_bulk3 or other

    dtav   = dtav_glob
    timeav = timeav_glob
    ! reading the namelist NAMBULKMICROSTAT
    if(myid==0)then
      open (ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMBULKMICROSTAT,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMBULKMICROSTAT'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMBULKMICROSTAT'
      endif
      write(6,NAMBULKMICROSTAT)
      close(ifnamopt)
    endif

    ! adjust these number when needed
    nbaspout = 0 ! 23           ! - to describe the number of basic mphys outputs
    ncolout  = nbaspout + 81    ! - to describe number of microphysical outpus up to collisions
    nmphyout = ncolout  + 21    ! - to describe number of microphysical outputs before total ones
    nvar_mphys = nmphyout + 25  ! then no adjustment of nvar_mphys needed

    ! transmitting setting to other processors
    call MPI_BCAST(lmicrostat  ,1,MPI_LOGICAL  ,0,comm3d,mpierr)
    call MPI_BCAST(dtav        ,1,MY_REAL      ,0,comm3d,mpierr)
    call MPI_BCAST(timeav      ,1,MY_REAL      ,0,comm3d,mpierr)
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   + btime
    tnextwrite = itimeav + btime
    nsamples   = itimeav / idtav

    if (.not. lmicrostat) return

    if (abs(timeav/dtav - nsamples) > 1e-4) then
      stop 'timeav must be an integer multiple of dtav (NAMBULKMICROSTAT)'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax - nint(dtav/dtmax)) > 1e-4) then
      stop 'dtav must be an integer multiple of dtmax (NAMBULKMICROSTAT)'
    end if

    ! allocating array for statistics that should be calculated
    allocate( Npav    (k1, nrfields)  &
             ,Npmn    (k1, nrfields)  &
             ,qlpav   (k1, nrfields)  &
             ,qlpmn   (k1, nrfields)  &
             ,qtpav   (k1, nrfields)  &
             ,qtpmn   (k1, nrfields)  )
    allocate( precmn        (k1)      &
             ,preccountmn   (k1)      &
             ,prec_prcmn    (k1)      &
             ,cloudcountavl (k1)      &
             ,cloudcountav  (k1)      &
             ,cloudcountmn  (k1)      &
             ,raincountavl  (k1)      &
             ,raincountav   (k1)      &
             ,raincountmn   (k1)      &
             ,Nrrainavl     (k1)      &
             ,Nrrainav      (k1)      &
             ,Nrrainmn      (k1)      &
             ,qravl         (k1)      &
             ,qrav          (k1)      &
             ,qrmn          (k1)      &
             ,Dvravl        (k1)      &
             ,Dvrav         (k1)      &
             ,Dvrmn         (k1)      )

    allocate( cl_countavl   (k1)      &
             ,ci_countavl   (k1)      &
             ,hr_countavl   (k1)      &
             ,hs_countavl   (k1)      &
             ,hg_countavl   (k1)      )

    allocate( cl_countav    (k1)      &
             ,ci_countav    (k1)      &
             ,hr_countav    (k1)      &
             ,hs_countav    (k1)      &
             ,hg_countav    (k1)      &
             ,cl_countmn    (k1)      &
             ,ci_countmn    (k1)      &
             ,hr_countmn    (k1)      &
             ,hs_countmn    (k1)      &
             ,hg_countmn    (k1)      &
             ,prec_r_mn     (k1)      &
             ,prec_i_mn     (k1)      &
             ,prec_s_mn     (k1)      &
             ,prec_g_mn     (k1)      )

    Npmn        = 0.0
    qlpmn       = 0.0
    qtpmn       = 0.0
    precmn      = 0.0
    preccountmn = 0.0
    prec_prcmn  = 0.0
    cloudcountmn= 0.0
    raincountmn = 0.0
    Nrrainmn    = 0.0
    qrmn        = 0.0
    Dvrmn       = 0.0

    cl_countmn  = 0.0
    ci_countmn  = 0.0
    hr_countmn  = 0.0
    hs_countmn  = 0.0
    hg_countmn  = 0.0
    prec_r_mn   = 0.0
    prec_i_mn   = 0.0
    prec_s_mn   = 0.0
    prec_g_mn   = 0.0

    ! opening output files
    if (myid == 0 .and. .not. lwarmstart) then
      open (ifoutput,file = 'precep.'//cexpnr ,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'nptend.'//cexpnr ,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'qlptend.'//cexpnr,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'qtptend.'//cexpnr,status = 'replace')
      close(ifoutput)
    end if
    if (lnetcdf) then
      idtav      = idtav_prof
      itimeav    = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples   = itimeav/idtav
      if (myid==0) then
        call ncinfo(tncname(1,:),'time','Time','s','time')
        call ncinfo(ncname( 1,:),'cfrac','Cloud fraction','-','tt')
        call ncinfo(ncname( 2,:),'rainrate','Echo Rain Rate','W/m^2','tt')
        call ncinfo(ncname( 3,:),'preccount','Preccount','W/m^2','tt')
        call ncinfo(ncname( 4,:),'nrrain','nrrain','W/m^2','tt')
        call ncinfo(ncname( 5,:),'raincount','raincount','W/m^2','tt')
        call ncinfo(ncname( 6,:),'precmn','precmn','W/m^2','tt')
        call ncinfo(ncname( 7,:),'dvrmn','dvrmn','W/m^2','tt')
        call ncinfo(ncname( 8,:),'qrmn','qrmn','W/m^2','tt')
        call ncinfo(ncname( 9,:),'npauto','Autoconversion rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(10,:),'npaccr','Accretion rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(11,:),'npsed','Sedimentation rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(12,:),'npevap','Evaporation rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(13,:),'nptot','Total rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(14,:),'qrpauto','Autoconversion rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(15,:),'qrpaccr','Accretion rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(16,:),'qrpsed','Sedimentation rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(17,:),'qrpevap','Evaporation rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(18,:),'qrptot','Total rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(19,:),'qtpauto','Autoconversion total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(20,:),'qtpaccr','Accretion total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(21,:),'qtpsed','Sedimentation total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(22,:),'qtpevap','Evaporation total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(23,:),'qtptot','Total water content tendency','kg/kg/s','tt')

        call define_nc(ncid_prof, nvar , ncname)

        ! and these shall go in a new profile file
        fname_mphys(15:17) = cexpnr

        allocate(ncmname(nvar_mphys,4))
        ! and now for bulk ice microphysics
        ! - cloud formation
        call nctiminfo(tncmname(1,:))
        call ncinfo(ncmname(nbaspout+1,:),'dn_c_nuc',  'cloud droplet nucleation tendency',                     '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+2,:),'dq_c_nuc',  'cloud liquid water content nucleation tendency',       'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+3,:),'dn_i_inuc', 'cloud ice particle nucleation tendency',                '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+4,:),'dq_i_inuc', 'cloud ice content nucleation tendency',                'kg/kg/s','tt')
        ! - deposition and sublimation
        call ncinfo(ncmname(nbaspout+5,:),'dq_i_dep',  'cloud ice content deposition tendency',                'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+6,:),'dq_s_dep',  'snow content deposition tendency',                     'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+7,:),'dq_g_dep',  'graupel content deposition tendency',                  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+8,:),'dq_i_sub',  'cloud ice content sublimation tendency',               'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+9,:),'dq_s_sub',  'snow content sublimation tendency',                    'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+10,:),'dq_g_sub', 'graupel content sublimation  tendency',                'kg/kg/s','tt')
        ! - freezing
        ! -- homogeneous freezing
        call ncinfo(ncmname(nbaspout+11,:),'dq_i_hom',  'homogeneous freezing: cloud ice content tendency',    'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+12,:),'dn_i_hom',  'homogeneous freezing: cloud ice number tendency',      '#/kg/s','tt')
        ! -- heterogeneous freezing
        call ncinfo(ncmname(nbaspout+13,:),'dq_i_het',  'heterogeneous freezing: cloud ice content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+14,:),'dn_i_het',  'heterogeneous freezing: cloud ice number tendency',    '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+15,:),'dq_g_het',  'heterogeneous freezing: graupel content tendency',    'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+16,:),'dn_g_het',  'heterogeneous freezing: graupel number tendency',      '#/kg/s','tt')
        ! - warm microphysics
        ! - saturation aqdjustment tendency
        call ncinfo(ncmname(nbaspout+17,:),'dq_c_cond', 'saturation adjustment: cloud water condensation tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+18,:),'dq_c_ev',   'saturation adjustment: cloud water evaporation tendency',    'kg/kg/s','tt')
        ! -- cloud water self-collection
        call ncinfo(ncmname(nbaspout+19,:),'dn_c_sc',   'cloud water number self-collection tendency',          '#/kg/s','tt')
        ! -- cloud water autoconversion
        call ncinfo(ncmname(nbaspout+20,:),'dq_c_au',   'autoconversion: cloud water content tendency',        'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+21,:),'dn_c_au',   'autoconversion: cloud water number tendency',          '#/kg/s','tt')
        ! -- cloud water accretion
        call ncinfo(ncmname(nbaspout+22,:),'dq_c_ac',   'accretion: cloud water content tendency',             'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+23,:),'dn_c_ac',   'accretion: cloud water number tendency',               '#/kg/s','tt')
        ! -- rain selfcollection and breakup
        call ncinfo(ncmname(nbaspout+24,:),'dn_r_sc',   'rain water number self-collection tendency',           '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+25,:),'dn_r_br',   'rain water number break-up tendency',                  '#/kg/s','tt')
        ! -- rain evaporation
        call ncinfo(ncmname(nbaspout+26,:),'dq_r_ev',   'evaporation: rain water content tendency',            'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+27,:),'dn_r_ev',   'evaporation: raindrop number tendency',                '#/kg/s','tt')
        ! - cold aggregation and self-collection
        ! -- ice aggregation
        call ncinfo(ncmname(nbaspout+28,:),'dq_i_agg',  'ice aggregation: cloud ice water content tendency',   'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+29,:),'dn_i_agg',  'ice aggregation: cloud ice number tendency',           '#/kg/s','tt')
        ! -- snow aggregation
        call ncinfo(ncmname(nbaspout+30,:),'dn_s_sc',   'snow number tendency due to snow selfcollection',      '#/kg/s','tt')
        ! - loss of liquid water due to riming
        call ncinfo(ncmname(nbaspout+31,:),'dq_c_rime', 'riming: cloud liquid water content tendency',         'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+32,:),'dn_c_rime', 'riming:  cloud liquid water number tendency',          '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+33,:),'dq_r_rime', 'riming: rain water content tendency',                 'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+34,:),'dn_r_rime', 'riming: raindrop number tendency',                     '#/kg/s','tt')
        ! - growth by riming
        call ncinfo(ncmname(nbaspout+35,:),'dq_i_rime_ic', 'cloud ice riming by cloud liquid water tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+36,:),'dq_s_rime_sc', 'snow riming by cloud liquid water tendency',       'kg/kg/s','tt')
       !call ncinfo(ncmname(nbaspout+37,:),'dq_s_rime_sr', 'snow riming by rain water tendency',               'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+37,:),'dq_g_rime_gc', 'graupel riming by cloud liquid water tendency',    'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+38,:),'dq_g_rime_gr', 'graupel riming by rain water tendency',            'kg/kg/s','tt')
        ! - collisions
        ! -- collision conversions
        call ncinfo(ncmname(nbaspout+39,:),'dq_i_col_rig', 'collision of r+i: graupel content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+40,:),'dn_i_col_rig', 'collision of r+i: graupel number tendency',   '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+41,:),'dq_r_col_rig', 'collision of r+i: rain water content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+42,:),'dn_r_col_rig', 'collision of r+i: raindrop number tendency',   '#/kg/s','tt')
        ! -- collision conversions
        call ncinfo(ncmname(nbaspout+43,:),'dq_s_col_rsg', 'collision of r+s: snow content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+44,:),'dn_s_col_rsg', 'collision of r+s: snow number tendency',   '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+45,:),'dq_r_col_rsg', 'collision of r+s: rain water content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+46,:),'dn_r_col_rsg', 'collision of r+s: raindrop number tendency',   '#/kg/s','tt')
        ! -- collision collection
        call ncinfo(ncmname(nbaspout+47,:),'dq_i_col_si',  'collection by snow: ice content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+48,:),'dn_i_col_si',  'collection by snow: ice number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+49,:),'dq_s_col_gsg', 'collection by graupel: snow content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+50,:),'dn_s_col_gsg', 'collection by graupel: snow number tendency',  '#/kg/s','tt')
        ! -- enhanced melting
        call ncinfo(ncmname(nbaspout+51,:),'dq_i_eme_ic',  'enhanced melting of ice by liquid clouds: ice content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+52,:),'dn_i_eme_ic',  'enhanced melting of ice by liquid clouds: ice number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+53,:),'dq_i_eme_ri',  'enhanced melting of ice by rain: ice content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+54,:),'dn_i_eme_ri',  'enhanced melting of ice by rain: ice number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+55,:),'dq_s_eme_sc',  'enhanced melting of snow by liquid clouds: snow content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+56,:),'dn_s_eme_sc',  'enhanced melting of snow by liquid clouds: snow number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+57,:),'dq_s_eme_rs',  'enhanced melting of ice by rain: snow content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+58,:),'dn_s_eme_rs',  'enhanced melting of ice by rain: snow number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+59,:),'dq_g_eme_gc',  'enhanced melting of graupel by liquid clouds: graupel content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+60,:),'dn_g_eme_gc',  'enhanced melting of graupel by liquid clouds: graupel number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+61,:),'dq_g_eme_gr',  'enhanced melting of graupel by rain: graupel content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+62,:),'dn_g_eme_gr',  'enhanced melting of graupel by rain: graupel number tendency',  '#/kg/s','tt')
        !
        ! - evaporation and melting
        ! -- melting
        call ncinfo(ncmname(nbaspout+63,:),'dq_i_me', 'partial melting of ice: ice content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+64,:),'dn_i_me', 'partial melting of ice: ice number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+65,:),'dq_s_me', 'partial melting of snow: snow content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+66,:),'dn_s_me', 'partial melting of snow: snow number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+67,:),'dq_g_me', 'partial melting of graupel: graupel content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+68,:),'dn_g_me', 'partial melting of graupel: graupel number tendency',  '#/kg/s','tt')
        ! -- evaporation
        call ncinfo(ncmname(nbaspout+69,:),'dq_i_ev', 'evaporation of ice: ice content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+70,:),'dn_i_ev', 'evaporation of ice: ice number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+71,:),'dq_s_ev', 'evaporation of snow: snow content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+72,:),'dn_s_ev', 'evaporation of snow: snow number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+73,:),'dq_g_ev', 'evaporation of graupel: graupel content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+74,:),'dn_g_ev', 'evaporation of graupel: graupel number tendency',  '#/kg/s','tt')
        !
        ! - conversions and other processes
        ! -- partial conversion to graupel
        call ncinfo(ncmname(nbaspout+75,:),'dq_i_cv_ig','partial conversion to graupel: ice content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+76,:),'dn_i_cv_ig','partial conversion to graupel: ice number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+77,:),'dq_s_cv_sg','partial conversion to graupel: snow content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+78,:),'dn_s_cv_sg','partial conversion to graupel: snow number tendency',  '#/kg/s','tt')
        !
        ! -- ice multiplication
        call ncinfo(ncmname(nbaspout+79,:),'dq_i_mul','ice multiplication: ice content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+80,:),'dn_i_mul','ice multiplication: ice number tendency',  '#/kg/s','tt')
        !
        ! -- impacty of saturation adjustment
        call ncinfo(ncmname(nbaspout+81,:),'dn_c_sadj', 'saturation adjustment: cloud water number tendency',  '#/kg/s','tt')
        !
        !-> add extra processes here
        !   -> then adjust nmphyout
        !
        ! - sedimentation tendency
        call ncinfo(ncmname(ncolout+1,:), 'dq_c_sed',  'cloud liquid water content sedimentation tendency','kg/kg/s','tt')
        call ncinfo(ncmname(ncolout+2,:), 'dq_r_sed',  'rain water content sedimentation tendency',        'kg/kg/s','tt')
        call ncinfo(ncmname(ncolout+3,:), 'dq_i_sed',  'cloud ice content sedimentation tendency',         'kg/kg/s','tt')
        call ncinfo(ncmname(ncolout+4,:), 'dq_s_sed',  'snow content sedimentation tendency',              'kg/kg/s','tt')
        call ncinfo(ncmname(ncolout+5,:), 'dq_g_sed',  'graupel content sedimentation tendency',           'kg/kg/s','tt')
        ! -- and for numbers
        call ncinfo(ncmname(ncolout+6,:), 'dn_c_sed',  'cloud droplet number sedimentation tendency',    '#/kg/s','tt')
        call ncinfo(ncmname(ncolout+7,:), 'dn_r_sed',  'raindrop number sedimentation tendency',         '#/kg/s','tt')
        call ncinfo(ncmname(ncolout+8,:), 'dn_i_sed',  'cloud ice number sedimentation tendency',        '#/kg/s','tt')
        call ncinfo(ncmname(ncolout+9,:), 'dn_s_sed',  'snow number sedimentation tendency',             '#/kg/s','tt')
        call ncinfo(ncmname(ncolout+10,:),'dn_g_sed',  'graupel number sedimentation tendency',          '#/kg/s','tt')
        !
        ! - total tendency
        ! - content of hydrometeors
        call ncinfo(ncmname(ncolout+11,:), 'dq_c_tot', 'total cloud liquid water content tendency','kg/kg/s','tt')
        call ncinfo(ncmname(ncolout+12,:), 'dq_r_tot', 'total rain water content tendency',        'kg/kg/s','tt')
        call ncinfo(ncmname(ncolout+13,:), 'dq_i_tot', 'total cloud ice content tendency',         'kg/kg/s','tt')
        call ncinfo(ncmname(ncolout+14,:), 'dq_s_tot', 'total snow content tendency',              'kg/kg/s','tt')
        call ncinfo(ncmname(ncolout+15,:), 'dq_g_tot', 'total graupel content tendency',           'kg/kg/s','tt')
        ! -- numbers of hydrometeors
        call ncinfo(ncmname(ncolout+16,:), 'dn_c_tot', 'total cloud droplet number tendency',    '#/kg/s','tt')
        call ncinfo(ncmname(ncolout+17,:), 'dn_r_tot', 'total raindrop number tendency',         '#/kg/s','tt')
        call ncinfo(ncmname(ncolout+18,:), 'dn_i_tot', 'total cloud ice number tendency',        '#/kg/s','tt')
        call ncinfo(ncmname(ncolout+19,:), 'dn_s_tot', 'total snow number tendency',             '#/kg/s','tt')
        call ncinfo(ncmname(ncolout+20,:), 'dn_g_tot', 'total graupel number tendency',          '#/kg/s','tt')
        ! -- numbers of CCN
        call ncinfo(ncmname(ncolout+21,:),'dn_ccn_tot','total ccn number tendency',              '#/kg/s','tt')
        ! -- temperature tendency
        call ncinfo(ncmname(nmphyout+1,:),'dthl_mphys','total microphysics liquid water temperature tendency', 'K/kg/s','tt')
        call ncinfo(ncmname(nmphyout+2,:),'dth_mphys', 'total microphysics potential temperature tendency',     'K/kg/s','tt')
        call ncinfo(ncmname(nmphyout+3,:),'dth_freeze','potential temperature tendency due to freezing',     'K/kg/s','tt')
        call ncinfo(ncmname(nmphyout+4,:),'dth_melt',  'potential temperature tendency due to melting',     'K/kg/s','tt')
        call ncinfo(ncmname(nmphyout+5,:),'dth_cond',  'potential temperature tendency due to condensation',     'K/kg/s','tt')
        call ncinfo(ncmname(nmphyout+6,:),'dth_ev',    'potential temperature tendency due to evaporation',     'K/kg/s','tt')
        call ncinfo(ncmname(nmphyout+7,:),'dth_dep',   'potential temperature tendency due to depostion',     'K/kg/s','tt')
        call ncinfo(ncmname(nmphyout+8,:),'dth_sub',   'potential temperature tendency due to sublimation',     'K/kg/s','tt')
        ! -- cloud fraction
        call ncinfo(ncmname(nmphyout+9,:),'cfrac_l',    'liquid cloud fraction','-','tt')
        call ncinfo(ncmname(nmphyout+10,:),'cfrac_i',   'ice cloud fraction','-','tt')
        call ncinfo(ncmname(nmphyout+11,:),'cfrac_tot', 'total cloud fraction','-','tt')
        ! -- precipitation
        call ncinfo(ncmname(nmphyout+12,:),'rain_rate', 'rain fall rate',   'kg/m2','tt')
        call ncinfo(ncmname(nmphyout+13,:),'ice_rate',  'ice-fall  rate',    'kg/m2','tt')
        call ncinfo(ncmname(nmphyout+14,:),'snow_rate', 'snow fall rate',    'kg/m2','tt')
        call ncinfo(ncmname(nmphyout+15,:),'graupel_rate', 'graupel fall rate',    'kg/m2','tt')

        call ncinfo(ncmname(nmphyout+16,:), 'dq_c_totf', 'total cloud liquid water content tendency','kg/kg/s','tt')
        call ncinfo(ncmname(nmphyout+17,:), 'dq_r_totf', 'total rain water content tendency',        'kg/kg/s','tt')
        call ncinfo(ncmname(nmphyout+18,:), 'dq_i_totf', 'total cloud ice content tendency',         'kg/kg/s','tt')
        call ncinfo(ncmname(nmphyout+19,:), 'dq_s_totf', 'total snow content tendency',              'kg/kg/s','tt')
        call ncinfo(ncmname(nmphyout+20,:), 'dq_g_totf', 'total graupel content tendency',           'kg/kg/s','tt')
        ! -- numbers of hydrometeors
        call ncinfo(ncmname(nmphyout+21,:), 'dn_c_totf', 'total cloud droplet number tendency',    '#/kg/s','tt')
        call ncinfo(ncmname(nmphyout+22,:), 'dn_r_totf', 'total raindrop number tendency',         '#/kg/s','tt')
        call ncinfo(ncmname(nmphyout+23,:), 'dn_i_totf', 'total cloud ice number tendency',        '#/kg/s','tt')
        call ncinfo(ncmname(nmphyout+24,:), 'dn_s_totf', 'total snow number tendency',             '#/kg/s','tt')
        call ncinfo(ncmname(nmphyout+25,:), 'dn_g_totf', 'total graupel number tendency',          '#/kg/s','tt')
        ! -- finish
        ! adding dimensions
        call open_nc(fname_mphys, ncid_mphys,nrec_mphys,n3=kmax)
        if (nrec_mphys == 0) then
          call define_nc( ncid_mphys, 1, tncmname)
          call writestat_dims_nc(ncid_mphys)
        end if
        call define_nc( ncid_mphys, NVar_mphys, ncmname)
      end if
   end if

  end subroutine initbulkmicrostat3


!------------------------------------------------------------------------------!
!> General routine, does the timekeeping
  subroutine bulkmicrostat3
    use modglobal,     only  : rk3step, timee, dt_lim
    implicit none

    if (.not. lmicrostat)  return
    if (rk3step /= 3)  return
    if (timee == 0)    return
    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if
    if (timee >= tnext) then
      tnext = tnext + idtav
      call dobulkmicrostat3   ! #sb3
    end if
    if (timee >= tnextwrite) then
      tnextwrite = tnextwrite + itimeav
      call writebulkmicrostat3   ! #sb3
    end if
  end subroutine bulkmicrostat3


!------------------------------------------------------------------------------!
!> Performs the calculations for rainrate etc.
  subroutine dobulkmicrostat3
    use modmpi,        only : my_real, mpi_sum, comm3d, mpierr
    use modglobal,     only : i1, j1, k1, ijtot
    use modmicrodata,  only : Dvr, epscloud, epsqr, epsprec
    use modmicrodata3, only : eps_hprec
    use modmicrodata3, only : precep_l,precep_i                       &
                             ,in_hr,in_cl,in_ci,in_hs,in_hg,in_cc     &
                             ,iq_hr,iq_cl,iq_ci,iq_hs,iq_hg
    use modfields,  only    : ql0, rhof, sv0
    implicit none

    integer      :: k

    cloudcountav  = 0.0
    raincountav   = 0.0
    Nrrainav      = 0.0
    qrav          = 0.0
    Dvrav         = 0.0

    cl_countavl   = 0.0
    ci_countavl   = 0.0
    hr_countavl   = 0.0
    hs_countavl   = 0.0
    hg_countavl   = 0.0

    ! and for large-scale variables
    cl_countav    = 0.0
    ci_countav    = 0.0
    hr_countav    = 0.0
    hs_countav    = 0.0
    hg_countav    = 0.0

    ! calculating the values
    do k = 1,k1
      cloudcountavl(k)  = count(ql0(2:i1,2:j1,k      ) > epscloud)
      raincountavl (k)  = count(sv0(2:i1,2:j1,k,iq_hr) > epsqr)
      Nrrainavl    (k)  = sum  (sv0(2:i1,2:j1,k,in_hr)) * rhof(k)
      qravl        (k)  = sum  (sv0(2:i1,2:j1,k,iq_hr))
      Dvravl       (k)  = sum  (Dvr(2:i1,2:j1,k), sv0(2:i1,2:j1,k,iq_hr) > epsqr)
    end do

    do k = 1,k1
      ! cloud column counts
      cl_countavl(k)  = count(sv0(2:i1,2:j1,k,iq_cl) > epscloud)
      ci_countavl(k)  = count(sv0(2:i1,2:j1,k,iq_ci) > epscloud)
      hr_countavl(k)  = count(sv0(2:i1,2:j1,k,iq_hr) > eps_hprec)
      hs_countavl(k)  = count(sv0(2:i1,2:j1,k,iq_hs) > eps_hprec)
      hg_countavl(k)  = count(sv0(2:i1,2:j1,k,iq_hg) > eps_hprec)
    end do

    ! transmitting to other processors
    call MPI_ALLREDUCE(cloudcountavl   ,cloudcountav,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(raincountavl    ,raincountav ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(Dvravl          ,Dvrav       ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(Nrrainavl       ,Nrrainav    ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qravl           ,qrav        ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(cl_countavl     ,cl_countav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(ci_countavl     ,ci_countav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(hr_countavl     ,hr_countav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(hs_countavl     ,hs_countav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(hg_countavl     ,hg_countav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)

    ! and normalising over all processors
    cloudcountmn = cloudcountmn  +  cloudcountav/ijtot
    raincountmn  = raincountmn   +  raincountav /ijtot
    Dvrmn        = Dvrmn         +  Dvrav       /ijtot
    Nrrainmn     = Nrrainmn      +  Nrrainav    /ijtot
    qrmn         = qrmn          +  qrav        /ijtot

    ! extra stats normalising over all processors
    cl_countmn = cl_countmn + cl_countav /ijtot
    ci_countmn = ci_countmn + ci_countav /ijtot
    hr_countmn = hr_countmn + hr_countav /ijtot
    hs_countmn = hs_countmn + hs_countav /ijtot
    hg_countmn = hg_countmn + hg_countav /ijtot
  end subroutine dobulkmicrostat3


!------------------------------------------------------------------------------!
!> Performs the calculations for the tendencies etc.
  subroutine bulkmicrotend3
  end subroutine bulkmicrotend3

!------------------------------------------------------------------------------!
!> Write the stats to file
  subroutine writebulkmicrostat3
!    use modmpi,        only : myid
!    use modglobal,     only : rtimee, ifoutput, cexpnr, k1,kmax,rlv,zf,cp
!    use modfields,     only : presf,rhof,exnf
!    use modstat_nc,    only : lnetcdf, writestat_nc
!    use modgenstat,    only : ncid_prof=>ncid, nrec_prof=>nrec
!    use modmicrodata3, only : rlme, rlvi
!
!    implicit none
!    real,dimension(k1,nvar) :: vars
!    real,allocatable,dimension(:,:) :: vars_mphys
!
!    integer    :: nsecs, nhrs, nminut
!    integer    :: k
!
!    ! allocate outut field
!    allocate(vars_mphys(k1,nvar_mphys))
!    allocate( cfrac_l  (k1)             &
!             ,cfrac_i  (k1)             &
!             ,cfrac_tot(k1)             )
!
!    cfrac_l      = 0.0
!    cfrac_i      = 0.0
!    cfrac_tot    = 0.0
!
!    ! TODO: these are all in the statistics_patch from modmicrodata3
!    dth_mphys    = 0.0
!    dth_freeze   = 0.0
!    dth_melt     = 0.0
!    dth_ev       = 0.0
!    dth_cond     = 0.0
!    dth_dep      = 0.0
!    dth_sub      = 0.0
!
!    nsecs  = nint(rtimee)
!    nhrs   = int (nsecs/3600)
!    nminut = int (nsecs/60)-nhrs*60
!    nsecs  = mod (nsecs,60)
!
!    cloudcountmn = cloudcountmn  /nsamples
!    raincountmn  = raincountmn   /nsamples
!    preccountmn  = preccountmn   /nsamples
!    prec_prcmn   = prec_prcmn    /nsamples
!    Dvrmn        = Dvrmn         /nsamples
!    Nrrainmn     = Nrrainmn      /nsamples
!    precmn       = precmn        /nsamples
!    qrmn         = qrmn          /nsamples
!
!    where (raincountmn > 0.)
!      Dvrmn        = Dvrmn / raincountmn
!    elsewhere
!      Dvrmn = 0.0
!    end where
!    where (preccountmn > 0.)
!      prec_prcmn = prec_prcmn/preccountmn
!    elsewhere
!      prec_prcmn = 0.0
!    end where
!
!    !  for cloud columns
!    cl_countmn = cl_countmn  /nsamples
!    ci_countmn = ci_countmn  /nsamples
!    hr_countmn = hr_countmn  /nsamples
!    hs_countmn = hs_countmn  /nsamples
!    hg_countmn = hg_countmn  /nsamples
!
!    ! for precipitation
!    prec_r_mn  = prec_r_mn   /nsamples
!    prec_i_mn  = prec_i_mn   /nsamples
!    prec_s_mn  = prec_s_mn   /nsamples
!    prec_g_mn  = prec_g_mn   /nsamples
!
!    !------------------------------------------
!    ! calculation of aggregated statistics
!    !------------------------------------------
!    ! cloud fractions
!    ! - liquid cloud fraction
!    cfrac_l = cl_countmn
!    ! - ice cloud fraction
!    cfrac_i = ci_countmn
!    ! - total cloud fraction
!    do k=1,k1
!      cfrac_tot(k) = max(cl_countmn(k),ci_countmn(k))
!      ! this is just an approximation
!      !-> later improve
!    enddo
!
!    ! process contributions
!
!    ! TODO: slabsum on statistics
!
!    if (myid == 0) then
!      !-----------------
!      ! text outputs
!      !--------------------
!      open (ifoutput,file='precep.'//cexpnr,position='append')
!      write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
!        '#-------------------------------------------------------------'   &
!        ,'---------------------------------)'           &
!        ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
!        ,nhrs,':',nminut,':',nsecs             &
!        ,'   HRS:MIN:SEC AFTER INITIALIZATION '
!      write (ifoutput,'(2A/A/A/2A/2A/2A)')             &
!        '#------------------------------------------------------------'     &
!        ,'------------'                 &
!        ,'#               --------   PRECIPITATION ------    '       &
!        ,'#                                                           '     &
!        ,'# LEV HEIGHT   RHO(k)  PRES  |CLOUDCOVER  ECHORAINRATE  PRECCOUNT '   &
!        ,'    NRRAIN      RAINCOUNT     PREC(k)     <Dvr(k)>     <qr(k)>'   &
!        ,'#      (M)             (MB)  |----------  ---W/M2----   --------- '   &
!        ,'    ------      ---------     -------     --------    ---------'   &
!        ,'#-----------------------------------------------------------------'   &
!        ,'---------------------------------------------------------------'
!      write(ifoutput,'(I4,F10.2,F8.3,F7.1,8E13.5)') &
!        (k                           , &
!        zf          (k)              , &
!        rhof        (k)              , &
!        presf       (k)/100.         , &
!        cloudcountmn(k)              , &
!        prec_prcmn  (k)*rhof(k)*rlv  , &
!        preccountmn (k)              , &
!        Nrrainmn    (k)              , &
!        raincountmn (k)              , &
!        precmn      (k)*rhof(k)*rlv  , &
!        Dvrmn       (k)              , &
!        qrmn        (k)              , &
!        k=1,kmax)
!      close(ifoutput)
!
!      open (ifoutput,file='nptend.'//cexpnr,position='append')
!      write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
!        '#-------------------------------------------------------------'   &
!        ,'---------------------------------)'           &
!        ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
!        ,nhrs,':',nminut,':',nsecs             &
!        ,'   HRS:MIN:SEC AFTER INITIALIZATION '
!      write (ifoutput,'(2A/A/A/2A/A/A)')             &
!        '#------------------------------------------------------------'     &
!        , '------------'               &
!        ,'#               --------   T E N D E N C I E S NRAIN ------    '     &
!        ,'#                                                           '     &
!        ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
!        ,'     EVAP         TOT '             &
!        ,'#      (M)   (MB)  |  ---------   (#/M3/S)      ----------'     &
!        ,'#-----------------------------------------------------------'
!      write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
!        (k                      , &
!        zf      (k)             , &
!        presf   (k)/100.        , &
!        Npmn    (k,iauto)       , &
!        Npmn    (k,iaccr)       , &
!        Npmn    (k,ised)        , &
!        Npmn    (k,ievap)       , &
!        sum(Npmn(k,2:nrfields)) , &
!        k=1,kmax)
!      close(ifoutput)
!
!      open (ifoutput,file='qlptend.'//cexpnr,position='append')
!      write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
!        '#-------------------------------------------------------------'   &
!        ,'---------------------------------)'           &
!        ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
!        ,nhrs,':',nminut,':',nsecs             &
!        ,'   HRS:MIN:SEC AFTER INITIALIZATION '
!      write (ifoutput,'(2A/A/A/2A/A/A)')             &
!        '#------------------------------------------------------------'     &
!        , '------------'               &
!        ,'#               --------   T E N D E N C I E S QRAIN ------    '   &
!        ,'#                                                           '     &
!        ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
!        ,'     EVAP         TOT '             &
!        ,'#      (M)   (MB)  |  ---------   (KG/KG/S)      ----------'     &
!        ,'#-----------------------------------------------------------'
!      write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
!        (k                         , &
!        zf       (k)               , &
!        presf    (k)/100.          , &
!        qlpmn    (k,iauto)         , &
!        qlpmn    (k,iaccr)         , &
!        qlpmn    (k,ised)          , &
!        qlpmn    (k,ievap)         , &
!        sum(qlpmn  (k,2:nrfields)) , &
!        k=1,kmax)
!      close(ifoutput)
!
!      open (ifoutput,file='qtptend.'//cexpnr,position='append')
!      write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
!        '#-------------------------------------------------------------'   &
!        ,'---------------------------------)'           &
!        ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
!        ,nhrs,':',nminut,':',nsecs             &
!        ,'   HRS:MIN:SEC AFTER INITIALIZATION '
!      write (ifoutput,'(2A/A/A/2A/A/A)')             &
!        '#------------------------------------------------------------'     &
!        , '------------'               &
!        ,'#               --------   T E N D E N C I E S QTP ------    '   &
!        ,'#                                                           '     &
!        ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
!        ,'     EVAP         TOT '             &
!        ,'#      (M)   (MB)  |  ---------   (KG/KG/S)      ----------'     &
!        ,'#-----------------------------------------------------------'
!      write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
!        (k                           , &
!        zf    (k)                    , &
!        presf (k)/100.               , &
!        qtpmn (k,iauto)              , &
!        qtpmn (k,iaccr)              , &
!        qtpmn (k,ised)               , &
!        qtpmn (k,ievap)              , &
!        sum   (qtpmn(k,2:nrfields))  , &
!        k=1,kmax)
!      close(ifoutput)
!
!      !------------------------------
!      ! NetCDF output
!      !------------------------------
!      if (lnetcdf) then
!        vars(:, 1) = cloudcountmn
!        vars(:, 2) = prec_prcmn*rhof(:)*rlv
!        vars(:, 3) = preccountmn 
!        vars(:, 4) = Nrrainmn    
!        vars(:, 5) = raincountmn 
!        vars(:, 6) = precmn*rhof(:)*rlv
!        vars(:, 7) = Dvrmn
!        vars(:, 8) = qrmn
!        vars(:, 9) = Npmn(:,iauto)
!        vars(:,10) = Npmn(:,iaccr)
!        vars(:,11) = Npmn(:,ised)
!        vars(:,12) = Npmn(:,ievap)
!        do k=1,k1
!        vars(k,13) =sum(Npmn  (k,2:nrfields))
!        enddo
!        vars(:,14) =qlpmn(:,iauto)
!        vars(:,15) =qlpmn(:,iaccr)
!        vars(:,16) =qlpmn(:,ised)
!        vars(:,17) =qlpmn(:,ievap)
!        do k=1,k1
!        vars(k,18) =sum(qlpmn  (k,2:nrfields))
!        enddo
!        vars(:,19) =qtpmn(:,iauto)
!        vars(:,20) =qtpmn(:,iaccr)
!        vars(:,21) =qtpmn(:,ised)
!        vars(:,22) =qtpmn(:,ievap)
!        do k=1,k1
!        vars(k,23) =sum(qtpmn  (k,2:nrfields))
!        enddo
!        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof,kmax)
!
!        ! recording added statistics
!        vars_mphys (:,nbaspout+1)  = dnc_nuc
!        vars_mphys (:,nbaspout+2)  = dqc_nuc
!        vars_mphys (:,nbaspout+3)  = dni_inuc
!        vars_mphys (:,nbaspout+4)  = dqi_inuc
!
!        ! - deposition and sublimation
!        vars_mphys (:,nbaspout+5)  = dqi_dep
!        vars_mphys (:,nbaspout+6)  = dqs_dep
!        vars_mphys (:,nbaspout+7)  = dqg_dep
!        vars_mphys (:,nbaspout+8)  = dqi_sub
!        vars_mphys (:,nbaspout+9)  = dqs_sub
!        vars_mphys (:,nbaspout+10) = dqg_sub
!
!        ! - freezing
!        ! -- homogeneous freezing
!        vars_mphys (:,nbaspout+11) = -dqc_hom ! dqi_hom
!        vars_mphys (:,nbaspout+12) = -dnc_hom ! dni_hom
!
!        ! -- heterogeneous freezing
!        vars_mphys (:,nbaspout+13) = -dqc_het ! dqi_het
!        vars_mphys (:,nbaspout+14) = -dnc_het ! dni_het
!        vars_mphys (:,nbaspout+15) = -dqr_het ! dqg_het
!        vars_mphys (:,nbaspout+16) = -dnr_het ! dng_het
!
!        ! - warm microphysics
!        ! - saturation aqdjustment tendency
!        vars_mphys (:,nbaspout+17) = dqc_sadjpl
!        vars_mphys (:,nbaspout+18) = dqc_sadjneg
!
!        ! -- cloud water self-collection
!        vars_mphys (:,nbaspout+19) = dnc_sc
!
!        ! -- cloud water autoconversion
!        vars_mphys (:,nbaspout+20) = -dqr_au ! dqc_au
!        vars_mphys (:,nbaspout+21) = dnc_au
!
!        ! -- cloud water accretion
!        vars_mphys (:,nbaspout+22) = -dqr_ac !dqc_ac
!        vars_mphys (:,nbaspout+23) = dnc_ac
!
!        ! -- rain selfcollection and breakup
!        vars_mphys (:,nbaspout+24) = dnr_sc
!        vars_mphys (:,nbaspout+25) = dnr_br
!
!        ! -- rain evaporation
!        vars_mphys (:,nbaspout+26) = dqr_ev
!        vars_mphys (:,nbaspout+27) = dnr_ev
!
!        ! - cold aggregation and self-collection
!        ! -- ice aggregation
!        vars_mphys (:,nbaspout+28) = dqi_agg
!        vars_mphys (:,nbaspout+29) = dni_agg
!
!        ! -- snow self-collection
!        vars_mphys (:,nbaspout+30) = dns_sc
!
!        ! - loss of liquid water due to riming
!        vars_mphys (:,nbaspout+31) = dqc_rime
!        vars_mphys (:,nbaspout+32) = dnc_rime
!        vars_mphys (:,nbaspout+33) = dqr_rime
!        vars_mphys (:,nbaspout+34) = dnr_rime
!
!        ! - growth by riming
!        vars_mphys (:,nbaspout+35) = dqi_rime_ic
!        vars_mphys (:,nbaspout+36) = dqs_rime_sc
!
!        !vars_mphys (:,nbaspout+37) = dqs_rime_sr
!        vars_mphys (:,nbaspout+37) = dqg_rime_gc
!        vars_mphys (:,nbaspout+38) = dqg_rime_gr
!
!        ! - collisions
!        ! -- collision conversions
!        vars_mphys (:,nbaspout+39) = dqi_col_rig
!        vars_mphys (:,nbaspout+40) = dni_col_rig
!        vars_mphys (:,nbaspout+41) = dqr_col_rig
!        vars_mphys (:,nbaspout+42) = dnr_col_rig
!
!        vars_mphys (:,nbaspout+43) = dqs_col_rsg
!        vars_mphys (:,nbaspout+44) = dns_col_rsg
!        vars_mphys (:,nbaspout+45) = dqr_col_rsg
!        vars_mphys (:,nbaspout+46) = dnr_col_rsg
!
!        ! -- collision collection
!        vars_mphys (:,nbaspout+47) = -dqs_col_si
!        vars_mphys (:,nbaspout+48) = dni_col_si
!        vars_mphys (:,nbaspout+49) = dqs_col_gsg
!        vars_mphys (:,nbaspout+50) = dns_col_gsg
!
!        ! -- enhanced melting
!        vars_mphys (:,nbaspout+51) = dqi_eme_ic
!        vars_mphys (:,nbaspout+52) = dni_eme_ic
!        vars_mphys (:,nbaspout+53) = dqi_eme_ri
!        vars_mphys (:,nbaspout+54) = dni_eme_ri
!        vars_mphys (:,nbaspout+55) = dqs_eme_sc
!        vars_mphys (:,nbaspout+56) = dns_eme_sc
!        vars_mphys (:,nbaspout+57) = dqs_eme_rs
!        vars_mphys (:,nbaspout+58) = dns_eme_rs
!        vars_mphys (:,nbaspout+59) = dqg_eme_gc
!        vars_mphys (:,nbaspout+60) = dng_eme_gc
!        vars_mphys (:,nbaspout+61) = dqg_eme_gr
!        vars_mphys (:,nbaspout+62) = dng_eme_gr
!
!        ! - evaporation and melting
!        ! -- melting
!        vars_mphys (:,nbaspout+63) = dqi_me
!        vars_mphys (:,nbaspout+64) = dni_me
!        vars_mphys (:,nbaspout+65) = dqs_me
!        vars_mphys (:,nbaspout+66) = dns_me
!        vars_mphys (:,nbaspout+67) = dqg_me
!        vars_mphys (:,nbaspout+68) = dng_me
!
!        ! -- evaporation
!        vars_mphys (:,nbaspout+69) = dqi_ev
!        vars_mphys (:,nbaspout+70) = dni_ev
!        vars_mphys (:,nbaspout+71) = dqs_ev
!        vars_mphys (:,nbaspout+72) = dns_ev
!        vars_mphys (:,nbaspout+73) = dqg_ev
!        vars_mphys (:,nbaspout+74) = dng_ev
!
!        ! - conversions and other processes
!        vars_mphys (:,nbaspout+75) = dqi_cv_ig
!        vars_mphys (:,nbaspout+76) = dni_cv_ig
!        vars_mphys (:,nbaspout+77) = dqs_cv_sg
!        vars_mphys (:,nbaspout+78) = dns_cv_sg
!
!        ! -- ice multiplication
!        vars_mphys (:,nbaspout+79) = dqi_mul
!        vars_mphys (:,nbaspout+80) = dni_mul
!
!        ! -- satuartion adjustments
!        vars_mphys (:,nbaspout+81) = dnc_sadj
!
!        !-> add extra processes here
!        !   then adjust nmphyout
!
!        ! - sedimentation tendency
!        vars_mphys (:,ncolout+1)  = dqc_sed
!        vars_mphys (:,ncolout+2)  = dqr_sed
!        vars_mphys (:,ncolout+3)  = dqi_sed
!        vars_mphys (:,ncolout+4)  = dqs_sed
!        vars_mphys (:,ncolout+5)  = dqg_sed
!
!        ! -- and for numbers
!        vars_mphys (:,ncolout+6)  = dnc_sed
!        vars_mphys (:,ncolout+7)  = dnr_sed
!        vars_mphys (:,ncolout+8)  = dni_sed
!        vars_mphys (:,ncolout+9)  = dns_sed
!        vars_mphys (:,ncolout+10) = dng_sed
!
!        ! - total tendency
!        ! - content of hydrometeors
!        vars_mphys (:,ncolout+11) = dqc_tot
!        vars_mphys (:,ncolout+12) = dqr_tot
!        vars_mphys (:,ncolout+13) = dqi_tot
!        vars_mphys (:,ncolout+14) = dqs_tot
!        vars_mphys (:,ncolout+15) = dqg_tot
!
!        ! -- numbers of hydrometeors
!        vars_mphys (:,ncolout+16) = dnc_tot
!        vars_mphys (:,ncolout+17) = dnr_tot
!        vars_mphys (:,ncolout+18) = dni_tot
!        vars_mphys (:,ncolout+19) = dns_tot
!        vars_mphys (:,ncolout+20) = dng_tot
!
!        ! -- numbers of CCN
!        vars_mphys (:,ncolout+21) = dn_ccn_tot
!
!        ! -- temperature tendency
!        vars_mphys (:,nmphyout+1) = dthl_mphys
!        vars_mphys (:,nmphyout+2) = dth_mphys
!        vars_mphys (:,nmphyout+3) = dth_freeze
!        vars_mphys (:,nmphyout+4) = dth_melt
!        vars_mphys (:,nmphyout+5) = dth_cond
!        vars_mphys (:,nmphyout+6) = dth_ev
!        vars_mphys (:,nmphyout+7) = dth_dep
!        vars_mphys (:,nmphyout+8) = dth_sub
!
!        ! -- cloud fraction
!        vars_mphys (:,nmphyout+9) = cfrac_l
!        vars_mphys (:,nmphyout+10) = cfrac_i
!        vars_mphys (:,nmphyout+11) = cfrac_tot
!
!        ! -- precipitation
!        vars_mphys (:,nmphyout+12) = prec_r_mn ! rain_rate
!        vars_mphys (:,nmphyout+13) = prec_i_mn ! ice_rate
!        vars_mphys (:,nmphyout+14) = prec_s_mn ! snow_rate
!        vars_mphys (:,nmphyout+15) = prec_g_mn ! graupel_rate
!
!        ! - content of hydrometeors
!        vars_mphys (:,nmphyout+16) = dqc_totf
!        vars_mphys (:,nmphyout+17) = dqr_totf
!        vars_mphys (:,nmphyout+18) = dqi_totf
!        vars_mphys (:,nmphyout+19) = dqs_totf
!        vars_mphys (:,nmphyout+20) = dqg_totf
!
!        ! -- numbers of hydrometeors
!        vars_mphys (:,nmphyout+21) = dnc_totf
!        vars_mphys (:,nmphyout+22) = dnr_totf
!        vars_mphys (:,nmphyout+23) = dni_totf
!        vars_mphys (:,nmphyout+24) = dns_totf
!        vars_mphys (:,nmphyout+25) = dng_totf
!
!        call writestat_nc(ncid_mphys,1,tncmname,(/rtimee/),nrec_mphys,.true.)
!        call writestat_nc(ncid_mphys,nvar_mphys,ncmname,vars_mphys(1:kmax,:),nrec_mphys,kmax)
!      end if ! lnetcdf
!
!    end if ! myid == 0
!
!    !------------------------------
!    ! reseting variables
!    !------------------------------
!
!    cloudcountmn= 0.0
!    raincountmn = 0.0
!    preccountmn = 0.0
!    prec_prcmn  = 0.0
!    Dvrmn       = 0.0
!    Nrrainmn    = 0.0
!    precmn      = 0.0
!    qrmn        = 0.0
!    Npmn        = 0.0
!    qlpmn       = 0.0
!    qtpmn       = 0.0
!
!    cl_countmn  = 0.0
!    ci_countmn  = 0.0
!    hr_countmn  = 0.0
!    hs_countmn  = 0.0
!    hg_countmn  = 0.0
!    prec_r_mn   = 0.0
!    prec_i_mn   = 0.0
!    prec_s_mn   = 0.0
!    prec_g_mn   = 0.0
!
!    ! and for variables describing processes
!    ! TODO: zero the tendencies
!
!    deallocate( vars_mphys )
!    deallocate( cfrac_l,cfrac_i,cfrac_tot                    &
!               ,dqc_rime,dqr_rime                            &
!               ,dnc_rime,dnr_rime                            &
!               ,dth_mphys,dth_freeze,dth_melt                &
!               ,dth_ev,dth_cond,dth_dep,dth_sub              )
!
  end subroutine writebulkmicrostat3

!------------------------------------------------------------------------------!

  subroutine exitbulkmicrostat3
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none
    if (.not. lmicrostat)  return

    if(lnetcdf .and. myid==0) call exitstat_nc(ncid_mphys)
    ! profile file closed separately

    deallocate(Npav     , &
               Npmn     , &
               qlpav    , &
               qlpmn    , &
               qtpav    , &
               qtpmn      )
    deallocate(precmn       , &
               preccountmn  , &
               prec_prcmn   , &
               cloudcountavl, &
               cloudcountav , &
               cloudcountmn , &
               raincountavl , &
               raincountav  , &
               raincountmn  , &
               Nrrainavl    , &
               Nrrainav     , &
               Nrrainmn     , &
               qravl        , &
               qrav         , &
               qrmn         , &
               Dvravl       , &
               Dvrav        , &
               Dvrmn)

    ! deallocating new variables
    deallocate( cl_countmn,hr_countmn                          &
               ,ci_countmn,hs_countmn,hg_countmn               &
               ,prec_r_mn,prec_i_mn,prec_s_mn,prec_g_mn        )
    deallocate( cl_countavl, hr_countavl                       &
               ,ci_countavl,hs_countavl,hg_countavl            )
    deallocate( cl_countav,hr_countav                          &
               ,ci_countav,hs_countav,hg_countav               )

  end subroutine exitbulkmicrostat3



!------------------------------------------------------------------------------

end module
