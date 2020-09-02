!> \file modbulkmicrostat3.f90
!!  Calculates profiles coming from the bulkmicrophysics
!>
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
  integer,parameter :: hdump_nvar = 10 ! 8
  integer :: nvar_mphys = 0 ! adjusted later
  integer :: ncid_mphys = 0
  integer :: nrec_mphys = 0
  integer :: hdump_ncid,hdump_nrec = 0
  character(80) :: fname_mphys = 'mphysprofiles.xxx.nc'
  character(80) :: hdump_fname = 'hydrodump.xxx.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  character(80),allocatable, dimension(:,:) :: ncmname
  character(80),dimension(1,4) :: tncmname
  character(80),dimension(hdump_nvar,4) :: hdump_ncname
  character(80),dimension(1,4) :: hdump_tncname
  integer      ,dimension(hdump_nvar) :: isv_vec
  real          :: dtav, timeav,hdump_dtav
  integer(kind=longint):: idtav, itimeav, tnext, tnextwrite
  integer(kind=longint):: hdump_idtav, hdump_tnext
  integer          :: nsamples, hdump_klow, hdump_khigh
  logical          :: lmicrostat = .false.
  integer, parameter      :: nrfields = 5  , &
                 iauto    = 2 , &
                  iaccr    = 3 , &
               ievap    = 4 , &
               ised      = 5
  integer          :: nmphyout,nbaspout,ncolout   ! #sb3
  real, allocatable, dimension(:,:)  :: Npav    , &
               Npmn    , &
               qlpav  , &
               qlpmn  , &
               qtpav  , &
               qtpmn
  real, allocatable, dimension(:)    :: precavl  , &
               precav  , &
               precmn  , &
               preccountavl  , &
               preccountav  , &
               preccountmn  , &
               prec_prcavl  , &
               prec_prcav  , &
               prec_prcmn  , &
               cloudcountavl, &
               cloudcountav  , &
               cloudcountmn  , &
               raincountavl  , &
               raincountav  , &
               raincountmn  , &
               Nrrainavl  , &
               Nrrainav  , &
               Nrrainmn  , &
               qravl  , &
               qrav    , &
               qrmn    , &
               Dvravl  , &
               Dvrav  , &
               Dvrmn  
    ! #sb3 start extend here onwards 
    real, allocatable, dimension(:) :: dnc_nuc ,dqc_nuc       &
              ,dni_inuc ,dqi_inuc                             &
              ,dqi_dep ,dqi_sub ,dqs_dep ,dqs_sub             &
              ,dqg_dep ,dqg_sub                               &
              ,dnc_hom ,dqc_hom                               &
              ,dnc_het ,dqc_het ,dnr_het ,dqr_het             &
              ,dnc_sadj ,dqc_sadjpl, dqc_sadjneg              &
              ,dnc_sc                                         &
              ,dnc_au ,dnr_au  ,dqr_au                        &
              ,dqr_ac , dnc_ac                                &
              ,dnr_sc ,dnr_br                                 &
              ,dqr_ev ,dnr_ev                             &
              ,dqi_agg ,dni_agg                               &
              ,dns_sc                                         &
              ,dqi_rime_ic ,dqs_rime_sc ,dqg_rime_gc          &
              ,dnc_rime_ic ,dnc_rime_sc ,dnc_rime_gc          &
              ,dqg_rime_gr,dnr_rime_gr                        &
              ,dqi_col_rig ,dqr_col_rig                       &
              ,dni_col_rig ,dnr_col_rig                       &
              ,dqs_col_rsg ,dqr_col_rsg                       &
              ,dns_col_rsg ,dnr_col_rsg                       &              
              ,dqs_col_si ,dni_col_si                         &
              ,dqs_col_gsg ,dns_col_gsg                       &
              ,dqi_cv_ig ,dni_cv_ig ,dqs_cv_sg ,dns_cv_sg     &
              ,dqc_sed,dqr_sed,dqi_sed,dqs_sed,dqg_sed        &
              ,dnc_sed,dnr_sed,dni_sed,dns_sed,dng_sed        &
              ,dqc_tot,dqi_tot,dqr_tot,dqs_tot,dqg_tot        &
              ,dnc_tot,dni_tot,dnr_tot,dns_tot, dng_tot       &
              ,dqc_totf,dqi_totf,dqr_totf,dqs_totf,dqg_totf   & ! #
              ,dnc_totf,dni_totf,dnr_totf,dns_totf, dng_totf  & ! #            
              ,dn_ccn_tot                                     &
              ,dthl_mphys,dth_mphys,dth_freeze,dth_melt       &
              ,dth_ev,dth_cond,dth_dep,dth_sub               &
              ,dqi_eme_ic,dqs_eme_sc,dqg_eme_gc               &
              ,dni_eme_ic,dns_eme_sc,dng_eme_gc               &
              ,dqi_eme_ri,dqs_eme_rs,dqg_eme_gr               &
              ,dni_eme_ri,dns_eme_rs,dng_eme_gr               &
              ,dqi_me ,dqs_me ,dqg_me                         & 
              ,dni_me ,dns_me ,dng_me                         &
              ,dqi_ev,dqs_ev,dqg_ev                        & 
              ,dni_ev,dns_ev,dng_ev                        &
              ,dni_mul,dqi_mul                                &
              ,cl_countmn,hr_countmn                          &
              ,ci_countmn,hs_countmn,hg_countmn               &
              ,prec_r_mn,prec_i_mn,prec_s_mn,prec_g_mn
              !-> continue here if needed           
              ! ,cfrac_l,cfrac_i,cfrac_tot    <- in writebulkmicrostat3
              ! ,rain_rate,ice_rate,snow_rate,graupel_rate      &
      real, allocatable, dimension(:) :: dnc_nuc_av,dqc_nuc_av &
              ,dni_inuc_av,dqi_inuc_av                         &
              ,dqi_dep_av,dqi_sub_av,dqs_dep_av,dqs_sub_av     &
              ,dqg_dep_av,dqg_sub_av                           &
              ,dnc_hom_av,dqc_hom_av                           &
              ,dnc_het_av,dqc_het_av,dnr_het_av,dqr_het_av     &
              ,dnc_sadj_av,dqc_sadjpl_av, dqc_sadjneg_av       &
              ,dnc_sc_av                                       &
              ,dnc_au_av,dnr_au_av,dqr_au_av                   &
              ,dqr_ac_av,dnc_ac_av                             &
              ,dnr_sc_av,dnr_br_av                             &
              ,dqr_ev_av,dnr_ev_av                         &
              ,dqi_agg_av,dni_agg_av                           &
              ,dns_sc_av                                       &
              ,dqi_rime_icav,dqs_rime_scav,dqg_rime_gcav       &
              ,dnc_rime_icav,dnc_rime_scav,dnc_rime_gcav       &
              ,dqg_rime_grav,dnr_rime_grav                     &
              ,dqi_col_rigav,dqr_col_rigav                     &
              ,dni_col_rigav,dnr_col_rigav                     &
              ,dqs_col_rsgav,dqr_col_rsgav                     &
              ,dns_col_rsgav,dnr_col_rsgav                     &              
              ,dqs_col_siav, dni_col_siav                      &
              ,dqs_col_gsgav,dns_col_gsgav                     &
              ,dqi_cv_igav,dni_cv_igav,dqs_cv_sgav,dns_cv_sgav &
              ,dn_ccn_av                                       &
              ,dthl_mphys_av                                   &
              ,dqc_sedav,dqr_sedav                             &
              ,dqi_sedav,dqs_sedav,dqg_sedav                   &
              ,dnc_sedav,dnr_sedav                             &
              ,dni_sedav,dns_sedav,dng_sedav                   &
              ,dqc_totav,dqr_totav                             &
              ,dqi_totav,dqs_totav,dqg_totav                   &
              ,dnc_totav, dnr_totav                            &
              ,dni_totav,dns_totav, dng_totav                  &
              ,dqc_totfav,dqr_totfav                             & ! #
              ,dqi_totfav,dqs_totfav,dqg_totfav                   & ! #
              ,dnc_totfav, dnr_totfav                            & ! #
              ,dni_totfav,dns_totfav, dng_totfav                  & ! #             
              ,dqi_eme_icav,dqs_eme_scav,dqg_eme_gcav          &
              ,dni_eme_icav,dns_eme_scav,dng_eme_gcav          &
              ,dqi_eme_riav,dqs_eme_rsav,dqg_eme_grav          &
              ,dni_eme_riav,dns_eme_rsav,dng_eme_grav          &
              ,dqi_me_av ,dqs_me_av ,dqg_me_av                 & 
              ,dni_me_av ,dns_me_av ,dng_me_av                 &
              ,dqi_ev_av,dqs_ev_av,dqg_ev_av                & 
              ,dni_ev_av,dns_ev_av,dng_ev_av                &
              ,dni_mul_av,dqi_mul_av                             &
              ,cl_countav,hr_countav                           &
              ,ci_countav,hs_countav,hg_countav                &
              ,prec_r_av,prec_i_av,prec_s_av,prec_g_av         
      real, allocatable, dimension(:) :: cl_countavl           &
            ,hr_countavl,ci_countavl,hs_countavl,hg_countavl   &
            ,prec_r_avl,prec_i_avl,prec_s_avl,prec_g_avl   
      real, allocatable, dimension(:) ::  cfrac_tot            &
                   ,cfrac_l,cfrac_i                            &
                   ,dqc_rime,dqr_rime,dnc_rime,dnr_rime        
                 !  ,dqi_col_si                                 &
                 !  ,dqi_rime,dqs_rime,dqg_rime                 &
                 !  ,dqi_col,dqs_col ,dqg_col                   &
                 !  ,dni_col,dns_col ,dng_col                   &
                 !  ,dqi_cv,dqs_cv,dqg_cv                       &
                 !  ,dni_cv,dns_cv,dng_cv    
     ! #sb3 END                 
               

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
    use modmicrodata3, only : l_hdump ! #sb3
    implicit none
    integer      :: ierr

    ! read nambulkmicrostat instead of initbulkmicrostat in modbulkmicro.f90
    namelist/NAMBULKMICROSTAT/ &
    lmicrostat, dtav, timeav
    

    ! if ((imicro /=imicro_bulk) .and. (imicro /= imicro_sice)) return #sb3
    ! - since we want to generate statistics even for imicro_bulk3 or other

    dtav  = dtav_glob
    timeav  = timeav_glob
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
      
      ! whether to write hydrodumps 
      ! if (lfielddump) then ! <- not possible since not public
      l_hdump = .true.
      ! endif
    endif 
    
    ! #sb3 START  number of basic mphys outputs
    ! adjust these number when needed 
    nbaspout = 0 ! 23               ! - to describe the number of basic mphys outputs
    ncolout  = nbaspout +81     ! - to describe number of microphysical outpus up to collisions 
    nmphyout = ncolout+21       !  = ncolout ! = ncolout+36       ! - to describe number of microphysical outputs before total ones
    nvar_mphys = nmphyout+25   ! then no adjustment of nvar_mphys needed
    ! #sb3 END 
    

    ! transmitting setting to other processors
    call MPI_BCAST(lmicrostat  ,1,MPI_LOGICAL  ,0,comm3d,mpierr)
    call MPI_BCAST(l_hdump  ,1,MPI_LOGICAL  ,0,comm3d,mpierr)
    call MPI_BCAST(dtav    ,1,MY_REAL  ,0,comm3d,mpierr)
    call MPI_BCAST(timeav    ,1,MY_REAL  ,0,comm3d,mpierr)
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav


    if (.not. lmicrostat) return
    if (abs(timeav/dtav - nsamples) > 1e-4) then
      stop 'timeav must be an integer multiple of dtav (NAMBULKMICROSTAT)'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax - nint(dtav/dtmax)) > 1e-4) then
      stop 'dtav must be an integer multiple of dtmax (NAMBULKMICROSTAT)'
    end if
    
    ! calling initialisation of hydrometeor dump
    if(l_hdump) then 
      call inithydrodump
    endif
    
!#t 
    ! allocating array for statistics that should be calculated
    allocate(Npav    (k1, nrfields)  , &
       Npmn    (k1, nrfields)  , &
       qlpav    (k1, nrfields)  , &
       qlpmn    (k1, nrfields)  , &
       qtpav    (k1, nrfields)  , &
       qtpmn    (k1, nrfields)  )
    allocate(precavl  (k1)    , &
       precav    (k1)    , &
       precmn    (k1)    , &
       preccountavl  (k1)    , &
       preccountav  (k1)    , &
       preccountmn  (k1)    , &
       prec_prcavl  (k1)    , &
       prec_prcav  (k1)    , &
       prec_prcmn  (k1)    , &
       cloudcountavl  (k1)    , &
       cloudcountav  (k1)    , &
       cloudcountmn  (k1)    , &
       raincountavl  (k1)    , &
       raincountav  (k1)    , &
       raincountmn  (k1)    , &
       Nrrainavl  (k1)    , &
       Nrrainav  (k1)    , &
       Nrrainmn  (k1)    , &
       qravl    (k1)    , &
       qrav    (k1)    , &
       qrmn    (k1)    , &
       Dvravl    (k1)    , &
       Dvrav    (k1)    , &
       Dvrmn    (k1))
    ! #sb3 START
     allocate( cl_countavl (k1)                      &               
              ,ci_countavl (k1)                      &
              ,hr_countavl (k1)                      &
              ,hs_countavl (k1)                      &
              ,hg_countavl (k1)                      &
              ,prec_r_avl  (k1)                      &
              ,prec_i_avl  (k1)                      &
              ,prec_s_avl  (k1)                      &
              ,prec_g_avl  (k1)                      &
             )
     allocate(cl_countav (k1)        &
             ,ci_countav (k1)        &
             ,hr_countav (k1)        &
             ,hs_countav (k1)        &
             ,hg_countav (k1)        &
             ,prec_r_av  (k1)        &
             ,prec_i_av  (k1)        &
             ,prec_s_av  (k1)        &
             ,prec_g_av  (k1)        &
             ,cl_countmn (k1)        &
             ,ci_countmn (k1)        &
             ,hr_countmn (k1)        &
             ,hs_countmn (k1)        &
             ,hg_countmn (k1)        &
             ,prec_r_mn  (k1)        &
             ,prec_i_mn  (k1)        &
             ,prec_s_mn  (k1)        &
             ,prec_g_mn  (k1)        &
            )
     allocate(dnc_nuc_av    (k1)     & 
             ,dqc_nuc_av    (k1)     &
             ,dni_inuc_av   (k1)     &
             ,dqi_inuc_av   (k1)     &
             ,dqi_dep_av    (k1)     &
             ,dqi_sub_av    (k1)     &
             ,dqs_dep_av    (k1)     &
             ,dqs_sub_av    (k1)     &
             ,dqg_dep_av    (k1)     &
             ,dqg_sub_av    (k1)     &
             ,dnc_hom_av    (k1)     &
             ,dqc_hom_av    (k1)     &
             ,dnc_het_av    (k1)     &
             ,dqc_het_av    (k1)     &
             ,dnr_het_av    (k1)     &
             ,dqr_het_av    (k1)     &
             ,dnc_sadj_av   (k1)     &
             ,dqc_sadjpl_av (k1)     &
             ,dqc_sadjneg_av(k1)     &
             ,dnc_sc_av     (k1)     &
             ,dnc_au_av     (k1)     &
             ,dnr_au_av     (k1)     &
             ,dqr_au_av     (k1)     &
             ,dqr_ac_av     (k1)     &
             ,dnc_ac_av     (k1)     &
             ,dnr_sc_av     (k1)     &
             ,dnr_br_av     (k1)     &
             ,dqr_ev_av   (k1)     &
             ,dnr_ev_av   (k1)     &
             ,dqi_agg_av    (k1)     &
             ,dni_agg_av    (k1)     &
             ,dns_sc_av     (k1)     &
             ,dqi_rime_icav (k1)     &
             ,dqs_rime_scav (k1)     &
             ,dqg_rime_gcav (k1)     &
             ,dnc_rime_icav (k1)     &
             ,dnc_rime_scav (k1)     &
             ,dnc_rime_gcav (k1)     &
             ,dqg_rime_grav (k1)     &
             ,dnr_rime_grav (k1)     &
             ,dqi_col_rigav (k1)     &
             ,dqr_col_rigav (k1)     &
             ,dni_col_rigav (k1)     &
             ,dnr_col_rigav (k1)     &
             ,dqs_col_rsgav (k1)     &
             ,dqr_col_rsgav (k1)     &
             ,dns_col_rsgav (k1)     &
             ,dnr_col_rsgav (k1)     &             
             ,dqs_col_siav  (k1)     &
             ,dni_col_siav  (k1)     &
             ,dqs_col_gsgav (k1)     & 
             ,dns_col_gsgav (k1)     & 
             ,dqi_cv_igav   (k1)     &
             ,dni_cv_igav   (k1)     &
             ,dqs_cv_sgav   (k1)     &
             ,dns_cv_sgav   (k1)     &
             ,dqc_sedav     (k1)     &
             ,dqi_sedav     (k1)     &
             ,dqr_sedav     (k1)     &
             ,dqs_sedav     (k1)     &
             ,dqg_sedav     (k1)     &
             ,dnc_sedav     (k1)     &
             ,dni_sedav     (k1)     &
             ,dnr_sedav     (k1)     &
             ,dns_sedav     (k1)     &
             ,dng_sedav     (k1)     & 
             ,dqc_totav     (k1)     & 
             ,dqi_totav     (k1)     & 
             ,dqr_totav     (k1)     & 
             ,dqs_totav     (k1)     & 
             ,dqg_totav     (k1)     & 
             ,dnc_totav     (k1)     & 
             ,dni_totav     (k1)     & 
             ,dnr_totav     (k1)     & 
             ,dns_totav     (k1)     & 
             ,dng_totav     (k1)     & 
             ,dn_ccn_av     (k1)     &
             ,dthl_mphys_av (k1)     &
             ,dqi_eme_icav  (k1)     &
             ,dqs_eme_scav  (k1)     &
             ,dqg_eme_gcav  (k1)     &
             ,dni_eme_icav  (k1)     &
             ,dns_eme_scav  (k1)     &
             ,dng_eme_gcav  (k1)     &
             ,dqi_eme_riav  (k1)     &
             ,dqs_eme_rsav  (k1)     &
             ,dqg_eme_grav  (k1)     &
             ,dni_eme_riav  (k1)     &
             ,dns_eme_rsav  (k1)     &
             ,dng_eme_grav  (k1)     &
             ,dqi_me_av     (k1)     &
             ,dqs_me_av     (k1)     &
             ,dqg_me_av     (k1)     &
             ,dni_me_av     (k1)     &
             ,dns_me_av     (k1)     &
             ,dng_me_av     (k1)     &
             ,dqi_ev_av    (k1)     &
             ,dqs_ev_av    (k1)     &
             ,dqg_ev_av    (k1)     &
             ,dni_ev_av    (k1)     &
             ,dns_ev_av    (k1)     &
             ,dng_ev_av    (k1)     &
             ,dni_mul_av    (k1)     &
             ,dqi_mul_av    (k1)     &
             ,dqc_totfav     (k1)     & 
             ,dqi_totfav     (k1)     & 
             ,dqr_totfav     (k1)     & 
             ,dqs_totfav     (k1)     & 
             ,dqg_totfav     (k1)     & 
             ,dnc_totfav     (k1)     & 
             ,dni_totfav     (k1)     & 
             ,dnr_totfav     (k1)     & 
             ,dns_totfav     (k1)     & 
             ,dng_totfav     (k1)     &             
             
            )
     allocate(  dnc_nuc     (k1)     & 
               ,dqc_nuc     (k1)     &
               ,dni_inuc    (k1)     &
               ,dqi_inuc    (k1)     &
               ,dqi_dep     (k1)     &
               ,dqi_sub     (k1)     &
               ,dqs_dep     (k1)     &
               ,dqs_sub     (k1)     &
               ,dqg_dep     (k1)     &
               ,dqg_sub     (k1)     &
               ,dnc_hom     (k1)     &
               ,dqc_hom     (k1)     &
               ,dnc_het     (k1)     &
               ,dqc_het     (k1)     &
               ,dnr_het     (k1)     &
               ,dqr_het     (k1)     &
               ,dnc_sadj    (k1)     &
               ,dqc_sadjpl  (k1)     &
               ,dqc_sadjneg (k1)     &
               ,dnc_sc      (k1)     &
               ,dnc_au      (k1)     &
               ,dnr_au      (k1)     &
               ,dqr_au      (k1)     &
               ,dqr_ac      (k1)     &
               ,dnc_ac      (k1)     &
               ,dnr_sc      (k1)     &
               ,dnr_br      (k1)     &
               ,dqr_ev    (k1)     &
               ,dnr_ev    (k1)     &
               ,dqi_agg     (k1)     &
               ,dni_agg     (k1)     &
               ,dns_sc      (k1)     &
               ,dqi_rime_ic (k1)     &
               ,dqs_rime_sc (k1)     &
               ,dqg_rime_gc (k1)     &
               ,dnc_rime_ic (k1)     &
               ,dnc_rime_sc (k1)     &
               ,dnc_rime_gc (k1)     &
               ,dqg_rime_gr (k1)     &
               ,dnr_rime_gr (k1)     &
               ,dqi_col_rig (k1)     &
               ,dqr_col_rig (k1)     &
               ,dni_col_rig (k1)     &
               ,dnr_col_rig (k1)     &
               ,dqs_col_rsg (k1)     &
               ,dqr_col_rsg (k1)     &
               ,dns_col_rsg (k1)     &
               ,dnr_col_rsg (k1)     &               
               ,dqs_col_si  (k1)     &
               ,dni_col_si  (k1)     &
               ,dqs_col_gsg (k1)     &
               ,dns_col_gsg (k1)     &
               ,dqi_cv_ig   (k1)     &
               ,dni_cv_ig   (k1)     &
               ,dqs_cv_sg   (k1)     &
               ,dns_cv_sg   (k1)     &
               ,dqc_sed     (k1)     &
               ,dqi_sed     (k1)     &
               ,dqr_sed     (k1)     &
               ,dqs_sed     (k1)     &
               ,dqg_sed     (k1)     &
               ,dnc_sed     (k1)     &
               ,dni_sed     (k1)     &
               ,dnr_sed     (k1)     &
               ,dns_sed     (k1)     &
               ,dng_sed     (k1)     &   
               ,dqc_tot     (k1)     & 
               ,dqi_tot     (k1)     & 
               ,dqr_tot     (k1)     & 
               ,dqs_tot     (k1)     & 
               ,dqg_tot     (k1)     & 
               ,dnc_tot     (k1)     & 
               ,dni_tot     (k1)     & 
               ,dnr_tot     (k1)     & 
               ,dns_tot     (k1)     & 
               ,dng_tot     (k1)     &                  
               ,dn_ccn_tot  (k1)     &
               ,dthl_mphys  (k1)     &
               ,dqi_eme_ic  (k1)     &
               ,dqs_eme_sc  (k1)     &
               ,dqg_eme_gc  (k1)     &
               ,dni_eme_ic  (k1)     &
               ,dns_eme_sc  (k1)     &
               ,dng_eme_gc  (k1)     &
               ,dqi_eme_ri  (k1)     &
               ,dqs_eme_rs  (k1)     &
               ,dqg_eme_gr  (k1)     &
               ,dni_eme_ri  (k1)     &
               ,dns_eme_rs  (k1)     &
               ,dng_eme_gr  (k1)     &
               ,dqi_me      (k1)     &
               ,dqs_me      (k1)     &
               ,dqg_me      (k1)     & 
               ,dni_me      (k1)     &
               ,dns_me      (k1)     &
               ,dng_me      (k1)     &
               ,dqi_ev     (k1)     &
               ,dqs_ev     (k1)     &
               ,dqg_ev     (k1)     & 
               ,dni_ev     (k1)     &
               ,dns_ev     (k1)     &
               ,dng_ev     (k1)     &
               ,dni_mul     (k1)     &
               ,dqi_mul     (k1)     &
               ,dqc_totf     (k1)     & 
               ,dqi_totf     (k1)     & 
               ,dqr_totf     (k1)     & 
               ,dqs_totf     (k1)     & 
               ,dqg_totf     (k1)     & 
               ,dnc_totf     (k1)     & 
               ,dni_totf     (k1)     & 
               ,dnr_totf     (k1)     & 
               ,dns_totf     (k1)     & 
               ,dng_totf     (k1)     &              
            )   
    ! #sb3 END
    

    
    Npmn    = 0.0
    qlpmn    = 0.0
    qtpmn    = 0.0
    precmn    = 0.0
    preccountmn  = 0.0
    prec_prcmn  = 0.0
    cloudcountmn  = 0.0
    raincountmn  = 0.0
    Nrrainmn  = 0.0
    qrmn    = 0.0
    Dvrmn    = 0.0
    ! #sb3 START
     cl_countmn  = 0.0
     ci_countmn  = 0.0
     hr_countmn  = 0.0
     hs_countmn  = 0.0
     hg_countmn  = 0.0
     prec_r_mn   = 0.0
     prec_i_mn   = 0.0
     prec_s_mn   = 0.0
     prec_g_mn   = 0.0   
     ! and process variables
     dnc_nuc      = 0.0 
     dqc_nuc      = 0.0
     dni_inuc     = 0.0
     dqi_inuc     = 0.0
     dqi_dep      = 0.0
     dqi_sub      = 0.0
     dqs_dep      = 0.0
     dqs_sub      = 0.0
     dqg_dep      = 0.0
     dqg_sub      = 0.0
     dnc_hom      = 0.0
     dqc_hom      = 0.0
     dnc_het      = 0.0
     dqc_het      = 0.0
     dnr_het      = 0.0
     dqr_het      = 0.0
     dnc_sadj     = 0.0
     dqc_sadjpl   = 0.0
     dqc_sadjneg  = 0.0
     dnc_sc       = 0.0 
     dnc_au       = 0.0
     dnr_au       = 0.0
     dqr_au       = 0.0
     dqr_ac       = 0.0
     dnc_ac       = 0.0
     dnr_sc       = 0.0
     dnr_br       = 0.0
     dqr_ev     = 0.0
     dnr_ev     = 0.0
     dqi_agg      = 0.0
     dni_agg      = 0.0
     dns_sc       = 0.0
     dqi_rime_ic  = 0.0
     dqs_rime_sc  = 0.0
     dqg_rime_gc  = 0.0
     dnc_rime_ic  = 0.0
     dnc_rime_sc  = 0.0
     dnc_rime_gc  = 0.0
     !dqs_rime_sr  = 0.0
     dqg_rime_gr  = 0.0
     ! dnr_rime_sr  = 0.0
     dnr_rime_gr  = 0.0
     dqi_col_rig = 0.0
     dqr_col_rig = 0.0
     dni_col_rig = 0.0
     dnr_col_rig = 0.0
     dqs_col_rsg = 0.0
     dqr_col_rsg = 0.0
     dns_col_rsg = 0.0
     dnr_col_rsg = 0.0    
     dqs_col_si   = 0.0
     dni_col_si   = 0.0
     dqs_col_gsg  = 0.0
     dns_col_gsg  = 0.0
     dqi_cv_ig    = 0.0
     dni_cv_ig    = 0.0
     dqs_cv_sg    = 0.0
     dns_cv_sg    = 0.0
     dqc_sed      = 0.0
     dqi_sed      = 0.0
     dqr_sed      = 0.0
     dqs_sed      = 0.0
     dqg_sed      = 0.0
     dnc_sed      = 0.0
     dni_sed      = 0.0
     dnr_sed      = 0.0
     dns_sed      = 0.0
     dng_sed      = 0.0 
     dqc_tot      = 0.0
     dqi_tot      = 0.0
     dqr_tot      = 0.0
     dqs_tot      = 0.0
     dqg_tot      = 0.0
     dnc_tot      = 0.0
     dni_tot      = 0.0
     dnr_tot      = 0.0
     dns_tot      = 0.0
     dng_tot      = 0.0  
     dn_ccn_tot   = 0.0
     dthl_mphys   = 0.0 
     dqi_eme_ic   = 0.0
     dqs_eme_sc   = 0.0
     dqg_eme_gc   = 0.0
     dni_eme_ic   = 0.0
     dns_eme_sc   = 0.0
     dng_eme_gc   = 0.0
     dqi_eme_ri   = 0.0
     dqs_eme_rs   = 0.0
     dqg_eme_gr   = 0.0
     dni_eme_ri   = 0.0
     dns_eme_rs   = 0.0
     dng_eme_gr   = 0.0
     dqi_me       = 0.0
     dqs_me       = 0.0
     dqg_me       = 0.0 
     dni_me       = 0.0
     dns_me       = 0.0
     dng_me       = 0.0
     dqi_ev      = 0.0
     dqs_ev      = 0.0
     dqg_ev      = 0.0 
     dni_ev      = 0.0
     dns_ev      = 0.0
     dng_ev      = 0.0
     dni_mul      = 0.0
     dqi_mul      = 0.0  
     dqc_totf      = 0.0
     dqi_totf      = 0.0
     dqr_totf      = 0.0
     dqs_totf      = 0.0
     dqg_totf      = 0.0
     dnc_totf      = 0.0
     dni_totf      = 0.0
     dnr_totf      = 0.0
     dns_totf      = 0.0
     dng_totf      = 0.0  
     
    ! #sb3 END    

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
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = itimeav/idtav
      if (myid==0) then  ! #t if (myid==10) then  ! if (myid==0) then  ! #t
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
        ! end of original profiles
        call define_nc( ncid_prof, nvar , ncname)
        !    
        ! and these shall go in a new profile file 
        fname_mphys(15:17) = cexpnr
        ! nvar_mphys= ! --> do
        allocate(ncmname(nvar_mphys,4))
        ! and now for bulk ice microphysics
        ! - cloud formation
        call nctiminfo(tncmname(1,:))
        ! call ncinfo(tncname(1,:),'time','Time','s','time
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
        call ncinfo(ncmname(nbaspout+18,:),'dq_c_ev', 'saturation adjustment: cloud water evaporation tendency',    'kg/kg/s','tt')
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
        call ncinfo(ncmname(nbaspout+26,:),'dq_r_ev', 'evaporation: rain water content tendency',            'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+27,:),'dn_r_ev', 'evaporation: raindrop number tendency',                '#/kg/s','tt')
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
        !call ncinfo(ncmname(nbaspout+37,:),'dq_s_rime_sr', 'snow riming  by rain water tendency',              'kg/kg/s','tt') 
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
        call ncinfo(ncmname(nbaspout+49,:),'dq_s_col_gsg',  'collection by graupel: snow content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+50,:),'dn_s_col_gsg',  'collection by graupel: snow number tendency',  '#/kg/s','tt')        
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
        call ncinfo(ncmname(nbaspout+63,:),'dq_i_me',  'partial melting of ice: ice content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+64,:),'dn_i_me',  'partial melting of ice: ice number tendency',  '#/kg/s','tt') 
        call ncinfo(ncmname(nbaspout+65,:),'dq_s_me',  'partial melting of snow: snow content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+66,:),'dn_s_me',  'partial melting of snow: snow number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+67,:),'dq_g_me',  'partial melting of graupel: graupel content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+68,:),'dn_g_me',  'partial melting of graupel: graupel number tendency',  '#/kg/s','tt')
        ! -- evaporation
        call ncinfo(ncmname(nbaspout+69,:),'dq_i_ev',  'evaporation of ice: ice content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+70,:),'dn_i_ev',  'evaporation of ice: ice number tendency',  '#/kg/s','tt') 
        call ncinfo(ncmname(nbaspout+71,:),'dq_s_ev',  'evaporation of snow: snow content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+72,:),'dn_s_ev',  'evaporation of snow: snow number tendency',  '#/kg/s','tt')
        call ncinfo(ncmname(nbaspout+73,:),'dq_g_ev',  'evaporation of graupel: graupel content tendency',  'kg/kg/s','tt')
        call ncinfo(ncmname(nbaspout+74,:),'dn_g_ev',  'evaporation of graupel: graupel number tendency',  '#/kg/s','tt')        
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
        call ncinfo(ncmname(nmphyout+6,:),'dth_ev',  'potential temperature tendency due to evaporation',     'K/kg/s','tt') 
        call ncinfo(ncmname(nmphyout+7,:),'dth_dep',  'potential temperature tendency due to depostion',     'K/kg/s','tt')
        call ncinfo(ncmname(nmphyout+8,:),'dth_sub',  'potential temperature tendency due to sublimation',     'K/kg/s','tt')        
        ! -- cloud fraction
        call ncinfo(ncmname(nmphyout+9,:),'cfrac_l',   'liquid cloud fraction','-','tt')
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
    use modmicrodata3, only  : l_hdump 
    implicit none
    
!     write(6,*) "bulkmicrostat3"     ! #debug
!     write(6,*) 'timee= ', timee ! #debug
!     write(6,*) 'tnext= ', tnext ! #debug
!     write(6,*) 'tnextwrite= ', tnextwrite ! #debug
!     write(6,*) 'dt_lim =', dt_lim ! #debug
!     write(6,*) 'rk3step =', rk3step ! #debug
!     write(6,*) 'lmicrostat =' , lmicrostat ! #debug
!     write(6,*) 'idtav = ', idtav ! #debug
!     write(6,*) 'itimeav = ', itimeav ! #debug
    
    if (.not. lmicrostat)  return
    if (rk3step /= 3)  return
    if (timee == 0)    return
    ! write(6,*) "going to following blocks" ! #debug
    if (timee < tnext .and. timee < tnextwrite) then
      ! write(6,*) "entering dt_lim block"
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if
    if (timee >= tnext) then
      ! write(6,*) "entering dobulkmicrostat3 block" ! #debug
      tnext = tnext + idtav
      ! write(6,*) "calling dobulkmicrostat3" ! #debug
      call dobulkmicrostat3   ! #sb3
      ! write(6,*) "leaving dobulkmicrostat3" ! #debug
    end if
    if (timee >= tnextwrite) then
      ! write(6,*) "entering writebulkmicrostat3 block" ! #debug
      tnextwrite = tnextwrite + itimeav
      call writebulkmicrostat3   ! #sb3
    end if
    
    ! call hydrodump
    if (l_hdump) then
      call hydrodump
    endif
    
  end subroutine bulkmicrostat3

!------------------------------------------------------------------------------!
!> Performs the calculations for rainrate etc.
  subroutine dobulkmicrostat3
    use modmpi,    only  : my_real, mpi_sum, comm3d, mpierr
    use modglobal,     only : i1, j1, k1, ijtot
    use modmicrodata,  only : qr,  Dvr, Nr, epscloud, epsqr       & ! precep
                             ,epsprec,imicro                      &
                             ,imicro_bulk,imicro_bulk3                               
    use modmicrodata3, only : eps_hprec                           &
                             ,in_hr,in_cl,in_ci,in_hs,in_hg,in_cc &
                             ,iq_hr,iq_cl,iq_ci,iq_hs,iq_hg       &
                             ,precep_hr,precep_ci                 &
                             ,precep_hs,precep_hg                 &
                             ,precep_l,precep_i                   &
                             ,n_cl,n_ci,n_hr,n_hs,n_hg            &
                             ,q_cl,q_ci,q_hr,q_hs,q_hg                           
    use modfields,  only  : ql0, rhof
    implicit none

    integer      :: k
 
    ! write(6,*) "in dobulkmicrostat3" ! #debug
    precav      = 0.0
    preccountav    = 0.0
    prec_prcav    = 0.0
    cloudcountav    = 0.0
    raincountav    = 0.0
    Nrrainav    = 0.0
    qrav      = 0.0
    Dvrav      = 0.0
    ! #sb3 START
     cl_countavl   = 0.0
     ci_countavl   = 0.0
     hr_countavl   = 0.0
     hs_countavl   = 0.0
     hg_countavl   = 0.0
     prec_r_avl    = 0.0
     prec_i_avl    = 0.0
     prec_s_avl    = 0.0
     prec_g_avl    = 0.0   
     ! and for large-scale variables
     cl_countav    = 0.0
     ci_countav    = 0.0
     hr_countav    = 0.0
     hs_countav    = 0.0
     hg_countav    = 0.0
     prec_r_av     = 0.0
     prec_i_av     = 0.0
     prec_s_av     = 0.0
     prec_g_av     = 0.0
    ! #sb3 END  

    ! calculating the values
    ! write(6,*) "calculating the values of old statistics " ! #debug
    do k = 1,k1
      cloudcountavl(k)  = count(ql0     (2:i1,2:j1,k) > epscloud)
      raincountavl (k)  = count(q_hr    (2:i1,2:j1,k) > epsqr)                          ! #sb3
      ! o: = count(qr      (2:i1,2:j1,k) > epsqr) #sb3
      preccountavl (k)  = count(precep_l(2:i1,2:j1,k) > epsprec)                        ! #sb3
      !o: preccountavl (k)  = count(precep (2:i1,2:j1,k) > epsprec)
      prec_prcavl  (k)  = sum  (precep_l(2:i1,2:j1,k), precep_l(2:i1,2:j1,k) > epsprec) !#sb3
      !o: prec_prcavl  (k)  = sum  (precep  (2:i1,2:j1,k)  , precep(2:i1,2:j1,k) > epsprec)
      Nrrainavl    (k)  = rhof(k)*sum  (n_hr      (2:i1,2:j1,k))
      !o: = sum  (Nr      (2:i1,2:j1,k)) #sb3
      precavl      (k)  = sum  (precep_l (2:i1,2:j1,k))                                 ! #sb3
      ! precavl      (k)  = sum  (precep  (2:i1,2:j1,k))
      qravl        (k)  = sum  (q_hr  (2:i1,2:j1,k))                                    ! #sb3
      !o: = sum  (qr  (2:i1,2:j1,k))
      ! if (imicro==imicro_bulk) then
      Dvravl     (k)  = sum  (Dvr  (2:i1,2:j1,k)  , q_hr  (2:i1,2:j1,k) > epsqr)
        !o: Dvravl     (k)  = sum  (Dvr  (2:i1,2:j1,k)  , qr  (2:i1,2:j1,k) > epsqr) #sb3
      ! end if
    end do
    
    ! #sb3 START
    ! additional statistics
    ! write(6,*) "calculating the values of additional statistics " ! #debug
    do k = 1,k1
      ! cloud column counts 
      ! write(6,*) " q_cl" ! #debug
      cl_countavl(k)  = count( q_cl(2:i1,2:j1,k) > epscloud)
      ! write(6,*) " q_ci" ! #debug
      ci_countavl(k)  = count( q_ci(2:i1,2:j1,k) > epscloud)
      ! write(6,*) " q_hr" ! #debug
      hr_countavl(k)  = count( q_hr(2:i1,2:j1,k) > eps_hprec)
      ! write(6,*) " q_hs" ! #debug
      hs_countavl(k)  = count( q_hs(2:i1,2:j1,k) > eps_hprec)
      ! write(6,*) " q_hg" ! #debug
      hg_countavl(k)  = count( q_hg(2:i1,2:j1,k) > eps_hprec)
      ! -- mixed clouds columns
      ! cmx_countavl(k) = count((q_cl(2:i1,2:j1,k) > epscloud).and.(q_ci(2:i1,2:j1,k) > epscloud))
      ! -- miced mphys columns 
      ! mx_countavl(k) = count(((q_cl(2:i1,2:j1,k)+                 & 
      !    q_hr(2:i1,2:j1,k))> eps_hprec).and.((q_ci(2:i1,2:j1,k) + &
      !     q_hs(2:i1,2:j1,k)+q_hg(2:i1,2:j1,k) )> eps_hprec))
      ! precipitation
      ! write(6,*) " precep_hr" ! #debug
      ! write(6,*) " prec_r_avl(k) = ", prec_r_avl (k)! #debug
      !write(6,*) " precep_hr  (2,2,k) = ", precep_hr (2,2,k)! #debug
      prec_r_avl (k)  = sum  (precep_hr  (2:i1,2:j1,k))
      ! write(6,*) " precep_ci" ! #debug
      prec_i_avl (k)  = sum  (precep_ci  (2:i1,2:j1,k))
      ! write(6,*) " precep_hs" ! #debug
      prec_s_avl (k)  = sum  (precep_hs  (2:i1,2:j1,k)) 
      ! write(6,*) " precep_hg" ! #debug
      prec_g_avl (k)  = sum  (precep_hg  (2:i1,2:j1,k))       
    end do
    ! #sb3 END

    ! write(6,*) "transmitting  " ! #debug
    ! transmitting to other processors
    call MPI_ALLREDUCE(cloudcountavl,cloudcountav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(raincountavl ,raincountav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(preccountavl  ,preccountav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(prec_prcavl  ,prec_prcav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(Dvravl  ,Dvrav    ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(Nrrainavl  ,Nrrainav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(precavl  ,precav    ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qravl  ,qrav    ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    ! transmitting to other processors #sb3 START
    call MPI_ALLREDUCE(cl_countavl ,cl_countav ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(ci_countavl ,ci_countav ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(hr_countavl ,hr_countav ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(hs_countavl ,hs_countav ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)  
    call MPI_ALLREDUCE(hg_countavl ,hg_countav ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(prec_r_avl  ,prec_r_av  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(prec_i_avl  ,prec_i_av  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(prec_s_avl  ,prec_s_av  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(prec_g_avl  ,prec_g_av  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    ! #sb3 END
    
    ! and normalising over all processors
    cloudcountmn  = cloudcountmn  +  cloudcountav  /ijtot
    raincountmn  = raincountmn  +  raincountav  /ijtot
    preccountmn  = preccountmn  +  preccountav  /ijtot
    prec_prcmn  = prec_prcmn  +  prec_prcav  /ijtot
    Dvrmn    = Dvrmn    +  Dvrav    /ijtot
    Nrrainmn  = Nrrainmn  +  Nrrainav  /ijtot
    precmn    = precmn  +  precav    /ijtot
    qrmn    = qrmn    +  qrav    /ijtot
    ! #sb3 START
    ! extra stats normalising over all processors  
    cl_countmn = cl_countmn + cl_countav /ijtot
    ci_countmn = ci_countmn + ci_countav /ijtot
    hr_countmn = hr_countmn + hr_countav /ijtot
    hs_countmn = hs_countmn + hs_countav /ijtot
    hg_countmn = hg_countmn + hg_countav /ijtot
    ! and precipitation 
    prec_r_mn  = prec_r_mn  + prec_r_av /ijtot
    prec_i_mn  = prec_i_mn  + prec_i_av /ijtot
    prec_s_mn  = prec_s_mn  + prec_s_av /ijtot
    prec_g_mn  = prec_g_mn  + prec_g_av /ijtot
    ! #sb3 END

  end subroutine dobulkmicrostat3

!------------------------------------------------------------------------------!
!> Performs the calculations for the tendencies etc.
  subroutine bulkmicrotend3
    use modmpi,    only  : slabsum
    use modglobal,    only  : rk3step, timee, dt_lim, k1, ih, i1, jh, j1, ijtot
    use modfields,    only  : qtp, rhof, svp    ! # test svp added
    use modmicrodata,  only : qtpmcr,thlpmcr
    use modmicrodata3, only : q_hrp, n_hrp                                  & !o: qrp, Nrp
                ,q_clp,q_cip,q_hsp,q_hgp                                    &
                ,n_clp,n_cip,n_hsp,n_hgp                                    &
                ,n_ccp                                                      &
                ,x_cnuc,x_inuc                                              &
                ,in_cl, in_ci, in_hr, in_hs, in_hg, in_cc                   &
                ,iq_cl, iq_ci, iq_hr, iq_hs, iq_hg                          &
                ,dn_cl_nu ,dn_ci_inu, dn_cl_au                              & 
                ,dn_hr_au , dn_cl_ac, dq_hr_au, dq_hr_ac                    &
                ,dn_hr_sc, dn_hr_br                                         & 
                ,dq_hr_ev, dn_hr_ev                                         &
                ,dq_ci_dep,dq_hs_dep  ,dq_hg_dep                            &
                ,dq_ci_rime ,dn_cl_rime_ci                                  & 
                ,dq_hs_rime ,dn_cl_rime_hs                                  &
                ,dq_hg_rime ,dn_cl_rime_hg                                  &
                ,dq_hghr_rime,dn_hr_rime_hg                                 &
                ,dq_hr_col_ri,dq_ci_col_ri,dn_ci_col_ri,dn_hr_col_ri        & 
                ,dq_hr_col_rs,dq_hs_col_rs,dn_hs_col_rs,dn_hr_col_rs        & 
                ,dq_cl_het,dn_cl_het,dq_hr_het,dn_hr_het                    &
                ,dq_cl_hom,dn_cl_hom                                        &
                ,dq_ci_col_iis,dn_ci_col_iis,dn_hs_col_sss                  &
                ,dq_hsci_col ,dn_ci_col_hs                                  & 
                ,dq_hghs_col, dn_hs_col_hg                                  &
                ,dq_ci_cv,dn_ci_cv,dq_hs_cv,dn_hs_cv                        &
                ,dn_cl_sc                                                   &
                ,dn_ci_mul, dq_ci_mul                                       &
                ,dq_ci_me,dq_hs_me, dq_hg_me                                & 
                ,dn_ci_me,dn_hs_me, dn_hg_me                                &
                ,dq_ci_ev,dq_hs_ev, dq_hg_ev                             & 
                ,dn_ci_ev,dn_hs_ev, dn_hg_ev                             &
                ,dn_ci_eme_ic,dq_ci_eme_ic,dn_ci_eme_ri,dq_ci_eme_ri        &  
                ,dn_hs_eme_sc,dq_hs_eme_sc,dn_hs_eme_rs,dq_hs_eme_rs        &
                ,dn_hg_eme_gc,dq_hg_eme_gc,dn_hg_eme_gr,dq_hg_eme_gr        &
                ,dn_cl_se,dq_cl_se,dn_ci_se,dq_ci_se ,dn_hr_se              &
                ,dq_hr_se,dn_hs_se,dq_hs_se,dn_hg_se,dq_hg_se               & 
                ,dq_cl_sa,dn_cl_sa,ret_cc
    
    implicit none

    real, dimension(:), allocatable  :: avfield
    real, dimension(:,:,:), allocatable  :: depos, subli,sadjpl,sadjneg, termtest    ! #sb3 # 
    integer        :: ifield = 0
    integer        :: i,j,k      ! loop variables 

    if (.not. lmicrostat)  return
    if (rk3step /= 3)  return
    if (timee == 0)    return
    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if
!    tnext = tnext+dtav

    allocate( avfield(k1))
    ! #sb3 START - allocate temporary variables
    allocate( depos(2-ih:i1+ih,2-jh:j1+jh,k1)                 &
             ,subli(2-ih:i1+ih,2-jh:j1+jh,k1)                 &
             ,sadjpl(2-ih:i1+ih,2-jh:j1+jh,k1)                &
             ,sadjneg(2-ih:i1+ih,2-jh:j1+jh,k1)               &
             ,termtest(2-ih:i1+ih,2-jh:j1+jh,k1)              &  !#
            )
    ! and initial values 
     depos = 0.0
     subli = 0.0
     sadjpl  = 0.0
     sadjneg = 0.0
     termtest = 0.0
     !
     dnc_nuc_av= 0.0 
     dqc_nuc_av= 0.0
     dni_inuc_av= 0.0
     dqi_inuc_av= 0.0
     dqi_dep_av= 0.0
     dqi_sub_av= 0.0
     dqs_dep_av= 0.0
     dqs_sub_av= 0.0
     dqg_dep_av= 0.0
     dqg_sub_av= 0.0
     dnc_hom_av= 0.0
     dqc_hom_av= 0.0
     dnc_het_av= 0.0
     dqc_het_av= 0.0
     dnr_het_av= 0.0
     dqr_het_av= 0.0
     dnc_sadj_av= 0.0
     dqc_sadjpl_av= 0.0
     dqc_sadjneg_av= 0.0
     dnc_sc_av  = 0.0
     dnc_au_av= 0.0
     dnr_au_av= 0.0
     dqr_au_av= 0.0
     dqr_ac_av= 0.0
     dnc_ac_av= 0.0
     dnr_sc_av= 0.0
     dnr_br_av= 0.0
     dqr_ev_av= 0.0
     dnr_ev_av= 0.0
     dqi_agg_av= 0.0
     dni_agg_av= 0.0
     dns_sc_av= 0.0
     dqi_rime_icav= 0.0
     dqs_rime_scav= 0.0
     dqg_rime_gcav= 0.0
     dnc_rime_icav= 0.0
     dnc_rime_scav= 0.0
     dnc_rime_gcav= 0.0
     ! dqs_rime_srav= 0.0
     dqg_rime_grav= 0.0
     ! dnr_rime_srav= 0.0
     dnr_rime_grav = 0.0
     dqi_col_rigav = 0.0
     dqr_col_rigav = 0.0
     dni_col_rigav = 0.0
     dnr_col_rigav = 0.0
     dqs_col_rsgav = 0.0
     dqr_col_rsgav = 0.0
     dns_col_rsgav = 0.0
     dnr_col_rsgav = 0.0     
     dqs_col_siav  = 0.0
     dni_col_siav  = 0.0
     dqs_col_gsgav = 0.0
     dns_col_gsgav = 0.0
     dqi_cv_igav   = 0.0
     dni_cv_igav   = 0.0
     dqs_cv_sgav   = 0.0
     dns_cv_sgav   = 0.0
     dn_ccn_av     = 0.0
     dqc_sedav     = 0.0
     dqi_sedav     = 0.0
     dqr_sedav     = 0.0
     dqs_sedav     = 0.0
     dqg_sedav     = 0.0
     dnc_sedav     = 0.0
     dni_sedav     = 0.0
     dnr_sedav     = 0.0
     dns_sedav     = 0.0
     dng_sedav     = 0.0 
     dqc_totav     = 0.0
     dqi_totav     = 0.0
     dqr_totav     = 0.0
     dqs_totav     = 0.0
     dqg_totav     = 0.0
     dnc_totav     = 0.0
     dni_totav     = 0.0
     dnr_totav     = 0.0
     dns_totav     = 0.0
     dng_totav     = 0.0 
     dqc_totfav     = 0.0
     dqi_totfav     = 0.0
     dqr_totfav     = 0.0
     dqs_totfav     = 0.0
     dqg_totfav     = 0.0
     dnc_totfav     = 0.0
     dni_totfav     = 0.0
     dnr_totfav     = 0.0
     dns_totfav     = 0.0
     dng_totfav     = 0.0      
     dthl_mphys_av = 0.0
     dqi_eme_icav  = 0.0
     dqs_eme_scav  = 0.0
     dqg_eme_gcav  = 0.0
     dni_eme_icav  = 0.0
     dns_eme_scav  = 0.0
     dng_eme_gcav  = 0.0
     dqi_eme_riav  = 0.0
     dqs_eme_rsav  = 0.0
     dqg_eme_grav  = 0.0
     dni_eme_riav  = 0.0
     dns_eme_rsav  = 0.0
     dng_eme_grav  = 0.0
     dqi_me_av      = 0.0
     dqs_me_av      = 0.0
     dqg_me_av      = 0.0 
     dni_me_av      = 0.0
     dns_me_av      = 0.0
     dng_me_av      = 0.0
     dqi_ev_av     = 0.0
     dqs_ev_av     = 0.0
     dqg_ev_av     = 0.0 
     dni_ev_av     = 0.0
     dns_ev_av     = 0.0
     dng_ev_av     = 0.0
     dni_mul_av     = 0.0
     dqi_mul_av     = 0.0  
    ! #sb3 END       
    !o: ifield    = mod(ifield, nrfields) + 1
    
    ifield = iauto ! autoconversion
    
    avfield    = 0.0
    call slabsum(avfield  ,1,k1,dn_hr_au  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) ! #sb3
    !o: call slabsum(avfield  ,1,k1,Nrp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    do k=1,k1
      Npav(k,ifield)  = rhof(k)*avfield(k)  ! since n_hrp is in #/kg, Npav in #/m3 
    end do
    !o: Npav(:,ifield)  = avfield - sum(Npav  (:,1:ifield-1),2)

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,dq_hr_au  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) ! #sb3
    !o: call slabsum(avfield  ,1,k1,qrp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qlpav(:,ifield) = avfield
    !o: qlpav(:,ifield) = avfield - sum(qlpav  (:,1:ifield-1),2)
    !
    ! avfield    = 0.0
    ! call slabsum(avfield  ,1,k1, dq_hr_au ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qtpav(:,ifield) = - avfield
    !o: qtpav(:,ifield) = avfield - sum(qtpav  (:,1:ifield-1),2)
    
    ! #sb3 START
    ! and repeat it for other processes
    
    ifield = iaccr ! accretion
    
    ! avfield    = 0.0
    ! call slabsum(avfield  ,1,k1,dn_hr_ac  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    do k=1,k1
      Npav(k,ifield)  = 0.0  ! since accretion does not change n_hr
    end do

    avfield    = 0.0
    ! call slabsum(avfield  ,1,k1, ac  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    ! qlpav(:,ifield) = - avfield
    !
    ! avfield    = 0.0
    ! call slabsum(avfield  ,1,k1, ac ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qtpav(:,ifield) = avfield 
    
    ifield = ised ! sedimentation
    
    avfield    = 0.0
    call slabsum(avfield  ,1,k1,dn_hr_se  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    do k=1,k1
      Npav(k,ifield)  = rhof(k)*avfield(k)  ! since n_hrp is in #/kg, Npav in #/m3 
    end do

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,dq_hr_se  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qlpav(:,ifield) = avfield

    avfield    = 0.0
    ! rain sedimentation does not change qt
    qtpav(:,ifield) = avfield     
    
    ifield = ievap ! evaporation
    
    avfield    = 0.0
    call slabsum(avfield  ,1,k1,dn_hr_ev ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    do k=1,k1
      Npav(k,ifield)  = rhof(k)*avfield(k)  ! since n_hrp is in #/kg, Npav in #/m3 
    end do

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,dq_hr_ev ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qlpav(:,ifield) = avfield

    ! since evaporation leads direcly to increase in qtp
    qtpav(:,ifield) = - avfield     
    ! #sb3 END
    
    ! #sb3 START
    ! and other tendencies
    ! avfield    = 0.0
     !-nucelation
     call slabsum(dnc_nuc_av,1,k1, dn_cl_nu  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
       dqc_nuc_av(k) = x_cnuc*dnc_nuc_av(k)
      end do
     call slabsum(dni_inuc_av,1,k1, dn_ci_inu  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
       dqi_inuc_av(k) =x_inuc*dni_inuc_av(k) 
      end do 
     !- deposition, sublimation
      !-  separating on depostion and sublimation
      do j=2,j1
      do i=2,i1
      do k=1,k1       
        depos(i,j,k) = max(0.0,dq_ci_dep(i,j,k))
        subli(i,j,k) = min(0.0,dq_ci_dep(i,j,k))
        
      enddo
      enddo
      enddo
      call slabsum(dqi_dep_av,1,k1,depos,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqi_sub_av,1,k1,subli,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     
      do j=2,j1
      do i=2,i1
      do k=1,k1       
        depos(i,j,k) = max(0.0,dq_hs_dep(i,j,k))
        subli(i,j,k) = min(0.0,dq_hs_dep(i,j,k))
      enddo
      enddo
      enddo
      call slabsum(dqs_dep_av,1,k1,depos,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqs_sub_av,1,k1,subli,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
     
      do j=2,j1
      do i=2,i1
      do k=1,k1       
        depos(i,j,k) = max(0.0,dq_hg_dep(i,j,k))
        subli(i,j,k) = min(0.0,dq_hg_dep(i,j,k))
      enddo
      enddo
      enddo
      call slabsum(dqg_dep_av,1,k1,depos,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
      call slabsum(dqg_sub_av,1,k1,subli,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    ! - freezing
     ! -- homogeneous freezing
      call slabsum(dnc_hom_av,1,k1,dn_cl_hom,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqc_hom_av,1,k1,dq_cl_hom,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! -- heterogeneous freezing 
      call slabsum(dnc_het_av,1,k1,dn_cl_het,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqc_het_av,1,k1,dq_cl_het,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)   
      call slabsum(dnr_het_av,1,k1,dn_hr_het,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqr_het_av,1,k1,dq_hr_het,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    ! - saturation adjustment
     ! -- separating onto positive and negative
      do j=2,j1
      do i=2,i1
      do k=1,k1       
        sadjpl(i,j,k)  = max(0.0,dq_cl_sa (i,j,k))
        sadjneg(i,j,k) = min(0.0,dq_cl_sa (i,j,k))
      enddo
      enddo
      enddo
      call slabsum(dqc_sadjpl_av,1,k1,sadjpl ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqc_sadjneg_av,1,k1,sadjneg,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)    
      !o call slabsum(dqc_sadj_av,1,k1,dq_cl_sa,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dnc_sadj_av,1,k1,dn_cl_sa,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    ! - warm microphysics interactions
     ! -- self-collection
      call slabsum(dnc_sc_av,1,k1,dn_cl_sc,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! -- autoconversion
      call slabsum(dqr_au_av,1,k1,dq_hr_au,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
      ! dqc_au_av = - dqr_au_av ! since change in cloud water is oposite
      call slabsum(dnc_au_av,1,k1,dn_cl_au,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
      call slabsum(dnr_au_av,1,k1,dn_hr_au,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
     ! -- accretion      
      call slabsum(dqr_ac_av,1,k1,dq_hr_ac,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
      ! dqc_ac_av=-dqr_ac_av
      call slabsum(dnc_ac_av,1,k1,dn_cl_ac,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)      
     ! -- rain selfcollection and breakup 
      call slabsum(dnr_sc_av,1,k1,dn_hr_sc,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
      call slabsum(dnr_br_av,1,k1,dn_hr_br,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! -- rain evaporation 
      call slabsum(dqr_ev_av,1,k1,dq_hr_ev,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
      call slabsum(dnr_ev_av,1,k1,dn_hr_ev,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! - cold aggregation and self-collection
     ! -- ice aggregation 
      call slabsum(dqi_agg_av,1,k1,dq_ci_col_iis,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
      call slabsum(dni_agg_av,1,k1,dn_ci_col_iis,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! -- snow self-collection
      call slabsum(dns_sc_av,1,k1,dn_hs_col_sss,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! - riming 
     ! -- riming by cloud water
      call slabsum(dqi_rime_icav,1,k1,dq_ci_rime,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqs_rime_scav,1,k1,dq_hs_rime,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqg_rime_gcav,1,k1,dq_hg_rime,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dnc_rime_icav,1,k1,dn_cl_rime_ci,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dnc_rime_scav,1,k1,dn_cl_rime_hs,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dnc_rime_gcav,1,k1,dn_cl_rime_hg,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! -- riming by rain water
      ! call slabsum(dqs_rime_srav,1,k1,dq_hshr_rime,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqg_rime_grav,1,k1,dq_hghr_rime,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      ! call slabsum(dnr_rime_srav,1,k1,dn_hr_rime_hs,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dnr_rime_grav,1,k1,dn_hr_rime_hg,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)   
     ! - collisions
     ! --  collision conversions
      call slabsum(dqi_col_rigav,1,k1, dq_ci_col_ri,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqr_col_rigav,1,k1, dq_hr_col_ri,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dni_col_rigav,1,k1, dn_ci_col_ri,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dnr_col_rigav,1,k1, dn_hr_col_ri,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqs_col_rsgav,1,k1, dq_hs_col_rs,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqr_col_rsgav,1,k1, dq_hr_col_rs,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dns_col_rsgav,1,k1, dn_hs_col_rs,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dnr_col_rsgav,1,k1, dn_hr_col_rs,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)      
     ! -- collision collection
      call slabsum(dqs_col_siav,1,k1,dq_hsci_col ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dni_col_siav,1,k1,dn_ci_col_hs ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqs_col_gsgav,1,k1,-dq_hghs_col ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dns_col_gsgav,1,k1,dn_hs_col_hg ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! -- enhanced melting
      call slabsum(dqi_eme_icav,1,k1,dq_ci_eme_ic,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dni_eme_icav,1,k1,dn_ci_eme_ic,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqi_eme_riav,1,k1,dq_ci_eme_ri,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dni_eme_riav,1,k1,dn_ci_eme_ri,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqs_eme_scav,1,k1,dq_hs_eme_sc,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dns_eme_scav,1,k1,dn_hs_eme_sc,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqs_eme_rsav,1,k1,dq_hs_eme_rs,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dns_eme_rsav,1,k1,dn_hs_eme_rs,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqg_eme_gcav,1,k1,dq_hg_eme_gc,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dng_eme_gcav,1,k1,dn_hg_eme_gc,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqg_eme_grav,1,k1,dq_hg_eme_gr,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dng_eme_grav,1,k1,dn_hg_eme_gr,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)      
     ! - melting and evaporation of ice
     ! -- partial melting
     call slabsum(dqi_me_av,1,k1,dq_ci_me,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dni_me_av,1,k1,dn_ci_me,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dqs_me_av,1,k1,dq_hs_me,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dns_me_av,1,k1,dn_hs_me,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dqg_me_av,1,k1,dq_hg_me,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dng_me_av,1,k1,dn_hg_me,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! -- evaporation
     call slabsum(dqi_ev_av,1,k1,dq_ci_ev,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dni_ev_av,1,k1,dn_ci_ev,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dqs_ev_av,1,k1,dq_hs_ev,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dns_ev_av,1,k1,dn_hs_ev,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dqg_ev_av,1,k1,dq_hg_ev,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dng_ev_av,1,k1,dn_hg_ev,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)    
     !
     ! - conversions and other processes
     ! -- partial conversion
      call slabsum(dqi_cv_igav,1,k1,dq_ci_cv,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dni_cv_igav,1,k1,dn_ci_cv,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqs_cv_sgav,1,k1,dq_hs_cv,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dns_cv_sgav,1,k1,dn_hs_cv,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! -- ice multiplication
     call slabsum(dqi_mul_av,1,k1,dq_ci_mul,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     call slabsum(dni_mul_av,1,k1,dn_ci_mul,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     !
     ! sedimentation
      call slabsum(dqc_sedav,1,k1,dq_cl_se,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqi_sedav,1,k1,dq_ci_se,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqr_sedav,1,k1,dq_hr_se,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqs_sedav,1,k1,dq_hs_se,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)      
      call slabsum(dqg_sedav,1,k1,dq_hg_se,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dnc_sedav,1,k1,dn_cl_se,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dni_sedav,1,k1,dn_ci_se,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dnr_sedav,1,k1,dn_hr_se,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dns_sedav,1,k1,dn_hs_se,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)      
      call slabsum(dng_sedav,1,k1,dn_hg_se,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! total tendencies
     ! - total tendencies in water species
      call slabsum(dnc_totav,1,k1,n_clp,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dni_totav,1,k1,n_cip,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dnr_totav,1,k1,n_hrp,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dns_totav,1,k1,n_hsp,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dng_totav,1,k1,n_hgp,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqc_totav,1,k1,q_clp,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqi_totav,1,k1,q_cip,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqr_totav,1,k1,q_hrp,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqs_totav,1,k1,q_hsp,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(dqg_totav,1,k1,q_hgp,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)     
      ! -total tendencies 
      !#
      do k=1,k1
      do j=2,j1
      do i=2,i1
        termtest(i,j,k) = svp(i,j,k,in_cl)
      enddo
      enddo
      enddo
      call slabsum(dnc_totfav,1,k1, termtest,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
      do j=2,j1
      do i=2,i1
        termtest(i,j,k) = svp(i,j,k,in_ci)
      enddo
      enddo     
      enddo
      call slabsum(dni_totfav,1,k1, termtest,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
      do j=2,j1
      do i=2,i1
        termtest(i,j,k) = svp(i,j,k,in_hr)
      enddo
      enddo     
      enddo      
      call slabsum(dnr_totfav,1,k1, termtest,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
      do j=2,j1
      do i=2,i1
        termtest(i,j,k) = svp(i,j,k,in_hs)
      enddo
      enddo     
      enddo      
      call slabsum(dns_totfav,1,k1,termtest ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
      do j=2,j1
      do i=2,i1
        termtest(i,j,k) = svp(i,j,k,in_hg)
      enddo
      enddo     
      enddo      
      call slabsum(dng_totfav,1,k1,termtest ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
      do j=2,j1
      do i=2,i1
        termtest(i,j,k) = svp(i,j,k,iq_cl)
      enddo
      enddo     
      enddo      
      call slabsum(dqc_totfav,1,k1,termtest ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
      do j=2,j1
      do i=2,i1
        termtest(i,j,k) = svp(i,j,k,iq_ci)
      enddo
      enddo     
      enddo      
      call slabsum(dqi_totfav,1,k1,termtest ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
      do j=2,j1
      do i=2,i1
        termtest(i,j,k) = svp(i,j,k,iq_hr)
      enddo
      enddo     
      enddo      
      call slabsum(dqr_totfav,1,k1,termtest ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
      do j=2,j1
      do i=2,i1
        termtest(i,j,k) = svp(i,j,k,iq_hs)
      enddo
      enddo     
      enddo      
      call slabsum(dqs_totfav,1,k1,termtest ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do k=1,k1
      do j=2,j1
      do i=2,i1
        termtest(i,j,k) = svp(i,j,k,iq_hg)
      enddo
      enddo     
      enddo      
      call slabsum(dqg_totfav,1,k1,termtest ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)      
      
     ! - tendency in  ccn      
     ! -- total tendency in ccn
      call slabsum(dn_ccn_av,1,k1,n_ccp,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! -- recovery of ccn 
     ! call slabsum(dn_ccn_retav,1,k1,n_ccp,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
     ! - total tendenty in thlpmcr
      call slabsum(dthl_mphys_av,1,k1,thlpmcr,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
      ! call slabsum(qtp_totav,1,k1,qtpmcr,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) 
    ! #sb3 END
    
    
    !o: if (ifield == nrfields) then
    !o:   Npmn    = Npmn  + Npav  /nsamples/ijtot
    !o:   qlpmn    = qlpmn  + qlpav /nsamples/ijtot
    !o:   qtpmn    = qtpmn  + qtpav /nsamples/ijtot
    !o:   Npav    = 0.0
    !o:   qlpav    = 0.0
    !o:   qtpav    = 0.0
    !end if
    ! #sb3 START
      Npmn    = Npmn  + Npav  /nsamples/ijtot
      qlpmn    = qlpmn  + qlpav /nsamples/ijtot
      qtpmn    = qtpmn  + qtpav /nsamples/ijtot
      Npav    = 0.0
      qlpav    = 0.0
      qtpav    = 0.0
    ! #sb3 END
    
    ! #sb3 START
    ! and processing other tendencies
     dnc_nuc    = dnc_nuc    + dnc_nuc_av    /nsamples/ijtot
     dqc_nuc    = dqc_nuc    + dqc_nuc_av    /nsamples/ijtot
     dni_inuc   = dni_inuc   + dni_inuc_av   /nsamples/ijtot
     dqi_inuc   = dqi_inuc   + dqi_inuc_av   /nsamples/ijtot  
     !
     dqi_dep    = dqi_dep    + dqi_dep_av    /nsamples/ijtot
     dqi_sub    = dqi_sub    + dqi_sub_av    /nsamples/ijtot
     dqs_dep    = dqs_dep    + dqs_dep_av    /nsamples/ijtot
     dqs_sub    = dqs_sub    + dqs_sub_av    /nsamples/ijtot
     dqg_dep    = dqs_dep    + dqs_dep_av    /nsamples/ijtot
     dqg_sub    = dqs_sub    + dqs_sub_av    /nsamples/ijtot    
     !
     dnc_hom    = dnc_hom    + dnc_hom_av    /nsamples/ijtot 
     dqc_hom    = dqc_hom    + dqc_hom_av    /nsamples/ijtot  
     dnc_het    = dnc_het    + dnc_het_av    /nsamples/ijtot 
     dqc_het    = dqc_het    + dqc_het_av    /nsamples/ijtot 
     dnr_het    = dnr_het    + dnr_het_av    /nsamples/ijtot
     dqr_het    = dqr_het    + dqr_het_av    /nsamples/ijtot
     !
     dnc_sadj   = dnc_sadj   + dnc_sadj_av   /nsamples/ijtot
     dqc_sadjpl = dqc_sadjpl + dqc_sadjpl_av /nsamples/ijtot
     dqc_sadjneg= dqc_sadjneg+ dqc_sadjneg_av/nsamples/ijtot
     ! 
     dnc_au     = dnc_au     + dnc_au_av     /nsamples/ijtot
     dnr_au     = dnr_au     + dnr_au_av     /nsamples/ijtot
     dqr_au     = dqr_au     + dqr_au_av     /nsamples/ijtot
     !
     dqr_ac     = dqr_ac     + dqr_ac_av     /nsamples/ijtot
     dnc_ac     = dnc_ac     + dnc_ac_av     /nsamples/ijtot
     !
     dnr_sc     = dnr_sc     + dnr_sc_av     /nsamples/ijtot
     dnr_br     = dnr_br     + dnr_br_av     /nsamples/ijtot
     !
     dqr_ev   = dqr_ev   + dqr_ev_av   /nsamples/ijtot
     dnr_ev   = dnr_ev   + dnr_ev_av   /nsamples/ijtot
     !
     dqi_agg    = dqi_agg    + dqi_agg_av    /nsamples/ijtot
     dni_agg    = dni_agg    + dni_agg_av    /nsamples/ijtot
     !
     dns_sc     = dns_sc     + dns_sc_av     /nsamples/ijtot
     dnc_sc     = dnc_sc     + dnc_sc_av     /nsamples/ijtot
     !
     dqi_rime_ic= dqi_rime_ic+ dqi_rime_icav /nsamples/ijtot
     dqs_rime_sc= dqs_rime_sc+ dqs_rime_scav /nsamples/ijtot
     dqg_rime_gc= dqg_rime_gc+ dqg_rime_gcav /nsamples/ijtot
     dnc_rime_ic= dnc_rime_ic+ dnc_rime_icav /nsamples/ijtot
     dnc_rime_sc= dnc_rime_sc+ dnc_rime_scav /nsamples/ijtot
     dnc_rime_gc= dnc_rime_gc+ dnc_rime_gcav /nsamples/ijtot
     !
     ! dqs_rime_sr= dqs_rime_sr+ dqs_rime_srav /nsamples/ijtot
     dqg_rime_gr= dqg_rime_gr+ dqg_rime_grav /nsamples/ijtot
     ! dnr_rime_sr= dnr_rime_sr+ dnr_rime_srav /nsamples/ijtot
     dnr_rime_gr= dnr_rime_gr+ dnr_rime_grav /nsamples/ijtot
     !
     dqi_col_rig = dqi_col_rig + dqi_col_rigav /nsamples/ijtot
     dqr_col_rig = dqr_col_rig + dqr_col_rigav /nsamples/ijtot
     dni_col_rig = dni_col_rig + dni_col_rigav /nsamples/ijtot
     dnr_col_rig = dnr_col_rig + dnr_col_rigav /nsamples/ijtot
     !
     dqs_col_rsg = dqs_col_rsg + dqs_col_rsgav /nsamples/ijtot
     dqr_col_rsg = dqr_col_rsg + dqr_col_rsgav /nsamples/ijtot
     dns_col_rsg = dns_col_rsg + dns_col_rsgav /nsamples/ijtot
     dnr_col_rsg = dnr_col_rsg + dnr_col_rsgav /nsamples/ijtot     
     !
     dqs_col_si  = dqs_col_si  + dqs_col_siav  /nsamples/ijtot
     dni_col_si  = dni_col_si  + dni_col_siav  /nsamples/ijtot
     dqs_col_gsg = dqs_col_gsg + dqs_col_gsgav /nsamples/ijtot
     dns_col_gsg = dns_col_gsg + dns_col_gsgav /nsamples/ijtot
     !
     dqi_eme_ic = dqi_eme_ic   + dqi_eme_icav  /nsamples/ijtot
     dni_eme_ic = dni_eme_ic   + dni_eme_icav  /nsamples/ijtot
     dqi_eme_ri = dqi_eme_ri   + dqi_eme_riav  /nsamples/ijtot
     dni_eme_ri = dni_eme_ri   + dni_eme_riav  /nsamples/ijtot
     dqs_eme_sc = dqs_eme_sc   + dqs_eme_scav  /nsamples/ijtot
     dns_eme_sc = dns_eme_sc   + dns_eme_scav  /nsamples/ijtot
     dqs_eme_rs = dqs_eme_rs   + dqs_eme_rsav  /nsamples/ijtot
     dns_eme_rs = dns_eme_rs   + dns_eme_rsav  /nsamples/ijtot
     dqg_eme_gc = dqg_eme_gc   + dqg_eme_gcav  /nsamples/ijtot
     dng_eme_gc = dng_eme_gc   + dng_eme_gcav  /nsamples/ijtot
     dqg_eme_gr = dqg_eme_gr   + dqg_eme_grav  /nsamples/ijtot
     dng_eme_gr = dng_eme_gr   + dng_eme_grav  /nsamples/ijtot
     !
     dqi_me     = dqi_me       + dqi_me_av     /nsamples/ijtot
     dni_me     = dni_me       + dni_me_av     /nsamples/ijtot
     dqs_me     = dqs_me       + dqs_me_av     /nsamples/ijtot
     dns_me     = dns_me       + dns_me_av     /nsamples/ijtot
     dqg_me     = dqg_me       + dqg_me_av     /nsamples/ijtot
     dng_me     = dng_me       + dng_me_av     /nsamples/ijtot
     !
     dqi_ev    = dqi_ev      + dqi_ev_av    /nsamples/ijtot
     dni_ev    = dni_ev      + dni_ev_av    /nsamples/ijtot
     dqs_ev    = dqs_ev      + dqs_ev_av    /nsamples/ijtot
     dns_ev    = dns_ev      + dns_ev_av    /nsamples/ijtot
     dqg_ev    = dqg_ev      + dqg_ev_av    /nsamples/ijtot
     dng_ev    = dng_ev      + dng_ev_av    /nsamples/ijtot     
     !
     dqi_cv_ig  = dqi_cv_ig    + dqi_cv_igav   /nsamples/ijtot
     dni_cv_ig  = dni_cv_ig    + dni_cv_igav   /nsamples/ijtot
     dqs_cv_sg  = dqs_cv_sg    + dqs_cv_sgav   /nsamples/ijtot
     dns_cv_sg  = dns_cv_sg    + dns_cv_sgav   /nsamples/ijtot
     !
     dqi_cv_ig  = dqi_cv_ig    + dqi_cv_igav   /nsamples/ijtot
     dni_cv_ig  = dni_cv_ig    + dni_cv_igav   /nsamples/ijtot
     dqs_cv_sg  = dqs_cv_sg    + dqs_cv_sgav   /nsamples/ijtot
     dns_cv_sg  = dns_cv_sg    + dns_cv_sgav   /nsamples/ijtot
     !
     dqi_mul    = dqi_mul      + dqi_mul_av    /nsamples/ijtot
     dni_mul    = dni_mul      + dni_mul_av    /nsamples/ijtot
     ! 
     ! sedimentation
     dqc_sed    = dqc_sed      + dqc_sedav     /nsamples/ijtot
     dqi_sed    = dqi_sed      + dqi_sedav     /nsamples/ijtot
     dqr_sed    = dqr_sed      + dqr_sedav     /nsamples/ijtot
     dqs_sed    = dqs_sed      + dqs_sedav     /nsamples/ijtot
     dqg_sed    = dqg_sed      + dqg_sedav     /nsamples/ijtot
     dnc_sed    = dnc_sed      + dnc_sedav     /nsamples/ijtot
     dni_sed    = dni_sed      + dni_sedav     /nsamples/ijtot
     dnr_sed    = dnr_sed      + dnr_sedav     /nsamples/ijtot
     dns_sed    = dns_sed      + dns_sedav     /nsamples/ijtot
     dng_sed    = dng_sed      + dng_sedav     /nsamples/ijtot
     
     ! and total tendencies
     dqc_tot    = dqc_tot + dqc_totav        /nsamples/ijtot
     dqr_tot    = dqr_tot + dqr_totav        /nsamples/ijtot
     dqi_tot    = dqi_tot + dqi_totav        /nsamples/ijtot
     dqs_tot    = dqs_tot + dqs_totav        /nsamples/ijtot
     dqg_tot    = dqg_tot + dqg_totav        /nsamples/ijtot
     ! 
     dnc_tot    = dnc_tot + dnc_totav        /nsamples/ijtot
     dnr_tot    = dnr_tot + dnr_totav        /nsamples/ijtot
     dni_tot    = dni_tot + dni_totav        /nsamples/ijtot
     dns_tot    = dns_tot + dns_totav        /nsamples/ijtot
     dng_tot    = dng_tot + dng_totav        /nsamples/ijtot
     
      ! and total tendencies
     dqc_totf    = dqc_totf + dqc_totfav        /nsamples/ijtot
     dqr_totf    = dqr_totf + dqr_totfav        /nsamples/ijtot
     dqi_totf    = dqi_totf + dqi_totfav        /nsamples/ijtot
     dqs_totf    = dqs_totf + dqs_totfav        /nsamples/ijtot
     dqg_totf    = dqg_totf + dqg_totfav        /nsamples/ijtot
     ! 
     dnc_totf    = dnc_totf + dnc_totfav        /nsamples/ijtot
     dnr_totf    = dnr_totf + dnr_totfav        /nsamples/ijtot
     dni_totf    = dni_totf + dni_totfav        /nsamples/ijtot
     dns_totf    = dns_totf + dns_totfav        /nsamples/ijtot
     dng_totf    = dng_totf + dng_totfav        /nsamples/ijtot    
     !
     dn_ccn_tot = dn_ccn_tot + dn_ccn_av     /nsamples/ijtot
     !
     ! thl and qtp
     dthl_mphys = dthl_mphys + dthl_mphys_av /nsamples/ijtot
     ! qtp        = qtp_totav
      
     ! and the aggregated statistics
     !- done in: writebulkmicrostat3
     ! #sb3 START 

    deallocate(avfield)
    deallocate(depos,subli,sadjpl,sadjneg)
    deallocate(termtest)

  end subroutine bulkmicrotend3

!------------------------------------------------------------------------------!
!> Write the stats to file
  subroutine writebulkmicrostat3
    use modmpi,    only  : myid
    use modglobal,    only  : rtimee, ifoutput, cexpnr, k1,kmax  &
                             ,rlv, zf                            &
                             ,cp   ! #sb3
    use modfields,    only  : presf,rhof                         &
                             ,exnf ! #sb3
    use modstat_nc, only: lnetcdf, writestat_nc
    use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec
    use modmicrodata3, only: rlme, rlvi   ! #sb3

    implicit none
     real,dimension(k1,nvar) :: vars
     real,allocatable,dimension(:,:) :: vars_mphys

     integer    :: nsecs, nhrs, nminut
     integer    :: k

      ! #sb3 START
      allocate( vars_mphys (k1,nvar_mphys) )   ! allocate ouptut field
      allocate(  cfrac_l     (k1)                             &
                ,cfrac_i     (k1)                             &
                ,cfrac_tot   (k1)                             &
                ,dqc_rime    (k1)                             &
                ,dqr_rime    (k1)                             &
                ,dnc_rime    (k1)                             &
                ,dnr_rime    (k1)                             &        
                ,dth_mphys   (k1)                             &
                ,dth_freeze  (k1)                             &
                ,dth_melt    (k1)                             &
                ,dth_ev     (k1)                             &
                ,dth_cond    (k1)                             &
                ,dth_dep     (k1)                             &
                ,dth_sub     (k1)                             & 
              )
              
!                 ,dqi_rime  (k1)                                  &
!                 ,dqs_rime  (k1)                                  &
!                 ,dqg_rime  (k1)                                  &  
!                 ,dqi_hom   (k1)                                  &
!                 ,dqi_het   (k1)                                  &
!                 ,dqg_het   (k1)                                  &
!                 ,dni_hom   (k1)                                  &
!                 ,dni_het   (k1)                                  &
!                 ,dng_het   (k1)                                  &
!                 ,dqi_col_si(k1)                                  &
!                 ,dqi_col   (k1)                                  &
!                 ,dqs_col   (k1)                                  &
!                 ,dqg_col   (k1)                                  &
!                 ,dni_col   (k1)                                  &
!                 ,dns_col   (k1)                                  &
!                 ,dng_col   (k1)                                  &
!                 ,dqi_cv    (k1)                                  &
!                 ,dqs_cv    (k1)                                  &
!                 ,dqg_cv    (k1)                                  &
!                 ,dni_cv    (k1)                                  &
!                 ,dns_cv    (k1)                                  &
!                 ,dng_cv    (k1)                                  &
!                )
      cfrac_l      = 0.0
      cfrac_i      = 0.0
      cfrac_tot    = 0.0 
      dth_mphys    = 0.0
      dth_freeze   = 0.0
      dth_melt     = 0.0
      dth_ev      = 0.0
      dth_cond     = 0.0
      dth_dep      = 0.0
      dth_sub      = 0.0
      ! #sb3 END
    
    nsecs    = nint(rtimee)
    nhrs    = int (nsecs/3600)
    nminut    = int (nsecs/60)-nhrs*60
    nsecs    = mod (nsecs,60)

    cloudcountmn    = cloudcountmn  /nsamples
                raincountmn     = raincountmn   /nsamples
                preccountmn     = preccountmn   /nsamples
                prec_prcmn      = prec_prcmn    /nsamples
                Dvrmn           = Dvrmn         /nsamples
                Nrrainmn        = Nrrainmn      /nsamples
                precmn          = precmn        /nsamples
                qrmn            = qrmn          /nsamples

    where (raincountmn > 0.)
      Dvrmn        = Dvrmn / raincountmn
    elsewhere
      Dvrmn = 0.0
    end where
    where (preccountmn > 0.)
      prec_prcmn = prec_prcmn/preccountmn
    elsewhere
      prec_prcmn = 0.0
    end where
    
    ! #sb3 START
    !  for cloud columns
    cl_countmn = cl_countmn  /nsamples 
    ci_countmn = ci_countmn  /nsamples 
    hr_countmn = hr_countmn  /nsamples 
    hs_countmn = hs_countmn  /nsamples 
    hg_countmn = hg_countmn  /nsamples  
    ! for precipuitation
    prec_r_mn  = prec_r_mn   /nsamples 
    prec_i_mn  = prec_i_mn   /nsamples 
    prec_s_mn  = prec_s_mn   /nsamples 
    prec_g_mn  = prec_g_mn   /nsamples  
    ! #sb3 END
    
    ! #sb3 START
    !------------------------------------------
    ! calculation of aggregated statistics
    !------------------------------------------
    ! cloud fractions
    ! - liquid cloud fraction
    cfrac_l = cl_countmn
    ! - ice cloud fraction
    cfrac_i = ci_countmn
    ! - total cloud fraction
    do k=1,k1
     cfrac_tot(k) = max(cl_countmn(k),ci_countmn(k))
     ! this is just an approximation
     !-> later improve
    enddo
    
    ! process contributions
     ! - freezing
     !o dqi_hom  = -dqc_hom
     !o dqi_het  = -dqc_het
     !o dqg_het  = -dqr_het
     !o dni_hom  = -dnc_hom
     !o dni_het  = -dnc_het
     !o dng_het  = -dnr_het
     ! - collection of ice by snow
     !o dqi_col_si = - dqs_col_si
      
    ! aggregated process contribution
     ! - growth of ice hydrometeors
      !o dqi_rime = dqi_rime_ic
      !o dqs_rime = dqs_rime_sc! +dqs_rime_sr
      !o dqg_rime = dqg_rime_gc+dqg_rime_gr
      !o dqi_col  = 0.0  ! ice crystals after collection are not ice crystals anymore
      !o  dqs_col  = -dqi_agg+dqs_col_si
      !o dqg_col  = -(dqi_col_rig+dqr_col_rig)
      !o dni_col  = 0.0  ! ice crystals after collection are not ice crystals anymore 
      !o dns_col  = -dni_agg +dns_sc
      !o dng_col  = -dni_col_rig-dns_col_rsg
      !o dqi_cv   = -dqc_hom-dqc_het
      !o dqs_cv   =  0.0 ! currently no conversion to snow
      !o dqg_cv   = -dqr_het-dqi_cv_ig-dqs_cv_sg
      !o dni_cv   = -dnc_hom-dnc_het
      !o dns_cv   =  0.0 ! currently no conversion to snow
      !o dng_cv   = -dnr_het-dni_cv_ig-dns_cv_sg
      
     ! - impact on liquid hydrometeors 
      dqc_rime = -dqi_rime_ic-dqs_rime_sc-dqg_rime_gc
      dqr_rime = -dqg_rime_gr+dqr_col_rig+dqr_col_rsg
      dnc_rime = dnc_rime_ic + dnc_rime_sc+dnc_rime_gc
      dnr_rime = dnr_rime_gr+dnr_col_rig+dnr_col_rsg 
    
      
    
    ! thermodynamics 
    do k=1,k1
     ! - change in th due to freezing
     dth_freeze(k)=  (rlme/(cp*exnf(k)))*(-dqc_hom(k)           & ! homogeneous frezing
         -dqc_het(k)-dqr_het(k)                                 & ! heterogeneous freezing
         +dqi_rime_ic(k)+dqs_rime_sc(k)+dqg_rime_gc(k)          & ! riming by cloud
         +dqg_rime_gr(k)                                        & ! riming by rain
         -dqr_col_rsg(k)-dqr_col_rig(k)                       )   ! collection r+i, s+i
     !
     ! - change in th due to melting
     dth_melt (k) = ( rlme/(cp*exnf(k)))*(                      &
         dqi_me (k)+dqs_me(k)+dqg_me(k)                         & ! partial melting 
         +dqi_ev (k)+dqs_ev(k)+dqg_ev(k)                     & ! partial followed by evaporation 
         +dqi_eme_ic(k)+dqs_eme_sc(k)+dqg_eme_gc(k)             & ! enhanced melting by cloud water
         +dqi_eme_ri(k)+dqs_eme_rs(k)+dqg_eme_gr(k)           )   ! enhanced melting by rain
     !
     ! - change in th due to condensation
     dth_cond (k) = ( rlv/(cp*exnf(k)))*( dqc_nuc_av(k)         & ! cloud nucleation
                      +dqc_sadjpl (k)                        )    ! cloud growth
     !
     ! - change in th due to evaporation 
     dth_ev (k) = ( rlv/(cp*exnf(k)))* ( dqc_sadjneg (k)       & ! evaporation of clouds
        +dqr_ev(k)+dqi_ev(k)+dqs_ev(k)+dqg_ev(k)         )   ! evaporation of hydrometeors
     !
     ! - change in th due to deposition
     dth_dep (k) = (rlvi/(cp*exnf(k)))* ( dqi_inuc (k)          & ! cloud nucleation
        + dqi_dep (k)+dqs_dep(k)+dqg_dep(k)                   )   ! ice deposition on existing particles
     !
     ! - change in th due to sublimation
     dth_sub (k) =(rlvi/(cp*exnf(k)))*(dqi_sub (k)+dqs_sub(k)+dqg_sub(k)) ! sublimation of ice particles
     !
     ! - change in th due to mphys in total
     dth_mphys(k) = dth_freeze(k)+dth_melt(k)+dth_cond(k)+dth_ev(k)+dth_sub(k)
     ! 
    enddo
    ! #sb3 END
    
    
    !-----------------
    ! text outputs 
    !--------------------
    
    if (myid == 0) then
    open (ifoutput,file='precep.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/2A/2A)')             &
      '#------------------------------------------------------------'     &
      ,'------------'                 &
      ,'#               --------   PRECIPITATION ------    '       &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   RHO(k)  PRES  |CLOUDCOVER  ECHORAINRATE  PRECCOUNT '   &
      ,'    NRRAIN      RAINCOUNT     PREC(k)     <Dvr(k)>     <qr(k)>'   &
      ,'#      (M)             (MB)  |----------  ---W/M2----   --------- '   &
      ,'    ------      ---------     -------     --------    ---------'   &
      ,'#-----------------------------------------------------------------'   &
      ,'---------------------------------------------------------------'
    write(ifoutput,'(I4,F10.2,F8.3,F7.1,8E13.5)') &
      (k          , &
      zf    (k)      , &
      rhof    (k)      , &
      presf    (k)/100.    , &
      cloudcountmn  (k)      , &
      prec_prcmn  (k)*rhof(k)*rlv  , &
      preccountmn  (k)      , &
      Nrrainmn  (k)      , &
      raincountmn  (k)      , &
      precmn    (k)*rhof(k)*rlv  , &
      Dvrmn    (k)      , &
      qrmn    (k)      , &
      k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='nptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S NRAIN ------    '     &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (#/M3/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      Npmn    (k,iauto)    , &
      Npmn    (k,iaccr)    , &
      Npmn    (k,ised)    , &
      Npmn    (k,ievap)    , &
      sum(Npmn  (k,2:nrfields))    , &
      k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='qlptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S QRAIN ------    '   &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (KG/KG/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      qlpmn    (k,iauto)    , &
      qlpmn    (k,iaccr)    , &
      qlpmn    (k,ised)    , &
      qlpmn    (k,ievap)    , &
      sum(qlpmn  (k,2:nrfields))    , &
                        k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='qtptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S QTP ------    '   &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (KG/KG/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      qtpmn    (k,iauto)    , &
      qtpmn    (k,iaccr)    , &
      qtpmn    (k,ised)    , &
      qtpmn    (k,ievap)    , &
      sum    (qtpmn(k,2:nrfields))  , &
      k=1,kmax)
      close(ifoutput)
      
    !------------------------------
    ! NetCDF output
    !------------------------------  
      if (lnetcdf) then
        vars(:, 1) = cloudcountmn
        vars(:, 2) = prec_prcmn  (:)*rhof(:)*rlv
        vars(:, 3) = preccountmn  (:)
        vars(:, 4) = Nrrainmn  (:)
        vars(:, 5) = raincountmn  (:)
        vars(:, 6) = precmn    (:)*rhof(:)*rlv
        vars(:, 7) = Dvrmn    (:)
        vars(:, 8) = qrmn    (:)
        vars(:, 9) =Npmn    (:,iauto)
        vars(:,10) =Npmn    (:,iaccr)
        vars(:,11) =Npmn    (:,ised)
        vars(:,12) =Npmn    (:,ievap)
        do k=1,k1
        vars(k,13) =sum(Npmn  (k,2:nrfields))
        enddo
        vars(:,14) =qlpmn    (:,iauto)
        vars(:,15) =qlpmn    (:,iaccr)
        vars(:,16) =qlpmn    (:,ised)
        vars(:,17) =qlpmn    (:,ievap)
        do k=1,k1
        vars(k,18) =sum(qlpmn  (k,2:nrfields))
        enddo
        vars(:,19) =qtpmn    (:,iauto)
        vars(:,20) =qtpmn    (:,iaccr)
        vars(:,21) =qtpmn    (:,ised)
        vars(:,22) =qtpmn    (:,ievap)
        do k=1,k1
        vars(k,23) =sum(qtpmn  (k,2:nrfields))
        enddo
        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof,kmax)
        ! #sb3 START
        ! recording added statistics
        vars_mphys (:,nbaspout+1)  = dnc_nuc
        vars_mphys (:,nbaspout+2)  = dqc_nuc
        vars_mphys (:,nbaspout+3)  = dni_inuc
        vars_mphys (:,nbaspout+4)  = dqi_inuc
        ! - deposition and sublimation 
        vars_mphys (:,nbaspout+5)  = dqi_dep
        vars_mphys (:,nbaspout+6)  = dqs_dep
        vars_mphys (:,nbaspout+7)  = dqg_dep
        vars_mphys (:,nbaspout+8)  = dqi_sub
        vars_mphys (:,nbaspout+9)  = dqs_sub
        vars_mphys (:,nbaspout+10) = dqg_sub 
        ! - freezing 
        ! -- homogeneous freezing 
        vars_mphys (:,nbaspout+11) = -dqc_hom ! dqi_hom
        vars_mphys (:,nbaspout+12) = -dnc_hom ! dni_hom
        ! -- heterogeneous freezing 
        vars_mphys (:,nbaspout+13) = -dqc_het ! dqi_het
        vars_mphys (:,nbaspout+14) = -dnc_het ! dni_het        
        vars_mphys (:,nbaspout+15) = -dqr_het ! dqg_het
        vars_mphys (:,nbaspout+16) = -dnr_het ! dng_het      
        ! - warm microphysics
        ! - saturation aqdjustment tendency
        vars_mphys (:,nbaspout+17) = dqc_sadjpl
        vars_mphys (:,nbaspout+18) = dqc_sadjneg
        ! -- cloud water self-collection 
        vars_mphys (:,nbaspout+19) = dnc_sc
        ! -- cloud water autoconversion
        vars_mphys (:,nbaspout+20) = -dqr_au ! dqc_au
        vars_mphys (:,nbaspout+21) = dnc_au
        ! -- cloud water accretion
        vars_mphys (:,nbaspout+22) = -dqr_ac !dqc_ac
        vars_mphys (:,nbaspout+23) = dnc_ac
        ! -- rain selfcollection and breakup
        vars_mphys (:,nbaspout+24) = dnr_sc
        vars_mphys (:,nbaspout+25) = dnr_br
        ! -- rain evaporation 
        vars_mphys (:,nbaspout+26) = dqr_ev
        vars_mphys (:,nbaspout+27) = dnr_ev
        ! - cold aggregation and self-collection
        ! -- ice aggregation 
        vars_mphys (:,nbaspout+28) = dqi_agg
        vars_mphys (:,nbaspout+29) = dni_agg
        ! -- snow self-collection
        vars_mphys (:,nbaspout+30) = dns_sc
        ! - loss of liquid water due to riming
        vars_mphys (:,nbaspout+31) = dqc_rime
        vars_mphys (:,nbaspout+32) = dnc_rime
        vars_mphys (:,nbaspout+33) = dqr_rime
        vars_mphys (:,nbaspout+34) = dnr_rime    
        ! - growth by riming
        vars_mphys (:,nbaspout+35) = dqi_rime_ic
        vars_mphys (:,nbaspout+36) = dqs_rime_sc
        !vars_mphys (:,nbaspout+37) = dqs_rime_sr
        vars_mphys (:,nbaspout+37) = dqg_rime_gc
        vars_mphys (:,nbaspout+38) = dqg_rime_gr
        ! - collisions
        ! -- collision conversions 
        vars_mphys (:,nbaspout+39) = dqi_col_rig
        vars_mphys (:,nbaspout+40) = dni_col_rig
        vars_mphys (:,nbaspout+41) = dqr_col_rig
        vars_mphys (:,nbaspout+42) = dnr_col_rig 
        !
        vars_mphys (:,nbaspout+43) = dqs_col_rsg
        vars_mphys (:,nbaspout+44) = dns_col_rsg
        vars_mphys (:,nbaspout+45) = dqr_col_rsg
        vars_mphys (:,nbaspout+46) = dnr_col_rsg        
        ! -- collision collection 
        vars_mphys (:,nbaspout+47) = -dqs_col_si
        vars_mphys (:,nbaspout+48) = dni_col_si
        vars_mphys (:,nbaspout+49) = dqs_col_gsg
        vars_mphys (:,nbaspout+50) = dns_col_gsg
        ! -- enhanced melting
        vars_mphys (:,nbaspout+51) = dqi_eme_ic
        vars_mphys (:,nbaspout+52) = dni_eme_ic  
        vars_mphys (:,nbaspout+53) = dqi_eme_ri
        vars_mphys (:,nbaspout+54) = dni_eme_ri  
        vars_mphys (:,nbaspout+55) = dqs_eme_sc
        vars_mphys (:,nbaspout+56) = dns_eme_sc  
        vars_mphys (:,nbaspout+57) = dqs_eme_rs
        vars_mphys (:,nbaspout+58) = dns_eme_rs        
        vars_mphys (:,nbaspout+59) = dqg_eme_gc
        vars_mphys (:,nbaspout+60) = dng_eme_gc  
        vars_mphys (:,nbaspout+61) = dqg_eme_gr
        vars_mphys (:,nbaspout+62) = dng_eme_gr         
        !
        ! - evaporation and melting
        ! -- melting
        vars_mphys (:,nbaspout+63) = dqi_me
        vars_mphys (:,nbaspout+64) = dni_me
        vars_mphys (:,nbaspout+65) = dqs_me
        vars_mphys (:,nbaspout+66) = dns_me
        vars_mphys (:,nbaspout+67) = dqg_me
        vars_mphys (:,nbaspout+68) = dng_me
        ! -- evaporation
        vars_mphys (:,nbaspout+69) = dqi_ev
        vars_mphys (:,nbaspout+70) = dni_ev
        vars_mphys (:,nbaspout+71) = dqs_ev
        vars_mphys (:,nbaspout+72) = dns_ev
        vars_mphys (:,nbaspout+73) = dqg_ev
        vars_mphys (:,nbaspout+74) = dng_ev      
        !
        ! - conversions and other processes
        vars_mphys (:,nbaspout+75) = dqi_cv_ig
        vars_mphys (:,nbaspout+76) = dni_cv_ig
        vars_mphys (:,nbaspout+77) = dqs_cv_sg
        vars_mphys (:,nbaspout+78) = dns_cv_sg
        ! -- ice multiplication
        vars_mphys (:,nbaspout+79) = dqi_mul
        vars_mphys (:,nbaspout+80) = dni_mul 
        ! -- satuartion adjustments 
        vars_mphys (:,nbaspout+81) = dnc_sadj
        !
        !-> add extra processes here
        !
        !   then adjust nmphyout
        ! - sedimentation tendency 
        vars_mphys (:,ncolout+1)  = dqc_sed
        vars_mphys (:,ncolout+2)  = dqr_sed      
        vars_mphys (:,ncolout+3)  = dqi_sed
        vars_mphys (:,ncolout+4)  = dqs_sed
        vars_mphys (:,ncolout+5)  = dqg_sed
        ! -- and for numbers
        vars_mphys (:,ncolout+6)  = dnc_sed
        vars_mphys (:,ncolout+7)  = dnr_sed   
        vars_mphys (:,ncolout+8)  = dni_sed
        vars_mphys (:,ncolout+9)  = dns_sed
        vars_mphys (:,ncolout+10) = dng_sed
        ! - total tendency
        ! - content of hydrometeors
        vars_mphys (:,ncolout+11) = dqc_tot
        vars_mphys (:,ncolout+12) = dqr_tot       
        vars_mphys (:,ncolout+13) = dqi_tot
        vars_mphys (:,ncolout+14) = dqs_tot
        vars_mphys (:,ncolout+15) = dqg_tot
        ! -- numbers of hydrometeors
        vars_mphys (:,ncolout+16) = dnc_tot
        vars_mphys (:,ncolout+17) = dnr_tot
        vars_mphys (:,ncolout+18) = dni_tot
        vars_mphys (:,ncolout+19) = dns_tot
        vars_mphys (:,ncolout+20) = dng_tot
        ! -- numbers of CCN
        vars_mphys (:,ncolout+21) = dn_ccn_tot
        ! -- temperature tendency
        vars_mphys (:,nmphyout+1) = dthl_mphys
        vars_mphys (:,nmphyout+2) = dth_mphys
        vars_mphys (:,nmphyout+3) = dth_freeze
        vars_mphys (:,nmphyout+4) = dth_melt
        vars_mphys (:,nmphyout+5) = dth_cond
        vars_mphys (:,nmphyout+6) = dth_ev
        vars_mphys (:,nmphyout+7) = dth_dep
        vars_mphys (:,nmphyout+8) = dth_sub
        ! -- cloud fraction
        vars_mphys (:,nmphyout+9) = cfrac_l
        vars_mphys (:,nmphyout+10) = cfrac_i
        vars_mphys (:,nmphyout+11) = cfrac_tot
        ! -- precipitation 
        vars_mphys (:,nmphyout+12) = prec_r_mn ! rain_rate
        vars_mphys (:,nmphyout+13) = prec_i_mn ! ice_rate
        vars_mphys (:,nmphyout+14) = prec_s_mn ! snow_rate
        vars_mphys (:,nmphyout+15) = prec_g_mn ! graupel_rate
        !
        ! - content of hydrometeors
        vars_mphys (:,nmphyout+16) = dqc_totf
        vars_mphys (:,nmphyout+17) = dqr_totf       
        vars_mphys (:,nmphyout+18) = dqi_totf
        vars_mphys (:,nmphyout+19) = dqs_totf
        vars_mphys (:,nmphyout+20) = dqg_totf
        ! -- numbers of hydrometeors
        vars_mphys (:,nmphyout+21) = dnc_totf
        vars_mphys (:,nmphyout+22) = dnr_totf
        vars_mphys (:,nmphyout+23) = dni_totf
        vars_mphys (:,nmphyout+24) = dns_totf
        vars_mphys (:,nmphyout+25) = dng_totf       
        ! #sb3 END
        !
        ! call writestat_nc(ncid_prof,nvar,ncmname,vars(1:kmax,:),nrec_prof,kmax)
        call writestat_nc(ncid_mphys,1,tncmname,(/rtimee/),nrec_mphys,.true.)
        call writestat_nc(ncid_mphys,nvar_mphys,ncmname,vars_mphys(1:kmax,:),nrec_mphys,kmax)
      end if

    end if

    !------------------------------
    ! reseting variables
    !------------------------------  
    
     cloudcountmn    = 0.0
     raincountmn    = 0.0
     preccountmn    = 0.0
     prec_prcmn    = 0.0
     Dvrmn      = 0.0
     Nrrainmn    = 0.0
     precmn      = 0.0
     qrmn      = 0.0
     Npmn      = 0.0
     qlpmn      = 0.0
     qtpmn      = 0.0
    ! #sb3 START
     cl_countmn  = 0.0
     ci_countmn  = 0.0
     hr_countmn  = 0.0
     hs_countmn  = 0.0
     hg_countmn  = 0.0
     prec_r_mn   = 0.0
     prec_i_mn   = 0.0
     prec_s_mn   = 0.0
     prec_g_mn   = 0.0   
     ! and for variables describing processes
     dnc_nuc      = 0.0 
     dqc_nuc      = 0.0
     dni_inuc     = 0.0
     dqi_inuc     = 0.0
     dqi_dep      = 0.0
     dqi_sub      = 0.0
     dqs_dep      = 0.0
     dqs_sub      = 0.0
     dqg_dep      = 0.0
     dqg_sub      = 0.0
     dnc_hom      = 0.0
     dqc_hom      = 0.0
     dnc_het      = 0.0
     dqc_het      = 0.0
     dnr_het      = 0.0
     dqr_het      = 0.0
     dnc_sadj     = 0.0
     dqc_sadjpl   = 0.0
     dqc_sadjneg  = 0.0
     dnc_au       = 0.0
     dnr_au       = 0.0
     dqr_au       = 0.0
     dqr_ac       = 0.0
     dnc_ac       = 0.0
     dnr_sc       = 0.0
     dnr_br       = 0.0
     dqr_ev     = 0.0
     dnr_ev     = 0.0
     dqi_agg      = 0.0
     dni_agg      = 0.0
     dns_sc       = 0.0
     dqi_rime_ic  = 0.0
     dqs_rime_sc  = 0.0
     dqg_rime_gc  = 0.0
     dnc_rime_ic  = 0.0
     dnc_rime_sc  = 0.0
     dnc_rime_gc  = 0.0
     ! dqs_rime_sr  = 0.0
     dqg_rime_gr  = 0.0
     ! dnr_rime_sr  = 0.0
     dnr_rime_gr  = 0.0
     dqi_col_rig  = 0.0
     dqr_col_rig  = 0.0
     dni_col_rig  = 0.0
     dnr_col_rig  = 0.0
     dqs_col_rsg  = 0.0
     dqr_col_rsg  = 0.0
     dns_col_rsg  = 0.0
     dnr_col_rsg  = 0.0    
     dqs_col_si   = 0.0
     dni_col_si   = 0.0
     dqs_col_gsg  = 0.0
     dns_col_gsg  = 0.0
     dqi_cv_ig    = 0.0
     dni_cv_ig    = 0.0
     dqs_cv_sg    = 0.0
     dns_cv_sg    = 0.0
     dqc_sed      = 0.0
     dqi_sed      = 0.0
     dqr_sed      = 0.0
     dqs_sed      = 0.0
     dqg_sed      = 0.0
     dnc_sed      = 0.0
     dni_sed      = 0.0
     dnr_sed      = 0.0
     dns_sed      = 0.0
     dng_sed      = 0.0 
     dqc_tot      = 0.0
     dqi_tot      = 0.0
     dqr_tot      = 0.0
     dqs_tot      = 0.0
     dqg_tot      = 0.0
     dnc_totf     = 0.0
     dni_totf     = 0.0
     dnr_totf     = 0.0
     dns_totf     = 0.0
     dng_totf     = 0.0  
     dqc_totf     = 0.0
     dqi_totf     = 0.0
     dqr_totf     = 0.0
     dqs_totf     = 0.0
     dqg_totf     = 0.0
     dnc_totf     = 0.0
     dni_totf     = 0.0
     dnr_totf     = 0.0
     dns_totf     = 0.0
     dng_totf     = 0.0     
     dn_ccn_tot   = 0.0
     dthl_mphys   = 0.0 
     dqi_eme_ic   = 0.0
     dqs_eme_sc   = 0.0
     dqg_eme_gc   = 0.0
     dni_eme_ic   = 0.0
     dns_eme_sc   = 0.0
     dng_eme_gc   = 0.0
     dqi_eme_ri   = 0.0
     dqs_eme_rs   = 0.0
     dqg_eme_gr   = 0.0
     dni_eme_ri   = 0.0
     dns_eme_rs   = 0.0
     dng_eme_gr   = 0.0
     dqi_me       = 0.0
     dqs_me       = 0.0
     dqg_me       = 0.0 
     dni_me       = 0.0
     dns_me       = 0.0
     dng_me       = 0.0
     dqi_ev      = 0.0
     dqs_ev      = 0.0
     dqg_ev      = 0.0 
     dni_ev      = 0.0
     dns_ev      = 0.0
     dng_ev      = 0.0
     dni_mul      = 0.0
     dqi_mul      = 0.0  
    ! #sb3 END         
     
    ! #sb3 START 
    deallocate( vars_mphys )
    deallocate( cfrac_l,cfrac_i,cfrac_tot                    &
               ,dqc_rime,dqr_rime                            &
               ,dnc_rime,dnr_rime                            &
               ,dth_mphys,dth_freeze,dth_melt                &
               ,dth_ev,dth_cond,dth_dep,dth_sub             &
              )
           !    ,dqi_hom,dqi_het,dqg_het                      &
           !    ,dni_hom,dni_het,dng_het                      &
           !    ,dqi_col_si                                   &
           !    ,dqi_rime,dqs_rime,dqg_rime                   &
           !    ,dqi_col,dqs_col ,dqg_col                     &
           !    ,dni_col,dns_col ,dng_col                     &
           !    ,dqi_cv,dqs_cv,dqg_cv                         &
           !    ,dni_cv,dns_cv,dng_cv                         &
           !   )
    ! #sb3 END       

  end subroutine writebulkmicrostat3

!------------------------------------------------------------------------------!

  subroutine exitbulkmicrostat3
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none
    if (.not. lmicrostat)  return
    
    ! calling initialisation of hydrometeor dump
    !if(l_hdump) then 
      call exithydrodump
    ! endif
    if(lnetcdf .and. myid==0) call exitstat_nc(ncid_mphys)
    ! profile file closed separately
    
    deallocate(Npav      , &
         Npmn      , &
         qlpav    , &
         qlpmn    , &
         qtpav    , &
         qtpmn    )
    deallocate(precavl    , &
         precav    , &
         precmn    , &
         preccountavl    , &
         preccountav    , &
         preccountmn    , &
         prec_prcavl    , &
         prec_prcav    , &
         prec_prcmn    , &
         cloudcountavl  , &
         cloudcountav    , &
         cloudcountmn    , &
         raincountavl    , &
         raincountav    , &
         raincountmn    , &
         Nrrainavl    , &
         Nrrainav    , &
         Nrrainmn    , &
         qravl    , &
         qrav      , &
         qrmn      , &
         Dvravl    , &
         Dvrav    , &
         Dvrmn)
    ! #sb3 START
    ! deallocating new variables
     deallocate( dnc_nuc ,dqc_nuc                             &
              ,dni_inuc ,dqi_inuc                             &
              ,dqi_dep ,dqi_sub ,dqs_dep ,dqs_sub             &
              ,dqg_dep ,dqg_sub                               &
              ,dnc_hom ,dqc_hom                               &
              ,dnc_het ,dqc_het ,dnr_het ,dqr_het             &
              ,dnc_sadj ,dqc_sadjpl, dqc_sadjneg              &
              ,dnc_sc                                         &
              ,dnc_au ,dnr_au ,dqr_au                         &
              ,dqr_ac ,dnc_ac                                 &
              ,dnr_sc ,dnr_br                                 &
              ,dqr_ev ,dnr_ev                             &
              ,dqi_agg ,dni_agg                               &
              ,dns_sc                                         &
              ,dqi_rime_ic ,dqs_rime_sc ,dqg_rime_gc          &
              ,dnc_rime_ic ,dnc_rime_sc ,dnc_rime_gc          &
              ,dqg_rime_gr,dnr_rime_gr                        &
              ,dqi_col_rig ,dqr_col_rig                       &
              ,dni_col_rig ,dnr_col_rig                       &
              ,dqs_col_rsg ,dqr_col_rsg                       &
              ,dns_col_rsg ,dnr_col_rsg                       &              
              ,dqs_col_si ,dni_col_si                         &
              ,dqs_col_gsg ,dns_col_gsg                       &
              ,dqi_cv_ig ,dni_cv_ig ,dqs_cv_sg ,dns_cv_sg     &              
              ,dqc_sed,dqr_sed,dqi_sed,dqs_sed,dqg_sed        &
              ,dnc_sed,dnr_sed,dni_sed,dns_sed,dng_sed        &
              ,dqc_tot,dqi_tot,dqr_tot,dqs_tot,dqg_tot        &
              ,dnc_tot,dni_tot,dnr_tot,dns_tot, dng_tot       &
              ,dqc_totf,dqi_totf,dqr_totf,dqs_totf,dqg_totf   & !#
              ,dnc_totf,dni_totf,dnr_totf,dns_totf, dng_totf  & !#        
              ,dn_ccn_tot                                     &
              ,dthl_mphys                                     &
              ,dqi_eme_ic,dqs_eme_sc,dqg_eme_gc               &
              ,dni_eme_ic,dns_eme_sc,dng_eme_gc               &
              ,dqi_eme_ri,dqs_eme_rs,dqg_eme_gr               &
              ,dni_eme_ri,dns_eme_rs,dng_eme_gr               &
              ,dqi_me ,dqs_me ,dqg_me                         & 
              ,dni_me ,dns_me ,dng_me                         &
              ,dqi_ev,dqs_ev,dqg_ev                        & 
              ,dni_ev,dns_ev,dng_ev                        &
              ,dni_mul,dqi_mul                                &
             )
     deallocate( cl_countmn,hr_countmn                        &
              ,ci_countmn,hs_countmn,hg_countmn               &
              ,prec_r_mn,prec_i_mn,prec_s_mn,prec_g_mn        &
             )
             ! ,cfrac_l,cfrac_i,cfrac_tot                <- in writebulkmicrostat3
             !,rain_rate,ice_rate,snow_rate,graupel_rate      &
    !#sb3 END
    ! #sb3 START deallocation
     deallocate( cl_countavl, hr_countavl                     &
               ,ci_countavl,hs_countavl,hg_countavl           &
               ,prec_r_avl,prec_i_avl,prec_s_avl,prec_g_avl   &
              )
     deallocate(dnc_nuc_av,dqc_nuc_av,dni_inuc_av,dqi_inuc_av  &
              ,dqi_dep_av,dqi_sub_av,dqs_dep_av,dqs_sub_av     &
              ,dqg_dep_av,dqg_sub_av                           &
              ,dnc_hom_av,dqc_hom_av                           &
              ,dnc_het_av,dqc_het_av,dnr_het_av,dqr_het_av     &
              ,dnc_sadj_av,dqc_sadjpl_av,dqc_sadjneg_av        &
              ,dnc_sc_av                                       &
              ,dnc_au_av,dnr_au_av, dqr_au_av                  &
              ,dqr_ac_av, dnc_ac_av                            &
              ,dnr_sc_av,dnr_br_av                             &
              ,dqr_ev_av,dnr_ev_av                         &
              ,dqi_agg_av,dni_agg_av                           &
              ,dns_sc_av                                       &
              ,dqi_rime_icav,dqs_rime_scav,dqg_rime_gcav       &
              ,dnc_rime_icav,dnc_rime_scav,dnc_rime_gcav       &
              ,dqg_rime_grav,dnr_rime_grav                     &
              ,dqi_col_rigav,dqr_col_rigav                     &
              ,dni_col_rigav,dnr_col_rigav                     &
              ,dqs_col_rsgav,dqr_col_rsgav                     &
              ,dns_col_rsgav,dnr_col_rsgav                     &              
              ,dqs_col_siav,dni_col_siav                       &
              ,dqs_col_gsgav ,dns_col_gsgav                    &
              ,dqi_cv_igav,dni_cv_igav,dqs_cv_sgav,dns_cv_sgav &
              ,dn_ccn_av                                       &
              ,dqc_sedav,dqr_sedav                             &
              ,dqi_sedav,dqs_sedav,dqg_sedav                   &
              ,dnc_sedav,dnr_sedav                             &
              ,dni_sedav,dns_sedav,dng_sedav                   &
              ,dqc_totav,dqr_totav                             &
              ,dqi_totav,dqs_totav,dqg_totav                   &
              ,dnc_totav, dnr_totav                            &
              ,dni_totav,dns_totav, dng_totav                  &
              ,dqc_totfav,dqr_totfav                           & !#
              ,dqi_totfav,dqs_totfav,dqg_totfav                & !#             
              ,dnc_totfav, dnr_totfav                          & !#
              ,dni_totfav,dns_totfav, dng_totfav               & !#
              ,dqi_eme_icav,dqs_eme_scav,dqg_eme_gcav          &
              ,dni_eme_icav,dns_eme_scav,dng_eme_gcav          &
              ,dqi_eme_riav,dqs_eme_rsav,dqg_eme_grav          &
              ,dni_eme_riav,dns_eme_rsav,dng_eme_grav          &
              ,dqi_me_av ,dqs_me_av ,dqg_me_av                 & 
              ,dni_me_av ,dns_me_av ,dng_me_av                 &
              ,dqi_ev_av,dqs_ev_av,dqg_ev_av                & 
              ,dni_ev_av,dns_ev_av,dng_ev_av                &
              ,dni_mul_av,dqi_mul_av                             &
              )
      deallocate( cl_countav,hr_countav                        &
              ,ci_countav,hs_countav,hg_countav                &
              ,prec_r_av,prec_i_av,prec_s_av,prec_g_av         &
              )
     ! deallocate(dth_freeze,dth_melt,dth_mphys, dthl_mphys)       
    ! #sb3 END
    

  end subroutine exitbulkmicrostat3

!------------------------------------------------------------------------------!
! subroutines for field output
! - strongly based on fielddump

!> Initializing cloudddump.
!> based on namelist from fielddump initializing the variables
!> some variables taken from fielddump 
! 
  subroutine inithydrodump
    use modmpi,   only :myid,my_real,comm3d,mpi_logical,mpi_integer,myidx,myidy
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options   & 
                  ,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime,tres
    ! use modfielddump, only: dtav, idtav, tnext, klow, khigh 
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
    use modmicrodata3, only: l_hdump, def_hdump_dtav         &
                        ,iq_cl,iq_hr,iq_ci,iq_hs,iq_hg       &
                        ,in_cl,in_hr,in_ci,in_hs,in_hg       &
                        ,in_cc                            
    implicit none
    integer :: ierr
    
    ! --- setting ------------------------
    ! -> some of them to bae later read from a namelist 
    hdump_klow = 1 ! 2
    hdump_khigh = kmax-10 ! <- adjust later !! 
    hdump_dtav = def_hdump_dtav ! <- adjust later through namelist
    
    ! --- initialisation --------------------
    hdump_idtav = hdump_dtav/tres
    hdump_tnext      = hdump_idtav   +btime
    if(.not.(l_hdump)) return
    ! done by fielddump:
    dt_lim = min(dt_lim,hdump_tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if
    
    
    ! which hydrometeor specie are we outputting 
    isv_vec( 1) =   iq_cl
    isv_vec( 2) =   iq_hr
    isv_vec( 3) =   iq_ci
    isv_vec( 4) =   iq_hs    
    isv_vec( 5) =   iq_hg
    isv_vec( 6) =   in_cl
    isv_vec( 7) =   in_hr
    isv_vec( 8) =   in_ci
    isv_vec( 9) =   in_hs    
    isv_vec(10) =   in_hg   
     
    if (lnetcdf) then
      write(hdump_fname,'(A,i3.3,A,i3.3,A)') 'hydrodump.', myidx, '.', myidy, '.xxx.nc'! 'fielddump.', myidx, '.', myidy, '.xxx.nc'
      hdump_fname(19:21) = cexpnr
      call ncinfo(hdump_tncname(1,:),'time','Time','s','time')
      call ncinfo(hdump_ncname( 1,:),'q_c','cloud liquid water specific humidity','kg/kg','tttt')
      call ncinfo(hdump_ncname( 2,:),'q_r','rain water specific humidity'        ,'kg/kg','tttt')
      call ncinfo(hdump_ncname( 3,:),'q_i','cloud ice water specific humidity'   ,'kg/kg','tttt')
      call ncinfo(hdump_ncname( 4,:),'q_s','snow water specific humidity'        ,'kg/kg','tttt')     
      call ncinfo(hdump_ncname( 5,:),'q_g','graupel water specific humidity'     ,'kg/kg','tttt') 
      call ncinfo(hdump_ncname( 6,:),'n_c','cloud liquid water number'   ,'0e/kg','tttt')
      call ncinfo(hdump_ncname( 7,:),'n_r','rain water number','0e/kg','tttt')
      call ncinfo(hdump_ncname( 8,:),'n_i','cloud ice water number','0e/kg','tttt')
      call ncinfo(hdump_ncname( 9,:),'n_s','snow water number','0e/kg','tttt')
      call ncinfo(hdump_ncname(10,:),'n_g','graupel water number','0e/kg','tttt')       
      call open_nc(hdump_fname,hdump_ncid,hdump_nrec,n1=imax,n2=jmax,n3=hdump_khigh-hdump_klow+1)
      if (hdump_nrec==0) then
        call define_nc( hdump_ncid,1, hdump_tncname)
        call writestat_dims_nc(hdump_ncid)
      end if
     call define_nc(hdump_ncid, hdump_nvar,hdump_ncname)
    end if

  end subroutine inithydrodump
  
!> Clean up when leaving the run
  subroutine exithydrodump
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modmicrodata3, only: l_hdump 
    implicit none

    if(l_hdump .and. lnetcdf) call exitstat_nc(hdump_ncid)
  end subroutine exithydrodump  

  
!> Do hydrodump. Collect data to truncated (2 byte) integers, and write them to file
!> or use 4 byte integers?
  subroutine hydrodump
    use modfields, only : svm     ! um,vm,wm,thlm,qtm,ql0,
    ! use modsurfdata,only : thls,qts,thvs
    use modglobal, only : imax,i1,ih,jmax,j1,jh,k1,rk3step,&
                          timee,dt_lim,cexpnr,ifoutput,rtimee
    use modmpi,    only : myid,cmyidx, cmyidy
    use modstat_nc, only : lnetcdf, writestat_nc
    ! use modmicrodata,  only : imicro, imicro_none          
    use modmicrodata3, only : l_hdump, l_hbinary,l_hdiracc      &
                      ,iq_cl,iq_ci,iq_hr,iq_hs,iq_hg            &
                      ,in_cl,in_ci,in_hr,in_hs,in_hg            & ! &
                      ,in_cc
                      ! ,mul_fac, hdump_byte
    implicit none

    ! integer(KIND=selected_int_kind(4)), allocatable :: field(:,:,:),hyvars(:,:,:,:)
    real, allocatable :: hyfield(:,:,:),hyvars(:,:,:,:)
    integer i,j,k, isvout, iisv
    integer :: writecounter = 1
    integer :: reclength
    ! real    :: mul_fac  ! multiplication factor for variables


    
    if (.not. l_hdump) return
    if (rk3step/=3) return

    if(timee<hdump_tnext) then
      ! this is done in fielddump :
      ! dt_lim = min(dt_lim,tnext-timee)
      !
      return
    end if
    
    ! this is done in fielddump : 
    hdump_tnext = hdump_tnext+hdump_idtav
    dt_lim = minval((/dt_lim,hdump_tnext-timee/))
    !

    allocate(hyfield(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(hyvars(imax,jmax,hdump_khigh-hdump_klow+1,hdump_nvar))

    reclength = imax*jmax*(hdump_khigh-hdump_klow+1)*2

    ! loop over variable we want ot output
    do isvout =1,hdump_nvar   ! hdump_nvar
    
      ! get the index of variable in svm array
      iisv = isv_vec(isvout) 
      
    
     !   if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        ! field(i,j,k) = NINT(mul_fac*svm(i,j,k,iqr),2)
        hyfield(i,j,k) = svm(i,j,k,iisv)
      enddo
      enddo
      enddo
    ! else
    !  field = 0.
    ! endif

    if (lnetcdf) hyvars(:,:,:,isvout) = hyfield(2:i1,2:j1,hdump_klow:hdump_khigh)
    if (l_hbinary) then
      if (l_hdiracc) then
        open (ifoutput,file='wb_hdump.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) hyfield(2:i1,2:j1,hdump_klow:hdump_khigh)
      else
        open  (ifoutput,file='wb_hdump.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((hyfield(i,j,k),i=2,i1),j=2,j1),k=hdump_klow,hdump_khigh)
      end if
      close (ifoutput)
    endif
    
    enddo ! hdump_nvar

    if(lnetcdf) then
      call writestat_nc(hdump_ncid,1,hdump_tncname,(/rtimee/),hdump_nrec,.true.)
      call writestat_nc(hdump_ncid,hdump_nvar,hdump_ncname,hyvars,hdump_nrec,imax,jmax,hdump_khigh-hdump_klow+1)
    end if

    writecounter=writecounter+1

    deallocate(hyfield,hyvars)

  end subroutine hydrodump

!------------------------------------------------------------------------------

end module
