!> \file modbulkmicro3.f90
!!
!!  Bulk microphysics
!>
!! Calculates bulk microphysics using a two moment scheme of S&B, 2006
!! further sources:
!! \see  Seifert and Beheng (Atm. Res., 2001)
!! \see  Seifert and Beheng (Met Atm Phys, 2006)
!! \see  Stevens and Seifert (J. Meteorol. Soc. Japan, 2008)  (rain sedim, mur param)
!! \see  Seifert (J. Atm Sc., 2008) (rain evap)
!! \see  Khairoutdinov and Kogan (2000) (drizzle param : auto, accr, sedim, evap)
!!  \author Olivier Geoffroy, K.N.M.I.
!!  \author Margreet van Zanten, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Jan Chylik, IGM U.z.Koeln
!!  \par Revision list
!! \todo documentation
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


module modbulkmicro3
!*********************************************************************
!
!   Amount of liquid water is splitted into cloud water and precipitable
!   water (as such it is a two moment scheme). Cloud droplet number conc. is
!   fixed in place and time.
!
!   same rhof value used for some diagnostics calculation (in modbulkmicrostat, modtimestat)
!
!   Cond. sampled timeav averaged profiles are weighted with fraction of condition,
!   similarly as is done in sampling.f90
!
!   bulkmicro is called from *modmicrophysics*
!
!*********************************************************************
!  use modmicrodata
!  use modmicrodata3

  implicit none
  real :: gamma25
  real :: gamma3
  real :: gamma35
  contains

!> Initializes and allocates the arrays
  subroutine initbulkmicro3
    use modglobal, only : ih,i1,jh,j1,k1,lwarmstart,ifnamopt,fname_options
    use modmpi,    only : myid,my_real,comm3d,mpi_integer,mpi_logical
    implicit none
    integer :: ierr
    ! #sb3 START - namelist for setting of bulk microphysics
    ! set some initial values before loading namelist
    !
    namelist/NAMBULK3/  &
     l_sb_classic, l_sb_dumpall                              &
     ,l_sb_all_or, l_sb_dbg                                  &
     ,l_setclouds, l_setccn                                  & ! flag whether to set cloud
     ,l_corr_neg_qt                                          & ! flag whether to adjust qt and thlp in hydrometeor corrections
     ,l_sb_lim_aggr, l_sb_stickyice                          & ! ice and snow aggregation flags
     ,l_sb_conv_par                                          & ! conversion flags
     ,l_c_ccn,l_sb_sat_max                                   & ! flags for cloud nucleation
     ,l_sb_nuc_sat,l_sb_nuc_expl,l_sb_nuc_diff               & ! flags for cloud nucleation
     ,l_sb_inuc_sat,l_sb_inuc_expl, l_sb_reisner             & ! flags for ice nucleation
     ,N_inuc, n_i_max, tmp_inuc, x_inuc                      & !  parameters for ice nucleation
     ,N_inuc_R, c_inuc_R, a1_inuc_R, a2_inuc_R               & !  parameters for ice nucleation - Reisner correction
     ,c_ccn, n_clmax                                         & ! C_CCN parameter, used when l_c_ccn
     ,kappa_ccn, x_cnuc,sat_max                              & ! parameters for liquid cloud nucleation
     ,Nc0, xc0_min, Nccn0                                      !! setting of initial clouds


    if(myid==0) then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMBULK3,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMBULK3 '
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMBULK3 '
      endif
      write(6 ,NAMBULK3)
      ! check values
      if (xc0_min.GE.xc_bmax) then
        write(6,*)  'Warning: xc0_min is invalid'
        write(6,*)  '  xc0_min= ',xc0_min,' is larger than xc_bmax =',xc_bmax
        write(6,*)  '  initial min size of droplets cannot be larger than the max size'
        write(6,*)  '  setting xc0_min to default value xcmin =',xcmin
        xc0_min = xcmin
      endif
      if (xc0_min .LE. 0.0) then
        write(6,*)  'Warning: xc0_min is invalid'
        write(6,*)  '  xc0_min= ',xc0_min,' is negative '
        write(6,*)  '  setting xc0_min to default value xcmin =',xcmin
        xc0_min = xcmin
      endif
      if (Nc0 .LE. 0.0) then
        write(6,*)  'Warning: Nc0 is invalid'
        write(6,*)  '  Nc0= ',Nc0,' is negative '
        write(6,*)  '  setting Nc0 to default value Nc_0 = ',Nc_0
        Nc0 = Nc_0
      endif
      ! close the namelist
      close(ifnamopt)
    end if

     ! #sb3 START  - checking if cloud initialisation should be done
    if(myid.eq.0) then
     if(lwarmstart) then
      l_clouds_init     = .true. ! in warm start, clouds are already there
      l_ccn_init        = .true. ! in warm start, ccn are already there
     else
      if(l_setclouds) then   ! if namelist says to initialise clouds
        l_clouds_init     = .false.
      else
        l_clouds_init     = .true.
      endif
      if(l_setccn)    then   ! if namelist says to initialise constant CCN
        if(l_c_ccn)   then
          write(6,*) 'modbulkmicro3: l_c_ccn = .TRUE., l_setccn ignored'
          l_ccn_init        = .true.
        else
          l_ccn_init        = .false.
        endif
      else
        l_ccn_init        = .true.
      endif
     endif
    endif
     ! #sb3 END

    ! send values
     call MPI_BCAST(l_sb_classic,      1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_dumpall,      1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_all_or,       1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_dbg,          1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_clouds_init,     1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_ccn_init,        1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_corr_neg_qt,     1, MPI_LOGICAL ,0,comm3d,ierr)
     !
     call MPI_BCAST(l_sb_lim_aggr,     1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_stickyice,    1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_conv_par ,    1, MPI_LOGICAL ,0,comm3d,ierr)
     !
     call MPI_BCAST(l_c_ccn,           1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_sat_max,      1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_nuc_sat,      1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_nuc_expl,     1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_nuc_diff,     1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_inuc_sat,     1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_inuc_expl,    1, MPI_LOGICAL ,0,comm3d,ierr)
     call MPI_BCAST(l_sb_reisner,      1, MPI_LOGICAL ,0,comm3d,ierr)
     !
     call MPI_BCAST(l_sb_reisner,      1, MPI_LOGICAL ,0,comm3d,ierr)
     !
     call MPI_BCAST(N_inuc_R ,         1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(c_inuc_R ,         1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(a1_inuc_R ,        1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(a2_inuc_R ,        1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(n_i_max ,          1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(N_inuc ,           1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST( tmp_inuc,         1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(x_inuc ,           1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(c_ccn,             1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(n_clmax,           1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(kappa_ccn ,        1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(sat_max ,          1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(x_cnuc ,           1, MY_REAL     ,0,comm3d,ierr)
     !
     call MPI_BCAST(Nc0,               1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(xc0_min,           1, MY_REAL     ,0,comm3d,ierr)
     call MPI_BCAST(Nccn0,             1, MY_REAL     ,0,comm3d,ierr)



    allocate( qltot    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,nuc      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,thlpmcr  (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,qtpmcr   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,sedc     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,sed_qr   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,sed_Nr   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Dvc      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,xc       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Dvr      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,mur      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,lbdr     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,phi      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,tau      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,qr_spl   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Nr_spl   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,wfall_qr (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,wfall_Nr (2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate( n_cc      (2-ih:i1+ih,2-jh:j1+jh,k1)  &  ! N_{ccn} number content [ kg^{-1}] of cloud condensation nuclei
             ,n_ccp     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,n_cl      (2-ih:i1+ih,2-jh:j1+jh,k1)  &  ! N_{c,l}  number content [ kg^{-1}] for liquid cloud droplets,
             ,n_clp     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,n_ci      (2-ih:i1+ih,2-jh:j1+jh,k1)  &  ! N_{c,i}  number content [ kg^{-1}] for ice cloud droplets,
             ,n_cip     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,n_hr      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! N_{h,r} number content [ kg^{-1}] for rain
             ,n_hrp     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,n_hs      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! N_{h,s} number content [ kg^{-1}] for snow
             ,n_hsp     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,n_hg      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! N_{h,g} number content [ kg^{-1}] for graupel
             ,n_hgp     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,q_cl      (2-ih:i1+ih,2-jh:j1+jh,k1)  &  ! q_{c,l}  water content [kg kg^{-1}] for liquid cloud droplets,
             ,q_clp     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,q_ci      (2-ih:i1+ih,2-jh:j1+jh,k1)  &  ! q_{c,i}  water content [kg kg^{-1}] for ice cloud droplets,
             ,q_cip     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,q_hr      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! q_{h,r} water content [kg kg^{-1}] for rain
             ,q_hrp     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,q_hs      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! q_{h,s} water content [kg kg^{-1}] for snow
             ,q_hsp     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,q_hg      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! q_{h,g} water content [kg kg^{-1}] for graupel
             ,q_hgp     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,x_ci      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! mean cloud ice size
             ,x_cl      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! mean cloud water size
             ,x_hs      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! mean snow size
             ,x_hg      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! mean graupel size
             ,x_hr      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! mean raindrop size
            )

    allocate( precep_l    (2-ih:i1+ih,2-jh:j1+jh,k1)   &
             ,precep_i    (2-ih:i1+ih,2-jh:j1+jh,k1)   &
            )
    allocate( qcmask    (2-ih:i1+ih,2-jh:j1+jh,k1) )
    allocate( q_ci_mask    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,q_cl_mask    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,q_hr_mask    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,q_hs_mask    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,q_hg_mask    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
            )

  gamma25=lacz_gamma(2.5)
  gamma3=2.
  gamma35=lacz_gamma(3.5)

  ! adding calculation of the constant part for moment
  c_mmt_1cl = calc_cons_mmt (1, mu_cl_cst, nu_cl_cst)
  c_mmt_1hr = calc_cons_mmt (1, mu_hr_cst, nu_hr_cst)
  c_mmt_2cl = calc_cons_mmt (2, mu_cl_cst, nu_cl_cst)
  c_mmt_2hr = calc_cons_mmt (2, mu_hr_cst, nu_hr_cst)

  ! adding calculation for constants in terminal velocity calculation
  c_v_c0 = calc_cons_v (0, mu_cl_cst, nu_cl_cst, al_cl, be_cl)
  c_v_i0 = calc_cons_v (0, mu_ci_cst, nu_ci_cst, al_ci, be_ci)
  c_v_s0 = calc_cons_v (0, mu_hs_cst, nu_hs_cst, al_hs, be_hs)
  c_v_g0 = calc_cons_v (0, mu_hg_cst, nu_hg_cst, al_hg, be_hg)

  c_v_c1 = calc_cons_v (1, mu_cl_cst, nu_cl_cst, al_cl, be_cl)
  c_v_i1 = calc_cons_v (1, mu_ci_cst, nu_ci_cst, al_ci, be_ci)
  c_v_s1 = calc_cons_v (1, mu_hs_cst, nu_hs_cst, al_hs, be_hs)
  c_v_g1 = calc_cons_v (1, mu_hg_cst, nu_hg_cst, al_hg, be_hg)

  ! adding ventilation coefficients
  aven_0r = calc_avent(0, mu_hr_cst,nu_hr_cst, a_hr, b_hr, a_v_r)
  aven_0i = calc_avent(0, mu_ci_cst,nu_ci_cst, a_ci, b_ci, a_v_i)
  aven_0s = calc_avent(0, mu_hs_cst,nu_hs_cst, a_hs, b_hs, a_v_s)
  aven_0g = calc_avent(0, mu_hg_cst,nu_hg_cst, a_hg, b_hg, a_v_g)

  bven_0r = calc_bvent(0, mu_hr_cst,nu_hr_cst, a_hr, b_hr, be_hr, b_v_r)
  bven_0i = calc_bvent(0, mu_ci_cst,nu_ci_cst, a_ci, b_ci, be_ci, b_v_i)
  bven_0s = calc_bvent(0, mu_hs_cst,nu_hs_cst, a_hs, b_hs, be_hs, b_v_s)
  bven_0g = calc_bvent(0, mu_hg_cst,nu_hg_cst, a_hg, b_hg, be_hg, b_v_g)

  aven_1r = calc_avent(1, mu_hr_cst,nu_hr_cst, a_hr, b_hr, a_v_r)
  aven_1i = calc_avent(1, mu_ci_cst,nu_ci_cst, a_ci, b_ci, a_v_i)
  aven_1s = calc_avent(1, mu_hs_cst,nu_hs_cst, a_hs, b_hs, a_v_s)
  aven_1g = calc_avent(1, mu_hg_cst,nu_hg_cst, a_hg, b_hg, a_v_g)

  bven_1r = calc_bvent(1, mu_hr_cst,nu_hr_cst, a_hr, b_hr, be_hr, b_v_r)
  bven_1i = calc_bvent(1, mu_ci_cst,nu_ci_cst, a_ci, b_ci, be_ci, b_v_i)
  bven_1s = calc_bvent(1, mu_hs_cst,nu_hs_cst, a_hs, b_hs, be_hs, b_v_s)
  bven_1g = calc_bvent(1, mu_hg_cst,nu_hg_cst, a_hg, b_hg, be_hg, b_v_g)

  ! adding collision coefficients

  ! coefficient for collected particle b
  !   b  \in  {l, r, i, s, g}
  dlt_c0   = calc_delta_b (0, mu_cl_cst, nu_cl_cst, b_cl)
  dlt_c1   = calc_delta_b (1, mu_cl_cst, nu_cl_cst, b_cl)
  th_c0    = calc_th_b    (0, mu_cl_cst, nu_cl_cst, b_cl, be_cl)
  th_c1    = calc_th_b    (1, mu_cl_cst, nu_cl_cst, b_cl, be_cl)

  dlt_r0   = calc_delta_b (0, mu_hr_cst, nu_hr_cst, b_hr)
  dlt_r1   = calc_delta_b (1, mu_hr_cst, nu_hr_cst, b_hr)
  th_r0    = calc_th_b    (0, mu_hr_cst, nu_hr_cst, b_hr, be_hr)
  th_r1    = calc_th_b    (1, mu_hr_cst, nu_hr_cst, b_hr, be_hr)

  dlt_i0   = calc_delta_b (0, mu_ci_cst, nu_ci_cst, b_ci)
  dlt_i1   = calc_delta_b (1, mu_ci_cst, nu_ci_cst, b_ci)
  th_i0    = calc_th_b    (0, mu_ci_cst, nu_ci_cst, b_ci, be_ci)
  th_i1    = calc_th_b    (1, mu_ci_cst, nu_ci_cst, b_ci, be_ci)

  dlt_s0   = calc_delta_b (0, mu_hs_cst, nu_hs_cst, b_hs)
  dlt_s1   = calc_delta_b (1, mu_hs_cst, nu_hs_cst, b_hs)
  th_s0    = calc_th_b    (0, mu_hs_cst, nu_hs_cst, b_hs, be_hs)
  th_s1    = calc_th_b    (1, mu_hs_cst, nu_hs_cst, b_hs, be_hs)

  dlt_g0   = calc_delta_b (0, mu_hg_cst, nu_hg_cst, b_hg)
  dlt_g1   = calc_delta_b (1, mu_hg_cst, nu_hg_cst, b_hg)
  th_g0    = calc_th_b    (0, mu_hg_cst, nu_hg_cst, b_hg, be_hg)
  th_g1    = calc_th_b    (1, mu_hg_cst, nu_hg_cst, b_hg, be_hg)

  ! and for particle couple {a,b}
  !   a in {i,s,g}
  !   b in {l, r, i, s, g}

  dlt_i0c  = calc_delta_ab (0, mu_ci_cst, nu_ci_cst, b_ci, mu_cl_cst, nu_cl_cst, b_cl)
  dlt_i0r  = calc_delta_ab (0, mu_ci_cst, nu_ci_cst, b_ci, mu_hr_cst, nu_hr_cst, b_hr)
  dlt_i0i  = calc_delta_ab (0, mu_ci_cst, nu_ci_cst, b_ci, mu_ci_cst, nu_ci_cst, b_ci)

  dlt_r0i  = calc_delta_ab (0, mu_hr_cst, nu_hr_cst, b_hr, mu_ci_cst, nu_ci_cst, b_ci)
  dlt_r0s  = calc_delta_ab (0, mu_hr_cst, nu_hr_cst, b_hr, mu_hs_cst, nu_hs_cst, b_hs)
  dlt_r0g  = calc_delta_ab (0, mu_hr_cst, nu_hr_cst, b_hr, mu_hg_cst, nu_hg_cst, b_hg)

  dlt_s0c  = calc_delta_ab (0, mu_hs_cst, nu_hs_cst, b_hs, mu_cl_cst, nu_cl_cst, b_cl)
  dlt_s0r  = calc_delta_ab (0, mu_hs_cst, nu_hs_cst, b_hs, mu_hr_cst, nu_hr_cst, b_hr)
  dlt_s0i  = calc_delta_ab (0, mu_hs_cst, nu_hs_cst, b_hs, mu_ci_cst, nu_ci_cst, b_ci)
  dlt_s0s  = calc_delta_ab (0, mu_hs_cst, nu_hs_cst, b_hs, mu_hs_cst, nu_hs_cst, b_hs)

  dlt_g0c  = calc_delta_ab (0, mu_hg_cst, nu_hg_cst, b_hg, mu_cl_cst, nu_cl_cst, b_cl)
  dlt_g0r  = calc_delta_ab (0, mu_hg_cst, nu_hg_cst, b_hg, mu_hr_cst, nu_hr_cst, b_hr)
  dlt_g0i  = calc_delta_ab (0, mu_hg_cst, nu_hg_cst, b_hg, mu_ci_cst, nu_ci_cst, b_ci)
  dlt_g0s  = calc_delta_ab (0, mu_hg_cst, nu_hg_cst, b_hg, mu_hs_cst, nu_hs_cst, b_hs)
  dlt_g0g  = calc_delta_ab (0, mu_hg_cst, nu_hg_cst, b_hg, mu_hg_cst, nu_hg_cst, b_hg)


  th_i0c   = calc_th_ab (0, mu_ci_cst, nu_ci_cst, b_ci, be_ci, mu_cl_cst, nu_cl_cst, b_cl, be_cl)
  th_i0r   = calc_th_ab (0, mu_ci_cst, nu_ci_cst, b_ci, be_ci, mu_hr_cst, nu_hr_cst, b_hr, be_hr)
  th_i0i   = calc_th_ab (0, mu_ci_cst, nu_ci_cst, b_ci, be_ci, mu_ci_cst, nu_ci_cst, b_ci, be_ci)

  th_r0i   = calc_th_ab (0, mu_hr_cst, nu_hr_cst, b_hr, be_hr, mu_ci_cst, nu_ci_cst, b_ci, be_ci)
  th_r0s   = calc_th_ab (0, mu_hr_cst, nu_hr_cst, b_hr, be_hr, mu_hs_cst, nu_hs_cst, b_hs, be_hs)
  th_r0g   = calc_th_ab (0, mu_hr_cst, nu_hr_cst, b_hr, be_hr, mu_hg_cst, nu_hg_cst, b_hg, be_hg)

  th_s0c   = calc_th_ab (0, mu_hs_cst, nu_hs_cst, b_hs, be_hs, mu_cl_cst, nu_cl_cst, b_cl, be_cl)
  th_s0r   = calc_th_ab (0, mu_hs_cst, nu_hs_cst, b_hs, be_hs, mu_hr_cst, nu_hr_cst, b_hr, be_hr)
  th_s0i   = calc_th_ab (0, mu_hs_cst, nu_hs_cst, b_hs, be_hs, mu_ci_cst, nu_ci_cst, b_ci, be_ci)
  th_s0s   = calc_th_ab (0, mu_hs_cst, nu_hs_cst, b_hs, be_hs, mu_hs_cst, nu_hs_cst, b_hs, be_hs)

  th_g0c   = calc_th_ab (0, mu_hg_cst, nu_hg_cst, b_hg, be_hg, mu_cl_cst, nu_cl_cst, b_cl, be_cl)
  th_g0r   = calc_th_ab (0, mu_hg_cst, nu_hg_cst, b_hg, be_hg, mu_hr_cst, nu_hr_cst, b_hr, be_hr)
  th_g0i   = calc_th_ab (0, mu_hg_cst, nu_hg_cst, b_hg, be_hg, mu_ci_cst, nu_ci_cst, b_ci, be_ci)
  th_g0s   = calc_th_ab (0, mu_hg_cst, nu_hg_cst, b_hg, be_hg, mu_hs_cst, nu_hs_cst, b_hs, be_hs)
  th_g0g   = calc_th_ab (0, mu_hg_cst, nu_hg_cst, b_hg, be_hg, mu_hg_cst, nu_hg_cst, b_hg, be_hg)

  dlt_i1c  = calc_delta_ab (1, mu_ci_cst, nu_ci_cst, b_ci, mu_cl_cst, nu_cl_cst, b_cl)
  dlt_i1r  = calc_delta_ab (1, mu_ci_cst, nu_ci_cst, b_ci, mu_hr_cst, nu_hr_cst, b_hr)
  dlt_i1i  = calc_delta_ab (1, mu_ci_cst, nu_ci_cst, b_ci, mu_ci_cst, nu_ci_cst, b_ci)

  dlt_r1i  = calc_delta_ab (1, mu_hr_cst, nu_hr_cst, b_hr, mu_ci_cst, nu_ci_cst, b_ci)
  dlt_r1s  = calc_delta_ab (1, mu_hr_cst, nu_hr_cst, b_hr, mu_hs_cst, nu_hs_cst, b_hs)
  dlt_r1g  = calc_delta_ab (1, mu_hr_cst, nu_hr_cst, b_hr, mu_hg_cst, nu_hg_cst, b_hg)

  dlt_s1c  = calc_delta_ab (1, mu_hs_cst, nu_hs_cst, b_hs, mu_cl_cst, nu_cl_cst, b_cl)
  dlt_s1r  = calc_delta_ab (1, mu_hs_cst, nu_hs_cst, b_hs, mu_hr_cst, nu_hr_cst, b_hr)
  dlt_s1i  = calc_delta_ab (1, mu_hs_cst, nu_hs_cst, b_hs, mu_ci_cst, nu_ci_cst, b_ci)
  dlt_s1s  = calc_delta_ab (1, mu_hs_cst, nu_hs_cst, b_hs, mu_hs_cst, nu_hs_cst, b_hs)

  dlt_g1c  = calc_delta_ab (1, mu_hg_cst, nu_hg_cst, b_hg, mu_cl_cst, nu_cl_cst, b_cl)
  dlt_g1r  = calc_delta_ab (1, mu_hg_cst, nu_hg_cst, b_hg, mu_hr_cst, nu_hr_cst, b_hr)
  dlt_g1i  = calc_delta_ab (1, mu_hg_cst, nu_hg_cst, b_hg, mu_ci_cst, nu_ci_cst, b_ci)
  dlt_g1s  = calc_delta_ab (1, mu_hg_cst, nu_hg_cst, b_hg, mu_hs_cst, nu_hs_cst, b_hs)
  dlt_g1g  = calc_delta_ab (1, mu_hg_cst, nu_hg_cst, b_hg, mu_hg_cst, nu_hg_cst, b_hg)


  th_i1c   = calc_th_ab (1, mu_ci_cst, nu_ci_cst, b_ci, be_ci, mu_cl_cst, nu_cl_cst, b_cl, be_cl)
  th_i1r   = calc_th_ab (1, mu_ci_cst, nu_ci_cst, b_ci, be_ci, mu_hr_cst, nu_hr_cst, b_hr, be_hr)
  th_i1i   = calc_th_ab (1, mu_ci_cst, nu_ci_cst, b_ci, be_ci, mu_ci_cst, nu_ci_cst, b_ci, be_ci)

  th_r1i   = calc_th_ab (1, mu_hr_cst, nu_hr_cst, b_hr, be_hr, mu_ci_cst, nu_ci_cst, b_ci, be_ci)
  th_r1s   = calc_th_ab (1, mu_hr_cst, nu_hr_cst, b_hr, be_hr, mu_hs_cst, nu_hs_cst, b_hs, be_hs)
  th_r1g   = calc_th_ab (1, mu_hr_cst, nu_hr_cst, b_hr, be_hr, mu_hg_cst, nu_hg_cst, b_hg, be_hg)

  th_s1c   = calc_th_ab (1, mu_hs_cst, nu_hs_cst, b_hs, be_hs, mu_cl_cst, nu_cl_cst, b_cl, be_cl)
  th_s1r   = calc_th_ab (1, mu_hs_cst, nu_hs_cst, b_hs, be_hs, mu_hr_cst, nu_hr_cst, b_hr, be_hr)
  th_s1i   = calc_th_ab (1, mu_hs_cst, nu_hs_cst, b_hs, be_hs, mu_ci_cst, nu_ci_cst, b_ci, be_ci)
  th_s1s   = calc_th_ab (1, mu_hs_cst, nu_hs_cst, b_hs, be_hs, mu_hs_cst, nu_hs_cst, b_hs, be_hs)

  th_g1c   = calc_th_ab (1, mu_hg_cst, nu_hg_cst, b_hg, be_hg, mu_cl_cst, nu_cl_cst, b_cl, be_cl)
  th_g1r   = calc_th_ab (1, mu_hg_cst, nu_hg_cst, b_hg, be_hg, mu_hr_cst, nu_hr_cst, b_hr, be_hr)
  th_g1i   = calc_th_ab (1, mu_hg_cst, nu_hg_cst, b_hg, be_hg, mu_ci_cst, nu_ci_cst, b_ci, be_ci)
  th_g1s   = calc_th_ab (1, mu_hg_cst, nu_hg_cst, b_hg, be_hg, mu_hs_cst, nu_hs_cst, b_hs, be_hs)
  th_g1g   = calc_th_ab (1, mu_hg_cst, nu_hg_cst, b_hg, be_hg, mu_hg_cst, nu_hg_cst, b_hg, be_hg)

  ! calculate max sedim velocities
  !  both for
  wfallmax_hr = d_wfallmax_hr   ! use default value for rain
  wfallmax_cl = d_wfallmax_hr   ! wfallmax_cl =   max(c_v_c0*x_cl_bmax**be_cl, c_v_c1*x_cl_bmax**be_cl)
  wfallmax_ci = d_wfallmax_hr   ! wfallmax_ci =   max( c_v_i0*x_ci_bmax**be_ci, c_v_i1*x_cl_bmax**be_ci)
  wfallmax_hs = d_wfallmax_hr   ! wfallmax_hs =   max( c_v_s0*x_hs_bmax**be_hs, c_v_s1*x_hs_bmax**be_hs)
  !#b2t3 wfallmax_hg = d_wfallmax_hg   ! wfallmax_hg =   max( c_v_g0*x_hg_bmax**be_hg, c_v_g1*x_hg_bmax**be_hg)
  wfallmax_hg = max(c_v_g0, c_v_g1) ! #b2t4

  ! coefficient for lbdr calculation l_sb_classic
   c_lbdr = calc_cons_lbd (mu_hr_cst, nu_hs_cst)

  ! calculating eslt3 = esl(T_3)/T_3
  ! using thermodynamic routine
  eslt3 = esl_at_te(T_3)

  ! #sb3 END

  ! #sb3 START write settings
   if (myid .eq. 0) then
     write (6,*) ' === Settings of Mphys ===='
     write (6,*) 'imicro = ', imicro
     write (6,*) 'l_sb = ', l_sb
     write (6,*) 'l_sb_classic = ', l_sb_classic
     write (6,*) 'l_sb_all_or = ', l_sb_all_or
     write (6,*) 'l_sb_dumpall = ', l_sb_dumpall
     write (6,*) 'l_sb_lim_aggr = ',  l_sb_lim_aggr
     write (6,*) 'l_sedc = ', l_sedc
     write (6,*) 'l_rain = ', l_rain
     write (6,*) 'l_mur_cst  = ', l_mur_cst

     write (6,*)   ' --------------------------------------------------- '
     write (6,*)   '  '
     write (6,*)   ' some calculated Mphys integrals'
     write (6,*)   '    2*dlt_s0 + dlt_s0s = ',        2*dlt_s0 + dlt_s0s
     write (6,*)   '    2*dlt_i0 + dlt_i0i = ',        2*dlt_i0 + dlt_i0i

     write (6,*)   '    dlt_s0 + dlt_s1s+ dlt_s1 = ' , dlt_s0 + dlt_s1s+ dlt_s1
     write (6,*)   '    dlt_i0 + dlt_i1i + dlt_i1 = ', dlt_i0 + dlt_i1i + dlt_i1

     write (6,*)   '    2*th_s0 - th_s0s  = ', 2*th_s0 - th_s0s    ! from Seifert, 2002
     write (6,*)   '    2*th_i0 - th_i0i  = ', 2*th_i0 - th_i0i    ! from Seifert, 2002

     write (6,*)   '    th_s0 - th_s1s + th_s1', th_s0 - th_s1s + th_s1
     write (6,*)   '    th_i0 - th_i1i + th_i1', th_i0 - th_i1i + th_i1  ! th_i0i - th_i1i + th_i1 ! th_1ab    = th_i1i

     write (6,*)   '  '
     write (6,*)   ' some ventilation parameters'
     write (6,*)   '   aven_0r=', aven_0r,'aven_1r=',aven_1r,'bven_0r=', bven_0r,'bven_1r=',bven_1r
     write (6,*)   '   aven_0i=', aven_0i,'aven_1i=',aven_1i,'bven_0i=', bven_0i,'bven_1i=',bven_1i
     write (6,*)   '   aven_0s=', aven_0s,'aven_1s=',aven_1s,'bven_0s=', bven_0s,'bven_1s=',bven_1s
     write (6,*)   '   aven_0g=', aven_0g,'aven_1g=',aven_1g,'bven_0g=', bven_0g,'bven_1g=',bven_1g
     write (6,*)   '  '
     write (6,*)   ' some calculated fall coefficients'
     write (6,*)   '   c_v_c0=', c_v_c0,' c_v_c1=', c_v_c1
     write (6,*)   '   c_v_i0=', c_v_i0,' c_v_i1=', c_v_i1
     write (6,*)   '   c_v_s0=', c_v_s0,' c_v_s1=', c_v_s1
     write (6,*)   '   c_v_g0=', c_v_g0,' c_v_g1=', c_v_g1
     write (6,*)   '  '
     write (6,*)   '   c_mmt_1cl=',c_mmt_1cl,'c_mmt_2cl=',c_mmt_2cl
     write (6,*)   '   c_mmt_1hr=',c_mmt_1hr,'c_mmt_2hr=',c_mmt_2hr
     write (6,*)   '  '
     write (6,*)   ' therefore max sedim velocities '
     write (6,*)   '   w0_cl_max=', c_v_c0*x_cl_bmax**be_cl,' w1_cl_max=', c_v_c1*x_cl_bmax**be_cl
     write (6,*)   '   w0_ci_max=', c_v_i0*x_ci_bmax**be_ci,' w1_ci_max=', c_v_i1*x_cl_bmax**be_ci
     write (6,*)   '   w0_hs_max=', c_v_s0*x_hs_bmax**be_hs,' w1_hs_max=', c_v_s1*x_hs_bmax**be_hs
     write (6,*)   '   w0_hg_max=', c_v_g0*x_hg_bmax**be_hg,' w1_hg_max=', c_v_g1*x_hg_bmax**be_hg

     write (6,*)   ' --------------------------------------------------- '

     write (6,*)   ' ice and snow aggregation debug '
     write (6,*)   '   dlt_i0 =',dlt_i0
     write (6,*)   '   dlt_0aa   = 2*dlt_i0 + dlt_i0i = ',2*dlt_i0 + dlt_i0i
      ! dlt_0b    = dlt_0a
     write (6,*)   '   dlt_1a    = dlt_i1 = ',dlt_i1
     write (6,*)   '   dlt_1aa   = dlt_i0 + dlt_i1i + dlt_i1 = ',dlt_i0 + dlt_i1i + dlt_i1
      ! dlt_1b    = dlt_1a
     write (6,*)   '   th_0a     = th_i0 = ',th_i0
     write (6,*)   '   th_0aa    = 2*th_i0 - th_i0i = ', 2*th_i0 - th_i0i ! from Seifert, 2002
      ! th_0b     = th_0a
     write (6,*)   '  th_1a     = th_i1 = ', th_i1
     write (6,*)   '  th_1aa    = th_i0 - th_i1i + th_i1 = ', th_i0 - th_i1i + th_i1  ! th_i0i - th_i1i + th_i1 ! th_1ab    = th_i1i
      write (6,*)   ' --------------------------------------------------- '
      ! th_1b     = th_1a

   endif
  end subroutine initbulkmicro3

!> Cleaning up after the run
  subroutine exitbulkmicro3
  !*********************************************************************
  ! subroutine exitbulkmicro
  !*********************************************************************
    implicit none


    deallocate(sedc,sed_qr,sed_Nr,Dvc,xc,Dvr,mur,lbdr      &
               ,phi,tau,qr_spl,Nr_spl,wfall_qr,wfall_Nr)
    deallocate(qltot,nuc,thlpmcr,qtpmcr)
    deallocate(n_cc,n_ccp                                     &
      ,n_cl,n_clp,n_ci,n_cip,n_hr,n_hrp,n_hs,n_hsp,n_hg,n_hgp &
      ,q_cl,q_clp,q_ci,q_cip,q_hr,q_hrp,q_hs,q_hsp,q_hg,q_hgp &
      ,x_cl, x_ci, x_hr, x_hs, x_hg                           &
      )
    deallocate(qcmask,q_ci_mask,q_cl_mask, q_hr_mask, q_hs_mask, q_hg_mask )
    deallocate(precep_l, precep_i)

  end subroutine exitbulkmicro3

!> Calculates the microphysical source term.
  subroutine bulkmicro3
    use modglobal, only : i1,j1,k1,rdt,rk3step,timee,rlv,cp, ih, jh

    use modfields, only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof,qvsl,qvsi,qt0, thl0
    use modbulkmicrostat3, only : bulkmicrotend3,bulkmicrostat3 ! #sb3 #t7
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k
    integer :: insv, iqsv ! #sb3 - species number
    logical :: flag_warn_update, flag_warn_test, flag_warn_size,flag_nan   &
              ,flag_do_dbg,flag_do_dbg_up
    real :: qrtest,nrtest
    real :: qicew, qcliq ! , cogr_max, ql_res
    real :: temp_thlp,temp_thlp1
    real, allocatable :: xtest (:,:,:)


    allocate( xtest(2-ih:i1+ih,2-jh:j1+jh,k1))  ! #d
    xtest = 0.0

    allocate( dn_cl_nu      (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< droplet nucleation rate
             ,dn_ci_inu     (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< ice nucleation rate
             ,dn_cl_au      (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< change in number of cloud droplets due to autoconversion
             ,dq_hr_au      (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< change in mass of raindrops due to autoconversion
             ,dn_hr_au      (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< change in number of raindrops due to autoconversion
             ,dq_hr_ac      (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< change in mass of raindrops due to accretion
             ,dn_cl_ac      (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< change in number of cloud droplets due to accretion
             ,dn_hr_br      (2-ih:i1+ih,2-jh:j1+jh,k1) &
             ,dn_hr_sc      (2-ih:i1+ih,2-jh:j1+jh,k1) &
             ,dq_hr_ev      (2-ih:i1+ih,2-jh:j1+jh,k1) &
             ,dn_hr_ev      (2-ih:i1+ih,2-jh:j1+jh,k1) &
             ,dq_ci_rime    (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< riming growth of ice
             ,dn_cl_rime_ci (2-ih:i1+ih,2-jh:j1+jh,k1) &    !<  - and impact on n_cl
             ,dq_hs_rime    (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< riming growth of snow
             ,dn_cl_rime_hs (2-ih:i1+ih,2-jh:j1+jh,k1) &    !<  - and impact on n_cl
             ,dq_hg_rime    (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< riming growth for graupel
             ,dn_cl_rime_hg (2-ih:i1+ih,2-jh:j1+jh,k1) &    !<  - and impact on n_cl
             ,dq_hghr_rime  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< riming growth for graupel with rain
             ,dn_hr_rime_hg (2-ih:i1+ih,2-jh:j1+jh,k1) &    !<  - and impact on n_hr
             ,dq_hshr_rime  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< riming growth for snow with rain
             ,dn_hr_rime_hs (2-ih:i1+ih,2-jh:j1+jh,k1) &    !<  - and impact on n_hr
             ,dq_hr_rime_ri (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< rain loss from riming of ice+rain->gr
             ,dq_ci_rime_ri (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< ice loss from riming of ice+rain->gr
             ,dn_ci_rime_ri (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< ice number loss from riming of ice+rain->gr
             ,dn_hr_rime_ri (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< rain number loss from riming of ice+rain->gr
             ,dq_hr_col_rs  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< rain loss from riming of ice+snow->gr
             ,dq_hs_col_rs  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< rain number loss from riming of ice+snow->gr
             ,dn_hr_col_rs  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< snow loss from riming of ice+snow->gr
             ,dn_hs_col_rs  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< snow number loss from riming of ice+snow->gr
             ,dq_hr_col_ri  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< rain loss from riming of ice+rain->gr
             ,dq_ci_col_ri  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< ice loss from riming of ice+rain->gr
             ,dn_ci_col_ri  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< ice number loss from riming of ice+rain->gr
             ,dn_hr_col_ri  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< rain number loss from riming of ice+rain->gr
             ,dq_cl_het     (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< heterogeneou freezing of cloud water
             ,dq_hr_het     (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< heterogeneou freezing of raindrops
             ,dn_hr_het     (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< heterogeneou freezing of raindrops
             ,dq_cl_hom     (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< homogeneous freezing of cloud water
             ,dq_ci_col_iis (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< self-collection of cloud ice
             ,dn_ci_col_iis (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< self-collection of cloud ice
             ,dn_hs_col_sss (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< self-collection of snow
             ,dq_hsci_col   (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< collection s+i - trend in q_hs
             ,dn_ci_col_hs  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< collection s+i - trend in n_ci
             ,dq_hghs_col   (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< collection g+s - trend in q_hg
             ,dn_hs_col_hg  (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< collection g+s - trend in n_hs
             ,dq_ci_cv      (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< partial conversion ice -> graupel
             ,dn_ci_cv      (2-ih:i1+ih,2-jh:j1+jh,k1) &
             ,dq_hs_cv      (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< partial conversion snow-> graupel
             ,dn_hs_cv      (2-ih:i1+ih,2-jh:j1+jh,k1) &
             ,dn_cl_sc      (2-ih:i1+ih,2-jh:j1+jh,k1) &    !< cloud self-collection
             ,dn_ci_mul     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< ice multiplication
             ,dq_ci_mul     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< ice multiplication
             ,dn_ci_me      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency melting of cloud ice
             ,dq_ci_me      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency melting of cloud ice
             ,dn_hs_me      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency melting of snow
             ,dq_hs_me      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency melting of snow
             ,dn_hg_me      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency melting of graupel
             ,dq_hg_me      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency melting of graupel
             ,dn_ci_ev      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency evaporation of cloud ice
             ,dq_ci_ev      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency evaporation of cloud ice
             ,dn_hs_ev      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency evaporation of snow
             ,dq_hs_ev      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency evaporation of snow
             ,dn_hg_ev      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency evaporation of graupel
             ,dq_hg_ev      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency evaporation of graupel
             ,dn_ci_eme_ic  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency enhanced melting of cloud ice by cloud water
             ,dq_ci_eme_ic  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency enhanced melting of cloud ice by cloud water
             ,dn_ci_eme_ri  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency enhanced melting of cloud ice by rain
             ,dq_ci_eme_ri  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency enhanced melting of cloud ice  by rain
             ,dn_hs_eme_sc  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency enhanced melting of snow by cloud water
             ,dq_hs_eme_sc  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency enhanced melting of snow by cloud water
             ,dn_hs_eme_rs  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency enhanced melting of snow by rain
             ,dq_hs_eme_rs  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency enhanced melting of snow by rain
             ,dn_hg_eme_gc  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency enhanced melting of graupel
             ,dq_hg_eme_gc  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency enhanced melting of graupel
             ,dn_hg_eme_gr  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< number tendency enhanced melting of graupel by rain
             ,dq_hg_eme_gr  (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< mass tendency enhanced melting of graupel by rain
             ,dn_cl_se      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< sedimentation for clouds water - number
             ,dq_cl_se      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !<       -||-- mixing ration
             ,dn_ci_se      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< sedimentation for cloud ice - number
             ,dq_ci_se      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !<       -||-- mixing ration
             ,dn_hr_se      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< sedimentation for rain - number
             ,dq_hr_se      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !<       -||-- mixing ration
             ,dn_hs_se      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< sedimentation for snow - number
             ,dq_hs_se      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !<       -||-- mixing ration
             ,dn_hg_se      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< sedimentation for graupel - number
             ,dq_hg_se      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !<       -||-- mixing ration
             ,precep_hr     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of raindrops
             ,precep_ci     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of ice crystals
             ,precep_hs     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of snow
             ,precep_hg     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of graupel
             ,dq_cl_sa      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< saturation adjustment
             ,dn_cl_sa      (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< change in n_cl due to saturation adjustment
             ,ret_cc        (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< recovery of ccn
            )

    allocate( D_cl          (2-ih:i1+ih,2-jh:j1+jh,k1) &   ! \bar{D}_cl mean diameter for cloud water particle
             ,D_ci          (2-ih:i1+ih,2-jh:j1+jh,k1) &   ! \bar{D}_ci mean diameter for cloud ice particle
             ,D_hr          (2-ih:i1+ih,2-jh:j1+jh,k1) &   ! \bar{D}_hr mean diameter for raindrops
             ,D_hs          (2-ih:i1+ih,2-jh:j1+jh,k1) &   ! \bar{D}_hs mean diameter for snow particle
             ,D_hg          (2-ih:i1+ih,2-jh:j1+jh,k1) &   ! \bar{D}_hg mean diameter for graupel particle
             ,v_cl          (2-ih:i1+ih,2-jh:j1+jh,k1) &   ! \bar{v}_cl mean velocity for cloud water droplets
             ,v_ci          (2-ih:i1+ih,2-jh:j1+jh,k1) &   ! \bar{v}_ci mean velocity for cloud ice particle
             ,v_hr          (2-ih:i1+ih,2-jh:j1+jh,k1) &   ! \bar{v}_hr mean velocity for raindrops
             ,v_hs          (2-ih:i1+ih,2-jh:j1+jh,k1) &   ! \bar{v}_hs mean velocity for snow particle
             ,v_hg          (2-ih:i1+ih,2-jh:j1+jh,k1) &   ! \bar{v}_hg mean velocity for graupel particle
             )

    ! check if cloud were already initialised
    if (.not. l_ccn_init) then
      call initccn3       ! initialise clouds
      l_ccn_init = .true.
      if(myid.eq.0) then
        write(6,*) ' modbulkmicro3: ccn initialised by initccn3'
      endif
    endif
    if (.not. l_clouds_init) then
      call initclouds3       ! initialise clouds
      l_clouds_init = .true.
      if (myid.eq.0) then
        write(6,*) ' modbulkmicro3: clouds initialised by initcloud3'
      endif
    endif

    ! loading cloud water number and cloud water densities
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! aerosols:
      n_cc  (i,j,k) = sv0(i,j,k,in_cc)

      ! reading number densities:
      n_cl  (i,j,k) = sv0(i,j,k,in_cl)
      n_ci  (i,j,k) = sv0(i,j,k,in_ci)
      n_hr  (i,j,k) = sv0(i,j,k,in_hr)
      n_hs  (i,j,k) = sv0(i,j,k,in_hs)
      n_hg  (i,j,k) = sv0(i,j,k,in_hg)

      ! reading mass densities
      q_cl  (i,j,k) = sv0(i,j,k,iq_cl)
      q_ci  (i,j,k) = sv0(i,j,k,iq_ci)
      q_hr  (i,j,k) = sv0(i,j,k,iq_hr)
      q_hs  (i,j,k) = sv0(i,j,k,iq_hs)
      q_hg  (i,j,k) = sv0(i,j,k,iq_hg)
    enddo
    enddo
    enddo

    ! allocating fields for processes
    thlpmcr = 0.0
    qtpmcr  = 0.0

    ! 0 to values of the update
    n_ccp  = 0.0
    n_clp  = 0.0
    n_cip  = 0.0
    n_hrp  = 0.0
    n_hsp  = 0.0
    n_hgp  = 0.0
    q_hrp  = 0.0
    q_clp  = 0.0
    q_cip  = 0.0
    q_hrp  = 0.0
    q_hsp  = 0.0
    q_hgp  = 0.0
    x_ci   = 0.0
    x_cl   = 0.0
    x_hr   = 0.0
    x_hs   = 0.0
    x_hg   = 0.0
    !
    D_ci   = 0.0
    D_cl   = 0.0
    D_hr   = 0.0
    D_hs   = 0.0
    D_hg   = 0.0
    v_ci   = 0.0
    v_cl   = 0.0
    v_hr   = 0.0
    v_hs   = 0.0
    v_hg   = 0.0

    ! 0 to values of the update
    dn_cl_nu      = 0.0
    dn_ci_inu     = 0.0
    dn_cl_au      = 0.0
    dq_hr_au      = 0.0
    dn_hr_au      = 0.0
    dq_hr_ac      = 0.0
    dn_cl_ac      = 0.0
    dn_hr_br      = 0.0
    dn_hr_sc      = 0.0
    dq_hr_ev      = 0.0
    dn_hr_ev      = 0.0
    dq_ci_dep     = 0.0
    dq_hs_dep     = 0.0
    dq_hg_dep     = 0.0
    dq_ci_rime    = 0.0
    dn_cl_rime_ci = 0.0
    dq_hs_rime    = 0.0
    dn_cl_rime_hs = 0.0
    dq_hg_rime    = 0.0
    dn_cl_rime_hg = 0.0
    dq_hshr_rime  = 0.0
    dn_hr_rime_hs = 0.0
    dq_hr_col_rs  = 0.0
    dq_hs_col_rs  = 0.0
    dn_hr_col_rs  = 0.0
    dn_hs_col_rs  = 0.0
    dq_hr_col_ri  = 0.0
    dq_ci_col_ri  = 0.0
    dn_ci_col_ri  = 0.0
    dn_hr_col_ri  = 0.0
    dq_hghr_rime  = 0.0
    dn_hr_rime_hg = 0.0
    dq_cl_het     = 0.0
    dq_hr_het     = 0.0
    dn_hr_het     = 0.0
    dq_cl_hom     = 0.0
    dq_ci_col_iis = 0.0
    dn_ci_col_iis = 0.0
    dn_hs_col_sss = 0.0
    dq_hsci_col   = 0.0
    dn_ci_col_hs  = 0.0
    dq_hghs_col   = 0.0
    dn_hs_col_hg  = 0.0
    dq_ci_cv      = 0.0
    dn_ci_cv      = 0.0
    dq_hs_cv      = 0.0
    dn_hs_cv      = 0.0
    dn_cl_sc      = 0.0
    dn_ci_mul     = 0.0
    dq_ci_mul     = 0.0
    dn_ci_me      = 0.0
    dq_ci_me      = 0.0
    dn_hs_me      = 0.0
    dq_hs_me      = 0.0
    dn_hg_me      = 0.0
    dq_hg_me      = 0.0
    dn_ci_ev      = 0.0
    dq_ci_ev      = 0.0
    dn_hs_ev      = 0.0
    dq_hs_ev      = 0.0
    dn_hg_ev      = 0.0
    dq_hg_ev      = 0.0
    dn_ci_eme_ic  = 0.0
    dq_ci_eme_ic  = 0.0
    dn_ci_eme_ri  = 0.0
    dq_ci_eme_ri  = 0.0
    dn_hs_eme_sc  = 0.0
    dq_hs_eme_sc  = 0.0
    dn_hs_eme_rs  = 0.0
    dq_hs_eme_rs  = 0.0
    dn_hg_eme_gc  = 0.0
    dq_hg_eme_gc  = 0.0
    dn_hg_eme_gr  = 0.0
    dq_hg_eme_gr  = 0.0
    dn_cl_se      = 0.0
    dq_cl_se      = 0.0
    dn_ci_se      = 0.0
    dq_ci_se      = 0.0
    dn_hr_se      = 0.0
    dq_hr_se      = 0.0
    dn_hs_se      = 0.0
    dq_hs_se      = 0.0
    dn_hg_se      = 0.0
    dq_hg_se      = 0.0
    precep_hr     = 0.0
    precep_ci     = 0.0
    precep_hs     = 0.0
    precep_hg     = 0.0
    dq_cl_sa      = 0.0
    dn_cl_sa      = 0.0
    ret_cc        = 0.0

    precep_l      = 0.0
    precep_i      = 0.0

    ! reset flags
    flag_warn_update = .false.
    flag_warn_test   = .false.
    flag_warn_size   = .true.
    flag_nan         = .false.
    if (l_sb_dbg) then
      flag_do_dbg_up  = .true.
    endif

    delt = rdt/ (4. - dble(rk3step))

    if (timee .eq. 0. .and. rk3step .eq. 1 .and. myid .eq. 0) then
      write(*,*) 'l_lognormal',l_lognormal
      write(*,*) 'rhof(1)', rhof(1),' rhof(10)', rhof(10)
      write(*,*) 'l_mur_cst',l_mur_cst,' mur_cst',mur_cst
      write(*,*) 'nuc = param'
    endif

  !*********************************************************************
  ! remove negative values for hydrometereors and clouds
  !*********************************************************************
    if (l_rain) then
       if ( sum(q_hr, q_hr<0.) > 0.000001*sum(q_hr) ) then
         write(*,*)'amount of neg. q_hr thrown away is too high  ',timee,' sec'
       end if
       if ( sum(n_hr, n_hr<0.) > 0.000001*sum(n_hr) ) then
         write(*,*)'amount of neg. n_hr thrown away is too high  ',timee,' sec'
       end if
       if ((sum(q_hs,q_hs<0.)+sum(q_hg,q_hg<0.))>0.000001*(sum(q_hs)+sum(q_hg))) then
         write(*,*)'amount of neg. q_hs, q_hg thrown away is too high  ',timee,' sec'
       end if
       if ((sum(n_hs,n_hs<0.)+sum(n_hg,n_hg<0.))>0.000001*(sum(n_hs)+sum(n_hg))) then
          write(*,*)'amount of neg. n_hs, n_hg and n_hg thrown away is too high ',timee,' sec'
       end if
       if ((sum(q_cl,q_cl<0.)+sum(q_ci,q_ci<0.))>0.000001*(sum(q_cl)+sum(q_ci))) then
          write(*,*)'amount of neg. q_cl and q_ci thrown away is too high ',timee,' sec'
       end if
       if ((sum(n_cl,n_cl<0.)+sum(n_ci,n_ci<0.))>0.000001*(sum(n_cl)+sum(n_ci))) then
          write(*,*)'amount of neg. n_cl and n_ci thrown away is too high ',timee,' sec'
       end if
       if ((sum(n_cc,n_cc<0.))>0.000001*(sum(n_cc))) then
          write(*,*)'amount of neg. n_cc thrown away is too high ',timee,' sec'
       end if

       ! replacing negative values
       do k=1,k1
       do j=2,j1
       do i=2,i1
          if (n_hr(i,j,k) < 0.)  then
            n_hr(i,j,k) = 0.
          endif
          if (q_hr(i,j,k) < 0.)  then
            q_hr(i,j,k) = 0.
          endif
          !#sb3 for other species
          if (n_hs(i,j,k) < 0.)  then
            n_hs(i,j,k) = 0.
          endif
          if (n_hg(i,j,k) < 0.)  then
            n_hg(i,j,k) = 0.
          endif
          if (q_hs(i,j,k) < 0.)  then
            q_hs(i,j,k) = 0.
          endif
          if (q_hg(i,j,k) < 0.)  then
            q_hg(i,j,k) = 0.
          endif
          if (n_cl(i,j,k) < 0.)  then
            n_cl(i,j,k) = 0.
          endif
          if (n_ci(i,j,k) < 0.)  then
            n_ci(i,j,k) = 0.
          endif
          if (q_cl(i,j,k) < 0.)  then
            q_cl(i,j,k) = 0.
          endif
          if (q_ci(i,j,k) < 0.)  then
            q_ci(i,j,k) = 0.
          endif
          if (n_cc(i,j,k) < 0.)  then
            n_cc(i,j,k) = 0.
          endif
       enddo
       enddo
       enddo

    end if   ! l_rain

    !testing for a noticable amount of graupel and snow
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! rain :
      if ((q_hr(i,j,k) .gt. q_hr_min).and.(n_hr(i,j,k).gt.0.0))  then
         q_hr_mask (i,j,k) = .true.
      else
         q_hr_mask (i,j,k) = .false.
      endif
      ! snow :
      if ((q_hs(i,j,k) .gt. q_hs_min).and.(n_hs(i,j,k).gt.0.0))  then
         q_hs_mask (i,j,k) = .true.
      else
         q_hs_mask (i,j,k) = .false.
      endif
      ! graupel :
      if ((q_hg(i,j,k) .gt. q_hg_min).and.(n_hg(i,j,k).gt.0.0))  then
         q_hg_mask (i,j,k) = .true.
      else
         q_hg_mask (i,j,k) = .false.
      endif
    enddo
    enddo
    enddo


  !******************************************************
  ! calculate qltot and initialize cloud droplet number Nc
  !******************************************************
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qltot (i,j,k) = q_cl(i,j,k) + q_hr (i,j,k)
      ! instead of: ql0  (i,j,k) + q_hr (i,j,k)
      ! i.e - cloud water + rain areas

      ! calculating ice water content already existing
      qicew = q_ci(i,j,k) ! +q_hs(i,j,k)+q_hg(i,j,k)

      ! liquid clouds - if there is enough liquid in clouds
      qcliq =  q_cl(i,j,k) ! ql0(i,j,k)-q_ci(i,j,k)

      ! clouds- original definition
      if (ql0(i,j,k) .ge. qcmin)  then ! if (ql0(i,j,k) >qcmin)  then
         qcmask (i,j,k) = .true.
      else
         qcmask (i,j,k) = .false.
      end if

      ! liquid clouds
      if ((q_cl(i,j,k).ge. qcliqmin).and.(n_cl(i,j,k).gt.0.0))  then
         q_cl_mask (i,j,k) = .true.
      else
         q_cl_mask (i,j,k) = .false.
      end if

      ! ice clouds
      if ((qicew.ge. qicemin).and.(n_ci(i,j,k).gt.0.0))  then
         q_ci_mask (i,j,k) = .true.
      else
         q_ci_mask (i,j,k) = .false.
      endif
    enddo
    enddo
    enddo

    if(l_sb_dbg_extra) then
      write(6,*) 'total count of q_cl_mask ',count(q_cl_mask(2:i1,2:j1,1:k1))
      write(6,*) 'total count of q_ci_mask ',count(q_ci_mask(2:i1,2:j1,1:k1))
      write(6,*) 'total count of q_hr_mask ',count(q_hr_mask(2:i1,2:j1,1:k1))
      write(6,*) 'total count of q_hs_mask ',count(q_hs_mask(2:i1,2:j1,1:k1))
      write(6,*) 'total count of q_hg_mask ',count(q_hg_mask(2:i1,2:j1,1:k1))

      if (any((qt0(2:i1,2:j1,1:k1)-qvsl(2:i1,2:j1,1:k1)-q_cl(2:i1,2:j1,1:k1)-q_ci(2:i1,2:j1,1:k1)).lt. 0.)) then
        write(6,*) 'undersaturated clouds in gridpoints ', &
          count((qt0(2:i1,2:j1,1:k1)-qvsl(2:i1,2:j1,1:k1)-q_cl(2:i1,2:j1,1:k1)-q_ci(2:i1,2:j1,1:k1)).lt. 0.)
      endif

      do k=1,k1
      do j=2,j1
      do i=2,i1
        nrtest=n_cl(i,j,k) ! considering n_clp only
        ! in clouds only
        if (nrtest.le.0.0) then
          xtest(i,j,k) =0.0
        else
          !  calculating new droplet size
          xtest(i,j,k) =  q_cl(i,j,k) / ( n_cl(i,j,k)+eps0)
                  ! o: rhof(k) * (svm(i,j,k,iqsv) +  delt*q_clp(i,j,k)) / (svm(i,j,k,insv)+delt*n_clp(i,j,k))
        endif
      enddo
      enddo
      enddo

      ! warning if droplets too large
      if(any( xtest(2:i1,2:j1,1:k1) .gt. xc_bmax)) then
        write(6,*) 'warning: cloud droplets at the start of the step too large'
        write(6,*) '   in gridpoints ', count(xtest(2:i1,2:j1,1:k1).gt.xc_bmax), ' larger than ',xc_bmax
        write(6,*) '   max value ', maxval(xtest), ' in ', maxloc(xtest)
      endif
    endif ! l_sb_dbg_extra

  !*********************************************************************
  ! calculate Rain DSD integral properties & parameters xr, Dvr, lbdr, mur
  !*********************************************************************
    call integrals_bulk3
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'before.')! #d
    call nucleation3      ! cloud nucleation
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'nuclea.')! #d
    call icenucle3        ! ice nucleation  ! #Bb ! #iceout
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'ice nuc')! #d

    ! --------------------------------------------------------------
    ! freezing of water droplets
    ! --------------------------------------------------------------
    call hethomfreez3

    ! --------------------------------------------------------------
    ! deposition processes
    ! --------------------------------------------------------------
    call deposit3

    ! --------------------------------------------------------------
    ! snow aggregation and self-collection
    ! collision processes for snow
    ! --------------------------------------------------------------
    call ice_aggr_snow_self

    ! --------------------------------------------------------------
    ! - rimings
    ! --------------------------------------------------------------
    call coll3
    call coll_grg3 ! riming g+r -> g  !#t2 !#t4
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'g+r->g ')! #d

    ! --------------------------------------------------------------
    ! - raindrop freezing
    ! --------------------------------------------------------------
    call rainhetfreez3        ! heterogeneous freezing! #iceout
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'fre.rai')! #d

    ! --------------------------------------------------------------
    ! - collision with conversion
    ! --------------------------------------------------------------
    call coll_rig3 ! riming r+i -> g  !#t2 !#t4
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'r+i->g ')! #d
    call coll_rsg3 ! riming r+i -> g  !#t7
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'r+s->g ')! #d

    !--------------------------------------------------------------
    ! conversions
    !--------------------------------------------------------------
    call ice_multi3 ! ice multiplication of Hallet and Mossop
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'ice.mul')! #d

    if (l_sb_conv_par) then
      call conv_partial3 ! partial conversion
      if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'par.con')! #d
    endif

    ! --------------------------------------------------------------
    ! - melting and evaporation of ice particles
    ! --------------------------------------------------------------
    call evapmelting3 ! melting of ice particles
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'melting')! #d

    ! --------------------------------------------------------------
    ! - basic warm processes
    !  --------------------------------------------------------------
    call autoconversion3
    if(l_sb_dbg_extra)   call check_nan(flag_do_dbg,flag_nan,'autocon')! #d
    call cloud_self3
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'cl self')! #d
    call accretion3
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'accret.')! #d

    call evap_rain3 ! rain evaporation
    if(l_sb_dbg_extra) call check_nan(flag_do_dbg,flag_nan,'evap_r.')! #d

  ! =============================
  ! saturation adjustment
  ! =============================
    call satadj3

  !==================================
  ! sedimentation part
  !==================================
    call sedim_rain3
    call sedim_cl3
    call sedim_ice3
    call sedim_snow3
    call sedim_graupel3 !#b2t3

    ! recovery of ccn aerosols
    call recover_cc

  !*********************************************************************
  ! remove negative values and non physical low values
  !*********************************************************************
  if (l_corr_neg_qt) then
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iq_hr)+(svp(i,j,k,iq_hr)+q_hrp(i,j,k))*delt
      nrtest=svm(i,j,k,in_hr)+(svp(i,j,k,in_hr)+n_hrp(i,j,k))*delt
      if ((qrtest < q_hr_min) .or. (nrtest < 0.0) ) then
        ! correction, after Jerome's implementation in Gales
        qtp(i,j,k) = qtp(i,j,k) + svm(i,j,k,iq_hr)/delt + svp(i,j,k,iq_hr) + q_hrp(i,j,k)
        thlp(i,j,k) = thlp(i,j,k)  -  &
               (rlv/(cp*exnf(k)))*(svm(i,j,k,iq_hr)/delt + svp(i,j,k,iq_hr) + q_hrp(i,j,k))

        svp(i,j,k,iq_hr) = - svm(i,j,k,iq_hr)/delt
        svp(i,j,k,in_hr) = - svm(i,j,k,in_hr)/delt
      else
        svp(i,j,k,iq_hr)=svp(i,j,k,iq_hr)+q_hrp(i,j,k)
        svp(i,j,k,in_hr)=svp(i,j,k,in_hr)+n_hrp(i,j,k)
        ! adjust negative qr tendencies at the end of the time-step
      end if
    enddo
    enddo
    enddo

    ! == snow ==
    iqsv = iq_hs
    insv = in_hs
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iqsv)+(svp(i,j,k,iqsv)+q_hsp(i,j,k))*delt
      nrtest=svm(i,j,k,insv)+(svp(i,j,k,insv)+n_hsp(i,j,k))*delt
      if ((qrtest .lt. qsnowmin) .or. (nrtest .le. 0.0) ) then
        ! correction, after Jerome's implementation in Gales
        qtp(i,j,k) = qtp(i,j,k) + svm(i,j,k,iqsv)/delt + svp(i,j,k,iqsv) + q_hsp(i,j,k)
        thlp(i,j,k) = thlp(i,j,k)  -  &
               (rlvi/(cp*exnf(k)))*(svm(i,j,k,iqsv)/delt + svp(i,j,k,iqsv) + q_hsp(i,j,k))
        svp(i,j,k,iqsv) = - svm(i,j,k,iqsv)/delt
        svp(i,j,k,insv) = - svm(i,j,k,insv)/delt
      else
        svp(i,j,k,iqsv)=svp(i,j,k,iqsv)+q_hsp(i,j,k)
        svp(i,j,k,insv)=svp(i,j,k,insv)+n_hsp(i,j,k)
      end if
    enddo
    enddo
    enddo

    ! == graupel ==
    iqsv = iq_hg
    insv = in_hg
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iqsv)+(svp(i,j,k,iqsv)+q_hgp(i,j,k))*delt
      nrtest=svm(i,j,k,insv)+(svp(i,j,k,insv)+n_hgp(i,j,k))*delt
      if ((qrtest .lt. qgrmin) .or. (nrtest .le. 0.0) ) then
        ! correction, after Jerome's implementation in Gales
        qtp(i,j,k) = qtp(i,j,k)+svm(i,j,k,iqsv)/delt+svp(i,j,k,iqsv)+q_hgp(i,j,k)
        thlp(i,j,k) = thlp(i,j,k) -  &
               (rlvi/(cp*exnf(k)))*(svm(i,j,k,iqsv)/delt+svp(i,j,k,iqsv)+q_hgp(i,j,k))
        svp(i,j,k,iqsv) = - svm(i,j,k,iqsv)/delt
        svp(i,j,k,insv) = - svm(i,j,k,insv)/delt
      else
        svp(i,j,k,iqsv)=svp(i,j,k,iqsv)+q_hgp(i,j,k)
        svp(i,j,k,insv)=svp(i,j,k,insv)+n_hgp(i,j,k)
      endif
    enddo
    enddo
    enddo

    ! == cloud ice ==
    iqsv = iq_ci
    insv = in_ci
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iqsv)+(svp(i,j,k,iqsv)+q_cip(i,j,k))*delt
      nrtest=svm(i,j,k,insv)+(svp(i,j,k,insv)+n_cip(i,j,k))*delt
      if ((qrtest .lt. qicemin) .or. (nrtest .le. 0.0) ) then
        ! correction, after Jerome's implementation in Gales
        qtp(i,j,k) = qtp(i,j,k)+svm(i,j,k,iqsv)/delt+svp(i,j,k,iqsv)+q_cip(i,j,k)
        thlp(i,j,k) = thlp(i,j,k) -  &
               (rlvi/(cp*exnf(k)))*(svm(i,j,k,iqsv)/delt+svp(i,j,k,iqsv)+q_cip(i,j,k))

        svp(i,j,k,iqsv) = - svm(i,j,k,iqsv)/delt
        svp(i,j,k,insv) = - svm(i,j,k,insv)/delt
        ! recovery of aerosols
        ! ret_cc(i,j,k)   = ret_cc(i,j,k)+max(0.0,nrtest/delt)
      else
        svp(i,j,k,iqsv)=svp(i,j,k,iqsv)+q_cip(i,j,k)
        svp(i,j,k,insv)=svp(i,j,k,insv)+n_cip(i,j,k)
      endif
    enddo
    enddo
    enddo

    ! == cloud liquid water ==
    iqsv = iq_cl
    insv = in_cl
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iqsv)+(svp(i,j,k,iqsv)+q_clp(i,j,k))*delt
      nrtest=svm(i,j,k,insv)+(svp(i,j,k,insv)+n_clp(i,j,k))*delt
      if ((qrtest .lt. qcliqmin) .or. (nrtest .le. 0.0) ) then
        ! correction, after Jerome's implementation in Gales
        svp(i,j,k,iqsv) = - svm(i,j,k,iqsv)/delt
        svp(i,j,k,insv) = - svm(i,j,k,insv)/delt
      else
        svp(i,j,k,iqsv)=svp(i,j,k,iqsv)+q_clp(i,j,k)
        svp(i,j,k,insv)=svp(i,j,k,insv)+n_clp(i,j,k)
      endif
    enddo
    enddo
    enddo
  else  ! l_corr_neg_qt
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iq_hr)+(svp(i,j,k,iq_hr)+q_hrp(i,j,k))*delt
      nrtest=svm(i,j,k,in_hr)+(svp(i,j,k,in_hr)+n_hrp(i,j,k))*delt
      if ((qrtest < q_hr_min) .or. (nrtest < 0.0) ) then
        ! correction, after Jerome's implementation in Gales
        svp(i,j,k,iq_hr) = - svm(i,j,k,iq_hr)/delt
        svp(i,j,k,in_hr) = - svm(i,j,k,in_hr)/delt
      else
        svp(i,j,k,iq_hr)=svp(i,j,k,iq_hr)+q_hrp(i,j,k)
        svp(i,j,k,in_hr)=svp(i,j,k,in_hr)+n_hrp(i,j,k)
      end if
    enddo
    enddo
    enddo

    ! == snow ==
    iqsv = iq_hs
    insv = in_hs
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iqsv)+(svp(i,j,k,iqsv)+q_hsp(i,j,k))*delt
      nrtest=svm(i,j,k,insv)+(svp(i,j,k,insv)+n_hsp(i,j,k))*delt
      if ((qrtest .lt. qsnowmin) .or. (nrtest .le. 0.0) ) then
        ! correction, after Jerome's implementation in Gales
        svp(i,j,k,iqsv) = - svm(i,j,k,iqsv)/delt
        svp(i,j,k,insv) = - svm(i,j,k,insv)/delt
      else
       svp(i,j,k,iqsv)=svp(i,j,k,iqsv)+q_hsp(i,j,k)
       svp(i,j,k,insv)=svp(i,j,k,insv)+n_hsp(i,j,k)
      end if
    enddo
    enddo
    enddo

    ! == graupel ==
    iqsv = iq_hg
    insv = in_hg
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iqsv)+(svp(i,j,k,iqsv)+q_hgp(i,j,k))*delt
      nrtest=svm(i,j,k,insv)+(svp(i,j,k,insv)+n_hgp(i,j,k))*delt
      if ((qrtest .lt. qgrmin) .or. (nrtest .le. 0.0) ) then
        ! correction, after Jerome's implementation in Gales
        svp(i,j,k,iqsv) = - svm(i,j,k,iqsv)/delt
        svp(i,j,k,insv) = - svm(i,j,k,insv)/delt
      else
        svp(i,j,k,iqsv)=svp(i,j,k,iqsv)+q_hgp(i,j,k)
        svp(i,j,k,insv)=svp(i,j,k,insv)+n_hgp(i,j,k)
      end if
    enddo
    enddo
    enddo

    ! == cloud ice ==
    iqsv = iq_ci
    insv = in_ci
    do k=1,k1
    do j=2,j1
    do i=2,i1
       qrtest=svm(i,j,k,iqsv)+(svp(i,j,k,iqsv)+q_cip(i,j,k))*delt
       nrtest=svm(i,j,k,insv)+(svp(i,j,k,insv)+n_cip(i,j,k))*delt
      if ((qrtest .lt. qicemin) .or. (nrtest .le. 0.0) ) then
        ! (nrtest < 0.0)  correction, after Jerome's implementation in Gales
        svp(i,j,k,iqsv) = - svm(i,j,k,iqsv)/delt
        svp(i,j,k,insv) = - svm(i,j,k,insv)/delt
      else
        svp(i,j,k,iqsv)=svp(i,j,k,iqsv)+q_cip(i,j,k)
        svp(i,j,k,insv)=svp(i,j,k,insv)+n_cip(i,j,k)
      endif
    enddo
    enddo
    enddo

    ! == cloud liquid water ==
    iqsv = iq_cl
    insv = in_cl
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iqsv)+(svp(i,j,k,iqsv)+q_clp(i,j,k))*delt
      nrtest=svm(i,j,k,insv)+(svp(i,j,k,insv)+n_clp(i,j,k))*delt
      if ((qrtest .lt. qcliqmin) .or. (nrtest .le. 0.0) ) then
        ! (nrtest < n_c_min) ) then ! (nrtest < 0.0) ) correction, after Jerome's implementation in Gales
        svp(i,j,k,iqsv) = - svm(i,j,k,iqsv)/delt
        svp(i,j,k,insv) = - svm(i,j,k,insv)/delt
      else
        svp(i,j,k,iqsv)=svp(i,j,k,iqsv)+q_clp(i,j,k)
        svp(i,j,k,insv)=svp(i,j,k,insv)+n_clp(i,j,k)
      endif
    enddo
    enddo
    enddo

  endif ! l_corr__neg_qt

  ! == CCN checking ==
  if(.not.l_c_ccn) then
    insv = in_cc
    do k=1,k1
    do j=2,j1
    do i=2,i1
      nrtest=svm(i,j,k,insv)+(svp(i,j,k,insv)+n_ccp(i,j,k))*delt
      if (nrtest .le.0.0 ) then ! correction
        svp(i,j,k,insv) = - svm(i,j,k,insv)/delt
      else
        svp(i,j,k,insv)=svp(i,j,k,insv)+n_ccp(i,j,k)
        ! no change for qt and th
      end if
    enddo
    enddo
    enddo
  endif ! not l_c_ccn

  ! == update =====
  ! and updating main prognostic variables by contribution from mphys processes
  do k=1,k1
  do j=2,j1
  do i=2,i1
    thlp(i,j,k) = thlp(i,j,k)+thlpmcr(i,j,k)
    qtp(i,j,k) = qtp(i,j,k)+qtpmcr(i,j,k)
  enddo
  enddo
  enddo

  if(l_sb_dbg_extra) then
    ! estimating new xc
    insv = in_cl
    iqsv = iq_cl
    do k=1,k1
    do j=2,j1
    do i=2,i1
      nrtest=svm(i,j,k,insv)+delt*svp(i,j,k,insv) ! since n_clp already included
      ! in clouds only
      if (nrtest<0.0 ) then
        xtest(i,j,k) =0.0
      else
        !  calculating new droplet size
        xtest(i,j,k) = (svm(i,j,k,iqsv) + delt*svp(i,j,k,iqsv)) / (svm(i,j,k,insv)+delt*svp(i,j,k,insv)+eps0)
      endif
    enddo
    enddo
    enddo

    ! warning if droplets too large
    if(any(xtest(2:i1,2:j1,1:k1) .gt. xc_bmax)) then
      write(6,*) 'warning: cloud droplets at the end of the step too large'
      write(6,*) '   in ', count(xtest(2:i1,2:j1,1:k1).gt.xc_bmax), ' gridpoints larger than ',xc_bmax
      write(6,*) '   max value ', maxval(xtest), ' in ', maxloc(xtest)
    endif
  endif ! l_sb_dbg


  ! # statistics
  ! call microphysics statistics - just once per step
  call bulkmicrotend3 ! #t5
  call bulkmicrostat3 ! #t5

  if(l_sb_dbg_extra) call check_sizes(flag_do_dbg,flag_warn_size,'all_ups')
  call check_allupdates(flag_do_dbg_up)

  ! ending line
  if(l_sb_dbg_extra) then
      write(6,*) ' =========================='
  endif

  ! deallocating process variables
  deallocate( xtest)  ! #d
  !
  deallocate( dn_cl_nu,dn_ci_inu                                      &
             ,dq_hr_au,dn_cl_au,dn_hr_au,dq_hr_ac,dn_cl_ac            &
             ,dn_hr_br,dn_hr_sc                                       &
             ,dq_hr_ev, dn_hr_ev                                      &
             ,dq_ci_rime,dq_hs_rime,dq_hg_rime,dq_hghr_rime           &
             ,dq_hshr_rime ,dn_hr_rime_hs                             &
             ,dn_cl_rime_ci,dn_cl_rime_hs,dn_cl_rime_hg,dn_hr_rime_hg &
             ,dq_hr_rime_ri,dq_ci_rime_ri                             &
             ,dn_ci_rime_ri ,dn_hr_rime_ri                            &
             ,dq_hr_col_rs,dq_hs_col_rs                               &
             ,dn_hr_col_rs,dn_hs_col_rs                               &
             ,dq_hr_col_ri,dq_ci_col_ri                               &
             ,dn_ci_col_ri,dn_hr_col_ri                               &
             ,dq_cl_het,dq_hr_het,dn_hr_het                           &
             ,dq_cl_hom,                                              &
             ,dq_ci_col_iis,dn_ci_col_iis                             &
             ,dn_hs_col_sss,dq_hsci_col,dn_ci_col_hs                  &
             ,dq_hghs_col, dn_hs_col_hg                               &
             ,dq_ci_cv,dq_hs_cv, dn_ci_cv,dn_hs_cv                    &
             ,dn_cl_sc                                                &
             ,dn_ci_mul, dq_ci_mul                                    &
             ,dq_ci_me,dq_hs_me, dq_hg_me                             &
             ,dn_ci_me,dn_hs_me, dn_hg_me                             &
             ,dq_ci_ev,dq_hs_ev, dq_hg_ev                             &
             ,dn_ci_ev,dn_hs_ev, dn_hg_ev                             &
             ,dn_ci_eme_ic,dq_ci_eme_ic,dn_ci_eme_ri,dq_ci_eme_ri     &
             ,dn_hs_eme_sc,dq_hs_eme_sc,dn_hs_eme_rs,dq_hs_eme_rs     &
             ,dn_hg_eme_gc,dq_hg_eme_gc,dn_hg_eme_gr,dq_hg_eme_gr     &
             ,dq_cl_se,dq_ci_se,dq_hr_se,dq_hs_se,dq_hg_se            &
             ,dn_cl_se,dn_ci_se,dn_hr_se,dn_hs_se,dn_hg_se            &
             ,precep_hr, precep_ci, precep_hs,precep_hg               &
             ,dq_cl_sa ,dn_cl_sa                                      &
             ,ret_cc                                                  &
            )
   ! deallocate sizes
   deallocate( D_cl,D_ci,D_hr,D_hs,D_hg                               &
              ,v_cl,v_ci,v_hr,v_hs,v_hg                               &
             )

 end subroutine bulkmicro3

! ===============================
!  cloud initialisation
! ===============================
! initialises cloud water and cloud droplet number when called
!  called from : bulkmicro3 if flag l_cloud_init==.false.
! steps:
!  - turns all water above saturation to clouds
!  - sets reasonable cloud droplet number
!     - based on proposed values Nc0 [ m^-3]
!     - limit if below n_ccn
!     - limit if mean x_cl smaller than xcmin
!
subroutine initclouds3
  use modglobal, only : ih,i1,jh,j1,k1,kmax,zf
  use modfields, only : rhof, qt0, svm, sv0, svp, qvsl, qsat
  implicit none

  integer :: i,j,k
  real :: x_min, Nc_set, q_tocl,n_prop  ! availabel water and proposed size

  x_min = xc0_min  ! minimal size of droplets
  Nc_set =  Nc0    ! prescribed number of droplets

  ! initialisation itself
  if (l_c_ccn) then
    do k=1,k1
    do j=2,j1
    do i=2,i1
      q_tocl = qt0(i,j,k)-qvsl(i,j,k)           ! get amount of available water for liquid clouds
      if (q_tocl>0.0) then
        ! prepare number of droplet
        n_prop = Nc_set/rhof(k)                 ! to get number of droplets in kg^{-1}
        ! n_prop = min(n_prop,sv0(i,j,k,in_cc)) ! number has to be lower than number of available CCN
        n_prop = min(n_prop,q_tocl/x_min)       ! droplets smaller than minimal size
        ! save values to arrays
        sv0(i,j,k,in_cl) = n_prop
        svm(i,j,k,in_cl) = n_prop
        sv0(i,j,k,iq_cl) = q_tocl
        svm(i,j,k,iq_cl) = q_tocl
      endif
    enddo
    enddo
    enddo
  else
    do k=1,k1
    do j=2,j1
    do i=2,i1
      q_tocl = qt0(i,j,k)-qvsl(i,j,k)           ! get amount of available water for liquid clouds
      if (q_tocl>0.0) then
        ! prepare number of droplet
        n_prop = Nc_set/rhof(k)                 ! to get number of droplets in kg^{-1}
        n_prop = min(n_prop,sv0(i,j,k,in_cc))   ! number has to be lower than number of available CCN
        n_prop = min(n_prop,q_tocl/x_min)       ! droplets smaller than minimal size

        ! save values to arrays
        sv0(i,j,k,in_cl) = n_prop
        svm(i,j,k,in_cl) = n_prop
        sv0(i,j,k,iq_cl) = q_tocl
        svm(i,j,k,iq_cl) = q_tocl
      endif
    enddo
    enddo
    enddo
  endif

end subroutine initclouds3

! ===============================
!  ccn initialisation
! ===============================
! initialises cloud water and cloud droplet number when called
!  called from : bulkmicro3 if flag l_ccn_init==.false.
! steps:
!     - based on proposed values Nc0 [ m^-3]
!     - set ccn fields
!
subroutine initccn3
  use modglobal, only : ih,i1,jh,j1,k1,kmax
  use modfields, only : rhof, svm, sv0, svp
  implicit none

  integer :: i,j,k
  real ::   Nccn_set, n_prop            ! available water and proposed size

  Nccn_set = Nccn0                      ! prescribed number of ccn

  do k=1,k1
  do j=2,j1
  do i=2,i1
    ! prepare number of CCNs
    n_prop = Nccn_set/rhof(k)           ! to get number of droplets in kg^{-1}
    ! save values to arrays
    sv0(i,j,k,in_cc) = n_prop
    svm(i,j,k,in_cc) = n_prop
  enddo
  enddo
  enddo

end subroutine initccn3


! =============================
! saturation adjustment
! ============================
! Performs the saturation adjustment to cloud water
!
! Uses split operator:
!   - calculates how much water is nucleated and consumed
!   - calculates remaining amount of water available for condensation
!   - dumps remaining water into cloud
!
!   In addition, check the average size of cloud droplets
!   and evaporates droplets that proportionally
!
subroutine satadj3
  use modglobal, only : ih,i1,jh,j1,k1,kmax
  use modfields, only : rhof, qt0, svm, svp, qvsl
  implicit none

  integer :: i,j,k
  real :: ntest, n_bmax, cogr_max, ql_res

  ! initialisation of process updates
  dq_cl_sa = 0.0
  dn_cl_sa = 0.0

  ! calculation
  if (l_sb_all_or) then
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! note: threshold might be lower then threshold for cloud computations, obviously
      ntest=svm(i,j,k,in_cl)+n_clp(i,j,k)*delt
      if (ntest.gt.0.0) then
        !
        ! remaining water =
        !    + condesable water available
        !    - already condensed
        !    - newly condendsed
        !    - removed by mphys processed
        !
        !  calculating amount of available water
        ql_res=(qt0(i,j,k)-qvsl(i,j,k)-svm(i,j,k,iq_cl))+delt*(qtpmcr(i,j,k)-q_clp(i,j,k))
        dq_cl_sa(i,j,k) = ql_res/delt
        if ((svm(i,j,k,iq_cl)+delt*(q_clp(i,j,k)+dq_cl_sa(i,j,k))).lt.0.0) then
          dn_cl_sa(i,j,k) = -svm(i,j,k,in_cl)/delt-n_clp(i,j,k)
          dq_cl_sa(i,j,k) = -svm(i,j,k,iq_cl)/delt-q_clp(i,j,k)
        endif
      endif
    enddo
    enddo
    enddo
  else ! l_sb_all_or
    if (l_sb_dumpall) then
      ! dump all water
      do k=1,k1
      do j=2,j1
      do i=2,i1
        ! note: threshold might be lower then threshold for cloud computations, obviously
        if (q_cl_mask(i,j,k)) then
          ntest= svm(i,j,k,in_cl)+n_clp(i,j,k)*delt ! #t1
          if (ntest.gt.0.0) then
            !
            ! remaining water =
            !    + condesable water available
            !    - already condensed
            !    - newly condendsed
            !    - removed by mphys processed
            !

            !  calculating amount of available water
            ql_res=(qt0(i,j,k)-qvsl(i,j,k)-svm(i,j,k,iq_cl))+delt*(qtpmcr(i,j,k)-q_clp(i,j,k))

            ! limiting so it does not remove more water than in clouds
            dq_cl_sa(i,j,k) = max((-svm(i,j,k,iq_cl)/delt)-q_clp(i,j,k),ql_res/delt)

            ! adjusting number of cloud droplets
            ! - calculate min size with this amount of water
            n_bmax = (svm(i,j,k,iq_cl)+delt*(q_clp(i,j,k)+dq_cl_sa(i,j,k)))/(0.1*x_cl_bmin)
            n_bmax = max(n_bmax, 0.0)

            ! of course we do not want negative values - but that is alread sorted above
            ! - remove droplets so that mean size in not less than
            dn_cl_sa(i,j,k) = min(0.0, n_bmax-ntest)/delt

            ! limit change so not in negative numbers
            dn_cl_sa(i,j,k)=max((-svm(i,j,k,in_cl)/delt)-svp(i,j,k,in_cl)-n_clp(i,j,k),dn_cl_sa(i,j,k))
          endif
        endif
      enddo
      enddo
      enddo
    else ! l_sb_dumpall
      do k=1,k1
      do j=2,j1
      do i=2,i1
        ! note: threshold might be lower then threshold for cloud computations, obviously
        if (q_cl_mask(i,j,k)) then
          ntest=svm(i,j,k,in_cl)+n_clp(i,j,k)*delt !#t1
          if (ntest.gt.0.0) then
            !
            ! and now if we want to enforce limit on cloud droplet size
            !
            ql_res=(qt0(i,j,k)-qvsl(i,j,k)-svm(i,j,k,iq_cl))+delt*(qtpmcr(i,j,k)-q_clp(i,j,k))

            ! calculate maximal available update
            cogr_max =(1.0/delt)*(x_cogr_max*(svm(i,j,k,in_cl)+delt*n_clp(i,j,k))-svm(i,j,k,iq_cl))

            ! dump just what is below max size by condensation growth
            dq_cl_sa(i,j,k) = min(cogr_max-q_clp(i,j,k),ql_res/delt)

            ! ie. either whole amount, or only that much that droplet size will be: xc_cogr_max
            ! other possibility: require it to be larger than 0
            ! and prevent negative values of svm + delt *svp
            dq_cl_sa(i,j,k) = max((-svm(i,j,k,iq_cl)/delt)-q_clp(i,j,k),dq_cl_sa(i,j,k))

            ! adjusting number of cloud droplets
            ! - calculate min size with this amount of water
            n_bmax = (svm(i,j,k,iq_cl)+delt*(q_clp(i,j,k)+dq_cl_sa(i,j,k)))/(0.5*x_cl_bmin)
            n_bmax = max(n_bmax, 0.0)

            ! - remove droplets so that mean size in not by order of magnitude less than x_cl_bmin
            dn_cl_sa(i,j,k) = min(0.0, n_bmax-ntest)/delt

            ! limit change so not in negative numbers
            dn_cl_sa(i,j,k)=max((-svm(i,j,k,in_cl)/delt)-svp(i,j,k,in_cl)-n_clp(i,j,k),dn_cl_sa(i,j,k))

            !-------------------------------------
            !! cloud water mixing ratio
            ! + mphys changes in cloud water
            ! + remaining water:
            !    + condesable water available
            !    - already condensed
            !    - newly condendsed
            !      ( {change in cloud water} = {newly condendsed or deposited} + {removed by cloud processes}  )
            !      ( - {newly condendsed or deposited} =
            !          - ({change in liquid cloud water}+{change in ice cloud water})
            !          + {removed by cloud processes}
            !      )
            ! -----------------------------
          endif
        endif
      enddo
      enddo
      enddo
    endif
  endif ! l_sb_all_or

  ! and update
  do k=1,k1
  do j=2,j1
  do i=2,i1
    q_clp(i,j,k)  = q_clp(i,j,k) + dq_cl_sa(i,j,k)
    n_clp(i,j,k)  = n_clp(i,j,k) + dn_cl_sa(i,j,k)
  enddo
  enddo
  enddo
end subroutine satadj3


! calculating rain and other integrals
subroutine integrals_bulk3
  use modglobal, only : ih,i1,jh,j1,k1,kmax
  use modfields, only : rhof
  implicit none

  real , allocatable :: N_r0(:,:,:), lbdr_try(:,:,:)  ! #sb3
  integer :: i,j,k
  real :: xr, xr_try, N_r0_try

  allocate ( N_r0     (2-ih:i1+ih,2-jh:j1+jh,k1)      &
            ,lbdr_try (2-ih:i1+ih,2-jh:j1+jh,k1)      &  ! #sb3
           )

  N_r0     = 0.0
  lbdr_try = 0.0

  ! calculations
  if (l_rain) then
    Dvr  (2:i1,2:j1,1:k1) = 0.
    mur  (2:i1,2:j1,1:k1) = 30.
    lbdr (2:i1,2:j1,1:k1) = 0.

    if (l_sb) then
      if(l_sb_classic) then
        do k=1,k1
        do j=2,j1
        do i=2,i1
          if (q_hr_mask(i,j,k)) then
            ! limiting procedure (as per S&B)
            x_hr (i,j,k)   = q_hr(i,j,k)/(n_hr(i,j,k)+eps0) ! JvdD Added eps0 to avoid floating point exception
            xr_try         = max(xrmin,min(xrmax,x_hr(i,j,k)))  !
            x_hr (i,j,k)   = max(xrmin,min(xrmax,x_hr(i,j,k)))  ! same as above, for cold mphys
            !
            D_hr(i,j,k) = a_hr *x_hr(i,j,k)**b_hr
            v_hr(i,j,k) = al_hr*x_hr(i,j,k)**be_hr*(rho0s/rhof(k))**ga_hr
            !
            Dvr(i,j,k)     = (xr_try/pirhow)**(1./3.)
            N_r0_try       = rhof(k)*n_hr(i,j,k)/Dvr(i,j,k)
            N_r0(i,j,k)    = max(N_0min,min(N_0max,N_r0_try))

            lbdr_try(i,j,k)= (pirhow*N_r0(i,j,k)/(rhof(k)*q_hr(i,j,k)))**0.25  ! c_lbdr*x_hr(i,j,k)**(-mu_hr_cst)
            lbdr(i,j,k)    = max(lbdr_min, min(lbdr_max,lbdr_try(i,j,k)))
          endif
        enddo
        enddo
        enddo

      else ! l_sb_classic

        do k=1,k1
        do j=2,j1
        do i=2,i1
           ! #sb3 changing the variable name
           if (q_hr_mask(i,j,k)) then
             x_hr (i,j,k) = q_hr(i,j,k)/(n_hr(i,j,k)+eps0)
             ! TODO: limit x_hr?
             xr = min(max(x_hr(i,j,k),xrmin),xrmax) ! to ensure xr is within borders
             Dvr(i,j,k) = (xr/pirhow)**(1./3.)
             !
             D_hr(i,j,k) = a_hr *x_hr(i,j,k)**b_hr
             v_hr(i,j,k) = al_hr*x_hr(i,j,k)**be_hr*(rho0s/rhof(k))**ga_hr
           endif
        enddo
        enddo
        enddo

        if (l_mur_cst) then
          ! mur = cst
          do k=1,k1
          do j=2,j1
          do i=2,i1
            if (q_hr_mask(i,j,k)) then
              mur(i,j,k) = mur_cst
              lbdr(i,j,k) = ((mur(i,j,k)+3.)*(mur(i,j,k)+2.)*(mur(i,j,k)+1.))**(1./3.)/Dvr(i,j,k)
            endif
          enddo
          enddo
          enddo
        else
          ! mur = f(Dv)
          do k=1,k1
          do j=2,j1
          do i=2,i1
            if (q_hr_mask(i,j,k)) then
              mur(i,j,k) = min(mur0_G09b,- 1. + c_G09b/ (q_hr(i,j,k)*rhof(k))**exp_G09b)  ! G09b
              lbdr(i,j,k) = ((mur(i,j,k)+3.)*(mur(i,j,k)+2.)*(mur(i,j,k)+1.))**(1./3.)/Dvr(i,j,k)
            endif
          enddo
          enddo
          enddo
        endif
      endif !  l_sb_classic
    else !  l_sb
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (q_hr_mask(i,j,k).and.n_hr(i,j,k).gt.0.) then
          x_hr (i,j,k) = q_hr(i,j,k)/(n_hr(i,j,k)+eps0)
          ! TODO: limit x_hr?
          Dvr  (i,j,k) = (x_hr(i,j,k)/pirhow)**(1./3.)
          !
          D_hr(i,j,k) = a_hr *x_hr(i,j,k)**b_hr
          v_hr(i,j,k) = al_hr*x_hr(i,j,k)**be_hr*(rho0s/rhof(k))**ga_hr
        endif
      enddo
      enddo
      enddo
    endif ! l_sb
  endif   ! l_rain

  ! cloud water
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_cl_mask(i,j,k)) then
      x_cl (i,j,k) = q_cl(i,j,k)/(n_cl(i,j,k)+eps0)
      ! xc (i,j,k) = q_cl(i,j,k)/(n_cl(i,j,k)+eps0)
      xc (i,j,k) = min(max(x_cl(i,j,k),x_cl_bmin),x_cl_bmax) ! to ensure x is within borders
      x_cl (i,j,k) = min(max(x_cl(i,j,k),x_cl_bmin),x_cl_bmax) ! as well - to limit extrme autoconversion
      !
      D_cl(i,j,k) = a_cl *x_cl(i,j,k)**b_cl
      v_cl(i,j,k) = al_cl*x_cl(i,j,k)**be_cl*(rho0s/rhof(k))**ga_cl
    endif
  enddo
  enddo
  enddo

  ! cloud ice
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_ci_mask(i,j,k)) then
      x_ci (i,j,k) = q_ci(i,j,k)/(n_ci(i,j,k)+eps0)
      x_ci (i,j,k) = min(max(x_ci(i,j,k),x_ci_bmin),x_ci_bmax) ! to ensure x is within borders
      !
      D_ci(i,j,k) = a_ci *x_ci(i,j,k)**b_ci
      v_ci(i,j,k) = al_ci*x_ci(i,j,k)**be_ci*(rho0s/rhof(k))**ga_ci
    endif
  enddo
  enddo
  enddo

  ! snow
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_hs_mask(i,j,k)) then
      x_hs (i,j,k) = q_hs(i,j,k)/(n_hs(i,j,k)+eps0)
      x_hs (i,j,k) = min(max(x_hs(i,j,k),x_hs_bmin),x_hs_bmax) ! to ensure x is within borders
      !
      D_hs(i,j,k) = a_hs *x_hs(i,j,k)**b_hs
      v_hs(i,j,k) = al_hs*x_hs(i,j,k)**be_hs*(rho0s/rhof(k))**ga_hs
    endif
  enddo
  enddo
  enddo

  ! graupel
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_hg_mask(i,j,k)) then
      x_hg (i,j,k) = q_hg(i,j,k)/(n_hg(i,j,k)+eps0)
      x_hg (i,j,k) = min(max(x_hg(i,j,k),x_hg_bmin),x_hg_bmax) ! to ensure x is within borders
      !
      D_hg(i,j,k) = a_hg *x_hg(i,j,k)**b_hg
      v_hg(i,j,k) = al_hg*x_hg(i,j,k)**be_hg*(rho0s/rhof(k))**ga_hg
    endif
  enddo
  enddo
  enddo

  deallocate ( N_r0, lbdr_try)

end subroutine integrals_bulk3


!> Determine autoconversion rate and adjust qrp and Nrp accordingly
!!
!!   The autoconversion rate is formulated for f(x)=A*x**(nuc)*exp(-Bx),
!!   decaying exponentially for droplet mass x.
!!   It can easily be reformulated for f(x)=A*x**(nuc)*exp(-Bx**(mu)) and
!!   by chosing mu=1/3 one would get a gamma distribution in drop diameter
!!   -> faster rain formation. (Seifert)
subroutine autoconversion3
  use modglobal, only : i1,j1,k1,kmax,rlv,cp
  use modmpi,    only : myid
  use modfields, only : exnf,rhof,ql0, svm
  implicit none

  integer :: i,j,k
  real :: rem_cf

  if (l_sb) then
    !
    ! SB autoconversion
    !
    tau(2:i1,2:j1,1:k1) = 0.
    phi(2:i1,2:j1,1:k1) = 0.
    nuc(2:i1,2:j1,1:k1) = 0.
    ! calculating the constant in the front part of the term

    if (l_sb_classic) then  ! l_sb_classic - ie. S&B version

      ! autoconversion coefficient
      k_au = k_cc/(20.0*x_s)
      ! remain coefficient
      rem_cf = (1.0-rem_n_cl_min)/delt

      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (q_cl_mask(i,j,k)) then
          dq_hr_au(i,j,k) = k_au * ( nu_cl_cst+2.) * ( nu_cl_cst +4.) / ( nu_cl_cst +1.)**2.    &
                  * (q_cl(i,j,k) * x_cl(i,j,k))**2. * rho0s ! *rho**2/rho/rho (= 1)
          tau    (i,j,k) = 1.0 - q_cl(i,j,k) / qltot(i,j,k)

          ! phi_au computation
          phi    (i,j,k) = k_1 * tau(i,j,k)**k_2 * (1.0 -tau(i,j,k)**k_2)**3
          dq_hr_au(i,j,k) = dq_hr_au(i,j,k) * (1.0 + phi(i,j,k)/(1.0 -tau(i,j,k))**2)

          ! basic au correction
          dq_hr_au(i,j,k)  = min(dq_hr_au(i,j,k),max(0.0,svm(i,j,k,iq_cl)/delt+q_clp(i,j,k)))

          ! and calculate cloud droplet numbers
          dn_cl_au(i,j,k) = (-2.0/x_s)*dq_hr_au(i,j,k)                                          ! and the droplet number
          dn_cl_au(i,j,k) = max(dn_cl_au(i,j,k),min(0.0,-rem_cf*svm(i,j,k,in_cl)-n_clp(i,j,k))) ! correct droplet number
          dn_hr_au(i,j,k) = -0.5*dn_cl_au(i,j,k)                                                ! and the raindrop number

          ! and adding updates
          q_hrp (i,j,k) = q_hrp (i,j,k) + dq_hr_au (i,j,k)
          n_hrp (i,j,k) = n_hrp (i,j,k) + dn_hr_au(i,j,k)
          q_clp (i,j,k) = q_clp (i,j,k) - dq_hr_au (i,j,k)
          n_clp (i,j,k) = n_clp (i,j,k) + dn_cl_au(i,j,k)! #sb3

          thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*dq_hr_au(i,j,k)
          qtpmcr (i,j,k) = qtpmcr (i,j,k) - dq_hr_au (i,j,k)
        endif
      enddo
      enddo
      enddo
    else ! l_sb_classic = .false. - original version in DALES

      k_au = k_c/(20.0*x_s)

      do k=1,k1
      do j=2,j1
      do i=2,i1
         if (q_cl_mask(i,j,k)) then
            nuc    (i,j,k) = k1nuc*(rhof(k)*q_cl(i,j,k)*1000.) +k2nuc-1. ! #sb3 G09a
            dq_hr_au(i,j,k) = k_au * (nuc(i,j,k)+2.) * (nuc(i,j,k)+4.) / (nuc(i,j,k)+1.)**2.    &
                    * (q_cl(i,j,k) * x_cl(i,j,k))**2. * rho0s ! *rho**2/rho/rho (= 1)
            tau    (i,j,k) = 1.0 - q_cl(i,j,k) / qltot(i,j,k) ! #sb3
            phi    (i,j,k) = k_1 * tau(i,j,k)**k_2 * (1.0 -tau(i,j,k)**k_2)**3   ! phi_au computation
            dq_hr_au(i,j,k) = dq_hr_au(i,j,k) * (1.0 + phi(i,j,k)/(1.0 -tau(i,j,k))**2)

            ! cloud water numbers
            dn_hr_au(i,j,k) = dq_hr_au(i,j,k)/x_s
            dn_cl_au(i,j,k) =  (-2.0/x_s)*dq_hr_au(i,j,k)

            ! #sb3 START outputs
            q_hrp    (i,j,k) = q_hrp    (i,j,k) + dq_hr_au (i,j,k)
            n_hrp    (i,j,k) = n_hrp    (i,j,k) + dn_hr_au(i,j,k)
            q_clp (i,j,k) = q_clp (i,j,k) - dq_hr_au (i,j,k)
            n_clp (i,j,k) = n_clp (i,j,k) + dn_cl_au(i,j,k) ! o:  n_clp (i,j,k) - (2.0/x_s)*rhof(k)*au(i,j,k)

            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*dq_hr_au(i,j,k)
            qtpmcr (i,j,k) = qtpmcr (i,j,k) - dq_hr_au (i,j,k)
         endif
      enddo
      enddo
      enddo
    endif ! l_sb_classic

  endif ! l_sb

  if (l_sb_dbg) then
    if (any(q_cl(2:i1,2:j1,1:kmax)/delt - dq_hr_au(2:i1,2:j1,1:kmax) .lt. 0.)) then
      write(6,*) 'WARNING: autoconversion too high'
      write(6,*) '  removing more cloud water than available in ', count(q_cl(2:i1,2:j1,1:kmax)/delt - dq_hr_au(2:i1,2:j1,1:kmax) .lt. 0.)
    end if

    if (any(n_cl(2:i1,2:j1,1:kmax)/delt -(2.0/x_s)*dq_hr_au(2:i1,2:j1,1:kmax) .lt. 0.)) then
      write(6,*) 'WARNING: autoconversion too high'
      write(6,*) '  removing more droplets than available in' , count(n_cl(2:i1,2:j1,1:kmax)/delt - (2.0/x_s)*dq_hr_au(2:i1,2:j1,1:kmax) .lt. 0.)
    end if
  endif

end subroutine autoconversion3


subroutine accretion3
!*********************************************************************
! determine accr. + self coll. + br-up rate and adjust qrp and Nrp
! base don S&B
!*********************************************************************
  use modglobal, only : ih,i1,jh,j1,k1,kmax,rlv,cp
  use modfields, only : exnf,rhof,ql0,qt0, svm, qvsl, qvsi
  use modmpi,    only : myid
  implicit none

  real , allocatable :: phi_br(:,:,:), Dvrf(:,:,:)
  real :: rem_cf, rem_clcf
  integer :: i,j,k

  allocate ( phi_br(2-ih:i1+ih,2-jh:j1+jh,k1)      &
            ,Dvrf  (2-ih:i1+ih,2-jh:j1+jh,k1)      &  ! #sb3
           )

  phi_br=0.0   ! #sb3
  Dvrf =0.0  ! #sb3

  ! calculating coef for ratio for minimal remaing number of droplets #sb3
  rem_cf = (1.0-rem_n_hr_min)/delt
  ! remain coefficient for clouds #sb3
  rem_clcf = (1.0-rem_n_cl_min)/delt

  if (l_sb) then
    !
    ! SB accretion
    !
    ! #t write(6,*) 'l_sb accretion...'
    if (l_sb_classic) then
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (q_hr_mask(i,j,k) .and. q_cl_mask(i,j,k)) then
          ! since it is forming only where rain
          tau    (i,j,k) = 1.0 - q_cl(i,j,k)/(qltot(i,j,k))
          phi    (i,j,k) = (tau(i,j,k)/(tau(i,j,k) + k_l))**4.
          dq_hr_ac(i,j,k) = k_cr *rhof(k)*q_cl(i,j,k) * q_hr(i,j,k) * phi(i,j,k) * &
                           (rho0s/rhof(k))**0.5  ! rho*rho / rho  = rho

          ! basic ac correction
          dq_hr_ac     (i,j,k) = min(dq_hr_ac(i,j,k),q_cl(i,j,k)/delt) ! min(dq_hr_ac(i,j,k),svm(i,j,k,iq_cl)/delt)\

          ! number of cloud droplets
          dn_cl_ac(i,j,k) = -dq_hr_ac(i,j,k)/x_cl(i,j,k)
          dn_cl_ac(i,j,k) = max(dn_cl_ac(i,j,k),min(0.0,-rem_clcf*svm(i,j,k,in_cl)-n_clp(i,j,k)))

          ! update
          q_hrp  (i,j,k) = q_hrp (i,j,k)  + dq_hr_ac(i,j,k)

          ! no change n_hrp
          ! and changes in water number for clouds
          q_clp (i,j,k) = q_clp (i,j,k) - dq_hr_ac(i,j,k)
          n_clp (i,j,k) = n_clp (i,j,k) + dn_cl_ac(i,j,k)  !o: n_clp (i,j,k) - rhof(k)*ac(i,j,k)/xc(i,j,k)

          qtpmcr (i,j,k) = qtpmcr (i,j,k) - dq_hr_ac(i,j,k)
          thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*dq_hr_ac(i,j,k)
        endif
      enddo
      enddo
      enddo
   else ! l_sb_classic
     !*********************************************************************
     ! determine accr. rate and adjust qrp and Nrp
     ! accordingly. Break-up : Seifert (2007), a
     !*********************************************************************
     tau    (i,j,k) = 1.0 - ql0(i,j,k)/(qltot(i,j,k))
     phi    (i,j,k) = (tau(i,j,k)/(tau(i,j,k) + k_l))**4.

     dq_hr_ac     (i,j,k) = k_r *rhof(k)*ql0(i,j,k) * q_hr(i,j,k) * phi(i,j,k) * &
                      (1.225/rhof(k))**0.5
     q_hrp  (i,j,k) = q_hrp (i,j,k) + dq_hr_ac(i,j,k)

     ! no change n_hrp
     q_clp  (i,j,k) = q_clp (i,j,k) - dq_hr_ac(i,j,k)
     qtpmcr (i,j,k) = qtpmcr (i,j,k) - dq_hr_ac(i,j,k)
     thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*dq_hr_ac(i,j,k)

     ! and update on cloud water number
     dn_cl_ac(i,j,k)= -dq_hr_ac(i,j,k)/x_cl(i,j,k)
     n_clp (i,j,k)  = n_clp(i,j,k) + dn_cl_ac(i,j,k)
   endif

   !
   ! SB self-collection & Break-up
   !
   if (l_sb_classic) then
     do k=1,k1
     do j=2,j1
     do i=2,i1
       if (q_hr_mask(i,j,k)) then
          dn_hr_sc(i,j,k) = -k_rr *rhof(k)* q_hr(i,j,k) * n_hr(i,j,k)  &
            * (1.0 + kappa_r/lbdr(i,j,k))**(-9.)*(rho0s/rhof(k))**0.5

          ! and calculating size of droplets - adjusted
          Dvrf(i,j,k)= Dvr(i,j,k) ! for now leaving the same
          if (Dvrf(i,j,k) .gt. dvrlim) then
            if (Dvrf(i,j,k) .gt. dvrbiglim) then
              ! for big drops
              phi_br(i,j,k) = 2.0*exp(kappa_br*(Dvrf(i,j,k)-D_eq))-1.0
            else
              ! for smaller drops
              phi_br(i,j,k) = k_br * (Dvrf(i,j,k)-D_eq)
            endif
            dn_hr_br(i,j,k) = -(phi_br(i,j,k) + 1.) * dn_hr_sc(i,j,k)
          else
            dn_hr_br(i,j,k) = 0. ! (phi_br = -1)
          endif
       endif ! q_hr_mask
     enddo
     enddo
     enddo
   else ! l_sb_classic
     do k=1,k1
     do j=2,j1
     do i=2,i1
       ! adjusting variable names in accretion
       if (q_hr_mask(i,j,k)) then
          dn_hr_sc(i,j,k) = -k_rr *rhof(k)* q_hr(i,j,k) * n_hr(i,j,k)  &
            * (1.0 + kappa_r/lbdr(i,j,k)*pirhow**(1./3.))**(-9.)*(rho0s/rhof(k))**0.5
       endif
       if ((Dvr(i,j,k) .gt. dvrlim) .and. q_hr_mask(i,j,k)) then
         if (Dvr(i,j,k) .gt. dvrbiglim) then
           ! for big drops
           phi_br(i,j,k) = 2.0*exp(kappa_br*(Dvr(i,j,k)-D_eq))-1.0
         else
           ! for smaller drops
           phi_br(i,j,k) = k_br* (Dvr(i,j,k)-D_eq)
         endif
         dn_hr_br(i,j,k) = -(phi_br(i,j,k) + 1.) * dn_hr_sc(i,j,k)
       else
          dn_hr_br(i,j,k) = 0. ! (phi_br = -1)
       endif
     enddo
     enddo
     enddo
   endif ! l_sb_classic
  endif

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_hr_mask(i,j,k)) then
     ! correcting for low values
     n_hrp(i,j,k) = n_hrp(i,j,k)                         &
        +max(min(0.0,-rem_cf*svm(i,j,k,in_hr)-n_hrp(i,j,k)),dn_hr_sc(i,j,k)+dn_hr_br(i,j,k))
    endif
  enddo
  enddo
  enddo

  if (l_sb_dbg) then
    if( any((svm(2:i1,2:j1,1:k1,iq_cl)/delt-dq_hr_ac(2:i1,2:j1,1:k1)).lt. 0) ) then
      write(6,*) 'WARNING: accretion removing too much water'
      write(6,*) ' updated below 0 in gridpoints', count((svm(2:i1,2:j1,1:k1,iq_cl)/delt- dq_hr_ac(2:i1,2:j1,1:k1)).lt. 0.0)
    endif

    if( any((svm(2:i1,2:j1,1:k1,in_hr)/delt+dn_hr_sc(2:i1,2:j1,1:k1)+dn_hr_br(2:i1,2:j1,1:k1)).lt. 0) ) then
      write(6,*) 'WARNING: self-collection of rain too high'
      write(6,*) ' removing more n_hr than available in gridpoints', count((svm(2:i1,2:j1,1:k1,in_hr)/delt-dn_hr_sc(2:i1,2:j1,1:k1)+dn_hr_br(2:i1,2:j1,1:k1)).lt. 0.0 )
      write(6,*) ' n_hr updated below 0 in gridpoints', count((n_hr(2:i1,2:j1,1:k1)/delt-dn_hr_sc(2:i1,2:j1,1:k1)+dn_hr_br(2:i1,2:j1,1:k1)).lt. 0.0 )
    endif

    if (any((svm(2:i1,2:j1,1:kmax,in_cl)/delt - dq_hr_ac(2:i1,2:j1,1:kmax)/x_cl(i,j,k)).lt. 0.)) then
      write(6,*)'WARNING: ac too large, removing too many droplets'
      write(6,*)'   in', count((svm(2:i1,2:j1,1:kmax,in_cl)/delt    &
         - dq_hr_ac(2:i1,2:j1,1:kmax)/x_cl(i,j,k)).lt. 0.),myid
    endif
  endif ! l_sb_dbg

  deallocate (phi_br,Dvrf)

end subroutine accretion3


!> Self-collection of cloud droplets
!! written based on S&B to evaluate the self-collection of cloud droplets
subroutine cloud_self3
  use modglobal, only : ih,i1,jh,j1,k1,kmax,rlv,cp
  use modfields, only : exnf,rhof, svm
  use modmpi,    only : myid
  implicit none

  real    ::  k_clsc, rem_cf
  integer :: i,j,k

  ! initialize
  dn_cl_sc = 0.0

  ! calculate constant
  k_clsc = - k_cc*rho0s
  rem_cf = (1.0-rem_n_cl_min)/delt

  ! loop sb
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_cl_mask(i,j,k)) then
       dn_cl_sc(i,j,k) = k_clsc*(q_cl(i,j,k)**2)* &
          (nu_cl_cst+2.0)/(nu_cl_cst+1.0)-dn_cl_au(i,j,k)

       ! basic sc collection
       dn_cl_sc(i,j,k) = min(0.0,dn_cl_sc(i,j,k))
       dn_cl_sc(i,j,k) = max(dn_cl_sc(i,j,k),min(0.0,-rem_cf*svm(i,j,k,in_cl)-n_clp(i,j,k)))

       ! update
       n_clp(i,j,k) = n_clp(i,j,k)+dn_cl_sc(i,j,k)

       ! no change to q_cl
    endif
  enddo
  enddo
  enddo

  if (l_sb_dbg) then
    ! testing for too high values
    if (any(q_cl_mask(2:i1,2:j1,1:k1).and.((svm(2:i1,2:j1,1:k1,in_cl)/delt - dn_cl_sc(2:i1,2:j1,1:k1)).lt. 0.))) then
      write(6,*)'WARNING: cloud sc too high'
      write(6,*) '  getting to negative n_cl  in gridpoints  ',count((svm(2:i1,2:j1,1:k1,in_cl)/delt -dn_cl_sc(2:i1,2:j1,1:k1)).lt. 0.)
    endif
  endif

end subroutine cloud_self3

!> Cloud nucleation
!! Written to prognostically evaluate the cloud water number content [ kg^{-1}]
!! directly follows Seifert&Beheng scheme
subroutine nucleation3
  use modglobal, only : dzf,i1,j1,ih,jh,k1,kmax,rlv,cp
  use modfields, only : w0, rhof,exnf, qvsl, qt0,ql0,svm
  use modmpi,    only : myid

  implicit none
  integer :: i,j,k
  real  :: coef_ccn, n_act

  ! allocatable variables - supersaturation, derivation of supersaturation
  real, allocatable :: ssat(:,:,:), wdssatdz(:,:,:) !, dn_cl_nu(:,:,:)

  ! allocate
  allocate(  ssat (2-ih:i1+ih,2-jh:j1+jh,k1)      &
            ,wdssatdz (2-ih:i1+ih,2-jh:j1+jh,k1)  &
           )

  ssat = 0.0  ! note that supersaturation is
              ! not always how supersaturated is water vapour,
              ! depending on a flag, it can also include water already in droplets

  wdssatdz = 0.0
  dn_cl_nu = 0.0
  coef_ccn  = 1.0/sat_max**kappa_ccn ! 1.0
  ! allows to keep both definitions consistent

  ! calculating supersaturation:
  if (l_sb_nuc_sat) then
    ! loops to determine supersaturation
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! calculating supersaturation of water vapour only
      ssat(i,j,k)=(100./qvsl(i,j,k))*(qt0(i,j,k)-q_cl(i,j,k)-qvsl(i,j,k))
    enddo
    enddo
    enddo
  else ! l_sb_nuc_sat
    ! ie. cloud liquid water is also included supersaturation
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ssat(i,j,k) = (100./qvsl(i,j,k))*(qt0(i,j,k)-qvsl(i,j,k))
    enddo
    enddo
    enddo
  endif ! l_sb_nuc_sat

  if (l_sb_nuc_expl.and.l_sb_nuc_diff) then
    ! calculating the derivation - second order estimation?
    ! ? add switches for different derivation calculation? - so first the firt order only
    ! option A: central differences lax
    do k=2,k1-1
    do j=2,j1
    do i=2,i1
      ! if in cloud
      if ( ssat(i,j,k).gt.0.0 ) then
        ! calculating the derivation
        wdssatdz(i,j,k) = 0.5*(w0(i,j,k+1)+ w0(i,j,k))*(ssat(i,j,k+1)-ssat(i,j,k-1))/(dzf(k)+dzf(k-1))
      endif
    enddo
    enddo
    enddo

    ! now for separate cases k=1 and k = k1
    k=1
    do j=2,j1
    do i=2,i1
      if (ssat(i,j,k).gt.0.0) then
        ! just an approximation - same as for the second level
        wdssatdz(i,j,k) = wdssatdz(i,j,k+1)  ! or 0 to prevent condensation there
      endif
    enddo
    enddo

    k=k1
    do j=2,j1
    do i=2,i1
      if (ssat(i,j,k).gt.0.0) then
        ! first order approximation of the difference
        wdssatdz(i,j,k) = w0(i,j,k)* (ssat(i,j,k)-ssat(i,j,k-1)) /dzf(k-1)
      endif
    enddo
    enddo
  endif ! l_sb_nuc_expl.and.l_sb_nuc_diff

  ! calculation of explicit nucleation
  if (l_sb_nuc_expl) then
    if (l_c_ccn) then ! c_ccn is constant
      if (l_sb_nuc_diff) then
        if (l_sb_sat_max) then
          do k=1,k1
          do j=2,j1
          do i=2,i1
            ! of course only in cloud
            if(ssat(i,j,k).gt.sat_min) then
              ! if conditions for nucleation, calculate it
              if ((ssat(i,j,k).lt.sat_max).and.(wdssatdz(i,j,k).gt.0.0)) then
                dn_cl_nu(i,j,k)=c_ccn*kappa_ccn*wdssatdz(i,j,k)*ssat(i,j,k)**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
                ! note - written this way on purpose
                ! in case of further adjustment of nucleation subroutin
              endif
              ! basic limiting
              dn_cl_nu(i,j,k)=max(0.0,min(dn_cl_nu(i,j,k),(n_clmax -n_cl(i,j,k))/delt))
            endif
          enddo
          enddo
          enddo
        else ! NOT l_sb_sat_max  AND l_sb_nuc_diff
          do k=1,k1
          do j=2,j1
          do i=2,i1
            ! of course only in cloud
            if( ssat(i,j,k).gt.sat_min) then
              ! if conditions for nucleation, calculate it
              if (wdssatdz(i,j,k).gt.0.0) then
                dn_cl_nu(i,j,k)=c_ccn*kappa_ccn*wdssatdz(i,j,k)*ssat(i,j,k)**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
                ! note - written this way on purpose
                ! in case of further adjustment of nucleation subroutin
                ! basic limiting
                dn_cl_nu(i,j,k)=max(0.0,min(dn_cl_nu(i,j,k),(n_clmax -n_cl(i,j,k))/delt))
              endif
            endif
          enddo
          enddo
          enddo
        endif   ! l_sb_sat_max
      else ! NOT l_sb_nuc_diff
        if (l_sb_sat_max)  then ! l_sb_sat_max AND NOT l_sb_nuc_diff
          do k=1,k1
          do j=2,j1
          do i=2,i1
            ! of course only in cloud
            if(ssat(i,j,k).gt.sat_min) then
              ! if conditions for nucleation, calculate it
              if ((ssat(i,j,k).lt.sat_max).and.(w0(i,j,k).gt.0.0)) then
                dn_cl_nu(i,j,k)=(c_ccn/dzf(k-1))*max(0.0,w0(i,j,k)*      &
                   (ssat(i,j,k)**kappa_ccn-ssat(i,j,k-1)**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if ((ssat(i,j,k).lt.sat_max).and.(w0(i,j,k+1).lt.0.0)) then
                dn_cl_nu(i,j,k)=(c_ccn/dzf(k-1))*max(0.0,w0(i,j,k+1)*    &
                   (ssat(i,j,k)**kappa_ccn-ssat(i,j,k+1)**kappa_ccn))   ! (1/rhof) *rhof = 1
               endif
               ! basic limiting
               dn_cl_nu(i,j,k)=max(0.0,min(dn_cl_nu(i,j,k),(n_clmax -n_cl(i,j,k))/delt))
            endif
          enddo
          enddo
          enddo
        else ! NOT l_sb_sat_max AND NOT l_sb_nuc_diff
          do k=1,k1
          do j=2,j1
          do i=2,i1
            ! of course only in cloud
            if( ssat(i,j,k).gt.sat_min) then
              ! if conditions for nucleation, calculate it
              if (w0(i,j,k).gt.0.0) then
                dn_cl_nu(i,j,k)=(c_ccn/dzf(k-1))*max(0.0,w0(i,j,k)*      &
                   (ssat(i,j,k)**kappa_ccn-ssat(i,j,k-1)**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if (w0(i,j,k+1).lt.0.0) then
                dn_cl_nu(i,j,k)=(c_ccn/dzf(k-1))*max(0.0,w0(i,j,k+1)*    &
                   (ssat(i,j,k)**kappa_ccn-ssat(i,j,k+1)**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              ! basic limiting
              dn_cl_nu(i,j,k)=max(0.0,min(dn_cl_nu(i,j,k),(n_clmax -n_cl(i,j,k))/delt))
            endif
          enddo
          enddo
          enddo
        endif
      endif  ! NOT l_sb_sat_max  AND NOT l_sb_nuc_diff
    else ! l_c_ccn, ie. c_ccn is not constant, bun dependent on n_cc
      if (l_sb_nuc_diff) then   ! l_sb_nuc_diff
        if (l_sb_sat_max ) then ! l_sb_sat_max
          do k=1,k1
          do j=2,j1
          do i=2,i1
            ! of course only in cloud
            if( ssat(i,j,k).gt.sat_min) then
              ! if conditions for nucleation, calculate it
              if ((ssat(i,j,k).lt.sat_max).and.(wdssatdz(i,j,k).gt.0.0)) then
                dn_cl_nu(i,j,k)=coef_ccn*n_cc(i,j,k)*kappa_ccn*wdssatdz(i,j,k)*ssat(i,j,k)**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
                ! basic limiting
                dn_cl_nu(i,j,k)=max(0.0,min(dn_cl_nu(i,j,k),(n_cc(i,j,k)-n_cl(i,j,k))/delt))
              endif
            endif
          enddo
          enddo
          enddo
        else  ! l_sb_sat_max
          do k=1,k1
          do j=2,j1
          do i=2,i1
            ! of course only in cloud
            if( ssat(i,j,k).gt.sat_min) then
              ! if conditions for nucleation, calculate it
              if (wdssatdz(i,j,k).gt.0.0) then
                dn_cl_nu(i,j,k)=coef_ccn*n_cc(i,j,k)*kappa_ccn*wdssatdz(i,j,k)*ssat(i,j,k)**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
                ! basic limiting
                dn_cl_nu(i,j,k)=max(0.0,min(dn_cl_nu(i,j,k),(n_cc(i,j,k)-n_cl(i,j,k))/delt))
              endif
            endif
          enddo
          enddo
          enddo
        endif  ! l_sb_sat_max
      else  ! l_sb_nuc_diff
        if (l_sb_sat_max)  then ! l_sb_sat_max  AND NOT l_sb_nuc_diff
          do k=1,k1
          do j=2,j1
          do i=2,i1
            ! of course only in cloud
            if(ssat(i,j,k).gt.sat_min) then
              ! if conditions for nucleation, calculate it
              if ((ssat(i,j,k).lt.sat_max).and.(w0(i,j,k).gt.0.0)) then
                dn_cl_nu(i,j,k)=(coef_ccn*n_cc(i,j,k)/dzf(k-1))*        &
                   max(0.0,w0(i,j,k)*                                   &
                  (ssat(i,j,k)**kappa_ccn-ssat(i,j,k-1)**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if ((ssat(i,j,k).lt.sat_max).and.(w0(i,j,k+1).lt.0.0)) then
                dn_cl_nu(i,j,k)=(coef_ccn*n_cc(i,j,k)/dzf(k-1))*         &
                   max(0.0,w0(i,j,k+1)*                                  &
                   (ssat(i,j,k)**kappa_ccn-ssat(i,j,k+1)**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              ! basic limiting
              dn_cl_nu(i,j,k)=max(0.0,min(dn_cl_nu(i,j,k),(n_cc(i,j,k)-n_cl(i,j,k))/delt))
            endif
          enddo
          enddo
          enddo
        else ! NOT l_sb_sat_max  AND NOT l_sb_nuc_diff
          do k=1,k1
          do j=2,j1
          do i=2,i1
            if(ssat(i,j,k).gt.sat_min) then
              ! if conditions for nucleation, calculate it
              if (w0(i,j,k).gt.0.0) then
                dn_cl_nu(i,j,k)=(coef_ccn*n_cc(i,j,k)/dzf(k-1))*         &
                    max(0.0,w0(i,j,k)*                                   &
                    (ssat(i,j,k)**kappa_ccn-ssat(i,j,k-1)**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if (w0(i,j,k+1).lt.0.0) then
                dn_cl_nu(i,j,k)=(coef_ccn*n_cc(i,j,k)/dzf(k-1))*         &
                   max(0.0,w0(i,j,k+1)*                                  &
                   (ssat(i,j,k)**kappa_ccn-ssat(i,j,k+1)**kappa_ccn))    ! (1/rhof) *rhof = 1
              endif
              ! basic limiting
              dn_cl_nu(i,j,k)=max(0.0,min(dn_cl_nu(i,j,k),(n_cc(i,j,k)-n_cl(i,j,k))/delt))
            endif
          enddo
          enddo
          enddo
        endif  ! l_sb_sat_max
      endif  ! l_sb_nuc_diff
    endif ! l_c_ccn
  else ! l_sb_nuc_expl
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! of course only in cloud
      if( ssat(i,j,k).gt.0.0) then
        ! calculate number of activated n_ccn
        n_act = coef_ccn*n_cc(i,j,k)*min(sat_max,ssat(i,j,k))**kappa_ccn
        n_act = max(n_cc(i,j,k),n_act)
        !
        if (n_act.gt.n_cl(i,j,k)) then
          dn_cl_nu(i,j,k) = (n_act-n_cl(i,j,k))/delt
        else
          dn_cl_nu(i,j,k) = 0.0
        endif
        ! basic limiting - not needed in this case
        ! dn_cl_nu(i,j,k)=min(dn_cl_nu(i,j,k), (n_cc(i,j,k)-n_cl(i,j,k))/delt)
      endif
    enddo
    enddo
    enddo
  endif ! l_sb_nuc_expl

  ! update itself
  do k=1,k1
  do j=2,j1
  do i=2,i1
     ! basic limiting
     n_clp(i,j,k) = n_clp(i,j,k)+dn_cl_nu(i,j,k)  ! increase in cloud water number

     ! update water density [kg kg^{-1}]
     q_clp (i,j,k) = q_clp (i,j,k) + x_cnuc*dn_cl_nu(i,j,k)
  enddo
  enddo
  enddo

  if (l_sb_dbg) then
    ! warning if too high liquid water
    if (any((ssat(i,j,k).gt.0.0).and.((qt0(2:i1,2:j1,1:kmax)-qvsl(2:i1,2:j1,1:kmax)-svm(2:i1,2:j1,1:k1,iq_cl)-svm(2:i1,2:j1,1:k1,iq_ci))/delt - x_cnuc*dn_cl_nu(2:i1,2:j1,1:k1)).lt. 0.)) then
      write(6,*) 'WARNING: cloud nucleation too high'
      write(6,*) ' removing too much water in gridpoints ',count((ssat(i,j,k).gt.0.0)  &
             .and.((qt0(2:i1,2:j1,1:kmax)-qvsl(2:i1,2:j1,1:kmax)-svm(2:i1,2:j1,1:k1,iq_cl)-svm(2:i1,2:j1,1:k1,iq_ci))/delt - x_cnuc*dn_cl_nu(2:i1,2:j1,1:k1)).lt. 0.)
    end if
  endif !l_sb_dbg

  ! deallocate variables
  deallocate(ssat,wdssatdz)

end subroutine  nucleation3

! ****************************************************************
!  Limiting cloud dropler and cloud ice nucleation
!  - to prevent too low ccp
subroutine cor_nucl3
  use modglobal, only : dzf,i1,j1,ih,jh,k1,kmax,rlv, cp
  use modfields, only : rhof,exnf, qt0, svm
  use modmpi,    only : myid
  implicit none

  integer :: i,j,k
  real :: coef_mincc, cor_coef

  ! allocatable varibles - supersaturation, derivation of supersaturation
  real, allocatable :: cor_dn_cl(:,:,:), cor_dn_ci(:,:,:), n_ccm(:,:,:)

  ! ------
  ! allocate
  allocate(  cor_dn_cl (2-ih:i1+ih,2-jh:j1+jh,k1)           &
            ,cor_dn_ci (2-ih:i1+ih,2-jh:j1+jh,k1)           &
            ,n_ccm     (2-ih:i1+ih,2-jh:j1+jh,k1)           &
          )

  ! multiplicative coefficient
  coef_mincc = (1.0-cc_min_ratio)/delt

  ! insert values
  cor_dn_cl = 0.0
  cor_dn_ci = 0.0
  n_ccm(2-ih:i1+ih,2-jh:j1+jh,1:k1)=svm(2-ih:i1+ih,2-jh:j1+jh,1:k1,in_cc)

  ! test how many low values
  if(.not.l_c_ccn) then
    if ( any((coef_mincc *n_ccm(2:i1,2:j1,1:k1)-dn_cl_nu(2:i1,2:j1,1:k1)-dn_ci_inu(2:i1,2:j1,1:k1) ).lt. 0.)) then
      if(l_sb_dbg) then
        write(6,*) 'correction for n_cc, n_clp and n_cip needed in gridpoints ', count((coef_mincc   &
                  *n_ccm(2:i1,2:j1,1:k1)-dn_cl_nu(2:i1,2:j1,1:k1)-dn_ci_inu(2:i1,2:j1,1:k1)).lt. 0.)
      endif

      ! run only in case of low resulting n_cc
      do k=1,k1
      do j=2,j1
      do i=2,i1
        ! now correct those low n_ccp
        if (( coef_mincc *n_ccm(i,j,k) -dn_cl_nu(i,j,k)).lt. 0.)  then  ! #inp
          ! calculating adjustment
          cor_coef = coef_mincc*n_ccm(i,j,k)/(dn_cl_nu(i,j,k)+eps0)-1.0 ! #inp
          cor_coef = min(0.0, cor_coef)   ! so the correction does not increase number of particles
          cor_coef = max(-1.0, cor_coef)  ! so it does not send tendency to negative values

          ! calculating correctors
          cor_dn_cl(i,j,k) = cor_coef*dn_cl_nu(i,j,k)
          cor_dn_ci(i,j,k) = 0.0 ! #inp

          ! correcting
          ! adjusted nucleation rates
          dn_cl_nu(i,j,k)  = dn_cl_nu(i,j,k) + cor_dn_cl(i,j,k)
          dn_ci_inu(i,j,k) = dn_ci_inu(i,j,k)+ cor_dn_ci(i,j,k)

          ! update arrays for droplet number and ccn number
          n_clp(i,j,k) = n_clp(i,j,k)+cor_dn_cl(i,j,k)  ! i.e. decreasing cloud water number
          n_cip(i,j,k) = n_cip(i,j,k)+cor_dn_ci(i,j,k)

          ! update water density [kg kg^{-1}]
          q_clp (i,j,k) = q_clp (i,j,k) + x_cnuc*cor_dn_cl(i,j,k)
          q_cip (i,j,k) = q_cip (i,j,k) + x_inuc*cor_dn_ci(i,j,k)
          qtpmcr(i,j,k) = qtpmcr(i,j,k) - x_inuc*cor_dn_ci(i,j,k) ! #iceout

          ! update potential temperature
          thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlvi/(cp*exnf(k)))*x_inuc*cor_dn_ci(i,j,k)
        endif
      enddo
      enddo
      enddo
    endif
  endif

  ! cleaning up
  deallocate(  cor_dn_cl, cor_dn_ci, n_ccm)

end subroutine cor_nucl3


! ****************************************************************
!  Ice nucleation
! - based on Seifert & Beheng (2003), p. 53
!  ************************************************************
subroutine icenucle3
  use modglobal, only : dzf,i1,j1,k1,ih,jh,kmax,rlv,cp
  use modfields, only : w0, rhof,exnf, qvsi, qt0, ql0, svm, sv0,tmp0   ! <- later remove svm, sv0 - just for testing
  implicit none

  integer :: i,j,k

  ! allocatable varibles - supersaturation, derivation of supersaturation
  real, allocatable :: ssice(:,:,:) !,  dn_ci_inu(:,:,:)
  real:: n_in,  n_tid, dq_inuc

  ! preparing constant values
  ! rlfcp = rlfr/cp

   allocate( ssice (2-ih:i1+ih,2-jh:j1+jh,k1) ) !

   ssice = 0.0 ! not always how supersaturated is water vapour,
               ! depending on a flag,  it can also include water already in ice particles

  dn_ci_inu = 0.0

  ! calculate supersaturation with respect to ice
  if (l_sb_inuc_sat) then  ! l_sb_inuc_sat
    ! calculating supersaturation of water vapour only
    do k=1,k1
    do j=2,j1
    do i=2,i1
       ! calculating supersaturation
       ssice(i,j,k) = (qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k) -1.0
       ! ssice(i,j,k) = max(0.0, (qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k) -1.0)
    enddo
    enddo
    enddo
  else  ! l_sb_inuc_sat
    ! ie. cloud ice water is also included supersaturation
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! calculating supersaturation
      ssice(i,j,k) = (qt0(i,j,k)-q_cl(i,j,k)+q_ci(i,j,k))/qvsi(i,j,k) -1.0
      ! ssice(i,j,k) = max(0.0, (qt0(i,j,k)-q_cl(i,j,k)+q_ci(i,j,k))/qvsi(i,j,k) -1.0)
    enddo
    enddo
    enddo
  endif   ! l_sb_inuc_sat

  ! loops to calculate ice nucleation
  if (l_sb_inuc_expl) then ! l_sb_inuc_expl
    ! explicit ice nucleation - not yet included
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! -> add later
      !dn_ci_inu(i,j,k) =
    enddo
    enddo
    enddo
  else  ! l_sb_inuc_expl
    if (l_sb_reisner) then ! whether to apply reisnerr correction
      do k=1,k1
      do j=2,j1
      do i=2,i1
        ! not using qcmask condition, since air can be saturated with respect to ice
        if((tmp0(i,j,k).lt.tmp_inuc).and.(ssice(i,j,k).gt.ssice_min)) then
          ! if conditions for nucleation, calculate it
          ! Meyers et al. (1992)
          n_in = (1.0/rhof(k))*N_inuc*exp( a_M92              &
                +b_M92*min(ssice(i,j,k),ssice_lim))! o: N_inuc*exp( a_M92 + b_M92*ssice(i,j,k) )

          ! prepare Reisnerr (1998) correction
          ! w: n_tid = (1.0/rhof(k))*N_inuc_R*exp( - min(tmp0(i,j,k),c_inuc_R)-T_3)
          n_tid = (1.0/rhof(k))*N_inuc_R*exp(b_inuc_R*(T_3- max(tmp0(i,j,k),c_inuc_R)))

          ! performing reisner correction
          n_in = max(a1_inuc_R*n_tid,min(a2_inuc_R*n_tid,n_in))

          ! limiting n_in
          n_in = min(n_i_max, n_in)

          ! checking conditions if suitable for nucleation
          if (n_ci(i,j,k).lt.n_in) then ! condition intentionally left this way
            dn_ci_inu(i,j,k) = (n_in-n_ci(i,j,k))/delt
          else
            dn_ci_inu(i,j,k) = 0.0
            ! note - written this way on purpose
            ! in case of further adjustment of nucleation subroutine
          endif
        endif

        ! basic correction
        n_cip(i,j,k) = n_cip(i,j,k)+dn_ci_inu(i,j,k)  ! increase in cloud water number

        ! update water density [kg kg^{-1}]
        q_cip (i,j,k) = q_cip (i,j,k) + x_inuc*dn_ci_inu(i,j,k)

        ! update liquid water potential temperature
        !  - due to latent heat of melting (freezing in this case)
        qtpmcr(i,j,k) = qtpmcr(i,j,k) - x_inuc*dn_ci_inu(i,j,k)  ! #iceout
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlvi/(cp*exnf(k)))*x_inuc*dn_ci_inu(i,j,k)
      enddo
      enddo
      enddo
    else !l_sb_reisner
      ! i.e. when not applying reisnerr correction
      do k=1,k1
      do j=2,j1
      do i=2,i1
        ! not using qcmask condition,
        ! since air can be saturated with respect to ice
        if((tmp0(i,j,k).lt.tmp_inuc).and.(ssice(i,j,k).gt.ssice_min)) then
          ! if conditions for nucleation, calculate it
          ! Meyers et al. (1992)
          n_in = (1.0/rhof(k))*N_inuc*exp( a_M92              &
                +b_M92*min(ssice(i,j,k),ssice_lim))! o: N_inuc*exp( a_M92 + b_M92*ssice(i,j,k) )
          ! limiting n_in
          n_in = min(n_i_max, n_in)
          ! checking conditions if suitable for nucleation
          if (n_ci(i,j,k).lt.n_in) then ! condition intentionally left this way
            dn_ci_inu(i,j,k) = (n_in-n_ci(i,j,k) )/delt
          else
            dn_ci_inu(i,j,k) = 0.0
            ! note - written this way on purpose
            ! in case of further adjustment of nucleation subroutine
          endif
        endif
        ! basic correction  - not suitable

        ! update water numbers and
        n_cip(i,j,k) = n_cip(i,j,k)+dn_ci_inu(i,j,k)  ! increase in cloud water number

        ! no change n_ccp #inp

        ! update water density [kg kg^{-1}]
        q_cip (i,j,k) = q_cip (i,j,k) + x_inuc*dn_ci_inu(i,j,k)

        ! update liquid water potential temperature
        !  - due to latent heat of melting (freezing in this case)
        qtpmcr(i,j,k) = qtpmcr(i,j,k) - x_inuc*dn_ci_inu(i,j,k)  ! #iceout
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlvi/(cp*exnf(k)))*x_inuc*dn_ci_inu(i,j,k)
      enddo
      enddo
      enddo
    endif !l_sb_reisner
  endif  ! l_sb_inuc_expl

  if (l_sb_dbg) then
    if (any((ssice(i,j,k).gt.0.0).and.((qt0(2:i1,2:j1,1:kmax)-qvsi(2:i1,2:j1,1:kmax)-svm(2:i1,2:j1,1:k1,iq_cl)-svm(2:i1,2:j1,1:k1,iq_ci))/delt - x_inuc*dn_ci_inu(2:i1,2:j1,1:k1)).lt. 0.)) then
      write(6,*) 'WARNING: high ice nucleation'
      write(6,*) ' removing too much water in gridpoints ',count((ssice(i,j,k).gt.0.0).and.   &
        ((qt0(2:i1,2:j1,1:kmax)-qvsi(2:i1,2:j1,1:kmax)-svm(2:i1,2:j1,1:k1,iq_cl)-svm(2:i1,2:j1,1:k1,iq_ci))/delt - x_inuc*dn_ci_inu(2:i1,2:j1,1:k1)).lt. 0.)
    endif
  endif

  ! deallocation
  deallocate ( ssice )

end subroutine icenucle3


!> Sedimentaion of rain
!! sedimentation of drizzle water
!! - gen. gamma distr is assumed. Terminal velocities param according to
!!   Stevens & Seifert. Flux are calc. anal.
!! - l_lognormal =T : lognormal DSD is assumed with D_g and N known and
!!   sig_g assumed. Flux are calc. numerically with help of a
!!   polynomial function
  subroutine sedim_rain3
    use modglobal, only : ih,i1,jh,j1,k1,kmax,eps1,dzf
    use modfields, only : rhof
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop
    real    :: pwcont
    real    :: xr_try, N_r0_try   ! #sb3
    real, allocatable ::  wvar(:,:,:), xr_spl(:,:,:),Dvr_spl(:,:,:)  &
                         ,mur_spl(:,:,:),lbdr_spl(:,:,:),Dgr(:,:,:)  &
                         ,wvar0(:,:,:)
    real, allocatable ::  N_r0(:,:,:), lbdr_try(:,:,:)  ! #sb3
    real,save         ::  dt_spl,wfallmax

    ! allocation

    allocate(  wvar    (2-ih:i1+ih,2-jh:j1+jh,k1)    &!<  work variable
              ,xr_spl  (2-ih:i1+ih,2-jh:j1+jh,k1)    &!<  for time splitting
              ,Dvr_spl (2-ih:i1+ih,2-jh:j1+jh,k1)    &!<     -
              ,mur_spl (2-ih:i1+ih,2-jh:j1+jh,k1)    &!<     -
              ,lbdr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    &!<     -
              ,Dgr     (2-ih:i1+ih,2-jh:j1+jh,k1)    &!<  lognormal geometric diameter
              ,wvar0   (2-ih:i1+ih,2-jh:j1+jh,k1)    &!<  extra test
              ,N_r0     (2-ih:i1+ih,2-jh:j1+jh,k1)   &!< rain integral stuff
              ,lbdr_try (2-ih:i1+ih,2-jh:j1+jh,k1)   &!< rain integral stuff
            )

    ! an zero values
    wvar       = 0.0
    xr_spl     = 0.0
    Dvr_spl    = 0.0
    mur_spl    = 0.0
    lbdr_spl   = 0.0
    Dgr        = 0.0
    wvar0      = 0.0
    N_r0       = 0.0
    lbdr_try   = 0.0

    qr_spl(2:i1,2:j1,1:k1) = q_hr(2:i1,2:j1,1:k1)
    Nr_spl(2:i1,2:j1,1:k1)  = n_hr(2:i1,2:j1,1:k1)


    ! wfallmax = 9.9
    wfallmax = wfallmax_hr
    n_spl = ceiling(split_factor*wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)

    do jn = 1 , n_spl ! time splitting loop

      sed_qr(2:i1,2:j1,1:k1) = 0.
      sed_Nr(2:i1,2:j1,1:k1) = 0.

      if (l_sb ) then

      if (l_sb_classic) then
        do k=1,k1
        do j=2,j1
        do i=2,i1
          if (qr_spl(i,j,k) > qrmin) then
            ! limiting procedure (as per S&B)
            xr_spl (i,j,k) = qr_spl(i,j,k)/(Nr_spl(i,j,k)+eps0)
            xr_try         = max(xrmin,min(xrmax,x_hr(i,j,k)))
            Dvr_spl(i,j,k) = (xr_try/pirhow)**(1./3.)
            N_r0_try       = n_hr(i,j,k)/Dvr(i,j,k) ! rhof(k)*n_hr(i,j,k)/Dvr(i,j,k)
            N_r0(i,j,k)    = max(N_0min/rhof(k),min(N_0max/rhof(k),N_r0_try))

            lbdr_try(i,j,k)= (pirhow*N_r0(i,j,k)/(rhof(k)*qr_spl(i,j,k)))**0.25
            lbdr_spl(i,j,k)    = max(lbdr_min, min(lbdr_max,lbdr_try(i,j,k)))

            ! calculation of velocities
            wfall_qr(i,j,k) = max(0.,((rho0s/rhof(k))**0.5)*(a_tvsbc           &
                      -b_tvsbc*(1.+c_tvsbc/lbdr_spl(i,j,k))**(-4.0))) ! k=1
            wfall_Nr(i,j,k) = max(0.,((rho0s/rhof(k))**0.5)*(a_tvsbc           &
                      -b_tvsbc*(1.+c_tvsbc/lbdr_spl(i,j,k))**(-1.0))) ! k=0
            sed_qr  (i,j,k) = wfall_qr(i,j,k)*qr_spl(i,j,k)*rhof(k)
            sed_Nr  (i,j,k) = wfall_Nr(i,j,k)*Nr_spl(i,j,k)*rhof(k)
          endif
        enddo
        enddo
        enddo
      else  ! l_sb_classic

       if (l_lognormal) then
         do k = 1,kmax
         do j = 2,j1
         do i = 2,i1
           if (qr_spl(i,j,k) > qrmin) then
             ! correction for width of DSD
             Dgr(i,j,k) = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr_spl(i,j,k)

             sed_qr(i,j,k) = 1.*sed_flux3(Nr_spl(i,j,k),Dgr(i,j,k),log(sig_gr)**2,D_s,3)
             sed_Nr(i,j,k) = 1./pirhow*sed_flux3(Nr_spl(i,j,k),Dgr(i,j,k) ,log(sig_gr)**2,D_s,0)

             ! correction for the fact that pwcont .ne. qr_spl
             ! actually in this way for every grid box a fall velocity is determined
             pwcont = liq_cont3(Nr_spl(i,j,k),Dgr(i,j,k),log(sig_gr)**2,D_s,3)         ! note : kg m-3
             if (pwcont > eps1) then
               sed_qr(i,j,k) = (qr_spl(i,j,k)*rhof(k)/pwcont)*sed_qr(i,j,k)
             end if
           end if ! qr_spl
         end do
         end do
         end do

       else
         !
         ! SB rain sedimentation
         !
         if (l_mur_cst) then
           mur_spl(2:i1,2:j1,1:k1) = mur_cst
         else
           do k=1,k1
           do j=2,j1
           do i=2,i1
             if (qr_spl(i,j,k) > qrmin) then
!              mur_spl(i,j,k) = 10. * (1+tanh(1200.*(Dvr_spl(i,j,k)-0.0014))) ! SS08
               mur_spl(i,j,k) = min(30.,- 1. + 0.008/ (qr_spl(i,j,k)*rhof(k))**0.6)  ! G09b
             endif
           enddo
           enddo
           enddo
         endif
         do k=1,k1
         do j=2,j1
         do i=2,i1
           if (qr_spl(i,j,k) > qrmin) then
             lbdr_spl(i,j,k) = ((mur_spl(i,j,k)+3.)*(mur_spl(i,j,k)+2.)* &
                                (mur_spl(i,j,k)+1.))**(1./3.)/Dvr_spl(i,j,k)
             wfall_qr(i,j,k) = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl(i,j,k))**(-1.*(mur_spl(i,j,k)+4.))))
             wfall_Nr(i,j,k) = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl(i,j,k))**(-1.*(mur_spl(i,j,k)+1.))))
             sed_qr  (i,j,k) = wfall_qr(i,j,k)*qr_spl(i,j,k)*rhof(k)
             sed_Nr  (i,j,k) = wfall_Nr(i,j,k)*Nr_spl(i,j,k)*rhof(k)
           endif
         enddo
         enddo
         enddo
      endif ! l_lognormal
     endif ! l_sb_classic
    else
      !
      ! KK00 rain sedimentation
      !
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (qr_spl(i,j,k) > qrmin) then
           !JvdD added eps0 to avoid division by zero
           xr_spl(i,j,k) = rhof(k)*qr_spl(i,j,k)/(Nr_spl(i,j,k)+eps0)

           ! to ensure xr is within borders
           xr_spl(i,j,k) = min(xr_spl(i,j,k),xrmaxkk)

           Dvr_spl(i,j,k) = (xr_spl(i,j,k)/pirhow)**(1./3.)
           sed_qr(i,j,k) = max(0., 0.006*1.0E6*Dvr_spl(i,j,k)- 0.2) * qr_spl(i,j,k)*rhof(k)
           sed_Nr(i,j,k) = max(0.,0.0035*1.0E6*Dvr_spl(i,j,k)- 0.1) * Nr_spl(i,j,k)
        endif
      enddo
      enddo
      enddo

    end if !l_sb

    do k = 1,kmax
      do j=2,j1
      do i=2,i1
        wvar(i,j,k)  = qr_spl(i,j,k) + (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzf(k)*rhof(k))
        wvar0(i,j,k) = Nr_spl(i,j,k) + (sed_Nr(i,j,k+1) - sed_Nr(i,j,k))*dt_spl/(dzf(k)*rhof(k))
      enddo
      enddo
      if (any(wvar(2:i1,2:j1,k) .lt. 0.)) then
        write(6,*)'  rain sedim of q_r too large',count(wvar(2:i1,2:j1,k).lt.0.0),myid, minval(wvar), minloc(wvar)
      end if
      if (any(wvar0(2:i1,2:j1,k) .lt. 0.)) then
        write(6,*)'  rain sedim of N_r too large',count(wvar0(2:i1,2:j1,k).lt.0.0),myid, minval(wvar0), minloc(wvar0)
      end if
      do j=2,j1
      do i=2,i1
        Nr_spl(i,j,k) = Nr_spl(i,j,k) + &
                (sed_Nr(i,j,k+1) - sed_Nr(i,j,k))*dt_spl/(dzf(k)*rhof(k))
        qr_spl(i,j,k) = qr_spl(i,j,k) + &
                (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzf(k)*rhof(k))
      enddo
      enddo

      if ( jn == 1. ) then
      do j=2,j1
      do i=2,i1
        precep_hr(i,j,k) = sed_qr(i,j,k)/rhof(k) ! kg kg-1 m s-1
        precep_l(i,j,k) = precep_l(i,j,k) +precep_hr(i,j,k) ! kg kg-1 m s-1
      enddo
      enddo
      endif

     enddo  ! second k loop
!
    enddo ! time splitting loop

    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! tendencies
      dn_hr_se(i,j,k) = (Nr_spl(i,j,k) - n_hr(i,j,k))/delt
      dq_hr_se(i,j,k) = (qr_spl(i,j,k) - q_hr(i,j,k))/delt

      ! updates
      n_hrp(i,j,k)= n_hrp(i,j,k) + dn_hr_se(i,j,k)
      q_hrp(i,j,k)= q_hrp(i,j,k) + dq_hr_se(i,j,k)
    enddo
    enddo
    enddo

    deallocate (wvar, xr_spl,Dvr_spl,mur_spl,lbdr_spl,Dgr)
    deallocate (wvar0,N_r0, lbdr_try) ! #sb3
  end subroutine sedim_rain3


   !*********************************************************************
   ! sedimentation of snow
   !*********************************************************************
  subroutine sedim_snow3 ! sedim_ice3
    use modglobal, only : ih,i1,jh,j1,k1,kmax,eps1,dzf
    use modfields, only : rhof
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop
    real    :: pwcont, xpmin, xpmax, qip_min                       &
               ,c_v_0, c_v_1, be_ip, aip, bip
    real, allocatable,dimension(:,:,:)  :: qip_spl, nip_spl
    real, allocatable,dimension(:,:,:)  :: sed_qip, sed_nip        &
                                          ,wfall_nip, wfall_qip
    real, allocatable :: wvar(:,:,:),wvar0(:,:,:), xip_spl(:,:,:),Dvp_spl(:,:,:) &
                        ,mur_spl(:,:,:)  ! ,lbdr_spl(:,:,:),Dgr(:,:,:)
    real,save :: dt_spl,wfallmax

    ! --outer part of the code

    ! set constants
    xpmin = x_hs_bmin
    xpmax = x_hs_bmax
    qip_min = qsnowmin
    c_v_0 = c_v_s0
    c_v_1 = c_v_s1
    be_ip = be_hs
    aip   = a_hs
    bip   = b_hs
    wfallmax = wfallmax_hs

    ! allocate

    allocate(  sed_qip(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,sed_nip(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,qip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    &
              ,nip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    &
              ,wfall_qip(2-ih:i1+ih,2-jh:j1+jh,k1)   &
              ,wfall_nip(2-ih:i1+ih,2-jh:j1+jh,k1)   &
            )

    sed_qip = 0.0
    sed_nip = 0.0
    wfall_qip = 0.0
    wfall_nip = 0.0

    qip_spl(2:i1,2:j1,1:k1)  = q_hs(2:i1,2:j1,1:k1)
    nip_spl(2:i1,2:j1,1:k1)  = n_hs(2:i1,2:j1,1:k1)


    ! --- inner part of the code -------------------------

    allocate(  wvar(2-ih:i1+ih,2-jh:j1+jh,k1)       & !<  work variable
              ,wvar0(2-ih:i1+ih,2-jh:j1+jh,k1)       &
              ,xip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<  for time splitting
              ,Dvp_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<
             )
            !  ,mur_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<     -
            !  ,Dgr(2-ih:i1+ih,2-jh:j1+jh,k1)        & !<  lognormal geometric diameter
            ! )

      wvar    = 0.0
      wvar0   = 0.0
      xip_spl = 0.0
      Dvp_spl = 0.0



    ! inner part of the code
    n_spl = ceiling(split_factor*wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)


    do jn = 1 , n_spl ! time splitting loop
      sed_qip(2:i1,2:j1,1:k1) = 0.
      sed_nip(2:i1,2:j1,1:k1) = 0.

      do k=1,k1
      do j=2,j1
      do i=2,i1
        if ((qip_spl(i,j,k) > qip_min).and.(nip_spl(i,j,k) > 0.0)) then
          xip_spl (i,j,k) = qip_spl(i,j,k)/(nip_spl(i,j,k)+eps0) ! JvdD Added eps0 to avoid division by zero
          xip_spl (i,j,k) = min(max(xip_spl(i,j,k),xpmin),xpmax) ! to ensure xr is within borders
          ! Dvp_spl(i,j,k) = aip*xip_spl(i,j,k)**bip
        endif
      enddo
      enddo
      enddo

      ! terminal fall velocity
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if ( (qip_spl(i,j,k) > qip_min).and.(nip_spl(i,j,k) > 0.0) ) then
          wfall_qip(i,j,k) = max(0.0,c_v_1 * xip_spl (i,j,k)**be_ip)
          wfall_nip(i,j,k) = max(0.0,c_v_0 * xip_spl (i,j,k)**be_ip)
          sed_qip(i,j,k)   = wfall_qip(i,j,k)*qip_spl(i,j,k)*rhof(k)
          sed_nip(i,j,k)   = wfall_nip(i,j,k)*nip_spl(i,j,k)*rhof(k)
        endif
      enddo
      enddo
      enddo

      ! segmentation over levels
      do k = 1,kmax
        do j=2,j1
        do i=2,i1
           wvar(i,j,k) = qip_spl(i,j,k) + (sed_qip(i,j,k+1) - sed_qip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
           wvar0(i,j,k) = nip_spl(i,j,k) + (sed_nip(i,j,k+1) - sed_nip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
        enddo
        enddo
        if (any(wvar(2:i1,2:j1,k) .lt. 0.)) then
          write(6,*)'  snow sedim too large', count(wvar(2:i1,2:j1,k) .lt. 0.),myid, minval(wvar), minloc(wvar)
        end if
        if (any(wvar0(2:i1,2:j1,k) .lt. 0.)) then
          write(6,*)'  snow sedim too large', count(wvar0(2:i1,2:j1,k) .lt. 0.),myid, minval(wvar0), minloc(wvar0)
        end if
        do j=2,j1
        do i=2,i1
          nip_spl(i,j,k) = nip_spl(i,j,k) + &
                  (sed_nip(i,j,k+1) - sed_nip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
          qip_spl(i,j,k) = qip_spl(i,j,k) + &
                  (sed_qip(i,j,k+1) - sed_qip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
        enddo
        enddo

        if ( jn == 1. ) then
          do j=2,j1
          do i=2,i1
            precep_hs(i,j,k) = sed_qip(i,j,k)/rhof(k) ! kg kg-1 m s-1
            precep_i(i,j,k) = precep_i(i,j,k)+ precep_hs(i,j,k)  ! kg kg-1 m s-1
          enddo
          enddo
        endif
      enddo  ! second k loop
    enddo ! time splitting loop


    ! --- end of the inner part of the code

    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! tendencies
      dn_hs_se(i,j,k) = (nip_spl(i,j,k) - n_hs(i,j,k))/delt
      dq_hs_se(i,j,k) = (qip_spl(i,j,k) - q_hs(i,j,k))/delt

      ! updates
      n_hsp(i,j,k)= n_hsp(i,j,k) + dn_hs_se(i,j,k)
      q_hsp(i,j,k)= q_hsp(i,j,k) + dq_hs_se(i,j,k)
    enddo
    enddo
    enddo

    deallocate (nip_spl, qip_spl,wfall_nip, wfall_qip, sed_nip, sed_qip)
    deallocate (wvar, wvar0, xip_spl,Dvp_spl) ! mur_spl,lbdr_spl,Dgr)
  end subroutine sedim_snow3


!*********************************************************************
! sedimentation of graupel
!*********************************************************************
  subroutine sedim_graupel3 ! sedim_ice3
    use modglobal, only : ih,i1,jh,j1,k1,kmax,eps1,dzf
    use modfields, only : rhof
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop
    real    :: pwcont, xpmin, xpmax, c_v_0, c_v_1, be_ip, aip, bip
    real    :: qip_min
    real, allocatable,dimension(:,:,:)  :: qip_spl, nip_spl
    real, allocatable,dimension(:,:,:)  :: sed_qip, sed_nip        &
                                          ,wfall_nip, wfall_qip
    real, allocatable :: wvar(:,:,:),wvar0(:,:,:), xip_spl(:,:,:),Dvp_spl(:,:,:) &
                        ,mur_spl(:,:,:)  ! ,lbdr_spl(:,:,:),Dgr(:,:,:)
    real :: dt_spl,wfallmax ! real,save :: dt_spl,wfallmax


    ! --outer part of the code

    allocate(  sed_qip(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,sed_nip(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,qip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,nip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,wfall_qip(2-ih:i1+ih,2-jh:j1+jh,k1)   &
              ,wfall_nip(2-ih:i1+ih,2-jh:j1+jh,k1)   &
            )

    sed_qip = 0.0
    sed_nip = 0.0
    wfall_qip = 0.0
    wfall_nip = 0.0

    ! write(6,*)'starting sedimentation'

    qip_spl(2:i1,2:j1,1:k1)  = q_hg(2:i1,2:j1,1:k1)
    nip_spl(2:i1,2:j1,1:k1)  = n_hg(2:i1,2:j1,1:k1)

    ! set constants
    xpmin = x_hg_bmin
    xpmax = x_hg_bmax
    c_v_0 = c_v_g0
    c_v_1 = c_v_g1
    be_ip = be_hs ! be_hg
    aip   = a_hg
    bip   = b_hg
    qip_min= qgrmin
    wfallmax = wfallmax_hg

    ! --- inner part of the code -------------------------

    allocate( wvar(2-ih:i1+ih,2-jh:j1+jh,k1)        & !<  work variable
              ,wvar0(2-ih:i1+ih,2-jh:j1+jh,k1)       &
              ,xip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<  for time splitting
              ,Dvp_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<
             )
            !  ,mur_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<     -
            !  ,Dgr(2-ih:i1+ih,2-jh:j1+jh,k1)        & !<  lognormal geometric diameter
            ! )

    wvar    = 0.0
    wvar0   = 0.0
    xip_spl = 0.0
    Dvp_spl = 0.0

    ! inner part of the code
    n_spl = ceiling(split_factor*wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)


    do jn = 1 , n_spl ! time splitting loop
      sed_qip(2:i1,2:j1,1:k1) = 0.
      sed_nip(2:i1,2:j1,1:k1) = 0.

      do k=1,k1
      do j=2,j1
      do i=2,i1
       if ( (qip_spl(i,j,k) > qip_min).and.(nip_spl(i,j,k) > 0.0) ) then
         xip_spl (i,j,k) = qip_spl(i,j,k)/(nip_spl(i,j,k)+eps0) ! JvdD Added eps0 to avoid division by zero
         xip_spl (i,j,k) = min(max(xip_spl(i,j,k),xpmin),xpmax) ! to ensure xr is within borders
         ! Dvp_spl(i,j,k) = aip*xip_spl(i,j,k)**bip
       endif
      enddo
      enddo
      enddo

      ! terminal fall velocity
      do k=1,k1
      do j=2,j1
      do i=2,i1
         if ( (qip_spl(i,j,k) > qip_min).and.(nip_spl(i,j,k) > 0.0) ) then
            wfall_qip(i,j,k) = max(0.0,c_v_1 * xip_spl (i,j,k)**be_ip)
            wfall_nip(i,j,k) = max(0.0,c_v_0 * xip_spl (i,j,k)**be_ip)
            sed_qip(i,j,k)   = wfall_qip(i,j,k)*qip_spl(i,j,k)*rhof(k)
            sed_nip(i,j,k)   = wfall_nip(i,j,k)*nip_spl(i,j,k)*rhof(k)
         endif
      enddo
      enddo
      enddo

      ! segmentation over levels
      do k = 1,kmax
        do j=2,j1
        do i=2,i1
          wvar(i,j,k) = qip_spl(i,j,k) +  (sed_qip(i,j,k+1) - sed_qip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
          wvar0(i,j,k) = nip_spl(i,j,k) + (sed_nip(i,j,k+1) - sed_nip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
        enddo
        enddo
        if (any(wvar(2:i1,2:j1,k) .lt. 0.)) then
          write(6,*)'  graupel sedim too large', count(wvar(2:i1,2:j1,k) .lt. 0.),myid, minval(wvar), minloc(wvar)
        end if
        if (any(wvar0(2:i1,2:j1,k) .lt. 0.)) then
          write(6,*)'  graupel sedim too large', count(wvar0(2:i1,2:j1,k) .lt. 0.),myid, minval(wvar0), minloc(wvar0)
        end if
        do j=2,j1
        do i=2,i1
          nip_spl(i,j,k) = nip_spl(i,j,k) + &
                  (sed_nip(i,j,k+1) - sed_nip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
          qip_spl(i,j,k) = qip_spl(i,j,k) + &
                  (sed_qip(i,j,k+1) - sed_qip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
        enddo
        enddo

        ! -> check this part properly later
        if ( jn == 1. ) then
          do j=2,j1
          do i=2,i1
            precep_hg(i,j,k) = precep_hg(i,j,k)+ sed_qip(i,j,k)/rhof(k)   ! kg kg-1 m s-1
            precep_i(i,j,k) =  precep_i(i,j,k)+ precep_hg(i,j,k)  ! kg kg-1 m s-1
          enddo
          enddo
        endif
      enddo  ! second k loop
    enddo ! time splitting loop

    ! --- end of the inner part of the code

    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! tendencies
      dn_hg_se(i,j,k) = (nip_spl(i,j,k) - n_hg(i,j,k))/delt
      dq_hg_se(i,j,k) = (qip_spl(i,j,k) - q_hg(i,j,k))/delt

      ! updates
      n_hgp(i,j,k)= n_hgp(i,j,k) + dn_hg_se(i,j,k)
      q_hgp(i,j,k)= q_hgp(i,j,k) + dq_hg_se(i,j,k)
    enddo
    enddo
    enddo

    deallocate (nip_spl, qip_spl,wfall_nip, wfall_qip, sed_nip, sed_qip)
    deallocate (wvar, wvar0, xip_spl,Dvp_spl) ! mur_spl,lbdr_spl,Dgr)
  end subroutine sedim_graupel3

!*********************************************************************
! sedimentation of cloud ice
!*********************************************************************
  subroutine sedim_ice3 ! sedim_ice3
    use modglobal, only : ih,i1,jh,j1,k1,kmax,eps1,dzf,rlv,cp
    use modfields, only : rhof, exnf
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop
    real    :: pwcont, xpmin, xpmax, c_v_0, c_v_1, be_ip, aip, bip
    real    :: qip_min
    real, allocatable,dimension(:,:,:)  :: qip_spl, nip_spl
    real, allocatable,dimension(:,:,:)  :: sed_qip, sed_nip        &
                                          ,wfall_nip, wfall_qip
    real, allocatable :: wvar(:,:,:), xip_spl(:,:,:),Dvp_spl(:,:,:) &
                        ,mur_spl(:,:,:)  ! ,lbdr_spl(:,:,:),Dgr(:,:,:)
    real,save :: dt_spl,wfallmax


    ! --outer part of the code

    allocate(  sed_qip(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,sed_nip(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,qip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,nip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,wfall_qip(2-ih:i1+ih,2-jh:j1+jh,k1)   &
              ,wfall_nip(2-ih:i1+ih,2-jh:j1+jh,k1)   &
            )

    sed_qip = 0.0
    sed_nip = 0.0
    wfall_qip = 0.0
    wfall_nip = 0.0

    ! write(6,*)'starting sedimentation'

    qip_spl(2:i1,2:j1,1:k1)  = q_ci(2:i1,2:j1,1:k1)
    nip_spl(2:i1,2:j1,1:k1)  = n_ci(2:i1,2:j1,1:k1)

    ! set constants
    xpmin = x_ci_bmin
    xpmax = x_ci_bmax
    c_v_0 = c_v_i0
    c_v_1 = c_v_i1
    be_ip = be_ci
    aip   = a_ci
    bip   = b_ci
    qip_min= qicemin
    wfallmax = wfallmax_ci

    ! --- inner part of the code -------------------------

    allocate( wvar(2-ih:i1+ih,2-jh:j1+jh,k1)        & !<  work variable
              ,xip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<  for time splitting
              ,Dvp_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<
             )
            !  ,mur_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<     -
            !  ,Dgr(2-ih:i1+ih,2-jh:j1+jh,k1)        & !<  lognormal geometric diameter
            ! )

    wvar    = 0.0
    xip_spl = 0.0
    Dvp_spl = 0.0

    ! inner part of the code
    n_spl = ceiling(split_factor*wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)

    do jn = 1 , n_spl ! time splitting loop

      sed_qip(2:i1,2:j1,1:k1) = 0.
      sed_nip(2:i1,2:j1,1:k1) = 0.

      do k=1,k1
      do j=2,j1
      do i=2,i1
       if ( (qip_spl(i,j,k) > qip_min).and.(nip_spl(i,j,k) > 0.0) ) then
         xip_spl (i,j,k) = qip_spl(i,j,k)/(nip_spl(i,j,k)+eps0) ! JvdD Added eps0 to avoid division by zero
         xip_spl (i,j,k) = min(max(xip_spl(i,j,k),xpmin),xpmax) ! to ensure xr is within borders
       endif
      enddo
      enddo
      enddo

      ! terminal fall velocity
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if ( (qip_spl(i,j,k) > qip_min).and.(nip_spl(i,j,k) > 0.0) ) then
          wfall_qip(i,j,k) = max(0.0,c_v_1 * xip_spl (i,j,k)**be_ip)
          wfall_nip(i,j,k) = max(0.0,c_v_0 * xip_spl (i,j,k)**be_ip)
          sed_qip(i,j,k)   = wfall_qip(i,j,k)*qip_spl(i,j,k)*rhof(k)
          sed_nip(i,j,k)   = wfall_nip(i,j,k)*nip_spl(i,j,k)*rhof(k)
        endif
      enddo
      enddo
      enddo

      ! segmentation over levels
      do k = 1,kmax
        do j=2,j1
        do i=2,i1
          wvar(i,j,k) = qip_spl(i,j,k) + (sed_qip(i,j,k+1) - sed_qip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
        enddo
        enddo
        if (any(wvar(2:i1,2:j1,k) .lt. 0.)) then
          write(6,*)'ice sedim too large', count(wvar(2:i1,2:j1,k) .lt. 0.),myid, minval(wvar), minloc(wvar)
        end if
        do j=2,j1
        do i=2,i1
          nip_spl(i,j,k) = nip_spl(i,j,k) + &
                  (sed_nip(i,j,k+1) - sed_nip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
          qip_spl(i,j,k) = qip_spl(i,j,k) + &
                  (sed_qip(i,j,k+1) - sed_qip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
        enddo
        enddo

        if ( jn == 1. ) then
          do j=2,j1
          do i=2,i1
            precep_ci(i,j,k) = sed_qip(i,j,k)/rhof(k)   ! kg kg-1 m s-1
            precep_i(i,j,k) = precep_i(i,j,k)+ precep_ci(i,j,k)   ! kg kg-1 m s-1
          enddo
          enddo
        endif

    enddo ! time splitting loop

    ! --- end of the inner part of the code
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! tendencies
      dn_ci_se(i,j,k) = (nip_spl(i,j,k) - n_ci(i,j,k))/delt
      dq_ci_se(i,j,k) = (qip_spl(i,j,k) - q_ci(i,j,k))/delt

      ! updates
      n_cip(i,j,k)= n_cip(i,j,k) + dn_ci_se(i,j,k)
      q_cip(i,j,k)= q_cip(i,j,k) + dq_ci_se(i,j,k)

      ! also qtpmcr and thlpmcr change
      qtpmcr(i,j,k)  = qtpmcr(i,j,k) + 0.0
      thlpmcr(i,j,k) = thlpmcr(i,j,k) +0.0
    enddo
    enddo
    enddo

    deallocate (nip_spl, qip_spl,wfall_nip, wfall_qip, sed_nip, sed_qip)
    deallocate (wvar, xip_spl,Dvp_spl) ! mur_spl,lbdr_spl,Dgr)
  end subroutine sedim_ice3

!*********************************************************************
! sedimentation of cloud water
!*********************************************************************
  subroutine sedim_cl3 ! sedim_ice3
    use modglobal, only : ih,i1,jh,j1,k1,kmax,eps1,dzf,rlv,cp
    use modfields, only : rhof, exnf
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop
    real    :: pwcont, xpmin, xpmax, c_v_0, c_v_1, be_ip, aip, bip
    real    :: qip_min
    real, allocatable,dimension(:,:,:)  :: qip_spl, nip_spl
    real, allocatable,dimension(:,:,:)  :: sed_qip, sed_nip        &
                                          ,wfall_nip, wfall_qip
    real, allocatable :: wvar(:,:,:), xip_spl(:,:,:),Dvp_spl(:,:,:) &
                        ,mur_spl(:,:,:)  ! ,lbdr_spl(:,:,:),Dgr(:,:,:)
    real,save :: dt_spl,wfallmax


    ! --outer part of the code

    allocate(  sed_qip(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,sed_nip(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,qip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,nip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,wfall_qip(2-ih:i1+ih,2-jh:j1+jh,k1)   &
              ,wfall_nip(2-ih:i1+ih,2-jh:j1+jh,k1)   &
            )

    sed_qip = 0.0
    sed_nip = 0.0
    wfall_qip = 0.0
    wfall_nip = 0.0

    qip_spl(2:i1,2:j1,1:k1)  = q_cl(2:i1,2:j1,1:k1)
    nip_spl(2:i1,2:j1,1:k1)  = n_cl(2:i1,2:j1,1:k1)

    ! set constants
    xpmin = x_cl_bmin
    xpmax = x_cl_bmax
    c_v_0 = c_v_c0
    c_v_1 = c_v_c1
    be_ip = be_cl
    aip   = a_cl
    bip   = b_cl
    qip_min= qcliqmin
    wfallmax = wfallmax_cl

    ! --- inner part of the code -------------------------

    allocate( wvar(2-ih:i1+ih,2-jh:j1+jh,k1)        & !<  work variable
              ,xip_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<  for time splitting
              ,Dvp_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<
             )
            !  ,mur_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & !<     -
            !  ,Dgr(2-ih:i1+ih,2-jh:j1+jh,k1)        & !<  lognormal geometric diameter
            ! )

    wvar    = 0.0
    xip_spl = 0.0
    Dvp_spl = 0.0


    ! inner part of the code
    n_spl = ceiling(split_factor*wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)


    do jn = 1 , n_spl ! time splitting loop

      sed_qip(2:i1,2:j1,1:k1) = 0.
      sed_nip(2:i1,2:j1,1:k1) = 0.

      do k=1,k1
      do j=2,j1
      do i=2,i1
        if ((qip_spl(i,j,k) > qip_min).and.(nip_spl(i,j,k) > 0.0)) then
          xip_spl (i,j,k) = qip_spl(i,j,k)/(nip_spl(i,j,k)+eps0) ! JvdD Added eps0 to avoid division by zero
          xip_spl (i,j,k) = min(max(xip_spl(i,j,k),xpmin),xpmax) ! to ensure xr is within borders
        endif
      enddo
      enddo
      enddo

      ! terminal fall velocity
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if ( (qip_spl(i,j,k) > qip_min).and.(nip_spl(i,j,k) > 0.0) ) then
           wfall_qip(i,j,k) = max(0.0,c_v_1 * xip_spl (i,j,k)**be_ip)
           wfall_nip(i,j,k) = max(0.0,c_v_0 * xip_spl (i,j,k)**be_ip)
           sed_qip(i,j,k)   = wfall_qip(i,j,k)*qip_spl(i,j,k)*rhof(k)
           sed_nip(i,j,k)   = wfall_nip(i,j,k)*nip_spl(i,j,k)*rhof(k)
        endif
      enddo
      enddo
      enddo

      ! segmentation over levels
      do k = 1,kmax
        do j=2,j1
        do i=2,i1
          wvar(i,j,k) = qip_spl(i,j,k) + (sed_qip(i,j,k+1) - sed_qip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
        enddo
        enddo
        if (any(wvar(2:i1,2:j1,k) .lt. 0.)) then
          ! t
          write(6,*)'cloud sedim too large', count(wvar(2:i1,2:j1,k) .lt. 0.),myid, minval(wvar), minloc(wvar)
        end if
        do j=2,j1
        do i=2,i1
          nip_spl(i,j,k) = nip_spl(i,j,k) + &
                  (sed_nip(i,j,k+1) - sed_nip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
          qip_spl(i,j,k) = qip_spl(i,j,k) + &
                  (sed_qip(i,j,k+1) - sed_qip(i,j,k))*dt_spl/(dzf(k)*rhof(k))
        enddo
        enddo
      enddo  ! second k loop

    enddo ! time splitting loop


    ! --- end of the inner part of the code
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! tendencies
      dn_cl_se(i,j,k) = (nip_spl(i,j,k) - n_cl(i,j,k))/delt
      dq_cl_se(i,j,k) = (qip_spl(i,j,k) - q_cl(i,j,k))/delt

      ! updates
      n_clp(i,j,k)= n_clp(i,j,k) + dn_cl_se(i,j,k)
      q_clp(i,j,k)= q_clp(i,j,k) + dq_cl_se(i,j,k)

      ! also qtpmcr and thlpmcr change
      qtpmcr(i,j,k)  = qtpmcr(i,j,k) +  dq_cl_se(i,j,k)
      thlpmcr(i,j,k) = thlpmcr(i,j,k)-(rlv/(cp*exnf(k)))*dq_cl_se(i,j,k)
    enddo
    enddo
    enddo

    deallocate (nip_spl, qip_spl,wfall_nip, wfall_qip, sed_nip, sed_qip)
    deallocate (wvar, xip_spl,Dvp_spl) ! mur_spl,lbdr_spl,Dgr)
  end subroutine sedim_cl3


  subroutine evap_rain3
  !*********************************************************************
  ! Evaporation of prec. : Seifert & Beheng
  ! Cond. (S>0.) neglected (all water is condensed on cloud droplets)
  !*********************************************************************

    use modglobal, only : ih,i1,jh,j1,k1,rv,rlv,cp,pi,mygamma251,mygamma21,lacz_gamma
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,qvsl,rhof,exnf
    implicit none
    integer :: i,j,k
    real, allocatable, dimension(:,:,:) :: f0,f1,S,G,vihr,nrex, x_hrf
    integer :: numel

    allocate(  f0(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! ventilation factor - moment 0
              ,f1(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! ventilation factor - moment 1
              ,S(2-ih:i1+ih,2-jh:j1+jh,k1)     & ! super or undersaturation
              ,G(2-ih:i1+ih,2-jh:j1+jh,k1)     & ! cond/evap rate of a drop
              ,vihr(2-ih:i1+ih,2-jh:j1+jh,k1)  & ! mean terminal velocity
              ,nrex(2-ih:i1+ih,2-jh:j1+jh,k1)  & ! Reynolds number N_re(xr)
              ,x_hrf(2-ih:i1+ih,2-jh:j1+jh,k1) & ! full x_hr without bounds
             )

    S = 0.0
    f0 = 0.0
    f1 = 0.0
    G = 0.0
    vihr = 0.0
    nrex = 0.0
    x_hrf= 0.0

    do k=1,k1
    do j=2,j1
    do i=2,i1
      if (q_hr_mask(i,j,k)) then
        ! adjusting the calculation for saturation
        S (i,j,k) = min(0.,((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)- 1.0))
        G   (i,j,k) = (rv * tmp0(i,j,k)) / (Dv*esl(i,j,k)) + rlv/(Kt*tmp0(i,j,k))*(rlv/(rv*tmp0(i,j,k)) -1.)
        G   (i,j,k) = 1./G(i,j,k)

        ! terminal velocity  (from mixed scheme)
        vihr(i,j,k) = al_hr*((rho0s/rhof(k))**0.5)*x_hr(i,j,k)**be_hr  ! xr in this calculation !!

        ! calculating  N_re Reynolds number
        nrex(i,j,k) = Dvr(i,j,k)*vihr(i,j,k)/nu_a
        x_hrf(i,j,k) = q_hr(i,j,k)/(n_hr(i,j,k)+eps0)

        f0(i,j,k) = aven_0r+bven_0r*Sc_num**(1.0/3.0)*nrex(i,j,k)**0.5
        f1(i,j,k) = aven_1r+bven_1r*Sc_num**(1.0/3.0)*nrex(i,j,k)**0.5
        f0(i,j,k) = max(0.0,f0(i,j,k))

        dq_hr_ev(i,j,k) = 2*pi*n_hr(i,j,k)*G(i,j,k)*  &
             Dvr(i,j,k)*f1(i,j,k)*S(i,j,k)

        dn_hr_ev(i,j,k) = 2*pi*n_hr(i,j,k)*G(i,j,k)* &
             Dvr(i,j,k)*f0(i,j,k)*S(i,j,k)/x_hrf(i,j,k) ! x_hr here, not xr

        ! and limiting it
        dn_hr_ev(i,j,k) = min(dn_hr_ev(i,j,k), 0.0)
        dn_hr_ev(i,j,k) = max(dn_hr_ev(i,j,k),dq_hr_ev(i,j,k)/x_hrf(i,j,k))

        if (((dq_hr_ev(i,j,k) +svm(i,j,k,iq_hr)/delt).lt.0).or.((dn_hr_ev(i,j,k) +svm(i,j,k,in_hr)/delt).lt.0)) then
          dn_hr_ev(i,j,k) = - svm(i,j,k,in_hr)/delt
          dq_hr_ev(i,j,k) = - svm(i,j,k,iq_hr)/delt
        endif
      endif
    enddo
    enddo
    enddo

    do k=1,k1
    do j=2,j1
    do i=2,i1
      q_hrp(i,j,k) = q_hrp(i,j,k) + dq_hr_ev(i,j,k)
      n_hrp(i,j,k) = n_hrp(i,j,k) + dn_hr_ev(i,j,k)
      qtpmcr(i,j,k)  = qtpmcr(i,j,k) -dq_hr_ev(i,j,k)
      thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*dq_hr_ev(i,j,k)

      ! recovery of aerosols ?
      ret_cc(i,j,k) = ret_cc(i,j,k) - c_ccn_ev_r*min(0.0,dn_hr_ev(i,j,k))
    enddo
    enddo
    enddo

    deallocate (f0, f1 ,S,G, vihr, nrex,x_hrf)

  end subroutine evap_rain3

! ****************************************
! Depositional growth of cloud ice particles
!  later tunr into:
!
! Depositional growth of various ice particles
! Inner subroutine - takes input from a wrapper
!  following S&B
! ****************************************
subroutine deposit3
  use modglobal, only : i1,j1,k1,rv,rd,cp,pi
  use modfields, only : exnf,qt0,svm,tmp0,,qvsi,rhof,presf
  implicit none

  real :: esi

  real    :: tocon,precon,cond_cf
  real    :: cor_dqci_dep,cor_dqhs_dep,cor_dqhg_dep

  real :: F_ci    & ! ventilation factor
         ,q_avail & ! available water for deposition
         ,Si      & ! super or undersaturation with respect to ice
         ,G       & ! G_iv function
         ,Dvic_ci & ! size of ice parti
         ,viic_ci & ! v for ice cloud particles
         ,nrex_ci & ! reynolds number
         ,xip     & ! particle size
         ,nip     & ! particle number

  integer :: i,j,k

  do k=1,k1
  do j=2,j1
  do i=2,i1
    ! TODO: any of the masks q_ci_mask, q_hg_mask, q_hs_mask
    if (q_ci_mask(i,j,k).and.(tmp0(i,j,k).le.T_3)) then
      xip = x_ci(2:i1,2:j1,1:k1)
      nip = n_ci(2:i1,2:j1,1:k1)

      q_avail = qt0(i,j,k)-q_cl(i,j,k)- qvsi(i,j,k)
      Si = q_avail/qvsi(i,j,k)

      ! calculating G_iv
      esi = qvsi(i,j,k)*presf(k)/(rd/rv+(1.0-rd/rv)*qvsi(i,j,k))
      G = (rv * tmp0(i,j,k)) / (Dv*esi) + rlvi/(Kt*tmp0(i,j,k))*(rlvi/(rv*tmp0(i,j,k)) -1.)
      G = 1./G

      ! diameter
      Dvic_ci = a_ci*xip**b_ci
      Dvic_hs = a_hs*xip**b_hs
      Dvic_hg = a_hg*xip**b_hg

      ! terminal velocity
      viic_ci = al_ci*((rho0s/rhof(k))**0.5)*xip**be_ci
      viic_hs = al_hs*((rho0s/rhof(k))**0.5)*xip**be_hs
      viic_hg = al_hg*((rho0s/rhof(k))**0.5)*xip**be_hg

      ! N_re Reynolds number
      nrex_ci = Dvic_ci*viic_ci/nu_a
      nrex_hs = Dvic_hs*viic_hs/nu_a
      nrex_hg = Dvic_hg*viic_hg/nu_a

      ! calculating from prepared ventilation coefficients
      F_ci = avent_1i+bvent_1i*Sc_num**(1.0/3.0)*nrex_ci**0.5
      F_hs = aven_1s+bven_1s*Sc_num**(1.0/3.0)*nrex_hs**0.5
      F_hg = aven_1g+bven_1g*Sc_num**(1.0/3.0)*nrex_hg**0.5

      ! depositional growth
      ! k_depos = 4*pi/c  for spherical particles
      ! k_depos = 4*pi/c  for hexagonal plates - included
      dq_ci_dep = (4*pi/c_ci)*nip*G*Dvic_ci*F_ci*Si
      dq_hs_dep = (4*pi/c_hs)*nip*G*Dvic_hs*F_hs*Si
      dq_hg_dep = (4*pi/c_hg)*nip*G*Dvic_hg*F_hg*Si

      ! and limiting not to deposit more than available
      ! ie. allows any negative dq_ci_dep but positive only smaller than amount of available water
      dq_ci_dep = max(dq_ci_dep,min(0.0,-svm(i,j,k,iq_ci)/delt-q_cip(i,j,k)))
      dq_ci_dep = min(dq_ci_dep,max(0.0,q_avail/delt))

      dq_hs_dep = max(dq_hs_dep,min(0.0,-svm(i,j,k,iq_hs)/delt-q_hsp(i,j,k)))
      dq_hs_dep = min(dq_hs_dep,max(0.0,q_avail/delt))

      dq_hg_dep = max(dq_hg_dep,min(0.0,-svm(i,j,k,iq_hg)/delt-q_hgp))
      dq_hg_dep = min(dq_hg_dep,max(0.0,q_avail/delt))

      ! adds to outputs
      q_cip = q_cip + dq_ci_dep
      q_hsp = q_hsp + dq_hs_dep
      q_hgp = q_hgp + dq_hg_dep

      qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_ci_dep - dq_hs_dep - dq_hg_dep
      thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlvi/(cp*exnf(k)))*(dq_ci_dep+dq_hs_dep+dq_hg_dep)

      if (l_sb_dbg) then
        if((svm(i,j,k,iq_ci)+delt*dq_ci_dep).lt. 0.0) then
          write(6,*) 'WARNING: ice_deposit3 too low'
        endif
        if(((qt0(i,j,k)-qvsi(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci)-delt*dq_ci_dep).lt.0.0) then
          write(6,*) 'WARNING: ice_deposit3 too high depositing more ice than available',
          ! too high:    ((qt0(i,j,k)-qvsi(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci)-delt*dq_ci_dep).lt. 0.0)
          ! negative qt: ((qt0(i,j,k)-delt*dq_ci_dep).lt. 0.0)
        endif

        if((svm(i,j,k,iq_hs)+delt*dq_hs_dep).lt. 0.0) then
          write(6,*) 'WARNING: deposit_snow3 too low'
          write(6,*) '  sublimating more snow than available'
        endif
        if((((qt0(i,j,k)-qvsi(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci)-delt*dq_hs_dep)).lt. 0.0)) then
          write(6,*) 'WARNING: deposit_snow3 too high'
          write(6,*) '  depositing more water than available'
          ! count((qt0(i,j,k)-qvsi(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci)-delt*dq_hs_dep).lt. 0.0)
          ! getting negative q_t in count((qt0(i,j,k)-delt*dq_hs_dep).lt. 0.0)
        endif

        if((svm(i,j,k,iq_hg)+delt*dq_hg_dep).lt.0.) then
          write(6,*) 'WARNING: deposit_graupel3 too low'
          write(6,*) '  sublimating more snow than available'
          ! count((svm(i,j,k,iq_hg)+delt*dq_hg_dep).lt. 0.0)
        endif
        if(((qt0(i,j,k)-qvsi(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci)-delt*dq_hg_dep).lt. 0.0)) then
          write(6,*) 'WARNING: deposit_graupel3 too high'
          write(6,*) '  depositing more water than available'
          ! count(qt0(i,j,k)-qvsi(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci)-delt*dq_hg_dep .lt. 0.0)
          ! getting negative q_t in  count((qt0(i,j,k)-delt*dq_hg_dep).lt. 0.0)
        endif
      endif
    endif

    ! available water vapour for deposition
    tocon = (qt0(i,j,k)-svm(i,j,k,iq_cl)-qvsi(i,j,k))/delt

    ! consumption of water vapour calculated by nucleation and deposition processes
    precon = dq_ci_dep+dq_hs_dep+dq_hg_dep

    if ((precon.gt.0.0).and.(tocon-precon).lt.0.0) then ! run only in oversaturted conditions
      ! preparing additive correctors:
      cond_cf = (tocon/(precon+eps0) - 1.0)  ! later added eps0 to prevent 0 values
      cond_cf = max(min(0.0, cond_cf),-1.0)

      ! - corrector for deposition - only if positive deposition
      cor_dqci_dep = cond_cf*max(0.0, dq_ci_dep)
      cor_dqhs_dep = cond_cf*max(0.0, dq_hs_dep)
      cor_dqhg_dep = cond_cf*max(0.0, dq_hg_dep)

      ! and updating values:
      ! - cloud water content
      q_cip = q_cip+cor_dqci_dep

      ! - correction for hydrometeors
      q_hsp = q_hsp+cor_dqhs_dep
      q_hgp = q_hgp+cor_dqhg_dep

      ! - correction for total water
      qtpmcr(i,j,k) = qtpmcr(i,j,k)-cor_dqhs_dep-cor_dqhg_dep-cor_dqci_dep

      ! - and correcting for heat
      thlpmcr(i,j,k) = thlpmcr(i,j,k)                     &
            +(rlvi/(cp*exnf(k)))*cor_dqci_dep             &
            +(rlvi/(cp*exnf(k)))*cor_dqhs_dep             &
            +(rlvi/(cp*exnf(k)))*cor_dqhg_dep

      ! corrector for process values
      dq_ci_dep =dq_ci_dep+cor_dqci_dep
      dq_hs_dep =dq_hs_dep+cor_dqhs_dep
      dq_hg_dep =dq_hg_dep+cor_dqhg_dep
    endif
  enddo
  enddo
  enddo
end subroutine deposit3

! *******************************************************
! Heterogenous and Homogeneous freezing of cloud droplets
! *******************************************************
subroutine hethomfreez3
  use modglobal, only : i1,j1,k1,rlv,cp
  use modfields, only : svm,tmp0,exnf
  implicit none

  real    :: J_hom   ! freezing rate
  real    :: tmpj    ! adjusted temperature

  real    :: dq_cl_hom !< homogeneous freezing of cloud water
  real    :: dn_cl_hom !< homogeneous freezing of cloud water

  integer :: i,j,k

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_cl_mask(i,j,k).and.tmp0(i,j,k).lt.T_3) then

      ! homogeneous freezing of cloud droplets
      ! --------------------------------------

      ! inserting temperature values
      !  -- can be later adjusted to include the effect of chemicals in clouds
      tmpj = tmp0(i,j,k)

      if (tmpj>tmp_lim1_Cf02) then
          J_hom = exp(C_Cf02+B_Cf02*(tmpj+CC_Cf02))
      else
        if(tmpj<tmp_lim2_Cf02) then
          J_hom = exp(C_30_Cf02)
        else
          J_hom = exp(C_20_Cf02                 &
             + B_21_Cf02*(tmpj-offset_Cf02)     &
             + B_22_Cf02*(tmpj-offset_Cf02)**2  &
             + B_23_Cf02*(tmpj-offset_Cf02)**3  &
             + B_24_Cf02*(tmpj-offset_Cf02)**4  )
        endif ! quick freezing
      endif ! slow freezing

      dn_cl_hom        = -c_mmt_1cl *n_cl(i,j,k) * x_cl(i,j,k)* J_hom
      dq_cl_hom        = -c_mmt_2cl *q_cl(i,j,k) * x_cl(i,j,k)* J_hom
      dn_cl_hom        = max(dn_cl_hom,min(0.0,-svm(i,j,k,in_cl)/delt-n_clp(i,j,k)))
      dq_cl_hom        = max(dq_cl_hom,min(0.0,-svm(i,j,k,iq_cl)/delt-q_clp(i,j,k)))

      ! heterogeneous freezing of cloud droplets
      ! ---------------------------------------
      J_het = A_het *exp( B_het*(T_3-tmp0(i,j,k)) -1)

      dn_cl_het = -c_mmt_1cl *n_cl(i,j,k) * x_cl(i,j,k)* J_het
      dq_cl_het = -c_mmt_2cl *q_cl(i,j,k) * x_cl(i,j,k)* J_het

      ! basic correction
      dn_cl_het = max(dn_cl_het,min(0.0,-svm(i,j,k,in_cl)/delt-n_clp(i,j,k)))
      dq_cl_het = max(dq_cl_het,min(0.0,-svm(i,j,k,iq_cl)/delt-q_clp(i,j,k)))

      ! Apply changes
      ! -------------

      ! changes in cloud water
      n_clp(i,j,k) = n_clp(i,j,k) + dn_cl_hom + dn_cl_het
      q_clp(i,j,k) = q_clp(i,j,k) + dq_cl_hom + dq_cl_het


      ! increase in cloud ice
      n_cip(i,j,k) = n_cip(i,j,k) - dn_cl_hom - dn_cl_het
      q_cip(i,j,k) = q_cip(i,j,k) - dq_cl_hom - dq_cl_het

      ! changes in the total amount of water
      qtpmcr(i,j,k) = qtpmcr(i,j,k) + dq_cl_hom + dq_cl_het

      ! change in th_l due to freezing
      thlpmcr(i,j,k) = thlpmcr(i,j,k)-((rlv+rlme)/(cp*exnf(k)))*(dq_cl_hom+dq_cl_het)
    endif ! cl_mask
  enddo
  enddo
  enddo

end subroutine hethomfreez3

! ****************************************
! Heterogeneou freezing of rain
!
!
! ****************************************
    subroutine rainhetfreez3

    use modglobal, only : ih,i1,jh,j1,k1,rv,rd, rlv,cp,pi,lacz_gamma
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,qvsl,qvsi,rhof,exnf,presf
    implicit none
    integer :: i,j,k
    real, allocatable, dimension(:,:,:) :: J_het
    logical ,allocatable :: qfr_mask(:,:,:)

    ! --- inner part -------------------

    ! allocate fields and fill
     allocate( qfr_mask   (2-ih:i1+ih,2-jh:j1+jh,k1))
     allocate( J_het      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! freezing rate
             )
    ! - filling
    J_het  = 0.0
    dn_hr_het = 0.0
    dq_hr_het = 0.0

    ! putting together masks
    qfr_mask(2:i1,2:j1,1:k1) = q_hr_mask(2:i1,2:j1,1:k1)

    ! performing calculation for the freezing rate
    do j=2,j1
    do i=2,i1
    do k=1,k1
      ! calculating for all cells with the value
      if (qfr_mask(i,j,k).and.(tmp0(i,j,k).lt.T_3)) then
          ! maybe only for temperatures below T_3 ?
          J_het(i,j,k) = A_het *exp( B_het*(T_3-tmp0(i,j,k)) -1)
      endif
    enddo
    enddo
    enddo


    ! freezing of rain
    do k=1,k1
    do j=2,j1
    do i=2,i1
      if (q_hr_mask(i,j,k).and.(tmp0(i,j,k).lt.T_3)) then
        dn_hr_het(i,j,k) = -c_mmt_1hr *n_hr(i,j,k) * x_hr(i,j,k)* J_het(i,j,k)
        dq_hr_het(i,j,k) = -c_mmt_2hr *q_hr(i,j,k) * x_hr(i,j,k)* J_het(i,j,k)

        ! basic correction
        dn_hr_het(i,j,k) = max(dn_hr_het(i,j,k),min(0.0,-svm(i,j,k,in_hr)/delt-n_hrp(i,j,k)))
        dq_hr_het(i,j,k) = max(dq_hr_het(i,j,k),min(0.0,-svm(i,j,k,iq_hr)/delt-q_hrp(i,j,k)))

        ! decrease in raindrops
        n_hrp(i,j,k) = n_hrp(i,j,k) + dn_hr_het(i,j,k)
        q_hrp(i,j,k) = q_hrp(i,j,k) + dq_hr_het(i,j,k)

        ! increase in graupel
        n_hgp(i,j,k) = n_hgp(i,j,k) - dn_hr_het(i,j,k)
        q_hgp(i,j,k) = q_hgp(i,j,k) - dq_hr_het(i,j,k)

        ! and consumption of aerosols for heterogeneous freezing  ?
        ! n_ccp(i,j,k) = n_ccp(i,j,k) - dn_hr_het(i,j,k)

        ! no changes in the total amount of water
        ! qtpmcr(i,j,k) = qtpmcr(i,j,k)

        ! change in th_l due to freezing
        thlpmcr(i,j,k) = thlpmcr(i,j,k) - (rlme/(cp*exnf(k)))*dq_hr_het(i,j,k)
      endif
    enddo
    enddo
    enddo

    deallocate(J_het)
    deallocate(qfr_mask)

   end subroutine rainhetfreez3


! Collision/collection processes
! ****************************************
! ice selfcollection and aggregation to snow
!  - all collection of ice by ice treated as snow
!
! snow selfcollection
! snow collecting cloud ice
!
! follows Seifert, 2002
! ****************************************
subroutine ice_aggr_snow_self
  use modglobal, only : i1,j1,k1,rv,pi
  use modfields, only : qt0,svm,tmp0,rhof
  implicit none
  integer :: i,j,k

  real    :: dlt_0aa_i ,dlt_0aa_s          &
            ,dlt_1aa_i ,dlt_1aa_s          &
            ,th_0aa_i  ,th_0aa_s           &
            ,th_1aa_i  ,th_1aa_s           &
            ,rem_cf_i  ,rem_cf_s

  real :: dif_D_10, x_crit_ii, x_minagg_ii

  ! Temporary variables, value changes several times in this routine
  real :: E_ab, E_stick

  ! calculate constants
  dlt_0aa_i   = 2*dlt_i0 + dlt_i0i
  dlt_1aa_i   = dlt_i0 + dlt_i1i + dlt_i1
  th_0aa_i    = 2*th_i0 - th_i0i   ! from Seifert, 2002
  th_1aa_i    = th_i0 - th_i1i + th_i1

  dlt_0aa_s   = 2*dlt_s0 + dlt_s0s
  dlt_1aa_s   = dlt_s0 + dlt_s1s + dlt_s1
  th_0aa_s    = 2*th_s0 - th_s0s   ! from Seifert, 2002
  th_1aa_s    = th_s0 - th_s1s + th_s1 ! th_1ab    = th_i1i

  ! adjusting coefficient
  ! prepare coefficient for remaining water number
  rem_cf_i = (1.0-rem_n_ci_min)/delt
  rem_cf_s = (1.0-rem_n_hs_min)/delt

  ! and the minimal conversion size
  x_minagg_ii = (D_i_b/a_ci)**(1.0/b_ci)

  ! and the critical size for the start of conversion
  x_crit_ii   = (D_crit_ii/a_ci)**(1.0/b_ci)


  do k=1,k1
  do j=2,j1
  do i=2,i1

    ! Process ice
    ! -----------
    dif_D_10  = D_i_b-D_i_a  ! difference in ice cloud droplet intervals

    if(q_ci_mask.and.(x_ci.gt.x_crit_ii).and.(q_ci.gt.q_crit_ii)) then
      ! calculating sticking efficiency
      if (l_sb_stickyice) then
        E_stick = c_E_o_s*exp(B_stick *(tmp0(i,j,k)+stick_off))
        E_stick = min(c_E_o_s,E_stick)
      else ! l_sb_stickyice
        E_stick = exp(B_stick_ii*(tmp0(i,j,k)+stick_off)+C_stick_ii)
        E_stick = min(E_ii_maxst,E_stick)
      endif ! l_sb_stickyice

      ! collision efficiency
      if (l_sb_lim_aggr) then
        ! checking whether sufficient size
        if( D__ci.gt.D_i_a ) then
          if( D_ci.gt.D_i_b ) then
            E_ab = E_ee_m*E_stick
          else
            E_ab = (E_ee_m*E_stick /dif_D_10)* (D_ci - D_i_a)
          endif
        endif
      else ! l_sb_lim_aggr
        E_ab = E_ee_m*E_stick
      endif

      dn_ci_col_iis = -rhof(k)*(pi/4)*E_ab                      &
                      *n_ci(i,j,k)**2 *dlt_0aa_i*D_ci**2        &
                      *(th_0aa_i *v_ci**2+2*sigma_ci**2)**0.5

      dq_ci_col_iis = -rhof(k)*(pi/4)*E_ab                           &
                      *dlt_1aa_i*n_ci(i,j,k)*q_ci(i,j,k)*D_ci**2     &
                      *(th_1aa_i*v_ci**2+2*sigma_ci**2)**0.5


      ! limiting dq_ci  --   limit x_cv_ii as per ICON 2017
      dq_ci_col_iis = max(dq_ci_col_iis, (-svm(i,j,k,iq_ci)/delt-q_cip(i,j,k)))
      dn_ci_col_iis = max(dn_ci_col_iis, (-rem_cf_i*svm(i,j,k,in_ci)-n_cip(i,j,k)))
      dn_ci_col_iis = max(dn_ci_col_iis, dq_ci_col_iis/x_minagg_ii)
    else
      dq_ci_col_iis = 0.
      dn_ci_col_iis = 0.
    endif

    ! Process snow
    ! ------------
    if(q_hs_mask) then
      ! calculating sticking efficiency
      if (l_sb_stickyice) then
        E_stick = c_E_o_s*exp(B_stick *(tmp0(i,j,k)+stick_off))
        E_stick = min(c_E_o_s,E_stick)
        E_ab_s = E_ee_m*E_stick
      else ! l_sb_stickyice
        E_stick = exp(B_stick_ii*(tmp0(i,j,k)+stick_off)+C_stick_ii)
        E_stick = min(E_ss_maxst,E_stick)
        E_ab_s  = E_ee_m*E_stick
      endif ! l_sb_stickyice

      dn_hs_col_sss = -rhof(k)*(pi/4)*E_ab_s                  &
                      *n_hs(i,j,k)**2 *dlt_0aa_s*D_hs**2      &
                      *(th_0aa_s *v_hs**2+2*sigma_hs**2)**0.5

      dn_hs_col_sss = max(min(0.0,dn_hs_col_sss),-rem_cf_s*svm(i,j,k,in_hs)-n_hsp(i,j,k))
    else
      dn_hs_col_sss = 0.
    endif

    ! snow collecting cloud ice
    ! -------------------------
    if (q_hs_mask(i,j,k).and.q_ci_mask(i,j,k)) then
      ! calculating sticking efficiency
      E_stick = c_E_o_s*exp(B_stick *(tmp0(i,j,k)+stick_off))
      E_stick =min(c_E_o_s,E_stick)
      E_ab = E_ee_m*E_stick

      dq_hsci_col = (rhof(k)*pi/4)*E_ab*n_hs(i,j,k)                &
                    *q_ci(i,j,k)*(dlt_s0*D_hs**2                   &
                    +dlt_s1i*D_hs*D_ci+dlt_i1*D_ci**2)             &
                    *(th_s0*v_hs**2-th_s1i*v_ci*v_hs               &
                    +th_1i*v_ci**2+sigma_hs**2+sigma_ci**2)**0.5

      dn_ci_col_hs = -(rhof(k)*pi/4)*E_ab*n_hs(i,j,k)              &
                     *n_ci(i,j,k)*(dlt_s0*D_hs**2                  &
                     +dlt_s0i*D_hs*D_ci+dlt_i0*D_ci**2)            &
                     *(th_s0*v_hs**2-th_s0i*v_ci*v_hs              &
                     +th_i0*v_ci**2+sigma_hs**2+sigma_ci**2)**0.5

      dq_hsci_col = min(dq_hsci_col,max(0.0,svm(i,j,k,iq_ci)/delt+q_cip(i,j,k)))   ! following ICON, 2017
      dn_ci_col_hs = max(dn_ci_col_hs,min(0.0,-rem_cf_i*svm(i,j,k,in_ci)-n_cip(i,j,k)))
    else
      dq_hsci_col = 0.
      dn_ci_col_hs = 0.
    endif

    ! graupel collecting snow
    if (q_hg_mask.and.q_hs_mask) then
      ! calculating sticking efficiency
      E_stick = c_E_o_s*exp(B_stick *(tmp0(i,j,k)+stick_off))
      E_stick =min(c_E_o_s,E_stick)
      E_ab = E_ee_m*E_stick

      dq_hghs_col  = (rhof(k)*pi/4)*E_ab*n_hg(i,j,k)       &
           *q_hs(i,j,k)*(dlt_g0*D_hg**2                    &
             +dlt_g1s*D_hg*D_hs+dlt_s1*D_hs**2)            &
           *( th_g0*v_hg**2-th_g1s*v_hs*v_hg               &
             +th_s1*v_hs**2+sigma_hg**2+sigma_hs**2)**0.5

      dn_hs_col_hg = -(rhof(k)*pi/4)*E_ab*n_hg(i,j,k)      &
           *n_hs(i,j,k)*(dlt_g0*D_hg**2                    &
             +dlt_g0s*D_hg*D_hs+dlt_s0*D_hs**2)            &
           *( th_g0*v_hg**2-th_g0s*v_hs*v_hg               &
             +th_s0*v_hs**2+sigma_hg**2+sigma_hs**2)**0.5


      dq_hghs_col = min(dq_hghs_col,max(0.0,svm(i,j,k,iq_hs)/delt+q_hsp(i,j,k)))   ! following ICON, 2017
      dn_hs_col_hg = max(dn_hs_col_hg,min(0.0,-rem_cf_s*svm(i,j,k,in_hs)-n_hsp(i,j,k)))
    else
      dq_hghs_col = 0.
      dn_hs_col_hg = 0.
    endif

    n_cip    = n_cip + dn_ci_col_iis + dn_ci_col_hs
    q_cip    = q_cip + dq_ci_col_iis - dq_hsci_col
    n_hsp    = n_hsp - 0.5 * dn_ci_col_iis + dn_hs_col_sss + dn_hs_col_hg
    q_hsp    = q_hsp - dq_ci_col_iis + dq_hsci_col - dq_hghs_col
    q_hgp    = q_hgp + dq_hghs_col

    if (l_sb_dbg) then
      if((svm(i,j,k,iq_ci)+delt*dq_ci_col_iis .lt. 0.0)) then
        write(6,*) 'WARNING: ice_aggr3 too high'
        write(6,*) ' removing more ice than available'
        ! count((svm(i,j,k,iq_ci) +delt*dq_ci_col_iis).lt. 0.0)
        write(6,*) ' removing too much ice'
        ! count(( q_ci(i,j,k)+delt*dq_ci_col_iis).lt. 0.0 )
        write(6,*) ' getting negative q_t'
        ! count(( qt0(i,j,k)+delt*q_cip(i,j,k)).lt. 0.0 )
      endif

      if((svm(i,j,k,in_ci)+delt*dn_ci_col_iis .lt. 0.0)) then
        write(6,*) 'WARNING: ice_aggr3 too high'
        write(6,*) ' removing more ice particles then available'
        ! count((svm(i,j,k,in_ci) +delt*dn_ci_col_iis).lt. 0.0)
        write(6,*) ' removing too much ice particles'
        ! count((n_ci(i,j,k)+delt*dn_ci_col_iis).lt. 0.0)
      endif

      if(( svm(i,j,k,in_hs)+delt*dn_hs_col_sss .lt. 0.0 )) then
        write(6,*) 'WARNING: snow self-collection too high'
        write(6,*) ' decreasing number of snowflakes below 0'
        ! count((svm(i,j,k,in_hs)+delt*dn_hs_col_sss).lt. 0.0 )
        write(6,*) ' decreasing number of snowflakes too much'
        ! count((n_hs(i,j,k)+delt*dn_hs_col_sss).lt. 0.0 )
      endif

      if(svm(i,j,k,iq_ci)-delt*dq_hsci_col.lt. 0.0) then
        write(6,*) 'WARNING: coll_sis3 too high removing more ice than available'
        write(6,*) ' removing more ice than available'
        ! count((svm(i,j,k,iq_ci)-delt*dq_hsci_col).lt. 0.0)
        write(6,*) ' removing too much ice'
        ! count((q_ci(i,j,k)-delt*dq_hsci_col).lt. 0.0)
        write(6,*) ' getting negative q_t'
        ! count((qt0(i,j,k)-delt*q_cip).lt. 0.0 )
      endif

      if(svm(i,j,k,in_ci)+delt*dn_col_b.lt. 0.0) then
        write(6,*) 'WARNING: coll_sis3 too high'
        write(6,*) ' removing more ice particles then available in gridpoints'
        ! count(svm(i,j,k,in_ci)+delt*dn_ci_col_hs.lt. 0.0)
        write(6,*) ' removing too many ice particles in gridpoints '
        ! count((n_ci(i,j,k)+delt*dn_ci_col_hs).lt. 0.0)
      endif

      if((svm(i,j,k,iq_hs)-delt*dq_hghs_col).lt. 0.0) then
        write(6,*) 'WARNING: coll_gsg3 too high removing more ice than available'
        write(6,*) ' removing more ice than available'
        ! count(svm(i,j,k,iq_hs)-delt*dq_hghs_col.lt. 0.0)
        write(6,*) ' removing too much ice'
        ! count(q_hs-delt*dq_hghs_col.lt. 0.0)
        write(6,*) ' getting negative q_t'
        ! count(qt0(i,j,k)-delt*q_hsp.lt. 0.0)
      endif

      if(svm(i,j,k,in_hs)+delt*dn_hs_col_hg.lt. 0.0) then
        write(6,*) 'WARNING: coll_gsg3 too high'
        write(6,*) ' removing more ice particles then available'
        ! count(svm(i,j,k,in_hs)+delt*dn_hs_col_hg.lt. 0.0)
        write(6,*) ' removing too many ice particles'
        ! count(n_hs+delt*dn_hs_col_hg.lt. 0.0)
      endif

    endif
  enddo
  enddo
  enddo

end subroutine ice_aggr_snow_self


subroutine coll3
  use modglobal, only : i1,j1,k1,cp,pi
  use modfields, only : qt0,svm,tmp0,rhof,exnf

  implicit none
  integer :: i,j,k

  ! Temporary variables, value changes several times in this routine
  real :: E_ab
  real :: dif_D_10
  real :: k_enhm
  real :: rem_cf
  real :: dn_col_ici, dq_col_ici
  real :: dn_col_b, dq_col_a

  do k=1,k1
  do j=2,j1
  do i=2,i1

    ! Collection of cloud droplets by ice
    ! - based on Seifert & Beheng (2004)
    ! - resulting process is:
    !    - ice riming by cloud ice
    !    - enhanced melting
    ! ------------------------------------

    ! collision efficiency
    if (q_ci_mask.and.q_cl_mask) then
      ! checking whether sufficient size
      if( (D_cl.gt.D_c_a).and.(D_ci.gt.D_i0)) then
        if( D_cl.gt.D_c_b ) then
          E_ab = E_i_m
        else
          ! denominator in calculationg collision efficiency
          dif_D_10  = D_c_b-D_c_a

          E_ab = (E_i_m /dif_D_10)* (D_cl - D_c_a)
        endif
      endif
      dq_col_ici = (rhof(k)*pi/4)*E_ab*n_ci       &
           *q_cl*(dlt_i0*D_ci**2                &
             +dlt_i1c*D_ci*D_cl+dlt_c1*D_cl**2) &
           *( th_i0*v_ci**2-th_i1c*v_cl*v_ci    &
             +th_c1*v_cl**2+sigma_ci**2+sigma_cl**2)**0.5

      dn_col_ici = -(rhof(k)*pi/4)*E_ab*n_ci      &
           *n_cl(i,j,k)*(dlt_i0*D_ci**2         &
             +dlt_i0c*D_ci*D_cl+dlt_c0*D_cl**2) &
           *( th_i0*v_ci**2-th_i0c*v_cl*v_ci    &
             +th_c0*v_cl**2+sigma_ci**2+sigma_cl**2)**0.5

      if(tmp0(i,j,k).lt.T_3) then
        ! riming
        ! adjusting coefficient
        ! prepare coefficient for remaining water number
        rem_cf = (1.0-rem_n_cl_min)/delt

        dq_ci_rime = min(dq_col_ici,max(0.0,svm(i,j,k,iq_cl)/delt+q_clp))   ! following ICON, 2017
        dn_cl_rime_ci =max(dn_col_ici,min(0.0,-rem_cf*svm(i,j,k,in_cl)-n_clp))
      else
        ! scheding and enhanced melting
        k_enhm  = c_water/rlme

        ! calculating the melting
        dq_ci_eme_ic = -k_enhm*(tmp0(i,j,k)-T_3)*min(dq_col_ici,q_cl/delt)
        dq_ci_eme_ic = max(dq_ci_eme_ic,min(0.0,-svm(i,j,k,iq_ci)/delt-q_cip))

        ! calculating number of melted particles
        ! - expected to be proportional to melted mass
        dn_ci_eme_ic = dq_ci_eme_ic*n_ci/(q_ci+eps0)

        ! - but not more than number of interacting particles
        dn_ci_eme_ic = max(dn_ci_eme_ic, max(dn_col_ici,-n_cl(i,j,k)/delt))

        ! - and not more than total number of particles
        dn_ci_eme_ic = max(dn_ci_eme_ic,-svm(i,j,k,in_ci)/delt-n_cip)
      endif
    endif

    n_clp    = n_clp + dn_cl_rime_ci
    q_clp    = q_clp - dq_ci_rime - dq_ci_eme_ic
    n_cip    = n_cip + dn_ci_eme_ic
    q_cip    = qcip + dq_ci_rime + dq_ci_eme_ic

    qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_ci_rime - dq_ci_eme_ic
    thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlvi/(cp*exnf(k)))*(dq_ci_rime + dq_ci_eme_ic)

    if (l_sb_dbg) then
      if(svm(i,j,k,iq_cl)-delt*dq_col_ici.lt. 0.0) then
        write(6,*) 'WARNING: coll_ici3 too high'
        write(6,*) ' removing more cloud water then available'
        ! count(svm(i,j,k,iq_cl)-delt*dq_ci_rime.lt. 0.0)
        write(6,*) ' removing too much water'
        ! count(q_cl-delt*dq_ci_rime.lt. 0.0)
        write(6,*) ' getting negative q_t'
        ! count(qt0(i,j,k)+delt*q_cip.lt. 0.0)
      endif
      if(svm(i,j,k,in_cl)+delt*dn_col_ici.lt. 0.0) then
        write(6,*) 'WARNING: coll_ici3 too high'
        write(6,*) ' removing more droplets then available'
        ! count(svm(i,j,k,in_cl)+delt*dn_cl_rime_ci.lt. 0.0)
        write(6,*) ' removing too many droplets'
        ! count(n_cl+delt*dn_cl_rime_ci.lt. 0.0)
      endif
    endif

    ! riming of snow
    ! ---------------
    if (q_hs_mask.and.q_cl_mask) then

      ! checking whether sufficient size
      if( (D_cl.gt.D_c_a).and.(D_hs.gt.D_i0)) then
        if( D_cl.gt.D_c_b ) then
          E_ab = E_s_m
        else
          ! denominator in calculationg collision efficiency
          dif_D_10  = D_c_b-D_c_a
          E_ab = (E_s_m/dif_D_10)* (D_cl - D_c_a)
        endif
      endif

      dq_col_a = (rhof(k)*pi/4)*E_ab*n_hs(i,j,k)                  &
           *q_cl*(dlt_s0*D_hs**2                                  &
             +dlt_s1c*D_hs*D_cl+dlt_c1*D_cl**2)                   &
           *( th_s0*v_hs**2-th_s1c*v_cl(i,j,k)*v_hs(i,j,k)        &
             +th_c1*v_cl(i,j,k)**2+sigma_hs**2+sigma_cl**2)**0.5

      dn_col_b = -(rhof(k)*pi/4)*E_ab*n_hs(i,j,k)                 &
           *n_cl*(dlt_s0*D_hs**2                                  &
             +dlt_s0c*D_hs*D_cl+dlt_c0*D_cl**2)                   &
           *( th_s0*v_hs(i,j,k)**2-th_s0c*v_cl(i,j,k)*v_hs(i,j,k) &
             +th_c0*v_cl(i,j,k)**2+sigma_hs**2+sigma_cl**2)**0.5

      ! basic limiting
      ! limited by amount of water in droplets
      dq_col_a = min(dq_col_a,q_cl/delt)
      dq_col_a = min(dq_col_a,max(0.0,svm(i,j,k,iq_cl)/delt+q_clp))
      dn_col_b = max(dn_col_b,min(0.0,-svm(i,j,k,in_cl)/delt-n_clp))
      ! limited by number of droplets
      ! then based on temperature

      if(tmp0(i,jk).lt.T_3) then
        ! riming only
        dq_hs_rime = dq_col_a ! following ICON, 2017

        ! = min(dq_hs_rime,svm(i,j,k,iq_cl)/delt)
        dn_cl_rime_hs = dn_col_b ! min(dn_cl_rime_hs,-svm(i,j,k,in_cl)/delt)

        ! record the change
        ! o dq_hs_rime = dq_col_a
        ! change in the amount of cloud ice
        q_hsp = q_hsp + dq_hs_rime

        ! change in the amount of cloud water
        n_clp = n_clp + dn_cl_rime_hs
        q_clp = q_clp - dq_hs_rime

        ! change in q_t
        qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_hs_rime

        ! change in th_l - freezing and removal
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+ (rlvi/(cp*exnf(k)))*dq_hs_rime
      elseif(tmp0(i,jk).gt.T_3) then ! BUG: should be plain else
        ! not riming, but enhanced melting and scheding

        ! enhanced melting
        k_enhm  = c_water/rlme

        ! calculating the melting
        dq_hs_eme_sc = -k_enhm*(tmp0(i,j,k)-T_3)*dq_col_a
        dq_hs_eme_sc = max(dq_hs_eme_sc, &
                           min(0.0, &
                               -svm(i,j,k,iq_hs)/delt-q_hsp))

        ! calculating number of melted particles
        ! - expected to be proportional to melted mass
        dn_hs_eme_sc = dq_hs_eme_sc*n_hs(i,j,k)/(q_hs(i,j,k)+eps0) ! q_hs here is always some small positive number

        ! - not more than number of interacting particles
        ! dn_hs_eme_sc = max(dn_hs_eme_sc, max(dn_col_b,-n_cl/delt))
        ! - and not more than total number of particles
        dn_hs_eme_sc = max(dn_hs_eme_sc, &
                           -min((svm(i,j,k,in_cl)/delt+n_clp), &
                                (svm(i,j,k,in_hs)/delt+n_hsp(i,j,k))))

        ! updating tendencies
        ! updating rain
        ! melted snow is turning into rain  RH84
        q_hrp = q_hrp -  dq_hs_eme_sc + dq_col_a ! both mass of snow and droplets
        n_hrp = n_hrp -  dn_hs_eme_sc

        ! snow
        q_hsp = q_hsp + dq_hs_eme_sc
        n_hsp = n_hsp + dn_hs_eme_sc

        ! updating cloud water
        q_clp = q_clp - dq_col_a
        n_clp = n_clp + dn_col_b

        ! updating thermodynamic
        ! qtp : increased by melted water
        qtpmcr(i,j,k) = qtpmcr(i,j,k) -dq_col_a ! -dq_hs_eme_sc

        ! thlp : melting and adding liquid water
        ! thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlvi/(cp*exnf(k)))*dq_hs_eme_sc
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+            &
            (rlme/(cp*exnf(k)))*dq_hs_eme_sc        &
            +(rlvi/(cp*exnf(k)))*dq_col_a
      endif

      if (l_sb_dbg) then
        if((svm(i,j,k,iq_cl)-delt*dq_col_a).lt. 0.0) then
          write(6,*) 'WARNING: coll_scs too high'
          write(6,*) ' removing more cloud water than available'
          ! count((svm(i,j,k,iq_cl)-delt*dq_hs_rime).lt. 0.0)
          write(6,*) ' removing too much cloud water',
          ! count(q_cl-delt*dq_hs_rime.lt. 0.0)
          write(6,*) ' getting negative q_t'
          ! count(qt0(i,j,k)+delt*q_clp).lt. 0.0)
        endif

        if(svm(i,j,k,in_cl)+delt*dn_col_b.lt. 0.0) then
          write(6,*) 'WARNING: coll_scs too high'
          write(6,*) ' removing more droplets then available'
          ! count(svm(i,j,k,in_cl)+delt*dn_cl_rime_hs).lt. 0.0)
          write(6,*) ' removing too many droplets'
          ! count(n_cl+delt*dn_cl_rime_hs.lt. 0.0)
        endif
      endif
    endif

    ! riming of graupel by cloud droplets
    ! -----------------------------------
    if (q_hg_mask.and.q_cl_mask) then
      ! collision efficiency
      ! checking whether sufficient size
      if( (D_cl.gt.D_c_a).and.(D_hg.gt.D_i0)) then
        if( D_cl.gt.D_c_b ) then
          E_ab = E_g_m
        else
          dif_D_10 = D_c_b-D_c_a
          E_ab = (E_g_m/dif_D_10)* (D_cl- D_c_a)
        endif
      endif

      dq_col_a = (rhof(k)*pi/4)*E_ab*n_hg    &
           *q_cl*(dlt_g0*D_hg**2                     &
             +dlt_g1c*D_hg*D_cl+dlt_c1*D_cl**2) &
           *( th_g0*v_hg**2-th_g1c*v_cl*v_hg    &
             +th_c1*v_cl**2+sigma_hg**2+sigma_cl**2)**0.5

      dn_col_b = -(rhof(k)*pi/4)*E_ab*n_hg   &
           *n_cl*(dlt_g0*D_hg**2                     &
             +dlt_g0c*D_hg*D_cl+dlt_c0*D_cl**2) &
           *( th_g0*v_hg**2-th_g0c*v_cl*v_hg    &
             +th_c0*v_cl**2+sigma_hg**2+sigma_cl**2)**0.5


      ! initial correction based on amount of cloud water
      rem_cf = (1.0-rem_n_cl_min)/delt
      dq_col_a = min(dq_col_a,max(0.0,svm(i,j,k,iq_cl)/delt+q_clp)) ! following ICON, 2017
      dn_col_b = max(dn_col_b,min(0.0,-rem_cf*svm(i,j,k,in_cl)-n_clp))

      ! then based on temeprature
      if(tmp0(i,j,k).lt.T_3) then
        dq_hg_rime = dq_col_a    ! = min(dq_hg_rime,svm(i,j,k,iq_cl)/delt)
        dn_cl_rime_hg = dn_col_b ! = min(dn_cl_rime_hg,-svm(i,j,k,in_cl)/delt)

        ! change in the amount of graupel
        q_hgp = q_hgp + dq_hg_rime

        ! change in the amount of rain
        n_clp = n_clp + dn_cl_rime_hg
        q_clp = q_clp - dq_hg_rime

        ! no change in q_t
        qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_hg_rime

        ! change in th_l - heat release from freezing and removal
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+ (rlvi/(cp*exnf(k)))*dq_hg_rime
      else
        ! not riming,but enhanced melting and scheding

        ! enhanced melting
        k_enhm  = c_water/rlme

        ! calculating the melting
        dq_hg_eme_gc = -k_enhm*(tmp0(i,j,k)-T_3)*min(dq_col_a,q_cl/delt)
        dq_hg_eme_gc = max(dq_hg_eme_gc,min(0.0,-svm(i,j,k,iq_hg)/delt-q_hgp))

        ! calculating number of melted particles
        ! - expected to be proportional to melted mass
        dn_hg_eme_gc = dq_hg_eme_gc*n_hg/(q_hg+eps0) ! q_hg here is always some small positive number

        ! - but not more than number of interacting particles
        ! dn_hg_eme_gc = max(dn_hg_eme_gc, max(dn_col_b,-n_cl/delt))
        ! - and not more than total number of particles
        dn_hg_eme_gc = max(dn_hg_eme_gc, -min((svm(i,j,k,in_cl)/delt+n_clp),(svm(i,j,k,in_hg)/delt+n_hgp)))

        ! updating tendencies
        ! updating rain
        ! based on RH84
        q_hrp = q_hrp - dq_hg_eme_gc + dq_col_a
        n_hrp = n_hrp - dn_col_b

        ! cloud ice
        q_hgp = q_hgp + dq_hg_eme_gc
        n_hgp = n_hgp + dn_hg_eme_gc

        ! updating cloud water
        q_clp = q_clp - dq_col_a !   dq_hg_eme_gc
        n_clp = n_clp + dn_col_b ! + dn_cl_rime_hg

        ! updating thermodynamic
        ! qtp : removal of droplets
        qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_col_a ! dq_hg_eme_gc

        ! thlp : melting and adding liquid water
        ! thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlvi/(cp*exnf(k)))*dq_hg_eme_gc
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+                     &
          (rlme/(cp*exnf(k)))*dq_hg_eme_gc+(rlvi/(cp*exnf(k)))*dq_col_a
      endif

      if (l_sb_dbg) then
        if(svm(i,j,k,iq_cl)-delt*dq_hg_rime.lt. 0.0) then
          write(6,*) 'WARNING: coll_gcg3 too high'
          write(6,*) ' removing more cloud water then available'
          ! count(svm(i,j,k,iq_cl)-delt*dq_hg_rime).lt. 0.0)
          write(6,*) ' removing too much water'
          ! count(q_cl-delt*dq_hg_rime.lt. 0.0)
          write(6,*) ' getting negative q_t'
          ! count(qt0(i,j,k)+delt*q_clp).lt. 0.0)
        endif

        if(svm(i,j,k,in_cl)+delt*dn_cl_rime_hg.lt. 0.0) then
          write(6,*) 'WARNING: coll_gcg3 too high'
          write(6,*) ' removing more droplets then available'
          ! count(svm(i,j,k,in_cl)+delt*dn_cl_rime_hg.lt. 0.0)
          write(6,*) ' removing too many droplets'
          ! count(n_cl+delt*dn_cl_rime_hg.lt. 0.0)
        endif
      endif
    endif

endsubroutine coll3

 ! ***********************************
 !  riming of ice + rain to graupel
 !
 !
 ! ***********************************
   subroutine coll_rig3

    use modglobal, only : ih,i1,jh,j1,k1,rv,rd, rlv,cp,pi
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,qvsl,qvsi,rhof,exnf,presf
    implicit none
    integer :: i,j,k
    real    :: sigma_a, sigma_b ,E_coli              &
               ,a_a, b_a, al_a, be_a, ga_a           &
               ,a_b, b_b, al_b, be_b, ga_b           &
               ,dlt_0a, dlt_0ab, dlt_0b              &
               ,dlt_1ab, dlt_1b                      &
               ,th_0a, th_0ab, th_0b                 &
               ,th_1ab, th_1b                        &
               ,dlt_0ba, dlt_1ba, dlt_1a             &
               ,th_0ba, th_1ba, th_1a
    real    :: dif_D_10, ntest, qtest, rem_ci_cf, rem_hr_cf,k_enhm
    real, allocatable, dimension(:,:,:) :: D_a, D_b, v_a, v_b ! <- later move outside of this subroutine
    real, allocatable, dimension(:,:,:) :: E_ab, E_stick, dn_col_b, dq_col_a, dq_col_b
    logical ,allocatable :: qcol_mask(:,:,:)

    ! start of the code

    ! allocate fields and fill
     allocate( D_a     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,D_b     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,v_a     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,v_b     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
             )

      allocate( E_ab    (2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,E_stick (2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,dn_col_b(2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,dq_col_a(2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,dq_col_b(2-ih:i1+ih,2-jh:j1+jh,k1)    &
             )
      allocate( qcol_mask(2-ih:i1+ih,2-jh:j1+jh,k1)   )

      D_a      = 0.0
      D_b      = 0.0
      v_a      = 0.0
      v_b      = 0.0
      E_ab     = 0.0
      E_stick  = 0.0
      dn_col_b = 0.0
      dq_col_a = 0.0
      dq_col_b = 0.0

    ! set constants
      sigma_a   = sigma_ci
      sigma_b   = sigma_hr
      a_a       = a_ci
      b_a       = b_ci
      al_a      = al_ci
      be_a      = be_ci
      ga_a      = ga_ci
      a_b       = a_hr
      b_b       = b_hr
      al_b      = al_hr
      be_b      = be_hr
      ga_b      = ga_hr
      E_coli    = E_i_m ! collision efficienty
      dlt_0a    = dlt_i0 ! dlt_s0
      dlt_0ab   = dlt_i0r ! dlt_s0i
      dlt_0b    = dlt_r0
      dlt_1ab   = dlt_i1r
      dlt_1b    = dlt_r1
      th_0a     = th_i0
      th_0ab    = th_i0r
      th_0b     = th_r0
      th_1ab    = th_i1r
      th_1b     = th_r1
      !

    ! setting up extra coefficients
        ! remain coefficients
        rem_ci_cf = (1.0-rem_n_ci_min)/delt
        rem_hr_cf = (1.0-rem_n_hr_min)/delt
        ! enhanced melting
        k_enhm  = c_water/rlme

    ! setting up mask
    qcol_mask(2:i1,2:j1,1:k1) =   &
        q_ci_mask(2:i1,2:j1,1:k1).and.q_hr_mask(2:i1,2:j1,1:k1)

    ! setting up diameters and velocities
     D_a (2:i1,2:j1,1:k1) = D_ci (2:i1,2:j1,1:k1)
     v_a (2:i1,2:j1,1:k1) = v_ci (2:i1,2:j1,1:k1)
     D_b (2:i1,2:j1,1:k1) = D_hr (2:i1,2:j1,1:k1)
     v_b (2:i1,2:j1,1:k1) = v_hr (2:i1,2:j1,1:k1)


    ! -- inner part -------------------

    ! calculating
    do k=1,k1
    do j=2,j1
    do i=2,i1
      if ( qcol_mask(i,j,k)) then
       E_ab(i,j,k) =  E_coli
       !
       dq_col_b(i,j,k) = -(rhof(k)*pi/4)*E_ab(i,j,k)*n_ci(i,j,k)    &
            *q_hr(i,j,k)*(dlt_0a*D_a(i,j,k)**2                     &
              +dlt_1ab*D_a(i,j,k)*D_b(i,j,k)+dlt_1b*D_b(i,j,k)**2) &
            *( th_0a*v_a(i,j,k)**2-th_1ab*v_b(i,j,k)*v_a(i,j,k)    &
              +th_1b*v_b(i,j,k)**2+sigma_a**2+sigma_b**2)**0.5
       !
       dn_col_b(i,j,k) = -(rhof(k)*pi/4)*E_ab(i,j,k)*n_ci(i,j,k)   &
            *n_hr(i,j,k)*(dlt_0a*D_a(i,j,k)**2                     &
              +dlt_0ab*D_a(i,j,k)*D_b(i,j,k)+dlt_0b*D_b(i,j,k)**2) &
            *( th_0a*v_a(i,j,k)**2-th_0ab*v_b(i,j,k)*v_a(i,j,k)    &
              +th_0b*v_b(i,j,k)**2+sigma_a**2+sigma_b**2)**0.5
      endif
    enddo
    enddo
    enddo

    ! -- and the second part
      dlt_0ba   = dlt_r0i ! dlt_0ab   = dlt_i0r
      ! dlt_0b    = dlt_r0
      dlt_1ba   = dlt_r1i ! dlt_1ab   = dlt_i1r
      dlt_1a    = dlt_i1  ! dlt_1b    = dlt_r1
      th_1a     = th_i1 ! th_0a     = th_i0
      th_0ba    = th_r0i ! th_0ab    = th_i0r
      th_0b     = th_r0
      th_1ba    = th_r1i ! th_1ab    = th_i1r


    ! -- inner piece ----
    ! calculating
    do k=1,k1
    do j=2,j1
    do i=2,i1
      if ( qcol_mask(i,j,k)) then
       dq_col_a(i,j,k) = (rhof(k)*pi/4)*E_ab(i,j,k)*q_ci(i,j,k)    &
            *n_hr(i,j,k)*(dlt_1a*D_a(i,j,k)**2                     &
              +dlt_1ba*D_a(i,j,k)*D_b(i,j,k)+dlt_0b*D_b(i,j,k)**2) &
            *( th_1a*v_a(i,j,k)**2-th_1ba*v_b(i,j,k)*v_a(i,j,k)    &
              +th_0b*v_b(i,j,k)**2+sigma_a**2+sigma_b**2)**0.5
      endif
    enddo
    enddo
    enddo

    !---------------------

    ! -- outputs -----------------------
    dq_hr_col_ri (2:i1,2:j1,1:k1) = dq_col_b(2:i1,2:j1,1:k1)
    dn_ci_col_ri (2:i1,2:j1,1:k1) = dn_col_b(2:i1,2:j1,1:k1)
    dn_hr_col_ri (2:i1,2:j1,1:k1) = dn_col_b(2:i1,2:j1,1:k1)
    dq_ci_col_ri (2:i1,2:j1,1:k1) = -dq_col_a(2:i1,2:j1,1:k1)


    do k=1,k1
    do j=2,j1
    do i=2,i1
     if(qcol_mask(i,j,k)) then
        ! first adjustment
         dq_ci_col_ri(i,j,k) = max(dq_ci_col_ri(i,j,k),min(0.0,-svm(i,j,k,iq_ci)/delt-q_cip(i,j,k)))   ! following ICON, 2017
         dq_hr_col_ri(i,j,k) = max(dq_hr_col_ri(i,j,k),min(0.0,-svm(i,j,k,iq_hr)/delt-q_hrp(i,j,k)))   ! following ICON, 2017

        ! adjustment of numbers - both ice and water
         dn_ci_col_ri(i,j,k) = max(dn_ci_col_ri(i,j,k),min(0.0,-rem_ci_cf*svm(i,j,k,in_ci)-n_cip(i,j,k)))
         dn_ci_col_ri(i,j,k) = max(dn_ci_col_ri(i,j,k),min(0.0,-rem_hr_cf*svm(i,j,k,in_hr)-n_hrp(i,j,k)))
      if(tmp0(i,j,k).lt.T_3) then
        ! the collection is just riming
        ! decrease in numeber of raindrops same as decrease in number of ice
        dn_hr_col_ri(i,j,k) = dn_ci_col_ri(i,j,k)

        ! record the change in cloud ice
        q_cip (i,j,k) = q_cip(i,j,k) + dq_ci_col_ri (i,j,k)
        n_cip (i,j,k) = n_cip(i,j,k) + dn_ci_col_ri (i,j,k)

        ! change in rain
        q_hrp (i,j,k) = q_hrp(i,j,k) + dq_hr_col_ri (i,j,k)
        n_hrp (i,j,k) = n_hrp(i,j,k) + dn_hr_col_ri (i,j,k)

        ! and for graupel
        q_hgp (i,j,k) =q_hgp(i,j,k)-dq_ci_col_ri(i,j,k)-dq_hr_col_ri(i,j,k)
        n_hgp (i,j,k) =n_hgp(i,j,k)-dn_ci_col_ri(i,j,k)

        ! change in q_t - decrease in cloud ice
        qtpmcr(i,j,k) = qtpmcr(i,j,k) + 0.0 ! dq_ci_col_ri (i,j,k)

        ! change in th_l - release from freezing and removal of ice
        thlpmcr(i,j,k) = thlpmcr(i,j,k)                  &
            -(rlme/(cp*exnf(k)))*dq_hr_col_ri (i,j,k)
      else  ! tmp0(i,j,k).gt.T_3
        ! enhanced melting and graupel formation

        ! calculate the melting
        dq_ci_eme_ri(i,j,k) = k_enhm*(tmp0(i,j,k)-T_3)*dq_hr_col_ri(i,j,k) ! with + due to negative value of dq_hr here

        ! limit melting
        dq_ci_eme_ri(i,j,k) = max(dq_ci_eme_ri(i,j,k),min(0.0,-svm(i,j,k,iq_ci)/delt-q_cip(i,j,k)))

        ! calculate how many ice perticles melted
        ! q_ci here is always some small positive number
        dn_ci_eme_ri(i,j,k) = dq_ci_eme_ri(i,j,k)*n_ci(i,j,k)/(q_ci(i,j,k)+eps0)

        ! limit so it dos not melt more than interacting
        ! also limit so that new graupel not larger that max mean size of source ice ?
        dn_ci_eme_ri(i,j,k) =max(dn_ci_eme_ri(i,j,k),dn_ci_col_ri(i,j,k))

        ! update ice
        q_cip (i,j,k) = q_cip(i,j,k) + dq_ci_eme_ri(i,j,k) ! dq_ci_col_ri (i,j,k)
        n_cip (i,j,k) = n_cip(i,j,k) + dn_ci_eme_ri(i,j,k) ! dn_ci_col_ri (i,j,k)

        ! no collection of raindrops
        dq_hr_col_ri(i,j,k) = 0.0
        dq_ci_col_ri(i,j,k) = 0.0
        dn_hr_col_ri(i,j,k) = 0.0
        dn_ci_col_ri(i,j,k) = 0.0

        ! increase in rain mass
        q_hrp (i,j,k) = q_hrp(i,j,k) - dq_ci_eme_ri (i,j,k)

        ! change in thl - heat spent on melting
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlme/(cp*exnf(k)))*dq_ci_eme_ri (i,j,k)
      endif
     endif
    enddo
    enddo
    enddo


    ! deallocating
    deallocate (D_a, D_b, v_a, v_b, E_ab, E_stick, dn_col_b, dq_col_a, dq_col_b)
    deallocate (qcol_mask)

   end subroutine coll_rig3



 ! ***********************************
 !  riming of rain + snow to graupel
 !
 !
 ! ***********************************
   subroutine coll_rsg3

    use modglobal, only : ih,i1,jh,j1,k1,rv,rd, rlv,cp,pi
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,qvsl,qvsi,rhof,exnf,presf
    implicit none
    integer :: i,j,k
    real    :: sigma_a, sigma_b ,E_coli              &
               ,a_a, b_a, al_a, be_a, ga_a           &
               ,a_b, b_b, al_b, be_b, ga_b           &
               ,dlt_0a, dlt_0ab, dlt_0b              &
               ,dlt_1ab, dlt_1b                      &
               ,th_0a, th_0ab, th_0b                 &
               ,th_1ab, th_1b                        &
               ,dlt_0ba, dlt_1ba, dlt_1a             &
               ,th_0ba, th_1ba, th_1a
    real    :: dif_D_10, ntest, qtest, rem_hs_cf, rem_hr_cf, k_enhm
    real, allocatable, dimension(:,:,:) :: D_a, D_b, v_a, v_b ! <- later move outside of this subroutine
    real, allocatable, dimension(:,:,:) :: E_ab, E_stick, dn_col_b, dq_col_a, dq_col_b
    logical ,allocatable :: qcol_mask(:,:,:)

    ! start of the code

    ! allocate fields and fill
     allocate( D_a     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,D_b     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,v_a     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,v_b     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
             )

      allocate( E_ab    (2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,E_stick (2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,dn_col_b(2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,dq_col_a(2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,dq_col_b(2-ih:i1+ih,2-jh:j1+jh,k1)    &
             )
      allocate( qcol_mask(2-ih:i1+ih,2-jh:j1+jh,k1)   )

      D_a      = 0.0
      D_b      = 0.0
      v_a      = 0.0
      v_b      = 0.0
      E_ab     = 0.0
      E_stick  = 0.0
      dn_col_b = 0.0
      dq_col_a = 0.0
      dq_col_b = 0.0

    ! set constants
      sigma_a   = sigma_hs
      sigma_b   = sigma_hr
      a_a       = a_hs
      b_a       = b_hs
      al_a      = al_hs
      be_a      = be_hs
      ga_a      = ga_hs
      a_b       = a_hr
      b_b       = b_hr
      al_b      = al_hr
      be_b      = be_hr
      ga_b      = ga_hr
      E_coli    = E_s_m ! collision efficienty
      dlt_0a    = dlt_s0 ! dlt_s0
      dlt_0ab   = dlt_s0r ! dlt_s0i
      dlt_0b    = dlt_r0
      dlt_1ab   = dlt_s1r
      dlt_1b    = dlt_r1
      th_0a     = th_s0
      th_0ab    = th_s0r
      th_0b     = th_r0
      th_1ab    = th_s1r
      th_1b     = th_r1
      !

    ! setting up extra coefficients
        ! remain coefficients
        rem_hs_cf = (1.0-rem_n_hs_min)/delt
        rem_hr_cf = (1.0-rem_n_hr_min)/delt

        ! enhanced melting
        k_enhm  = c_water/rlme

    ! setting up mask
    qcol_mask(2:i1,2:j1,1:k1) =   &
        q_hs_mask(2:i1,2:j1,1:k1).and.q_hr_mask(2:i1,2:j1,1:k1)

    ! setting up diameters and velocities
     D_a (2:i1,2:j1,1:k1) = D_hs (2:i1,2:j1,1:k1)
     v_a (2:i1,2:j1,1:k1) = v_hs (2:i1,2:j1,1:k1)
     D_b (2:i1,2:j1,1:k1) = D_hr (2:i1,2:j1,1:k1)
     v_b (2:i1,2:j1,1:k1) = v_hr (2:i1,2:j1,1:k1)


    ! -- inner part -------------------

    ! calculating
    do k=1,k1
    do j=2,j1
    do i=2,i1
      if ( qcol_mask(i,j,k)) then
       E_ab(i,j,k) =  E_coli
       !
       dq_col_b(i,j,k) = -(rhof(k)*pi/4)*E_ab(i,j,k)*n_hs(i,j,k)    &
            *q_hr(i,j,k)*(dlt_0a*D_a(i,j,k)**2                     &
              +dlt_1ab*D_a(i,j,k)*D_b(i,j,k)+dlt_1b*D_b(i,j,k)**2) &
            *( th_0a*v_a(i,j,k)**2-th_1ab*v_b(i,j,k)*v_a(i,j,k)    &
              +th_1b*v_b(i,j,k)**2+sigma_a**2+sigma_b**2)**0.5
       !
       dn_col_b(i,j,k) = -(rhof(k)*pi/4)*E_ab(i,j,k)*n_hs(i,j,k)   &
            *n_hr(i,j,k)*(dlt_0a*D_a(i,j,k)**2                     &
              +dlt_0ab*D_a(i,j,k)*D_b(i,j,k)+dlt_0b*D_b(i,j,k)**2) &
            *( th_0a*v_a(i,j,k)**2-th_0ab*v_b(i,j,k)*v_a(i,j,k)    &
              +th_0b*v_b(i,j,k)**2+sigma_a**2+sigma_b**2)**0.5
      endif
    enddo
    enddo
    enddo

    ! -- and the second part
      dlt_0ba   = dlt_r0s ! dlt_0ab   = dlt_i0r
      ! dlt_0b    = dlt_r0
      dlt_1ba   = dlt_r1s ! dlt_1ab   = dlt_i1r
      dlt_1a    = dlt_s1  ! dlt_1b    = dlt_r1
      th_1a     = th_s1 ! th_0a     = th_i0
      th_0ba    = th_r0s ! th_0ab    = th_i0r
      th_0b     = th_r0
      th_1ba    = th_r1s ! th_1ab    = th_i1r


    ! -- inner piece ----
    ! calculating
    do k=1,k1
    do j=2,j1
    do i=2,i1
      if ( qcol_mask(i,j,k)) then
       !
       dq_col_a(i,j,k) = (rhof(k)*pi/4)*E_ab(i,j,k)*q_hs(i,j,k)    &
            *n_hr(i,j,k)*(dlt_1a*D_a(i,j,k)**2                     &
              +dlt_1ba*D_a(i,j,k)*D_b(i,j,k)+dlt_0b*D_b(i,j,k)**2) &
            *( th_1a*v_a(i,j,k)**2-th_1ba*v_b(i,j,k)*v_a(i,j,k)    &
              +th_0b*v_b(i,j,k)**2+sigma_a**2+sigma_b**2)**0.5
      endif
    enddo
    enddo
    enddo


    do k=1,k1
    do j=2,j1
    do i=2,i1
      if(qcol_mask(i,j,k)) then
        ! first adjustment
        dq_hs_col_rs(i,j,k) = max(-dq_col_a(i,j,k),min(0.0,-svm(i,j,k,iq_hs)/delt-q_hsp(i,j,k)))   ! following ICON, 2017
        dq_hr_col_rs(i,j,k) = max(dq_col_b(i,j,k),min(0.0,-svm(i,j,k,iq_hr)/delt-q_hrp(i,j,k)))   ! following ICON, 2017

        ! adjustment of numbers - both ice and snow
        dn_hs_col_rs(i,j,k) = max(dn_col_b(i,j,k) ,min(0.0,-rem_hs_cf*svm(i,j,k,in_hs)-n_hsp(i,j,k)))
        dn_hs_col_rs(i,j,k) = max(dn_col_b(i,j,k) ,min(0.0,-rem_hr_cf*svm(i,j,k,in_hr)-n_hrp(i,j,k)))

         if (tmp0(i,j,k).lt.T_3) then
           ! and copying it to the second one
           dn_hr_col_rs(i,j,k) = dn_hs_col_rs(i,j,k)

           ! record the change in cloud ice
           q_hsp (i,j,k) = q_hsp(i,j,k) + dq_hs_col_rs (i,j,k)
           n_hsp (i,j,k) = n_hsp(i,j,k) + dn_hs_col_rs (i,j,k)

           ! change in rain
           q_hrp (i,j,k) = q_hrp(i,j,k) + dq_hr_col_rs (i,j,k)
           n_hrp (i,j,k) = n_hrp(i,j,k) + dn_hr_col_rs (i,j,k)

           ! and for graupel
           q_hgp (i,j,k) =q_hgp(i,j,k)-dq_hs_col_rs(i,j,k)-dq_hr_col_rs(i,j,k)
           n_hgp (i,j,k) =n_hgp(i,j,k)-dn_hs_col_rs(i,j,k)

           ! change in th_l - release from freezing and removal of ice
           thlpmcr(i,j,k) = thlpmcr(i,j,k)                  &
                -(rlme/(cp*exnf(k)))*dq_hr_col_rs (i,j,k)
         else  ! tmp0(i,j,k).gt.T_3
           ! enhanced melting and graupel formation
           ! calculate the melting
           ! with + due to negative value of dq_hr here
           dq_hs_eme_rs(i,j,k) = k_enhm*(tmp0(i,j,k)-T_3)*dq_hr_col_rs(i,j,k)

           ! snow melting
           dq_hs_eme_rs(i,j,k) = max(dq_hs_eme_rs(i,j,k),min(0.0,-svm(i,j,k,iq_hs)/delt-q_hsp(i,j,k)))

           ! calculate how many snow perticles melted
           ! q_hs here is always some small positive number
           dn_hs_eme_rs(i,j,k) = dq_hs_eme_rs(i,j,k)*n_hs(i,j,k)/q_hs(i,j,k)

           ! limit so it dos not melt more than interacting
           ! also limit so that new graupel not larger that max mean size of source ice ?
           dn_hs_eme_rs(i,j,k) =max(dn_hs_eme_rs(i,j,k),dn_hs_col_rs(i,j,k))

           ! update snow
           q_hsp (i,j,k) = q_hsp(i,j,k) + dq_hs_eme_rs(i,j,k)
           n_hsp (i,j,k) = n_hsp(i,j,k) + dn_hs_eme_rs(i,j,k)

           ! no collection of raindrops
           dq_hr_col_rs(i,j,k) = 0.0
           dq_hs_col_rs(i,j,k) = 0.0
           dn_hr_col_rs(i,j,k) = 0.0
           dn_hs_col_rs(i,j,k) = 0.0

           ! increase in rain mass
           q_hrp (i,j,k) = q_hrp(i,j,k) - dq_hs_eme_rs (i,j,k)

           ! change in thl - heat spent on melting
           thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlme/(cp*exnf(k)))*dq_hs_eme_rs (i,j,k)
       endif
      endif
    enddo
    enddo
    enddo

    deallocate (D_a, D_b, v_a, v_b, E_ab, E_stick, dn_col_b, dq_col_a, dq_col_b)
    deallocate (qcol_mask)

   end subroutine coll_rsg3


! ****************************************
!  riming of graupel
!
!
! ****************************************
    subroutine coll_grg3

    use modglobal, only : ih,i1,jh,j1,k1,rv,rd, rlv,cp,pi
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,qvsl,qvsi,rhof,exnf,presf
    implicit none
    integer :: i,j,k
    real    :: sigma_a, sigma_b ,E_coli              &
               ,a_a, b_a, al_a, be_a, ga_a           &
               ,a_b, b_b, al_b, be_b, ga_b           &
               ,dlt_0a, dlt_0ab, dlt_0b              &
               ,dlt_1ab, dlt_1b                      &
               ,th_0a, th_0ab, th_0b                 &
               ,th_1ab, th_1b
    real    :: dif_D_10, ntest, qtest, rem_cf, k_enhm
    real, allocatable, dimension(:,:,:) :: D_a, D_b, v_a, v_b ! <- later move outside of this subroutine
    real, allocatable, dimension(:,:,:) :: E_ab, E_stick, dn_col_b, dq_col_a
    logical ,allocatable :: qcol_mask(:,:,:)

    ! start of the code

    ! allocate fields and fill
     allocate( D_a     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,D_b     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,v_a     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
              ,v_b     (2-ih:i1+ih,2-jh:j1+jh,k1)     &
             )

      allocate( E_ab    (2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,E_stick (2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,dn_col_b(2-ih:i1+ih,2-jh:j1+jh,k1)    &
               ,dq_col_a(2-ih:i1+ih,2-jh:j1+jh,k1)    &
             )
      allocate( qcol_mask(2-ih:i1+ih,2-jh:j1+jh,k1)   )

      D_a      = 0.0
      D_b      = 0.0
      v_a      = 0.0
      v_b      = 0.0
      E_ab     = 0.0
      E_stick  = 0.0
      dn_col_b = 0.0
      dq_col_a = 0.0

    ! set constants
      sigma_a   = sigma_hg
      sigma_b   = sigma_hr
      a_a       = a_hg
      b_a       = b_hg
      al_a      = al_hg
      be_a      = be_hg
      ga_a      = ga_hg
      a_b       = a_hr
      b_b       = b_hr
      al_b      = al_hr
      be_b      = be_hr
      ga_b      = ga_hr
      E_coli    = E_er_m ! collision efficienty
      dlt_0a    = dlt_g0 ! dlt_s0
      dlt_0ab   = dlt_g0r ! dlt_s0i
      dlt_0b    = dlt_r0
      dlt_1ab   = dlt_g1r
      dlt_1b    = dlt_r1
      th_0a     = th_g0
      th_0ab    = th_g0r
      th_0b     = th_r0
      th_1ab    = th_g1r
      th_1b     = th_r1
      !
      dif_D_10  = D_c_b-D_c_a   !< denominator in calculationg collision efficiency
      !
      ! remaining number of particles
       rem_cf = (1.0-rem_n_hr_min)/delt
      !
      ! enhanced melting
       k_enhm  = c_water/rlme


    ! setting up mask
    qcol_mask(2:i1,2:j1,1:k1) =   &
        q_hg_mask(2:i1,2:j1,1:k1).and.q_hr_mask(2:i1,2:j1,1:k1)

    ! setting up diameters and velocities
     D_a (2:i1,2:j1,1:k1) = D_hg (2:i1,2:j1,1:k1)
     v_a (2:i1,2:j1,1:k1) = v_hg (2:i1,2:j1,1:k1)
     D_b (2:i1,2:j1,1:k1) = D_hr (2:i1,2:j1,1:k1)
     v_b (2:i1,2:j1,1:k1) = v_hr (2:i1,2:j1,1:k1)



    ! -- inner part -------------------


    ! calculating
    do k=1,k1
    do j=2,j1
    do i=2,i1
      if ( qcol_mask(i,j,k)) then
       E_ab(i,j,k) =  E_coli
       !
       dq_col_a(i,j,k) = (rhof(k)*pi/4)*E_ab(i,j,k)*n_hg(i,j,k)    &
            *q_hr(i,j,k)*(dlt_0a*D_a(i,j,k)**2                     &
              +dlt_1ab*D_a(i,j,k)*D_b(i,j,k)+dlt_1b*D_b(i,j,k)**2) &
            *( th_0a*v_a(i,j,k)**2-th_1ab*v_b(i,j,k)*v_a(i,j,k)    &
              +th_1b*v_b(i,j,k)**2+sigma_a**2+sigma_b**2)**0.5
       !
       dn_col_b(i,j,k) = -(rhof(k)*pi/4)*E_ab(i,j,k)*n_hg(i,j,k)   &
            *n_hr(i,j,k)*(dlt_0a*D_a(i,j,k)**2                     &
              +dlt_0ab*D_a(i,j,k)*D_b(i,j,k)+dlt_0b*D_b(i,j,k)**2) &
            *( th_0a*v_a(i,j,k)**2-th_0ab*v_b(i,j,k)*v_a(i,j,k)    &
              +th_0b*v_b(i,j,k)**2+sigma_a**2+sigma_b**2)**0.5
      endif
    enddo
    enddo
    enddo

    ! -- outputs -----------------------
    ! dq_hghr_rime  (2:i1,2:j1,1:k1) = dq_col_a(2:i1,2:j1,1:k1)
    ! dn_hr_rime_hg (2:i1,2:j1,1:k1) = dn_col_b(2:i1,2:j1,1:k1)

    do k=1,k1
    do j=2,j1
    do i=2,i1
     if(qcol_mask(i,j,k)) then
      if (tmp0(i,j,k).lt.T_3) then
       ! riming only
         dq_hghr_rime(i,j,k) = min(dq_col_a(i,j,k),max(0.0,svm(i,j,k,iq_hr)/delt+q_hrp(i,j,k)))   ! following ICON, 2017
        ! = min(dq_hghr_rime(i,j,k),svm(i,j,k,iq_hr)/delt)
         dn_hr_rime_hg(i,j,k)= max(dn_col_b(i,j,k),min(0.0,-rem_cf*svm(i,j,k,in_hr)-n_hrp(i,j,k)))
        ! =min(dn_hr_rime_hg(i,j,k),-svm(i,j,k,in_hr)/delt)
        ! change in the amount of graupel
         q_hgp (i,j,k) = q_hgp(i,j,k) + dq_hghr_rime (i,j,k)
        ! change in the amount of rain
         n_hrp (i,j,k) = n_hrp(i,j,k) + dn_hr_rime_hg (i,j,k)
         q_hrp (i,j,k) = q_hrp(i,j,k) - dq_hghr_rime (i,j,k)
        !
        ! no change in q_t
        ! qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_col_a (i,j,k)
        ! change in th_l - just heat release from freezing
         thlpmcr(i,j,k) = thlpmcr(i,j,k)+ (rlme/(cp*exnf(k)))*dq_hghr_rime (i,j,k)
      else
       ! not riming,but enhanced melting and scheding
       ! calculating the melting
        dq_hg_eme_gr(i,j,k) = -k_enhm*(tmp0(i,j,k)-T_3)*min(dq_col_a(i,j,k),q_hr(i,j,k)/delt)
        dq_hg_eme_gr(i,j,k) = max(dq_hg_eme_gr(i,j,k),min(0.0,-svm(i,j,k,iq_hg)/delt-q_hgp(i,j,k)))
       ! calculating number of melted particles
        ! - expected to be proportional to melted mass
        dn_hg_eme_gr(i,j,k) = dq_hg_eme_gr(i,j,k)*n_hg(i,j,k)/(q_hg(i,j,k)+eps0) ! q_hg here is always some small positive number
        ! - but not more than number of interacting particles
        ! dn_hg_eme_gr(i,j,k) = max(dn_hg_eme_gr(i,j,k), max(dn_col_b(i,j,k),-n_hr(i,j,k)/delt))
        ! - and not more than total number of particles
        dn_hg_eme_gr(i,j,k) = max(dn_hg_eme_gr(i,j,k), -min((svm(i,j,k,in_hr)/delt+n_hrp(i,j,k)),(svm(i,j,k,in_hg)/delt+n_hgp(i,j,k))))
       ! updating tendencies
       ! updating rain
         q_hrp (i,j,k) = q_hrp(i,j,k) - dq_hg_eme_gr (i,j,k)
         ! n_hrp (i,j,k) = n_hrp(i,j,k) + dn_hr_rime_hg (i,j,k)
        ! graupel
         q_hgp (i,j,k) = q_hgp(i,j,k) + dq_hg_eme_gr (i,j,k)
         n_hgp (i,j,k) = n_hgp(i,j,k) + dn_hg_eme_gr (i,j,k)
        ! updating cloud water
         ! no change - all turned into rain
       ! updating thermodynamic
        ! qt : no change -  increased by melted water goes to rain
        ! qtpmcr(i,j,k) = qtpmcr(i,j,k) +0.0
        ! thl : melting
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlme/(cp*exnf(k)))*dq_hg_eme_gr(i,j,k)
      endif
     endif
    enddo
    enddo
    enddo

    ! #hh checking the sizes
    if (l_sb_dbg) then
     if(any(( svm(2:i1,2:j1,1:k1,iq_hr)-delt*dq_col_a(2:i1,2:j1,1:k1) ).lt. 0.0 )) then
      write(6,*) 'WARNING: coll_grg3 too high'
      write(6,*) ' removing more rain water then there in ', count(( svm(2:i1,2:j1,1:k1,iq_hr)-delt*dq_hghr_rime(2:i1,2:j1,1:k1) ).lt. 0.0 )
      write(6,*) ' removing too much rain water in ', count(( q_hr(2:i1,2:j1,1:k1)-delt*dq_hghr_rime(2:i1,2:j1,1:k1) ).lt. 0.0 )
      write(6,*) ' getting negative q_t in  ', count(( qt0(2:i1,2:j1,1:k1)+delt*q_hrp(2:i1,2:j1,1:k1) ).lt. 0.0 )
     endif
     if(any(( svm(2:i1,2:j1,1:k1,in_hr)+delt*dn_col_b(2:i1,2:j1,1:k1) ).lt. 0.0 )) then
      write(6,*) 'WARNING: coll_grg too high'
      write(6,*) ' removing more raindrops than available in gridpoints ', count(( svm(2:i1,2:j1,1:k1,in_hr)+delt*dn_hr_rime_hg(2:i1,2:j1,1:k1) ).lt. 0.0 )
      write(6,*) ' removing too many raindrops in gridpoints ', count(( n_hr(2:i1,2:j1,1:k1)+delt*dn_hr_rime_hg(2:i1,2:j1,1:k1) ).lt. 0.0 )
     endif
    endif


     ! deallocating
     deallocate (D_a, D_b, v_a, v_b, E_ab, E_stick, dn_col_b, dq_col_a)
     deallocate (qcol_mask)

   end subroutine coll_grg3

! ****************************************
!  partial conversion to graupel
!
!  - based on Seifgert and Beheng, 2004
!       - the implementation differs from the implementation in ICON
!
!
!  - have to be called after enhanced melting subroutine
!
! ****************************************
    subroutine conv_partial3

    use modglobal, only : ih,i1,jh,j1,k1,rv,rd, rlv,cp,pi,rhow
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,qvsl,qvsi,rhof,exnf,presf
    implicit none
    integer :: i,j,k
    real    ::pi6rhoe, cc_ci, cc_hs, rem_ci_cf, rem_hs_cf
    logical ,allocatable, dimension(:,:,:) :: qcol_mask_i      &
         ,qcol_mask_s ! , qcv_mask_i, qcv_mask_s


      allocate(  qcol_mask_i(2-ih:i1+ih,2-jh:j1+jh,k1) &
                ,qcol_mask_s(2-ih:i1+ih,2-jh:j1+jh,k1) &
              )
              !  ,qcv_mask_i (2-ih:i1+ih,2-jh:j1+jh,k1) &
              !  ,qcv_mask_s (2-ih:i1+ih,2-jh:j1+jh,k1) &
              !)

    ! calculate constants
    pi6rhoe = (pi/6.0)*rhoeps
    cc_ci   = al_0ice*rhow/rhoeps
    cc_hs   = al_0snow*rhow/rhoeps

    ! remain coefficients
    rem_ci_cf = (1.0-rem_n_min_cv)/delt
    rem_hs_cf = (1.0-rem_n_min_cv)/delt

    ! setting values
    ! D_ci = 0.0
    ! D_hs = 0.0
    ! G_ci = 0.0
    ! G_hs = 0.0
    dq_ci_cv = 0.0
    dq_hs_cv = 0.0

    ! set up masks
    qcol_mask_i(2:i1,2:j1,1:k1) =  &
        q_ci_mask(2:i1,2:j1,1:k1).and.q_cl_mask(2:i1,2:j1,1:k1)
    qcol_mask_s(2:i1,2:j1,1:k1) =  &
        q_hs_mask(2:i1,2:j1,1:k1).and.q_cl_mask(2:i1,2:j1,1:k1)


    ! inner term for ice conversion
    do k=1,k1
    do j=2,j1
    do i=2,i1
     if ( qcol_mask_i(i,j,k)) then
      if ((D_ci(i,j,k).gt.D_mincv_ci).and.(tmp0(i,j,k).lt.T_3)) then
       dq_ci_cv(i,j,k) = cc_ci*(pi6rhoe*D_ci(i,j,k)**3 /x_ci(i,j,k)-1.0)
        ! ? not to exceed conversion rate 1 ?
        ! G_ci(i,j,k) = min( 1.0, G_ci(i,j,k))
       dq_ci_cv(i,j,k) = -dq_ci_rime(i,j,k)/dq_ci_cv(i,j,k)
       dq_ci_cv(i,j,k) = max( dq_ci_cv(i,j,k),min(0.0,-svm(i,j,k,iq_ci)/delt-q_cip(i,j,k)))  ! based on ICON, 2017
       ! = max( dq_ci_cv(i,j,k),-svm(i,j,k,iq_ci)/delt)
       dn_ci_cv(i,j,k) = dq_ci_cv(i,j,k)/max(x_ci(i,j,k),x_ci_cvmin)
       ! = dq_ci_cv(i,j,k)/max(x_ci(i,j,k),x_ci_cvmin)
       dn_ci_cv(i,j,k) = max(dn_ci_cv(i,j,k),min(0.0,-rem_ci_cf*svm(i,j,k,in_ci)-n_cip(i,j,k)))
       ! = max(dn_hs_cv(i,j,k), -svm(i,j,k,in_ci)/delt)
       ! change in the amount of graupel
        n_hgp (i,j,k) = n_hgp(i,j,k)-dn_ci_cv (i,j,k)
        q_hgp (i,j,k) = q_hgp(i,j,k)-dq_ci_cv (i,j,k)
       ! change in the amount of cloud ice
        n_cip (i,j,k) = n_cip(i,j,k) + dn_ci_cv (i,j,k)
        q_cip (i,j,k) = q_cip(i,j,k) + dq_ci_cv (i,j,k)
       ! change in q_t -- in case of i->g
       !#iceout qtpmcr(i,j,k) = qtpmcr(i,j,k) + dq_ci_cv (i,j,k)
       ! change in th_l - in case of i->g
       !#iceout thlpmcr(i,j,k) = thlpmcr(i,j,k)- (rlv/(cp*exnf(k)))*dq_ci_cv(i,j,k)
      endif
     endif
    enddo
    enddo
    enddo

    ! inner term for snow conversion
    do k=1,k1
    do j=2,j1
    do i=2,i1
     if ( qcol_mask_s(i,j,k)) then
      if ((D_hs(i,j,k).gt.D_mincv_hs).and.(tmp0(i,j,k).lt.T_3)) then
       dq_hs_cv(i,j,k) = cc_hs*(pi6rhoe*D_hs(i,j,k)**3 /x_hs(i,j,k)-1)
       ! ? at the sam time, the value should be limited
       ! ? not to exceed conversion rate 1
       ! G_hs(i,j,k) = min( 1.0, G_hs(i,j,k))
       dq_hs_cv(i,j,k) = -dq_hs_rime(i,j,k)/dq_hs_cv(i,j,k)
       ! correction - not removing more than available
       ! dq_hs_cv(i,j,k) = max(dq_hs_cv(i,j,k),-q_hs(i,j,k)/delt )
       ! basic correction of the tendency
       dq_hs_cv(i,j,k) = max( dq_hs_cv(i,j,k),min(0.0,-svm(i,j,k,iq_hs)/delt-q_hsp(i,j,k)))
       ! dq_hs_cv(i,j,k) = max( dq_hs_cv(i,j,k),-svm(i,j,k,iq_hs)/delt)
       dn_hs_cv(i,j,k) = dq_hs_cv(i,j,k)/max(x_hs(i,j,k),x_hs_cvmin)
       ! dn_hs_cv(i,j,k) = dq_hs_cv(i,j,k)/max(x_hs(i,j,k),x_hs_cvmin)
       dn_hs_cv(i,j,k) = max(dn_hs_cv(i,j,k),min(0.0,-rem_hs_cf*svm(i,j,k,in_hs)-n_hsp(i,j,k)))
       ! dn_hs_cv(i,j,k) = max(dn_hs_cv(i,j,k), -svm(i,j,k,in_hs)/delt)
       ! and the second correction of the q tendency
      ! change in the amount of graupel
      n_hgp (i,j,k) = n_hgp(i,j,k)-dn_hs_cv (i,j,k)
      q_hgp (i,j,k) = q_hgp(i,j,k)-dq_hs_cv (i,j,k)
      ! change in the amount of snow
      n_hsp (i,j,k) = n_hsp(i,j,k) + dn_hs_cv (i,j,k)
      q_hsp (i,j,k) = q_hsp(i,j,k) + dq_hs_cv (i,j,k)
      ! change in q_t    -- none
      ! change in th_l   -- none
      endif
     endif
    enddo
    enddo
    enddo


    ! warnings
    if (l_sb_dbg) then
     if(any(( svm(2:i1,2:j1,1:k1,iq_ci)+delt*dq_ci_cv(2:i1,2:j1,1:k1) ).lt. 0 )) then
      write(6,*) 'WARNING: conv_partial3 removing too much ice'
      write(6,*) ' removing more ice than available in ', count((svm(2:i1,2:j1,1:k1,iq_ci) +delt*dq_ci_cv(2:i1,2:j1,1:k1) ).lt. 0 )
      write(6,*) ' removing too much ice in ', count(( q_ci(2:i1,2:j1,1:k1)+delt*dq_ci_cv(2:i1,2:j1,1:k1) ).lt. 0 )
      write(6,*) ' getting negative q_t in  ', count(( qt0(2:i1,2:j1,1:k1)+delt*dq_ci_cv(2:i1,2:j1,1:k1) ).lt. 0 )
     endif
     if(any(( svm(2:i1,2:j1,1:k1,iq_hs)+delt*dq_hs_cv(2:i1,2:j1,1:k1) ).lt. 0 )) then
      write(6,*) 'WARNING: conv_partial3 removing too much snow'
      write(6,*) ' removing more snow than available in ', count((svm(2:i1,2:j1,1:k1,iq_hs) +delt*dq_hs_cv(2:i1,2:j1,1:k1) ).lt. 0 )
      write(6,*) ' removing too much snow in ', count(( q_hs(2:i1,2:j1,1:k1)+delt*dq_hs_cv(2:i1,2:j1,1:k1) ).lt. 0 )
      ! write(6,*) ' getting negative q_t in  ', count(( qt0(2:i1,2:j1,1:k1)+delt*q_ci_cv(2:i1,2:j1,1:k1) ).lt. 0 )
     endif
    endif

     ! deallocating
     deallocate (qcol_mask_i, qcol_mask_s) ! deallocate (qcol_mask_i, qcol_mask_s, qcv_mask_i, qcv_mask_s)

   end subroutine conv_partial3

! ***************************************************************
! Melting and evaporation of ice particles
!
! - wrapper
!   - calls separately melting processes for each of species
!
!
! ***************************************************************
   subroutine evapmelting3

    use modglobal, only : ih,i1,jh,j1,k1,rv,rd, rlv,cp,pi,rhow
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,qvsl,qvsi,rhof,exnf,presf
    implicit none
    integer :: i,j,k
    real    ::pi6rhoe
    ! real, allocatable, dimension(:,:,:) :: dq_ci_spl, dn_ci_spl, dq_hg_temp
    ! logical ,allocatable, dimension(:,:,:) :: qcol_mask


    ! allocate( )

    ! ------------------------------------
    ! calling separate melting processes
    ! ------------------------------------

    ! ice
    call sb_evmelt3(aven_0i,aven_1i,bven_0i,bven_1i,x_ci_bmin        &
       ,n_ci,q_ci,x_ci,D_ci,v_ci,dq_ci_me,dn_ci_me,dq_ci_ev,dn_ci_ev)
     ! loss of ice
     n_cip (2:i1,2:j1,1:k1) = n_cip(2:i1,2:j1,1:k1)+ dn_ci_me(2:i1,2:j1,1:k1)+dn_ci_ev(2:i1,2:j1,1:k1)
     q_cip (2:i1,2:j1,1:k1) = q_cip(2:i1,2:j1,1:k1)+ dq_ci_me(2:i1,2:j1,1:k1)+dq_ci_ev(2:i1,2:j1,1:k1)
     ! transformed to rain or cloud
     ! for now into rain - improve possibly later
     n_hrp (2:i1,2:j1,1:k1) = n_hrp(2:i1,2:j1,1:k1)- dn_ci_me(2:i1,2:j1,1:k1)
     q_hrp (2:i1,2:j1,1:k1) = q_hrp(2:i1,2:j1,1:k1)- dq_ci_me(2:i1,2:j1,1:k1)
     ! transfomed to water vapour
     qtpmcr(2:i1,2:j1,1:k1) = qtpmcr(2:i1,2:j1,1:k1)-dq_ci_ev(2:i1,2:j1,1:k1)
     ! and heat production : heat spent on melting and evaporation - done lower

    ! snow
    call sb_evmelt3(aven_0s,aven_1s,bven_0s,bven_1s,x_hs_bmin        &
       ,n_hs,q_hs,x_hs,D_hs,v_hs,dq_hs_me,dn_hs_me,dq_hs_ev,dn_hs_ev)
     ! loss of snow
     n_hsp (2:i1,2:j1,1:k1) = n_hsp(2:i1,2:j1,1:k1)+ dn_hs_me(2:i1,2:j1,1:k1)+dn_hs_ev(2:i1,2:j1,1:k1)
     q_hsp (2:i1,2:j1,1:k1) = q_hsp(2:i1,2:j1,1:k1)+ dq_hs_me(2:i1,2:j1,1:k1)+dq_hs_ev(2:i1,2:j1,1:k1)
     ! transformed to rain
     n_hrp (2:i1,2:j1,1:k1) = n_hrp(2:i1,2:j1,1:k1)- dn_hs_me(2:i1,2:j1,1:k1)
     q_hrp (2:i1,2:j1,1:k1) = q_hrp(2:i1,2:j1,1:k1)- dq_hs_me(2:i1,2:j1,1:k1)
     ! transfomed to water vapour
     qtpmcr(2:i1,2:j1,1:k1) = qtpmcr(2:i1,2:j1,1:k1)-dq_hs_ev(2:i1,2:j1,1:k1)
     ! and heat production : heat spent on melting and evaporation - done lower

    !
    ! graupel
    call sb_evmelt3(aven_0g,aven_1g,bven_0g,bven_1g,x_hg_bmin         &
      ,n_hg,q_hg,x_hg,D_hg,v_hg,dq_hg_me,dn_hg_me,dq_hg_ev,dn_hg_ev)
     ! correction
     !write(6,*) '   n_hg(i,j,k), q_hg(i,j,k), dn_hg_ev(i,j,k), dn_hg_me (i,j,k)' ! #b2t17
     do k=1,k1
     do j=2,j1
     do i=2,i1
       if (dn_hg_me(i,j,k).lt.0.0) then
          ! limiting not to remove more water than available
          ! dn_hg_me (i,j,k) = max(dn_hg_me (i,j,k),(-svm(i,j,k,in_hg)/delt-n_hgp(i,j,k)) )
          dn_hg_ev(i,j,k) = max(dn_hg_ev(i,j,k),(-svm(i,j,k,in_hg)/delt-n_hgp(i,j,k)) )
          ! dq_hg_me (i,j,k) = max(dq_hg_me (i,j,k),(-svm(i,j,k,iq_hg)/delt-q_hgp(i,j,k)) )
          dq_hg_ev(i,j,k) = max(dq_hg_ev(i,j,k),(-svm(i,j,k,iq_hg)/delt-q_hgp(i,j,k)) )
          dn_hg_me (i,j,k) = max(dn_hg_me (i,j,k),(-svm(i,j,k,in_hg)/delt-n_hgp(i,j,k)) )
          dq_hg_me (i,j,k) = max(dq_hg_me (i,j,k),(-svm(i,j,k,iq_hg)/delt-q_hgp(i,j,k)) )
          !
       endif
     enddo
     enddo
     enddo

     ! loss of graupel
     n_hgp (2:i1,2:j1,1:k1) = n_hgp(2:i1,2:j1,1:k1)+ dn_hg_me(2:i1,2:j1,1:k1)+ dn_hg_ev(2:i1,2:j1,1:k1)
     q_hgp (2:i1,2:j1,1:k1) = q_hgp(2:i1,2:j1,1:k1)+ dq_hg_me(2:i1,2:j1,1:k1)+ dq_hg_ev(2:i1,2:j1,1:k1)
     ! transformed to rain
     !
     n_hrp (2:i1,2:j1,1:k1) = n_hrp(2:i1,2:j1,1:k1)- dn_hg_me(2:i1,2:j1,1:k1)
     q_hrp (2:i1,2:j1,1:k1) = q_hrp(2:i1,2:j1,1:k1)- dq_hg_me(2:i1,2:j1,1:k1)
     ! transfomed to water vapour
     qtpmcr(2:i1,2:j1,1:k1) = qtpmcr(2:i1,2:j1,1:k1)-dq_hg_ev(2:i1,2:j1,1:k1)
    !

     ! and heat production : heat spent on melting and evaporation - done lower
     !  - melting goes to rain
     !  - evaporation goes water vapour
     do k=1,k1
     do j=2,j1
     do i=2,i1
       thlpmcr(i,j,k) = thlpmcr(i,j,k)+                            &
          (rlme/(cp*exnf(k)))*(dq_ci_me(i,j,k)+dq_hs_me(i,j,k)+    &
            dq_hg_me(i,j,k))+                                      &
          ((rlv+rlme)/(cp*exnf(k)))*(dq_ci_ev(i,j,k)+             &
            dq_hs_ev(i,j,k)+ dq_hg_ev(i,j,k))
     enddo
     enddo
     enddo
    ! clean-up
    ! deallocate( )

   end subroutine evapmelting3
! ice_melt3(avent0,avent1,bvent0,bvent1,n_e,q_e,x_e,D_e,v_e,dq_e_me,dn_e_me)



! ***************************************************************
! Ice multiplication
!
! - wrapper
!   - calls separately H-M process for each of the riming processes
!
!
! ***************************************************************
   subroutine ice_multi3

    use modglobal, only : ih,i1,jh,j1,k1,rv,rd, rlv,cp,pi,rhow
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,qvsl,qvsi,rhof,exnf,presf
    implicit none
    integer :: i,j,k
    real    ::pi6rhoe
    real, allocatable, dimension(:,:,:) :: dq_ci_spl, dn_ci_spl, dq_hg_temp
    ! logical ,allocatable, dimension(:,:,:) :: qcol_mask


    ! allocate
    allocate( dq_hg_temp (2-ih:i1+ih,2-jh:j1+jh,k1)   &
             ,dq_ci_spl  (2-ih:i1+ih,2-jh:j1+jh,k1)   &
             ,dn_ci_spl  (2-ih:i1+ih,2-jh:j1+jh,k1)   &
            )

    dq_hg_temp = 0.0
    dq_ci_spl  = 0.0
    dn_ci_spl  = 0.0


    ! ------------------------------------
    ! calling separate H-M processes
    ! ------------------------------------


    ! i+l -> i
    call hallet_mossop3 (dq_ci_rime,q_ci,dq_ci_spl,dn_ci_spl)
     dn_ci_mul(2:i1,2:j1,1:k1) = dn_ci_mul(2:i1,2:j1,1:k1) + dn_ci_spl(2:i1,2:j1,1:k1)
     ! no change in ice content
     ! dq_ci_rime(2:i1,2:j1,1:k1) = dq_ci_rime(2:i1,2:j1,1:k1)-dq_ci_spl(2:i1,2:j1,1:k1)
     ! no change in ice content

    ! s+l -> s
    call hallet_mossop3 (dq_hs_rime,q_hs,dq_ci_spl,dn_ci_spl)
     dn_ci_mul (2:i1,2:j1,1:k1) = dn_ci_mul(2:i1,2:j1,1:k1)+ dn_ci_spl(2:i1,2:j1,1:k1)
     dq_ci_mul (2:i1,2:j1,1:k1) = dq_ci_mul(2:i1,2:j1,1:k1)+ dq_ci_spl(2:i1,2:j1,1:k1)
     ! dq_hs_rime(2:i1,2:j1,1:k1) = dq_hs_rime(2:i1,2:j1,1:k1)- dq_ci_spl(2:i1,2:j1,1:k1)
     q_hsp     (2:i1,2:j1,1:k1) = q_hsp(2:i1,2:j1,1:k1)    - dq_ci_spl(2:i1,2:j1,1:k1)  ! effect on snow

    ! g+l -> g
    call hallet_mossop3 (dq_hg_rime,q_hg,dq_ci_spl,dn_ci_spl)
     dn_ci_mul (2:i1,2:j1,1:k1) = dn_ci_mul(2:i1,2:j1,1:k1)+ dn_ci_spl(2:i1,2:j1,1:k1)
     dq_ci_mul (2:i1,2:j1,1:k1) = dq_ci_mul(2:i1,2:j1,1:k1)+ dq_ci_spl(2:i1,2:j1,1:k1)
     ! dq_hg_rime(2:i1,2:j1,1:k1) = dq_hg_rime(2:i1,2:j1,1:k1)    - dq_ci_spl(2:i1,2:j1,1:k1)
     q_hgp     (2:i1,2:j1,1:k1) = q_hgp(2:i1,2:j1,1:k1)    - dq_ci_spl(2:i1,2:j1,1:k1)  ! effect on snow

    ! g+r -> g
    call hallet_mossop3 (dq_hghr_rime,q_hg,dq_ci_spl,dn_ci_spl)
     dn_ci_mul(2:i1,2:j1,1:k1) = dn_ci_mul(2:i1,2:j1,1:k1)+ dn_ci_spl(2:i1,2:j1,1:k1)
     dq_ci_mul(2:i1,2:j1,1:k1) = dq_ci_mul(2:i1,2:j1,1:k1)+ dq_ci_spl(2:i1,2:j1,1:k1)
     ! dq_hghr_rime(2:i1,2:j1,1:k1) = dq_hghr_rime(2:i1,2:j1,1:k1)    - dq_ci_spl(2:i1,2:j1,1:k1)
     q_hgp    (2:i1,2:j1,1:k1) = q_hgp(2:i1,2:j1,1:k1)    - dq_ci_spl(2:i1,2:j1,1:k1)  ! effect on snow

    ! s+r -> g  -- does it make sense to include it?
    dq_hg_temp(2:i1,2:j1,1:k1) = -dq_hs_col_rs(2:i1,2:j1,1:k1)-dq_hr_col_rs(2:i1,2:j1,1:k1)
    call hallet_mossop3 (dq_hg_temp ,q_hg,dq_ci_spl,dn_ci_spl)
     dn_ci_mul(2:i1,2:j1,1:k1) = dn_ci_mul(2:i1,2:j1,1:k1)+ dn_ci_spl(2:i1,2:j1,1:k1)
     dq_ci_mul(2:i1,2:j1,1:k1) = dq_ci_mul(2:i1,2:j1,1:k1)+ dq_ci_spl(2:i1,2:j1,1:k1)
     ! dq_hghr_rime(2:i1,2:j1,1:k1) = dq_hghr_rime(2:i1,2:j1,1:k1)    - dq_ci_spl(2:i1,2:j1,1:k1)
     q_hgp    (2:i1,2:j1,1:k1) = q_hgp(2:i1,2:j1,1:k1)    - dq_ci_spl(2:i1,2:j1,1:k1)  ! effect on snow


    ! i+r-> g   -- does it make sense to include it?
    dq_hg_temp(2:i1,2:j1,1:k1) = -dq_ci_col_ri(2:i1,2:j1,1:k1)-dq_hr_col_ri(2:i1,2:j1,1:k1)
    call hallet_mossop3 (dq_hg_temp,q_hg,dq_ci_spl,dn_ci_spl)
     dn_ci_mul(2:i1,2:j1,1:k1) = dn_ci_mul(2:i1,2:j1,1:k1)+ dn_ci_spl(2:i1,2:j1,1:k1)
     dq_ci_mul(2:i1,2:j1,1:k1) = dq_ci_mul(2:i1,2:j1,1:k1)+ dq_ci_spl(2:i1,2:j1,1:k1)
     ! dq_hghr_rime(2:i1,2:j1,1:k1) = dq_hghr_rime(2:i1,2:j1,1:k1)    - dq_ci_spl(2:i1,2:j1,1:k1)
     q_hgp    (2:i1,2:j1,1:k1) = q_hgp(2:i1,2:j1,1:k1)    - dq_ci_spl(2:i1,2:j1,1:k1)  ! effect on snow

    ! ------------------------------------
    !   update
    ! ------------------------------------
    ! add updates to dq_ci, dn_ci
    n_cip(2:i1,2:j1,1:k1)=n_cip(2:i1,2:j1,1:k1)+dn_ci_mul(2:i1,2:j1,1:k1)
    q_cip(2:i1,2:j1,1:k1)=q_cip(2:i1,2:j1,1:k1)+dq_ci_mul(2:i1,2:j1,1:k1)


    ! clean-up
    deallocate(dq_hg_temp,dq_ci_spl, dn_ci_spl)

   end subroutine ice_multi3




! ***************************************************************
! Ice multiplication of Hallet and Mossop (1974)
!
! - written as described in Seifert (2002)
! - implementation similar to the one in ICON model
!
! - returns:
!    dq_i_hm  - change in moass content during the process
!    dn_i_hm  - number of newly produced ice particles
!
!
! ***************************************************************
   subroutine hallet_mossop3(dq_rime,q_e,dq_i_hm,dn_i_hm)

   use modglobal, only : ih,i1,jh,j1,k1,rv,rd, rlv,cp,pi,rhow
   use modfields, only : exnf,qt0,tmp0
   implicit none

   ! inputs
   real, intent(in)  :: dq_rime  (:,:,:)  ! riming rate
   real, intent(in)  :: q_e      (:,:,:)  ! amount of that ice phase
   ! outputs
   real, intent(out) :: dq_i_hm   (:,:,:) ! tendency in q_i by H-M process
   real, intent(out) :: dn_i_hm   (:,:,:) ! tendency in n_i by H-M process

   ! local variables
   integer :: i,j,k
   real    :: c_spl,c_1_hm,c_2_hm                        &  ! constants in calculation
             ,mult_1, mult_2, mint_1, mint_2             &  ! calculation variables
             ,dn_try, dq_try, rem_cf                           ! trial variables

   ! allocated variables
   ! real, allocatable, dimension(:,:,:) :: f_spl
   !
   ! allocating
   ! allocate( f_spl (2-ih:i1+ih,2-jh:j1+jh,k1) )

   ! f_spl = 0.0

   ! setting constants
   c_spl   = c_spl_hm74
   c_1_hm  = 1.0/(tmp_opt_hm74-tmp_min_hm74)
   c_2_hm  = 1.0/(tmp_opt_hm74-tmp_max_hm74)

   ! setting coefficient for reminder
   rem_cf  = (1.0-rem_q_e_hm)/delt


   ! calculation
    do k=1,k1
    do j=2,j1
    do i=2,i1
     if ((dq_rime(i,j,k).gt.0).and.(tmp0(i,j,k).lt.T_3)) then ! only if riming going on temperature below 0
        ! f_spl calculation following ICON
        mult_1 = c_1_hm *(tmp0(i,j,k)-tmp_min_hm74)   ! positive in the target interval
        mult_2 = c_2_hm *(tmp0(i,j,k)-tmp_max_hm74)   ! positive in the target interval
        ! now for intervals
        mint_1 = max(0.0,min(1.0,mult_1))             !  0 for T<T_min, 1 for T>T_opt
        mint_2 = max(0.0,min(1.0,mult_2))             !  0 for T>T_max, 1 for T<T_opt
        ! calculating prediction for the process
        dn_try = c_spl*mint_1*mint_2*dq_rime(i,j,k)
        dq_try = x_ci_spl*dn_try
        ! correcting
        dq_try = min(dq_try,rem_cf*q_e(i,j,k)+dq_rime(i,j,k))    !  limit splintering
        dq_try = max(dq_try,0.0)
        !
        ! prepare updates
        dq_i_hm(i,j,k) = dq_try
        dn_i_hm(i,j,k) = dq_try/x_ci_spl
        !
     else
        dq_i_hm(i,j,k) = 0.0
        dn_i_hm(i,j,k) = 0.0
     endif
    enddo
    enddo
    enddo

   ! deallocate(f_spl)

   end subroutine hallet_mossop3



! ***************************************************************
! Melting of ice particles
!
! - this is a inner subroutine called from a wrapper
! - written as described in Seifert&Beheng (2004)
!   - based on Pruppacher and Klett (1997)
! - implementation similar to the one in ICON model
!
! - returns:
!    dq_me  - mass content melting tendency
!    dn_me  - number content melting tendency
!    dq_ev  - mass content evaporation tendency
!    dn_ev  - number content evaporation tendency
!
! ***************************************************************
  subroutine sb_evmelt3(avent0,avent1,bvent0,bvent1,x_bmin,n_e   &
           ,q_e,x_e,D_e,v_e,dq_me,dn_me,dq_ev,dn_ev)
   ! vapour_deposit(aic,bic,alpaic,betaic,xic,dqdep)

    use modglobal, only : ih,i1,jh,j1,k1,rv,rd, rlv,cp,pi
    use modfields, only : exnf,qt0,svm,tmp0,esl,qvsl,rhof,exnf,presf
    implicit none

    ! inputs -------------------------------------------------------------
    real, intent(in)  :: n_e     (:,:,:)  ! number density of ice particles
    real, intent(in)  :: q_e     (:,:,:)  ! mass density of ice particles
    real, intent(in)  :: x_e     (:,:,:)  ! mean size of the ice particles
    real, intent(in)  :: D_e     (:,:,:)  ! mean diameter of particles
    real, intent(in)  :: v_e     (:,:,:)  ! mean terminal velocity of particles
    real, intent(in)  :: avent0           ! ventilation coefficient
    real, intent(in)  :: avent1           ! ventilation coefficient
    real, intent(in)  :: bvent0           ! ventilation coefficient
    real, intent(in)  :: bvent1           ! ventilation coefficient
    real, intent(in)  :: x_bmin           ! minimal size of the hydrometeor
    ! real, intent(in)  :: k_melt         ! depositional growth constant

    ! outputs  -----------------------------------------------------------
    real, intent(out) :: dq_me   (:,:,:) ! melting tendency in q_e
    real, intent(out) :: dn_me   (:,:,:) ! melting tendency in n_e
    real, intent(out) :: dq_ev   (:,:,:) ! evaporation tendency in q_e
    real, intent(out) :: dn_ev   (:,:,:) ! evaporation tendency in n_e

    ! loacal variables -------------------------------------------------

    integer :: i,j,k
    real                 :: k_melt, k_ev, ktdtodv,dvleorv, g_ev_const &
                            ,eslt3t3
    real, allocatable    :: S(:,:,:),nrex(:,:,:),f0(:,:,:),f1(:,:,:) &
                           ,g_me(:,:,:),g_ev(:,:,:),x_er(:,:,:)    &
                           ,me_q(:,:,:), me_n(:,:,:)
    logical ,allocatable :: melt_mask(:,:,:)

    ! preparing arrays and values  ----------------------------------


    !  allocating mask
    allocate(  melt_mask (2-ih:i1+ih,2-jh:j1+jh,k1))    ! mask
    !
    ! allocate fields and fill
    allocate(  S       (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! subsaturation
              ,g_me    (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! thermodynamic term for melting
              ,g_ev   (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! thermodynamic term for evaporation
              ,x_er    (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! mean size of particles
              ,f0      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! ventilation factor
              ,f1      (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! ventilation factor
              ,nrex    (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! reynolds number
              ,me_q    (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! basic melting rate in q
              ,me_n    (2-ih:i1+ih,2-jh:j1+jh,k1)  & ! basic melting rate in n
             )
    ! - filling
    S           =  0.0
    g_me        =  0.0
    g_ev       =  0.0
    x_er        =  0.0
    f0          =  0.0
    f1          =  0.0
    nrex        =  0.0
    me_q        =  0.0
    me_n        =  0.0
    !
    melt_mask   = .false.

    ! set constants
    k_melt      = 2*pi/rlme
    k_ev       = 2*pi
    ktdtodv     = Kt**2/(cp*rho0s*Dv)   ! |<  = K_T * D_T / D_v  and D_T = K_T/(c_p \rho_s)
    dvleorv     = Dv*rlvi/rv
    ! constant evaporation parameter for melting particles
    g_ev_const = (rv*T_3)/(Dv*eslt3)+rlv/(Kt*T_3)*(rlv/(rv*T_3) -1.)
    ! based on: G   (i,j,k) = (rv * tmp0(i,j,k)) / (Dv*esl(i,j,k)) + rlv/(Kt*tmp0(i,j,k))*(rlv/(rv*tmp0(i,j,k)) -1.)
    g_ev_const = 1.0 / g_ev_const
    ! additional calculations
    eslt3t3 = eslt3/T_3

    ! set tendencies to 0
    dq_me  = 0.0
    dn_me  = 0.0
    dq_ev = 0.0
    dn_ev = 0.0

    ! depositional growth constant
    ! k_depos = 4*pi/cip   ! for spherical particles


    ! -- preaparing the mask  ---------------------------------------
    do k=1,k1
    do j=2,j1
    do i=2,i1
     ! if temperature above 0 and enough ice particles
     if (( tmp0(i,j,k).gt.T_3 ).and.(q_e(i,j,k).gt.qicemin)) then
       melt_mask(i,j,k) = .true.
     endif
    enddo
    enddo
    enddo

    ! preparing calculation
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! calculating for all cells with the value
      if (melt_mask(i,j,k)) then
       ! calculation of the subsaturation
        S(i,j,k) = min(0.0,((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)- 1.0))
       ! calculating the thermodynamic term for evaporation
        ! g_ev(i,j,k) =(rv*tmp0(i,j,k))/(Dv*esl(i,j,k))+rlv/(Kt*tmp0(i,j,k))*(rlv/(rv*tmp0(i,j,k)) -1.)
        ! g_ev(i,j,k) =(rv*T_3)/(Dv*eslt3)+rlv/(Kt*T_3)*(rlv/(rv*T_3) -1.)
        g_ev(i,j,k) = g_ev_const ! 1.0/g_ev(i,j,k)
       ! calculation of the thermodynamic term for melting
        g_me(i,j,k)= - k_melt*(ktdtodv*(tmp0(i,j,k)-T_3)+          &
                    dvleorv*(esl(i,j,k)/tmp0(i,j,k)-eslt3t3))
       ! calculating real mean particle mass
        x_er(i,j,k) = q_e(i,j,k)/(n_e(i,j,k)+eps0)
       ! calculating N_re Reynolds number
        nrex(i,j,k)= D_e(i,j,k)*v_e(i,j,k)/nu_a
       ! calculating from prepared ventilation coefficients
        f0(i,j,k)  = avent0+bvent0*Sc_num**(1.0/3.0)*nrex(i,j,k)**0.5
        f1(i,j,k)  = avent1+bvent1*Sc_num**(1.0/3.0)*nrex(i,j,k)**0.5
       ! preapring updates for evaporation
        dn_ev(i,j,k) = k_ev*g_ev(i,j,k)*S(i,j,k)*n_e(i,j,k)*  &
                        D_e(i,j,k)*f0(i,j,k)/max(x_bmin, x_er(i,j,k))
        dq_ev(i,j,k) = k_ev*g_ev(i,j,k)*S(i,j,k)*n_e(i,j,k)*  &
                        D_e(i,j,k)*f1(i,j,k)
       ! preapring updates for melting
        me_n(i,j,k) = g_me(i,j,k)*n_e(i,j,k)*D_e(i,j,k)*         &
                    f0(i,j,k)/max(x_bmin, x_er(i,j,k))
        me_q(i,j,k) = g_me(i,j,k)*n_e(i,j,k)*D_e(i,j,k)*f1(i,j,k)
      ! and limiting so not removing more than available
        ! basic correction of melting rate
        me_q(i,j,k) = min(0.0,max(me_q(i,j,k),-q_e(i,j,k)/delt))
        ! basic correction for number tendency in melting
        me_n(i,j,k) = min(0.0,max(me_n(i,j,k),-n_e(i,j,k)/delt))
        ! and prevent melting of all particles while leaving mass ?
        !out : me_n(i,j,k) = max(me_q(i,j,k)/x_er(i,j,k),me_n(i,j,k))
        ! basic correction for mass tendency in evaporation
        dq_ev(i,j,k) = min(0.0,max(dq_ev(i,j,k),-q_e(i,j,k)/delt))
        ! basic correction for number tendency in evaporation
        dn_ev(i,j,k) = min(0.0,max(dn_ev(i,j,k),-n_e(i,j,k)/delt))
        ! prevent evaporation of all particles while leaving mass ?
        ! out : dn_ev(i,j,k) = max(dq_ev(i,j,k)/x_er(i,j,k),dn_ev(i,j,k))
        ! now what melts without evaporating
        dq_me(i,j,k) = me_q(i,j,k) ! min(0.0, me_q(i,j,k) - dq_ev(i,j,k))
        dn_me(i,j,k) = me_n(i,j,k) ! min(0.0, me_n(i,j,k) - dn_ev(i,j,k))
        ! basic correction for mass tendency in melting
        !d dq_me(i,j,k) = min(0.0,max(dq_me(i,j,k),-dq_ev(i,j,k)-q_e(i,j,k)/delt))
        ! basic correction for number tendency in melting
        !d dn_me(i,j,k) = min(0.0,max(dn_me(i,j,k),-dn_ev(i,j,k)-n_e(i,j,k)/delt))
        ! and prevent melting of all particles while leaving mass ?
        ! out: dn_me(i,j,k) = max(dq_me(i,j,k)/x_er(i,j,k),dn_me(i,j,k))  ! i.e. prevent increasing the mean size of particles during melting
        !
        ! #ICON - based on assumption that x_s is app. constant during melting
        ! dn_em(i,j,k) = (q_e(i,j,k)/delt+ dq_em(i,j,k))/x_e(i,j,k)-n_e(i,j,k)/delt
        ! dn_em(i,j,k) = min(0.0,max(-n_e(i,j,k)/delt,dn_em(i,j,k)))
      endif
    enddo
    enddo
    enddo


    ! --- end of the inner part ---------------------------

   ! deallocating fields
    deallocate (S,g_me,g_ev,x_er,f0, f1,nrex,me_q,me_n)
    deallocate (melt_mask)

   end subroutine sb_evmelt3


! ***************************************************************
! Enhanced Melting of ice particles
!
!  - not used, instead calculated separately in collection subroutine
!
! - this is a inner subroutine called from a wrapper
! - written as described in Seifert&Beheng (2004)
!   - based on Rutledge & Hobbs (1984)
! - implementation similar to the one in ICON model


! ***************************************************************
!    recovery of ccn
!
!   - to be later replaced based on advance literature
! ***************************************************************
    subroutine recover_cc

    use modglobal, only : ih,i1,j1,k1  ! ,jh,rv,rd, rlv,cp,pi,rhow
    ! use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,qvsl,qvsi,rhof,exnf,presf
    implicit none
    integer :: i,j,k

    ! allocations

    !   - to be later replaced based on advance literature
    !
    !   - so far just and easy recovery
    !     of ccn based on number of water particles that evaporated, sublimated
    !     or got removed with remaining positive n_
    if(.not.(l_c_ccn)) then ! ie. no change if l_c_ccn
     do k=1,k1
     do j=2,j1
     do i=2,i1
      ! decrease in total amount of potential CCN
      ! decrease by cloud processes:
      !   self-collection, sedim, precip, riming
      n_ccp(i,j,k)=n_ccp(i,j,k)+dn_cl_sc(i,j,k)+dn_cl_se(i,j,k)+   &
        dn_cl_au(i,j,k)+dn_cl_ac(i,j,k)+                           &
        dn_cl_hom(i,j,k)+dn_cl_het(i,j,k)+                         &
        dn_cl_rime_ci(i,j,k)+dn_cl_rime_hs(i,j,k)+dn_cl_rime_hg(i,j,k)
      ! recovery of potential CCN
      n_ccp (i,j,k) = n_ccp (i,j,k)+c_rec_cc*ret_cc(i,j,k)
     enddo
     enddo
     enddo
    endif

    !deallocations

    end subroutine recover_cc


! ****************************************
!  Debugging subroutines
!
! ****************************************

   ! checking for strange of incorrect resuls

    ! checking for strange of incorrect resuls
   subroutine check_sizes(flag_dbg,flag,proc_name)
     use modglobal, only : ih,i1,jh,j1,k1
     use modfields, only : sv0,svm,ql0,qvsl,qvsi,qt0, thl0, tmp0, thlp, qtp, rhof
   implicit none
    logical, intent (inout) ::flag_dbg,flag
    character (len=7), intent (in) :: proc_name
    real :: lim_x_c, lim_x_h , ssat, ssat_ice, uncond, nrtest, facmax
    real ,allocatable ,dimension(:,:,:) :: xt_cl, xt_ci, xt_hr, xt_hs, xt_hg  ! size testing
    logical ,allocatable ,dimension(:,:,:) :: ma_cl, ma_ci, ma_hr, ma_hs, ma_hg ! size testing
    integer:: i,j,k

   ! inner setting
    ! -- limit for the change of the size
    lim_x_c = 2.0e-7
    lim_x_h = 1.0e-5

    facmax = 100.0 ! how many times exceed standard max size


    if(flag_dbg) then
    ! allocate
        ! allocate fields and fill
     allocate(  xt_cl   (2-ih:i1+ih,2-jh:j1+jh,k1)     &
               ,xt_ci   (2-ih:i1+ih,2-jh:j1+jh,k1)     &
               ,xt_hr   (2-ih:i1+ih,2-jh:j1+jh,k1)     &
               ,xt_hs   (2-ih:i1+ih,2-jh:j1+jh,k1)     &
               ,xt_hg   (2-ih:i1+ih,2-jh:j1+jh,k1)     &
             )

      allocate(  ma_cl(2-ih:i1+ih,2-jh:j1+jh,k1)       &
                ,ma_ci(2-ih:i1+ih,2-jh:j1+jh,k1)       &
                ,ma_hr(2-ih:i1+ih,2-jh:j1+jh,k1)       &
                ,ma_hs(2-ih:i1+ih,2-jh:j1+jh,k1)       &
                ,ma_hg(2-ih:i1+ih,2-jh:j1+jh,k1)       &
              )

     xt_cl = 0.0
     xt_ci = 0.0
     xt_hr = 0.0
     xt_hs = 0.0
     xt_hg = 0.0

     ma_cl = .false.
     ma_ci = .false.
     ma_hr = .false.
     ma_hs = .false.
     ma_hg = .false.

    ! calculate sizes after update
    do k=1,k1
    do j=2,j1
    do i=2,i1
      !
      ! cl
      nrtest=svm(i,j,k,in_cl)+delt*n_clp(i,j,k) ! considering n_clp only
      if (nrtest.le.0.0) then
        xt_cl(i,j,k) =0.0
      else
        !  calculating new droplet size
        xt_cl(i,j,k) =  (svm(i,j,k,iq_cl) +    &
                delt*q_clp(i,j,k)) / (svm(i,j,k,in_cl)+delt*n_clp(i,j,k)+eps0)
                ! o: rhof(k) * (svm(i,j,k,iqsv) +  delt*q_clp(i,j,k)) / (svm(i,j,k,insv)+delt*n_clp(i,j,k))
      endif
      !
      ! ci
      nrtest=svm(i,j,k,in_cl)+delt*n_cip(i,j,k) ! considering n_clp only
      if (nrtest.le.0.0) then
        xt_ci(i,j,k) =0.0
      else
        !  calculating new droplet size
        xt_ci(i,j,k) =  (svm(i,j,k,iq_ci) +    &
                delt*q_cip(i,j,k)) / (svm(i,j,k,in_ci)+delt*n_cip(i,j,k)+eps0)
      endif
      !
      ! hr
      nrtest=svm(i,j,k,in_hr)+delt*n_hrp(i,j,k) ! considering n_clp only
      if (nrtest.le.0.0) then
        xt_hr(i,j,k) =0.0
      else
        !  calculating new droplet size
        xt_hr(i,j,k) =  (svm(i,j,k,iq_hr) +    &
                delt*q_hrp(i,j,k)) / (svm(i,j,k,in_hr)+delt*n_hrp(i,j,k)+eps0)
      endif
      !
      ! hs
      nrtest=svm(i,j,k,in_hs)+delt*n_hsp(i,j,k) ! considering n_clp only
      if (nrtest.le.0.0) then
        xt_hs(i,j,k) =0.0
      else
        !  calculating new droplet size
        xt_hs(i,j,k) =  (svm(i,j,k,iq_hs) +    &
                delt*q_hsp(i,j,k)) / (svm(i,j,k,in_hs)+delt*n_hsp(i,j,k)+eps0)
      endif
      !
      ! hg
      nrtest=svm(i,j,k,in_hg)+delt*n_hgp(i,j,k) ! considering n_clp only
      if (nrtest.le.0.0) then
        xt_hg(i,j,k) =0.0
      else
        !  calculating new droplet size
        xt_hg(i,j,k) =  (svm(i,j,k,iq_hg) +    &
                delt*q_hgp(i,j,k)) / (svm(i,j,k,in_hg)+delt*n_hgp(i,j,k)+eps0)
      endif
      !
    enddo
    enddo
    enddo

    ! masks
    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! clouds
      if ((xt_cl(i,j,k).gt.(facmax*x_cl_bmax)).or.((xt_cl(i,j,k)-x_cl(i,j,k)).gt. lim_x_c)) then
         ma_cl(i,j,k) = .true.
      endif
      if ((xt_cl(i,j,k).gt.(facmax*x_ci_bmax)).or.((xt_ci(i,j,k)-x_ci(i,j,k)).gt. lim_x_c)) then
         ma_ci(i,j,k) = .true.
      endif
      !
      ! and hydrometeors
      if ((xt_hr(i,j,k).gt.(facmax*x_hr_bmax)).or.((xt_hr(i,j,k)-x_hr(i,j,k)).gt. lim_x_h)) then
         ma_hr(i,j,k) = .true.
      endif
      if ((xt_hs(i,j,k).gt.(facmax*x_hs_bmax)).or.((xt_hs(i,j,k)-x_hs(i,j,k)).gt. lim_x_h)) then
         ma_hs(i,j,k) = .true.
      endif
      if ((xt_hg(i,j,k).gt.(facmax*x_hg_bmax)).or.((xt_hg(i,j,k)-x_hg(i,j,k)).gt. lim_x_h)) then
         ma_hg(i,j,k) = .true.
      endif
      !
    enddo
    enddo
    enddo



     if (any(ma_cl.or.ma_ci.or.ma_hr.or.ma_hs.or.ma_hg)) then ! problem detected
       write(6,*) 'WARNING: issue_sizes in: ',proc_name
       write(6,*) ' issue with size in no. of (cl, ci, hr, hs, hg): ',count(ma_cl), count(ma_ci), count(ma_hr),count(ma_hs), count(ma_hg)
      if (flag) then
       write(6,*)  ' most of them already listed '
      else ! flag
       flag = .true.
       do k=1,k1
       do j=2,j1
       do i=2,i1
        if (ma_cl(i,j,k).or.ma_ci(i,j,k).or.ma_hr(i,j,k).or.ma_hs(i,j,k).or.ma_hg(i,j,k)) then
          ! calculate
          uncond=ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k
          write(6,*) '   thlpmcr=', thlpmcr(i,j,k),'qtpmcr=', qtpmcr(i,j,k)
          write(6,*) '   thlp=',thlp(i,j,k),'qtp=', qtp(i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),'thl = ', thl0(i,j,k), 'tmp0 = ', tmp0(i,j,k)
          write(6,*) '   rhof=',rhof(k),' n_ccn = ', n_cc(i,j,k)
          write(6,*) '   x_cl =', x_cl(i,j,k), ' xt_cl =', xt_cl(i,j,k)
          write(6,*) '   x_ci =', x_ci(i,j,k), ' xt_ci =', xt_ci(i,j,k)
          write(6,*) '   x_hr =', x_hr(i,j,k), ' xt_hr =', xt_hr(i,j,k)
          write(6,*) '   x_hs =', x_hs(i,j,k), ' xt_hs =', xt_hs(i,j,k)
          write(6,*) '   x_hg =', x_hg(i,j,k),' xt_hg =', xt_hg (i,j,k)
          write(6,*) '   n_cl =', n_cl(i,j,k),' q_cl =', q_cl (i,j,k)
          write(6,*) '   n_ci =', n_ci(i,j,k),' q_ci =', q_ci (i,j,k)
          write(6,*) '   n_hr =', n_hr(i,j,k),' q_hr =', q_hr (i,j,k)
          write(6,*) '   n_hs =', n_hs(i,j,k),' q_hs =', q_hs (i,j,k)
          write(6,*) '   n_hg =', n_hg(i,j,k),' q_hg =', q_hg (i,j,k)
          ! updates  and past values
          write(6,*) '   q_clp =',q_clp(i,j,k),' n_clp =', n_clp(i,j,k),' n_cl_m =', svm(i,j,k,in_cl),' q_cl_m =',svm(i,j,k,iq_cl)
          write(6,*) '   q_cip =',q_cip(i,j,k),' n_cip =', n_cip(i,j,k),' n_ci_m =', svm(i,j,k,in_ci),' q_ci_m =',svm(i,j,k,iq_ci)
          write(6,*) '   q_hrp =',q_hrp(i,j,k),' n_hrp =', n_hrp(i,j,k),' n_hr_m =', svm(i,j,k,in_hr),' q_hr_m =',svm(i,j,k,iq_hr)
          write(6,*) '   q_hsp =',q_hsp(i,j,k),' n_hsp =', n_hsp(i,j,k),' n_hs_m =', svm(i,j,k,in_hs),' q_hs_m =',svm(i,j,k,iq_hs)
          write(6,*) '   q_hgp =',q_hgp(i,j,k),' n_hgp =', n_hgp(i,j,k),' n_hg_m =', svm(i,j,k,in_hg),' q_hg_m =',svm(i,j,k,iq_hg)
          ! and processes
          write(6,*) '   dn_cl_nu=',dn_cl_nu(i,j,k),' dn_ci_inu=',dn_ci_inu(i,j,k)
          write(6,*) '   dq_hr_au=',dq_hr_au(i,j,k), 'dq_hr_ac=', dq_hr_ac(i,j,k),'dq_hr_ev=' ,dq_hr_ev(i,j,k)
          write(6,*) '   dn_cl|au=',(2.0/x_s)*dq_hr_au(i,j,k), 'dn_cl|ac=', dq_hr_ac(i,j,k)/(x_cl(i,j,k)+eps0),'dn_hr_ev=' ,dn_hr_ev(i,j,k)
          write(6,*) '   dq_ci_dep=',dq_ci_dep(i,j,k),'dq_hs_dep=',dq_hs_dep(i,j,k),'dq_hg_dep=',dq_hg_dep    (i,j,k)
          write(6,*) '   dq_ci_rime=',dq_ci_rime (i,j,k),'dq_hs_rime =',dq_hs_rime (i,j,k),'dq_hg_rime =',dq_hg_rime   (i,j,k)
          write(6,*) '   dn_cl_rime_ci=',dn_cl_rime_ci(i,j,k),'dn_cl_rime_hs=',dn_cl_rime_hs(i,j,k),'dn_cl_rime_h=',dn_cl_rime_hg(i,j,k)
          write(6,*) '   dq_hshr_rime=',dq_hshr_rime (i,j,k),'dq_hghr_rime=',dq_hghr_rime (i,j,k)
          write(6,*) '   dn_hr_rime_hs=',dn_hr_rime_hs(i,j,k),'dn_hr_rime_hg=',dn_hr_rime_hg(i,j,k)
          write(6,*) '   dq_cl_het=' ,dq_cl_het    (i,j,k),'dq_hr_het=',dq_hr_het    (i,j,k), 'dq_cl_hom=',dq_cl_hom    (i,j,k)
          write(6,*) '   dn_cl_het=',dn_cl_het    (i,j,k),'dn_hr_het=',dn_hr_het    (i,j,k), 'dn_cl_hom=',dn_cl_hom    (i,j,k)
          write(6,*) '   dq_ci_col_iis=',dq_ci_col_iis    (i,j,k),'dq_hsci_col=',dq_hsci_col  (i,j,k)
          write(6,*) '   dn_ci_col_iis=',dn_ci_col_iis    (i,j,k), 'dn_hs_col_sss=' ,dn_hs_col_sss  (i,j,k), 'dn_ci_col_iis='  ,dn_ci_col_iis (i,j,k)
          write(6,*) '   dq_ci_cv=',dq_ci_cv (i,j,k), 'dq_hs_cv=' ,dq_hs_cv     (i,j,k)
          write(6,*) '   dn_ci_cv=',dn_ci_cv     (i,j,k), 'dn_hs_cv=' ,dn_hs_cv     (i,j,k)
          write(6,*) '   dq_cl_se=',dq_cl_se     (i,j,k),'dq_ci_se=',dq_ci_se     (i,j,k)
          write(6,*) '   dn_cl_sc=',dn_cl_sc     (i,j,k),'dn_cl_se=',dn_cl_se     (i,j,k),'dn_ci_se=',dn_ci_se (i,j,k)
          write(6,*) '   dq_hr_se=',dq_hr_se     (i,j,k), 'dq_hs_se=',dq_hs_se     (i,j,k),'dq_hg_se=',dq_hg_se     (i,j,k)
          write(6,*) '   dn_hr_se=',dn_hr_se     (i,j,k),'dn_hs_se=',dn_hs_se     (i,j,k) ,'dn_hg_se=',dn_hg_se     (i,j,k)
          write(6,*) '   precep_hs=',precep_hs    (i,j,k),'precep_hg=',precep_hg    (i,j,k)
          write(6,*) '   dq_cl_sa=',dq_cl_sa     (i,j,k) ,'ret_cc=',ret_cc       (i,j,k)
        endif
       enddo
       enddo
       enddo
      endif ! flag
    ! else ! no problem detected
    ! ! reset flag
    !  ! flag = .false.
    endif ! problem detected

    deallocate( xt_cl,xt_ci,xt_hr,xt_hs,xt_hg )

    deallocate(  ma_cl,ma_ci,ma_hr,ma_hs,ma_hg)
   endif !flag_dbg

   end subroutine check_sizes



   subroutine check_update(flag_dbg,flag,proc_name)
     use modglobal, only : ih,i1,jh,j1,k1
     use modfields, only : sv0,svm,ql0,qvsl,qvsi,qt0, thl0
   implicit none
    logical, intent (inout) ::flag_dbg,flag
    character (len=7), intent (in) :: proc_name
    real :: lim_thlp, lim_qtp, ssat, ssat_ice, uncond
    integer:: i,j,k

   ! inner setting
    lim_thlp = 4.0
    lim_qtp = 0.0005

   if(flag_dbg) then

   ! checking for large updates
   if(any(thlpmcr(2:i1,2:j1,1:k1).gt.lim_thlp).or.any(thlpmcr(2:i1,2:j1,1:k1).lt.(-lim_thlp))) then
    write(6,*) 'WARNING: thl problem in: ', proc_name
    if (flag) then
      write(6,*) 'problems in: ',count(thlpmcr(2:i1,2:j1,1:k1).gt.lim_thlp),'and',count(thlpmcr(2:i1,2:j1,1:k1).lt.(-lim_thlp))
    else ! flag
     flag = .true.
    do k=1,k1
    do j=2,j1
    do i=2,i1
        if (thlpmcr(i,j,k).gt.lim_thlp) then
          ! calculate
          uncond=ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k,'  thlpmcr = ', thlpmcr(i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),' thl = ', thl0(i,j,k), ' n_ccn = ', n_cc(i,j,k)
          write(6,*) '  x_cl =', x_cl(i,j,k), ' n_cl =', n_cl(i,j,k),' q_cl =', q_cl (i,j,k)
          write(6,*) '  x_ci =', x_ci(i,j,k), ' n_ci =', n_ci(i,j,k),' q_ci =', q_ci (i,j,k)
          write(6,*) '  x_hr =', x_hr(i,j,k), ' n_hr =', n_hr(i,j,k),' q_hr =', q_hr (i,j,k)
          write(6,*) '  x_hs =', x_hs(i,j,k), ' n_hs =', n_hs(i,j,k),' q_hs =', q_hs (i,j,k)
          write(6,*) '  x_hg =', x_hg(i,j,k), ' n_hg =', n_hg(i,j,k),' q_hg =', q_hg (i,j,k)
          write(6,*) '  q_clp =',q_clp(i,j,k),' n_clp =', n_clp(i,j,k),' n_cl_m =', svm(i,j,k,in_cl),' q_cl_m =',svm(i,j,k,iq_cl)
          write(6,*) '  q_cip =',q_cip(i,j,k),' n_cip =', n_cip(i,j,k),' n_ci_m =', svm(i,j,k,in_ci),' q_ci_m =',svm(i,j,k,iq_ci)
          write(6,*) '  q_hrp =',q_hrp(i,j,k),' n_hrp =', n_hrp(i,j,k),' n_hr_m =', svm(i,j,k,in_hr),' q_hr_m =',svm(i,j,k,iq_hr)
          write(6,*) '  q_hsp =',q_hsp(i,j,k),' n_hsp =', n_hsp(i,j,k),' n_hs_m =', svm(i,j,k,in_hs),' q_hs_m =',svm(i,j,k,iq_hs)
          write(6,*) '  q_hgp =',q_hgp(i,j,k),' n_hgp =', n_hgp(i,j,k),' n_hg_m =', svm(i,j,k,in_hg),' q_hg_m =',svm(i,j,k,iq_hg)
        endif
        if (thlpmcr(i,j,k).lt.(-lim_thlp)) then
          ! calculate
          uncond   = ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k,'  thlpmcr = ', thlpmcr(i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),' thl = ', thl0(i,j,k), ' n_ccn = ', n_cc(i,j,k)
          write(6,*) '  x_cl =', x_cl(i,j,k), ' n_cl =', n_cl(i,j,k),' q_cl =', q_cl (i,j,k)
          write(6,*) '  x_ci =', x_ci(i,j,k), ' n_ci =', n_ci(i,j,k),' q_ci =', q_ci (i,j,k)
          write(6,*) '  x_hr =', x_hr(i,j,k), ' n_hr =', n_hr(i,j,k),' q_hr =', q_hr (i,j,k)
          write(6,*) '  x_hs =', x_hs(i,j,k), ' n_hs =', n_hs(i,j,k),' q_hs =', q_hs (i,j,k)
          write(6,*) '  x_hg =', x_hg(i,j,k), ' n_hg =', n_hg(i,j,k),' q_hg =', q_hg (i,j,k)
          write(6,*) '  q_clp =',q_clp(i,j,k),' n_clp =', n_clp(i,j,k),' n_cl_m =', svm(i,j,k,in_cl),' q_cl_m =',svm(i,j,k,iq_cl)
          write(6,*) '  q_cip =',q_cip(i,j,k),' n_cip =', n_cip(i,j,k),' n_ci_m =', svm(i,j,k,in_ci),' q_ci_m =',svm(i,j,k,iq_ci)
          write(6,*) '  q_hrp =',q_hrp(i,j,k),' n_hrp =', n_hrp(i,j,k),' n_hr_m =', svm(i,j,k,in_hr),' q_hr_m =',svm(i,j,k,iq_hr)
          write(6,*) '  q_hsp =',q_hsp(i,j,k),' n_hsp =', n_hsp(i,j,k),' n_hs_m =', svm(i,j,k,in_hs),' q_hs_m =',svm(i,j,k,iq_hs)
          write(6,*) '  q_hgp =',q_hgp(i,j,k),' n_hgp =', n_hgp(i,j,k),' n_hg_m =', svm(i,j,k,in_hg),' q_hg_m =',svm(i,j,k,iq_hg)
        endif
    enddo
    enddo
    enddo
    endif ! flag
   endif

   if(any(qtpmcr(2:i1,2:j1,1:k1).lt.(-lim_qtp)).or.any(qtpmcr(2:i1,2:j1,1:k1).gt.lim_qtp)) then
    write(6,*) 'WARNING: qtp problem in: ', proc_name
    if (flag) then
      write(6,*) 'problems in: ',count(qtpmcr(2:i1,2:j1,1:k1).lt.(-lim_qtp)),'and',count(qtpmcr(2:i1,2:j1,1:k1).gt.lim_qtp)
    else ! flag
     flag = .true.
    do k=1,k1
    do j=2,j1
    do i=2,i1
        if (qtpmcr(i,j,k).gt.lim_qtp) then
          ! calculate
          uncond   = ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond   = ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k,'  qtpmcr =', qtpmcr (i,j,k)
          write(6,*) '   uncondensed =', uncond,'  ssat =',ssat,'% ','  ssat_i =',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),' thl = ', thl0(i,j,k), ' n_ccn =', n_cc(i,j,k)
          write(6,*) '  x_cl =', x_cl(i,j,k), ' n_cl =', n_cl(i,j,k),' q_cl =', q_cl (i,j,k)
          write(6,*) '  x_ci =', x_ci(i,j,k), ' n_ci =', n_ci(i,j,k),' q_ci =', q_ci (i,j,k)
          write(6,*) '  x_hr =', x_hr(i,j,k), ' n_hr =', n_hr(i,j,k),' q_hr =', q_hr (i,j,k)
          write(6,*) '  x_hs =', x_hs(i,j,k), ' n_hs =', n_hs(i,j,k),' q_hs =', q_hs (i,j,k)
          write(6,*) '  x_hg =', x_hg(i,j,k), ' n_hg =', n_hg(i,j,k),' q_hg =', q_hg (i,j,k)
          write(6,*) '  q_clp =',q_clp(i,j,k),' n_clp =', n_clp(i,j,k),' n_cl_m =', svm(i,j,k,in_cl),' q_cl_m =',svm(i,j,k,iq_cl)
          write(6,*) '  q_cip =',q_cip(i,j,k),' n_cip =', n_cip(i,j,k),' n_ci_m =', svm(i,j,k,in_ci),' q_ci_m =',svm(i,j,k,iq_ci)
          write(6,*) '  q_hrp =',q_hrp(i,j,k),' n_hrp =', n_hrp(i,j,k),' n_hr_m =', svm(i,j,k,in_hr),' q_hr_m =',svm(i,j,k,iq_hr)
          write(6,*) '  q_hsp =',q_hsp(i,j,k),' n_hsp =', n_hsp(i,j,k),' n_hs_m =', svm(i,j,k,in_hs),' q_hs_m =',svm(i,j,k,iq_hs)
          write(6,*) '  q_hgp =',q_hgp(i,j,k),' n_hgp =', n_hgp(i,j,k),' n_hg_m =', svm(i,j,k,in_hg),' q_hg_m =',svm(i,j,k,iq_hg)
        endif
        if (qtpmcr(i,j,k).lt.(-lim_qtp)) then
          ! calculate
          uncond   = ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond   = ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k,'  qtpmcr = ', qtpmcr (i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),' thl = ', thl0(i,j,k), ' n_ccn = ', n_cc(i,j,k)
          write(6,*) '  x_cl =', x_cl(i,j,k), ' n_cl =', n_cl(i,j,k),' q_cl =', q_cl (i,j,k)
          write(6,*) '  x_ci =', x_ci(i,j,k), ' n_ci =', n_ci(i,j,k),' q_ci =', q_ci (i,j,k)
          write(6,*) '  x_hr =', x_hr(i,j,k), ' n_hr =', n_hr(i,j,k),' q_hr =', q_hr (i,j,k)
          write(6,*) '  x_hs =', x_hs(i,j,k), ' n_hs =', n_hs(i,j,k),' q_hs =', q_hs (i,j,k)
          write(6,*) '  x_hg =', x_hg(i,j,k), ' n_hg =', n_hg(i,j,k),' q_hg =', q_hg (i,j,k)
          write(6,*) '  q_clp =',q_clp(i,j,k),' n_clp =', n_clp(i,j,k),' n_cl_m =', svm(i,j,k,in_cl),' q_cl_m =',svm(i,j,k,iq_cl)
          write(6,*) '  q_cip =',q_cip(i,j,k),' n_cip =', n_cip(i,j,k),' n_ci_m =', svm(i,j,k,in_ci),' q_ci_m =',svm(i,j,k,iq_ci)
          write(6,*) '  q_hrp =',q_hrp(i,j,k),' n_hrp =', n_hrp(i,j,k),' n_hr_m =', svm(i,j,k,in_hr),' q_hr_m =',svm(i,j,k,iq_hr)
          write(6,*) '  q_hsp =',q_hsp(i,j,k),' n_hsp =', n_hsp(i,j,k),' n_hs_m =', svm(i,j,k,in_hs),' q_hs_m =',svm(i,j,k,iq_hs)
          write(6,*) '  q_hgp =',q_hgp(i,j,k),' n_hgp =', n_hgp(i,j,k),' n_hg_m =', svm(i,j,k,in_hg),' q_hg_m =',svm(i,j,k,iq_hg)
        endif
    enddo
    enddo
    enddo
    endif ! flag
   endif
   endif ! flag_dbg

   end subroutine check_update

   subroutine check_allupdates(flag)
     use modglobal, only : ih,i1,jh,j1,k1
     use modfields, only : sv0,svm,svp,ql0,qvsl,qvsi,qt0   &
                           ,thl0,tmp0, qtp,qtm,thlp,rhof
   implicit none
    logical, intent (inout) ::flag
    logical     :: prob_thl, prob_qt
    ! character (len=7), intent (in) :: proc_name
    real :: lim_thlp, lim_qtp, ssat, ssat_ice, uncond
    integer:: i,j,k

   ! inner setting
    lim_thlp = 7.0
    lim_qtp = 0.0009
    prob_thl =  .false.
    prob_qt = .false.

   ! checking for large updates
   if(any(thlpmcr(2:i1,2:j1,1:k1).gt.lim_thlp).or.any(thlpmcr(2:i1,2:j1,1:k1).lt.(-lim_thlp))) then
    write(6,*) 'WARNING: thlpmcr problem'
    prob_thl = .true.
   endif
   if(any(thlp(2:i1,2:j1,1:k1).gt.lim_thlp).or.any(thlp(2:i1,2:j1,1:k1).lt.(-lim_thlp))) then
    write(6,*) 'WARNING: thlp problem'
    prob_thl = .true.
   endif
   if(any(qtpmcr(2:i1,2:j1,1:k1).gt.lim_qtp).or.any(qtpmcr(2:i1,2:j1,1:k1).lt.(-lim_qtp))) then
     write(6,*) 'WARNING: qtpmcr problem'
     prob_qt = .true.
   endif
   if(any(qtp(2:i1,2:j1,1:k1).gt.lim_qtp).or.any(qtp(2:i1,2:j1,1:k1).lt.(-lim_qtp))) then
     write(6,*) 'WARNING: qtp problem'
     prob_qt = .true.
   endif
   ! checking for update into negative numbers
   if(any((qtm(2:i1,2:j1,1:k1)+delt*qtpmcr(2:i1,2:j1,1:k1)).lt.0.0)) then
    write(6,*) 'WARNING: negative_qt problem after mphys'
    prob_qt = .true.
   endif
   if(any((qtm(2:i1,2:j1,1:k1)+delt*qtp(2:i1,2:j1,1:k1)).lt.0.0)) then
    write(6,*) 'WARNING: negative_qt problem'
    prob_qt= .true.
   endif
    if(any(n_hgp(2:i1,2:j1,1:k1).gt.100.0)) then
    write(6,*) 'WARNING: high n_hgp problem'
    prob_qt= .true.
   endif


   if (prob_thl .or. prob_qt)then
    if (flag) then
    do k=1,k1
    do j=2,j1
    do i=2,i1
        if ((thlpmcr(i,j,k).gt.lim_thlp).or.(thlpmcr(i,j,k).lt.(                   &
          -lim_thlp)).or.(qtpmcr(i,j,k).gt.lim_qtp).or.(qtpmcr(i,j,k).lt.(         &
          -lim_qtp)).or.(qtp(i,j,k).gt.lim_qtp).or.(qtp(i,j,k).lt.(                &
          -lim_qtp)).or.(thlp(i,j,k).gt.lim_thlp).or.(thlp(i,j,k).lt.(             &
          -lim_thlp)).or.((qtm(i,j,k)+delt*qtp(i,j,k)).lt.0.0).or.(n_hgp(i,j,k).gt.100))then
          ! calculate
          uncond   = ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond   = ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k
          write(6,*) '   thlpmcr=', thlpmcr(i,j,k),'qtpmcr=', qtpmcr(i,j,k)
          write(6,*) '   thlp=',thlp(i,j,k),'qtp=', qtp(i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),'thl = ', thl0(i,j,k), 'tmp0 = ', tmp0(i,j,k)
          write(6,*) '   rhof=',rhof(k),' n_ccn = ', n_cc(i,j,k)
          write(6,*) '   x_cl =', x_cl(i,j,k)
          write(6,*) '   x_ci =', x_ci(i,j,k)
          write(6,*) '   x_hr =', x_hr(i,j,k)
          write(6,*) '   x_hs =', x_hs(i,j,k)
          write(6,*) '   x_hg =', x_hg(i,j,k)
          ! outputs
          write(6,*) '   n_cl=', n_cl(i,j,k),'n_cl_m=', svm(i,j,k,in_cl),'n_cl_1=',svm(i,j,k,in_cl)+delt*svp(i,j,k,in_cl)
          write(6,*) '   n_ci=', n_ci(i,j,k),'n_ci_m=', svm(i,j,k,in_ci),'n_ci_1=',svm(i,j,k,in_ci)+delt*svp(i,j,k,in_ci)
          write(6,*) '   n_hr=', n_hr(i,j,k),'n_hr_m=', svm(i,j,k,in_hr),'n_hr_1=',svm(i,j,k,in_hr)+delt*svp(i,j,k,in_hr)
          write(6,*) '   n_hs=', n_hs(i,j,k),'n_hs_m=', svm(i,j,k,in_hs),'n_hs_1=',svm(i,j,k,in_hs)+delt*svp(i,j,k,in_hs)
          write(6,*) '   n_hg=', n_hs(i,j,k),'n_hg_m=', svm(i,j,k,in_hg),'n_hg_1=',svm(i,j,k,in_hg)+delt*svp(i,j,k,in_hg)

          write(6,*) '   q_cl=', q_cl (i,j,k), 'q_cl_m=',svm(i,j,k,iq_cl),'q_cl_1=',svm(i,j,k,iq_cl)+delt*svp(i,j,k,iq_cl)
          write(6,*) '   q_ci=', q_ci (i,j,k), 'q_ci_m=',svm(i,j,k,iq_ci),'q_ci_1=',svm(i,j,k,iq_ci)+delt*svp(i,j,k,iq_ci)
          write(6,*) '   q_hr=', q_hr (i,j,k), 'q_hr_m=',svm(i,j,k,iq_hr),'q_hr_1=',svm(i,j,k,iq_hr)+delt*svp(i,j,k,iq_hr)
          write(6,*) '   q_hs=', q_hs (i,j,k), 'q_hs_m=',svm(i,j,k,iq_hs),'q_hs_1=',svm(i,j,k,iq_hs)+delt*svp(i,j,k,iq_hs)
          write(6,*) '   q_hg=', q_hs (i,j,k), 'q_hg_m=',svm(i,j,k,iq_hg),'q_hg_1=',svm(i,j,k,iq_hg)+delt*svp(i,j,k,iq_hg)
          write(6,*) '   '
          ! updates  and past values
          write(6,*) '   n_clp =', n_clp(i,j,k),' svp n_cl=', svp(i,j,k,in_cl)
          write(6,*) '   n_cip =', n_cip(i,j,k),' svp n_ci=', svp(i,j,k,in_ci)
          write(6,*) '   n_hrp =', n_hrp(i,j,k),' svp n_hr=', svp(i,j,k,in_hr)
          write(6,*) '   n_hsp =', n_hsp(i,j,k),' svp n_hs=', svp(i,j,k,in_hs)
          write(6,*) '   n_hgp =', n_hgp(i,j,k),' svp n_hg=', svp(i,j,k,in_hg)
          write(6,*) '   q_clp =', q_clp(i,j,k),' svp q_cl=', svp(i,j,k,iq_cl)
          write(6,*) '   q_cip =', q_cip(i,j,k),' svp q_ci=', svp(i,j,k,iq_ci)
          write(6,*) '   q_hrp =', q_hrp(i,j,k),' svp q_hr=', svp(i,j,k,iq_hr)
          write(6,*) '   q_hsp =', q_hsp(i,j,k),' svp q_hs=', svp(i,j,k,iq_hs)
          write(6,*) '   q_hgp =', q_hgp(i,j,k),' svp q_hg=', svp(i,j,k,iq_hg)
          write(6,*) '   '
         ! and processes
          write(6,*) '   dn_cl_nu=',dn_cl_nu(i,j,k)
          write(6,*) '   dn_ci_inu=',dn_ci_inu(i,j,k)
          write(6,*) '   dq_hr_ev=' ,dq_hr_ev(i,j,k)
          write(6,*) '   dn_hr_ev=' ,dn_hr_ev(i,j,k)
          write(6,*) '   dq_hr_au=',dq_hr_au(i,j,k)
          write(6,*) '   dn_cl|au=',dn_cl_au(i,j,k)
          write(6,*) '   dq_hr_ac=', dq_hr_ac(i,j,k)
          write(6,*) '   dn_cl|ac=', dn_cl_ac(i,j,k)
          write(6,*) '   dn_hr_sc=',dn_hr_sc(i,j,k)
          write(6,*) '   dn_hr_br=',dn_hr_br(i,j,k)
          write(6,*) '   dq_ci_dep=',dq_ci_dep(i,j,k)
          write(6,*) '   dq_hs_dep=',dq_hs_dep(i,j,k)
          write(6,*) '   dq_hg_dep=',dq_hg_dep    (i,j,k)
          write(6,*) '   dq_ci_rime=',dq_ci_rime (i,j,k)
          write(6,*) '   dq_hs_rime =',dq_hs_rime (i,j,k)
          write(6,*) '   dq_hg_rime =',dq_hg_rime   (i,j,k)
          write(6,*) '   dn_cl_rime_ci=',dn_cl_rime_ci(i,j,k)
          write(6,*) '   dn_cl_rime_hs=',dn_cl_rime_hs(i,j,k)
          write(6,*) '   dn_cl_rime_hg=',dn_cl_rime_hg(i,j,k)
          write(6,*) '   dq_hshr_rime=',dq_hshr_rime (i,j,k)
          write(6,*) '   dq_hghr_rime=',dq_hghr_rime (i,j,k)
          write(6,*) '   dn_hr_rime_hs=',dn_hr_rime_hs(i,j,k)
          write(6,*) '   dn_hr_rime_hg=',dn_hr_rime_hg(i,j,k)
          write(6,*) '   dq_cl_hom=',dq_cl_hom    (i,j,k)
          write(6,*) '   dq_cl_het=' ,dq_cl_het    (i,j,k)
          write(6,*) '   dq_hr_het=',dq_hr_het    (i,j,k)
          write(6,*) '   dn_cl_hom=',dn_cl_hom    (i,j,k)
          write(6,*) '   dn_cl_het=',dn_cl_het    (i,j,k)
          write(6,*) '   dn_hr_het=',dn_hr_het    (i,j,k)
          write(6,*) '    '
          write(6,*) '   dq_ci_col_iis=',dq_ci_col_iis    (i,j,k)
          write(6,*) '   dq_hsci_col=',dq_hsci_col  (i,j,k)
          write(6,*) '   dn_hr_col_rs=', dn_hr_col_rs  (i,j,k)
          write(6,*) '   dn_hs_col_rs=', dn_hs_col_rs  (i,j,k)
          write(6,*) '   dn_ci_col_ri=', dn_ci_col_ri  (i,j,k)
          write(6,*) '   dn_hr_col_ri=', dn_hr_col_ri  (i,j,k)
          write(6,*) '   dn_hs_col_hg=', dn_hs_col_hg  (i,j,k)
          write(6,*) '   dn_ci_col_iis=',dn_ci_col_iis    (i,j,k)
          write(6,*) '   dn_hs_col_sss=' ,dn_hs_col_sss  (i,j,k)
          write(6,*) '   dn_ci_col_hs='  ,dn_ci_col_hs (i,j,k)
          write(6,*) '   dq_ci_cv=',dq_ci_cv (i,j,k)
          write(6,*) '   dq_hs_cv=' ,dq_hs_cv     (i,j,k)
          write(6,*) '   dn_ci_cv=',dn_ci_cv     (i,j,k)
          write(6,*) '   dn_hs_cv=' ,dn_hs_cv     (i,j,k)
          write(6,*) '    '
          write(6,*) '   dq_ci_me      ', dq_ci_me    (i,j,k)
          write(6,*) '   dq_hs_me      ', dq_hs_me    (i,j,k)
          write(6,*) '   dq_hg_me      ', dq_hg_me    (i,j,k)
          write(6,*) '   dq_ci_ev      ', dq_ci_ev   (i,j,k)
          write(6,*) '   dq_hs_ev     ', dq_hs_ev   (i,j,k)
          write(6,*) '   dq_hg_ev     ', dq_hg_ev   (i,j,k)
          write(6,*) '   dn_ci_me      ', dn_ci_me    (i,j,k)
          write(6,*) '   dn_hs_me      ', dn_hs_me    (i,j,k)
          write(6,*) '   dn_hg_me      ', dn_hg_me    (i,j,k)
          write(6,*) '   dn_ci_ev     ', dn_ci_ev   (i,j,k)
          write(6,*) '   dn_hs_ev     ', dn_hs_ev   (i,j,k)
          write(6,*) '   dn_hg_ev     ', dn_hg_ev   (i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_hs_eme_rs  ', dq_hs_eme_rs(i,j,k)
          write(6,*) '   dq_hs_eme_rs  ', dq_hs_eme_rs(i,j,k)
          write(6,*) '   dq_hg_eme_gc  ', dq_hg_eme_gc(i,j,k)
          write(6,*) '   dq_hg_eme_gr  ', dq_hg_eme_gr(i,j,k)
          write(6,*) '   dn_hs_eme_rs  ', dn_hs_eme_rs(i,j,k)
          write(6,*) '   dn_hs_eme_sc  ', dn_hs_eme_sc(i,j,k)
          write(6,*) '   dn_hg_eme_gc  ', dn_hg_eme_gc(i,j,k)
          write(6,*) '   dn_hg_eme_gr  ', dn_hg_eme_gr(i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_cl_se=',dq_cl_se     (i,j,k)
          write(6,*) '   dq_ci_se=',dq_ci_se     (i,j,k)
          write(6,*) '   dn_cl_sc=',dn_cl_sc     (i,j,k)
          write(6,*) '   dn_cl_se=',dn_cl_se     (i,j,k)
          write(6,*) '   dn_ci_se=',dn_ci_se (i,j,k)
          write(6,*) '   dq_hr_se=',dq_hr_se     (i,j,k)
          write(6,*) '   dq_hs_se=',dq_hs_se     (i,j,k)
          write(6,*) '   dq_hg_se=',dq_hg_se     (i,j,k)
          write(6,*) '   dn_hr_se=',dn_hr_se     (i,j,k)
          write(6,*) '   dn_hs_se=',dn_hs_se     (i,j,k)
          write(6,*) '   dn_hg_se=',dn_hg_se     (i,j,k)
          write(6,*) '    '
          write(6,*) '   precep_hs=',precep_hs    (i,j,k)
          write(6,*) '   precep_hg=',precep_hg    (i,j,k)
          write(6,*) '   dq_cl_sa=',dq_cl_sa     (i,j,k)
          write(6,*) '   ret_cc=',ret_cc       (i,j,k)
        endif
    enddo
    enddo
    enddo
    else ! flag
      write(6,*) 'mphys problems in: ',count(thlpmcr(2:i1,2:j1,1:k1).gt.lim_thlp),'and',count(thlpmcr(2:i1,2:j1,1:k1).lt.(-lim_thlp))
      write(6,*) 'mphys problems in: ',count(qtpmcr(2:i1,2:j1,1:k1).lt.(-lim_qtp)),'and',count(qtpmcr(2:i1,2:j1,1:k1).gt.lim_qtp)
    endif ! flag
   endif

   end subroutine check_allupdates


    ! checking for strange of incorrect resuls
   subroutine check_test(flag_dbg,flag, proc_name)
     use modglobal, only : ih,i1,jh,j1,k1
     use modfields, only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof,qvsl,qvsi,qt0, thl0
   implicit none
    logical, intent (inout) ::flag_dbg,flag
    character (len=7), intent (in) :: proc_name
    real :: lim_thlp, lim_qtp, ssat, ssat_ice, uncond
    integer:: i,j,k

   ! inner setting
    lim_thlp = 4.0
    lim_qtp = 0.0005

   if(flag_dbg) then

   ! checking for large updates
   if(any(thlp(2:i1,2:j1,1:k1).gt.lim_thlp).or.any(thlp(2:i1,2:j1,1:k1).lt.(-lim_thlp))) then
    write(6,*) 'WARNING: thl problem in: ', proc_name
    if (flag) then
      write(6,*) 'problems in: ',count(thlp(2:i1,2:j1,1:k1).lt.(-lim_qtp)),'and',count(thlp(2:i1,2:j1,1:k1).gt.lim_qtp)
    else ! flag
     flag = .true.
    do k=1,k1
    do j=2,j1
    do i=2,i1
        if (thlp(i,j,k).gt.lim_thlp) then
          ! calculate
          uncond=ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k,'  thlp = ', thlp(i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),' thl = ', thl0(i,j,k), ' n_ccn = ', n_cc(i,j,k)
          write(6,*) '  x_cl =', x_cl(i,j,k), ' n_cl =', n_cl(i,j,k),' q_cl =', q_cl (i,j,k)
          write(6,*) '  x_ci =', x_ci(i,j,k), ' n_ci =', n_ci(i,j,k),' q_ci =', q_ci (i,j,k)
          write(6,*) '  x_hr =', x_hr(i,j,k), ' n_hr =', n_hr(i,j,k),' q_hr =', q_hr (i,j,k)
          write(6,*) '  x_hs =', x_hs(i,j,k), ' n_hs =', n_hs(i,j,k),' q_hs =', q_hs (i,j,k)
          write(6,*) '  x_hg =', x_hg(i,j,k), ' n_hg =', n_hg(i,j,k),' q_hg =', q_hg (i,j,k)
          write(6,*) '  q_clp =',q_clp(i,j,k),' n_clp =', n_clp(i,j,k),' n_cl_m =', svm(i,j,k,in_cl),' q_cl_m =',svm(i,j,k,iq_cl)
          write(6,*) '  q_cip =',q_cip(i,j,k),' n_cip =', n_cip(i,j,k),' n_ci_m =', svm(i,j,k,in_ci),' q_ci_m =',svm(i,j,k,iq_ci)
          write(6,*) '  q_hrp =',q_hrp(i,j,k),' n_hrp =', n_hrp(i,j,k),' n_hr_m =', svm(i,j,k,in_hr),' q_hr_m =',svm(i,j,k,iq_hr)
          write(6,*) '  q_hsp =',q_hsp(i,j,k),' n_hsp =', n_hsp(i,j,k),' n_hs_m =', svm(i,j,k,in_hs),' q_hs_m =',svm(i,j,k,iq_hs)
          write(6,*) '  q_hgp =',q_hgp(i,j,k),' n_hgp =', n_hgp(i,j,k),' n_hg_m =', svm(i,j,k,in_hg),' q_hg_m =',svm(i,j,k,iq_hg)
        endif
        if (thlp(i,j,k).lt.(-lim_thlp)) then
          ! calculate
          uncond=ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k,'  thlp = ', thlp(i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),' thl = ', thl0(i,j,k), ' n_ccn = ', n_cc(i,j,k)
          write(6,*) '  x_cl =', x_cl(i,j,k), ' n_cl =', n_cl(i,j,k),' q_cl =', q_cl (i,j,k)
          write(6,*) '  x_ci =', x_ci(i,j,k), ' n_ci =', n_ci(i,j,k),' q_ci =', q_ci (i,j,k)
          write(6,*) '  x_hr =', x_hr(i,j,k), ' n_hr =', n_hr(i,j,k),' q_hr =', q_hr (i,j,k)
          write(6,*) '  x_hs =', x_hs(i,j,k), ' n_hs =', n_hs(i,j,k),' q_hs =', q_hs (i,j,k)
          write(6,*) '  x_hg =', x_hg(i,j,k), ' n_hg =', n_hg(i,j,k),' q_hg =', q_hg (i,j,k)
          write(6,*) '  q_clp =',q_clp(i,j,k),' n_clp =', n_clp(i,j,k),' n_cl_m =', svm(i,j,k,in_cl),' q_cl_m =',svm(i,j,k,iq_cl)
          write(6,*) '  q_cip =',q_cip(i,j,k),' n_cip =', n_cip(i,j,k),' n_ci_m =', svm(i,j,k,in_ci),' q_ci_m =',svm(i,j,k,iq_ci)
          write(6,*) '  q_hrp =',q_hrp(i,j,k),' n_hrp =', n_hrp(i,j,k),' n_hr_m =', svm(i,j,k,in_hr),' q_hr_m =',svm(i,j,k,iq_hr)
          write(6,*) '  q_hsp =',q_hsp(i,j,k),' n_hsp =', n_hsp(i,j,k),' n_hs_m =', svm(i,j,k,in_hs),' q_hs_m =',svm(i,j,k,iq_hs)
          write(6,*) '  q_hgp =',q_hgp(i,j,k),' n_hgp =', n_hgp(i,j,k),' n_hg_m =', svm(i,j,k,in_hg),' q_hg_m =',svm(i,j,k,iq_hg)
        endif
    enddo
    enddo
    enddo
    endif ! flag
   endif

   if(any(qtp(2:i1,2:j1,1:k1).lt.(-lim_qtp)).or.any(qtp(2:i1,2:j1,1:k1).gt.lim_qtp)) then
    write(6,*) 'WARNING: qtp problem in: ', proc_name
    if (flag) then
      write(6,*) 'problems in: ',count(qtp(2:i1,2:j1,1:k1).lt.(-lim_qtp)),'and',count(qtp(2:i1,2:j1,1:k1).gt.lim_qtp)
    else ! flag
     flag = .true.
    do k=1,k1
    do j=2,j1
    do i=2,i1
        if (qtp(i,j,k).gt.lim_qtp) then
          ! calculate
          uncond=ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k,'  qtp = ', qtp (i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),' thl = ', thl0(i,j,k), ' n_ccn = ', n_cc(i,j,k)
          write(6,*) '  x_cl =', x_cl(i,j,k), ' n_cl =', n_cl(i,j,k),' q_cl =', q_cl (i,j,k)
          write(6,*) '  x_ci =', x_ci(i,j,k), ' n_ci =', n_ci(i,j,k),' q_ci =', q_ci (i,j,k)
          write(6,*) '  x_hr =', x_hr(i,j,k), ' n_hr =', n_hr(i,j,k),' q_hr =', q_hr (i,j,k)
          write(6,*) '  x_hs =', x_hs(i,j,k), ' n_hs =', n_hs(i,j,k),' q_hs =', q_hs (i,j,k)
          write(6,*) '  x_hg =', x_hg(i,j,k), ' n_hg =', n_hg(i,j,k),' q_hg =', q_hg (i,j,k)
          write(6,*) '  q_clp =',q_clp(i,j,k),' n_clp =', n_clp(i,j,k),' n_cl_m =', svm(i,j,k,in_cl),' q_cl_m =',svm(i,j,k,iq_cl)
          write(6,*) '  q_cip =',q_cip(i,j,k),' n_cip =', n_cip(i,j,k),' n_ci_m =', svm(i,j,k,in_ci),' q_ci_m =',svm(i,j,k,iq_ci)
          write(6,*) '  q_hrp =',q_hrp(i,j,k),' n_hrp =', n_hrp(i,j,k),' n_hr_m =', svm(i,j,k,in_hr),' q_hr_m =',svm(i,j,k,iq_hr)
          write(6,*) '  q_hsp =',q_hsp(i,j,k),' n_hsp =', n_hsp(i,j,k),' n_hs_m =', svm(i,j,k,in_hs),' q_hs_m =',svm(i,j,k,iq_hs)
          write(6,*) '  q_hgp =',q_hgp(i,j,k),' n_hgp =', n_hgp(i,j,k),' n_hg_m =', svm(i,j,k,in_hg),' q_hg_m =',svm(i,j,k,iq_hg)
          !
        endif
        if (qtp(i,j,k).lt.(-lim_qtp)) then
            ! calculate
          uncond=ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k,'  qtp = ', qtp (i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),' thl = ', thl0(i,j,k), ' n_ccn = ', n_cc(i,j,k)
          write(6,*) '  x_cl =', x_cl(i,j,k), ' n_cl =', n_cl(i,j,k),' q_cl =', q_cl (i,j,k)
          write(6,*) '  x_ci =', x_ci(i,j,k), ' n_ci =', n_ci(i,j,k),' q_ci =', q_ci (i,j,k)
          write(6,*) '  x_hr =', x_hr(i,j,k), ' n_hr =', n_hr(i,j,k),' q_hr =', q_hr (i,j,k)
          write(6,*) '  x_hs =', x_hs(i,j,k), ' n_hs =', n_hs(i,j,k),' q_hs =', q_hs (i,j,k)
          write(6,*) '  x_hg =', x_hg(i,j,k), ' n_hg =', n_hg(i,j,k),' q_hg =', q_hg (i,j,k)
          write(6,*) '  q_clp =',q_clp(i,j,k),' n_clp =', n_clp(i,j,k),' n_cl_m =', svm(i,j,k,in_cl),' q_cl_m =',svm(i,j,k,iq_cl)
          write(6,*) '  q_cip =',q_cip(i,j,k),' n_cip =', n_cip(i,j,k),' n_ci_m =', svm(i,j,k,in_ci),' q_ci_m =',svm(i,j,k,iq_ci)
          write(6,*) '  q_hrp =',q_hrp(i,j,k),' n_hrp =', n_hrp(i,j,k),' n_hr_m =', svm(i,j,k,in_hr),' q_hr_m =',svm(i,j,k,iq_hr)
          write(6,*) '  q_hsp =',q_hsp(i,j,k),' n_hsp =', n_hsp(i,j,k),' n_hs_m =', svm(i,j,k,in_hs),' q_hs_m =',svm(i,j,k,iq_hs)
          write(6,*) '  q_hgp =',q_hgp(i,j,k),' n_hgp =', n_hgp(i,j,k),' n_hg_m =', svm(i,j,k,in_hg),' q_hg_m =',svm(i,j,k,iq_hg)
          !
        endif
    enddo
    enddo
    enddo
    endif
   endif
   endif ! flag_dbg

   end subroutine check_test


 !! point out- gives output in a point
     subroutine point_out(i,j,k, proc_name)
     use modglobal, only : ih,i1,jh,j1,k1
     use modfields, only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof,qvsl,qvsi,qt0, thl0
   implicit none
    ! logical, intent (inout) ::flag
    character (len=7), intent (in) :: proc_name
    real :: lim_thlp, lim_qtp, ssat, ssat_ice, uncond
    integer, intent (in) :: i,j,k

          ! calculate
          uncond=ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
       ! ssat(i,j,k) = (100./qvsl(i,j,k))*max( qt0(i,j,k)-qvsl(i,j,k),0.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k,' update in: ', proc_name
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),' thl = ', thl0(i,j,k), ' n_ccn = ', n_cc(i,j,k)
          write(6,*) '  x_cl =', x_cl(i,j,k), ' n_cl =', n_cl(i,j,k),' q_cl =', q_cl (i,j,k)
          write(6,*) '  x_ci =', x_ci(i,j,k), ' n_ci =', n_ci(i,j,k),' q_ci =', q_ci (i,j,k)
          write(6,*) '  x_hr =', x_hr(i,j,k), ' n_hr =', n_hr(i,j,k),' q_hr =', q_hr (i,j,k)
          write(6,*) '  x_hs =', x_hs(i,j,k), ' n_hs =', n_hs(i,j,k),' q_hs =', q_hs (i,j,k)
          write(6,*) '  x_hg =', x_hg(i,j,k), ' n_hg =', n_hg(i,j,k),' q_hg =', q_hg (i,j,k)
          write(6,*) '  q_clp =',q_clp(i,j,k),' n_clp =', n_clp(i,j,k),' n_cl_m =', svm(i,j,k,in_cl),' q_cl_m =',svm(i,j,k,iq_cl)
          write(6,*) '  q_cip =',q_cip(i,j,k),' n_cip =', n_cip(i,j,k),' n_ci_m =', svm(i,j,k,in_ci),' q_ci_m =',svm(i,j,k,iq_ci)
          write(6,*) '  q_hrp =',q_hrp(i,j,k),' n_hrp =', n_hrp(i,j,k),' n_hr_m =', svm(i,j,k,in_hr),' q_hr_m =',svm(i,j,k,iq_hr)
          write(6,*) '  q_hsp =',q_hsp(i,j,k),' n_hsp =', n_hsp(i,j,k),' n_hs_m =', svm(i,j,k,in_hs),' q_hs_m =',svm(i,j,k,iq_hs)
          write(6,*) '  q_hgp =',q_hgp(i,j,k),' n_hgp =', n_hgp(i,j,k),' n_hg_m =', svm(i,j,k,in_hg),' q_hg_m =',svm(i,j,k,iq_hg)


   end subroutine point_out


  subroutine check_ssat(flag_dbg)

    use modglobal, only : dzf,i1,j1,ih,jh,k1,kmax,rlv, cp  ! dzf,pi
    use modfields, only : w0, rhof,exnf, qvsl,svm,sv0,svp, qvsi, qt0,ql0,thl0,tmp0,thlp,qtp,rhof   ! ,exnf,ql0
    use modmpi,    only : myid
    ! -> add other needed variables
    implicit none
     ! character (len=10), intent (in) :: proc_name
     logical, intent (inout) ::flag_dbg ! ,flag
     integer :: i,j,k
     ! allocatable varibles - supersaturation, derivation of supersaturation
     real, allocatable :: ssat(:,:,:), ssice(:,:,:)
     real :: max_sat, uncond

    ! setting
    max_sat = 3.0 ! 9.0

    if (flag_dbg) then
    ! ------
    ! allocate
     allocate( ssat  (2-ih:i1+ih,2-jh:j1+jh,k1)      &
              ,ssice (2-ih:i1+ih,2-jh:j1+jh,k1)      &
             )
     ! allocate( dn_cl_nu (2-ih:i1+ih,2-jh:j1+jh,k1) )
      ssat = 0.0
      ssice = 0.0


    ! -----


    ! loops to determine supersaturation
     do k=1,k1
     do j=2,j1
     do i=2,i1
       ! calculating supersaturation
       ! q_4condense = (qt0 - q_condensed) -q_saturated(T)
       ssat(i,j,k) = (100./qvsl(i,j,k))* ((qt0(i,j,k)-q_cl(i,j,k))-qvsl(i,j,k))
       ssice(i,j,k) =(100.0/qvsi(i,j,k))*((qt0(i,j,k)-q_cl(i,j,k))-qvsi(i,j,k))
       !#iceout ssat(i,j,k) = (100./qvsl(i,j,k))*max( (qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)  )-qvsl(i,j,k),0.0)
       !#iceout ssice(i,j,k) =(100.0/qvsi(i,j,k))*max( (qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)  )-qvsi(i,j,k),0.0)
       ! ssat(i,j,k) = (100./qvsl(i,j,k))*max( qt0(i,j,k)-qvsl(i,j,k),0.0)
     enddo
     enddo
     enddo

     ! now check
    if(any(ssat(2:i1,2:j1,1:k1).gt.max_sat)) then
    ! if(any(ssat(2:i1,2:j1,1:k1).gt.max_sat).or.any(ssice(2:i1,2:j1,1:k1).gt.max_sat) ) then
    write(6,*) 'WARNING: saturation too high: '
    do k=1,k1
    do j=2,j1
    do i=2,i1
        if ( ssat(i,j,k).gt.max_sat) then
          ! calculate
         ! calculate
          uncond=qt0(i,j,k)-q_cl(i,j,k)-qvsl(i,j,k)
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          ! ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          ! ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k
          write(6,*) '   thlpmcr=', thlpmcr(i,j,k),'qtpmcr=', qtpmcr(i,j,k)
          write(6,*) '   thlp=',thlp(i,j,k),'qtp=', qtp(i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat(i,j,k),'% ','  ssat_i = ',ssice(i,j,k),'% '
          write(6,*) '   qt = ',qt0(i,j,k),'thl = ', thl0(i,j,k), 'tmp0 = ', tmp0(i,j,k)
          write(6,*) '   rhof=',rhof(k),' n_ccn = ', n_cc(i,j,k)
          write(6,*) '   x_cl =', x_cl(i,j,k)
          write(6,*) '   x_ci =', x_ci(i,j,k)
          write(6,*) '   x_hr =', x_hr(i,j,k)
          write(6,*) '   x_hs =', x_hs(i,j,k)
          write(6,*) '   x_hg =', x_hg(i,j,k)
          ! outputs
          write(6,*) '   n_cl=', n_cl(i,j,k),'n_cl_m=', svm(i,j,k,in_cl),'n_cl_1=',svm(i,j,k,in_cl)+delt*svp(i,j,k,in_cl)
          write(6,*) '   n_ci=', n_ci(i,j,k),'n_ci_m=', svm(i,j,k,in_ci),'n_ci_1=',svm(i,j,k,in_ci)+delt*svp(i,j,k,in_ci)
          write(6,*) '   n_hr=', n_hr(i,j,k),'n_hr_m=', svm(i,j,k,in_hr),'n_hr_1=',svm(i,j,k,in_hr)+delt*svp(i,j,k,in_hr)
          write(6,*) '   n_hs=', n_hs(i,j,k),'n_hs_m=', svm(i,j,k,in_hs),'n_hs_1=',svm(i,j,k,in_hs)+delt*svp(i,j,k,in_hs)
          write(6,*) '   n_hg=', n_hs(i,j,k),'n_hg_m=', svm(i,j,k,in_hg),'n_hg_1=',svm(i,j,k,in_hg)+delt*svp(i,j,k,in_hg)

          write(6,*) '   q_cl=', q_cl (i,j,k), 'q_cl_m=',svm(i,j,k,iq_cl),'q_cl_1=',svm(i,j,k,iq_cl)+delt*svp(i,j,k,iq_cl)
          write(6,*) '   q_ci=', q_ci (i,j,k), 'q_ci_m=',svm(i,j,k,iq_ci),'q_ci_1=',svm(i,j,k,iq_ci)+delt*svp(i,j,k,iq_ci)
          write(6,*) '   q_hr=', q_hr (i,j,k), 'q_hr_m=',svm(i,j,k,iq_hr),'q_hr_1=',svm(i,j,k,iq_hr)+delt*svp(i,j,k,iq_hr)
          write(6,*) '   q_hs=', q_hs (i,j,k), 'q_hs_m=',svm(i,j,k,iq_hs),'q_hs_1=',svm(i,j,k,iq_hs)+delt*svp(i,j,k,iq_hs)
          write(6,*) '   q_hg=', q_hs (i,j,k), 'q_hg_m=',svm(i,j,k,iq_hg),'q_hg_1=',svm(i,j,k,iq_hg)+delt*svp(i,j,k,iq_hg)
          ! updates  and past values
          write(6,*) '   n_clp =', n_clp(i,j,k),' svp q_cl=', svp(i,j,k,in_cl)
          write(6,*) '   n_cip =', n_cip(i,j,k),' svp q_ci=', svp(i,j,k,in_ci)
          write(6,*) '   n_hrp =', n_hrp(i,j,k),' svp q_hr=', svp(i,j,k,in_hr)
          write(6,*) '   n_hsp =', n_hsp(i,j,k),' svp q_hs=', svp(i,j,k,in_hs)
          write(6,*) '   n_hgp =', n_hgp(i,j,k),' svp q_hg=', svp(i,j,k,in_hg)
          write(6,*) '   q_clp =', q_clp(i,j,k),' svp q_cl=', svp(i,j,k,iq_cl)
          write(6,*) '   q_cip =', q_cip(i,j,k),' svp q_ci=', svp(i,j,k,iq_ci)
          write(6,*) '   q_hrp =', q_hrp(i,j,k),' svp q_hr=', svp(i,j,k,iq_hr)
          write(6,*) '   q_hsp =', q_hsp(i,j,k),' svp q_hs=', svp(i,j,k,iq_hs)
          write(6,*) '   q_hgp =', q_hgp(i,j,k),' svp q_hg=', svp(i,j,k,iq_hg)
         ! and processes
          write(6,*) '   dn_cl_nu=',dn_cl_nu(i,j,k)
          write(6,*) '   dn_ci_inu=',dn_ci_inu(i,j,k)
          write(6,*) '   dq_hr_ev=' ,dq_hr_ev(i,j,k)
          write(6,*) '   dn_hr_ev=' ,dn_hr_ev(i,j,k)
          write(6,*) '   dq_hr_au=',dq_hr_au(i,j,k)
          write(6,*) '   dn_cl|au=',(2.0/x_s)*dq_hr_au(i,j,k)
          write(6,*) '   dq_hr_ac=', dq_hr_ac(i,j,k)
          write(6,*) '   dn_cl|ac=', dq_hr_ac(i,j,k)/(x_cl(i,j,k)+eps0)
          write(6,*) '   dn_hr_sc=',dn_hr_sc(i,j,k)
          write(6,*) '   dn_hr_br=',dn_hr_br(i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_ci_dep=    ',dq_ci_dep(i,j,k)
          write(6,*) '   dq_hs_dep=    ',dq_hs_dep(i,j,k)
          write(6,*) '   dq_hg_dep=    ',dq_hg_dep    (i,j,k)
          write(6,*) '   dq_ci_rime=   ',dq_ci_rime (i,j,k)
          write(6,*) '   dq_hs_rime =  ',dq_hs_rime (i,j,k)
          write(6,*) '   dq_hg_rime =  ',dq_hg_rime   (i,j,k)
          write(6,*) '   dn_cl_rime_ci=',dn_cl_rime_ci(i,j,k)
          write(6,*) '   dn_cl_rime_hs=',dn_cl_rime_hs(i,j,k)
          write(6,*) '   dn_cl_rime_h= ',dn_cl_rime_hg(i,j,k)
          write(6,*) '   dq_hshr_rime= ',dq_hshr_rime (i,j,k)
          write(6,*) '   dq_hghr_rime= ',dq_hghr_rime (i,j,k)
          write(6,*) '   dn_hr_rime_hs=',dn_hr_rime_hs(i,j,k)
          write(6,*) '   dn_hr_rime_hg=',dn_hr_rime_hg(i,j,k)
          write(6,*) '   dq_hr_het=    ' ,dq_hr_het    (i,j,k)
          write(6,*) '   dq_cl_het=    ',dq_cl_het    (i,j,k)
          write(6,*) '   dq_cl_hom=    ',dq_cl_hom    (i,j,k)
          write(6,*) '   dn_hr_het=    ',dn_hr_het    (i,j,k)
          write(6,*) '   dn_cl_het=    ',dn_cl_het    (i,j,k)
          write(6,*) '   dn_cl_hom=    ',dn_cl_hom    (i,j,k)
          write(6,*) '   dq_ci_col_iis=    ',dq_ci_col_iis    (i,j,k)
          write(6,*) '   dq_hsci_col=  ',dq_hsci_col  (i,j,k)
          write(6,*) '   dn_ci_col_iis=    ',dn_ci_col_iis    (i,j,k)
          write(6,*) '   dn_hs_col_sss=    ',dn_hs_col_sss  (i,j,k)
          write(6,*) '   dn_ci_col_hs= ',dn_ci_col_hs (i,j,k)
          write(6,*) '   dq_ci_cv=     ',dq_ci_cv (i,j,k)
          write(6,*) '   dq_hs_cv=     ',dq_hs_cv     (i,j,k)
          write(6,*) '   dn_ci_cv=     ',dn_ci_cv     (i,j,k)
          write(6,*) '   dn_hs_cv=     ',dn_hs_cv     (i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_ci_me      ', dq_ci_me    (i,j,k)
          write(6,*) '   dq_hs_me      ', dq_hs_me    (i,j,k)
          write(6,*) '   dq_hg_me      ', dq_hg_me    (i,j,k)
          write(6,*) '   dq_ci_ev     ', dq_ci_ev   (i,j,k)
          write(6,*) '   dq_hs_ev     ', dq_hs_ev   (i,j,k)
          write(6,*) '   dq_hg_ev     ', dq_hg_ev   (i,j,k)
          write(6,*) '   dn_ci_me      ', dn_ci_me    (i,j,k)
          write(6,*) '   dn_hs_me      ', dn_hs_me    (i,j,k)
          write(6,*) '   dn_hg_me      ', dn_hg_me    (i,j,k)
          write(6,*) '   dn_ci_ev     ', dn_ci_ev   (i,j,k)
          write(6,*) '   dn_hs_ev     ', dn_hs_ev   (i,j,k)
          write(6,*) '   dn_hg_ev     ', dn_hg_ev   (i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_hs_eme_rs  ', dq_hs_eme_rs(i,j,k)
          write(6,*) '   dq_hs_eme_rs  ', dq_hs_eme_rs(i,j,k)
          write(6,*) '   dq_hg_eme_gc  ', dq_hg_eme_gc(i,j,k)
          write(6,*) '   dq_hg_eme_gr  ', dq_hg_eme_gr(i,j,k)
          write(6,*) '   dq_hs_eme_rs  ', dq_hs_eme_rs(i,j,k)
          write(6,*) '   dq_hs_eme_rs  ', dq_hs_eme_rs(i,j,k)
          write(6,*) '   dq_hg_eme_gc  ', dq_hg_eme_gc(i,j,k)
          write(6,*) '   dq_hg_eme_gr  ', dq_hg_eme_gr(i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_cl_se=     ',dq_cl_se     (i,j,k)
          write(6,*) '   dq_ci_se=     ',dq_ci_se     (i,j,k)
          write(6,*) '   dq_hr_se=     ',dq_hr_se     (i,j,k)
          write(6,*) '   dq_hs_se=     ',dq_hs_se     (i,j,k)
          write(6,*) '   dq_hg_se=     ',dq_hg_se     (i,j,k)
          write(6,*) '   dn_cl_se=     ',dn_cl_se     (i,j,k)
          write(6,*) '   dn_ci_se=     ',dn_ci_se     (i,j,k)
          write(6,*) '   dn_hr_se=     ',dn_hr_se     (i,j,k)
          write(6,*) '   dn_hs_se=     ',dn_hs_se     (i,j,k)
          write(6,*) '   dn_hg_se=     ',dn_hg_se     (i,j,k)
          write(6,*) '   '
          write(6,*) '   precep_hs=    ',precep_hs    (i,j,k)
          write(6,*) '   precep_hg=    ',precep_hg    (i,j,k)
          write(6,*) '   dq_cl_sa=     ',dq_cl_sa     (i,j,k)
          write(6,*) '   ret_cc=       ',ret_cc       (i,j,k)

          !
        endif
    enddo
    enddo
    enddo
   endif
   endif ! flag_dbg

   end subroutine check_ssat

   subroutine check_nan(flag_dbg,flag,proc_name)
     use modglobal, only : ih,i1,jh,j1,k1
     use modfields, only : sv0,svm,svp,ql0,qvsl,qvsi,qt0, thl0, qtp,thlp,tmp0,rhof
   implicit none
    logical, intent (inout) ::flag_dbg,flag
    character (len=7), intent (in) :: proc_name
    real :: lim_thlp, lim_qtp, ssat, ssat_ice, uncond
    logical ,allocatable ,dimension(:,:,:) :: nan_det ! ma_cl, ma_ci, ma_hr, ma_hs, ma_hg ! size testing
    integer:: i,j,k

    if(flag_dbg) then
    allocate (  nan_det(2-ih:i1+ih,2-jh:j1+jh,k1))
      ! allocate(  ma_cl(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             !   ,ma_ci(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             !   ,ma_hr(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             !   ,ma_hs(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             !   ,ma_hg(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             ! )
    ! fill with .false. at the start
    nan_det = .false.

    do k=1,k1
    do j=2,j1
    do i=2,i1
       if(isnan( n_clp(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( n_cip(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( n_hrp(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( n_hsp(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( n_hgp(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( q_clp(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( q_cip(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( q_hrp(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( q_hsp(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( q_hgp(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( thlpmcr(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( qtpmcr(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( thlp(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
       if(isnan( qtp(i,j,k))) then
          nan_det(i,j,k) = .true.
       endif
    enddo
    enddo
    enddo

    ! if any NaN values detected
    if(any(nan_det(2:i1,2:j1,1:k1))) then
    write(6,*) 'ERROR: NaN problem in: ', proc_name
    if (flag) then
      write(6,*) 'problems in: ',count(nan_det(2:i1,2:j1,1:k1))
    else ! flag
     flag = .true.
    do k=1,k1
    do j=2,j1
    do i=2,i1
        if ( nan_det(i,j,k)) then
          ! calculate
         ! calculate
          uncond=ql0(i,j,k)-q_cl(i,j,k)
          ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsl(i,j,k)-1.0)
          ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k))/qvsi(i,j,k)-1.0)
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k)
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci(i,j,k))/qvsi(i,j,k)-1.0)
          ! outputs
          write(6,*) ' (i,j,k) = ',i,j,k
          write(6,*) '   thlpmcr=', thlpmcr(i,j,k),'qtpmcr=', qtpmcr(i,j,k)
          write(6,*) '   thlp=',thlp(i,j,k),'qtp=', qtp(i,j,k)
          write(6,*) '   uncondensed = ', uncond,'  ssat = ',ssat,'% ','  ssat_i = ',ssat_ice,'% '
          write(6,*) '   qt = ',qt0(i,j,k),'thl = ', thl0(i,j,k), 'tmp0 = ', tmp0(i,j,k)
          write(6,*) '   rhof=',rhof(k),' n_ccn = ', n_cc(i,j,k)
          write(6,*) '   x_cl =', x_cl(i,j,k)
          write(6,*) '   x_ci =', x_ci(i,j,k)
          write(6,*) '   x_hr =', x_hr(i,j,k)
          write(6,*) '   x_hs =', x_hs(i,j,k)
          write(6,*) '   x_hg =', x_hg(i,j,k)
          ! outputs
          write(6,*) '   n_cl=', n_cl(i,j,k),'n_cl_m=', svm(i,j,k,in_cl),'n_cl_1=',svm(i,j,k,in_cl)+delt*svp(i,j,k,in_cl)
          write(6,*) '   n_ci=', n_ci(i,j,k),'n_ci_m=', svm(i,j,k,in_ci),'n_ci_1=',svm(i,j,k,in_ci)+delt*svp(i,j,k,in_ci)
          write(6,*) '   n_hr=', n_hr(i,j,k),'n_hr_m=', svm(i,j,k,in_hr),'n_hr_1=',svm(i,j,k,in_hr)+delt*svp(i,j,k,in_hr)
          write(6,*) '   n_hs=', n_hs(i,j,k),'n_hs_m=', svm(i,j,k,in_hs),'n_hs_1=',svm(i,j,k,in_hs)+delt*svp(i,j,k,in_hs)
          write(6,*) '   n_hg=', n_hs(i,j,k),'n_hg_m=', svm(i,j,k,in_hg),'n_hg_1=',svm(i,j,k,in_hg)+delt*svp(i,j,k,in_hg)

          write(6,*) '   q_cl=', q_cl (i,j,k), 'q_cl_m=',svm(i,j,k,iq_cl),'q_cl_1=',svm(i,j,k,iq_cl)+delt*svp(i,j,k,iq_cl)
          write(6,*) '   q_ci=', q_ci (i,j,k), 'q_ci_m=',svm(i,j,k,iq_ci),'q_ci_1=',svm(i,j,k,iq_ci)+delt*svp(i,j,k,iq_ci)
          write(6,*) '   q_hr=', q_hr (i,j,k), 'q_hr_m=',svm(i,j,k,iq_hr),'q_hr_1=',svm(i,j,k,iq_hr)+delt*svp(i,j,k,iq_hr)
          write(6,*) '   q_hs=', q_hs (i,j,k), 'q_hs_m=',svm(i,j,k,iq_hs),'q_hs_1=',svm(i,j,k,iq_hs)+delt*svp(i,j,k,iq_hs)
          write(6,*) '   q_hg=', q_hs (i,j,k), 'q_hg_m=',svm(i,j,k,iq_hg),'q_hg_1=',svm(i,j,k,iq_hg)+delt*svp(i,j,k,iq_hg)
          ! updates  and past values
          write(6,*) '   n_clp =', n_clp(i,j,k),' svp q_cl=', svp(i,j,k,in_cl)
          write(6,*) '   n_cip =', n_cip(i,j,k),' svp q_ci=', svp(i,j,k,in_ci)
          write(6,*) '   n_hrp =', n_hrp(i,j,k),' svp q_hr=', svp(i,j,k,in_hr)
          write(6,*) '   n_hsp =', n_hsp(i,j,k),' svp q_hs=', svp(i,j,k,in_hs)
          write(6,*) '   n_hgp =', n_hgp(i,j,k),' svp q_hg=', svp(i,j,k,in_hg)
          write(6,*) '   q_clp =', q_clp(i,j,k),' svp q_cl=', svp(i,j,k,iq_cl)
          write(6,*) '   q_cip =', q_cip(i,j,k),' svp q_ci=', svp(i,j,k,iq_ci)
          write(6,*) '   q_hrp =', q_hrp(i,j,k),' svp q_hr=', svp(i,j,k,iq_hr)
          write(6,*) '   q_hsp =', q_hsp(i,j,k),' svp q_hs=', svp(i,j,k,iq_hs)
          write(6,*) '   q_hgp =', q_hgp(i,j,k),' svp q_hg=', svp(i,j,k,iq_hg)
          write(6,*) '   '
          ! and processes
          write(6,*) '   dn_cl_nu=',dn_cl_nu(i,j,k)
          write(6,*) '   dn_ci_inu=',dn_ci_inu(i,j,k)
          write(6,*) '   dq_hr_ev=' ,dq_hr_ev(i,j,k)
          write(6,*) '   dn_hr_ev=' ,dn_hr_ev(i,j,k)
          write(6,*) '   dq_hr_au=',dq_hr_au(i,j,k)
          write(6,*) '   dn_cl|au=',(2.0/x_s)*dq_hr_au(i,j,k)
          write(6,*) '   dq_hr_ac=', dq_hr_ac(i,j,k)
          write(6,*) '   dn_cl|ac=', dq_hr_ac(i,j,k)/(x_cl(i,j,k)+eps0)
          write(6,*) '   dn_hr_sc=',dn_hr_sc(i,j,k)
          write(6,*) '   dn_hr_br=',dn_hr_br(i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_ci_dep=    ',dq_ci_dep(i,j,k)
          write(6,*) '   dq_hs_dep=    ',dq_hs_dep(i,j,k)
          write(6,*) '   dq_hg_dep=    ',dq_hg_dep    (i,j,k)
          write(6,*) '   dq_ci_rime=   ',dq_ci_rime (i,j,k)
          write(6,*) '   dq_hs_rime =  ',dq_hs_rime (i,j,k)
          write(6,*) '   dq_hg_rime =  ',dq_hg_rime   (i,j,k)
          write(6,*) '   dn_cl_rime_ci=',dn_cl_rime_ci(i,j,k)
          write(6,*) '   dn_cl_rime_hs=',dn_cl_rime_hs(i,j,k)
          write(6,*) '   dn_cl_rime_h= ',dn_cl_rime_hg(i,j,k)
          write(6,*) '   dq_hshr_rime= ',dq_hshr_rime (i,j,k)
          write(6,*) '   dq_hghr_rime= ',dq_hghr_rime (i,j,k)
          write(6,*) '   dn_hr_rime_hs=',dn_hr_rime_hs(i,j,k)
          write(6,*) '   dn_hr_rime_hg=',dn_hr_rime_hg(i,j,k)
          write(6,*) '   dq_hr_het=    ' ,dq_hr_het    (i,j,k)
          write(6,*) '   dq_cl_het=    ',dq_cl_het    (i,j,k)
          write(6,*) '   dq_cl_hom=    ',dq_cl_hom    (i,j,k)
          write(6,*) '   dn_hr_het=    ',dn_hr_het    (i,j,k)
          write(6,*) '   dn_cl_het=    ',dn_cl_het    (i,j,k)
          write(6,*) '   dn_cl_hom=    ',dn_cl_hom    (i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_ci_col_iis=    ',dq_ci_col_iis    (i,j,k)
          write(6,*) '   dq_hsci_col=  ',dq_hsci_col  (i,j,k)
          write(6,*) '   dn_ci_col_iis=    ',dn_ci_col_iis    (i,j,k)
          write(6,*) '   dn_hs_col_sss=    ',dn_hs_col_sss  (i,j,k)
          write(6,*) '   dn_ci_col_hs= ',dn_ci_col_hs (i,j,k)
          write(6,*) '   dq_ci_cv=     ',dq_ci_cv (i,j,k)
          write(6,*) '   dq_hs_cv=     ',dq_hs_cv     (i,j,k)
          write(6,*) '   dn_ci_cv=     ',dn_ci_cv     (i,j,k)
          write(6,*) '   dn_hs_cv=     ',dn_hs_cv     (i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_ci_me      ', dq_ci_me    (i,j,k)
          write(6,*) '   dq_hs_me      ', dq_hs_me    (i,j,k)
          write(6,*) '   dq_hg_me      ', dq_hg_me    (i,j,k)
          write(6,*) '   dq_ci_ev     ', dq_ci_ev   (i,j,k)
          write(6,*) '   dq_hs_ev     ', dq_hs_ev   (i,j,k)
          write(6,*) '   dq_hg_ev     ', dq_hg_ev   (i,j,k)
          write(6,*) '   dn_ci_me      ', dn_ci_me    (i,j,k)
          write(6,*) '   dn_hs_me      ', dn_hs_me    (i,j,k)
          write(6,*) '   dn_hg_me      ', dn_hg_me    (i,j,k)
          write(6,*) '   dn_ci_ev     ', dn_ci_ev   (i,j,k)
          write(6,*) '   dn_hs_ev     ', dn_hs_ev   (i,j,k)
          write(6,*) '   dn_hg_ev     ', dn_hg_ev   (i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_hs_eme_rs  ', dq_hs_eme_rs(i,j,k)
          write(6,*) '   dq_hs_eme_rs  ', dq_hs_eme_rs(i,j,k)
          write(6,*) '   dq_hg_eme_gc  ', dq_hg_eme_gc(i,j,k)
          write(6,*) '   dq_hg_eme_gr  ', dq_hg_eme_gr(i,j,k)
          write(6,*) '   dq_hs_eme_rs  ', dq_hs_eme_rs(i,j,k)
          write(6,*) '   dq_hs_eme_rs  ', dq_hs_eme_rs(i,j,k)
          write(6,*) '   dq_hg_eme_gc  ', dq_hg_eme_gc(i,j,k)
          write(6,*) '   dq_hg_eme_gr  ', dq_hg_eme_gr(i,j,k)
          write(6,*) '   '
          write(6,*) '   dq_cl_se=     ',dq_cl_se     (i,j,k)
          write(6,*) '   dq_ci_se=     ',dq_ci_se     (i,j,k)
          write(6,*) '   dq_hr_se=     ',dq_hr_se     (i,j,k)
          write(6,*) '   dq_hs_se=     ',dq_hs_se     (i,j,k)
          write(6,*) '   dq_hg_se=     ',dq_hg_se     (i,j,k)
          write(6,*) '   dn_cl_se=     ',dn_cl_se     (i,j,k)
          write(6,*) '   dn_ci_se=     ',dn_ci_se     (i,j,k)
          write(6,*) '   dn_hr_se=     ',dn_hr_se     (i,j,k)
          write(6,*) '   dn_hs_se=     ',dn_hs_se     (i,j,k)
          write(6,*) '   dn_hg_se=     ',dn_hg_se     (i,j,k)
          write(6,*) '   '
          write(6,*) '   precep_hs=    ',precep_hs    (i,j,k)
          write(6,*) '   precep_hg=    ',precep_hg    (i,j,k)
          write(6,*) '   dq_cl_sa=     ',dq_cl_sa     (i,j,k)
          write(6,*) '   ret_cc=       ',ret_cc       (i,j,k)
        endif
    enddo
    enddo
    enddo
    endif
    endif

   deallocate(  nan_det)
   endif

   end subroutine check_nan




  !*********************************************************************
  !*********************************************************************

  real function sed_flux3(Nin,Din,sig2,Ddiv,nnn)

  !*********************************************************************
  ! Function to calculate numerically the analytical solution of the
  ! sedimentation flux between Dmin and Dmax based on
  ! Feingold et al 1986 eq 17 -20.
  ! fall velocity is determined by alfa* D^beta with alfa+ beta taken as
  ! specified in Rogers and Yau 1989 Note here we work in D and in SI
  ! (in Roger+Yau in cm units + radius)
  ! flux is multiplied outside sed_flux with 1/rho_air to get proper
  ! kg/kg m/s units
  !
  ! M.C. van Zanten    August 2005
  !*********************************************************************
    use modglobal, only : pi,rhow
    implicit none

    real, intent(in) :: Nin, Din, sig2, Ddiv
    integer, intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real, parameter ::   C = rhow*pi/6.     &
                        ,D_intmin = 1e-6    &
                        ,D_intmax = 4.3e-3

    real ::  alfa         & ! constant in fall velocity relation
            ,beta         & ! power in fall vel. rel.
            ,D_min        & ! min integration limit
            ,D_max        & ! max integration limit
            ,flux           ![kg m^-2 s^-1]

    integer :: k

    flux = 0.0

    if (Din < Ddiv) then
      alfa = 3.e5*100  ![1/ms]
      beta = 2
      D_min = D_intmin
      D_max = Ddiv
      flux = C*Nin*alfa*erfint3(beta,Din,D_min,D_max,sig2,nnn)
    else
      do k = 1,3
        select case(k)
        case(1)        ! fall speed ~ D^2
          alfa = 3.e5*100 ![1/m 1/s]
          beta = 2
          D_min = Ddiv
          D_max = 133e-6
        case(2)        ! fall speed ~ D
          alfa = 4e3     ![1/s]
          beta = 1
          D_min = 133e-6
          D_max = 1.25e-3
        case default         ! fall speed ~ sqrt(D)
          alfa = 1.4e3 *0.1  ![m^.5 1/s]
          beta = .5
          D_min = 1.25e-3
          D_max = D_intmax
        end select
        flux = flux + C*Nin*alfa*erfint3(beta,Din,D_min,D_max,sig2,nnn)
      end do
    end if
      sed_flux3 = flux
  end function sed_flux3





  !*********************************************************************
  !*********************************************************************

  real function liq_cont3(Nin,Din,sig2,Ddiv,nnn)

  !*********************************************************************
  ! Function to calculate numerically the analytical solution of the
  ! liq. water content between Dmin and Dmax based on
  ! Feingold et al 1986 eq 17 -20.
  !
  ! M.C. van Zanten    September 2005
  !*********************************************************************
    use modglobal, only : pi,rhow
    implicit none

    real, intent(in) :: Nin, Din, sig2, Ddiv
    integer, intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real, parameter :: beta = 0           &
                      ,C = pi/6.*rhow     &
                      ,D_intmin = 80e-6    &   ! value of start of rain D
                      ,D_intmax = 3e-3         !4.3e-3    !  value is now max value for sqrt fall speed rel.

    real ::  D_min        & ! min integration limit
            ,D_max          ! max integration limit

    if (Din < Ddiv) then
    D_min = D_intmin
    D_max = Ddiv
    else
    D_min = Ddiv
    D_max = D_intmax
    end if

    liq_cont3 = C*Nin*erfint3(beta,Din,D_min,D_max,sig2,nnn)

  end function liq_cont3

  !*********************************************************************
  !*********************************************************************

  real function erfint3(beta, D, D_min, D_max, sig2,nnn )

  !*********************************************************************
  ! Function to calculate erf(x) approximated by a polynomial as
  ! specified in 7.1.27 in Abramowitz and Stegun
  ! NB phi(x) = 0.5(erf(0.707107*x)+1) but 1 disappears by substraction
  !
  !*********************************************************************
    implicit none
    real, intent(in) :: beta, D, D_min, D_max, sig2
    integer, intent(in) :: nnn

    real, parameter :: eps = 1e-10       &
                      ,a1 = 0.278393    & !a1 till a4 constants in polynomial fit to the error
                      ,a2 = 0.230389    & !function 7.1.27 in Abramowitz and Stegun
                      ,a3 = 0.000972    &
                      ,a4 = 0.078108
    real :: nn, ymin, ymax, erfymin, erfymax, D_inv

    D_inv = 1./(eps + D)
    nn = beta + nnn

    ymin = 0.707107*(log(D_min*D_inv) - nn*sig2)/(sqrt(sig2))
    ymax = 0.707107*(log(D_max*D_inv) - nn*sig2)/(sqrt(sig2))

    erfymin = 1.-1./((1.+a1*abs(ymin) + a2*abs(ymin)**2 + a3*abs(ymin)**3 +a4*abs(ymin)**4)**4)
    erfymax = 1.-1./((1.+a1*abs(ymax) + a2*abs(ymax)**2 + a3*abs(ymax)**3 +a4*abs(ymax)**4)**4)
    if (ymin < 0.) then
      erfymin = -1.*erfymin
    end if
    if (ymax < 0.) then
      erfymax = -1.*erfymax
    end if
    erfint3 = D**nn*exp(0.5*nn**2*sig2)*0.5*(erfymax-erfymin)
  !  if (erfint < 0.) write(*,*)'erfint neg'
    if (erfint3 < 0.) erfint3 = 0.
  end function erfint3

  ! #sb3 START
  real function calc_avent (nn, mu_a, nu_a, a_a,b_a, av)

  !*********************************************************************
  ! Function to calculate coefficient \bar{a}_vent,n
  ! for ventilation parameter
  ! specified in Appendix B in S&B
  !
  !*********************************************************************
    use modglobal, only : lacz_gamma ! LACZ_GAMMA
    implicit none
    real, intent(in) ::  mu_a, nu_a, a_a,b_a, av
    integer, intent(in) :: nn
    real :: arg11, arg12, arg21, arg22, mtt11, mtt12, mtt21,mtt22, expon02

    ! calculation itself
    arg11     = (nu_a+nn+b_a)/mu_a
    arg21     = (nu_a+1.0)/mu_a
    arg12     = (nu_a+1.0)/mu_a ! arg12     = (nu_a+nn+1.0)/mu_a
    arg22     = (nu_a+2.0)/mu_a
    expon02   = b_a+nn-1.0

    mtt11 = lacz_gamma(arg11)
    mtt21 = lacz_gamma(arg21)
    mtt12 = lacz_gamma(arg12)
    mtt22 = lacz_gamma(arg22)

    ! putting it together
    calc_avent = avf*(mtt11/mtt21)*(mtt12/mtt22)**expon02

  end function calc_avent

  real function calc_bvent (nn, mu_a, nu_a, a_a,b_a, beta_a, bv )

  !*********************************************************************
  ! Function to calculate coefficient \bar{a}_vent,n
  ! for ventilation parameter
  ! specified in Appendix B in S&B
  !
  !*********************************************************************
    use modglobal, only : lacz_gamma ! LACZ_GAMMA
    implicit none
    real, intent(in) ::  mu_a, nu_a, a_a,b_a, beta_a, bv
    integer, intent(in) :: nn
    real :: arg11, arg12, arg21, arg22, mtt11, mtt12, mtt21,mtt22, expon02

    ! calculation itself
    arg11     = (nu_a+nn+(3.0/2.0)*b_a+0.5*beta_a)/mu_a
    arg21     = (nu_a+1.0)/mu_a
    arg12     = (nu_a+1.0)/mu_a
    arg22     = (nu_a+2.0)/mu_a
    expon02   = (3.0/2.0)*b_a+0.5*beta_a+nn-1.0

    mtt11 = lacz_gamma(arg11)
    mtt21 = lacz_gamma(arg21)
    mtt12 = lacz_gamma(arg12)
    mtt22 = lacz_gamma(arg22)

    ! putting it together
    calc_bvent = bv*(mtt11/mtt21)*(mtt12/mtt22)**expon02

  end function calc_bvent
  ! #sb3 END


  ! #sb3 START
  real function calc_delta_b (kk, mu_b, nu_b, b_b)

  !*********************************************************************
  ! Function to calculate coefficient \delta^k_b
  !  for collision/collection
  ! specified in Appendix C in S&B
  !
  !*********************************************************************
    use modglobal, only : lacz_gamma ! LACZ_GAMMA
    implicit none
    real, intent(in) ::  mu_b, nu_b, b_b
    integer, intent(in) :: kk
    real :: arg11, arg12, arg21, arg22, mtt11, mtt12, mtt21,mtt22, expon02

    ! calculation itself
    arg11     = (2.0*b_b+nu_b+1+kk)/mu_b
    arg21     = (nu_b+1)/mu_b
    arg12     = (nu_b+1)/mu_b
    arg22     = (nu_b+2)/mu_b
    expon02   = 2.0*b_b+kk

    mtt11 = lacz_gamma(arg11)
    mtt21 = lacz_gamma(arg21)
    mtt12 = lacz_gamma(arg12)
    mtt22 = lacz_gamma(arg22)

    ! putting it together
    calc_delta_b = (mtt11/mtt21)*(mtt12/mtt22)**expon02

  end function calc_delta_b

  real function calc_delta_ab (kk, mu_a, nu_a, b_a, mu_b, nu_b, b_b)

  !*********************************************************************
  ! Function to calculate coefficient \delta^k_a,b
  ! for collision/collection
  ! specified in Appendix C in S&B
  !
  !*********************************************************************
    use modglobal, only : lacz_gamma ! LACZ_GAMMA
    implicit none
    real, intent(in) ::  mu_a, nu_a, b_a, mu_b, nu_b, b_b
    integer, intent(in) :: kk
    real :: arg11, arg12, arg21, arg22     &
           ,arg13, arg14, arg23, arg24     &
           ,mtt11, mtt12, mtt21,mtt22      &
           ,mtt13, mtt14, mtt23,mtt24      &
           ,exp03,exp04

    ! calculation itself
    arg11     = (b_a+nu_a+1+kk)/mu_a
    arg21     = (nu_a+1)/mu_a
    arg12     = (b_b+nu_b+1)/mu_b
    arg22     = (nu_b+1)/mu_b
    arg13     = arg21
    arg23     = (nu_a+2)/mu_a
    arg14     = arg22
    arg24     = (nu_b+2)/mu_b
    exp03     = b_a+kk
    exp04     = b_b

    mtt11 = lacz_gamma(arg11)
    mtt21 = lacz_gamma(arg21)
    mtt12 = lacz_gamma(arg12)
    mtt22 = lacz_gamma(arg22)
    mtt13 = mtt21
    mtt23 = lacz_gamma(arg23)
    mtt14 = mtt22
    mtt24 = lacz_gamma(arg24)

    ! putting it together
    calc_delta_ab = 2.0*(mtt11/mtt21)*(mtt12/mtt22)*(mtt13/mtt23)**exp03*(mtt14/mtt24)**exp04

  end function calc_delta_ab


  real function calc_th_b (kk, mu_b, nu_b, b_b, beta_b)

  !*********************************************************************
  ! Function to calculate coefficient \vartheta^k_a,b
  ! for collision/collection
  ! specified in Appendix C in S&B
  !
  !*********************************************************************
    use modglobal, only : lacz_gamma ! LACZ_GAMMA
    implicit none
    real, intent(in) ::  mu_b, nu_b, b_b, beta_b
    integer, intent(in) :: kk
    real :: arg11, arg12, arg21, arg22     &
           ,mtt11, mtt12, mtt21,mtt22      &
           ,exp02

    ! calculation itself
    arg11     = (2.0*beta_b+2.0*b_b+nu_b+1.0+kk)/mu_b
    arg21     = (2.0*b_b+nu_b+1.0+kk )/mu_b
    arg12     = (nu_b+1.0)/mu_b
    arg22     = (nu_b+2.0)/mu_b
    exp02     = 2.0*beta_b

    mtt11 = lacz_gamma(arg11)
    mtt21 = lacz_gamma(arg21)
    mtt12 = lacz_gamma(arg12)
    mtt22 = lacz_gamma(arg22)

    ! putting it together
    calc_th_b = (mtt11/mtt21)*(mtt12/mtt22)**exp02

  end function calc_th_b

  real function calc_th_ab (kk, mu_a, nu_a, b_a, beta_a, mu_b, nu_b, b_b, beta_b)

  !*********************************************************************
  ! Function to calculate coefficient \vartheta^k_a,b
  ! for collision/collection
  ! specified in Appendix C in S&B
  !
  !*********************************************************************
    use modglobal, only : lacz_gamma ! LACZ_GAMMA
    implicit none
    real, intent(in) ::  mu_a, nu_a, b_a, beta_a, mu_b, nu_b, b_b, beta_b
    integer, intent(in) :: kk
    real :: arg11, arg12, arg21, arg22     &
           ,arg13, arg14, arg23, arg24     &
           ,mtt11, mtt12, mtt21,mtt22      &
           ,mtt13, mtt14, mtt23,mtt24      &
           ,exp03,exp04

    ! calculation itself
    arg11     = (beta_a+b_a+nu_a+1+kk)/mu_a
    arg21     = (b_a+nu_a+1+kk)/mu_a
    arg12     = (beta_b+b_b+nu_b+1)/mu_b
    arg22     = (b_b+nu_b+1)/mu_b
    arg13     = (nu_a+1)/mu_a
    arg23     = (nu_a+2)/mu_a
    arg14     = (nu_b+1)/mu_b
    arg24     = (nu_b+2)/mu_b
    exp03     = beta_a
    exp04     = beta_b

    mtt11 = lacz_gamma(arg11)
    mtt21 = lacz_gamma(arg21)
    mtt12 = lacz_gamma(arg12)
    mtt22 = lacz_gamma(arg22)
    mtt13 = lacz_gamma(arg13)
    mtt23 = lacz_gamma(arg23)
    mtt14 = lacz_gamma(arg14)
    mtt24 = lacz_gamma(arg24)

    ! putting it together
    calc_th_ab = 2.0*(mtt11/mtt21)*(mtt12/mtt22)*(mtt13/mtt23)**exp03*(mtt14/mtt24)**exp04

  end function calc_th_ab

  real function calc_cons_mmt (kk, mu_a, nu_a)

  !*********************************************************************
  ! Function to calculate the constatnt part of the momement
  ! used in some constants
  ! specified in Appendix C in S&B
  !
  !*********************************************************************
    use modglobal, only : lacz_gamma ! LACZ_GAMMA
    implicit none
    real, intent(in) ::    mu_a, nu_a
    integer, intent(in) :: kk
    real :: arg11, arg12, arg21, arg22     &
           ,mtt11, mtt12, mtt21,mtt22      &
           ,exp02

    ! calculation itself
    arg11     = (nu_a+1+kk)/mu_a
    arg21     = (nu_a+1)/mu_a
    arg12     = (nu_a+1)/mu_a
    arg22     = (nu_a+2)/mu_a
    exp02     =  kk

    mtt11 = lacz_gamma(arg11)
    mtt21 = lacz_gamma(arg21)
    mtt12 = lacz_gamma(arg12)
    mtt22 = lacz_gamma(arg22)

    ! putting it together
    calc_cons_mmt = (mtt11/mtt21)*(mtt12/mtt22)**exp02

  end function calc_cons_mmt

  real function calc_cons_v (kk, mu_a, nu_a, al_a, be_a)

  !*********************************************************************
  ! Function to calculate the constatnt part of the momement
  ! used in some constants
  ! specified in Appendix C in S&B
  !
  !*********************************************************************
    use modglobal, only : lacz_gamma ! LACZ_GAMMA
    implicit none
    real, intent(in) ::    mu_a, nu_a, al_a, be_a
    integer, intent(in) :: kk
    real :: arg11, arg12, arg21, arg22     &
           ,mtt11, mtt12, mtt21,mtt22      &
           ,exp02

    ! calculation itself
    arg11     = (nu_a+be_a+1+kk)/mu_a
    arg21     = (nu_a+1+kk)/mu_a
    arg12     = (nu_a+1)/mu_a
    arg22     = (nu_a+2)/mu_a
    exp02     =  be_a

    mtt11 = lacz_gamma(arg11)
    mtt21 = lacz_gamma(arg21)
    mtt12 = lacz_gamma(arg12)
    mtt22 = lacz_gamma(arg22)

    ! putting it together
    calc_cons_v = al_a*(mtt11/mtt21)*(mtt12/mtt22)**exp02

  end function calc_cons_v


  !*********************************************************************
  ! Function to calculate the constatnt part of the momement
  ! used in some constants
  ! specified in Appendix C in S&B
  !
  !*********************************************************************
  real function calc_cons_lbd (mu_a, nu_a)


    use modglobal, only : lacz_gamma ! LACZ_GAMMA
    implicit none
    real, intent(in) ::    mu_a, nu_a
    ! integer, intent(in) :: kk
    real :: arg11, arg21                   &
           ,mtt11, mtt21                   &
           ,exp01

    ! calculation itself
    arg11     = (nu_a+1.0)/mu_a
    arg21     = (nu_a+2.0)/mu_a
    exp01     =  -mu_a

    mtt11 = lacz_gamma(arg11)
    mtt21 = lacz_gamma(arg21)

    ! putting it together
    calc_cons_lbd = (mtt11/mtt21)**exp01

  end function calc_cons_lbd


  !*********************************************************************
  ! Function to calculate the saturation pressure at a specific temperature
  !
  ! copied from subroutine initglobal in modglobal
  !
  !*********************************************************************
  real function  esl_at_te(te_a)
  real, intent(in) ::    te_a
  real             :: esl_try

   esl_try=exp(54.842763-6763.22/te_a-4.21*log(te_a)+         &
         0.000367*te_a+tanh(0.0415*(te_a-218.8))*                &
         (53.878-1331.22/te_a-9.44523*log(te_a)+ 0.014025*te_a))

   esl_at_te = max(0.0,esl_try)

  end function esl_at_te

end module modbulkmicro3
