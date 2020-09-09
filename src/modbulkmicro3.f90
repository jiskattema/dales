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
     call MPI_BCAST(tmp_inuc,          1, MY_REAL     ,0,comm3d,ierr)
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


     allocate(thlpmcr  (2-ih:i1+ih,2-jh:j1+jh,k1) &
             ,qtpmcr   (2-ih:i1+ih,2-jh:j1+jh,k1) &
             ,qr_spl   (2-ih:i1+ih,2-jh:j1+jh,k1) &
             ,Nr_spl   (2-ih:i1+ih,2-jh:j1+jh,k1)




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


    deallocate( &
              ,qr_spl,Nr_spl)
    deallocate(qltot,thlpmcr,qtpmcr)
    deallocate(n_cc,n_ccp                                     &
      ,n_cl,n_clp,n_ci,n_cip,n_hr,n_hrp,n_hs,n_hsp,n_hg,n_hgp &
      ,q_cl,q_clp,q_ci,q_cip,q_hr,q_hrp,q_hs,q_hsp,q_hg,q_hgp &
      ,x_cl, x_ci, x_hr, x_hs, x_hg                           &
      )

  end subroutine exitbulkmicro3

!> Calculates the microphysical source term.
  subroutine bulkmicro3
    use modglobal, only : i1,j1,k1,rdt,rk3step,timee,rlv,cp,ih,jh

    use modfields, only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof,qvsl,qvsi,qt0
    use modbulkmicrostat3, only : bulkmicrotend3,bulkmicrostat3 ! #sb3 #t7
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k
    integer :: insv, iqsv ! #sb3 - species number
    logical :: flag_warn_size,flag_do_dbg_up
    real :: qrtest,nrtest
    real, allocatable :: xtest (:,:,:)

    allocate(xtest(2-ih:i1+ih,2-jh:j1+jh,k1))  ! #d
    xtest = 0.0

    real ::   dn_cl_nu       &    !< droplet nucleation rate
             ,dn_ci_inu      &    !< ice nucleation rate
             ,dn_cl_au       &    !< change in number of cloud droplets due to autoconversion
             ,dq_hr_au       &    !< change in mass of raindrops due to autoconversion
             ,dn_hr_au       &    !< change in number of raindrops due to autoconversion
             ,dq_hr_ac       &    !< change in mass of raindrops due to accretion
             ,dn_cl_ac       &    !< change in number of cloud droplets due to accretion
             ,dn_hr_br       &
             ,dn_hr_sc       &
             ,dq_hr_ev       &
             ,dn_hr_ev       &
             ,dq_ci_rime     &    !< riming growth of ice
             ,dn_cl_rime_ci  &    !<  - and impact on n_cl
             ,dq_hs_rime     &    !< riming growth of snow
             ,dn_cl_rime_hs  &    !<  - and impact on n_cl
             ,dq_hg_rime     &    !< riming growth for graupel
             ,dn_cl_rime_hg  &    !<  - and impact on n_cl
             ,dq_hghr_rime   &    !< riming growth for graupel with rain
             ,dn_hr_rime_hg  &    !<  - and impact on n_hr
             ,dq_hshr_rime   &    !< riming growth for snow with rain
             ,dn_hr_rime_hs  &    !<  - and impact on n_hr
             ,dq_hr_rime_ri  &    !< rain loss from riming of ice+rain->gr
             ,dq_ci_rime_ri  &    !< ice loss from riming of ice+rain->gr
             ,dn_ci_rime_ri  &    !< ice number loss from riming of ice+rain->gr
             ,dn_hr_rime_ri  &    !< rain number loss from riming of ice+rain->gr
             ,dq_hr_col_rs   &    !< rain loss from riming of ice+snow->gr
             ,dq_hs_col_rs   &    !< rain number loss from riming of ice+snow->gr
             ,dn_hr_col_rs   &    !< snow loss from riming of ice+snow->gr
             ,dn_hs_col_rs   &    !< snow number loss from riming of ice+snow->gr
             ,dq_hr_col_ri   &    !< rain loss from riming of ice+rain->gr
             ,dq_ci_col_ri   &    !< ice loss from riming of ice+rain->gr
             ,dn_ci_col_ri   &    !< ice number loss from riming of ice+rain->gr
             ,dn_hr_col_ri   &    !< rain number loss from riming of ice+rain->gr
             ,dq_cl_het      &    !< heterogeneou freezing of cloud water
             ,dq_hr_het      &    !< heterogeneou freezing of raindrops
             ,dn_hr_het      &    !< heterogeneou freezing of raindrops
             ,dq_cl_hom      &    !< homogeneous freezing of cloud water
             ,dq_ci_col_iis  &    !< self-collection of cloud ice
             ,dn_ci_col_iis  &    !< self-collection of cloud ice
             ,dn_hs_col_sss  &    !< self-collection of snow
             ,dq_hsci_col    &    !< collection s+i - trend in q_hs
             ,dn_ci_col_hs   &    !< collection s+i - trend in n_ci
             ,dq_hghs_col    &    !< collection g+s - trend in q_hg
             ,dn_hs_col_hg   &    !< collection g+s - trend in n_hs
             ,dq_ci_cv       &    !< partial conversion ice -> graupel
             ,dn_ci_cv       &
             ,dq_hs_cv       &    !< partial conversion snow-> graupel
             ,dn_hs_cv       &
             ,dn_cl_sc       &    !< cloud self-collection
             ,dn_ci_mul      &      !< ice multiplication
             ,dq_ci_mul      &      !< ice multiplication
             ,dn_ci_me       &      !< number tendency melting of cloud ice
             ,dq_ci_me       &      !< mass tendency melting of cloud ice
             ,dn_hs_me       &      !< number tendency melting of snow
             ,dq_hs_me       &      !< mass tendency melting of snow
             ,dn_hg_me       &      !< number tendency melting of graupel
             ,dq_hg_me       &      !< mass tendency melting of graupel
             ,dn_ci_ev       &      !< number tendency evaporation of cloud ice
             ,dq_ci_ev       &      !< mass tendency evaporation of cloud ice
             ,dn_hs_ev       &      !< number tendency evaporation of snow
             ,dq_hs_ev       &      !< mass tendency evaporation of snow
             ,dn_hg_ev       &      !< number tendency evaporation of graupel
             ,dq_hg_ev       &      !< mass tendency evaporation of graupel
             ,dn_ci_eme_ic   &      !< number tendency enhanced melting of cloud ice by cloud water
             ,dq_ci_eme_ic   &      !< mass tendency enhanced melting of cloud ice by cloud water
             ,dn_ci_eme_ri   &      !< number tendency enhanced melting of cloud ice by rain
             ,dq_ci_eme_ri   &      !< mass tendency enhanced melting of cloud ice  by rain
             ,dn_hs_eme_sc   &      !< number tendency enhanced melting of snow by cloud water
             ,dq_hs_eme_sc   &      !< mass tendency enhanced melting of snow by cloud water
             ,dn_hs_eme_rs   &      !< number tendency enhanced melting of snow by rain
             ,dq_hs_eme_rs   &      !< mass tendency enhanced melting of snow by rain
             ,dn_hg_eme_gc   &      !< number tendency enhanced melting of graupel
             ,dq_hg_eme_gc   &      !< mass tendency enhanced melting of graupel
             ,dn_hg_eme_gr   &      !< number tendency enhanced melting of graupel by rain
             ,dq_hg_eme_gr   &      !< mass tendency enhanced melting of graupel by rain
             ,dn_cl_se       &      !< sedimentation for clouds water - number
             ,dq_cl_se       &      !<       -||-- mixing ration
             ,dn_ci_se       &      !< sedimentation for cloud ice - number
             ,dq_ci_se       &      !<       -||-- mixing ration
             ,dn_hr_se       &      !< sedimentation for rain - number
             ,dq_hr_se       &      !<       -||-- mixing ration
             ,dn_hs_se       &      !< sedimentation for snow - number
             ,dq_hs_se       &      !<       -||-- mixing ration
             ,dn_hg_se       &      !< sedimentation for graupel - number
             ,dq_hg_se              !<       -||-- mixing ration
             ,dq_cl_sa       &      !< saturation adjustment
             ,dn_cl_sa       &      !< change in n_cl due to saturation adjustment

    allocate (precep_hr     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of raindrops
             ,precep_ci     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of ice crystals
             ,precep_hs     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of snow
             ,precep_hg     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of graupel
             ,ret_cc        (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< recovery of ccn
            )

    real ::  D_cl  ! \bar{D}_cl mean diameter for cloud water particle
            ,D_ci  ! \bar{D}_ci mean diameter for cloud ice particle
            ,D_hr  ! \bar{D}_hr mean diameter for raindrops
            ,D_hs  ! \bar{D}_hs mean diameter for snow particle
            ,D_hg  ! \bar{D}_hg mean diameter for graupel particle
            ,v_cl  ! \bar{v}_cl mean velocity for cloud water droplets
            ,v_ci  ! \bar{v}_ci mean velocity for cloud ice particle
            ,v_hr  ! \bar{v}_hr mean velocity for raindrops
            ,v_hs  ! \bar{v}_hs mean velocity for snow particle
            ,v_hg  ! \bar{v}_hg mean velocity for graupel particle

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
      n_cl = sv0(i,j,k,in_cl)
      n_ci = sv0(i,j,k,in_ci)
      n_hr = sv0(i,j,k,in_hr)
      n_hs = sv0(i,j,k,in_hs)
      n_hg = sv0(i,j,k,in_hg)

      ! reading mass densities
      q_cl = sv0(i,j,k,iq_cl)
      q_ci = sv0(i,j,k,iq_ci)
      q_hr = sv0(i,j,k,iq_hr)
      q_hs = sv0(i,j,k,iq_hs)
      q_hg = sv0(i,j,k,iq_hg)
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

    ! reset flags
    flag_warn_size   = .true.
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
          if (q_hr < 0.)  then
            q_hr = 0.
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
          if (q_ci < 0.)  then
            q_ci = 0.
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
      if ((q_hr .gt. q_hr_min).and.(n_hr(i,j,k).gt.0.0))  then
         q_hr_mask = .true.
      else
         q_hr_mask = .false.
      endif
      ! snow :
      if ((q_hs(i,j,k) .gt. q_hs_min).and.(n_hs(i,j,k).gt.0.0))  then
         q_hs_mask = .true.
      else
         q_hs_mask = .false.
      endif
      ! graupel :
      if ((q_hg(i,j,k) .gt. q_hg_min).and.(n_hg(i,j,k).gt.0.0))  then
         q_hg_mask = .true.
      else
         q_hg_mask = .false.
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
      qltot = q_cl(i,j,k) + q_hr
      ! instead of: ql0  (i,j,k) + q_hr
      ! i.e - cloud water + rain areas

      ! clouds- original definition
      if (ql0(i,j,k) .ge. qcmin)  then
         qcmask = .true.
      else
         qcmask = .false.
      end if

      ! liquid clouds
      if ((q_cl(i,j,k).ge. qcliqmin).and.(n_cl(i,j,k).gt.0.0))  then
         q_cl_mask = .true.
      else
         q_cl_mask = .false.
      end if

      ! ice clouds

      ! calculating ice water content already existing
      ! qicew = q_ci+q_hs(i,j,k)+q_hg(i,j,k)

      if ((qci .ge. qicemin).and.(n_ci(i,j,k).gt.0.0))  then
         q_ci_mask = .true.
      else
         q_ci_mask = .false.
      endif
    enddo
    enddo
    enddo

    !*********************************************************************
    ! calculate Rain DSD integral properties & parameters xr, Dvr, lbdr
    !*********************************************************************
    call integrals_bulk3
    call nucleation3      ! cloud nucleation
    call icenucle3        ! ice nucleation  ! #Bb ! #iceout

    ! --------------------------------------------------------------
    ! freezing of water droplets
    ! --------------------------------------------------------------
    call homfreez3        ! homogeneous freezing of cloud droplets
    call hetfreez3        ! heterogeneous freezing

    ! --------------------------------------------------------------
    ! deposition processes
    ! --------------------------------------------------------------
    call deposit_ice3      ! deposition of vapour to cloud ice
    call deposit_snow3     ! deposition of vapour to snow 
    call deposit_graupel3  ! deposition of vapour to graupel 
    call cor_deposit3      ! correction for deposition

    ! --------------------------------------------------------------
    ! snow aggregation and self-collection
    ! collision processes for snow
    ! --------------------------------------------------------------
    call coll_sis3  ! snow selfcollection s+i -> s
    call coll_gsg3  ! snow selfcollection g+s -> g

    ! --------------------------------------------------------------
    ! - rimings
    ! --------------------------------------------------------------
    call coll_ici3 ! riming i+c -> i
    call coll_scs3 ! riming s+c -> s
    call coll_gcg3 ! riming g+c -> g
    call coll_grg3 ! riming g+r -> g

    ! --------------------------------------------------------------
    ! - raindrop freezing
    ! --------------------------------------------------------------
    call rainhetfreez3        ! heterogeneous freezing! #iceout

    ! --------------------------------------------------------------
    ! - collision with conversion
    ! --------------------------------------------------------------
    call coll_rig3 ! riming r+i -> g  !#t2 !#t4
    call coll_rsg3 ! riming r+i -> g  !#t7

    !--------------------------------------------------------------
    ! conversions
    !--------------------------------------------------------------
    call ice_multi3 ! ice multiplication of Hallet and Mossop

    if (l_sb_conv_par) then
      call conv_partial3 ! partial conversion
    endif

    ! --------------------------------------------------------------
    ! - melting and evaporation of ice particles
    ! --------------------------------------------------------------
    call evapmelting3 ! melting of ice particles

    ! --------------------------------------------------------------
    ! - basic warm processes
    !  --------------------------------------------------------------
    call autoconversion3
    call cloud_self3
    call accretion3
    call evap_rain3 ! rain evaporation

    ! =============================
    ! saturation adjustment
    ! =============================
    call satadj3

    !==================================
    ! sedimentation part
    !==================================
    call sedim_rain3 !jja sed_qr sed_Nr
    call sedim_cl3 !jja from routine local, sed_qip, sed_nip
    call sedim_ice3 !jja from routine local, sed_qip, sed_nip
    call sedim_snow3 !jja from routine local, sed_qip, sed_nip
    call sedim_graupel3 !#b2t3 !jja from routine local, sed_qip, sed_nip

    ! recovery of ccn aerosols
    call recover_cc

  !*********************************************************************
  ! remove negative values and non physical low values
  ! JJA: uses
  ! q_hrp, n_hrp, qtp, thlp,
  ! q_hsp, n_hsp, qtp, thlp,
  ! q_hgp, n_hgp, qtp, thlp,
  ! q_cip, n_cip, qtp, thlp,
  ! q_clp, n_clp
  ! thlpmcr, qtpmcr
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

  if(l_sb_dbg_extra) call check_sizes(flag_warn_size,'all_ups')
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


! ****************************************
!  Debugging subroutines
!
! ****************************************

   ! checking for strange of incorrect resuls

    ! checking for strange of incorrect resuls
   subroutine check_sizes(flag,proc_name)
     use modglobal, only : ih,i1,jh,j1,k1
     use modfields, only : sv0,svm,ql0,qvsl,qvsi,qt0, thl0, tmp0, thlp, qtp, rhof
   implicit none
    logical, intent (inout) ::flag
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


    if(.true.) then
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
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl-q_ci)/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl-q_ci)/qvsi(i,j,k)-1.0)
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
          write(6,*) '   n_cl =', n_cl(i,j,k),' q_cl =', q_cl
          write(6,*) '   n_ci =', n_ci(i,j,k),' q_ci =', q_ci
          write(6,*) '   n_hr =', n_hr(i,j,k),' q_hr =', q_hr
          write(6,*) '   n_hs =', n_hs(i,j,k),' q_hs =', q_hs
          write(6,*) '   n_hg =', n_hg(i,j,k),' q_hg =', q_hg
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
   endif

   end subroutine check_sizes

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
          !#iceout uncond   = ql0(i,j,k)-q_cl(i,j,k)-q_ci
          !#iceout ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci)/qvsl(i,j,k)-1.0)
          !#iceout ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci)/qvsi(i,j,k)-1.0)
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

          write(6,*) '   q_cl=', q_cl, 'q_cl_m=',svm(i,j,k,iq_cl),'q_cl_1=',svm(i,j,k,iq_cl)+delt*svp(i,j,k,iq_cl)
          write(6,*) '   q_ci=', q_ci, 'q_ci_m=',svm(i,j,k,iq_ci),'q_ci_1=',svm(i,j,k,iq_ci)+delt*svp(i,j,k,iq_ci)
          write(6,*) '   q_hr=', q_hr, 'q_hr_m=',svm(i,j,k,iq_hr),'q_hr_1=',svm(i,j,k,iq_hr)+delt*svp(i,j,k,iq_hr)
          write(6,*) '   q_hs=', q_hs, 'q_hs_m=',svm(i,j,k,iq_hs),'q_hs_1=',svm(i,j,k,iq_hs)+delt*svp(i,j,k,iq_hs)
          write(6,*) '   q_hg=', q_hs, 'q_hg_m=',svm(i,j,k,iq_hg),'q_hg_1=',svm(i,j,k,iq_hg)+delt*svp(i,j,k,iq_hg)
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
       !#iceout ssat(i,j,k) = (100./qvsl(i,j,k))*max( (qt0(i,j,k)-q_cl(i,j,k)-q_ci  )-qvsl(i,j,k),0.0)
       !#iceout ssice(i,j,k) =(100.0/qvsi(i,j,k))*max( (qt0(i,j,k)-q_cl(i,j,k)-q_ci  )-qvsi(i,j,k),0.0)
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
          !#iceout uncond=ql0(i,j,k)-q_cl(i,j,k)-q_ci
          ! ssat     = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci)/qvsl(i,j,k)-1.0)
          ! ssat_ice = 100.0*((qt0(i,j,k)-q_cl(i,j,k)-q_ci)/qvsi(i,j,k)-1.0)
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

          write(6,*) '   q_cl=', q_cl, 'q_cl_m=',svm(i,j,k,iq_cl),'q_cl_1=',svm(i,j,k,iq_cl)+delt*svp(i,j,k,iq_cl)
          write(6,*) '   q_ci=', q_ci, 'q_ci_m=',svm(i,j,k,iq_ci),'q_ci_1=',svm(i,j,k,iq_ci)+delt*svp(i,j,k,iq_ci)
          write(6,*) '   q_hr=', q_hr, 'q_hr_m=',svm(i,j,k,iq_hr),'q_hr_1=',svm(i,j,k,iq_hr)+delt*svp(i,j,k,iq_hr)
          write(6,*) '   q_hs=', q_hs, 'q_hs_m=',svm(i,j,k,iq_hs),'q_hs_1=',svm(i,j,k,iq_hs)+delt*svp(i,j,k,iq_hs)
          write(6,*) '   q_hg=', q_hs, 'q_hg_m=',svm(i,j,k,iq_hg),'q_hg_1=',svm(i,j,k,iq_hg)+delt*svp(i,j,k,iq_hg)
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
