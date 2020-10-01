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
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI, NLeSC
!

! TODO:
! subroutine check_sizes
! subroutine check_updates
!     compare updates to limit
! subroutine check_ssat(flag_dbg)
!     ssat(i,j,k) = (100.0/qvsl(i,j,k))*((qt0(i,j,k)-q_cl(i,j,k))-qvsl(i,j,k)) > max_sat
!     ssice(i,j,k) =(100.0/qvsi(i,j,k))*((qt0(i,j,k)-q_cl(i,j,k))-qvsi(i,j,k)) > max_sat


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
  use modmicrodata
  use modmicrodata3

  implicit none
  private
  public initbulkmicro3, exitbulkmicro3, bulkmicro3

  contains

!> Initializes and allocates the arrays
  subroutine initbulkmicro3
    use modglobal, only : lwarmstart,ifnamopt,fname_options,i1,ih,j1,jh,k1
    use modmpi,    only : myid,my_real,comm3d,mpi_logical
    implicit none
    integer :: ierr

    ! set some initial values before loading namelist
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
     ,Nc0, xc0_min, Nccn0                                    & ! setting of initial clouds
     ,l_statistics, l_tendencies                               ! output

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
  wfallmax_ci = d_wfallmax_hr   ! wfallmax_ci =   max(c_v_i0*x_ci_bmax**be_ci, c_v_i1*x_cl_bmax**be_ci)
  wfallmax_hs = d_wfallmax_hr   ! wfallmax_hs =   max(c_v_s0*x_hs_bmax**be_hs, c_v_s1*x_hs_bmax**be_hs)
  wfallmax_hg = max(c_v_g0, c_v_g1)

  ! coefficient for lbdr calculation l_sb_classic
  c_lbdr = calc_cons_lbd (mu_hr_cst, nu_hs_cst) ! TODO: Not used

  ! calculating eslt3 = esl(T_3)/T_3
  ! using thermodynamic routine
  eslt3 = esl_at_te(T_3)

   if (myid .eq. 0) then
     write (6,*) ' === Settings of Mphys ===='
     write (6,*) 'imicro = ', imicro
     write (6,*) 'l_sb = ', l_sb
     write (6,*) 'l_lognormal', l_lognormal
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
     write (6,*)   '   dlt_i0    =', dlt_i0
     write (6,*)   '   dlt_0aa   = 2*dlt_i0 + dlt_i0i = ',2*dlt_i0 + dlt_i0i
     write (6,*)   '   dlt_1a    = dlt_i1 = ',dlt_i1
     write (6,*)   '   dlt_1aa   = dlt_i0 + dlt_i1i + dlt_i1 = ',dlt_i0 + dlt_i1i + dlt_i1
     write (6,*)   '   th_0a     = th_i0 = ',th_i0
     write (6,*)   '   th_0aa    = 2*th_i0 - th_i0i = ', 2*th_i0 - th_i0i ! from Seifert, 2002
     write (6,*)   '   th_1a     = th_i1 = ', th_i1
     write (6,*)   '   th_1aa    = th_i0 - th_i1i + th_i1 = ', th_i0 - th_i1i + th_i1  ! th_i0i - th_i1i + th_i1 ! th_1ab    = th_i1i
     write (6,*)   ' --------------------------------------------------- '

   endif

   allocate(precep_l     (2-ih:i1+ih,2-jh:j1+jh) &
           ,precep_i     (2-ih:i1+ih,2-jh:j1+jh) )
   allocate(statistics_patch(nmphys,k1) &
           ,tend_patch(ntends,k1)       )
  end subroutine initbulkmicro3


!> Cleaning up after the run
! ---------------------------------------------------------------------------------
subroutine exitbulkmicro3
  implicit none
  deallocate(precep_l, precep_i)
  deallocate(statistics_patch,tend_patch)
end subroutine exitbulkmicro3


!> Calculates the microphysical source term.
! ---------------------------------------------------------------------------------
subroutine bulkmicro3
  use modglobal, only : i1,j1,k1,rdt,rk3step
  use modfields, only : exnf,rhof,presf
  use modbulkmicrostat3, only : bulkmicrotend3, bulkmicrostat3
  use modmicrodata3, only : in_cl
  use modmpi,    only : myid
  use modbulkmicro3_point, only : point_processes
  use modbulkmicro3_column, only : nucleation3, column_processes
  implicit none
  integer :: i,j,k,ks,ke

  integer :: k_low(ncols)    & ! lowest k with non-zero values for svp(:,:,k,i)
            ,k_high(ncols)     ! highest k with non-zero values for svp(:,:,k,i)

  ! transposed arrays with optimized memory layout
  real :: sv0_t   (ncols,k1,i1+3,j1) &
         ,svp_t   (ncols,k1,i1+3,j1) &
         ,svm_t   (ncols,k1,i1+3,j1) &
         ,prg_t   (nprgs,k1,i1+3,j1) &
         ,thlp_t  (      k1,i1  ,j1) &
         ,qtp_t   (      k1,i1  ,j1)

  ! column variables
  real :: tend_col      (ntends, k1)  &
         ,statistics_col(nmphys, k1)

  ! and variables as the surface (k=1) of the column
  real :: precep_hr,precep_ci,precep_hs,precep_hg


  ! check if ccn and clouds were already initialised
  ! ------------------------------------------------
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

  ! Zero the summed statistics and tendencies
  ! no need to zero:
  ! thlp_t, qtp_t      : they are set in point_processes
  ! precep_i, precep_l : they are set at the end of the column loop
  if (l_tendencies) then
    tend_patch = 0.
    tend_col = 0.
  endif
  if (l_statistics) then
    statistics_patch = 0.
    statistics_col = 0.
  endif

  delt = rdt / (4. - dble(rk3step))

  ! Re-order the data in column format
  ! NOTE: svp is also adjusted by advection and canopy routines,
  !       so we cannot assume they are zero here.
  !       For qtp and thlp this is not an issue, as the microphysics scheme
  !       does not use those, but only adds to them.
  !
  ! BUG: how to treat existing values in svp, wrt. limiting and removing low values?
  !      Current code works on *total* tendencies, which means we include interactions
  !      with fi. the canopy.
  ! --------------------------------------------------------------------
  call transpose_svs(sv0_t, svm_t, svp_t, prg_t)

  ! Prevent doing unnecessary calculations (updating tendencies etc.)
  ! for levels without any non-zero values for svp.
  ! Reset the values here to all-levels-zero, and update them when
  ! removing negative and low values, below.
  k_low = k1
  k_high = 1

  ! loop over all (i,j) columns
  do i=2,i1
  do j=2,j1

  ! Column processes
  ! ------------------------------------------------------------------
    call nucleation3(prg_t(n_qt0,:,i,j), prg_t(n_qvsl,:,i,j), prg_t(n_w0,:,i,j)  &
                    ,sv0_t(iq_cl,:,i,j), svp_t(iq_cl,:,i,j)                      &
                    ,sv0_t(in_cl,:,i,j), svp_t(in_cl,:,i,j)                      &
                    ,sv0_t(in_cc,:,i,j)                                          &
                    ,statistics_col, tend_col                                    )

  !  - Point processes at k-point
  ! ------------------------------------------------------------------
    if (rk3step /= 1) then
      do k=1,k1
        call point_processes(prg_t(:,k,i,j),exnf(k),rhof(k),presf(k)         &
                            ,sv0_t(:,k,i,j),svp_t(:,k,i,j),svm_t(:,k,i,j)    &
                            ,thlp_t(k,i,j),qtp_t(k,i,j)                      &
                            ,statistics_col(:,k), tend_col(:,k)              )
      enddo
    else ! rk3step
      ! use svm == sv0
      do k=1,k1
        call point_processes(prg_t(:,k,i,j),exnf(k),rhof(k),presf(k)         &
                            ,sv0_t(:,k,i,j),svp_t(:,k,i,j),sv0_t(:,k,i,j)    &
                            ,thlp_t(k,i,j),qtp_t(k,i,j)                      &
                            ,statistics_col(:,k), tend_col(:,k)              )
      enddo
    endif ! rk3step


  ! Column processes
  ! ------------------------------------------------------------------
    call column_processes(sv0_t(:,:,i,j),svp_t(:,:,i,j)            &
                         ,thlp_t(:,i,j),qtp_t(:,i,j)               &
                         ,precep_hr,precep_ci,precep_hs,precep_hg  &
                         ,tend_col                                 )

  ! Remove negative values and non physical low values
  ! -----------------------------------------------------------------
  if (rk3step /= 1)  then
    call correct_neg_qt(svp_t(:,:,i,j),svm_t(:,:,i,j)            &
                       ,thlp_t(:,i,j),qtp_t(:,i,j),k_low,k_high  )
  else
    ! use svm == sv0
    call correct_neg_qt(svp_t(:,:,i,j),sv0_t(:,:,i,j)            &
                       ,thlp_t(:,i,j),qtp_t(:,i,j),k_low,k_high  )
  endif ! rk3step

  ! Keep track of output
  ! -----------------------------------------------------------------
    precep_l(i,j) = precep_hr
    precep_i(i,j) = precep_ci + precep_hs + precep_hg

    ! Accumulate non-zero tendencies and statistics, and reset to zero.
    ! NOTE BUG: can we have tendencies while all svp() are exactly zero,
    !       but that seems extremely unlikely.
    if (l_tendencies) then
      ks = minval(k_low(:))
      ke = maxval(k_high(:))
      if (ks <= ke) then
        do k=ks,ke
          tend_patch(:,k) = tend_patch(:,k) + tend_col(:,k)
        enddo
        tend_col(:,ks:ke) = 0.
      endif
    endif
    if (l_statistics) then
      ks = minval(k_low(:))
      ke = maxval(k_high(:))
      if (ks <= ke) then
        do k=ks,ke
          statistics_patch(:,k) = statistics_patch(:,k) + statistics_col(:,k)
        enddo
        statistics_col(:,ks:ke) = 0.
      endif
    endif

  enddo ! loop over i
  enddo ! loop over j


  ! update main prognostic variables by contribution from mphys processes
  ! ------------------------------------------------------------------
  call untranspose_svs(svp_t,thlp_t,qtp_t,k_low,k_high)

  ! microphysics statistics - just once per step
  ! ------------------------------------------------------------------
  !call bulkmicrotend3 ! #t5
  !call bulkmicrostat3 ! #t5
end subroutine bulkmicro3


! Re-order the data to column layout
!
! When the arrays get larger, they will no longer fit in the cache and data
! needs to be fetched from a next cache, and finally the main memory.
! This is very slow.
! In this subroutine we transpose the data such that all elements necessary
! to calculate a grid point are close together (ie. on the same memory page).
! Ideally, the microphysics can then run purely on the cache.
! The transpose itself also accesses a lot of memory, but the total amount of
! cache misses is reduced because of the cache-friendly(er) access pattern.
!
! tmp0,qt0,ql0,esl,qvsl,qvsi,w0  (i,j,k)  => prg_t (iprg,k,j)
!
! sv0,svm,svp0 (i,j,k,isv)  =>  sv0_t, svm_t, svp_t (isv,k,i,j)
!
! NOTE: Loop ordering and partial unrolling to improve performance
!       Total performance is quite sensitive, please test before making any changes
!       We cant do the dont-copy-zeros trick, as we need to initialize the arrays anyways
! ----------------------------------------------
subroutine transpose_svs(sv0_t, svm_t, svp_t, prg_t)
  use modfields, only : tmp0, qt0, ql0, esl, qvsl, qvsi, w0
  use modfields, only : sv0, svm, svp
  use modglobal, only : i1,j1,k1,rk3step
  implicit none
  real, intent(out) ::  sv0_t   (ncols,k1,i1+3,j1) &
                       ,svp_t   (ncols,k1,i1+3,j1) &
                       ,svm_t   (ncols,k1,i1+3,j1) &
                       ,prg_t   (nprgs,k1,i1+3,j1)

  integer :: i,j,k,isv

  do j=2,j1
  do k=1,k1
  do i=2,i1
    prg_t(n_tmp0,k,i,j) = tmp0(i,j,k)
    prg_t(n_qt0 ,k,i,j) = qt0 (i,j,k)
    prg_t(n_ql0 ,k,i,j) = ql0 (i,j,k)
    prg_t(n_esl ,k,i,j) = esl (i,j,k)
    prg_t(n_qvsl,k,i,j) = qvsl(i,j,k)
    prg_t(n_qvsi,k,i,j) = qvsi(i,j,k)
    prg_t(n_w0  ,k,i,j) = w0  (i,j,k)
  enddo
  enddo
  enddo

  ! BUG: check rk3step details
  ! for rk3step == 1 we have sv0 = svm, so no need for to transpose it too
  if (rk3step /= 1) then
  do j=2,j1
  do k=1,k1
  do i=2,i1,4
  do isv=1,ncols
    svm_t(isv,k,i+0,j) = max(svm(i+0,j,k,isv), 0.)
    svm_t(isv,k,i+1,j) = max(svm(i+1,j,k,isv), 0.)
    svm_t(isv,k,i+2,j) = max(svm(i+2,j,k,isv), 0.)
    svm_t(isv,k,i+3,j) = max(svm(i+3,j,k,isv), 0.)
  enddo
  enddo
  enddo
  enddo
  endif

  do j=2,j1
  do k=1,k1
  do i=2,i1,4
  do isv=1,ncols
    sv0_t(isv,k,i+0,j) = max(sv0(i+0,j,k,isv), 0.)
    sv0_t(isv,k,i+1,j) = max(sv0(i+1,j,k,isv), 0.)
    sv0_t(isv,k,i+2,j) = max(sv0(i+2,j,k,isv), 0.)
    sv0_t(isv,k,i+3,j) = max(sv0(i+3,j,k,isv), 0.)
  enddo
  enddo
  enddo
  enddo

  do j=2,j1
  do k=1,k1
  do i=2,i1,4
  do isv=1,ncols
    svp_t(isv,k,i+0,j) = svp(i+0,j,k,isv)
    svp_t(isv,k,i+1,j) = svp(i+1,j,k,isv)
    svp_t(isv,k,i+2,j) = svp(i+2,j,k,isv)
    svp_t(isv,k,i+3,j) = svp(i+3,j,k,isv)
  enddo
  enddo
  enddo
  enddo

endsubroutine transpose_svs


! copying from the columns to the full 3d fields
!
! NOTE: this will be slow when the arrays are large,
!       as we have a memory-bandwidth bottleneck
! ----------------------------------------------
subroutine untranspose_svs(svp_t,thlp_t,qtp_t,k_low,k_high)
  use modfields, only : svp,thlp,qtp
  use modglobal, only : i1,j1,k1
  implicit none
  real, intent(in) ::  svp_t(ncols,k1,i1+3,j1) &
                      ,thlp_t(k1,i1,j1)        &
                      ,qtp_t (k1,i1,j1)
  integer, intent(in) :: k_low(ncols), k_high(ncols)

  integer :: i,j,k,isv,ks,ke

  ! We initialized svp_t from svp, so now we can replace the svp values.
  ! To be a bit more efficient, go over the isv in pairs,
  ! which should be the in_XX and iq_XX fields.

  do isv=1,ncols,2 ! NOTE: ncols should be even!

  ks = min(k_low(isv),k_low(isv+1))
  ke = max(k_high(isv),k_high(isv+1))
  if(ks.le.ke) then
    do k=ks,ke
    do j=2,j1
    do i=2,i1
       svp(i,j,k,isv+0) = svp_t(isv+0,k,i,j)
       svp(i,j,k,isv+1) = svp_t(isv+1,k,i,j)
    enddo ! j
    enddo ! i
    if (ks.gt. 1) svp(:,:,1:ks ,isv:isv+1) = 0.
    if (ke.lt.k1) svp(:,:,ke:k1,isv:isv+1) = 0.
    enddo ! k
  else
    svp(:,:,:,isv:isv+1) = 0.
  endif
  enddo ! isv

  ! Add the microphysics tendencies to the total tendencies
  ks = minval(k_low)
  ke = maxval(k_high)
  if (ks.le.ke) then
    do j=1,j1
    do i=2,i1 - 1,2 ! careful not to go out-of-array
    do k=ks,ke
      thlp(i+0,j,k) = thlp(i+0,j,k) + thlp_t(k,i+0,j)
      thlp(i+1,j,k) = thlp(i+1,j,k) + thlp_t(k,i+1,j)

      qtp (i+0,j,k) = qtp (i+0,j,k) + qtp_t (k,i+0,j)
      qtp (i+1,j,k) = qtp (i+1,j,k) + qtp_t (k,i+1,j)
    enddo ! k
    enddo ! i

    ! Do corner case i = i1 and i1 is odd
    if (mod(i1,2) == 1) then
      do k=ks,ke
        thlp(i1,j,k) = thlp(i1,j,k) + thlp_t(k,i1,j)
        qtp (i1,j,k) = qtp (i1,j,k) + qtp_t (k,i1,j)
      enddo ! k
    endif
    enddo ! j
  endif
endsubroutine untranspose_svs


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
! ----------------------------------------------------------------------------
subroutine initclouds3
use modglobal, only : i1,j1,k1
use modfields, only : rhof, qt0, svm, sv0, qvsl
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


!  ccn initialisation
! ===============================
! initialises cloud water and cloud droplet number when called
!  called from : bulkmicro3 if flag l_ccn_init==.false.
! steps:
!     - based on proposed values Nc0 [ m^-3]
!     - set ccn fields
! ------------------------------------------------------------------------------
subroutine initccn3
use modglobal, only : i1,j1,k1
use modfields, only : rhof, svm, sv0
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


! remove negative values and non physical low values
! also set k_low/k_high to contain only the levels with non-zero svp() values
! --------------------------------------------------
subroutine correct_neg_qt(svp_col,svm_col,thlp_col,qtp_col,k_low,k_high)
  use modglobal, only : cp, rlv, k1
  use modfields, only : exnf
  implicit none
  real, intent(in)    :: svm_col(ncols,k1)
  real, intent(inout) :: svp_col(ncols,k1),thlp_col(k1),qtp_col(k1)
  integer, intent(inout) :: k_low(ncols), k_high(ncols)

  integer :: k
  real :: nrtest, qrtest

  ! correction, after Jerome's implementation in Gales
  real :: qtp_cor(k1), thlp_cor(k1)

  qtp_cor = 0.
  thlp_cor = 0.

  do k=1,k1
    ! == rain ==
    qrtest=svm_col(iq_hr,k)+svp_col(iq_hr,k)*delt
    nrtest=svm_col(in_hr,k)+svp_col(in_hr,k)*delt
    if ((qrtest < q_hr_min) .or. (nrtest < 0.0)) then
      qtp_cor(k) = qtp_cor(k) + qrtest
      thlp_cor(k) = thlp_cor(k) + (rlv/(cp*exnf(k))) * qrtest

      svp_col(iq_hr,k) = - svm_col(iq_hr,k)/delt
      svp_col(in_hr,k) = - svm_col(in_hr,k)/delt
    end if
    if ((svp_col(iq_hr,k).ne.0.).or.(svp_col(in_hr,k).ne.0.)) then
      k_low(iq_hr) = max(k_low(iq_hr),k)
      k_low(in_hr) = max(k_low(in_hr),k)

      k_high(iq_hr) = min(k_high(iq_hr),k)
      k_high(in_hr) = min(k_high(in_hr),k)
    endif

    ! == snow ==
    qrtest=svm_col(iq_hs,k)+svp_col(iq_hs,k)*delt
    nrtest=svm_col(in_hs,k)+svp_col(in_hs,k)*delt
    if ((qrtest .lt. qsnowmin) .or. (nrtest < 0.0) ) then
      qtp_cor(k) = qtp_cor(k) + qrtest
      thlp_cor(k) = thlp_cor(k) - (rlvi/(cp*exnf(k))) * qrtest

      svp_col(iq_hs,k) = - svm_col(iq_hs,k)/delt
      svp_col(in_hs,k) = - svm_col(in_hs,k)/delt
    endif
    if ((svp_col(iq_hs,k).ne.0.).or.(svp_col(in_hs,k).ne.0.)) then
      k_low(iq_hs) = max(k_low(iq_hs),k)
      k_low(in_hs) = max(k_low(in_hs),k)

      k_high(iq_hs) = min(k_high(iq_hs),k)
      k_high(in_hs) = min(k_high(in_hs),k)
    endif

    ! == graupel ==
    qrtest=svm_col(iq_hg,k)+svp_col(iq_hg,k)*delt
    nrtest=svm_col(in_hg,k)+svp_col(in_hg,k)*delt
    if ((qrtest .lt. qgrmin) .or. (nrtest < 0.0) ) then
      qtp_cor(k) = qtp_cor(k) + qrtest
      thlp_cor(k) = thlp_cor(k) - (rlvi/(cp*exnf(k))) * qrtest

      svp_col(iq_hg,k) = - svm_col(iq_hg,k)/delt
      svp_col(in_hg,k) = - svm_col(in_hg,k)/delt
    endif
    if ((svp_col(iq_hg,k).ne.0.).or.(svp_col(in_hg,k).ne.0.)) then
      k_low(iq_hg) = max(k_low(iq_hg),k)
      k_low(in_hg) = max(k_low(in_hg),k)

      k_high(iq_hg) = min(k_high(iq_hg),k)
      k_high(in_hg) = min(k_high(in_hg),k)
    endif

    ! == cloud ice ==
    qrtest=svm_col(iq_ci,k)+svp_col(iq_ci,k)*delt
    nrtest=svm_col(in_ci,k)+svp_col(in_ci,k)*delt
    if ((qrtest .lt. qicemin) .or. (nrtest < 0.0) ) then
      qtp_cor(k)  = qtp_cor(k) + qrtest
      thlp_cor(k) = thlp_cor(k) - (rlvi/(cp*exnf(k))) * qrtest

      svp_col(iq_ci,k) = - svm_col(iq_ci,k)/delt
      svp_col(in_ci,k) = - svm_col(in_ci,k)/delt
    endif
    if ((svp_col(iq_ci,k).ne.0.).or.(svp_col(in_ci,k).ne.0.)) then
      k_low(iq_ci) = max(k_low(iq_ci),k)
      k_low(in_ci) = max(k_low(in_ci),k)

      k_high(iq_ci) = min(k_high(iq_ci),k)
      k_high(in_ci) = min(k_high(in_ci),k)
    endif

    ! == cloud liquid water ==
    qrtest=svm_col(iq_cl,k)+svp_col(iq_cl,k)*delt
    nrtest=svm_col(in_cl,k)+svp_col(in_cl,k)*delt
    if ((qrtest .lt. qcliqmin) .or. (nrtest .le. 0.0) ) then
      svp_col(iq_cl,k) = - svm_col(iq_cl,k)/delt
      svp_col(in_cl,k) = - svm_col(in_cl,k)/delt
    endif
    if ((svp_col(iq_cl,k).ne.0.).or.(svp_col(in_cl,k).ne.0.)) then
      k_low(iq_cl) = max(k_low(iq_cl),k)
      k_low(in_cl) = max(k_low(in_cl),k)

      k_high(iq_cl) = min(k_high(iq_cl),k)
      k_high(in_cl) = min(k_high(in_cl),k)
    endif
  enddo

  ! apply correction
  if (l_corr_neg_qt) then
    qtp_col = qtp_col + qtp_cor / delt
    thlp_col = thlp_col + thlp_cor / delt
  endif

  ! == CCN checking ==
  if(.not.l_c_ccn) then
    do k=1,k1
      nrtest=svm_col(in_cc,k)+svp_col(in_cc,k)*delt
      if (nrtest < 0.0) then
        svp_col(in_cc,k) = - svm_col(in_cc,k)/delt
      end if
    enddo
    if (svp_col(in_cc,k).ne.0.) then
      k_low(in_cc) = max(k_low(in_cc),k)
      k_high(in_cc) = min(k_high(in_cc),k)
    endif
  endif ! not l_c_ccn
end subroutine correct_neg_qt


!*********************************************************************
! Function to calculate coefficient \bar{a}_vent,n
! for ventilation parameter
! specified in Appendix B in S&B
!
!*********************************************************************
real function calc_avent (nn, mu_a, nu_a, a_a,b_a, av)
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


!*********************************************************************
! Function to calculate coefficient \bar{a}_vent,n
! for ventilation parameter
! specified in Appendix B in S&B
!*********************************************************************
real function calc_bvent (nn, mu_a, nu_a, a_a, b_a, beta_a, bv)
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


!*********************************************************************
! Function to calculate coefficient \delta^k_b
!  for collision/collection
! specified in Appendix C in S&B
!*********************************************************************
real function calc_delta_b (kk, mu_b, nu_b, b_b)
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


!*********************************************************************
! Function to calculate coefficient \delta^k_a,b
! for collision/collection
! specified in Appendix C in S&B
!
!*********************************************************************
real function calc_delta_ab (kk, mu_a, nu_a, b_a, mu_b, nu_b, b_b)
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


!*********************************************************************
! Function to calculate coefficient \vartheta^k_a,b
! for collision/collection
! specified in Appendix C in S&B
!
!*********************************************************************
real function calc_th_b (kk, mu_b, nu_b, b_b, beta_b)
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


!*********************************************************************
! Function to calculate coefficient \vartheta^k_a,b
! for collision/collection
! specified in Appendix C in S&B
!
!*********************************************************************
real function calc_th_ab (kk, mu_a, nu_a, b_a, beta_a, mu_b, nu_b, b_b, beta_b)
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


!*********************************************************************
! Function to calculate the constatnt part of the momement
! used in some constants
! specified in Appendix C in S&B
!
!*********************************************************************
real function calc_cons_mmt (kk, mu_a, nu_a)
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


!*********************************************************************
! Function to calculate the constatnt part of the momement
! used in some constants
! specified in Appendix C in S&B
!
!*********************************************************************
real function calc_cons_v (kk, mu_a, nu_a, al_a, be_a)
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
