module modbulkmicro_column
  use modglobal, only : k1
  implicit none
  real ::  
          ,qt0(k1)  &
          ,qvsl(k1) &
          ,w0(k1)

          ,n_cc(k1), n_ccp(k1) &     ! N_{ccn} nr content [ kg^{-1}] of cloud condensation nuclei
          ,n_cl(k1), n_clp(k1) &     ! N_{c,l} nr content [ kg^{-1}] for liquid cloud droplets,
          ,n_ci(k1), n_cip(k1) &     ! N_{c,i} nr content [ kg^{-1}] for ice cloud droplets,
          ,n_hr(k1), n_hrp(k1) &     ! N_{h,r} nr content [ kg^{-1}] for rain
          ,n_hs(k1), n_hsp(k1) &     ! N_{h,s} nr content [ kg^{-1}] for snow
          ,n_hg(k1), n_hgp(k1) &     ! N_{h,g} nr content [ kg^{-1}] for graupel
          ,q_cl(k1), q_clp(k1) &     ! q_{c,l} water content [kg/kg] for liquid cloud droplets,
          ,q_ci(k1), q_cip(k1) &     ! q_{c,i} water content [kg/kg] for ice cloud droplets,
          ,q_hr(k1), q_hrp(k1) &     ! q_{h,r} water content [kg/kg] for rain
          ,q_hs(k1), q_hsp(k1) &     ! q_{h,s} water content [kg/kg] for snow
          ,q_hg(k1), q_hgp(k1)       ! q_{h,g} water content [kg/kg] for graupel

! TODO: who uses these?
! (precep_hr     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of raindrops
! ,precep_ci     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of ice crystals
! ,precep_hs     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of snow
! ,precep_hg     (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< precipitation of graupel
! ,ret_cc        (2-ih:i1+ih,2-jh:j1+jh,k1) &      !< recovery of ccn
! )

contains

subroutine column_processes(i,j)
  implicit none
  integer, intent(in) :: i,j

  ! sedimentation
  ! -----------------------------------------------------------------
  call sedim_rain3
  call sedim_cl3
  call sedim_ice3
  call sedim_snow3
  call sedim_graupel3

  ! recovery of ccn
  !   - to be later replaced based on advance literature
  ! -----------------------------------------------------------------
  call recover_cc

  ! -----------------------------------------------------------------
  ! NEXT:
  ! remove negative values and non physical low values svp
  ! integrate svp = svp + ..._p
  ! prognostic variables
  ! thlp(i,j,k) = thlp(i,j,k)+thlpmcr(i,j,k)
  ! qtp(i,j,k) = qtp(i,j,k)+qtpmcr(i,j,k)
  ! -----------------------------------------------------------------
end subroutine column_processes


!> Cloud nucleation
!! Written to prognostically evaluate the cloud water number content [ kg^{-1}]
!! directly follows Seifert&Beheng scheme
subroutine nucleation3
  use modglobal, only : dzf,k1
  implicit none

  integer :: k
  real  :: coef_ccn, n_act

  ! note that supersaturation is
  ! not always how supersaturated is water vapour,
  ! depending on a flag, it can also include water already in droplets
  real :: ssat_u    ! supersaturation at (...,k+1)
         ,ssat(k1)  !                 at (...,k)
         ,ssat_d    !                 at (...,k-1)
         ,wdssatdz  ! derivation of supersaturation

  real :: ssat(k1)

  dn_cl_nu = 0.0

  coef_ccn  = 1.0/sat_max**kappa_ccn ! 1.0
  ! allows to keep both definitions consistent

  ! calculating supersaturation
  if (l_sb_nuc_sat) then
    ! calculating supersaturation of water vapour only
    ssat = (100./qvsl)*(qt0-q_cl-qvsl)
  else ! l_sb_nuc_sat
    ! ie. cloud liquid water is also included supersaturation
    ssat = (100./qvsl)*(qt0-qvsl)
  endif ! l_sb_nuc_sat

  ! calculating the derivation - second order estimation?
  ! ? add switches for different derivation calculation? - so first the first order only
  ! option A: central differences lax

  ! just an approximation - same as for the second level,
  ! or 0 to prevent condensation there
  wdssatdz(1) = w0(2)*(ssat(2)-ssat(1))/dzf(2)
  do k=2,k1 - 1
      wdssatdz = 0.5*(w0(k+1)+w0(k))*(ssat(k+1) - ssat(k-1))/(dzf(k)+dzf(k-1))
  enddo
  ! first order approximation of the difference
  wdssatdz(k1) = w0(k1)*(ssat(k1) - ssat(k1-1)) / dzf(k1-1)

  do k=1,k1
    ! BUG: original code went out of bounds for k=1 and k=k1
    !      for now, take d ssat / dz = 0 at boundaries
    if (k.eq.1) then
      ssat_d = ssat(1)
    else
      ssat_d = ssat(k-1)
    endif
    if (k.eq.k1) then
      ssat_u = ssat(k1)
    else
      ssat_u = ssat(k+1)
    endif

    if (l_sb_nuc_expl) then ! calculation of explicit nucleation

      if(ssat(k).gt.sat_min) then ! of course only in cloud

        if (l_c_ccn) then ! c_ccn is constant

          if (l_sb_nuc_diff) then

            if (l_sb_sat_max) then

              if ((ssat(k).lt.sat_max).and.(wdssatdz(k).gt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = c_ccn*kappa_ccn*wdssatdz(k)*ssat(k)**(kappa_ccn-1.0) ! (1/rhof) *rhof = 1
              endif

            else ! NOT l_sb_sat_max  AND l_sb_nuc_diff

              if (wdssatdz(k).gt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = c_ccn*kappa_ccn*wdssatdz(k)*ssat(k)**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
              endif

            endif ! l_sb_sat_max

          else ! NOT l_sb_nuc_diff

            if (l_sb_sat_max) then ! l_sb_sat_max AND NOT l_sb_nuc_diff

              if ((ssat(k).lt.sat_max).and.(w0(k).gt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = (c_ccn/dzf(k-1))*max(0.0,w0(k)*      &
                   (ssat(k)**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if ((ssat(k).lt.sat_max).and.(w0(k+1).lt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = (c_ccn/dzf(k-1))*max(0.0,w0(k+1)*    &
                   (ssat(k)**kappa_ccn-ssat_u**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif

            else ! NOT l_sb_sat_max AND NOT l_sb_nuc_diff

              if (w0(k).gt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = (c_ccn/dzf(k-1))*max(0.0,w0(k)*      &
                   (ssat(k)**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if (w0(k+1).lt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = (c_ccn/dzf(k-1))*max(0.0,w0(k+1)*    &
                   (ssat(k)**kappa_ccn-ssat_u**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif

            endif  ! NOT l_sb_sat_max  AND NOT l_sb_nuc_diff

          endif ! NOT l_sb_nuc_diff

          ! basic limiting
          dn_cl_nu(k) = max(0.0,min(dn_cl_nu(k),(n_clmax-n_cl(k))/delt))

        else ! c_ccn is not constant, but dependent on n_cc

          if (l_sb_nuc_diff) then

            if (l_sb_sat_max) then ! l_sb_sat_max

              if ((ssat(k).lt.sat_max).and.(wdssatdz(k).gt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = coef_ccn*n_cc(k)*kappa_ccn*wdssatdz(k)*ssat(k)**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
              endif

            else  ! l_sb_sat_max

              if (wdssatdz(k).gt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = coef_ccn*n_cc(k)*kappa_ccn*wdssatdz(k)*ssat(k)**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
              endif

            endif  ! l_sb_sat_max

          else  ! l_sb_nuc_diff

            if (l_sb_sat_max)  then ! l_sb_sat_max  AND NOT l_sb_nuc_diff

              if ((ssat(k).lt.sat_max).and.(w0(k).gt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = (coef_ccn*n_cc(k)/dzf(k-1))*max(0.0,w0(k)*   &
                  (ssat(k)**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if ((ssat(k).lt.sat_max).and.(w0(k+1).lt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = (coef_ccn*n_cc(k)/dzf(k-1))*max(0.0,w0(k+1)* &
                   (ssat(k)**kappa_ccn-ssat_u**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif

            else ! NOT l_sb_sat_max  AND NOT l_sb_nuc_diff

              if (w0(k).gt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = (coef_ccn*n_cc(k)/dzf(k-1))*max(0.0,w0(k)*   &
                    (ssat(k)**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if (w0(k+1).lt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = (coef_ccn*n_cc(k)/dzf(k-1))*max(0.0,w0(k+1)* &
                   (ssat(k)**kappa_ccn-ssat_u**kappa_ccn))    ! (1/rhof) *rhof = 1
              endif

            endif  ! l_sb_sat_max

          endif  ! l_sb_nuc_diff

          ! basic limiting
          dn_cl_nu(k) = max(0.0,min(dn_cl_nu(k),(n_cc(k)-n_cl(k))/delt))

        endif ! l_c_ccn

      endif ! ssat.gt.sat_min

    else ! l_sb_nuc_expl

      if(ssat.gt.0.0) then

        ! calculate number of activated n_ccn
        n_act = coef_ccn*n_cc(k)*min(sat_max,ssat(k))**kappa_ccn
        n_act = max(n_cc(k),n_act(k))

        if (n_act.gt.n_cl(k)) then
          dn_cl_nu(k) = (n_act-n_cl(k))/delt
        endif

        ! basic limiting - not needed in this case
        ! dn_cl_nu = min(dn_cl_nu, (n_cc - n_cl))/delt)

      endif ! ssat.gt.0.0

    endif ! l_sb_nuc_expl

    !if (l_sb_dbg) then
    !  ! warning if too high liquid water
    !  if (ssat.gt.0.0 .and.( &
    !  (qt0-qvsl-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci))/delt-x_cnuc*dn_cl_nu.lt. 0.)
    !  ) then
    !    write(6,*) 'WARNING: cloud nucleation too high'
    !    write(6,*) ' removing too much water'
    !    ! count(ssat.gt.0.0 .and.( &
    !    ! (qt0-qvsl-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci))/delt-x_cnuc*dn_cl_nu.lt.0.))
    !  end if
    !endif ! l_sb_dbg
  enddo

  ! increase in cloud water number
  n_clp = n_clp + dn_cl_nu

  ! update water density [kg kg^{-1}]
  q_clp = q_clp + x_cnuc * dn_cl_nu

end subroutine  nucleation3


!> Sedimentaion of rain
!! sedimentation of drizzle water
!! - gen. gamma distr is assumed. Terminal velocities param according to
!!   Stevens & Seifert. Flux are calc. anal.
!! - l_lognormal =T : lognormal DSD is assumed with D_g and N known and
!!   sig_g assumed. Flux are calc. numerically with help of a
!!   polynomial function
subroutine sedim_rain3
  use modglobal, only : k1,kmax,eps1,dzf
  use modfields, only : rhof
  implicit none

  integer :: k,jn,n_spl

  real :: wvar        &!< work variable
         ,wvar0       &!< extra test
         ,dt_spl      &!<
         ,wfallmax    &!<
         ,xr_spl      &!< for time splitting
         ,xr_try      &!<
         ,Dvr_spl     &!<     -
         ,mur_spl     &!<     -
         ,lbdr_spl    &!<     -
         ,Dgr         &!< lognormal geometric diameter
         ,N_r0        &!< rain integral stuff
         ,N_r0_try    &!<
         ,lbdr_try    &!< rain integral stuff
         ,pwcont       !<

  real, allocatable :: sed_qr(:), sed_Nr(:)
  allocate(sed_qr(1:k1), sed_Nr(1:k1))

  ! wfallmax = 9.9
  wfallmax = wfallmax_hr
  n_spl = ceiling(split_factor*wfallmax*delt/(minval(dzf)))
  dt_spl = delt/real(n_spl)

  qr_spl(1:k1) = q_hr(1:k1)
  Nr_spl(1:k1) = n_hr(1:k1)

  do jn = 1 , n_spl ! time splitting loop

    sed_qr(1:k1) = 0.
    sed_Nr(1:k1) = 0.

    if (l_sb) then
      if (l_sb_classic) then
        do k=1,k1
          if (qr_spl(k) > qrmin) then
            ! limiting procedure (as per S&B)
            xr_spl     = qr_spl(k)/(Nr_spl(k)+eps0)
            xr_try     = max(xrmin,min(xrmax,x_hr(k))) ! BUG: x_hr should be xr_spl?
            Dvr_spl    = (xr_try/pirhow)**(1./3.)
            N_r0_try   = n_hr(k)/Dvr(k) ! rhof(k)*n_hr(k)/Dvr(k)
            N_r0       = max(N_0min/rhof(k),min(N_0max/rhof(k),N_r0_try))

            lbdr_try  = (pirhow*N_r0/(rhof(k)*qr_spl(k)))**0.25
            lbdr_spl  = max(lbdr_min, min(lbdr_max,lbdr_try))

            ! calculation of velocities
            wfall_qr        = max(0.,((rho0s/rhof(k))**0.5)*(a_tvsbc  &
                      -b_tvsbc*(1.+c_tvsbc/lbdr_spl)**(-4.0))) ! k=1
            wfall_Nr        = max(0.,((rho0s/rhof(k))**0.5)*(a_tvsbc  &
                      -b_tvsbc*(1.+c_tvsbc/lbdr_spl)**(-1.0))) ! k=0
            sed_qr  (k) = wfall_qr       *qr_spl(k)*rhof(k)
            sed_Nr  (k) = wfall_Nr       *Nr_spl(k)*rhof(k)
          endif
        enddo
      else  ! l_sb_classic
        if (l_lognormal) then
          do k=1, kmax ! BUG: kmax should be k1?
            if (qr_spl(k) > qrmin) then
              ! correction for width of DSD
              ! BUG: Dvr_spl unset!
              Dgr = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr_spl

              sed_qr(k) = 1.*sed_flux3(Nr_spl(k),Dgr,log(sig_gr)**2,D_s,3)
              sed_Nr(k) = 1./pirhow*sed_flux3(Nr_spl(k),Dgr,log(sig_gr)**2,D_s,0)

              ! correction for the fact that pwcont .ne. qr_spl
              ! actually in this way for every grid box a fall velocity is determined
              pwcont = liq_cont3(Nr_spl(k),Dgr,log(sig_gr)**2,D_s,3)         ! note : kg m-3
              if (pwcont > eps1) then
                sed_qr(k) = (qr_spl(k)*rhof(k)/pwcont)*sed_qr(k)
              end if
            end if ! qr_spl
          enddo
        else ! l_lognormal
          !
          ! SB rain sedimentation
          !
          do k=1,k1
            if (l_mur_cst) then
              mur_spl = mur_cst
            else
              if (qr_spl(k) > qrmin) then
                ! SS08
                ! BUG: Dvr_spl is unset
                ! mur_spl = 10. * (1+tanh(1200.*(Dvr_spl(k)-0.0014)))

                ! G09b
                mur_spl = min(30.,- 1. + 0.008/ (qr_spl(k)*rhof(k))**0.6)
              endif
            endif ! l_mur_cst

            if (qr_spl(k) > qrmin) then
              lbdr_spl  = ((mur_spl+3.)*(mur_spl+2.)*(mur_spl+1.))**(1./3.)/Dvr_spl ! BUG: Dvr_spl is unset
              wfall_qr  = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl)**(-1.*(mur_spl+4.))))
              wfall_Nr  = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl)**(-1.*(mur_spl+1.))))
              sed_qr(k) = wfall_qr * qr_spl(k) * rhof(k)
              sed_Nr(k) = wfall_Nr * Nr_spl(k) * rhof(k)
            endif
          enddo
        endif ! l_lognormal
      endif ! l_sb_classic
    else ! l_sb
      !
      ! KK00 rain sedimentation
      !
      do k=1,k1
        if (qr_spl(k) > qrmin) then
          !JvdD added eps0 to avoid division by zero
          xr_spl = rhof(k)*qr_spl(k)/(Nr_spl(k)+eps0)

          ! to ensure xr is within borders
          xr_spl = min(xr_spl,xrmaxkk)

          Dvr_spl = (xr_spl/pirhow)**(1./3.)
          sed_qr(k) = max(0., 0.006*1.0E6*Dvr_spl - 0.2) * qr_spl(k)*rhof(k)
          sed_Nr(k) = max(0.,0.0035*1.0E6*Dvr_spl - 0.1) * Nr_spl(k)
        endif
      enddo
    end if ! l_sb

    do k = 1,kmax ! second k loop
      wvar  = qr_spl(k) + (sed_qr(k+1) - sed_qr(k))*dt_spl/(dzf(k)*rhof(k))
      wvar0 = Nr_spl(k) + (sed_Nr(k+1) - sed_Nr(k))*dt_spl/(dzf(k)*rhof(k))
      if (wvar.lt.0.) then
        write(6,*)'  rain sedim of q_r too large'
      end if
      if (wvar0.lt.0.) then
        write(6,*)'  rain sedim of N_r too large'
      end if

      Nr_spl(k) = Nr_spl(k) + (sed_Nr(k+1) - sed_Nr(k))*dt_spl/(dzf(k)*rhof(k))
      qr_spl(k) = qr_spl(k) + (sed_qr(k+1) - sed_qr(k))*dt_spl/(dzf(k)*rhof(k))

      if (jn == 1.) then
        precep_hr(k) = sed_qr(k)/rhof(k) ! kg kg-1 m s-1
      endif
    enddo  ! second k loop
  enddo ! time splitting loop

  do k=1,k1
    ! tendencies
    dn_hr_se(k) = (Nr_spl(k) - n_hr(k))/delt
    dq_hr_se(k) = (qr_spl(k) - q_hr(k))/delt

    ! updates
    n_hrp(k) = n_hrp(k) + dn_hr_se(k)
    q_hrp(k) = q_hrp(k) + dq_hr_se(k)
  enddo
end subroutine sedim_rain3
  

! sedimentation of snow
! ---------------------
subroutine sedim_snow3
  use modglobal, only : k1,kmax,dzf
  use modfields, only : rhof
  implicit none

  integer :: k,jn
  integer :: n_spl      !<  sedimentation time splitting loop
  real    :: pwcont, xpmin, xpmax, qip_min                       &
            ,c_v_0, c_v_1, be_ip, aip, bip

  real, allocatable,dimension(:)  :: qip_spl, nip_spl
  real, allocatable,dimension(:)  :: sed_qip, sed_nip

  real :: wfall_qip, wfall_nip, xip_spl, wvar
  real :: dt_spl,wfallmax

  ! set constants 
  xpmin = x_hs_bmin
  xpmax = x_hs_bmax
  qip_min = qsnowmin
  c_v_0 = c_v_s0
  c_v_1 = c_v_s1
  be_ip = be_hs
  aip   = a_hs
  bip   = b_hs

  ! wfallmax = 9.9   ! <- replace with a highest terminal velocity for particles
  wfallmax = wfallmax_hs

  ! allocate 
  allocate( sed_qip(k1)    &
           ,sed_nip(k1)    &
           ,qip_spl(k1)    &
           ,nip_spl(k1)    )             
     
  sed_qip = 0.0
  sed_nip = 0.0  
  wfall_qip = 0.0
  wfall_nip = 0.0

  qip_spl(1:k1)  = q_hs(1:k1)
  nip_spl(1:k1)  = n_hs(1:k1)

  n_spl = ceiling(split_factor*wfallmax*delt/(minval(dzf)))
  dt_spl = delt/real(n_spl)

  do jn = 1 , n_spl ! time splitting loop
    sed_qip = 0.
    sed_nip = 0.

    do k=1,k1
      ! terminal fall velocity
      if ((qip_spl(k) > qip_min).and.(nip_spl(k) > 0.0)) then
        xip_spl = qip_spl(k)/(nip_spl(k)+eps0) ! JvdD Added eps0 to avoid division by zero
        xip_spl = min(max(xip_spl,xpmin),xpmax) ! to ensure xr is within borders
        ! Dvp_spl = aip*xip_spl**bip          
        wfall_qip = max(0.0, c_v_1 * xip_spl**be_ip)
        wfall_nip = max(0.0, c_v_0 * xip_spl**be_ip)
        sed_qip(k) = wfall_qip*qip_spl(k)*rhof(k)
        sed_nip(k) = wfall_nip*nip_spl(k)*rhof(k)
      endif
    enddo

    ! segmentation over levels
    do k = 1,kmax

      wvar = qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k))
      if (wvar.lt. 0.) then
        write(6,*)'  snow sedim too large'
      end if

      wvar = nip_spl(k) + (sed_nip(k+1) - sed_nip(k))*dt_spl/(dzf(k)*rhof(k))
      if (wvar.lt. 0.) then
        write(6,*)'  snow sedim too large'
      end if    

      nip_spl(k) = nip_spl(k) + (sed_nip(k+1) - sed_nip(k))*dt_spl/(dzf(k)*rhof(k))
      qip_spl(k) = qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k))

      ! -> check this part properly later
      if ( jn == 1. ) then
        precep_hs(k) = sed_qip(k)/rhof(k)       ! kg kg-1 m s-1
        precep_i(k) = precep_i(k)+ precep_hs(k) ! kg kg-1 m s-1       
      endif
    enddo  ! second k loop
  enddo ! time splitting loop

  do k=1,k1      
    ! tendencies
    dn_hs_se(k) = (nip_spl(k) - n_hs(k))/delt
    dq_hs_se(k) = (qip_spl(k) - q_hs(k))/delt

    ! updates 
    n_hsp(k)= n_hsp(k) + dn_hs_se(k) 
    q_hsp(k)= q_hsp(k) + dq_hs_se(k)
  enddo       

  deallocate (nip_spl, qip_spl,sed_nip, sed_qip)
end subroutine sedim_snow3
  
  
! sedimentation of graupel
! ------------------------
subroutine sedim_graupel3
  use modglobal, only : k1,kmax,dzf
  use modfields, only : rhof
  implicit none

  integer :: k,jn
  integer :: n_spl      !<  sedimentation time splitting loop
  real    :: pwcont, xpmin, xpmax, c_v_0, c_v_1, be_ip, aip, bip
  real    :: qip_min

  real, allocatable,dimension(:)  :: qip_spl, nip_spl
  real, allocatable,dimension(:)  :: sed_qip, sed_nip

  real :: wvar,xip_spl,Dvp_spl,mur_spl

  real :: dt_spl,wfallmax, wfall_nip, wfall_qip

  allocate( sed_qip(k1) &
           ,sed_nip(k1) &
           ,qip_spl(k1) &
           ,nip_spl(k1) )

  sed_qip = 0.0
  sed_nip = 0.0  
  wfall_qip = 0.0
  wfall_nip = 0.0

  qip_spl(1:k1)  = q_hg(1:k1)
  nip_spl(1:k1)  = n_hg(1:k1)

  ! set constants 
  xpmin = x_hg_bmin
  xpmax = x_hg_bmax
  c_v_0 = c_v_g0
  c_v_1 = c_v_g1
  be_ip = be_hs ! be_hg
  aip   = a_hg
  bip   = b_hg
  qip_min= qgrmin

  ! wfallmax = 15.9 ! 9.9   ! <- replace with a highest terminal velocity for particles
  wfallmax = wfallmax_hg

  n_spl = ceiling(split_factor*wfallmax*delt/(minval(dzf)))
  dt_spl = delt/real(n_spl)

  do jn = 1 , n_spl ! time splitting loop
    sed_qip = 0.
    sed_nip = 0.

    do k=1,k1
      if ((qip_spl(k) > qip_min).and.(nip_spl(k) > 0.0) ) then
        xip_spl = qip_spl(k)/(nip_spl(k)+eps0) ! JvdD Added eps0 to avoid division by zero
        xip_spl = min(max(xip_spl,xpmin),xpmax) ! to ensure xr is within borders
        ! Dvp_spl = aip*xip_spl**bip          

        wfall_qip(k) = max(0.0,c_v_1 * xip_spl**be_ip)
        wfall_nip(k) = max(0.0,c_v_0 * xip_spl**be_ip)
        sed_qip(k)   = wfall_qip(k)*qip_spl(k)*rhof(k)
        sed_nip(k)   = wfall_nip(k)*nip_spl(k)*rhof(k)
      endif
    enddo
    enddo
    enddo

    ! segmentation over levels
    do k = 1,kmax
      wvar = qip_spl(k) +  (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k))
      if (wvar.lt.0.) then
        write(6,*)'  graupel sedim too large'
      end if

      wvar(k) = nip_spl(k) + (sed_nip(k+1) - sed_nip(k))*dt_spl/(dzf(k)*rhof(k))
      if (wvar.lt.0.) then
        write(6,*)'  graupel sedim too large'
      end if

      nip_spl(k) = nip_spl(k) + (sed_nip(k+1) - sed_nip(k))*dt_spl/(dzf(k)*rhof(k))
      qip_spl(k) = qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k))

      ! -> check this part properly later
      if ( jn == 1. ) then
        precep_hg(k) = precep_hg(k) + sed_qip(k) /rhof(k) ! kg kg-1 m s-1
        precep_i(k) = precep_i(k) + precep_hg(k) ! kg kg-1 m s-1
      endif
    enddo  ! second k loop
  enddo ! time splitting loop


  do k=1,k1   
    ! tendencies
    dn_hg_se(k) = (nip_spl(k) - n_hg(k))/delt
    dq_hg_se(k) = (qip_spl(k) - q_hg(k))/delt

    ! updates 
    n_hgp(k)= n_hgp(k) + dn_hg_se(k) 
    q_hgp(k)= q_hgp(k) + dq_hg_se(k)
  enddo 

  deallocate (nip_spl, qip_spl,sed_nip,sed_qip)
end subroutine sedim_graupel3
 

! sedimentation of cloud ice 
! --------------------------
subroutine sedim_ice3
  use modglobal, only : k1,kmax,dzf
  use modfields, only : rhof
  implicit none

  integer :: k,jn,n_spl

  real, allocatable,dimension(:)  :: qip_spl, nip_spl
  real, allocatable,dimension(:)  :: sed_qip, sed_nip

  real :: dt_spl, xip_spl, wvar, wfall

  allocate( sed_qip(k1) &
           ,sed_nip(k1) &
           ,qip_spl(k1) &
           ,nip_spl(k1) )
   

  n_spl = ceiling(split_factor*wfallmax_ci*delt/(minval(dzf)))
  dt_spl = delt/real(n_spl)

  qip_spl(1:k1)  = q_ci(1:k1)
  nip_spl(1:k1)  = n_ci(1:k1)

  do jn = 1 , n_spl ! time splitting loop
    sed_qip(1:k1) = 0.
    sed_nip(1:k1) = 0.

    do k=1,k1
      if ( (qip_spl(k) > qicemin).and.(nip_spl(k) > 0.0) ) then
        xip_spl = qip_spl(k)/(nip_spl(k)+eps0) ! JvdD Added eps0 to avoid division by zero
        xip_spl = min(max(xip_spl,x_ci_bmin),x_ci_bmax) ! to ensure xr is within borders
        ! Dvp_spl = a_ci*xip_spl**b_ci         

        ! terminal fall velocity
        wfall = max(0.0,c_v_i * xip_spl**be_ci)
        sed_qip(k) = wfall*qip_spl(k)*rhof(k)

        wfall = max(0.0,c_v_i0 * xip_spl**be_ci)
        sed_nip(k) = wfall*nip_spl(k)*rhof(k)
      endif
    enddo

    ! segmentation over levels
    do k = 1,kmax
      wvar = qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k))
      if (wvar.lt. 0.) then
        write(6,*) 'ice sedim too large'
      end if

      nip_spl(k) = nip_spl(k) + (sed_nip(k+1) - sed_nip(k))*dt_spl/(dzf(k)*rhof(k))
      qip_spl(k) = qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k))

      !d -> check this part properly later
      if ( jn == 1. ) then
        precep_ci(k) = sed_qip(k)/rhof(k) ! kg kg-1 m s-1
        precep_i(k) = precep_i(k) + precep_ci(k) ! kg kg-1 m s-1
      endif

    enddo  ! second k loop
  enddo ! time splitting loop

  do k=1,k1
    ! tendencies
    dn_ci_se(k) = (nip_spl(k) - n_ci(k))/delt
    dq_ci_se(k) = (qip_spl(k) - q_ci(k))/delt

    ! updates 
    n_cip(k)= n_cip(k) + dn_ci_se(k) 
    q_cip(k)= q_cip(k) + dq_ci_se(k)

    ! BUG: also qtpmcr and thlpmcr change?
    ! qtpmcr(k) = qtpmcr + 0.0
    ! thlpmcr(k) = thlpmcr + 0.0
  enddo

  deallocate(nip_spl,qip_spl,sed_nip,sed_qip)
end subroutine sedim_ice3


!*********************************************************************
! sedimentation of cloud water
!*********************************************************************
subroutine sedim_cl3 ! sedim_ice3
  use modglobal, only : k1,kmax,dzf,rlv,cp
  use modfields, only : rhof, exnf
  implicit none

  integer :: k, jn, n_spl
  real, allocatable,dimension(:,:,:)  :: qip_spl, nip_spl
  real, allocatable,dimension(:,:,:)  :: sed_qip, sed_nip
  real, allocatable :: wvar, xip_spl
  real,save :: dt_spl,wfall


  allocate( sed_qip(k1)   &
           ,sed_nip(k1)   &
           ,qip_spl(k1)   &
           ,nip_spl(k1)   &  
   
  qip_spl(1:k1)  = q_cl(1:k1)
  nip_spl(1:k1)  = n_cl(1:k1)

  n_spl = ceiling(split_factor*wfallmax_cl*delt/(minval(dzf)))
  dt_spl = delt/real(n_spl)

  do jn = 1 , n_spl ! time splitting loop

    sed_qip(1:k1) = 0.
    sed_nip(1:k1) = 0.

    do k=1,k1
      if ((qip_spl(k) > qcliqmin).and.(nip_spl(k) > 0.0)) then
        xip_spl = qip_spl(k)/(nip_spl(k)+eps0) ! JvdD Added eps0 to avoid division by zero
        xip_spl = min(max(xip_spl,x_cl_bmin),x_cl_bmax) ! to ensure xr is within borders
        ! Dvp_spl = a_cl*xip_spl**b_cl         

        ! terminal fall velocity
        wfall = max(0.0,c_v_c_1 * xip_spl**be_cl)
        sed_qip(k)   = wfall*qip_spl(k)*rhof(k)

        wfall = max(0.0,c_v_c_0 * xip_spl**be_cl)
        sed_nip(k) = wfall*nip_spl(k)*rhof(k)
      endif
    enddo

    ! segmentation over levels
    do k = 1,kmax
      wvar = qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k))
      if (wvar.lt. 0.) then
        write(6,*)'cloud sedim too large'
      end if
      nip_spl(k) = nip_spl(k) + (sed_nip(k+1) - sed_nip(k))*dt_spl/(dzf(k)*rhof(k))
      qip_spl(k) = qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k))

      ! BUG: no precep?
    enddo  ! second k loop

  enddo ! time splitting loop

  do k=1,k1
    ! tendencies
    dn_cl_se(k) = (nip_spl(k) - n_cl(k))/delt
    dq_cl_se(k) = (qip_spl(k) - q_cl(k))/delt

    ! updates
    n_clp(k) = n_clp(k) + dn_cl_se(k) 
    q_clp(k) = q_clp(k) + dq_cl_se(k)

    ! also qtpmcr and thlpmcr change      
    qtpmcr(k) = qtpmcr(k) + dq_cl_se(k)
    thlpmcr(k) = thlpmcr(k)-(rlv/(cp*exnf(k)))*dq_cl_se(k)
  enddo         

  deallocate (nip_spl,qip_spl,sed_nip, sed_qip)
end subroutine sedim_cl3


!    recovery of ccn
!
!   - to be later replaced based on advance literature
!   - so far just and easy recovery
!     of ccn based on number of water particles that evaporated, sublimated
!     or got removed with remaining positive n_
! -------------------------------------------------------------------------
subroutine recover_cc
  use modglobal, only : k1
  implicit none

  integer :: k

  if(.not.(l_c_ccn)) then
    do k=1,k1
        ! decrease in total amount of potential CCN
        n_ccp = n_ccp         &
              + dn_cl_sc      &
              + dn_cl_se      &
              + dn_cl_au      &
              + dn_cl_ac      &
              + dn_cl_hom     &
              + dn_cl_het     &
              + dn_cl_rime_ci &
              + dn_cl_rime_hs &
              + dn_cl_rime_hg

        n_ccp = n_ccp + c_rec_cc*ret_cc
    enddo
  endif
end subroutine recover_cc


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
! ---------------------------------------------------------------------
real function sed_flux3(Nin,Din,sig2,Ddiv,nnn)
  use modglobal, only : pi,rhow
  implicit none

  real, intent(in)    :: Nin, Din, sig2, Ddiv
  integer, intent(in) :: nnn

  ! para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
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

  if (Din < Ddiv) then
    alfa = 3.e5*100  ![1/ms]
    beta = 2
    D_min = D_intmin
    D_max = Ddiv
    flux = C*Nin*alfa*erfint3(beta,Din,D_min,D_max,sig2,nnn)
  else
    ! fall speed ~ D^2
    alfa = 3.e5*100 ![1/m 1/s]
    beta = 2
    D_min = Ddiv
    D_max = 133e-6
    flux = flux + C*Nin*alfa*erfint3(beta,Din,D_min,D_max,sig2,nnn)

    ! fall speed ~ D
    alfa = 4e3     ![1/s]
    beta = 1
    D_min = 133e-6
    D_max = 1.25e-3
    flux = flux + C*Nin*alfa*erfint3(beta,Din,D_min,D_max,sig2,nnn)

    ! fall speed ~ sqrt(D)
    alfa = 1.4e3 *0.1  ![m^.5 1/s]
    beta = .5
    D_min = 1.25e-3
    D_max = D_intmax
    flux = flux + C*Nin*alfa*erfint3(beta,Din,D_min,D_max,sig2,nnn)
    end do
  end if
  sed_flux3 = flux
end function sed_flux3


! Function to calculate numerically the analytical solution of the
! liq. water content between Dmin and Dmax based on
! Feingold et al 1986 eq 17 -20.
!
! M.C. van Zanten    September 2005
! -----------------------------------------------------------------
real function liq_cont3(Nin,Din,sig2,Ddiv,nnn)
  use modglobal, only : pi,rhow
  implicit none

  real, intent(in)    :: Nin, Din, sig2, Ddiv
  integer, intent(in) :: nnn

  ! para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
  ! ,power of of D in integral
  real, parameter :: beta = 0           &
                    ,C = pi/6.*rhow     &
                    ,D_intmin = 80e-6   &   ! value of start of rain D
                    ,D_intmax = 3e-3        !4.3e-3    ! value is now max value for sqrt fall speed rel.

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


! Function to calculate erf(x) approximated by a polynomial as
! specified in 7.1.27 in Abramowitz and Stegun
! NB phi(x) = 0.5(erf(0.707107*x)+1) but 1 disappears by substraction
! -------------------------------------------------------------------
real function erfint3(beta, D, D_min, D_max, sig2,nnn )
  implicit none
  real, intent(in)    :: beta, D, D_min, D_max, sig2
  integer, intent(in) :: nnn

  real, parameter :: eps = 1e-10      &
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
  if (erfint3 < 0.) erfint3 = 0.
end function erfint3

end module modbulkmicro_column
