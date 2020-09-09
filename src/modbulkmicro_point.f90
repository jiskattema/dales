module modbulkmicro_point
  use modmicrodata3, only : in_hr,iq_hr,in_cl,in_tr1,iq_cl,in_cc, &
                            in_ci,iq_ci,in_hs,iq_hs,in_hg,iq_hg

  implicit none

  logical :: q_hr_mask,q_hs_mask,q_hg_mask

  real    :: delt

  real ::  qltot    &
          ,n_cc     & ! N_{ccn} number content [ kg^{-1}] of cloud condensation nuclei
          ,n_ccp    &
          ,n_cl     & ! N_{c,l}  number content [ kg^{-1}] for liquid cloud droplets,
          ,n_clp    &
          ,n_ci     & ! N_{c,i}  number content [ kg^{-1}] for ice cloud droplets,
          ,n_cip    &
          ,n_hr     & ! N_{h,r} number content [ kg^{-1}] for rain
          ,n_hrp    &
          ,n_hs     & ! N_{h,s} number content [ kg^{-1}] for snow
          ,n_hsp    &
          ,n_hg     & ! N_{h,g} number content [ kg^{-1}] for graupel
          ,n_hgp    &
          ,q_cl     & ! q_{c,l}  water content [kg kg^{-1}] for liquid cloud droplets,
          ,q_clp    &
          ,q_ci     & ! q_{c,i}  water content [kg kg^{-1}] for ice cloud droplets,
          ,q_cip    &
          ,q_hr     & ! q_{h,r} water content [kg kg^{-1}] for rain
          ,q_hrp    &
          ,q_hs     & ! q_{h,s} water content [kg kg^{-1}] for snow
          ,q_hsp    &
          ,q_hg     & ! q_{h,g} water content [kg kg^{-1}] for graupel
          ,q_hgp    &
          ,x_ci     & ! mean cloud ice size
          ,x_cl     & ! mean cloud water size
          ,x_hs     & ! mean snow size
          ,x_hg     & ! mean graupel size
          ,x_hr     & ! mean raindrop size
          ,D_ci     & ! \bar{D}_ci mean diameter for cloud ice particle
          ,D_cl     & ! \bar{D}_cl mean diameter for cloud water particle
          ,D_hr     & ! \bar{D}_hr mean diameter for raindrops
          ,D_hs     & ! \bar{D}_hs mean diameter for snow particle
          ,D_hg     & ! \bar{D}_hg mean diameter for graupel particle
          ,v_ci     & ! \bar{v}_ci mean velocity for cloud ice particle
          ,v_cl     & ! \bar{v}_cl mean velocity for cloud water droplets
          ,v_hr     & ! \bar{v}_hr mean velocity for raindrops
          ,v_hs     & ! \bar{v}_hs mean velocity for snow particle
          ,v_hg     & ! \bar{v}_hg mean velocity for graupel particle
          ,xc       &
          ,Dvr      &
          ,lbdr

  real ::     dn_cl_nu       &    !< droplet nucleation rate
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
             ,dq_cl_sa       &      !< saturation adjustment
             ,dn_cl_sa              !< change in n_cl due to saturation adjustment

contains

! calculating rain and other integrals
subroutine integrals_bulk3
  use modglobal, only : i1,j1,k1,kmax
  use modfields, only : rhof
  implicit none

  integer :: i,j,k
  real :: xr, N_r0, xr_try, libdr_try, N_r0_try, mur

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (l_rain) then
      Dvr  = 0.
      mur  = 30.
      lbdr = 0.

      if (l_sb) then
        if (q_hr_mask) then
          if(l_sb_classic) then
            ! limiting procedure (as per S&B)
            x_hr    = q_hr/(n_hr+eps0) ! JvdD Added eps0 to avoid floating point exception
            xr_try  = max(xrmin,min(xrmax,x_hr))  !
            x_hr    = max(xrmin,min(xrmax,x_hr))  ! same as above, for cold mphys

            D_hr = a_hr *x_hr**b_hr
            v_hr = al_hr*x_hr**be_hr*(rho0s/rhof(k))**ga_hr

            Dvr      = (xr_try/pirhow)**(1./3.)
            N_r0_try = rhof(k)*n_hr/Dvr
            N_r0     = max(N_0min,min(N_0max,N_r0_try))

            lbdr_try = (pirhow*N_r0/(rhof(k)*q_hr))**0.25  ! c_lbdr*x_hr**(-mu_hr_cst)
            lbdr     = max(lbdr_min, min(lbdr_max,lbdr_try))
          else ! l_sb_classic
            ! #sb3 changing the variable name
            x_hr = q_hr/(n_hr+eps0)
            ! TODO: limit x_hr?
            xr = min(max(x_hr,xrmin),xrmax) ! to ensure xr is within borders
            Dvr = (xr/pirhow)**(1./3.)
            !
            D_hr = a_hr *x_hr**b_hr
            v_hr = al_hr*x_hr**be_hr*(rho0s/rhof(k))**ga_hr

            if (l_mur_cst) then
              ! mur = cst
              mur  = mur_cst
              lbdr = ((mur+3.)*(mur+2.)*(mur+1.))**(1./3.)/Dvr
            else
              ! mur = f(Dv)
              mur  = min(mur0_G09b,- 1. + c_G09b/ (q_hr*rhof(k))**exp_G09b)  ! G09b
              lbdr = ((mur+3.)*(mur+2.)*(mur+1.))**(1./3.)/Dvr
            endif ! l_mur_cst
          endif !  l_sb_classic
        endif ! q_hr_mask
      else !  l_sb
        if (q_hr_mask.and.n_hr.gt.0.) then
          x_hr = q_hr/(n_hr+eps0)
          ! TODO: limit x_hr?
          Dvr = (x_hr/pirhow)**(1./3.)

          D_hr = a_hr *x_hr**b_hr
          v_hr = al_hr*x_hr**be_hr*(rho0s/rhof(k))**ga_hr
        endif
      endif ! l_sb
    endif   ! l_rain

    ! cloud water
    if (q_cl_mask) then
      x_cl = q_cl/(n_cl+eps0)
      ! xc = q_cl/(n_cl+eps0)
      xc   = min(max(x_cl,x_cl_bmin),x_cl_bmax) ! to ensure x is within borders
      x_cl = min(max(x_cl,x_cl_bmin),x_cl_bmax) ! as well - to limit extrme autoconversion

      D_cl = a_cl *x_cl**b_cl
      v_cl = al_cl*x_cl**be_cl*(rho0s/rhof(k))**ga_cl
    endif

    ! cloud ice
    if (q_ci_mask) then
      x_ci = q_ci/(n_ci+eps0)
      x_ci = min(max(x_ci,x_ci_bmin),x_ci_bmax) ! to ensure x is within borders

      D_ci = a_ci *x_ci**b_ci
      v_ci = al_ci*x_ci**be_ci*(rho0s/rhof(k))**ga_ci
    endif

    ! snow
    if (q_hs_mask) then
      x_hs = q_hs/(n_hs+eps0)
      x_hs = min(max(x_hs,x_hs_bmin),x_hs_bmax) ! to ensure x is within borders

      D_hs = a_hs *x_hs**b_hs
      v_hs = al_hs*x_hs**be_hs*(rho0s/rhof(k))**ga_hs
    endif

    ! graupel
    if (q_hg_mask) then
      x_hg = q_hg/(n_hg+eps0)
      x_hg = min(max(x_hg,x_hg_bmin),x_hg_bmax) ! to ensure x is within borders

      D_hg = a_hg *x_hg**b_hg
      v_hg = al_hg*x_hg**be_hg*(rho0s/rhof(k))**ga_hg
    endif
  enddo
  enddo
  enddo
end subroutine integrals_bulk3

!> Cloud nucleation
!! Written to prognostically evaluate the cloud water number content [ kg^{-1}]
!! directly follows Seifert&Beheng scheme
subroutine nucleation3
  use modglobal, only : dzf,i1,j1,k1
  use modfields, only : w0,rhof,qvsl,qt0,svm
  implicit none

  integer :: i,j,k
  real  :: coef_ccn, n_act

  ! note that supersaturation is
  ! not always how supersaturated is water vapour,
  ! depending on a flag, it can also include water already in droplets
  real :: ssat_u    ! supersaturation at (...,k+1)
         ,ssat      !                 at (...,k)
         ,ssat_d    !                 at (...,k-1)
         ,wdssatdz  ! derivation of supersaturation

  dn_cl_nu = 0.0

  coef_ccn  = 1.0/sat_max**kappa_ccn ! 1.0
  ! allows to keep both definitions consistent

  do k=1,k1
  do j=2,j1
  do i=2,i1
    ! calculating supersaturation
    ! calculating the derivation - second order estimation?
    ! ? add switches for different derivation calculation? - so first the first order only
    ! option A: central differences lax
    if (l_sb_nuc_sat) then
      ! calculating supersaturation of water vapour only
      ! NOTE:  we'll need q_cl at points above and below our gridpoint
      !        as q_cl = sv0(i,j,k,iq_cl), we use sv0 for those points instead.
      if (k>1) then
        ssat_d =(100./qvsl(i,j,k-1))*(qt0(i,j,k-1)-sv0(i,j,k-1,iq_cl)-qvsl(i,j,k-1))
      endif
      ssat =(100./qvsl(i,j,k))*(qt0(i,j,k)-q_cl-qvsl(i,j,k))
      if (k<k1) then
        ssat_u =(100./qvsl(i,j,k+1))*(qt0(i,j,k+1)-sv0(i,j,k+1,iq_cl)-qvsl(i,j,k+1))
      else
        ssat_u == ssat
      endif
    else ! l_sb_nuc_sat
      ! ie. cloud liquid water is also included supersaturation
      if (k>1) then
        ssat_d = (100./qvsl(i,j,k-1))*(qt0(i,j,k-1)-qvsl(i,j,k-1))
      endif
      ssat = (100./qvsl(i,j,k))*(qt0(i,j,k)-qvsl(i,j,k))
      if (k<k1) then
        ssat_u = (100./qvsl(i,j,k+1))*(qt0(i,j,k+1)-qvsl(i,j,k+1))
      endif
    endif ! l_sb_nuc_sat

    if (ssat.gt.sat_min) then ! only used above this threshold
      if (k.eq.1) then
        ! just an approximation - same as for the second level, or 0 to prevent condensation there
        wdssatdz = w0(i,j,2)*(ssat_u-ssat)/dzf(2)
        ssat_d = ssat ! BUG: div0? original code went out of array
      elseif (k.eq.k1) then
        ! first order approximation of the difference
        wdssatdz = w0(i,j,k1)*(ssat - ssat_d) / dzf(k1-1)
        ssat_u = ssat + wdssatdz * dzf(k1) / w0(i,j,k1) ! BUG: div0? original code went out of array
      else
        wdssatdz = 0.5*(w0(i,j,k+1)+ w0(i,j,k))*(ssat_u - ssat_d)/(dzf(k)+dzf(k-1))
      endif
    endif

    if (l_sb_nuc_expl) then ! calculation of explicit nucleation

      if(ssat.gt.sat_min) then ! of course only in cloud

        if (l_c_ccn) then ! c_ccn is constant

          if (l_sb_nuc_diff) then

            if (l_sb_sat_max) then

              if ((ssat.lt.sat_max).and.(wdssatdz.gt.0.0)) then ! condition for nucleation
                dn_cl_nu = c_ccn*kappa_ccn*wdssatdz*ssat**(kappa_ccn-1.0) ! (1/rhof) *rhof = 1
              endif

            else ! NOT l_sb_sat_max  AND l_sb_nuc_diff

              if (wdssatdz.gt.0.0) then ! condition for nucleation
                dn_cl_nu = c_ccn*kappa_ccn*wdssatdz*ssat**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
              endif

            endif ! l_sb_sat_max

          else ! NOT l_sb_nuc_diff

            if (l_sb_sat_max) then ! l_sb_sat_max AND NOT l_sb_nuc_diff

              if ((ssat.lt.sat_max).and.(w0(i,j,k).gt.0.0)) then ! condition for nucleation
                dn_cl_nu = (c_ccn/dzf(k-1))*max(0.0,w0(i,j,k)*      &
                   (ssat**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if ((ssat.lt.sat_max).and.(w0(i,j,k+1).lt.0.0)) then ! condition for nucleation
                dn_cl_nu = (c_ccn/dzf(k-1))*max(0.0,w0(i,j,k+1)*    &
                   (ssat**kappa_ccn-ssat_u**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif

            else ! NOT l_sb_sat_max AND NOT l_sb_nuc_diff

              if (w0(i,j,k).gt.0.0) then ! condition for nucleation
                dn_cl_nu = (c_ccn/dzf(k-1))*max(0.0,w0(i,j,k)*      &
                   (ssat**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if (w0(i,j,k+1).lt.0.0) then ! condition for nucleation
                dn_cl_nu = (c_ccn/dzf(k-1))*max(0.0,w0(i,j,k+1)*    &
                   (ssat**kappa_ccn-ssat_u**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif

            endif  ! NOT l_sb_sat_max  AND NOT l_sb_nuc_diff

          endif ! NOT l_sb_nuc_diff

          ! basic limiting
          dn_cl_nu = max(0.0,min(dn_cl_nu,(n_clmax-n_cl)/delt))

        else ! c_ccn is not constant, but dependent on n_cc

          if (l_sb_nuc_diff) then

            if (l_sb_sat_max) then ! l_sb_sat_max

              if ((ssat.lt.sat_max).and.(wdssatdz.gt.0.0)) then ! condition for nucleation
                dn_cl_nu = coef_ccn*n_cc*kappa_ccn*wdssatdz*ssat**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
              endif

            else  ! l_sb_sat_max

              if (wdssatdz.gt.0.0) then ! condition for nucleation
                dn_cl_nu = coef_ccn*n_cc*kappa_ccn*wdssatdz*ssat**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
              endif

            endif  ! l_sb_sat_max

          else  ! l_sb_nuc_diff

            if (l_sb_sat_max)  then ! l_sb_sat_max  AND NOT l_sb_nuc_diff

              if ((ssat.lt.sat_max).and.(w0(i,j,k).gt.0.0)) then ! condition for nucleation
                dn_cl_nu = (coef_ccn*n_cc/dzf(k-1))*max(0.0,w0(i,j,k)*   &
                  (ssat**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if ((ssat.lt.sat_max).and.(w0(i,j,k+1).lt.0.0)) then ! condition for nucleation
                dn_cl_nu = (coef_ccn*n_cc/dzf(k-1))*max(0.0,w0(i,j,k+1)* &
                   (ssat**kappa_ccn-ssat_u**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif

            else ! NOT l_sb_sat_max  AND NOT l_sb_nuc_diff

              if (w0(i,j,k).gt.0.0) then ! condition for nucleation
                dn_cl_nu = (coef_ccn*n_cc/dzf(k-1))*max(0.0,w0(i,j,k)*   &
                    (ssat**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if (w0(i,j,k+1).lt.0.0) then ! condition for nucleation
                dn_cl_nu = (coef_ccn*n_cc/dzf(k-1))*max(0.0,w0(i,j,k+1)* &
                   (ssat**kappa_ccn-ssat_u**kappa_ccn))    ! (1/rhof) *rhof = 1
              endif

            endif  ! l_sb_sat_max

          endif  ! l_sb_nuc_diff

          ! basic limiting
          dn_cl_nu = max(0.0,min(dn_cl_nu,(n_cc-n_cl)/delt))

        endif ! l_c_ccn

      endif ! ssat.gt.sat_min

    else ! l_sb_nuc_expl
      if(ssat.gt.0.0) then
        ! calculate number of activated n_ccn
        n_act = coef_ccn*n_cc*min(sat_max,ssat)**kappa_ccn
        n_act = max(n_cc,n_act)

        if (n_act.gt.n_cl) then
          dn_cl_nu = (n_act-n_cl)/delt
        endif

        ! basic limiting - not needed in this case
        ! dn_cl_nu = min(dn_cl_nu, (n_cc - n_cl))/delt)
      endif ! ssat.gt.0.0
    endif ! l_sb_nuc_expl

    ! increase in cloud water number
    n_clp = n_clp + dn_cl_nu

    ! update water density [kg kg^{-1}]
    q_clp = q_clp + x_cnuc * dn_cl_nu

    if (l_sb_dbg) then
      ! warning if too high liquid water
      if (ssat.gt.0.0 .and.( &
      (qt0(i,j,k)-qvsl(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci))/delt-x_cnuc*dn_cl_nu.lt. 0.)
      ) then
        write(6,*) 'WARNING: cloud nucleation too high'
        write(6,*) ' removing too much water'
        ! count(ssat.gt.0.0 .and.( &
        ! (qt0(i,j,k)-qvsl(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci))/delt-x_cnuc*dn_cl_nu.lt.0.))
      end if
    endif ! l_sb_dbg
  enddo
  enddo
  enddo
end subroutine  nucleation3


! ****************************************************************
! Ice nucleation
!  - based on Seifert & Beheng (2003), p. 53
! ****************************************************************
subroutine icenucle3
  use modglobal, only : i1,j1,k1,kmax,cp
  use modfields, only : rhof,exnf,qvsi,qt0,svm,tmp0   ! <- later remove svm - just for testing
  implicit none

  integer :: i,j,k

  real :: ssice
  real::  n_in,n_tid

  ! not always how supersaturated is water vapour, depending on a flag,
  ! it can also include water already in ice particles
  ssice = 0.0

  dn_ci_inu = 0.0

  do k=1,k1
  do j=2,j1
  do i=2,i1
    ! NOTE: not using qcmask condition, since air can be saturated with respect to ice
    if (tmp0(i,j,k).lt.tmp_inuc) then

      ! calculate supersaturation with respect to ice
      if (l_sb_inuc_sat) then  ! l_sb_inuc_sat
        ! calculating supersaturation of water vapour only
        ssice = (qt0(i,j,k)-q_cl/qvsi(i,j,k) -1.0
        ! ssice = max(0.0, (qt0(i,j,k)-q_cl)/qvsi(i,j,k) -1.0)
      else  ! l_sb_inuc_sat
        ! ie. cloud ice water is also included supersaturation
        ssice = (qt0(i,j,k)-q_cl+q_ci)/qvsi(i,j,k) -1.0
        ! ssice = max(0.0, (qt0(i,j,k)-q_cl+q_ci)/qvsi(i,j,k) -1.0)
      endif ! l_sb_inuc_sat

      if (ssice.gt.ssice_min) then ! condition for nucleation
        if (l_sb_inuc_expl) then ! l_sb_inuc_expl
          ! explicit ice nucleation - not yet included
          dn_ci_inu = 0. ! BUG: should give not-implemented error
        else  ! l_sb_inuc_expl
          ! Meyers et al. (1992)
          n_in = (1.0/rhof(k))*N_inuc*exp( a_M92              &
              +b_M92*min(ssice,ssice_lim))! o: N_inuc*exp( a_M92 + b_M92*ssice)

          ! whether to apply reisnerr correction
          if (l_sb_reisner) then
            ! prepare Reisnerr (1998) correction
            ! w: n_tid = (1.0/rhof(k))*N_inuc_R*exp( - min(tmp0(i,j,k),c_inuc_R)-T_3)
            n_tid = (1.0/rhof(k))*N_inuc_R*exp(b_inuc_R*(T_3- max(tmp0(i,j,k),c_inuc_R)))

            ! performing reisner correction
            n_in = max(a1_inuc_R*n_tid,min(a2_inuc_R*n_tid,n_in))
          endif ! l_sb_reisner

          ! limiting n_in
          n_in = min(n_i_max, n_in)

          ! checking conditions if suitable for nucleation
          if (n_ci.lt.n_in) then ! condition intentionally left this way
            dn_ci_inu = (n_in-n_ci)/delt
          else
            dn_ci_inu = 0.0
            ! note - written this way on purpose
            ! in case of further adjustment of nucleation subroutine
          endif
        endif ! l_sb_inuc_expl

        ! basic correction  - not suitable

        ! update cloud water number
        n_cip = n_cip + dn_ci_inu

        ! update water density [kg kg^{-1}]
        q_cip = q_cip + x_inuc*dn_ci_inu

        ! update liquid water potential temperature
        !  - due to latent heat of melting (freezing in this case)
        qtpmcr(i,j,k) = qtpmcr(i,j,k) - x_inuc*dn_ci_inu
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlvi/(cp*exnf(k)))*x_inuc*dn_ci_inu

        if (l_sb_dbg) then
          if ((qt0(i,j,k)-qvsi(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci))/delt - x_inuc*dn_ci_inu.lt. 0.) then
            write(6,*) 'WARNING: high ice nucleation'
            write(6,*) ' removing too much water'
            ! count((qt0(i,j,k)-qvsi(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci))/delt - x_inuc*dn_ci_inu.lt. 0.)
          endif
        endif ! l_sb_dbg
      endif ! ssice.gt.ssice_min
    endif ! tmp0(i,j,k).lt.tmp_inuc
  enddo
  enddo
  enddo
end subroutine icenucle3


! ****************************************
! Homogeneous freezing
! of cloud droplets
! 
! ****************************************     
subroutine homfreez3
  use modglobal, only : i1,j1,k1,rlv,cp
  use modfields, only : svm,tmp0,exnf
  implicit none

  real    :: J_hom   ! freezing rate
  real    :: tmpj    ! adjusted temperature

  real    :: dq_cl_hom !< homogeneous freezing of cloud water
  real    :: dn_cl_hom !< homogeneous freezing of cloud water

  ! calculate constants
  real, parameter :: expC_30 = exp(C_30_Cf02)

  integer :: i,j,k

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_cl_mask.and.tmp0(i,j,k).lt.T_3) then

      ! homogeneous freezing of cloud droplets
      ! --------------------------------------

      ! inserting temperature values
      !  -- can be later adjusted to include the effect of chemicals in clouds
      tmpj = tmp0(i,j,k)

      if (tmpj>tmp_lim1_Cf02) then
        J_hom = exp(C_Cf02+B_Cf02*(tmpj+CC_Cf02))
      else if(tmpj<tmp_lim2_Cf02) then
        J_hom = expC_30
      else
        J_hom = exp(C_20_Cf02                 &
           + B_21_Cf02*(tmpj-offset_Cf02)     &
           + B_22_Cf02*(tmpj-offset_Cf02)**2  &
           + B_23_Cf02*(tmpj-offset_Cf02)**3  &
           + B_24_Cf02*(tmpj-offset_Cf02)**4  )
      endif

      dn_cl_hom        = -c_mmt_1cl *n_cl(i,j,k) * x_cl(i,j,k)* J_hom
      dq_cl_hom        = -c_mmt_2cl *q_cl(i,j,k) * x_cl(i,j,k)* J_hom
      dn_cl_hom        = max(dn_cl_hom,min(0.0,-svm(i,j,k,in_cl)/delt-n_clp(i,j,k)))
      dq_cl_hom        = max(dq_cl_hom,min(0.0,-svm(i,j,k,iq_cl)/delt-q_clp(i,j,k)))

      ! changes in cloud water
      n_clp = n_clp + dn_cl_hom
      q_clp = q_clp + dq_cl_hom

      ! increase in cloud ice
      n_cip = n_cip - dn_cl_hom
      q_cip = q_cip - dq_cl_hom

      ! changes in the total amount of water
      qtpmcr(i,j,k) = qtpmcr(i,j,k) + dq_cl_hom

      ! change in th_l due to freezing
      thlpmcr(i,j,k) = thlpmcr(i,j,k)-((rlv+rlme)/(cp*exnf(k)))*dq_cl_hom
    endif ! cl_mask
  enddo
  enddo
  enddo
 end subroutine homfreez3      


! ****************************************
! Heterogeneous freezing of cloud droplets
! ****************************************     
subroutine hetfreez3
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
    if (q_cl_mask.and.tmp0(i,j,k).lt.T_3) then
      J_het = A_het *exp( B_het*(T_3-tmp0(i,j,k)) -1)

      dn_cl_het = -c_mmt_1cl *n_cl(i,j,k) * x_cl(i,j,k)* J_het
      dq_cl_het = -c_mmt_2cl *q_cl(i,j,k) * x_cl(i,j,k)* J_het

      ! basic correction
      dn_cl_het = max(dn_cl_het,min(0.0,-svm(i,j,k,in_cl)/delt-n_clp(i,j,k)))
      dq_cl_het = max(dq_cl_het,min(0.0,-svm(i,j,k,iq_cl)/delt-q_clp(i,j,k)))

      ! changes in cloud water
      n_clp = n_clp + dn_cl_het
      q_clp = q_clp + dq_cl_het

      ! increase in cloud ice
      n_cip = n_cip - dn_cl_het
      q_cip = q_cip - dq_cl_het

      ! changes in the total amount of water
      qtpmcr(i,j,k) = qtpmcr(i,j,k) + dq_cl_het

      ! change in th_l due to freezing
      thlpmcr(i,j,k) = thlpmcr(i,j,k)-((rlv+rlme)/(cp*exnf(k)))*(dq_cl_het)
    endif ! cl_mask
  enddo
  enddo
  enddo
end subroutine hetfreez3   


! ****************************************
! Depositional growth of cloud ice particles
!  later tunr into: 
!  
! Depositional growth of various ice particles
! Inner subroutine - takes input from a wrapper
!  following S&B 
! ****************************************
subroutine deposit_ice3
  use modglobal, only : i1,j1,k1,rv,rd,cp,pi
  use modfields, only : exnf,qt0,svm,tmp0,,qvsi,rhof,presf
  implicit none

  real :: esi

  real :: F_ci    & ! ventilation factor
         ,q_avail & ! available water for deposition
         ,Si      & ! super or undersaturation with respect to ice
         ,G       & ! G_iv function
         ,Dvic_ci & ! size of ice parti
         ,viic_ci & ! v for ice cloud particles
         ,nrex_ci   ! reynolds number

  integer :: i,j,k

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_ci_mask.AND.(tmp0(i,j,k).le.T_3)) then
      q_avail = qt0(i,j,k)-q_cl(i,j,k)- qvsi(i,j,k)
      Si = q_avail/qvsi(i,j,k)

      ! calculating G_iv
      esi = qvsi(i,j,k)*presf(k)/(rd/rv+(1.0-rd/rv)*qvsi(i,j,k))
      G = (rv * tmp0(i,j,k)) / (Dv*esi) + rlvi/(Kt*tmp0(i,j,k))*(rlvi/(rv*tmp0(i,j,k)) -1.)
      G = 1./G

      ! diameter
      Dvic_ci = a_ci*x_ci**b_ci

      ! terminal velocity
      viic_ci = al_ci*((rho0s/rhof(k))**0.5)*x_ci**be_ci

      ! N_re Reynolds number
      nrex_ci = Dvic_ci*viic_ci/nu_a

      ! calculating from prepared ventilation coefficients
      F_ci = avent_1i+bvent_1i*Sc_num**(1.0/3.0)*nrex_ci**0.5

      ! depositional growth
      ! k_depos = 4*pi/c  for spherical particles
      ! k_depos = 4*pi/c  for hexagonal plates - included
      dq_ci_dep = (4*pi/c_ci)*n_ci*G*Dvic_ci*F_ci*Si

      ! and limiting not to deposit more than available
      ! ie. allows any negative dq_ci_dep but positive only smaller than amount of available water
      dq_ci_dep = max(dq_ci_dep,min(0.0,-svm(i,j,k,iq_ci)/delt-q_cip(i,j,k)))
      dq_ci_dep = min(dq_ci_dep,max(0.0,q_avail/delt))

      ! adds to outputs
      q_cip = q_cip + dq_ci_dep

      qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_ci_dep
      thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlvi/(cp*exnf(k)))*(dq_ci_dep)

      if (l_sb_dbg) then
        if((svm(i,j,k,iq_ci)+delt*dq_ci_dep).lt. 0.0) then
          write(6,*) 'WARNING: ice_deposit3 too low'
        endif
        if(((qt0(i,j,k)-qvsi(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci)-delt*dq_ci_dep).lt.0.0) then
          write(6,*) 'WARNING: ice_deposit3 too high depositing more ice than available',
          ! too high:    ((qt0(i,j,k)-qvsi(i,j,k)-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci)-delt*dq_ci_dep).lt. 0.0)
          ! negative qt: ((qt0(i,j,k)-delt*dq_ci_dep).lt. 0.0)
        endif
      endif
    endif
  enddo
  enddo
  enddo
end subroutine deposit_ice3


! ****************************************
! Depositional growth of snow particles
!  later turn into: 
!  
! Depositional growth of various ice particles
! Inner subroutine - takes input from a wrapper
!  following S&B 
! ****************************************  
subroutine deposit_snow3
  use modglobal, only : i1,j1,k1,rv,rd,cp,pi
  use modfields, only : exnf,qt0,svm,tmp0,,qvsi,rhof,presf
  implicit none

  real :: esi

  real :: F_hs    & ! ventilation factor
         ,q_avail & ! available water for deposition
         ,Si      & ! super or undersaturation with respect to ice
         ,G       & ! G_iv function
         ,Dvic_hs & ! size of ice parti
         ,viic_hs & ! v for ice cloud particles
         ,nrex_hs   ! reynolds number

  integer :: i,j,k

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_hs_mask.AND.(tmp0(i,j,k).le.T_3)) then
      q_avail = qt0(i,j,k)-q_cl(i,j,k)- qvsi(i,j,k)
      Si = q_avail/qvsi(i,j,k)

      ! calculating G_iv
      esi = qvsi(i,j,k)*presf(k)/(rd/rv+(1.0-rd/rv)*qvsi(i,j,k))
      G = (rv * tmp0(i,j,k)) / (Dv*esi) + rlvi/(Kt*tmp0(i,j,k))*(rlvi/(rv*tmp0(i,j,k)) -1.)
      G = 1./G

      ! diameter
      Dvic_hs = a_hs*x_hs**b_hs

      ! terminal velocity
      viic_hs = al_hs*((rho0s/rhof(k))**0.5)*x_hs**be_hs

      ! N_re Reynolds number
      nrex_hs = Dvic_hs*viic_hs/nu_a

      ! calculating from prepared ventilation coefficients
      F_hs = aven_1s+bven_1s*Sc_num**(1.0/3.0)*nrex_hs**0.5

      ! depositional growth
      ! k_depos = 4*pi/c  for spherical particles
      ! k_depos = 4*pi/c  for hexagonal plates - included
      dq_hs_dep = (4*pi/c_hs)*n_hs*G*Dvic_hs*F_hs*Si

      ! and limiting not to deposit more than available
      dq_hs_dep = max(dq_hs_dep,min(0.0,-svm(i,j,k,iq_hs)/delt-q_hsp(i,j,k)))
      dq_hs_dep = min(dq_hs_dep,max(0.0,q_avail/delt))

      ! adds to outputs
      q_hsp = q_hsp + dq_hs_dep

      qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_hs_dep
      thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlvi/(cp*exnf(k)))*dq_hs_dep

      if (l_sb_dbg) then
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
      endif
    endif
  enddo
  enddo
  enddo
end subroutine deposit_snow3


! ****************************************
! Depositional growth of graupel particles
! as well as sublimation
!  later turn into: 
!  
! Depositional growth of various ice particles
! Inner subroutine - takes input from a wrapper
!  following S&B 
! ****************************************  
subroutine deposit_graupel3
  use modglobal, only : i1,j1,k1,rv,rd,cp,pi
  use modfields, only : exnf,qt0,svm,tmp0,,qvsi,rhof,presf
  implicit none

  real :: esi

  real :: cor_dqhg_dep

  real :: F_hg    & ! ventilation factor
         ,q_avail & ! available water for deposition
         ,Si      & ! super or undersaturation with respect to ice
         ,G       & ! G_iv function
         ,Dvic_hg & ! size of particles
         ,viic_hg & ! v for particles
         ,nrex_hg   ! reynolds number

  integer :: i,j,k

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_hg_mask.AND.(tmp0(i,j,k).le.T_3)) then
      q_avail = qt0(i,j,k)-q_cl(i,j,k)- qvsi(i,j,k)
      Si = q_avail/qvsi(i,j,k)

      ! calculating G_iv
      esi = qvsi(i,j,k)*presf(k)/(rd/rv+(1.0-rd/rv)*qvsi(i,j,k))
      G = (rv * tmp0(i,j,k)) / (Dv*esi) + rlvi/(Kt*tmp0(i,j,k))*(rlvi/(rv*tmp0(i,j,k)) -1.)
      G = 1./G

      ! diameter
      Dvic_hg = a_hg*x_hg**b_hg

      ! terminal velocity
      viic_hg = al_hg*((rho0s/rhof(k))**0.5)*x_hg**be_hg

      ! N_re Reynolds number
      nrex_hg = Dvic_hg*viic_hg/nu_a

      ! calculating from prepared ventilation coefficients
      F_hg = aven_1g+bven_1g*Sc_num**(1.0/3.0)*nrex_hg**0.5

      ! depositional growth
      ! k_depos = 4*pi/c  for spherical particles
      ! k_depos = 4*pi/c  for hexagonal plates - included
      dq_hg_dep = (4*pi/c_hg)*n_hg*G*Dvic_hg*F_hg*Si

      dq_hg_dep = max(dq_hg_dep,min(0.0,-svm(i,j,k,iq_hg)/delt-q_hgp))
      dq_hg_dep = min(dq_hg_dep,max(0.0,q_avail/delt))

      ! adds to outputs
      q_hgp = q_hgp + dq_hg_dep

      qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_hg_dep
      thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlvi/(cp*exnf(k)))*dq_hg_dep

      if (l_sb_dbg) then
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
  enddo
  enddo
  enddo
end subroutine deposit_graupel3   


!! ****************************************************************
!!  Limiting condensation and deposition
!!  - to prevent negative values
!! 
!!  - call should be located:
!!      - after depositions
!!      - after nucleation correction
!!      - before heterogeneous freezing 
!! 
!!  ************************************************************
subroutine cor_deposit3
  use modglobal, only : i1,j1,k1,rv,rd,cp,pi
  use modfields, only : exnf,qt0,svm,tmp0,,qvsi,rhof,presf
  implicit none

  real    :: tocon,precon,cond_cf
  real    :: cor_dqci_dep,cor_dqhs_dep,cor_dqhg_dep
  integer :: i,j,k

  do k=1,k1
  do j=2,j1
  do i=2,i1
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
      q_cip = q_cip + cor_dqci_dep
      q_hsp = q_hsp + cor_dqhs_dep
      q_hgp = q_hgp + cor_dqhg_dep

      ! - correction for total water
      qtpmcr(i,j,k) = qtpmcr(i,j,k) &
                      -cor_dqhs_dep &
                      -cor_dqhg_dep &
                      -cor_dqci_dep

      ! - and correcting for heat
      thlpmcr(i,j,k) = thlpmcr(i,j,k)                     &
            +(rlvi/(cp*exnf(k)))*cor_dqci_dep             &
            +(rlvi/(cp*exnf(k)))*cor_dqhs_dep             &
            +(rlvi/(cp*exnf(k)))*cor_dqhg_dep

      ! corrector for process values
      dq_ci_dep = dq_ci_dep + cor_dqci_dep
      dq_hs_dep = dq_hs_dep + cor_dqhs_dep
      dq_hg_dep = dq_hg_dep + cor_dqhg_dep
    endif
  enddo
  enddo
  enddo
end subroutine cor_deposit3


! ****************************************
! ice selfcollection and aggregation to snow
!  - all collection of ice by ice treated as snow 
! 
! follows Seifert, 2002
! ****************************************     
subroutine ice_aggr3
  use modglobal, only : i1,j1,k1,rv,pi
  use modfields, only : qt0,svm,tmp0,rhof
  implicit none
  integer :: i,j,k

  real    :: dlt_0aa_i &
            ,dlt_1aa_i &
            ,th_0aa_i  &
            ,th_1aa_i

  real :: dif_D_10, x_crit_ii, x_minagg_ii, rem_cf_i

  ! Temporary variables, value changes several times in this routine
  real :: E_ab, E_stick

  ! calculate constants
  dlt_0aa_i   = 2*dlt_i0 + dlt_i0i
  dlt_1aa_i   = dlt_i0 + dlt_i1i + dlt_i1
  th_0aa_i    = 2*th_i0 - th_i0i   ! from Seifert, 2002
  th_1aa_i    = th_i0 - th_i1i + th_i1

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if(q_ci_mask.and.(x_ci.gt.x_crit_ii).and.(q_ci.gt.q_crit_ii)) then
      ! prepare coefficient for remaining water number
      rem_cf_i = (1.0-rem_n_ci_min)/delt

      ! and the minimal conversion size
      x_minagg_ii = (D_i_b/a_ci)**(1.0/b_ci)

      ! and the critical size for the start of conversion
      x_crit_ii   = (D_crit_ii/a_ci)**(1.0/b_ci)

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
            dif_D_10  = D_i_b-D_i_a  ! difference in ice cloud droplet intervals
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
                      *dlt_1aa_i*n_ci(i,j,k)*q_ci*D_ci**2     &
                      *(th_1aa_i*v_ci**2+2*sigma_ci**2)**0.5


      ! limiting dq_ci  --   limit x_cv_ii as per ICON 2017
      dq_ci_col_iis = max(dq_ci_col_iis, (-svm(i,j,k,iq_ci)/delt-q_cip(i,j,k)))
      dn_ci_col_iis = max(dn_ci_col_iis, (-rem_cf_i*svm(i,j,k,in_ci)-n_cip(i,j,k)))
      dn_ci_col_iis = max(dn_ci_col_iis, dq_ci_col_iis/x_minagg_ii)
    endif

    n_cip    = n_cip + dn_ci_col_iis
    q_cip    = q_cip + dq_ci_col_iis
    n_hsp    = n_hsp - dn_ci_col_iis
    q_hsp    = q_hsp - dq_ci_col_iis

    if (l_sb_dbg) then
      if((svm(i,j,k,iq_ci)+delt*dq_ci_col_iis .lt. 0.0)) then
        write(6,*) 'WARNING: ice_aggr3 too high'
        write(6,*) ' removing more ice than available'
        ! count((svm(i,j,k,iq_ci) +delt*dq_ci_col_iis).lt. 0.0)
        write(6,*) ' removing too much ice'
        ! count(( q_ci+delt*dq_ci_col_iis).lt. 0.0 )
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
    endif
  enddo
  enddo
  enddo
end subroutine ice_aggr3


! *************************************************
! snow selfcollection
! 
! *************************************************
subroutine snow_self3
  use modglobal, only : i1,j1,k1,rv,pi
  use modfields, only : qt0,svm,tmp0,rhof
  implicit none
  integer :: i,j,k

  real :: dif_D_10, rem_cf_s

  ! Temporary variables, value changes several times in this routine
  real :: E_ab, E_stick

  ! adjusting coefficient
  ! prepare coefficient for remaining water number
  rem_cf_s = (1.0-rem_n_hs_min)/delt

  do k=1,k1
  do j=2,j1
  do i=2,i1
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
    endif

    n_hsp    = n_hsp + dn_hs_col_sss

    if (l_sb_dbg) then
      if(( svm(i,j,k,in_hs)+delt*dn_hs_col_sss .lt. 0.0 )) then
        write(6,*) 'WARNING: snow self-collection too high'
        write(6,*) ' decreasing number of snowflakes below 0'
        ! count((svm(i,j,k,in_hs)+delt*dn_hs_col_sss).lt. 0.0 )
        write(6,*) ' decreasing number of snowflakes too much'
        ! count((n_hs(i,j,k)+delt*dn_hs_col_sss).lt. 0.0 )
      endif
  enddo
  enddo
  enddo
end subroutine snow_self3


! ****************************************
! snow collecting cloud ice
! 
! s + i -> s
! ****************************************     
subroutine coll_sis3
  use modglobal, only : i1,j1,k1,rv,pi
  use modfields, only : qt0,svm,tmp0,rhof
  implicit none
  integer :: i,j,k

  ! Temporary variables, value changes several times in this routine
  real :: E_ab, E_stick, rem_cf_i

  ! adjusting coefficient
  ! prepare coefficient for remaining water number
  rem_cf_i = (1.0-rem_n_ci_min)/delt

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_hs_mask.and.q_ci_mask) then
      ! calculating sticking efficiency
      E_stick = c_E_o_s*exp(B_stick *(tmp0(i,j,k)+stick_off))
      E_stick =min(c_E_o_s,E_stick)
      E_ab = E_ee_m*E_stick

      dq_hsci_col = (rhof(k)*pi/4)*E_ab*n_hs(i,j,k)                &
                    *q_ci*(dlt_s0*D_hs**2                   &
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
    endif

    n_cip    = n_cip + dn_ci_col_hs
    q_cip    = q_cip + dq_hsci_col
    q_hsp    = q_hsp - dq_hsci_col

    if (l_sb_dbg) then
      if(svm(i,j,k,iq_ci)-delt*dq_hsci_col.lt. 0.0) then
        write(6,*) 'WARNING: coll_sis3 too high removing more ice than available'
        write(6,*) ' removing more ice than available'
        ! count((svm(i,j,k,iq_ci)-delt*dq_hsci_col).lt. 0.0)
        write(6,*) ' removing too much ice'
        ! count((q_ci-delt*dq_hsci_col).lt. 0.0)
        write(6,*) ' getting negative q_t'
        ! count((qt0(i,j,k)-delt*q_cip).lt. 0.0 )
      endif

      if(svm(i,j,k,in_ci)+delt*dn_ci_col_hs.lt. 0.0) then
        write(6,*) 'WARNING: coll_sis3 too high'
        write(6,*) ' removing more ice particles then available in gridpoints'
        ! count(svm(i,j,k,in_ci)+delt*dn_ci_col_hs.lt. 0.0)
        write(6,*) ' removing too many ice particles in gridpoints '
        ! count((n_ci(i,j,k)+delt*dn_ci_col_hs).lt. 0.0)
      endif
    endif
  enddo
  enddo
  enddo
end subroutine coll_sis3    


! ****************************************
! snow collecting cloud ice
! 
! g+s -> g 
! ****************************************     
subroutine coll_gsg3
  use modglobal, only : i1,j1,k1,rv,pi
  use modfields, only : qt0,svm,tmp0,rhof
  implicit none
  integer :: i,j,k

  real :: E_ab, E_stick, rem_cf_s

  ! adjusting coefficient
  ! prepare coefficient for remaining water number
  rem_cf_s = (1.0-rem_n_hs_min)/delt

  do k=1,k1
  do j=2,j1
  do i=2,i1
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

      dq_hghs_col = min(dq_hghs_col,&
                        max(0.0,    &
                        svm(i,j,k,iq_hs)/delt+q_hsp(i,j,k)))   ! following ICON, 2017
      dn_hs_col_hg = max(dn_hs_col_hg,&
                         min(0.0,     &
                         -rem_cf_s*svm(i,j,k,in_hs)-n_hsp(i,j,k)))
    endif

    n_hsp    = n_hsp + dn_hs_col_hg
    q_hsp    = q_hsp - dq_hghs_col
    q_hgp    = q_hgp + dq_hghs_col

    if (l_sb_dbg) then
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
end subroutine coll_gsg3    


!****************************************
! Collection of clpoud droplets by ice
! 
! - based on Seifert & Beheng (2004)
! - resulting process is:
!    - ice riming by cloud ice
!    - enhanced melting 
!
!****************************************
subroutine coll_ici3
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

  do k=1,k1
  do j=2,j1
  do i=2,i1

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
        dn_cl_rime_ci = max(dn_col_ici,min(0.0,-rem_cf*svm(i,j,k,in_cl)-n_clp))
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
  enddo
  enddo
  enddo
end subroutine coll_ici3  


! ****************************************
!  riming of snow
! ****************************************     
subroutine coll_scs3
  use modglobal, only : i1,j1,k1,cp,pi
  use modfields, only : qt0,svm,tmp0,rhof,exnf

  implicit none
  integer :: i,j,k

  real :: E_ab
  real :: dif_D_10
  real :: k_enhm
  real :: rem_cf
  real :: dn_col_b, dq_col_a

  do k=1,k1
  do j=2,j1
  do i=2,i1
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
  enddo
  enddo
  enddo
end subroutine coll_scs3  


! ****************************************
!  riming of graupel by cloud droplets
! 
! 
! ****************************************     
subroutine coll_gcg3
  use modglobal, only : i1,j1,k1,cp,pi
  use modfields, only : qt0,svm,tmp0,rhof,exnf

  implicit none
  integer :: i,j,k

  real :: E_ab
  real :: dif_D_10
  real :: k_enhm
  real :: rem_cf
  real :: dn_col_b, dq_col_a

  do k=1,k1
  do j=2,j1
  do i=2,i1
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
  enddo
  enddo
  enddo
end subroutine coll_gcg3    


! ****************************************
!  riming of graupel
! 
! 
! ****************************************     
subroutine coll_grg3
  use modglobal, only : i1,j1,k1,cp,pi
  use modfields, only : qt0,svm,tmp0,rhof,exnf

  implicit none
  integer :: i,j,k

  real :: E_ab
  real :: dif_D_10
  real :: k_enhm
  real :: rem_cf
  real :: dn_col_ici, dq_col_ici
  real :: dn_col_b, dq_col_a

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_hg_mask.and.q_hr_mask) then
      E_ab =  E_er_m

      dq_col_a = (rhof(k)*pi/4)*E_ab*n_hg    &
            *q_hr*(dlt_g0*D_hg**2                     &
              +dlt_g1r*D_hg*D_hr+dlt_r1*D_hr**2) &
            *( th_g0*v_hg**2-th_g1r*v_hr*v_hg    &
              +th_r1*v_hr**2+sigma_hg**2+sigma_hr**2)**0.5

      dn_col_b = -(rhof(k)*pi/4)*E_ab*n_hg   &
            *n_hr*(dlt_g0*D_hg**2                     &
              +dlt_g0r*D_hg*D_hr+dlt_r0*D_hr**2) &
            *( th_g0*v_hg**2-th_g0r*v_hr*v_hg    &
              +th_r0*v_hr**2+sigma_hg**2+sigma_hr**2)**0.5

      if (tmp0(i,j,k).lt.T_3) then
        ! riming only
        ! remaining number of particles
        rem_cf = (1.0-rem_n_hr_min)/delt

        dq_hghr_rime = min(dq_col_a,&
                           max(0.0,&
                           svm(i,j,k,iq_hr)/delt+q_hrp))   ! following ICON, 2017
        ! = min(dq_hghr_rime,svm(i,j,k,iq_hr)/delt)

        dn_hr_rime_hg = max(dn_col_b,&
                            min(0.0,&
                            -rem_cf*svm(i,j,k,in_hr)-n_hrp))
        ! = min(dn_hr_rime_hg,-svm(i,j,k,in_hr)/delt)

        ! change in the amount of graupel
        q_hgp = q_hgp + dq_hghr_rime

        ! change in the amount of rain
        n_hrp = n_hrp + dn_hr_rime_hg
        q_hrp = q_hrp - dq_hghr_rime

        ! no change in q_t
        ! qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_col_a

        ! change in th_l - just heat release from freezing
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+ (rlme/(cp*exnf(k)))*dq_hghr_rime
      else
        ! not riming,but enhanced melting and scheding

        ! enhanced melting
        k_enhm  = c_water/rlme

        ! calculating the melting
        dq_hg_eme_gr = -k_enhm*(tmp0(i,j,k)-T_3)*min(dq_col_a,q_hr/delt)
        dq_hg_eme_gr = max(dq_hg_eme_gr,min(0.0,-svm(i,j,k,iq_hg)/delt-q_hgp))

        ! calculating number of melted particles
        ! - expected to be proportional to melted mass
        dn_hg_eme_gr = dq_hg_eme_gr*n_hg/(q_hg+eps0) ! q_hg here is always some small positive number

        ! - but not more than number of interacting particles
        ! dn_hg_eme_gr = max(dn_hg_eme_gr, max(dn_col_b,-n_hr/delt))
        ! - and not more than total number of particles
        dn_hg_eme_gr = max(dn_hg_eme_gr,&
                           -min((svm(i,j,k,in_hr)/delt+n_hrp),&
                                (svm(i,j,k,in_hg)/delt+n_hgp(i,j,k))))

        ! updating tendencies
        ! updating rain
        q_hrp = q_hrp - dq_hg_eme_gr
        ! n_hrp = n_hrp + dn_hr_rime_hg

        ! graupel
        q_hgp = q_hgp + dq_hg_eme_gr
        n_hgp = n_hgp + dn_hg_eme_gr

        ! updating cloud water
        ! no change - all turned into rain
        ! updating thermodynamic
        ! qt : no change -  increased by melted water goes to rain
        ! qtpmcr(i,j,k) = qtpmcr(i,j,k) +0.0
        ! thl : melting
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlme/(cp*exnf(k)))*dq_hg_eme_gr
      endif

      if (l_sb_dbg) then
        if(svm(i,j,k,iq_hr)-delt*dq_col_a.lt. 0.0) then
          write(6,*) 'WARNING: coll_grg3 too high'
          write(6,*) ' removing more rain water than available'
          ! count(svm(i,j,k,iq_hr)-delt*dq_hghr_rime.lt. 0.0)
          write(6,*) ' removing too much rain water '
          ! count(q_hr-delt*dq_hghr_rime.lt. 0.0)
          write(6,*) ' getting negative q_t in  '
          ! count(qt0(i,j,k)+delt*q_hrp).lt. 0.0)
        endif
        if(svm(i,j,k,in_hr)+delt*dn_col_b.lt. 0.0) then
          write(6,*) 'WARNING: coll_grg too high'
          write(6,*) ' removing more raindrops than available '
          ! count(svm(i,j,k,in_hr)+delt*dn_hr_rime_hg.lt. 0.0)
          write(6,*) ' removing too many raindrops'
          ! count(n_hr+delt*dn_hr_rime_hg.lt. 0.0)
        endif
      endif
    endif
end subroutine coll_grg3


! ****************************************
! Heterogeneos freezing of rain
! ****************************************
subroutine rainhetfreez3
  use modglobal, only : i1,jh,j1,k1,pi
  use modfields, only : svm,tmp0,exnf
  implicit none

  integer :: i,j,k
  real    :: J_het

  dn_hr_het = 0.0
  dq_hr_het = 0.0

  ! performing calculation for the freezing rate
  do j=2,j1
  do i=2,i1
  do k=1,k1
    ! Heterogeneos freezing of rain
    ! -----------------------------
    if (q_hr_mask.and.(tmp0(i,j,k).lt.T_3)) then
      ! maybe only for temperatures below T_3 ?
      J_het = A_het *exp( B_het*(T_3-tmp0(i,j,k)) -1)

      dn_hr_het = -c_mmt_1hr * n_hr * x_hr * J_het
      dq_hr_het = -c_mmt_2hr * q_hr * x_hr * J_het

      ! basic correction
      dn_hr_het = max(dn_hr_het,min(0.0,-svm(i,j,k,in_hr)/delt-n_hrp))
      dq_hr_het = max(dq_hr_het,min(0.0,-svm(i,j,k,iq_hr)/delt-q_hrp))

      ! decrease in raindrops
      n_hrp = n_hrp + dn_hr_het
      q_hrp = q_hrp + dq_hr_het

      ! increase in graupel
      n_hgp = n_hgp - dn_hr_het
      q_hgp = q_hgp - dq_hr_het

      ! and consumption of aerosols for heterogeneous freezing  ?
      ! n_ccp = n_ccp - dn_hr_het

      ! no changes in the total amount of water
      ! qtpmcr(i,j,k) = qtpmcr(i,j,k)

      ! change in th_l due to freezing
      thlpmcr(i,j,k) = thlpmcr(i,j,k) - (rlme/(cp*exnf(k)))*dq_hr_het
    endif
  enddo
  enddo
  enddo
end subroutine rainhetfreez3


! riming of ice + rain to graupel
! -------------------------------
subroutine coll_rig3
  use modglobal, only : i1,j1,k1,cp,pi
  use modfields, only : exnf,svm,tmp0,rhof
  implicit none
  integer :: i,j,k
  real    :: rem_ci_cf,rem_hr_cf,k_enhm
  real    :: E_ab, dn_col_b, dq_col_a, dq_col_b

  ! setting up extra coefficients
  ! remain coefficients
  rem_ci_cf = (1.0-rem_n_ci_min)/delt
  rem_hr_cf = (1.0-rem_n_hr_min)/delt

  ! enhanced melting
  k_enhm  = c_water/rlme

  ! calculating
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_ci_mask.and.q_hr_mask) then
      E_ab =  E_i_m

      dq_col_b = -(rhof(k)*pi/4)*E_ab*n_ci           &
           *q_hr*(dlt_i0*D_ci**2                     &
             +dlt_i1r*D_ci*D_hr+dlt_r1*D_hr**2)      &
           *( th_i0*v_ci**2-th_i1r*v_hr*v_ci         &
             +th_r1*v_hr**2+sigma_ci**2+sigma_hr**2)**0.5

      dn_col_b = -(rhof(k)*pi/4)*E_ab*n_ci           &
           *n_hr*(dlt_i0*D_ci**2                     &
             +dlt_i0r*D_ci*D_hr+dlt_r0*D_hr**2)      &
           *( th_i0*v_ci**2-th_i0r*v_hr*v_ci         &
             +th_r0*v_hr**2+sigma_ci**2+sigma_hr**2)**0.5

      dq_col_a = (rhof(k)*pi/4)*E_ab*q_ci            &
           *n_hr*(dlt_i1*D_ci**2                     &
             +dlt_r1i*D_ci*D_hr+dlt_r0*D_hr**2)      &
           *( th_i1*v_ci**2-th_r1i*v_hr*v_ci         &
             +th_r0*v_hr**2+sigma_ci**2+sigma_hr**2)**0.5

      dq_hr_col_ri =  dq_col_b
      dn_ci_col_ri =  dn_col_b
      dn_hr_col_ri =  dn_col_b
      dq_ci_col_ri = -dq_col_a


      ! first adjustment
      dq_ci_col_ri = max(dq_ci_col_ri,min(0.0,-svm(i,j,k,iq_ci)/delt-q_cip))   ! following ICON, 2017
      dq_hr_col_ri = max(dq_hr_col_ri,min(0.0,-svm(i,j,k,iq_hr)/delt-q_hrp))   ! following ICON, 2017

      ! adjustment of numbers - both ice and water
      dn_ci_col_ri = max(dn_ci_col_ri,min(0.0,-rem_ci_cf*svm(i,j,k,in_ci)-n_cip))
      dn_ci_col_ri = max(dn_ci_col_ri,min(0.0,-rem_hr_cf*svm(i,j,k,in_hr)-n_hrp))

      if(tmp0(i,j,k).lt.T_3) then
        ! the collection is just riming
        ! decrease in numeber of raindrops same as decrease in number of ice
        dn_hr_col_ri = dn_ci_col_ri

        ! record the change in cloud ice
        q_cip = q_cip + dq_ci_col_ri
        n_cip = n_cip + dn_ci_col_ri

        ! change in rain
        q_hrp = q_hrp + dq_hr_col_ri
        n_hrp = n_hrp + dn_hr_col_ri

        ! and for graupel
        q_hgp = q_hgp - dq_ci_col_ri - dq_hr_col_ri
        n_hgp = n_hgp - dn_ci_col_ri

        ! change in q_t - decrease in cloud ice
        qtpmcr(i,j,k) = qtpmcr(i,j,k) + 0.0 ! dq_ci_col_ri

        ! change in th_l - release from freezing and removal of ice
        thlpmcr(i,j,k) = thlpmcr(i,j,k)                  &
            -(rlme/(cp*exnf(k)))*dq_hr_col_ri
      else  ! tmp0(i,j,k).gt.T_3
        ! enhanced melting and graupel formation

        ! calculate the melting
        dq_ci_eme_ri = k_enhm*(tmp0(i,j,k)-T_3)*dq_hr_col_ri ! with + due to negative value of dq_hr here

        ! limit melting
        dq_ci_eme_ri = max(dq_ci_eme_ri,min(0.0,-svm(i,j,k,iq_ci)/delt-q_cip))

        ! calculate how many ice perticles melted
        ! q_ci here is always some small positive number
        dn_ci_eme_ri = dq_ci_eme_ri*n_ci/(q_ci+eps0)

        ! limit so it dos not melt more than interacting
        ! also limit so that new graupel not larger that max mean size of source ice ?
        dn_ci_eme_ri =max(dn_ci_eme_ri,dn_ci_col_ri)

        ! update ice
        q_cip = q_cip + dq_ci_eme_ri ! dq_ci_col_ri
        n_cip = n_cip + dn_ci_eme_ri ! dn_ci_col_ri

        ! no collection of raindrops
        dq_hr_col_ri = 0.0
        dq_ci_col_ri = 0.0
        dn_hr_col_ri = 0.0
        dn_ci_col_ri = 0.0

        ! increase in rain mass
        q_hrp = q_hrp - dq_ci_eme_ri

        ! change in thl - heat spent on melting
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlme/(cp*exnf(k)))*dq_ci_eme_ri
      endif
   endif
  enddo
  enddo
  enddo
end subroutine coll_rig3


! riming of rain + snow to graupel
! --------------------------------
subroutine coll_rsg3

  use modglobal, only : i1,j1,k1,cp,pi
  use modfields, only : exnf,svm,tmp0,rhof
  implicit none

  integer :: i,j,k
  real    ::
  real    :: rem_cf, k_enhm
  real    :: E_ab, dn_col_b, dq_col_a, dq_col_b

  dn_col_b = 0.0
  dq_col_a = 0.0
  dq_col_b = 0.0

  ! calculating
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_hs_mask.and.q_hr_mask) then
      E_ab =  E_s_m

      dq_col_b = -(rhof(k)*pi/4)*E_ab*n_hs           &
           *q_hr*(dlt_s0*D_hs**2                     &
             +dlt_s1r*D_hs*D_hr+dlt_r1*D_hr**2)      &
           *( th_s0*v_hs**2-th_s1r*v_hr*v_hs         &
             +th_r1*v_hr**2+sigma_hs**2+sigma_hr**2)**0.5

      dn_col_b = -(rhof(k)*pi/4)*E_ab*n_hs           &
           *n_hr*(dlt_s0*D_hs**2                     &
             +dlt_s0r*D_hs*D_hr+dlt_r0*D_hr**2)      &
           *( th_s0*v_hs**2-th_s0r*v_hr*v_hs         &
             +th_r0*v_hr**2+sigma_hs**2+sigma_hr**2)**0.5

      dq_col_a = (rhof(k)*pi/4)*E_ab*q_hs            &
           *n_hr*(dlt_s1*D_hs**2                     &
             +dlt_r1s*D_hs*D_hr+dlt_r0*D_hr**2)      &
           *( th_s1*v_hs**2-th_r1s*v_hr*v_hs         &
             +th_r0*v_hr**2+sigma_hs**2+sigma_hr**2)**0.5

      ! first adjustment following ICON, 2017
      dq_hs_col_rs = max(-dq_col_a,min(0.0,-svm(i,j,k,iq_hs)/delt-q_hsp))
      dq_hr_col_rs = max( dq_col_b,min(0.0,-svm(i,j,k,iq_hr)/delt-q_hrp))

      ! adjustment of numbers - both ice and snow
      rem_cf = (1.0-rem_n_hs_min)/delt
      dn_hs_col_rs = max(dn_col_b,min(0.0,-rem_cf*svm(i,j,k,in_hs)-n_hsp))

      rem_cf = (1.0-rem_n_hr_min)/delt
      dn_hs_col_rs = max(dn_col_b,min(0.0,-rem_cf*svm(i,j,k,in_hr)-n_hrp))

      if (tmp0(i,j,k).lt.T_3) then
        ! and copying it to the second one
        dn_hr_col_rs = dn_hs_col_rs

        ! record the change in cloud ice
        q_hsp = q_hsp + dq_hs_col_rs
        n_hsp = n_hsp + dn_hs_col_rs

        ! change in rain
        q_hrp = q_hrp + dq_hr_col_rs
        n_hrp = n_hrp + dn_hr_col_rs

        ! and for graupel
        q_hgp = q_hgp - dq_hs_col_rs - dq_hr_col_rs
        n_hgp = n_hgp - dn_hs_col_rs

        ! change in th_l - release from freezing and removal of ice
        thlpmcr(i,j,k) = thlpmcr(i,j,k) - (rlme/(cp*exnf(k)))*dq_hr_col_rs
      else  ! tmp0(i,j,k).gt.T_3
        ! enhanced melting and graupel formation
        k_enhm  = c_water/rlme

        ! calculate the melting
        ! with + due to negative value of dq_hr here
        dq_hs_eme_rs = k_enhm*(tmp0(i,j,k)-T_3)*dq_hr_col_rs

        ! snow melting
        dq_hs_eme_rs = max(dq_hs_eme_rs,min(0.0,-svm(i,j,k,iq_hs)/delt-q_hsp))

        ! calculate how many snow perticles melted
        ! q_hs here is always some small positive number
        dn_hs_eme_rs = dq_hs_eme_rs*n_hs/q_hs

        ! limit so it dos not melt more than interacting
        ! also limit so that new graupel not larger that max mean size of source ice ?
        dn_hs_eme_rs = max(dn_hs_eme_rs,dn_hs_col_rs)

        ! update snow
        q_hsp = q_hsp + dq_hs_eme_rs
        n_hsp = n_hsp + dn_hs_eme_rs

        ! no collection of raindrops
        dq_hr_col_rs = 0.0
        dq_hs_col_rs = 0.0
        dn_hr_col_rs = 0.0
        dn_hs_col_rs = 0.0

        ! increase in rain mass
        q_hrp = q_hrp - dq_hs_eme_rs

        ! change in thl - heat spent on melting
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlme/(cp*exnf(k)))*dq_hs_eme_rs
      endif
   endif
  enddo
  enddo
  enddo

end subroutine coll_rsg3


! Ice multiplication
!   - calls separately H-M process for each of the riming processes
! ------------------
subroutine ice_multi3
  implicit none

  integer :: i,j,k
  real    :: dq_ci_spl, dn_ci_spl, dq_hg_temp

  dq_hg_temp = 0.0
  dq_ci_spl  = 0.0
  dn_ci_spl  = 0.0

  ! ------------------------------------
  ! calling separate H-M processes
  ! ------------------------------------

  ! i+l -> i
  call hallet_mossop3 (tmp0(i,j,k),dq_ci_rime,q_ci,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl

  ! s+l -> s
  call hallet_mossop3 (tmp0(i,j,k),dq_hs_rime,q_hs,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl
  dq_ci_mul = dq_ci_mul + dq_ci_spl
  q_hsp     = q_hsp - dq_ci_spl ! effect on snow

  ! g+l -> g
  call hallet_mossop3 (tmp0(i,j,k),dq_hg_rime,q_hg,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl
  dq_ci_mul = dq_ci_mul + dq_ci_spl
  q_hgp     = q_hgp - dq_ci_spl ! effect on snow

  ! g+r -> g
  call hallet_mossop3 (tmp0(i,j,k),dq_hghr_rime,q_hg,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl
  dq_ci_mul = dq_ci_mul + dq_ci_spl
  q_hgp     = q_hgp - dq_ci_spl  ! effect on snow

  ! s+r -> g  -- does it make sense to include it?
  dq_hg_temp = -dq_hs_col_rs - dq_hr_col_rs
  call hallet_mossop3 (tmp0(i,j,k),dq_hg_temp ,q_hg,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl
  dq_ci_mul = dq_ci_mul + dq_ci_spl
  q_hgp     = q_hgp - dq_ci_spl ! effect on snow

  ! i+r-> g   -- does it make sense to include it?
  dq_hg_temp = -dq_ci_col_ri - dq_hr_col_ri
  call hallet_mossop3 (tmp0(i,j,k),dq_hg_temp,q_hg,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl
  dq_ci_mul = dq_ci_mul + dq_ci_spl
  q_hgp     = q_hgp - dq_ci_spl  ! effect on snow

  ! ------------------------------------
  !   update
  ! ------------------------------------
  ! add updates to dq_ci, dn_ci
  n_cip = n_cip + dn_ci_mul
  q_cip = q_cip + dq_ci_mul

end subroutine ice_multi3


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
  use modmicrodata3, only : rhoeps , al_0ice, al_0snow

  implicit none

  real, parameter  :: pi6rhoe = (pi/6.0)*rhoeps
                     ,cc_ci = al_0ice*rhow/rhoeps
                     ,cc_hs = al_0snow*rhow/rhoeps
  integer :: i,j,k

  real :: rem_ci_cf, rem_hs_cf

  ! remain coefficients
  rem_ci_cf = (1.0-rem_n_min_cv)/delt
  rem_hs_cf = (1.0-rem_n_min_cv)/delt

  ! setting values
  dq_ci_cv = 0.0
  dq_hs_cv = 0.0

  ! inner term for ice conversion
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_cl_mask) then
      if (q_ci_mask) then
        if ((D_ci.gt.D_mincv_ci).and.(tmp0(i,j,k).lt.T_3)) then
          ! remain coefficient
          rem_ci_cf = (1.0-rem_n_min_cv)/delt

          dq_ci_cv = cc_ci*(pi6rhoe*D_ci**3 /x_ci-1.0)
          ! ? not to exceed conversion rate 1 ?
          ! G_ci = min(1.0, G_ci)

          dq_ci_cv = -dq_ci_rime/dq_ci_cv
          dq_ci_cv = max(dq_ci_cv,&
                         min(0.0, &
                         -svm(i,j,k,iq_ci)/delt-q_cip))  ! based on ICON, 2017
          ! = max(dq_ci_cv,-svm(i,j,k,iq_ci)/delt)

          dn_ci_cv = dq_ci_cv/max(x_ci,x_ci_cvmin)
          ! = dq_ci_cv/max(x_ci,x_ci_cvmin)

          dn_ci_cv = max(dn_ci_cv,&
                         min(0.0,        &
                         -rem_ci_cf*svm(i,j,k,in_ci)-n_cip(i,j,k)))
          ! = max(dn_hs_cv, -svm(i,j,k,in_ci)/delt)

          ! change in the amount of graupel
          n_hgp = n_hgp - dn_ci_cv
          q_hgp = q_hgp - dq_ci_cv

          ! change in the amount of cloud ice
          n_cip = n_cip + dn_ci_cv
          q_cip = q_cip + dq_ci_cv
        endif

      ! term for snow conversion
      if (q_hs_mask) then
        if ((D_hs.gt.D_mincv_hs).and.(tmp0(i,j,k).lt.T_3)) then
          ! remain coefficient
          rem_hs_cf = (1.0-rem_n_min_cv)/delt

          dq_hs_cv = cc_hs*(pi6rhoe*D_hs**3 /x_hs-1)
          ! ? at the sam time, the value should be limited
          ! ? not to exceed conversion rate 1
          ! G_hs = min(1.0, G_hs)
          dq_hs_cv = -dq_hs_rime/dq_hs_cv

          ! correction - not removing more than available
          ! dq_hs_cv = max(dq_hs_cv,-q_hs/delt )

          ! basic correction of the tendency
          dq_hs_cv = max(dq_hs_cv,&
                         min(0.0, &
                         -svm(i,j,k,iq_hs)/delt-q_hsp))
          ! dq_hs_cv = max( dq_hs_cv,-svm(i,j,k,iq_hs)/delt)

          dn_hs_cv = dq_hs_cv/max(x_hs,x_hs_cvmin)
          ! dn_hs_cv = dq_hs_cv/max(x_hs,x_hs_cvmin)

          dn_hs_cv = max(dn_hs_cv,min(0.0,-rem_hs_cf*svm(i,j,k,in_hs)-n_hsp))
          ! dn_hs_cv = max(dn_hs_cv, -svm(i,j,k,in_hs)/delt)

          ! and the second correction of the q tendency
          ! change in the amount of graupel
          n_hgp = n_hgp - dn_hs_cv
          q_hgp = q_hgp - dq_hs_cv

          ! change in the amount of snow
          n_hsp = n_hsp + dn_hs_cv
          q_hsp = q_hsp + dq_hs_cv
        endif
      endif

      ! warnings
      if (l_sb_dbg) then
        if(svm(i,j,k,iq_ci)+delt*dq_ci_cv.lt. 0) then
          write(6,*) 'WARNING: conv_partial3 removing too much ice'
          write(6,*) ' removing more ice than available'
          ! count(svm(i,j,k,iq_ci)+delt*dq_ci_cv.lt. 0)
          write(6,*) ' removing too much ice'
          ! count(q_ci+delt*dq_ci_cv.lt. 0)
          write(6,*) ' getting negative q_t'
          ! count(qt0(i,j,k)+delt*dq_ci_cv.lt. 0)
        endif
        if(svm(i,j,k,iq_hs)+delt*dq_hs_cv.lt. 0) then
          write(6,*) 'WARNING: conv_partial3 removing too much snow'
          write(6,*) ' removing more snow than available'
          ! count(svm(i,j,k,iq_hs)+delt*dq_hs_cv).lt. 0)
          write(6,*) ' removing too much snow'
          ! count(q_hs+delt*dq_hs_cv.lt. 0)
        endif
      endif
    endif
  enddo
  enddo
  enddo
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
  use modglobal, only : rlv,cp
  use modfields, only : exnf,svm
  implicit none

  integer :: i,j,k

  ! ------------------------------------
  ! calling separate melting processes
  ! ------------------------------------

  ! ice
  call sb_evmelt3(aven_0i,aven_1i,bven_0i,bven_1i,x_ci_bmin        &
                 ,n_ci,q_ci,x_ci,D_ci,v_ci,dq_ci_me,dn_ci_me,dq_ci_ev,dn_ci_ev)

  ! loss of ice
  n_cip  = n_cip + dn_ci_me + dn_ci_ev
  q_cip  = q_cip + dq_ci_me + dq_ci_ev

  ! transformed to rain or cloud

  ! for now into rain - improve possibly later
  n_hrp = n_hrp - dn_ci_me
  q_hrp = q_hrp - dq_ci_me

  ! transfomed to water vapour
  qtpmcr(i,j,k) = qtpmcr(i,j,k) - dq_ci_ev
  ! and heat production : heat spent on melting and evaporation - done lower

  ! snow
  call sb_evmelt3(aven_0s,aven_1s,bven_0s,bven_1s,x_hs_bmin        &
                 ,n_hs,q_hs,x_hs,D_hs,v_hs,dq_hs_me,dn_hs_me,dq_hs_ev,dn_hs_ev)

  ! loss of snow
  n_hsp = n_hsp + dn_hs_me +dn_hs_ev
  q_hsp = q_hsp + dq_hs_me +dq_hs_ev

  ! transformed to rain
  n_hrp = n_hrp - dn_hs_me
  q_hrp = q_hrp - dq_hs_me

  ! transfomed to water vapour
  qtpmcr(i,j,k) = qtpmcr(i,j,k)-dq_hs_ev
  ! and heat production : heat spent on melting and evaporation - done lower

  ! graupel
  call sb_evmelt3(aven_0g,aven_1g,bven_0g,bven_1g,x_hg_bmin         &
                 ,n_hg,q_hg,x_hg,D_hg,v_hg,dq_hg_me,dn_hg_me,dq_hg_ev,dn_hg_ev)

  ! limiting not to remove more water than available
  if (dn_hg_me.lt.0.0) then
    dn_hg_ev = max(dn_hg_ev,-svm(i,j,k,in_hg)/delt-n_hgp)
    dq_hg_ev = max(dq_hg_ev,-svm(i,j,k,iq_hg)/delt-q_hgp)
    dn_hg_me = max(dn_hg_me,-svm(i,j,k,in_hg)/delt-n_hgp)
    dq_hg_me = max(dq_hg_me,-svm(i,j,k,iq_hg)/delt-q_hgp)
  endif

  ! loss of graupel
  n_hgp = n_hgp + dn_hg_me + dn_hg_ev
  q_hgp = q_hgp + dq_hg_me + dq_hg_ev

  ! transformed to rain
  n_hrp = n_hrp - dn_hg_me
  q_hrp = q_hrp - dq_hg_me

  ! transfomed to water vapour
  qtpmcr(i,j,k) = qtpmcr(i,j,k)-dq_hg_ev

  ! and heat production : heat spent on melting and evaporation - done lower
  !  - melting goes to rain
  !  - evaporation goes water vapour
  thlpmcr(i,j,k) = thlpmcr(i,j,k)+                             &
            (rlme/(cp*exnf(k)))*(dq_ci_me+dq_hs_me+dq_hg_me) + &
      ((rlv+rlme)/(cp*exnf(k)))*(dq_ci_ev+dq_hs_ev+dq_hg_ev)
end subroutine evapmelting3


!> Determine autoconversion rate and adjust qrp and Nrp accordingly
!!
!!   The autoconversion rate is formulated for f(x)=A*x**(nuc)*exp(-Bx),
!!   decaying exponentially for droplet mass x.
!!   It can easily be reformulated for f(x)=A*x**(nuc)*exp(-Bx**(mu)) and
!!   by chosing mu=1/3 one would get a gamma distribution in drop diameter
!!   -> faster rain formation. (Seifert)
subroutine autoconversion3
  use modglobal, only : i1,j1,k1,kmax,rlv,cp
  use modfields, only : exnf,rhof,svm
  implicit none

  integer :: i,j,k
  real :: rem_cf
  real :: tau, phi, nuc

  if (l_sb) then
    !
    ! SB autoconversion
    !
    if (l_sb_classic) then  ! l_sb_classic - ie. S&B version

      ! autoconversion coefficient
      k_au = k_cc/(20.0*x_s)

      ! remain coefficient
      rem_cf = (1.0-rem_n_cl_min)/delt

      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (q_cl_mask) then
          dq_hr_au = k_au * ( nu_cl_cst+2.) * ( nu_cl_cst +4.) / ( nu_cl_cst +1.)**2.    &
                  * (q_cl * x_cl)**2. * rho0s ! *rho**2/rho/rho (= 1)
          tau = 1.0 - q_cl / qltot

          ! phi_au computation
          phi = k_1 * tau**k_2 * (1.0 -tau**k_2)**3
          dq_hr_au = dq_hr_au * (1.0 + phi/(1.0 -tau)**2)

          ! basic au correction
          dq_hr_au  = min(dq_hr_au,max(0.0,svm(i,j,k,iq_cl)/delt+q_clp))

          ! and calculate cloud droplet numbers
          dn_cl_au = (-2.0/x_s)*dq_hr_au                                   ! and the droplet number
          dn_cl_au = max(dn_cl_au,min(0.0,-rem_cf*svm(i,j,k,in_cl)-n_clp)) ! correct droplet number
          dn_hr_au = -0.5*dn_cl_au                                         ! and the raindrop number

          ! and adding updates
          q_hrp = q_hrp + dq_hr_au
          n_hrp = n_hrp + dn_hr_au
          q_clp = q_clp - dq_hr_au
          n_clp = n_clp + dn_cl_au

          thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*dq_hr_au
          qtpmcr (i,j,k) = qtpmcr (i,j,k) - dq_hr_au
        endif
      enddo
      enddo
      enddo
    else ! l_sb_classic = .false. - original version in DALES

      k_au = k_c/(20.0*x_s)

      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (q_cl_mask) then
          nuc = k1nuc*(rhof(k)*q_cl*1000.) +k2nuc-1. ! #sb3 G09a
          dq_hr_au = k_au * (nuc+2.) * (nuc+4.) / (nuc+1.)**2.    &
                  * (q_cl * x_cl)**2. * rho0s ! *rho**2/rho/rho (= 1)
          tau            = 1.0 - q_cl / qltot ! #sb3
          phi = k_1 * tau**k_2 * (1.0 -tau**k_2)**3   ! phi_au computation
          dq_hr_au = dq_hr_au * (1.0 + phi/(1.0 -tau)**2)

          ! cloud water numbers
          dn_hr_au = dq_hr_au/x_s
          dn_cl_au = (-2.0/x_s)*dq_hr_au

          ! #sb3 START outputs
          q_hrp = q_hrp + dq_hr_au
          n_hrp = n_hrp + dn_hr_au
          q_clp = q_clp - dq_hr_au
          n_clp = n_clp + dn_cl_au  ! o:  n_clp - (2.0/x_s)*rhof(k)*au

          thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*dq_hr_au
          qtpmcr (i,j,k) = qtpmcr (i,j,k) - dq_hr_au
        endif
      enddo
      enddo
      enddo
    endif ! l_sb_classic

  endif ! l_sb

  if (l_sb_dbg) then
    if (any(q_cl/delt - dq_hr_au.lt. 0.)) then
      write(6,*) 'WARNING: autoconversion too high'
      write(6,*) '  removing more cloud water than available in ', count(q_cl/delt - dq_hr_au.lt. 0.)
    end if

    if (any(n_cl(2:i1,2:j1,1:kmax)/delt -(2.0/x_s)*dq_hr_au.lt. 0.)) then
      write(6,*) 'WARNING: autoconversion too high'
      write(6,*) '  removing more droplets than available in' , count(n_cl(2:i1,2:j1,1:kmax)/delt - (2.0/x_s)*dq_hr_au.lt. 0.)
    end if
  endif

end subroutine autoconversion3


!> Self-collection of cloud droplets
!! written based on S&B to evaluate the self-collection of cloud droplets
subroutine cloud_self3
  use modmicrodata3, only : k_cc,rho0s
  use modglobal, only : i1,j1,k1
  use modfields, only : svm
  implicit none

  real, parameter :: k_clsc = - k_cc*rho0s

  real    :: rem_cf
  integer :: i,j,k

  ! initialize
  dn_cl_sc = 0.0

  ! calculate constant
  rem_cf = (1.0-rem_n_cl_min)/delt

  ! loop sb
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_cl_mask) then
      dn_cl_sc = k_clsc*(q_cl(i,j,k)**2)* &
        (nu_cl_cst+2.0)/(nu_cl_cst+1.0)-dn_cl_au(i,j,k)

      ! basic sc collection
      dn_cl_sc = min(0.0,dn_cl_sc)
      dn_cl_sc = max(dn_cl_sc,min(0.0,-rem_cf*svm(i,j,k,in_cl)-n_clp(i,j,k)))

      ! update
      n_clp(i,j,k) = n_clp(i,j,k)+dn_cl_sc

      ! no change to q_cl

      if (l_sb_dbg) then
        ! testing for too high values
        if (svm(i,j,k,in_cl)/delt - dn_cl_sc.lt. 0.) then
          write(6,*)'WARNING: cloud sc too high'
          write(6,*) '  getting to negative n_cl'
          ! count(svm(i,j,k,in_cl)/delt -dn_cl_sc).lt. 0.)
        endif
      endif
    endif
  enddo
  enddo
  enddo
end subroutine cloud_self3


!*********************************************************************
! determine accr. + self coll. + br-up rate and adjust qrp and Nrp
! base don S&B
!*********************************************************************
subroutine accretion3
  use modglobal, only : i1,j1,k1,rlv,cp
  use modfields, only : exnf,rhof,ql0,svm
  implicit none

  real :: phi_br, Dvrf, tau
  real :: rem_cf
  integer :: i,j,k

  if (l_sb) then
    !
    ! SB accretion
    !
    ! #t write(6,*) 'l_sb accretion...'
    if (l_sb_classic) then
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (q_hr_mask.and.q_cl_mask) then
          ! since it is forming only where rain
          tau = 1.0 - q_cl/qltot
          phi = (tau/(tau + k_l))**4.
          dq_hr_ac = k_cr *rhof(k)*q_cl*q_hr * phi * &
                           (rho0s/rhof(k))**0.5  ! rho*rho / rho  = rho

          ! basic ac correction
          dq_hr_ac = min(dq_hr_ac,q_cl/delt) ! min(dq_hr_ac,svm(i,j,k,iq_cl)/delt)

          ! number of cloud droplets
          ! remain coefficient for clouds #sb3
          rem_cf = (1.0-rem_n_cl_min)/delt

          dn_cl_ac = -dq_hr_ac/x_cl
          dn_cl_ac = max(dn_cl_ac,&
                         min(0.0, &
                         -rem_cf*svm(i,j,k,in_cl)-n_clp))

          ! update
          q_hrp = q_hrp + dq_hr_ac

          ! no change n_hrp
          ! and changes in water number for clouds
          q_clp = q_clp - dq_hr_ac
          n_clp = n_clp + dn_cl_ac !o: n_clp - rhof(k)*ac/xc

          qtpmcr (i,j,k) = qtpmcr (i,j,k) - dq_hr_ac
          thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*dq_hr_ac
        endif
      enddo
      enddo
      enddo
    else ! l_sb_classic
      !*********************************************************************
      ! determine accr. rate and adjust qrp and Nrp
      ! accordingly. Break-up : Seifert (2007), a
      !*********************************************************************
      tau = 1.0 - ql0(i,j,k)/qltot
      phi = (tau/(tau + k_l))**4.

      dq_hr_ac = k_r *rhof(k)*ql0(i,j,k) * q_hr * phi * &
                  (1.225/rhof(k))**0.5
      q_hrp = q_hrp + dq_hr_ac

      ! no change n_hrp
      q_clp  = q_clp - dq_hr_ac
      qtpmcr (i,j,k) = qtpmcr (i,j,k) - dq_hr_ac
      thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*dq_hr_ac

      ! and update on cloud water number
      dn_cl_ac = -dq_hr_ac/x_cl
      n_clp = n_clp + dn_cl_ac
    endif

    !
    ! SB self-collection & Break-up
    !
    if (l_sb_classic) then
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (q_hr_mask) then
          dn_hr_sc = -k_rr *rhof(k)* q_hr * n_hr  &
            * (1.0 + kappa_r/lbdr)**(-9.)*(rho0s/rhof(k))**0.5

          ! and calculating size of droplets - adjusted
          Dvrf = Dvr ! for now leaving the same
          if (Dvrf.gt.dvrlim) then
            if (Dvrf.gt.dvrbiglim) then
              ! for big drops
              phi_br = 2.0*exp(kappa_br*(Dvrf-D_eq))-1.0
            else
              ! for smaller drops
              phi_br = k_br * (Dvrf-D_eq)
            endif
            dn_hr_br = -(phi_br + 1.) * dn_hr_sc
          else
            dn_hr_br = 0. ! (phi_br = -1)
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
        if (q_hr_mask) then
          dn_hr_sc = -k_rr *rhof(k)* q_hr * n_hr  &
            * (1.0 + kappa_r/lbdr*pirhow**(1./3.))**(-9.)*(rho0s/rhof(k))**0.5
        endif
        if ((Dvr .gt. dvrlim) .and. q_hr_mask) then
          if (Dvr .gt. dvrbiglim) then
            ! for big drops
            phi_br = 2.0*exp(kappa_br*(Dvr-D_eq))-1.0
          else
            ! for smaller drops
            phi_br = k_br* (Dvr-D_eq)
          endif
          dn_hr_br = -(phi_br + 1.) * dn_hr_sc
        else
          dn_hr_br = 0. ! (phi_br = -1)
        endif
      enddo
      enddo
      enddo
  endif ! l_sb_classic
  endif

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_hr_mask) then
      ! correcting for low values
      ! calculating coef for ratio for minimal remaing number of droplets #sb3
      rem_cf = (1.0-rem_n_hr_min)/delt

      n_hrp = n_hrp +max(min(0.0,                             &
                             -rem_cf*svm(i,j,k,in_hr)-n_hrp), &
                         dn_hr_sc+dn_hr_br)
    endif
  enddo
  enddo
  enddo

  if (l_sb_dbg) then
    if(svm(i,j,k,iq_cl)/delt-dq_hr_ac.lt. 0.) then
      write(6,*) 'WARNING: accretion removing too much water'
      write(6,*) ' updated below 0'
      ! count(svm(i,j,k,iq_cl)/delt-dq_hr_ac.lt. 0.0)
    endif

    if(svm(i,j,k,in_hr)/delt+dn_hr_sc+dn_hr_br.lt. 0.) then
      write(6,*) 'WARNING: self-collection of rain too high'
      write(6,*) ' removing more n_hr than available'
      ! count(svm(i,j,k,in_hr)/delt-dn_hr_sc+dn_hr_br.lt. 0.0)
      write(6,*) ' n_hr updated below 0'
      ! count(n_hr/delt-dn_hr_sc+dn_hr_br.lt. 0.0)
    endif

    if (svm(i,j,k,in_cl)/delt - dq_hr_ac/x_cl.lt. 0.) then
      write(6,*)'WARNING: ac too large, removing too many droplets'
      write(6,*)'   in'
      ! count(svm(i,j,k,in_cl)/delt-dq_hr_ac/x_cl.lt. 0.)
    endif
  endif ! l_sb_dbg

end subroutine accretion3


!*********************************************************************
! Evaporation of prec. : Seifert & Beheng
! Cond. (S>0.) neglected (all water is condensed on cloud droplets)
!*********************************************************************
subroutine evap_rain3
  use modglobal, only : i1,j1,k1,rv,rlv,cp,pi
  use modfields, only : exnf,qt0,svm,qvsl,tmp0,rhof
  implicit none

  integer :: i,j,k

  real :: f0    & ! ventilation factor - moment 0
         ,f1    & ! ventilation factor - moment 1
         ,S     & ! super or undersaturation
         ,G     & ! cond/evap rate of a drop
         ,vihr  & ! mean terminal velocity
         ,nrex  & ! Reynolds number N_re(xr)
         ,x_hrf   ! full x_hr without bounds

  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (q_hr_mask) then
      ! adjusting the calculation for saturation
      S = min(0.,((qt0(i,j,k)-q_cl)/qvsl(i,j,k)- 1.0))
      G = (rv * tmp0(i,j,k)) / (Dv*esl(i,j,k)) &
        + rlv/(Kt*tmp0(i,j,k))*(rlv/(rv*tmp0(i,j,k)) -1.)
      G = 1./G

      ! terminal velocity  (from mixed scheme)
      ! xr in this calculation !!
      vihr = al_hr*((rho0s/rhof(k))**0.5)*x_hr(i,j,k)**be_hr

      ! calculating  N_re Reynolds number
      nrex = Dvr*vihr/nu_a
      x_hrf = q_hr/(n_hr+eps0)

      f0 = aven_0r+bven_0r*Sc_num**(1.0/3.0)*nrex**0.5
      f1 = aven_1r+bven_1r*Sc_num**(1.0/3.0)*nrex**0.5
      f0 = max(0.0,f0)

      dq_hr_ev = 2*pi*n_hr*G*Dvr*f1*S

      dn_hr_ev = 2*pi*n_hr*G*Dvr*f0*S/x_hrf ! x_hr here, not xr

      ! and limiting it
      dn_hr_ev = min(dn_hr_ev, 0.0)
      dn_hr_ev = max(dn_hr_ev,dq_hr_ev/x_hrf)

      if ((dq_hr_ev+svm(i,j,k,iq_hr)/delt.lt.0).or.&
          (dn_hr_ev+svm(i,j,k,in_hr)/delt.lt.0)) then
        dn_hr_ev = - svm(i,j,k,in_hr)/delt
        dq_hr_ev = - svm(i,j,k,iq_hr)/delt
      endif

      q_hrp = q_hrp + dq_hr_ev
      n_hrp = n_hrp + dn_hr_ev
      qtpmcr(i,j,k)  = qtpmcr(i,j,k) -dq_hr_ev
      thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*dq_hr_ev

      ! recovery of aerosols ?
      ret_cc = ret_cc - c_ccn_ev_r*min(0.0,dn_hr_ev)

    endif
  enddo
  enddo
  enddo
end subroutine evap_rain3



! TOOD: removed warnings when throwing away too much negative
! moisture

  subroutine point_processes(i,j,k)
    use modfields, only : sv0, svm, svp
    implicit none
    integer, intent(in) :: i,j,k

    ! aerosols:
    n_cc = max(sv0(i,j,k,in_cc),0.)

    ! reading number densities:
    n_cl = max(sv0(i,j,k,in_cl),0.)
    n_ci = max(sv0(i,j,k,in_ci),0.)
    n_hr = max(sv0(i,j,k,in_hr),0.)
    n_hs = max(sv0(i,j,k,in_hs),0.)
    n_hg = max(sv0(i,j,k,in_hg),0.)

    ! reading mass densities
    q_cl = max(sv0(i,j,k,iq_cl),0.)
    q_ci = max(sv0(i,j,k,iq_ci),0.)
    q_hr = max(sv0(i,j,k,iq_hr),0.)
    q_hs = max(sv0(i,j,k,iq_hs),0.)
    q_hg = max(sv0(i,j,k,iq_hg),0.)

    ! Reset all values
    ! ----------------

    ! 0 bulk integrals
    xc = 0.
    Dvr = 0.
    lbdr = 0.

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
    dq_cl_sa      = 0.0
    dn_cl_sa      = 0.0

    delt = rdt / (4. - dble(rk3step))  ! TODO: get from modbulkmicro3


    ! Testing for a noticable amount of rain graupel and snow
    ! -------------------------------------------------------

    ! rain :
    q_hr_mask = (q_hr.gt.q_hr_min).and.(n_hr.gt.0.0)

    ! snow :
    q_hs_mask = (q_hs.gt.q_hs_min).and.(n_hs.gt.0.0)

    ! graupel :
    q_hg_mask = (q_hg.gt.q_hg_min).and.(n_hg.gt.0.0)


    ! calculate qltot and initialize cloud droplet number Nc
    ! ------------------------------------------------------
    ! instead of: ql0  (i,j,k) + q_hr
    ! i.e - cloud water + rain areas
    qltot = q_cl + q_hr

    ! clouds- original definition
    qcmask = (ql0.gt.qcmin)

    ! liquid clouds
    q_cl_mask = (q_cl(i,j,k).ge.qcliqmin).and.(n_cl(i,j,k).gt.0.0)

    ! ice clouds

    q_ci_mask = (qci .ge. qicemin).and.(n_ci(i,j,k).gt.0.0)


    ! calculate Rain DSD integral properties & parameters xr, Dvr, lbdr
    ! -----------------------------------------------------------------
    call integrals_bulk3
    call nucleation3      ! cloud nucleation
    call icenucle3        ! ice nucleation


    ! freezing of water droplets
    ! -----------------------------------------------------------------
    call homfreez3        ! homogeneous freezing of cloud droplets
    call hetfreez3        ! heterogeneous freezing


    ! deposition processes
    ! -----------------------------------------------------------------
    call deposit_ice3      ! deposition of vapour to cloud ice
    call deposit_snow3     ! deposition of vapour to snow
    call deposit_graupel3  ! deposition of vapour to graupel
    call cor_deposit3      ! correction for deposition


    ! snow aggregation and self-collection
    ! -----------------------------------------------------------------
    call ice_aggr3         ! ice selfcollection
    call snow_self3        ! snow selfcollection - tendency only


    ! collision processes for snow
    ! -----------------------------------------------------------------
    call coll_sis3         ! snow selfcollection s+i -> s
    call coll_gsg3         ! snow selfcollection g+s -> g
    call coll_ici3         ! riming i+c -> i
    call coll_scs3         ! riming s+c -> s
    call coll_gcg3         ! riming g+c -> g
    call coll_grg3         ! riming g+r -> g


    ! raindrop freezing
    ! -----------------------------------------------------------------
    call rainhetfreez3     ! heterogeneous freezing


    ! collision with conversion
    ! -----------------------------------------------------------------
    call coll_rig3         ! riming r+i -> g  !#t2 !#t4
    call coll_rsg3         ! riming r+i -> g  !#t7


    ! conversions
    ! -----------------------------------------------------------------
    call ice_multi3        ! ice multiplication of Hallet and Mossop

    if (l_sb_conv_par) then
      call conv_partial3 ! partial conversion
    endif


    ! melting and evaporation of ice particles
    ! -----------------------------------------------------------------
    call evapmelting3      ! melting of ice particles


    ! basic warm processes
    ! -----------------------------------------------------------------
    call autoconversion3
    call cloud_self3
    call accretion3
    call evap_rain3        ! rain evaporation

    ! -----------------------------------------------------------------
end subroutine point_processes


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
! Ice multiplication of Hallet and Mossop (1974)
!
! - written as described in Seifert (2002)
! - implementation similar to the one in ICON model
!
! - returns:
!    dq_i_hm  - change in moass content during the process
!    dn_i_hm  - number of newly produced ice particles
!
! ***************************************************************
subroutine hallet_mossop3(t0,dq_rime,q_e,dq_i_hm,dn_i_hm)
  use modmicrodata3, only : c_spl_hm74, rem_q_e_hm, tmp_opt_hm74 &
                           ,tmp_min_hm74,tmp_max_hm74

  implicit none

  ! inputs
  real, intent(in)  :: t0        ! tmp0 at this gridpoint
  real, intent(in)  :: dq_rime   ! riming rate
  real, intent(in)  :: q_e       ! amount of that ice phase

  ! outputs
  real, intent(out) :: dq_i_hm   ! tendency in q_i by H-M process
  real, intent(out) :: dn_i_hm   ! tendency in n_i by H-M process

  ! local variables

! constants in calculation
  real, parameter :: c_spl  = c_spl_hm74 &
                    ,c_1_hm = 1.0/(tmp_opt_hm74-tmp_min_hm74) &
                    ,c_2_hm = 1.0/(tmp_opt_hm74-tmp_max_hm74)

  real :: mult_1,mult_2,mint_1,mint_2    &  ! calculation variables
         ,dn_try,dq_try,rem_cf              ! trial variables

  ! setting coefficient for reminder
  rem_cf  = (1.0-rem_q_e_hm)/delt

  ! only if riming going on temperature below 0
  if ((dq_rime.gt.0).and.(t0.lt.T_3)) then

     ! f_spl calculation following ICON
     mult_1 = c_1_hm * (t0-tmp_min_hm74)   ! positive in the target interval
     mult_2 = c_2_hm * (t0-tmp_max_hm74)   ! positive in the target interval

     ! now for intervals
     mint_1 = max(0.0,min(1.0,mult_1))     !  0 for T<T_min, 1 for T>T_opt
     mint_2 = max(0.0,min(1.0,mult_2))     !  0 for T>T_max, 1 for T<T_opt

     ! calculating prediction for the process
     dn_try = c_spl*mint_1*mint_2*dq_rime
     dq_try = x_ci_spl*dn_try

     ! correcting
     dq_try = min(dq_try,rem_cf*q_e+dq_rime)    !  limit splintering
     dq_try = max(dq_try,0.0)

     ! prepare updates
     dq_i_hm = dq_try
     dn_i_hm = dq_try/x_ci_spl
     !
  else
     dq_i_hm = 0.0
     dn_i_hm = 0.0
  endif
end subroutine hallet_mossop3

end module modbulkmicro_point
