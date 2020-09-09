module modbulkmicro_column

  implicit none
contains


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
  use modfields, only : rhof,qt0,svm,svp,qvsl
  implicit none

  integer :: i,j,k
  real :: ntest, n_bmax, cogr_max, ql_res

  ! initialisation of process updates
  dq_cl_sa = 0.0
  dn_cl_sa = 0.0

  ! calculation
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if (l_sb_all_or) then
      ! note: threshold might be lower then threshold for cloud computations, obviously
      ntest=svm(i,j,k,in_cl)+n_clp*delt
      if (ntest.gt.0.0) then
        !
        ! remaining water =
        !    + condesable water available
        !    - already condensed
        !    - newly condendsed
        !    - removed by mphys processed
        !
        !  calculating amount of available water
        ql_res=(qt0(i,j,k)-qvsl(i,j,k)-svm(i,j,k,iq_cl))+delt*(qtpmcr(i,j,k)-q_clp)
        dq_cl_sa = ql_res/delt
        if ((svm(i,j,k,iq_cl)+delt*(q_clp+dq_cl_sa)).lt.0.0) then
          dn_cl_sa = -svm(i,j,k,in_cl)/delt-n_clp
          dq_cl_sa = -svm(i,j,k,iq_cl)/delt-q_clp
        endif
      endif
    else ! l_sb_all_or
      if (l_sb_dumpall) then
        ! dump all water
        ! note: threshold might be lower then threshold for cloud computations, obviously
        if (q_cl_mask) then
          ntest = svm(i,j,k,in_cl)+n_clp*delt ! #t1
          if (ntest.gt.0.0) then
            !
            ! remaining water =
            !    + condesable water available
            !    - already condensed
            !    - newly condendsed
            !    - removed by mphys processed
            !

            !  calculating amount of available water
            ql_res=(qt0(i,j,k)-qvsl(i,j,k)-svm(i,j,k,iq_cl))+delt*(qtpmcr(i,j,k)-q_clp)

            ! limiting so it does not remove more water than in clouds
            dq_cl_sa = max((-svm(i,j,k,iq_cl)/delt)-q_clp,ql_res/delt)

            ! adjusting number of cloud droplets
            ! - calculate min size with this amount of water
            n_bmax = (svm(i,j,k,iq_cl)+delt*(q_clp+dq_cl_sa))/(0.1*x_cl_bmin)
            n_bmax = max(n_bmax, 0.0)

            ! of course we do not want negative values - but that is alread sorted above
            ! - remove droplets so that mean size in not less than
            dn_cl_sa = min(0.0, n_bmax-ntest)/delt

            ! limit change so not in negative numbers
            dn_cl_sa = max((-svm(i,j,k,in_cl)/delt)-svp(i,j,k,in_cl)-n_clp,dn_cl_sa)
          endif
        endif
      else ! l_sb_dumpall
        ! note: threshold might be lower then threshold for cloud computations, obviously
        if (q_cl_mask) then
          ntest=svm(i,j,k,in_cl)+n_clp*delt !#t1
          if (ntest.gt.0.0) then
            !
            ! and now if we want to enforce limit on cloud droplet size
            !
            ql_res=(qt0(i,j,k)-qvsl(i,j,k)-svm(i,j,k,iq_cl))+delt*(qtpmcr(i,j,k)-q_clp)

            ! calculate maximal available update
            cogr_max =(1.0/delt)*(x_cogr_max*(svm(i,j,k,in_cl)+delt*n_clp)-svm(i,j,k,iq_cl))

            ! dump just what is below max size by condensation growth
            dq_cl_sa = min(cogr_max-q_clp,ql_res/delt)

            ! ie. either whole amount, or only that much that droplet size will be: xc_cogr_max
            ! other possibility: require it to be larger than 0
            ! and prevent negative values of svm + delt *svp
            dq_cl_sa = max((-svm(i,j,k,iq_cl)/delt)-q_clp,dq_cl_sa)

            ! adjusting number of cloud droplets
            ! - calculate min size with this amount of water
            n_bmax = (svm(i,j,k,iq_cl)+delt*(q_clp+dq_cl_sa))/(0.5*x_cl_bmin)
            n_bmax = max(n_bmax, 0.0)

            ! - remove droplets so that mean size in not by order of magnitude less than x_cl_bmin
            dn_cl_sa = min(0.0, n_bmax-ntest)/delt

            ! limit change so not in negative numbers
            dn_cl_sa = max((-svm(i,j,k,in_cl)/delt)-svp(i,j,k,in_cl)-n_clp,dn_cl_sa)

            !-------------------------------------
            !! cloud water mixing ratio
            ! + mphys changes in cloud water
            ! + remaining water:
            !    + condesable water available
            !    - already condensed
            !    - newly condendsed
            !      ( {change in cloud water} = {newly condendsed or deposited}
            !                                         {removed by cloud processes}  )
            !      ( - {newly condendsed or deposited} =
            !          - ({change in liquid cloud water}+{change in ice cloud water})
            !          + {removed by cloud processes}
            !      )
            ! -----------------------------
          endif
        endif
      endif
    endif ! l_sb_all_or

    ! and update
    q_clp  = q_clp + dq_cl_sa
    n_clp  = n_clp + dn_cl_sa
  enddo
  enddo
  enddo
end subroutine satadj3


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
            wfall_qr        = max(0.,((rho0s/rhof(k))**0.5)*(a_tvsbc           &
                      -b_tvsbc*(1.+c_tvsbc/lbdr_spl(i,j,k))**(-4.0))) ! k=1
            wfall_Nr        = max(0.,((rho0s/rhof(k))**0.5)*(a_tvsbc           &
                      -b_tvsbc*(1.+c_tvsbc/lbdr_spl(i,j,k))**(-1.0))) ! k=0
            sed_qr  (i,j,k) = wfall_qr       *qr_spl(i,j,k)*rhof(k)
            sed_Nr  (i,j,k) = wfall_Nr       *Nr_spl(i,j,k)*rhof(k)
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
             wfall_qr        = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl(i,j,k))**(-1.*(mur_spl(i,j,k)+4.))))
             wfall_Nr        = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl(i,j,k))**(-1.*(mur_spl(i,j,k)+1.))))
             sed_qr  (i,j,k) = wfall_qr       *qr_spl(i,j,k)*rhof(k)
             sed_Nr  (i,j,k) = wfall_Nr       *Nr_spl(i,j,k)*rhof(k)
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


! ***************************************************************
!    recovery of ccn
!
!   - to be later replaced based on advance literature
! ***************************************************************
subroutine recover_cc
  use modglobal, only : ih,i1,j1,k1
  implicit none

  integer :: i,j,k

  !   - to be later replaced based on advance literature
  !
  !   - so far just and easy recovery
  !     of ccn based on number of water particles that evaporated, sublimated
  !     or got removed with remaining positive n_
  do k=1,k1
  do j=2,j1
  do i=2,i1
    if(.not.(l_c_ccn)) then ! ie. no change if constant ccn
        ! decrease in total amount of potential CCN
        n_ccp = n_ccp    &
              + dn_cl_sc &
              + dn_cl_se &
              + dn_cl_au &
              + dn_cl_ac &
              + dn_cl_hom  &
              + dn_cl_het  &
              + dn_cl_rime_ci &
              + dn_cl_rime_hs &
              + dn_cl_rime_hg

        ! recovery of potential CCN
        n_ccp = n_ccp + c_rec_cc*ret_cc
    endif
  enddo
  enddo
  enddo
end subroutine recover_cc


subroutine column_processes(i,j)
  implicit none
  integer, intent(in) :: i,j

    ! -----------------------------------------------------------------
    ! column processes:
    ! saturation adjustment
    ! sedimentation part
    ! recovery of ccn aerosols
    ! remove negative values and non physical low values svp
    ! integrate svp = svp + ..._p
    ! prognostic variables
    ! thlp(i,j,k) = thlp(i,j,k)+thlpmcr(i,j,k)
    ! qtp(i,j,k) = qtp(i,j,k)+qtpmcr(i,j,k)
    ! -----------------------------------------------------------------
    !call satadj3
    !call sedim_rain3 !jja sed_qr sed_Nr
    !call sedim_cl3 !jja from routine local, sed_qip, sed_nip
    !call sedim_ice3 !jja from routine local, sed_qip, sed_nip
    !call sedim_snow3 !jja from routine local, sed_qip, sed_nip
    !call sedim_graupel3 !#b2t3 !jja from routine local, sed_qip, sed_nip
    !call recover_cc
    ! -----------------------------------------------------------------

end subroutine column_processes

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

end module modbulkmicro_column
