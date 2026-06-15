module flux_atmocn_bulk_mod
  !-------------------------------------------------------------------------------
  ! PURPOSE:
  !   computes atm/ocn surface fluxes using Large and Pond
  !
  ! NOTES:
  !   o all fluxes are positive downward
  !   o net heat flux = net sw + lw up + lw down + sen + lat
  !   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
  !   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
  !
  ! ASSUMPTIONS:
  !  Large:
  !   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
  !   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
  !                                 ctn = .0180 sqrt(cdn), stable
  !   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
  !   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
  !-------------------------------------------------------------------------------

  use ESMF,         only: ESMF_FAILURE, ESMF_SUCCESS
  use ufs_kind_mod, only: R8=>SHR_KIND_R8, IN=>SHR_KIND_IN
  use ufs_flux_mod, only: alpha, maxscl, td0
  use ufs_flux_mod, only: flux_con_tol, flux_con_max_iter
  use ufs_flux_mod, only: loc_cpdair, loc_cpvir, loc_karman, loc_g, loc_zvir
  use ufs_flux_mod, only: loc_latvap, loc_stebol, use_coldair_outbreak_mod

  implicit none

  private
  public :: flux_atmocn_bulk

contains
  subroutine flux_atmocn_bulk(logunit, nMax, mask,                   &
       zbot, ubot, vbot, qbot, rbot, tbot, ts, us, vs, thbot, spval, &
       sen, lat, lwup, taux, tauy, evap, tref, qref, duu10n, rc)

    !--- input arguments --------------------------------
    integer,     intent(in) :: logunit
    integer(IN), intent(in) ::       nMax  ! data vector length
    integer(IN), intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
    real(R8),    intent(in) :: zbot (nMax) ! atm level height           (m)
    real(R8),    intent(in) :: ubot (nMax) ! atm u wind               (m/s)
    real(R8),    intent(in) :: vbot (nMax) ! atm v wind               (m/s)
    real(R8),    intent(in) :: qbot (nMax) ! atm specific humidity  (kg/kg)
    real(R8),    intent(in) :: rbot (nMax) ! atm air density       (kg/m^3)
    real(R8),    intent(in) :: tbot (nMax) ! atm T                      (K)
    real(R8),    intent(in) :: ts   (nMax) ! ocn temperature            (K)
    real(R8),    intent(in) :: us   (nMax) ! ocn u-velocity           (m/s)
    real(R8),    intent(in) :: vs   (nMax) ! ocn v-velocity           (m/s)
    real(R8),    intent(in) :: thbot(nMax) ! atm potential T            (K)
    real(R8),    intent(in) :: spval       ! masked value

    !--- output arguments -------------------------------
    real(R8),    intent(out) ::  sen  (nMax)    ! heat flux: sensible      (W/m^2)
    real(R8),    intent(out) ::  lat  (nMax)    ! heat flux: latent        (W/m^2)
    real(R8),    intent(out) ::  lwup (nMax)    ! heat flux: lw upward     (W/m^2)
    real(R8),    intent(out) ::  taux (nMax)    ! surface stress, zonal        (N)
    real(R8),    intent(out) ::  tauy (nMax)    ! surface stress, maridional   (N)
    real(R8),    intent(out) ::  evap (nMax)    ! water flux: evap    ((kg/s)/m^2)
    real(R8),    intent(out) ::  tref (nMax)    ! diag:  2m ref height T       (K)
    real(R8),    intent(out) ::  qref (nMax)    ! diag:  2m ref humidity   (kg/kg)
    real(R8),    intent(out) :: duu10n(nMax)    ! diag: 10m wind speed squared (m/s)^2
    integer,     intent(out) :: rc

    !--- local constants --------------------------------
    real(R8), parameter :: umin  =  0.5_R8 ! minimum wind speed       (m/s)
    real(R8), parameter :: zref  = 10.0_R8 ! reference height           (m)
    real(R8), parameter :: ztref =  2.0_R8 ! reference height for air T (m)

    !--- local variables --------------------------------
    integer     :: n      ! vector loop index
    integer     :: iter
    real(R8)    :: vmag   ! surface wind magnitude   (m/s)
    real(R8)    :: ssq    ! sea surface humidity     (kg/kg)
    real(R8)    :: delt   ! potential T difference   (K)
    real(R8)    :: delq   ! humidity difference      (kg/kg)
    real(R8)    :: stable ! stability factor
    real(R8)    :: rdn    ! sqrt of neutral exchange coeff (momentum)
    real(R8)    :: rhn    ! sqrt of neutral exchange coeff (heat)
    real(R8)    :: ren    ! sqrt of neutral exchange coeff (water)
    real(R8)    :: rd     ! sqrt of exchange coefficient (momentum)
    real(R8)    :: rh     ! sqrt of exchange coefficient (heat)
    real(R8)    :: re     ! sqrt of exchange coefficient (water)
    real(R8)    :: ustar  ! ustar
    real(r8)    :: ustar_prev
    real(R8)    :: qstar  ! qstar
    real(R8)    :: tstar  ! tstar
    real(R8)    :: hol    ! H (at zbot) over L
    real(R8)    :: xsq    ! ?
    real(R8)    :: xqq    ! ?
    real(R8)    :: psimh  ! stability function at zbot (momentum)
    real(R8)    :: psixh  ! stability function at zbot (heat and water)
    real(R8)    :: psix2  ! stability function at ztref reference height
    real(R8)    :: alz    ! ln(zbot/zref)
    real(R8)    :: al2    ! ln(zref/ztref)
    real(R8)    :: u10n   ! 10m neutral wind
    real(R8)    :: tau    ! stress at zbot
    real(R8)    :: cp     ! specific heat of moist air
    real(R8)    :: fac    ! vertical interpolation factor

    !--- local functions --------------------------------
    real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
    !!++ Large only (formula v*=[c4/U10+c5+c6*U10]*U10 in Large et al. 1994)
    real(R8)    :: cdn    ! function: neutral drag coeff at 10m
    !!++ Large only (stability functions)
    real(R8)    :: psimhu ! function: unstable part of psimh
    real(R8)    :: psixhu ! function: unstable part of psimx
    real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
    real(R8)    :: Tk     ! dummy arg ~ temperature (K)
    real(R8)    :: xd     ! dummy arg ~ ?
    !--- for cold air outbreak calc --------------------------------
    real(R8)    :: tdiff(nMax)               ! tbot - ts
    real(R8)    :: vscl

    qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)
    ! Large and Pond
    cdn(Umps)  =   0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps
    psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
    psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)

    u10n = spval
    rh = spval
    psixh = spval
    hol=spval

    rc = ESMF_SUCCESS

    !--- for cold air outbreak calc --------------------------------
    tdiff= tbot - ts

    al2 = log(zref/ztref)
    do n=1,nMax
       if (mask(n) /= 0) then

          !--- compute some needed quantities ---
          vmag   = max(umin, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )
          if (use_coldair_outbreak_mod) then
             ! Cold Air Outbreak Modification:
             ! Increase windspeed for negative tbot-ts
             ! based on Mahrt & Sun 1995,MWR

             if (tdiff(n).lt.td0) then
                vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag))),maxscl)
                vmag=vmag*vscl
             endif
          endif
          ssq    = 0.98_R8 * qsat(ts(n)) / rbot(n)   ! sea surf hum (kg/kg)
          delt   = thbot(n) - ts(n)                  ! pot temp diff (K)
          delq   = qbot(n) - ssq                     ! spec hum dif (kg/kg)
          alz    = log(zbot(n)/zref)
          cp     = loc_cpdair*(1.0_R8 + loc_cpvir*ssq)

          !------------------------------------------------------------
          ! first estimate of Z/L and ustar, tstar and qstar
          !------------------------------------------------------------
          !--- neutral coefficients, z/L = 0.0 ---
          stable = 0.5_R8 + sign(0.5_R8 , delt)
          rdn    = sqrt(cdn(vmag))
          rhn    = (1.0_R8-stable) * 0.0327_R8 + stable * 0.018_R8
          !(1.0_R8-stable) * chxcdu + stable * chxcds
          ren    = 0.0346_R8 !cexcd

          !--- ustar, tstar, qstar ---
          ustar = rdn * vmag
          tstar = rhn * delt
          qstar = ren * delq
          ustar_prev = ustar*2.0_R8
          iter = 0
          do while( abs((ustar - ustar_prev)/ustar) > flux_con_tol .and. iter < flux_con_max_iter)
             iter = iter + 1
             ustar_prev = ustar
             !--- compute stability & evaluate all stability functions ---
             hol  = loc_karman*loc_g*zbot(n)*  &
                  (tstar/thbot(n)+qstar/(1.0_R8/loc_zvir+qbot(n)))/ustar**2
             hol  = sign( min(abs(hol),10.0_R8), hol )
             stable = 0.5_R8 + sign(0.5_R8 , hol)
             xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
             xqq    = sqrt(xsq)
             psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
             psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)

             !--- shift wind speed using old coefficient ---
             rd   = rdn / (1.0_R8 + rdn/loc_karman*(alz-psimh))
             u10n = vmag * rd / rdn

             !--- update transfer coeffs at 10m and neutral stability ---
             rdn = sqrt(cdn(u10n))
             ren = 0.0346_R8 !cexcd
             rhn = (1.0_R8-stable)*0.0327_R8 + stable * 0.018_R8
             !(1.0_R8-stable) * chxcdu + stable * chxcds

             !--- shift all coeffs to measurement height and stability ---
             rd = rdn / (1.0_R8 + rdn/loc_karman*(alz-psimh))
             rh = rhn / (1.0_R8 + rhn/loc_karman*(alz-psixh))
             re = ren / (1.0_R8 + ren/loc_karman*(alz-psixh))

             !--- update ustar, tstar, qstar using updated, shifted coeffs --
             ustar = rd * vmag
             tstar = rh * delt
             qstar = re * delq
          enddo
          if (iter < 1) then
             write(logunit,*) 'iter<1 ',ustar,ustar_prev,flux_con_tol,flux_con_max_iter
             rc = ESMF_FAILURE
             return
          end if

          !------------------------------------------------------------
          ! compute the fluxes
          !------------------------------------------------------------

          tau = rbot(n) * ustar * ustar

          !--- momentum flux ---
          taux(n) = tau * (ubot(n)-us(n)) / vmag
          tauy(n) = tau * (vbot(n)-vs(n)) / vmag

          !--- heat flux ---
          sen (n) =          cp * tau * tstar / ustar
          lat (n) =  loc_latvap * tau * qstar / ustar
          lwup(n) = -loc_stebol * ts(n)**4

          !--- water flux ---
          evap(n) = lat(n)/loc_latvap

          !------------------------------------------------------------
          ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
          !------------------------------------------------------------
          hol = hol*ztref/zbot(n)
          xsq = max( 1.0_R8, sqrt(abs(1.0_R8-16.0_R8*hol)) )
          xqq = sqrt(xsq)
          psix2   = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
          fac     = (rh/loc_karman) * (alz + al2 - psixh + psix2 )
          tref(n) = thbot(n) - delt*fac
          tref(n) = tref(n) - 0.01_R8*ztref   ! pot temp to temp correction
          fac     = (re/loc_karman) * (alz + al2 - psixh + psix2 )
          qref(n) =  qbot(n) - delq*fac

          duu10n(n) = u10n*u10n ! 10m wind speed squared

       else
          !------------------------------------------------------------
          ! no valid data here -- out of domain
          !------------------------------------------------------------
          sen   (n) = spval  ! sensible         heat flux  (W/m^2)
          lat   (n) = spval  ! latent           heat flux  (W/m^2)
          lwup  (n) = spval  ! long-wave upward heat flux  (W/m^2)
          evap  (n) = spval  ! evaporative water flux ((kg/s)/m^2)
          taux  (n) = spval  ! x surface stress (N)
          tauy  (n) = spval  ! y surface stress (N)
          tref  (n) = spval  !  2m reference height temperature (K)
          qref  (n) = spval  !  2m reference height humidity (kg/kg)
          duu10n(n) = spval  ! 10m wind speed squared (m/s)^2

       endif
    enddo
  end subroutine flux_atmocn_bulk
end module flux_atmocn_bulk_mod
