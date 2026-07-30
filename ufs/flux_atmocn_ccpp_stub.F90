module flux_atmocn_ccpp_mod

  use ESMF,            only : ESMF_GridComp, ESMF_FAILURE
  use med_kind_mod,    only : R8=>SHR_KIND_R8

  implicit none

  private ! default private

  public :: flux_atmocn_ccpp ! computes atm/ocn fluxes
contains

  subroutine flux_atmocn_ccpp(gcomp, garea, maintask, logunit, nMax, mask, &
       pbot, tbot, qbot, zbot, ubot, vbot, rbot, ts, usfc, vsfc,           &
       psfc, lwdn, spval, sen, lat, lwup, evap, taux, tauy, tref, qref,    &
       duu10n, ustar_sv, re_sv, ssq_sv, rc)

    implicit none

    !--- input arguments --------------------------------
    type(ESMF_GridComp), intent(in) :: gcomp       ! gridded component
    real(r8), intent(in)  :: garea(nMax) ! grid area                      (m^2)
    logical , intent(in)  :: maintask    ! main task
    integer , intent(in)  :: logunit     ! log file unit number
    integer , intent(in)  :: nMax        ! data vector length
    integer , intent(in)  :: mask (nMax) ! ocn domain mask
    real(r8), intent(in)  :: pbot(nMax)  ! atm P (bottom)                 (Pa)
    real(r8), intent(in)  :: tbot(nMax)  ! atm T (bottom)                 (K)
    real(r8), intent(in)  :: qbot(nMax)  ! atm specific humidity (bottom) (kg/kg)
    real(r8), intent(in)  :: zbot(nMax)  ! atm level height               (m)
    real(r8), intent(in)  :: ubot(nMax)  ! atm u wind (bottom)            (m/s)
    real(r8), intent(in)  :: vbot(nMax)  ! atm v wind (bottom)            (m/s)
    real(r8), intent(in)  :: rbot(nMax)  ! atm density                    (kg/m^3)
    real(r8), intent(in)  :: ts(nMax)    ! ocn surface temperature        (K)
    real(r8), intent(in)  :: usfc(nMax)  ! atm u wind (surface)           (m/s)
    real(r8), intent(in)  :: vsfc(nMax)  ! atm v wind (surface)           (m/s)
    real(r8), intent(in)  :: psfc(nMax)  ! atm P (surface)                (Pa)
    real(r8), intent(in)  :: lwdn(nMax)  ! atm lw downward                (W/m^2)
    real(r8), intent(in)  :: spval       ! masked value

    !--- output arguments -------------------------------
    real(r8), intent(out) :: sen(nMax)      ! heat flux: sensible            (W/m^2)
    real(r8), intent(out) :: lat(nMax)      ! heat flux: latent              (W/m^2)
    real(r8), intent(out) :: lwup(nMax)     ! heat flux: lw upward           (W/m^2)
    real(r8), intent(out) :: evap(nMax)     ! heat flux: evap                ((kg/s)
    real(r8), intent(out) :: taux(nMax)     ! surface stress, zonal          (N)
    real(r8), intent(out) :: tauy(nMax)     ! surface stress, maridional     (N)
    real(r8), intent(out) :: tref (nMax)    ! diag: 2m ref height T          (K)
    real(r8), intent(out) :: qref(nMax)     ! diag: 2m ref humidity          (kg/kg)
    real(r8), intent(out) :: duu10n(nMax)   ! diag: 10m wind speed squared (m/s)^2
    real(r8), intent(out) :: ustar_sv(nMax) ! diag: ustar
    real(r8), intent(out) :: re_sv (nMax)   ! diag: sqrt of exchange coefficient (water)
    real(r8), intent(out) :: ssq_sv(nMax)   ! diag: sea surface humidity (kg/kg)
    integer,  intent(out) :: rc             ! return code

    ! provide explicit values for intent(out)
    sen = spval
    lat = spval
    lwup = spval
    evap = spval
    taux = spval
    tauy = spval
    tref = spval
    qref = spval
    duu10n = spval
    ustar_sv = spval
    re_sv = spval
    ssq_sv = spval

    if (maintask) then
       write(logunit,*) 'ERROR: ocn_surface_flux_scheme=1 and CMEPS_AOFLUX=OFF '
    end if
    rc = ESMF_FAILURE

  end subroutine flux_atmocn_ccpp
end module flux_atmocn_ccpp_mod
