module esmFldsExchange_coastal_mod

  use ESMF
  use NUOPC
  use med_utils_mod         , only : chkerr => med_utils_chkerr
  use med_kind_mod          , only : CX=>SHR_KIND_CX
  use med_kind_mod          , only : CS=>SHR_KIND_CS
  use med_kind_mod          , only : CL=>SHR_KIND_CL
  use med_kind_mod          , only : R8=>SHR_KIND_R8
  use med_internalstate_mod , only : compmed
  use med_internalstate_mod , only : compatm
  use med_internalstate_mod , only : compocn
  use med_internalstate_mod , only : compwav
  use med_internalstate_mod , only : compice
  use med_internalstate_mod , only : ncomps
  use med_internalstate_mod , only : coupling_mode

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  implicit none
  public

  public :: esmFldsExchange_coastal

  character(*), parameter :: u_FILE_u = &
       __FILE__

  type gcomp_attr
    character(len=CX) :: atm2ocn_fmap = 'unset'
    character(len=CX) :: atm2ocn_smap = 'unset'
    character(len=CX) :: atm2ocn_vmap = 'unset'
    character(len=CX) :: atm2wav_smap = 'unset'
    character(len=CX) :: ocn2atm_fmap = 'unset'
    character(len=CX) :: ocn2atm_smap = 'unset'
    character(len=CX) :: ocn2wav_smap = 'unset'
    character(len=CX) :: wav2ocn_smap = 'unset'
    character(len=CX) :: wav2atm_smap = 'unset'

    character(len=CX) :: atm2ice_smap = 'unset'
    character(len=CX) :: ocn2ice_smap = 'unset'
    character(len=CX) :: ice2ocn_smap = 'unset'
    character(len=CX) :: ice2atm_smap = 'unset'

    character(len=CS) :: mapnorm      = 'one'
    logical           :: atm_present  = .false.
    logical           :: ocn_present  = .false.
    logical           :: wav_present  = .false.
    logical           :: ice_present  = .false.
  end type

!===============================================================================
contains
!===============================================================================

  subroutine esmFldsExchange_coastal(gcomp, phase, rc)

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    character(len=*) , parameter   :: subname='(esmFldsExchange_coastal)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    if (phase == 'advertise') then
      call esmFldsExchange_coastal_advt(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    elseif (phase == 'fieldcheck') then
      call esmFldsExchange_coastal_fchk(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    elseif (phase == 'initialize') then
      call esmFldsExchange_coastal_init(gcomp, phase, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogSetError(ESMF_FAILURE, &
         msg=trim(subname)//": Phase is set to "//trim(phase), &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_coastal

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_coastal_advt(gcomp, phase, rc)

    use esmFlds, only : addfld_to => med_fldList_addfld_to
    use esmFlds, only : addfld_from => med_fldList_addfld_from

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    integer             :: n
    logical             :: isPresent
    character(len=CL)   :: cvalue
    character(len=CS)   :: fldname
    character(len=CS)   :: fldname1, fldname2
    type(gcomp_attr)    :: coastal_attr
    character(len=CS), allocatable :: S_flds(:)
    character(len=CS), allocatable :: F_flds(:,:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_coastal_advt)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    !=====================================================================
    ! scalar information
    !=====================================================================

    call NUOPC_CompAttributeGet(gcomp, name='ScalarFieldName', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", &
          value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,ncomps
          call addfld_from(n, trim(cvalue))
          call addfld_to(n, trim(cvalue))
       end do
    end if

    !=====================================================================
    ! attribute settings
    !=====================================================================
    call esmFldsExchange_coastal_attr(gcomp, coastal_attr, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !=====================================================================
    ! FIELDS TO MEDIATOR component (for fractions and atm/ocn flux calculation)
    !=====================================================================

    !----------------------------------------------------------
    ! to med: masks from components
    !----------------------------------------------------------
    call addfld_from(compocn, 'So_omask')
    call addfld_from(compice, 'Si_imask')

    !----------------------------------------------------------
    ! to med: frac from components
    !----------------------------------------------------------
    call addfld_to(compatm, 'So_ofrac')
          
    if (coastal_attr%ice_present .and. coastal_attr%atm_present) then
         call addfld_from(compice , 'Si_ifrac')
         call addfld_to(compatm   , 'Si_ifrac')

    end if

    if (coastal_attr%ice_present .and. coastal_attr%ocn_present) then
         call addfld_from(compice , 'Si_ifrac')
         call addfld_to(compocn   , 'Si_ifrac')
    end if

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ocn: state fields
    ! ---------------------------------------------------------------------
    if (coastal_attr%atm_present .and. coastal_attr%ocn_present) then
      allocate(S_flds(9))
      S_flds = (/'Sa_u10m   ', & ! inst_zonal_wind_height10m
                 'Sa_v10m   ', & ! inst_merid_wind_height10m
                 'Sa_pslv   ', & ! inst_pres_height_surface!
                 'Sa_t2m    ', & ! inst_temp_height2m
                 'Sa_q2m    ', & ! inst_spec_humid_height2m
                 'Faxa_lwdn ', & ! mean_down_lw_flx
                 'Faxa_swdn ', & ! mean_down_sw_flx
                 'Faxa_swnet', & ! mean_net_sw_flx
                 'Faxa_rain ' /) ! mean_prec_rate
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compatm, trim(fldname))
         call addfld_to(compocn, trim(fldname))
      end do
      deallocate(S_flds)
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: wave fields
    ! ---------------------------------------------------------------------
    if (coastal_attr%wav_present .and. coastal_attr%ocn_present) then
      allocate(S_flds(16))
      S_flds = (/'Sw_hs     ', & ! significant wave height
                 'Sw_bhd    ', & ! Bernoulli head (J term)
                 'Sw_tauox  ', & ! wave to ocean momentum flux x
                 'Sw_tauoy  ', & ! wave to ocean momentum flux y
                 'Sw_taubblx', & ! momentum flux due to bottom friction x
                 'Sw_taubbly', & ! momentum flux due to bottom friction y
                 'Sw_ubrx   ', & ! near bottom rms wave velocities x
                 'Sw_ubry   ', & ! near bottom rms wave velocities y
                 'Sw_thm    ', & ! mean wave direction
                 'Sw_t0m1   ', & ! mean wave period
                 'Sw_wnmean ', & ! mean wave number
                 'Sw_ustokes', & ! eastward_surface_stokes_drift_current
                 'Sw_vstokes', & ! northward_surface_stokes_drift_current
                 'Sw_wavsuu ', & ! eastward_wave_radiation_stress 
                 'Sw_wavsuv ', & ! eastward_northward_wave_radiation_stress
                 'Sw_wavsvv ' /) ! northward_wave_radiation_stress
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compwav, trim(fldname))
         call addfld_to(compocn, trim(fldname))
      end do
      deallocate(S_flds)
    end if

    !=====================================================================
    ! FIELDS TO WAVE (compwav)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to wav: state fields
    ! ---------------------------------------------------------------------
    if (coastal_attr%atm_present .and. coastal_attr%wav_present) then
      allocate(S_flds(2))
      S_flds = (/'Sa_u10m', & ! inst_zonal_wind_height10m
                 'Sa_v10m' /) ! inst_merid_wind_height10m
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compatm, trim(fldname))
         call addfld_to(compwav, trim(fldname))
      end do
      deallocate(S_flds)
    end if

    ! ---------------------------------------------------------------------
    ! to wav: ocean fields 
    ! ---------------------------------------------------------------------
    if (coastal_attr%ocn_present .and. coastal_attr%wav_present) then
      allocate(S_flds(3))
      S_flds = (/'So_h', & ! sea_surface_height_above_sea_level
                 'So_u', & ! ocn_current_zonal
                 'So_v' /) ! ocn_current_merid
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compocn, trim(fldname))
         call addfld_to(compwav, trim(fldname))
      end do
      deallocate(S_flds)
    end if

    ! ---------------------------------------------------------------------
    ! to ice: atm fields
    ! ---------------------------------------------------------------------
    if (coastal_attr%atm_present .and. coastal_attr%ice_present) then
      allocate(S_flds(14))
      S_flds = (/'Sa_u      ', & ! inst_zonalv_wind_height10m
                 'Sa_v      ', & ! inst_merid_wind_height10m                                     
                 'Sa_z      ', & ! lowest atmospheric height                                     
                 'Sa_tbot   ', & ! tenperature at lowest atm level                               
                 'Sa_shum   ', & ! Specific Humidity                                             
                 'Sa_pbot   ', & ! pressure bottom  
                 'Faxa_snow ', & ! Snow
                 'Faxa_rain ', & ! Rain
                 'Faxa_swvdr', & ! Short-wave
                 'Faxa_swvdf', & ! Short-wave
                 'Faxa_swndr', & ! Short-wave
                 'Faxa_swndf', & ! Short-wave
                 'Faxa_swnet', & ! Short-wave net
                 'Faxa_lwdn ' /) ! Long-wave 
                 
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compatm, trim(fldname))
         call addfld_to(compice, trim(fldname))
      end do
      deallocate(S_flds)
    end if

    ! ---------------------------------------------------------------------
    ! to ice: ocean fields 
    ! ---------------------------------------------------------------------
    if (coastal_attr%ocn_present .and. coastal_attr%ice_present) then
      allocate(S_flds(8))
      S_flds = (/'So_u   ',  & ! ocn_current_zonal
                 'So_v   ',  & ! ocn_current_merid
                 'So_t   ',  & ! SST
                 'So_s   ',  & ! SSS 
                 'So_h   ',  & ! depth of ml
                 'So_dhdx',  & ! tilt-x 
                 'So_dhdy',  & ! tilt-y 
                 'Fioo_q '/)    ! heat flux 
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compocn, trim(fldname))
         call addfld_to(compice, trim(fldname))
      end do
      deallocate(S_flds)
    end if 

    ! ---------------------------------------------------------------------
    ! to ocean: ice fields 
    ! ---------------------------------------------------------------------
    if (coastal_attr%ocn_present .and. coastal_attr%ice_present) then
      allocate(S_flds(13))
      S_flds = (/'Si_ifrac  ',    & ! Ice area fraction aice 
                 'Si_vice   ',    & ! Ice volume 
                 'Si_vsno   ',    & ! Snow volume 
                 'Si_uvel   ',    & ! Ice u-velocity
                 'Si_vvel   ',    & ! Ice v-velocity
                 'Fioi_meltw',    & ! Fresh water flux due to melting/freeze
                 'Fioi_melth',    & ! Ice-to-ocean heat flux at base of ice
                 'Fioi_taux ',    & ! Ice-to-ocean surface stress in x-dir
                 'Fioi_tauy ',    & ! Ice-to-ocean surface stress in y-dir
                 'Fioi_salt ',    & ! Salinity flux due to melting/freeze
                 'Fioi_swpen',    & ! Penatrive SW radiation flux 
                 'Si_frzmlt ',    & ! Energy used to create ice (all energy sst-tfrz)
                 'Si_CdnIO  ' /)    ! ice drag coeff

      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         call addfld_from(compice, trim(fldname))
         call addfld_to(compocn, trim(fldname))
      end do
      deallocate(S_flds)
    end if



    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_coastal_advt

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_coastal_fchk(gcomp, phase, rc)

    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    character(len=*) , parameter   :: subname='(esmFldsExchange_coastal_fchk)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (fldchk(is_local%wrap%FBImp(compocn,compocn),'So_omask',rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": Field connected "//"So_omask", &
          ESMF_LOGMSG_INFO)
    else
       call ESMF_LogSetError(ESMF_FAILURE, &
          msg=trim(subname)//": Field is not connected "//"So_omask", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_coastal_fchk

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_coastal_init(gcomp, phase, rc)

    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState
    use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch
    use med_internalstate_mod , only : mapfcopy, mapnstod, mapnstod_consd
    use med_internalstate_mod , only : mapfillv_bilnr
    use med_internalstate_mod , only : mapnstod_consf
    use med_internalstate_mod , only : mapbilnr_nstod
    use esmFlds               , only : addmap_from => med_fldList_addmap_from
    use esmFlds               , only : addmrg_to   => med_fldList_addmrg_to

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: n
    character(len=CS)   :: fldname
    character(len=CS)   :: fldname1, fldname2
    type(gcomp_attr)    :: coastal_attr
    character(len=CS), allocatable :: S_flds(:)
    character(len=CS), allocatable :: F_flds(:,:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_coastal_init)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------------
    ! Merging arguments:
    ! mrg_fromN = source component index that for the field to be merged
    ! mrg_fldN  = souce field name to be merged
    ! mrg_typeN = merge type ('copy', 'copy_with_weights', 'sum',
    !                         'sum_with_weights', 'merge')
    ! NOTE:
    ! mrg_from(compmed) can either be for mediator computed fields for atm/ocn
    ! fluxes or for ocn albedos
    !
    ! NOTE:
    ! FBMed_aoflux_o only refer to output fields to the atm/ocn that computed in
    ! the atm/ocn flux calculations. Input fields required from either the atm
    ! or the ocn for these computation will use the logical 'use_med_aoflux'
    ! below. This is used to determine mappings between the atm and ocn needed
    ! for these computations.
    !--------------------------------------

    !=====================================================================
    ! attribute settings
    !=====================================================================
    call esmFldsExchange_coastal_attr(gcomp, coastal_attr, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ocn: state fields
    ! ---------------------------------------------------------------------
    if (coastal_attr%atm_present .and. coastal_attr%ocn_present) then
      allocate(S_flds(9))
      S_flds = (/'Sa_u10m   ', & ! inst_zonal_wind_height10m
                 'Sa_v10m   ', & ! inst_merid_wind_height10m
                 'Sa_pslv   ', & ! inst_pres_height_surface!
                 'Sa_t2m    ', & ! inst_temp_height2m
                 'Sa_q2m    ', & ! inst_spec_humid_height2m
                 'Faxa_lwdn ', & ! mean_down_lw_flx
                 'Faxa_swdn ', & ! mean_down_sw_flx
                 'Faxa_swnet', & ! mean_net_sw_flx
                 'Faxa_rain ' /) ! mean_prec_rate
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname),rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
            ) then
            call addmap_from(compatm, trim(fldname), compocn, &
                 mapbilnr_nstod, coastal_attr%mapnorm, coastal_attr%atm2ocn_smap)
            call addmrg_to(compocn, trim(fldname), &
                 mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
         end if
      end do
      deallocate(S_flds)
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: wave fields
    ! ---------------------------------------------------------------------
    if (coastal_attr%wav_present .and. coastal_attr%ocn_present) then
      allocate(S_flds(16))
      S_flds = (/'Sw_hs     ', & ! significant wave height
                 'Sw_bhd    ', & ! Bernoulli head (J term)
                 'Sw_tauox  ', & ! wave to ocean momentum flux x
                 'Sw_tauoy  ', & ! wave to ocean momentum flux y
                 'Sw_taubblx', & ! momentum flux due to bottom friction x
                 'Sw_taubbly', & ! momentum flux due to bottom friction y
                 'Sw_ubrx   ', & ! near bottom rms wave velocities x
                 'Sw_ubry   ', & ! near bottom rms wave velocities y
                 'Sw_thm    ', & ! mean wave direction
                 'Sw_t0m1   ', & ! mean wave period
                 'Sw_wnmean ', & ! mean wave number
                 'Sw_ustokes', & ! eastward_surface_stokes_drift_current
                 'Sw_vstokes', & ! northward_surface_stokes_drift_current
                 'Sw_wavsuu ', & ! eastward_wave_radiation_stress 
                 'Sw_wavsuv ', & ! eastward_northward_wave_radiation_stress
                 'Sw_wavsvv ' /) ! northward_wave_radiation_stress
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname),rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compwav,compwav),trim(fldname),rc=rc) &
            ) then
            call addmap_from(compwav, trim(fldname), compocn, &
                 mapbilnr_nstod, coastal_attr%mapnorm, coastal_attr%wav2ocn_smap)
            call addmrg_to(compocn, trim(fldname), &
                 mrg_from=compwav, mrg_fld=trim(fldname), mrg_type='copy')
         end if
      end do
      deallocate(S_flds)
    end if

    !=====================================================================
    ! FIELDS TO WAVE (compwav)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to wav: state fields 
    ! ---------------------------------------------------------------------
    if (coastal_attr%atm_present .and. coastal_attr%wav_present) then
      allocate(S_flds(2))
      S_flds = (/'Sa_u10m', & ! inst_zonal_wind_height10m
                 'Sa_v10m' /) ! inst_merid_wind_height10m
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compwav),trim(fldname),rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
            ) then
            call addmap_from(compatm, trim(fldname), compwav, &
                 mapbilnr_nstod, coastal_attr%mapnorm, coastal_attr%atm2wav_smap)
            call addmrg_to(compwav, trim(fldname), &
                 mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
         end if
      end do
      deallocate(S_flds)
    end if

    ! ---------------------------------------------------------------------
    ! to wav: ocean fields 
    ! ---------------------------------------------------------------------
    if (coastal_attr%atm_present .and. coastal_attr%wav_present) then
      allocate(S_flds(3))
      S_flds = (/'So_h', & ! sea_surface_height_above_sea_level
                 'So_u', & ! ocn_current_zonal
                 'So_v' /) ! ocn_current_merid
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compwav),trim(fldname),rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compocn,compocn),trim(fldname),rc=rc) &
            ) then
            call addmap_from(compocn, trim(fldname), compwav, &
                 mapbilnr_nstod, coastal_attr%mapnorm, coastal_attr%ocn2wav_smap)
            call addmrg_to(compwav, trim(fldname), &
                 mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
         end if
      end do
      deallocate(S_flds)
    end if
    
    
    !======================================================================
    ! BEGIN FIELDS TO ICE (compice)
    !======================================================================
    ! ---------------------------------------------------------------------                      
    ! To ice from atm
    ! ---------------------------------------------------------------------                      
    ! ----------------------------------------  
    ! | Atmospheric Fields:                  |
    ! ----------------------------------------
    !   - zonal wind height 10m     [m/s]
    !   - merid wind height 10m     [m/s]
    !   - lowest atmospheric height [m] (defualt in CICE is 10m see model comp)
    !   - Temperature at 10m        [k]
    !   - Specific Humidity         [kg/kg]
    !   - Pressure at 10m           [Pa]
    !   - Rainfall rate             [kg/m/m/s]
    !   - Snowfall rate             [kg/m/m/s]
    !   - Net shortwave rad         [W/m/m] Other bands are derived from this
    !   - Net longwave  rad         [W/m/m] 
    ! ---------------------------------------
    if (coastal_attr%atm_present .and. coastal_attr%ice_present) then                            
      allocate(S_flds(14))                                                                       
      S_flds = (/'Sa_u      ', & ! inst_zonal_wind_height10m
                 'Sa_v      ', & ! inst_merid_wind_height10m                                     
                 'Sa_z      ', & ! lowest atmospheric height                                     
                 'Sa_tbot   ', & ! tenperature at lowest atm level                               
                 'Sa_shum   ', & ! Specific Humidity                                             
                 'Sa_pbot   ', & ! pressure bottom                                               
                 'Faxa_rain ', & ! Rain
                 'Faxa_snow ', & ! Snow
                 'Faxa_swvdr', & ! Short-wave
                 'Faxa_swvdf', & ! Short-wave                                                    
                 'Faxa_swndr', & ! Short-wave                                                    
                 'Faxa_swndf', & ! Short-wave
                 'Faxa_swnet', & ! Short-wave net
                 'Faxa_lwdn '/)  ! long-wave                                                     
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compice),trim(fldname),rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm,compatm),trim(fldname),rc=rc) &
            ) then
            call addmap_from(compatm, trim(fldname), compice, &
                 mapnstod_consf, coastal_attr%mapnorm, coastal_attr%atm2ice_smap) !!mapnstod_consf conservative interp 
            call addmrg_to(compice, trim(fldname), &
                 mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
         end if
      end do
      deallocate(S_flds)
    end if

    ! ---------------------------------------------------------------------                      
    ! To ice from ocn fields
    ! --------------------------------------------------------------------- 
    ! ----------------------------------------
    ! | Oceanic fields:                      |
    ! ----------------------------------------
    !   - zonal current             [m/s]
    !   - merid current             [m/s]
    !   - Sea Surface temperature   [K]
    !   - Sea surface salinity      [psu. (g/kg)]
    !   - Mixed layer depth         [m]
    !   - Zonal sea surface tilt    [1]
    !   - Merid sea surface tilt    [1]
    !   - deep ocean heatflux       [W/m/m]
    ! ---------------------------------------

    if (coastal_attr%ocn_present .and. coastal_attr%ice_present) then
      allocate(S_flds(8))
      S_flds = (/'So_u   ',  & ! ocn_current_zonal
                 'So_v   ',  & ! ocn_current_merid
                 'So_t   ',  & ! SST
                 'So_s   ',  & ! SSS 
                 'So_h   ',  & ! depth of ml
                 'So_dhdx',  & ! tilt-x 
                 'So_dhdy',  & ! tilt-y 
                 'Fioo_qi'/)   ! heat flux  

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compice),trim(fldname),rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compocn,compocn),trim(fldname),rc=rc) &
            ) then
            call addmap_from(compocn, trim(fldname), compice, &
                 mapbilnr_nstod, coastal_attr%mapnorm, coastal_attr%ocn2ice_smap)
            call addmrg_to(compice, trim(fldname), &
                 mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
         end if
      end do
      deallocate(S_flds)
    end if

    !======================================================================
    ! END FIELDS TO ICE (compice)
    !======================================================================
    
    !======================================================================
    ! BEGIN FIELDS FROM ICE (compice)
    !======================================================================
    ! ---------------------------------------------------------------------
    ! To ocn from ice fields
    ! ---------------------------------------------------------------------
    ! ----------------------------------------
    ! | Ice fields:                          |
    ! ----------------------------------------
    ! - Ice area fraction        [1]
    ! - Ice volume               [m^3]
    ! - Snow volume              [m^3]
    ! - Ice merid vel.           [m/s]
    ! - Ice Zonal vel.           [m/s]
    ! - Fresh water flux         [kg/m/m/s]
    ! - Heat flux at base of ice [W/m/m]
    ! - Ice-to-Ocn merid stress  [N/m/m]
    ! - Ice-to-Ocn zonal stress  [N/m/m]
    ! - Salinity flux            [kg/m/m/s]
    ! - Shortwave pen. rad       [W/m/m]
    ! - Freeze melt potential    [W/m/m] (Internal energy mixed layer)
    ! - Ice ocean drag           [1] (Computed by ice component)

    if (coastal_attr%ocn_present .and. coastal_attr%ice_present) then
      allocate(S_flds(13))
      S_flds = (/'Si_ifrac  ',    & ! Ice area fraction aice
                 'Si_vice   ',    & ! Ice volume
                 'Si_vsno   ',    & ! Snow volume
                 'Si_uvel   ',    & ! Ice u-velocity
                 'Si_vvel   ',    & ! Ice v-velocity
                 'Fioi_meltw',    & ! Fresh water flux due to melting/freeze
                 'Fioi_melth',    & ! Ice-to-ocean heat flux at base of ice
                 'Fioi_taux ',    & ! Ice-to-ocean surface stress in x-dir
                 'Fioi_tauy ',    & ! Ice-to-ocean surface stress in y-dir
                 'Fioi_salt ',    & ! Salinity flux due to melting/freeze
                 'Fioi_swpen',    & ! Penatrive SW radiation flux 
                 'Si_frzmlt ',    & ! Energy used to create ice (all energy sst-tfrz)
                 'Si_CdnIO  ' /)      ! ice drag
      
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
      do n = 1,size(S_flds)
         fldname = trim(S_flds(n))
         if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname),rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compice,compice),trim(fldname),rc=rc) &
            ) then
            call addmap_from(compice, trim(fldname), compocn, &
                 mapnstod_consf, coastal_attr%mapnorm, coastal_attr%ice2ocn_smap)
            call addmrg_to(compocn, trim(fldname), &
                 mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy')
         end if
      end do
      deallocate(S_flds)
    end if

    !======================================================================
    ! END FIELDS FROM ICE (compice)
    !======================================================================

  end subroutine esmFldsExchange_coastal_init

  !-----------------------------------------------------------------------------

  subroutine esmFldsExchange_coastal_attr(gcomp, coastal_attr, rc)

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    type(gcomp_attr) , intent(inout) :: coastal_attr
    integer          , intent(inout) :: rc

    ! local variables:
    character(32)       :: cname
    integer             :: verbosity, diagnostic
    character(len=CL)   :: cvalue
    logical             :: isPresent, isSet
    character(len=*) , parameter   :: subname='(esmFldsExchange_coastal_attr)'
    !--------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call NUOPC_CompGet(gcomp, name=cname, verbosity=verbosity, &
      diagnostic=diagnostic, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------
    ! Component active or not?
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='ATM_model', &
       value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'satm') coastal_attr%atm_present = .true.
    end if

    call NUOPC_CompAttributeGet(gcomp, name='OCN_model', &
       value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'socn') coastal_attr%ocn_present = .true.
    end if

    call NUOPC_CompAttributeGet(gcomp, name='WAV_model', &
       value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'swav') coastal_attr%wav_present = .true.
    end if

    call NUOPC_CompAttributeGet(gcomp, name='ICE_model', &
       value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'sice') coastal_attr%ice_present = .true.
    end if

    !----------------------------------------------------------
    ! Normalization type
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='normalization', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='normalization', &
          value=coastal_attr%mapnorm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------------------------------------------
    ! Initialize mapping file names
    !----------------------------------------------------------

    ! to atm
    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', &
          value=coastal_attr%ocn2atm_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', &
          value=coastal_attr%ocn2atm_fmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! to ocn
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', &
          value=coastal_attr%atm2ocn_fmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', &
       value=coastal_attr%atm2ocn_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', &
          value=coastal_attr%atm2ocn_vmap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! to wav
    call NUOPC_CompAttributeGet(gcomp, name='atm2wav_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2wav_smapname', &
          value=coastal_attr%atm2wav_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ocn2wav_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2wav_smapname', &
          value=coastal_attr%ocn2wav_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! to ice
    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='atm2ice_smapname', &
          value=coastal_attr%atm2ice_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ocn2ice_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ocn2ice_smapname', &
          value=coastal_attr%ocn2ice_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    call NUOPC_CompAttributeGet(gcomp, name='ice2ocn_smapname', &
       isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='ice2ocn_smapname', &
          value=coastal_attr%ice2ocn_smap, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Log Attribute Settings
    if (btest(verbosity,16)) then
       write(cvalue,"(I0)") verbosity
       call ESMF_LogWrite(trim(subname)//': Verbosity        = '// &
          trim(cvalue), ESMF_LOGMSG_INFO)
       write(cvalue,"(I0)") diagnostic
       call ESMF_LogWrite(trim(subname)//': Diagnostic       = '// &
          trim(cvalue), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': normalization    = '// &
          trim(coastal_attr%mapnorm), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': ocn2atm_smapname = '// &
          trim(coastal_attr%ocn2atm_smap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': ocn2atm_fmapname = '// &
          trim(coastal_attr%ocn2atm_fmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_fmapname = '// &
          trim(coastal_attr%atm2ocn_fmap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_smapname = '// &
          trim(coastal_attr%atm2ocn_smap), ESMF_LOGMSG_INFO)
       call ESMF_LogWrite(trim(subname)//': atm2ocn_vmapname = '// &
          trim(coastal_attr%atm2ocn_vmap), ESMF_LOGMSG_INFO)
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine esmFldsExchange_coastal_attr

end module esmFldsExchange_coastal_mod
