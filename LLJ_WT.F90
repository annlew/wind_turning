subroutine calc_wind_turn (state, wturn, mflux, h, pblh, clubb_called)
    !
    !---------------------------------------------------------------
    !
    ! Purpose:  Calculate wind turning diagnostic
    !
    !---------------------------------------------------------------

    use cam_history,        only: addfld, add_default, horiz_only
    use ppgrid,             only: pcols, pver
    use hycoef,             only: hyam, hybm, ps0

    !
    ! Arguments
    !

    type(physics_state), intent(in) :: state  
    
    !
    !---------------------------Local workspace-----------------------------
    !

    ! Initialize local variables
    integer :: ncol
    integer :: pbl_idx(pcols)
    integer :: p_layer_850
    integer :: lchnk
    integer :: h_i
    integer :: t_i
    integer :: ixq

    logical, intent(in), optional :: clubb_called                            ! Whether CLUBB has been called
  
    real(r8) :: gas_constant_dry
    real(r8) :: gas_constant_dry_heat_capacity
    real(r8) :: g
    real(r8) :: pi
    real(r8) :: m(pcols, pver)
    real(r8) :: m0(pcols)
    real(r8) :: t0(pcols)
    real(r8) :: q0(pcols)
    real(r8) :: rho(pcols, pver)
    real(r8) :: rho_0(pcols)
    real(r8) :: rho_i(pcols, pver)
    real(r8) :: dp(pcols, pver)
    real(r8) :: dh(pcols, pver)
    real(r8) :: h(pcols, pver)
    real(r8) :: theta(pcols, pver)
    real(r8) :: theta_0(pcols)
    real(r8) :: theta_v(pcols, pver)
    real(r8) :: theta_v0(pcols)
    real(r8) :: Ri_b(pcols, pver)
    real(r8) :: k(pcols)
    real(r8) :: u_pbl(pcols)
    real(r8) :: v_pbl(pcols)
    !real(r8) :: Ri_b_850(pcols)
    real(r8) :: wdir_0(pcols)
    real(r8) :: wdir_pbl(pcols)
    real(r8) :: wdir(pcols, pver)
    real(r8) :: rotang
    real(r8) :: agwind(pcols, pver)
    real(r8) :: wturn(pcols, pver)
    real(r8) :: mflux(pcols, pver)
    real(r8) :: pblh(pcols)
    real(r8) :: t(pcols, pver)
    real(r8) :: q(pcols, pver)
    real(r8) :: u(pcols, pver)
    real(r8) :: v(pcols, pver)
    real(r8) :: p(pcols, pver)
    real(r8) :: ps(pcols)

    !
    !-----------------------------------------------------------------------
    !

    ! Load state variables
    call cnst_get_ind('Q', ixq)
    ncol = state%ncol
    lchnk = state%lchnk
    t = state%t(:ncol, :pver)
    q = state%q(:ncol, :pver, ixq)
    u = state%u(:ncol, :pver)
    v = state%v(:ncol, :pver)
    ps = state%ps(:ncol)

    q0 = q(:, pver)
    t0 = t(:, pver)

    ! Define constants
    pi = 3.14159265359
    gas_constant_dry = 2.87e2
    gas_constant_dry_heat_capacity = 2.86e-1
    g = 9.81

    ! Calculate pressure
    do t_i = 1, ncol
      do h_i = 1, pver
        p(t_i, h_i) = (hyam(h_i) * ps0) + (hybm(h_i) * ps(t_i))
      enddo
    enddo
  
    ! Calculate specific humidity from mixing ratio
    m = q / (1 - q)
    m0 = q0 / (1 - q0)
  
    ! Calculate rho
    rho = p / (gas_constant_dry * t)
    rho_0 = ps / (gas_constant_dry * t0)
  
    ! Calculate rho between model level interfaces
    rho_i(:, 1 : pver - 1) = ((rho(:, 1 : pver - 1) * p(:, 1 : pver - 1)) + (rho(:, 2 : pver) * p(:, 2 : pver))) / (p(:, 1 : pver - 1) + p(:, 2 : pver))
    rho_i(:, pver) = ((rho_0 * ps) + (rho(:, pver) * p(:, pver))) / (ps + p(:, pver))
  
    ! Calculate P differences between model levels
    dp(:, 1 : pver - 1) = p(:, 2 : pver) - p(:, 1 : pver - 1)

    do t_i = 1, ncol
      dp(t_i, pver) = ps(t_i) - p(t_i, pver)
    enddo
  
    ! Calculate height differences between model levels
    dh = dp / (g * rho_i)
  
    ! Calculate height
    do t_i = 1, ncol
     do h_i = 1, pver
       if (h_i == pver) then
         h(t_i, h_i) = dh(t_i, pver)
       else
         h(t_i, h_i) = sum(dh(t_i, h_i : pver))
       endif
     enddo
    enddo

    ! Calculate potential T
    theta = t * (1e5 / p) ** gas_constant_dry_heat_capacity
    theta_0 = t0 * (1e5 / ps) ** gas_constant_dry_heat_capacity
  
    ! Calculate virtual potential T
    theta_v = theta * (1 + (0.61 * m))
    theta_v0 = theta_0 * (1 + (0.61 * m0))
  
    ! Calculate bulk Richardson numbers
    do t_i = 1, ncol
      do h_i = 1, pver
        Ri_b(t_i, h_i) = (g / theta_v0(t_i)) * (theta_v(t_i, h_i) - theta_v0(t_i)) * h(t_i, h_i) / (u(t_i, h_i) ** 2 + v(t_i, h_i) ** 2)
      enddo
    enddo
  
    ! Calculate boundary layer index
    ! and interpolation coefficient
    oloop: do t_i = 1, ncol
      iloop: do h_i = 1, pver
        if (Ri_b(t_i, pver - h_i + 1) > 0.25) then
            pbl_idx(t_i) = pver - h_i + 1
            k(t_i) = h(t_i, pbl_idx(t_i)) - h(t_i, pbl_idx(t_i) + 1) / (Ri_b(t_i, pbl_idx(t_i)) - Ri_b(t_i, pbl_idx(t_i) + 1))
            exit iloop
        end if
      enddo iloop
    enddo oloop

    !do h_i = 1 : size(h)
    !  if (p(h_i) < 85000) then
    !      p_layer_850 = h_i
    !      exit
    !  end if
    !enddo
  
    do t_i = 1, ncol
      v_pbl(t_i) = v(t_i, pbl_idx(t_i))
      u_pbl(t_i) = u(t_i, pbl_idx(t_i))
    enddo
  
    ! Calculate wind directions
    wdir = 270 - ((180 / pi) * atan2(v, u))
    do t_i = 1, ncol
      do h_i = pver, 1
        if (wdir_0(t_i) > 360) then
          wdir(t_i, h_i) = wdir(t_i, h_i) - 360
        end if
      enddo
    enddo
  
    ! Calculate wind-turning over boundary layer and set it to fall between -180 and 180 degrees
    do t_i = 1, ncol
      do h_i = pbl_idx(t_i), pver
        wturn(t_i, h_i) = wdir(t_i, h_i) - wdir(t_i, pver)
        if (wturn(t_i, h_i) < -180) then
          wturn(t_i, h_i) = wturn(t_i, h_i) + 360
        else if (wturn(t_i, h_i) > 180) then
          wturn(t_i, h_i) = wturn(t_i, h_i) - 360
        endif
      enddo
      wturn(t_i, 1 : pbl_idx(t_i) - 2) = -9999
      pblh(t_i) = h(t_i, pbl_idx(t_i) + 1) + k(t_i) * (0.25 - Ri_b(t_i, pbl_idx(t_i) + 1))
    enddo

    ! Calculate ageostrophic wind and mass flux
    do t_i = 1, ncol
      ! Calculate rotation angle
      rotang = (270 - wdir(t_i, pbl_idx(t_i))) * (pi / 180)

      ! Calculate ageostrophic wind
      ! To do: calculate geostropic wind at each height using pressure gradient force
      do h_i = 1, pver
        agwind(t_i, h_i) = (u(t_i, h_i) - u(t_i, pbl_idx(t_i))) * sin(rotang) + (v(t_i, h_i) - v(t_i, pbl_idx(t_i))) * cos(rotang)
      enddo

      ! Use ageostrophic wind to calculate mass fluxes and set mass fluxes below PBLH to a large negative number
      do h_i = 1, pver
        mflux(t_i, h_i) = sum(rho(t_i, h_i : pver) * agwind(t_i, h_i : pver) * dh(t_i, h_i : pver))
      enddo
      
      mflux(t_i, 1 : pbl_idx(t_i) - 2) = -999999999 
    enddo

    call find_lljs(state, u, v, h, clubb_called)
  
    if (clubb_called == .FALSE.) then
      if (hist_fld_active('WTURN_BCL')) then
        call outfld('WTURN_BCL', wturn, pcols, lchnk)
        call outfld('MFLUX_BCL', mflux, pcols, lchnk)
        call outfld('PBLH_RI_BCL', pblh, pcols, lchnk)
        call outfld('VSTATE_BCL', state%v, pcols, lchnk)
        call outfld('LHEIGHT_BCL', h, pcols, lchnk)
        call outfld('USTATE_BCL', state%u, pcols, lchnk)
        call outfld('QSTATE_BCL', state%q, pcols, lchnk)
        call outfld('TSTATE_BCL', state%t, pcols, lchnk)
        call outfld('PSSTATE_BCL', state%ps, pcols, lchnk)
      end if
    else
      if (hist_fld_active('WTURN_ACL')) then
        call outfld('WTURN_ACL', wturn, pcols, lchnk)
        call outfld('MFLUX_ACL', mflux, pcols, lchnk)
        call outfld('PBLH_RI_ACL', pblh, pcols, lchnk)
        call outfld('VSTATE_ACL', state%v, pcols, lchnk)
        call outfld('LHEIGHT_ACL', h, pcols, lchnk)
        call outfld('USTATE_ACL', state%u, pcols, lchnk)
        call outfld('QSTATE_ACL', state%q, pcols, lchnk)
        call outfld('TSTATE_ACL', state%t, pcols, lchnk)
        call outfld('PSSTATE_ACL', state%ps, pcols, lchnk)
      end if
    end if

  end subroutine calc_wind_turn

  subroutine find_lljs (state, u, v, h, clubb_called)
    !
    !---------------------------------------------------------------
    !
    ! Purpose:  Find LLJ statistics
    !
    !---------------------------------------------------------------
    
    use cam_history,        only: addfld, add_default, horiz_only
    use ppgrid,             only: pcols, pver
    use hycoef,             only: hyam, hybm, ps0
    
    !
    ! Arguments
    !
    
    type(physics_state), intent(in) :: state  
        
    !
    !---------------------------Local workspace-----------------------------
    !
    
    ! Initialize local variables
    integer :: ncol
    integer :: pbl_idx(pcols)
    integer :: lchnk
    integer :: h_i
    integer :: t_i
    integer :: ixq
    integer :: max_idx_nr
    integer :: max_idx(pcols, 8)
    integer :: min_idx(pcols, 8)
    integer :: llj_idx(pcols)
    integer :: llj_pos_idx
    
    logical, intent(in), optional :: clubb_called                            ! Whether CLUBB has been called
    logical :: llj_found(pcols)

    real(r8) :: llj_u(pcols)
    real(r8) :: llj_v(pcols)
    real(r8) :: llj_h(pcols)
    real(r8) :: llj_wspeed(pcols)

    real(r8) :: u(pcols, pver)
    real(r8) :: v(pcols, pver)
    real(r8) :: h(pcols, pver)
    real(r8) :: wspeed(pcols, pver)

    real(r8) :: max_h(pcols, 8)
    real(r8) :: max_w(pcols, 8)
    real(r8) :: max_u(pcols, 8)
    real(r8) :: max_v(pcols, 8)

    real(r8) :: min_h(pcols, 8)
    real(r8) :: min_w(pcols, 8)

    u = state%u(:ncol, :pver)
    v = state%v(:ncol, :pver)

    wspeed = (u ** 2 + v ** 2) ** 0.5

    ! Give a default values for the min- and max arrays
    min_idx(:, :) = -9999
    max_idx(:, :) = -9999

    min_h(:, :) = -9999.
    min_w(:, :) = -9999.
    max_h(:, :) = -9999.
    max_w(:, :) = -9999.
    max_u(:, :) = -9999.
    max_v(:, :) = -9999.
    llj_found(:) = .FALSE.

    ! Check wind at each height to find maximums and minimums
    do t_i = 1, ncol
      do h_i = 1, pver
        if (h_i == 0) then 
          if (wspeed(t_i, pver - h_i + 1) > wspeed(t_i, h_i)) then
            max_h(t_i, h_i) = h(t_i, pver - h_i + 1)
            max_w(t_i, h_i) = wspeed(t_i, pver - h_i + 1)
            max_idx(t_i, h_i) = pver - h_i + 1
          endif
        else if (h_i == pver) then
          cycle
        else if (h(t_i, pver - h_i + 1) > 1500) then
          !if ((len(max_h) > 0) .AND. (len(min_h) == len(max_h))) then
          min_h(t_i, h_i) = h(t_i, pver - h_i + 2)
          min_w(t_i, h_i) = wspeed(t_i, pver - h_i + 2)
          min_idx(t_i, h_i) = pver - h_i + 1
          exit
        else if ((wspeed(t_i, pver - h_i + 1) > wspeed(t_i, pver - h_i + 2)) .AND. (wspeed(t_i, pver - h_i + 1) > wspeed(t_i, pver - h_i))) then
          max_h = h(t_i, pver - h_i + 1)
          max_w = wspeed(t_i, pver - h_i + 1)
          max_idx(t_i, h_i) = pver - h_i + 1
        else if ((h_i < pver) .AND. (wspeed(t_i, pver - h_i + 1) < wspeed(t_i, pver - h_i + 2)) .AND. (wspeed(t_i, pver - h_i + 1) < wspeed(t_i, pver - h_i))) then
          min_h(t_i, h_i) = h(t_i, pver - h_i + 1)
          min_w(t_i, h_i) = wspeed(t_i, pver - h_i + 1)
          min_idx(t_i, h_i) = pver - h_i + 1
        endif
      enddo
  
      ## Check if maximums are LLJs
      if (size(max_h) == 0) then
        cycle
      else
        do max_idx_nr = 1, 8
          if((max_w(t_i, max_idx_nr) >= min_w(t_i, max_idx_nr) + 2) .AND. (max_w(t_i, max_idx_nr) >= min_w(t_i, max_idx_nr + 1) + 2)) then
            if ((max_w(t_i, max_idx_nr) >= min_w(t_i, max_idx_nr) * 1.25) .AND. (max_w(t_i, max_idx_nr) >= min_w(t_i, max_idx_nr + 1) * 1.25)) then
              llj_pos_idx = max_idx_nr
              llj_found(t_i) = .TRUE.
              exit
            endif
          endif
        enddo
      endif
  
      ## If LLJ is found, add to LLJ array
      if (llj_found(t_i) == .TRUE.) then
        llj_h(t_i) = max_h(t_i, llj_pos_idx)
        llj_wspeed(t_i) = max_w(t_i, llj_pos_idx)
        llj_u(t_i) = max_u(t_i, llj_pos_idx)
        llj_v(t_i) = max_v(t_i, llj_pos_idx)
        llj_idx(t_i) = max_idx(t_i, llj_pos_idx)
      endif

    enddo

    if (clubb_called == .FALSE.) then
      if (hist_fld_active('WTURN_BCL')) then
        call outfld('LLJ_U_BCL', llj_u, pcols, lchnk)
        call outfld('LLJ_V_BCL', llj_v, pcols, lchnk)
        call outfld('LLJ_H_BCL', llj_h, pcols, lchnk)
        call outfld('LLJ_WSPEED_BCL', llj_wspeed, pcols, lchnk)
      end if
    else
      if (hist_fld_active('WTURN_ACL')) then
        call outfld('LLJ_U_ACL', llj_u, pcols, lchnk)
        call outfld('LLJ_V_ACL', llj_v, pcols, lchnk)
        call outfld('LLJ_H_ACL', llj_h, pcols, lchnk)
        call outfld('LLJ_WSPEED_ACL', llj_wspeed, pcols, lchnk)
      end if
    end if
  end subroutine find_lljs