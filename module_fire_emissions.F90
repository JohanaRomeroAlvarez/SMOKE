!> Debugged and corrected module_fire_emissions

module module_fire_emissions
    use mpas_kind_types
    use mpas_smoke_init, only : p_smoke_fine, p_smoke_coarse
    use mpas_smoke_config
  
    implicit none
    private
 
    public :: compute_emission_factors, calculate_smoke_emissions, diurnal_cycle
  

    REAL(RKIND), PARAMETER :: secsperhour = 3600.0_RKIND
    LOGICAL, PARAMETER :: use_classic_pm25_ef = .false.

    !>-- Emissions factors parameters:frst,hwd,mxd,shrb,shrb_grs,grs,org_soil_mix,org_soil_woody
    real(RKIND), dimension(8), parameter :: EF_FLM = &
      (/19.0_RKIND, 9.4_RKIND, 14.6_RKIND, 9.3_RKIND, &
        10.7_RKIND, 13.3_RKIND, 10.8_RKIND, 11.8_RKIND/)

    real(RKIND), dimension(8), parameter :: EF_SML = &
      (/28.0_RKIND, 37.7_RKIND, 17.6_RKIND, 36.6_RKIND, &
        36.7_RKIND, 38.4_RKIND, 20.7_RKIND, 23.3_RKIND/)                               
    
    !>-- Emissions factors parameters GEUFE:IGBP 10, 12 14
    real(RKIND), dimension(3), parameter :: EF_EAST = &
      (/9.29_RKIND, 11.91_RKIND, 8.20_RKIND/)
    !>-- Biomass Combustion Coefficient (Beta Factor, β)
    real(RKIND), parameter :: BETA_WOOSTER = 0.38_RKIND
    !>-- Biomass Combustion Coefficient: SA,SAOS,AG,AGOS,EF,EFOS,PEAT
    real(RKIND), dimension(7), parameter :: BETA_KAISER = &
      (/0.78_RKIND, 0.26_RKIND, 0.29_RKIND, 0.13_RKIND, &
        0.49_RKIND, 1.55_RKIND, 5.87_RKIND/)

    ! Andreae (2019) biomass-burning categories.
    ! I remapped it as much as I can
    integer, parameter :: cat_none      = 0
    integer, parameter :: cat_savanna   = 1
    integer, parameter :: cat_tropical  = 2
    integer, parameter :: cat_temperate = 3
    integer, parameter :: cat_boreal    = 4
    integer, parameter :: cat_peat      = 5
    integer, parameter :: cat_agric     = 6
    integer, parameter :: num_fire_cat  = 6

    ! Classic category-average PM2.5 emission factors in g kg-1 dry matter.
    ! Order is savanna/grassland, tropical forest, temperate forest,
    ! boreal forest, peat, and agricultural residue.
    ! This is recommendation while we are doing retro test for the
    ! New emission based on smoldering and flaming.
    real(RKIND), dimension(num_fire_cat), parameter :: EF_PM25_CATEGORY = &
      (/ 6.7_RKIND, 8.3_RKIND, 18.5_RKIND, 18.7_RKIND, &
         18.9_RKIND, 8.2_RKIND /)
  contains
  
  !-------------------------------------------------------------------------------
  ! Subroutine: compute_emission_factors
  ! Description:
  !   Computes a 2D emission factor map (EFs_map) based on vegetation fractions
  !   and ecoregions. Emission factors are calculated using a combination of
  !   flaming/smoldering coefficients and beta scaling. If bb_beta = 1, beta_kaiser
  !   is used only for Canadian forest/shrubland types; otherwise, beta_wooster is applied.
  !   EF_east factors are applied to Grasslands and Croplands based on region masks.
  !-------------------------------------------------------------------------------
    subroutine compute_emission_factors(EFs_map, landusef, bb_beta,                      &
                                        eco_id, efs_smold, efs_flam, efs_rsmold,         &
                                        fmc_avg, wetness,                                &
                                        ef_diurnal_opt, nlcat, dt, gmt, xlong,           & !JR added xlong,ef_diurnal_opt
                                        ids, ide, jds, jde, kds, kde,                    &
                                        ims, ime, jms, jme, kms, kme,                    &
                                        its, ite, jts, jte, kts, kte)
  
      implicit none

      integer, intent(in) :: bb_beta, nlcat
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte

      logical, intent(in) :: ef_diurnal_opt

      integer, intent(in), dimension(ims:ime,jms:jme) :: eco_id

      real(RKIND), intent(in) :: dt, gmt
      real(RKIND), intent(in), dimension(ims:ime,jms:jme) :: &
              efs_smold, efs_flam, efs_rsmold, fmc_avg, wetness, xlong
      real(RKIND), intent(in), dimension(ims:ime,nlcat,jms:jme) :: landusef
      real(RKIND), intent(out), dimension(ims:ime,jms:jme) :: EFs_map

      integer :: i, j, n
      real(RKIND) :: class_factor

      EFs_map = 0.0_RKIND
      
      do j = jts, jte
         do i = its, ite
            do n = 1, min(20,nlcat)

               if (use_classic_pm25_ef) then
                  call get_classic_pm25_factor( &
                       n, eco_id(i,j), bb_beta, class_factor)    
               else
                  call get_land_class_factor( &
                       n, eco_id(i,j), bb_beta, &
                       ef_diurnal_opt, gmt, xlong(i,j), &
                       class_factor) 
               endif     
               EFs_map(i,j) = EFs_map(i,j) + max(0.0_RKIND,landusef(i,n,j)) * class_factor

            enddo
         enddo
      enddo  
  
    end subroutine compute_emission_factors
 

      !--------------------------------------------------------------------------
   ! Return the emission-factor and beta product for one IGBP land class.
   !--------------------------------------------------------------------------
   subroutine get_land_class_factor(land_class, eco, bb_beta,                  &
                                    ef_diurnal_opt, gmt, xlong, class_factor)

      implicit none

      integer, intent(in) :: land_class, eco, bb_beta
      logical, intent(in) :: ef_diurnal_opt
      real(RKIND), intent(in) :: gmt, xlong
      real(RKIND), intent(out) :: class_factor

      real(RKIND) :: base_flam_frac
      real(RKIND) :: ef, beta_val

      ef             = 0.0_RKIND
      beta_val       = 0.0_RKIND
      base_flam_frac = 0.0_RKIND

      select case (land_class)

      case (1)
         base_flam_frac = 0.75_RKIND
         call mixed_combustion_ef(base_flam_frac, EF_FLM(1), EF_SML(1),        &
                                  ef_diurnal_opt, gmt, xlong, ef)
         if (bb_beta == 1 .and. (eco == 2 .or. eco == 3 .or. eco == 4)) then
            beta_val = BETA_KAISER(6)
         else
            beta_val = BETA_WOOSTER
         end if
      case (2)
         base_flam_frac = 0.75_RKIND
         call mixed_combustion_ef(base_flam_frac, EF_FLM(1), EF_SML(1),        &
                                  ef_diurnal_opt, gmt, xlong, ef)
         if (bb_beta == 1 .and. eco == 2) then
            beta_val = BETA_KAISER(6)
         else if (bb_beta == 1 .and. eco == 5) then
            beta_val = BETA_KAISER(5)
         else
            beta_val = BETA_WOOSTER
         end if

      case (3,4)
         base_flam_frac = 0.80_RKIND
         call mixed_combustion_ef(base_flam_frac, EF_FLM(2), EF_SML(2),        &
                                  ef_diurnal_opt, gmt, xlong, ef)
         if (bb_beta == 1 .and. eco == 2) then
            beta_val = BETA_KAISER(6)
         else if (bb_beta == 1 .and. eco == 5) then
            beta_val = BETA_KAISER(5)
         else
            beta_val = BETA_WOOSTER
         end if

      case (5)
         base_flam_frac = 0.85_RKIND
         call mixed_combustion_ef(base_flam_frac, EF_FLM(3), EF_SML(3),        &
                                  ef_diurnal_opt, gmt, xlong, ef)
         if (bb_beta == 1 .and. eco == 2) then
            beta_val = BETA_KAISER(6)
         else if (bb_beta == 1 .and. eco == 5) then
            beta_val = BETA_KAISER(5)
         else
            beta_val = BETA_WOOSTER
         end if

      case (6,7)
         base_flam_frac = 0.95_RKIND
         call mixed_combustion_ef(base_flam_frac, EF_FLM(4), EF_SML(4),        &
                                  ef_diurnal_opt, gmt, xlong, ef)
         if (bb_beta == 1 .and. (eco == 2 .or. eco == 3 .or. eco == 4)) then
            beta_val = BETA_KAISER(2)
         else
            beta_val = BETA_WOOSTER
         end if

      case (8)
         base_flam_frac = 0.95_RKIND
         call mixed_combustion_ef(base_flam_frac, EF_FLM(5), EF_SML(5),        &
                                  ef_diurnal_opt, gmt, xlong, ef)
         beta_val = BETA_WOOSTER

      case (9)
         base_flam_frac = 0.95_RKIND
         call mixed_combustion_ef(base_flam_frac, EF_FLM(6), EF_SML(6),        &
                                  ef_diurnal_opt, gmt, xlong, ef)
         beta_val = BETA_WOOSTER

      case (10)
         if (eco == 10) then
            ef = EF_EAST(1)
         else
            ef = EF_FLM(6)
         end if
         beta_val = BETA_WOOSTER

      case (11)
         ef       = 18.9_RKIND
         beta_val = BETA_WOOSTER

      case (12)
         if (eco == 8 .or. eco == 9) then
            ef = EF_EAST(2)
         else
            ef = 8.2_RKIND
         end if
         beta_val = BETA_WOOSTER

      case (14)
         if (eco == 8 .or. eco == 9) then
            ef = EF_EAST(3)
         else
            ef = 8.2_RKIND
         end if
         beta_val = BETA_WOOSTER

      case (19)
         base_flam_frac = 0.70_RKIND
         call mixed_combustion_ef(base_flam_frac, EF_FLM(7), EF_SML(7),        &
                                  ef_diurnal_opt, gmt, xlong, ef)
         beta_val = BETA_WOOSTER

      case default
         ef       = 0.0_RKIND
         beta_val = 0.0_RKIND

      end select

      class_factor = ef * beta_val

   end subroutine get_land_class_factor

   !--------------------------------------------------------------------------
   ! Return the classic category-average PM2.5 EF and beta product for one
   ! IGBP land class. The land-class/ecoregion mapping is shared with the
   ! category-average gas-phase calculation.
   !--------------------------------------------------------------------------
   subroutine get_classic_pm25_factor(land_class, eco, bb_beta, class_factor)

      implicit none

      integer, intent(in) :: land_class, eco, bb_beta
      real(RKIND), intent(out) :: class_factor

      integer :: fire_category
      real(RKIND) :: beta_val

      class_factor = 0.0_RKIND

      fire_category = get_fire_category(land_class, eco)
      if (fire_category == cat_none) return

      call get_gas_beta_factor(land_class, eco, bb_beta, beta_val)
      if (beta_val <= 0.0_RKIND) return

      class_factor = EF_PM25_CATEGORY(fire_category) * beta_val

   end subroutine get_classic_pm25_factor

   !--------------------------------------------------------------------------
   ! Map fractional IGBP land cover and EPA Level 1 ecoregion to the
   ! Andreae (2019) combustion categories.
   !--------------------------------------------------------------------------
   integer function get_fire_category(land_class, eco) result(category)

      implicit none

      integer, intent(in) :: land_class, eco

      select case (land_class)

      case (1:5)
         select case (eco)
         case (1:5)
            category = cat_boreal
         case (14,15)
            category = cat_tropical
         case (6:13)
            category = cat_temperate
         case default
            if (land_class == 2) then
               category = cat_tropical
            else if (land_class == 3) then
               category = cat_boreal
            else
               category = cat_temperate
            end if
         end select
  
      case (6:10)
         category = cat_savanna

      case (11)
         ! Permanent wetland is represented by the peat-fire category.
         category = cat_peat

      case (12,14)
         category = cat_agric

      case (18,19)
         category = cat_boreal

      case default
         ! Urban, snow and ice, barren land, water, and barren tundra.
         category = cat_none

      end select

   end function get_fire_category

   !--------------------------------------------------------------------------
   ! Reproduce the beta selection used by the current PM2.5 calculation.
   ! This duplication intentionally avoids changing the working PM2.5 code.
   !--------------------------------------------------------------------------
   subroutine get_gas_beta_factor(land_class, eco, bb_beta, beta_val)

      implicit none

      integer, intent(in) :: land_class, eco, bb_beta
      real(RKIND), intent(out) :: beta_val

      beta_val = 0.0_RKIND

      select case (land_class)

      case (1)
         if (bb_beta == 1 .and. &
             (eco == 2 .or. eco == 3 .or. eco == 4)) then
            beta_val = BETA_KAISER(6)
         else
            beta_val = BETA_WOOSTER
         end if

      case (2:5)
         if (bb_beta == 1 .and. eco == 2) then
            beta_val = BETA_KAISER(6)
         else if (bb_beta == 1 .and. eco == 5) then
            beta_val = BETA_KAISER(5)
         else
            beta_val = BETA_WOOSTER
         end if

      case (6,7)
         if (bb_beta == 1 .and. &
             (eco == 2 .or. eco == 3 .or. eco == 4)) then
            beta_val = BETA_KAISER(2)
         else
            beta_val = BETA_WOOSTER
         end if

      case (8:12,14,18,19)
         beta_val = BETA_WOOSTER

      case default
         beta_val = 0.0_RKIND

      end select

   end subroutine get_gas_beta_factor

   !--------------------------------------------------------------------------
   ! Blend flaming and smoldering PM2.5 emission factors.
   !--------------------------------------------------------------------------
   subroutine mixed_combustion_ef(base_flam_frac, ef_flam, ef_smold,           &
                                  ef_diurnal_opt, gmt, xlong, ef)

      implicit none

      real(RKIND), intent(in) :: base_flam_frac, ef_flam, ef_smold
      real(RKIND), intent(in) :: gmt, xlong
      logical, intent(in) :: ef_diurnal_opt
      real(RKIND), intent(out) :: ef

      real(RKIND) :: flam_frac, smold_frac

      if (ef_diurnal_opt) then
         flam_frac = get_flam_frac(base_flam_frac, gmt, xlong)
      else
         flam_frac = base_flam_frac
      endif

      smold_frac = 1.0_RKIND - flam_frac
      ef = flam_frac * ef_flam + smold_frac * ef_smold

   end subroutine mixed_combustion_ef

   !--------------------------------------------------------------------------
   ! Select the EBB2 time block from elapsed model time.
   !--------------------------------------------------------------------------
   integer function get_fre_block(ktau, dt, bb_input_prevh, nblocks) result(blk)

      implicit none

      integer, intent(in) :: ktau, bb_input_prevh, nblocks
      real(RKIND), intent(in) :: dt

      integer :: hour_int

      if (nblocks <= 1 .or. bb_input_prevh <= 0) then
         blk = 1
         return
      end if

      hour_int = floor(real(ktau, RKIND) * dt / secsperhour)
      blk = floor(real(mod(hour_int, 24), RKIND) / &
                  real(bb_input_prevh, RKIND)) + 1
      blk = min(max(blk, 1), nblocks)

   end function get_fre_block


   !--------------------------------------------------------------------------
   ! Local solar-hour helper for the optional flaming-fraction cycle.
   !--------------------------------------------------------------------------
   real(RKIND) function get_local_hour(gmt, xlong) result(local_hour)

      implicit none

      real(RKIND), intent(in) :: gmt, xlong
      real(RKIND) :: lon_deg

      lon_deg = xlong
      if (lon_deg > 180.0_RKIND) lon_deg = lon_deg - 360.0_RKIND

      local_hour = modulo(gmt + lon_deg / 15.0_RKIND, 24.0_RKIND)

   end function get_local_hour


   !--------------------------------------------------------------------------
   ! Smooth local-time variation of the flaming fraction.
   !--------------------------------------------------------------------------
   real(RKIND) function get_flam_frac(base_flam_frac, gmt, xlong) result(flam_frac)

      implicit none

      real(RKIND), intent(in) :: base_flam_frac, gmt, xlong
      ! -----------------------------------------------------------
      ! future namelist options
      ! ef_peak_local_hour
      ! ef_diurnal_amplitude
      ! ef_percent_change
      ! ef_flam_floor
      ! ef_flam_ceiling
      ! -------------------------------------------------------

      real(RKIND), parameter :: PI = 3.1415_RKIND
      real(RKIND), parameter :: EF_PEAK_LOCAL_HOUR = 15.0_RKIND
      real(RKIND), parameter :: EF_DIURNAL_AMPLITUDE = 0.20_RKIND
      real(RKIND), parameter :: EF_PERCENT_CHANGE    = 0.0_RKIND
      real(RKIND), parameter :: EF_FLAM_FLOOR   = 0.05_RKIND
      real(RKIND), parameter :: EF_FLAM_CEILING = 0.98_RKIND

      real(RKIND) :: local_hour, phase, adjusted_base

      local_hour = get_local_hour(gmt, xlong)
      ! Apply a relative percentage adjustment to the original
      ! vegetation-dependent flaming fraction
      adjusted_base = base_flam_frac * &
                    (1.0_RKIND + EF_PERCENT_CHANGE / 100.0_RKIND)
      adjusted_base = min(max(adjusted_base, EF_FLAM_FLOOR), &
                     EF_FLAM_CEILING)

       ! Smooth 24-hour cycle, maximum at EF_PEAK_LOCAL_HOUR.
      phase = cos(2.0_RKIND * PI * &
                  (local_hour - EF_PEAK_LOCAL_HOUR) / 24.0_RKIND)

      ! Apply the diurnal perturbation.
      !flam_frac = base_flam_frac + AMPLITUDE * phase
      !flam_frac = min(max(flam_frac, FLAM_MIN), FLAM_MAX)
      flam_frac = adjusted_base + EF_DIURNAL_AMPLITUDE * phase
      flam_frac = min(max(flam_frac, EF_FLAM_FLOOR), &
              EF_FLAM_CEILING)
   end function get_flam_frac


  !-------------------------------------------------------------------------------
  ! Subroutine: calculate_smoke_emissions
  ! Description:
  !   Calculates smoke emissions (ebu) based on Fire Radiative Energy (fre_avg),
  !   emission factors (EFs_map), and target grid area. Emissions are estimated 
  !   for each temporal block, in units of micrograms per second per square meter.
  !   Only grid points with positive FRP are considered; emissions for zero/non-burning
  !   grid points are set to zero. 
  !   fre_avg declared always as a 3-D dummy
  !-------------------------------------------------------------------------------
    subroutine calculate_smoke_emissions(dt, julday, nlcat, EFs_map, fre_avg,       &
                                         ebb_dcycle, EBB, area, nblocks, ktau,      & !JR added EBB
                                         bb_input_prevh, ebu, num_e_bb_in,          &
                                         num_bb_pt,                                 &
                                         myrank,bb_pt_rank,bb_pt_local_cell_idx,         &
                                         ids, ide, jds, jde, kds, kde,              &
                                         ims, ime, jms, jme, kms, kme,              &
                                         its, ite, jts, jte, kts, kte)
  
      implicit none
      INTEGER, INTENT(IN) :: EBB !JR added
      INTEGER, INTENT(IN) :: julday, nlcat, num_e_bb_in, num_bb_pt, myrank
      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN) :: its, ite, jts, jte, kts, kte
       ! Timestep, day, constants
      REAL(RKIND), INTENT(IN) :: dt
      integer,intent(in):: ktau
      INTEGER, INTENT(IN) :: ebb_dcycle, nblocks, bb_input_prevh

      INTEGER, INTENT(IN), DIMENSION(1:num_bb_pt) :: bb_pt_rank,bb_pt_local_cell_idx
      REAL(RKIND), INTENT(IN),   DIMENSION(1:num_bb_pt,nblocks)         :: fre_avg
      REAL(RKIND), INTENT(IN),   DIMENSION(ims:ime, jms:jme)            :: EFs_map
      REAL(RKIND), INTENT(IN),   DIMENSION(ims:ime, jms:jme)            :: area !we need it for first v level 
      REAL(RKIND), INTENT(INOUT),DIMENSION(1:num_bb_pt, kms:kme, num_e_bb_in)   :: ebu

      !Local variables
      !INTEGER :: i, j, blk, hour_int, nv, ii
      INTEGER :: i, j, blk, nv, ii
      REAL(RKIND) ::  hour_tmp
      REAL(RKIND), PARAMETER :: fg_to_ug = 1.0e6
      REAL(RKIND), PARAMETER :: to_s = 3600.0
       
     ! for EBB = 2, blk > 1
      !JR st -> future forecast block selection
      !IF ( ebb_dcycle .eq. 2 ) then
      if (ebb_dcycle == -1 .and. (EBB == 2 .or. EBB == 3)) then
         blk = get_fre_block(ktau,dt,bb_input_prevh,nblocks)
      else
         blk = 1
      endif        
      !JR ends        
    
      !if ebb_dcycle == 2 then ebb_dc1 do not use blocks? should be 24
      j = 1
      DO nv = 1, num_e_bb_in
      DO ii = 1, num_bb_pt
         ! Skip if we are on the wrong rank for this EGU
         if (myrank /= bb_pt_rank(ii)) cycle
         ! Otherwise the index is just the local index
         i = bb_pt_local_cell_idx(ii)
         ! Final check in case the cell wasn't actually on any rank
         if ( i .le. 0 .or. i .gt. ite ) cycle
         IF (fre_avg(ii,blk) > 0.0_RKIND) THEN
            ebu(ii,kts,nv) = (fre_avg(ii,blk) * EFs_map(i,j) * fg_to_ug) / ( area(i,j) * to_s) !we need the area of the grids
         ELSE
            ebu(ii,kts,nv) = 0.0_RKIND
         END IF
      END DO
      END DO
  
    end subroutine calculate_smoke_emissions
  
  !peak_hr only needed at first time step -> from wrapper?
  !-------------------------------------------------------------------------------
  ! Subroutine: diurnal_cycle
  ! Apply fire diurnal cycle emission scaling depending on fire type and time.
  ! Handles multiple fire age "blocks" (nblocks) using coef_bb_dc and fire_hist.
  !-------------------------------------------------------------------------------
  ! Apply diurnal cycle depending on fire type
subroutine diurnal_cycle(  dtstep,dz8w,rho_phy,pi,ebb_min,          &
                           chem,num_chem,julday,gmt,xlat,xlong,     &
                           fire_end_hr,peak_hr,time_int,coef_bb_dc, &
                           fire_hist,hwp,hwp_avg,hwp_prev_day,      &  !I think fire_hist replaced sc_factor
                           vegfrac, eco_id, nblocks,                &
                           lu_nofire, lu_qfire, lu_sfire,           &
                           swdown,ebb_dcycle,ebu,                   &
                           num_e_bb_in,num_bb_pt,                   &
                           myrank,bb_pt_rank,bb_pt_local_cell_idx,  &
                           index_e_bb_in_smoke_fine,                &
                           fire_type,  q_vap, add_fire_moist_flux,  &
                           plume_beta_qv, hwp_alpha,                &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte                )

  implicit none

  INTEGER,      INTENT(IN   ) ::  julday, num_chem, num_e_bb_in,    &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,        &
                                  its,ite, jts,jte, kts,kte
  INTEGER, INTENT(IN) :: myrank, num_bb_pt
  INTEGER, INTENT(IN) :: index_e_bb_in_smoke_fine

  INTEGER, DIMENSION(1:num_bb_pt), INTENT(IN) :: bb_pt_rank,bb_pt_local_cell_idx

  REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem   ! shall we set num_chem=1 here?

  REAL(RKIND), DIMENSION( 1:num_bb_pt, kms:kme, num_e_bb_in ),                 &
         INTENT(INOUT ) ::                                   ebu
  REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: q_vap ! SRB: added q_vap
  INTEGER,      INTENT(IN), DIMENSION(ims:ime,jms:jme) :: eco_id
  REAL(RKIND),  INTENT(IN), DIMENSION(ims:ime,jms:jme,20) :: vegfrac      
 
  REAL(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: hwp_prev_day ! JLS --- TODO
  REAL(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(INOUT)    :: fire_hist
 
  real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: xlat,xlong, swdown
  real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: hwp, peak_hr, fire_end_hr !RAR: Shall we make fire_end integer?
  real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(OUT)    :: lu_nofire, lu_qfire, lu_sfire, coef_bb_dc
  real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: hwp_avg
  real(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: dz8w,rho_phy  !,rel_hum
  real(RKIND), INTENT(IN) ::  dtstep, gmt
  real(RKIND), INTENT(IN) ::  time_int, pi, ebb_min       ! RAR: time in seconds since start of simulation
  INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(OUT) :: fire_type
  integer, INTENT(IN) ::  ebb_dcycle     ! RAR: this is going to be namelist dependent, ebb_dcycle=means 
  real(RKIND), INTENT(IN) :: plume_beta_qv ! SRB: namelist scaling factor for fire qv fluxes

  INTEGER, INTENT(IN) ::  nblocks ! e.g. 1, 4, or 24
  
  !>--local 
  logical, intent(in)  :: add_fire_moist_flux
  integer :: i,j,k,n,m, blk, ii
  integer :: icall=0
  real(RKIND) :: conv_rho, conv, dm_smoke, dc_hwp, dc_gp, dc_fn !daero_num_wfa, daero_num_ifa !, lu_sum1_5, lu_sum12_14

  INTEGER, PARAMETER :: kfire_max=51    ! max vertical level for BB plume rise
  real(RKIND), PARAMETER :: ef_h2o=324.22  ! Emission factor for water vapor
  ! Constants for the fire diurnal cycle calculation ! JLS - needs to be
  ! defined below due to intent(in) of pi
  real(RKIND) :: coef_con
  real(RKIND) :: timeq, fire_age, age_hr, dt1,dt2,dtm         ! For BB emis. diurnal cycle calculation
  
  ! For Gaussian diurnal cycle
  real(RKIND), INTENT(IN) :: hwp_alpha  ! to scale up the wildfire emissions, Jordan please make this a namelist option
  real(RKIND), PARAMETER :: rinti=2.1813936e-8, ax2=3400., const2=130., &
                           coef2=10.6712963e-4, cx2=7200., timeq_max=3600.*24.
  !>-- Fire parameters: Fores west, Forest east, Shrubland, Savannas, Grassland, Cropland
  real(RKIND), dimension(1:5), parameter :: avg_fire_dur   = (/8.9, 4.2, 3.3, 3.0, 1.4/)
  real(RKIND), dimension(1:5), parameter :: sigma_fire_dur = (/8.7, 6.0, 5.5, 5.2, 2.4/)
  ! For fire diurnal cycle calculation
  !real(kind_phys), parameter :: avgx1=-2.0, sigmx1=0.7, C1=0.083  ! Ag fires
  !real(kind_phys), parameter :: avgx2=-0.1, sigmx2=0.8, C2=0.55   ! Grass fires, slash burns
  real(RKIND), parameter :: avgx1=0.,  sigmx1=2.2, C1=0.2  ! Ag fires
  real(RKIND), parameter :: avgx2=0.5, sigmx2=0.8, C2=1.1   ! Grass fires, slash burns

  timeq= gmt*3600._RKIND + real(time_int,4)
  timeq= mod(timeq,timeq_max)
  coef_con=1._RKIND/((2._RKIND*pi)**0.5)

  !peak_hr     = 0.
  fire_type   = 0
  lu_qfire    = 0.
  lu_sfire    = 0.
  lu_nofire   = 0.

  !RAR: change this to the fractional LU type; fire_type: 0- no fires, 1- Ag
  ! or urban fires, 2- prescribed fires in wooded area, 3- wildfires
  j = 1
  !JR st -> future forecast block selection
  !if (ebb_dcycle==2) then
  if (ebb_dcycle == -1) then       !.and. (EBB == 2 .or. EBB == 3)) then
  !JR ends
     do ii = 1, num_bb_pt
        ! Skip if we are on the wrong rank for this EGU
        if (myrank /= bb_pt_rank(ii)) cycle
        ! Otherwise the index is just the local index
        i = bb_pt_local_cell_idx(ii)
        ! Final check in case the cell wasn't actually on any rank
        if ( i .le. 0 .or. i .gt. ite ) cycle
        if (ebu(ii,kts,index_e_bb_in_smoke_fine)<ebb_min) then
           fire_type(i,j) = 0
           lu_nofire(i,j) = 1.0
        else
          ! cropland, urban, cropland/natural mosaic, barren and sparsely
          ! vegetated and non-vegetation areas:
          lu_qfire(i,j) = lu_nofire(i,j) + vegfrac(i,12,j) + vegfrac(i,13,j) + vegfrac(i,14,j) + vegfrac(i,16,j)
          ! Savannas and grassland fires, these fires last longer than the Ag fires:
          lu_sfire(i,j) = lu_qfire(i,j) + vegfrac(i,8,j) + vegfrac(i,9,j) + vegfrac(i,10,j)
          if (lu_nofire(i,j)>0.95) then ! no fires
            fire_type(i,j) = 0
          else if (lu_qfire(i,j)>0.9) then   ! Ag. and urban fires
            fire_type(i,j) = 1
          else if (eco_id(i,j)==8) then
            fire_type(i,j) = 2    ! slash burn and wildfires in the east, eastern temperate forest ecosystem
            ! SRB: Eastern wildland fires if FRE>1E6MJ and eco region is 8 and the fire is not older than 7 hrs
            if ( ebu(ii,kts,index_e_bb_in_smoke_fine)*3600*(1/0.416) .ge. 1.E6 .and. fire_end_hr(i,j) .le. 8 ) then
              fire_type(i,j) = 4
            endif 
          else if (lu_sfire(i,j)>0.8) then
            fire_type(i,j) = 3    ! savanna and grassland fires
          else
            fire_type(i,j) = 4    ! potential wildfires     
          end if
        end if
     end do
  endif


  !JR st -> future forecast block selection
  !if (ebb_dcycle==2) then 
  if (ebb_dcycle == -1) then   !.and. (EBB == 2 .or. EBB == 3)) then
  !JR ends        
    !do blk = 1, nblocks
     do ii = 1, num_bb_pt
        ! Skip if we are on the wrong rank for this EGU
        if (myrank /= bb_pt_rank(ii)) cycle
        ! Otherwise the index is just the local index
        i = bb_pt_local_cell_idx(ii)
        ! Final check in case the cell wasn't actually on any rank
        if ( i .le. 0 .or. i .gt. ite ) cycle
        fire_age= MAX(0.01_RKIND,time_int/3600. + (fire_end_hr(i,j)-1.0))  !One hour delay is due to the latency of the RAVE files, hours

          SELECT CASE ( fire_type(i,j) )   !Ag, urban fires, bare land etc.
          CASE (1)
             ! these fires will have exponentially decreasing diurnal cycle,
             !coef_bb_dc(i,j) = coef_con*1._kind_phys/(sigma_fire_dur(5) *fire_age) *                          &
             !                exp(- ( log(fire_age) - avg_fire_dur(5))**2 /(2._kind_phys*sigma_fire_dur(5)**2 ))
            coef_bb_dc(i,j)= C1/(sigmx1* fire_age)* exp(- (log(fire_age) - avgx1)**2 /(2.*sigmx1**2 ) )

             IF ( dbg_opt .AND. time_int<5000.) then
               WRITE(6,*) 'i,j,peak_hr(i,j),fire_type(i,j) ',i,j,peak_hr(i,j),fire_type(i,j)
               WRITE(6,*) 'coef_bb_dc(i,j) ',coef_bb_dc(i,j)
             END IF

          CASE (2)    ! Savanna and grassland fires, or fires in the eastern US 
            !  coef_bb_dc(i,j) = coef_con*1._kind_phys/(sigma_fire_dur(4) *fire_age) *                          &
            !                  exp(- ( log(fire_age) - avg_fire_dur(4))**2 /(2._kind_phys*sigma_fire_dur(4)**2 ))
            coef_bb_dc(i,j)= C2/(sigmx2* fire_age)* exp(- (log(fire_age) - avgx2)**2 /(2.*sigmx2**2 ) )

              IF ( dbg_opt .AND. time_int<5000.) then
                WRITE(6,*) 'i,j,peak_hr(i,j),fire_type(i,j) ',i,j,peak_hr(i,j),fire_type(i,j)
                WRITE(6,*) 'coef_bb_dc(i,j) ',coef_bb_dc(i,j)
              END IF

          CASE (3,4)    ! wildfires 
             IF (swdown(i,j)<.1 .AND. fire_age> 6. .AND. fire_hist(i,j)>0.75) THEN
                 fire_hist(i,j)= 0.75_RKIND
             ENDIF
             IF (swdown(i,j)<.1 .AND. fire_age> 18. .AND. fire_hist(i,j)>0.5) THEN
                 fire_hist(i,j)= 0.5_RKIND
             ENDIF
             IF (swdown(i,j)<.1 .AND. fire_age> 36. .AND. fire_hist(i,j)>0.25) THEN
                 fire_hist(i,j)= 0.25_RKIND
             ENDIF
             ! this is based on hwp, hourly or instantenous TBD
             dc_hwp= hwp_alpha * hwp(i,j)/ MAX(10._RKIND,hwp_prev_day(i,j)) + (1-hwp_alpha) * (hwp(i,j)/MAX(10._RKIND,hwp_avg(i,j)))
             dc_hwp= MAX(0._RKIND,dc_hwp)
             dc_hwp= MIN(20._RKIND,dc_hwp)

       
             coef_bb_dc(i,j) = fire_hist(i,j)* dc_hwp     ! RAR: scaling factor is applied to the forest fires only, except the eastern US

             IF ( dbg_opt .AND. time_int<5000.) then
               WRITE(6,*) 'i,j,fire_hist(i,j),peak_hr(i,j) ', i,j,fire_hist(i,j),peak_hr(i,j)
               WRITE(6,*) 'dc_gp,dc_hwp,dc_fn,coef_bb_dc(i,j) ',dc_gp,dc_hwp,dc_fn,coef_bb_dc(i,j)
             END IF

          CASE DEFAULT
          END SELECT
       !enddo
     enddo
  endif
end subroutine diurnal_cycle

end module module_fire_emissions
