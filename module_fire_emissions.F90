!> \file module_fire_emissions.F90
!! Online biomass-burning PM2.5 emission calculation.
!! Note for developer
!!   1. Compute a grid-cell PM2.5 emission factor from land-use fractions
!!      and ecoregion information.
!!   2. Convert hourly FRE to a surface smoke_fine emission flux.
!!   3. Add vertically distributed biomass-burning emissions to chemistry
!!      tracers.
!!
!! Plume-rise allocation remains in module_plumerise and is called from
!! mpas_smoke_wrapper between emission-flux calculation and tracer addition.
!! Code written by Johana Romero-Alvarez, CIRES/NOAA GSL
!! johana.romero-alvarez@noaa.gov, smoke.gsl@noaa.gov

module module_fire_emissions

   use mpas_kind_types, only : RKIND
   use mpas_smoke_config, only : epsilc
   use mpas_smoke_init, only : p_smoke_ultrafine, p_smoke_fine,     &
                               p_smoke_coarse,                       &
                               p_no3_a_fine, p_so4_a_fine,          &
                               p_nh4_a_fine,                         &
                               p_co, p_nox, p_ch4, p_so2, p_nh3,    &
                               p_bbvoc

   implicit none

   private

   public :: compute_emission_factors
   public :: compute_gas_emission_factors
   public :: calculate_smoke_emissions_research
   public :: calculate_smoke_emissions_fcst
   public :: calculate_gas_emissions_research
   public :: calculate_gas_emissions_fcst
   public :: add_emis_burn

   public :: num_fire_gas
   public :: id_co, id_ch4, id_nox, id_nh3, id_so2

   real(RKIND), parameter :: gram2ugram = 1.0e6_RKIND
   real(RKIND), parameter :: secsperhour  = 3600.0_RKIND

   ! PM2.5 emission-factor method.
   !   .true.  : category-average Andreae (2019) emission factors
   !   .false. : current flaming/smoldering emission-factor formulation
   ! This is a temporary module-level switch.
   ! I recommend turn this on while we are testing new emission factor
   ! using smoldering and flaming
   logical, parameter :: use_classic_pm25_ef = .false.

   ! PM2.5 emission factors for flaming and smoldering combustion.
   ! Land classes are forest, hardwood, mixed forest, shrubland,
   ! shrub/grass, grass, organic-soil mix, and organic-soil woody.
   real(RKIND), dimension(8), parameter :: EF_FLM = &
      (/19.0_RKIND, 9.4_RKIND, 14.6_RKIND, 9.3_RKIND, &
        10.7_RKIND, 13.3_RKIND, 10.8_RKIND, 11.8_RKIND/)

   real(RKIND), dimension(8), parameter :: EF_SML = &
      (/28.0_RKIND, 37.7_RKIND, 17.6_RKIND, 36.6_RKIND, &
        36.7_RKIND, 38.4_RKIND, 20.7_RKIND, 23.3_RKIND/)

   ! GEUFE emission factors for eastern-U.S. grasslands and croplands.
   real(RKIND), dimension(3), parameter :: EF_EAST = &
      (/9.29_RKIND, 11.91_RKIND, 8.20_RKIND/)

   ! Two main beta factor
   ! Biomass-combustion coefficients.
   ! Order is SA, SAOS, AG, AGOS, EF, EFOS, and peat.
   real(RKIND), dimension(7), parameter :: BETA_KAISER = &
      (/0.78_RKIND, 0.26_RKIND, 0.29_RKIND, 0.13_RKIND, &
        0.49_RKIND, 1.55_RKIND, 5.87_RKIND/)

   real(RKIND), parameter :: BETA_WOOSTER = 0.38_RKIND

   ! Gas-species IDs used by gas_EFs_map and gas_ebu_index.
   integer, parameter :: id_co   = 1
   integer, parameter :: id_ch4  = 2
   integer, parameter :: id_nox  = 3
   integer, parameter :: id_nh3  = 4
   integer, parameter :: id_so2  = 5
   integer, parameter :: num_fire_gas = 5

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

   ! Category-average gas emission factors in g kg-1 dry matter burned.
   ! First dimension is CO, CH4, NOx as NO, NH3, and SO2.
   ! Second dimension is savanna/grassland, tropical forest,
   ! temperate forest, boreal forest, peat, and agricultural residue.
   ! This is from Andrea as well
   real(RKIND), dimension(num_fire_gas, num_fire_cat), parameter :: &
      EF_GAS = reshape((/ &
       69.0_RKIND, 2.7_RKIND, 2.50_RKIND, 0.89_RKIND, 0.47_RKIND, &
      104.0_RKIND, 6.5_RKIND, 2.80_RKIND, 1.33_RKIND, 0.77_RKIND, &
      113.0_RKIND, 5.2_RKIND, 3.00_RKIND, 0.98_RKIND, 0.70_RKIND, &
      121.0_RKIND, 5.5_RKIND, 1.18_RKIND, 2.50_RKIND, 0.75_RKIND, &
      260.0_RKIND, 9.1_RKIND, 1.20_RKIND, 4.20_RKIND, 4.30_RKIND, &
       76.0_RKIND, 5.7_RKIND, 2.40_RKIND, 0.99_RKIND, 0.80_RKIND  &
      /), (/num_fire_gas, num_fire_cat/))

   ! Molecular weights in g mol-1. Andreae reports NOx as NO.
   real(RKIND), dimension(num_fire_gas), parameter :: gas_mw = &
      (/28.010_RKIND, 16.043_RKIND, 30.006_RKIND, &
        17.031_RKIND, 64.066_RKIND/)

   real(RKIND), parameter :: joule2megajoule = 1.0e-6_RKIND

contains

   !--------------------------------------------------------------------------
   ! Compute the grid-cell PM2.5 emission factor.
   !
   ! EFs_map is the vegetation-fraction-weighted product of the PM2.5
   ! emission factor and biomass-combustion coefficient.
   !
   ! efs_smold, efs_flam, efs_rsmold, fmc_avg, wetness, and dt are retained
   ! in the interface for compatibility with the current wrapper. The present
   ! hard-coded emission-factor formulation does not use them.
   ! Note that, we are no longer using the "vegfrac", now we are using landuse
   ! fraction (20 category USGS)
   ! map available ask Minsu Choi
   !--------------------------------------------------------------------------
   subroutine compute_emission_factors(EFs_map, landusef, bb_beta,             &
                                       eco_id, efs_smold, efs_flam,           &
                                       efs_rsmold, fmc_avg, wetness,          &
                                       ef_diurnal_opt, nlcat, dt, gmt,        &
                                       xlong,                                 &
                                       ids, ide, jds, jde, kds, kde,          &
                                       ims, ime, jms, jme, kms, kme,          &
                                       its, ite, jts, jte, kts, kte)

      implicit none

      integer, intent(in) :: bb_beta, nlcat
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte

      logical, intent(in) :: ef_diurnal_opt

      integer, intent(in), dimension(ims:ime, jms:jme) :: eco_id

      real(RKIND), intent(in) :: dt, gmt
      real(RKIND), intent(in), dimension(ims:ime, jms:jme) :: &
         efs_smold, efs_flam, efs_rsmold, fmc_avg, wetness, xlong
      real(RKIND), intent(in), dimension(ims:ime, nlcat, jms:jme) :: &
         landusef

      real(RKIND), intent(out), dimension(ims:ime, jms:jme) :: EFs_map

      integer :: i, j, n
      real(RKIND) :: class_factor

      EFs_map = 0.0_RKIND

      do j = jts, jte
        do i = its, ite
          do n = 1, min(20, nlcat)
            if (use_classic_pm25_ef) then
               call get_classic_pm25_factor(n, eco_id(i,j), bb_beta,        &
                                             class_factor)
            else
               call get_land_class_factor(n, eco_id(i,j), bb_beta,          &
                                          ef_diurnal_opt, gmt, xlong(i,j),  &
                                          class_factor)
            end if
            EFs_map(i,j) = EFs_map(i,j) + &
                           max(0.0_RKIND, landusef(i,n,j)) * class_factor
          end do
          ! Debugging....
!          do n = 1, min(20, nlcat)
!             call get_land_class_factor( n, eco_id(i,j), bb_beta,      &
!                                         .false., gmt, xlong(i,j),      &
!                                         class_factor )
          
!             if (landusef(i,n,j) > 0.0_RKIND) then
!                write(*,*) 'EBB1 EF input ', i, j, n,                  &
!                           eco_id(i,j), landusef(i,n,j), class_factor
!             endif
!             EFs_map(i,j) = EFs_map(i,j) +                  &
!                                 max(0.0_RKIND, landusef(i,n,j)) *      &
!                                 class_factor
!          enddo
!          if (EFs_map(i,j) > 0.0_RKIND) then
!             write(*,*) 'EBB1 EFs_map', i, j, EFs_map(i,j)
!          endif
        end do
      end do

   end subroutine compute_emission_factors


   !--------------------------------------------------------------------------
   ! Compute category-average gas emission-factor maps.
   !
   ! gas_EFs_map has units of g species MJ-1. The existing PM2.5 EFs_map and
   ! its flaming/smoldering formulation are not modified by this routine.
   !--------------------------------------------------------------------------
   subroutine compute_gas_emission_factors(gas_EFs_map, landusef, bb_beta,     &
                                           eco_id, nlcat,                     &
                                           ids, ide, jds, jde, kds, kde,      &
                                           ims, ime, jms, jme, kms, kme,      &
                                           its, ite, jts, jte, kts, kte)

      implicit none

      integer, intent(in) :: bb_beta, nlcat
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte

      integer, intent(in), dimension(ims:ime, jms:jme) :: eco_id
      real(RKIND), intent(in), dimension(ims:ime, nlcat, jms:jme) :: &
         landusef
      real(RKIND), intent(out), &
         dimension(ims:ime, jms:jme, num_fire_gas) :: gas_EFs_map

      integer :: i, j, n
      real(RKIND) :: veg_fraction
      real(RKIND), dimension(num_fire_gas) :: class_factor

      gas_EFs_map = 0.0_RKIND

      do j = jts, jte
         do i = its, ite
            if (eco_id(i,j) == 0) cycle

            do n = 1, min(20, nlcat)
               veg_fraction = max(0.0_RKIND, landusef(i,n,j))
               if (veg_fraction <= 0.0_RKIND) cycle

               call get_gas_class_factor(n, eco_id(i,j), bb_beta,             &
                                         class_factor)

               gas_EFs_map(i,j,:) = gas_EFs_map(i,j,:) + &
                                    veg_fraction * class_factor(:)
            end do
         end do
      end do

   end subroutine compute_gas_emission_factors


   !--------------------------------------------------------------------------
   ! EBB1
   ! Convert the current hourly FRE field to a surface smoke_fine flux.
   ! The resulting ebu has units of micrograms m-2 s-1 and is placed at kts.
   !--------------------------------------------------------------------------
   subroutine calculate_smoke_emissions_research(dt, julday, nlcat, EFs_map,       &
                                             fre_in, area, ktau,               &
                                             bb_input_prevh, ebu,              &
                                             ids, ide, jds, jde, kds, kde,     &
                                             ims, ime, jms, jme, kms, kme,     &
                                             its, ite, jts, jte, kts, kte)

      implicit none

      integer, intent(in) :: julday, nlcat, ktau, bb_input_prevh
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte

      real(RKIND), intent(in) :: dt
      real(RKIND), intent(in), dimension(ims:ime, jms:jme) :: &
         EFs_map, fre_in, area
      real(RKIND), intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: &
         ebu

      call set_surface_smoke_flux(fre_in, EFs_map, area, ebu,                  &
                                  ims, ime, jms, jme, kms, kme,                &
                                  its, ite, jts, jte, kts)

   end subroutine calculate_smoke_emissions_research


   !--------------------------------------------------------------------------
   ! EBB2
   ! Select the current FRE time block and convert it to a surface smoke_fine
   ! flux. The resulting ebu has units of micrograms m-2 s-1 at kts.
   !--------------------------------------------------------------------------
   subroutine calculate_smoke_emissions_fcst(dt, julday, nlcat, EFs_map,       &
                                             fre_avg, area, nblocks, ktau,     &
                                             bb_input_prevh, ebu,              &
                                             ids, ide, jds, jde, kds, kde,     &
                                             ims, ime, jms, jme, kms, kme,     &
                                             its, ite, jts, jte, kts, kte)

      implicit none

      integer, intent(in) :: julday, nlcat, nblocks, ktau, bb_input_prevh
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte

      real(RKIND), intent(in) :: dt
      real(RKIND), intent(in), dimension(ims:ime, jms:jme) :: EFs_map, area
      real(RKIND), intent(in), dimension(ims:ime, jms:jme, nblocks) :: fre_avg
      real(RKIND), intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: &
         ebu

      integer :: blk

      blk = get_fre_block(ktau, dt, bb_input_prevh, nblocks)

      call set_surface_smoke_flux(fre_avg(:,:,blk), EFs_map, area, ebu,        &
                                  ims, ime, jms, jme, kms, kme,                &
                                  its, ite, jts, jte, kts)

   end subroutine calculate_smoke_emissions_fcst


   !--------------------------------------------------------------------------
   ! EBB1 gas emissions.
   ! Convert current-hour FRE to mol m-2 s-1 for the supported gas species.
   !--------------------------------------------------------------------------
   subroutine calculate_gas_emissions_research(dt, julday, gas_EFs_map, fre_in,   &
                                           area, ktau, bb_input_prevh,        &
                                           ebu, num_e_bb_in, gas_ebu_index,   &
                                           ids, ide, jds, jde, kds, kde,      &
                                           ims, ime, jms, jme, kms, kme,      &
                                           its, ite, jts, jte, kts, kte)

      implicit none

      integer, intent(in) :: julday, ktau, bb_input_prevh, num_e_bb_in
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte
      integer, intent(in), dimension(num_fire_gas) :: gas_ebu_index

      real(RKIND), intent(in) :: dt
      real(RKIND), intent(in), dimension(ims:ime, jms:jme) :: fre_in, area
      real(RKIND), intent(in), &
         dimension(ims:ime, jms:jme, num_fire_gas) :: gas_EFs_map
      real(RKIND), intent(inout), &
         dimension(ims:ime, kms:kme, jms:jme, num_e_bb_in) :: ebu

      call set_surface_gas_flux(fre_in, gas_EFs_map, area, ebu,               &
                                num_e_bb_in, gas_ebu_index,                   &
                                ims, ime, jms, jme, kms, kme,                 &
                                its, ite, jts, jte, kts)

   end subroutine calculate_gas_emissions_research


   !--------------------------------------------------------------------------
   ! EBB2 gas emissions. This follows the existing PM2.5 block selection.
   !--------------------------------------------------------------------------
   subroutine calculate_gas_emissions_fcst(dt, julday, gas_EFs_map, fre_avg, &
                                           area, nblocks, ktau,               &
                                           bb_input_prevh, ebu,               &
                                           num_e_bb_in, gas_ebu_index,        &
                                           ids, ide, jds, jde, kds, kde,      &
                                           ims, ime, jms, jme, kms, kme,      &
                                           its, ite, jts, jte, kts, kte)

      implicit none

      integer, intent(in) :: julday, nblocks, ktau, bb_input_prevh
      integer, intent(in) :: num_e_bb_in
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte
      integer, intent(in), dimension(num_fire_gas) :: gas_ebu_index

      real(RKIND), intent(in) :: dt
      real(RKIND), intent(in), dimension(ims:ime, jms:jme) :: area
      real(RKIND), intent(in), dimension(ims:ime, jms:jme, nblocks) :: fre_avg
      real(RKIND), intent(in), &
         dimension(ims:ime, jms:jme, num_fire_gas) :: gas_EFs_map
      real(RKIND), intent(inout), &
         dimension(ims:ime, kms:kme, jms:jme, num_e_bb_in) :: ebu

      integer :: blk

      blk = get_fre_block(ktau, dt, bb_input_prevh, nblocks)

      call set_surface_gas_flux(fre_avg(:,:,blk), gas_EFs_map, area, ebu,     &
                                num_e_bb_in, gas_ebu_index,                   &
                                ims, ime, jms, jme, kms, kme,                 &
                                its, ite, jts, jte, kts)

   end subroutine calculate_gas_emissions_fcst


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

      case (18)
         base_flam_frac = 0.70_RKIND
         call mixed_combustion_ef(base_flam_frac, EF_FLM(8), EF_SML(8),        &
                                  ef_diurnal_opt, gmt, xlong, ef)
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
   ! Return gas emission-factor and beta products for one IGBP land class.
   ! Gas EFs use category-average values and do not use the PM2.5
   ! flaming/smoldering fractions.
   !--------------------------------------------------------------------------
   subroutine get_gas_class_factor(land_class, eco, bb_beta, class_factor)

      implicit none

      integer, intent(in) :: land_class, eco, bb_beta
      real(RKIND), intent(out), dimension(num_fire_gas) :: class_factor

      integer :: fire_category
      real(RKIND) :: beta_val

      class_factor = 0.0_RKIND
      fire_category = get_fire_category(land_class, eco)
      if (fire_category == cat_none) return

      call get_gas_beta_factor(land_class, eco, bb_beta, beta_val)
      if (beta_val <= 0.0_RKIND) return

      class_factor(:) = EF_GAS(:,fire_category) * beta_val

   end subroutine get_gas_class_factor


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
   ! Convert an hourly FRE field to the surface smoke_fine emission flux.
   !--------------------------------------------------------------------------
   subroutine set_surface_smoke_flux(fre, EFs_map, area, ebu,                  &
                                     ims, ime, jms, jme, kms, kme,             &
                                     its, ite, jts, jte, kts)

      implicit none

      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts

      real(RKIND), intent(in), dimension(ims:ime, jms:jme) :: &
         fre, EFs_map, area
      real(RKIND), intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: &
         ebu

      integer :: i, j

      do j = jts, jte
         do i = its, ite
            if (fre(i,j) > 0.0_RKIND .and. area(i,j) > 0.0_RKIND) then
               ebu(i,kts,j) = fre(i,j) * EFs_map(i,j) * 1.0e-6_RKIND * &
                              gram2ugram / &
                              (area(i,j) * secsperhour)
            else
               ebu(i,kts,j) = 0.0_RKIND
            end if
!            if (fre(i,j) > 0.0_RKIND) then
!              write(*,*) 'EBB1 flux ', i, j,     &
!              fre(i,j), EFs_map(i,j), area(i,j), &
!              ebu(i,kts,j)
!            endif
         end do
      end do

   end subroutine set_surface_smoke_flux


   !--------------------------------------------------------------------------
   ! Convert hourly FRE to surface gas emission fluxes in mol m-2 s-1.
   ! FRE is supplied in J and gas_EFs_map in g MJ-1.
   !--------------------------------------------------------------------------
   subroutine set_surface_gas_flux(fre, gas_EFs_map, area, ebu,               &
                                   num_e_bb_in, gas_ebu_index,                &
                                   ims, ime, jms, jme, kms, kme,              &
                                   its, ite, jts, jte, kts)

      implicit none

      integer, intent(in) :: num_e_bb_in
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts
      integer, intent(in), dimension(num_fire_gas) :: gas_ebu_index

      real(RKIND), intent(in), dimension(ims:ime, jms:jme) :: fre, area
      real(RKIND), intent(in), &
         dimension(ims:ime, jms:jme, num_fire_gas) :: gas_EFs_map
      real(RKIND), intent(inout), &
         dimension(ims:ime, kms:kme, jms:jme, num_e_bb_in) :: ebu

      integer :: i, j, gas_id, nv

      do gas_id = 1, num_fire_gas
         nv = gas_ebu_index(gas_id)
         if (nv < 1 .or. nv > num_e_bb_in) cycle

         do j = jts, jte
            do i = its, ite
               ebu(i,kts,j,nv) = 0.0_RKIND

               if (fre(i,j) > 0.0_RKIND .and. area(i,j) > 0.0_RKIND) then
                  ebu(i,kts,j,nv) =                                  &
                     fre(i,j) * gas_EFs_map(i,j,gas_id) *             &
                     joule2megajoule / gas_mw(gas_id) /               &
                     (area(i,j) * secsperhour)
               end if
            end do
         end do
      end do

   end subroutine set_surface_gas_flux


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
   ! Optional smooth local-time variation of the flaming fraction.
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


   !--------------------------------------------------------------------------
   ! Entire module move to subroutine 
   !--------------------------------------------------------------------------
   subroutine add_emis_burn(dtstep, dz8w, rho_phy, pi, ebb_min,               &
                            chem, num_chem, julday, gmt, xlat, xlong,         &
                            fire_end_hr, peak_hr, time_int,                   &
                            coef_bb_dc, fire_hist, hwp, hwp_prevd,            &
                            swdown, ebb_dcycle, EBB, ebu, num_e_bb_in,        &
                            fire_type, q_vap, add_fire_moist_flux,            &
                            do_mpas_sna,                                      &
                            sc_factor,                                        &
                            index_e_bb_in_smoke_ultrafine,                    &
                            index_e_bb_in_smoke_fine,                         &
                            index_e_bb_in_smoke_coarse,                       &
                            index_e_bb_in_co, index_e_bb_in_nh3,              &
                            index_e_bb_in_ch4,                                &
                            index_e_bb_in_nox, index_e_bb_in_so2,             &
                            index_e_bb_in_voc,                                &
                            ids, ide, jds, jde, kds, kde,                     &
                            ims, ime, jms, jme, kms, kme,                     &
                            its, ite, jts, jte, kts, kte                      )

      implicit none

      integer, intent(in) :: julday, num_chem, num_e_bb_in
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte

      integer, intent(in) :: index_e_bb_in_smoke_ultrafine
      integer, intent(in) :: index_e_bb_in_smoke_fine
      integer, intent(in) :: index_e_bb_in_smoke_coarse
      integer, intent(in) :: index_e_bb_in_co, index_e_bb_in_nh3
      integer, intent(in) :: index_e_bb_in_ch4
      integer, intent(in) :: index_e_bb_in_nox, index_e_bb_in_so2
      integer, intent(in) :: index_e_bb_in_voc

      real(RKIND), intent(inout), dimension(ims:ime, kms:kme, jms:jme,        &
                                             1:num_chem) :: chem
      real(RKIND), intent(inout), dimension(ims:ime, kms:kme, jms:jme,        &
                                             num_e_bb_in) :: ebu
      real(RKIND), intent(inout), dimension(ims:ime, kms:kme, jms:jme) ::     &
         q_vap

      real(RKIND), intent(in), dimension(ims:ime, jms:jme) ::                 &
         xlat, xlong, swdown
      real(RKIND), intent(in), dimension(ims:ime, jms:jme) ::                 &
         hwp, peak_hr, fire_end_hr
      real(RKIND), intent(inout), dimension(ims:ime, jms:jme) :: coef_bb_dc
      real(RKIND), intent(in), dimension(ims:ime, jms:jme) :: hwp_prevd
      real(RKIND), intent(in), dimension(ims:ime, kms:kme, jms:jme) ::        &
         dz8w, rho_phy
      integer, intent(in), dimension(ims:ime, jms:jme) :: fire_type
      real(RKIND), intent(inout), dimension(ims:ime, jms:jme) :: fire_hist

      real(RKIND), intent(in) :: dtstep, gmt
      real(RKIND), intent(in) :: time_int, pi, ebb_min
      real(RKIND), intent(in) :: sc_factor
      integer, intent(in) :: ebb_dcycle
      integer, intent(in) :: EBB
      logical, intent(in) :: add_fire_moist_flux, do_mpas_sna

      logical :: do_bb_sna_speciation
      integer :: i, j, k
      integer :: icall = 0
      real(RKIND) :: conv, conv_gas
      real(RKIND) :: dm_smoke_fine, dm_smoke_ultrafine
      real(RKIND) :: dm_smoke_coarse
      real(RKIND) :: dm_co, dm_nox, dm_so2, dm_nh3, dm_ch4, dm_bbvoc
      real(RKIND) :: ext2

      integer, parameter :: kfire_max = 51
      real(RKIND), parameter :: ef_h2o = 324.22_RKIND
      real(RKIND), parameter :: sc_me = 4.0_RKIND
      real(RKIND), parameter :: ab_me = 0.5_RKIND
      real(RKIND), parameter :: frac_bb_no3 = 0.03_RKIND
      real(RKIND), parameter :: frac_bb_so4 = 0.01_RKIND
      real(RKIND), parameter :: frac_bb_nh4 = 0.03_RKIND
      real(RKIND), parameter :: frac_bb_smoke =                         &
         1.0_RKIND - frac_bb_no3 - frac_bb_so4 - frac_bb_nh4

      if (mod(int(time_int), 1800) == 0) then
         icall = 0
      end if

      do_bb_sna_speciation = do_mpas_sna .and.                         &
                              p_no3_a_fine > 0 .and.                    &
                              p_so4_a_fine > 0 .and.                    &
                              p_nh4_a_fine > 0

      ext2 = sc_me + ab_me

      do j = jts, jte
         do i = its, ite
            do k = kts, min(kte, kfire_max)

               if (ebu(i,k,j,index_e_bb_in_smoke_fine) < ebb_min) cycle

               !JR st --> ebb_dcycle=2 no longer exists
               if (ebb_dcycle == 1 .or. ebb_dcycle == -1) then
               
                  conv = dtstep / (rho_phy(i,k,j) * dz8w(i,k,j))
                  conv_gas = dtstep * 0.02897_RKIND /                    &
                             (rho_phy(i,k,j) * dz8w(i,k,j))
               !JR st      
               !else if (ebb_dcycle == 2) then
               elseif (ebb_dcycle == -1 .and. (EBB == 2 .or. EBB == 3)) then
               !JR ends        
                  conv = coef_bb_dc(i,j) * dtstep /                      &
                         (rho_phy(i,k,j) * dz8w(i,k,j))
                  conv_gas = coef_bb_dc(i,j) * dtstep * 0.02897_RKIND /  &
                             (rho_phy(i,k,j) * dz8w(i,k,j))
               end if

               conv = conv * sc_factor
               conv_gas = conv_gas * sc_factor

               dm_smoke_fine =                                         &
                  conv * ebu(i,k,j,index_e_bb_in_smoke_fine)

               if (do_bb_sna_speciation) then
                  chem(i,k,j,p_smoke_fine) =                            &
                     chem(i,k,j,p_smoke_fine) +                         &
                     frac_bb_smoke * dm_smoke_fine

                  chem(i,k,j,p_no3_a_fine) =                            &
                     chem(i,k,j,p_no3_a_fine) +                         &
                     frac_bb_no3 * dm_smoke_fine

                  chem(i,k,j,p_so4_a_fine) =                            &
                     chem(i,k,j,p_so4_a_fine) +                         &
                     frac_bb_so4 * dm_smoke_fine

                  chem(i,k,j,p_nh4_a_fine) =                            &
                     chem(i,k,j,p_nh4_a_fine) +                         &
                     frac_bb_nh4 * dm_smoke_fine

                  chem(i,k,j,p_no3_a_fine) =                            &
                     min(max(chem(i,k,j,p_no3_a_fine), epsilc),         &
                         5.0e3_RKIND)
                  chem(i,k,j,p_so4_a_fine) =                            &
                     min(max(chem(i,k,j,p_so4_a_fine), epsilc),         &
                         5.0e3_RKIND)
                  chem(i,k,j,p_nh4_a_fine) =                            &
                     min(max(chem(i,k,j,p_nh4_a_fine), epsilc),         &
                         5.0e3_RKIND)
               else
                  chem(i,k,j,p_smoke_fine) =                            &
                     chem(i,k,j,p_smoke_fine) + dm_smoke_fine
               end if

               chem(i,k,j,p_smoke_fine) =                              &
                  min(max(chem(i,k,j,p_smoke_fine), epsilc),            &
                      5.0e3_RKIND)

!               write(*,*) 'EBB1 smoke_fine ', i, j, k,                 &
!               ebu(i,k,j,index_e_bb_in_smoke_fine),         &
!               dm_smoke_fine,                               &
!               chem(i,k,j,p_smoke_fine)

               if (p_smoke_ultrafine > 0) then
                  dm_smoke_ultrafine =                                  &
                     conv * ebu(i,k,j,index_e_bb_in_smoke_ultrafine)
                  chem(i,k,j,p_smoke_ultrafine) =                       &
                     chem(i,k,j,p_smoke_ultrafine) + dm_smoke_ultrafine
                  chem(i,k,j,p_smoke_ultrafine) =                       &
                     min(max(chem(i,k,j,p_smoke_ultrafine), epsilc),    &
                         5.0e3_RKIND)
               end if

               if (p_smoke_coarse > 0) then
                  dm_smoke_coarse =                                     &
                     conv * ebu(i,k,j,index_e_bb_in_smoke_coarse)
                  chem(i,k,j,p_smoke_coarse) =                          &
                     chem(i,k,j,p_smoke_coarse) + dm_smoke_coarse
                  chem(i,k,j,p_smoke_coarse) =                          &
                     min(max(chem(i,k,j,p_smoke_coarse), epsilc),       &
                         5.0e3_RKIND)
               end if

               if (p_co > 0) then
                  dm_co = conv_gas * ebu(i,k,j,index_e_bb_in_co)
                  chem(i,k,j,p_co) = chem(i,k,j,p_co) + dm_co
                  chem(i,k,j,p_co) =                                    &
                     min(max(chem(i,k,j,p_co), epsilc), 5.0e3_RKIND)
               end if

               if (p_nox > 0) then
                  dm_nox = conv_gas * ebu(i,k,j,index_e_bb_in_nox)
                  chem(i,k,j,p_nox) = chem(i,k,j,p_nox) + dm_nox
                  chem(i,k,j,p_nox) =                                   &
                     min(max(chem(i,k,j,p_nox), epsilc), 5.0e3_RKIND)
               end if

               if (p_ch4 > 0) then
                  dm_ch4 = conv_gas * ebu(i,k,j,index_e_bb_in_ch4)
                  chem(i,k,j,p_ch4) = chem(i,k,j,p_ch4) + dm_ch4
                  chem(i,k,j,p_ch4) =                                   &
                     min(max(chem(i,k,j,p_ch4), epsilc), 5.0e3_RKIND)
               end if

               if (p_so2 > 0) then
                  dm_so2 = conv_gas * ebu(i,k,j,index_e_bb_in_so2)
                  chem(i,k,j,p_so2) = chem(i,k,j,p_so2) + dm_so2
                  chem(i,k,j,p_so2) =                                   &
                     min(max(chem(i,k,j,p_so2), epsilc), 5.0e3_RKIND)
               end if

               if (p_nh3 > 0) then
                  dm_nh3 = conv_gas * ebu(i,k,j,index_e_bb_in_nh3)
                  chem(i,k,j,p_nh3) = chem(i,k,j,p_nh3) + dm_nh3
                  chem(i,k,j,p_nh3) =                                   &
                     min(max(chem(i,k,j,p_nh3), epsilc), 5.0e3_RKIND)
               end if

               if (p_bbvoc > 0) then
                  dm_bbvoc = conv_gas * ebu(i,k,j,index_e_bb_in_voc)
                  chem(i,k,j,p_bbvoc) = chem(i,k,j,p_bbvoc) + dm_bbvoc
                  chem(i,k,j,p_bbvoc) =                                 &
                     min(max(chem(i,k,j,p_bbvoc), epsilc), 5.0e3_RKIND)
               end if

               if (add_fire_moist_flux) then
                  q_vap(i,k,j) = q_vap(i,k,j) +                         &
                     dm_smoke_fine * ef_h2o * 1.0e-9_RKIND
                  q_vap(i,k,j) =                                        &
                     min(max(q_vap(i,k,j), 0.0_RKIND), 1.0e3_RKIND)
               end if

            end do

            icall = icall + 1

         end do
      end do

   end subroutine add_emis_burn

end module module_fire_emissions
