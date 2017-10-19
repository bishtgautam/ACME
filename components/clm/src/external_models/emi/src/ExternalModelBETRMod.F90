module ExternalModelBETRMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides wrapper for BeTR in ALM
  !
  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use ExternalModelInterfaceDataMod, only : emi_data_list, emi_data
  use ExternalModelBaseType        , only : em_base_type
  use decompMod                    , only : bounds_type
  use BeTRSimulationALM            , only : betr_simulation_alm_type

  !
  implicit none
  !

  type, public, extends(em_base_type) :: em_betr_type

     ! ----------------------------------------------------------------------
     ! Indicies required during the initialization
     ! ----------------------------------------------------------------------
     integer :: index_l2e_init_max_patch_per_col

     integer :: index_l2e_init_patch_active
     integer :: index_l2e_init_patch_itype
     integer :: index_l2e_init_patch_wt_col

     integer :: index_l2e_init_col_active
     integer :: index_l2e_init_col_landunit
     integer :: index_l2e_init_col_snl
     integer :: index_l2e_init_col_zi
     integer :: index_l2e_init_col_dz
     integer :: index_l2e_init_col_z
     integer :: index_l2e_init_col_pfti
     integer :: index_l2e_init_col_pftf
     integer :: index_l2e_init_col_npft

     integer :: index_l2e_init_landunit_itype

     integer :: index_l2e_init_finundated_col
     integer :: index_l2e_init_frac_h2osfc_col

     integer :: index_l2e_init_h2osoi_liq_col
     integer :: index_l2e_init_h2osoi_ice_col
     integer :: index_l2e_init_h2osoi_liqvol_col
     integer :: index_l2e_init_h2osoi_icevol_col
     integer :: index_l2e_init_h2osoi_vol_col
     integer :: index_l2e_init_air_vol_col
     integer :: index_l2e_init_rho_vap_col
     integer :: index_l2e_init_rhvap_soi_col
     integer :: index_l2e_init_smp_l_col
     
     ! ----------------------------------------------------------------------
     ! Indicies required during timestepping
     ! ----------------------------------------------------------------------

     integer :: index_l2e_col_snl
     integer :: index_l2e_col_zi
     integer :: index_l2e_col_dz
     integer :: index_l2e_col_z
     integer :: index_l2e_col_pfti
     integer :: index_l2e_col_pftf
     integer :: index_l2e_col_npft

     integer :: index_l2e_state_frac_h2osfc
     integer :: index_l2e_state_finundated
     integer :: index_l2e_state_h2osoi_liq
     integer :: index_l2e_state_h2osoi_ice
     integer :: index_l2e_state_h2osoi_liqvol
     integer :: index_l2e_state_h2osoi_icevol
     integer :: index_l2e_state_h2osoi_vol
     integer :: index_l2e_state_air_vol
     integer :: index_l2e_state_rho_vap
     integer :: index_l2e_state_rhvap_soi
     integer :: index_l2e_state_smp_l

     integer :: index_l2e_annsum_npp
     integer :: index_l2e_agnpp
     integer :: index_l2e_bgnpp

     integer :: index_l2e_flux_infl
     integer :: index_l2e_flux_totdrain
     integer :: index_l2e_flux_gross_evap_soil
     integer :: index_l2e_flux_gross_infl_soil
     integer :: index_l2e_flux_surf
     integer :: index_l2e_flux_dew_grnd
     integer :: index_l2e_flux_dew_snow
     integer :: index_l2e_flux_sub_snow_vol
     integer :: index_l2e_flux_sub_snow
     integer :: index_l2e_flux_h2osfc2topsoi
     integer :: index_l2e_flux_snow2topsoi
     integer :: index_l2e_flux_rootsoi
     integer :: index_l2e_flux_adv
     integer :: index_l2e_flux_drain_vr
     integer :: index_l2e_flux_tran_veg
     integer :: index_l2e_flux_rootsoi_frac

     integer :: index_l2e_state_t_soi10cm
     integer :: index_l2e_state_t_soil_nlevsoi
     integer :: index_l2e_state_t_veg

     integer :: index_l2e_state_qcharge
     integer :: index_l2e_state_fracice

     integer :: index_l2e_state_altmax
     integer :: index_l2e_state_altmax_lastyear
     integer :: index_l2e_state_lbl_rsc_h2o
     integer :: index_l2e_state_elai

     integer :: index_l2e_state_forc_pbot_downscaled
     integer :: index_l2e_state_forc_t_downscaled

     integer :: index_l2e_state_soil_ph

     integer :: index_l2e_parameter_cellorg
     integer :: index_l2e_parameter_cellclay
     integer :: index_l2e_parameter_cellsand
     integer :: index_l2e_parameter_bd
     integer :: index_l2e_parameter_watfc
     integer :: index_l2e_parameter_rootfr_patch
     integer :: index_l2e_parameter_watsatc
     integer :: index_l2e_parameter_bswc
     integer :: index_l2e_parameter_sucsatc
     integer :: index_l2e_parameter_effporosityc

     integer :: index_l2e_filter_nolakec
     integer :: index_l2e_filter_num_nolakec

     class(betr_simulation_alm_type), pointer            :: betr

   contains
     procedure, public :: Populate_L2E_Init_List  => EM_BeTR_Populate_L2E_Init_List
     !procedure, public :: Populate_E2L_Init_List  => EM_BeTR_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EM_BeTR_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_BeTR_Populate_E2L_List
     procedure, public :: Init                    => EM_BeTR_Init
     procedure, public :: Solve                   => EM_Betr_Solve
  end type em_betr_type

contains

  !------------------------------------------------------------------------
  subroutine EM_BeTR_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by VSFM from ALM
    !
    ! !USES:
    use ExternalModelConstants , only : EM_INITIALIZATION_STAGE
    use ExternalModelConstants , only : L2E_VAR_MAX_PATCH_PER_COL

    use ExternalModelConstants , only : L2E_PATCH_ACTIVE
    use ExternalModelConstants , only : L2E_PATCH_TYPE
    use ExternalModelConstants , only : L2E_PATCH_WT_COL

    use ExternalModelConstants , only : L2E_COLUMN_ACTIVE
    use ExternalModelConstants , only : L2E_COLUMN_LANDUNIT_INDEX
    use ExternalModelConstants , only : L2E_COLUMN_NUM_SNOW_LAYERS
    use ExternalModelConstants , only : L2E_COLUMN_ZI
    use ExternalModelConstants , only : L2E_COLUMN_DZ
    use ExternalModelConstants , only : L2E_COLUMN_Z
    use ExternalModelConstants , only : L2E_COLUMN_PATCH_INDEX_BEGIN
    use ExternalModelConstants , only : L2E_COLUMN_PATCH_INDEX_END
    use ExternalModelConstants , only : L2E_COLUMN_NUM_PATCH

    use ExternalModelConstants , only : L2E_LANDUNIT_TYPE

    use ExternalModelConstants , only : L2E_STATE_FRAC_INUNDATED
    use ExternalModelConstants , only : L2E_STATE_FRAC_H2OSFC
    use ExternalModelConstants , only : L2E_STATE_H2OSOI_LIQ_NLEVSOI
    use ExternalModelConstants , only : L2E_STATE_H2OSOI_ICE_NLEVSOI
    use ExternalModelConstants , only : L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI
    use ExternalModelConstants , only : L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI
    use ExternalModelConstants , only : L2E_STATE_H2OSOI_VOL_NLEVSOI
    use ExternalModelConstants , only : L2E_STATE_AIR_VOL_NLEVSOI
    use ExternalModelConstants , only : L2E_STATE_RHO_VAP_NLEVSOI
    use ExternalModelConstants , only : L2E_STATE_RHVAP_SOI_NLEVSOI
    use ExternalModelConstants , only : L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI

    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    ! ALMVar
    id = L2E_VAR_MAX_PATCH_PER_COL
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_max_patch_per_col = index

    ! Patch-level
    id = L2E_PATCH_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_patch_active = index

    id = L2E_PATCH_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_patch_itype = index

    id = L2E_PATCH_WT_COL
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_patch_wt_col = index

    ! Column
    id = L2E_COLUMN_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_active = index

    id = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_landunit = index

    id = L2E_COLUMN_NUM_SNOW_LAYERS
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_snl = index

    id = L2E_COLUMN_ZI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_zi = index

    id = L2E_COLUMN_DZ
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_dz = index

    id = L2E_COLUMN_Z
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_z = index

    id = L2E_COLUMN_PATCH_INDEX_BEGIN
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_pfti = index

    id = L2E_COLUMN_PATCH_INDEX_END
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_pftf = index

    id = L2E_COLUMN_NUM_PATCH
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_npft = index

    ! Landunit-level
    id = L2E_LANDUNIT_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_itype = index

    ! Waterstate
    id = L2E_STATE_FRAC_INUNDATED
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_finundated_col = index

    id = L2E_STATE_FRAC_H2OSFC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_frac_h2osfc_col = index

    id = L2E_STATE_H2OSOI_LIQ_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_liq_col = index

    id = L2E_STATE_H2OSOI_ICE_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_ice_col = index

    id = L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_liqvol_col = index

    id = L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_icevol_col = index

    id = L2E_STATE_H2OSOI_VOL_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_vol_col = index

    id = L2E_STATE_AIR_VOL_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_air_vol_col = index

    id = L2E_STATE_RHO_VAP_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_rho_vap_col = index

    id = L2E_STATE_RHVAP_SOI_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_rhvap_soi_col = index

    id = L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_smp_l_col = index

    deallocate(em_stages)

  end subroutine EM_BeTR_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list

    call EM_BETR_Populate_L2E_List_WaterState_Vars(this, l2e_list)
    call EM_BETR_Populate_L2E_List_Column_Vars(this, l2e_list)
    call EM_BETR_Populate_L2E_List_WaterFlux_Vars(this, l2e_list)
    call EM_BETR_Populate_L2E_List_Temperature_Vars(this, l2e_list)
    call EM_BETR_Populate_L2E_List_SoilHydrology_Vars(this, l2e_list)
    call EM_BETR_Populate_L2E_List_CanopyState_Vars(this, l2e_list)
    call EM_BETR_Populate_L2E_List_Atm2Lnd_Vars(this, l2e_list)
    call EM_BETR_Populate_L2E_List_ChemState_Vars(this, l2e_list)
    call EM_BETR_Populate_L2E_List_SoilState_Vars(this, l2e_list)


  end subroutine EM_BETR_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List_WaterState_Vars(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    use ExternalModelConstants    , only : EM_BETR_BEGIN_MASS_BALANCE_STAGE
    use ExternalModelConstants    , only : EM_BETR_PRE_DIAG_WATER_FLUX_STAGE
    use ExternalModelConstants    , only : EM_BETR_PRE_DIAG_DTRACER_FREEZE_THAW_STAGE
    use ExternalModelConstants    , only : L2E_STATE_FRAC_H2OSFC
    use ExternalModelConstants    , only : L2E_STATE_FRAC_INUNDATED
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_LIQ_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_ICE_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_AIR_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_RHO_VAP_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_RHVAP_SOI_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
    use ExternalModelConstants    , only : L2E_FILTER_NOLAKEC
    use ExternalModelConstants    , only : L2E_FILTER_NUM_NOLAKEC
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 2
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_PRE_DIAG_WATER_FLUX_STAGE
    em_stages(2) = EM_BETR_PRE_DIAG_DTRACER_FREEZE_THAW_STAGE

    id                            = L2E_STATE_FRAC_H2OSFC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_frac_h2osfc   = index

    id                            = L2E_STATE_FRAC_INUNDATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_finundated    = index

    id                            = L2E_STATE_H2OSOI_LIQ_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq    = index

    id                            = L2E_STATE_H2OSOI_ICE_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice    = index

    id                            = L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liqvol = index

    id                            = L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_icevol = index

    id                            = L2E_STATE_H2OSOI_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_vol    = index

    id                            = L2E_STATE_AIR_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_air_vol       = index

    id                            = L2E_STATE_RHO_VAP_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_rho_vap       = index

    id                            = L2E_STATE_RHVAP_SOI_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_rhvap_soi     = index

    id                            = L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_smp_l         = index

    id                            = L2E_FILTER_NOLAKEC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_nolakec      = index

    id                            = L2E_FILTER_NUM_NOLAKEC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num_nolakec  = index

    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List_WaterState_Vars

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List_Column_Vars(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    use ExternalModelConstants    , only : EM_BETR_PRE_DIAG_DTRACER_FREEZE_THAW_STAGE
    use ExternalModelConstants    , only : L2E_COLUMN_NUM_SNOW_LAYERS
    use ExternalModelConstants    , only : L2E_COLUMN_ZI
    use ExternalModelConstants    , only : L2E_COLUMN_DZ
    use ExternalModelConstants    , only : L2E_COLUMN_Z
    use ExternalModelConstants    , only : L2E_COLUMN_PATCH_INDEX_BEGIN
    use ExternalModelConstants    , only : L2E_COLUMN_PATCH_INDEX_END
    use ExternalModelConstants    , only : L2E_COLUMN_NUM_PATCH
    use ExternalModelConstants    , only : L2E_COLUMN_NUM_PATCH
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_PRE_DIAG_DTRACER_FREEZE_THAW_STAGE

    id = L2E_COLUMN_NUM_SNOW_LAYERS
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_snl = index

    id = L2E_COLUMN_ZI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_zi = index

    id = L2E_COLUMN_DZ
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_dz = index

    id = L2E_COLUMN_Z
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_z = index

    id = L2E_COLUMN_PATCH_INDEX_BEGIN
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_pfti = index

    id = L2E_COLUMN_PATCH_INDEX_END
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_pftf = index

    id = L2E_COLUMN_NUM_PATCH
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_npft = index

    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List_Column_Vars

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List_WaterFlux_Vars(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    use ExternalModelConstants    , only : EM_BETR_STEP_WITHOUT_DRAINGE_STAGE
    use ExternalModelConstants    , only : L2E_FLUX_INFL
    use ExternalModelConstants    , only : L2E_FLUX_TOTDRAIN
    use ExternalModelConstants    , only : L2E_FLUX_GROSS_EVAP_SOIL
    use ExternalModelConstants    , only : L2E_FLUX_GROSS_INFL_SOIL
    use ExternalModelConstants    , only : L2E_FLUX_SURF
    use ExternalModelConstants    , only : L2E_FLUX_DEW_GRND
    use ExternalModelConstants    , only : L2E_FLUX_DEW_SNOW
    use ExternalModelConstants    , only : L2E_FLUX_SUB_SNOW_VOL
    use ExternalModelConstants    , only : L2E_FLUX_SUB_SNOW
    use ExternalModelConstants    , only : L2E_FLUX_H2OSFC2TOPSOI
    use ExternalModelConstants    , only : L2E_FLUX_SNOW2TOPSOI
    use ExternalModelConstants    , only : L2E_FLUX_ROOTSOI
    use ExternalModelConstants    , only : L2E_FLUX_ADV
    use ExternalModelConstants    , only : L2E_FLUX_DRAIN_VR
    use ExternalModelConstants    , only : L2E_FLUX_TRAN_VEG
    use ExternalModelConstants    , only : L2E_FLUX_ROOTSOI_FRAC
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_STEP_WITHOUT_DRAINGE_STAGE

    id = L2E_FLUX_INFL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_infl = index

    id = L2E_FLUX_TOTDRAIN
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_totdrain = index

    id = L2E_FLUX_GROSS_EVAP_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_gross_evap_soil = index

    id = L2E_FLUX_GROSS_INFL_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_gross_infl_soil = index

    id = L2E_FLUX_SURF
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_surf = index

    id = L2E_FLUX_DEW_GRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_dew_grnd = index

    id = L2E_FLUX_DEW_SNOW
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_dew_snow = index

    id = L2E_FLUX_SUB_SNOW_VOL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_sub_snow_vol = index

    id = L2E_FLUX_SUB_SNOW
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_sub_snow = index

    id = L2E_FLUX_H2OSFC2TOPSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_h2osfc2topsoi = index

    id = L2E_FLUX_SNOW2TOPSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snow2topsoi = index

    id = L2E_FLUX_ROOTSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_rootsoi = index

    id = L2E_FLUX_ADV
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_adv = index

    id = L2E_FLUX_DRAIN_VR
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_drain_vr = index

    id = L2E_FLUX_TRAN_VEG
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_tran_veg = index

    id = L2E_FLUX_ROOTSOI_FRAC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_rootsoi_frac = index


    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List_WaterFlux_Vars

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List_Temperature_Vars(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    use ExternalModelConstants    , only : EM_BETR_STEP_WITHOUT_DRAINGE_STAGE
    use ExternalModelConstants    , only : L2E_STATE_T_SOI10CM
    use ExternalModelConstants    , only : L2E_STATE_T_SOIL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_T_VEG
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_STEP_WITHOUT_DRAINGE_STAGE

    id = L2E_STATE_T_SOI10CM
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_t_soi10cm = index

    id = L2E_STATE_T_SOIL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_t_soil_nlevsoi = index

    id = L2E_STATE_T_VEG
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_t_veg = index


    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List_Temperature_Vars

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List_SoilHydrology_Vars(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    use ExternalModelConstants    , only : EM_BETR_STEP_WITHOUT_DRAINGE_STAGE
    use ExternalModelConstants    , only : L2E_STATE_QCHARGE
    use ExternalModelConstants    , only : L2E_STATE_FRACICE
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_STEP_WITHOUT_DRAINGE_STAGE

    id = L2E_STATE_QCHARGE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_qcharge = index

    id = L2E_STATE_FRACICE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_fracice = index


    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List_SoilHydrology_Vars

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List_CanopyState_Vars(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    use ExternalModelConstants    , only : EM_BETR_STEP_WITHOUT_DRAINGE_STAGE
    use ExternalModelConstants    , only : L2E_STATE_ALTMAX
    use ExternalModelConstants    , only : L2E_STATE_ALTMAX_LASTYEAR
    use ExternalModelConstants    , only : L2E_STATE_LBL_RSC_H2O
    use ExternalModelConstants    , only : L2E_STATE_ELAI
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_STEP_WITHOUT_DRAINGE_STAGE

    id = L2E_STATE_ALTMAX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_altmax = index

    id = L2E_STATE_ALTMAX_LASTYEAR
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_altmax_lastyear = index

    id = L2E_STATE_LBL_RSC_H2O
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_lbl_rsc_h2o = index

    id = L2E_STATE_ELAI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_elai = index


    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List_CanopyState_Vars

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List_Atm2Lnd_Vars(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    use ExternalModelConstants    , only : EM_BETR_STEP_WITHOUT_DRAINGE_STAGE
    use ExternalModelConstants    , only : L2E_STATE_FORC_PBOT_DOWNSCALED
    use ExternalModelConstants    , only : L2E_STATE_FORC_T_DOWNSCALED
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_STEP_WITHOUT_DRAINGE_STAGE

    id = L2E_STATE_FORC_PBOT_DOWNSCALED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_forc_pbot_downscaled = index

    id = L2E_STATE_FORC_T_DOWNSCALED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_forc_t_downscaled = index


    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List_Atm2Lnd_Vars

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List_ChemState_Vars(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    use ExternalModelConstants    , only : EM_BETR_STEP_WITHOUT_DRAINGE_STAGE
    use ExternalModelConstants    , only : L2E_STATE_SOIL_PH
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_STEP_WITHOUT_DRAINGE_STAGE

    id = L2E_STATE_SOIL_PH
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_soil_ph = index


    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List_ChemState_Vars

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List_SoilState_Vars(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    use ExternalModelConstants    , only : EM_BETR_STEP_WITHOUT_DRAINGE_STAGE
    use ExternalModelConstants    , only : L2E_PARAMETER_CELLORG
    use ExternalModelConstants    , only : L2E_PARAMETER_CELLCLAY
    use ExternalModelConstants    , only : L2E_PARAMETER_CELLSAND
    use ExternalModelConstants    , only : L2E_PARAMETER_BD
    use ExternalModelConstants    , only : L2E_PARAMETER_WATFC
    use ExternalModelConstants    , only : L2E_PARAMETER_ROOTFR_PATCH
    use ExternalModelConstants    , only : L2E_PARAMETER_WATSATC
    use ExternalModelConstants    , only : L2E_PARAMETER_BSWC
    use ExternalModelConstants    , only : L2E_PARAMETER_SUCSATC
    use ExternalModelConstants    , only : L2E_PARAMETER_EFFPOROSITYC
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_STEP_WITHOUT_DRAINGE_STAGE

    id = L2E_PARAMETER_CELLORG
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_cellorg = index

    id = L2E_PARAMETER_CELLCLAY
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_cellclay = index

    id = L2E_PARAMETER_CELLSAND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_cellsand = index

    id = L2E_PARAMETER_BD
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_bd = index

    id = L2E_PARAMETER_WATFC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_watfc = index

    id = L2E_PARAMETER_ROOTFR_PATCH
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_rootfr_patch = index

    id = L2E_PARAMETER_WATSATC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_watsatc = index

    id = L2E_PARAMETER_BSWC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_bswc = index

    id = L2E_PARAMETER_SUCSATC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_sucsatc = index

    id = L2E_PARAMETER_EFFPOROSITYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_effporosityc = index


    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List_SoilState_Vars

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_E2L_List(this, e2l_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables returned by BeTR to ALM
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages


  end subroutine EM_BETR_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_BeTR_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)

    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ALMbetrNLMod      , only : betr_namelist_buffer
    use BeTRSimulationALM , only : create_betr_simulation_alm
    use BeTR_decompMod    , only : betr_bounds_type
    use tracer_varcon     , only : is_active_betr_bgc    
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                          :: this
    class(emi_data_list) , intent(in)            :: l2e_init_list
    class(emi_data_list) , intent(inout)         :: e2l_init_list
    integer              , intent(in)            :: iam
    type(bounds_type)    , intent(in)            :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer                :: l2e_max_patch_per_col

    integer  , pointer     :: l2e_init_patch_active      (:)
    integer  , pointer     :: l2e_init_patch_itype       (:)
    real(r8) , pointer     :: l2e_init_patch_wt_col      (:)

    integer  , pointer     :: l2e_init_col_active        (:)
    integer  , pointer     :: l2e_init_col_landunit      (:)
    integer  , pointer     :: l2e_init_col_snl           (:)
    real(r8) , pointer     :: l2e_init_col_zi            (:,:)
    real(r8) , pointer     :: l2e_init_col_dz            (:,:)
    real(r8) , pointer     :: l2e_init_col_z             (:,:)
    integer  , pointer     :: l2e_init_col_pfti          (:)
    integer  , pointer     :: l2e_init_col_pftf          (:)
    integer  , pointer     :: l2e_init_col_npft          (:)

    integer  , pointer     :: l2e_init_landunit_itype    (:)
    real(r8) , pointer     :: l2e_init_finundated_col    (:)
    real(r8) , pointer     :: l2e_init_frac_h2osfc_col   (:)
    real(r8) , pointer     :: l2e_init_h2osoi_liq_col    (:,:)
    real(r8) , pointer     :: l2e_init_h2osoi_ice_col    (:,:)
    real(r8) , pointer     :: l2e_init_h2osoi_liqvol_col (:,:)
    real(r8) , pointer     :: l2e_init_h2osoi_icevol_col (:,:)
    real(r8) , pointer     :: l2e_init_h2osoi_vol_col    (:,:)
    real(r8) , pointer     :: l2e_init_air_vol_col       (:,:)
    real(r8) , pointer     :: l2e_init_rho_vap_col       (:,:)
    real(r8) , pointer     :: l2e_init_rhvap_soi_col     (:,:)
    real(r8) , pointer     :: l2e_init_smp_l_col         (:,:)

    logical                :: masterproc
    logical  , pointer     :: column_active              (:)
    logical  , pointer     :: patch_active               (:)
    integer                :: begc, endc, c, p
    type(betr_bounds_type) :: betr_bounds

    this%betr => create_betr_simulation_alm()

    call l2e_init_list%GetPointerToInt1D  ( this%index_l2e_init_patch_active      , l2e_init_patch_active      )
    call l2e_init_list%GetPointerToInt1D  ( this%index_l2e_init_patch_itype       , l2e_init_patch_itype       )
    call l2e_init_list%GetPointerToReal1D ( this%index_l2e_init_patch_wt_col      , l2e_init_patch_wt_col      )

    call l2e_init_list%GetPointerToInt1D  ( this%index_l2e_init_col_active        , l2e_init_col_active        )
    call l2e_init_list%GetPointerToInt1D  ( this%index_l2e_init_col_landunit      , l2e_init_col_landunit      )
    call l2e_init_list%GetPointerToInt1D  ( this%index_l2e_init_col_snl           , l2e_init_col_snl           )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_col_zi            , l2e_init_col_zi            )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_col_dz            , l2e_init_col_dz            )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_col_z             , l2e_init_col_z             )
    call l2e_init_list%GetPointerToInt1D  ( this%index_l2e_init_col_pfti          , l2e_init_col_pfti          )
    call l2e_init_list%GetPointerToInt1D  ( this%index_l2e_init_col_pftf          , l2e_init_col_pftf          )
    call l2e_init_list%GetPointerToInt1D  ( this%index_l2e_init_col_npft          , l2e_init_col_npft          )
    
    call l2e_init_list%GetPointerToInt1D  ( this%index_l2e_init_landunit_itype    , l2e_init_landunit_itype    )

    call l2e_init_list%GetPointerToReal1D ( this%index_l2e_init_finundated_col    , l2e_init_finundated_col    )
    call l2e_init_list%GetPointerToReal1D ( this%index_l2e_init_frac_h2osfc_col   , l2e_init_frac_h2osfc_col   )

    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_h2osoi_liq_col    , l2e_init_h2osoi_liq_col    )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_h2osoi_ice_col    , l2e_init_h2osoi_ice_col    )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_h2osoi_liqvol_col , l2e_init_h2osoi_liqvol_col )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_h2osoi_icevol_col , l2e_init_h2osoi_icevol_col )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_h2osoi_vol_col    , l2e_init_h2osoi_vol_col    )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_air_vol_col       , l2e_init_air_vol_col       )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_rho_vap_col       , l2e_init_rho_vap_col       )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_rhvap_soi_col     , l2e_init_rhvap_soi_col     )
    call l2e_init_list%GetPointerToReal2D ( this%index_l2e_init_smp_l_col         , l2e_init_smp_l_col         )

    call l2e_init_list%GetIntValue        ( this%index_l2e_init_max_patch_per_col , l2e_max_patch_per_col      )

    begc = bounds_clump%begc
    endc = bounds_clump%endc

    if (iam==0) then
       masterproc = .true.
    else
       masterproc = .false.
    end if

    allocate(patch_active(bounds_clump%begp:bounds_clump%endp))

    do p = bounds_clump%begp, bounds_clump%endp
       if (l2e_init_patch_active(p) == 1) then
          patch_active(p) = .true.
       else
          patch_active(p) = .false.
       end if
    end do

    allocate(column_active(bounds_clump%begc:bounds_clump%endc))
    do c = begc, endc
       if (l2e_init_col_active(c) == 1) then
          column_active(c) = .true.
       else
          column_active(c) = .false.
       end if
    end do

    call this%betr%BeTRSetFilter(maxpft_per_col=l2e_max_patch_per_col, boffline=.false.)
    
    call this%betr%betr_time%Init(betr_namelist_buffer, masterproc)

    call this%betr%SetBaseFilename()

    call this%betr%betr_time%Init(betr_namelist_buffer, masterproc)

    !allocate memory
    call this%betr%BeTRInitAllocateMemory(begc, endc)

    !grid horizontal bounds
    call this%betr%BeTRSetBounds(betr_bounds)

    call this%betr%BeTRInitializeBiophysForc (betr_bounds, begc, endc)

    call this%betr%BeTRInitializeCols        (betr_bounds, begc, endc)

    call this%betr%BeTRInitializePFT         (betr_bounds, begc, endc)

    call this%betr%BeTRInitializeActiveCols  (betr_bounds, begc, endc, &
         l2e_init_col_landunit, l2e_init_landunit_itype)

    call this%betr%BeTRSetcpsCol(begc, endc,                                 &
         l2e_init_col_snl, l2e_init_col_zi, l2e_init_col_dz, l2e_init_col_z, &
         l2e_init_col_pfti, l2e_init_col_pftf, l2e_init_col_npft)

    call this%betr%BeTRSetcpsPFT(begc, endc,                                 &
         l2e_init_col_pfti, l2e_init_col_npft, patch_active,                 &
         l2e_init_patch_itype, l2e_init_patch_wt_col)

    call this%betr%BeTRSetBiophysForcingWaterStateVars( &
         begc, endc,                                    &
         betr_bounds%lbj, betr_bounds%ubj,              &
         l2e_init_finundated_col,                       &
         l2e_init_frac_h2osfc_col,                      &
         l2e_init_h2osoi_liq_col,                       &
         l2e_init_h2osoi_ice_col,                       &
         l2e_init_h2osoi_liqvol_col,                    &
         l2e_init_h2osoi_icevol_col,                    &
         l2e_init_h2osoi_vol_col,                       &
         l2e_init_air_vol_col,                          &
         l2e_init_rho_vap_col,                          &
         l2e_init_rhvap_soi_col,                        &
         l2e_init_smp_l_col)

    call this%betr%BeTRInitialize_1(betr_bounds, begc, &
         endc, betr_namelist_buffer)

    !call this%betr%BeTRCreateHistory(bounds_clump)

    !call this%betr%RestAlloc(begc, endc)

    call this%betr%set_active(begc, endc, column_active)

    is_active_betr_bgc = this%betr%do_soibgc()

    deallocate(patch_active)
    deallocate(column_active)

  end subroutine EM_BeTR_Init

    !------------------------------------------------------------------------
  subroutine EM_BETR_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    ! 
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use ExternalModelConstants    , only : EM_BETR_BEGIN_MASS_BALANCE_STAGE
    use ExternalModelConstants    , only : EM_BETR_PRE_DIAG_WATER_FLUX_STAGE
    use ExternalModelConstants    , only : EM_BETR_PRE_DIAG_DTRACER_FREEZE_THAW_STAGE
    use clm_varctl                , only : iulog
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                  :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent (in)   :: bounds_clump

    select case(em_stage)

    case (EM_BETR_BEGIN_MASS_BALANCE_STAGE)
       call EM_BETR_BeginMassBalance_Solve(this, dt, nstep, bounds_clump, l2e_list, e2l_list)

    case (EM_BETR_PRE_DIAG_WATER_FLUX_STAGE)
       call EM_BETR_PreDiagSoilColWaterFlux_Solve(this, bounds_clump, l2e_list, e2l_list)

    case (EM_BETR_PRE_DIAG_DTRACER_FREEZE_THAW_STAGE)
       call EM_BETR_PreDiagDtracerFreezeThaw_Solve(this, bounds_clump, l2e_list, e2l_list)

    case default
       write(iulog,*)'EM_BETR_Solve: Unknown em_stage.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine EM_BETR_Solve

    !------------------------------------------------------------------------
  subroutine EM_BETR_BeginMassBalance_Solve(this, dt, nstep, bounds_clump, l2e_list, e2l_list)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use ExternalModelConstants    , only : EM_BETR_BEGIN_MASS_BALANCE_STAGE
    use ExternalModelConstants    , only : EM_BETR_PRE_DIAG_WATER_FLUX_STAGE
    use clm_varctl                , only : iulog
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                  :: this
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    type(bounds_type)    , intent(in)    :: bounds_clump
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list

    call this%betr%SetClock(dtime = dt, nelapstep = nstep)

    call this%betr%BeginMassBalanceCheck(bounds_clump)


  end subroutine EM_BETR_BeginMassBalance_Solve

    !------------------------------------------------------------------------
  subroutine EM_BETR_PreDiagSoilColWaterFlux_Solve(this, bounds, l2e_list, e2l_list)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use ExternalModelConstants    , only : EM_BETR_BEGIN_MASS_BALANCE_STAGE
    use ExternalModelConstants    , only : EM_BETR_PRE_DIAG_WATER_FLUX_STAGE
    use clm_varctl                , only : iulog
    use clm_varpar                , only : nlevsoi
    use BeTR_decompMod            , only : betr_bounds_type
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                  :: this
    type(bounds_type)    , intent(in)    :: bounds
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer      :: l2e_frac_h2osfc_col(:)
    real(r8), pointer      :: l2e_finundated_col(:)
    real(r8), pointer      :: l2e_h2osoi_liq_col(:,:)
    real(r8), pointer      :: l2e_h2osoi_ice_col(:,:)
    real(r8), pointer      :: l2e_h2osoi_liqvol_col(:,:)
    real(r8), pointer      :: l2e_h2osoi_icevol_col(:,:)
    real(r8), pointer      :: l2e_h2osoi_vol_col(:,:)
    real(r8), pointer      :: l2e_air_vol_col(:,:)
    real(r8), pointer      :: l2e_rho_vap_col(:,:)
    real(r8), pointer      :: l2e_rhvap_soi_col(:,:)
    real(r8), pointer      :: l2e_smp_l_col(:,:)
    integer, pointer       :: l2e_filter_nolakec(:)
    integer                :: l2e_num_nolakec
    integer                :: cc, c, fc, lbj, ubj

    type(betr_bounds_type) :: betr_bounds

    call l2e_list%GetPointerToReal1D(this%index_l2e_state_frac_h2osfc  , l2e_frac_h2osfc_col    )
    call l2e_list%GetPointerToReal1D(this%index_l2e_state_finundated   , l2e_finundated_col     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq   , l2e_h2osoi_liq_col     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice   , l2e_h2osoi_ice_col     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liqvol, l2e_h2osoi_liqvol_col  )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_icevol, l2e_h2osoi_icevol_col  )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_vol   , l2e_h2osoi_vol_col     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_air_vol      , l2e_air_vol_col        )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_rho_vap      , l2e_rho_vap_col        )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_rhvap_soi    , l2e_rhvap_soi_col      )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_smp_l        , l2e_smp_l_col          )

    call l2e_list%GetPointerToInt1D(this%index_l2e_filter_nolakec      , l2e_filter_nolakec )
    call l2e_list%GetIntValue(this%index_l2e_filter_num_nolakec        , l2e_num_nolakec    )

    call this%betr%BeTRSetBounds(betr_bounds)

    call this%betr%BeTRSetBiophysForcingWaterstateVars( &
         bounds%begc, bounds%endc,                      &
         betr_bounds%lbj, betr_bounds%ubj,              &
         l2e_finundated_col,                            &
         l2e_frac_h2osfc_col,                           &
         l2e_h2osoi_liq_col,                            &
         l2e_h2osoi_ice_col,                            &
         l2e_h2osoi_liqvol_col,                         &
         l2e_h2osoi_icevol_col,                         &
         l2e_h2osoi_vol_col,                            &
         l2e_air_vol_col,                               &
         l2e_rho_vap_col,                               &
         l2e_rhvap_soi_col,                             &
         l2e_smp_l_col)


    call this%betr%PreDiagSoilColWaterFlux(l2e_num_nolakec , l2e_filter_nolakec)

  end subroutine EM_BETR_PreDiagSoilColWaterFlux_Solve

    !------------------------------------------------------------------------
  subroutine EM_BETR_PreDiagDtracerFreezeThaw_Solve(this, bounds, l2e_list, e2l_list)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use ExternalModelConstants    , only : EM_BETR_BEGIN_MASS_BALANCE_STAGE
    use ExternalModelConstants    , only : EM_BETR_PRE_DIAG_WATER_FLUX_STAGE
    use clm_varctl                , only : iulog
    use clm_varpar                , only : nlevsoi
    use BeTR_decompMod            , only : betr_bounds_type
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                  :: this
    type(bounds_type)    , intent(in)    :: bounds
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    integer  , pointer     :: l2e_col_snl           (:)
    real(r8) , pointer     :: l2e_col_zi            (:,:)
    real(r8) , pointer     :: l2e_col_dz            (:,:)
    real(r8) , pointer     :: l2e_col_z             (:,:)
    integer  , pointer     :: l2e_col_pfti          (:)
    integer  , pointer     :: l2e_col_pftf          (:)
    integer  , pointer     :: l2e_col_npft          (:)
    !
    real(r8) , pointer     :: l2e_frac_h2osfc_col   (:)
    real(r8) , pointer     :: l2e_finundated_col    (:)
    real(r8) , pointer     :: l2e_h2osoi_liq_col    (:,:)
    real(r8) , pointer     :: l2e_h2osoi_ice_col    (:,:)
    real(r8) , pointer     :: l2e_h2osoi_liqvol_col (:,:)
    real(r8) , pointer     :: l2e_h2osoi_icevol_col (:,:)
    real(r8) , pointer     :: l2e_h2osoi_vol_col    (:,:)
    real(r8) , pointer     :: l2e_air_vol_col       (:,:)
    real(r8) , pointer     :: l2e_rho_vap_col       (:,:)
    real(r8) , pointer     :: l2e_rhvap_soi_col     (:,:)
    real(r8) , pointer     :: l2e_smp_l_col         (:,:)
    integer  , pointer     :: l2e_filter_nolakec    (:)
    !
    integer                :: l2e_num_nolakec

    type(betr_bounds_type) :: betr_bounds

    call l2e_list%GetPointerToInt1D  ( this%index_l2e_col_snl           , l2e_col_snl           )
    call l2e_list%GetPointerToReal2D ( this%index_l2e_col_zi            , l2e_col_zi            )
    call l2e_list%GetPointerToReal2D ( this%index_l2e_col_dz            , l2e_col_dz            )
    call l2e_list%GetPointerToReal2D ( this%index_l2e_col_z             , l2e_col_z             )
    call l2e_list%GetPointerToInt1D  ( this%index_l2e_col_pfti          , l2e_col_pfti          )
    call l2e_list%GetPointerToInt1D  ( this%index_l2e_col_pftf          , l2e_col_pftf          )
    call l2e_list%GetPointerToInt1D  ( this%index_l2e_col_npft          , l2e_col_npft          )

    call l2e_list%GetPointerToReal1D(this%index_l2e_state_frac_h2osfc  , l2e_frac_h2osfc_col    )
    call l2e_list%GetPointerToReal1D(this%index_l2e_state_finundated   , l2e_finundated_col     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq   , l2e_h2osoi_liq_col     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice   , l2e_h2osoi_ice_col     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liqvol, l2e_h2osoi_liqvol_col  )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_icevol, l2e_h2osoi_icevol_col  )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_vol   , l2e_h2osoi_vol_col     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_air_vol      , l2e_air_vol_col        )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_rho_vap      , l2e_rho_vap_col        )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_rhvap_soi    , l2e_rhvap_soi_col      )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_smp_l        , l2e_smp_l_col          )

    call l2e_list%GetPointerToInt1D(this%index_l2e_filter_nolakec      , l2e_filter_nolakec )
    call l2e_list%GetIntValue(this%index_l2e_filter_num_nolakec        , l2e_num_nolakec    )

    call this%betr%BeTRSetBounds(betr_bounds)

    call this%betr%BeTRSetBiophysForcingWaterstateVars( &
         bounds%begc, bounds%endc,                      &
         betr_bounds%lbj, betr_bounds%ubj,              &
         l2e_finundated_col,                            &
         l2e_frac_h2osfc_col,                           &
         l2e_h2osoi_liq_col,                            &
         l2e_h2osoi_ice_col,                            &
         l2e_h2osoi_liqvol_col,                         &
         l2e_h2osoi_icevol_col,                         &
         l2e_h2osoi_vol_col,                            &
         l2e_air_vol_col,                               &
         l2e_rho_vap_col,                               &
         l2e_rhvap_soi_col,                             &
         l2e_smp_l_col)

    call this%betr%BeTRSetcpsCol(bounds%begc, bounds%endc, &
         l2e_col_snl, l2e_col_zi, l2e_col_dz, l2e_col_z, &
         l2e_col_pfti, l2e_col_pftf, l2e_col_npft)

    call this%betr%DiagnoseDtracerFreezeThaw(bounds, l2e_num_nolakec , l2e_filter_nolakec)

  end subroutine EM_BETR_PreDiagDtracerFreezeThaw_Solve

end module ExternalModelBETRMod
