module ExternalModelVSFMSPACMod

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides
  !
  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use ExternalModelInterfaceDataMod, only : emi_data_list, emi_data
  use mpp_varctl                   , only : iulog
  use ExternalModelBaseType        , only : em_base_type
  use ExternalModelVSFMMod         , only : em_vsfm_type
  use MultiPhysicsProbVSFM         , only : mpp_vsfm_type
  use decompMod                    , only : bounds_type
  use petscsys
  !
  implicit none
  !

  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x_column = 1.d0
  PetscReal , parameter :: y_column = 1.d0
  PetscReal , parameter :: z_column = 1.d0
  PetscInt              :: nz
  PetscInt              :: ncells_ghost

  PetscInt, parameter  :: nx_soil        =   1
  PetscInt, parameter  :: ny_soil        =   1
  PetscInt, parameter  :: nz_soil        =  50
  PetscInt, parameter  :: nz_soil_top    =   3
  PetscInt, parameter  :: nz_soil_mid    =   3
  PetscInt, parameter  :: nz_soil_bot    =  44

  PetscInt, parameter  :: nx_root        =   1
  PetscInt, parameter  :: ny_root        =   1
  PetscInt, parameter  :: nz_root        =  30

  PetscInt, parameter  :: nx_xylem       =   1
  PetscInt, parameter  :: ny_xylem       =   1
  PetscInt, parameter  :: nz_xylem       = 170

  PetscReal, parameter  :: dx            = 1.d0         ! [m]
  PetscReal, parameter  :: dy            = 1.d0         ! [m]
  PetscReal, parameter  :: dz            = 0.1d0        ! [m]

  PetscReal, parameter  :: dx_xylem      = 0.25d0       ! [m]
  PetscReal, parameter  :: dy_xylem      = 0.25d0       ! [m]
  PetscReal, parameter  :: xylem_height  = 17.d0        ! [m]
  PetscReal, parameter  :: xylem_LAI     = 5.d0         ! [leaf area m^2/ soil area m^2]

  PetscReal, parameter  :: scalez        = 0.025d0

  ! Parameters for root length density = length-of-root/volume-of-soil  [m_root/m^3_soil]
  PetscReal, parameter :: qz             = 9.d0          ! [-]
  PetscReal, parameter :: d              = 3.2d0         ! [m]
  PetscReal, parameter :: rld0           = 4.d4          ! [m/m^3]
  PetscReal, parameter :: root_radius    = 2.d-3         ! [m]

  PetscReal, parameter :: pi              = 3.1416d0
  PetscBool :: single_pde_formulation

  PetscReal, parameter :: perm_xy_top    = 6.83d-11      ! [m^2]
  PetscReal, parameter :: perm_z_top     = 6.83d-11      ! [m^2]
  PetscReal, parameter :: sat_res_top    = 0.06d0        ! [-]
  PetscReal, parameter :: alpha_top      = 0.00005d0     ! [Pa^{-1}]
  PetscReal, parameter :: vg_m_top       = 0.33d0        ! [-]
  PetscReal, parameter :: por_top        = 0.5d0         ! [-]

  PetscReal, parameter :: perm_xy_mid    = 6.83d-11      ! [m^2]
  PetscReal, parameter :: perm_z_mid     = 6.83d-11      ! [m^2]
  PetscReal, parameter :: sat_res_mid    = 0.06d0        ! [-]
  PetscReal, parameter :: alpha_mid      = 0.00005d0     ! [Pa^{-1}]
  PetscReal, parameter :: vg_m_mid       = 0.33d0        ! [-]
  PetscReal, parameter :: por_mid        = 0.5d0         ! [-]

  PetscReal, parameter :: perm_xy_bot    = 6.83d-11      ! [m^2]
  PetscReal, parameter :: perm_z_bot     = 6.83d-11      ! [m^2]
  PetscReal, parameter :: sat_res_bot    = 0.06d0        ! [-]
  PetscReal, parameter :: alpha_bot      = 0.00005d0     ! [Pa^{-1}]
  PetscReal, parameter :: vg_m_bot       = 0.33d0        ! [-]
  PetscReal, parameter :: por_bot        = 0.5d0         ! [-]

  PetscReal, parameter :: perm_root       = 6.83d-11      ! [m^2]
  PetscReal, parameter :: sat_res_root    = 0.06d0        ! [-]
  PetscReal, parameter :: alpha_root      = 0.00005d0     ! [Pa^{-1}]
  PetscReal, parameter :: vg_m_root       = 0.33d0        ! [-]
  PetscReal, parameter :: por_root        = 0.5d0         ! [-]

  PetscReal, parameter :: perm_xylem      = 6.83d-11      ! [m^2]
  PetscReal, parameter :: sat_res_xylem   = 0.06d0        ! [-]
  PetscReal, parameter :: alpha_xylem     = 0.00005d0     ! [Pa^{-1}]
  PetscReal, parameter :: vg_m_xylem      = 0.33d0        ! [-]
  PetscReal, parameter :: por_xylem       = 0.5d0         ! [-]

  PetscReal, parameter :: press_initial  = 3.5355d3      ! [Pa]

  type, public, extends(em_vsfm_type) :: em_vsfm_spac_type

     integer :: index_e2l_init_state_h2oroot_liq
     integer :: index_e2l_init_state_h2oxylem_liq

     integer :: index_e2l_state_h2oroot_liq
     integer :: index_e2l_state_h2oxylem_liq
     integer :: index_e2l_state_xylemp

   contains
     !procedure, public :: Populate_L2E_Init_List  => EM_VSFM_SPAC_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EM_VSFM_SPAC_Populate_E2L_Init_List
     !procedure, public :: Populate_L2E_List       => EM_VSFM_SPAC_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_VSFM_SPAC_Populate_E2L_List
     procedure, public :: Init                    => EM_VSFM_SPAC_Init
     procedure, public :: Solve                   => EM_VSFM_SPAC_Solve
  end type em_vsfm_spac_type

contains

  !------------------------------------------------------------------------
  subroutine EM_VSFM_SPAC_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)

    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use mpp_varctl                , only : vsfm_use_dynamic_linesearch
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use SystemOfEquationsVSFMType , only : VSFMSOEUpdateConnections
    use petscsnes
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_spac_type)             :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent(in)    :: bounds_clump

    single_pde_formulation = PETSC_FALSE

    !
    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp(this%vsfm_mpp, iam)

    ! 2. Add all meshes needed for the MPP
    call add_meshes(this%vsfm_mpp)

    ! 3. Add all governing equations
    call add_goveqns(this%vsfm_mpp)

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns(this%vsfm_mpp)

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars(this%vsfm_mpp)

    ! 6. Setup the MPP
    !call vsfm_mpp%SetupProblem(vsfm_use_dynamic_linesearch)
    call setup_petsc_snes(this%vsfm_mpp)

    ! 7. Add material properities associated with all governing equations
    call set_material_properties(this%vsfm_mpp)

    ! 8. Set initial conditions
    call set_initial_conditions(this%vsfm_mpp)

    ! 9. Determine IDs for various source-sink condition
    call determine_condition_ids(this)

    !10.
    call VSFMSOEUpdateConnections(this%vsfm_mpp%sysofeqns, MPP_VSFM_SNES_CLM)

    !11.
    call extract_data_for_alm(this, l2e_init_list, e2l_init_list, bounds_clump)

  end subroutine EM_VSFM_SPAC_Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp(vsfm_mpp, iam)
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    !
    implicit none
    !
    ! !ARGUMENTS
    type(mpp_vsfm_type) :: vsfm_mpp
    integer, intent(in) :: iam
    !
    ! Set up the multi-physics problem
    !
    call vsfm_mpp%Init       ()
    call vsfm_mpp%SetName    ('Variably-Saturated-Flow-Model')
    call vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
    call vsfm_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes(vsfm_mpp)
    !
    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_AGAINST_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants , only : MESH_SPAC_XYLEM_COL
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : VAR_VOLUME
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    !
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    !
    PetscInt            :: imesh, ii, jj, kk
    PetscInt            :: ncells_soil, ncells_root, ncells_xylem
    PetscInt            :: ncells_srx
    PetscInt            :: count
    PetscInt            :: iconn, nconn

    PetscReal , pointer :: soil_xc(:)           ! x-position of grid cell [m]
    PetscReal , pointer :: soil_yc(:)           ! y-position of grid cell [m]
    PetscReal , pointer :: soil_zc(:)           ! z-position of grid cell [m]
    PetscReal , pointer :: soil_dx(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: soil_dy(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: soil_area(:)         ! area of grid cell [m^2]
    PetscInt  , pointer :: soil_filter(:)       ! 
    PetscReal , pointer :: soil_vol(:)

    PetscReal , pointer :: root_xc(:)           ! x-position of grid cell [m]
    PetscReal , pointer :: root_yc(:)           ! y-position of grid cell [m]
    PetscReal , pointer :: root_zc(:)           ! z-position of grid cell [m]
    PetscReal , pointer :: root_dx(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: root_dy(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: root_dz(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: root_area(:)         ! area of grid cell [m^2]
    PetscInt  , pointer :: root_filter(:)       ! 
    PetscReal , pointer :: root_vol(:)

    PetscReal , pointer :: xylem_xc(:)           ! x-position of grid cell [m]
    PetscReal , pointer :: xylem_yc(:)           ! y-position of grid cell [m]
    PetscReal , pointer :: xylem_zc(:)           ! z-position of grid cell [m]
    PetscReal , pointer :: xylem_dx(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: xylem_dy(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: xylem_dz(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: xylem_area(:)         ! area of grid cell [m^2]
    PetscInt  , pointer :: xylem_filter(:)       ! 
    PetscReal , pointer :: xylem_vol(:)

    PetscReal , pointer :: srx_xc(:)           ! x-position of grid cell [m]
    PetscReal , pointer :: srx_yc(:)           ! y-position of grid cell [m]
    PetscReal , pointer :: srx_zc(:)           ! z-position of grid cell [m]
    PetscReal , pointer :: srx_dx(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: srx_dy(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: srx_dz(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: srx_area(:)         ! area of grid cell [m^2]
    PetscInt  , pointer :: srx_filter(:)       ! 
    PetscReal , pointer :: srx_vol(:)

    PetscReal , pointer :: root_len_den(:)
    PetscReal , pointer :: root_surf_area(:)
    PetscReal           :: root_len

    PetscInt  , pointer :: soil_conn_id_up(:)   !
    PetscInt  , pointer :: soil_conn_id_dn(:)   !
    PetscReal , pointer :: soil_conn_dist_up(:) !
    PetscReal , pointer :: soil_conn_dist_dn(:) !
    PetscReal , pointer :: soil_conn_area(:)    !
    PetscInt  , pointer :: soil_conn_type(:)    !

    PetscInt  , pointer :: root_conn_id_up(:)   !
    PetscInt  , pointer :: root_conn_id_dn(:)   !
    PetscReal , pointer :: root_conn_dist_up(:) !
    PetscReal , pointer :: root_conn_dist_dn(:) !
    PetscReal , pointer :: root_conn_area(:)    !
    PetscInt  , pointer :: root_conn_type(:)    !

    PetscInt  , pointer :: xylem_conn_id_up(:)   !
    PetscInt  , pointer :: xylem_conn_id_dn(:)   !
    PetscReal , pointer :: xylem_conn_dist_up(:) !
    PetscReal , pointer :: xylem_conn_dist_dn(:) !
    PetscReal , pointer :: xylem_conn_area(:)    !
    PetscInt  , pointer :: xylem_conn_type(:)    !

    PetscInt  , pointer :: srx_conn_id_up(:)   !
    PetscInt  , pointer :: srx_conn_id_dn(:)   !
    PetscReal , pointer :: srx_conn_dist_up(:) !
    PetscReal , pointer :: srx_conn_dist_dn(:) !
    PetscReal , pointer :: srx_conn_area(:)    !
    PetscInt  , pointer :: srx_conn_type(:)    !

    PetscReal , pointer :: zv_x(:)
    PetscReal , pointer :: zv_y(:)
    PetscReal , pointer :: xv_2d(:,:)
    PetscReal , pointer :: yv_2d(:,:)
    PetscReal , pointer :: zv_2d(:,:)
    PetscReal , pointer :: xc_3d(:,:,:)
    PetscReal , pointer :: yc_3d(:,:,:)
    PetscReal , pointer :: zc_3d(:,:,:)
    PetscInt  , pointer :: id_3d(:,:,:)
    PetscReal           :: dist_x, dist_y, dist_z, dist

    PetscErrorCode :: ierr

    call mpp_varpar_set_nlevsoi(nz_soil)
    call mpp_varpar_set_nlevgrnd(nz_soil)

    !
    ! Set up the meshes
    !    
    call vsfm_mpp%SetNumMeshes(3)

    ncells_ghost = 0
    ncells_soil  = nx_soil *ny_soil *nz_soil
    ncells_root  = nx_root *ny_root *nz_root
    ncells_xylem = nx_xylem*ny_xylem*nz_xylem
    ncells_srx   = ncells_soil + ncells_root + ncells_xylem

    allocate(soil_xc            (ncells_soil))
    allocate(soil_yc            (ncells_soil))
    allocate(soil_zc            (ncells_soil))
    allocate(soil_dx            (ncells_soil))
    allocate(soil_dy            (ncells_soil))
    allocate(soil_dz            (ncells_soil))
    allocate(soil_area          (ncells_soil))
    allocate(soil_filter        (ncells_soil))
    allocate(soil_vol           (ncells_soil))

    allocate(root_xc            (ncells_root))
    allocate(root_yc            (ncells_root))
    allocate(root_zc            (ncells_root))
    allocate(root_dx            (ncells_root))
    allocate(root_dy            (ncells_root))
    allocate(root_dz            (ncells_root))
    allocate(root_area          (ncells_root))
    allocate(root_filter        (ncells_root))
    allocate(root_vol           (ncells_root))
    allocate(root_len_den       (ncells_root))
    allocate(root_surf_area     (ncells_root))

    allocate(xylem_xc            (ncells_xylem))
    allocate(xylem_yc            (ncells_xylem))
    allocate(xylem_zc            (ncells_xylem))
    allocate(xylem_dx            (ncells_xylem))
    allocate(xylem_dy            (ncells_xylem))
    allocate(xylem_dz            (ncells_xylem))
    allocate(xylem_area          (ncells_xylem))
    allocate(xylem_filter        (ncells_xylem))
    allocate(xylem_vol           (ncells_xylem))

    allocate(srx_xc            (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_yc            (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_zc            (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_dx            (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_dy            (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_dz            (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_area          (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_filter        (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_vol           (ncells_soil + ncells_root + ncells_xylem))

    allocate(soil_conn_id_up    (ncells_soil-1))
    allocate(soil_conn_id_dn    (ncells_soil-1))
    allocate(soil_conn_dist_up  (ncells_soil-1))
    allocate(soil_conn_dist_dn  (ncells_soil-1))
    allocate(soil_conn_area     (ncells_soil-1))
    allocate(soil_conn_type     (ncells_soil-1))

    allocate(root_conn_id_up    (ncells_root-1))
    allocate(root_conn_id_dn    (ncells_root-1))
    allocate(root_conn_dist_up  (ncells_root-1))
    allocate(root_conn_dist_dn  (ncells_root-1))
    allocate(root_conn_area     (ncells_root-1))
    allocate(root_conn_type     (ncells_root-1))

    allocate(xylem_conn_id_up   (ncells_xylem-1))
    allocate(xylem_conn_id_dn   (ncells_xylem-1))
    allocate(xylem_conn_dist_up (ncells_xylem-1))
    allocate(xylem_conn_dist_dn (ncells_xylem-1))
    allocate(xylem_conn_area    (ncells_xylem-1))
    allocate(xylem_conn_type    (ncells_xylem-1))


    allocate(srx_conn_id_up   (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    allocate(srx_conn_id_dn   (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    allocate(srx_conn_dist_up (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    allocate(srx_conn_dist_dn (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    allocate(srx_conn_area    (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    allocate(srx_conn_type    (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    
    ! Soil mesh
    soil_filter (:) = 1
    soil_area   (:) = dx*dy
    soil_dx     (:) = dx
    soil_dy     (:) = dy
    soil_dz     (:) = dz
    soil_xc     (:) = dx/2.d0
    soil_yc     (:) = dy/2.d0
    soil_zc     (:) = 0.d0
    soil_zc     (1) = -dz/2.d0
    do kk = 2,nz_soil
       soil_zc(kk) = soil_zc(kk-1) - dz
    enddo

    root_filter (:) = 1
    root_area   (:) = pi*(root_radius**2.d0)
    root_dx     (:) = dx
    root_dy     (:) = dy
    root_dz     (:) = dz
    root_xc     (:) = dx/2.d0 + root_radius
    root_yc     (:) = dy/2.d0
    root_zc     (:) = 0.d0
    do kk = 1,nz_root
       root_zc(kk) = soil_zc(kk)
    enddo

    xylem_filter (:) = 1
    xylem_area   (:) = dx_xylem*dy_xylem
    xylem_dx     (:) = dx_xylem
    xylem_dy     (:) = dy_xylem
    xylem_dz     (:) = dz
    xylem_xc     (:) = dx/2.d0 + root_radius
    xylem_yc     (:) = dy/2.d0
    xylem_zc     (1) = nz_xylem*dz - dz/2.d0
    kk = 1;
    do kk = 2, nz_xylem
       xylem_zc(kk) = xylem_zc(kk-1) - dz
    enddo

    do kk = 1, nz_soil

       if (kk <= nz_root) then
          root_len_den(kk) = rld0 * (1.d0 - abs(root_zc(kk))/d)*exp(-qz*abs(root_zc(kk))/d)

          root_len = root_len_den(kk)*     &       ! [m/m^3]
               (dx*dy*soil_dz(kk))           ! [m^3]

          root_surf_area(kk) = 2.d0*pi*root_radius*root_len   ! [m^2]

          root_vol(kk) = pi*(root_radius**2.d0)*root_len      ! [m^3]

          soil_vol(kk) = dx * dy * soil_dz(kk) - root_vol(kk) ! [m^3]
       else
          soil_vol(kk) = dx * dy * soil_dz(kk)
       endif

    end do
    xylem_vol(:) = dx*dy*dz
    soil_vol(:)  = dx*dy*dz
    root_vol(:)  = dx*dy*dz
    root_area(:) = dx*dy
    xylem_area(:)= dx*dy

    srx_filter (:) = 1
    srx_area   (:) = dx*dy
    srx_dx     (:) = dx
    srx_dy     (:) = dy
    srx_dz     (:) = dz
    srx_vol    (:) = dx*dy*dz

    srx_xc     (1                        :ncells_soil                         ) = soil_xc(1:ncells_soil )
    srx_xc     (ncells_soil+1            :ncells_soil+ncells_root             ) = root_xc(1:ncells_root )
    srx_xc     (ncells_soil+ncells_root+1:ncells_soil+ncells_root+ncells_xylem) = xylem_xc(1:ncells_xylem)

    srx_yc     (1                        :ncells_soil                         ) = soil_yc(1:ncells_soil )
    srx_yc     (ncells_soil+1            :ncells_soil+ncells_root             ) = root_yc(1:ncells_root )
    srx_yc     (ncells_soil+ncells_root+1:ncells_soil+ncells_root+ncells_xylem) = xylem_yc(1:ncells_xylem)

    srx_zc     (1                        :ncells_soil                         ) = soil_zc(1:ncells_soil )
    srx_zc     (ncells_soil+1            :ncells_soil+ncells_root             ) = root_zc(1:ncells_root )
    srx_zc     (ncells_soil+ncells_root+1:ncells_soil+ncells_root+ncells_xylem) = xylem_zc(1:ncells_xylem)

    ! *********** Soil mesh ***********

    if (.not. single_pde_formulation) then

       imesh        = 1

       call vsfm_mpp%MeshSetName                (imesh, 'Soil mesh')
       call vsfm_mpp%MeshSetOrientation         (imesh, MESH_ALONG_GRAVITY)
       call vsfm_mpp%MeshSetID                  (imesh, MESH_CLM_SOIL_COL)
       call vsfm_mpp%MeshSetDimensions          (imesh, ncells_soil, ncells_ghost, nz_soil)

       call vsfm_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME  , soil_vol)
       
       !
       ! Vertical connections
       !
       iconn = 0
       do kk = 1, nz_soil-1

          iconn               = iconn + 1
          soil_conn_id_up(iconn)   = kk
          soil_conn_id_dn(iconn)   = kk + 1
          soil_conn_dist_up(iconn) = 0.5d0*soil_dz(kk  )
          soil_conn_dist_dn(iconn) = 0.5d0*soil_dz(kk+1)
          soil_conn_area(iconn)    = soil_area(kk)
          soil_conn_type(iconn)    = CONN_VERTICAL
       end do
       nconn = iconn

       call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
            nconn,  soil_conn_id_up, soil_conn_id_dn,               &
            soil_conn_dist_up, soil_conn_dist_dn, soil_conn_area,   &
            soil_conn_type)

       ! *********** Root mesh *********** 

       imesh        = 2

       ncells_ghost = 0

       call vsfm_mpp%MeshSetName                (imesh, 'Root mesh')
       call vsfm_mpp%MeshSetOrientation         (imesh, MESH_ALONG_GRAVITY)
       call vsfm_mpp%MeshSetID                  (imesh, MESH_SPAC_ROOT_COL)
       call vsfm_mpp%MeshSetDimensions          (imesh, ncells_root, ncells_ghost, nz_root)
       call vsfm_mpp%MeshSetGridCellFilter      (imesh, root_filter)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , root_xc)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , root_yc)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , root_zc)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , root_dx)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , root_dy)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , root_dz)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , root_area)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME  , root_vol)

       !
       ! Vertical connections
       !
       iconn = 0
       do kk = 1, nz_root-1
          iconn                    = iconn + 1
          root_conn_id_up(iconn)   = kk
          root_conn_id_dn(iconn)   = kk + 1
          root_conn_dist_up(iconn) = 0.5d0*root_dz(kk  )
          root_conn_dist_dn(iconn) = 0.5d0*root_dz(kk+1)
          root_conn_area(iconn)    = root_area(kk)
          root_conn_type(iconn)    = CONN_VERTICAL
       end do
       nconn = iconn

       call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
            nconn,  root_conn_id_up, root_conn_id_dn,               &
            root_conn_dist_up, root_conn_dist_dn, root_conn_area,   &
            root_conn_type)

       ! *********** Xylem mesh *********** 
       imesh        = 3

       ncells_ghost = 0

       call vsfm_mpp%MeshSetName                (imesh , 'Xylem mesh'                         )
       call vsfm_mpp%MeshSetOrientation         (imesh , MESH_ALONG_GRAVITY                   )
       call vsfm_mpp%MeshSetID                  (imesh , MESH_SPAC_XYLEM_COL)
       call vsfm_mpp%MeshSetDimensions          (imesh , ncells_xylem, ncells_ghost, nz_xylem )
       call vsfm_mpp%MeshSetGridCellFilter      (imesh , xylem_filter                         )
       call vsfm_mpp%MeshSetGeometricAttributes (imesh , VAR_XC     , xylem_xc                )
       call vsfm_mpp%MeshSetGeometricAttributes (imesh , VAR_YC     , xylem_yc                )
       call vsfm_mpp%MeshSetGeometricAttributes (imesh , VAR_ZC     , xylem_zc                )
       call vsfm_mpp%MeshSetGeometricAttributes (imesh , VAR_DX     , xylem_dx                )
       call vsfm_mpp%MeshSetGeometricAttributes (imesh , VAR_DY     , xylem_dy                )
       call vsfm_mpp%MeshSetGeometricAttributes (imesh , VAR_DZ     , xylem_dz                )
       call vsfm_mpp%MeshSetGeometricAttributes (imesh , VAR_AREA   , xylem_area              )
       call vsfm_mpp%MeshSetGeometricAttributes (imesh , VAR_VOLUME , xylem_vol               )

       !
       ! Vertical connections
       !
       iconn = 0
       do kk = 1, nz_xylem-1

          iconn                     = iconn + 1
          xylem_conn_id_up(iconn)   = kk
          xylem_conn_id_dn(iconn)   = kk + 1
          xylem_conn_dist_up(iconn) = 0.5d0*xylem_dz(kk  )
          xylem_conn_dist_dn(iconn) = 0.5d0*xylem_dz(kk+1)
          xylem_conn_area(iconn)    = xylem_area(kk)
          xylem_conn_type(iconn)    = CONN_VERTICAL
       end do
       nconn = iconn

       call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL,  &
            nconn,  xylem_conn_id_up, xylem_conn_id_dn,              &
            xylem_conn_dist_up, xylem_conn_dist_dn, xylem_conn_area, &
            xylem_conn_type)

    else

       imesh        = 1

       call vsfm_mpp%MeshSetName                (imesh, 'Soil+Root+Xylem mesh')
       call vsfm_mpp%MeshSetOrientation         (imesh, MESH_ALONG_GRAVITY)
       call vsfm_mpp%MeshSetID                  (imesh, MESH_CLM_SOIL_COL)
       call vsfm_mpp%MeshSetDimensions          (imesh, ncells_srx, ncells_ghost, ncells_soil)

       call vsfm_mpp%MeshSetGridCellFilter      (imesh, srx_filter)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , srx_xc)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , srx_yc)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , srx_zc)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , srx_dx)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , srx_dy)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , srx_dz)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , srx_area)
       call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME  , srx_vol)

       !
       ! Vertical connections
       !
       iconn = 0
       do kk = 1, nz_soil-1

          iconn                   = iconn + 1
          srx_conn_id_up(iconn)   = kk
          srx_conn_id_dn(iconn)   = kk + 1
          srx_conn_dist_up(iconn) = 0.5d0*soil_dz(kk  )
          srx_conn_dist_dn(iconn) = 0.5d0*soil_dz(kk+1)
          srx_conn_area(iconn)    = soil_area(kk)
          srx_conn_type(iconn)    = CONN_VERTICAL
       end do

       do kk = 1, nz_root-1
          iconn                    = iconn + 1
          srx_conn_id_up(iconn)   = kk + ncells_soil
          srx_conn_id_dn(iconn)   = kk + 1 + ncells_soil
          srx_conn_dist_up(iconn) = 0.5d0*root_dz(kk  )
          srx_conn_dist_dn(iconn) = 0.5d0*root_dz(kk+1)
          srx_conn_area(iconn)    = root_area(kk)
          srx_conn_type(iconn)    = CONN_VERTICAL
       end do

       do kk = 1, nz_xylem-1

          iconn                   = iconn + 1
          srx_conn_id_up(iconn)   = kk + ncells_soil + ncells_root
          srx_conn_id_dn(iconn)   = kk + 1 + ncells_soil + ncells_root
          srx_conn_dist_up(iconn) = 0.5d0*xylem_dz(kk  )
          srx_conn_dist_dn(iconn) = 0.5d0*xylem_dz(kk+1)
          srx_conn_area(iconn)    = xylem_area(kk)
          srx_conn_type(iconn)    = CONN_VERTICAL
       end do

       do kk = 1, nz_root
          iconn                   = iconn + 1
          srx_conn_id_up(iconn)   = kk
          srx_conn_id_dn(iconn)   = kk + nz_soil
          srx_conn_dist_up(iconn) = 0.5d0*root_radius
          srx_conn_dist_dn(iconn) = 0.5d0*root_radius
          srx_conn_area(iconn)    = soil_area(kk)
          srx_conn_type(iconn)    = CONN_VERTICAL
       end do

       iconn                   = iconn + 1
       srx_conn_id_up(iconn)   = nz_xylem + nz_root + nz_soil
       srx_conn_id_dn(iconn)   = nz_soil + 1
       srx_conn_dist_up(iconn) = 0.5d0*xylem_dz(nz_xylem)
       srx_conn_dist_dn(iconn) = 0.5d0*root_dz(1)
       srx_conn_area(iconn)    = xylem_area(kk)
       srx_conn_type(iconn)    = CONN_VERTICAL

       nconn = iconn

       call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
            nconn,  srx_conn_id_up, srx_conn_id_dn,                 &
            srx_conn_dist_up, srx_conn_dist_dn, srx_conn_area,      &
            srx_conn_type)
    endif

    deallocate(soil_xc            )
    deallocate(soil_yc            )
    deallocate(soil_zc            )
    deallocate(soil_dx            )
    deallocate(soil_dy            )
    deallocate(soil_dz            )
    deallocate(soil_area          )
    deallocate(soil_filter        )
    deallocate(soil_vol           )

    deallocate(root_xc            )
    deallocate(root_yc            )
    deallocate(root_zc            )
    deallocate(root_dx            )
    deallocate(root_dy            )
    deallocate(root_dz            )
    deallocate(root_area          )
    deallocate(root_filter        )
    deallocate(root_vol           )
    deallocate(root_len_den       )
    deallocate(root_surf_area     )

    deallocate(xylem_xc            )
    deallocate(xylem_yc            )
    deallocate(xylem_zc            )
    deallocate(xylem_dx            )
    deallocate(xylem_dy            )
    deallocate(xylem_dz            )
    deallocate(xylem_area          )
    deallocate(xylem_filter        )
    deallocate(xylem_vol           )

    deallocate(soil_conn_id_up    )
    deallocate(soil_conn_id_dn    )
    deallocate(soil_conn_dist_up  )
    deallocate(soil_conn_dist_dn  )
    deallocate(soil_conn_area     )
    deallocate(soil_conn_type     )
    deallocate(root_conn_id_up    )
    deallocate(root_conn_id_dn    )
    deallocate(root_conn_dist_up  )
    deallocate(root_conn_dist_dn  )
    deallocate(root_conn_area     )
    deallocate(root_conn_type     )
    deallocate(xylem_conn_id_up   )
    deallocate(xylem_conn_id_dn   )
    deallocate(xylem_conn_dist_up )
    deallocate(xylem_conn_dist_dn )
    deallocate(xylem_conn_area    )
    deallocate(xylem_conn_type    )

  end subroutine add_meshes

  !------------------------------------------------------------------------
  
  subroutine add_goveqns(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_RE
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_XYLEM_COL
    !
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp

    if (.not. single_pde_formulation) then
       call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Soil', MESH_CLM_SOIL_COL)
       call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Root', MESH_SPAC_ROOT_COL)
       call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Xylem', MESH_SPAC_XYLEM_COL)
    else
       call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Soil+Root+Xylem', MESH_CLM_SOIL_COL)
    endif

    call vsfm_mpp%SetMeshesOfGoveqns()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_MASS_RATE
    use MultiPhysicsProbConstants , only : SOIL_CELLS
    use ConnectionSetType         , only : connection_set_type
    use ConnectionSetType         , only : ConnectionSetDestroy
    use MeshType                  , only : MeshCreateConnectionSet
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscInt                            :: ieqn, ieqn_other
    PetscInt                            :: kk
    PetscInt                            :: nconn
    PetscInt                            :: ncells_root
    PetscReal                           :: root_len
    PetscReal                 , pointer :: root_len_den(:)
    PetscReal                 , pointer :: root_vol(:)
    PetscReal                 , pointer :: root_surf_area(:)
    PetscInt                  , pointer :: id_up(:)
    PetscInt                  , pointer :: id_dn(:)
    PetscReal                 , pointer :: dist_up(:)
    PetscReal                 , pointer :: dist_dn(:)
    PetscReal                 , pointer :: root_zc(:)           ! z-position of grid cell [m]
    PetscReal                 , pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
    PetscReal                 , pointer :: area(:)
    PetscInt                  , pointer :: itype(:)
    PetscReal                 , pointer :: unit_vec(:,:)
    type(connection_set_type) , pointer :: conn_set

    if (single_pde_formulation) return

    ncells_root  = nx_root*ny_root*nz_root

    allocate(soil_dz            (ncells_root))
    allocate(root_zc            (ncells_root))
    allocate(root_len_den       (ncells_root))
    allocate(root_surf_area     (ncells_root))
    
    soil_dz     (:) = dz
    root_zc     (:) = 0.d0
    do kk = 2,nz_root
       root_zc(kk) = root_zc(kk-1) - dz
    enddo

    do kk = 1, nz_root

       root_len_den(kk) = rld0 * (1.d0 - &
            abs(root_zc(kk))/d)*exp(-qz*abs(root_zc(kk))/d)

       root_len = root_len_den(kk)*     &       ! [m/m^3]
            (dx*dy*soil_dz(kk))           ! [m^3]
       root_surf_area(kk) = 2.d0*pi*root_radius*root_len   ! [m^2]

    end do
    root_surf_area(:) = dx*dy

    nconn         = nz_root

    allocate(id_up    (nconn   ))
    allocate(id_dn    (nconn   ))
    allocate(dist_up  (nconn   ))
    allocate(dist_dn  (nconn   ))
    allocate(area     (nconn   ))
    allocate(itype    (nconn   ))
    allocate(unit_vec (nconn,3 ))

    do kk = 1, nz_root
       id_up(kk)      = 0
       id_dn(kk)      = kk
       dist_up(kk)    = 0.d0
       dist_dn(kk)    = root_radius/2.d0
       area(kk)       = root_surf_area(kk)
       unit_vec(kk,1) = -1.d0
       unit_vec(kk,2) = 0.d0
       unit_vec(kk,3) = 0.d0
       itype(kk)      = CONN_VERTICAL
    enddo

    allocate(conn_set)

    ieqn       = 1
    ieqn_other = 2

    call MeshCreateConnectionSet(vsfm_mpp%meshes(1), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%GovEqnAddCondition(ieqn, ss_or_bc_type=COND_BC,   &
         name='Root BC in soil equation', unit='Pa', cond_type=COND_DIRICHLET_FRM_OTR_GOVEQ, &
         region_type=SOIL_TOP_CELLS, id_of_other_goveq=ieqn_other, &
         conn_set=conn_set)

    call ConnectionSetDestroy(conn_set)

    ieqn       = 2
    ieqn_other = 1

    unit_vec(:,1) = 1.d0

    call MeshCreateConnectionSet(vsfm_mpp%meshes(2), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%GovEqnAddCondition(ieqn, ss_or_bc_type=COND_BC,   &
         name='Soil BC in root equation', unit='Pa', cond_type=COND_DIRICHLET_FRM_OTR_GOVEQ, &
         region_type=SOIL_TOP_CELLS, id_of_other_goveq=ieqn_other, &
         conn_set=conn_set)

    call ConnectionSetDestroy(conn_set)

    ieqn       = 2
    ieqn_other = 3

    call vsfm_mpp%GovEqnAddCondition(ieqn, ss_or_bc_type=COND_BC,   &
         name='Xylem BC in root equation', unit='Pa', cond_type=COND_DIRICHLET_FRM_OTR_GOVEQ, &
         region_type=SOIL_TOP_CELLS, id_of_other_goveq=ieqn_other)

    ieqn       = 3
    ieqn_other = 2

    call vsfm_mpp%GovEqnAddCondition(ieqn, ss_or_bc_type=COND_BC,   &
         name='Root BC in xylem equation', unit='Pa', cond_type=COND_DIRICHLET_FRM_OTR_GOVEQ, &
         region_type=SOIL_BOTTOM_CELLS, id_of_other_goveq=ieqn_other)


    ! Coupling ALM with VSFM-SPAC

    ieqn       = 1
    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Infiltration_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Dew_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Drainage_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_CELLS)

    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Snow_Disappearance_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Sublimation_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    ieqn       = 3
    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Evapotranspiration_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_CELLS)


    deallocate(soil_dz        )
    deallocate(root_zc        )
    deallocate(root_len_den   )
    deallocate(root_surf_area )
    deallocate(id_up          )
    deallocate(id_dn          )
    deallocate(dist_up        )
    deallocate(dist_dn        )
    deallocate(area           )
    deallocate(itype          )
    deallocate(unit_vec       )

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    !
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    integer           :: ieqn
    integer           :: nvars_for_coupling
    integer, pointer  :: var_ids_for_coupling(:)
    integer, pointer  :: goveqn_ids_for_coupling(:)
    integer, pointer  :: is_bc(:)

    if (.not. single_pde_formulation) then
       !
       ! SOIL <---> ROOT
       !
       ieqn               = 1
       nvars_for_coupling = 1

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_PRESSURE
       goveqn_ids_for_coupling (1) = 2

       call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)

       !
       ! ROOT <---> SOIL
       ! ROOT <---> XYLEM
       !
       ieqn               = 2
       nvars_for_coupling = 2

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_PRESSURE
       var_ids_for_coupling    (2) = VAR_PRESSURE
       goveqn_ids_for_coupling (1) = 1
       goveqn_ids_for_coupling (2) = 3

       call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)

       !
       ! XYLEM <---> ROOT
       !
       ieqn               = 3
       nvars_for_coupling = 1

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_PRESSURE
       goveqn_ids_for_coupling (1) = 2

       call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)
    endif

    !
    ! Allocate auxvars
    !
    call vsfm_mpp%AllocateAuxVars()

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine set_material_properties(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoils
    use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use mpp_varcon                , only : denh2o
    use mpp_varcon                , only : grav
    use SystemOfEquationsBasePointerType,only : SOEIFunction, SOEIJacobian
    use GoverningEquationBaseType
    use GoveqnRichardsODEPressureType
    use RichardsODEPressureAuxType
    use SaturationFunction, only : SatFunc_Set_VG
    use SaturationFunction, only : SatFunc_Set_Weibull_RelPerm
    use PorosityFunctionMod, only : PorosityFunctionSetConstantModel
    use ConditionType, only : condition_type
    use ConnectionSetType, only : connection_set_type
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants , only : MESH_SPAC_XYLEM_COL
    !
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscReal        , parameter :: porosity            = 0.368d0
    PetscReal        , parameter :: lambda              = 0.5d0
    PetscReal        , parameter :: alpha               = 3.4257d-4
    PetscReal        , parameter :: perm                = 8.3913d-12
    !
    PetscReal        , pointer   :: vsfm_watsat(:,:)
    PetscReal        , pointer   :: vsfm_hksat(:,:)
    PetscReal        , pointer   :: vsfm_bsw(:,:)
    PetscReal        , pointer   :: vsfm_sucsat(:,:)
    PetscReal        , pointer   :: vsfm_eff_porosity(:,:)
    PetscReal        , pointer   :: vsfm_residual_sat(:,:)
    PetscReal        , parameter :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC
    PetscInt                     :: begc , endc
    integer          , pointer   :: vsfm_filter(:)
    character(len=32)            :: satfunc_type
    !
    class(goveqn_base_type),pointer                     :: cur_goveq
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_in(:)  !!< Internal state.
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_bc(:)  !!< Boundary conditions.
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_ss(:)  !!< Source-sink.
    type(condition_type),pointer                        :: cur_cond
    type(connection_set_type), pointer                  :: cur_conn_set
    PetscInt                                            :: ghosted_id
    PetscInt                                            :: sum_conn
    PetscInt                                            :: iconn
    PetscInt                                            :: ncells
    PetscInt                                            :: neqn
    !-----------------------------------------------------------------------

    begc = 1
    endc = nx*ny

    satfunc_type = 'van_genuchten'


    cur_goveq => vsfm_mpp%sysofeqns%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
          class is (goveqn_richards_ode_pressure_type)

             ! Set properties for internal auxvars
          aux_vars_in => cur_goveq%aux_vars_in

          select case(cur_goveq%mesh_itype)
          case (MESH_CLM_SOIL_COL)
             do ghosted_id = 1,cur_goveq%mesh%ncells_local

                if (ghosted_id <= nx_soil*ny_soil*nz_soil_bot) then

                   aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_bot
                   aux_vars_in(ghosted_id)%perm(3)              = perm_z_bot
                   aux_vars_in(ghosted_id)%por                  = por_bot

                   call PorosityFunctionSetConstantModel(    &
                        aux_vars_in(ghosted_id)%porParams, &
                        por_bot)

                   call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                        sat_res_bot,                        &
                        alpha_bot,                          &
                        vg_m_bot)

                else if (ghosted_id <= nx_soil*ny_soil*(nz_soil_bot+nz_soil_mid)) then

                   aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_mid
                   aux_vars_in(ghosted_id)%perm(3)              = perm_z_mid
                   aux_vars_in(ghosted_id)%por                  = por_mid

                   call PorosityFunctionSetConstantModel(    &
                        aux_vars_in(ghosted_id)%porParams, &
                        por_mid)

                   call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                        sat_res_mid,                        &
                        alpha_mid,                          &
                        vg_m_mid)
                else

                   aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_top
                   aux_vars_in(ghosted_id)%perm(3)              = perm_z_top
                   aux_vars_in(ghosted_id)%por                  = por_top

                   call PorosityFunctionSetConstantModel(  &
                        aux_vars_in(ghosted_id)%porParams, &
                        por_top)

                   call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                        sat_res_top,                                       &
                        alpha_top,                                         &
                        vg_m_top)
                end if

                aux_vars_in(ghosted_id)%pressure_prev        = press_initial
             enddo
          case (MESH_SPAC_ROOT_COL)

             do ghosted_id = 1,cur_goveq%mesh%ncells_local
                !                                                [1/s] *[m]*[Ns/m^2]/[kg/m^3]/[m/s^2]
                !                aux_vars_in(ghosted_id)%perm(1)              = cond_plant*dx*8.9d-4/denh2o/9.8068d0
                !                aux_vars_in(ghosted_id)%perm(2)              = cond_plant*dy*8.9d-4/denh2o/9.8068d0
                !                aux_vars_in(ghosted_id)%perm(3)              = cond_plant*dz*8.9d-4/denh2o/9.8068d0
                !                aux_vars_in(ghosted_id)%por                  = por_plant

                aux_vars_in(ghosted_id)%perm(1:2)            = perm_root
                aux_vars_in(ghosted_id)%perm(3)              = perm_root
                aux_vars_in(ghosted_id)%por                  = por_root

                call PorosityFunctionSetConstantModel(  &
                     aux_vars_in(ghosted_id)%porParams, &
                     por_root)

                call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                     sat_res_root,                                      &
                     alpha_root,                                        &
                     vg_m_root)

                ! Use Weibull relative perm model
                !call SatFunc_Set_Weibull_RelPerm(aux_vars_in(ghosted_id)%satParams, &
                !       weibull_d, weibull_c)

                aux_vars_in(ghosted_id)%pressure_prev        = press_initial
             enddo

          case (MESH_SPAC_XYLEM_COL)

             do ghosted_id = 1,cur_goveq%mesh%ncells_local
                !                                                [1/s] *[m]*[Ns/m^2]/[kg/m^3]/[m/s^2]
                !                aux_vars_in(ghosted_id)%perm(1)              = cond_plant*dx*8.9d-4/denh2o/9.8068d0
                !                aux_vars_in(ghosted_id)%perm(2)              = cond_plant*dy*8.9d-4/denh2o/9.8068d0
                !                aux_vars_in(ghosted_id)%perm(3)              = cond_plant*dz*8.9d-4/denh2o/9.8068d0
                !                aux_vars_in(ghosted_id)%por                  = por_plant

                aux_vars_in(ghosted_id)%perm(1:2)            = perm_xylem
                aux_vars_in(ghosted_id)%perm(3)              = perm_xylem
                aux_vars_in(ghosted_id)%por                  = por_xylem

                call PorosityFunctionSetConstantModel(  &
                     aux_vars_in(ghosted_id)%porParams, &
                     por_xylem)

                call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                     sat_res_top,                                       &
                     alpha_xylem,                                       &
                     vg_m_xylem)

                ! Use Weibull relative perm model
                !call SatFunc_Set_Weibull_RelPerm(aux_vars_in(ghosted_id)%satParams, &
                !       weibull_d, weibull_c)

                aux_vars_in(ghosted_id)%pressure_prev        = press_initial
             enddo
          case default
             write(*,*)'SPACMPPSetupPetscSNESSetup: Unknown mesh_itype.'
             stop
          end select

          ! Set hydraulic properties for boundary-condition auxvars
          aux_vars_bc => cur_goveq%aux_vars_bc
          sum_conn = 0
          cur_cond => cur_goveq%boundary_conditions%first
          do
             if (.not.associated(cur_cond)) exit
             cur_conn_set => cur_cond%conn_set

             do iconn = 1, cur_conn_set%num_connections
                sum_conn = sum_conn + 1
                ghosted_id = cur_conn_set%id_dn(iconn)

                aux_vars_bc(sum_conn)%perm(:)             = 6.83d-21 !aux_vars_in(ghosted_id)%perm(:)
                aux_vars_bc(sum_conn)%por                 = aux_vars_in(ghosted_id)%por
                aux_vars_bc(sum_conn)%satParams           = aux_vars_in(ghosted_id)%satParams
                aux_vars_bc(sum_conn)%porParams           = aux_vars_in(ghosted_id)%porParams

             enddo
             cur_cond => cur_cond%next
          enddo

          ! Set hydraulic properties for source-sink auxvars
          aux_vars_ss => cur_goveq%aux_vars_ss
          sum_conn = 0
          cur_cond => cur_goveq%source_sinks%first
          do
             if (.not.associated(cur_cond)) exit
             cur_conn_set => cur_cond%conn_set

             do iconn = 1, cur_conn_set%num_connections
                sum_conn = sum_conn + 1
                ghosted_id = cur_conn_set%id_dn(iconn)

                aux_vars_ss(sum_conn)%perm(:)             = aux_vars_in(ghosted_id)%perm(:)
                aux_vars_ss(sum_conn)%por                 = aux_vars_in(ghosted_id)%por
                aux_vars_ss(sum_conn)%satParams           = aux_vars_in(ghosted_id)%satParams
                aux_vars_ss(sum_conn)%porParams           = aux_vars_in(ghosted_id)%porParams

             enddo
             cur_cond => cur_cond%next
          enddo

          class default
          write(*,*)'SPACMPPSetupPetscSNESSetup: Unknown goveq type.'
          stop
       end select

       cur_goveq => cur_goveq%next
    enddo

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscErrorCode                                    :: ierr

    !
    call VecSet (vsfm_mpp%sysofeqns%soln, press_initial, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_mpp%sysofeqns%soln, vsfm_mpp%sysofeqns%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_mpp%sysofeqns%soln, vsfm_mpp%sysofeqns%soln_prev_clm, ierr); CHKERRQ(ierr)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine setup_petsc_snes(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    use mpp_varctl                       , only : iulog
    use mpp_abortutils                   , only : endrun
    use mpp_shr_log_mod                  , only : errMsg => shr_log_errMsg
    use GoverningEquationBaseType        , only : goveqn_base_type
    use GoveqnRichardsODEPressureType    , only : goveqn_richards_ode_pressure_type
    use SystemOfEquationsVSFMType        , only : sysofeqns_vsfm_type
    use SystemOfEquationsBasePointerType , only : sysofeqns_base_pointer_type
    use SystemOfEquationsBasePointerType , only : SOEResidual
    use SystemOfEquationsBasePointerType , only : SOEJacobian
    use SystemOfEquationsVSFMType        , only : VSFMSOESetAuxVars
    use MultiPhysicsProbConstants        , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants        , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants        , only : MPP_VSFM_SNES_CLM
    use MultiPhysicsProbConstants        , only : SOE_RE_ODE
    use mpp_abortutils                   , only : endrun
    use petscmat
    use petscdm
    use petscdmda
    !
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscInt, pointer                                 :: mesh_size(:)
    PetscInt :: igoveqn
    DM                                                :: dm_s                ! DM object for soil equation
    DM                                                :: dm_r                ! DM object for root equation
    DM                                                :: dm_x                ! DM object for xylem equation
    PetscErrorCode                                    :: ierr
    class(goveqn_base_type),pointer                   :: cur_goveq
    class (goveqn_richards_ode_pressure_type),pointer :: goveq_richards_pres
    class(sysofeqns_vsfm_type),pointer                :: vsfm_soe
    PetscReal, parameter                              :: atol    = PETSC_DEFAULT_REAL
    PetscReal, parameter                              :: rtol    = PETSC_DEFAULT_REAL
    PetscReal, parameter                              :: stol    = 1.d-10
    PetscInt, parameter                               :: max_it  = PETSC_DEFAULT_INTEGER
    PetscInt, parameter                               :: max_f   = PETSC_DEFAULT_INTEGER
    !

    vsfm_soe => vsfm_mpp%sysofeqns
    vsfm_mpp%sysofeqns_ptr%ptr => vsfm_mpp%sysofeqns

    allocate(mesh_size(vsfm_soe%ngoveqns))

    ! Get pointers to governing-equations

    igoveqn = 0
    cur_goveq => vsfm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
          class is (goveqn_richards_ode_pressure_type)
          igoveqn = igoveqn + 1
          mesh_size(igoveqn) = cur_goveq%mesh%ncells_local
       end select

       cur_goveq => cur_goveq%next
    enddo

    ! DM-Composite approach

    ! Create PETSc DM for pressure-equation
    !size = goveq_richards_pres%mesh%ncells_local

    select case(vsfm_mpp%id)
    case (MPP_VSFM_SNES_CLM)
       call DMDACreate1d(PETSC_COMM_SELF,     &
            DM_BOUNDARY_NONE,    &
            mesh_size(1),        &
            1,                   &
            1,                   &
            PETSC_NULL_INTEGER,  &
            dm_s,                &
            ierr); CHKERRQ(ierr)

       if (.not.single_pde_formulation) then
          call DMDACreate1d(PETSC_COMM_SELF,     &
               DM_BOUNDARY_NONE,    &
               mesh_size(2),        &
               1,                   &
               1,                   &
               PETSC_NULL_INTEGER,  &
               dm_r,                &
               ierr); CHKERRQ(ierr)

          call DMDACreate1d(PETSC_COMM_SELF,     &
               DM_BOUNDARY_NONE,    &
               mesh_size(3),        &
               1,                   &
               1,                   &
               PETSC_NULL_INTEGER,  &
               dm_x,                &
               ierr); CHKERRQ(ierr)
       endif
    case default
       write(iulog,*)'VSFMMPPSetupPetscSNESSetup: Unknown vsfm_mpp%id'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    call DMSetOptionsPrefix (dm_s , "fs_" , ierr           ); CHKERRQ(ierr)
    call DMSetFromOptions   (dm_s , ierr                   ); CHKERRQ(ierr)
    call DMSetUp            (dm_s , ierr                       ); CHKERRQ(ierr)
    call DMDASetFieldName   (dm_s , 0     , "soil_" , ierr ); CHKERRQ(ierr)

    if (.not.single_pde_formulation) then
       call DMSetOptionsPrefix (dm_r , "fr_" , ierr           ); CHKERRQ(ierr)
       call DMSetFromOptions   (dm_r , ierr                   ); CHKERRQ(ierr)
       call DMSetUp            (dm_r , ierr                   ); CHKERRQ(ierr)
       call DMDASetFieldName   (dm_r , 0     , "root_" , ierr ); CHKERRQ(ierr)

       call DMSetOptionsPrefix (dm_x , "fx_" , ierr           ); CHKERRQ(ierr)
       call DMSetFromOptions   (dm_x , ierr                   ); CHKERRQ(ierr)
       call DMSetUp            (dm_x , ierr                   ); CHKERRQ(ierr)
       call DMDASetFieldName   (dm_x , 0     , "xylem_" , ierr ); CHKERRQ(ierr)
    endif

    ! Create DMComposite: pressure
    call DMCompositeCreate  (PETSC_COMM_SELF , vsfm_soe%dm , ierr ); CHKERRQ(ierr)
    call DMSetOptionsPrefix (vsfm_soe%dm     , "pressure_" , ierr ); CHKERRQ(ierr)

    call DMCompositeAddDM   (vsfm_soe%dm     , dm_s        , ierr ); CHKERRQ(ierr)
    if (.not.single_pde_formulation) then
       call DMCompositeAddDM   (vsfm_soe%dm     , dm_r        , ierr ); CHKERRQ(ierr)
       call DMCompositeAddDM   (vsfm_soe%dm     , dm_x        , ierr ); CHKERRQ(ierr)
    endif

    call DMDestroy          (dm_s            , ierr               ); CHKERRQ(ierr)
    if (.not.single_pde_formulation) then
       call DMDestroy          (dm_r            , ierr               ); CHKERRQ(ierr)
       call DMDestroy          (dm_x            , ierr               ); CHKERRQ(ierr)
    endif

    call DMSetUp            (vsfm_soe%dm     , ierr               ); CHKERRQ(ierr)

    call DMCreateMatrix     (vsfm_soe%dm     , vsfm_soe%jac, ierr ); CHKERRQ(ierr)
    call MatSetOption       (vsfm_soe%jac, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
    call MatSetOption       (vsfm_soe%jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)

    ! Create solution vector
    call DMCreateGlobalVector(vsfm_soe%dm , vsfm_soe%soln          , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(vsfm_soe%dm , vsfm_soe%soln_prev     , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(vsfm_soe%dm , vsfm_soe%soln_prev_clm , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(vsfm_soe%dm , vsfm_soe%res          , ierr); CHKERRQ(ierr)

    call VecZeroEntries(vsfm_soe%soln         , ierr); CHKERRQ(ierr)
    call VecZeroEntries(vsfm_soe%soln_prev    , ierr); CHKERRQ(ierr)
    call VecZeroEntries(vsfm_soe%soln_prev_clm, ierr); CHKERRQ(ierr)
    call VecZeroEntries(vsfm_soe%res          , ierr); CHKERRQ(ierr)

    ! SNES
    call SNESCreate(PETSC_COMM_SELF, vsfm_soe%snes, ierr); CHKERRQ(ierr)
    call SNESSetTolerances(vsfm_soe%snes, atol, rtol, stol, &
         max_it, max_f, ierr); CHKERRQ(ierr)

    call SNESSetFunction(vsfm_soe%snes, vsfm_soe%res, SOEResidual, &
         vsfm_mpp%sysofeqns_ptr, ierr); CHKERRQ(ierr)

    call SNESSetJacobian(vsfm_soe%snes, vsfm_soe%jac, vsfm_soe%jac,     &
         SOEJacobian, vsfm_mpp%sysofeqns_ptr, ierr); CHKERRQ(ierr)

    call SNESSetFromOptions(vsfm_soe%snes, ierr); CHKERRQ(ierr)

    ! Get pointers to governing-equations
    call vsfm_soe%CreateVectorsForGovEqn()

    vsfm_mpp%sysofeqns%solver_type = vsfm_mpp%solver_type
    vsfm_mpp%sysofeqns%itype       = SOE_RE_ODE

  end subroutine setup_petsc_snes

  !-----------------------------------------------------------------------
  subroutine determine_condition_ids(this)
    !
    !DESCRIPTION
    !  Determines the IDs of various source-sink conditions in VSFM
    !
    use mpp_varctl                       , only : vsfm_lateral_model_type
    use MultiPhysicsProbConstants        , only : COND_SS
    use MultiPhysicsProbConstants        , only : COND_NULL
    use mpp_varctl                       , only : iulog
    use abortutils                       , only : endrun
    use shr_log_mod                      , only : errMsg => shr_log_errMsg
    ! !ARGUMENTS:
    implicit none
    !
    ! !ARGUMENTS
    class(em_vsfm_type)          :: this
    integer                      :: ier ! error status
    !
    ! !LOCAL VARIABLES:
    character (len=256), pointer :: cond_names(:)
    integer                      :: num_conds
    integer                      :: num_conds_expected
    integer                      :: nn
    integer                      :: kk
    character (len=256)          :: cond_name
    !------------------------------------------------------------------------------
    
    this%vsfm_cond_id_for_infil        = -1
    this%vsfm_cond_id_for_et           = -1
    this%vsfm_cond_id_for_dew          = -1
    this%vsfm_cond_id_for_drainage     = -1
    this%vsfm_cond_id_for_snow         = -1
    this%vsfm_cond_id_for_sublimation  = -1
    this%vsfm_cond_id_for_lateral_flux = -1

    num_conds_expected = 6

    if (vsfm_lateral_model_type == 'source_sink' ) then
       num_conds_expected = num_conds_expected + 1
    end if

    ! Get the number of conditions
    call this%vsfm_mpp%sysofeqns%GetConditionNames(COND_SS, COND_NULL, num_conds, cond_names)

    if (num_conds /= num_conds_expected) then
       write(iulog,*)'In init_vsfm_condition_ids: Source-sink conditions /= ', num_conds_expected
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do nn = 1, num_conds
       select case(trim(cond_names(nn)))

       case ("Infiltration_Flux")
          this%vsfm_cond_id_for_infil        = nn

       case ("Evapotranspiration_Flux")
          this%vsfm_cond_id_for_et           = nn

       case ("Dew_Flux")
          this%vsfm_cond_id_for_dew          = nn

       case ("Drainage_Flux")
          this%vsfm_cond_id_for_drainage     = nn

       case ("Snow_Disappearance_Flux")
          this%vsfm_cond_id_for_snow         = nn

       case ("Sublimation_Flux")
          this%vsfm_cond_id_for_sublimation  = nn

       case ("Lateral_flux")
          this%vsfm_cond_id_for_lateral_flux = nn

       case default
          write(iulog,*) trim(cond_names(nn))
          write(iulog,*)'In init_vsfm_condition_ids: Unknown flux.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    enddo

    if (this%vsfm_cond_id_for_infil == -1) then
       write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_infil not defined.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%vsfm_cond_id_for_et == -1) then
       write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_et not defined.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%vsfm_cond_id_for_dew == -1) then
       write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_dew not defined.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%vsfm_cond_id_for_drainage == -1) then
       write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_drainage not defined.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%vsfm_cond_id_for_snow == -1) then
       write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_snow not defined.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (this%vsfm_cond_id_for_sublimation == -1) then
       write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_sublimation not defined.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_lateral_model_type == 'source_sink') then
       if (this%vsfm_cond_id_for_lateral_flux == -1) then
          write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_lateral_flux not defined.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

    endif

    deallocate(cond_names)

  end subroutine determine_condition_ids

  !-----------------------------------------------------------------------
  subroutine extract_data_for_alm(this, l2e_init_list, e2l_init_list, bounds_clump)
    !
    !DESCRIPTION
    !  Saves
    !
    use mpp_varctl                , only : restart_vsfm
    use mpp_bounds                , only : bounds_proc_begc, bounds_proc_endc
    use mpp_varpar                , only : nlevgrnd
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_vsfm_spac_type)             :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! !LOCAL VARIABLES:
    integer                              :: p,c,fc,j,g  ! do loop indices

    real(r8)  , pointer                  :: l2e_col_zi(:,:)

    real(r8)  , pointer                  :: l2e_soilp(:,:)
    real(r8)  , pointer                  :: l2e_mflx_snowlyr_col(:)

    real(r8)  , pointer                  :: vsfm_soilp_col_1d(:)
    real(r8)  , pointer                  :: vsfm_mass_col_1d(:)
    real(r8)  , pointer                  :: vsfm_smpl_col_1d(:)
    real(r8)  , pointer                  :: e2l_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: e2l_smp_l(:,:)
    real(r8)  , pointer                  :: e2l_zwt(:)
    real(r8)  , pointer                  :: e2l_mflx_snowlyr_col(:)
    integer                              :: jwt
    integer                              :: idx
    integer                              :: soe_auxvar_id
    real(r8)                             :: z_up, z_dn
    integer                              :: num
    real(r8)  , pointer                  :: e2l_h2oroot_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2oxylem_liq(:,:)
!    integer :: bounds_proc_begc, bounds_proc_endc

    !-----------------------------------------------------------------------
    call l2e_init_list%GetPointerToReal1D(this%index_l2e_init_flux_mflx_snowlyr_col , l2e_mflx_snowlyr_col )

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_state_soilp           , l2e_soilp            )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_zi                , l2e_col_zi           )

    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_state_wtd             , e2l_zwt              )
    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_flux_mflx_snowlyr_col , e2l_mflx_snowlyr_col )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_liq      , e2l_h2osoi_liq       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_ice      , e2l_h2osoi_ice       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_smp             , e2l_smp_l            )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2oroot_liq     , e2l_h2oroot_liq      )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2oxylem_liq    , e2l_h2oxylem_liq     )

    !bounds_proc_begc     = bounds_clump%begc
    !bounds_proc_endc     = bounds_clump%endc

    ! PreSolve: Allows saturation value to be computed based on ICs and stored
    !           in GE auxvar
    call this%vsfm_mpp%sysofeqns%SetDtime(1.d0)
    call this%vsfm_mpp%sysofeqns%PreSolve()

    ! PostSolve: Allows saturation value stored in GE auxvar to be copied into
    !            SoE auxvar
    call this%vsfm_mpp%sysofeqns%PostSolve()

    num = this%vsfm_mpp%sysofeqns%num_auxvars_in
    allocate(vsfm_soilp_col_1d(num))
    allocate(vsfm_mass_col_1d (num))
    allocate(vsfm_smpl_col_1d (num))

    if (restart_vsfm) then

       write(*,*)'In restart_vsfm not support in VSFM-SPAC'
       stop
       ! Save 1D array for VSFM (vsfm_soilp_col_1d) and
       ! set initial value of mflx_snowlyr_col for ALM
       do c = bounds_proc_begc, bounds_proc_endc
          do j = 1, nlevgrnd
             idx = (c - bounds_proc_begc)*nlevgrnd + j
             vsfm_soilp_col_1d(idx) = l2e_soilp(c,j)
          end do
          idx = c-bounds_proc_begc+1
          e2l_mflx_snowlyr_col(c) = l2e_mflx_snowlyr_col(c)
       end do

       ! Set the initial conditions
       call this%vsfm_mpp%Restart(vsfm_soilp_col_1d)

       ! PreSolve: Allows saturation value to be computed based on ICs and stored
       !           in GE auxvar
       call this%vsfm_mpp%sysofeqns%SetDtime(1.d0)
       call this%vsfm_mpp%sysofeqns%PreSolve()

       ! PostSolve: Allows saturation value stored in GE auxvar to be copied into
       !            SoE auxvar
       call this%vsfm_mpp%sysofeqns%PostSolve()
       
    else
       ! Set initial value of mflx_snowlyr_col for ALM
       e2l_mflx_snowlyr_col(:) = 0._r8
    end if

    ! Get total mass
    soe_auxvar_id = 1;
    call this%vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL,   &
         VAR_MASS,          &
         soe_auxvar_id,     &
         vsfm_mass_col_1d)

    ! Get liquid soil matrix potential
    soe_auxvar_id = 1;
    call this%vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL,       &
         VAR_SOIL_MATRIX_POT,   &
         soe_auxvar_id,         &
         vsfm_smpl_col_1d)

    do c = bounds_proc_begc, bounds_proc_endc
       ! initialization
       jwt = -1

       ! Loops in decreasing j so WTD can be computed in the same loop
       do j = nlevgrnd, 1, -1
          idx = (c-bounds_proc_begc)*nlevgrnd + j

          e2l_h2osoi_liq(c,j) = vsfm_mass_col_1d(idx)
          e2l_h2osoi_ice(c,j) = 0.d0
          e2l_smp_l(c,j)      = vsfm_smpl_col_1d(idx)*1000._r8      ! [m] --> [mm]

          if (jwt == -1) then
             ! Find the first soil that is unsaturated
             if (e2l_smp_l(c,j) < 0._r8) jwt = j
          end if

       end do

       if (jwt == -1 .or. jwt == nlevgrnd) then
          ! Water table below or in the last layer
          e2l_zwt(c) = l2e_col_zi(c,nlevgrnd)
       else
          z_dn = (l2e_col_zi(c,jwt-1) + l2e_col_zi(c,jwt  ))/2._r8
          z_up = (l2e_col_zi(c,jwt ) + l2e_col_zi(c,jwt+1))/2._r8
          e2l_zwt(c) = (0._r8 - e2l_smp_l(c,jwt))/(e2l_smp_l(c,jwt) - e2l_smp_l(c,jwt+1))*(z_dn - z_up) + z_dn
       endif

       do j = nz_soil + 1, nz_soil + nz_root
          idx = (c-bounds_proc_begc)*(nz_soil + nz_root + nz_xylem) + j
          e2l_h2oroot_liq(c,j - nz_soil) = vsfm_mass_col_1d(idx)
       enddo

       do j = nz_soil + nz_root + 1, nz_soil + nz_root + nz_xylem
          idx = (c-bounds_proc_begc)*(nz_soil + nz_root + nz_xylem) + j
          e2l_h2oxylem_liq(c,j - nz_soil - nz_root) = vsfm_mass_col_1d(idx)
      enddo

    enddo

  end subroutine extract_data_for_alm

  !------------------------------------------------------------------------
  subroutine EM_VSFM_SPAC_Populate_E2L_Init_List(this, e2l_init_list)
    !
    !
    ! !DESCRIPTION:
    ! Create a list of all variables to be returned by VSFM from ALM
    !
    ! !USES:
    use ExternalModelVSFMMod  , only : EM_VSFM_Populate_E2L_Init_List
    use ExternalModelConstants, only : EM_INITIALIZATION_STAGE
    use ExternalModelConstants, only : E2L_STATE_H2OROOT_LIQ
    use ExternalModelConstants, only : E2L_STATE_H2OXYLEM_LIQ
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_spac_type)             :: this
    class(emi_data_list) , intent(inout) :: e2l_init_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages

    call EM_VSFM_Populate_E2L_Init_List(this, e2l_init_list)

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    call e2l_init_list%AddDataByID(E2L_STATE_H2OROOT_LIQ, number_em_stages, &
                              em_stages, this%index_e2l_init_state_h2oroot_liq)
    call e2l_init_list%AddDataByID(E2L_STATE_H2OXYLEM_LIQ, number_em_stages, &
                              em_stages, this%index_e2l_init_state_h2oxylem_liq)

    deallocate(em_stages)

end subroutine EM_VSFM_SPAC_Populate_E2L_Init_List


  !------------------------------------------------------------------------
  subroutine EM_VSFM_SPAC_Populate_E2L_List(this, e2l_list)
    !
    !
    ! !DESCRIPTION:
    ! Create a list of all variables to be returned by VSFM from ALM
    !
    ! !USES:
    use ExternalModelVSFMMod  , only : EM_VSFM_Populate_E2L_List
    use ExternalModelConstants, only : EM_VSFM_SOIL_HYDRO_STAGE
    use ExternalModelConstants, only : E2L_STATE_H2OROOT_LIQ
    use ExternalModelConstants, only : E2L_STATE_H2OXYLEM_LIQ
    use ExternalModelConstants, only : E2L_STATE_XYLEM_MATRIC_POTENTIAL
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_spac_type)             :: this
    class(emi_data_list) , intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    call EM_VSFM_Populate_E2L_List(this, e2l_list)

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_VSFM_SOIL_HYDRO_STAGE

    call e2l_list%AddDataByID(E2L_STATE_H2OROOT_LIQ, number_em_stages, &
                              em_stages, this%index_e2l_state_h2oroot_liq)
    call e2l_list%AddDataByID(E2L_STATE_H2OXYLEM_LIQ, number_em_stages, &
                              em_stages, this%index_e2l_state_h2oxylem_liq)

    id = E2L_STATE_XYLEM_MATRIC_POTENTIAL
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_xylemp = index

    deallocate(em_stages)

end subroutine EM_VSFM_SPAC_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_VSFM_SPAC_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
      bounds_clump)
    !
    ! !DESCRIPTION:
    ! The VSFM dirver subroutine
    !
    ! !USES:
    use ExternalModelConstants , only : EM_VSFM_SOIL_HYDRO_STAGE
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_vsfm_spac_type)             :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    select case (em_stage)
    case (EM_VSFM_SOIL_HYDRO_STAGE)
       call EM_VSFM_SPAC_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, bounds_clump)
    case default
       write(iulog,*)'EM_FATES_Solve: Unknown em_stage.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine EM_VSFM_SPAC_Solve

  !------------------------------------------------------------------------
  subroutine EM_VSFM_SPAC_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, bounds_clump)
    !
    ! !DESCRIPTION:
    ! Solve the Variably Saturated Flow Model (VSFM) in soil columns.
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_FRAC_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use MultiPhysicsProbConstants , only : VAR_LATERAL_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : VAR_BC_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use mpp_varpar                , only : nlevgrnd
    use mpp_bounds                , only : bounds_proc_begc, bounds_proc_endc
    use petscsnes
    !
    implicit none
    !
    !
    ! !ARGUMENTS:
    class(em_vsfm_spac_type)             :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer                              :: p,c,fc,j,g                                                       ! do loop indices
    integer                              :: pi                                                               ! pft index
    real(r8)                             :: dzsum                                                            ! summation of dzmm of layers below water table (mm)
    real(r8)                             :: dtime

    real(r8)  , pointer                  :: mflx_et_col_1d         (:)
    real(r8)  , pointer                  :: mflx_infl_col_1d       (:)
    real(r8)  , pointer                  :: mflx_dew_col_1d        (:)
    real(r8)  , pointer                  :: mflx_drain_col_1d      (:)
    real(r8)  , pointer                  :: mflx_sub_snow_col_1d   (:)
    real(r8)  , pointer                  :: mflx_snowlyr_col_1d    (:)
    real(r8)  , pointer                  :: t_soil_col_1d          (:)

    real(r8)  , pointer                  :: vsfm_fliq_col_1d       (:)
    real(r8)  , pointer                  :: vsfm_mass_col_1d       (:)
    real(r8)  , pointer                  :: vsfm_smpl_col_1d       (:)
    real(r8)  , pointer                  :: vsfm_soilp_col_1d      (:)
    real(r8)  , pointer                  :: vsfm_sat_col_1d        (:)

    real(r8)  , pointer                  :: frac_ice                    (:,:) ! fraction of ice
    real(r8)  , pointer                  :: total_mass_flux_col         (:)            ! Sum of all source-sinks conditions for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_et_col      (:)            ! ET sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_infl_col    (:)            ! Infiltration source for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_dew_col     (:)            ! Dew source for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_drain_col   (:)            ! Drainage sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_snowlyr_col (:)            ! Flux due to disappearance of snow for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_sub_col     (:)            ! Sublimation sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_lateral_col (:)            ! Lateral flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_seepage_col (:)            ! Seepage flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: qflx_seepage                (:)            ! Seepage flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: vsfm_mass_prev_col          (:,:) ! Mass of water before a VSFM solve
    real(r8)  , pointer                  :: vsfm_dmass_col              (:)            ! Change in mass of water after a VSFM solve
    real(r8)  , pointer                  :: mass_beg_col                (:)            ! Total mass before a VSFM solve
    real(r8)  , pointer                  :: mass_end_col                (:)            ! Total mass after a VSFM solve
    integer                              :: ier                                                              ! error status

    integer                              :: begc, endc
    integer                              :: idx, idx_c
    real(r8)                             :: area

    PetscInt                             :: soe_auxvar_id                                                    ! Index of system-of-equation's (SoE's) auxvar
    PetscErrorCode                       :: ierr                                                             ! PETSc return error code

    PetscBool                            :: converged                                                        ! Did VSFM solver converge to a solution with given PETSc SNES tolerances
    PetscInt                             :: converged_reason                                                 ! SNES converged due to which criteria
    PetscReal                            :: atol_default                                                     ! Default SNES absolute convergance tolerance
    PetscReal                            :: rtol_default                                                     ! Default SNES relative convergance tolerance
    PetscReal                            :: stol_default                                                     ! Default SNES solution convergance tolerance
    PetscInt                             :: max_it_default                                                   ! Default SNES maximum number of iteration
    PetscInt                             :: max_f_default                                                    ! Default SNES maximum number of function evaluation
    PetscReal                            :: stol                                                             ! solution convergance tolerance
    PetscReal                            :: rtol                                                             ! relative convergance tolerance
    PetscReal,parameter                  :: stol_alternate = 1.d-10                                          ! Alternate solution convergance tolerance

    PetscReal                            :: mass_beg                                                         ! Sum of mass of water for all active soil columns before VSFM is called
    PetscReal                            :: mass_end                                                         ! Sum of mass of water for all active soil columns after VSFM is called
    PetscReal                            :: total_mass_flux_et                                               ! Sum of mass ET mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_infl                                             ! Sum of mass infiltration mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_dew                                              ! Sum of mass dew mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_drain                                            ! Sum of mass drainage mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_snowlyr                                          ! Sum of mass snow layer disappearance mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_sub                                              ! Sum of mass sublimation mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_lateral                                          ! Sum of lateral mass flux for all active soil columns
    PetscReal                            :: total_mass_flux                                                  ! Sum of mass ALL mass flux of water for all active soil columns
    PetscInt                             :: iter_count                                                       ! How many times VSFM solver is called

    PetscInt, parameter                  :: max_iter_count = 10                                              ! Maximum number of times VSFM can be called
    PetscInt                             :: diverged_count                                                   ! Number of time VSFM solver diverged
    PetscInt                             :: mass_bal_err_count                                               ! Number of time VSFM solver returns a solution that isn't within acceptable mass balance error threshold
    PetscReal                            :: abs_mass_error_col                                               ! Maximum absolute error for any active soil column
    PetscReal, parameter                 :: max_abs_mass_error_col  = 1.e-5                                  ! Acceptable mass balance error
    PetscBool                            :: successful_step                                                  ! Is the solution return by VSFM acceptable
    PetscReal , pointer                  :: vsfm_soilp_col_ghosted_1d(:)
    PetscReal , pointer                  :: vsfm_fliq_col_ghosted_1d(:)
    PetscReal , pointer                  :: mflx_lateral_col_1d(:)
    PetscReal , pointer                  :: lat_mass_exc_col_1d(:)
    PetscReal , pointer                  :: seepage_mass_exc_col_1d(:)
    PetscReal , pointer                  :: seepage_press_1d(:)

    integer                              :: jwt
    real(r8)                             :: z_dn, z_up

    real(r8)  , pointer                  :: l2e_mflux_infil(:)
    real(r8)  , pointer                  :: l2e_mflux_dew(:)
    real(r8)  , pointer                  :: l2e_mflux_sub_snow(:)
    real(r8)  , pointer                  :: l2e_mflux_snowlyr(:)
    real(r8)  , pointer                  :: l2e_mflux_et(:,:)
    real(r8)  , pointer                  :: l2e_mflux_drain(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: l2e_zi(:,:)
    integer   , pointer                  :: l2e_filter_hydrologyc(:)
    integer                              :: l2e_num_hydrologyc

    real(r8)  , pointer                  :: e2l_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: e2l_smp(:,:)
    real(r8)  , pointer                  :: e2l_wtd(:)
    real(r8)  , pointer                  :: e2l_soilp(:,:)
    real(r8)  , pointer                  :: e2l_xylemp(:,:)
    real(r8)  , pointer                  :: e2l_qrecharge(:)
    real(r8)                             :: et_integrated_over_root
    real(r8)  , pointer                  :: e2l_h2oroot_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2oxylem_liq(:,:)

    !integer                              :: bounds_proc_begc, bounds_proc_endc

    !-----------------------------------------------------------------------

    ! Get time step

    dtime = dt

    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_infil       , l2e_mflux_infil       )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_dew         , l2e_mflux_dew         )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_snow_sub    , l2e_mflux_sub_snow    )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_snowlyr     , l2e_mflux_snowlyr     )

    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_et          , l2e_mflux_et          )
    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_drainage    , l2e_mflux_drain       )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq , l2e_h2osoi_liq        )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice , l2e_h2osoi_ice        )

    call l2e_list%GetPointerToInt1D(this%index_l2e_filter_hydrologyc , l2e_filter_hydrologyc )
    call l2e_list%GetIntValue(this%index_l2e_filter_num_hydrologyc   , l2e_num_hydrologyc    )

    call l2e_list%GetPointerToReal2D(this%index_l2e_column_zi        , l2e_zi                )

    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2oroot_liq , e2l_h2oroot_liq      )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2oxylem_liq, e2l_h2oxylem_liq     )

    call e2l_list%GetPointerToReal1D(this%index_e2l_state_wtd        , e2l_wtd               )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_liq , e2l_h2osoi_liq        )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_ice , e2l_h2osoi_ice        )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_smp        , e2l_smp               )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_soilp      , e2l_soilp             )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_xylemp     , e2l_xylemp            )

    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_qrecharge   , e2l_qrecharge         )

    !bounds_proc_begc     = bounds_clump%begc
    !bounds_proc_endc     = bounds_clump%endc

    begc = bounds_proc_begc
    endc = bounds_proc_endc

    allocate(frac_ice                    (begc:endc,1:nlevgrnd))
    allocate(total_mass_flux_col         (begc:endc))
    allocate(total_mass_flux_et_col      (begc:endc))
    allocate(total_mass_flux_infl_col    (begc:endc))
    allocate(total_mass_flux_dew_col     (begc:endc))
    allocate(total_mass_flux_drain_col   (begc:endc))
    allocate(total_mass_flux_snowlyr_col (begc:endc))
    allocate(total_mass_flux_sub_col     (begc:endc))
    allocate(total_mass_flux_lateral_col (begc:endc))
    allocate(total_mass_flux_seepage_col (begc:endc))
    allocate(qflx_seepage                (begc:endc))
    allocate(vsfm_mass_prev_col          (begc:endc,1:nlevgrnd))
    allocate(vsfm_dmass_col              (begc:endc))
    allocate(mass_beg_col                (begc:endc))
    allocate(mass_end_col                (begc:endc))

    allocate(mflx_et_col_1d              (nz_xylem))
    allocate(t_soil_col_1d               (nz_soil + nz_root + nz_xylem))

    allocate(mflx_drain_col_1d           ((endc-begc+1)*nlevgrnd))
    allocate(mflx_infl_col_1d            (endc-begc+1))
    allocate(mflx_dew_col_1d             (endc-begc+1))
    allocate(mflx_sub_snow_col_1d        (endc-begc+1))
    allocate(mflx_snowlyr_col_1d         (endc-begc+1))

    allocate(vsfm_mass_col_1d            (nz_soil+nz_root+nz_xylem))
    allocate(vsfm_fliq_col_1d            (nz_soil+nz_root+nz_xylem))
    allocate(vsfm_smpl_col_1d            (nz_soil+nz_root+nz_xylem))
    allocate(vsfm_soilp_col_1d           (nz_soil+nz_root+nz_xylem))
    allocate(vsfm_sat_col_1d             (nz_soil+nz_root+nz_xylem))

    !allocate(vsfm_sat_col_1d             ((endc-begc+1)*nlevgrnd))

    ! initialize

    mflx_et_col_1d(:)                = 0.d0
    mflx_infl_col_1d(:)              = 0.d0
    mflx_dew_col_1d(:)               = 0.d0
    mflx_drain_col_1d(:)             = 0.d0
    mflx_sub_snow_col_1d(:)          = 0.d0
    mflx_snowlyr_col_1d(:)           = 0.d0
    t_soil_col_1d(:)                 = 298.15d0

    mass_beg                         = 0.d0
    mass_end                         = 0.d0
    total_mass_flux                  = 0.d0
    total_mass_flux_et               = 0.d0
    total_mass_flux_infl             = 0.d0
    total_mass_flux_dew              = 0.d0
    total_mass_flux_drain            = 0.d0
    total_mass_flux_snowlyr          = 0.d0
    total_mass_flux_sub              = 0.d0
    total_mass_flux_lateral          = 0.d0

    mass_beg_col(:)                  = 0.d0
    mass_end_col(:)                  = 0.d0
    total_mass_flux_col(:)           = 0.d0
    total_mass_flux_et_col(:)        = 0.d0
    total_mass_flux_infl_col(:)      = 0.d0
    total_mass_flux_dew_col(:)       = 0.d0
    total_mass_flux_drain_col(:)     = 0.d0
    total_mass_flux_snowlyr_col(:)   = 0.d0
    total_mass_flux_sub_col(:)       = 0.d0
    total_mass_flux_lateral_col(:)   = 0.d0

    vsfm_mass_prev_col(:,:)          = 0.d0
    vsfm_dmass_col(:)                = 0.d0

    ! Get total mass
    soe_auxvar_id = 1;
    call this%vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL ,       &
         VAR_MASS        ,       &
         soe_auxvar_id   ,       &
         vsfm_mass_col_1d        &
         )

    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)

       et_integrated_over_root = 0.d0

       do j = 1, nlevgrnd

          idx = (c - begc)*nlevgrnd + j

          et_integrated_over_root      = et_integrated_over_root      + l2e_mflux_et(c,j)
          total_mass_flux_et           = total_mass_flux_et           + l2e_mflux_et(c,j)
          total_mass_flux_et_col(c)    = total_mass_flux_et_col(c)    + l2e_mflux_et(c,j)

          mflx_drain_col_1d(idx)       = l2e_mflux_drain(c,j)
          total_mass_flux_drain        = total_mass_flux_drain        + mflx_drain_col_1d(idx)
          total_mass_flux_drain_col(c) = total_mass_flux_drain_col(c) + mflx_drain_col_1d(idx)

          mass_beg                     = mass_beg                     + vsfm_mass_col_1d(idx)
          mass_beg_col(c)              = mass_beg_col(c)              + vsfm_mass_col_1d(idx)
          vsfm_mass_prev_col(c,j)      = vsfm_mass_col_1d(idx)
       end do

       do j = nz_soil + 1, nz_soil + nz_root + nz_xylem
          idx = (c-begc)*nlevgrnd + j
          mass_beg_col(c) = mass_beg_col(c) + vsfm_mass_col_1d(idx)
          mass_beg        = mass_beg        + vsfm_mass_col_1d(idx)
       enddo

       do j = 1, nz_xylem
          idx = (c - begc)*nz_xylem + j
          mflx_et_col_1d(idx) = et_integrated_over_root/nz_xylem
       enddo
       
       idx = c - begc+1

       mflx_dew_col_1d(idx)           = l2e_mflux_dew(c)
       mflx_infl_col_1d(idx)          = l2e_mflux_infil(c)
       mflx_snowlyr_col_1d(idx)       = l2e_mflux_snowlyr(c)
       mflx_sub_snow_col_1d(idx)      = l2e_mflux_sub_snow(c)

       total_mass_flux_dew            = total_mass_flux_dew            + mflx_dew_col_1d(idx)
       total_mass_flux_dew_col(c)     = total_mass_flux_dew_col(c)     + mflx_dew_col_1d(idx)

       total_mass_flux_infl           = total_mass_flux_infl           + mflx_infl_col_1d(idx)
       total_mass_flux_infl_col(c)    = total_mass_flux_infl_col(c)    + mflx_infl_col_1d(idx)

       total_mass_flux_snowlyr        = total_mass_flux_snowlyr        + mflx_snowlyr_col_1d(idx)
       total_mass_flux_snowlyr_col(c) = total_mass_flux_snowlyr_col(c) + mflx_snowlyr_col_1d(idx)

       total_mass_flux_sub            = total_mass_flux_sub            + mflx_sub_snow_col_1d(idx)
       total_mass_flux_sub_col(c)     = total_mass_flux_sub_col(c)     + mflx_sub_snow_col_1d(idx)

       total_mass_flux_col(c) = total_mass_flux_et_col(c)      + &
            total_mass_flux_infl_col(c)    + &
            total_mass_flux_dew_col(c)     + &
            total_mass_flux_drain_col(c)   + &
            total_mass_flux_snowlyr_col(c) + &
            total_mass_flux_sub_col(c)     + &
            total_mass_flux_lateral_col(c)
    end do
    total_mass_flux        = total_mass_flux_et        + &
         total_mass_flux_infl      + &
         total_mass_flux_dew       + &
         total_mass_flux_drain     + &
         total_mass_flux_snowlyr   + &
         total_mass_flux_sub       + &
         total_mass_flux_lateral

    ! Set temperature
    soe_auxvar_id = 1;
    call this%vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_INTERNAL ,      &
         VAR_TEMPERATURE ,      &
         soe_auxvar_id   ,      &
         t_soil_col_1d          &
         )
    ! Set Infiltration
    soe_auxvar_id = this%vsfm_cond_id_for_infil;
    call this%vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS           ,  &
         VAR_BC_SS_CONDITION ,  &
         soe_auxvar_id       ,  &
         mflx_infl_col_1d       &
         )
    ! Set ET
    soe_auxvar_id = this%vsfm_cond_id_for_et;
    call this%vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS           ,  &
         VAR_BC_SS_CONDITION ,  &
         soe_auxvar_id       ,  &
         mflx_et_col_1d         &
         )
    ! Set Dew
    soe_auxvar_id = this%vsfm_cond_id_for_dew;
    call this%vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS           ,  &
         VAR_BC_SS_CONDITION ,  &
         soe_auxvar_id       ,  &
         mflx_dew_col_1d        &
         )
    ! Set Drainage sink
    soe_auxvar_id = this%vsfm_cond_id_for_drainage;
    call this%vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS           ,  &
         VAR_BC_SS_CONDITION ,  &
         soe_auxvar_id       ,  &
         mflx_drain_col_1d      &
         )
    ! Set mass flux associated with disappearance of snow layer
    ! from last time step
    soe_auxvar_id = this%vsfm_cond_id_for_snow;
    call this%vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS           ,  &
         VAR_BC_SS_CONDITION ,  &
         soe_auxvar_id       ,  &
         mflx_snowlyr_col_1d    &
         )
    ! Set mass flux associated with sublimation of snow
    soe_auxvar_id = this%vsfm_cond_id_for_sublimation;
    call this%vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS            , &
         VAR_BC_SS_CONDITION  , &
         soe_auxvar_id        , &
         mflx_sub_snow_col_1d   &
         )

    frac_ice(:,:)       = 0.d0
    vsfm_fliq_col_1d(:) = 1.d0
    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)
       do j = 1, nlevgrnd

          frac_ice(c,j) = l2e_h2osoi_ice(c,j)/(l2e_h2osoi_liq(c,j) + l2e_h2osoi_ice(c,j))

          idx = (c - begc)*nlevgrnd + j
          vsfm_fliq_col_1d(idx) = 1._r8 - frac_ice(c,j)
       end do
    end do

    ! Set frac_liq
    soe_auxvar_id = 1;
    call this%vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_INTERNAL  , &
         VAR_FRAC_LIQ_SAT , &
         soe_auxvar_id    , &
         vsfm_fliq_col_1d   &
         )

    ! Preform Pre-StepDT operations
    call this%vsfm_mpp%sysofeqns%PreStepDT()

    ! Get default SNES settings
    call SNESGetTolerances(this%vsfm_mpp%sysofeqns%snes , &
         atol_default            , &
         rtol_default            , &
         stol_default            , &
         max_it_default          , &
         max_f_default           , &
         ierr                      &
         )
    CHKERRQ(ierr)

    stol = stol_default
    rtol = rtol_default

    !
    ! Solve the VSFM.
    !
    iter_count           = 0
    diverged_count       = 0
    mass_bal_err_count   = 0
    abs_mass_error_col   = 0.d0
    successful_step      = PETSC_FALSE

    do

       iter_count = iter_count + 1

       call SNESSetTolerances(this%vsfm_mpp%sysofeqns%snes , &
            atol_default            , &
            rtol                    , &
            stol                    , &
            max_it_default          , &
            max_f_default           , &
            ierr                      &
            );
       CHKERRQ(ierr)

       call this%vsfm_mpp%sysofeqns%StepDT(dtime, nstep, &
            converged, converged_reason, ierr); CHKERRQ(ierr)


       if (.not. converged) then

          ! VSFM solver did not converge, so let's try again with different
          ! solver settings.

          stol             = stol_alternate
          diverged_count   = diverged_count + 1
          successful_step  = PETSC_FALSE

          ! Reduce total run length time by the amount VSFM ran successfully
          ! with previous solver settings
          dtime = dtime - this%vsfm_mpp%sysofeqns%time

          if (diverged_count > 1) then
             ! Set frac_liq
             vsfm_fliq_col_1d(:) = 1.d0
             soe_auxvar_id = 1;

             call this%vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_INTERNAL  , &
                  VAR_FRAC_LIQ_SAT , &
                  soe_auxvar_id    , &
                  vsfm_fliq_col_1d   &
                  )
          end if
       else

          ! Solver converged, so let's copy data from VSFM model to
          ! CLM's data structure.

          ! Get Liquid saturation
          soe_auxvar_id = 1;
          call this%vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL , &
               VAR_LIQ_SAT     , &
               soe_auxvar_id   , &
               vsfm_sat_col_1d   &
               )

          ! Get total mass
          soe_auxvar_id = 1;
          call this%vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL  , &
               VAR_MASS         , &
               soe_auxvar_id    , &
               vsfm_mass_col_1d   &
               )

          ! Get liquid soil matrix potential
          soe_auxvar_id = 1;
          call this%vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL     , &
               VAR_SOIL_MATRIX_POT , &
               soe_auxvar_id       , &
               vsfm_smpl_col_1d      &
               )

          ! Get soil liquid pressure. This is the prognostic state of VSFM
          ! and needs to be saved in the restart file.
          soe_auxvar_id = 1;
          call this%vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL   , &
               VAR_PRESSURE      , &
               soe_auxvar_id     , &
               vsfm_soilp_col_1d   &
               )

          qflx_seepage(:) = 0._r8

          ! Put the data in CLM's data structure
          mass_end        = 0.d0
          area            = 1.d0 ! [m^2]

          do fc = 1, l2e_num_hydrologyc
             c = l2e_filter_hydrologyc(fc)

             !if (lateral_connectivity) then
             !  g    = col%gridCell(c)
             !   area = ldomain_lateral%ugrid%areaGrid_ghosted(g)
             !endif

             ! initialization
             jwt = -1

             ! Loops in decreasing j so WTD can be computed in the same loop
             do j = nlevgrnd, 1, -1
                idx = (c-begc)*nlevgrnd + j

                e2l_h2osoi_liq(c,j) = (1.d0 - frac_ice(c,j))*vsfm_mass_col_1d(idx)/area
                e2l_h2osoi_ice(c,j) = frac_ice(c,j)         *vsfm_mass_col_1d(idx)/area

                mass_end        = mass_end        + vsfm_mass_col_1d(idx)
                mass_end_col(c) = mass_end_col(c) + vsfm_mass_col_1d(idx)

                vsfm_dmass_col(c) = vsfm_dmass_col(c) + &
                     (vsfm_mass_col_1d(idx)-vsfm_mass_prev_col(c,j))

                e2l_smp(c,j)    = vsfm_smpl_col_1d(idx)*1000.0_r8      ! [m] --> [mm]

                if (jwt == -1) then
                   ! Find the first soil that is unsaturated
                   if (e2l_smp(c,j) < 0._r8) jwt = j
                end if

             end do

             do j = nz_soil + 1, nz_soil + nz_root + nz_xylem
                idx = (c-begc)*nlevgrnd + j
                mass_end_col(c) = mass_end_col(c) + vsfm_mass_col_1d(idx)
                mass_end        = mass_end        + vsfm_mass_col_1d(idx)
             enddo

             idx_c                   = (c-begc) + 1
             e2l_h2oroot_liq(idx_c,:)  = 0.d0
             e2l_h2oxylem_liq(idx_c,:) = 0.d0

             do j = nz_soil + 1, nz_soil + nz_root
                idx = (c-begc)*nlevgrnd + j
                e2l_h2oroot_liq(idx_c,j-nz_soil) = vsfm_mass_col_1d(idx)
             enddo

             do j = nz_soil + nz_root + 1, nz_soil + nz_root + nz_xylem
                idx = (c-begc)*nlevgrnd + j
                e2l_h2oxylem_liq(idx_c, j - nz_soil - nz_root) = vsfm_mass_col_1d(idx)
             enddo

             ! Find maximum water balance error over the column
             abs_mass_error_col = max(abs_mass_error_col,                     &
                  abs(mass_beg_col(c) - mass_end_col(c) + &
                  total_mass_flux_col(c)*dt))
             e2l_qrecharge(c) = 0._r8

             if (jwt == -1 .or. jwt == nlevgrnd) then
                ! Water table below or in the last layer
                e2l_wtd(c) = l2e_zi(c,nlevgrnd)
             else
                z_dn = (l2e_zi(c,jwt-1) + l2e_zi(c,jwt  ))/2._r8
                z_up = (l2e_zi(c,jwt ) + l2e_zi(c,jwt+1))/2._r8
                e2l_wtd(c) = (0._r8 - e2l_smp(c,jwt))/(e2l_smp(c,jwt) - e2l_smp(c,jwt+1))*(z_dn - z_up) + z_dn
             endif
          end do

          ! Save soil liquid pressure from VSFM for all (active+nonactive) cells.
          ! soilp_col is used for restarting VSFM.
          do c = begc, endc
             do j = 1, nlevgrnd
                idx = (c - begc)*nlevgrnd + j
                e2l_soilp(c,j) = vsfm_soilp_col_1d(idx)
             end do
          end do

          ! Save xylem matric potential [mm]
          do c = begc, endc
             do j = 1, 170
                idx = (c - begc)*nlevgrnd + j + nz_soil + nz_root
                e2l_xylemp(c,j) = vsfm_smpl_col_1d(idx)*1000.0_r8 ! [m] --> [mm]
             end do
          end do

          ! For the solution that did converge, is the mass error acceptable?

          if (abs_mass_error_col >= max_abs_mass_error_col) then

             ! For the solution that converged, the mass error
             ! is unacceptable. So let's try again with tighter
             ! solution tolerance (stol) for SNES.

             mass_bal_err_count  = mass_bal_err_count + 1

             if (converged_reason == SNES_CONVERGED_FNORM_RELATIVE) then
                rtol = rtol/10._r8
             else if (converged_reason == SNES_CONVERGED_SNORM_RELATIVE) then
                stol = stol/10._r8
             endif

             dtime               = dt
             successful_step     = PETSC_FALSE
             abs_mass_error_col  = 0._r8
             mass_end_col(:)     = 0._r8

             ! Perform Pre-StepDT operations
             call this%vsfm_mpp%sysofeqns%PreStepDT()

          else

             successful_step  = PETSC_TRUE

          endif

       endif

       if (successful_step) exit

       if (iter_count >= max_iter_count) then
          write(iulog,*)'In soilwater_vsfm: VSFM failed to converge after multiple attempts.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

    end do

#if 0
    ! Add seepage flux from VSFM to surface runoff
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)

       qflx_surf(c) = qflx_surf(c) + qflx_seepage(c)
    enddo
#endif

    call SNESSetTolerances(this%vsfm_mpp%sysofeqns%snes, atol_default, rtol_default, stol_default, &
         max_it_default, max_f_default, ierr); CHKERRQ(ierr)

    call this%vsfm_mpp%sysofeqns%PostStepDT()

#if VSFM_DEBUG
    write(iulog,*)'VSFM-DEBUG: nstep                      = ',get_nstep()
    write(iulog,*)'VSFM-DEBUG: dtime                      = ',dt
    write(iulog,*)'VSFM-DEBUG: change in mass between dt  = ',-(mass_beg - mass_end)
    write(iulog,*)'VSFM-DEBUG: change in mass due to flux = ',total_mass_flux*dt
    write(iulog,*)'VSFM-DEBUG: Error in mass conservation = ',mass_beg - mass_end + total_mass_flux*dt
    write(iulog,*)'VSFM-DEBUG: et_flux    * dtime         = ',total_mass_flux_et*dt
    write(iulog,*)'VSFM-DEBUG: infil_flux * dtime         = ',total_mass_flux_infl*dt
    write(iulog,*)'VSFM-DEBUG: dew_flux   * dtime         = ',total_mass_flux_dew*dt
    write(iulog,*)'VSFM-DEBUG: drain_flux * dtime         = ',total_mass_flux_drain*dt
    write(iulog,*)'VSFM-DEBUG: snow_flux  * dtime         = ',total_mass_flux_snowlyr*dt
    write(iulog,*)'VSFM-DEBUG: sub_flux   * dtime         = ',total_mass_flux_sub*dt
    write(iulog,*)'VSFM-DEBUG: lat_flux   * dtime         = ',total_mass_flux_lateral*dt
    write(iulog,*)'VSFM-DEBUG: total_mass_flux            = ',total_mass_flux!/flux_unit_conversion
    write(iulog,*)'VSFM-DEBUG: et_flux                    = ',total_mass_flux_et
    write(iulog,*)'VSFM-DEBUG: infil_flux                 = ',total_mass_flux_infl
    write(iulog,*)'VSFM-DEBUG: dew_flux                   = ',total_mass_flux_dew
    write(iulog,*)'VSFM-DEBUG: drain_flux                 = ',total_mass_flux_drain
    write(iulog,*)'VSFM-DEBUG: snow_flux                  = ',total_mass_flux_snowlyr
    write(iulog,*)'VSFM-DEBUG: sub_flux                   = ',total_mass_flux_sub
    write(iulog,*)''
#endif

    deallocate(frac_ice                    )
    deallocate(total_mass_flux_col         )
    deallocate(total_mass_flux_et_col      )
    deallocate(total_mass_flux_infl_col    )
    deallocate(total_mass_flux_dew_col     )
    deallocate(total_mass_flux_drain_col   )
    deallocate(total_mass_flux_snowlyr_col )
    deallocate(total_mass_flux_sub_col     )
    deallocate(total_mass_flux_lateral_col )
    deallocate(total_mass_flux_seepage_col )
    deallocate(qflx_seepage                )
    deallocate(vsfm_mass_prev_col          )
    deallocate(vsfm_dmass_col              )
    deallocate(mass_beg_col                )
    deallocate(mass_end_col                )

    deallocate(mflx_et_col_1d              )
    deallocate(mflx_drain_col_1d           )
    deallocate(mflx_infl_col_1d            )
    deallocate(mflx_dew_col_1d             )
    deallocate(mflx_sub_snow_col_1d        )
    deallocate(mflx_snowlyr_col_1d         )
    deallocate(t_soil_col_1d               )

    deallocate(vsfm_mass_col_1d            )
    deallocate(vsfm_fliq_col_1d            )
    deallocate(vsfm_smpl_col_1d            )
    deallocate(vsfm_soilp_col_1d           )
    deallocate(vsfm_sat_col_1d             )

  end subroutine EM_VSFM_SPAC_Solve_Soil_Hydro

#endif

end module ExternalModelVSFMSPACMod
