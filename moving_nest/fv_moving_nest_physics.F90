!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!***********************************************************************
!> @file
!! @brief Provides Moving Nest functionality for physics and surface variables
!! @author W. Ramstrom.  Collaboration with Bin Liu and Chunxi Zhang, EMC
!! @email William.Ramstrom@noaa.gov
! =======================================================================!


! =======================================================================!
!
! Notes
!
!------------------------------------------------------------------------
! Moving Nest Subroutine Naming Convention
!-----------------------------------------------------------------------
!
! mn_meta_* subroutines perform moving nest operations for FV3 metadata.
!               These routines will run only once per nest move.
!
! mn_var_*  subroutines perform moving nest operations for an individual FV3 variable.
!               These routines will run many times per nest move.
!
! mn_prog_* subroutines perform moving nest operations for the list of prognostic fields.
!               These routines will run only once per nest move.
!
! mn_phys_* subroutines perform moving nest operations for the list of physics fields.
!               These routines will run only once per nest move.
!
! =======================================================================!

module fv_moving_nest_physics_mod

  use block_control_mod,      only: block_control_type
  use mpp_mod,                only: mpp_pe, mpp_sync, mpp_sync_self, mpp_send, mpp_error, NOTE, FATAL
  use mpp_domains_mod,        only: mpp_update_domains, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,        only: mpp_define_nest_domains, mpp_shift_nest_domains, nest_domain_type, domain2d
  use mpp_domains_mod,        only: mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod,        only: mpp_get_F2C_index, mpp_update_nest_coarse
  use mpp_domains_mod,        only: NORTH, SOUTH, EAST, WEST, CORNER, CENTER
  use mpp_domains_mod,        only: NUPDATE, SUPDATE, EUPDATE, WUPDATE, DGRID_NE

  use GFS_typedefs,           only: GFS_data_type, GFS_control_type, kind_phys
  use GFS_init,               only: GFS_grid_populate

  use boundary_mod,           only: update_coarse_grid, update_coarse_grid_mpp
#ifdef OVERLOAD_R4
  use constantsR4_mod,        only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks, hlv
#else
  use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks, hlv
#endif
  use field_manager_mod,      only: MODEL_ATMOS
  use fv_arrays_mod,          only: fv_atmos_type, fv_nest_type, fv_grid_type, R_GRID
  use fv_moving_nest_types_mod,   only: fv_moving_nest_prog_type, fv_moving_nest_physics_type, mn_surface_grids, fv_moving_nest_type
  use fv_arrays_mod,          only: allocate_fv_nest_bc_type, deallocate_fv_nest_bc_type
  use fv_grid_tools_mod,      only: init_grid
  use fv_grid_utils_mod,      only: grid_utils_init, ptop_min, dist2side_latlon
  use fv_mapz_mod,            only: Lagrangian_to_Eulerian, moist_cv, compute_total_energy
  use fv_nesting_mod,         only: dealloc_nested_buffers
  use fv_nwp_nudge_mod,       only: do_adiabatic_init
  use init_hydro_mod,         only: p_var
  use tracer_manager_mod,     only: get_tracer_index, get_tracer_names
  use fv_moving_nest_utils_mod,  only: alloc_halo_buffer, grid_geometry, output_grid_to_nc
  use fv_moving_nest_utils_mod,  only: fill_nest_from_buffer, fill_nest_from_buffer_cell_center, fill_nest_from_buffer_nearest_neighbor
  use fv_moving_nest_utils_mod,  only: fill_nest_halos_from_parent, fill_grid_from_supergrid, fill_weight_grid
  use fv_moving_nest_utils_mod,  only: alloc_read_data
  use fv_moving_nest_utils_mod,  only: fill_nest_from_buffer_cell_center_masked
  use fv_moving_nest_utils_mod,  only: fill_nest_halos_from_parent_masked

  use fv_moving_nest_mod,     only: mn_var_fill_intern_nest_halos, mn_var_dump_to_netcdf, mn_var_shift_data, calc_nest_alignment
  use fv_moving_nest_types_mod, only: Moving_nest
  implicit none

#ifdef NO_QUAD_PRECISION
  ! 64-bit precision (kind=8)
  integer, parameter:: f_p = selected_real_kind(15)
#else
  ! Higher precision (kind=16) for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

#ifdef OVERLOAD_R4
  real, parameter:: kind_dyn_snan=x'FFBFFFFF'
#else
  real, parameter:: kind_dyn_snan=x'FFF7FFFFFFFFFFFF'
#endif

#ifdef CCPP_32BIT
  real(kind_phys), parameter:: kind_phys_snan=x'FFBFFFFF'
#else
  real(kind_phys), parameter:: kind_phys_snan=x'FFF7FFFFFFFFFFFF'
#endif

  logical :: debug_log = .false.
  logical :: move_physics = .true.       ! Always true, unless developer sets move_physics to .False. here for debugging.
  logical :: move_nsst = .true.          ! Value is reset in fv_moving_nest_main.F90 from namelist options

  ! FIXME: Put move_hafs_subset in a namelist

  integer, parameter, private :: interp_A_grid = 1        ! cell-centered A-grid
  integer, parameter, private :: interp_A_grid_lmask = 7  ! land mask, cell-centered A-grid

  ! Persistent variables to enable debug printing after range warnings.
  type (fv_atmos_type), pointer                 :: save_Atm_n
  type (block_control_type), pointer            :: save_Atm_block
  type(GFS_control_type), pointer               :: save_GFS_Control
  type(GFS_data_type), pointer                  :: save_GFS_Data(:)

  integer, parameter :: DO_FILL_NEST_HALOS_FROM_PARENT = 1
  integer, parameter :: DO_FILL_INTERN_NEST_HALOS = 2
  integer, parameter :: DO_SHIFT_DATA = 3
  integer, parameter :: DO_COPY_TO_BLOCK_ARRAYS = 4
  integer, parameter :: DO_COPY_FROM_BLOCK_ARRAYS = 5

  type movement_info
    integer :: action = -1 !< One of the five "DO_" parameters in this module
    type(fv_moving_nest_physics_type), pointer :: mn_phys => null() !< Physics variables with halos
    type(fv_atmos_type), pointer :: ChildGrid => null() !< Grid data for the nest domain
    type(nest_domain_type), pointer :: nest_domain => null() !< Nest domain for FMS
    type(domain2d), pointer :: domain_fine => null() !< The nest domain
    logical :: is_fine_pe = .false. !< Is this the nest in the parent-nest communication?
    integer :: delta_i_c = -1 !< Movement in X direction
    integer :: delta_j_c = -1 !< MOvement in Y direction
    integer :: x_refine = -1 !< Nest refinement ratio, X direction
    integer :: y_refine = -1 !< Nest refinement ratio, Y direction
    integer :: isd = -1 !< i starting index of allocated area with halo
    integer :: jsd = -1 !< j starting index of allocated area with halo
    integer, pointer :: ii(:) => null() !< Physics block i index in domain
    integer, pointer :: jj(:) => null() !< Physics block j index in domain
    real(kind=kind_phys), pointer :: slmsk(:) => null() !< Physics block sea-land-ice mask: sea=0, land=1, ice=2
  end type movement_info

  interface mover
     module procedure mover_r4_2d, mover_r8_2d, mover_phys_3d, mover_phys_4d, mover_int_2d
  end interface mover

  interface block_copy
     module procedure block_copy_phys_4D, block_copy_phys_3D, block_copy_r4_2D, block_copy_r8_2D, block_copy_sfcfsw, block_copy_sfcflw, block_copy_int_2d
  end interface block_copy

contains

  !>@brief The subroutine 'mn_phys_reset_sfc_props' sets the static surface parameters from the high-resolution input file data
  !>@details This subroutine relies on earlier code reading the data from files into the mn_static data structure
  !!  This subroutine does not yet handle ice points or frac_grid - fractional landfrac/oceanfrac values
  subroutine mn_phys_reset_sfc_props(Atm, n, mn_static, Atm_block, GFS_data, ioffset, joffset, refine)
    type(fv_atmos_type), intent(inout),allocatable   :: Atm(:)              !< Array of atmospheric data
    integer, intent(in)                              :: n                   !< Current grid number
    type(mn_surface_grids), intent(in)               :: mn_static           !< Static surface data
    type(block_control_type), intent(in)             :: Atm_block           !< Physics block layout
    type(GFS_data_type), intent(inout)               :: GFS_data(:)         !< Physics variable data
    integer, intent(in)                              :: ioffset, joffset    !< Current nest offset in i,j direction
    integer, intent(in)                              :: refine              !< Nest refinement ratio

    ! For iterating through physics/surface vector data
    integer                 :: nb, blen, ix, i_pe, j_pe, i_idx, j_idx
    real(kind=kind_phys)    :: phys_oro

    ! Setup local land sea mask grid for masked interpolations
    do i_pe = Atm(n)%bd%isd, Atm(n)%bd%ied
      do j_pe = Atm(n)%bd%jsd, Atm(n)%bd%jed
        i_idx = (ioffset-1)*refine + i_pe
        j_idx = (joffset-1)*refine + j_pe

        Moving_nest(n)%mn_phys%slmsk(i_pe, j_pe) = mn_static%ls_mask_grid(i_idx, j_idx)
      enddo
    enddo

    !  Reset the variables from the fix_sfc files
    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
        i_pe = Atm_block%index(nb)%ii(ix)
        j_pe = Atm_block%index(nb)%jj(ix)

        i_idx = (ioffset-1)*refine + i_pe
        j_idx = (joffset-1)*refine + j_pe

        ! Reset the land sea mask from the hires parent data
        GFS_data(nb)%Sfcprop%slmsk(ix) = mn_static%ls_mask_grid(i_idx, j_idx)

        !  IFD values are 0 for land, and 1 for oceans/lakes -- reverse of the land sea mask
        !  Land Sea Mask has values of 0 for oceans/lakes, 1 for land, 2 for sea ice
        !  TODO figure out what ifd should be for sea ice
        if (mn_static%ls_mask_grid(i_idx, j_idx) .eq. 1 ) then
          if (move_nsst) GFS_data(nb)%Sfcprop%ifd(ix) = 0         ! Land
          GFS_data(nb)%Sfcprop%oceanfrac(ix) = 0   ! Land -- TODO permit fractions
          GFS_data(nb)%Sfcprop%landfrac(ix) = 1    ! Land -- TODO permit fractions
        else
          if (move_nsst) GFS_data(nb)%Sfcprop%ifd(ix) = 1         ! Ocean
          GFS_data(nb)%Sfcprop%oceanfrac(ix) = 1   ! Ocean -- TODO permit fractions
          GFS_data(nb)%Sfcprop%landfrac(ix) = 0    ! Ocean -- TODO permit fractions
        endif

        GFS_data(nb)%Sfcprop%tg3(ix) = mn_static%deep_soil_temp_grid(i_idx, j_idx)

        ! Follow logic from FV3/io/fv3atm_sfc_io.F90
        ! TODO this will need to be more complicated if we support frac_grid
        !if (nint(mn_static%soil_type_grid(i_idx, j_idx)) == 14 .or. int(mn_static%soil_type_grid(i_idx, j_idx)+0.5) <= 0) then
        !if (nint(mn_static%soil_type_grid(i_idx, j_idx)) == 14 .or.

        !if ( (mn_static%ls_mask_grid(i_idx, j_idx) .eq. 1 .and. nint(mn_static%land_frac_grid(i_idx, j_idx)) == 0) .or. &
        !    mn_static%soil_type_grid(i_idx, j_idx) < 0.5) then
        if (mn_static%ls_mask_grid(i_idx, j_idx) .eq. 1 .and. nint(mn_static%land_frac_grid(i_idx, j_idx)) == 0 ) then
          ! Water soil type == lake, etc. -- override the other variables and make this water
          !!print '("mn_phys_reset_sfc_props LAKE SOIL npe=",I0," x,y=",I0,",",I0," lat=",F10.3," lon=",F10.3)', mpp_pe(), i_idx, j_idx, GFS_data(nb)%Grid%xlat_d(ix), GFS_data(nb)%Grid%xlon_d(ix)-360.0

          if (move_nsst) GFS_data(nb)%Sfcprop%ifd(ix) = 1         ! Ocean
          GFS_data(nb)%Sfcprop%oceanfrac(ix) = 1   ! Ocean -- TODO permit fractions
          GFS_data(nb)%Sfcprop%landfrac(ix) = 0    ! Ocean -- TODO permit fractions

          GFS_data(nb)%Sfcprop%stype(ix) = 0
          GFS_data(nb)%Sfcprop%slmsk(ix) = 0
        else
          GFS_data(nb)%Sfcprop%stype(ix) = nint(mn_static%soil_type_grid(i_idx, j_idx))
        endif

        !GFS_data(nb)%Sfcprop%vfrac(ix) = mn_static%veg_frac_grid(i_idx, j_idx)
        GFS_data(nb)%Sfcprop%vtype(ix) = nint(mn_static%veg_type_grid(i_idx, j_idx))
        GFS_data(nb)%Sfcprop%slope(ix) = nint(mn_static%slope_type_grid(i_idx, j_idx))
        GFS_data(nb)%Sfcprop%snoalb(ix) = mn_static%max_snow_alb_grid(i_idx, j_idx)

        GFS_data(nb)%Sfcprop%facsf(ix) = mn_static%facsf_grid(i_idx, j_idx)
        GFS_data(nb)%Sfcprop%facwf(ix) = mn_static%facwf_grid(i_idx, j_idx)

        GFS_data(nb)%Sfcprop%alvsf(ix) = mn_static%alvsf_grid(i_idx, j_idx)
        GFS_data(nb)%Sfcprop%alvwf(ix) = mn_static%alvwf_grid(i_idx, j_idx)
        GFS_data(nb)%Sfcprop%alnsf(ix) = mn_static%alnsf_grid(i_idx, j_idx)
        GFS_data(nb)%Sfcprop%alnwf(ix) = mn_static%alnwf_grid(i_idx, j_idx)

        ! Reset the orography in the physics arrays, using the smoothed values from above
        phys_oro =  Atm(n)%phis(i_pe, j_pe) / grav
        GFS_data(nb)%Sfcprop%oro(ix) = phys_oro
        GFS_data(nb)%Sfcprop%oro_uf(ix) = phys_oro

      enddo
    enddo

  end subroutine mn_phys_reset_sfc_props

  !>@brief The subroutine 'mn_phys_reset_phys_latlon' sets the lat/lons from the high-resolution input file data
  !>@details This subroutine sets lat/lons of the moved nest, then recalculates all the derived quantities (dx,dy,etc.)
  subroutine mn_reset_phys_latlon(Atm, n, tile_geo, fp_super_tile_geo, Atm_block, GFS_control, GFS_data)
    type(fv_atmos_type), allocatable, intent(in)      :: Atm(:)               !< Array of atmospheric data
    integer, intent(in)                  :: n                    !< Current grid number
    type(grid_geometry), intent(in)      :: tile_geo             !< Bounds of this grid
    type(grid_geometry), intent(in)      :: fp_super_tile_geo    !< Bounds of high-resolution parent grid
    type(block_control_type), intent(in) :: Atm_block            !< Physics block layout
    type(GFS_control_type), intent(in)   :: GFS_control          !< Physics metadata
    type(GFS_data_type), intent(inout)   :: GFS_data(:)          !< Physics variable data

    integer :: isc, jsc, iec, jec
    integer :: x, y, fp_i, fp_j
    integer :: nest_x, nest_y, parent_x, parent_y
    integer :: this_pe

    real(kind=kind_phys), allocatable :: lats(:,:), lons(:,:), area(:,:)

    this_pe = mpp_pe()

    isc = Atm(n)%bd%isc
    jsc = Atm(n)%bd%jsc
    iec = Atm(n)%bd%iec
    jec = Atm(n)%bd%jec

    allocate(lats(isc:iec, jsc:jec))
    allocate(lons(isc:iec, jsc:jec))
    allocate(area(isc:iec, jsc:jec))

    call calc_nest_alignment(Atm, n, nest_x, nest_y, parent_x, parent_y)

    do x = isc, iec
      do y = jsc, jec
        fp_i = (x - nest_x) * 2 + parent_x
        fp_j = (y - nest_y) * 2 + parent_y

        lons(x,y) = fp_super_tile_geo%lons(fp_i, fp_j)
        lats(x,y) = fp_super_tile_geo%lats(fp_i, fp_j)

        ! Need to add the areas from 4 squares, because the netCDF file has areas calculated for the supergrid cells
        !  We need the area of the whole center of the cell.
        !  Example dimensions for C288_grid.tile6.nc
        !   longitude -- x(577,577)
        !   latitude  -- y(577,577)
        !   area      -- x(576,576)

        !  Extracting lat/lon/area from Supergrid
        !
        !   1,1----2,1----3,1
        !    |      |      |
        !    | a1,1 | a2,1 |
        !    |      |      |
        !   1,2----2,2----3,2
        !    |      |      |
        !    | a1,2 | a2,2 |
        !    |      |      |
        !   1,3----2,3----3,3
        !
        !  The model A-grid cell 1,1 is centered at supergrid location 2,2
        !    The area of the A-grid cell is the sum of the 4 supergrid areas   A = a(1,1) + a(1,2) + a(2,1) + a(2,2)

        area(x,y) = fp_super_tile_geo%area(fp_i - 1, fp_j - 1) + fp_super_tile_geo%area(fp_i - 1, fp_j) + &
            fp_super_tile_geo%area(fp_i, fp_j - 1) + fp_super_tile_geo%area(fp_i, fp_j)   ! TODO make sure these offsets are correct.
      enddo
    enddo

    call GFS_grid_populate(GFS_data%Grid, lons, lats, area)

    deallocate(lats)
    deallocate(lons)
    deallocate(area)

  end subroutine mn_reset_phys_latlon

  !>@brief The subroutine 'mn_phys_fill_temp_variables' extracts 1D physics data into a 2D array for nest motion
  !>@details This subroutine fills in the mn_phys structure on the Atm object with 2D arrays of physics/surface variables.
  !!  Note that ice variables are not yet handled.
  subroutine mn_phys_fill_temp_variables(Atm, Atm_block, GFS_Control, GFS_Data, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)            !< Array of atmospheric data
    type (block_control_type), target, intent(inout)         :: Atm_block         !< Physics block layout
    type(GFS_control_type), target, intent(in)               :: GFS_Control       !< Physics metadata
    type(GFS_data_type), target, intent(inout)               :: GFS_Data(:)       !< Physics variable data
    integer, intent(in)                                      :: n, child_grid_num !< Current grid number, child grid number
    logical, intent(in)                                      :: is_fine_pe        !< Is this a nest PE?
    integer, intent(in)                                      :: npz               !< Number of vertical levels

    save_Atm_n => Atm(n)
    save_Atm_block => Atm_block
    save_GFS_Control => GFS_Control
    save_GFS_Data => GFS_Data

    call mn_phys_copy_all_blocks(Moving_nest(n)%mn_phys, .false., Atm, Atm_block, GFS_Control, GFS_Data, n)

  end subroutine mn_phys_fill_temp_variables

  pure subroutine block_copy_sfcfsw(mi, mn_phys, block_array)
    use module_radsw_parameters,  only: sfcfsw_type
    implicit none
    type(movement_info), intent(in) :: mi
    type(fv_moving_nest_physics_type), intent(inout) :: mn_phys !< Physics variables with halos
    type(sfcfsw_type), intent(inout) :: block_array(:)
    integer :: ix
    if(mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
       do ix = 1, size(block_array,1)
          block_array(ix)%upfxc = mn_phys%sfcfsw_upfxc(mi%ii(ix),mi%jj(ix))
          block_array(ix)%upfx0 = mn_phys%sfcfsw_upfx0(mi%ii(ix),mi%jj(ix))
          block_array(ix)%dnfxc = mn_phys%sfcfsw_dnfxc(mi%ii(ix),mi%jj(ix))
          block_array(ix)%dnfx0 = mn_phys%sfcfsw_dnfx0(mi%ii(ix),mi%jj(ix))
       enddo
    else
       do ix = 1, size(block_array,1)
          mn_phys%sfcfsw_upfxc(mi%ii(ix),mi%jj(ix)) = block_array(ix)%upfxc
          mn_phys%sfcfsw_upfx0(mi%ii(ix),mi%jj(ix)) = block_array(ix)%upfx0
          mn_phys%sfcfsw_dnfxc(mi%ii(ix),mi%jj(ix)) = block_array(ix)%dnfxc
          mn_phys%sfcfsw_dnfx0(mi%ii(ix),mi%jj(ix)) = block_array(ix)%dnfx0
       enddo
    endif
  end subroutine block_copy_sfcfsw

  pure subroutine block_copy_sfcflw(mi, mn_phys, block_array)
    use module_radlw_parameters,  only: sfcflw_type
    implicit none
    type(movement_info), intent(in) :: mi
    type(fv_moving_nest_physics_type), intent(inout) :: mn_phys !< Physics variables with halos
    type(sfcflw_type), intent(inout) :: block_array(:)
    integer :: ix
    if(mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
       do ix = 1, size(block_array,1)
          block_array(ix)%upfxc = mn_phys%sfcflw_upfxc(mi%ii(ix),mi%jj(ix))
          block_array(ix)%upfx0 = mn_phys%sfcflw_upfx0(mi%ii(ix),mi%jj(ix))
          block_array(ix)%dnfxc = mn_phys%sfcflw_dnfxc(mi%ii(ix),mi%jj(ix))
          block_array(ix)%dnfx0 = mn_phys%sfcflw_dnfx0(mi%ii(ix),mi%jj(ix))
       enddo
    else
       do ix = 1, size(block_array,1)
          mn_phys%sfcflw_upfxc(mi%ii(ix),mi%jj(ix)) = block_array(ix)%upfxc
          mn_phys%sfcflw_upfx0(mi%ii(ix),mi%jj(ix)) = block_array(ix)%upfx0
          mn_phys%sfcflw_dnfxc(mi%ii(ix),mi%jj(ix)) = block_array(ix)%dnfxc
          mn_phys%sfcflw_dnfx0(mi%ii(ix),mi%jj(ix)) = block_array(ix)%dnfx0
       enddo
    endif
  end subroutine block_copy_sfcflw

  pure subroutine block_copy_phys_4D(mi, work_array, block_array)
    implicit none
    type(movement_info), intent(in) :: mi
    real(kind_phys), intent(inout) :: work_array(mi%isd:,mi%jsd:,:,:)
    real(kind_phys), intent(inout) :: block_array(:,:,:)
    integer :: m, k, ix
    if(mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
       do m = lbound(block_array,3), ubound(block_array,3)
          do k = lbound(block_array,2), ubound(block_array,2)
             do ix = 1, size(block_array,1)
                block_array(ix,k,m) = work_array(mi%ii(ix),mi%jj(ix),k,m)
             enddo
          enddo
       enddo
    else
       do m = lbound(block_array,3), ubound(block_array,3)
          do k = lbound(block_array,2), ubound(block_array,2)
             do ix = 1, size(block_array,1)
                work_array(mi%ii(ix),mi%jj(ix),k,m) = block_array(ix,k,m)
             enddo
          enddo
       enddo
    endif
  end subroutine block_copy_phys_4D

  pure subroutine block_copy_phys_3D(mi, work_array, block_array)
    implicit none
    type(movement_info), intent(in) :: mi
    real(kind_phys), intent(inout) :: work_array(mi%isd:,mi%jsd:,:)
    real(kind_phys), intent(inout) :: block_array(:,:)
    integer :: k, ix
    if(mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
       do k = lbound(block_array,2), ubound(block_array,2)
          do ix = 1, size(block_array,1)
             block_array(ix,k) = work_array(mi%ii(ix),mi%jj(ix),k)
          enddo
       enddo
    else
       do k = lbound(block_array,2), ubound(block_array,2)
          do ix = 1, size(block_array,1)
             work_array(mi%ii(ix),mi%jj(ix),k) = block_array(ix,k)
          enddo
       enddo
    endif
  end subroutine block_copy_phys_3D

  pure subroutine block_copy_r4_2D(mi, work_array, block_array, &
       if_missing, if_missing_on_sea, if_missing_on_land, if_negative)
    implicit none
    type(movement_info), intent(in) :: mi
    real(kind=4), intent(inout) :: work_array(mi%isd:,mi%jsd:)
    real(kind=4), intent(inout) :: block_array(:)
    real(kind=4), optional, intent(in) :: if_missing, if_missing_on_sea, if_missing_on_land, if_negative
    integer :: ix
    if(mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
       do ix = 1, size(block_array,1)
          block_array(ix) = work_array(mi%ii(ix),mi%jj(ix))
       enddo

       ! Handle missing values. Generally, only one of these blocks will be reached,
       ! but the logic should work for any combination.
       ! Note that sea/land filters require slmsk or they'll crash when accessing it.

       if(present(if_missing)) then
          where(block_array>1e6)                          block_array = if_missing
       endif
       if(present(if_missing_on_sea)) then
          where(nint(mi%slmsk)==0 .and. block_array>1e6)  block_array = if_missing_on_sea
       endif
       if(present(if_missing_on_land)) then
          where(nint(mi%slmsk)==1 .and. block_array>1e6)  block_array = if_missing_on_land
       endif
       if(present(if_negative)) then
          where(block_array<0.0)                          block_array = if_negative
       endif
    else
       do ix = 1, size(block_array,1)
          work_array(mi%ii(ix),mi%jj(ix)) = block_array(ix)
       enddo
    endif
  end subroutine block_copy_r4_2D

  pure subroutine block_copy_r8_2D(mi, work_array, block_array, &
       if_missing, if_missing_on_sea, if_missing_on_land, if_negative)
    implicit none
    type(movement_info), intent(in) :: mi
    real(kind=8), intent(inout) :: work_array(mi%isd:,mi%jsd:)
    real(kind=8), intent(inout) :: block_array(:)
    real(kind=8), optional, intent(in) :: if_missing, if_missing_on_sea, if_missing_on_land, if_negative
    integer :: ix
    if(mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
       do ix = 1, size(block_array,1)
          block_array(ix) = work_array(mi%ii(ix),mi%jj(ix))
       enddo

       ! Handle missing values. Generally, only one of these blocks will be reached,
       ! but the logic should work for any combination.

       if(present(if_missing)) then
          where(block_array>1e6)                          block_array = if_missing
       endif
       if(present(if_missing_on_sea)) then
          where(nint(mi%slmsk)==0 .and. block_array>1e6)  block_array = if_missing_on_sea
       endif
       if(present(if_missing_on_land)) then
          where(nint(mi%slmsk)==1 .and. block_array>1e6)  block_array = if_missing_on_land
       endif
       if(present(if_negative)) then
          where(block_array<0.0)                          block_array = if_negative
       endif
    else
       do ix = 1, size(block_array,1)
          work_array(mi%ii(ix),mi%jj(ix)) = block_array(ix)
       enddo
    endif
  end subroutine block_copy_r8_2D

  pure subroutine block_copy_int_2D(mi, work_array, block_array, &
       if_missing, if_missing_on_sea, if_missing_on_land, if_negative)
    implicit none
    type(movement_info), intent(in) :: mi
    integer, intent(inout) :: work_array(mi%isd:,mi%jsd:)
    integer, intent(inout) :: block_array(:)
    integer, optional, intent(in) :: if_missing, if_missing_on_sea, if_missing_on_land, if_negative
    integer :: ix
    if(mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
       do ix = 1, size(block_array,1)
          block_array(ix) = work_array(mi%ii(ix),mi%jj(ix))
       enddo

       ! Handle missing values. Generally, only one of these blocks will be reached,
       ! but the logic should work for any combination.

       if(present(if_missing)) then
          where(block_array>1e6)                          block_array = if_missing
       endif
       if(present(if_missing_on_sea)) then
          where(nint(mi%slmsk)==0 .and. block_array>1e6)  block_array = if_missing_on_sea
       endif
       if(present(if_missing_on_land)) then
          where(nint(mi%slmsk)==1 .and. block_array>1e6)  block_array = if_missing_on_land
       endif
       if(present(if_negative)) then
          where(block_array<0.0)                          block_array = if_negative
       endif
    else
       do ix = 1, size(block_array,1)
          work_array(mi%ii(ix),mi%jj(ix)) = block_array(ix)
       enddo
    endif
  end subroutine block_copy_int_2D


  !>@brief The subroutine 'mn_phys_apply_temp_variables' copies moved 2D data back into 1D physics arryas for nest motion
  !>@details This subroutine fills the 1D physics arrays from the mn_phys structure on the Atm object
  !!  Note that ice variables are not yet handled.
  subroutine mn_phys_apply_temp_variables(Atm, Atm_block, GFS_Control, GFS_Data, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), target, intent(inout)               :: Atm(:)            !< Array of atmospheric data
    type (block_control_type), target, intent(inout)         :: Atm_block         !< Physics block layout
    type(GFS_control_type), intent(in)                       :: GFS_Control       !< Physics metadata
    type(GFS_data_type), intent(inout)                       :: GFS_Data(:)       !< Physics variable data
    integer, intent(in)                                      :: n, child_grid_num !< Current grid number, child grid number
    logical, intent(in)                                      :: is_fine_pe        !< Is this a nest PE?
    integer, intent(in)                                      :: npz               !< Number of vertical levels

    !  Needed to fill the local grids for parent and nest PEs in order to transmit/interpolate data from parent to nest
    !  But only the nest PE's have changed the values with nest motion, so they are the only ones that need to update the original arrays
    if (is_fine_pe) then
      call mn_phys_copy_all_blocks(Moving_nest(n)%mn_phys, .true., Atm, Atm_block, GFS_Control, GFS_Data, n)
    endif
  end subroutine mn_phys_apply_temp_variables

  subroutine mn_phys_copy_all_blocks(mn_phys, to_block, Atm, Atm_block, GFS_Control, GFS_Data, n)
    type(fv_moving_nest_physics_type), target, intent(inout) :: mn_phys
    logical, intent(in)                                      :: to_block
    type(fv_atmos_type), intent(inout)                       :: Atm(:)            !< Array of atmospheric data
    type (block_control_type), target, intent(inout)         :: Atm_block         !< Physics block layout
    type(GFS_control_type), intent(in)                       :: GFS_Control       !< Physics metadata
    type(GFS_data_type), intent(inout)                       :: GFS_Data(:)       !< Physics variable data
    integer, intent(in)                                      :: n                 !< Current grid number

    integer :: nb, blen, i, j, ix
    integer :: is, ie, js, je
    type(movement_info) :: mi

    is = Atm(n)%bd%is
    ie = Atm(n)%bd%ie
    js = Atm(n)%bd%js
    je = Atm(n)%bd%je

    ! SST directly in Atm structure
    Atm(n)%ts(is:ie, js:je) =  mn_phys%ts(is:ie, js:je)

    block_loop: do nb = 1, Atm_block%nblks
       if(to_block) then
          mi%action = DO_COPY_TO_BLOCK_ARRAYS
       else
          mi%action = DO_COPY_FROM_BLOCK_ARRAYS
       endif
       mi%isd = mn_phys%isd
       mi%jsd = mn_phys%jsd
       mi%ii => Atm_block%index(nb)%ii
       mi%jj => Atm_block%index(nb)%jj
       mi%slmsk => GFS_Data(nb)%Sfcprop%slmsk
       mi%mn_phys => mn_phys

       call mn_phys_mover(mi, GFS_Data(nb), GFS_Control)
    enddo block_loop
  end subroutine mn_phys_copy_all_blocks

  subroutine mn_phys_mover(mi, GFS_Data, GFS_Control, wt_h)
    implicit none
    type(movement_info), intent(inout)                       :: mi
    type(GFS_Control_type), intent(in)                       :: GFS_Control
    type(GFS_data_type), intent(inout)                       :: GFS_Data          !< Physics variable data for the current block
    real, allocatable, intent(in), optional :: wt_h(:,:,:)
    type(fv_moving_nest_physics_type), pointer :: mn_phys

    integer :: ix

    mn_phys => mi%mn_phys

    if_move_physics: if (move_physics) then

       lsm_choice: if (GFS_Control%lsm == GFS_Control%lsm_noah .or. GFS_Control%lsm == GFS_Control%lsm_noahmp) then
         call mover_phys_3d(mi, 'smc', mn_phys%smc, GFS_Data%Sfcprop%smc, wt_h=wt_h)
         call mover_phys_3d(mi, 'stc', mn_phys%stc, GFS_Data%Sfcprop%stc, wt_h=wt_h)
         call mover_phys_3d(mi, 'slc', mn_phys%slc, GFS_Data%Sfcprop%slc, wt_h=wt_h)
       elseif (GFS_Control%lsm == GFS_Control%lsm_ruc) then
         call mover_phys_3d(mi, 'sh2o', mn_phys%sh2o, GFS_Data%Sfcprop%sh2o, wt_h=wt_h)
         call mover_phys_3d(mi, 'smois', mn_phys%smois, GFS_Data%Sfcprop%smois, wt_h=wt_h)
         call mover_phys_3d(mi, 'tslb', mn_phys%tslb, GFS_Data%Sfcprop%tslb, wt_h=wt_h)
       endif lsm_choice

       ! EMIS PATCH - When copying back to GFS_Data, force to positive at all locations.
       call mover_r8_2d(mi, 'emis_lnd', mn_phys%emis_lnd, GFS_Data%Sfcprop%emis_lnd, wt_h=wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover_r8_2d(mi, 'emis_ice', mn_phys%emis_ice, GFS_Data%Sfcprop%emis_ice, wt_h=wt_h, &
            halo_ice_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover_r8_2d(mi, 'emis_wat', mn_phys%emis_wat, GFS_Data%Sfcprop%emis_wat, wt_h=wt_h, &
            halo_sea_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)

       ! When copying back to GFS_Data, set albedo values to physically reasonable values if they have negative fill values.
       call mover_r8_2d(mi, 'albdirvis_lnd ', mn_phys%albdirvis_lnd , GFS_Data%Sfcprop%albdirvis_lnd, wt_h=wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover_r8_2d(mi, 'albdirnir_lnd ', mn_phys%albdirnir_lnd , GFS_Data%Sfcprop%albdirnir_lnd, wt_h=wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover_r8_2d(mi, 'albdifvis_lnd ', mn_phys%albdifvis_lnd , GFS_Data%Sfcprop%albdifvis_lnd, wt_h=wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover_r8_2d(mi, 'albdifnir_lnd ', mn_phys%albdifnir_lnd , GFS_Data%Sfcprop%albdifnir_lnd, wt_h=wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)

       !call mover_r8_2d(mi, 'sfalb_lnd', mn_phys%sfalb_lnd, GFS_Data%Sfcprop%sfalb_lnd, wt_h=wt_h)
       !call mover_r8_2d(mi, 'sfalb_lnd_bck', mn_phys%sfalb_lnd_bck, GFS_Data%Sfcprop%sfalb_lnd_bck, wt_h=wt_h)
       !call mover_r8_2d(mi, 'semis', mn_phys%semis, GFS_Data%Radtend%semis, wt_h=wt_h)
       !call mover_r8_2d(mi, 'semisbase', mn_phys%semisbase, GFS_Data%Sfcprop%semisbase, wt_h=wt_h)
       !call mover_r8_2d(mi, 'sfalb', mn_phys%sfalb, GFS_Data%Radtend%sfalb, wt_h=wt_h)

       call mover_r8_2d(mi, 'u10m', mn_phys%u10m, GFS_Data%IntDiag%u10m, wt_h=wt_h)
       call mover_r8_2d(mi, 'v10m', mn_phys%v10m, GFS_Data%IntDiag%v10m, wt_h=wt_h)
       call mover_r8_2d(mi, 'tprcp', mn_phys%tprcp, GFS_Data%Sfcprop%tprcp, wt_h=wt_h)

       call mover_phys_3d(mi, 'hprime', mn_phys%hprime, GFS_Data%Sfcprop%hprime, wt_h=wt_h)

       call mover_r8_2d(mi, 'lakefrac', mn_phys%lakefrac, GFS_Data%Sfcprop%lakefrac, wt_h=wt_h)
       call mover_r8_2d(mi, 'lakedepth', mn_phys%lakedepth, GFS_Data%Sfcprop%lakedepth, wt_h=wt_h)

       call mover_r8_2d(mi, 'canopy', mn_phys%canopy, GFS_Data%Sfcprop%canopy, wt_h=wt_h)
       call mover_r8_2d(mi, 'vegfrac', mn_phys%vegfrac, GFS_Data%Sfcprop%vfrac, wt_h=wt_h)
       call mover_r8_2d(mi, 'uustar', mn_phys%uustar, GFS_Data%Sfcprop%uustar, wt_h=wt_h)
       call mover_r8_2d(mi, 'shdmin', mn_phys%shdmin, GFS_Data%Sfcprop%shdmin, wt_h=wt_h)
       call mover_r8_2d(mi, 'shdmax', mn_phys%shdmax, GFS_Data%Sfcprop%shdmax, wt_h=wt_h)

       ! When copying back to GFS_Data, set roughness lengths to physically reasonable values if they have
       ! fill value (possible at coastline) sea/land mask array (sea:0,land:1,sea-ice:2)
       call mover_r8_2d(mi, 'zorll', mn_phys%zorll, GFS_Data%Sfcprop%zorll, wt_h=wt_h, &
            halo_land_mask_fill=86.0_kind_phys, if_missing_on_land=82.0_kind_phys)
       call mover_r8_2d(mi, 'zorlw', mn_phys%zorlw, GFS_Data%Sfcprop%zorlw, wt_h=wt_h, &
            halo_land_mask_fill=78.0_kind_phys, if_missing_on_sea=83.0_kind_phys)
       call mover_r8_2d(mi, 'zorlwav', mn_phys%zorlwav, GFS_Data%Sfcprop%zorlwav, wt_h=wt_h, &
            halo_land_mask_fill=77.0_kind_phys, if_missing_on_sea=84.0_kind_phys)
       call mover_r8_2d(mi, 'usfco', mn_phys%usfco, GFS_Data%Sfcprop%usfco, wt_h=wt_h, &
            halo_land_mask_fill=0.0_kind_phys, if_missing_on_sea=0.0_kind_phys)
       call mover_r8_2d(mi, 'vsfco', mn_phys%vsfco, GFS_Data%Sfcprop%vsfco, wt_h=wt_h, &
            halo_land_mask_fill=0.0_kind_phys, if_missing_on_sea=0.0_kind_phys)

       call mover_r8_2d(mi, 'zorl', mn_phys%zorl, GFS_Data%Sfcprop%zorl, wt_h=wt_h, &
            if_missing=85.0_kind_phys)

       call mover_r8_2d(mi, 'tsfco', mn_phys%tsfco, GFS_Data%Sfcprop%tsfco, wt_h=wt_h)
       call mover_r8_2d(mi, 'tsfcl', mn_phys%tsfcl, GFS_Data%Sfcprop%tsfcl, wt_h=wt_h)
       call mover_r8_2d(mi, 'tsfc', mn_phys%tsfc, GFS_Data%Sfcprop%tsfc, wt_h=wt_h)

       call mover_phys_3d(mi, 'phy_f2d', mn_phys%phy_f2d, GFS_Data%Tbd%phy_f2d, wt_h=wt_h)
       call mover_phys_4d(mi, 'phy_f3d', mn_phys%phy_f3d, GFS_Data%Tbd%phy_f3d, wt_h=wt_h)

       call mover_r8_2d(mi, 'cv', mn_phys%cv, GFS_Data%Cldprop%cv, wt_h=wt_h)
       call mover_r8_2d(mi, 'cvt', mn_phys%cvt, GFS_Data%Cldprop%cvt, wt_h=wt_h)
       call mover_r8_2d(mi, 'cvb', mn_phys%cvb, GFS_Data%Cldprop%cvb, wt_h=wt_h)


       ! --------------------------------------------------------------------------------
       ! New stuff

       if(.false.) then
          ! FIXME: This should be in mv_static
          call mover_int_2d(mi, 'scolor', mn_phys%scolor, GFS_Data%Sfcprop%scolor, wt_h)
       endif


      ! Flake
      if (GFS_Control%lkm > 0 .and. GFS_Control%iopt_lake==GFS_Control%iopt_lake_flake) then
        call mover(mi, 'T_snow', mn_phys%T_snow, GFS_Data%Sfcprop%T_snow, wt_h=wt_h)
        call mover(mi, 'T_ice', mn_phys%T_ice, GFS_Data%Sfcprop%T_ice, wt_h=wt_h)
        call mover(mi, 'h_ML', mn_phys%h_ML, GFS_Data%Sfcprop%h_ML, wt_h=wt_h)
        call mover(mi, 't_ML', mn_phys%t_ML, GFS_Data%Sfcprop%t_ML, wt_h=wt_h)
        call mover(mi, 't_mnw', mn_phys%t_mnw, GFS_Data%Sfcprop%t_mnw, wt_h=wt_h)
        call mover(mi, 'h_talb', mn_phys%h_talb, GFS_Data%Sfcprop%h_talb, wt_h=wt_h)
        call mover(mi, 't_talb', mn_phys%t_talb, GFS_Data%Sfcprop%t_talb, wt_h=wt_h)
        call mover(mi, 't_bot1', mn_phys%t_bot1, GFS_Data%Sfcprop%t_bot1, wt_h=wt_h)
        call mover(mi, 't_bot2', mn_phys%t_bot2, GFS_Data%Sfcprop%t_bot2, wt_h=wt_h)
        call mover(mi, 'c_t', mn_phys%c_t, GFS_Data%Sfcprop%c_t, wt_h=wt_h)
      endif

      if (GFS_Control%lsm == GFS_Control%lsm_ruc) then
        call mover(mi, 'wetness', mn_phys%wetness, GFS_Data%Sfcprop%wetness, wt_h)
        call mover(mi, 'clw_surf_land', mn_phys%clw_surf_land, GFS_Data%Sfcprop%clw_surf_land, wt_h)
        call mover(mi, 'clw_surf_ice', mn_phys%clw_surf_ice, GFS_Data%Sfcprop%clw_surf_ice, wt_h)
        call mover(mi, 'qwv_surf_land', mn_phys%qwv_surf_land, GFS_Data%Sfcprop%qwv_surf_land, wt_h)
        call mover(mi, 'qwv_surf_ice', mn_phys%qwv_surf_ice, GFS_Data%Sfcprop%qwv_surf_ice, wt_h)
        call mover(mi, 'tsnow_land', mn_phys%tsnow_land, GFS_Data%Sfcprop%tsnow_land, wt_h)
        call mover(mi, 'tsnow_ice', mn_phys%tsnow_ice, GFS_Data%Sfcprop%tsnow_ice, wt_h)
        call mover(mi, 'snowfallac_land', mn_phys%snowfallac_land, GFS_Data%Sfcprop%snowfallac_land, wt_h)
        call mover(mi, 'snowfallac_ice', mn_phys%snowfallac_ice, GFS_Data%Sfcprop%snowfallac_ice, wt_h)
        call mover(mi, 'sfalb_lnd', mn_phys%sfalb_lnd, GFS_Data%Sfcprop%sfalb_lnd, wt_h)
        call mover(mi, 'sfalb_lnd_bck', mn_phys%sfalb_lnd_bck, GFS_Data%Sfcprop%sfalb_lnd_bck, wt_h)
        call mover(mi, 'sfalb_ice', mn_phys%sfalb_ice, GFS_Data%Sfcprop%sfalb_ice, wt_h)
        if (GFS_Control%rdlai) then
          call mover(mi, 'xlaixy', mn_phys%xlaixy, GFS_Data%Sfcprop%xlaixy, wt_h)
        endif
      endif

      !--------------------------------------------------------------------------------
      ! Everything up to this point does not change the HAFS.
      !--------------------------------------------------------------------------------
      
      if (GFS_Control%use_cice_alb .or. GFS_Control%lsm == GFS_Control%lsm_ruc) then
       ! When copying back to GFS_Data, set albedo values to physically reasonable values if they have negative fill values.
         call mover_r8_2d(mi, 'albdirvis_ice', mn_phys%albdirvis_ice , GFS_Data%Sfcprop%albdirvis_ice, wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
         call mover_r8_2d(mi, 'albdirnir_ice', mn_phys%albdirnir_ice , GFS_Data%Sfcprop%albdirnir_ice, wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
         call mover_r8_2d(mi, 'albdifvis_ice', mn_phys%albdifvis_ice , GFS_Data%Sfcprop%albdifvis_ice, wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
         call mover_r8_2d(mi, 'albdifnir_ice', mn_phys%albdifnir_ice , GFS_Data%Sfcprop%albdifnir_ice, wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
      endif
      
      ! FIXME: These should be enabled.
      disable_the_new_stuff: if(.false.) then

         call mover(mi, 'pgr', mn_phys%pgr, GFS_Data%Statein%pgr, wt_h)
         call mover(mi, 'tisfc', mn_phys%tisfc, GFS_Data%Sfcprop%tisfc, wt_h)
         call mover(mi, 'snowd', mn_phys%snowd, GFS_Data%Sfcprop%snowd, wt_h)
         call mover(mi, 'fice', mn_phys%fice, GFS_Data%Sfcprop%fice, wt_h)
         call mover(mi, 'sncovr', mn_phys%sncovr, GFS_Data%Sfcprop%sncovr, wt_h)
         call mover(mi, 'hice', mn_phys%hice, GFS_Data%Sfcprop%hice, wt_h)
         call mover(mi, 'weasd', mn_phys%weasd, GFS_Data%Sfcprop%weasd, wt_h)
         call mover(mi, 'ffmm', mn_phys%ffmm, GFS_Data%Sfcprop%ffmm, wt_h)
         call mover(mi, 'ffhh', mn_phys%ffhh, GFS_Data%Sfcprop%ffhh, wt_h)
         call mover(mi, 'f10m', mn_phys%f10m, GFS_Data%Sfcprop%f10m, wt_h)
         call mover(mi, 'srflag', mn_phys%srflag, GFS_Data%Sfcprop%srflag, wt_h)
         call mover(mi, 't2m', mn_phys%t2m, GFS_Data%Sfcprop%t2m, wt_h=wt_h)
         call mover(mi, 'q2m', mn_phys%q2m, GFS_Data%Sfcprop%q2m, wt_h=wt_h)
         call mover(mi, 'nirbmdi', mn_phys%nirbmdi, GFS_Data%Coupling%nirbmdi, wt_h=wt_h)
         call mover(mi, 'nirdfdi', mn_phys%nirdfdi, GFS_Data%Coupling%nirdfdi, wt_h=wt_h)
         call mover(mi, 'visbmdi', mn_phys%visbmdi, GFS_Data%Coupling%visbmdi, wt_h=wt_h)
         call mover(mi, 'visdfdi', mn_phys%visdfdi, GFS_Data%Coupling%visdfdi, wt_h=wt_h)
         call mover(mi, 'nirbmui', mn_phys%nirbmui, GFS_Data%Coupling%nirbmui, wt_h=wt_h)
         call mover(mi, 'nirdfui', mn_phys%nirdfui, GFS_Data%Coupling%nirdfui, wt_h=wt_h)
         call mover(mi, 'visbmui', mn_phys%visbmui, GFS_Data%Coupling%visbmui, wt_h=wt_h)
         call mover(mi, 'visdfui', mn_phys%visdfui, GFS_Data%Coupling%visdfui, wt_h=wt_h)
         call mover(mi, 'sfcdsw', mn_phys%sfcdsw, GFS_Data%Coupling%sfcdsw, wt_h=wt_h)
         call mover(mi, 'sfcnsw', mn_phys%sfcnsw, GFS_Data%Coupling%sfcnsw, wt_h=wt_h)
         call mover(mi, 'sfcdlw', mn_phys%sfcdlw, GFS_Data%Coupling%sfcdlw, wt_h=wt_h)
         call mover(mi, 'sfalb', mn_phys%sfalb, GFS_Data%Radtend%sfalb, wt_h=wt_h)
         call mover(mi, 'coszen', mn_phys%coszen, GFS_Data%Radtend%coszen, wt_h=wt_h)
         call mover(mi, 'tsflw', mn_phys%tsflw, GFS_Data%Radtend%tsflw, wt_h=wt_h)
         call mover(mi, 'semis', mn_phys%semis, GFS_Data%Radtend%semis, wt_h=wt_h)
         call mover(mi, 'coszdg', mn_phys%coszdg, GFS_Data%Radtend%coszdg, wt_h=wt_h)

         call mover_phys_3d(mi, 'htrsw', mn_phys%htrsw, GFS_Data%Radtend%htrsw, wt_h=wt_h)
         call mover_phys_3d(mi, 'htrlw', mn_phys%htrlw, GFS_Data%Radtend%htrlw, wt_h=wt_h)
         call mover_phys_3d(mi, 'swhc', mn_phys%swhc, GFS_Data%Radtend%swhc, wt_h=wt_h)
         call mover_phys_3d(mi, 'lwhc', mn_phys%lwhc, GFS_Data%Radtend%lwhc, wt_h=wt_h)

         rrtmg_types: if(mi%action == DO_COPY_FROM_BLOCK_ARRAYS .or. mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
            ! Copy between arrays of derived type (with no halos) to arrays of reals (with halos).
            call block_copy(mi, mn_phys, GFS_Data%Radtend%sfcfsw)
            call block_copy(mi, mn_phys, GFS_Data%Radtend%sfcflw)
         else
            ! While shifting data, everything is already in intrinsic type arrays with halos.
            call mover(mi, 'sfcflw_upfxc', mn_phys%sfcflw_upfxc, wt_h=wt_h)
            call mover(mi, 'sfcflw_upfx0', mn_phys%sfcflw_upfx0, wt_h=wt_h)
            call mover(mi, 'sfcflw_dnfxc', mn_phys%sfcflw_dnfxc, wt_h=wt_h)
            call mover(mi, 'sfcflw_dnfx0', mn_phys%sfcflw_dnfx0, wt_h=wt_h)
            call mover(mi, 'sfcflw_upfxc', mn_phys%sfcfsw_upfxc, wt_h=wt_h)
            call mover(mi, 'sfcflw_upfx0', mn_phys%sfcfsw_upfx0, wt_h=wt_h)
            call mover(mi, 'sfcflw_dnfxc', mn_phys%sfcfsw_dnfxc, wt_h=wt_h)
            call mover(mi, 'sfcflw_dnfx0', mn_phys%sfcfsw_dnfx0, wt_h=wt_h)
         endif rrtmg_types

         call mover(mi, 'tiice', mn_phys%tiice, GFS_Data%Sfcprop%tiice, wt_h=wt_h)
         call mover(mi, 'sncovr_ice', mn_phys%sncovr_ice, GFS_Data%Sfcprop%sncovr_ice, wt_h=wt_h)

         lsm_choice_2: if (GFS_Control%lsm == GFS_Control%lsm_noahmp) then
            call mover(mi, 'snowxy', mn_phys%snowxy, GFS_Data%Sfcprop%snowxy, wt_h=wt_h)
            call mover(mi, 'tvxy', mn_phys%tvxy, GFS_Data%Sfcprop%tvxy, wt_h=wt_h)
            call mover(mi, 'tgxy', mn_phys%tgxy, GFS_Data%Sfcprop%tgxy, wt_h=wt_h)
            call mover(mi, 'canicexy', mn_phys%canicexy, GFS_Data%Sfcprop%canicexy, wt_h=wt_h)
            call mover(mi, 'canliqxy', mn_phys%canliqxy, GFS_Data%Sfcprop%canliqxy, wt_h=wt_h)
            call mover(mi, 'eahxy', mn_phys%eahxy, GFS_Data%Sfcprop%eahxy, wt_h=wt_h)
            call mover(mi, 'tahxy', mn_phys%tahxy, GFS_Data%Sfcprop%tahxy, wt_h=wt_h)
            call mover(mi, 'cmxy', mn_phys%cmxy, GFS_Data%Sfcprop%cmxy, wt_h=wt_h)
            call mover(mi, 'chxy', mn_phys%chxy, GFS_Data%Sfcprop%chxy, wt_h=wt_h)
            call mover(mi, 'fwetxy', mn_phys%fwetxy, GFS_Data%Sfcprop%fwetxy, wt_h=wt_h)
            call mover(mi, 'sneqvoxy', mn_phys%sneqvoxy, GFS_Data%Sfcprop%sneqvoxy, wt_h=wt_h)
            call mover(mi, 'alboldxy', mn_phys%alboldxy, GFS_Data%Sfcprop%alboldxy, wt_h=wt_h)
            call mover(mi, 'qsnowxy', mn_phys%qsnowxy, GFS_Data%Sfcprop%qsnowxy, wt_h=wt_h)
            call mover(mi, 'wslakexy', mn_phys%wslakexy, GFS_Data%Sfcprop%wslakexy, wt_h=wt_h)
            call mover(mi, 'zwtxy', mn_phys%zwtxy, GFS_Data%Sfcprop%zwtxy, wt_h=wt_h)
            call mover(mi, 'waxy', mn_phys%waxy, GFS_Data%Sfcprop%waxy, wt_h=wt_h)
            call mover(mi, 'wtxy', mn_phys%wtxy, GFS_Data%Sfcprop%wtxy, wt_h=wt_h)
            call mover(mi, 'lfmassxy', mn_phys%lfmassxy, GFS_Data%Sfcprop%lfmassxy, wt_h=wt_h)
            call mover(mi, 'rtmassxy', mn_phys%rtmassxy, GFS_Data%Sfcprop%rtmassxy, wt_h=wt_h)
            call mover(mi, 'stmassxy', mn_phys%stmassxy, GFS_Data%Sfcprop%stmassxy, wt_h=wt_h)
            call mover(mi, 'woodxy', mn_phys%woodxy, GFS_Data%Sfcprop%woodxy, wt_h=wt_h)
            call mover(mi, 'stblcpxy', mn_phys%stblcpxy, GFS_Data%Sfcprop%stblcpxy, wt_h=wt_h)
            call mover(mi, 'fastcpxy', mn_phys%fastcpxy, GFS_Data%Sfcprop%fastcpxy, wt_h=wt_h)
            call mover(mi, 'xsaixy', mn_phys%xsaixy, GFS_Data%Sfcprop%xsaixy, wt_h=wt_h)
            call mover(mi, 'xlaixy', mn_phys%xlaixy, GFS_Data%Sfcprop%xlaixy, wt_h=wt_h)
            call mover(mi, 'taussxy', mn_phys%taussxy, GFS_Data%Sfcprop%taussxy, wt_h=wt_h)
            call mover(mi, 'smcwtdxy', mn_phys%smcwtdxy, GFS_Data%Sfcprop%smcwtdxy, wt_h=wt_h)
            call mover(mi, 'deeprechxy', mn_phys%deeprechxy, GFS_Data%Sfcprop%deeprechxy, wt_h=wt_h)
            call mover(mi, 'rechxy', mn_phys%rechxy, GFS_Data%Sfcprop%rechxy, wt_h=wt_h)

            ! These five arrays use bizarre indexing.
            !do k=-2,0
            call mover(mi, 'snicexy', mn_phys%snicexy, GFS_Data%Sfcprop%snicexy, wt_h)
            !do k=-2,0
            call mover(mi, 'snliqxy', mn_phys%snliqxy, GFS_Data%Sfcprop%snliqxy, wt_h)
            !do k=-2,0
            call mover(mi, 'tsnoxy', mn_phys%tsnoxy, GFS_Data%Sfcprop%tsnoxy, wt_h)
            !do k=1,4
            call mover(mi, 'smoiseq', mn_phys%smoiseq, GFS_Data%Sfcprop%smoiseq, wt_h)
            !do k=-2,4
            call mover(mi, 'zsnsoxy', mn_phys%zsnsoxy, GFS_Data%Sfcprop%zsnsoxy, wt_h)
         elseif (GFS_Control%lsm == GFS_Control%lsm_ruc) then
            ! Already moved up
         endif lsm_choice_2

         if (GFS_Control%nctp > 0 .and. GFS_Control%cscnv) then
            call mover(mi, 'phy_fctd', mn_phys%phy_fctd, GFS_Data%Tbd%phy_fctd, wt_h=wt_h)
         endif

         ! Do we need to move Statein and Stateout variables?
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3dlevsp1, GFS_Data%Statein%phii)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3dlevsp1, GFS_Data%Statein%prsi)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3dlevsp1, GFS_Data%Statein%prsik)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Statein%phil)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Statein%prsl)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Statein%prslk)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Statein%ugrs)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Statein%vgrs)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Statein%vvl)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Statein%tgrs)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Stateout%gu0)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Stateout%gv0)
         ! call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Stateout%gt0)
         ! do k = 1,ntr
         !   call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Statein%qgrs(:,:,k))
         !   call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data%Stateout%gq0(:,:,k))
         ! enddo
      end if disable_the_new_stuff

      ! End of new stuff
      ! --------------------------------------------------------------------------------

       check_stype_vtype: if(mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
          ! Check if stype and vtype are properly set for land points.  Set to reasonable values if they have fill values.
          do ix = 1, size(GFS_Data%Sfcprop%slmsk)
             if ( (int(GFS_data%Sfcprop%slmsk(ix)) .eq. 1) )  then

                if (GFS_data%Sfcprop%vtype(ix) .lt. 0.5) then
                   GFS_data%Sfcprop%vtype(ix) = 7    ! Force to grassland
                endif

                if (GFS_data%Sfcprop%stype(ix) .lt. 0.5) then
                   GFS_data%Sfcprop%stype(ix) = 3    ! Force to sandy loam
                endif

                if (GFS_data%Sfcprop%vtype_save(ix) .lt. 0.5) then
                   GFS_data%Sfcprop%vtype_save(ix) = 7    ! Force to grassland
                endif
                if (GFS_data%Sfcprop%stype_save(ix) .lt. 0.5) then
                   GFS_data%Sfcprop%stype_save(ix) = 3    ! Force to sandy loam
                endif

             endif
          enddo
       endif check_stype_vtype
    endif if_move_physics

    if_move_nsst: if (move_nsst) then
       call mover_r8_2d(mi, 'tref', mn_phys%tref, GFS_Data%Sfcprop%tref, wt_h)
       call mover_r8_2d(mi, 'z_c', mn_phys%z_c, GFS_Data%Sfcprop%z_c, wt_h)
       call mover_r8_2d(mi, 'c_0', mn_phys%c_0, GFS_Data%Sfcprop%c_0, wt_h)
       call mover_r8_2d(mi, 'c_d', mn_phys%c_d, GFS_Data%Sfcprop%c_d, wt_h)
       call mover_r8_2d(mi, 'w_0', mn_phys%w_0, GFS_Data%Sfcprop%w_0, wt_h)
       call mover_r8_2d(mi, 'w_d', mn_phys%w_d, GFS_Data%Sfcprop%w_d, wt_h)
       call mover_r8_2d(mi, 'xt', mn_phys%xt, GFS_Data%Sfcprop%xt, wt_h)
       call mover_r8_2d(mi, 'xs', mn_phys%xs, GFS_Data%Sfcprop%xs, wt_h)
       call mover_r8_2d(mi, 'xu', mn_phys%xu, GFS_Data%Sfcprop%xu, wt_h)
       call mover_r8_2d(mi, 'xv', mn_phys%xv, GFS_Data%Sfcprop%xv, wt_h)
       call mover_r8_2d(mi, 'xz', mn_phys%xz, GFS_Data%Sfcprop%xz, wt_h)
       call mover_r8_2d(mi, 'zm', mn_phys%zm, GFS_Data%Sfcprop%zm, wt_h)
       call mover_r8_2d(mi, 'xtts', mn_phys%xtts, GFS_Data%Sfcprop%xtts, wt_h)
       call mover_r8_2d(mi, 'xzts', mn_phys%xzts, GFS_Data%Sfcprop%xzts, wt_h)
       call mover_r8_2d(mi, 'd_conv', mn_phys%d_conv, GFS_Data%Sfcprop%d_conv, wt_h)
       call mover_r8_2d(mi, 'dt_cool', mn_phys%dt_cool, GFS_Data%Sfcprop%dt_cool, wt_h)
       call mover_r8_2d(mi, 'qrain', mn_phys%qrain, GFS_Data%Sfcprop%qrain, wt_h)
    endif if_move_nsst
    !endif move_other_physics
  end subroutine mn_phys_mover

  subroutine mover_r8_2d(mi, name, var, block_array, wt_h, &
       halo_sea_mask_fill, halo_land_mask_fill, halo_ice_mask_fill, &
       if_missing, if_missing_on_sea, if_missing_on_land, if_negative)
    implicit none
    type(movement_info), intent(in) :: mi
    real(kind=8), optional, intent(inout) :: block_array(:)
    real(kind=8), optional, intent(in) :: halo_sea_mask_fill, halo_land_mask_fill, halo_ice_mask_fill
    real(kind=8), optional, intent(in) :: if_missing, if_missing_on_sea, if_missing_on_land, if_negative
    character(len=*), intent(in) :: name
    real(kind=8), allocatable, intent(inout) :: var(:,:)
    real, allocatable, intent(in), optional :: wt_h(:,:,:)

    integer :: halo_interp, halo_slmsk
    real(kind=8) :: halo_default

    halo_slmsk=-1
    halo_default=0

    if(mi%action == DO_FILL_NEST_HALOS_FROM_PARENT) then
       if(present(halo_sea_mask_fill)) then
          halo_slmsk = 0
          halo_default = halo_sea_mask_fill
       else if(present(halo_land_mask_fill)) then
          halo_slmsk = 1
          halo_default = halo_land_mask_fill
       else if(present(halo_ice_mask_fill)) then
          halo_slmsk = 2
          halo_default = halo_ice_mask_fill
       endif
       if(halo_slmsk >= 0) then
          call fill_nest_halos_from_parent_masked(name, var, interp_A_grid_lmask, mi%ChildGrid%neststruct%wt_h, & ! mover_r8_2d_maksed
               mi%ChildGrid%neststruct%ind_h, &
               mi%ChildGrid%neststruct%refinement, mi%ChildGrid%neststruct%refinement, &
               mi%is_fine_pe, mi%nest_domain, CENTER, mi%mn_phys%slmsk, halo_slmsk, halo_default)
       else
          call fill_nest_halos_from_parent(name, var, interp_A_grid, mi%ChildGrid%neststruct%wt_h, & ! mover_r8_2d
               mi%ChildGrid%neststruct%ind_h, &
               mi%ChildGrid%neststruct%refinement, mi%ChildGrid%neststruct%refinement, &
               mi%is_fine_pe, mi%nest_domain, CENTER)
       endif
    else if(mi%action == DO_FILL_INTERN_NEST_HALOS) then
       call mn_var_fill_intern_nest_halos(var, mi%domain_fine, mi%is_fine_pe)
    else if(mi%action == DO_SHIFT_DATA) then
       call mn_var_shift_data(var, 1, wt_h, mi%ChildGrid%neststruct%ind_h, & ! mover_r8_2d
            mi%delta_i_c, mi%delta_j_c, mi%x_refine, mi%y_refine, mi%is_fine_pe, mi%nest_domain, CENTER)
    else
       call block_copy(mi, var, block_array, if_missing, if_missing_on_sea, if_missing_on_land, if_negative)
    endif
  end subroutine mover_r8_2d

  subroutine mover_r4_2d(mi, name, var, block_array, wt_h, &
       if_missing, if_missing_on_sea, if_missing_on_land, if_negative)
    implicit none
    type(movement_info), intent(in) :: mi
    real(kind=4), optional, intent(inout) :: block_array(:)
    real(kind=4), optional, intent(in) :: if_missing, if_missing_on_sea, if_missing_on_land, if_negative
    character(len=*), intent(in) :: name
    real(kind=4), allocatable, intent(inout) :: var(:,:)
    real, allocatable, intent(in), optional :: wt_h(:,:,:)

    if(mi%action == DO_FILL_NEST_HALOS_FROM_PARENT) then
       call fill_nest_halos_from_parent(name, var, interp_A_grid, mi%ChildGrid%neststruct%wt_h, & ! mover_r8_2d
            mi%ChildGrid%neststruct%ind_h, &
            mi%ChildGrid%neststruct%refinement, mi%ChildGrid%neststruct%refinement, &
            mi%is_fine_pe, mi%nest_domain, CENTER)
    else if(mi%action == DO_FILL_INTERN_NEST_HALOS) then
       call mn_var_fill_intern_nest_halos(var, mi%domain_fine, mi%is_fine_pe)
    else if(mi%action == DO_SHIFT_DATA) then
       call mn_var_shift_data(var, 1, wt_h, mi%ChildGrid%neststruct%ind_h, & ! mover_r4_2d
            mi%delta_i_c, mi%delta_j_c, mi%x_refine, mi%y_refine, mi%is_fine_pe, mi%nest_domain, CENTER)
    else
       call block_copy(mi, var, block_array, if_missing, if_missing_on_sea, if_missing_on_land, if_negative)
    endif
  end subroutine mover_r4_2d

  subroutine mover_int_2d(mi, name, var, block_array, wt_h, &
       if_missing, if_missing_on_sea, if_missing_on_land, if_negative)
    implicit none
    type(movement_info), intent(in) :: mi
    integer, optional, intent(inout) :: block_array(:)
    integer, optional, intent(in) :: if_missing, if_missing_on_sea, if_missing_on_land, if_negative
    character(len=*), intent(in) :: name
    integer, allocatable, intent(inout) :: var(:,:)
    real, allocatable, intent(in), optional :: wt_h(:,:,:)

    if(mi%action == DO_FILL_NEST_HALOS_FROM_PARENT) then
       call fill_nest_halos_from_parent(name, var, interp_A_grid, mi%ChildGrid%neststruct%wt_h, & ! mover_r8_2d
            mi%ChildGrid%neststruct%ind_h, &
            mi%ChildGrid%neststruct%refinement, mi%ChildGrid%neststruct%refinement, &
            mi%is_fine_pe, mi%nest_domain, CENTER)
    else if(mi%action == DO_FILL_INTERN_NEST_HALOS) then
       call mn_var_fill_intern_nest_halos(var, mi%domain_fine, mi%is_fine_pe)
    else if(mi%action == DO_SHIFT_DATA) then
       call mn_var_shift_data(var, 1, wt_h, mi%ChildGrid%neststruct%ind_h, & ! mover_r4_2d
            mi%delta_i_c, mi%delta_j_c, mi%x_refine, mi%y_refine, mi%is_fine_pe, mi%nest_domain, CENTER)
    else
       call block_copy(mi, var, block_array, if_missing, if_missing_on_sea, if_missing_on_land, if_negative)
    endif
  end subroutine mover_int_2d

  subroutine mover_phys_3d(mi, name, var, block_array, wt_h)
    implicit none
    type(movement_info), intent(in) :: mi
    character(len=*), intent(in) :: name
    real(kind_phys), allocatable, intent(inout) :: var(:,:,:)
    real(kind_phys), intent(inout) :: block_array(:,:)
    real, allocatable, intent(in), optional :: wt_h(:,:,:)

    if(mi%action == DO_FILL_NEST_HALOS_FROM_PARENT) then
       call fill_nest_halos_from_parent(name, var, interp_A_grid, mi%ChildGrid%neststruct%wt_h, & ! mover_r8_3d
            mi%ChildGrid%neststruct%ind_h, &
            mi%ChildGrid%neststruct%refinement, mi%ChildGrid%neststruct%refinement, &
            mi%is_fine_pe, mi%nest_domain, CENTER, size(var,3))
    else if(mi%action == DO_FILL_INTERN_NEST_HALOS) then
       call mn_var_fill_intern_nest_halos(var, mi%domain_fine, mi%is_fine_pe)
    else if(mi%action == DO_SHIFT_DATA) then
       call mn_var_shift_data(var, 1, wt_h, mi%ChildGrid%neststruct%ind_h, & ! mover_phys_3d
            mi%delta_i_c, mi%delta_j_c, mi%x_refine, mi%y_refine, mi%is_fine_pe, mi%nest_domain, CENTER, size(var,3))
    else
       call block_copy(mi, var, block_array)
    endif
  end subroutine mover_phys_3d

  subroutine mover_phys_4d(mi, name, var, block_array, wt_h)
    implicit none
    type(movement_info), intent(in) :: mi
    real(kind_phys), intent(inout) :: block_array(:,:,:)
    character(len=*), intent(in) :: name
    real(kind=kind_phys), allocatable, intent(inout) :: var(:,:,:,:)
    real, allocatable, intent(in), optional :: wt_h(:,:,:)

    if(mi%action == DO_FILL_NEST_HALOS_FROM_PARENT) then
       call fill_nest_halos_from_parent(name, var, interp_A_grid, mi%ChildGrid%neststruct%wt_h, & ! mover_r4_4d
            mi%ChildGrid%neststruct%ind_h, &
            mi%ChildGrid%neststruct%refinement, mi%ChildGrid%neststruct%refinement, &
            mi%is_fine_pe, mi%nest_domain, CENTER, size(var,3))
    else if(mi%action == DO_FILL_INTERN_NEST_HALOS) then
       call mn_var_fill_intern_nest_halos(var, mi%domain_fine, mi%is_fine_pe)
    else if(mi%action == DO_SHIFT_DATA) then
       call mn_var_shift_data(var, 1, wt_h, mi%ChildGrid%neststruct%ind_h, & ! mover_phys_4d
            mi%delta_i_c, mi%delta_j_c, mi%x_refine, mi%y_refine, mi%is_fine_pe, mi%nest_domain, CENTER, size(var,3))
    else
       call block_copy(mi, var, block_array)
    endif
  end subroutine mover_phys_4d

  !>@brief The subroutine 'mn_physfill_nest_halos_from_parent' transfers data from the coarse grid to the nest edge
  !>@details This subroutine must run on parent and nest PEs to complete the data transfers
  subroutine mn_phys_fill_nest_halos_from_parent(Atm, GFS_Control, GFS_Data, mn_static, n, child_grid_num, is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)            !< Array of atmospheric data
    type(GFS_control_type), intent(in)                       :: GFS_Control       !< Physics metadata
    type(GFS_data_type), intent(inout)                       :: GFS_Data(:)       !< Physics variable data
    type(mn_surface_grids), intent(in)                       :: mn_static         !< Static data
    integer, intent(in)                                      :: n, child_grid_num !< Current grid number, child grid number
    logical, intent(in)                                      :: is_fine_pe        !< Is this a nest PE?
    type(nest_domain_type), target, intent(inout)            :: nest_domain       !< Nest domain for FMS
    integer, intent(in)                                      :: nz                !< Number of vertical levels

    type(movement_info) :: mi

    mi%action = DO_FILL_NEST_HALOS_FROM_PARENT
    mi%mn_phys => Moving_nest(n)%mn_phys
    mi%ChildGrid => Atm(child_grid_num)
    mi%nest_domain => nest_domain
    mi%domain_fine => null()
    mi%is_fine_pe = is_fine_pe
    
    !  Fill centered-grid variables

    call mn_phys_mover(mi, GFS_Data(1), GFS_Control)

  end subroutine mn_phys_fill_nest_halos_from_parent

  !>@brief The subroutine 'mn_phys_fill_intern_nest_halos' fills the intenal nest halos for the physics variables
  !>@details This subroutine is only called for the nest PEs.
  subroutine mn_phys_fill_intern_nest_halos(moving_nest, GFS_Control, GFS_Data, domain_fine, is_fine_pe)
    type(fv_moving_nest_type), target, intent(inout) :: moving_nest         !< Single instance of moving nest data
    type(GFS_control_type), intent(in)               :: GFS_Control         !< Physics metadata
    type(GFS_data_type), intent(inout)               :: GFS_Data(:)         !< Physics variable data
    type(domain2d), target, intent(inout)            :: domain_fine         !< Domain structure for this nest
    logical, intent(in)                              :: is_fine_pe          !< Is nest PE - should be True.  Argument is redundant.

    type(movement_info) :: mi

    mi%action = DO_FILL_INTERN_NEST_HALOS
    mi%mn_phys => moving_nest%mn_phys
    mi%ChildGrid => null()
    mi%nest_domain => null()
    mi%domain_fine => domain_fine
    mi%is_fine_pe = is_fine_pe

    call mn_phys_mover(mi, GFS_Data(1), GFS_Control)

  end subroutine mn_phys_fill_intern_nest_halos

  !>@brief The subroutine 'mn_phys_shift_data' shifts the variable in the nest, including interpolating at the leading edge
  !>@details This subroutine is called for the nest and parent PEs.
  subroutine mn_phys_shift_data(Atm, GFS_Control, GFS_Data, n, child_grid_num, wt_h, wt_u, wt_v, &
      delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)                  !< Array of atmospheric data
    type(GFS_control_type), intent(in)                       :: GFS_Control             !< Physics metadata
    type(GFS_data_type), intent(inout)                       :: GFS_Data(:)             !< Physics variable data
    integer, intent(in)                                      :: n, child_grid_num       !< Current grid number, child grid number
    real, allocatable, intent(in)                            :: wt_h(:,:,:), wt_u(:,:,:), wt_v(:,:,:) !< Interpolation weights
    integer, intent(in)                                      :: delta_i_c, delta_j_c    !< Nest motion in i,j direction
    integer, intent(in)                                      :: x_refine, y_refine      !< Nest refinement
    logical, intent(in)                                      :: is_fine_pe              !< Is this the nest PE?
    type(nest_domain_type), intent(inout), target            :: nest_domain             !< Nest domain structure
    integer, intent(in)                                      :: nz                      !< Number of vertical levels

    type(movement_info) :: mi

    mi%action = DO_SHIFT_DATA
    mi%mn_phys => Moving_nest(n)%mn_phys
    mi%ChildGrid => Atm(child_grid_num)
    mi%nest_domain => nest_domain
    mi%domain_fine => null()
    mi%is_fine_pe = is_fine_pe
    mi%delta_i_c = delta_i_c
    mi%delta_j_c = delta_j_c
    mi%x_refine = x_refine
    mi%y_refine = y_refine

    call mn_phys_mover(mi, GFS_Data(1), GFS_Control, wt_h=wt_h)

  end subroutine mn_phys_shift_data

  !>@brief The subroutine 'mn_phys_dump_to_netcdf' dumps physics variables to debugging netCDF files
  !>@details This subroutine is called for the nest and parent PEs.
  subroutine mn_phys_dump_to_netcdf(Atm, Atm_block, GFS_Control, GFS_Data, time_val, file_prefix, is_fine_pe, domain_coarse, domain_fine, nz)
    type(fv_atmos_type), intent(in)            :: Atm                           !< Single instance of atmospheric data
    type (block_control_type), intent(in)      :: Atm_block                     !< Physics block layout
    type(GFS_control_type), intent(in)         :: GFS_Control                   !< Physics metadata
    type(GFS_data_type), intent(in)            :: GFS_Data(:)                   !< Physics variable data
    integer, intent(in)                        :: time_val                      !< Timestep number for filename
    character(len=*), intent(in)               :: file_prefix                   !< Prefix for output netCDF filenames
    logical, intent(in)                        :: is_fine_pe                    !< Is this the nest PE?
    type(domain2d), intent(in)                 :: domain_coarse, domain_fine    !< Domain structures for parent and nest
    integer, intent(in)                        :: nz                            !< Number of vertical levels

    integer :: is, ie, js, je
    integer :: nb, blen, i, j, k, ix, nv
    integer :: this_pe

    integer            :: n_moist
    character(len=16)  :: out_var_name, phys_var_name
    integer            :: position = CENTER

    ! Coerce the double precision variables from physics into single precision for debugging netCDF output
    ! Does not affect values used in calculations.
    ! TODO do we want to dump these as double precision??
    real, allocatable :: smc_pr_local (:,:,:)  !< soil moisture content
    real, allocatable :: stc_pr_local (:,:,:)  !< soil temperature
    real, allocatable :: slc_pr_local (:,:,:)  !< soil liquid water content
    real, allocatable, dimension(:,:) :: sealand_pr_local, deep_soil_t_pr_local, soil_type_pr_local, veg_type_pr_local, slope_type_pr_local, max_snow_alb_pr_local
    real, allocatable, dimension(:,:) :: tsfco_pr_local, tsfcl_pr_local, tsfc_pr_local, vegfrac_pr_local
    real, allocatable, dimension(:,:) :: tref_pr_local, c_0_pr_local, xt_pr_local,  xu_pr_local,  xv_pr_local, ifd_pr_local
    real, allocatable, dimension(:,:) :: facsf_pr_local, facwf_pr_local
    real, allocatable, dimension(:,:) :: alvsf_pr_local, alvwf_pr_local, alnsf_pr_local, alnwf_pr_local
    real, allocatable, dimension(:,:) :: zorl_pr_local, zorll_pr_local, zorlw_pr_local, zorli_pr_local
    real, allocatable, dimension(:,:) :: usfco_pr_local, vsfco_pr_local
    real, allocatable :: phy_f2d_pr_local (:,:,:)
    real, allocatable :: phy_f3d_pr_local (:,:,:,:)
    real, allocatable :: lakefrac_pr_local (:,:)  !< lake fraction
    real, allocatable :: landfrac_pr_local (:,:)  !< land fraction
    real, allocatable :: emis_lnd_pr_local (:,:)  !< emissivity land

    this_pe = mpp_pe()

    !  Skin temp/SST
    call mn_var_dump_to_netcdf(Atm%ts, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SSTK")
    !  Terrain height == phis / grav
    call mn_var_dump_to_netcdf(Atm%phis / grav, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "orog")

    ! sgh and oro were only fully allocated if fv_land is True
    !      if false, oro is (1,1), and sgh is not allocated
    if ( Atm%flagstruct%fv_land ) then
      ! land frac --  called oro in fv_array.F90
      call mn_var_dump_to_netcdf(Atm%oro, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "LFRAC")
      ! terrain standard deviation --  called sgh in fv_array.F90
      call mn_var_dump_to_netcdf(Atm%sgh, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "STDDEV")
    endif

    is = Atm%bd%is
    ie = Atm%bd%ie
    js = Atm%bd%js
    je = Atm%bd%je

    ! Just allocate compute domain size here for outputs;  the nest moving code also has halos added, but we don't need them here.
    if (move_physics) then
      allocate ( smc_pr_local(is:ie, js:je, GFS_Control%lsoil) )
      allocate ( stc_pr_local(is:ie, js:je, GFS_Control%lsoil) )
      allocate ( slc_pr_local(is:ie, js:je, GFS_Control%lsoil) )
      allocate ( sealand_pr_local(is:ie, js:je) )
      allocate ( lakefrac_pr_local(is:ie, js:je) )
      allocate ( landfrac_pr_local(is:ie, js:je) )
      allocate ( emis_lnd_pr_local(is:ie, js:je) )
      allocate ( phy_f2d_pr_local(is:ie, js:je, GFS_Control%ntot2d) )
      allocate ( phy_f3d_pr_local(is:ie, js:je, GFS_Control%levs, GFS_Control%ntot3d) )
      allocate ( tsfco_pr_local(is:ie, js:je) )
      allocate ( tsfcl_pr_local(is:ie, js:je) )
      allocate ( tsfc_pr_local(is:ie, js:je) )
      allocate ( vegfrac_pr_local(is:ie, js:je) )
      allocate ( alvsf_pr_local(is:ie, js:je) )
      allocate ( alvwf_pr_local(is:ie, js:je) )
      allocate ( alnsf_pr_local(is:ie, js:je) )
      allocate ( alnwf_pr_local(is:ie, js:je) )
      allocate ( deep_soil_t_pr_local(is:ie, js:je) )
      allocate ( soil_type_pr_local(is:ie, js:je) )
      !allocate ( veg_frac_pr_local(is:ie, js:je) )
      allocate ( veg_type_pr_local(is:ie, js:je) )
      allocate ( slope_type_pr_local(is:ie, js:je) )
      allocate ( max_snow_alb_pr_local(is:ie, js:je) )
      allocate ( facsf_pr_local(is:ie, js:je) )
      allocate ( facwf_pr_local(is:ie, js:je) )
      allocate ( zorl_pr_local(is:ie, js:je) )
      allocate ( zorll_pr_local(is:ie, js:je) )
      allocate ( zorlw_pr_local(is:ie, js:je) )
      allocate ( zorli_pr_local(is:ie, js:je) )
      allocate ( usfco_pr_local(is:ie, js:je) )
      allocate ( vsfco_pr_local(is:ie, js:je) )
    endif

    if (move_nsst) then
      allocate ( tref_pr_local(is:ie, js:je) )
      allocate ( c_0_pr_local(is:ie, js:je) )
      allocate ( xt_pr_local(is:ie, js:je) )
      allocate ( xu_pr_local(is:ie, js:je) )
      allocate ( xv_pr_local(is:ie, js:je) )
      allocate ( ifd_pr_local(is:ie, js:je) )
    endif

    if (move_physics) then
      smc_pr_local = +99999.9
      stc_pr_local = +99999.9
      slc_pr_local = +99999.9
      sealand_pr_local = +99999.9
      lakefrac_pr_local = +99999.9
      landfrac_pr_local = +99999.9
      emis_lnd_pr_local = +99999.9
      phy_f2d_pr_local = +99999.9
      phy_f3d_pr_local = +99999.9
      tsfco_pr_local = +99999.9
      tsfcl_pr_local = +99999.9
      tsfc_pr_local = +99999.9
      vegfrac_pr_local = +99999.9
      alvsf_pr_local = +99999.9
      alvwf_pr_local = +99999.9
      alnsf_pr_local = +99999.9
      alnwf_pr_local = +99999.9
    endif
    if (move_nsst) then
      tref_pr_local = +99999.9
      c_0_pr_local = +99999.9
      xt_pr_local = +99999.9
      xu_pr_local = +99999.9
      xv_pr_local = +99999.9
      ifd_pr_local = +99999.9
    endif

    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
        i = Atm_block%index(nb)%ii(ix)
        j = Atm_block%index(nb)%jj(ix)

        if (move_physics) then
          do k = 1, GFS_Control%lsoil
            ! Use real() to lower the precision
            smc_pr_local(i,j,k) = real(GFS_Data(nb)%Sfcprop%smc(ix,k))
            stc_pr_local(i,j,k) = real(GFS_Data(nb)%Sfcprop%stc(ix,k))
            slc_pr_local(i,j,k) = real(GFS_Data(nb)%Sfcprop%slc(ix,k))
          enddo

          sealand_pr_local(i,j) = real(GFS_Data(nb)%Sfcprop%slmsk(ix))
          lakefrac_pr_local(i,j) = real(GFS_Data(nb)%Sfcprop%lakefrac(ix))
          landfrac_pr_local(i,j) = real(GFS_Data(nb)%Sfcprop%landfrac(ix))
          emis_lnd_pr_local(i,j) = real(GFS_Data(nb)%Sfcprop%emis_lnd(ix))
          deep_soil_t_pr_local(i, j) = GFS_data(nb)%Sfcprop%tg3(ix)
          soil_type_pr_local(i, j) = GFS_data(nb)%Sfcprop%stype(ix)
          !veg_frac_pr_local(i, j) = GFS_data(nb)%Sfcprop%vfrac(ix)
          veg_type_pr_local(i, j) = GFS_data(nb)%Sfcprop%vtype(ix)
          slope_type_pr_local(i, j) = GFS_data(nb)%Sfcprop%slope(ix)
          facsf_pr_local(i, j) = GFS_data(nb)%Sfcprop%facsf(ix)
          facwf_pr_local(i, j) = GFS_data(nb)%Sfcprop%facwf(ix)
          zorl_pr_local(i, j) = GFS_data(nb)%Sfcprop%zorl(ix)
          zorlw_pr_local(i, j) = GFS_data(nb)%Sfcprop%zorlw(ix)
          zorll_pr_local(i, j) = GFS_data(nb)%Sfcprop%zorll(ix)
          zorli_pr_local(i, j) = GFS_data(nb)%Sfcprop%zorli(ix)
          usfco_pr_local(i, j) = GFS_data(nb)%Sfcprop%usfco(ix)
          vsfco_pr_local(i, j) = GFS_data(nb)%Sfcprop%vsfco(ix)
          max_snow_alb_pr_local(i, j) = GFS_data(nb)%Sfcprop%snoalb(ix)
          tsfco_pr_local(i, j) = GFS_data(nb)%Sfcprop%tsfco(ix)
          tsfcl_pr_local(i, j) = GFS_data(nb)%Sfcprop%tsfcl(ix)
          tsfc_pr_local(i, j)  = GFS_data(nb)%Sfcprop%tsfc(ix)
          vegfrac_pr_local(i, j) = GFS_data(nb)%Sfcprop%vfrac(ix)
          alvsf_pr_local(i, j) = GFS_data(nb)%Sfcprop%alvsf(ix)
          alvwf_pr_local(i, j) = GFS_data(nb)%Sfcprop%alvwf(ix)
          alnsf_pr_local(i, j) = GFS_data(nb)%Sfcprop%alnsf(ix)
          alnwf_pr_local(i, j) = GFS_data(nb)%Sfcprop%alnwf(ix)

          do nv = 1, GFS_Control%ntot2d
            ! Use real() to lower the precision
            phy_f2d_pr_local(i,j,nv) = real(GFS_Data(nb)%Tbd%phy_f2d(ix, nv))
          enddo

          do k = 1, GFS_Control%levs
            do nv = 1, GFS_Control%ntot3d
              ! Use real() to lower the precision
              phy_f3d_pr_local(i,j,k,nv) = real(GFS_Data(nb)%Tbd%phy_f3d(ix, k, nv))
            enddo
          enddo
        endif

        if (move_nsst) then
          tref_pr_local(i,j) = GFS_data(nb)%Sfcprop%tref(ix)
          c_0_pr_local(i,j) = GFS_data(nb)%Sfcprop%c_0(ix)
          xt_pr_local(i,j) = GFS_data(nb)%Sfcprop%xt(ix)
          xu_pr_local(i,j) = GFS_data(nb)%Sfcprop%xu(ix)
          xv_pr_local(i,j) = GFS_data(nb)%Sfcprop%xv(ix)
          ifd_pr_local(i,j) = GFS_data(nb)%Sfcprop%ifd(ix)
        endif
      enddo
    enddo

    if (move_physics) then
      !call mn_var_dump_to_netcdf(stc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, GFS_Control%lsoil, time_val, Atm%global_tile, file_prefix, "SOILT")
      !call mn_var_dump_to_netcdf(smc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, GFS_Control%lsoil, time_val, Atm%global_tile, file_prefix, "SOILM")
      !call mn_var_dump_to_netcdf(slc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, GFS_Control%lsoil, time_val, Atm%global_tile, file_prefix, "SOILL")
      call mn_var_dump_to_netcdf(sealand_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "LMASK")
      call mn_var_dump_to_netcdf(lakefrac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "LAKEFRAC")
      call mn_var_dump_to_netcdf(landfrac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "LANDFRAC")
      call mn_var_dump_to_netcdf(emis_lnd_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "EMISLAND")
      call mn_var_dump_to_netcdf(deep_soil_t_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "DEEPSOIL")
      call mn_var_dump_to_netcdf(soil_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SOILTP")
      !call mn_var_dump_to_netcdf(veg_frac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "VEGFRAC")
      call mn_var_dump_to_netcdf(veg_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "VEGTYPE")
      call mn_var_dump_to_netcdf(slope_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SLOPE")
      call mn_var_dump_to_netcdf(max_snow_alb_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SNOWALB")
      call mn_var_dump_to_netcdf(tsfco_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "TSFCO")
      call mn_var_dump_to_netcdf(tsfcl_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "TSFCL")
      call mn_var_dump_to_netcdf(tsfc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "TSFC")
      call mn_var_dump_to_netcdf(vegfrac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "VEGFRAC")
      call mn_var_dump_to_netcdf(alvsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ALVSF")
      call mn_var_dump_to_netcdf(alvwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ALVWF")
      call mn_var_dump_to_netcdf(alnsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ALNSF")
      call mn_var_dump_to_netcdf(alnwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ALNWF")
      call mn_var_dump_to_netcdf(facsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "FACSF")
      call mn_var_dump_to_netcdf(facwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "FACWF")
      call mn_var_dump_to_netcdf(zorl_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ZORL")
      call mn_var_dump_to_netcdf(zorlw_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ZORLW")
      call mn_var_dump_to_netcdf(zorll_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ZORLL")
      call mn_var_dump_to_netcdf(zorli_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ZORLI")
      call mn_var_dump_to_netcdf(usfco_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SSU")
      call mn_var_dump_to_netcdf(vsfco_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SSV")

      do nv = 1, GFS_Control%ntot2d
        write (phys_var_name, "(A4,I0.3)")  'PH2D', nv
        !call mn_var_dump_to_netcdf(phy_f2d_pr_local(:,:,nv), is_fine_pe, domain_coarse, domain_fine, position, 1, &
        !    time_val, Atm%global_tile, file_prefix, phys_var_name)
      enddo

      do nv = 1, GFS_Control%ntot3d
        write (phys_var_name, "(A4,I0.3)")  'PH3D', nv
        !call mn_var_dump_to_netcdf(phy_f3d_pr_local(:,:,:,nv), is_fine_pe, domain_coarse, domain_fine, position, GFS_Control%levs, &
        !    time_val, Atm%global_tile, file_prefix, phys_var_name)
      enddo
    endif

    if (move_nsst) then
      call mn_var_dump_to_netcdf(tref_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "TREF")
      call mn_var_dump_to_netcdf(c_0_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "C_0")
      call mn_var_dump_to_netcdf(xt_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "XT")
      call mn_var_dump_to_netcdf(xu_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "XU")
      call mn_var_dump_to_netcdf(xv_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "XV")
      call mn_var_dump_to_netcdf(ifd_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "IFD")
    endif

    if (move_physics) then
      deallocate(smc_pr_local)
      deallocate(stc_pr_local)
      deallocate(slc_pr_local)
      deallocate(lakefrac_pr_local)
      deallocate(landfrac_pr_local)
      deallocate(emis_lnd_pr_local)
      deallocate(sealand_pr_local, deep_soil_t_pr_local, soil_type_pr_local, veg_type_pr_local, max_snow_alb_pr_local)
      deallocate(tsfco_pr_local, tsfcl_pr_local, tsfc_pr_local, vegfrac_pr_local)
      deallocate(alvsf_pr_local, alvwf_pr_local, alnsf_pr_local, alnwf_pr_local)
      deallocate(facsf_pr_local, facwf_pr_local)
      deallocate(zorl_pr_local, zorlw_pr_local, zorll_pr_local, zorli_pr_local)
      deallocate(usfco_pr_local, vsfco_pr_local)
      deallocate(phy_f2d_pr_local)
      deallocate(phy_f3d_pr_local)
    endif

    if (move_nsst) deallocate(tref_pr_local, c_0_pr_local, xt_pr_local,  xu_pr_local,  xv_pr_local, ifd_pr_local)

  end subroutine mn_phys_dump_to_netcdf

end module fv_moving_nest_physics_mod
