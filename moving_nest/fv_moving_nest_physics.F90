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

  logical, parameter, private :: omit_hafs_bugs = .true. ! Set to .false. to get original HAFS results which forgot to move some variables

  ! FIXME: Put move_hafs_subset in a namelist

  integer, parameter, private :: interp_A_grid = 1        ! cell-centered A-grid
  integer, parameter, private :: interp_A_grid_lmask = 7  ! land mask, cell-centered A-grid

  ! Persistent variables to enable debug printing after range warnings.
  type (fv_atmos_type), pointer                 :: save_Atm_n
  type (block_control_type), pointer            :: save_Atm_block
  type(GFS_control_type), pointer               :: save_GFS_Control
  type(GFS_data_type), pointer                  :: save_GFS_Data(:)

  integer, parameter, private :: DO_FILL_NEST_HALOS_FROM_PARENT = 1
  integer, parameter, private :: DO_FILL_INTERN_NEST_HALOS = 2
  integer, parameter, private :: DO_SHIFT_DATA = 3
  integer, parameter, private :: DO_COPY_TO_BLOCK_ARRAYS = 4
  integer, parameter, private :: DO_COPY_FROM_BLOCK_ARRAYS = 5

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
     module procedure mover_r4_2d, mover_r8_2d, mover_phys_3d, mover_phys_4d, mover_int_2d, mover_int_3d
  end interface mover

  interface block_copy
     module procedure block_copy_phys_4D, block_copy_phys_3D, block_copy_r4_2D, block_copy_r8_2D, block_copy_sfcfsw, block_copy_sfcflw, block_copy_int_2d, block_copy_int_3d
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
        GFS_data(nb)%Sfcprop%slope_save(ix) = nint(mn_static%slope_type_grid(i_idx, j_idx))
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

  pure subroutine block_copy_int_3D(mi, work_array, block_array, &
       if_missing, if_missing_on_sea, if_missing_on_land, if_negative)
    implicit none
    type(movement_info), intent(in) :: mi
    integer, intent(inout) :: work_array(mi%isd:,mi%jsd:,:)
    integer, intent(inout) :: block_array(:,:)
    integer, optional, intent(in) :: if_missing, if_missing_on_sea, if_missing_on_land, if_negative
    integer :: ix, k
    if(mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
       do k = lbound(block_array,2), ubound(block_array,2)
          do ix = 1, size(block_array,1)
             block_array(ix, k) = work_array(mi%ii(ix),mi%jj(ix),k)
          enddo
       enddo

       ! Handle missing values. Generally, only one of these blocks will be reached,
       ! but the logic should work for any combination.

       if(present(if_missing)) then
          where(block_array>1e6)
             block_array = if_missing
          end where
       endif
       if(present(if_missing_on_sea)) then
          do k = lbound(block_array,2), ubound(block_array,2)
             do ix = 1, size(block_array,1)
                if(nint(mi%slmsk(ix))==0 .and. block_array(ix,k)>1e6) then
                   block_array(ix,k) = if_missing_on_sea
                endif
             enddo
          enddo
       endif
       if(present(if_missing_on_land)) then
          do k = lbound(block_array,2), ubound(block_array,2)
             do ix = 1, size(block_array,1)
                if(nint(mi%slmsk(ix))==1 .and. block_array(ix,k)>1e6) then
                   block_array(ix,k) = if_missing_on_land
                endif
             enddo
          enddo
       endif
       if(present(if_negative)) then
          where(block_array<0.0)
             block_array = if_negative
          end where
       endif
    else
       do k = lbound(block_array,2), ubound(block_array,2)
          do ix = 1, size(block_array,1)
             work_array(mi%ii(ix),mi%jj(ix),k) = block_array(ix,k)
          enddo
       enddo
    endif
  end subroutine block_copy_int_3D


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

    integer :: ix, k

    mn_phys => mi%mn_phys

    if_move_physics: if (move_physics) then

       if (GFS_Control%lsm == GFS_Control%lsm_noah .or. GFS_Control%lsm == GFS_Control%lsm_noahmp .or. GFS_Control%lsm == GFS_Control%lsm_ruc) then
         call mover(mi, 'smc', mn_phys%smc, GFS_Data%Sfcprop%smc, wt_h=wt_h)
         call mover(mi, 'stc', mn_phys%stc, GFS_Data%Sfcprop%stc, wt_h=wt_h)
         call mover(mi, 'slc', mn_phys%slc, GFS_Data%Sfcprop%slc, wt_h=wt_h)
       endif
       if (GFS_Control%lsm == GFS_Control%lsm_ruc) then
         call mover(mi, 'sh2o', mn_phys%sh2o, GFS_Data%Sfcprop%sh2o, wt_h=wt_h)
         call mover(mi, 'smois', mn_phys%smois, GFS_Data%Sfcprop%smois, wt_h=wt_h)
         ! FIXME: BOUND TSLB
         call mover(mi, 'tslb', mn_phys%tslb, GFS_Data%Sfcprop%tslb, wt_h=wt_h)
         call mover(mi, 'keepsmfr', mn_phys%keepsmfr, GFS_Data%Sfcprop%keepsmfr, wt_h=wt_h)
         call mover(mi, 'flag_frsoil', mn_phys%flag_frsoil, GFS_Data%Sfcprop%flag_frsoil, wt_h=wt_h)
         call mover(mi, 'rhofr', mn_phys%rhofr, GFS_Data%Sfcprop%rhofr, wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
         call mover(mi, 'fire_heat_flux', mn_phys%fire_heat_flux, GFS_Data%Sfcprop%fire_heat_flux, wt_h=wt_h)
         call mover(mi, 'frac_grid_burned', mn_phys%frac_grid_burned, GFS_Data%Sfcprop%frac_grid_burned, wt_h=wt_h)
       endif

       if (GFS_Control%do_myjsfc .or. GFS_Control%do_myjpbl) then
          call mover(mi, 'z0base', mn_phys%z0base, GFS_Data%Sfcprop%z0base, wt_h=wt_h)
       endif

       ! EMIS PATCH - When copying back to GFS_Data, force to positive at all locations.
       call mover(mi, 'emis_lnd', mn_phys%emis_lnd, GFS_Data%Sfcprop%emis_lnd, wt_h=wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover(mi, 'emis_ice', mn_phys%emis_ice, GFS_Data%Sfcprop%emis_ice, wt_h=wt_h, &
            halo_ice_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover(mi, 'emis_wat', mn_phys%emis_wat, GFS_Data%Sfcprop%emis_wat, wt_h=wt_h, &
            halo_sea_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)

       ! When copying back to GFS_Data, set albedo values to physically reasonable values if they have negative fill values.
       call mover(mi, 'albdirvis_lnd ', mn_phys%albdirvis_lnd , GFS_Data%Sfcprop%albdirvis_lnd, wt_h=wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover(mi, 'albdirnir_lnd ', mn_phys%albdirnir_lnd , GFS_Data%Sfcprop%albdirnir_lnd, wt_h=wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover(mi, 'albdifvis_lnd ', mn_phys%albdifvis_lnd , GFS_Data%Sfcprop%albdifvis_lnd, wt_h=wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover(mi, 'albdifnir_lnd ', mn_phys%albdifnir_lnd , GFS_Data%Sfcprop%albdifnir_lnd, wt_h=wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)

       call mover(mi, 'semisbase', mn_phys%semisbase, GFS_Data%Sfcprop%semisbase, wt_h=wt_h, &
               halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
       call mover(mi, 'u10m', mn_phys%u10m, GFS_Data%IntDiag%u10m, wt_h=wt_h)
       call mover(mi, 'v10m', mn_phys%v10m, GFS_Data%IntDiag%v10m, wt_h=wt_h)
       call mover(mi, 'tprcp', mn_phys%tprcp, GFS_Data%Sfcprop%tprcp, wt_h=wt_h)

       call mover(mi, 'hprime', mn_phys%hprime, GFS_Data%Sfcprop%hprime, wt_h=wt_h)

       call mover(mi, 'lakefrac', mn_phys%lakefrac, GFS_Data%Sfcprop%lakefrac, wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
       call mover(mi, 'lakedepth', mn_phys%lakedepth, GFS_Data%Sfcprop%lakedepth, wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=10.0_kind_phys)

       call mover(mi, 'canopy', mn_phys%canopy, GFS_Data%Sfcprop%canopy, wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
       call mover(mi, 'vegfrac', mn_phys%vegfrac, GFS_Data%Sfcprop%vfrac, wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
       call mover(mi, 'uustar', mn_phys%uustar, GFS_Data%Sfcprop%uustar, wt_h=wt_h)
       call mover(mi, 'shdmin', mn_phys%shdmin, GFS_Data%Sfcprop%shdmin, wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
       call mover(mi, 'shdmax', mn_phys%shdmax, GFS_Data%Sfcprop%shdmax, wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)

       ! When copying back to GFS_Data, set roughness lengths to physically reasonable values if they have
       ! fill value (possible at coastline) sea/land mask array (sea:0,land:1,sea-ice:2)
       call mover(mi, 'zorll', mn_phys%zorll, GFS_Data%Sfcprop%zorll, wt_h=wt_h, &
            halo_land_mask_fill=86.0_kind_phys, if_missing_on_land=82.0_kind_phys)
       call mover(mi, 'zorlw', mn_phys%zorlw, GFS_Data%Sfcprop%zorlw, wt_h=wt_h, &
            halo_land_mask_fill=78.0_kind_phys, if_missing_on_sea=83.0_kind_phys)
       call mover(mi, 'zorlwav', mn_phys%zorlwav, GFS_Data%Sfcprop%zorlwav, wt_h=wt_h, &
            halo_land_mask_fill=77.0_kind_phys, if_missing_on_sea=84.0_kind_phys)
       call mover(mi, 'usfco', mn_phys%usfco, GFS_Data%Sfcprop%usfco, wt_h=wt_h, &
            halo_land_mask_fill=0.0_kind_phys, if_missing_on_sea=0.0_kind_phys)
       call mover(mi, 'vsfco', mn_phys%vsfco, GFS_Data%Sfcprop%vsfco, wt_h=wt_h, &
            halo_land_mask_fill=0.0_kind_phys, if_missing_on_sea=0.0_kind_phys)

       call mover(mi, 'zorl', mn_phys%zorl, GFS_Data%Sfcprop%zorl, wt_h=wt_h, &
            if_missing=85.0_kind_phys)

       call mover(mi, 'tsfco', mn_phys%tsfco, GFS_Data%Sfcprop%tsfco, wt_h=wt_h, &
               halo_sea_mask_fill=290.0_kind_phys, if_negative=290.0_kind_phys)
       call mover(mi, 'tsfcl', mn_phys%tsfcl, GFS_Data%Sfcprop%tsfcl, wt_h=wt_h, &
               halo_land_mask_fill=290.0_kind_phys, if_negative=290.0_kind_phys)
       call mover(mi, 'tsfc', mn_phys%tsfc, GFS_Data%Sfcprop%tsfc, wt_h=wt_h)

       call mover(mi, 'phy_f2d', mn_phys%phy_f2d, GFS_Data%Tbd%phy_f2d, wt_h=wt_h)
       call mover(mi, 'phy_f3d', mn_phys%phy_f3d, GFS_Data%Tbd%phy_f3d, wt_h=wt_h)

       call mover(mi, 'cv', mn_phys%cv, GFS_Data%Cldprop%cv, wt_h=wt_h)
       call mover(mi, 'cvt', mn_phys%cvt, GFS_Data%Cldprop%cvt, wt_h=wt_h)
       call mover(mi, 'cvb', mn_phys%cvb, GFS_Data%Cldprop%cvb, wt_h=wt_h)

       call mover(mi, 'rmol', mn_phys%rmol, GFS_Data%Sfcprop%rmol, wt_h=wt_h)

       if (GFS_Control%cpllnd .and. GFS_Control%cpllnd2atm) then
          call mover(mi, 'sncovr1_lnd', mn_phys%sncovr1_lnd, GFS_Data%Coupling%sncovr1_lnd, wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 'qsurf_lnd', mn_phys%qsurf_lnd  , GFS_Data%Coupling%qsurf_lnd  , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 'evap_lnd', mn_phys%evap_lnd   , GFS_Data%Coupling%evap_lnd   , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 'hflx_lnd', mn_phys%hflx_lnd   , GFS_Data%Coupling%hflx_lnd   , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 'ep_lnd', mn_phys%ep_lnd     , GFS_Data%Coupling%ep_lnd     , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 't2mmp_lnd', mn_phys%t2mmp_lnd  , GFS_Data%Coupling%t2mmp_lnd  , wt_h=wt_h, &
               halo_land_mask_fill=290.0_kind_phys, if_negative=290.0_kind_phys)
          call mover(mi, 'q2mp_lnd', mn_phys%q2mp_lnd   , GFS_Data%Coupling%q2mp_lnd   , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 'gflux_lnd', mn_phys%gflux_lnd  , GFS_Data%Coupling%gflux_lnd  , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 'runoff_lnd', mn_phys%runoff_lnd , GFS_Data%Coupling%runoff_lnd , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 'drain_lnd', mn_phys%drain_lnd  , GFS_Data%Coupling%drain_lnd  , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 'cmm_lnd', mn_phys%cmm_lnd    , GFS_Data%Coupling%cmm_lnd    , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 'chh_lnd', mn_phys%chh_lnd    , GFS_Data%Coupling%chh_lnd    , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
          call mover(mi, 'zvfun_lnd', mn_phys%zvfun_lnd  , GFS_Data%Coupling%zvfun_lnd  , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
       endif

       ! --------------------------------------------------------------------------------
       ! New stuff

       if(.false.) then
          ! FIXME: This should be in mv_static
          call mover(mi, 'scolor', mn_phys%scolor, GFS_Data%Sfcprop%scolor, wt_h)
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
        call mover(mi, 'wetness', mn_phys%wetness, GFS_Data%Sfcprop%wetness, wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
        call mover(mi, 'clw_surf_land', mn_phys%clw_surf_land, GFS_Data%Sfcprop%clw_surf_land, wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
        call mover(mi, 'clw_surf_ice', mn_phys%clw_surf_ice, GFS_Data%Sfcprop%clw_surf_ice, wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
        call mover(mi, 'qwv_surf_land', mn_phys%qwv_surf_land, GFS_Data%Sfcprop%qwv_surf_land, wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
        call mover(mi, 'qwv_surf_ice', mn_phys%qwv_surf_ice, GFS_Data%Sfcprop%qwv_surf_ice, wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
        call mover(mi, 'tsnow_land', mn_phys%tsnow_land, GFS_Data%Sfcprop%tsnow_land, wt_h, &
               halo_land_mask_fill=273.0_kind_phys, if_negative=273.0_kind_phys)
        call mover(mi, 'tsnow_ice', mn_phys%tsnow_ice, GFS_Data%Sfcprop%tsnow_ice, wt_h, &
               halo_land_mask_fill=273.0_kind_phys, if_negative=273.0_kind_phys)
        call mover(mi, 'snowfallac_land', mn_phys%snowfallac_land, GFS_Data%Sfcprop%snowfallac_land, wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
        call mover(mi, 'snowfallac_ice', mn_phys%snowfallac_ice, GFS_Data%Sfcprop%snowfallac_ice, wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
        call mover(mi, 'sfalb_lnd', mn_phys%sfalb_lnd, GFS_Data%Sfcprop%sfalb_lnd, wt_h, &
               halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
        call mover(mi, 'sfalb_lnd_bck', mn_phys%sfalb_lnd_bck, GFS_Data%Sfcprop%sfalb_lnd_bck, wt_h, &
               halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
        call mover(mi, 'sfalb_ice', mn_phys%sfalb_ice, GFS_Data%Sfcprop%sfalb_ice, wt_h, &
               halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
        if (GFS_Control%rdlai) then
          call mover(mi, 'xlaixy', mn_phys%xlaixy, GFS_Data%Sfcprop%xlaixy, wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
        endif
      endif
      
      if (GFS_Control%imfdeepcnv == GFS_Control%imfdeepcnv_gf .or. GFS_Control%imfdeepcnv == GFS_Control%imfdeepcnv_c3) then
         call mover(mi, 'maxupmf', mn_phys%maxupmf, GFS_Data%Sfcprop%maxupmf, wt_h=wt_h)
      endif

      if (GFS_Control%use_cice_alb .or. GFS_Control%lsm == GFS_Control%lsm_ruc) then
       ! When copying back to GFS_Data, set albedo values to physically reasonable values if they have negative fill values.
         call mover(mi, 'albdirvis_ice', mn_phys%albdirvis_ice , GFS_Data%Sfcprop%albdirvis_ice, wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
         call mover(mi, 'albdirnir_ice', mn_phys%albdirnir_ice , GFS_Data%Sfcprop%albdirnir_ice, wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
         call mover(mi, 'albdifvis_ice', mn_phys%albdifvis_ice , GFS_Data%Sfcprop%albdifvis_ice, wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
         call mover(mi, 'albdifnir_ice', mn_phys%albdifnir_ice , GFS_Data%Sfcprop%albdifnir_ice, wt_h, &
            halo_land_mask_fill=0.5_kind_phys, if_negative=0.5_kind_phys)
      endif

      lsm_choice_2: if (GFS_Control%lsm == GFS_Control%lsm_noahmp) then
         call mover(mi, 'rca', mn_phys%rca, GFS_Data%Sfcprop%rca, wt_h=wt_h)
         call mover(mi, 'alboldxy', mn_phys%alboldxy, GFS_Data%Sfcprop%alboldxy, wt_h=wt_h)
         call mover(mi, 'canicexy', mn_phys%canicexy, GFS_Data%Sfcprop%canicexy, wt_h=wt_h)
         call mover(mi, 'canliqxy', mn_phys%canliqxy, GFS_Data%Sfcprop%canliqxy, wt_h=wt_h)
         call mover(mi, 'chxy', mn_phys%chxy, GFS_Data%Sfcprop%chxy, wt_h=wt_h)
         call mover(mi, 'cmxy', mn_phys%cmxy, GFS_Data%Sfcprop%cmxy, wt_h=wt_h)
         call mover(mi, 'deeprechxy', mn_phys%deeprechxy, GFS_Data%Sfcprop%deeprechxy, wt_h=wt_h)
         call mover(mi, 'eahxy', mn_phys%eahxy, GFS_Data%Sfcprop%eahxy, wt_h=wt_h)
         call mover(mi, 'fastcpxy', mn_phys%fastcpxy, GFS_Data%Sfcprop%fastcpxy, wt_h=wt_h)
         call mover(mi, 'fwetxy', mn_phys%fwetxy, GFS_Data%Sfcprop%fwetxy, wt_h=wt_h)
         call mover(mi, 'lfmassxy', mn_phys%lfmassxy, GFS_Data%Sfcprop%lfmassxy, wt_h=wt_h)
         call mover(mi, 'qsnowxy', mn_phys%qsnowxy, GFS_Data%Sfcprop%qsnowxy, wt_h=wt_h)
         call mover(mi, 'rechxy', mn_phys%rechxy, GFS_Data%Sfcprop%rechxy, wt_h=wt_h)
         call mover(mi, 'rtmassxy', mn_phys%rtmassxy, GFS_Data%Sfcprop%rtmassxy, wt_h=wt_h)
         call mover(mi, 'smcwtdxy', mn_phys%smcwtdxy, GFS_Data%Sfcprop%smcwtdxy, wt_h=wt_h)
         call mover(mi, 'smoiseq', mn_phys%smoiseq, GFS_Data%Sfcprop%smoiseq, wt_h)
         call mover(mi, 'sneqvoxy', mn_phys%sneqvoxy, GFS_Data%Sfcprop%sneqvoxy, wt_h=wt_h)
         call mover(mi, 'snicexy', mn_phys%snicexy, GFS_Data%Sfcprop%snicexy, wt_h)
         call mover(mi, 'snliqxy', mn_phys%snliqxy, GFS_Data%Sfcprop%snliqxy, wt_h)
         call mover(mi, 'snowxy', mn_phys%snowxy, GFS_Data%Sfcprop%snowxy, wt_h=wt_h)
         call mover(mi, 'stblcpxy', mn_phys%stblcpxy, GFS_Data%Sfcprop%stblcpxy, wt_h=wt_h)
         call mover(mi, 'stmassxy', mn_phys%stmassxy, GFS_Data%Sfcprop%stmassxy, wt_h=wt_h)
         call mover(mi, 'tahxy', mn_phys%tahxy, GFS_Data%Sfcprop%tahxy, wt_h=wt_h)
         call mover(mi, 'taussxy', mn_phys%taussxy, GFS_Data%Sfcprop%taussxy, wt_h=wt_h)
         call mover(mi, 'tgxy', mn_phys%tgxy, GFS_Data%Sfcprop%tgxy, wt_h=wt_h)
         call mover(mi, 'tsnoxy', mn_phys%tsnoxy, GFS_Data%Sfcprop%tsnoxy, wt_h)
         call mover(mi, 'tvxy', mn_phys%tvxy, GFS_Data%Sfcprop%tvxy, wt_h=wt_h)
         call mover(mi, 'waxy', mn_phys%waxy, GFS_Data%Sfcprop%waxy, wt_h=wt_h)
         call mover(mi, 'woodxy', mn_phys%woodxy, GFS_Data%Sfcprop%woodxy, wt_h=wt_h)
         call mover(mi, 'wslakexy', mn_phys%wslakexy, GFS_Data%Sfcprop%wslakexy, wt_h=wt_h)
         call mover(mi, 'wtxy', mn_phys%wtxy, GFS_Data%Sfcprop%wtxy, wt_h=wt_h)
         call mover(mi, 'xlaixy', mn_phys%xlaixy, GFS_Data%Sfcprop%xlaixy, wt_h=wt_h)
         call mover(mi, 'xsaixy', mn_phys%xsaixy, GFS_Data%Sfcprop%xsaixy, wt_h=wt_h)
         call mover(mi, 'zsnsoxy', mn_phys%zsnsoxy, GFS_Data%Sfcprop%zsnsoxy, wt_h)
         call mover(mi, 'zwtxy', mn_phys%zwtxy, GFS_Data%Sfcprop%zwtxy, wt_h=wt_h)
      endif lsm_choice_2

      if (GFS_Control%nctp > 0 .and. GFS_Control%cscnv) then
         call mover(mi, 'phy_fctd', mn_phys%phy_fctd, GFS_Data%Tbd%phy_fctd, wt_h=wt_h)
      endif


         rrtmg_types: if(mi%action == DO_COPY_FROM_BLOCK_ARRAYS .or. mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
            ! Copy between arrays of derived type (with no halos) to arrays of reals (with halos).
            call block_copy(mi, mn_phys, GFS_Data%Radtend%sfcfsw)
            call block_copy(mi, mn_phys, GFS_Data%Radtend%sfcflw)
         else
            ! While shifting data, everything is already in intrinsic type arrays with halos.
            call mover(mi, 'sfcflw_dnfx0', mn_phys%sfcflw_dnfx0, wt_h=wt_h)
            call mover(mi, 'sfcflw_dnfx0', mn_phys%sfcfsw_dnfx0, wt_h=wt_h)
            call mover(mi, 'sfcflw_dnfxc', mn_phys%sfcflw_dnfxc, wt_h=wt_h)
            call mover(mi, 'sfcflw_dnfxc', mn_phys%sfcfsw_dnfxc, wt_h=wt_h)
            call mover(mi, 'sfcflw_upfx0', mn_phys%sfcflw_upfx0, wt_h=wt_h)
            call mover(mi, 'sfcflw_upfx0', mn_phys%sfcfsw_upfx0, wt_h=wt_h)
            call mover(mi, 'sfcflw_upfxc', mn_phys%sfcflw_upfxc, wt_h=wt_h)
            call mover(mi, 'sfcflw_upfxc', mn_phys%sfcfsw_upfxc, wt_h=wt_h)
         endif rrtmg_types

         call mover(mi, 'pgr', mn_phys%pgr, GFS_Data%Statein%pgr, wt_h)
         call mover(mi, 'tisfc', mn_phys%tisfc, GFS_Data%Sfcprop%tisfc, wt_h)
         call mover(mi, 'snowd', mn_phys%snowd, GFS_Data%Sfcprop%snowd, wt_h)
         call mover(mi, 'fice', mn_phys%fice, GFS_Data%Sfcprop%fice, wt_h)
         call mover(mi, 'hice', mn_phys%hice, GFS_Data%Sfcprop%hice, wt_h)
         call mover(mi, 'f10m', mn_phys%f10m, GFS_Data%Sfcprop%f10m, wt_h)
         call mover(mi, 'ffmm', mn_phys%ffmm, GFS_Data%Sfcprop%ffmm, wt_h)
         call mover(mi, 'ffhh', mn_phys%ffhh, GFS_Data%Sfcprop%ffhh, wt_h)
         call mover(mi, 't2m', mn_phys%t2m, GFS_Data%Sfcprop%t2m, wt_h=wt_h)
         call mover(mi, 'q2m', mn_phys%q2m, GFS_Data%Sfcprop%q2m, wt_h=wt_h)
         call mover(mi, 'semis', mn_phys%semis, GFS_Data%Radtend%semis, wt_h=wt_h)
         call mover(mi, 'coszdg', mn_phys%coszdg, GFS_Data%Radtend%coszdg, wt_h=wt_h)
         call mover(mi, 'sncovr_ice', mn_phys%sncovr_ice, GFS_Data%Sfcprop%sncovr_ice, wt_h=wt_h)

         call mover(mi, 'snodl', mn_phys%snodl, GFS_Data%Sfcprop%snodl, wt_h=wt_h)
         call mover(mi, 'weasdl', mn_phys%weasdl, GFS_Data%Sfcprop%weasdl, wt_h=wt_h)
         call mover(mi, 'snodi', mn_phys%snodi, GFS_Data%Sfcprop%snodi, wt_h=wt_h)
         call mover(mi, 'weasdi', mn_phys%weasdi, GFS_Data%Sfcprop%weasdi, wt_h=wt_h)
         call mover(mi, 'acsnow_land', mn_phys%acsnow_land, GFS_Data%Sfcprop%acsnow_land, wt_h=wt_h)
         call mover(mi, 'acsnow_ice', mn_phys%acsnow_ice, GFS_Data%Sfcprop%acsnow_ice, wt_h=wt_h)
         call mover(mi, 'th2m', mn_phys%th2m, GFS_Data%Sfcprop%th2m, wt_h=wt_h)


         call mover(mi, 'nirbmui', mn_phys%nirbmui, GFS_Data%Coupling%nirbmui, wt_h=wt_h)
         call mover(mi, 'nirdfui', mn_phys%nirdfui, GFS_Data%Coupling%nirdfui, wt_h=wt_h)
         call mover(mi, 'visbmui', mn_phys%visbmui, GFS_Data%Coupling%visbmui, wt_h=wt_h)
         call mover(mi, 'visdfui', mn_phys%visdfui, GFS_Data%Coupling%visdfui, wt_h=wt_h)

    if(GFS_Control%print_diff_pgr) then
       call mover(mi, 'old_pgr', mn_phys%old_pgr, GFS_Data%Intdiag%old_pgr, wt_h=wt_h)
    endif

    if(GFS_Control%lightning_threat) then
       call mover(mi, 'ltg1_max', mn_phys%ltg1_max, GFS_Data%Intdiag%ltg1_max, wt_h=wt_h)
       call mover(mi, 'ltg2_max', mn_phys%ltg2_max, GFS_Data%Intdiag%ltg2_max, wt_h=wt_h)
       call mover(mi, 'ltg3_max', mn_phys%ltg3_max, GFS_Data%Intdiag%ltg3_max, wt_h=wt_h)
    endif

        !--- Radiation
    call mover(mi, 'fluxr   ', mn_phys%fluxr   , GFS_Data%Intdiag%fluxr   , wt_h=wt_h)
!--- Physics
!--- In/Out
    call mover(mi, 'srunoff ', mn_phys%srunoff , GFS_Data%Intdiag%srunoff , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'evbsa   ', mn_phys%evbsa   , GFS_Data%Intdiag%evbsa   , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'evcwa   ', mn_phys%evcwa   , GFS_Data%Intdiag%evcwa   , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'snohfa  ', mn_phys%snohfa  , GFS_Data%Intdiag%snohfa  , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'transa  ', mn_phys%transa  , GFS_Data%Intdiag%transa  , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'sbsnoa  ', mn_phys%sbsnoa  , GFS_Data%Intdiag%sbsnoa  , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'snowca  ', mn_phys%snowca  , GFS_Data%Intdiag%snowca  , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'evbs    ', mn_phys%evbs    , GFS_Data%Intdiag%evbs    , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'evcw    ', mn_phys%evcw    , GFS_Data%Intdiag%evcw    , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'sbsno   ', mn_phys%sbsno   , GFS_Data%Intdiag%sbsno   , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'trans   ', mn_phys%trans   , GFS_Data%Intdiag%trans   , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'snowmt_land ', mn_phys%snowmt_land , GFS_Data%Intdiag%snowmt_land , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'snowmt_ice  ', mn_phys%snowmt_ice  , GFS_Data%Intdiag%snowmt_ice  , wt_h=wt_h)
    call mover(mi, 'soilm   ', mn_phys%soilm   , GFS_Data%Intdiag%soilm   , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'tmpmin  ', mn_phys%tmpmin  , GFS_Data%Intdiag%tmpmin  , wt_h=wt_h)
    call mover(mi, 'tmpmax  ', mn_phys%tmpmax  , GFS_Data%Intdiag%tmpmax  , wt_h=wt_h)
    call mover(mi, 'dusfc   ', mn_phys%dusfc   , GFS_Data%Intdiag%dusfc   , wt_h=wt_h)
    call mover(mi, 'dvsfc   ', mn_phys%dvsfc   , GFS_Data%Intdiag%dvsfc   , wt_h=wt_h)
    call mover(mi, 'dtsfc   ', mn_phys%dtsfc   , GFS_Data%Intdiag%dtsfc   , wt_h=wt_h)
    call mover(mi, 'dqsfc   ', mn_phys%dqsfc   , GFS_Data%Intdiag%dqsfc   , wt_h=wt_h)
    call mover(mi, 'totprcp ', mn_phys%totprcp , GFS_Data%Intdiag%totprcp , wt_h=wt_h)
    call mover(mi, 'totprcpb', mn_phys%totprcpb, GFS_Data%Intdiag%totprcpb, wt_h=wt_h)
    call mover(mi, 'gflux   ', mn_phys%gflux   , GFS_Data%Intdiag%gflux   , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'dlwsfc  ', mn_phys%dlwsfc  , GFS_Data%Intdiag%dlwsfc  , wt_h=wt_h)
    call mover(mi, 'ulwsfc  ', mn_phys%ulwsfc  , GFS_Data%Intdiag%ulwsfc  , wt_h=wt_h)
    call mover(mi, 'suntim  ', mn_phys%suntim  , GFS_Data%Intdiag%suntim  , wt_h=wt_h)
    call mover(mi, 'runoff  ', mn_phys%runoff  , GFS_Data%Intdiag%runoff  , wt_h=wt_h)
    call mover(mi, 'tecan   ', mn_phys%tecan   , GFS_Data%Intdiag%tecan   , wt_h=wt_h)
    call mover(mi, 'tetran  ', mn_phys%tetran  , GFS_Data%Intdiag%tetran  , wt_h=wt_h)
    call mover(mi, 'tedir   ', mn_phys%tedir   , GFS_Data%Intdiag%tedir   , wt_h=wt_h)
    call mover(mi, 'ep      ', mn_phys%ep      , GFS_Data%Intdiag%ep      , wt_h=wt_h)
    call mover(mi, 'cldwrk  ', mn_phys%cldwrk  , GFS_Data%Intdiag%cldwrk  , wt_h=wt_h)
    call mover(mi, 'dugwd   ', mn_phys%dugwd   , GFS_Data%Intdiag%dugwd   , wt_h=wt_h)
    call mover(mi, 'dvgwd   ', mn_phys%dvgwd   , GFS_Data%Intdiag%dvgwd   , wt_h=wt_h)
    call mover(mi, 'psmean  ', mn_phys%psmean  , GFS_Data%Intdiag%psmean  , wt_h=wt_h)
    call mover(mi, 'cnvprcp ', mn_phys%cnvprcp , GFS_Data%Intdiag%cnvprcp , wt_h=wt_h)
    call mover(mi, 'cnvprcpb', mn_phys%cnvprcpb, GFS_Data%Intdiag%cnvprcpb, wt_h=wt_h)
    call mover(mi, 'spfhmin ', mn_phys%spfhmin , GFS_Data%Intdiag%spfhmin , wt_h=wt_h)
    call mover(mi, 'spfhmax ', mn_phys%spfhmax , GFS_Data%Intdiag%spfhmax , wt_h=wt_h)
    call mover(mi, 'u10mmax ', mn_phys%u10mmax , GFS_Data%Intdiag%u10mmax , wt_h=wt_h)
    call mover(mi, 'v10mmax ', mn_phys%v10mmax , GFS_Data%Intdiag%v10mmax , wt_h=wt_h)
    call mover(mi, 'wind10mmax ', mn_phys%wind10mmax , GFS_Data%Intdiag%wind10mmax , wt_h=wt_h)
    call mover(mi, 'u10max ', mn_phys%u10max , GFS_Data%Intdiag%u10max , wt_h=wt_h)
    call mover(mi, 'v10max ', mn_phys%v10max , GFS_Data%Intdiag%v10max , wt_h=wt_h)
    call mover(mi, 'spd10max ', mn_phys%spd10max , GFS_Data%Intdiag%spd10max , wt_h=wt_h)
    call mover(mi, 'rain    ', mn_phys%rain    , GFS_Data%Intdiag%rain    , wt_h=wt_h)
    call mover(mi, 'rainc   ', mn_phys%rainc   , GFS_Data%Intdiag%rainc   , wt_h=wt_h)
    call mover(mi, 'ice     ', mn_phys%ice     , GFS_Data%Intdiag%ice     , wt_h=wt_h)
    call mover(mi, 'snow    ', mn_phys%snow    , GFS_Data%Intdiag%snow    , wt_h=wt_h)
    call mover(mi, 'graupel ', mn_phys%graupel , GFS_Data%Intdiag%graupel , wt_h=wt_h)
    call mover(mi, 'totice  ', mn_phys%totice  , GFS_Data%Intdiag%totice  , wt_h=wt_h)
    call mover(mi, 'totsnw  ', mn_phys%totsnw  , GFS_Data%Intdiag%totsnw  , wt_h=wt_h)
    call mover(mi, 'totgrp  ', mn_phys%totgrp  , GFS_Data%Intdiag%totgrp  , wt_h=wt_h)
    call mover(mi, 'toticeb ', mn_phys%toticeb , GFS_Data%Intdiag%toticeb , wt_h=wt_h)
    call mover(mi, 'totsnwb ', mn_phys%totsnwb , GFS_Data%Intdiag%totsnwb , wt_h=wt_h)
    call mover(mi, 'totgrpb ', mn_phys%totgrpb , GFS_Data%Intdiag%totgrpb , wt_h=wt_h)
    call mover(mi, 'u10m    ', mn_phys%u10m    , GFS_Data%Intdiag%u10m    , wt_h=wt_h)
    call mover(mi, 'v10m    ', mn_phys%v10m    , GFS_Data%Intdiag%v10m    , wt_h=wt_h)
    call mover(mi, 'dpt2m   ', mn_phys%dpt2m   , GFS_Data%Intdiag%dpt2m   , wt_h=wt_h)
    call mover(mi, 'zlvl    ', mn_phys%zlvl    , GFS_Data%Intdiag%zlvl    , wt_h=wt_h)
    call mover(mi, 'psurf   ', mn_phys%psurf   , GFS_Data%Intdiag%psurf   , wt_h=wt_h)
    call mover(mi, 'pwat    ', mn_phys%pwat    , GFS_Data%Intdiag%pwat    , wt_h=wt_h)
    call mover(mi, 't1      ', mn_phys%t1      , GFS_Data%Intdiag%t1      , wt_h=wt_h)
    call mover(mi, 'q1      ', mn_phys%q1      , GFS_Data%Intdiag%q1      , wt_h=wt_h)
    call mover(mi, 'u1      ', mn_phys%u1      , GFS_Data%Intdiag%u1      , wt_h=wt_h)
    call mover(mi, 'v1      ', mn_phys%v1      , GFS_Data%Intdiag%v1      , wt_h=wt_h)
    call mover(mi, 'chh     ', mn_phys%chh     , GFS_Data%Intdiag%chh     , wt_h=wt_h)
    call mover(mi, 'cmm     ', mn_phys%cmm     , GFS_Data%Intdiag%cmm     , wt_h=wt_h)
    call mover(mi, 'dlwsfci ', mn_phys%dlwsfci , GFS_Data%Intdiag%dlwsfci , wt_h=wt_h)
    call mover(mi, 'ulwsfci ', mn_phys%ulwsfci , GFS_Data%Intdiag%ulwsfci , wt_h=wt_h)
    call mover(mi, 'dswsfci ', mn_phys%dswsfci , GFS_Data%Intdiag%dswsfci , wt_h=wt_h)
    call mover(mi, 'nswsfci ', mn_phys%nswsfci , GFS_Data%Intdiag%nswsfci , wt_h=wt_h)
    call mover(mi, 'uswsfci ', mn_phys%uswsfci , GFS_Data%Intdiag%uswsfci , wt_h=wt_h)
    call mover(mi, 'dusfci  ', mn_phys%dusfci  , GFS_Data%Intdiag%dusfci  , wt_h=wt_h)
    call mover(mi, 'dvsfci  ', mn_phys%dvsfci  , GFS_Data%Intdiag%dvsfci  , wt_h=wt_h)
    call mover(mi, 'dtsfci  ', mn_phys%dtsfci  , GFS_Data%Intdiag%dtsfci  , wt_h=wt_h)
    call mover(mi, 'dqsfci  ', mn_phys%dqsfci  , GFS_Data%Intdiag%dqsfci  , wt_h=wt_h)
    call mover(mi, 'gfluxi  ', mn_phys%gfluxi  , GFS_Data%Intdiag%gfluxi  , wt_h=wt_h)
    call mover(mi, 'epi     ', mn_phys%epi     , GFS_Data%Intdiag%epi     , wt_h=wt_h)
    call mover(mi, 'smcwlt2 ', mn_phys%smcwlt2 , GFS_Data%Intdiag%smcwlt2 , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'smcref2 ', mn_phys%smcref2 , GFS_Data%Intdiag%smcref2 , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    call mover(mi, 'rhonewsn1 ', mn_phys%rhonewsn1 , GFS_Data%Intdiag%rhonewsn1 , wt_h=wt_h)
    call mover(mi, 'frzr    ', mn_phys%frzr    , GFS_Data%Intdiag%frzr    , wt_h=wt_h)
    call mover(mi, 'frzrb   ', mn_phys%frzrb   , GFS_Data%Intdiag%frzrb   , wt_h=wt_h)
    call mover(mi, 'frozr   ', mn_phys%frozr   , GFS_Data%Intdiag%frozr   , wt_h=wt_h)
    call mover(mi, 'frozrb  ', mn_phys%frozrb  , GFS_Data%Intdiag%frozrb  , wt_h=wt_h)
    call mover(mi, 'tsnowp  ', mn_phys%tsnowp  , GFS_Data%Intdiag%tsnowp  , wt_h=wt_h)
    call mover(mi, 'tsnowpb ', mn_phys%tsnowpb , GFS_Data%Intdiag%tsnowpb , wt_h=wt_h)
    if (.not. GFS_Control%lsm == GFS_Control%lsm_ruc) then
       call mover(mi, 'wet1    ', mn_phys%wet1    , GFS_Data%Intdiag%wet1    , wt_h=wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
    end if
    call mover(mi, 'sr       ', mn_phys%sr       , GFS_Data%Intdiag%sr       , wt_h=wt_h)
    call mover(mi, 'tdomr    ', mn_phys%tdomr    , GFS_Data%Intdiag%tdomr    , wt_h=wt_h)
    call mover(mi, 'tdomzr   ', mn_phys%tdomzr   , GFS_Data%Intdiag%tdomzr   , wt_h=wt_h)
    call mover(mi, 'tdomip   ', mn_phys%tdomip   , GFS_Data%Intdiag%tdomip   , wt_h=wt_h)
    call mover(mi, 'tdoms    ', mn_phys%tdoms    , GFS_Data%Intdiag%tdoms    , wt_h=wt_h)
    call mover(mi, 'zmtnblck ', mn_phys%zmtnblck , GFS_Data%Intdiag%zmtnblck , wt_h=wt_h)

    if(GFS_Control%lsm == GFS_Control%lsm_noahmp) then
       call mover(mi, 'paha    ', mn_phys%paha    , GFS_Data%Intdiag%paha    , wt_h=wt_h)
       call mover(mi, 'twa     ', mn_phys%twa     , GFS_Data%Intdiag%twa     , wt_h=wt_h)
       call mover(mi, 'pahi    ', mn_phys%pahi    , GFS_Data%Intdiag%pahi    , wt_h=wt_h)
    endif

    
    ! F-A MP scheme
    if (GFS_Control%imp_physics == GFS_Control%imp_physics_fer_hires) then
       call mover(mi, 'train     ', mn_phys%train     , GFS_Data%Intdiag%train     , wt_h=wt_h)
    end if
    call mover(mi, 'cldfra     ', mn_phys%cldfra     , GFS_Data%Intdiag%cldfra     , wt_h=wt_h)
    call mover(mi, 'cldfra2d   ', mn_phys%cldfra2d   , GFS_Data%Intdiag%cldfra2d   , wt_h=wt_h)
    call mover(mi, 'total_albedo ', mn_phys%total_albedo , GFS_Data%Intdiag%total_albedo , wt_h=wt_h)
    call mover(mi, 'lwp_ex ', mn_phys%lwp_ex , GFS_Data%Intdiag%lwp_ex , wt_h=wt_h)
    call mover(mi, 'iwp_ex ', mn_phys%iwp_ex , GFS_Data%Intdiag%iwp_ex , wt_h=wt_h)
    call mover(mi, 'lwp_fc ', mn_phys%lwp_fc , GFS_Data%Intdiag%lwp_fc , wt_h=wt_h)
    call mover(mi, 'iwp_fc ', mn_phys%iwp_fc , GFS_Data%Intdiag%iwp_fc , wt_h=wt_h)

    !--- 3D diagnostics
    if (GFS_Control%ldiag3d) then
       call mover(mi, 'dtend', mn_phys%dtend, GFS_Data%Intdiag%dtend, wt_h=wt_h)
      if (GFS_Control%qdiag3d) then
         call mover(mi, 'upd_mf ', mn_phys%upd_mf , GFS_Data%Intdiag%upd_mf , wt_h=wt_h)
         call mover(mi, 'dwn_mf ', mn_phys%dwn_mf , GFS_Data%Intdiag%dwn_mf , wt_h=wt_h)
         call mover(mi, 'det_mf ', mn_phys%det_mf , GFS_Data%Intdiag%det_mf , wt_h=wt_h)
      endif
      if (GFS_Control%oz_phys_2015) then
         call mover(mi, 'do3_dt_prd', mn_phys%do3_dt_prd, GFS_Data%Intdiag%do3_dt_prd, wt_h=wt_h)
         call mover(mi, 'do3_dt_ozmx', mn_phys%do3_dt_ozmx, GFS_Data%Intdiag%do3_dt_ozmx, wt_h=wt_h)
         call mover(mi, 'do3_dt_temp', mn_phys%do3_dt_temp, GFS_Data%Intdiag%do3_dt_temp, wt_h=wt_h)
         call mover(mi, 'do3_dt_ohoz', mn_phys%do3_dt_ohoz, GFS_Data%Intdiag%do3_dt_ohoz, wt_h=wt_h)
      endif
    endif

! UGWP
    call mover(mi, 'zmtb      ', mn_phys%zmtb      , GFS_Data%Intdiag%zmtb      , wt_h=wt_h)
    call mover(mi, 'zogw      ', mn_phys%zogw      , GFS_Data%Intdiag%zogw      , wt_h=wt_h)
    call mover(mi, 'zlwb      ', mn_phys%zlwb      , GFS_Data%Intdiag%zlwb      , wt_h=wt_h)
    call mover(mi, 'tau_ogw   ', mn_phys%tau_ogw   , GFS_Data%Intdiag%tau_ogw   , wt_h=wt_h)
    call mover(mi, 'tau_ngw   ', mn_phys%tau_ngw   , GFS_Data%Intdiag%tau_ngw   , wt_h=wt_h)
    call mover(mi, 'tau_mtb   ', mn_phys%tau_mtb   , GFS_Data%Intdiag%tau_mtb   , wt_h=wt_h)
    call mover(mi, 'tau_tofd  ', mn_phys%tau_tofd  , GFS_Data%Intdiag%tau_tofd  , wt_h=wt_h)
    call mover(mi, 'dudt_gw   ', mn_phys%dudt_gw   , GFS_Data%Intdiag%dudt_gw   , wt_h=wt_h)
    call mover(mi, 'dvdt_gw   ', mn_phys%dvdt_gw   , GFS_Data%Intdiag%dvdt_gw   , wt_h=wt_h)
    call mover(mi, 'dtdt_gw   ', mn_phys%dtdt_gw   , GFS_Data%Intdiag%dtdt_gw   , wt_h=wt_h)
    call mover(mi, 'kdis_gw   ', mn_phys%kdis_gw   , GFS_Data%Intdiag%kdis_gw   , wt_h=wt_h)

    if (GFS_Control%ldiag_ugwp) then
       call mover(mi, 'du3dt_dyn  ', mn_phys%du3dt_dyn  , GFS_Data%Intdiag%du3dt_dyn  , wt_h=wt_h)
       call mover(mi, 'du3dt_pbl  ', mn_phys%du3dt_pbl  , GFS_Data%Intdiag%du3dt_pbl  , wt_h=wt_h)
       call mover(mi, 'dv3dt_pbl  ', mn_phys%dv3dt_pbl  , GFS_Data%Intdiag%dv3dt_pbl  , wt_h=wt_h)
       call mover(mi, 'dt3dt_pbl  ', mn_phys%dt3dt_pbl  , GFS_Data%Intdiag%dt3dt_pbl  , wt_h=wt_h)
       call mover(mi, 'du3dt_ogw  ', mn_phys%du3dt_ogw  , GFS_Data%Intdiag%du3dt_ogw  , wt_h=wt_h)
       call mover(mi, 'du3dt_mtb  ', mn_phys%du3dt_mtb  , GFS_Data%Intdiag%du3dt_mtb  , wt_h=wt_h)
       call mover(mi, 'du3dt_tms  ', mn_phys%du3dt_tms  , GFS_Data%Intdiag%du3dt_tms  , wt_h=wt_h)
       call mover(mi, 'du3dt_ngw  ', mn_phys%du3dt_ngw  , GFS_Data%Intdiag%du3dt_ngw  , wt_h=wt_h)
       call mover(mi, 'dv3dt_ngw  ', mn_phys%dv3dt_ngw  , GFS_Data%Intdiag%dv3dt_ngw  , wt_h=wt_h)
       call mover(mi, 'dudt_tot  ', mn_phys%dudt_tot  , GFS_Data%Intdiag%dudt_tot  , wt_h=wt_h)
       call mover(mi, 'dvdt_tot  ', mn_phys%dvdt_tot  , GFS_Data%Intdiag%dvdt_tot  , wt_h=wt_h)
       call mover(mi, 'dtdt_tot  ', mn_phys%dtdt_tot  , GFS_Data%Intdiag%dtdt_tot  , wt_h=wt_h)
       call mover(mi, 'uav_ugwp  ', mn_phys%uav_ugwp  , GFS_Data%Intdiag%uav_ugwp  , wt_h=wt_h)
       call mover(mi, 'tav_ugwp  ', mn_phys%tav_ugwp  , GFS_Data%Intdiag%tav_ugwp  , wt_h=wt_h)
       call mover(mi, 'dws3dt_ogw ', mn_phys%dws3dt_ogw , GFS_Data%Intdiag%dws3dt_ogw , wt_h=wt_h)
       call mover(mi, 'dws3dt_obl ', mn_phys%dws3dt_obl , GFS_Data%Intdiag%dws3dt_obl , wt_h=wt_h)
       call mover(mi, 'dws3dt_oss ', mn_phys%dws3dt_oss , GFS_Data%Intdiag%dws3dt_oss , wt_h=wt_h)
       call mover(mi, 'dws3dt_ofd ', mn_phys%dws3dt_ofd , GFS_Data%Intdiag%dws3dt_ofd , wt_h=wt_h)
       call mover(mi, 'ldu3dt_ogw  ', mn_phys%ldu3dt_ogw  , GFS_Data%Intdiag%ldu3dt_ogw  , wt_h=wt_h)
       call mover(mi, 'ldu3dt_obl  ', mn_phys%ldu3dt_obl  , GFS_Data%Intdiag%ldu3dt_obl  , wt_h=wt_h)
       call mover(mi, 'ldu3dt_oss  ', mn_phys%ldu3dt_oss  , GFS_Data%Intdiag%ldu3dt_oss  , wt_h=wt_h)
       call mover(mi, 'ldu3dt_ofd  ', mn_phys%ldu3dt_ofd  , GFS_Data%Intdiag%ldu3dt_ofd  , wt_h=wt_h)
       call mover(mi, 'ldu3dt_ngw ', mn_phys%ldu3dt_ngw , GFS_Data%Intdiag%ldu3dt_ngw , wt_h=wt_h)
       call mover(mi, 'ldv3dt_ngw ', mn_phys%ldv3dt_ngw , GFS_Data%Intdiag%ldv3dt_ngw , wt_h=wt_h)
       call mover(mi, 'ldt3dt_ngw ', mn_phys%ldt3dt_ngw , GFS_Data%Intdiag%ldt3dt_ngw , wt_h=wt_h)
    endif

    if (GFS_Control%do_ugwp_v1 .or. GFS_Control%ldiag_ugwp) then
       call mover(mi, 'dudt_ogw  ', mn_phys%dudt_ogw  , GFS_Data%Intdiag%dudt_ogw  , wt_h=wt_h)
       call mover(mi, 'dvdt_ogw  ', mn_phys%dvdt_ogw  , GFS_Data%Intdiag%dvdt_ogw  , wt_h=wt_h)
       call mover(mi, 'dudt_obl  ', mn_phys%dudt_obl  , GFS_Data%Intdiag%dudt_obl  , wt_h=wt_h)
       call mover(mi, 'dvdt_obl  ', mn_phys%dvdt_obl  , GFS_Data%Intdiag%dvdt_obl  , wt_h=wt_h)
       call mover(mi, 'dudt_oss  ', mn_phys%dudt_oss  , GFS_Data%Intdiag%dudt_oss  , wt_h=wt_h)
       call mover(mi, 'dvdt_oss  ', mn_phys%dvdt_oss  , GFS_Data%Intdiag%dvdt_oss  , wt_h=wt_h)
       call mover(mi, 'dudt_ofd  ', mn_phys%dudt_ofd  , GFS_Data%Intdiag%dudt_ofd  , wt_h=wt_h)
       call mover(mi, 'dvdt_ofd  ', mn_phys%dvdt_ofd  , GFS_Data%Intdiag%dvdt_ofd  , wt_h=wt_h)
       call mover(mi, 'du_ogwcol ', mn_phys%du_ogwcol , GFS_Data%Intdiag%du_ogwcol , wt_h=wt_h)
       call mover(mi, 'dv_ogwcol ', mn_phys%dv_ogwcol , GFS_Data%Intdiag%dv_ogwcol , wt_h=wt_h)
       call mover(mi, 'du_oblcol ', mn_phys%du_oblcol , GFS_Data%Intdiag%du_oblcol , wt_h=wt_h)
       call mover(mi, 'dv_oblcol ', mn_phys%dv_oblcol , GFS_Data%Intdiag%dv_oblcol , wt_h=wt_h)
       call mover(mi, 'du_osscol ', mn_phys%du_osscol , GFS_Data%Intdiag%du_osscol , wt_h=wt_h)
       call mover(mi, 'dv_osscol ', mn_phys%dv_osscol , GFS_Data%Intdiag%dv_osscol , wt_h=wt_h)
       call mover(mi, 'du_ofdcol ', mn_phys%du_ofdcol , GFS_Data%Intdiag%du_ofdcol , wt_h=wt_h)
       call mover(mi, 'dv_ofdcol ', mn_phys%dv_ofdcol , GFS_Data%Intdiag%dv_ofdcol , wt_h=wt_h)
       call mover(mi, 'du3_ogwcol ', mn_phys%du3_ogwcol , GFS_Data%Intdiag%du3_ogwcol , wt_h=wt_h)
       call mover(mi, 'dv3_ogwcol ', mn_phys%dv3_ogwcol , GFS_Data%Intdiag%dv3_ogwcol , wt_h=wt_h)
       call mover(mi, 'du3_oblcol ', mn_phys%du3_oblcol , GFS_Data%Intdiag%du3_oblcol , wt_h=wt_h)
       call mover(mi, 'dv3_oblcol ', mn_phys%dv3_oblcol , GFS_Data%Intdiag%dv3_oblcol , wt_h=wt_h)
       call mover(mi, 'du3_osscol ', mn_phys%du3_osscol , GFS_Data%Intdiag%du3_osscol , wt_h=wt_h)
       call mover(mi, 'dv3_osscol ', mn_phys%dv3_osscol , GFS_Data%Intdiag%dv3_osscol , wt_h=wt_h)
       call mover(mi, 'du3_ofdcol ', mn_phys%du3_ofdcol , GFS_Data%Intdiag%du3_ofdcol , wt_h=wt_h)
       call mover(mi, 'dv3_ofdcol ', mn_phys%dv3_ofdcol , GFS_Data%Intdiag%dv3_ofdcol , wt_h=wt_h)
    else
       call mover(mi, 'dudt_ogw  ', mn_phys%dudt_ogw  , GFS_Data%Intdiag%dudt_ogw  , wt_h=wt_h)
    endif

    !--- 3D diagnostics for Thompson MP / GFDL MP
    call mover(mi, 'refl_10cm', mn_phys%refl_10cm, GFS_Data%Intdiag%refl_10cm, wt_h=wt_h)
    call mover(mi, 'max_hail_diam_sfc', mn_phys%max_hail_diam_sfc, GFS_Data%Intdiag%max_hail_diam_sfc, wt_h=wt_h)

    !--- New PBL Diagnostics
    call mover(mi, 'dkt', mn_phys%dkt, GFS_Data%Intdiag%dkt, wt_h=wt_h)
    call mover(mi, 'dku', mn_phys%dku, GFS_Data%Intdiag%dku, wt_h=wt_h)

    !--  New max hourly diag.
    call mover(mi, 'refdmax', mn_phys%refdmax, GFS_Data%Intdiag%refdmax, wt_h=wt_h)
    call mover(mi, 'refdmax263k', mn_phys%refdmax263k, GFS_Data%Intdiag%refdmax263k, wt_h=wt_h)
    call mover(mi, 't02max', mn_phys%t02max, GFS_Data%Intdiag%t02max, wt_h=wt_h)
    call mover(mi, 't02min', mn_phys%t02min, GFS_Data%Intdiag%t02min, wt_h=wt_h)
    call mover(mi, 'rh02max', mn_phys%rh02max, GFS_Data%Intdiag%rh02max, wt_h=wt_h)
    call mover(mi, 'rh02min', mn_phys%rh02min, GFS_Data%Intdiag%rh02min, wt_h=wt_h)
    call mover(mi, 'pratemax', mn_phys%pratemax, GFS_Data%Intdiag%pratemax, wt_h=wt_h)

    !--- MYNN variables:
    if (GFS_Control%do_mynnedmf) then
      if (GFS_Control%bl_mynn_output .ne. 0) then
         call mover(mi, 'edmf_a    ', mn_phys%edmf_a    , GFS_Data%Intdiag%edmf_a    , wt_h=wt_h)
         call mover(mi, 'edmf_w    ', mn_phys%edmf_w    , GFS_Data%Intdiag%edmf_w    , wt_h=wt_h)
         call mover(mi, 'edmf_qt   ', mn_phys%edmf_qt   , GFS_Data%Intdiag%edmf_qt   , wt_h=wt_h)
         call mover(mi, 'edmf_thl  ', mn_phys%edmf_thl  , GFS_Data%Intdiag%edmf_thl  , wt_h=wt_h)
         call mover(mi, 'edmf_ent  ', mn_phys%edmf_ent  , GFS_Data%Intdiag%edmf_ent  , wt_h=wt_h)
         call mover(mi, 'edmf_qc   ', mn_phys%edmf_qc   , GFS_Data%Intdiag%edmf_qc   , wt_h=wt_h)
         call mover(mi, 'sub_thl   ', mn_phys%sub_thl   , GFS_Data%Intdiag%sub_thl   , wt_h=wt_h)
         call mover(mi, 'sub_sqv   ', mn_phys%sub_sqv   , GFS_Data%Intdiag%sub_sqv   , wt_h=wt_h)
         call mover(mi, 'det_thl   ', mn_phys%det_thl   , GFS_Data%Intdiag%det_thl   , wt_h=wt_h)
         call mover(mi, 'det_sqv   ', mn_phys%det_sqv   , GFS_Data%Intdiag%det_sqv   , wt_h=wt_h)
      endif
      if (GFS_Control%tke_budget .gt. 0) then
         call mover(mi, 'dqke      ', mn_phys%dqke      , GFS_Data%Intdiag%dqke      , wt_h=wt_h)
         call mover(mi, 'qwt       ', mn_phys%qwt       , GFS_Data%Intdiag%qwt       , wt_h=wt_h)
         call mover(mi, 'qshear    ', mn_phys%qshear    , GFS_Data%Intdiag%qshear    , wt_h=wt_h)
         call mover(mi, 'qbuoy     ', mn_phys%qbuoy     , GFS_Data%Intdiag%qbuoy     , wt_h=wt_h)
         call mover(mi, 'qdiss     ', mn_phys%qdiss     , GFS_Data%Intdiag%qdiss     , wt_h=wt_h)
      endif
      call mover(mi, 'maxwidth  ', mn_phys%maxwidth  , GFS_Data%Intdiag%maxwidth  , wt_h=wt_h)
      call mover(mi, 'maxmf     ', mn_phys%maxmf     , GFS_Data%Intdiag%maxmf     , wt_h=wt_h)
      call mover(mi, 'ztop_plume', mn_phys%ztop_plume, GFS_Data%Intdiag%ztop_plume, wt_h=wt_h)
      call mover(mi, 'ktop_plume', mn_phys%ktop_plume, GFS_Data%Intdiag%ktop_plume, wt_h=wt_h)
      call mover(mi, 'exch_h    ', mn_phys%exch_h    , GFS_Data%Intdiag%exch_h    , wt_h=wt_h)
      call mover(mi, 'exch_m    ', mn_phys%exch_m    , GFS_Data%Intdiag%exch_m    , wt_h=wt_h)
    endif

    ! Extended diagnostics for Thompson MP
    if (GFS_Control%ext_diag_thompson) then
       call mover(mi, 'thompson_ext_diag3d', mn_phys%thompson_ext_diag3d, GFS_Data%Intdiag%thompson_ext_diag3d, wt_h=wt_h)
    endif

    ! Air quality diagnostics
    ! -- initialize diagnostic variables
    if (GFS_Control%cplaqm) then
       call mover(mi, 'aod', mn_phys%aod, GFS_Data%Intdiag%aod, wt_h=wt_h)
    end if

    ! Auxiliary arrays in output for debugging
    if (GFS_Control%naux2d>0) then
       call mover(mi, 'aux2d', mn_phys%aux2d, GFS_Data%Intdiag%aux2d, wt_h=wt_h)
    endif
    if (GFS_Control%naux3d>0) then
       call mover(mi, 'aux3d', mn_phys%aux3d, GFS_Data%Intdiag%aux3d, wt_h=wt_h)
    endif


         
      if(omit_hafs_bugs) then
         !--------------------------------------------------------------------------------
         ! Everything in this section DOES change the HAFS results.
         ! They should be enabled anyway, since they're needed by the physics.
         !--------------------------------------------------------------------------------
         call mover(mi, 'sncovr', mn_phys%sncovr, GFS_Data%Sfcprop%sncovr, wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
         call mover(mi, 'weasd', mn_phys%weasd, GFS_Data%Sfcprop%weasd, wt_h, &
               halo_land_mask_fill=0.0_kind_phys, if_negative=0.0_kind_phys)
         call mover(mi, 'srflag', mn_phys%srflag, GFS_Data%Sfcprop%srflag, wt_h)

         call mover(mi, 'nirbmdi', mn_phys%nirbmdi, GFS_Data%Coupling%nirbmdi, wt_h=wt_h)
         call mover(mi, 'nirdfdi', mn_phys%nirdfdi, GFS_Data%Coupling%nirdfdi, wt_h=wt_h)
         call mover(mi, 'visbmdi', mn_phys%visbmdi, GFS_Data%Coupling%visbmdi, wt_h=wt_h)
         call mover(mi, 'visdfdi', mn_phys%visdfdi, GFS_Data%Coupling%visdfdi, wt_h=wt_h)

         call mover(mi, 'sfcdsw', mn_phys%sfcdsw, GFS_Data%Coupling%sfcdsw, wt_h=wt_h)
         call mover(mi, 'sfcdlw', mn_phys%sfcdlw, GFS_Data%Coupling%sfcdlw, wt_h=wt_h)
         call mover(mi, 'sfcnsw', mn_phys%sfcnsw, GFS_Data%Coupling%sfcnsw, wt_h=wt_h)

         call mover(mi, 'sfalb', mn_phys%sfalb, GFS_Data%Radtend%sfalb, wt_h=wt_h)
         call mover(mi, 'coszen', mn_phys%coszen, GFS_Data%Radtend%coszen, wt_h=wt_h)
         call mover(mi, 'tsflw', mn_phys%tsflw, GFS_Data%Radtend%tsflw, wt_h=wt_h)


         call mover(mi, 'swhc', mn_phys%swhc, GFS_Data%Radtend%swhc, wt_h=wt_h)
         call mover(mi, 'lwhc', mn_phys%lwhc, GFS_Data%Radtend%lwhc, wt_h=wt_h)
         call mover(mi, 'htrsw', mn_phys%htrsw, GFS_Data%Radtend%htrsw, wt_h=wt_h)
         call mover(mi, 'htrlw', mn_phys%htrlw, GFS_Data%Radtend%htrlw, wt_h=wt_h)

         call mover(mi, 'tiice', mn_phys%tiice, GFS_Data%Sfcprop%tiice, wt_h=wt_h)
      endif


      ! Coupling
       if (GFS_control%do_RRTMGP) then
          call mover(mi, 'fluxlwUP_radtime', mn_phys%fluxlwUP_radtime, GFS_Data%Coupling%fluxlwUP_radtime, wt_h=wt_h)
          call mover(mi, 'fluxlwDOWN_radtime', mn_phys%fluxlwDOWN_radtime, GFS_Data%Coupling%fluxlwDOWN_radtime, wt_h=wt_h)
          call mover(mi, 'fluxlwUP_jac', mn_phys%fluxlwUP_jac, GFS_Data%Coupling%fluxlwUP_jac, wt_h=wt_h)
          call mover(mi, 'tsfc_radtime', mn_phys%tsfc_radtime, GFS_Data%Coupling%tsfc_radtime, wt_h=wt_h)
       endif

! BEGINNING OF CONFLICTS------------------------------------------------------------------------------------
! This section was from an older attempt to get rrfs physics working. It is untested in the latest version.
       UNTESTED_IN_LATEST_VERSION: IF(.FALSE.) THEN

       if (GFS_control%cplflx .or. GFS_control%cplchm .or. GFS_control%cplwav) then
          !--- instantaneous quantities
          call mover(mi, 'u10mi_cpl', mn_phys%u10mi_cpl, GFS_Data%Coupling%u10mi_cpl, wt_h=wt_h)
          call mover(mi, 'v10mi_cpl', mn_phys%v10mi_cpl, GFS_Data%Coupling%v10mi_cpl, wt_h=wt_h)
       endif

       if (GFS_control%cplflx .or. GFS_control%cplchm .or. GFS_control%cpllnd) then
          !--- instantaneous quantities
          call mover(mi, 'tsfci_cpl', mn_phys%tsfci_cpl, GFS_Data%Coupling%tsfci_cpl, wt_h=wt_h)
       endif

       if (GFS_control%cplflx .or. GFS_control%cpllnd) then
          call mover(mi, 'dlwsfci_cpl', mn_phys%dlwsfci_cpl, GFS_Data%Coupling%dlwsfci_cpl, wt_h=wt_h)
          call mover(mi, 'dswsfci_cpl', mn_phys%dswsfci_cpl, GFS_Data%Coupling%dswsfci_cpl, wt_h=wt_h)
          call mover(mi, 'dlwsfc_cpl', mn_phys%dlwsfc_cpl, GFS_Data%Coupling%dlwsfc_cpl, wt_h=wt_h)
          call mover(mi, 'dswsfc_cpl', mn_phys%dswsfc_cpl, GFS_Data%Coupling%dswsfc_cpl, wt_h=wt_h)
          call mover(mi, 'psurfi_cpl', mn_phys%psurfi_cpl, GFS_Data%Coupling%psurfi_cpl, wt_h=wt_h)
          call mover(mi, 'nswsfc_cpl', mn_phys%nswsfc_cpl, GFS_Data%Coupling%nswsfc_cpl, wt_h=wt_h)
          call mover(mi, 'nswsfci_cpl', mn_phys%nswsfci_cpl, GFS_Data%Coupling%nswsfci_cpl, wt_h=wt_h)
          call mover(mi, 'nnirbmi_cpl', mn_phys%nnirbmi_cpl, GFS_Data%Coupling%nnirbmi_cpl, wt_h=wt_h)
          call mover(mi, 'nnirdfi_cpl', mn_phys%nnirdfi_cpl, GFS_Data%Coupling%nnirdfi_cpl, wt_h=wt_h)
          call mover(mi, 'nvisbmi_cpl', mn_phys%nvisbmi_cpl, GFS_Data%Coupling%nvisbmi_cpl, wt_h=wt_h)
          call mover(mi, 'nvisdfi_cpl', mn_phys%nvisdfi_cpl, GFS_Data%Coupling%nvisdfi_cpl, wt_h=wt_h)
          call mover(mi, 'nnirbm_cpl', mn_phys%nnirbm_cpl, GFS_Data%Coupling%nnirbm_cpl, wt_h=wt_h)
          call mover(mi, 'nnirdf_cpl', mn_phys%nnirdf_cpl, GFS_Data%Coupling%nnirdf_cpl, wt_h=wt_h)
          call mover(mi, 'nvisbm_cpl', mn_phys%nvisbm_cpl, GFS_Data%Coupling%nvisbm_cpl, wt_h=wt_h)
          call mover(mi, 'nvisdf_cpl', mn_phys%nvisdf_cpl, GFS_Data%Coupling%nvisdf_cpl, wt_h=wt_h)
       endif

    if (GFS_control%cplflx) then
      !--- incoming quantities
       call mover(mi, 'slimskin_cpl', mn_phys%slimskin_cpl, GFS_Data%Coupling%slimskin_cpl, wt_h=wt_h)
       call mover(mi, 'dusfcin_cpl', mn_phys%dusfcin_cpl, GFS_Data%Coupling%dusfcin_cpl, wt_h=wt_h)
       call mover(mi, 'dvsfcin_cpl', mn_phys%dvsfcin_cpl, GFS_Data%Coupling%dvsfcin_cpl, wt_h=wt_h)
       call mover(mi, 'dtsfcin_cpl', mn_phys%dtsfcin_cpl, GFS_Data%Coupling%dtsfcin_cpl, wt_h=wt_h)
       call mover(mi, 'dqsfcin_cpl', mn_phys%dqsfcin_cpl, GFS_Data%Coupling%dqsfcin_cpl, wt_h=wt_h)
       call mover(mi, 'ulwsfcin_cpl', mn_phys%ulwsfcin_cpl, GFS_Data%Coupling%ulwsfcin_cpl, wt_h=wt_h)
       call mover(mi, 'hsnoin_cpl', mn_phys%hsnoin_cpl, GFS_Data%Coupling%hsnoin_cpl, wt_h=wt_h)

      ! -- Coupling options to retrive atmosphere-ocean fluxes from mediator
      if (GFS_control%use_med_flux) then
         call mover(mi, 'dusfcin_med', mn_phys%dusfcin_med, GFS_Data%Coupling%dusfcin_med, wt_h=wt_h)
         call mover(mi, 'dvsfcin_med', mn_phys%dvsfcin_med, GFS_Data%Coupling%dvsfcin_med, wt_h=wt_h)
         call mover(mi, 'dtsfcin_med', mn_phys%dtsfcin_med, GFS_Data%Coupling%dtsfcin_med, wt_h=wt_h)
         call mover(mi, 'dqsfcin_med', mn_phys%dqsfcin_med, GFS_Data%Coupling%dqsfcin_med, wt_h=wt_h)
         call mover(mi, 'ulwsfcin_med', mn_phys%ulwsfcin_med, GFS_Data%Coupling%ulwsfcin_med, wt_h=wt_h)
      end if

      !--- accumulated quantities
      call mover(mi, 'dusfc_cpl', mn_phys%dusfc_cpl, GFS_Data%Coupling%dusfc_cpl, wt_h=wt_h)
      call mover(mi, 'dvsfc_cpl', mn_phys%dvsfc_cpl, GFS_Data%Coupling%dvsfc_cpl, wt_h=wt_h)
      call mover(mi, 'dtsfc_cpl', mn_phys%dtsfc_cpl, GFS_Data%Coupling%dtsfc_cpl, wt_h=wt_h)
      call mover(mi, 'dqsfc_cpl', mn_phys%dqsfc_cpl, GFS_Data%Coupling%dqsfc_cpl, wt_h=wt_h)
      call mover(mi, 'dnirbm_cpl', mn_phys%dnirbm_cpl, GFS_Data%Coupling%dnirbm_cpl, wt_h=wt_h)
      call mover(mi, 'dnirdf_cpl', mn_phys%dnirdf_cpl, GFS_Data%Coupling%dnirdf_cpl, wt_h=wt_h)
      call mover(mi, 'dvisbm_cpl', mn_phys%dvisbm_cpl, GFS_Data%Coupling%dvisbm_cpl, wt_h=wt_h)
      call mover(mi, 'dvisdf_cpl', mn_phys%dvisdf_cpl, GFS_Data%Coupling%dvisdf_cpl, wt_h=wt_h)
      call mover(mi, 'nlwsfc_cpl', mn_phys%nlwsfc_cpl, GFS_Data%Coupling%nlwsfc_cpl, wt_h=wt_h)

      !--- instantaneous quantities
      call mover(mi, 'dusfci_cpl', mn_phys%dusfci_cpl, GFS_Data%Coupling%dusfci_cpl, wt_h=wt_h)
      call mover(mi, 'dvsfci_cpl', mn_phys%dvsfci_cpl, GFS_Data%Coupling%dvsfci_cpl, wt_h=wt_h)
      call mover(mi, 'dtsfci_cpl', mn_phys%dtsfci_cpl, GFS_Data%Coupling%dtsfci_cpl, wt_h=wt_h)
      call mover(mi, 'dqsfci_cpl', mn_phys%dqsfci_cpl, GFS_Data%Coupling%dqsfci_cpl, wt_h=wt_h)
      call mover(mi, 'dnirbmi_cpl', mn_phys%dnirbmi_cpl, GFS_Data%Coupling%dnirbmi_cpl, wt_h=wt_h)
      call mover(mi, 'dnirdfi_cpl', mn_phys%dnirdfi_cpl, GFS_Data%Coupling%dnirdfi_cpl, wt_h=wt_h)
      call mover(mi, 'dvisbmi_cpl', mn_phys%dvisbmi_cpl, GFS_Data%Coupling%dvisbmi_cpl, wt_h=wt_h)
      call mover(mi, 'dvisdfi_cpl', mn_phys%dvisdfi_cpl, GFS_Data%Coupling%dvisdfi_cpl, wt_h=wt_h)
      call mover(mi, 'nlwsfci_cpl', mn_phys%nlwsfci_cpl, GFS_Data%Coupling%nlwsfci_cpl, wt_h=wt_h)
      call mover(mi, 't2mi_cpl', mn_phys%t2mi_cpl, GFS_Data%Coupling%t2mi_cpl, wt_h=wt_h)
      call mover(mi, 'q2mi_cpl', mn_phys%q2mi_cpl, GFS_Data%Coupling%q2mi_cpl, wt_h=wt_h)
      call mover(mi, 'oro_cpl', mn_phys%oro_cpl, GFS_Data%Coupling%oro_cpl, wt_h=wt_h)
      call mover(mi, 'slmsk_cpl', mn_phys%slmsk_cpl, GFS_Data%Coupling%slmsk_cpl, wt_h=wt_h)
   endif


      ! ------------------------------------------------------------
      ! Radtend

      call mover(mi, 'lwhd', mn_phys%lwhc, GFS_Data%Radtend%lwhc, wt_h=wt_h)

      ! ------------------------------------------------------------
      ! Tbd

       call mover(mi, 'rann', mn_phys%rann, GFS_data%Tbd%rann, wt_h=wt_h)

       call mover(mi, 'acv', mn_phys%acv, GFS_data%Tbd%acv, wt_h=wt_h)
       call mover(mi, 'acvb', mn_phys%acvb, GFS_data%Tbd%acvb, wt_h=wt_h)
       call mover(mi, 'acvt', mn_phys%acvt, GFS_data%Tbd%acvt, wt_h=wt_h)

       if (GFS_control%cplflx .or. GFS_control%cplchm .or. GFS_control%cpllnd) then
          call mover(mi, 'drain_cpl', mn_phys%drain_cpl, GFS_data%Tbd%drain_cpl, wt_h=wt_h)
          call mover(mi, 'dsnow_cpl', mn_phys%dsnow_cpl, GFS_data%Tbd%dsnow_cpl, wt_h=wt_h)
       endif

       if (GFS_control%imfdeepcnv == GFS_control%imfdeepcnv_gf .or. GFS_control%imfdeepcnv == GFS_control%imfdeepcnv_ntiedtke .or.  GFS_control%imfdeepcnv == GFS_control%imfdeepcnv_c3) then
          call mover(mi, 'forcet', mn_phys%forcet, GFS_Data%Tbd%forcet, wt_h=wt_h)
          call mover(mi, 'forceq', mn_phys%forceq, GFS_Data%Tbd%forceq, wt_h=wt_h)
       end if

       if (GFS_control%imfdeepcnv == GFS_control%imfdeepcnv_gf .or.  GFS_control%imfdeepcnv == GFS_control%imfdeepcnv_c3) then
          call mover(mi, 'cactiv', mn_phys%cactiv, GFS_Data%Tbd%cactiv, wt_h=wt_h)
          call mover(mi, 'cactiv_m', mn_phys%cactiv_m, GFS_Data%Tbd%cactiv_m, wt_h=wt_h)
       endif

       ! MYJ variables
       if (GFS_control%do_myjsfc.or.GFS_control%do_myjpbl) then
          !print*,"Allocating all MYJ surface variables:"
          call mover(mi, 'phy_myj_qsfc', mn_phys%phy_myj_qsfc, GFS_data%Tbd%phy_myj_qsfc, wt_h=wt_h)
          call mover(mi, 'phy_myj_thz0', mn_phys%phy_myj_thz0, GFS_data%Tbd%phy_myj_thz0, wt_h=wt_h)
          call mover(mi, 'phy_myj_qz0', mn_phys%phy_myj_qz0, GFS_data%Tbd%phy_myj_qz0, wt_h=wt_h)
          call mover(mi, 'phy_myj_uz0', mn_phys%phy_myj_uz0, GFS_data%Tbd%phy_myj_uz0, wt_h=wt_h)
          call mover(mi, 'phy_myj_vz0', mn_phys%phy_myj_vz0, GFS_data%Tbd%phy_myj_vz0, wt_h=wt_h)
          call mover(mi, 'phy_myj_akhs', mn_phys%phy_myj_akhs, GFS_data%Tbd%phy_myj_akhs, wt_h=wt_h)
          call mover(mi, 'phy_myj_akms', mn_phys%phy_myj_akms, GFS_data%Tbd%phy_myj_akms, wt_h=wt_h)
          call mover(mi, 'phy_myj_chkqlm', mn_phys%phy_myj_chkqlm, GFS_data%Tbd%phy_myj_chkqlm, wt_h=wt_h)
          call mover(mi, 'phy_myj_elflx', mn_phys%phy_myj_elflx, GFS_data%Tbd%phy_myj_elflx, wt_h=wt_h)
          call mover(mi, 'phy_myj_a1u', mn_phys%phy_myj_a1u, GFS_data%Tbd%phy_myj_a1u, wt_h=wt_h)
          call mover(mi, 'phy_myj_a1t', mn_phys%phy_myj_a1t, GFS_data%Tbd%phy_myj_a1t, wt_h=wt_h)
          call mover(mi, 'phy_myj_a1q', mn_phys%phy_myj_a1q, GFS_data%Tbd%phy_myj_a1q, wt_h=wt_h)
       endif
    ENDIF UNTESTED_IN_LATEST_VERSION
! END OF CONFLICTS----------------------------------------------------------------------------------------

       check_stype_vtype: if(mi%action == DO_COPY_TO_BLOCK_ARRAYS) then
          if (GFS_Control%lsm == GFS_Control%lsm_ruc) then
             where(GFS_Data%Sfcprop%tsnow_land < 200)
                GFS_Data%Sfcprop%tsnow_land = 280
             end where
             where(GFS_Data%Sfcprop%tsnow_ice < 200)
                GFS_Data%Sfcprop%tsnow_ice = 280
             end where
             where(GFS_Data%Sfcprop%tslb < 200)
                GFS_Data%Sfcprop%tslb = 280
             end where
             GFS_data%Sfcprop%flag_frsoil = nint(GFS_data%Sfcprop%flag_frsoil)
          endif
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

    if_move_nsst: if (GFS_control%nstf_name(1) > 0) then
       call mover(mi, 'tref', mn_phys%tref, GFS_Data%Sfcprop%tref, wt_h)
       call mover(mi, 'z_c', mn_phys%z_c, GFS_Data%Sfcprop%z_c, wt_h)
       call mover(mi, 'c_0', mn_phys%c_0, GFS_Data%Sfcprop%c_0, wt_h)
       call mover(mi, 'c_d', mn_phys%c_d, GFS_Data%Sfcprop%c_d, wt_h)
       call mover(mi, 'w_0', mn_phys%w_0, GFS_Data%Sfcprop%w_0, wt_h)
       call mover(mi, 'w_d', mn_phys%w_d, GFS_Data%Sfcprop%w_d, wt_h)
       call mover(mi, 'xt', mn_phys%xt, GFS_Data%Sfcprop%xt, wt_h)
       call mover(mi, 'xs', mn_phys%xs, GFS_Data%Sfcprop%xs, wt_h)
       call mover(mi, 'xu', mn_phys%xu, GFS_Data%Sfcprop%xu, wt_h)
       call mover(mi, 'xv', mn_phys%xv, GFS_Data%Sfcprop%xv, wt_h)
       call mover(mi, 'xz', mn_phys%xz, GFS_Data%Sfcprop%xz, wt_h)
       call mover(mi, 'zm', mn_phys%zm, GFS_Data%Sfcprop%zm, wt_h)
       call mover(mi, 'xtts', mn_phys%xtts, GFS_Data%Sfcprop%xtts, wt_h)
       call mover(mi, 'xzts', mn_phys%xzts, GFS_Data%Sfcprop%xzts, wt_h)
       call mover(mi, 'd_conv', mn_phys%d_conv, GFS_Data%Sfcprop%d_conv, wt_h)
       call mover(mi, 'dt_cool', mn_phys%dt_cool, GFS_Data%Sfcprop%dt_cool, wt_h)
       call mover(mi, 'qrain', mn_phys%qrain, GFS_Data%Sfcprop%qrain, wt_h)
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

  subroutine mover_int_3d(mi, name, var, block_array, wt_h)
    implicit none
    type(movement_info), intent(in) :: mi
    character(len=*), intent(in) :: name
    integer, allocatable, intent(inout) :: var(:,:,:)
    integer, intent(inout) :: block_array(:,:)
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
  end subroutine mover_int_3d

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
