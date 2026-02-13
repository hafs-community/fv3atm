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
!! @brief   Provides subroutines to read subsets of static surface data
!! @author W. Ramstrom, AOML/HRD   07/08/2024
!! @email William.Ramstrom@noaa.gov
! =======================================================================!

module fv_moving_nest_tiled_read_mod

  use mpp_mod,           only: FATAL, WARNING, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod,           only: mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_set_warn_level
  use mpp_mod,           only: mpp_declare_pelist, mpp_set_current_pelist, mpp_sync, mpp_sync_self
  use mpp_mod,           only: input_nml_file
  use mpp_mod,           only: mpp_get_current_pelist, mpp_broadcast
  use mpp_domains_mod,   only: GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, CGRID_NE, DGRID_NE, AGRID
  use mpp_parameter_mod, only: AGRID_PARAM=>AGRID,CGRID_NE_PARAM=>CGRID_NE,SCALAR_PAIR
  use mpp_domains_mod,   only: FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE, FOLD_WEST_EDGE, FOLD_EAST_EDGE
  use mpp_domains_mod,   only: MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE, XUPDATE, YUPDATE, SCALAR_PAIR
  use mpp_domains_mod,   only: domain1D, domain2D, DomainCommunicator2D, BITWISE_EFP_SUM
  use mpp_domains_mod,   only: mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod,   only: mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
  use mpp_domains_mod,   only: mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
  use mpp_domains_mod,   only: mpp_update_domains, mpp_check_field, mpp_redistribute, mpp_get_memory_domain
  use mpp_domains_mod,   only: mpp_define_layout, mpp_define_domains, mpp_modify_domain
  use mpp_domains_mod,   only: mpp_define_io_domain
  use mpp_domains_mod,   only: mpp_get_neighbor_pe, mpp_define_mosaic, mpp_nullify_domain_list
  use mpp_domains_mod,   only: NORTH, NORTH_EAST, EAST, SOUTH_EAST, CORNER, CENTER
  use mpp_domains_mod,   only: SOUTH, SOUTH_WEST, WEST, NORTH_WEST, mpp_define_mosaic_pelist
  use mpp_domains_mod,   only: mpp_get_global_domain, ZERO, NINETY, MINUS_NINETY
  use mpp_domains_mod,   only: mpp_get_boundary, mpp_start_update_domains, mpp_complete_update_domains
  use mpp_domains_mod,   only: mpp_define_nest_domains, nest_domain_type
  use mpp_domains_mod,   only: mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod,   only: mpp_get_F2C_index, mpp_update_nest_coarse
  use mpp_domains_mod,   only: mpp_get_domain_shift, EDGEUPDATE, mpp_deallocate_domain
  use mpp_domains_mod,   only: mpp_group_update_type, mpp_create_group_update
  use mpp_domains_mod,   only: mpp_do_group_update, mpp_clear_group_update
  use mpp_domains_mod,   only: mpp_start_group_update, mpp_complete_group_update
  use mpp_domains_mod,   only: WUPDATE, SUPDATE, mpp_get_compute_domains, NONSYMEDGEUPDATE
  use mpp_domains_mod,   only: domainUG, mpp_define_unstruct_domain, mpp_get_UG_domain_tile_id
  use mpp_domains_mod,   only: mpp_get_UG_compute_domain, mpp_pass_SG_to_UG, mpp_pass_UG_to_SG
  use mpp_domains_mod,   only: mpp_get_ug_global_domain, mpp_global_field_ug
  use mpp_memutils_mod,  only: mpp_memuse_begin, mpp_memuse_end

#ifdef GFS_TYPES
  use GFS_typedefs,      only: kind_phys
#else
  use IPD_typedefs,      only: kind_phys => IPD_kind_phys
#endif

#ifdef OVERLOAD_R4
  use constantsR4_mod,  only: grav
#else
  use constants_mod,     only: grav
#endif
  use boundary_mod,      only: update_coarse_grid, update_coarse_grid_mpp
  use bounding_box_mod,  only: bbox, bbox_get_C2F_index, fill_bbox
  use fms2_io_mod,       only: read_data, write_data, open_file, close_file, register_axis, register_field
  use fms2_io_mod,       only: FmsNetcdfDomainFile_t, FmsNetcdfFile_t
  use fms2_io_mod,       only: is_dimension_registered, get_variable_size, get_variable_num_dimensions, get_dimension_names
  use fms_mod,           only: mpp_clock_begin, mpp_clock_end

  use fv_arrays_mod,     only: R_GRID
  use fv_arrays_mod,     only: fv_grid_type, fv_nest_type, fv_atmos_type
  use fv_surf_map_mod,   only: FV3_zs_filter
  use fv_moving_nest_types_mod, only: grid_geometry, fv_moving_nest_type, mn_surface_grids
  use ifport,            only: getcwd

  implicit none

#ifdef NO_QUAD_PRECISION
  ! 64-bit precision (kind=8)
  integer, parameter:: f_p = selected_real_kind(15)
#else
  ! Higher precision (kind=16) for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

  integer, parameter :: UWIND = 1
  integer, parameter :: VWIND = 2

  logical :: debug_log = .false.


  interface mn_static_read_tiled_hires
    module procedure  mn_static_read_tiled_hires_r4
    module procedure  mn_static_read_tiled_hires_r8
  end interface mn_static_read_tiled_hires

  interface alloc_read_tiled_data
#ifdef OVERLOAD_R8
    module procedure alloc_read_tiled_data_r4_2d
#endif
    module procedure alloc_read_tiled_data_r8_2d
  end interface alloc_read_tiled_data


#include <fms_platform.h>

contains

  !>@brief The subroutine 'mn_static_filename' generates the full pathname for a static file for each run
  !>@details Constructs the full pathname for a variable and refinement level and tests whether it exists
  subroutine mn_static_filename(surface_dir, tile_num, tag, refine, grid_filename)
    character(len=*), intent(in)       :: surface_dir     !< Directory
    character(len=*), intent(in)       :: tag             !< Variable name
    integer, intent(in)                :: tile_num        !< Tile number
    integer, intent(in)                :: refine          !< Nest refinement
    character(len=*), intent(out)      :: grid_filename   !< Output pathname to netCDF file

    character(len=256) :: refine_str, parent_str
    character(len=1)   :: divider
    logical            :: file_exists
    integer            :: this_pe

    this_pe = mpp_pe()

    write(parent_str, '(I0)'), tile_num

    if (refine .eq. 1 .and. (trim(tag) .eq. 'grid' .or. trim(tag) .eq. 'oro_data')) then
      ! For 1x files in INPUT directory; go at the symbolic link
      grid_filename = trim(trim(surface_dir) // '/' // trim(tag) // '.tile' // trim(parent_str) // '.nc')
    else
      if (refine .eq. 1) then
        grid_filename = trim(trim(surface_dir) // '/' // trim(tag) // '.tile' // trim(parent_str) // '.nc')
      else
        write(refine_str, '(I0,A1)'), refine, 'x'
        grid_filename = trim(trim(surface_dir) // '/' // trim(tag) // '.tile' // trim(parent_str) // '.' // trim(refine_str) // '.nc')
      endif
    endif

    grid_filename = trim(grid_filename)

    inquire(FILE=grid_filename, EXIST=file_exists)
    if (.not. file_exists) then
      !call mpp_error(FATAL, 'mn_static_filename DOES NOT EXIST '//trim(grid_filename))
      print '("[ERROR] WDR mn_static_filename DOES NOT EXIST npe=",I0," tile_num=",I0," tag=",A16," refine=",I0," grid_filename=",A120)', this_pe, tile_num, tag, refine, grid_filename
      print '("[ERROR] WDR mn_static_filename DNE npe=",I0," grid_filename=",A120)', this_pe, grid_filename
    endif

  end subroutine mn_static_filename


  subroutine compare_tile_grids(tile_grid, full_grid, var_name)
    real, _ALLOCATABLE, intent(in)   :: tile_grid(:,:)  !< 2D grid of data
    real, _ALLOCATABLE, intent(in)   :: full_grid(:,:)  !< 2D grid of data
    character(len=*), intent(in)     :: var_name                   !< Variable name in netCDF file

    integer :: tis, tie, tjs, tje
    integer :: fis, fie, fjs, fje
    integer :: i,j, num_matches, num_mismatches
    integer :: this_pe

    this_pe = mpp_pe()

    if (.not. allocated(tile_grid)) then
      print '("[ERROR] WDR compare_tile_grids npe=",I0," ",A32," tile_grid not allocated.")', this_pe, var_name
      return
    endif
    if (.not. allocated(full_grid)) then
      print '("[ERROR] WDR compare_tile_grids npe=",I0," ",A32," full_grid not allocated.")', this_pe, var_name
      return
    endif

    tis = lbound(tile_grid,1)
    tie = ubound(tile_grid,1)
    tjs = lbound(tile_grid,2)
    tje = ubound(tile_grid,2)

    fis = lbound(full_grid,1)
    fie = ubound(full_grid,1)
    fjs = lbound(full_grid,2)
    fje = ubound(full_grid,2)

    !print '("[WDR] compare_tile_grids npe=",I0," TILE BOUNDS ",A32,"(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, var_name, tis, tie, tjs, tje
    !print '("[WDR] compare_tile_grids npe=",I0," FULL BOUNDS ",A32,"(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, var_name, fis, fie, fjs, fje

    num_mismatches = 0
    num_matches = 0

    if (tis .lt. fis .or. tie .gt. fie .or. tjs .lt. fjs .or. tje .gt. fje) then
      print '("[ERROR] WDR compare_tile_grids npe=",I0," ",A32," tile_grid outside of full_grid.")', this_pe, var_name
      return
    endif

    do i=tis, tie
      do j=tjs, tje
        if (tile_grid(i,j) .eq. full_grid(i,j)) then
          num_matches = num_matches + 1
          !print '("[INFO] WDR compare_tile_grids npe=",I0," MATCH ",A32," tile_grid(",I0,",",I0,")=",F10.4," full_grid=",F10.4)', this_pe, var_name, i, j, tile_grid(i,j), full_grid(i,j)
        else
          num_mismatches = num_mismatches + 1
          print '("[ERROR] WDR compare_tile_grids npe=",I0," MISMATCH ",A32," tile_grid(",I0,",",I0,")=",F10.4," full_grid=",F10.4)', this_pe, var_name, i, j, tile_grid(i,j), full_grid(i,j)
        endif
      enddo
    enddo

    print '("[INFO] WDR compare_tile_grids npe=",I0," total ",A32," num_matches=",I0," num_mismatches=",I0)', this_pe, var_name, num_matches, num_mismatches

  end subroutine compare_tile_grids

  !>@brief The subroutine 'mn_replace_low_values' replaces low values with a default value.
  subroutine mn_replace_low_values(data_grid, low_value, new_value)
    real, _ALLOCATABLE, intent(inout)   :: data_grid(:,:)  !< 2D grid of data
    real, intent(in)                    :: low_value       !< Low value to check for; e.g. negative or fill value
    real, intent(in)                    :: new_value       !< Value to replace low value with

    integer :: i, j

    do i=lbound(data_grid,1),ubound(data_grid,1)
      do j=lbound(data_grid,2),ubound(data_grid,2)
        if (data_grid(i,j) .le. low_value) data_grid(i,j) = new_value
      enddo
    enddo
  end subroutine mn_replace_low_values


  subroutine initialize_static_tile_bounds(st, npx, npy, refine, nest_nx, nest_ny, ratio)
    type(mn_surface_grids), intent(inout) :: st
    integer, intent(in)                   :: npx, npy
    integer, intent(in)                   :: refine
    integer, intent(in)                   :: nest_nx, nest_ny
    real, intent(in)                      :: ratio

    ! This sets the sizes of the parent and tile for the static surface datasets
    !   These will not change during the model run
    !   This routine does not set the tile offsets, as they will change during the run

    call calc_fp_bounds(npx, npy, refine, st%fp_nx, st%fp_ny)

    call calc_tile_bounds(st%fp_nx, st%fp_ny, nest_nx, nest_ny, ratio, st%tile_nx, st%tile_ny)

    st%num_reads = 0

  end subroutine initialize_static_tile_bounds


  !>@brief The subroutine calc_fp_bounds calculates the full panel size,
  !>@where the full panel is the parent grid at nest resolution
  subroutine calc_fp_bounds(npx, npy, refine, fp_nx, fp_ny)
    integer, intent(in)   :: npx, npy
    integer, intent(in)   :: refine
    integer, intent(out)  :: fp_nx, fp_ny

    integer :: nx, ny

    nx = npx - 1
    ny = npy - 1

    fp_nx = nx * refine
    fp_ny = ny * refine

  end subroutine calc_fp_bounds


  subroutine calc_tile_bounds(fp_nx, fp_ny, nest_nx, nest_ny, ratio, tile_nx, tile_ny)
    integer, intent(in)   :: fp_nx, fp_ny
    integer, intent(in)   :: nest_nx, nest_ny
    real, intent(in)      :: ratio
    integer, intent(out)  :: tile_nx, tile_ny

    tile_nx = nest_nx + nint( (fp_nx - nest_nx) * ratio)
    tile_ny = nest_ny + nint( (fp_ny - nest_ny) * ratio)

    ! Safeguard -- don't let the tile be larger than the high-resolution parent.
    if (tile_nx .gt. fp_nx) tile_nx = fp_nx
    if (tile_ny .gt. fp_ny) tile_ny = fp_ny


  end subroutine calc_tile_bounds


  logical function is_grid_inside_tile(nest_nx, nest_ny, ioffset, joffset, refine, halo, tile_nx, tile_ny, tile_ioffset, tile_joffset)
    integer, intent(in)   :: nest_nx, nest_ny
    integer, intent(in)   :: ioffset, joffset
    integer, intent(in)   :: refine, halo
    integer, intent(in)   :: tile_nx, tile_ny
    integer, intent(in)   :: tile_ioffset, tile_joffset

    integer :: nsx, nex, nsy, ney
    integer :: tsx, tex, tsy, tey

    integer :: this_pe

    this_pe = mpp_pe()

    ! Nest runs in x direction from ioffset to ioffset + nest_nx
    nsx = (ioffset -1 ) * refine - halo
    nex = (ioffset -1 ) * refine + nest_nx + halo
    ! Nest runs in y direction from joffset to joffset + nest_ny
    nsy = (joffset - 1 ) * refine - halo
    ney = (joffset - 1 ) * refine + nest_ny + halo

    ! Tile runs in x direction from tile_ioffset to tile_ioffset + tile_nx
    tsx = tile_ioffset
    tex = tile_ioffset + tile_nx - 1
    ! Tile runs in y direction from tile_joffset to tile_joffset + tile_ny
    tsy = tile_joffset
    tey = tile_joffset + tile_ny - 1

    is_grid_inside_tile = .False.

    if (nsx .ge. tsx .and. nex .le. tex) then
      if (nsy .ge. tsy .and. ney .le. tey) then
        is_grid_inside_tile = .True.
      endif
    endif

    !if (this_pe .eq. 119) then
    !  print '("[INFO] WDR is_grid_inside_tile npe=",I0," tile(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, tsx, tex, tsy, tey
    !  print '("[INFO] WDR is_grid_inside_tile npe=",I0," nest(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, nsx, nex, nsy, ney
    !  print '("[INFO] WDR is_grid_inside_tile npe=",I0," nest nx=",I0," ny=",I0," ioffset=",I0," joffset=",I0)', this_pe, nest_ny, nest_ny, ioffset, joffset
    !  print '("[INFO] WDR is_grid_inside_tile npe=",I0," tile nx=",I0," ny=",I0," ioffset=",I0," joffset=",I0)', this_pe, tile_ny, tile_ny, tile_ioffset, tile_joffset
    !endif

  end function is_grid_inside_tile

  logical function is_grid_inside_full_tile(nest_nx, nest_ny, ioffset, joffset, tile_nx, tile_ny, tile_ioffset, tile_joffset)
    integer, intent(in)   :: nest_nx, nest_ny
    integer, intent(in)   :: ioffset, joffset
    integer, intent(in)   :: tile_nx, tile_ny
    integer, intent(in)   :: tile_ioffset, tile_joffset

    integer :: nsx, nex, nsy, ney
    integer :: tsx, tex, tsy, tey

    integer :: this_pe

    this_pe = mpp_pe()

    ! Nest runs in x direction from ioffset to ioffset + nest_nx
    nsx = ioffset
    nex = ioffset + nest_nx
    ! Nest runs in y direction from joffset to joffset + nest_ny
    nsy = joffset
    ney = joffset + nest_ny

    ! Tile runs in x direction from tile_ioffset to tile_ioffset + tile_nx
    tsx = tile_ioffset
    tex = tile_ioffset + tile_nx
    ! Tile runs in y direction from tile_joffset to tile_joffset + tile_ny
    tsy = tile_joffset
    tey = tile_joffset + tile_ny

    is_grid_inside_full_tile = .False.

    if (nsx .ge. tsx .and. nex .le. tex) then
      if (nsy .ge. tsy .and. ney .le. tey) then
        is_grid_inside_full_tile = .True.
      endif
    endif

    !if (this_pe .eq. 119) then
    !  print '("[INFO] WDR is_grid_inside_full_tile npe=",I0," tile(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, tsx, tex, tsy, tey
    !  print '("[INFO] WDR is_grid_inside_full_tile npe=",I0," nest(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, nsx, nex, nsy, ney
    !  print '("[INFO] WDR is_grid_inside_full_tile npe=",I0," nest nx=",I0," ny=",I0," ioffset=",I0," joffset=",I0)', this_pe, nest_ny, nest_ny, ioffset, joffset
    !  print '("[INFO] WDR is_grid_inside_full_tile npe=",I0," tile nx=",I0," ny=",I0," ioffset=",I0," joffset=",I0)', this_pe, tile_ny, tile_ny, tile_ioffset, tile_joffset
    !endif

  end function is_grid_inside_full_tile


  subroutine get_nest_grid_center(nx, ny, ioffset, joffset, refine, center_x, center_y)
    integer, intent(in)  :: nx, ny, ioffset, joffset, refine
    integer, intent(out) :: center_x, center_y

    center_x = nx/2 + ( ioffset - 1 ) * refine
    center_y = ny/2 + ( joffset - 1 ) * refine

  end subroutine get_nest_grid_center

  subroutine get_grid_edges_from_center(nx, ny, center_x, center_y, min_x, max_x, min_y, max_y)
    integer, intent(in)  :: nx, ny, center_x, center_y
    integer, intent(out) :: min_x, max_x, min_y, max_y


    integer :: half_nx, half_ny

!    integer :: this_pe
!    this_pe = mpp_pe()

    half_nx = nx/2
    half_ny = ny/2

    min_x = center_x - half_nx
    max_x = center_x + half_nx

    min_y = center_y - half_ny
    max_y = center_y + half_ny

!    if (this_pe .eq. 1199) print '("[INFO] WDR FX grid_edges npe=",I0," nx=",I0," ny=",I0," half_nx=",I0," half_ny=",I0," min_x=",I0," max_x=",I0," min_y=",I0," max_y=",I0)',this_pe, nx, ny, half_nx, half_ny, min_x, max_x, min_y, max_y

  end subroutine get_grid_edges_from_center

  subroutine get_grid_offset_from_center(nx, ny, center_x, center_y, ioffset, joffset)
    integer, intent(in)  :: nx, ny, center_x, center_y
    integer, intent(out) :: ioffset, joffset

    integer :: half_nx, half_ny

!    integer :: this_pe
!    this_pe = mpp_pe()

    half_nx = nx/2
    half_ny = ny/2

    ioffset = center_x - half_nx
    joffset = center_y - half_ny

!    if (this_pe .eq. 1199) print '("[INFO] WDR FX grid_offsets npe=",I0," nx=",I0," ny=",I0," half_nx=",I0," half_ny=",I0," ioffset=",I0," joffset=",I0)',this_pe, nx, ny, half_nx, half_ny, ioffset, joffset

  end subroutine get_grid_offset_from_center




  !>@brief The subroutine trigger_reread_static_data checks if the static data should be reread
  !>@ if so, it returns the new tile offset values that should be read.  Assume for now that tile_nx and tile_ny will be the same through a full run
  subroutine trigger_reread_static_data(fp_nx, fp_ny, nest_nx, nest_ny, tile_nx, tile_ny, nest_ioffset, nest_joffset, tile_ioffset, tile_joffset, refine, halo, first_call, allow_early_read, do_read, new_tile_ioffset, new_tile_joffset)
    integer, intent(in)   :: fp_nx, fp_ny
    integer, intent(in)   :: nest_nx, nest_ny
    integer, intent(in)   :: tile_nx, tile_ny
    integer, intent(in)   :: refine, halo
    integer, intent(in)   :: nest_ioffset, nest_joffset
    integer, intent(in)   :: tile_ioffset, tile_joffset
    logical, intent(in)   :: first_call
    logical, intent(in)   :: allow_early_read    ! Not yet implemented;  both values will produce same effect for now

    logical, intent(out)  :: do_read
    integer, intent(out)  :: new_tile_ioffset, new_tile_joffset

    logical :: is_inside
    integer :: nest_center_x, nest_center_y
    integer :: tile_center_x, tile_center_y
    integer :: x_edge, y_edge
    logical :: is_tile_inside_fp
    integer :: this_pe

    ! allow_early_read set to true will cause the tile to be reread if the center has shifted at all
    !  the idea is that other processors are already busy moving the nest, so we may want to take advantage of the
    !  ununsed cycles on this processor to realign the tile.

    this_pe = mpp_pe()

    if (first_call) then
      is_inside = .false.
    else
      do_read = .False.
      is_inside = is_grid_inside_tile(nest_nx, nest_ny, nest_ioffset, nest_joffset, refine, halo, tile_nx, tile_ny, tile_ioffset, tile_joffset)
    endif

    if (.not. is_inside) then

      do_read = .True.

      ! Move the tile area to be centered on the new nest location.
      !  This may be a big jump of location, many points -- we don't usually want to reread each nest move
      !  Choose new location, but then back off if it is over the edge of the parent domain

      call get_nest_grid_center(nest_nx, nest_ny, nest_ioffset, nest_joffset, refine, nest_center_x, nest_center_y)

      tile_center_x = nest_center_x
      tile_center_y = nest_center_y

      call get_grid_offset_from_center(tile_nx, tile_ny, tile_center_x, tile_center_y, new_tile_ioffset, new_tile_joffset)

      ! Adjust if the tile runs off the edge of the full panel (parent high-resolution grid)

      if (new_tile_ioffset .lt. 1) new_tile_ioffset = 1
      if (new_tile_joffset .lt. 1) new_tile_joffset = 1

      x_edge = (new_tile_ioffset + tile_nx) - fp_nx
      if (x_edge .gt. 0) new_tile_ioffset = new_tile_ioffset - x_edge

      y_edge = (new_tile_joffset + tile_ny) - fp_ny
      if (y_edge .gt. 0) new_tile_joffset = new_tile_joffset - y_edge


      is_tile_inside_fp = is_grid_inside_full_tile(tile_nx, tile_ny, new_tile_ioffset, new_tile_joffset, fp_nx, fp_ny, 0, 0)

      if (.not. is_tile_inside_fp) then
        print '("[ERROR] WDR trigger_reread_static_data out of bounds npe=",I0," new_tile_ioffset=",I0," new_tile_joffset=",I0)', this_pe, new_tile_ioffset, new_tile_joffset

      endif

    endif

  end subroutine trigger_reread_static_data


  subroutine check_update_static_tile_data(fp_nx, fp_ny, nest_nx, nest_ny, ioffset, joffset, refine, a_step, st, surface_dir, pelist, parent_tile, month, use_timers, id_movnest_readstatic, do_read)
    integer, intent(in)                               :: fp_nx, fp_ny, nest_nx, nest_ny, ioffset, joffset, refine, a_step
    type(mn_surface_grids), intent(inout) :: st
    character(len=*), intent(in)       :: surface_dir
    integer, intent(in)                :: pelist(:)
    integer, intent(in)                :: parent_tile
    integer, intent(in)                :: month
    logical, intent(in)                :: use_timers
    integer, intent(in)                :: id_movnest_readstatic
    logical, intent(out)               :: do_read

    logical :: allow_early_read = .False.

    integer :: tile_nx, tile_ny, tile_ioffset, tile_joffset
    integer :: new_tile_ioffset, new_tile_joffset
    integer :: num_bytes, total_bytes
    integer :: i,j
    integer :: this_pe
    integer :: halo = 3

    logical, save     :: first_call = .true.

    this_pe = mpp_pe()

    total_bytes = 0

    !type(fv_moving_nest_type), target ,intent(inout)  :: child_moving_nest
    !type(mn_surface_grids), pointer :: st
    !surface_dir = trim(child_moving_nest%mn_flag%surface_dir)
    !st => child_moving_nest%mn_static

    tile_nx = st%tile_nx
    tile_ny = st%tile_ny
    tile_ioffset = st%tile_ioffset
    tile_joffset = st%tile_joffset

    !print '("[INFO] WDR TC1 check_update_static_tile_data npe=",I0)', this_pe

    !print '("[INFO] WDR check_update_static_tile_data npe=",I0," tile_nx=",I0," tile_ny=",I0," tile_ioffset=",I0," tile_joffset=",I0," first_call=",L1)',this_pe, tile_nx, tile_ny, tile_ioffset, tile_joffset, first_call

    call trigger_reread_static_data(fp_nx, fp_ny, nest_nx, nest_ny, tile_nx, tile_ny, ioffset, joffset, tile_ioffset, tile_joffset, refine, halo, first_call, allow_early_read, do_read, new_tile_ioffset, new_tile_joffset)
    if (first_call) first_call=.False.

    !print '("[INFO] WDR check_update_static_tile_data npe=",I0," do_read=",L1)', this_pe, do_read
    !print '("[INFO] WDR TC2 check_update_static_tile_data npe=",I0)', this_pe

    if (do_read) then
      !print '("[INFO] WDR TC3 check_update_static_tile_data npe=",I0,L1,L1,L1,L1)', this_pe, allocated(st%fp_fix%deep_soil_temp_grid), allocated(st%fp_fix%slope_type_grid), allocated(st%fp_ls%soil_type_grid), allocated(st%nest_fix%deep_soil_temp_grid)
      !print '("[INFO] WDR TILE READ check_update_static_tile_data npe=",I0," a_step=",I0," tile_nx=",I0," tile_ny=",I0," new_tile_ioffset=",I0," new_tile_joffset=",I0)',this_pe, a_step, tile_nx, tile_ny, new_tile_ioffset, new_tile_joffset
      if (use_timers) call mpp_clock_begin (id_movnest_readstatic)

      !print '("[INFO] WDR TILE READ check npe=",I0," trim(surface_dir)=",A120)', this_pe, trim(surface_dir)
      !print '("[INFO] WDR TILE READ check npe=",I0," surface_dir=",A120)', this_pe, surface_dir
      !print '("[INFO] WDR TILE READ check npe=",I0," alloc deep_soil=",L1)', this_pe, allocated(st%fp_fix%deep_soil_temp_grid)
      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "substrate_temperature", "substrate_temperature", st%fp_fix%deep_soil_temp_grid, num_bytes, parent_tile)
      total_bytes = total_bytes + num_bytes
      ! set any -999s to +4C
      call mn_replace_low_values(st%fp_fix%deep_soil_temp_grid, -100.0, 277.0)

      !print '("[INFO] WDR TC4 check_update_static_tile_data npe=",I0)', this_pe

      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "substrate_temperature", "geolat", st%fp_fix%deep_lat, num_bytes, parent_tile)
      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "substrate_temperature", "geolon", st%fp_fix%deep_lon, num_bytes, parent_tile)

      !print '("[INFO] WDR TC5 check_update_static_tile_data npe=",I0)', this_pe

!      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "soil_type", "soil_type", st%fp_ls%soil_type_grid, num_bytes, parent_tile)
      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "soil_type", "soil_type", st%fp_ls%soil_type_grid, num_bytes, parent_tile)
      total_bytes = total_bytes + num_bytes
      !print '("[INFO] WDR TC5.1 check_update_static_tile_data npe=",I0)', this_pe
      ! To match initialization behavior, set any -999s to 0 in soil_type
      call mn_replace_low_values(st%fp_ls%soil_type_grid, -100.0, 0.0)
      !print '("[INFO] WDR TC5.2 check_update_static_tile_data npe=",I0)', this_pe


      !! TODO investigate reading high-resolution veg_frac and veg_greenness
      !call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "", mn_static%veg_frac_grid)

      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "vegetation_type", "vegetation_type", st%fp_fix%veg_type_grid, num_bytes, parent_tile)
      total_bytes = total_bytes + num_bytes

      ! To match initialization behavior, set any -999s to 0 in veg_type
      call mn_replace_low_values(st%fp_fix%veg_type_grid, -100.0, 0.0)


      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "slope_type", "slope_type", st%fp_fix%slope_type_grid, num_bytes, parent_tile)
      total_bytes = total_bytes + num_bytes

      ! To match initialization behavior, set any -999s to 0 in slope_type
      call mn_replace_low_values(st%fp_fix%slope_type_grid, -100.0, 0.0)


      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "maximum_snow_albedo", "maximum_snow_albedo", st%fp_fix%max_snow_alb_grid, num_bytes, parent_tile)
      total_bytes = total_bytes + num_bytes

      ! Set any -999s to 0.5
      call mn_replace_low_values(st%fp_fix%max_snow_alb_grid, -100.0, 0.5)

      ! Albedo fraction -- read and calculate
      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "facsf", "facsf", st%fp_fix%facsf_grid, num_bytes, parent_tile)
      total_bytes = total_bytes + num_bytes

      if (allocated(st%fp_fix%facwf_grid)) deallocate(st%fp_fix%facwf_grid)

      allocate(st%fp_fix%facwf_grid(lbound(st%fp_fix%facsf_grid,1):ubound(st%fp_fix%facsf_grid,1),lbound(st%fp_fix%facsf_grid,2):ubound(st%fp_fix%facsf_grid,2)))
      total_bytes = total_bytes + num_bytes

      ! For land points, set facwf = 1.0 - facsf
      ! To match initialization behavior, set any -999s to 0
      do i=lbound(st%fp_fix%facsf_grid,1),ubound(st%fp_fix%facsf_grid,1)
        do j=lbound(st%fp_fix%facsf_grid,2),ubound(st%fp_fix%facsf_grid,2)
          if (st%fp_fix%facsf_grid(i,j) .lt. -100) then
            st%fp_fix%facsf_grid(i,j) = 0
            st%fp_fix%facwf_grid(i,j) = 0
          else
            st%fp_fix%facwf_grid(i,j) = 1.0 - st%fp_fix%facsf_grid(i,j)
          endif
        enddo
      enddo

      ! Additional albedo variables
      !  black sky = strong cosz -- direct sunlight
      !  white sky = weak cosz -- diffuse light

      ! alvsf = visible strong cosz = visible_black_sky_albedo
      ! alvwf = visible weak cosz = visible_white_sky_albedo
      ! alnsf = near IR strong cosz = near_IR_black_sky_albedo
      ! alnwf = near IR weak cosz = near_IR_white_sky_albedo

      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "snowfree_albedo", "visible_black_sky_albedo", st%fp_fix%alvsf_grid, num_bytes, parent_tile, time=month)
      total_bytes = total_bytes + num_bytes

      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "snowfree_albedo", "visible_white_sky_albedo", st%fp_fix%alvwf_grid, num_bytes, parent_tile, time=month)
      total_bytes = total_bytes + num_bytes

      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "snowfree_albedo", "near_IR_black_sky_albedo", st%fp_fix%alnsf_grid, num_bytes, parent_tile, time=month)
      total_bytes = total_bytes + num_bytes

      call mn_static_read_tiled_hires_r4(new_tile_ioffset, new_tile_joffset, tile_nx, tile_ny, refine, pelist, surface_dir, "snowfree_albedo", "near_IR_white_sky_albedo", st%fp_fix%alnwf_grid, num_bytes, parent_tile, time=month)
      total_bytes = total_bytes + num_bytes

      ! Set the -999s to small value of 0.06, matching initialization code in chgres

      call mn_replace_low_values(st%fp_fix%alvsf_grid, -100.0, 0.06)
      call mn_replace_low_values(st%fp_fix%alvwf_grid, -100.0, 0.06)
      call mn_replace_low_values(st%fp_fix%alnsf_grid, -100.0, 0.06)
      call mn_replace_low_values(st%fp_fix%alnwf_grid, -100.0, 0.06)

      ! Complete the handling of a tile refresh

      st%num_reads = st%num_reads + 1

      st%tile_ioffset = new_tile_ioffset
      st%tile_joffset = new_tile_joffset

      print '("[INFO] WDR check_update_static_tile_data READ ",I0," COMPLETE npe=",I0," complete bytes=",I0," MB=",F10.3)', st%num_reads, this_pe, total_bytes, total_bytes/1024.0/1024.0

      if (use_timers) call mpp_clock_end (id_movnest_readstatic)

    endif

!    print '("[INFO] WDR TILE A4a check_update_static_tile_data npe=",I0," allocated(st%deep_soil_temp_grid)=",L1," allocated(child_moving_nest%mn_static%deep_soil_temp_grid)=",L1)', this_pe, allocated(st%deep_soil_temp_grid), allocated(child_moving_nest%mn_static%deep_soil_temp_grid)
!    print '("[INFO] WDR TILE A4a check_update_static_tile_data npe=",I0," allocated(st%deep_soil_temp_grid)=",L1)', this_pe, allocated(st%deep_soil_temp_grid)

  end subroutine check_update_static_tile_data







  !>@brief The subroutine 'mn_static_read_tiled_hires_r4' loads high resolution data from netCDF
  !>@details Gathers a single variable from the netCDF file
  subroutine mn_static_read_tiled_hires_r4(tile_ioffset, tile_joffset, x_size, y_size, refine, pelist, surface_dir, file_prefix, var_name, data_grid, num_bytes, parent_tile, time)
    integer, intent(in)                :: tile_ioffset, tile_joffset
    integer, intent(in)                :: x_size, y_size, refine
    integer, intent(in)                :: pelist(:)                  !< PE list for fms2_io
    character(len=*), intent(in)       :: surface_dir, file_prefix   !< Surface directory and file tag
    character(len=*), intent(in)       :: var_name                   !< Variable name in netCDF file
    real*4, allocatable, intent(out)   :: data_grid(:,:)             !< Output data grid
    integer, intent(out)               :: num_bytes
    integer, intent(in)                :: parent_tile                !< Parent tile number
    integer, intent(in), optional      :: time                       !< Optional month number for time-varying parameters

    character(len=512) :: nc_filename
    integer            :: this_pe

    this_pe = mpp_pe()

    !call mpp_sync(pelist)
!    if (present(time)) then
!      print '("[INFO] WDR SRT0.1 npe=",I0," time PRESENT time=",I0, " var_name=",A)', this_pe, time, var_name
!    else
!      print '("[INFO] WDR SRT0.1 npe=",I0," time NOT PRESENT var_name=",A)', this_pe, var_name
!    endif


    if (allocated(data_grid)) then
      !print '("[INFO] WDR mn_static_read_tiled_hires_r4 npe=",I0," start ",A24," with bounds ",I0,"-",I0,",",I0,"-",I0)', this_pe, var_name, lbound(data_grid,1) , ubound(data_grid,1) , lbound(data_grid,2) , ubound(data_grid,2)
      deallocate(data_grid)
    !else
    !  print '("[INFO] WDR mn_static_read_tiled_hires_r4 npe=",I0," start ",A24," unallocated.")', this_pe, var_name
    endif

    call mn_static_filename(surface_dir, parent_tile, file_prefix, refine, nc_filename)

    if (present(time)) then
      call alloc_read_tiled_data(nc_filename, var_name, tile_ioffset, tile_joffset, x_size, y_size, data_grid, pelist, time)
    else
      call alloc_read_tiled_data(nc_filename, var_name, tile_ioffset, tile_joffset, x_size, y_size, data_grid, pelist)
    endif

    num_bytes = sizeof(data_grid)

    !print '("[INFO] WDR mn_static_read_tiled_hires_r4 npe=",I0," complete ",A24," with bounds ",I0,"-",I0,",",I0,"-",I0," bytes=",I0)', this_pe, var_name, lbound(data_grid,1) , ubound(data_grid,1) , lbound(data_grid,2) , ubound(data_grid,2), num_bytes
  end subroutine mn_static_read_tiled_hires_r4

  subroutine mn_static_read_tiled_hires_r8(tile_ioffset, tile_joffset, x_size, y_size, refine, pelist, surface_dir, file_prefix, var_name, data_grid, num_bytes, parent_tile, time)
    integer, intent(in)                :: tile_ioffset, tile_joffset
    integer, intent(in)                :: x_size, y_size, refine
    integer, allocatable, intent(in)   :: pelist(:)                  !< PE list for fms2_io
    character(len=*), intent(in)       :: surface_dir, file_prefix   !< Surface directory and file tag
    character(len=*), intent(in)       :: var_name                   !< Variable name in netCDF file
    real*8, allocatable, intent(out)   :: data_grid(:,:)             !< Output data grid
    integer, intent(out)               :: num_bytes
    integer, intent(in)                :: parent_tile                !< Parent tile number
    integer, intent(in), optional      :: time                       !< Optional month number for time-varying parameters

    character(len=512) :: nc_filename
    integer            :: this_pe

    this_pe = mpp_pe()


    if (allocated(data_grid)) then
      !print '("[INFO] WDR mn_static_read_tiled_hires_r8 npe=",I0," start ",A24," with bounds ",I0,"-",I0,",",I0,"-",I0)', this_pe, var_name, lbound(data_grid,1) , ubound(data_grid,1) , lbound(data_grid,2) , ubound(data_grid,2)
      deallocate(data_grid)
    !else
    !  print '("[INFO] WDR mn_static_read_tiled_hires_r8 npe=",I0," start ",A24," unallocated.")', this_pe, var_name
    endif

    print '("[INFO] WDR call mn_static_filename R2 npe=",I0)', mpp_pe()

    call mn_static_filename(surface_dir, parent_tile, file_prefix, refine, nc_filename)

    if (present(time)) then
      call alloc_read_tiled_data(nc_filename, var_name, tile_ioffset, tile_joffset, x_size, y_size, data_grid, pelist, time)
    else
      call alloc_read_tiled_data(nc_filename, var_name, tile_ioffset, tile_joffset, x_size, y_size, data_grid, pelist)
    endif

    num_bytes = sizeof(data_grid)

    !print '("[INFO] WDR mn_static_read_tiled_hires_r8 npe=",I0," complete ",A24," with bounds ",I0,"-",I0,",",I0,"-",I0," bytes=",I0," MB=",F10.3)', this_pe, var_name, lbound(data_grid,1) , ubound(data_grid,1) , lbound(data_grid,2) , ubound(data_grid,2), num_bytes, num_bytes/1024.0/1024.0

  end subroutine mn_static_read_tiled_hires_r8



#ifdef OVERLOAD_R8
  subroutine alloc_read_tiled_data_r4_2d(nc_filename, var_name, tile_ioffset, tile_joffset, x_size, y_size, data_array, pes, time)
    character(len=*), intent(in)           :: nc_filename, var_name
    integer, intent(in)                    :: tile_ioffset, tile_joffset
    integer, intent(in)                    :: x_size, y_size
    real*4, allocatable, intent(inout)     :: data_array(:,:)
    integer, intent(in)                    :: pes(:)
    integer, intent(in), optional          :: time

    type(FmsNetcdfFile_t)        :: fileobj        !< Fms2_io fileobj
    real*4, allocatable          :: time_array(:,:,:)

    integer, dimension(3)        :: corner
    integer, dimension(3)        :: edge_lengths
    integer, dimension(3)        :: dim_sizes
    integer :: num_dims
    character(len=32), dimension(3) :: dim_names

    integer                      :: this_pe

    !real*4, allocatable :: local_data_array(:,:,:)

    ! Allocate data_array to match the expected data size, then read in the data
    ! This subroutine consolidates the allocation and reading of data to ensure consistency of data sizing and simplify code
    ! Could later extend this function to determine data size based on netCDF file metadata

    this_pe = mpp_pe()
    dim_sizes = -1

    corner(1) = tile_ioffset
    corner(2) = tile_joffset
    corner(3) = 0

    edge_lengths(1) = x_size
    edge_lengths(2) = y_size
    edge_lengths(3) = 1

    !allocate(data_array(x_size, y_size))
    ! Seems like we can allocate with lower and upper bounds, and return it with bounds as it is an allocatable argument
    !  But passed into read_data, without the allocatable designation, it reverts to lower bound of 0, which FMS expects
    !allocate(local_data_array(1, tile_ioffset:tile_ioffset+x_size, tile_joffset:tile_joffset+y_size))
    !local_data_array = -9999.9

    allocate(data_array(tile_ioffset:tile_ioffset+x_size-1, tile_joffset:tile_joffset+y_size-1))
    data_array = -9999.9

    if (present(time)) then
      corner(3) = time
      edge_lengths(3) = 1

      !allocate(time_array(x_size, y_size, 12)) ! assume monthly data; allocate 12 slots
      allocate(time_array(time:time,tile_ioffset:tile_ioffset+x_size-1, tile_joffset:tile_joffset+y_size-1))
      time_array = 0.00001

      if (open_file(fileobj, nc_filename, "read", pelist=pes, is_restart=.false.)) then
        !print '("[INFO] WDR alloc_read_tiled_data_r4_2d DDTIME npe=",I0)', this_pe
        !call get_dimension_names(fileobj, dim_names)
        !print '("[INFO] WDR alloc_read_tiled_data_r4_2d DIMNAMES npe=",I0," ",A32,"(",A8,",",A8,",",A8,")")', this_pe, var_name, dim_names(1), dim_names(2), dim_names(3)
        !call get_variable_size(fileobj, var_name, dim_sizes)
        !print '("[INFO] WDR alloc_read_tiled_data_r4_2d DIMSIZES npe=",I0," ",A32,"(",I0,",",I0,",",I0,")")', this_pe, var_name, dim_sizes(1), dim_sizes(2), dim_sizes(3)

        call read_data(fileobj, var_name, time_array, corner=corner, edge_lengths=edge_lengths)
        call close_file(fileobj)
      endif

      data_array = time_array(time,:,:)
      deallocate(time_array)
    else
      ! Following transition documents at https://github.com/NOAA-GFDL/FMS/tree/2021.03.01/fms2_io
      if (open_file(fileobj, nc_filename, "read", pelist=pes, is_restart=.false.)) then
        !num_dims = get_variable_num_dimensions(fileobj, var_name)
        !call get_variable_size(fileobj, var_name, dim_sizes(1:num_dims))
        !print '("[INFO] WDR alloc_read_tiled_data_r4_2d DIMSIZES npe=",I0," ",A32,"(",I0,",",I0,",",I0,")")', this_pe, var_name, dim_sizes(1), dim_sizes(2), dim_sizes(3)

        call read_data(fileobj, var_name, data_array, corner=corner(1:2), edge_lengths=edge_lengths(1:2))
        call close_file(fileobj)
      endif
    endif

    !data_array(:,:) = local_data_array(1,:,:)
    !print '("[INFO] WDR alloc_read_tiled_data_r4_2d ZZ npe=",I0)', this_pe

  end subroutine alloc_read_tiled_data_r4_2d
#endif

  subroutine alloc_read_tiled_data_r8_2d(nc_filename, var_name, tile_ioffset, tile_joffset, x_size, y_size, data_array, pes, time)
    character(len=*), intent(in)           :: nc_filename, var_name
    integer, intent(in)                    :: tile_ioffset, tile_joffset
    integer, intent(in)                    :: x_size, y_size
    real*8, allocatable, intent(inout)     :: data_array(:,:)
    integer, intent(in)                    :: pes(:)
    integer, intent(in), optional          :: time

    type(FmsNetcdfFile_t)        :: fileobj        !< Fms2_io fileobj
    real*4, allocatable          :: time_array(:,:,:)

    integer, dimension(3)        :: corner
    integer, dimension(3)        :: edge_lengths


    integer                      :: this_pe

    ! Allocate data_array to match the expected data size, then read in the data
    ! This subroutine consolidates the allocation and reading of data to ensure consistency of data sizing and simplify code
    ! Could later extend this function to determine data size based on netCDF file metadata

    this_pe = mpp_pe()

    corner(1) = tile_ioffset
    corner(2) = tile_joffset
    corner(3) = 0

    edge_lengths(1) = x_size
    edge_lengths(2) = y_size
    edge_lengths(3) = 0

    !allocate(data_array(x_size, y_size))
    ! Seems like we can allocate with lower and upper bounds, and return it with bounds as it is an allocatable argument
    !  But passed into read_data, without the allocatable designation, it reverts to lower bound of 0, which FMS expects
    ! TODO match this with the r4 expression
    allocate(data_array(tile_ioffset:tile_ioffset+x_size, tile_joffset:tile_joffset+y_size))
    data_array = -9999.9

    if (present(time)) then
      corner(3) = time
      edge_lengths(3) = 1

      !allocate(time_array(x_size, y_size, 12)) ! assume monthly data; allocate 12 slots
      if (open_file(fileobj, nc_filename, "read", pelist=pes, is_restart=.false.)) then
        call read_data(fileobj, var_name, data_array, corner=corner, edge_lengths=edge_lengths)
        call close_file(fileobj)
      endif

      !data_array = time_array(:,:,time)
      !deallocate(time_array)
    else
      ! Following transition documents at https://github.com/NOAA-GFDL/FMS/tree/2021.03.01/fms2_io
      if (open_file(fileobj, nc_filename, "read", pelist=pes, is_restart=.false.)) then
        call read_data(fileobj, var_name, data_array, corner=corner(1:2), edge_lengths=edge_lengths(1:2))
        call close_file(fileobj)
      endif
    endif

  end subroutine alloc_read_tiled_data_r8_2d



end module fv_moving_nest_tiled_read_mod
