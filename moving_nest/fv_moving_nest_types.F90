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
!! @brief Provides data structures for moving nest functionality
!! @author W. Ramstrom, AOML/HRD   03/24/2022
!! @email William.Ramstrom@noaa.gov
! =======================================================================!

module fv_moving_nest_types_mod

#include <fms_platform.h>

#ifdef GFS_TYPES
  use GFS_typedefs,           only: kind_phys
#else
  use IPD_typedefs,           only: kind_phys => IPD_kind_phys
#endif

  use fms_mod,         only: check_nml_error
  use fv_arrays_mod,   only: fv_atmos_type
  use fv_mp_mod,       only: MAX_NNEST
  use mpp_mod,         only: input_nml_file, mpp_pe, read_input_nml
  use GFS_typedefs,    only: GFS_control_type

  implicit none

  type fv_moving_nest_flag_type
    ! Moving Nest Namelist Variables
    logical               :: is_moving_nest = .false.
    character(len=120)    :: surface_dir = "INPUT/moving_nest"
    integer               :: terrain_smoother = 4
    integer               :: vortex_tracker = 0
    integer               :: ntrack = 1
    integer               :: corral_x = 5
    integer               :: corral_y = 5

    integer               :: outatcf_lun = 600

    ! Moving nest related variables
    integer               :: move_cd_x = 0
    integer               :: move_cd_y = 0
    logical               :: do_move = .false.
  end type fv_moving_nest_flag_type

  ! Encapsulates the grid definition data, such as read from the netCDF files
  type grid_geometry
    integer   :: nx, ny, nxp, nyp

    real(kind=kind_phys), allocatable  :: lats(:,:)
    real(kind=kind_phys), allocatable  :: lons(:,:)

    !real, allocatable  :: dx(:,:)
    !real, allocatable  :: dy(:,:)
    real(kind=kind_phys), allocatable  :: area(:,:)
  end type grid_geometry

  type fv_moving_nest_prog_type
    real, _ALLOCATABLE                  :: delz(:,:,:)      _NULL   !< layer thickness (meters)
  end type fv_moving_nest_prog_type

  ! TODO deallocate these at end of model run.  They are only allocated once, at first nest move, inside mn_static_read_hires().
  !  Note these are only 32 bits for now; matching the precision of the input netCDF files
  !  though the model generally handles physics variables with 64 bit precision
  type mn_surface_grids
    real, allocatable  :: orog_grid(:,:)               _NULL  ! orography -- raw or filtered depending on namelist option, in meters
    real, allocatable  :: orog_std_grid(:,:)           _NULL  ! terrain standard deviation for gravity wave drag, in meters (?)
    real, allocatable  :: ls_mask_grid(:,:)            _NULL  ! land sea mask -- 0 for ocean/lakes, 1, for land.  Perhaps 2 for sea ice.
    real, allocatable  :: land_frac_grid(:,:)          _NULL  ! Continuous land fraction - 0.0 ocean, 0.5 half of each, 1.0 all land

    real, allocatable  :: parent_orog_grid(:,:)        _NULL  ! parent orography -- only used for terrain_smoother=1.
    !     raw or filtered depending on namelist option,in meters

    ! Soil variables
    real, allocatable  :: deep_soil_temp_grid(:,:)     _NULL  ! deep soil temperature at 5m, in degrees K
    real, allocatable  :: soil_type_grid(:,:)          _NULL  ! STATSGO soil type

    ! Vegetation variables
    real, allocatable  :: veg_frac_grid(:,:)           _NULL  ! vegetation fraction
    real, allocatable  :: veg_type_grid(:,:)           _NULL  ! IGBP vegetation type
    real, allocatable  :: veg_greenness_grid(:,:)      _NULL  ! NESDIS vegetation greenness; netCDF file has monthly values

    ! Orography variables
    real, allocatable  :: slope_type_grid(:,:)         _NULL  ! legacy 1 degree GFS slope type

    ! Albedo variables
    real, allocatable  :: max_snow_alb_grid(:,:)       _NULL  ! max snow albedo
    real, allocatable  :: facsf_grid(:,:)              _NULL  ! fractional coverage with strong cosz dependency
    real, allocatable  :: facwf_grid(:,:)              _NULL  ! fractional coverage with weak cosz dependency

    ! Snow free albedo
    !   strong cosz angle dependence = black sky
    !   weak cosz angle dependence = white sky
    !  From the chgres code in static_data.F90, we see the linkage of variable names:
    !   type(esmf_field), public           :: alvsf_target_grid !< visible black sky albedo
    !   type(esmf_field), public           :: alvwf_target_grid !< visible white sky albedo
    !   type(esmf_field), public           :: alnsf_target_grid !< near ir black sky albedo
    !   type(esmf_field), public           :: alnwf_target_grid !< near ir white sky albedo

    real, allocatable  :: alvsf_grid(:,:)              _NULL  ! Visible black sky albedo; netCDF file has monthly values
    real, allocatable  :: alvwf_grid(:,:)              _NULL  ! Visible white sky albedo; netCDF file has monthly values
    real, allocatable  :: alnsf_grid(:,:)              _NULL  ! Near IR black sky albedo; netCDF file has monthly values
    real, allocatable  :: alnwf_grid(:,:)              _NULL  ! Near IR white sky albedo; netCDF file has monthly values

  end type mn_surface_grids

  type fv_moving_nest_physics_type
    integer                             :: isd, ied, jsd, jed, npz
    logical                             :: move_physics, move_nsst
    integer                             :: lsoil, nmtvr, levs, ntot2d, ntot3d

    real, _ALLOCATABLE                  :: ts(:,:)          _NULL   !< 2D skin temperature/SST
    real, _ALLOCATABLE                  :: slmsk(:,:)       _NULL   !< land sea mask -- 0 for ocean/lakes, 1, for land.  Perhaps 2 for sea ice.
    real (kind=kind_phys), _ALLOCATABLE :: smc (:,:,:)      _NULL   !< soil moisture content
    real (kind=kind_phys), _ALLOCATABLE :: stc (:,:,:)      _NULL   !< soil temperature
    real (kind=kind_phys), _ALLOCATABLE :: slc (:,:,:)      _NULL   !< soil liquid water content

    real (kind=kind_phys), _ALLOCATABLE :: u10m (:,:)       _NULL   !< 10m u wind (a-grid?)
    real (kind=kind_phys), _ALLOCATABLE :: v10m (:,:)       _NULL   !< 10m v wind (a-grid?)
    real (kind=kind_phys), _ALLOCATABLE :: hprime (:,:,:)   _NULL   !< orographic metrics (maybe standard deviation?)

    real (kind=kind_phys), _ALLOCATABLE :: tprcp (:,:)      _NULL   !< total (of all precip types) precipitation rate

    real (kind=kind_phys), _ALLOCATABLE :: zorl (:,:)       _NULL   !< roughness length
    real (kind=kind_phys), _ALLOCATABLE :: zorll (:,:)      _NULL   !< land roughness length
    !real (kind=kind_phys), _ALLOCATABLE :: zorli (:,:)     _NULL   !< ice surface roughness length ! TODO do we need this?
    real (kind=kind_phys), _ALLOCATABLE :: zorlw (:,:)      _NULL   !< wave surface roughness length
    real (kind=kind_phys), _ALLOCATABLE :: zorlwav (:,:)    _NULL   !< wave surface roughness in cm derived from wave model

    real (kind=kind_phys), _ALLOCATABLE :: usfco (:,:)      _NULL   !< sea surface current (zonal)
    real (kind=kind_phys), _ALLOCATABLE :: vsfco (:,:)      _NULL   !< sea surface current (meridional)

    real (kind=kind_phys), _ALLOCATABLE :: sfalb_lnd(:,:)   _NULL   !< surface albedo over land for LSM
    real (kind=kind_phys), _ALLOCATABLE :: emis_lnd(:,:)    _NULL   !< surface emissivity over land for LSM
    real (kind=kind_phys), _ALLOCATABLE :: emis_ice(:,:)    _NULL   !< surface emissivity over ice for LSM
    real (kind=kind_phys), _ALLOCATABLE :: emis_wat(:,:)    _NULL   !< surface emissivity over water for LSM
    real (kind=kind_phys), _ALLOCATABLE :: sfalb_lnd_bck(:,:) _NULL !< snow-free albedo over land

    !real (kind=kind_phys), _ALLOCATABLE :: semis(:,:)       _NULL   !< surface lw emissivity in fraction
    !real (kind=kind_phys), _ALLOCATABLE :: semisbase(:,:)   _NULL   !< background surface emissivity
    !real (kind=kind_phys), _ALLOCATABLE :: sfalb(:,:)       _NULL   !< mean surface diffused sw albedo

    real (kind=kind_phys), _ALLOCATABLE :: alvsf(:,:)       _NULL   !< visible black sky albedo
    real (kind=kind_phys), _ALLOCATABLE :: alvwf(:,:)       _NULL   !< visible white sky albedo
    real (kind=kind_phys), _ALLOCATABLE :: alnsf(:,:)       _NULL   !< near IR black sky albedo
    real (kind=kind_phys), _ALLOCATABLE :: alnwf(:,:)       _NULL   !< near IR white sky albedo

    real (kind=kind_phys), _ALLOCATABLE :: albdirvis_lnd(:,:)       _NULL   !<
    real (kind=kind_phys), _ALLOCATABLE :: albdirnir_lnd(:,:)       _NULL   !<
    real (kind=kind_phys), _ALLOCATABLE :: albdifvis_lnd(:,:)       _NULL   !<
    real (kind=kind_phys), _ALLOCATABLE :: albdifnir_lnd(:,:)       _NULL   !<

    real (kind=kind_phys), _ALLOCATABLE :: facsf(:,:)       _NULL   !< fractional coverage for strong zenith angle albedo
    real (kind=kind_phys), _ALLOCATABLE :: facwf(:,:)       _NULL   !< fractional coverage for strong zenith angle albedo

    real (kind=kind_phys), _ALLOCATABLE :: lakefrac (:,:)   _NULL   !< lake  fraction [0:1]
    real (kind=kind_phys), _ALLOCATABLE :: lakedepth (:,:)  _NULL   !< lake  depth [ m ]

    real (kind=kind_phys), _ALLOCATABLE :: canopy (:,:)     _NULL   !< canopy water content
    real (kind=kind_phys), _ALLOCATABLE :: vegfrac (:,:)    _NULL   !< vegetation fraction
    real (kind=kind_phys), _ALLOCATABLE :: uustar (:,:)     _NULL   !< u* wind in similarity theory
    real (kind=kind_phys), _ALLOCATABLE :: shdmin (:,:)     _NULL   !< min fractional coverage of green vegetation
    real (kind=kind_phys), _ALLOCATABLE :: shdmax (:,:)     _NULL   !< max fractional coverage of green vegetation
    real (kind=kind_phys), _ALLOCATABLE :: tsfco (:,:)      _NULL   !< surface temperature ocean
    real (kind=kind_phys), _ALLOCATABLE :: tsfcl (:,:)      _NULL   !< surface temperature land
    real (kind=kind_phys), _ALLOCATABLE :: tsfc (:,:)       _NULL   !< surface temperature
    !real (kind=kind_phys), _ALLOCATABLE :: tsfc_radtime (:,:) _NULL !< surface temperature on radiative timestep

    real (kind=kind_phys), _ALLOCATABLE :: cv  (:,:)        _NULL   !< fraction of convective cloud
    real (kind=kind_phys), _ALLOCATABLE :: cvt (:,:)        _NULL   !< convective cloud top pressure
    real (kind=kind_phys), _ALLOCATABLE :: cvb (:,:)        _NULL   !< convective cloud bottom pressure

    real (kind=kind_phys), _ALLOCATABLE :: phy_f2d (:,:,:)  _NULL   !< 2D physics variables
    real (kind=kind_phys), _ALLOCATABLE :: phy_f3d(:,:,:,:) _NULL   !< 3D physics variables

    ! NSST Variables

    real (kind=kind_phys), _ALLOCATABLE :: tref (:,:)       _NULL   !< reference temperature for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: z_c (:,:)        _NULL   !< coefficient for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: c_0 (:,:)        _NULL   !< coefficient for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: c_d (:,:)        _NULL   !< coefficient for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: w_0 (:,:)        _NULL   !< coefficient for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: w_d (:,:)        _NULL   !< coefficient for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xt (:,:)         _NULL   !< heat content for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xs (:,:)         _NULL   !< salinity for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xu (:,:)         _NULL   !< u current constant for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xv (:,:)         _NULL   !< v current constant for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xz (:,:)         _NULL   !< DTL thickness for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: zm (:,:)         _NULL   !< MXL for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xtts (:,:)       _NULL   !< d(xt)/d(ts) for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xzts (:,:)       _NULL   !< d(xz)/d(ts) for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: d_conv (:,:)     _NULL   !< think of free convection layer for NSSTM
    ! real (kind=kind_phys), _ALLOCATABLE :: ifd (:,:)      _NULL   !< index to start DTM run  for NSSTM   ! TODO Probably can't interpolate an index.
    !  IFD values are 0 for land, and 1 for oceans/lakes -- reverse of the land sea mask
    !  Land Sea Mask has values of 0 for oceans/lakes, 1 for land, 2 for sea ice
    real (kind=kind_phys), _ALLOCATABLE :: dt_cool (:,:)    _NULL   !< sub-layer cooling amount for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: qrain (:,:)      _NULL   !< sensible heat flux due to rainfall for NSSTM

  contains

    procedure, public :: alloc_dealloc => fv_moving_nest_physics_alloc_dealloc

  end type fv_moving_nest_physics_type

  type fv_moving_nest_type
    type(fv_moving_nest_flag_type)    :: mn_flag   ! Mostly namelist variables
    type(mn_surface_grids)            :: mn_static
    type(fv_moving_nest_prog_type)    :: mn_prog
    type(fv_moving_nest_physics_type) :: mn_phys

    type(grid_geometry)               :: parent_geo
    type(grid_geometry)               :: fp_super_tile_geo
  end type fv_moving_nest_type

  ! Moving Nest Namelist Variables
  logical, dimension(MAX_NNEST) :: is_moving_nest = .False.
  character(len=120)            :: surface_dir = "INPUT/moving_nest"
  integer, dimension(MAX_NNEST) :: terrain_smoother = 4  ! 0 -- all high-resolution data, 1 - static nest smoothing algorithm with blending zone of 5 points, 2 - blending zone of 10 points, 5 - 5 point smoother, 9 - 9 point smoother
  integer, dimension(MAX_NNEST) :: vortex_tracker = 0 ! 0 - not a moving nest, tracker not needed
  ! 1 - prescribed nest moving
  ! 2 - following child domain center
  ! 3 - tracking Min MSLP
  ! 6 - simplified version of GFDL tracker, adopted from HWRF's internal vortex tracker.
  ! 7 - nearly the full storm tracking algorithm from GFDL vortex tracker. The only part that is missing is the part that gives up when the storm dissipates, which is left out intentionally. Adopted from HWRF's internal vortex tracker.
  integer, dimension(MAX_NNEST) :: ntrack = 1 ! number of dt_atmos steps to call the vortex tracker, tracker time step = ntrack*dt_atmos
  integer, dimension(MAX_NNEST) :: move_cd_x = 0 ! the number of parent domain grid cells to move in i direction
  integer, dimension(MAX_NNEST) :: move_cd_y = 0 ! the number of parent domain grid cells to move in j direction
  ! used to control prescribed nest moving, when vortex_tracker=1
  ! the move happens every ntrack*dt_atmos seconds
  ! positive is to move in increasing i and j direction, and
  ! negative is to move in decreasing i and j direction.
  ! 0 means no move. The limitation is to move only 1 grid cell at each move.
  integer, dimension(MAX_NNEST) :: corral_x = 5 ! Minimum parent gridpoints on each side of nest in i direction
  integer, dimension(MAX_NNEST) :: corral_y = 5 ! Minimum parent gridpoints on each side of nest in j direction

  integer, dimension(MAX_NNEST) :: outatcf_lun = 600  ! base fortran unit number to write out the partial atcfunix file from the internal tracker

  type(fv_moving_nest_type), _ALLOCATABLE, target    :: Moving_nest(:)

contains

  subroutine fv_moving_nest_init(Atm, this_grid)
    type(fv_atmos_type), allocatable, intent(in) :: Atm(:)
    integer, intent(in)                          :: this_grid

    integer :: n, ngrids

    ! Allocate the array of fv_moving_nest_type structures of the proper length
    allocate(Moving_nest(size(Atm)))

    ! Configure namelist variables

    ngrids = size(Atm)

    call read_input_nml(Atm(1)%nml_filename) !re-reads top level file into internal namelist

    ! Read in namelist

    call read_namelist_moving_nest_nml

    do n=1,ngrids
      if (Atm(n)%neststruct%nested) then
        Moving_nest(n)%mn_flag%is_moving_nest         = is_moving_nest(n)
        Moving_nest(n)%mn_flag%surface_dir            = trim(surface_dir)
        Moving_nest(n)%mn_flag%terrain_smoother       = terrain_smoother(n)
        Moving_nest(n)%mn_flag%vortex_tracker         = vortex_tracker(n)
        Moving_nest(n)%mn_flag%ntrack                 = ntrack(n)
        Moving_nest(n)%mn_flag%move_cd_x              = move_cd_x(n)
        Moving_nest(n)%mn_flag%move_cd_y              = move_cd_y(n)
        Moving_nest(n)%mn_flag%corral_x               = corral_x(n)
        Moving_nest(n)%mn_flag%corral_y               = corral_y(n)
        Moving_nest(n)%mn_flag%outatcf_lun            = outatcf_lun(n)
      else
        Moving_nest(n)%mn_flag%is_moving_nest         = .false.
        Moving_nest(n)%mn_flag%vortex_tracker         = 0
        Moving_nest(n)%mn_flag%ntrack                 = 1
        Moving_nest(n)%mn_flag%move_cd_x              = 0
        Moving_nest(n)%mn_flag%move_cd_y              = 0
        Moving_nest(n)%mn_flag%corral_x               = 5
        Moving_nest(n)%mn_flag%corral_y               = 5
        Moving_nest(n)%mn_flag%outatcf_lun            = 600
      endif
    enddo


    call read_input_nml(Atm(this_grid)%nml_filename) !re-reads into internal namelist


  end subroutine fv_moving_nest_init

  subroutine read_namelist_moving_nest_nml
    integer :: f_unit, ios, ierr
    namelist /fv_moving_nest_nml/ surface_dir, is_moving_nest, terrain_smoother, &
        vortex_tracker, ntrack, move_cd_x, move_cd_y, corral_x, corral_y, outatcf_lun

#ifdef INTERNAL_FILE_NML
    read (input_nml_file,fv_moving_nest_nml,iostat=ios)
    ierr = check_nml_error(ios,'fv_moving_nest_nml')
#else
    f_unit=open_namelist_file()
    rewind (f_unit)
    read (f_unit,fv_moving_nest_nml,iostat=ios)
    ierr = check_nml_error(ios,'fv_moving_nest_nml')
    call close_file(f_unit)
#endif

  end subroutine read_namelist_moving_nest_nml

  subroutine deallocate_fv_moving_nests(n)
    integer, intent(in)   :: n

    integer :: i

    do i=1,n
      call deallocate_fv_moving_nest(i)
    enddo
    deallocate(Moving_nest)
  end subroutine deallocate_fv_moving_nests

  subroutine deallocate_fv_moving_nest(n)
    integer, intent(in)   :: n

    call deallocate_fv_moving_nest_prog_type(Moving_nest(n)%mn_prog)
    call deallocate_fv_moving_nest_physics_type(Moving_nest(n)%mn_phys)

  end subroutine deallocate_fv_moving_nest


  subroutine  allocate_fv_moving_nest_prog_type(isd, ied, jsd, jed, npz, mn_prog)
    integer, intent(in)                           :: isd, ied, jsd, jed, npz
    type(fv_moving_nest_prog_type), intent(inout) :: mn_prog

    allocate ( mn_prog%delz(isd:ied, jsd:jed, 1:npz) )
    mn_prog%delz = +99999.9

  end subroutine allocate_fv_moving_nest_prog_type

  subroutine  deallocate_fv_moving_nest_prog_type(mn_prog)
    type(fv_moving_nest_prog_type), intent(inout) :: mn_prog

    if (allocated(mn_prog%delz)) deallocate(mn_prog%delz)

  end subroutine deallocate_fv_moving_nest_prog_type

  subroutine  deallocate_fv_moving_nest_physics_type(mn_phys)
    implicit none
    type(fv_moving_nest_physics_type), intent(inout) :: mn_phys

    call mn_phys%alloc_dealloc(.false.)
  end subroutine deallocate_fv_moving_nest_physics_type

  subroutine  allocate_fv_moving_nest_physics_type(isd, ied, jsd, jed, npz, GFS_Control, move_physics, move_nsst, mn_phys)
    implicit none
    type(GFS_control_type), intent(in)               :: GFS_Control   !< Physics metadata
    integer, intent(in)                              :: npz           !< vertical dimension
    integer, intent(in)                              :: isd, ied, jsd, jed !< horizontal dimensions with halos added
    logical, intent(in)                              :: move_physics  !< Flag for moving physics variables. Always true.
    logical, intent(in)                              :: move_nsst     !< Flag for moving nsst variables
    type(fv_moving_nest_physics_type), target        :: mn_phys       !< Storage for physics fields with halos

    mn_phys%isd = isd
    mn_phys%ied = ied
    mn_phys%jsd = jsd
    mn_phys%jed = jed
    mn_phys%npz = npz
    mn_phys%move_physics = move_physics
    mn_phys%move_nsst = move_nsst
    mn_phys%lsoil = GFS_Control%lsoil
    mn_phys%nmtvr = GFS_Control%nmtvr
    mn_phys%levs = GFS_Control%levs
    mn_phys%ntot2d = GFS_Control%ntot2d
    mn_phys%ntot3d = GFS_Control%ntot3d

    call mn_phys%alloc_dealloc(.true.)
  end subroutine allocate_fv_moving_nest_physics_type

  subroutine fv_moving_nest_physics_alloc_dealloc(mn_phys, to_alloc)
    implicit none
    class(fv_moving_nest_physics_type) :: mn_phys
    logical, intent(in) :: to_alloc

    integer :: isd, ied, jsd, jed, npz
    logical :: move_physics, move_nsst
    integer :: lsoil, nmtvr, levs, ntot2d, ntot3d

    interface alloc_dealloc
       procedure alloc_dealloc_2d_r4, alloc_dealloc_2d_r8, alloc_dealloc_3d, alloc_dealloc_4d
    end interface alloc_dealloc

    ! Copy these to locals to reduce typing
    isd = mn_phys%isd
    ied = mn_phys%ied
    jsd = mn_phys%jsd
    jed = mn_phys%jed
    npz = mn_phys%npz
    move_physics = mn_phys%move_physics
    move_nsst = mn_phys%move_nsst
    lsoil = mn_phys%lsoil
    nmtvr = mn_phys%nmtvr
    levs = mn_phys%levs
    ntot2d = mn_phys%ntot2d
    ntot3d = mn_phys%ntot3d

    ! The local/temporary variables need to be allocated to the larger data (compute + halos) domain so that the nest motion code has halos to use
    call alloc_dealloc(mn_phys%ts)

    if (move_physics) then
       call alloc_dealloc(mn_phys%slmsk)
       call alloc_dealloc(mn_phys%smc, lsoil)
       call alloc_dealloc(mn_phys%stc, lsoil)
       call alloc_dealloc(mn_phys%slc, lsoil)

       call alloc_dealloc(mn_phys%sfalb_lnd)
       call alloc_dealloc(mn_phys%emis_lnd)
       call alloc_dealloc(mn_phys%emis_ice)
       call alloc_dealloc(mn_phys%emis_wat)
       call alloc_dealloc(mn_phys%sfalb_lnd_bck)

       !call alloc_dealloc(mn_phys%semis)
       !call alloc_dealloc(mn_phys%semisbase)
       !call alloc_dealloc(mn_phys%sfalb)

       call alloc_dealloc(mn_phys%u10m)
       call alloc_dealloc(mn_phys%v10m)
       call alloc_dealloc(mn_phys%tprcp)

       call alloc_dealloc(mn_phys%hprime, nmtvr)

       call alloc_dealloc(mn_phys%zorl)
       call alloc_dealloc(mn_phys%zorll)
       call alloc_dealloc(mn_phys%zorlwav)
       call alloc_dealloc(mn_phys%zorlw)

       call alloc_dealloc(mn_phys%usfco)
       call alloc_dealloc(mn_phys%vsfco)

       call alloc_dealloc(mn_phys%alvsf)
       call alloc_dealloc(mn_phys%alvwf)
       call alloc_dealloc(mn_phys%alnsf)
       call alloc_dealloc(mn_phys%alnwf)

       call alloc_dealloc(mn_phys%facsf)
       call alloc_dealloc(mn_phys%facwf)

       call alloc_dealloc(mn_phys%lakefrac)
       call alloc_dealloc(mn_phys%lakedepth)

       call alloc_dealloc(mn_phys%canopy)
       call alloc_dealloc(mn_phys%vegfrac)
       call alloc_dealloc(mn_phys%uustar)
       call alloc_dealloc(mn_phys%shdmin)
       call alloc_dealloc(mn_phys%shdmax)
       call alloc_dealloc(mn_phys%tsfco)
       call alloc_dealloc(mn_phys%tsfcl)
       call alloc_dealloc(mn_phys%tsfc)
       !call alloc_dealloc(mn_phys%tsfc_radtime)

       call alloc_dealloc(mn_phys%albdirvis_lnd)
       call alloc_dealloc(mn_phys%albdirnir_lnd)
       call alloc_dealloc(mn_phys%albdifvis_lnd)
       call alloc_dealloc(mn_phys%albdifnir_lnd)

       call alloc_dealloc(mn_phys%cv)
       call alloc_dealloc(mn_phys%cvt)
       call alloc_dealloc(mn_phys%cvb)

       call alloc_dealloc(mn_phys%phy_f2d, ntot2d)
       call alloc_dealloc(mn_phys%phy_f3d, levs, ntot3d)
    end if

    if (move_nsst) then
       call alloc_dealloc(mn_phys%tref)
       call alloc_dealloc(mn_phys%z_c)
       call alloc_dealloc(mn_phys%c_0)
       call alloc_dealloc(mn_phys%c_d)
       call alloc_dealloc(mn_phys%w_0)
       call alloc_dealloc(mn_phys%w_d)
       call alloc_dealloc(mn_phys%xt)
       call alloc_dealloc(mn_phys%xs)
       call alloc_dealloc(mn_phys%xu)
       call alloc_dealloc(mn_phys%xv)
       call alloc_dealloc(mn_phys%xz)
       call alloc_dealloc(mn_phys%zm)
       call alloc_dealloc(mn_phys%xtts)
       call alloc_dealloc(mn_phys%xzts)
       call alloc_dealloc(mn_phys%d_conv)
       !call alloc_dealloc(mn_phys%ifd)
       call alloc_dealloc(mn_phys%dt_cool)
       call alloc_dealloc(mn_phys%qrain)
    end if

  contains

    subroutine alloc_dealloc_2d_r8(var)
      implicit none
      real(kind=8), allocatable :: var(:,:)

      if(allocated(var)) then
         deallocate(var)
      endif
      if(to_alloc) then
         allocate(var(isd:ied, jsd:jed))
         var = +99999.9
      endif
    end subroutine alloc_dealloc_2d_r8

    subroutine alloc_dealloc_2d_r4(var)
      implicit none
      real(kind=4), allocatable :: var(:,:)

      if(allocated(var)) then
         deallocate(var)
      endif
      if(to_alloc) then
         allocate(var(isd:ied, jsd:jed))
         var = +99999.9
      endif
    end subroutine alloc_dealloc_2d_r4

    subroutine alloc_dealloc_3d(var, kdim)
      implicit none
      real(kind_phys), allocatable :: var(:,:,:)
      integer, intent(in) :: kdim

      if(allocated(var)) then
         deallocate(var)
      endif
      if(to_alloc) then
         allocate(var(isd:ied, jsd:jed, kdim))
         var = +99999.9
      endif
    end subroutine alloc_dealloc_3d

    subroutine alloc_dealloc_4d(var, kdim, mdim)
      implicit none
      real(kind_phys), allocatable :: var(:,:,:,:)
      integer, intent(in) :: kdim, mdim

      if(allocated(var)) then
         deallocate(var)
      endif
      if(to_alloc) then
         allocate(var(isd:ied, jsd:jed, kdim, mdim))
         var = +99999.9
      endif
    end subroutine alloc_dealloc_4d

  end subroutine fv_moving_nest_physics_alloc_dealloc


end module fv_moving_nest_types_mod
