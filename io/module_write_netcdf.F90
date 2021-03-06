#define ESMF_ERR_RETURN(rc) if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#define NC_ERR_STOP(status) \
    if (status /= nf90_noerr) write(0,*) "line ", __LINE__, trim(nf90_strerror(status)); \
    if (status /= nf90_noerr) call ESMF_Finalize(endflag=ESMF_END_ABORT)

module module_write_netcdf

  use esmf
  use netcdf
  use module_fv3_io_def,only : ideflate, nbits, &
                               output_grid,dx,dy,lon1,lat1,lon2,lat2

  implicit none
  private
  public write_netcdf
  public write_grid_netcdf

  contains

!----------------------------------------------------------------------------------------
  subroutine write_netcdf(fieldbundle, wrtfb, filename, mpi_comm, mype, im, jm, ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d, rc)
!
    type(ESMF_FieldBundle), intent(in) :: fieldbundle
    type(ESMF_FieldBundle), intent(in) :: wrtfb
    character(*), intent(in)           :: filename
    integer, intent(in)                :: mpi_comm
    integer, intent(in)                :: mype
    integer, intent(in)                :: im, jm
    integer, intent(in)                :: ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d
    integer, optional,intent(out)      :: rc
!
!** local vars
    integer :: i,j,m,n,k
    integer :: lm

    integer, dimension(:), allocatable     :: fldlev
    real(4), dimension(:,:), allocatable   :: arrayr4
    real(8), dimension(:,:), allocatable   :: arrayr8
    real(4), dimension(:,:,:), allocatable :: arrayr4_3d,arrayr4_3d_save
    real(8), dimension(:,:,:), allocatable :: arrayr8_3d

    real(8) x(im),y(jm)
    integer :: fieldCount, fieldDimCount, gridDimCount
    integer, dimension(:), allocatable   :: ungriddedLBound, ungriddedUBound

    type(ESMF_Field), allocatable        :: fcstField(:)
    type(ESMF_TypeKind_Flag)             :: typekind
    type(ESMF_TypeKind_Flag)             :: attTypeKind
    type(ESMF_Grid)                      :: wrtgrid
    type(ESMF_Array)                     :: array

    integer :: attcount
    character(len=ESMF_MAXSTR) :: attName, fldName

    integer :: varival
    real(4) :: varr4val, scale_fact, offset, dataMin, dataMax
    real(4), allocatable, dimension(:) :: compress_err
    real(8) :: varr8val
    character(len=ESMF_MAXSTR) :: varcval

    character(128) :: time_units

    integer :: ncerr
    integer :: ncid
    integer :: oldMode
    integer :: im_dimid, jm_dimid, pfull_dimid, phalf_dimid, time_dimid
    integer :: im_varid, jm_varid, lm_varid, time_varid, lon_varid, lat_varid
    integer, dimension(:), allocatable :: varids
    logical shuffle

    call ESMF_FieldBundleGet(fieldbundle, fieldCount=fieldCount, rc=rc); ESMF_ERR_RETURN(rc)

    allocate(compress_err(fieldCount)); compress_err=-999.
    allocate(fldlev(fieldCount)) ; fldlev = 0
    allocate(fcstField(fieldCount))
    allocate(varids(fieldCount))

    call ESMF_FieldBundleGet(fieldbundle, fieldList=fcstField, grid=wrtGrid, &
!                             itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                             rc=rc); ESMF_ERR_RETURN(rc)

    call ESMF_GridGet(wrtgrid, dimCount=gridDimCount, rc=rc); ESMF_ERR_RETURN(rc)

    do i=1,fieldCount
       call ESMF_FieldGet(fcstField(i), dimCount=fieldDimCount, rc=rc); ESMF_ERR_RETURN(rc)
       if (fieldDimCount > 3) then
          write(0,*)"write_netcdf: Only 2D and 3D fields are supported!"
          stop
       end if
       if (fieldDimCount > gridDimCount) then
         allocate(ungriddedLBound(fieldDimCount-gridDimCount))
         allocate(ungriddedUBound(fieldDimCount-gridDimCount))
         call ESMF_FieldGet(fcstField(i), &
                            ungriddedLBound=ungriddedLBound, &
                            ungriddedUBound=ungriddedUBound, rc=rc); ESMF_ERR_RETURN(rc)
         fldlev(i) = ungriddedUBound(fieldDimCount-gridDimCount) - &
                     ungriddedLBound(fieldDimCount-gridDimCount) + 1
         deallocate(ungriddedLBound)
         deallocate(ungriddedUBound)
       else if (fieldDimCount == 2) then
         fldlev(i) = 1
       end if
    end do

    lm = maxval(fldlev(:))

    allocate(arrayr4(im,jm))
    allocate(arrayr8(im,jm))
    allocate(arrayr4_3d(im,jm,lm),arrayr4_3d_save(im,jm,lm))
    allocate(arrayr8_3d(im,jm,lm))

! create netcdf file and enter define mode
    if (mype==0) then

    ncerr = nf90_create(trim(filename),&
            cmode=IOR(IOR(NF90_CLOBBER,NF90_NETCDF4),NF90_CLASSIC_MODEL),&
            ncid=ncid); NC_ERR_STOP(ncerr)
    ncerr = nf90_set_fill(ncid, NF90_NOFILL, oldMode); NC_ERR_STOP(ncerr)

    ! define dimensions
    ncerr = nf90_def_dim(ncid, "grid_xt", im, im_dimid); NC_ERR_STOP(ncerr)
    ncerr = nf90_def_dim(ncid, "grid_yt", jm, jm_dimid); NC_ERR_STOP(ncerr)
    ! define coordinate variables
    ncerr = nf90_def_var(ncid, "grid_xt", NF90_DOUBLE, im_dimid, im_varid); NC_ERR_STOP(ncerr)
    ncerr = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/im_dimid,jm_dimid/), lon_varid); NC_ERR_STOP(ncerr)
    ncerr = nf90_put_att(ncid, lon_varid, "long_name", "T-cell longitude"); NC_ERR_STOP(ncerr)
    ncerr = nf90_put_att(ncid, lon_varid, "units", "degrees_E"); NC_ERR_STOP(ncerr)
    ncerr = nf90_put_att(ncid, im_varid, "cartesian_axis", "X"); NC_ERR_STOP(ncerr)
    ncerr = nf90_def_var(ncid, "grid_yt", NF90_DOUBLE, jm_dimid, jm_varid); NC_ERR_STOP(ncerr)
    ncerr = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/im_dimid,jm_dimid/), lat_varid); NC_ERR_STOP(ncerr)
    ncerr = nf90_put_att(ncid, lat_varid, "long_name", "T-cell latitude"); NC_ERR_STOP(ncerr)
    ncerr = nf90_put_att(ncid, lat_varid, "units", "degrees_N"); NC_ERR_STOP(ncerr)
    ncerr = nf90_put_att(ncid, jm_varid, "cartesian_axis", "Y"); NC_ERR_STOP(ncerr)

    if (lm > 1) then
      call add_dim(ncid, "pfull", pfull_dimid, wrtgrid, rc)
      call add_dim(ncid, "phalf", phalf_dimid, wrtgrid, rc)
    end if

    call add_dim(ncid, "time", time_dimid, wrtgrid, rc)

    call get_global_attr(wrtfb, ncid, rc)

    do i=1, fieldCount
      call ESMF_FieldGet(fcstField(i), name=fldName, typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)

      ! define variables
      if (fldlev(i) == 1) then
        if (typekind == ESMF_TYPEKIND_R4) then
          if (ideflate > 0) then
             if (ichunk2d < 0 .or. jchunk2d < 0) then
                ! let netcdf lib choose chunksize
                ! shuffle filter on for 2d fields (lossless compression)
                ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                        (/im_dimid,jm_dimid,time_dimid/), varids(i), &
                        shuffle=.true.,deflate_level=ideflate); NC_ERR_STOP(ncerr)
             else
                ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                        (/im_dimid,jm_dimid,time_dimid/), varids(i), &
                        shuffle=.true.,deflate_level=ideflate,&
                        chunksizes=(/ichunk2d,jchunk2d,1/),cache_size=40*im*jm); NC_ERR_STOP(ncerr)
             endif
          else
             ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
             (/im_dimid,jm_dimid,time_dimid/), varids(i)); NC_ERR_STOP(ncerr)
          endif
        else if (typekind == ESMF_TYPEKIND_R8) then
           ncerr = nf90_def_var(ncid, trim(fldName), NF90_DOUBLE, &
                               (/im_dimid,jm_dimid,time_dimid/), varids(i)); NC_ERR_STOP(ncerr)
        else
           write(0,*)'Unsupported typekind ', typekind
           stop
        end if
      else if (fldlev(i) > 1) then
         if (typekind == ESMF_TYPEKIND_R4) then
           if (ideflate > 0) then
             ! shuffle filter off for 3d fields using lossy compression
             if (nbits > 0) then
                shuffle=.false.
             else
                shuffle=.true.
             endif
             if (ichunk3d < 0 .or. jchunk3d < 0 .or. kchunk3d < 0) then
                ! let netcdf lib choose chunksize
                ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                        (/im_dimid,jm_dimid,pfull_dimid,time_dimid/), varids(i), &
                        shuffle=shuffle,deflate_level=ideflate); NC_ERR_STOP(ncerr)
             else
                ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                       (/im_dimid,jm_dimid,pfull_dimid,time_dimid/), varids(i), &
                       shuffle=shuffle,deflate_level=ideflate,&
                       chunksizes=(/ichunk3d,jchunk3d,kchunk3d,1/)); NC_ERR_STOP(ncerr)
             endif
           else
             ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                     (/im_dimid,jm_dimid,pfull_dimid,time_dimid/), varids(i)); NC_ERR_STOP(ncerr)
           endif 
         else if (typekind == ESMF_TYPEKIND_R8) then
           ncerr = nf90_def_var(ncid, trim(fldName), NF90_DOUBLE, &
                                (/im_dimid,jm_dimid,pfull_dimid,time_dimid/), varids(i)); NC_ERR_STOP(ncerr)
        else
          write(0,*)'Unsupported typekind ', typekind
          stop
        end if
      end if

      ! define variable attributes
      call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, Count=attcount, &
                             rc=rc); ESMF_ERR_RETURN(rc)

      do j=1,attCount
        call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                               attnestflag=ESMF_ATTNEST_OFF, attributeIndex=j, &
                               name=attName, typekind=attTypeKind, itemCount=n, &
                               rc=rc); ESMF_ERR_RETURN(rc)

        if ( index(trim(attName),"ESMF") /= 0 ) then
           cycle
        endif

        if (attTypeKind==ESMF_TYPEKIND_I4) then
           call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                  name=trim(attName), value=varival, &
                                  rc=rc); ESMF_ERR_RETURN(rc)
           ncerr = nf90_put_att(ncid, varids(i), trim(attName), varival); NC_ERR_STOP(ncerr)

        else if (attTypeKind==ESMF_TYPEKIND_R4) then
           call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                  name=trim(attName), value=varr4val, &
                                  rc=rc); ESMF_ERR_RETURN(rc)
           ncerr = nf90_put_att(ncid, varids(i), trim(attName), varr4val); NC_ERR_STOP(ncerr)

        else if (attTypeKind==ESMF_TYPEKIND_R8) then
           call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                  name=trim(attName), value=varr8val, &
                                  rc=rc); ESMF_ERR_RETURN(rc)
           if (trim(attName) /= '_FillValue') then
              ! FIXME:  _FillValue must be cast to var type for recent versions of netcdf
              ncerr = nf90_put_att(ncid, varids(i), trim(attName), varr8val); NC_ERR_STOP(ncerr)
           endif

        else if (attTypeKind==ESMF_TYPEKIND_CHARACTER) then
           call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                  name=trim(attName), value=varcval, &
                                  rc=rc); ESMF_ERR_RETURN(rc)
           ncerr = nf90_put_att(ncid, varids(i), trim(attName), trim(varcval)); NC_ERR_STOP(ncerr)

        end if

      end do ! j=1,attCount

    end do   ! i=1,fieldCount

    ! write grid_xt, grid_yt attributes
    if (trim(output_grid) == 'gaussian_grid' .or. &
        trim(output_grid) == 'global_latlon' .or. &
        trim(output_grid) == 'regional_latlon') then
       ncerr = nf90_put_att(ncid, im_varid, "long_name", "T-cell longitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, im_varid, "units", "degrees_E"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "long_name", "T-cell latiitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "units", "degrees_N"); NC_ERR_STOP(ncerr)
    else if (trim(output_grid) == 'rotated_latlon') then
       ncerr = nf90_put_att(ncid, im_varid, "long_name", "rotated T-cell longiitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, im_varid, "units", "degrees"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "long_name", "rotated T-cell latiitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "units", "degrees"); NC_ERR_STOP(ncerr)
    else if (trim(output_grid) == 'lambert_conformal') then
       ncerr = nf90_put_att(ncid, im_varid, "long_name", "x-coordinate of projection"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, im_varid, "units", "meters"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "long_name", "y-coordinate of projection"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "units", "meters"); NC_ERR_STOP(ncerr)
    endif

    ncerr = nf90_enddef(ncid); NC_ERR_STOP(ncerr)
    end if

! end of define mode

    ! write grid_xt, grid_yt values
    call ESMF_GridGetCoord(wrtGrid, coordDim=1, array=array, rc=rc); ESMF_ERR_RETURN(rc)
    call ESMF_ArrayGather(array, arrayr8, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
    if (mype==0) then
       if (trim(output_grid) == 'gaussian_grid' .or. &
           trim(output_grid) == 'global_latlon' .or. &
           trim(output_grid) == 'regional_latlon') then
          ncerr = nf90_put_var(ncid, im_varid, values=arrayr8(:,1)  ); NC_ERR_STOP(ncerr)
       else if (trim(output_grid) == 'rotated_latlon') then
          do i=1,im
             x(i) = lon1 + (lon2-lon1)/(im-1) * (i-1)
          enddo
          ncerr = nf90_put_var(ncid, im_varid, values=x  ); NC_ERR_STOP(ncerr)
       else if (trim(output_grid) == 'lambert_conformal') then
          do i=1,im
             x(i) = dx * (i-1)
          enddo
          ncerr = nf90_put_var(ncid, im_varid, values=x  ); NC_ERR_STOP(ncerr)
       endif
       ncerr = nf90_put_var(ncid, lon_varid, values=arrayr8 ); NC_ERR_STOP(ncerr)
    endif

    call ESMF_GridGetCoord(wrtGrid, coordDim=2, array=array, rc=rc); ESMF_ERR_RETURN(rc)
    call ESMF_ArrayGather(array, arrayr8, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
    if (mype==0) then
       if (trim(output_grid) == 'gaussian_grid' .or. &
           trim(output_grid) == 'global_latlon' .or. &
           trim(output_grid) == 'regional_latlon') then
          ncerr = nf90_put_var(ncid, jm_varid, values=arrayr8(1,:)  ); NC_ERR_STOP(ncerr)
       else if (trim(output_grid) == 'rotated_latlon') then
          do j=1,jm
             y(j) = lat1 + (lat2-lat1)/(jm-1) * (j-1)
          enddo
          ncerr = nf90_put_var(ncid, jm_varid, values=y  ); NC_ERR_STOP(ncerr)
       else if (trim(output_grid) == 'lambert_conformal') then
          do j=1,jm
             y(j) = dy * (j-1)
          enddo
          ncerr = nf90_put_var(ncid, jm_varid, values=y  ); NC_ERR_STOP(ncerr)
       endif
       ncerr = nf90_put_var(ncid, lat_varid, values=arrayr8 ); NC_ERR_STOP(ncerr)
    endif

    do i=1, fieldCount

       call ESMF_FieldGet(fcstField(i),name=fldName,typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)

       if (fldlev(i) == 1) then
         if (typekind == ESMF_TYPEKIND_R4) then
           call ESMF_FieldGather(fcstField(i), arrayr4, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
           if (mype==0) then
             ncerr = nf90_put_var(ncid, varids(i), values=arrayr4, start=(/1,1,1/),count=(/im,jm,1/) ); NC_ERR_STOP(ncerr)
           end if
         else if (typekind == ESMF_TYPEKIND_R8) then
           call ESMF_FieldGather(fcstField(i), arrayr8, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
           if (mype==0) then
             ncerr = nf90_put_var(ncid, varids(i), values=arrayr8, start=(/1,1,1/),count=(/im,jm,1/) ); NC_ERR_STOP(ncerr)
           end if
         end if
      else if (fldlev(i) > 1) then
         if (typekind == ESMF_TYPEKIND_R4) then
           call ESMF_FieldGather(fcstField(i), arrayr4_3d, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
           if (mype==0) then
             if (ideflate > 0 .and. nbits > 0) then
                ! Lossy compression if nbits>0.
                ! The floating point data is quantized to improve compression
                ! See doi:10.5194/gmd-10-413-2017.  The method employed
                ! here is identical to the 'scaled linear packing' method in
                ! that paper, except that the data are scaling into an arbitrary
                ! range (2**nbits-1 not just 2**16-1) and are stored as
                ! re-scaled floats instead of short integers.
                ! The zlib algorithm does almost as
                ! well packing the re-scaled floats as it does the scaled
                ! integers, and this avoids the need for the client to apply the
                ! rescaling (plus it allows the ability to adjust the packing
                ! range).
                arrayr4_3d_save = arrayr4_3d
                dataMax = maxval(arrayr4_3d); dataMin = minval(arrayr4_3d)
                arrayr4_3d = quantized(arrayr4_3d_save, nbits, dataMin, dataMax)
                ! compute max abs compression error, save as a variable
                ! attribute.
                compress_err(i) = maxval(abs(arrayr4_3d_save-arrayr4_3d))
             endif
             ncerr = nf90_put_var(ncid, varids(i), values=arrayr4_3d, start=(/1,1,1/),count=(/im,jm,lm,1/) ); NC_ERR_STOP(ncerr)
           end if
         else if (typekind == ESMF_TYPEKIND_R8) then
           call ESMF_FieldGather(fcstField(i), arrayr8_3d, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
           if (mype==0) then
             ncerr = nf90_put_var(ncid, varids(i), values=arrayr8_3d, start=(/1,1,1/),count=(/im,jm,lm,1/) ); NC_ERR_STOP(ncerr)
           end if
         end if

      end if

    end do

    if (ideflate > 0 .and. nbits > 0 .and. mype == 0) then
       ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
       do i=1, fieldCount
          if (compress_err(i) > 0) then
             ncerr = nf90_put_att(ncid, varids(i), 'max_abs_compression_error', compress_err(i)); NC_ERR_STOP(ncerr)
             ncerr = nf90_put_att(ncid, varids(i), 'nbits', nbits); NC_ERR_STOP(ncerr)
          endif
       enddo
       ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
    endif

    deallocate(arrayr4)
    deallocate(arrayr8)
    deallocate(arrayr4_3d,arrayr4_3d_save)
    deallocate(arrayr8_3d)

    deallocate(fcstField)
    deallocate(varids)
    deallocate(compress_err)

    if (mype==0) then
    ncerr = nf90_close(ncid=ncid); NC_ERR_STOP(ncerr)
    end if

  end subroutine write_netcdf

!----------------------------------------------------------------------------------------
  subroutine write_grid_netcdf(grid, fileName, overwrite, status, &
    timeslice, iofmt, relaxedflag, regridArea, rc)
!
    type(ESMF_Grid),            intent(in)            :: grid
    character(len=*),           intent(in),  optional :: fileName
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    type(ESMF_IOFmt_Flag),      intent(in),  optional :: iofmt
    logical,                    intent(in),  optional :: relaxedflag
    logical,                    intent(in),  optional :: regridArea
    integer,                    intent(out)           :: rc
!
!** local vars
    logical                 :: ioCapable
    logical                 :: doItFlag
    character(len=64)       :: lfileName
    character(len=64)       :: gridName
    type(ESMF_Array)        :: array
    type(ESMF_ArrayBundle)  :: arraybundle
    logical                 :: isPresent
    integer                 :: dimCount
    integer                 :: dimIndex
    integer,allocatable     :: coordDimCount(:)
    integer                 :: coordDimMax
    integer                 :: stat
    logical                 :: lnclScript
    logical                 :: hasCorners
    logical                     :: lRegridArea
    type(ESMF_Field)            :: areaField
    type(ESMF_FieldStatus_Flag) :: areaFieldStatus
!
!!
!
    ioCapable = (ESMF_IO_PIO_PRESENT .and. &
      (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT))

    doItFlag = .true. ! default
    if (present(relaxedFlag)) then
      doItFlag = .not.relaxedflag .or. (relaxedflag.and.ioCapable)
    endif

    if (doItFlag) then

      if (present(fileName)) then
        lfileName = trim(fileName)
      else
        call ESMF_GridGet(grid, name=gridName, rc=rc); ESMF_ERR_RETURN(rc)
        lfileName = trim(gridName)//".nc"
      endif

      if (present(regridArea)) then
        lRegridArea = regridArea
      else
        lRegridArea = .FALSE.
      endif

      arraybundle = ESMF_ArrayBundleCreate(rc=rc); ESMF_ERR_RETURN(rc)

      ! -- centers --

      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        isPresent=isPresent, rc=rc); ESMF_ERR_RETURN(rc)
      if (isPresent) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, &
          rc=rc); ESMF_ERR_RETURN(rc)
        call ESMF_ArraySet(array, name="lon_center", rc=rc); ESMF_ERR_RETURN(rc)
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), &
          rc=rc); ESMF_ERR_RETURN(rc)
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, &
          rc=rc); ESMF_ERR_RETURN(rc)
        call ESMF_ArraySet(array, name="lat_center", rc=rc); ESMF_ERR_RETURN(rc)
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), &
          rc=rc); ESMF_ERR_RETURN(rc)
      endif

      ! -- corners --

      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        isPresent=hasCorners, rc=rc); ESMF_ERR_RETURN(rc)
      if (hasCorners) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
          call ESMF_ArraySet(array, name="lon_corner", &
            rc=rc); ESMF_ERR_RETURN(rc)
          call ESMF_ArrayBundleAdd(arraybundle,(/array/), &
            rc=rc); ESMF_ERR_RETURN(rc)
        endif
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
          call ESMF_ArraySet(array, name="lat_corner", &
            rc=rc); ESMF_ERR_RETURN(rc)
          call ESMF_ArrayBundleAdd(arraybundle,(/array/), &
            rc=rc); ESMF_ERR_RETURN(rc)
        endif
        if (lRegridArea) then
          areaField = ESMF_FieldCreate(grid=grid, typekind=ESMF_TYPEKIND_R8, &
            rc=rc); ESMF_ERR_RETURN(rc)
          call ESMF_FieldRegridGetArea(areaField, rc=rc); ESMF_ERR_RETURN(rc)
          call ESMF_FieldGet(areaField, array=array, rc=rc)
          if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
            call ESMF_ArraySet(array, name="regrid_area", &
              rc=rc); ESMF_ERR_RETURN(rc)
            call ESMF_ArrayBundleAdd(arraybundle,(/array/), &
              rc=rc); ESMF_ERR_RETURN(rc)
          endif
        endif
      endif

      ! -- mask --

      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, &
        rc=rc); ESMF_ERR_RETURN(rc)
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc); ESMF_ERR_RETURN(rc)
        call ESMF_ArraySet(array, name="mask", rc=rc); ESMF_ERR_RETURN(rc)
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), &
          rc=rc); ESMF_ERR_RETURN(rc)
      endif

      ! -- area --

      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, &
        rc=rc); ESMF_ERR_RETURN(rc)
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc); ESMF_ERR_RETURN(rc)
        call ESMF_ArraySet(array, name="area", rc=rc); ESMF_ERR_RETURN(rc)
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), &
          rc=rc); ESMF_ERR_RETURN(rc)
      endif

      call ESMF_ArrayBundleWrite(arraybundle, &
        fileName=trim(lfileName),rc=rc); ESMF_ERR_RETURN(rc)

      if (lRegridArea) then
        call ESMF_FieldGet(areaField, status=areaFieldStatus, &
          rc=rc); ESMF_ERR_RETURN(rc)
        if (areaFieldStatus.eq.ESMF_FIELDSTATUS_COMPLETE) then
          call ESMF_FieldDestroy(areaField, rc=rc); ESMF_ERR_RETURN(rc)
        endif
      endif

      call ESMF_ArrayBundleDestroy(arraybundle,rc=rc); ESMF_ERR_RETURN(rc)
    endif

  end subroutine write_grid_netcdf

!----------------------------------------------------------------------------------------
  subroutine get_global_attr(fldbundle, ncid, rc)
    type(ESMF_FieldBundle), intent(in) :: fldbundle
    integer, intent(in)                :: ncid
    integer, intent(out)               :: rc

! local variable
    integer :: i, attcount
    integer :: ncerr
    character(len=ESMF_MAXSTR) :: attName
    type(ESMF_TypeKind_Flag)   :: typekind

    integer :: varival
    real(ESMF_KIND_R4) :: varr4val
    real(ESMF_KIND_R4), dimension(:), allocatable :: varr4list
    real(ESMF_KIND_R8) :: varr8val
    real(ESMF_KIND_R8), dimension(:), allocatable :: varr8list
    integer :: itemCount
    character(len=ESMF_MAXSTR) :: varcval
!
    call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, Count=attcount, &
                           rc=rc); ESMF_ERR_RETURN(rc)

    do i=1,attCount

      call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
                             typekind=typekind, itemCount=itemCount, rc=rc); ESMF_ERR_RETURN(rc)

      if (typekind==ESMF_TYPEKIND_I4) then
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), value=varival, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), varival); NC_ERR_STOP(ncerr)

      else if (typekind==ESMF_TYPEKIND_R4) then
         allocate (varr4list(itemCount))
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), valueList=varr4list, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), varr4list); NC_ERR_STOP(ncerr)
         deallocate(varr4list)

      else if (typekind==ESMF_TYPEKIND_R8) then
         allocate (varr8list(itemCount))
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), valueList=varr8list, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), varr8list); NC_ERR_STOP(ncerr)
         deallocate(varr8list)

      else if (typekind==ESMF_TYPEKIND_CHARACTER) then
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), value=varcval, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), trim(varcval)); NC_ERR_STOP(ncerr)

      end if

    end do

  end subroutine get_global_attr
!
!----------------------------------------------------------------------------------------
  subroutine get_grid_attr(grid, prefix, ncid, varid, rc)
    type(ESMF_Grid), intent(in)  :: grid
    character(len=*), intent(in) :: prefix
    integer, intent(in)          :: ncid
    integer, intent(in)          :: varid
    integer, intent(out)         :: rc

! local variable
    integer :: i, attcount, n, ind
    integer :: ncerr
    character(len=ESMF_MAXSTR) :: attName
    type(ESMF_TypeKind_Flag)   :: typekind

    integer :: varival
    real(ESMF_KIND_R4) :: varr4val
    real(ESMF_KIND_R8) :: varr8val
    character(len=ESMF_MAXSTR) :: varcval
!
    call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, Count=attcount, &
                           rc=rc); ESMF_ERR_RETURN(rc)

    !write(0,*)'grid attcount = ', attcount
    do i=1,attCount

      call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
                             typekind=typekind, itemCount=n, rc=rc); ESMF_ERR_RETURN(rc)
      !write(0,*)'grid att = ',i,trim(attName), ' itemCount = ' , n

      if (index(trim(attName), trim(prefix)//":")==1) then
         ind = len(trim(prefix)//":")

         if (typekind==ESMF_TYPEKIND_I4) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varival, rc=rc); ESMF_ERR_RETURN(rc)
            ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), varival); NC_ERR_STOP(ncerr)

         else if (typekind==ESMF_TYPEKIND_R4) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varr4val, rc=rc); ESMF_ERR_RETURN(rc)
            ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), varr4val); NC_ERR_STOP(ncerr)

         else if (typekind==ESMF_TYPEKIND_R8) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varr8val, rc=rc); ESMF_ERR_RETURN(rc)
            if (trim(attName) /= '_FillValue') then
              ! FIXME:  _FillValue must be cast to var type for recent versions
              ! of netcdf
              ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), varr8val); NC_ERR_STOP(ncerr)
            endif

         else if (typekind==ESMF_TYPEKIND_CHARACTER) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varcval, rc=rc); ESMF_ERR_RETURN(rc)
            ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), trim(varcval)); NC_ERR_STOP(ncerr)

         end if

      end if

    end do

  end subroutine get_grid_attr

  subroutine add_dim(ncid, dim_name, dimid, grid, rc)
    integer, intent(in)             :: ncid
    character(len=*), intent(in)    :: dim_name
    integer, intent(inout) :: dimid
    type(ESMF_Grid), intent(in)     :: grid
    integer, intent(out)            :: rc

! local variable
    integer :: i, attcount, n, dim_varid
    integer :: ncerr
    character(len=ESMF_MAXSTR) :: attName
    type(ESMF_TypeKind_Flag)   :: typekind

    integer, allocatable  :: valueListI(:)
    real(ESMF_KIND_R4), allocatable  :: valueListR4(:)
    real(ESMF_KIND_R8), allocatable  :: valueListR8(:)
    character(len=ESMF_MAXSTR), allocatable  :: valueListC(:)
!
    call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, name=dim_name, &
                           typekind=typekind, itemCount=n, rc=rc); ESMF_ERR_RETURN(rc)

    if ( trim(dim_name) == "time" ) then
    ncerr = nf90_def_dim(ncid, trim(dim_name), NF90_UNLIMITED, dimid); NC_ERR_STOP(ncerr)
    else
    ncerr = nf90_def_dim(ncid, trim(dim_name), n, dimid); NC_ERR_STOP(ncerr)
    end if

    if (typekind==ESMF_TYPEKIND_R8) then
       ncerr = nf90_def_var(ncid, dim_name, NF90_REAL8, dimids=(/dimid/), varid=dim_varid); NC_ERR_STOP(ncerr)
       allocate(valueListR8(n))
       call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dim_name), valueList=valueListR8, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_var(ncid, dim_varid, values=valueListR8 ); NC_ERR_STOP(ncerr)
       ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
       deallocate(valueListR8)
     else if (typekind==ESMF_TYPEKIND_R4) then
       ncerr = nf90_def_var(ncid, dim_name, NF90_REAL4, dimids=(/dimid/), varid=dim_varid); NC_ERR_STOP(ncerr)
       allocate(valueListR4(n))
       call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dim_name), valueList=valueListR4, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_var(ncid, dim_varid, values=valueListR4 ); NC_ERR_STOP(ncerr)
       ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
       deallocate(valueListR4)
     else
        write(0,*)'Error in module_write_netcdf.F90(add_dim) unknown typekind for ',trim(dim_name)
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    call get_grid_attr(grid, dim_name, ncid, dim_varid, rc)

  end subroutine add_dim
!
!----------------------------------------------------------------------------------------
  subroutine nccheck(status)
    use netcdf
    implicit none
    integer, intent (in) :: status

    if (status /= nf90_noerr) then
      write(0,*) status, trim(nf90_strerror(status))
      stop "stopped"
    end if
  end subroutine nccheck
 
  elemental real function quantized(dataIn, nbits, dataMin, dataMax)
    integer, intent(in) :: nbits
    real(4), intent(in) :: dataIn, dataMin, dataMax
    real(4) offset, scale_fact
    ! convert data to 32 bit integers in range 0 to 2**nbits-1, then cast
    ! cast back to 32 bit floats (data is then quantized in steps
    ! proportional to 2**nbits so last 32-nbits in floating
    ! point representation should be zero for efficient zlib compression).
    scale_fact = (dataMax - dataMin) / (2**nbits-1); offset = dataMin
    quantized = scale_fact*(nint((dataIn - offset) / scale_fact)) + offset
  end function quantized

end module module_write_netcdf
