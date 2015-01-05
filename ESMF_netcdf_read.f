! $Id$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This module reads in grids stored in the netCDF form required by SCRIP
!
!-----------------------------------------------------------------------
      module ESMF_netcdf_read

!-----------------------------------------------------------------------

      use ESMF    

      implicit none

      include 'netcdf.inc'


!***********************************************************************

      contains

!***********************************************************************

      subroutine netcdf_read_grid_meta(
     &                                 grid_file,
     &                                 grid_title,
     &                                 grid_rank,
     &                                 grid_size,
     &                                 grid_corners )

!------------------------------------------------------------------------
!     call arguments
!------------------------------------------------------------------------

      character (ESMF_MAXSTR), intent(in) ::
     &        grid_file             ! filename of source grid file

      character (ESMF_MAXSTR), intent(out) ::
     &        grid_title            ! filename of source grid file

      integer, intent(out) ::
     &       grid_size,             ! total points on each grid
     &       grid_rank,             ! rank of each grid
     &       grid_corners           ! number of corners

!------------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------------

      integer  ::
     &       ncstat,                ! netCDF status variable
     &       nc_grid_id,            ! netCDF grid file id
     &       nc_gridsize_id,        ! netCDF grid size dim id
     &       nc_gridcorn_id,        ! netCDF grid corner dim id
     &       nc_gridrank_id         ! netCDF grid rank dim id

!-----------------------------------------------------------------------
!     open grid files and read  grid meta data
!-----------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! open netcdf file
      !-----------------------------------------------------------------

      ncstat = nf_open(grid_file, NF_NOWRITE, nc_grid_id)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! get ID of dimension and length
      !-----------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_grid_id, 'grid_size', nc_gridsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_grid_id, nc_gridsize_id, grid_size)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! get ID of 'grid_rank' dimension and length
      !-----------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_grid_id, 'grid_rank', nc_gridrank_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_grid_id, nc_gridrank_id, grid_rank)
      call netcdf_error_handler(ncstat)
      !-----------------------------------------------------------------
      ! get ID of 'grid_corners' dimension and length
      !-----------------------------------------------------------------
      ncstat = nf_inq_dimid(nc_grid_id,'grid_corners',nc_gridcorn_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_grid_id,nc_gridcorn_id,grid_corners)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! get 'title' of grid
      !-----------------------------------------------------------------
      ncstat = nf_get_att_text(nc_grid_id, nf_global,'title',grid_title)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------
!     close input file     
!------------------------------------------------------------------------

      ncstat = nf_close(nc_grid_id)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------


      end subroutine netcdf_read_grid_meta


!=======================================================================

!=======================================================================

      subroutine netcdf_read_grid_data(
     &                                 grid_file,
     &                                 grid_rank, 
     &                                 grid_size, 
     &                                 grid_corners,
     &                                 grid_dims,
     &                                 grid_corner_lat,
     &                                 grid_corner_lon,
     &                                 grid_center_lat,
     &                                 grid_center_lon,
     &                                 grid_center_msk )

!------------------------------------------------------------------------
!     call arguments
!------------------------------------------------------------------------

      character (ESMF_MAXSTR), intent(in) ::
     &        grid_file             ! filename of source grid file

      integer, intent(in) ::
     &       grid_size,             ! total points on each grid
     &       grid_rank,             ! rank of each grid
     &       grid_corners           ! number of corners

      integer, intent(out), dimension(grid_rank) ::
     &       grid_dims              ! total points on each grid


      real (ESMF_KIND_R8), intent(out), dimension(grid_size) ::
     &       grid_center_lat,
     &       grid_center_lon 

      real (ESMF_KIND_R8),intent(out), dimension(grid_size,grid_rank) ::
     &       grid_corner_lat,
     &       grid_corner_lon 

      integer (ESMF_KIND_I4), intent(out), dimension(grid_size) ::
     &       grid_center_msk

!------------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------------

      character (ESMF_MAXSTR) :: 
     &             grid_units       ! units for grid coords (degs/radians)

      integer  ::
     &       ncstat,                ! netCDF status variable
     &       nc_grid_id,            ! netCDF grid file id
     &       nc_griddims_id,        ! netCDF grid dims id
     &       nc_grdcrnrlat_id,      ! netCDF grid corner lat var id
     &       nc_grdcrnrlon_id,      ! netCDF grid corner lon var id
     &       nc_grdcntrlat_id,      ! netCDF grid center lat var id
     &       nc_grdcntrlon_id,      ! netCDF grid center lon var id
     &       nc_grdcntrmsk_id       ! netCDF grid center lon var id

!-----------------------------------------------------------------------
! local constants
!-----------------------------------------------------------------------

      real (ESMF_KIND_R8), parameter :: 
     &             zero     = 0.0,
     &             two      = 2.0,
     &             half     = 0.5,
     &             pi       = 3.14159265359,
     &             pi2      = two*pi,
     &             pih      = half*pi,
     &             deg2rad  = pi/180.    ! conversion for deg to rads

!-----------------------------------------------------------------------
!     open grid files and read grid data
!-----------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! open netcdf file
      !-----------------------------------------------------------------

      ncstat = nf_open(grid_file, NF_NOWRITE, nc_grid_id)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------
!     Read grid dimension data
!------------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! get variable ID for 'grid_dims' and extract
      !-----------------------------------------------------------------

      ncstat = nf_inq_varid(nc_grid_id, 'grid_dims', nc_griddims_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_int(nc_grid_id, nc_griddims_id, grid_dims)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------
!     Read grid fields
!------------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! get variable ID for 'grid_center_lat' and extract
      !-----------------------------------------------------------------

      ncstat = nf_inq_varid(nc_grid_id, 'grid_center_lat',
     &                                   nc_grdcntrlat_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_double(nc_grid_id, nc_grdcntrlat_id,
     &                                       grid_center_lat)
      call netcdf_error_handler(ncstat)


      !-----------------------------------------------------------------
      ! get variable ID for 'grid_center_lon' and extract
      !-----------------------------------------------------------------

      ncstat = nf_inq_varid(nc_grid_id, 'grid_center_lon',
     &                                   nc_grdcntrlon_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_double(nc_grid_id, nc_grdcntrlon_id,
     &                                       grid_center_lon)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! get variable ID for 'grid_imask' and extract
      !-----------------------------------------------------------------

      ncstat = nf_inq_varid(nc_grid_id, 'grid_imask',
     &                                   nc_grdcntrmsk_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_int(nc_grid_id, nc_grdcntrmsk_id,
     &                                    grid_center_msk)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! get variable ID for 'grid_corner_lat'
      !-----------------------------------------------------------------

      ncstat = nf_inq_varid(nc_grid_id, 'grid_corner_lat',
     &                                   nc_grdcrnrlat_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_double(nc_grid_id, nc_grdcrnrlat_id,
     &                                       grid_corner_lat)
      call netcdf_error_handler(ncstat)


      !-----------------------------------------------------------------
      ! get variable ID for 'grid_corner_lon'
      !-----------------------------------------------------------------

      ncstat = nf_inq_varid(nc_grid_id, 'grid_corner_lon',
     &                                   nc_grdcrnrlon_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_double(nc_grid_id, nc_grdcrnrlon_id,
     &                                       grid_corner_lon)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------
!     convert lat/lon to radians if units are not already
!------------------------------------------------------------------------
      grid_units = ' '
      ncstat = nf_get_att_text(nc_grid_id, nc_grdcntrlat_id, 'units',
     &                         grid_units)
      call netcdf_error_handler(ncstat)

      select case (grid_units(1:7))
      case ('degrees')

        !grid_center_lat = grid_center_lat*deg2rad
        !grid_center_lon = grid_center_lon*deg2rad

      case ('radians')

      !-----------------------------------------------------------------
      !*** no conversion necessary
      !-----------------------------------------------------------------

      case default

        print *,'unknown units supplied for grid center lat/lon: '
        print *,'proceeding assuming radians'

      end select

      grid_units = ' '
      ncstat = nf_get_att_text(nc_grid_id, nc_grdcrnrlat_id, 'units',
     &                         grid_units)
      call netcdf_error_handler(ncstat)

      select case (grid_units(1:7))
      case ('degrees')

        !grid_corner_lat = grid_corner_lat*deg2rad
        !grid_corner_lon = grid_corner_lon*deg2rad

      case ('radians')

      !-----------------------------------------------------------------
      !*** no conversion necessary
      !-----------------------------------------------------------------

      case default

        print *,'unknown units supplied for grid corner lat/lon: '
        print *,'proceeding assuming radians'

      end select

      !-----------------------------------------------------------------
      ! convert longitudes to 0,2pi interval
      !-----------------------------------------------------------------

!      where (grid_center_lon > pi2)  grid_center_lon =
!     &                                   grid_center_lon - pi2
!      where (grid_center_lon < zero) grid_center_lon =
!     &                                   grid_center_lon + pi2
!      where (grid_corner_lon > pi2)  grid_corner_lon =
!     &                                   grid_corner_lon - pi2
!      where (grid_corner_lon < zero) grid_corner_lon =
!     &                                   grid_corner_lon + pi2

      !-----------------------------------------------------------------
      ! make sure input latitude range is within the machine values
      ! for +/- pi/2
      !-----------------------------------------------------------------

!      where (grid_center_lat >  pih) grid_center_lat =  pih
!      where (grid_corner_lat >  pih) grid_corner_lat =  pih
!      where (grid_center_lat < -pih) grid_center_lat = -pih
!      where (grid_corner_lat < -pih) grid_corner_lat = -pih

!------------------------------------------------------------------------
!     close input file
!------------------------------------------------------------------------

      ncstat = nf_close(nc_grid_id)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------


      end subroutine netcdf_read_grid_data


!=======================================================================

!=======================================================================

      subroutine netcdf_read_grid_dims(
     &                                 grid_file,
     &                                 grid_rank, 
     &                                 grid_dims)

!------------------------------------------------------------------------
!     call arguments
!------------------------------------------------------------------------

      character (ESMF_MAXSTR), intent(in) ::
     &        grid_file             ! filename of source grid file

      integer, intent(in) ::
     &       grid_rank              ! rank of each grid

      integer, intent(out), dimension(grid_rank) ::
     &       grid_dims              ! total points on each grid

!------------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------------

      integer  ::
     &       ncstat,                ! netCDF status variable
     &       nc_grid_id,            ! netCDF grid file id
     &       nc_griddims_id         ! netCDF grid dims id

!-----------------------------------------------------------------------
!     open grid files and read grid data
!-----------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! open netcdf file
      !-----------------------------------------------------------------

      ncstat = nf_open(grid_file, NF_NOWRITE, nc_grid_id)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------
!     Read grid dimension data
!------------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! get variable ID for 'grid_dims' and extract
      !-----------------------------------------------------------------

      ncstat = nf_inq_varid(nc_grid_id, 'grid_dims', nc_griddims_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_int(nc_grid_id, nc_griddims_id, grid_dims)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------
!     close input file
!------------------------------------------------------------------------

      ncstat = nf_close(nc_grid_id)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------


      end subroutine netcdf_read_grid_dims


!=======================================================================

!=======================================================================

      subroutine netcdf_read_remap_meta(
     &                  remap_file,
     &                  num_links, 
     &                  num_wts,
     &                  src_grid_file_name,
     &                  dst_grid_file_name,
     &                  src_grid_meta,
     &                  dst_grid_meta )


!------------------------------------------------------------------------
!     call arguments
!------------------------------------------------------------------------

      character (ESMF_MAXSTR), intent(in) ::
     &       remap_file            ! filename of reampping weights file

      integer, intent(out) ::
     &       num_links,
     &       num_wts

      character (ESMF_MAXSTR), intent(out) ::
     &       src_grid_file_name,   ! filename of source grid file
     &       dst_grid_file_name   ! filename of destination grid file

      integer, intent(out), dimension(10) ::
     &       src_grid_meta,         ! 
     &       dst_grid_meta          ! 

!------------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------------

      integer ::
     &       ncstat,                ! netCDF status variable
     &       nc_file_id,            ! netCDF reamp file id
     &       nc_numlinks_id,        ! netCDF number of links id
     &       nc_numwgts_id,         ! netCDF number of wts id
     &       nc_gridsize_id,        ! netCDF grid size dim id
     &       nc_gridcorn_id,        ! netCDF grid corner dim id
     &       nc_gridrank_id,        ! netCDF grid rank dim id
     &       nc_griddims_id        ! netCDF grid dimension size id

      integer ::
     &       grid_rank,
     &       grid_size,
     &       grid_corners

      integer, allocatable ::
     &       grid_dims(:)

!-----------------------------------------------------------------------
!     open remap file and read meta data
!-----------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! open netcdf file
      !-----------------------------------------------------------------

      ncstat = nf_open(remap_file, NF_NOWRITE, nc_file_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
! read source grid meta data for consistency check
!-----------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! number of address pairs in the remappings
      !-----------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'n_s', nc_numlinks_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numlinks_id, num_links)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! number of weights per point/order of interpolation method
      !-----------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'num_wgts', nc_numwgts_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numwgts_id, num_wts)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
! read source grid meta data for consistency check
!-----------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! Source; Grid Size
      !-----------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'n_a',nc_gridsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id,nc_gridsize_id, grid_size)
      call netcdf_error_handler(ncstat)
      src_grid_meta(1) = grid_size

      !-----------------------------------------------------------------
      ! Source; number of Grid Corners
      !-----------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'nv_a', 
     &                                                 nc_gridcorn_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_gridcorn_id, grid_corners)
      call netcdf_error_handler(ncstat)
      src_grid_meta(2) = grid_corners

      !-----------------------------------------------------------------
      ! Source; Grid Rank
      !-----------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'src_grid_rank', nc_gridrank_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_gridrank_id, grid_rank)
      call netcdf_error_handler(ncstat)
      src_grid_meta(3) = grid_rank

      !-----------------------------------------------------------------
      ! Source; Grid Dims
      !-----------------------------------------------------------------

      allocate(grid_dims(grid_rank) )
      ncstat = nf_inq_varid(nc_file_id, 'src_grid_dims',nc_griddims_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_int(nc_file_id,nc_griddims_id,grid_dims)
      call netcdf_error_handler(ncstat)
      src_grid_meta(4:4+grid_rank-1) = grid_dims(1:grid_rank)
      deallocate(grid_dims )

      !-----------------------------------------------------------------
      ! get source grid name
      !-----------------------------------------------------------------

      ncstat = nf_get_att_text(nc_file_id, nf_global,'grid_file_src',
     &                                              src_grid_file_name)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
! read destination grid meta data for consistency check
!-----------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! Destination; Grid Size
      !-----------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'n_b',nc_gridsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id,nc_gridsize_id, grid_size)
      call netcdf_error_handler(ncstat)
      dst_grid_meta(1) = grid_size

      !-----------------------------------------------------------------
      ! Destination; number of Grid Corners
      !-----------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'nv_b', 
     &                                                 nc_gridcorn_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_gridcorn_id, grid_corners)
      call netcdf_error_handler(ncstat)
      dst_grid_meta(2) = grid_corners

      !-----------------------------------------------------------------
      ! Destination; Grid Rank
      !-----------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'dst_grid_rank', nc_gridrank_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_gridrank_id, grid_rank)
      call netcdf_error_handler(ncstat)
      dst_grid_meta(3) = grid_rank


      !-----------------------------------------------------------------
      ! Destination; Grid Dims
      !-----------------------------------------------------------------

      allocate(grid_dims(grid_rank) )
      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_dims',nc_griddims_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_int(nc_file_id,nc_griddims_id,grid_dims)
      call netcdf_error_handler(ncstat)
      dst_grid_meta(4:4+grid_rank-1) = grid_dims(1:grid_rank)
      deallocate(grid_dims )

      !-----------------------------------------------------------------
      ! get destination grid name
      !-----------------------------------------------------------------

      ncstat = nf_get_att_text(nc_file_id, nf_global, 'grid_file_dst',
     &                                              dst_grid_file_name)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------
!     close input file
!------------------------------------------------------------------------

      ncstat = nf_close(nc_file_id)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------

      end subroutine netcdf_read_remap_meta

!=======================================================================

!=======================================================================

      subroutine netcdf_read_remap_wts(
     &                             remap_file,
     &                             num_links,
     &                             num_wts,
     &                             map_wts,
     &                             mapIndxList)
!------------------------------------------------------------------------
!     call arguments
!------------------------------------------------------------------------

      character (ESMF_MAXSTR), intent(in) ::
     &       remap_file            ! filename of reampping weights file

      integer, intent(in) ::
     &       num_links,
     &       num_wts

      real (ESMF_KIND_R8), intent(out) ::
     &       map_wts(num_wts,num_links)

      integer, intent(out) ::
     &       mapIndxList(2,num_links)

!------------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------------

      integer  ::
     &       ncstat,                ! netCDF status variable
     &       nc_file_id,            ! netCDF remap file id
     &       nc_dstgrdadd_id,       ! netCDF remap address id
     &       nc_srcgrdadd_id,       ! netCDF reamp address id
     &       nc_rmpmatrix_id        ! netCDF weights id

      integer, allocatable  ::
     &       address(:)

!-----------------------------------------------------------------------
!     open remap file and read remap weights
!-----------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! open netcdf file
      !-----------------------------------------------------------------

      ncstat = nf_open(remap_file, NF_NOWRITE, nc_file_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
! get remapping weights and addresses
!-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! source addresses for weights
      !-----------------------------------------------------------------

      allocate( address(num_links) )
      ncstat = nf_inq_varid(nc_file_id, 'col', nc_srcgrdadd_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_int(nc_file_id, nc_srcgrdadd_id, address)
      call netcdf_error_handler(ncstat)
      mapIndxList(1,:) = address

      !-----------------------------------------------------------------
      ! destination addresss for weights 
      !-----------------------------------------------------------------

      ncstat = nf_inq_varid(nc_file_id, 'row', nc_dstgrdadd_id)
      ncstat = nf_get_var_int(nc_file_id, nc_dstgrdadd_id, address)
      call netcdf_error_handler(ncstat)
      mapIndxList(2,:) = address
      deallocate( address )

      !-----------------------------------------------------------------
      !     read all variables
      !-----------------------------------------------------------------

      ncstat = nf_inq_varid(nc_file_id, 'S', nc_rmpmatrix_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_get_var_double(nc_file_id, nc_rmpmatrix_id, map_wts)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------
!     close input file
!------------------------------------------------------------------------

      ncstat = nf_close(nc_file_id)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------


      end subroutine netcdf_read_remap_wts

!=======================================================================

!=======================================================================

      subroutine netcdf_write_test(
     &                  grid_rank,
     &                  grid_dims,
     &                  array1, array2)
!------------------------------------------------------------------------
!     call arguments
!------------------------------------------------------------------------
      integer, intent(in) ::
     &    grid_rank

      integer, dimension(:), intent(in) ::
     &    grid_dims
 
      real (ESMF_KIND_R8), dimension(:,:), intent(in) ::
     &    array1, array2

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
      character (ESMF_MAXSTR) ::
     &        map_name,      ! extra stuff
     &        output_file,   ! filename for test results
     &        dim_name       ! netCDF dimension name

      integer :: n
     
      integer, save ::    ! netCDF ids for files and arrays
     &        ncstat, nc_outfile_id,
     &        nc_srcgrdcntrlat_id, nc_srcgrdcntrlon_id,
     &        nc_dstgrdcntrlat_id, nc_dstgrdcntrlon_id,
     &        nc_forerror_id, nc_forrelerror_id,
     &        nc_backerror_id, nc_backrelerror_id,
     &        nc_srcgrdfield_id, nc_dstgrdfield_id
     
      integer, dimension(:), allocatable ::
     &        nc_grid1size_id, nc_grid2size_id

      integer, save :: i=0

      map_name = "Regrid Error"
      output_file="ESMF_Regrid_Results.nc"
!-----------------------------------------------------------------------
      i = i+1 !------------------------first time counter---------------
!-----------------------------------------------------------------------
 
      if (i .eq. 1) then 
!-----------------------------------------------------------------------
!     setup a NetCDF file for output
!-----------------------------------------------------------------------
      !-----------------------------------------------------------------
      ! create netCDF dataset
      !-----------------------------------------------------------------
      ncstat = nf_create (output_file, NF_CLOBBER, nc_outfile_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, NF_GLOBAL, 'title',
     &                          len_trim(map_name), map_name)
      call netcdf_error_handler(ncstat)

      ncstat = nf_close(nc_outfile_id)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! define grid size dimensions
      !-----------------------------------------------------------------
      ncstat = nf_open(output_file, NF_WRITE, nc_outfile_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_redef(nc_outfile_id)
      call netcdf_error_handler(ncstat)

      allocate( nc_grid1size_id(grid_rank) )

      do n=1,grid_rank
        write(dim_name,1000) 'grid1_dim',n
        ncstat = nf_def_dim (nc_outfile_id, dim_name,
     &                       grid_dims(n), nc_grid1size_id(n))
        call netcdf_error_handler(ncstat)
      end do
      !-----------------------------------------------------------------
      ! define grid center latitude array
      !-----------------------------------------------------------------
      ncstat = nf_def_var (nc_outfile_id, 'src_grid_center_lat',
     &                     NF_DOUBLE, grid_rank, nc_grid1size_id,
     &                     nc_srcgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, nc_srcgrdcntrlat_id,
     &                          'units', 7, 'radians')
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! define grid center longitude array
      !-----------------------------------------------------------------
      ncstat = nf_def_var (nc_outfile_id, 'src_grid_center_lon',
     &                     NF_DOUBLE, grid_rank, nc_grid1size_id,
     &                     nc_srcgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, nc_srcgrdcntrlon_id,
     &                          'units', 7, 'radians')
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! define grid 1 field array
      !-----------------------------------------------------------------
      ncstat = nf_def_var (nc_outfile_id, 'src_grid_field',
     &                     NF_DOUBLE, grid_rank, nc_grid1size_id,
     &                     nc_srcgrdfield_id)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! define error arrays
      !-----------------------------------------------------------------
      ncstat = nf_def_var (nc_outfile_id, 'forward_error',
     &                     NF_DOUBLE, grid_rank, nc_grid1size_id,
     &                     nc_backerror_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_outfile_id, 'forward_relative_error',
     &                     NF_DOUBLE, grid_rank, nc_grid1size_id,
     &                     nc_backrelerror_id)
      call netcdf_error_handler(ncstat)

      !------------------------------------------------------------------------
      !     cleanup
      !------------------------------------------------------------------------
      deallocate( nc_grid1size_id )

      ncstat = nf_close(nc_outfile_id)
      call netcdf_error_handler(ncstat)
      endif

      if (i .eq. 3) then
      ncstat = nf_open(output_file, NF_WRITE, nc_outfile_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_redef(nc_outfile_id)
      call netcdf_error_handler(ncstat)

      allocate(nc_grid2size_id(grid_rank))
      do n=1,grid_rank
        write(dim_name,1000) 'grid2_dim',n
        ncstat = nf_def_dim (nc_outfile_id, dim_name,
     &                       grid_dims(n), nc_grid2size_id(n))
        call netcdf_error_handler(ncstat)
      end do

      !-----------------------------------------------------------------
      ! define grid center latitude array
      !-----------------------------------------------------------------
      ncstat = nf_def_var (nc_outfile_id, 'dst_grid_center_lat',
     &                     NF_DOUBLE, grid_rank, nc_grid2size_id,
     &                     nc_dstgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, nc_dstgrdcntrlat_id,
     &                          'units', 7, 'radians')
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! define grid center longitude array
      !-----------------------------------------------------------------
      ncstat = nf_def_var (nc_outfile_id, 'dst_grid_center_lon',
     &                     NF_DOUBLE, grid_rank, nc_grid2size_id,
     &                     nc_dstgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, nc_dstgrdcntrlon_id,
     &                          'units', 7, 'radians')
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! define grid 2 field array
      !-----------------------------------------------------------------
      ncstat = nf_def_var (nc_outfile_id, 'dst_grid_field',
     &                     NF_DOUBLE, grid_rank, nc_grid2size_id,
     &                     nc_dstgrdfield_id)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! define error arrays
      !-----------------------------------------------------------------
      ncstat = nf_def_var (nc_outfile_id, 'backward_error',
     &                     NF_DOUBLE, grid_rank, nc_grid2size_id,
     &                     nc_forerror_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_outfile_id, 'backward_relative_error',
     &                     NF_DOUBLE, grid_rank, nc_grid2size_id,
     &                     nc_forrelerror_id)
      call netcdf_error_handler(ncstat)

      !------------------------------------------------------------------------
      !     cleanup
      !------------------------------------------------------------------------
      deallocate( nc_grid2size_id )

      ncstat = nf_close(nc_outfile_id)
      call netcdf_error_handler(ncstat)
      endif

 ! write format
 1000 format(a9,i1)

!-----------------------------------------------------------------------
!     write coordinate info
!-----------------------------------------------------------------------
      if (i .eq. 1) then
      ncstat = nf_open(output_file, NF_WRITE, nc_outfile_id)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! write grid center latitude array
      !-----------------------------------------------------------------
      ncstat = nf_put_var_double(nc_outfile_id, nc_srcgrdcntrlat_id,
     &                           array1)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! write grid center longitude array
      !-----------------------------------------------------------------
      ncstat = nf_put_var_double(nc_outfile_id, nc_srcgrdcntrlon_id,
     &                           array2)
      call netcdf_error_handler(ncstat)

      ncstat = nf_close(nc_outfile_id)
      call netcdf_error_handler(ncstat)

      endif
      if (i .eq. 2) then
      ncstat = nf_open(output_file, NF_WRITE, nc_outfile_id)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! write grid 1 field array
      !-----------------------------------------------------------------
      ncstat = nf_put_var_double(nc_outfile_id, nc_srcgrdfield_id,
     &                           array1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_close(nc_outfile_id)
      call netcdf_error_handler(ncstat)

      endif
      if (i .eq. 3) then
      ncstat = nf_open(output_file, NF_WRITE, nc_outfile_id)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! write grid center latitude array
      !-----------------------------------------------------------------
      ncstat = nf_put_var_double(nc_outfile_id, nc_dstgrdcntrlat_id,
     &                           array1)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! write grid center longitude array
      !-----------------------------------------------------------------
      ncstat = nf_put_var_double(nc_outfile_id, nc_dstgrdcntrlon_id,
     &                            array2)
      call netcdf_error_handler(ncstat)

      ncstat = nf_close(nc_outfile_id)
      call netcdf_error_handler(ncstat)

      endif
      if (i .eq. 4) then
      ncstat = nf_open(output_file, NF_WRITE, nc_outfile_id)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! write grid 2 field array
      !-----------------------------------------------------------------
      ncstat = nf_put_var_double(nc_outfile_id, nc_dstgrdfield_id,
     &                           array1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_close(nc_outfile_id)
      call netcdf_error_handler(ncstat)

      endif

!-----------------------------------------------------------------------
!     write test results to a NetCDF file for output
!-----------------------------------------------------------------------
      if (i .eq. 5) then
      ncstat = nf_open(output_file, NF_WRITE, nc_outfile_id)
      call netcdf_error_handler(ncstat)

      !-----------------------------------------------------------------
      ! create netCDF dataset
      !-----------------------------------------------------------------
      ncstat = nf_put_var_double(nc_outfile_id, nc_forerror_id, array1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_forrelerror_id,
     &                           array2)
      call netcdf_error_handler(ncstat)

      ncstat = nf_close(nc_outfile_id)
      call netcdf_error_handler(ncstat)

      endif
      if (i .eq. 6) then
      ncstat = nf_open(output_file, NF_WRITE, nc_outfile_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_backerror_id, array1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_backrelerror_id,
     &                           array2)
      call netcdf_error_handler(ncstat)

!------------------------------------------------------------------------
!     close input file
!------------------------------------------------------------------------

      ncstat = nf_close(nc_outfile_id)
      call netcdf_error_handler(ncstat)

      endif

      end subroutine netcdf_write_test

!=======================================================================

!=======================================================================

      subroutine netcdf_error_handler(istat)

!-----------------------------------------------------------------------
!     This routine provides a simple interface to netCDF error message
!     routine.
!-----------------------------------------------------------------------

      integer, intent(in) ::
     &       istat   ! integer status returned by netCDF function call

!-----------------------------------------------------------------------

      if (istat /= NF_NOERR) then
        print *,'Error in netCDF: ',nf_strerror(istat)
        stop
      endif

!-----------------------------------------------------------------------


      end subroutine netcdf_error_handler


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module ESMF_netcdf_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
