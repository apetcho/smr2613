      module user_model1
      use ESMF
      use ESMF_netcdf_read
!
      implicit none
      public userm1_register

      private
      integer, save :: rootPet=0, coordX=1, coordY=2
      type(ESMF_Grid), save :: grid
      type(ESMF_Array), save  :: centers_array_lat, centers_array_lon, &
                                 corners_array_lat, corners_array_lon, &
                                 centers_array_msk
      type(ESMF_Field), save :: field1, field2
!
      real (ESMF_KIND_R8), parameter ::                                 &
     &             zero     = 0.0,                                      &
     &             two      = 2.0,                                      &
     &             half     = 0.5,                                      &
     &             pi       = 3.14159265359,                            &
     &             pi2      = two*pi,                                   &
     &             pih      = half*pi,                                  &
     &             deg2rad  = pi/180.    ! conversion for deg to rads
!
      contains
!
!-----------------------------------------------------------------------
!     Set register routines
!-----------------------------------------------------------------------
!
      subroutine userm1_register(comp, rc)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp) :: comp
      integer, intent(out) :: rc
!
      rc = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!     Register the callback routines.
!-----------------------------------------------------------------------
!
      call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, user_init, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, user_run, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, user_final, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end subroutine

      subroutine user_init(comp, importState, exportState, clock, rc)
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, k, localDECount
      integer :: srcsize, corners, rank, gaddr(4)
      integer :: iy, jx, cpus_per_dim(2), petCount, myPet
      integer, allocatable :: dims(:)
      character(ESMF_MAXSTR) :: grid1_fname, grid1_title
      real(ESMF_KIND_R8), allocatable :: global_grid_center_lat(:)
      real(ESMF_KIND_R8), allocatable :: global_grid_center_lon(:)
      real(ESMF_KIND_R8), allocatable :: global_grid_corner_lat(:,:)
      real(ESMF_KIND_R8), allocatable :: global_grid_corner_lon(:,:)
      integer(ESMF_KIND_I4), allocatable :: global_mask_center(:)
      real(ESMF_KIND_R8), pointer :: reshaped_grid_center_lat(:,:)
      real(ESMF_KIND_R8), pointer :: reshaped_grid_center_lon(:,:)
      real(ESMF_KIND_R8), pointer :: reshaped_grid_corner_lat(:,:)
      real(ESMF_KIND_R8), pointer :: reshaped_grid_corner_lon(:,:)
      integer(ESMF_KIND_I4), pointer :: reshaped_mask_center(:,:)
      real(ESMF_KIND_R8), dimension(:,:), pointer :: ptr1, ptr2, latPtr, lonPtr
      real(ESMF_KIND_R8) :: MISSING_R8 = 1.0d20
      integer(ESMF_KIND_I4) :: tile(2)
      logical :: file_exists
!
      type(ESMF_Decomp_Flag) :: decompflag(2)
      type(ESMF_DistGrid) :: distGrid
      type(ESMF_VM) :: vm
      type(ESMF_ArraySpec) :: arraySpec
      type(ESMF_Config) :: cf
!
      rc = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!     Query component 
!-----------------------------------------------------------------------
!
      call ESMF_GridCompGet(comp, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_VMGet(vm, localPet=myPet, petCount=petCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Query namelist.rc to get tile distribution and SCRIP grid file 
!-----------------------------------------------------------------------
!
      inquire(file='namelist.rc', exist=file_exists)
      if (file_exists) then
        cf = ESMF_ConfigCreate(rc=rc)
        call ESMF_ConfigLoadFile(cf, 'namelist.rc', rc=rc)
        call ESMF_ConfigFindLabel(cf, 'gcomp1_tile:', rc=rc)
        do i = 1, 2
          call ESMF_ConfigGetAttribute(cf, tile(i), rc=rc)
          if (myPet == rootPet) then
            print*, 'tile(', i, ') = ', tile(i)
          end if
        end do
        call ESMF_ConfigFindLabel(cf, 'gcomp1_grid:', rc=rc)
        call ESMF_ConfigGetAttribute(cf, grid1_fname, rc=rc)
      end if
!
!-----------------------------------------------------------------------
!     Read grid information 
!-----------------------------------------------------------------------
!
      call netcdf_read_grid_meta(trim(grid1_fname), grid1_title,              &
                                 rank, srcsize, corners)
      if (myPet == rootPet) print*, 'grid 1:', trim(grid1_fname)
!
      if (.not. allocated(dims)) allocate(dims(rank))
      call netcdf_read_grid_dims(trim(grid1_fname), rank, dims)
      if (myPet == rootPet) print*, 'grid 1 dims:', dims
!
      gaddr(1) = corners
      gaddr(2) = 2*corners
      gaddr(3) = 2*corners+1
      gaddr(4) = 2*corners+2
!
      allocate(global_mask_center(srcsize))
      allocate(global_grid_center_lat(srcsize))
      allocate(global_grid_center_lon(srcsize))
      allocate(global_grid_corner_lat(srcsize,corners))
      allocate(global_grid_corner_lon(srcsize,corners))
      allocate(reshaped_mask_center(dims(1), dims(2)))
      allocate(reshaped_grid_center_lat(dims(1), dims(2)))
      allocate(reshaped_grid_center_lon(dims(1), dims(2)))
      allocate(reshaped_grid_corner_lat(dims(1), dims(2)))
      allocate(reshaped_grid_corner_lon(dims(1), dims(2)))
!
      if (myPet == rootPet) then
        call netcdf_read_grid_data(trim(grid1_fname), rank, srcsize, corners, dims, &
        global_grid_corner_lat, global_grid_corner_lon, &
        global_grid_center_lat, global_grid_center_lon, global_mask_center)
      endif
!
      if (myPet == rootPet) then
        do j=1,dims(2)
          do i=1,dims(1)
            k = (j-1)*dims(1) + i
            reshaped_grid_center_lat(i,j) = global_grid_center_lat(k)
            reshaped_grid_center_lon(i,j) = global_grid_center_lon(k)
            reshaped_mask_center(i,j) = global_mask_center(k)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!     Create Grid
!-----------------------------------------------------------------------
!
      decompflag = (/ ESMF_DECOMP_RESTLAST, ESMF_DECOMP_RESTLAST /)
!
      distGrid = ESMF_DistGridCreate(minIndex=(/ 1, 1 /),               &
                                     maxIndex=(/ dims(1),dims(2)/),             &
                                     regDecomp=(/tile(1),tile(2)/), &
                                     decompflag=decompflag,             &
                                     rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      grid = ESMF_GridCreate(distgrid=distGrid,          &
                                            indexflag=ESMF_INDEX_GLOBAL,&
                                            name="atm_grid",            &
                                            rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Allocate coordinates 
!-----------------------------------------------------------------------
!
      call ESMF_GridAddCoord(grid, staggerEdgeLWidth=(/0,0/), staggerEdgeUWidth=(/0,0/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_GridAddItem (grid, itemflag=ESMF_GRIDITEM_MASK, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get pointers and set coordinates for the grid 
!-----------------------------------------------------------------------
! 
      call ESMF_GridGetCoord(grid, coordDim=coordX, array=centers_array_lat, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_GridGetCoord(grid, coordDim=coordY, array=centers_array_lon, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_GridGetItem (grid, itemflag=ESMF_GRIDITEM_MASK, array=centers_array_msk, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Scatter the reshaped grid  
!-----------------------------------------------------------------------
!
      call ESMF_ArrayScatter(centers_array_lat, reshaped_grid_center_lat, &
                             rootPet=rootPet, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_ArrayScatter(centers_array_lon, reshaped_grid_center_lon, &
                             rootPet=rootPet, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_ArrayScatter(centers_array_msk, reshaped_mask_center,     &
                             rootPet=rootPet, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set array descriptor
!-----------------------------------------------------------------------
!
      call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8,      &
                             rank=2, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create import and export field and add it to the state 
!-----------------------------------------------------------------------
!
      field1 = ESMF_FieldCreate(grid, arraySpec, name='field', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_FieldGet(field1, localDe=0, farrayPtr=ptr1, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      ptr1 = MISSING_R8
!
      if (associated(ptr1)) then
        nullify(ptr1)
      end if
!
      field2 = ESMF_FieldCreate(grid, arraySpec, name='field', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_FieldGet(field2, localDe=0, farrayPtr=ptr2, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      ptr2 = MISSING_R8
!
      call ESMF_ArrayGet(centers_array_lat, localDe=0, farrayPtr=latPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_ArrayGet(centers_array_lon, localDe=0, farrayPtr=lonPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!     pseudo-spherical harmonic l=32, m=16
!
      do j = lbound(ptr2,2), ubound(ptr2,2)
        do i = lbound(ptr2,1), ubound(ptr2,1)
          ptr2(i,j) = 2.0 + sin(2.0*latPtr(i,j)*deg2rad)**16.0 * cos(16.0*lonPtr(i,j)*deg2rad)
        end do
      end do
!
      if (associated(ptr2)) then
        nullify(ptr2)
      end if
!
      call ESMF_FieldWrite(field1, 'gcomp1_data_imp.nc',   &
                           variableName='data', overwrite=.true.,   &
                           rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_FieldWrite(field2, 'gcomp1_data_exp.nc',   &
                           variableName='data', overwrite=.true.,   &
                           rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Add fields to state
!-----------------------------------------------------------------------
!
      call ESMF_StateAdd(importState, (/field1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_StateAdd(exportState, (/field2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Debug: write out component grid in VTK format 
!-----------------------------------------------------------------------
!
      call ESMF_GridWriteVTK(grid, filename="atm_grd", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      end subroutine user_init

      subroutine user_run(comp, importState, exportState, clock, rc)
      implicit none
!     
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: petCount, myPet
!
      type(ESMF_VM) :: vm
!
      rc = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!     Query component 
!-----------------------------------------------------------------------
!
      call ESMF_GridCompGet(comp, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_VMGet(vm, localPet=myPet, petCount=petCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)

      end subroutine user_run
!
      subroutine user_final(comp, importState, exportState, clock, rc)
      implicit none
!     
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: petCount, myPet
!
      type(ESMF_VM) :: vm
!
      rc = ESMF_SUCCESS

      end subroutine user_final
      end module user_model1
