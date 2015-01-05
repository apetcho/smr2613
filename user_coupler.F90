      module user_coupler
      use ESMF
      use ESMF_netcdf_read
!
      implicit none
      public usercpl_register

      private
      integer, save :: rootPet=0
!
      type(ESMF_RouteHandle), save :: routeHandleF1, routeHandleF2
      type(ESMF_RouteHandle), save :: routeHandleB1, routeHandleB2
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
      real*8, parameter :: MISSING_R8 = 1.0d20
      real  , parameter :: MISSING_R4 = 1.0e20
      real*8, parameter :: TOL_R8 = MISSING_R8/2.0d0
      real  , parameter :: TOL_R4 = MISSING_R4/2.0
!
      real(ESMF_KIND_R8), parameter :: ZERO_R8 = 0.0d0
      real(ESMF_KIND_R8), parameter :: ONE_R8 = 1.0d0
!
      integer(ESMF_KIND_I4), parameter :: MAPPED_MASK = 99
      integer(ESMF_KIND_I4), parameter :: UNMAPPED_MASK = 98
!
      character(len=ESMF_MAXSTR), save :: forward_init='forward_init',        &
                                    forward_run='forward_run',          &
                                    backward_init='backward_init',      &
                                    backward_run='backward_run'
!
      contains
!
!-----------------------------------------------------------------------
!     Set register routines
!-----------------------------------------------------------------------
!
      subroutine usercpl_register(comp, rc)
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
      integer :: petCount, myPet
      integer(ESMF_KIND_I4) :: f_i=0, b_i=0
!
      type(ESMF_VM) :: vm
      type(ESMF_UnmappedAction_Flag) :: unmap
      type(ESMF_RegridMethod_Flag) :: regridMethod
      type(ESMF_Field) :: srcField, dstField, tmpField
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
!     Reconcile state objects
!-----------------------------------------------------------------------
!
!      call ESMF_StateReconcile(importState, vm=vm, rc=rc)
!      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
!          line=__LINE__, file=__FILE__))                                &
!          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!      call ESMF_StateReconcile(exportState, vm=vm, rc=rc)
!      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
!          line=__LINE__, file=__FILE__))                                &
!          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Query direction 
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet(exportState, name=forward_init, rc=f_i)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_AttributeGet(exportState, name=backward_init, rc=b_i)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      if (myPet == rootPet) then
        if (b_i == ESMF_SUCCESS) then
          print*, "**** init coupler in backward direction ****"
        end if
        if (f_i == ESMF_SUCCESS) then
          print*, "**** init coupler in forward  direction ****"
        end if
      end if
!
!-----------------------------------------------------------------------
!     Get source and destination fields
!-----------------------------------------------------------------------
!
      call ESMF_StateGet(importState, "field", srcfield, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_StateGet(exportState, "field", dstfield, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create 1st routehandle
!-----------------------------------------------------------------------
!
      unmap = ESMF_UNMAPPEDACTION_IGNORE
      regridMethod = ESMF_REGRIDMETHOD_BILINEAR
!
      if (f_i == ESMF_SUCCESS) then
      call ESMF_FieldRegridStore(srcField=srcField,                     &
                                 dstField=dstField,                     &
                                 srcMaskValues=(/0/),                   &
                                 dstMaskValues=(/0/),                   &
                                 unmappedaction=unmap,                  &
                                 routeHandle=routeHandleF1,             &
                                 regridmethod=regridMethod,             &
                                 rc=rc)
      else
      call ESMF_FieldRegridStore(srcField=srcField,                     &
                                 dstField=dstField,                     &
                                 srcMaskValues=(/0/),                   &
                                 dstMaskValues=(/0/),                   &
                                 unmappedaction=unmap,                  &
                                 routeHandle=routeHandleB1,             &
                                 regridmethod=regridMethod,             &
                                 rc=rc)
      end if
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Modify grid mask to split masked and unmasked grid cells    
!-----------------------------------------------------------------------
!
      if (f_i == ESMF_SUCCESS) then
        call UTIL_FindUnmapped(srcField, dstField, 0, 0, .true., rc)
      else
        call UTIL_FindUnmapped(srcField, dstField, 0, 0, .false., rc)
      end if
!
!-----------------------------------------------------------------------
!     Create temporary field in destination grid
!-----------------------------------------------------------------------
!
      tmpField = UTIL_FieldCreate(dstField, 'field', ONE_R8, -1, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create 2nd routehandle
!-----------------------------------------------------------------------
!
      unmap = ESMF_UNMAPPEDACTION_IGNORE
      regridMethod = ESMF_REGRIDMETHOD_NEAREST_STOD
!
      if (f_i == ESMF_SUCCESS) then
      call ESMF_FieldRegridStore(srcField=tmpField,                     &
                                 dstField=dstField,                     &
                                 srcMaskValues=(/ 0, UNMAPPED_MASK /),  &
                                 dstMaskValues=(/ 0, MAPPED_MASK /),    &
                                 !dstMaskValues=(/ MAPPED_MASK /),       &
                                 unmappedaction=unmap,                  &
                                 routeHandle=routeHandleF2,             &
                                 regridmethod=regridMethod,             &
                                 rc=rc)
      else
      call ESMF_FieldRegridStore(srcField=tmpField,                     &
                                 dstField=dstField,                     &
                                 srcMaskValues=(/ 0, UNMAPPED_MASK /),  &
                                 dstMaskValues=(/ 0, MAPPED_MASK /),    &
                                 !dstMaskValues=(/ MAPPED_MASK /),       &
                                 unmappedaction=unmap,                  &
                                 routeHandle=routeHandleB2,             &
                                 regridmethod=regridMethod,             &
                                 rc=rc)
      end if
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Delete temporary field    
!-----------------------------------------------------------------------
!
      call ESMF_FieldDestroy(tmpField, rc=rc)
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
      integer(ESMF_KIND_I4) :: f_r=0, b_r=0
!
      type(ESMF_VM) :: vm
      type(ESMF_Field) :: srcField, dstField, tmpField
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
!     Query direction 
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet(exportState, name=forward_run, rc=f_r)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_AttributeGet(exportState, name=backward_run, rc=b_r)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
       
      if (myPet == rootPet) then
        if (b_r == ESMF_SUCCESS) then
          print*, "**** run coupler in backward direction ****"
        end if
        if (f_r == ESMF_SUCCESS) then
          print*, "**** run coupler in forward  direction ****"
        end if
      end if
!
!-----------------------------------------------------------------------
!     Get source and destination fields
!-----------------------------------------------------------------------
!
      call ESMF_StateGet(importState, "field", srcfield, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_StateGet(exportState, "field", dstfield, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create temporary field in destination grid
!-----------------------------------------------------------------------
!
      tmpField = UTIL_FieldCreate(dstField, "field", MISSING_R8, -1, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Perform 1st regrid operation
!-----------------------------------------------------------------------
!
      if (f_r == ESMF_SUCCESS) then      
        call ESMF_FieldRegrid(srcField, tmpField, routeHandleF1,        &
                              zeroregion=ESMF_REGION_SELECT, rc=rc)
      else
        call ESMF_FieldRegrid(srcField, tmpField, routeHandleB1,        &
                              zeroregion=ESMF_REGION_SELECT, rc=rc)
      end if
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Write result of regrid operation
!-----------------------------------------------------------------------
! 
      if (f_r == ESMF_SUCCESS) then
      call ESMF_FieldWrite(tmpField, 'remap_1_forward.nc',              &
                           variableName='data', overwrite=.true.,       &
                           rc=rc)
      else
      call ESMF_FieldWrite(tmpField, 'remap_1_backward.nc',            &
                           variableName='data', overwrite=.true.,       &
                           rc=rc)
      end if
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Copy content from temporary field to destination field 
!-----------------------------------------------------------------------
!
      call ESMF_FieldCopy(dstField, tmpField, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Perform 2nd regrid operation
!-----------------------------------------------------------------------
!
      if (f_r == ESMF_SUCCESS) then
        call ESMF_FieldRegrid(tmpField, dstField, routeHandleF2,        &
                              zeroregion=ESMF_REGION_SELECT, rc=rc)
      else
        call ESMF_FieldRegrid(tmpField, dstField, routeHandleB2,        &
                              zeroregion=ESMF_REGION_SELECT, rc=rc)
      end if
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Write result of regrid operation
!-----------------------------------------------------------------------
! 
      if (f_r == ESMF_SUCCESS) then
      call ESMF_FieldWrite(dstField, 'remap_2_forward.nc',              &
                           variableName='data', overwrite=.true.,       &
                           rc=rc)
      else
      call ESMF_FieldWrite(dstField, 'remap_2_backward.nc',             &
                           variableName='data', overwrite=.true.,       &
                           rc=rc)
      end if
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Deallocate temporary fields 
!-----------------------------------------------------------------------
!
      call ESMF_FieldDestroy(tmpField, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
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
!
      function UTIL_FieldCreate(field, fname, initVal, dstLandMask, rc)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_Field) :: UTIL_FieldCreate 
! 
      type(ESMF_Field), intent(in) :: field
      character(*), intent(in) :: fname
      real(ESMF_KIND_R8), intent(in) :: initVal
      integer(ESMF_KIND_I4), intent(in) :: dstLandMask      
      integer, intent(out) :: rc 
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, k, localDECount
      integer :: cLbnd(2), cUbnd(2)
      real(ESMF_KIND_R8), dimension(:,:), pointer :: ptr2d
      integer(ESMF_KIND_I4), dimension(:,:), pointer :: msk2d
      integer(ESMF_KIND_I4), allocatable, dimension(:,:) :: tlw, tuw
!
      type(ESMF_Grid) :: grid
      type(ESMF_DistGrid) :: distGrid
      type(ESMF_ArraySpec) :: arraySpec
      type(ESMF_StaggerLoc) :: staggerLoc      
!
      rc = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!     Query field
!-----------------------------------------------------------------------
!
      call ESMF_FieldGet(field, arrayspec=arraySpec,                    &
                         grid=grid, staggerloc=staggerLoc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Query grid
!-----------------------------------------------------------------------
!
      call ESMF_GridGet(grid, distgrid=distGrid,                        &
                        localDECount=localDECount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Allocate arrays for totalLWidth, totalUWidth and query field 
!-----------------------------------------------------------------------
!
      if (.not. allocated(tlw)) then
        allocate(tlw(2,localDECount))
        allocate(tuw(2,localDECount))
      end if
!
      call ESMF_FieldGet(field, totalLWidth=tlw, totalUWidth=tuw, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create field from base field attributes 
!-----------------------------------------------------------------------
!
      UTIL_FieldCreate = ESMF_FieldCreate(grid, arraySpec,              &
                                          staggerloc=staggerLoc,        &
                                          totalLWidth=tlw(:,1),         &
                                          totalUWidth=tuw(:,1),         &
                                          name=trim(fname), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      do k = 0, localDECount-1
!
!-----------------------------------------------------------------------
!     Get pointer from field 
!-----------------------------------------------------------------------
!
      call ESMF_FieldGet(UTIL_FieldCreate, localDe=k, farrayPtr=ptr2d,  &
                         computationalLBound=cLbnd,                     &
                         computationalUBound=cUbnd, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get pointer from grid (mask item) 
!-----------------------------------------------------------------------
!
      call ESMF_GridGetItem(grid, ESMF_GRIDITEM_MASK,                   &
                            staggerloc=staggerLoc,                      &
                            localDe=k, farrayPtr=msk2d, rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Initialize pointer 
!-----------------------------------------------------------------------
!
      do i = cLbnd(1), cUbnd(1)
        do j = cLbnd(2), cUbnd(2)
          if (msk2d(i,j) /= dstLandMask) then
            ptr2d(i,j) = initVal
          else
            ptr2d(i,j) = MISSING_R8
          end if
        end do
      end do
!
!-----------------------------------------------------------------------
!     Nullify pointer to make sure that it does not point on a random 
!     part in the memory 
!-----------------------------------------------------------------------
!
      if (associated(ptr2d)) then
        nullify(ptr2d)
      end if
      if (associated(msk2d)) then
        nullify(msk2d)
      end if
!
      end do
!
!-----------------------------------------------------------------------
!     Deallocate temporary fields
!-----------------------------------------------------------------------
!
      if (allocated(tlw)) then
        deallocate(tlw)
        deallocate(tuw)
      end if
!
!-----------------------------------------------------------------------
!     Check consistency of the created field 
!-----------------------------------------------------------------------
!
      call ESMF_FieldValidate(UTIL_FieldCreate, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!      
      end function UTIL_FieldCreate
!
      subroutine UTIL_FindUnmapped(srcField, dstField,                  &
                                   srcLandMask, dstLandMask,            &
                                   isForward, rc)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_Field), intent(in) :: srcField 
      type(ESMF_Field), intent(in) :: dstField
      integer, intent(in) :: srcLandMask 
      integer, intent(in) :: dstLandMask
      logical, intent(in) :: isForward 
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: cLbnd(2), cUbnd(2)
      integer :: i, j, k, p, localDECount
      character(ESMF_MAXSTR) :: fname
      real(ESMF_KIND_R8), dimension(:,:), pointer :: ptr2d, bdy2d
      integer(ESMF_KIND_I4), dimension(:,:), pointer :: msk2d
!
      type(ESMF_Grid) :: grid
      type(ESMF_Field) :: aField, bField, cField, mField
      type(ESMF_UnmappedAction_Flag) :: unmap
      type(ESMF_RegridMethod_Flag) :: regridMethod
      type(ESMF_RouteHandle) :: routeHandle
      type(ESMF_StaggerLoc) :: sLoc
!
      rc = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!     Create dummy fields 
!-----------------------------------------------------------------------
!
      fname = 'const_1'
      aField = UTIL_FieldCreate(srcField, fname, ONE_R8,                &
                                srcLandMask, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      fname = 'const_2'
      bField = UTIL_FieldCreate(dstField, fname, MISSING_R8,            &
                                dstLandMask, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      fname = 'const_3'
      cField = UTIL_FieldCreate(dstField, fname, ZERO_R8, -1, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create 1st routehandle 
!     Used to find the boundary of the destination grid
!-----------------------------------------------------------------------
!
      unmap = ESMF_UNMAPPEDACTION_IGNORE
      if (isForward) then
        regridMethod = ESMF_REGRIDMETHOD_NEAREST_STOD
      else
        regridMethod = ESMF_REGRIDMETHOD_NEAREST_DTOS
      end if
!
      call ESMF_FieldRegridStore(srcField=aField,                       &
                                 dstField=bField,                       &
                                 srcMaskValues=(/srcLandMask/),         &
                                 dstMaskValues=(/dstLandMask/),         &
                                 unmappedaction=unmap,                  &
                                 routeHandle=routeHandle,               &
                                 regridmethod=regridMethod,             &
                                 rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Perform regrid using 1st routehandle 
!-----------------------------------------------------------------------
!
      call ESMF_FieldRegrid(aField, bField, routeHandle,                &
                            zeroregion=ESMF_REGION_EMPTY, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create 2nd routehandle 
!     Used to find the unmapped grid cells
!-----------------------------------------------------------------------
!
      unmap = ESMF_UNMAPPEDACTION_IGNORE
      regridMethod = ESMF_REGRIDMETHOD_BILINEAR
!
      call ESMF_FieldRegridStore(srcField=aField,                       &
                                 dstField=cField,                       &
                                 srcMaskValues=(/srcLandMask/),         &
                                 dstMaskValues=(/dstLandMask/),         &
                                 unmappedaction=unmap,                  &
                                 routeHandle=routeHandle,               &
                                 regridmethod=regridMethod,             &
                                 rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Perform regrid using 2nd routehandle
!-----------------------------------------------------------------------
!
      call ESMF_FieldRegrid(aField, cField, routeHandle,                &
                            zeroregion=ESMF_REGION_TOTAL, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Query result field
!-----------------------------------------------------------------------
!
      call ESMF_FieldGet(cField, grid=grid, staggerloc=sLoc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of local DEs in the grid
!-----------------------------------------------------------------------
!
      call ESMF_GridGet(grid, localDECount=localDECount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      do k = 0, localDECount-1
!
!-----------------------------------------------------------------------
!     Get pointer from fields 
!-----------------------------------------------------------------------
!
      call ESMF_FieldGet(bField, localDe=k, farrayPtr=bdy2d, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_FieldGet(cField, localDe=k, farrayPtr=ptr2d,            &
                         computationalLBound=cLbnd,                     &
                         computationalUBound=cUbnd, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get pointer from grid (mask item) 
!-----------------------------------------------------------------------
!
      call ESMF_GridGetItem(grid, ESMF_GRIDITEM_MASK, staggerloc=sLoc,  &
                            localDe=k, farrayPtr=msk2d, rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Modify masking to split mapped and unmapped grid cells    
!-----------------------------------------------------------------------
!
      do i = cLbnd(1), cUbnd(1)
      do j = cLbnd(2), cUbnd(2)
        if ((bdy2d(i,j) < TOL_R8).and.(msk2d(i,j) /= dstLandMask)) then
          if (ptr2d(i,j) < ONE_R8/2.0d0) then
            msk2d(i,j) = UNMAPPED_MASK
          else
            msk2d(i,j) = MAPPED_MASK
          end if   
        end if
      end do
      end do
!
!-----------------------------------------------------------------------
!     Debug: Write mask to file    
!-----------------------------------------------------------------------
!
      mField = ESMF_FieldCreate(grid, msk2d, name='mask', rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      if (isForward) then
        call ESMF_FieldWrite(mField, 'mask_forward.nc', rc=rc)
      else
        call ESMF_FieldWrite(mField, 'mask_backward.nc', rc=rc)
      end if
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT) 
!
!-----------------------------------------------------------------------
!     Nullify pointer to make sure that it does not point on a random 
!     part in the memory 
!-----------------------------------------------------------------------
!
      if (associated(ptr2d)) then
        nullify(ptr2d)
      end if
      if (associated(bdy2d)) then
        nullify(bdy2d)
      end if
      if (associated(msk2d)) then
        nullify(msk2d)
      end if
!
!-----------------------------------------------------------------------
!     Remove temporary fields
!-----------------------------------------------------------------------
!
      call ESMF_FieldDestroy(aField, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_FieldDestroy(bField, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_FieldDestroy(cField, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      end do
! 
      end subroutine UTIL_FindUnmapped
      end module user_coupler
