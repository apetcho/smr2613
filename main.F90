      program main
      use ESMF
!
      use user_model1 , only : userm1_register
      use user_model2 , only : userm2_register
      use user_coupler, only : usercpl_register

      implicit none
!
!-----------------------------------------------------------------------
!     Local variables 
!-----------------------------------------------------------------------
!
      integer :: localPet, petCount, rc
      character(len=ESMF_MAXSTR) :: gname1, gname2, cname
      character(len=ESMF_MAXSTR) :: forward_init='forward_init',        &
                                    forward_run='forward_run',          &
                                    backward_init='backward_init',      &
                                    backward_run='backward_run'
      integer(ESMF_KIND_I4) :: f_i=1, f_r=1, b_i=1, b_r=1
!
      type(ESMF_VM):: vm
      type(ESMF_GridComp) :: gcomp1, gcomp2
      type(ESMF_CplComp)  :: ccomp
      type(ESMF_State) :: c1imp, c1exp, c2imp, c2exp
!
!-----------------------------------------------------------------------
!     Initialize framework and get back default global VM
!-----------------------------------------------------------------------
!
      call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of PETs we are running with
!-----------------------------------------------------------------------
!
      call ESMF_VMGet(vm, petCount=petCount, localPet=localPet, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create the 2 model components and coupler
!-----------------------------------------------------------------------
!
      gname1 = "model1"
      gcomp1 = ESMF_GridCompCreate(name=gname1, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      gname2 = "model2"
      gcomp2 = ESMF_GridCompCreate(name=gname2, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      cname = "coupler"
      ccomp = ESMF_CplCompCreate(name=cname, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Register components 
!-----------------------------------------------------------------------
!
      call ESMF_GridCompSetServices(gcomp1, userm1_register, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompSetServices(gcomp2, userm2_register, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_CplCompSetServices(ccomp, usercpl_register, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Init section (component 1) 
!-----------------------------------------------------------------------
!
      c1imp = ESMF_StateCreate(name="comp1 import",                     &
                               stateintent=ESMF_STATEINTENT_IMPORT,     &
                               rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      c1exp = ESMF_StateCreate(name="comp1 export",                     &
                               stateintent=ESMF_STATEINTENT_EXPORT,     &
                               rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompInitialize(gcomp1, importState=c1imp,           &
                                   exportState=c1exp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Init section (component 2) 
!-----------------------------------------------------------------------
!
      c2imp = ESMF_StateCreate(name="comp2 import",                     &
                               stateintent=ESMF_STATEINTENT_IMPORT,     &
                               rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      c2exp = ESMF_StateCreate(name="comp2 export",                     &
                               stateintent=ESMF_STATEINTENT_EXPORT,     &
                               rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompInitialize(gcomp2, importState=c2imp,           &
                                   exportState=c2exp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Init section (coupler) 
!-----------------------------------------------------------------------
!      
      call ESMF_AttributeSet(c2imp, forward_init, f_i, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_CplCompInitialize(ccomp, importState=c1exp, exportState=c2imp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_AttributeSet(c1imp, backward_init, b_i, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_CplCompInitialize(ccomp, importState=c2exp, exportState=c1imp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Run section 
!-----------------------------------------------------------------------
!
      call ESMF_GridCompRun(gcomp1, importState=c1imp, exportState=c1exp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_AttributeSet(c2imp, forward_run, f_r, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_CplCompRun(ccomp, importState=c1exp, exportState=c2imp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompRun(gcomp2, importState=c2imp, exportState=c2exp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_AttributeSet(c1imp, backward_run, b_r, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_CplCompRun(ccomp, importState=c2exp, exportState=c1imp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompRun(gcomp1, importState=c1imp, exportState=c1exp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Finalize 
!-----------------------------------------------------------------------
!
      call ESMF_GridCompFinalize(gcomp1, importState=c1imp, exportState=c1exp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompFinalize(gcomp2, importState=c2imp, exportState=c2exp, rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_CplCompFinalize(ccomp, importState=c1exp, exportState=c2imp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Destroy 
!-----------------------------------------------------------------------
!
      call ESMF_StateDestroy(c1imp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_StateDestroy(c1exp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_StateDestroy(c2imp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_StateDestroy(c2exp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompDestroy(gcomp1, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_GridCompDestroy(gcomp2, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_CplCompDestroy(ccomp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
          line=__LINE__, file=__FILE__))                                &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Finalize 
!-----------------------------------------------------------------------
!
      call ESMF_Finalize(rc=rc)
      end program main
