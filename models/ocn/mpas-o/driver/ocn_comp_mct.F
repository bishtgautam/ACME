module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ocn_comp_mct
! !INTERFACE:

! !DESCRIPTION:
!  This is the main driver for the Model for Predication Across Scales Ocean Model (MPAS-O).
!
! !REVISION HISTORY:
!  SVN:$Id:
!
! !USES:
   use mct_mod
   use esmf
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use seq_comm_mct,      only : seq_comm_suffix, seq_comm_inst, seq_comm_name
   use shr_file_mod 
   use shr_cal_mod,       only : shr_cal_date2ymd
   use shr_sys_mod
   use shr_pio_mod
   use perf_mod
   use ocn_communicator,  only: mpi_communicator_ocn

   use mpaso_cpl_indices
   use mpaso_mct_vars

   use mpas_framework
   use mpas_core
   use mpas_kind_types
   use mpas_io_units
   use ocn_constants, only : T0_Kelvin

!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct
  SAVE
  private                              ! By default make data private

!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Doug Jacobsen
!
!EOP
! !PRIVATE MODULE FUNCTIONS:
  private :: ocn_export_mct
  private :: ocn_import_mct
  private :: ocn_SetGSMap_mct
  private :: ocn_domain_mct
!
! !PRIVATE MODULE VARIABLES

  logical :: ldiag_cpl = .false.

  integer, private ::   &
      cpl_write_restart,   &! flag id for write restart
      cpl_write_history,   &! flag id for write history
      cpl_write_tavg,      &! flag id for write tavg      
      cpl_diag_global,     &! flag id for computing diagnostics
      cpl_diag_transp,     &! flag id for computing diagnostics
      my_task

   real (kind=RKIND) :: tlast_coupled

   integer  :: nsend, nrecv

   character(len=StrKIND) :: runtype, coupleTimeStamp

   type(seq_infodata_type), pointer :: infodata   
   type (iosystem_desc_t), pointer :: io_system 

   !! MPAS-O Datatypes
   type (dm_info), pointer :: dminfo
   type (domain_type), pointer :: domain
   type (io_output_object), save :: output_obj
   integer :: output_frame, itimestep, ocn_cpl_dt

!=======================================================================

contains

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_init_mct
!
! !INTERFACE:
  subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )!{{{
!
! !DESCRIPTION:
! Initialize MPAS-O 
!

    use ocn_time_average_coupled
! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(inout) :: EClock
    type(seq_cdata), intent(inout) :: cdata_o
    type(mct_aVect), intent(inout) :: x2o_o, o2x_o
    character(len=*), optional, intent(in) :: NLFilename ! Namelist filename
!
! !REVISION HISTORY:
! Author: Doug Jacobsen
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer :: OCNID, mpicom_o, lsize, start_ymd, start_tod, start_year, start_day,   &
       start_month, start_hour, iyear, mpas_o_cpl_dt, shrloglev, shrlogunit, stdout_shr, pio_iotype

    type(mct_gsMap), pointer :: gsMap_o

    type(mct_gGrid), pointer :: dom_o

    integer :: errorCode  ! error code

    integer :: nThreads

    real (kind=RKIND) :: precadj

    character(len=StrKIND) :: cesm_calendar
    character(len=StrKIND) :: calendar_name
    integer :: calendar_id
    type (MPAS_Time_Type) :: currTime
    integer :: iam,ierr 
    integer :: iyear0, imonth0
    character(len=StrKIND)  :: starttype          ! infodata start type
    character(len=StrKIND)  :: timeStamp
    character(len=StrKIND)  :: nml_filename
    character(len=16) :: inst_suffix
    integer :: lbnum

    type (MPAS_Time_Type) :: alarmStartTime
    type (MPAS_TimeInterval_Type) :: alarmTimeStep
    type (block_type), pointer :: block_ptr
    logical :: exists, tight_coupling
!-----------------------------------------------------------------------
!
!  set cdata pointers
!
!-----------------------------------------------------------------------

   errorCode = 0 

   call seq_cdata_setptrs(cdata_o, ID=OCNID, mpicom=mpicom_o, &
        gsMap=gsMap_o, dom=dom_o, infodata=infodata)

   MPASO_MCT_OCNID =  OCNID
   MPASO_MCT_gsMap_o => gsMap_o
   MPASO_MCT_dom_o   => dom_o

   call seq_timemgr_EClockGetData(EClock, calendar=cesm_calendar)
   if( trim(cesm_calendar) == trim(seq_timemgr_noleap)) then
      write(calendar_name, *) 'gregorian_noleap'
   else if( trim(cesm_calendar) == trim(seq_timemgr_gregorian)) then
      write(calendar_name, *) 'gregorian'
   end if

#if (defined _MEMTRACE)
    call MPI_comm_rank(mpicom_o,iam,ierr)
    if(iam == 0) then
        lbnum=1
        call memmon_dump_fort('memmon.out','ocn_init_mct:start::',lbnum) 
    endif
#endif


    ! The following communicator module variable will be utilized in init_communicate that
    ! is called by initial - this is done to make the code backwards compatible

    mpi_communicator_ocn = mpicom_o

!-----------------------------------------------------------------------
!
!  initialize the model run 
!
!-----------------------------------------------------------------------

    call mpaso_cpl_indices_set()

    call seq_infodata_GetData( infodata, start_type=starttype, ocean_tight_coupling=tight_coupling)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       call mpas_dmpar_global_abort(' ocn_comp_mct ERROR: unknown starttype')
    end if
   
!-----------------------------------------------------------------------
!
!  first initializaiton phase of mpas-o
!  initialize mpas-o because grid information is needed for
!  creation of GSMap_ocn.
!  call mpas-o initialization routines
!
!-----------------------------------------------------------------------
!  inst_name   = seq_comm_name(OCNID)
!  inst_index  = seq_comm_inst(OCNID)
   inst_suffix = seq_comm_suffix(OCNID)

   call t_startf('mpas-o_init')
   nml_filename = "mpaso.in"

   io_system => shr_pio_getiosys(ocnid)

   pio_iotype = shr_pio_getiotype(ocnid)
   call MPAS_io_set_iotype(pio_iotype)

   stdout_shr = shr_file_getUnit()

   call mpas_framework_init(dminfo, domain, mpi_communicator_ocn, nml_filename, io_system, calendar_name, stdout_shr, stdout_shr)

   my_task = dminfo % my_proc_id

   if(my_task==0) then
      inquire(file='ocn_modelio.nml'//trim(inst_suffix),exist=exists)
      if(exists) then
         call shr_file_setio('ocn_modelio.nml'//trim(inst_suffix),stdout_shr)
      endif
   endif
   call shr_file_getLogUnit (shrlogunit)
   call shr_file_getLogLevel(shrloglev)
   call shr_file_setLogUnit (stdout_shr)


   call mpas_timer_start("total time")

   call mpas_core_setup_packages(ierr)
   call mpas_input_state_for_domain(domain)

   call mpas_core_init(domain, timeStamp)

!-----------------------------------------------------------------------
!
!  initialize time-stamp information
!
!-----------------------------------------------------------------------
   call t_stopf ('mpas-o_init')

!----------------------------------------------------------------------------
!
! reset shr logging to my log file
!
!----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout_shr)

   
!-----------------------------------------------------------------------
!
!  check for consistency of mpas-o and sync clock initial time
!
!-----------------------------------------------------------------------

  if (runtype == 'initial' .or. runtype == 'continue') then
     call seq_timemgr_EClockGetData(EClock, ECurrTime=currTime % t)
     call mpas_set_clock_time(clock, currTime, MPAS_START_TIME, ierr)
     call mpas_set_clock_time(clock, currTime, MPAS_NOW, ierr)

     currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
     call mpas_get_time(curr_time=currTime, YYYY=iyear0, MM=imonth0, ierr=ierr)
     call seq_timemgr_EClockGetData(EClock, start_ymd=start_ymd, start_tod=start_tod)
     call shr_cal_date2ymd(start_ymd,start_year,start_month,start_day)

     if (iyear0 /= start_year) then
        call mpas_dmpar_global_abort(' iyear0 does not match start_year')
     end if
     if (imonth0 /= start_month) then
        call mpas_dmpar_global_abort(' imonth0 does not match start_year')
     end if
  end if

!-----------------------------------------------------------------------
!
!  initialize MCT attribute vectors and indices
!
!-----------------------------------------------------------------------

   call t_startf ('mpas-o_mct_init')

   call ocn_SetGSMap_mct( mpicom_o, OCNID, GSMap_o )
   lsize = mct_gsMap_lsize(gsMap_o, mpicom_o)

   ! Initialize mct ocn domain (needs ocn initialization info)
   call ocn_domain_mct( lsize, gsMap_o, dom_o )
   
   ! Inialize mct attribute vectors
   
   call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
   call mct_aVect_zero(x2o_o)
   
   call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize) 
   call mct_aVect_zero(o2x_o)
   
   nsend = mct_avect_nRattr(o2x_o)
   nrecv = mct_avect_nRattr(x2o_o)

!-----------------------------------------------------------------------
!
!   initialize necessary  coupling info
!
!-----------------------------------------------------------------------

    call seq_timemgr_EClockGetData(EClock, dtime=ocn_cpl_dt)
    call convert_seconds_to_timestamp(ocn_cpl_dt, coupleTimeStamp)

    ! set coupling alarm
    call mpas_set_timeInterval(alarmTimeStep, timeString=coupleTimeStamp, ierr=ierr)
    alarmStartTime = currTime
    call mpas_add_clock_alarm(clock, coupleAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)

    ! Ocean clock starts one coupling interval after everything else.
    if (runtype == 'initial' .and. .not. tight_coupling) then
       call mpas_advance_clock(clock, alarmTimeStep)
    end if

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call ocn_time_average_coupled_init(block_ptr % forcing)
      call ocn_time_average_coupled_accumulate(block_ptr % diagnostics, block_ptr % forcing)
      block_ptr => block_ptr % next
    end do

!-----------------------------------------------------------------------
!
!  send intial state to driver
!
!-----------------------------------------------------------------------

   call ocn_export_mct(o2x_o, errorCode)  
   if (errorCode /= 0) then
      call mpas_dmpar_global_abort('ERROR in ocn_export_mct')
   endif

   call t_stopf ('mpas-o_mct_init')

   call seq_infodata_PutData( infodata, ocn_prognostic=.true., ocnrof_prognostic=.true.)

!-----------------------------------------------------------------------
!
!  get intial state from driver
!
!-----------------------------------------------------------------------

   call ocn_import_mct(x2o_o, errorCode)  
   if (errorCode /= 0) then
      call mpas_dmpar_global_abort('ERROR in ocn_import_mct')
   endif

   ! Build forcing arrays.
   block_ptr => domain % blocklist
   do while(associated(block_ptr))
       call ocn_forcing_build_arrays(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % forcing, ierr)
       call ocn_forcing_build_transmission_array(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % forcing, ierr)
       block_ptr => block_ptr % next
   end do


!----------------------------------------------------------------------------
!
! Reset shr logging to original values
!
!----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(iam  == 0) then
        lbnum=1
        call memmon_dump_fort('memmon.out','ocn_init_mct:end::',lbnum) 
        call memmon_reset_addr()
    endif
#endif

      !
      ! Set up output streams to be written to by the MPAS core
      !
      output_frame = 1

      if(config_frames_per_outfile > 0) then
         call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp))
      else
         call mpas_output_state_init(output_obj, domain, "OUTPUT")         
      end if

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      write(shrlogunit,*) 'Initial time ', trim(timeStamp)

      call ocn_write_output_frame(output_obj, output_frame, domain)
      itimestep = 0

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call ocn_time_average_init(block_ptr % average)
        block_ptr => block_ptr % next
      end do

      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)


!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_init_mct!}}}

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_run_mct
!
! !INTERFACE:
  subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)!{{{

    use mpas_kind_types
    use mpas_grid_types
    use mpas_io_output
    use mpas_timer
    use ocn_forcing
    use ocn_time_average_coupled

    implicit none
!
! !DESCRIPTION:
! Run MPAS-O for one coupling interval
!
! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o

!! !REVISION HISTORY:
!! Author: Doug Jacobsen
!!EOP
!!-----------------------------------------------------------------------
!!
!!  local variables
!!
!!-----------------------------------------------------------------------
      integer :: ymd, tod, ihour, iminute, isecond
      integer :: iyear, imonth, iday, curr_ymd, curr_tod
      integer :: shrloglev, shrlogunit
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime, alarmStartTime
      type (MPAS_TimeInterval_Type) :: alarmTimeStep
      character(len=StrKIND) :: timeStamp
      integer :: ierr

      ! Setup log information.
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (stdoutunit)

      dt = config_dt
      call mpas_reset_clock_alarm(clock, coupleAlarmID, ierr=ierr)

      ! Import state from coupler
      call ocn_import_mct(x2o_o, ierr)

      ! Initialize time average fields
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call ocn_time_average_coupled_init(block_ptr % forcing)
        block_ptr => block_ptr % next
      end do

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      ! This integration loop continues for a single coupling interval.
      do while (.not. mpas_is_alarm_ringing(clock,coupleAlarmID, ierr=ierr))
         itimestep = itimestep + 1
         call mpas_advance_clock(clock)

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         write(shrlogunit,*) 'Doing timestep ', trim(timeStamp)

         ! Build forcing arrays.
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
           call ocn_forcing_build_arrays(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % forcing, ierr)
           call ocn_forcing_build_transmission_array(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % forcing, ierr)
           block_ptr => block_ptr % next
         end do

         call mpas_timer_start("time integration", .false., timeIntTimer)
         call mpas_timestep(domain, itimestep, dt, timeStamp)
         call mpas_timer_stop("time integration", timeIntTimer)
   
         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
            call mpas_shift_time_levels_state(block_ptr % state)
            block_ptr => block_ptr % next
         end do
      
         ! Check if output file should be written
         if (mpas_is_alarm_ringing(clock, outputAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, outputAlarmID, ierr=ierr)
            ! output_frame will always be > 1 here unless it was reset after the maximum number of frames per outfile was reached
            if(output_frame == 1) then
               call mpas_output_state_finalize(output_obj, domain % dminfo)
               call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp))
            end if

            block_ptr => domain % blocklist
            do while (associated(block_ptr))
                call ocn_time_average_normalize(block_ptr % average)
                block_ptr => block_ptr % next
            end do

            call ocn_write_output_frame(output_obj, output_frame, domain)

            block_ptr => domain % blocklist
            do while (associated(block_ptr))
                call ocn_time_average_init(block_ptr % average)
                block_ptr => block_ptr % next
            end do
         end if
      end do

      ! Check if coupler wants us to write a restart file.
      if (seq_timemgr_RestartAlarmIsOn(EClock)) then
         call mpas_reset_clock_alarm(clock, restartAlarmID, ierr=ierr)

         ! Write one restart time per file
         call mpas_output_state_init(restart_obj, domain, "RESTART", trim(timeStamp))
         call mpas_output_state_for_domain(restart_obj, domain, 1)
         call mpas_output_state_finalize(restart_obj, domain % dminfo)
      end if

      ! Normalize time averaged fields
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call ocn_time_average_coupled_normalize(block_ptr % forcing)
        block_ptr => block_ptr % next
      end do

      ! Export state to coupler
      call ocn_export_mct(o2x_o, ierr)


      ! Check if clocks are in sync
      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, YYYY=iyear, MM=imonth, DD=iday, H=ihour, M=iminute, S=isecond, ierr=ierr)
      call seq_timemgr_EClockGetData(EClock, curr_ymd=curr_ymd, curr_tod=curr_tod)
 
      ymd = iyear * 10000 + imonth*100 + iday
      tod = ihour * 3600 + iminute * 60 + isecond
      if (.not. seq_timemgr_EClockDateInSync( EClock, ymd, tod)) then
         write(shrlogunit,*) 'MPAS ymd=',ymd,' MPAS tod=', tod
         write(shrlogunit,*) 'sync ymd=',curr_ymd,' sync tod=', curr_tod
         write(shrlogunit,*) 'Internal mpas clock not in sync with sync clock'
      end if

      ! Reset I/O logs
      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_run_mct!}}}

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_final_mct
!
! !INTERFACE:
  subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)!{{{
!
! !DESCRIPTION:
! Finalize MPAS-O
!
! !USES:
! !ARGUMENTS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
!
! !REVISION HISTORY:
! Author: Doug Jacobsen
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer :: shrloglev, shrlogunit
   integer :: &
      errorCode         ! error code

!-----------------------------------------------------------------------

    ! Setup I/O logs
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdoutunit)

    ! Finalize MPAS-O
    call mpas_output_state_finalize(output_obj, domain % dminfo)
    call mpas_timer_stop("total time")
    call mpas_timer_write()
    call mpas_core_finalize(domain)
    call mpas_framework_finalize(dminfo, domain, io_system)

    ! Reset I/O logs
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

  end subroutine ocn_final_mct!}}}

!***********************************************************************
!BOP
!IROUTINE: ocn_SetGSMap_mct
! !INTERFACE:

 subroutine ocn_SetGSMap_mct( mpicom_ocn, OCNID, gsMap_ocn )!{{{

   use mpas_dmpar

! !DESCRIPTION:
!  This routine sets up the MPAS-O grid numbering for MCT
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: mpicom_ocn
    integer        , intent(in)    :: OCNID
    type(mct_gsMap), intent(inout) :: gsMap_ocn

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer,allocatable :: &
      gindex(:)

    integer ::   &
      i,j, k, n, iblock, &
      lsize, gsize,   &
      ier

    type (block_type), pointer :: block_ptr 

    ! Loop over all cells in all blocks to determine total number.
    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      n = n + block_ptr % mesh % nCellsSolve
      block_ptr => block_ptr % next
    end do

    ! Determine total number of cells across all processors
    lsize = n
    call mpas_dmpar_sum_int(dminfo, lsize, gsize)
    allocate(gindex(lsize),stat=ier)

    ! Setup the mapping (gindex)
    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      do i = 1, block_ptr % mesh % nCellsSolve
        n = n + 1
        gindex(n) = block_ptr % mesh % indexToCellID % array(i)
      end do
      block_ptr => block_ptr % next
    end do

    ! Init the gsMap with gindex
    call mct_gsMap_init( gsMap_ocn, gindex, mpicom_ocn, OCNID, lsize, gsize )

    deallocate(gindex)

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_SetGSMap_mct!}}}

!***********************************************************************
!BOP
! !IROUTINE: ocn_domain_mct
! !INTERFACE:

 subroutine ocn_domain_mct( lsize, gsMap_o, dom_o )!{{{

! !DESCRIPTION:
!  This routine sets up the MCT domain for MPAS-O
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_o
    type(mct_ggrid), intent(inout) :: dom_o     

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer, pointer :: idata(:)

    real(kind=RKIND), pointer :: data(:)
    real(kind=RKIND) :: r2d

    integer :: i,j, k, n, ier

    type (block_type), pointer :: block_ptr

    r2d = 180.0/pii

!-------------------------------------------------------------------
!
!  initialize mct domain type
!
!-------------------------------------------------------------------

    call mct_gGrid_init( GGrid=dom_o, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_o%data)
    allocate(data(lsize))

!-------------------------------------------------------------------
!
! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
!
!-------------------------------------------------------------------

    call mct_gsMap_orderedPoints(gsMap_o, my_task, idata)
    call mct_gGrid_importIAttr(dom_o,'GlobGridNum',idata,lsize)

!-------------------------------------------------------------------
!
! Determine domain 
! Initialize attribute vector with special value
!
!-------------------------------------------------------------------

    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_o,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"aream",data,lsize) 
    data(:) = 1.0_R8     
    call mct_gGrid_importRAttr(dom_o,"mask",data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"frac",data,lsize) 

!-------------------------------------------------------------------
!
! Fill in correct values for domain components
!
!-------------------------------------------------------------------

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      do i = 1, block_ptr % mesh % nCellsSolve
        n = n + 1
        data(n) = block_ptr % mesh % lonCell % array(i) * r2d
      end do
      
      block_ptr => block_ptr % next
    end do
    call mct_gGrid_importRattr(dom_o,"lon",data,lsize) 

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      do i = 1, block_ptr % mesh % nCellsSolve
        n = n + 1
        data(n) = block_ptr % mesh % latCell % array(i) * r2d
      end do
      block_ptr => block_ptr % next
    end do
    call mct_gGrid_importRattr(dom_o,"lat",data,lsize) 

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      do i = 1, block_ptr % mesh % nCellsSolve
        n = n + 1
        data(n) = block_ptr % mesh % areaCell % array(i) / (block_ptr % mesh % sphere_radius * block_ptr % mesh % sphere_radius)
      end do
      block_ptr => block_ptr % next
    end do
    call mct_gGrid_importRattr(dom_o,"area",data,lsize) 

    data(:) = 1.0_RKIND ! No land cells in MPAS-O, only Ocean cells
    call mct_gGrid_importRattr(dom_o,"mask",data,lsize) 
    call mct_gGrid_importRattr(dom_o,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_domain_mct!}}}


!***********************************************************************
!BOP
! !IROUTINE: ocn_import_mct
! !INTERFACE:

 subroutine ocn_import_mct(x2o_o, errorCode)!{{{

    use ocn_constants

! !DESCRIPTION:
!-----------------------------------------------------------------------
!  This routine receives message from cpl7 driver
!
!    The following fields are always received from the coupler:
! 
!    o  taux   -- zonal wind stress (taux)                 (W/m2   )
!    o  tauy   -- meridonal wind stress (tauy)             (W/m2   )
!    o  snow   -- water flux due to snow                   (kg/m2/s)
!    o  rain   -- water flux due to rain                   (kg/m2/s)
!    o  evap   -- evaporation flux                         (kg/m2/s)
!    o  meltw  -- snow melt flux                           (kg/m2/s)
!    o  salt   -- salt                                     (kg(salt)/m2/s)
!    o  swnet  -- net short-wave heat flux                 (W/m2   )
!    o  sen    -- sensible heat flux                       (W/m2   )
!    o  lwup   -- longwave radiation (up)                  (W/m2   )
!    o  lwdn   -- longwave radiation (down)                (W/m2   )
!    o  melth  -- heat flux from snow&ice melt             (W/m2   )
!    o  ifrac  -- ice fraction                             (%)
!    o  rofl   -- river runoff flux                        (kg/m2/s)
!    o  rofi   -- ice runoff flux                          (kg/m2/s)
! 
!    The following fields are sometimes received from the coupler,
!      depending on model options:
! 
!    o  pslv   -- sea-level pressure                       (Pa)
!    o  duu10n -- 10m wind speed squared                   (m^2/s^2)
!    o  co2prog-- bottom atm level prognostic co2
!    o  co2diag-- bottom atm level diagnostic co2
! 
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    type(mct_aVect)   , intent(inout) :: x2o_o

! !OUTPUT PARAMETERS:

   integer, intent(out) :: &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (len=StrKIND) ::   &
      label,                 &
      message
 
   integer ::  &
      i,n

   type (block_type), pointer :: block_ptr

!-----------------------------------------------------------------------
!
!  zero out padded cells 
!
!-----------------------------------------------------------------------

   errorCode = 0

!-----------------------------------------------------------------------
!
!  unpack and distribute wind stress, then convert to correct units
!  and rotate components to local coordinates
!
!-----------------------------------------------------------------------

   n = 0
   block_ptr => domain % blocklist
   do while(associated(block_ptr))
      do i = 1, block_ptr % mesh % nCellsSolve
        n = n + 1
        block_ptr % forcing % windStressZonal % array(i) = x2o_o % rAttr(index_x2o_Foxx_taux, n)
        block_ptr % forcing % windStressMeridional % array(i) = x2o_o % rAttr(index_x2o_Foxx_tauy, n)

        block_ptr % forcing % latentHeatFlux % array(i) = x2o_o % rAttr(index_x2o_Foxx_lat, n)
        block_ptr % forcing % sensibleHeatFlux % array(i) = x2o_o % rAttr(index_x2o_Foxx_sen, n)
        block_ptr % forcing % longWaveHeatFluxUp % array(i) = x2o_o % rAttr(index_x2o_Foxx_lwup, n)
        block_ptr % forcing % longWaveHeatFluxDown % array(i) = x2o_o % rAttr(index_x2o_Faxa_lwdn, n)
        block_ptr % forcing % evaporationFlux % array(i) = x2o_o % rAttr(index_x2o_Foxx_evap, n)
        block_ptr % forcing % seaIceHeatFlux % array(i) = x2o_o % rAttr(index_x2o_Fioi_melth, n)
        block_ptr % forcing % snowFlux % array(i) = x2o_o % rAttr(index_x2o_Faxa_snow, n)

        block_ptr % forcing % seaIceFreshWaterFlux % array(i) = x2o_o % rAttr(index_x2o_Fioi_meltw, n)
        block_ptr % forcing % seaIceSalinityFlux % array(i) = x2o_o % rAttr(index_x2o_Fioi_salt, n)
        block_ptr % forcing % riverRunoffFlux % array(i) = x2o_o % rAttr(index_x2o_Foxx_rofl, n)
        block_ptr % forcing % iceRunoffFlux % array(i) = x2o_o % rAttr(index_x2o_Foxx_rofi, n)

        block_ptr % forcing % shortWaveHeatFlux % array(i) = max(x2o_o % rAttr(index_x2o_Foxx_swnet, n), 0.0_RKIND)

        block_ptr % forcing % rainFlux % array(i) = x2o_o % rAttr(index_x2o_Faxa_rain, n)
        block_ptr % forcing % seaSurfacePressure % array(i) = x2o_o % rAttr(index_x2o_Sa_pslv, n)
        block_ptr % forcing % iceFraction % array(i) = x2o_o % rAttr(index_x2o_Si_ifrac, n)

! mrp causes out-of-bounds index error in debug mode
!       block_ptr % forcing % prognosticCO2 % array(i) = x2o_o % rAttr(index_x2o_Sa_co2prog, n)
!       block_ptr % forcing % diagnosticCO2 % array(i) = x2o_o % rAttr(index_x2o_Sa_co2diag, n)
        block_ptr % forcing % squaredWindSpeed10Meter % array(i) = x2o_o % rAttr(index_x2o_So_duu10n, n)

        if(block_ptr % forcing % iceRunoffFlux % array(n) < 0.0_RKIND) then
            call shr_sys_abort ('Error: incoming rofi_F is negative')
        end if
      end do

      block_ptr => block_ptr % next
   end do

   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % windStressMeridional)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % windStressZonal)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % latentHeatFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % sensibleHeatFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % longWaveHeatFluxUp)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % longWaveHeatFluxDown)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % evaporationFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % seaIceHeatFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % snowFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % seaIceFreshWaterFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % seaIceSalinityFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % riverRunoffFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % iceRunoffFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % shortWaveHeatFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % rainFlux)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % seaSurfacePressure)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % iceFraction)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % prognosticCO2)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % diagnosticCO2)
   call mpas_dmpar_exch_halo_field(domain % blocklist % forcing % squaredWindSpeed10Meter)

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_import_mct!}}}
!***********************************************************************
!BOP
! !IROUTINE: ocn_export_mct
! !INTERFACE:

 subroutine ocn_export_mct(o2x_o, errorCode)   !{{{

! !DESCRIPTION:
!  This routine calls the routines necessary to send MPAS-O fields to
!  the CCSM cpl7 driver
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)   , intent(inout) :: o2x_o

! !OUTPUT PARAMETERS:

   integer, intent(out) :: &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer :: i,n, index_temperature, index_salinity, index_zonal_vel, index_meridional_vel, index_zonal_ssh, index_meridional_ssh

   type (block_type), pointer :: block_ptr

   n = 0
   block_ptr => domain % blocklist
   do while(associated(block_ptr))
     do i = 1, block_ptr % mesh % nCellsSolve
       n = n + 1

       index_temperature = block_ptr % forcing % index_temperatureSurfaceValue
       index_salinity = block_ptr % forcing % index_salinitySurfaceValue
       index_zonal_vel = block_ptr % forcing % index_avgZonalSurfaceVelocity
       index_meridional_vel = block_ptr % forcing % index_avgMeridionalSurfaceVelocity
       index_zonal_ssh = block_ptr % forcing % index_avgZonalSSHGradient
       index_meridional_ssh = block_ptr % forcing % index_avgMeridionalSSHGradient

       o2x_o % rAttr(index_o2x_So_t, n) = block_ptr % forcing % avgTracersSurfaceValue % array(index_temperature, i)
       o2x_o % rAttr(index_o2x_So_s, n) = block_ptr % forcing % avgTracersSurfaceValue % array(index_salinity, i)
       o2x_o % rAttr(index_o2x_So_u, n) = block_ptr % forcing % avgSurfaceVelocity % array(index_zonal_vel, i)
       o2x_o % rAttr(index_o2x_So_v, n) = block_ptr % forcing % avgSurfaceVelocity % array(index_meridional_vel, i)
       o2x_o % rAttr(index_o2x_So_dhdx, n) = block_ptr % forcing % avgSSHGradient % array(index_zonal_ssh, i)
       o2x_o % rAttr(index_o2x_So_dhdy, n) = block_ptr % forcing % avgSSHGradient % array(index_meridional_ssh, i)
       o2x_o % rAttr(index_o2x_Fioo_q, n) = block_ptr % forcing % seaIceEnergy % array(i) / ocn_cpl_dt

       ! Reset SeaIce Energy
       block_ptr % forcing % seaIceEnergy % array(i) = 0.0_RKIND

!      o2x_o % rAttr(index_o2x_Faoo_fco2_ocn, n) = block_ptr % forcing % CO2Flux % array(i)
!      o2x_o % rAttr(index_o2x_Faoo_fdms_ocn, n) = block_ptr % forcing % DMSFlux % array(i)
!      o2x_o % rAttr(index_o2x_Faoo_fco2_ocn, n) = block_ptr % forcing % surfaceUpwardCO2Flux % array(i)
     end do

     block_ptr => block_ptr % next
   end do


!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_export_mct!}}}

 subroutine convert_seconds_to_timestamp(seconds, timeStamp)!{{{
   integer, intent(in) :: seconds
   character (len=StrKIND), intent(out) :: timeStamp
   real (kind=RKIND) :: secondsPerHour, secondsPerMinute, remaining
   integer :: minutes, hours, secondsLeft

   secondsPerHour = 3600
   secondsPerMinute = 60

   if(seconds < 0 .or. seconds > 86400) then
     secondsLeft = 00
     minutes = 00
     hours = 00
   else
     hours = int(seconds/secondsPerHour)
     remaining = seconds - real(hours) * secondsPerHour

     minutes = int(remaining/secondsPerMinute)
     remaining = remaining - real(minutes) * secondsPerMinute

     secondsLeft = int(remaining)
   end if

   write(timeStamp,"(a,i2.2,a,i2.2,a,i2.2)") "00_",hours,":",minutes,":",secondsLeft
   timeStamp = trim(timeStamp)

 end subroutine convert_seconds_to_timestamp!}}}

end module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
