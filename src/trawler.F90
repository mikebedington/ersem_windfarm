#include "fabm_driver.h"

module ersem_trawler

   use fabm_types
   use fabm_particle
   use fabm_expressions
   use fabm_builtin_models

   use ersem_shared

   implicit none

   private

   type,extends(type_particle_model),public :: type_ersem_trawler
      type (type_model_id),         allocatable,dimension(:) :: id_trawl_target,id_trawl_dest
      type (type_horizontal_dependency_id),    allocatable,dimension(:) :: id_trawl_targetc,id_trawl_targetn,id_trawl_targetp,id_trawl_targets
      type (type_bottom_state_variable_id), allocatable,dimension(:) :: id_trawl_destc,id_trawl_destn,id_trawl_destp,id_trawl_dests
      type (type_state_variable_id), allocatable,dimension(:) :: id_trawl_destpelc,id_trawl_destpeln,id_trawl_destpelp,id_trawl_destpels
      type (type_horizontal_dependency_id) :: id_sed_type,id_trawl           ! porosity (currently one single value is used across all layers)
      
!      type (type_bottom_state_variable_id) :: id_trawl

      type (type_horizontal_diagnostic_variable_id), allocatable,dimension(:) :: id_ftrawl_C
      type (type_horizontal_diagnostic_variable_id) :: id_effort
      integer  :: ntrawl_targets
     
      real(rk) :: trawl_hours, sand_code, mud_code, gravel_code
      real(rk),allocatable :: sand(:),mud(:),gravel(:)
      logical,allocatable  :: trawl_destispel(:)

   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_trawler),intent(inout),target :: self
      integer,                         intent(in)           :: configunit
      character(len=16) :: index
      integer           :: icount

      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are given in d-1.
      self%dt = 86400._rk
      !call self%register_state_variable(self%id_trawl,'trawling_activity','-','trawling factor',0._rk,minimum=0._rk)      
     
       call self%register_dependency(self%id_trawl,'trawling_activity','-','trawling factor')      
      
      ! Register parameters
      call self%get_parameter(self%ntrawl_targets,'ntrawl_targets','-','numbers of benthic variables affected by trawling')
      call self%get_parameter(self%trawl_hours,'trawling_hours','-','numbers of trawling hours per day',default=24._rk)
      call self%get_parameter(self%gravel_code,'gravel_code','-','code for gravel sediment type',default=0._rk)
      call self%get_parameter(self%sand_code,'sand_code','-','code for sandy sediment type',default=1._rk)
      call self%get_parameter(self%mud_code,'mud_code','-','code for muddy sediment type',default=2._rk)
      
      ! allocate and obtain impact of trawling
      allocate(self%sand(self%ntrawl_targets))
      allocate(self%mud(self%ntrawl_targets))
      allocate(self%gravel(self%ntrawl_targets))
      do icount= 1,self%ntrawl_targets
         write (index,'(i0)') icount
         call self%get_parameter(self%sand(icount),'sand'//trim(index),'-','relative impact of the trawler on target '//trim(index)//' on sand')
         call self%get_parameter(self%mud(icount),'mud'//trim(index),'-','relative impact of the trawler on target '//trim(index)//' on mud')
         call self%get_parameter(self%gravel(icount),'gravel'//trim(index),'-','relative impact of the trawler on target '//trim(index)//' on gravel')
      end do

      ! allocate and couple targets and destination modules
      allocate(self%id_trawl_target(self%ntrawl_targets))
      allocate(self%id_trawl_targetc(self%ntrawl_targets))
      allocate(self%id_trawl_targetn(self%ntrawl_targets))
      allocate(self%id_trawl_targetp(self%ntrawl_targets))
      allocate(self%id_trawl_targets(self%ntrawl_targets))
      allocate(self%id_trawl_dest(self%ntrawl_targets))
      allocate(self%id_trawl_destc(self%ntrawl_targets))
      allocate(self%id_trawl_destn(self%ntrawl_targets))
      allocate(self%id_trawl_destp(self%ntrawl_targets))
      allocate(self%id_trawl_dests(self%ntrawl_targets))
      allocate(self%id_trawl_destpelc(self%ntrawl_targets))
      allocate(self%id_trawl_destpeln(self%ntrawl_targets))
      allocate(self%id_trawl_destpelp(self%ntrawl_targets))
      allocate(self%id_trawl_destpels(self%ntrawl_targets))
      allocate(self%trawl_destispel(self%ntrawl_targets))
      
      allocate(self%id_ftrawl_C(self%ntrawl_targets))

      do icount= 1,self%ntrawl_targets
         write (index,'(i0)') icount
         call self%register_model_dependency(self%id_trawl_target(icount),'trawl_target'//trim(index))
         call self%register_dependency(self%id_trawl_targetn(icount),'trawl_target'//trim(index)//'n','mmol N/m^2','trawl target '//trim(index)//' nitrogen')
         call self%register_dependency(self%id_trawl_targetc(icount),'trawl_target'//trim(index)//'c','mmol C/m^2','trawl target '//trim(index)//' carbon') 
         call self%register_dependency(self%id_trawl_targetp(icount),'trawl_target'//trim(index)//'p','mmol P/m^2','trawl target '//trim(index)//' phosphorus')
         call self%register_dependency(self%id_trawl_targets(icount),'trawl_target'//trim(index)//'s','mmol Si/m^2','trawl target '//trim(index)//' silicate')
         call self%request_coupling_to_model(self%id_trawl_targetc(icount),self%id_trawl_target(icount),standard_variables%total_carbon)
         call self%request_coupling_to_model(self%id_trawl_targetn(icount),self%id_trawl_target(icount),standard_variables%total_nitrogen)
         call self%request_coupling_to_model(self%id_trawl_targetp(icount),self%id_trawl_target(icount),standard_variables%total_phosphorus)
         call self%request_coupling_to_model(self%id_trawl_targets(icount),self%id_trawl_target(icount),standard_variables%total_silicate)
         
         call self%register_model_dependency(self%id_trawl_dest(icount),'trawl_dest'//trim(index))
         call self%get_parameter(self%trawl_destispel(icount),'trawl_dest'//trim(index)//'ispel','','trawl_destination '//trim(index)//' is pelagic',default=.false.)
         if (self%trawl_destispel(icount)) then
            call self%register_state_dependency(self%id_trawl_destpelc(icount),'trawl_dest'//trim(index)//'c','mmol C/m^3','trawl_destination '//trim(index)//' carbon') 
            call self%register_state_dependency(self%id_trawl_destpeln(icount),'trawl_dest'//trim(index)//'n','mmol N/m^3','trawl_destination '//trim(index)//' nitrogen')
            call self%register_state_dependency(self%id_trawl_destpelp(icount),'trawl_dest'//trim(index)//'p','mmol P/m^3','trawl_destination '//trim(index)//' phosphorus')
            call self%register_state_dependency(self%id_trawl_destpels(icount),'trawl_dest'//trim(index)//'s','mmol Si/m^3','trawl_destination '//trim(index)//' silicate')
            call self%request_coupling_to_model(self%id_trawl_destpelc(icount),self%id_trawl_dest(icount),standard_variables%total_carbon)
            call self%request_coupling_to_model(self%id_trawl_destpeln(icount),self%id_trawl_dest(icount),standard_variables%total_nitrogen)
            call self%request_coupling_to_model(self%id_trawl_destpelp(icount),self%id_trawl_dest(icount),standard_variables%total_phosphorus)
            call self%request_coupling_to_model(self%id_trawl_destpels(icount),self%id_trawl_dest(icount),standard_variables%total_silicate)
         else
            call self%register_state_dependency(self%id_trawl_destc(icount),'trawl_dest'//trim(index)//'c','mmol C/m^2','trawl target '//trim(index)//' carbon') 
            call self%register_state_dependency(self%id_trawl_destn(icount),'trawl_dest'//trim(index)//'n','mmol N/m^2','trawl target '//trim(index)//' nitrogen')
            call self%register_state_dependency(self%id_trawl_destp(icount),'trawl_dest'//trim(index)//'p','mmol P/m^2','trawl target '//trim(index)//' phosphorus')
            call self%register_state_dependency(self%id_trawl_dests(icount),'trawl_dest'//trim(index)//'s','mmol Si/m^2','trawl target '//trim(index)//' silicate')
            call self%request_coupling_to_model(self%id_trawl_destc(icount),self%id_trawl_dest(icount),standard_variables%total_carbon)
            call self%request_coupling_to_model(self%id_trawl_destn(icount),self%id_trawl_dest(icount),standard_variables%total_nitrogen)
            call self%request_coupling_to_model(self%id_trawl_destp(icount),self%id_trawl_dest(icount),standard_variables%total_phosphorus)
            call self%request_coupling_to_model(self%id_trawl_dests(icount),self%id_trawl_dest(icount),standard_variables%total_silicate)
         end if
      call self%register_diagnostic_variable(self%id_ftrawl_C(icount),'ftrawl'//trim(index)//'C','mg C/m^2/d','trawling C',output=output_time_step_averaged,domain=domain_bottom,source=source_do_bottom)
      end do
      
      call self%register_dependency(self%id_sed_type,sediment_type)
      call self%register_diagnostic_variable(self%id_effort,'effort','-','effort',output=output_time_step_averaged,domain=domain_bottom,source=source_do_bottom) 
      end subroutine initialize
      
    subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)     
      
       class (type_ersem_trawler),intent(in) :: self           
       _DECLARE_ARGUMENTS_DO_BOTTOM_        

       integer  :: icount,istate
       real(rk) :: sedtype, rate, trawl
       real(rk) :: trawl_targetcP,trawl_targetnP,trawl_targetpP,trawl_targetsP, trawl_targetP
       real(rk) :: fluxc,fluxn,fluxp,fluxs

       character(len=16) :: index
       _HORIZONTAL_LOOP_BEGIN_
           
       _GET_HORIZONTAL_(self%id_sed_type,sedtype)
       _GET_HORIZONTAL_(self%id_trawl,trawl)
       _SET_HORIZONTAL_DIAGNOSTIC_ (self%id_effort,trawl)
       
       if (self%trawl_hours>0._rk) then
         ! Get target concentrations: benthic 
         do icount=1,self%ntrawl_targets
            _GET_HORIZONTAL_(self%id_trawl_targetc(icount),trawl_targetcP)
            _GET_HORIZONTAL_(self%id_trawl_targetn(icount),trawl_targetnP)
            _GET_HORIZONTAL_(self%id_trawl_targetp(icount),trawl_targetpP)
            _GET_HORIZONTAL_(self%id_trawl_targets(icount),trawl_targetsP)
            trawl_targetcP=trawl_targetcP*CMass

            ! first start by looking the type of sediment
            ! TODO: move as input file or in itialisaton to save the IF at evrytime step?
            if (sedtype==self%gravel_code) then
               rate=self%gravel(icount)
            elseif (sedtype==self%sand_code) then
               rate=self%sand(icount)
            elseif (sedtype==self%mud_code) then
               rate=self%mud(icount)
            else    ! everything with not classified type is considered not affected by trawling
               rate=0._rk
            endif

            
            ! neded to transform absolute reduction in rate based on first order kinetic
            ! calculated assuming that the "rate" reduction is achieved over "trawl"
            ! fraction of the grid cell in the "trawl_hours" interval
            ! reduction is capped to 95% reduction to avoid instability
            ! trawl_hours is cap between 1 and 24 to protect from instability and typos
            rate=log(1._rk+max(-0.95_rk,rate*trawl))*24._rk/min(max(1._rk,self%trawl_hours),24._rk)

            ! flux is positive leaving the source
            ! the sign is in the rate
            fluxc=rate*trawl_targetcP
            fluxn=rate*trawl_targetnP
            fluxp=rate*trawl_targetpP
            fluxs=rate*trawl_targetsP
            
            do istate=1,size(self%id_trawl_target(icount)%bottom_state)
               _GET_HORIZONTAL_(self%id_trawl_target(icount)%bottom_state(istate),trawl_targetP)
               _SET_BOTTOM_ODE_(self%id_trawl_target(icount)%bottom_state(istate),rate*trawl_targetP)
            end do            
            
            if (self%trawl_destispel(icount)) then
               if (_VARIABLE_REGISTERED_(self%id_trawl_destpelc(icount))) _SET_BOTTOM_EXCHANGE_(self%id_trawl_destpelc(icount),-fluxc/CMass)
               if (_VARIABLE_REGISTERED_(self%id_trawl_destpeln(icount))) _SET_BOTTOM_EXCHANGE_(self%id_trawl_destpeln(icount),-fluxn)
               if (_VARIABLE_REGISTERED_(self%id_trawl_destpelp(icount))) _SET_BOTTOM_EXCHANGE_(self%id_trawl_destpelp(icount),-fluxp)
               if (_VARIABLE_REGISTERED_(self%id_trawl_destpels(icount))) _SET_BOTTOM_EXCHANGE_(self%id_trawl_destpels(icount),-fluxs)
            else
               if (_VARIABLE_REGISTERED_(self%id_trawl_destc(icount)))  _SET_BOTTOM_ODE_(self%id_trawl_destc(icount),-fluxc/CMass)
               if (_VARIABLE_REGISTERED_(self%id_trawl_destn(icount)))  _SET_BOTTOM_ODE_(self%id_trawl_destn(icount),-fluxn)
               if (_VARIABLE_REGISTERED_(self%id_trawl_destp(icount)))  _SET_BOTTOM_ODE_(self%id_trawl_destp(icount),-fluxp)
               if (_VARIABLE_REGISTERED_(self%id_trawl_dests(icount)))  _SET_BOTTOM_ODE_(self%id_trawl_dests(icount),-fluxs)
            endif
            _SET_HORIZONTAL_DIAGNOSTIC_ (self%id_ftrawl_C(icount),-fluxc)
            
         end do
      end if

      _HORIZONTAL_LOOP_END_
    end subroutine do_bottom
 end module
