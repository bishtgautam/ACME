
module openacc_mod
#if USE_OPENACC
  use openacc
  use kinds, only              : real_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize, ntrac, nc, nep, nelemd
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, cp
  use derivative_mod, only     : gradient, vorticity, gradient_wk, derivative_t, divergence, &
                                 gradient_sphere
  use element_mod, only        : element_t
  use fvm_control_volume_mod, only        : fvm_struct
  use spelt_mod, only          : spelt_struct
  use filter_mod, only         : filter_t, filter_P
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, smooth, TimeLevel_Qdp
  use prim_si_mod, only        : preq_pressure
  use diffusion_mod, only      : scalar_diffusion, diffusion_init
  use control_mod, only        : integration, test_case, filter_freq_advection,  hypervis_order, &
        statefreq, moisture, TRACERADV_TOTAL_DIVERGENCE, TRACERADV_UGRADQ, &
        prescribed_wind, nu_q, nu_p, limiter_option, hypervis_subcycle_q, rsplit
  use edge_mod, only           : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, initedgebuffer, edgevunpackmin, ghostbuffer3D_t
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use viscosity_mod, only      : biharmonic_wk_scalar, biharmonic_wk_scalar_minmax, neighbor_minmax
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only   : abortmp
  implicit none
  private
  save

  public :: euler_step_oacc
  public :: openacc_init

  type (EdgeBuffer_t) :: edgeAdv, edgeAdvQ3, edgeAdv_p1, edgeAdvQ2, edgeAdv1,  edgeveloc
  type (ghostBuffer3D_t)   :: ghostbuf_tr

  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  real(kind=real_kind), allocatable :: qmin            (:,:,:)
  real(kind=real_kind), allocatable :: qmax            (:,:,:)
  real(kind=real_kind), allocatable :: Vstar           (:,:,:,:,:)
  real(kind=real_kind), allocatable :: Qtens           (:,:,:,:,:)
  real(kind=real_kind), allocatable :: dp              (:,:,:,:)
  real(kind=real_kind), allocatable :: Qtens_biharmonic(:,:,:,:,:)
  real(kind=real_kind), allocatable :: new_dinv        (:,:,:,:,:)



contains



  subroutine openacc_init(elem)
    use dimensions_mod, only : nlev, qsize, nelemd
    type (element_t), intent(inout) :: elem(:)

    ! Shared buffer pointers.
    ! Using "=> null()" in a subroutine is usually bad, because it makes
    ! the variable have an implicit "save", and therefore shared between
    ! threads. But in this case we want shared pointers.
    real(kind=real_kind), pointer :: buf_ptr(:) => null()
    real(kind=real_kind), pointer :: receive_ptr(:) => null()
    integer :: i , j , ie

    ! this might be called with qsize=0
    ! allocate largest one first
    ! Currently this is never freed. If it was, only this first one should
    ! be freed, as only it knows the true size of the buffer.
    call initEdgeBuffer(edgeAdvQ3,max(nlev,qsize*nlev*3), buf_ptr, receive_ptr)  ! Qtens,Qmin, Qmax

    ! remaining edge buffers can share %buf and %receive with edgeAdvQ3
    ! (This is done through the optional 1D pointer arguments.)
    call initEdgeBuffer(edgeAdv1,nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(edgeAdv,qsize*nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(edgeAdv_p1,qsize*nlev + nlev,buf_ptr,receive_ptr) 
    call initEdgeBuffer(edgeAdvQ2,qsize*nlev*2,buf_ptr,receive_ptr)  ! Qtens,Qmin, Qmax

    ! Don't actually want these saved, if this is ever called twice.
    nullify(buf_ptr)
    nullify(receive_ptr)

    !$OMP BARRIER
    !$OMP MASTER

    ! this static array is shared by all threads, so dimension for all threads (nelemd), not nets:nete:
    allocate( qmin            (      nlev  ,qsize,nelemd) )
    allocate( qmax            (      nlev  ,qsize,nelemd) )
    allocate( Vstar           (np,np,nlev,2      ,nelemd) )
    allocate( Qtens           (np,np,nlev  ,qsize,nelemd) )
    allocate( dp              (np,np,nlev        ,nelemd) )
    allocate( Qtens_biharmonic(np,np,nlev  ,qsize,nelemd) )
    allocate( new_dinv        (np,np,2,2         ,nelemd) )
    do ie = 1 , nelemd
      do j = 1 , np
        do i = 1 , np
          new_dinv(i,j,:,:,ie) = elem(ie)%Dinv(:,:,i,j)
        enddo
      enddo
    enddo

    !$OMP END MASTER
    !$OMP BARRIER
  end subroutine openacc_init



  subroutine euler_step_oacc( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.  
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds          , only : real_kind
  use dimensions_mod , only : np, npdg, nlev
  use hybrid_mod     , only : hybrid_t
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t, gradient_sphere, vorticity_sphere
  use edge_mod       , only : edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchangev
  use hybvcoord_mod  , only : hvcoord_t
  implicit none
  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=real_kind), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier
  ! local
  real(kind=real_kind), dimension(np,np                       ) :: divdp, dpdiss
  real(kind=real_kind), pointer, dimension(:,:,:)               :: DSSvar
  real(kind=real_kind) :: vstar_tmp(np,np,nlev,2)
  real(kind=real_kind) :: gradQ  (np,np,nlev,2)
  real(kind=real_kind) :: dp_star(np,np,nlev  )
  real(kind=real_kind) :: dp0
  real(kind=real_kind) :: dptmp
  integer :: ie,q,i,j,k
  integer :: rhs_viss = 0
  integer(kind=8) :: tc1,tc2,tr,tm
  logical, save :: first_time = .true.
! call t_barrierf('sync_euler_step', hybrid%par%comm)
!   call t_startf('euler_step')
!pw++
  call t_startf('euler_step_openacc')
!pw--

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss = 0
  if ( limiter_option == 8  ) then
    ! when running lim8, we also need to limit the biharmonic, so that term needs
    ! to be included in each euler step.  three possible algorithms here:
    ! 1) most expensive:
    !     compute biharmonic (which also computes qmin/qmax) during all 3 stages
    !     be sure to set rhs_viss=1
    !     cost:  3 biharmonic steps with 3 DSS
    !
    ! 2) cheapest:
    !     compute biharmonic (which also computes qmin/qmax) only on first stage
    !     be sure to set rhs_viss=3
    !     reuse qmin/qmax for all following stages (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps with 1 DSS
    !     main concern:  viscosity 
    !     
    ! 3)  compromise:
    !     compute biharmonic (which also computes qmin/qmax) only on last stage
    !     be sure to set rhs_viss=3
    !     compute qmin/qmax directly on first stage
    !     reuse qmin/qmax for 2nd stage stage (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps, 2 DSS
    !
    !  NOTE  when nu_p=0 (no dissipation applied in dynamics to dp equation), we should
    !        apply dissipation to Q (not Qdp) to preserve Q=1
    !        i.e.  laplace(Qdp) ~  dp0 laplace(Q)                
    !        for nu_p=nu_q>0, we need to apply dissipation to Q * diffusion_dp
    !
    ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
    do ie = nets , nete
      ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q)
#endif
      do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
        dp(:,:,k,ie) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,k) 
        do q = 1 , qsize
          Qtens_biharmonic(:,:,k,q,ie) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp)/dp(:,:,k,ie)
        enddo
      enddo
    enddo

    ! compute element qmin/qmax
    if ( rhs_multiplier == 0 ) then
      do ie = nets , nete
        do k = 1 , nlev    
          do q = 1 , qsize
            qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
            qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
            qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
          enddo
        enddo
      enddo
      ! update qmin/qmax based on neighbor data for lim8
      call neighbor_minmax(elem,hybrid,edgeAdvQ2,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
    endif

    ! lets just reuse the old neighbor min/max, but update based on local data
    if ( rhs_multiplier == 1 ) then
      do ie = nets , nete
        do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
          do q = 1 , qsize
            qmin(k,q,ie)=min(qmin(k,q,ie),minval(Qtens_biharmonic(:,:,k,q,ie)))
            qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
            qmax(k,q,ie)=max(qmax(k,q,ie),maxval(Qtens_biharmonic(:,:,k,q,ie)))
          enddo
        enddo
      enddo
    endif

    ! get niew min/max values, and also compute biharmonic mixing term
    if ( rhs_multiplier == 2 ) then
      rhs_viss = 3
      ! compute element qmin/qmax  
      do ie = nets , nete
        do k = 1  ,nlev    
          do q = 1 , qsize
            qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
            qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
            qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
          enddo
        enddo
      enddo
      ! two scalings depending on nu_p:
      ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
      ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
      if ( nu_p > 0 ) then
        do ie = nets , nete
#if (defined ELEMENT_OPENMP)
          !$omp parallel do private(k, q, dp0, dpdiss)
#endif
          do k = 1 , nlev    
            dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
#if 0
            dpdiss(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%derived%psdiss_ave(:,:)
#else
            dpdiss(:,:) = elem(ie)%derived%dpdiss_ave(:,:,k)
#endif
            do q = 1 , qsize
              ! NOTE: divide by dp0 since we multiply by dp0 below
              Qtens_biharmonic(:,:,k,q,ie)=Qtens_biharmonic(:,:,k,q,ie)*dpdiss(:,:)/dp0
            enddo
          enddo
        enddo
      endif
      call biharmonic_wk_scalar_minmax( elem , qtens_biharmonic , deriv , edgeAdvQ3 , hybrid , nets , nete , qmin(:,:,nets:nete) , qmax(:,:,nets:nete) )
      do ie = nets , nete
#if (defined ELEMENT_OPENMP)
        !$omp parallel do private(k, q, dp0)
#endif
        do k = 1 , nlev    !  Loop inversion (AAM)
          dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
          do q = 1 , qsize
            ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
            qtens_biharmonic(:,:,k,q,ie) = -rhs_viss*dt*nu_q*dp0*Qtens_biharmonic(:,:,k,q,ie) / elem(ie)%spheremp(:,:)
          enddo
        enddo
      enddo
    endif
  endif  ! compute biharmonic mixing term and qmin/qmax



!$OMP BARRIER
if (hybrid%ithr == 0) then   !!!!!!!!!!!!!!!!!!!!!!!!! OMP MASTER !!!!!!!!!!!!!!!!!!!!!!!!!
  !$acc data present_or_create( nelemd,qsize,n0_qdp,elem,deriv,dt,qtens,hvcoord,new_dinv,rhs_multiplier )
if (first_time) then
  !$acc update          device( nelemd,qsize,n0_qdp,elem,deriv,dt,      hvcoord         ,rhs_multiplier )
  first_time = .false.
else
  !$acc update          device(              n0_qdp,elem,      dt                       ,rhs_multiplier )
endif
  !$acc wait

  if (hybrid%masterthread) call system_clock(tc1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$acc parallel loop gang vector collapse(4) private(dptmp) async(1)
  do ie = 1 , nelemd
    do k = 1 , nlev
      do j = 1 , np
        do i = 1 , np
          dptmp = elem(ie)%derived%dp(i,j,k) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(i,j,k) 
          Vstar(i,j,k,1,ie) = elem(ie)%derived%vn0(i,j,1,k) / dptmp
          Vstar(i,j,k,2,ie) = elem(ie)%derived%vn0(i,j,2,k) / dptmp
        enddo
      enddo
    enddo
  enddo

  !$acc parallel loop gang vector collapse(5) async(1)
  do ie = 1 , nelemd
    do q = 1 , qsize
      do k = 1 , nlev
        do j = 1 , np
          do i = 1 , np 
            Qtens(i,j,k,q,ie) = elem(ie)%state%Qdp(i,j,k,q,n0_qdp) - dt * divergence_sphere_2( Vstar(1,1,1,1,ie) , elem(ie)%state%Qdp(1,1,1,q,n0_qdp) , deriv , elem(ie) , i , j , k )
          enddo
        enddo
      enddo
    enddo
  enddo

! if ( limiter_option == 8 ) then
!   do ie = 1 , nelemd
!     do q = 1 , qsize
!       do k = 1 , nlev  ! Loop index added (AAM)
!         ! UN-DSS'ed dp at timelevel n0+1:  
!         dp_star(:,:,k) = dp(:,:,k,ie) - dt * elem(ie)%derived%divdp(:,:,k)  
!         if ( nu_p > 0 .and. rhs_viss /= 0 ) then
!           ! add contribution from UN-DSS'ed PS dissipation
!           dpdiss(:,:) = elem(ie)%derived%dpdiss_biharmonic(:,:,k)
!           dp_star(:,:,k) = dp_star(:,:,k) - rhs_viss * dt * nu_q * dpdiss(:,:) / elem(ie)%spheremp(:,:)
!         endif
!       enddo
!       ! apply limiter to Q = Qtens / dp_star 
!       call limiter_optim_iter_full( Qtens(:,:,:,q,ie) , elem(ie)%spheremp(:,:) , qmin(:,q,ie) , qmax(:,q,ie) , dp_star(:,:,:) )
!     enddo
!   enddo
! endif

  ! apply mass matrix, overwrite np1 with solution:
  ! dont do this earlier, since we allow np1_qdp == n0_qdp 
  ! and we dont want to overwrite n0_qdp until we are done using it
  !$acc parallel loop gang vector collapse(5) async(1)
  do ie = 1 , nelemd
    do q = 1 , qsize
      do k = 1 , nlev
        do j = 1 , np
          do i = 1 , np
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%spheremp(i,j) * Qtens(i,j,k,q,ie)
          enddo
        enddo
      enddo
    enddo
  enddo

  !$acc parallel loop gang vector collapse(3) async(1)
  do ie = 1 , nelemd
    do q = 1 , qsize
      do k = nlev , 1 , -1
        if ( limiter_option == 4 ) then
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          ! sign-preserving limiter, applied after mass matrix
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          call limiter2d_zero( elem(ie)%state%Qdp(:,:,k,q,np1_qdp) , hvcoord )
        endif
      enddo
    enddo
  enddo

  !$acc wait(1)
  if (hybrid%masterthread) call system_clock(tc2,tr)
  if (hybrid%masterthread) write(*,*) 'MYTIMER: ', dble(tc2-tc1)/tr

  !$acc update host( elem )
  !$acc wait
  !$acc end data
endif   !!!!!!!!!!!!!!!!!!!!!!!!! OMP END MASTER !!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP BARRIER



  do ie = nets , nete
    if ( DSSopt == DSSno_var ) then
      call edgeVpack(edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
    else
      call edgeVpack(edgeAdv_p1 , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
      if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
      if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
      if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
      do k = 1 , nlev   ! also DSS extra field
        DSSvar(:,:,k) = elem(ie)%spheremp(:,:) * DSSvar(:,:,k) 
      enddo
      call edgeVpack( edgeAdv_p1 , DSSvar(:,:,1:nlev) , nlev , nlev*qsize , elem(ie)%desc )
    endif
  enddo
!pw++
  call t_startf('eus_bexchV')
!pw--
  if ( DSSopt == DSSno_var ) then
    call bndry_exchangeV( hybrid , edgeAdv    )
  else
    call bndry_exchangeV( hybrid , edgeAdv_p1 )
  endif
!pw++
  call t_stopf('eus_bexchV')
!pw--
  do ie = nets , nete
    if ( DSSopt == DSSno_var ) then
      call edgeVunpack( edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
      do q = 1 , qsize
        do k = 1 , nlev    !  Potential loop inversion (AAM)
          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
        enddo
      enddo
    else
      call edgeVunpack( edgeAdv_p1 , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
      do q = 1 , qsize
        do k = 1 , nlev    !  Potential loop inversion (AAM)
          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
        enddo
      enddo
      if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
      if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
      if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
      call edgeVunpack( edgeAdv_p1 , DSSvar(:,:,1:nlev) , nlev , qsize*nlev , elem(ie)%desc )
      do k = 1 , nlev
        DSSvar(:,:,k) = DSSvar(:,:,k) * elem(ie)%rspheremp(:,:)
      enddo
    endif
  enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
!pw++
  call t_stopf('euler_step_openacc')
!pw--
!   call t_stopf('euler_step')
  end subroutine euler_step_oacc



  function divergence_sphere_2(vstar,qdp,deriv,elem,i,j,k) result(div)
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v
    implicit none
    real(kind=real_kind), intent(in)        :: vstar(np,np,nlev,2)
    real(kind=real_kind), intent(in)        :: qdp(np,np,nlev)
    integer             , intent(in), value :: i,j,k
    type (derivative_t)                     :: deriv
    type (element_t)                        :: elem
    real(kind=real_kind)                    :: div
    ! Local
    integer :: i, j, s
    real(kind=real_kind) ::  dudx00
    real(kind=real_kind) ::  dvdy00
    dudx00=0.0d0
    dvdy00=0.0d0
    do s=1,np
      dudx00 = dudx00 + deriv%Dvv(s,i)*( elem%metdet(s,j)*(elem%Dinv(1,1,s,j)*vstar(s,j,k,1) + elem%Dinv(1,2,s,j)*vstar(s,j,k,2))*qdp(s,j,k) )
      dvdy00 = dvdy00 + deriv%Dvv(s,j)*( elem%metdet(i,s)*(elem%Dinv(2,1,i,s)*vstar(i,s,k,1) + elem%Dinv(2,2,i,s)*vstar(i,s,k,2))*qdp(i,s,k) )
    enddo
    div=(dudx00+dvdy00)*(elem%rmetdet(i,j)*rrearth)
  end function divergence_sphere_2



  subroutine limiter2d_zero(Q,hvcoord)
  ! mass conserving zero limiter (2D only).  to be called just before DSS
  !
  ! this routine is called inside a DSS loop, and so Q had already
  ! been multiplied by the mass matrix.  Thus dont include the mass
  ! matrix when computing the mass = integral of Q over the element
  !
  ! ps is only used when advecting Q instead of Qdp
  ! so ps should be at one timelevel behind Q
  implicit none
  real (kind=real_kind), intent(inout) :: Q(np,np)
  type (hvcoord_t)     , intent(in   ) :: hvcoord
  ! local
  real (kind=real_kind) :: mass,mass_new
  integer i,j
  mass = 0
  !$acc loop reduction(+:mass) collapse(2) seq
  do j = 1 , np
    do i = 1 , np
      mass = mass + Q(i,j)
    enddo
  enddo

  ! negative mass.  so reduce all postive values to zero 
  ! then increase negative values as much as possible
  if ( mass < 0 ) then
    !$acc loop collapse(2) seq
    do j = 1 , np
      do i = 1 , np
        Q(i,j) = -Q(i,j) 
      enddo
    enddo
  endif
  mass_new = 0
  !$acc loop reduction(+:mass_new) collapse(2) seq
  do j = 1 , np
    do i = 1 , np
      if ( Q(i,j) < 0 ) then
        Q(i,j) = 0
      else
        mass_new = mass_new + Q(i,j)
      endif
    enddo
  enddo

  ! now scale the all positive values to restore mass
  if ( mass_new > 0 ) then
    !$acc loop collapse(2) seq
    do j = 1 , np
      do i = 1 , np
        Q(i,j) = Q(i,j) * abs(mass) / mass_new
      enddo
    enddo
  endif
  if ( mass     < 0 ) then
    !$acc loop collapse(2) seq
    do j = 1 , np
      do i = 1 , np
        Q(i,j) = -Q(i,j) 
      enddo
    enddo
  endif
  end subroutine limiter2d_zero



  subroutine limiter2d_zero_vertical(Q,hvcoord)
  ! mass conserving zero limiter (2D only).  to be called just before DSS
  !
  ! this routine is called inside a DSS loop, and so Q had already
  ! been multiplied by the mass matrix.  Thus dont include the mass
  ! matrix when computing the mass = integral of Q over the element
  !
  ! ps is only used when advecting Q instead of Qdp
  ! so ps should be at one timelevel behind Q
  implicit none
  real (kind=real_kind), intent(inout) :: Q(np,np,nlev)
  type (hvcoord_t)     , intent(in   ) :: hvcoord
  ! local
  real (kind=real_kind) :: dp(np,np)
  real (kind=real_kind) :: mass,mass_new
  integer i,j,k
  mass = 0
  !$acc loop vector reduction(+:mass) collapse(3)
  do k = nlev , 1 , -1
    do j = 1 , np
      do i = 1 , np
        mass = mass + Q(i,j,k)
      enddo
    enddo
  enddo

  ! negative mass.  so reduce all postive values to zero 
  ! then increase negative values as much as possible
  if ( mass < 0 ) then
    !$acc loop vector collapse(3)
    do k = nlev , 1 , -1
      do j = 1 , np
        do i = 1 , np
          Q(i,j,k) = -Q(i,j,k) 
        enddo
      enddo
    enddo
  endif

  mass_new = 0
  !$acc loop vector reduction(+:mass_new) collapse(3)
  do k = nlev , 1 , -1
    do j = 1 , np
      do i = 1 , np
        if ( Q(i,j,k) < 0 ) then
          Q(i,j,k) = 0
        else
          mass_new = mass_new + Q(i,j,k)
        endif
      enddo
    enddo
  enddo

  ! now scale the all positive values to restore mass
  if ( mass_new > 0 ) then
    !$acc loop vector collapse(3)
    do k = nlev , 1 , -1
      do j = 1 , np
        do i = 1 , np
          Q(i,j,k) = Q(i,j,k) * abs(mass) / mass_new
        enddo
      enddo
    enddo
  endif
  if ( mass     < 0 ) then
    !$acc loop vector collapse(3)
    do k = nlev , 1 , -1
      do j = 1 , np
        do i = 1 , np
          Q(i,j,k) = -Q(i,j,k) 
        enddo
      enddo
    enddo
  endif
  end subroutine limiter2d_zero_vertical



  subroutine edgeVpack(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    type (EdgeBuffer_t)                      :: edge
    integer,              intent(in)   :: vlyr
    real (kind=real_kind),intent(in)   :: v(np,np,vlyr)
    integer,              intent(in)   :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc
    ! Local variables
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,k,ir,ll
    integer :: is,ie,in,iw
    call t_startf('edge_pack')
    is = desc%putmapP(south)
    ie = desc%putmapP(east)
    in = desc%putmapP(north)
    iw = desc%putmapP(west)
    if(MODULO(np,2) == 0 .and. UseUnroll) then 
      do k=1,vlyr
        do i=1,np,2
          edge%buf(kptr+k,is+i  ) = v(i  ,1  ,k)
          edge%buf(kptr+k,is+i+1) = v(i+1,1  ,k)
          edge%buf(kptr+k,ie+i  ) = v(np ,i  ,k)
          edge%buf(kptr+k,ie+i+1) = v(np ,i+1,k)
          edge%buf(kptr+k,in+i  ) = v(i  ,np ,k)
          edge%buf(kptr+k,in+i+1) = v(i+1,np ,k)
          edge%buf(kptr+k,iw+i  ) = v(1  ,i  ,k)
          edge%buf(kptr+k,iw+i+1) = v(1  ,i+1,k)
        enddo
      enddo
    else
      do k=1,vlyr
        do i=1,np
          edge%buf(kptr+k,is+i)   = v(i  ,1 ,k)
          edge%buf(kptr+k,ie+i)   = v(np ,i ,k)
          edge%buf(kptr+k,in+i)   = v(i  ,np,k)
          edge%buf(kptr+k,iw+i)   = v(1  ,i ,k)
        enddo
      enddo
    endif
    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    if(desc%reverse(south)) then
      is = desc%putmapP(south)
      do k=1,vlyr
        do i=1,np
          ir = np-i+1
          edge%buf(kptr+k,is+ir)=v(i,1,k)
        enddo
      enddo
    endif
    if(desc%reverse(east)) then
      ie = desc%putmapP(east)
      do k=1,vlyr
        do i=1,np
          ir = np-i+1
          edge%buf(kptr+k,ie+ir)=v(np,i,k)
        enddo
      enddo
    endif
    if(desc%reverse(north)) then
      in = desc%putmapP(north)
      do k=1,vlyr
        do i=1,np
          ir = np-i+1
          edge%buf(kptr+k,in+ir)=v(i,np,k)
        enddo
      enddo
    endif
    if(desc%reverse(west)) then
      iw = desc%putmapP(west)
      do k=1,vlyr
        do i=1,np
          ir = np-i+1
          edge%buf(kptr+k,iw+ir)=v(1,i,k)
        enddo
      enddo
    endif
    ! SWEST
    do ll=swest,swest+max_corner_elem-1
      if (desc%putmapP(ll) /= -1) then
        do k=1,vlyr
          edge%buf(kptr+k,desc%putmapP(ll)+1)=v(1  ,1 ,k)
        enddo
      endif
    enddo
    ! SEAST
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if (desc%putmapP(ll) /= -1) then
        do k=1,vlyr
          edge%buf(kptr+k,desc%putmapP(ll)+1)=v(np ,1 ,k)
        enddo
      endif
    enddo
    ! NEAST
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if (desc%putmapP(ll) /= -1) then
        do k=1,vlyr
          edge%buf(kptr+k,desc%putmapP(ll)+1)=v(np ,np,k)
        enddo
      endif
    enddo
    ! NWEST
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if (desc%putmapP(ll) /= -1) then
        do k=1,vlyr
          edge%buf(kptr+k,desc%putmapP(ll)+1)=v(1  ,np,k)
        enddo
      endif
    enddo
    call t_stopf('edge_pack')
  end subroutine edgeVpack






#endif
end module openacc_mod




