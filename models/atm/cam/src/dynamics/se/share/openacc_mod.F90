#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module openacc_mod
#if USE_OPENACC
  use openacc
  use kinds, only              : real_kind, int_kind, log_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize, ntrac, nc, nep, nelemd, max_neigh_edges, max_corner_elem
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
  use edge_mod, only           : EdgeBuffer_t, edgerotate, edgevunpack, initedgebuffer, edgevunpackmin, ghostbuffer3D_t
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use viscosity_mod, only      : biharmonic_wk_scalar, biharmonic_wk_scalar_minmax, neighbor_minmax
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only       : abortmp
  use schedule_mod, only       : Cycle_t
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
  integer :: nbuf

  real(kind=real_kind)   , allocatable :: qmin            (:,:,:)
  real(kind=real_kind)   , allocatable :: qmax            (:,:,:)
  real(kind=real_kind)   , allocatable :: Vstar           (:,:,:,:,:)
  real(kind=real_kind)   , allocatable :: Qtens           (:,:,:,:,:)
  real(kind=real_kind)   , allocatable :: dp              (:,:,:,:)
  real(kind=real_kind)   , allocatable :: Qtens_biharmonic(:,:,:,:,:)
  real(kind=real_kind)   , allocatable :: new_dinv        (:,:,:,:,:)
  real(kind=real_kind)   , allocatable :: edgebuf         (:,:)
  real(kind=real_kind)   , allocatable :: edgerecv        (:,:)
  real(kind=real_kind)   , allocatable :: sendbuf         (:,:)
  real(kind=real_kind)   , allocatable :: recvbuf         (:,:)
  integer                , allocatable :: sendbuf_cycbeg  (:)
  integer                , allocatable :: recvbuf_cycbeg  (:)
  integer                , allocatable :: send_ptrP       (:)
  integer                , allocatable :: recv_ptrP       (:)
  integer                , allocatable :: send_lengthP    (:)
  integer                , allocatable :: recv_lengthP    (:)
  integer(kind=int_kind) , allocatable :: putmapP         (:,:)
  integer(kind=int_kind) , allocatable :: getmapP         (:,:)
  logical(kind=log_kind) , allocatable :: reverse         (:,:)



contains



  subroutine openacc_init(elem)
    use dimensions_mod, only : nlev, qsize, nelemd
    use schedule_mod  , only: schedule_t, cycle_t, schedule
    type (element_t), intent(inout) :: elem(:)
    real(kind=real_kind), pointer :: buf_ptr(:) => null()
    real(kind=real_kind), pointer :: receive_ptr(:) => null()
    integer :: i , j , ie
    integer :: nSendCycles, nRecvCycles, tot_send_len, tot_recv_len, icycle
    type(Schedule_t),pointer :: pSchedule

    call initEdgeBuffer( edgeAdvQ3  , max(nlev,qsize*nlev*3) , buf_ptr , receive_ptr )  ! Qtens , Qmin , Qmax
    call initEdgeBuffer( edgeAdv1   , nlev                   , buf_ptr , receive_ptr )
    call initEdgeBuffer( edgeAdv    , qsize*nlev             , buf_ptr , receive_ptr )
    call initEdgeBuffer( edgeAdv_p1 , qsize*nlev + nlev      , buf_ptr , receive_ptr )
    call initEdgeBuffer( edgeAdvQ2  , qsize*nlev*2           , buf_ptr , receive_ptr )  ! Qtens , Qmin , Qmax

    nullify(buf_ptr)
    nullify(receive_ptr)

    !$OMP BARRIER
    !$OMP MASTER

#ifdef _PREDICT
    pSchedule => Schedule(iam)
#else
    pSchedule => Schedule(1)
#endif
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles
    tot_send_len = 0
    do icycle = 1 , nSendCycles
      tot_send_len = tot_send_len + pSchedule%SendCycle(icycle)%lengthP
    enddo
    tot_recv_len = 0
    do icycle = 1 , nRecvCycles
      tot_recv_len = tot_recv_len + pSchedule%RecvCycle(icycle)%lengthP 
    enddo
    nbuf = 4*(np+max_corner_elem)*nelemd

    allocate( putmapP         (max_neigh_edges   ,nelemd) )
    allocate( getmapP         (max_neigh_edges   ,nelemd) )
    allocate( reverse         (max_neigh_edges   ,nelemd) )
    allocate( edgebuf         (qsize*nlev,nbuf          ) )
    allocate( edgerecv        (qsize*nlev,nbuf          ) )
    allocate( qmin            (      nlev  ,qsize,nelemd) )
    allocate( qmax            (      nlev  ,qsize,nelemd) )
    allocate( Vstar           (np,np,nlev,2      ,nelemd) )
    allocate( Qtens           (np,np,nlev  ,qsize,nelemd) )
    allocate( dp              (np,np,nlev        ,nelemd) )
    allocate( Qtens_biharmonic(np,np,nlev  ,qsize,nelemd) )
    allocate( new_dinv        (np,np,2,2         ,nelemd) )
    allocate( sendbuf         (nlev*qsize,tot_send_len)   )
    allocate( recvbuf         (nlev*qsize,tot_send_len)   )
    if (nSendCycles > 0) then
      allocate( send_ptrP     (nSendCycles) )
      allocate( send_lengthP  (nSendCycles) )
      allocate( sendbuf_cycbeg(nSendCycles) )
    endif
    if (nRecvCycles > 0) then
      allocate( recv_ptrP     (nRecvCycles) )
      allocate( recv_lengthP  (nRecvCycles) )
      allocate( recvbuf_cycbeg(nRecvCycles) )
    endif
    do ie = 1 , nelemd
      do j = 1 , np
        do i = 1 , np
          new_dinv(i,j,:,:,ie) = elem(ie)%Dinv(:,:,i,j)
        enddo
      enddo
      putmapP(:,ie) = elem(ie)%desc%putmapP(:)
      getmapP(:,ie) = elem(ie)%desc%getmapP(:)
      reverse(:,ie) = elem(ie)%desc%reverse(:)
    enddo

    do icycle = 1 , nSendCycles
      send_ptrP   (icycle) = pSchedule%SendCycle(icycle)%ptrP
      send_lengthP(icycle) = pSchedule%SendCycle(icycle)%lengthP
    enddo
    do icycle = 1 , nRecvCycles
      recv_ptrP   (icycle) = pSchedule%RecvCycle(icycle)%ptrP
      recv_lengthP(icycle) = pSchedule%RecvCycle(icycle)%lengthP
    enddo

    !For the packing of MPI data, determine where each cycle begins in the packed buffer
    sendbuf_cycbeg(1) = 1
    do icycle = 1 , nSendCycles-1
      sendbuf_cycbeg(icycle+1) = sendbuf_cycbeg(icycle) + pSchedule%SendCycle(icycle)%lengthP
    enddo
    recvbuf_cycbeg(1) = 1
    do icycle = 1 , nRecvCycles-1
      recvbuf_cycbeg(icycle+1) = recvbuf_cycbeg(icycle) + pSchedule%RecvCycle(icycle)%lengthP
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
  !$acc data  pcreate( nelemd,qsize,n0_qdp,elem,deriv,dt,qtens,hvcoord,new_dinv,rhs_multiplier,edgebuf,putmapP,getmapP,reverse,Vstar,sendbuf,recvbuf,send_lengthP,recv_lengthP,send_ptrP,recv_ptrP,sendbuf_cycbeg,recvbuf_cycbeg )
if (first_time) then
  !$acc update device( nelemd,qsize,n0_qdp,elem,deriv,dt,qtens,hvcoord,new_dinv,rhs_multiplier,edgebuf,putmapP,getmapP,reverse,Vstar,sendbuf,recvbuf,send_lengthP,recv_lengthP,send_ptrP,recv_ptrP,sendbuf_cycbeg,recvbuf_cycbeg )
  first_time = .false.
else
  !$acc update device(              n0_qdp,elem,      dt                       ,rhs_multiplier                                                                                                                                   )
endif
  !$acc wait
  call t_startf('euler_step_openacc')

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

  !$acc parallel loop gang vector collapse(5) async(1) vector_length(256)
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
  !$acc parallel loop gang vector collapse(5) async(1) vector_length(256)
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

  !$acc parallel loop gang vector collapse(3) async(1) vector_length(64)
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

  call edgeVpack_qdp  (elem,edgebuf,nlev*qsize,0,np1_qdp,putmapP,reverse)
  call bndry_exchangeV_openacc(hybrid,edgebuf,edgerecv,nlev*qsize)
  call edgeVunpack_qdp(elem,edgebuf,nlev*qsize,0,np1_qdp,getmapP)

  !$acc parallel loop gang vector collapse(5) async(1) vector_length(128)
  do ie = 1 , nelemd
    do q = 1 , qsize
      do k = 1 , nlev
        do j = 1 , np
          do i = 1 , np
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%rspheremp(i,j) * elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
          enddo
        enddo
      enddo
    enddo
  enddo

  !$acc wait(1)
  call t_stopf('euler_step_openacc')
  !$acc update host( elem )
  !$acc wait
  !$acc end data
endif   !!!!!!!!!!!!!!!!!!!!!!!!! OMP END MASTER !!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP BARRIER



  if ( DSSopt /= DSSno_var ) then
    do ie = nets , nete
      if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
      if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
      if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
      do k = 1 , nlev
        DSSvar(:,:,k) = elem(ie)%spheremp(:,:) * DSSvar(:,:,k) 
      enddo
      call edgeVpack( edgeAdv1 , DSSvar(:,:,1:nlev) , nlev , 0 , elem(ie)%desc )
    enddo
    call bndry_exchangeV( hybrid , edgeAdv1 )
    do ie = nets , nete
      if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
      if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
      if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
      call edgeVunpack( edgeAdv1 , DSSvar(:,:,1:nlev) , nlev , 0 , elem(ie)%desc )
      do k = 1 , nlev
        DSSvar(:,:,k) = DSSvar(:,:,k) * elem(ie)%rspheremp(:,:)
      enddo
    enddo
  endif



#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
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



  subroutine edgeVpack_qdp(elem,edgebuf,nlyr,kptr,nt,putmapP,reverse)
    use edge_mod, only: EdgeDescriptor_t, EdgeBuffer_t
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    implicit none
    type (element_t)       ,intent(in   ) :: elem(:)
    real(kind=real_kind)   ,intent(  out) :: edgebuf(nlyr,nbuf)
    integer                ,intent(in   ) :: nlyr
    integer                ,intent(in   ) :: kptr
    integer                ,intent(in   ) :: nt
    integer(kind=int_kind) ,intent(in   ) :: putmapP(max_neigh_edges,nelemd)
    logical(kind=log_kind) ,intent(in   ) :: reverse(max_neigh_edges,nelemd)
    ! Local variables
    integer :: i,k,ir,ll,kq,ie,q
    !$acc parallel loop gang vector collapse(4) private(kq) async(1) vector_length(64)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          do i = 1 , np
            kq = (q-1)*nlev+k
            edgebuf(kptr+kq,putmapP(south,ie)+i) = elem(ie)%state%Qdp(i ,1 ,k,q,nt)
            edgebuf(kptr+kq,putmapP(east ,ie)+i) = elem(ie)%state%Qdp(np,i ,k,q,nt)
            edgebuf(kptr+kq,putmapP(north,ie)+i) = elem(ie)%state%Qdp(i ,np,k,q,nt)
            edgebuf(kptr+kq,putmapP(west ,ie)+i) = elem(ie)%state%Qdp(1 ,i ,k,q,nt)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(4) private(kq,ir) async(1) vector_length(128)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          do i = 1 , np
            kq = (q-1)*nlev+k
            ir = np-i+1
            if(reverse(south,ie)) edgebuf(kptr+kq,putmapP(south,ie)+ir) = elem(ie)%state%Qdp(i ,1 ,k,q,nt)
            if(reverse(east ,ie)) edgebuf(kptr+kq,putmapP(east ,ie)+ir) = elem(ie)%state%Qdp(np,i ,k,q,nt)
            if(reverse(north,ie)) edgebuf(kptr+kq,putmapP(north,ie)+ir) = elem(ie)%state%Qdp(i ,np,k,q,nt)
            if(reverse(west ,ie)) edgebuf(kptr+kq,putmapP(west ,ie)+ir) = elem(ie)%state%Qdp(1 ,i ,k,q,nt)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(4) private(kq,ll) async(1) vector_length(64)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          do i = 1 , max_corner_elem
            kq = (q-1)*nlev+k
            ll = swest+0*max_corner_elem+i-1
            if (putmapP(ll,ie) /= -1) edgebuf(kptr+kq,putmapP(ll,ie)+1) = elem(ie)%state%Qdp(1 ,1 ,k,q,nt)    ! SWEST
            ll = swest+1*max_corner_elem+i-1
            if (putmapP(ll,ie) /= -1) edgebuf(kptr+kq,putmapP(ll,ie)+1) = elem(ie)%state%Qdp(np,1 ,k,q,nt)    ! SEAST
            ll = swest+2*max_corner_elem+i-1
            if (putmapP(ll,ie) /= -1) edgebuf(kptr+kq,putmapP(ll,ie)+1) = elem(ie)%state%Qdp(1 ,np,k,q,nt)    ! NWEST
            ll = swest+3*max_corner_elem+i-1
            if (putmapP(ll,ie) /= -1) edgebuf(kptr+kq,putmapP(ll,ie)+1) = elem(ie)%state%Qdp(np,np,k,q,nt)    ! NEAST
          enddo
        enddo
      enddo
    enddo
  end subroutine edgeVpack_qdp



  subroutine edgeVunpack_qdp(elem,edgebuf,vlyr,kptr,nt,getmapP)
    use edge_mod, only: EdgeDescriptor_t, EdgeBuffer_t
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    implicit none
    type(element_t)        ,intent(  out) :: elem(:)
    integer                ,intent(in   ) :: vlyr
    real(kind=real_kind)   ,intent(in   ) :: edgebuf(vlyr,nbuf)
    integer                ,intent(in   ) :: kptr
    integer                ,intent(in   ) :: nt
    integer(kind=int_kind) ,intent(in   ) :: getmapP(max_neigh_edges,nelemd)
    ! Local
    integer :: i,k,ll,q,ie,kq
    !$acc parallel loop gang vector collapse(4) private(kq) async(1) vector_length(64)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          do i = 1 , np
            kq = (q-1)*nlev+k
            elem(ie)%state%Qdp(i ,1 ,k,q,nt) = elem(ie)%state%Qdp(i ,1 ,k,q,nt) + edgebuf(kptr+kq,getmapP(south,ie)+i)
            elem(ie)%state%Qdp(i ,np,k,q,nt) = elem(ie)%state%Qdp(i ,np,k,q,nt) + edgebuf(kptr+kq,getmapP(north,ie)+i)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(4) private(kq) async(1) vector_length(128)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          do i = 1 , np
            kq = (q-1)*nlev+k
            elem(ie)%state%Qdp(1 ,i ,k,q,nt) = elem(ie)%state%Qdp(1 ,i ,k,q,nt) + edgebuf(kptr+kq,getmapP(west ,ie)+i)
            elem(ie)%state%Qdp(np,i ,k,q,nt) = elem(ie)%state%Qdp(np,i ,k,q,nt) + edgebuf(kptr+kq,getmapP(east ,ie)+i)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(4) private(kq,ll) async(1) vector_length(32)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          do i = 1 , max_corner_elem
            kq = (q-1)*nlev+k
            ll = swest+0*max_corner_elem+i-1
            if(getmapP(ll,ie) /= -1) elem(ie)%state%Qdp(1 ,1 ,k,q,nt) = elem(ie)%state%Qdp(1 ,1 ,k,q,nt) + edgebuf(kptr+kq,getmapP(ll,ie)+1)    ! SWEST
            ll = swest+1*max_corner_elem+i-1
            if(getmapP(ll,ie) /= -1) elem(ie)%state%Qdp(np,1 ,k,q,nt) = elem(ie)%state%Qdp(np,1 ,k,q,nt) + edgebuf(kptr+kq,getmapP(ll,ie)+1)    ! SEAST
            ll = swest+2*max_corner_elem+i-1
            if(getmapP(ll,ie) /= -1) elem(ie)%state%Qdp(1 ,np,k,q,nt) = elem(ie)%state%Qdp(1 ,np,k,q,nt) + edgebuf(kptr+kq,getmapP(ll,ie)+1)    ! NWEST
            ll = swest+3*max_corner_elem+i-1
            if(getmapP(ll,ie) /= -1) elem(ie)%state%Qdp(np,np,k,q,nt) = elem(ie)%state%Qdp(np,np,k,q,nt) + edgebuf(kptr+kq,getmapP(ll,ie)+1)    ! NEAST
          enddo
        enddo
      enddo
    enddo
  end subroutine edgeVunpack_qdp



  subroutine bndry_exchangeV_openacc(hybrid,edgebuf,edgerecv,nlyr)
    use openacc
    use hybrid_mod, only : hybrid_t
    use kinds, only : log_kind
    use edge_mod, only : Edgebuffer_t
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd, np
#ifdef _MPI
    use parallel_mod, only : abortmp, status, srequest, rrequest, mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : abortmp
#endif
    implicit none
    type (hybrid_t)       , intent(in   ) :: hybrid
    real(kind=real_kind)  , intent(inout) :: edgebuf (nlyr,nbuf)
    real(kind=real_kind)  , intent(inout) :: edgerecv(nlyr,nbuf)
    integer(kind=int_kind), intent(in   ) :: nlyr
    type (Schedule_t),pointer        :: pSchedule
    type (Cycle_t),pointer           :: pCycle
    integer                          :: dest,length,tag
    integer                          :: icycle,ierr
    integer                          :: iptr,source
    integer                          :: nSendCycles,nRecvCycles
    integer                          :: errorcode,errorlen
    character(len=80)                :: errorstring
    integer                          :: i, j, max_lengthP
    logical(kind=log_kind),parameter :: Debug = .FALSE.
#ifdef _MPI
    !Setup the pointer to proper Schedule
#ifdef _PREDICT
    pSchedule => Schedule(iam)
#else
    pSchedule => Schedule(1)
#endif
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    !==================================================
    !  Post the Receives 
    !==================================================
    do icycle = 1 , nRecvCycles
      pCycle => pSchedule%RecvCycle(icycle)
      source =  pCycle%source - 1
      length =  nlyr * pCycle%lengthP
      tag    =  pCycle%tag
      iptr   =  pCycle%ptrP
      call MPI_Irecv(recvbuf(1,recvbuf_cycbeg(icycle)),length,MPIreal_t,source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
      endif
    enddo

    !Pack the send arrays into one array for PCI-e transfer
    max_lengthP = maxval(send_lengthP)
    !$acc parallel loop gang vector collapse(3) async(1)
    do icycle = 1 , nSendCycles
      do j = 0 , max_lengthP
        do i = 1 , nlyr
          if ( j < send_lengthP(icycle) ) sendbuf( i , sendbuf_cycbeg(icycle)+j ) = edgebuf( i , send_ptrP(icycle)+j  )
        enddo
      enddo
    enddo
    !$acc update host(sendbuf) async(1)
    !$acc wait(1)

    !==================================================
    !  Fire off the sends
    !==================================================
    do icycle = 1 , nSendCycles
      pCycle => pSchedule%SendCycle(icycle)
      dest   =  pCycle%dest - 1
      length =  nlyr * pCycle%lengthP
      tag    =  pCycle%tag
      iptr   =  pCycle%ptrP
      call MPI_Isend(sendbuf(1,sendbuf_cycbeg(icycle)),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
      endif
    enddo
    !==================================================
    !  Wait for all the receives to complete
    !==================================================
    call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
    !$acc update device(recvbuf) async(1)

    !==================================================
    !  Wait for all the sends to complete
    !==================================================
    call MPI_Waitall(nSendCycles,Srequest,status,ierr)

    !Unpack the send arrays into one array for PCI-e transfer
    max_lengthP = maxval(recv_lengthP)
    !$acc parallel loop gang vector collapse(3) async(1)
    do icycle = 1 , nRecvCycles
      do j = 0 , max_lengthP
        do i = 1 , nlyr
          if (j < recv_lengthP(icycle) ) edgebuf( i , recv_ptrP(icycle)+j  ) = recvbuf( i , recvbuf_cycbeg(icycle)+j )
        enddo
      enddo
    enddo
#endif
  end subroutine bndry_exchangeV_openacc






#endif
end module openacc_mod




