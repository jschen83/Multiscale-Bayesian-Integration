!****************************************************************************************
! Purpose:
!   Combine the prior and likelihoods from car and airborne surveys.
!
! Author:
!   Jinsong Chen and Haruko Wainwright
!   Lawrence Berkeley National Lab
!   Earth and Environmental Sciences Area
!   Berkeley, CA 94720
!
! Written: May of 2016
!
!----------------------------------------------------------------------------------------
! Copyright Notices:
!
! "Bayesian Integration of multiscale environmental data” Copyright (c) 2016, The Regents
! of the University of California, through Lawrence Berkeley National Laboratory (subject
! to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!   (1) Redistributions of source code must retain the above copyright notice, this list
!       of conditions and the following disclaimer.
!   (2) Redistributions in binary form must reproduce the above copyright notice, this
!       list of conditions and the following disclaimer in the documentation and/or other
!       materials provided with the distribution.
!   (3) Neither the name of the University of California, Lawrence Berkeley National 
!       Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used
!       to endorse or promote products derived from this software without specific prior
!       written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
! OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
! SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
! BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
! WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades
! to the features, functionality or performance of the source code ("Enhancements") to
! anyone; however, if you choose to make your Enhancements available either publicly, or
! directly to Lawrence Berkeley National Laboratory, without imposing a separate written
! license agreement for such Enhancements, then you hereby grant the following license:
! a non-exclusive, royalty-free perpetual license to install, use, modify, prepare 
! derivative works, incorporate into other computer software, distribute, and sublicense
! such enhancements or derivative works thereof, in binary and source code form.
!****************************************************************************************
!
module class_Post
use class_Global
use class_Subfun
use class_Prior
!
implicit none
!
  private
  public :: Post_Free,Post_InputPara,Post_CalcMeanCova
  !
  integer :: nx,ny   !Domain size 
  integer :: ngrid   !# of grid points 
  integer :: nobs    !# of walking survey points 
  integer :: ncar    !# of car survey data points 
  integer :: nair    !# of air survey data points 
  integer :: nsmooth !# of smoothing points 
  real    :: dx,dy   !discretization
  real    :: gmu     !overall mean values
  !
  type(geoModel)                       :: gstat(3)
  real, dimension(:,:), allocatable    :: smooth
  integer, dimension(:,:), allocatable :: lusemat
  !
contains
  !
  subroutine Post_Free
  !
  if(allocated(smooth)) deallocate(smooth)
  if(allocated(lusemat)) deallocate(lusemat)
  call Prior_Free
  end subroutine Post_Free

  subroutine Post_InputPara(grid,gmod,insmooth,inluse,walkdata)
  type(gridDomain), intent(in) :: grid     !Hold general grid information
  type(geoModel), intent(in)   :: gmod(3)  !Variogram model information
  type(surveyData), intent(in) :: walkdata !Walking survey data
  real, intent(in)             :: insmooth(:,:)
  integer, intent(in)          :: inluse(:,:)
  !
  ! Input domain, variogram models, and walking survey data
  !
  integer            :: ierr(6)
  integer, parameter :: locDebug=0
  !
  nx=grid%nx
  ny=grid%ny
  nsmooth=grid%nsmooth
  dx=grid%dx
  dy=grid%dy
  gmu=grid%gmean
  ngrid=nx*ny
  nobs=walkdata%num
  !
  ! Allocate memory
  !
  allocate(smooth(nsmooth,2),STAT=ierr(1))
  allocate(lusemat(ngrid,3),STAT=ierr(2))
  if(any(ierr(1:2)>0)) then
    write(ilog,*) 'Post_InputPara: Errors in allocating.'
    stop
  endif
  smooth(1:nsmooth,1:2)=insmooth(1:nsmooth,1:2)
  lusemat(1:ngrid,1:3)=inluse(1:ngrid,1:3)
  !
  ! Initilize Prior module
  !
  call Prior_InputPara(grid,gmod,walkdata)
  end subroutine Post_InputPara

  subroutine Post_CalcMeanCova(cardata,airdata,rtnmean,rtnstd)
  type(surveyData), intent(in) :: cardata       !Car survey data
  type(surveyData), intent(in) :: airdata       !Airborne survey data
  real, intent(inout)          :: rtnmean(:,:)  !Prior and posterior mean
  real, intent(inout)          :: rtnstd(:,:)   !Prior and posterior STD
  !
  ! Compute posterior distribution of dose given car and airborne data
  !
  integer      :: i,j,ierr(10)
  real         :: start_time,end_time                 
  real(kind=q) :: alpha,beta                 
  !
  real(kind=q), dimension(:), allocatable   :: zadose,zcdose
  real(kind=q), dimension(:), allocatable   :: tmpmu,primu,primean
  real(kind=q), dimension(:,:), allocatable :: pricova,invpricova
  real(kind=q), dimension(:,:), allocatable :: Amat,Cmat,invDa,invDc
  real(kind=q), dimension(:,:), allocatable :: Qmat,tmpMat1,tmpMat2
  real(kind=q), dimension(:,:), allocatable :: tAmat,tCmat
  integer, parameter                        :: locDebug=0
  !
  call cpu_time(start_time)
  !
  ncar=cardata%num
  nair=airdata%num
  if(locDebug==1) then
    write(idbg,*) 'ncar=',ncar
    write(idbg,*) 'nair=',nair
    write(*,*) 'ncar=',ncar
    write(*,*) 'nair=',nair
  endif
  !
  ! Prepare Amat and Cmat matrices
  !
  allocate(Amat(nair,ngrid),STAT=ierr(1))
  allocate(Cmat(ncar,ngrid),STAT=ierr(2))
  if(any(ierr(1:2)>0)) then
    write(ilog,*) 'Post_CalcMeanCova: Errors in allocating-1.'
    stop
  endif
  call Assemble_Amat(airdata,Amat)
  call Assemble_Cmat(cardata,Cmat)
  !
  ! Prepare airborne and car survey data using normal scores
  !
  allocate(invDa(nair,nair),STAT=ierr(1))
  allocate(invDc(ncar,ncar),STAT=ierr(2))
  allocate(zadose(nair),STAT=ierr(3))
  allocate(zcdose(ncar),STAT=ierr(4))
  if(any(ierr(1:4)>0)) then
    write(ilog,*) 'Post_CalcMeanCova: Errors in allocating-2.'
    stop
  endif
  call Prepare_Zadose_invDa(airdata,zadose,invDa)
  call Prepare_Zcdose_invDc(cardata,zcdose,invDc)
  !
  ! Compute t(Amat)*invDa and save the result to tmpMat1
  !
  allocate(tAmat(ngrid,nair),STAT=ierr(1))
  allocate(tmpMat1(ngrid,nair),STAT=ierr(2))
  if(any(ierr(1:2)>0)) then
    write(ilog,*) 'Post_CalcMeanCova: Errors in allocating-3.'
    stop
  endif
  tAmat=transpose(Amat) !Need for a large size of matrix
  call Subfun_mklMatmul(ngrid,nair,nair,tAmat,invDa,tmpMat1)
  if(allocated(tAmat)) deallocate(tAmat)
  if(allocated(invDa)) deallocate(invDa)
  !
  ! Compute t(Cmat)*invDc and save the result to tmpMat2
  !
  allocate(tCmat(ngrid,ncar),STAT=ierr(1))
  allocate(tmpMat2(ngrid,ncar),STAT=ierr(2))
  if(any(ierr(1:2)>0)) then
    write(ilog,*) 'Post_CalcMeanCova: Errors in allocating-4.'
    stop
  endif
  tCmat=transpose(Cmat)  !Need for a large size of matrix
  call Subfun_mklMatmul(ngrid,ncar,ncar,tCmat,invDc,tmpMat2)
  if(allocated(tCmat)) deallocate(tCmat)
  if(allocated(invDc)) deallocate(invDc)
  call cpu_time(end_time)
  write(*,*) 'CPU time for preparing matrices:', (end_time-start_time) 
  !
  ! Call Prior module to calculate prior mean and covariance matrix
  !
  call cpu_time(start_time)
  allocate(primean(ngrid),STAT=ierr(1))
  allocate(pricova(ngrid,ngrid),STAT=ierr(2))
  if(any(ierr(1:2)>0)) then
    write(ilog,*) 'Post_CalcMeanCova: Errors in allocating-5.'
    stop
  endif
  call Prior_CalcMeanCova(lusemat,primean,pricova)
  !
  ! Save prior mean and std for output
  !
  rtnmean(1:ngrid,1)=real(primean(1:ngrid))
  do i=1,ngrid
     rtnstd(i,1)=real(sqrt(pricova(i,i)))
  enddo
  call cpu_time(end_time)
  write(*,*) 'CPU time for calculating prior:', (end_time-start_time) 
  !
  ! Compute posterior covariance by conditioning to airborne and car data
  !
  call cpu_time(start_time)
  allocate(tmpmu(ngrid),STAT=ierr(1))
  allocate(invpricova(ngrid,ngrid),STAT=ierr(2))
  allocate(Qmat(ngrid,ngrid),STAT=ierr(3))
  if(any(ierr(1:3)>0)) then
    write(ilog,*) 'Post_CalcMeanCova: Errors in allocating-6.'
    stop
  endif
  invpricova(1:ngrid,1:ngrid)=pricova(1:ngrid,1:ngrid)
  if(allocated(pricova)) deallocate(pricova)
  !
  ! Invert prior covariance 
  !
  call Subfun_invert_pd_matrix_lapack2(invpricova,ngrid,ierr(1))
  if(ierr(1)>0) then
    write(ilog,*) 'Post_CalcMeanCova: invert prior covariance matrix.'
    stop
  endif
  !
  ! Compute intermediate combined mean 
  !
  tmpmu(1:ngrid)=matmul(invpricova,primean)+matmul(tmpMat1,zadose)+matmul(tmpMat2,zcdose)
  if(allocated(primean)) deallocate(primean)
  if(allocated(zadose)) deallocate(zadose)
  if(allocated(zcdose)) deallocate(zcdose)
  !
  ! Compute intermediate combined covariance matrix 
  !
  alpha=real(1.0,kind=q)
  beta=real(1.0,kind=q)
  Qmat(1:ngrid,1:ngrid)=invpricova(1:ngrid,1:ngrid)
  call Subfun_mklMatmul(ngrid,nair,ngrid,alpha,beta,tmpMat1,Amat,Qmat)
  call Subfun_mklMatmul(ngrid,ncar,ngrid,alpha,beta,tmpMat2,Cmat,Qmat)
  if(allocated(Amat)) deallocate(Amat)
  if(allocated(Cmat)) deallocate(Cmat)
  if(allocated(tmpMat1)) deallocate(tmpMat1)
  if(allocated(tmpMat2)) deallocate(tmpMat2)
  if(allocated(invpricova)) deallocate(invpricova)
  !
  ! Compute posterior mean and covariance by first inverting Qmat
  !
  call cpu_time(start_time)
  call Subfun_invert_pd_matrix_lapack2(Qmat,ngrid,ierr(1))
  call cpu_time(end_time)
  write(*,*) 'CPU time for inverting Qmat:', (end_time-start_time) 
  !
  rtnmean(1:ngrid,2)=real(matmul(Qmat,tmpmu))  !Qmat is inverted
  do i=1,ngrid
     rtnstd(i,2)=real(sqrt(Qmat(i,i)))
  enddo
  if(allocated(tmpmu)) deallocate(tmpmu)
  if(allocated(Qmat)) deallocate(Qmat)
  !
  if(locDebug==1) then
    open(22,file='post_mean_std.out',action='Write',iostat=ierr(1))
    do i=1,ngrid
       write(22,*) rtnmean(i,2),rtnstd(i,2)
    enddo
    close(22)
  endif
  end subroutine Post_CalcMeanCova

  subroutine Prepare_Zadose_invDa(airdata,zadose,invDa)
  type(surveyData), intent(in) :: airdata 
  real(kind=q), intent(inout)  :: zadose(:)
  real(kind=q), intent(inout)  :: invDa(:,:)
  !
  ! Convert airborne data to normal scores and calculate inverse error matrix invDa
  !
  integer :: i,k
  !
  invDa(1:nair,1:nair)=real(0.0,kind=q)
  do i=1,nair
     k=airdata%landuse(i)
     zadose(i)=real((airdata%datamat(i,3)-airdata%muvec(k))/airdata%stdvec(k),kind=q)
     invDa(i,i)=real(1.0/airdata%varvec(k),kind=q)
  enddo
  end subroutine Prepare_Zadose_invDa

  subroutine Prepare_Zcdose_invDc(cardata,zcdose,invDc)
  type(surveyData), intent(in) :: cardata 
  real(kind=q), intent(inout)  :: zcdose(:)
  real(kind=q), intent(inout)  :: invDc(:,:)
  !
  ! Convert car survey data to normal scores and calculate inverse error matrix invDc
  !
  integer :: i,k
  !
  invDc(1:ncar,1:ncar)=real(0.0,kind=q)
  do i=1,ncar
     k=cardata%landuse(i)
     zcdose(i)=real((cardata%datamat(i,3)-cardata%muvec(k))/cardata%stdvec(k),kind=q)
     invDc(i,i)=real(1.0/cardata%varvec(k),kind=q)
  enddo
  end subroutine Prepare_Zcdose_invDc

  subroutine Assemble_Amat(airdata,Amat)
  type(surveyData), intent(in) :: airdata 
  real(kind=q), intent(inout)  :: Amat(:,:)
  !
  ! Construct smoothing matrix-A from airborne survey data and weight matrix.
  !
  integer      :: i,k,ix,iy,mycount
  real(kind=q) :: x0,y0,x1,y1,tmpdx,tmpdy,tmpdist,tmpval
  !
  nair=airdata%num
  ! 
  Amat(1:nair,1:ngrid)=real(0.0,kind=q)
  do i=1,nair
     x0=real(airdata%datamat(i,1),kind=q)
     y0=real(airdata%datamat(i,2),kind=q)
     do k=1,ngrid
        call grid_location(k,ix,iy)
        tmpdx=real(x0-(ix-1)*dx,kind=q)
        tmpdy=real(y0-(iy-1)*dy,kind=q)
        tmpdist=sqrt(tmpdx*tmpdx+tmpdy*tmpdy)
        if(tmpdist<smooth(nsmooth,1)) then !Make sure within the largest distance
          tmpval=Subfun_linear_interp1d(nsmooth,real(smooth(1:nsmooth,1),kind=q), &
                                        real(smooth(1:nsmooth,2),kind=q),tmpdist)
        else
          tmpval=real(0.0,kind=q)
        endif
        Amat(i,k)=tmpval
      enddo
      mycount=count(Amat(i,1:ngrid).gt.0)
      if(mycount==0) then
        write(ilog,*) 'Assemble_Amatrix: Errors in constructing Amat'
        stop
      endif
      Amat(i,1:ngrid)=Amat(i,1:ngrid)/sum(Amat(i,1:ngrid))
  enddo
  end subroutine Assemble_Amat

  subroutine Assemble_Cmat(cardata,Cmat)
  type(surveyData), intent(in) :: cardata 
  real(kind=q), intent(inout)  :: Cmat(:,:)
  !
  real(kind=q) :: x0,y0,x1,y1,mydist,mysum,myradius
  integer      :: i,k,ii,jj
  !
  myradius=real(100.0,kind=q) !Effective range in meter
  Cmat(1:ncar,1:ngrid)=real(0.0,kind=q)
  do i=1,ncar
     x0=real(cardata%datamat(i,1),kind=q)
     y0=real(cardata%datamat(i,2),kind=q)
     do k=1,ngrid
        jj=mod(k,ny)
        if(jj==0) then
          jj=ny
          ii=k/ny
        else
          ii=(k-jj)/ny+1
        endif 
        x1=real((ii-1)*dx,kind=q)
        y1=real((jj-1)*dy,kind=q)
        mydist=sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
        if(mydist<=myradius) then
          Cmat(i,k)=real(1.0,kind=q)
        endif
     enddo
     !
     ! Normalizing
     !
     mysum=sum(Cmat(i,1:ngrid))
     if(mysum>1e-6) then !Nonzero
       Cmat(i,1:ngrid)=Cmat(i,1:ngrid)/mysum
     else
       write(ilog,*) 'Assemble_Cmat: the data point is not within the domain.'
       write(ilog,*) 'Please remove it from car survey data file.'
       stop
     endif
  enddo
  end subroutine Assemble_Cmat

  subroutine grid_location(k,ix,iy)
  integer, intent(in)  :: k
  integer, intent(out) :: ix,iy
  !
  if(mod(k,ny)==0) then !Right edge
    ix=k/ny
    iy=ny
  else !General case
    ix=k/ny+1
    iy=k-(ix-1)*ny
  endif
  end subroutine grid_location
end module class_Post
