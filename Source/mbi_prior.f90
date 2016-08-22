!****************************************************************************************
! Purpose:
!   Perform spatial interpolation using kring (geostatistical method). These are 
! transformations of the MATLAB code developed by Haruko Wainwright.
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
module class_Prior
  use class_Global
  use class_Subfun
  !
  implicit none
  !
  private
  !
  public :: Prior_Free, Prior_InputPara, Prior_CalcMeanCova
  !
  integer :: nx,ny !Domain size 
  integer :: ngrid !Total # of grid points 
  integer :: nobs  !Total # of walking survey points 
  !
  real               :: dx,dy !discretization
  real               :: gmu   !overall mean values
  real, dimension(3) :: muvec,stdvec,varvec
  !
  real, dimension(:), allocatable    :: xloc
  real, dimension(:), allocatable    :: yloc
  real, dimension(:), allocatable    :: dose
  integer, dimension(:), allocatable :: landuse
  type(geoModel)                     :: gstat(3)  !Variogram model information
  !
contains
  !
  subroutine Prior_Free
  !
  ! Deallocate memory
  !
  if(allocated(xloc)) deallocate(xloc)
  if(allocated(yloc)) deallocate(yloc)
  if(allocated(dose)) deallocate(dose)
  if(allocated(landuse)) deallocate(landuse)
  end subroutine Prior_Free

  subroutine Prior_InputPara(grid,gmod,walkdata)
  type(gridDomain) :: grid     !Hold general grid information
  type(geoModel)   :: gmod(3)  !Variogram model information
  type(surveyData) :: walkdata !Walking survey data
  !
  ! Input domain, variogram models, and walking survey data
  !
  integer            :: ierr(6)
  integer, parameter :: locDebug=0
  !
  nx=grid%nx
  ny=grid%ny
  dx=grid%dx
  dy=grid%dy
  gmu=grid%gmean
  ngrid=nx*ny
  !
  ! Allocate memory
  !
  nobs=walkdata%num
  allocate(xloc(nobs),STAT=ierr(1))
  allocate(yloc(nobs),STAT=ierr(2))
  allocate(dose(nobs),STAT=ierr(3))
  allocate(landuse(nobs),STAT=ierr(4))
  if(any(ierr(1:4)>0)) then
    write(ilog,*) 'Prior_InputPara: Errors in allocating'
    stop
  endif
  gstat(1:3)=gmod(1:3)
  muvec(1:3)=walkdata%muvec(1:3)
  stdvec(1:3)=walkdata%stdvec(1:3)
  varvec(1:3)=walkdata%varvec(1:3)
  xloc(1:nobs)=walkdata%datamat(1:nobs,1)
  yloc(1:nobs)=walkdata%datamat(1:nobs,2)
  dose(1:nobs)=walkdata%datamat(1:nobs,3)
  landuse(1:nobs)=walkdata%landuse(1:nobs)
  end subroutine Prior_InputPara

  subroutine Prior_CalcMeanCova(lusemat,rtnmean,rtncova)
  integer, intent(in)         :: lusemat(:,:)
  real(kind=q), intent(inout) :: rtnmean(:)
  real(kind=q), intent(inout) :: rtncova(:,:)
  !
  ! Calculate mean and standard deviation using kring by conditioning 
  ! to walking survey data.
  !
  type(geoModel) :: mymodel
  integer        :: i,j,k,ng,nw,ii,jj
  integer        :: ierr(10)
  real           :: start_time,end_time
  !
  integer, dimension(:), allocatable        :: myindex1,myindex2
  real(kind=q)                              :: distx,disty
  real(kind=q), dimension(:), allocatable   :: beta,tmpvec1,tmpvec2
  real(kind=q), dimension(:,:), allocatable :: S11,S21,S22,invS11
  real(kind=q), dimension(:,:), allocatable :: R11,R21,R22,invR11
  real(kind=q), dimension(:,:), allocatable :: tR21,tmpMat1,tmpMat2
  !
  integer, parameter :: locDebug=0
  ! 
  if(locDebug==1) then
    write(idbg,*) '***START of Prior_CalcMeanStd***'
    write(*,*) '***START of Prior_CalcMeanStd***'
  endif
  !
  allocate(S11(nobs,nobs),STAT=ierr(1))
  allocate(S21(ngrid,nobs),STAT=ierr(2))
  allocate(S22(ngrid,ngrid),STAT=ierr(3))
  if(any(ierr(1:3)>0)) then
    write(ilog,*) 'Prior_CalcMeanStd: Errors in allocating memory'
    stop
  endif
  ! 
  ! Initializing matrices with zero values
  !
  S11(1:nobs,1:nobs)=real(0.0,kind=q)
  S21(1:ngrid,1:nobs)=real(0.0,kind=q)
  S22(1:ngrid,1:ngrid)=real(0.0,kind=q)
  !
  ! Main loop for constructing S11, S21, and S22 matrices
  !
  if(locDebug==1) call CPU_Time(start_time)
  !
  do k=1,3
     ng=count(lusemat(:,3)==k)
     nw=count(landuse==k)
     !
     ! If no the given land-use type, it should not have observations of the type.
     ! If it does, there is inconsistency between input landuse-id and walking data.
     !
     if(ng==0) then
       if(nw>0) then
         write(ilog,*) 'Prior_CalcMeanCova: Landuse and observation types are inconsistent.'
         write(ilog,*) 'Please check and re-run.'
         stop
       else
         cycle
       endif
     endif
     !
     allocate(myindex1(nw),STAT=ierr(1))
     allocate(myindex2(ng),STAT=ierr(2))
     allocate(R11(nw,nw),STAT=ierr(3))
     allocate(R21(ng,nw),STAT=ierr(4))
     allocate(R22(ng,ng),STAT=ierr(5))
     allocate(tR21(nw,ng),STAT=ierr(6))
     allocate(tmpMat1(nw,ng),STAT=ierr(7))
     allocate(tmpMat2(ng,ng),STAT=ierr(8))
     allocate(invR11(nw,nw),STAT=ierr(9))
     if(any(ierr(1:9)>0)) then
       write(ilog,*) 'Prior_CalcMeanStd: Errors in allocating memory-1.'
       stop
     endif
     !
     ! Find indices for a given land-use type
     !
     call GetIndex(landuse,nobs,k,myindex1)       !In walking data points
     call GetIndex(lusemat(:,3),ngrid,k,myindex2) !In overall grids
     !
     ! Do conditioning
     !
     mymodel=gstat(k)
     do i=1,ng
     do j=1,ng
        ii=myindex2(i)
        jj=myindex2(j)
        distx=real((lusemat(ii,1)-lusemat(jj,1))*dx,kind=q)
        disty=real((lusemat(ii,2)-lusemat(jj,2))*dy,kind=q)
        !
        ! Consider nugget effects 
        !
        R22(i,j)=CorrCoeff(distx,disty,mymodel)*real(1.0-mymodel%nugget,kind=q)
        if(i==j) R22(i,j)=R22(i,j)+real(mymodel%nugget,kind=q)
     enddo
     enddo
     !
     if(nw>0) then
       do i=1,nw
       do j=1,nw
          ii=myindex1(i)
          jj=myindex1(j)
          distx=real(xloc(ii)-xloc(jj),kind=q)
          disty=real(yloc(ii)-yloc(jj),kind=q)
          R11(i,j)=CorrCoeff(distx,disty,mymodel)*real(1.0-mymodel%nugget,kind=q)
          if(i==j) R11(i,j)=R11(i,j)+real(mymodel%nugget,kind=q)
          !
          ! Construct global S11 matrix
          !
          S11(ii,jj)=R11(i,j)*real(mymodel%variance,kind=q)
       enddo
       enddo
       !
       ! Do inversion
       !
       invR11(1:nw,1:nw)=R11(1:nw,1:nw)
       call Subfun_invert_pd_matrix_lapack2(invR11,nw,ierr(1))
       if(ierr(1)>0) then
         write(ilog,*) 'Prior_CalcMeanCova: Errors in inverting R11 matrix.'
         stop
       endif
       !
       ! Cross-covariance
       !
       do i=1,ng
       do j=1,nw
          ii=myindex2(i)
          jj=myindex1(j)
          distx=real((lusemat(ii,1)-1)*dx-xloc(jj),kind=q)
          disty=real((lusemat(ii,2)-1)*dy-yloc(jj),kind=q)
          R21(i,j)=CrossCorrCoeff(distx,disty,mymodel)
          S21(ii,jj)=R21(i,j)*real(mymodel%variance,kind=q)
       enddo
       enddo
       tR21(1:nw,1:ng)=transpose(R21)
       !
       ! Update R22 by conditioning to walking data whenthey are available
       !
       call Subfun_mklMatmul(nw,nw,ng,invR11,tR21,tmpMat1)
       call Subfun_mklMatmul(ng,nw,ng,R21,tmpMat1,tmpMat2)
       R22(1:ng,1:ng)=R22(1:ng,1:ng)-tmpMat2(1:ng,1:ng)
     endif
     !
     ! Construct global S22 matrix
     !
     do i=1,ng
     do j=1,ng
          ii=myindex2(i)
          jj=myindex2(j)
          S22(ii,jj)=R22(i,j)*real(mymodel%variance,kind=q)
     enddo
     enddo
     !
     ! Deallocate memory
     !
     if(allocated(myindex1)) deallocate(myindex1)
     if(allocated(myindex2)) deallocate(myindex2)
     if(allocated(R11)) deallocate(R11)
     if(allocated(R21)) deallocate(R21)
     if(allocated(R22)) deallocate(R22)
     if(allocated(invR11)) deallocate(invR11)
     if(allocated(tR21)) deallocate(tR21)
     if(allocated(tmpMat1)) deallocate(tmpMat1)
     if(allocated(tmpMat2)) deallocate(tmpMat2)
  enddo
  if(locDebug==1) then
    call CPU_Time(end_time)
    write(*,*) 'CPU time for S11, S21, and S22 matrices:',(end_time-start_time)
  endif
  !
  ! Calculate conditional mean and covariance
  !
  if(locDebug==1) call CPU_Time(start_time)
  !
  allocate(tmpvec1(nobs),STAT=ierr(1))
  allocate(tmpvec2(nobs),STAT=ierr(2))
  allocate(beta(ngrid),STAT=ierr(3))
  allocate(invS11(nobs,nobs),STAT=ierr(4))
  if(any(ierr(1:4)>0)) then
     write(ilog,*) 'Prior_CalcMeanStd: Errors in allocating memory-2.'
     stop
  endif
  !
  ! Assign walking data to the nearest grid points and compute the 
  ! difference between measurements and overall mean.
  !
  beta(1:ngrid)=gmu  !Overall mean
  do i=1,nobs
     ii=floor(xloc(i)/dx)+1
     jj=floor(yloc(i)/dy)+1
     k=(ii-1)*ny+jj
     tmpvec1(i)=dose(i)-beta(k)
  enddo
  !
  ! Invert S11 matrix
  !
  invS11(1:nobs,1:nobs)=S11(1:nobs,1:nobs)
  call Subfun_invert_pd_matrix_lapack2(invS11,nobs,ierr(1))
  if(ierr(1)>0) then
     write(ilog,*) 'Prior_CalcMeanStd: Errors in inverting S11 matrix.'
     stop
  endif
  tmpvec2=matmul(invS11,tmpvec1)
  if(allocated(S11)) deallocate(S11)
  if(allocated(invS11)) deallocate(invS11)
  !
  ! Return conditional mean
  !
  if(nobs>0) then
    rtnmean(1:ngrid)=beta(1:ngrid)+matmul(S21,tmpvec2)
  else
    rtnmean(1:ngrid)=beta(1:ngrid)
  endif
  if(allocated(S21)) deallocate(S21)
  if(allocated(beta)) deallocate(beta)
  if(allocated(tmpvec1)) deallocate(tmpvec1)
  if(allocated(tmpvec2)) deallocate(tmpvec2)
  !
  ! Return conditional covariance 
  !
  rtncova(1:ngrid,1:ngrid)=S22(1:ngrid,1:ngrid)
  !
  if(locDebug==1) then
    call CPU_Time(end_time)
    write(*,*) 'CPU time for conditional mean and covariance:',(end_time-start_time)
  endif
  !
  ! Output prior mean and std for debuging
  !
  if(locDebug==1) then
    open(22,file='prior_mean_std.out',action='Write',iostat=ierr(1))
    do i=1,ngrid
       write(22,*) real(rtnmean(i)),real(sqrt(S22(i,i)))
    enddo 
    close(22)
  endif
  if(allocated(S22)) deallocate(S22)
  !
  if(locDebug==1) then
    write(*,*) '***END of Prior_CalcMeanCova***'
    write(idbg,*) '***END of Prior_CalcMeanCova***'
  endif
  end subroutine Prior_CalcMeanCova

  function CorrCoeff(distx,disty,mymodel) result(rtnval)
  real(kind=q), intent(in)   :: distx
  real(kind=q), intent(in)   :: disty
  type(geoModel), intent(in) :: mymodel
  real(kind=q)               :: rtnval
  !
  ! This function is just for correlation coefficient function calculation
  !
  real(kind=q) :: tmpx,tmpy,tmpdist
  !
  tmpx=distx/real(mymodel%crlengthx,kind=q)
  tmpy=disty/real(mymodel%crlengthy,kind=q)
  tmpdist=sqrt(tmpx*tmpx+tmpy*tmpy)
  rtnval=exp(-tmpdist)
  end function CorrCoeff

  function CrossCorrCoeff(distx,disty,mymodel) result(rtnval)
  real(kind=q), intent(in)   :: distx
  real(kind=q), intent(in)   :: disty
  type(geoModel), intent(in) :: mymodel
  real(kind=q)               :: rtnval
  !
  ! This function is just for cross-correlation coefficient function calculation
  !
  real(kind=q) :: tmpx,tmpy,tmpdist
  !
  tmpx=distx/real(mymodel%crlengthx,kind=q)
  tmpy=disty/real(mymodel%crlengthy,kind=q)
  tmpdist=sqrt(tmpx*tmpx+tmpy*tmpy)
  rtnval=exp(-tmpdist)*real(1.0-mymodel%nugget,kind=q)
  end function CrossCorrCoeff

  subroutine GetIndex(lusevec,ndim,landtype,rtnvec)
  integer, intent(in)    :: lusevec(:)
  integer, intent(in)    :: ndim
  integer, intent(in)    :: landtype
  integer, intent(inout) :: rtnvec(:)
  !
  integer :: i,k
  !
  k=1
  do i=1,ndim
     if(lusevec(i)==landtype) then
       rtnvec(k)=i
       k=k+1
     endif
  enddo 
  end subroutine GetIndex
end module class_Prior
