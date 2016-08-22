!****************************************************************************************
! Purpose: 
!   Provide basic functions and subroutines for multiscale Bayesian integration 
!
! Author:
!   Jinsong Chen and Haruko Wainwright
!   Lawrence Berkeley National Lab
!   Earth Sciences Division
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
module class_Subfun
  use class_Global
  !
  implicit none
  !
  private
  public :: Subfun_cholesky_lower_lapack,Subfun_invert_pd_matrix_lapack, &
            Subfun_linear_interp1d,Subfun_mklMatmul, &
            Subfun_invert_pd_matrix_lapack2
  ! 
  integer, parameter :: dp=kind(1d0)
  !
  interface Subfun_mklMatmul
    module procedure Subfun_mklMatmul1
    module procedure Subfun_mklMatmul2
  end interface Subfun_mklMatmul
  !
Contains
  !  
  subroutine Subfun_mklMatmul1(m,p,n,tmpAmat,tmpBmat,rtnCmat)
  integer, intent(in)         :: m,p,n        !dimension size
  real(kind=q), intent(in)    :: tmpAmat(:,:) !dimension of mxp
  real(kind=q), intent(in)    :: tmpBmat(:,:) !dimension of pxn
  real(kind=q), intent(inout) :: rtnCmat(:,:) !dimension of mxn
  !
  integer                                    :: ierr(3)
  real(kind=dp), parameter                   :: alpha=1.0_dp
  real(kind=dp), parameter                   :: beta=0.0_dp
  real(kind=dp), dimension(:,:), allocatable :: Amat,Bmat,Cmat
  !
  allocate(Amat(m,p),STAT=ierr(1))
  allocate(Bmat(p,n),STAT=ierr(2))
  allocate(Cmat(m,n),STAT=ierr(3))
  if(any(ierr(1:3)>0)) then
    write(ilog,*) 'Subfun_mklMatmul2: Errors in allocating memory.'
    stop
  endif
  !
  Amat=real(tmpAmat,kind=dp)
  Bmat=real(tmpBmat,kind=dp)
  Cmat(1:m,1:n)=real(0.0,kind=dp)
  !
  ! Call MKL subroutine DGEMM
  !
  call dgemm('N','N',m,n,p,alpha,Amat,m,Bmat,p,beta,Cmat,m)
  !
  rtnCmat=real(Cmat,kind=q)
  if(allocated(Amat)) deallocate(Amat)
  if(allocated(Bmat)) deallocate(Bmat)
  if(allocated(Cmat)) deallocate(Cmat)
  end subroutine Subfun_mklMatmul1

  subroutine Subfun_mklMatmul2(m,p,n,tmpalpha,tmpbeta,tmpAmat,tmpBmat,tmpCmat)
  integer, intent(in)         :: m,p,n        !dimension size
  real(kind=q), intent(in)    :: tmpalpha,tmpbeta   !for dgemm
  real(kind=q), intent(in)    :: tmpAmat(:,:) !dimension of mxp
  real(kind=q), intent(in)    :: tmpBmat(:,:) !dimension of pxn
  real(kind=q), intent(inout) :: tmpCmat(:,:) !dimension of mxn
  !
  integer                                    :: ierr(3)
  real(kind=dp)                              :: alpha,beta
  real(kind=dp), dimension(:,:), allocatable :: Amat,Bmat,Cmat
  !
  allocate(Amat(m,p),STAT=ierr(1))
  allocate(Bmat(p,n),STAT=ierr(2))
  allocate(Cmat(m,n),STAT=ierr(3))
  if(any(ierr(1:3)>0)) then
    write(ilog,*) 'Subfun_mklMatmul2: Errors in allocating memory.'
    stop
  endif
  !
  Amat=real(tmpAmat,kind=dp)
  Bmat=real(tmpBmat,kind=dp)
  Cmat=real(tmpCmat,kind=dp)
  alpha=real(tmpalpha,kind=dp)
  beta=real(tmpbeta,kind=dp)
  !
  ! Call MKL subroutine DGEMM
  !
  call dgemm('N','N',m,n,p,alpha,Amat,m,Bmat,p,beta,Cmat,m)
  !
  tmpCmat=real(Cmat,kind=q)
  if(allocated(Amat)) deallocate(Amat)
  if(allocated(Bmat)) deallocate(Bmat)
  if(allocated(Cmat)) deallocate(Cmat)
  end subroutine Subfun_mklMatmul2

  function Subfun_linear_interp1d(ndim,xvec,yvec,xin) result (yout)
  integer, intent(in)      :: ndim       !Dimension
  real(kind=q), intent(in) :: xvec(:)
  real(kind=q), intent(in) :: yvec(:)
  real(kind=q), intent(in) :: xin
  real(kind=q)             :: yout
  !
  ! Assume xvec is in the increasing order
  !
  integer      :: i,k
  real(kind=q) :: dx,dy
  !
  if(xin<xvec(1)) then
    yout=yvec(1)
    return
  elseif(xin>xvec(ndim)) then
    yout=yvec(ndim)
    return
  else
    k=0
    do i=1,ndim
       if(xvec(i)<xin) k=k+1 
    enddo
    dx=xvec(k+1)-xvec(k)
    dy=yvec(k+1)-yvec(k)
    if(abs(dx)<1e-8) then
      yout=yvec(k)
      return
    else
      yout=yvec(k)+(dy/dx)*(xin-xvec(k))
      return
    endif
  endif
  end function Subfun_linear_interp1d

  subroutine Subfun_cholesky_lower_lapack(a,n,ierr)
  real, intent(inout)  :: a(:,:)
  integer, intent(in)  :: n
  integer, intent(out) :: ierr
  !
  ! Get a lower triangular matrix from Cholesky decomposition of a given positive
  ! definite matrix using LAPACK subroutines.
  !
  external                  :: dpptrf !From LAPACK
  integer                   :: i,j,k,locerr
  real(kind=q), allocatable :: packvect(:)
  !
  allocate(packvect(n*(n+1)/2),STAT=locerr)
  if(locerr>0) then
    write(6,*) 'cholesky_lower_lapack: Errors in allocating memory.'
    ierr=1
    return
  endif
  !
  ! Store the lower part of a 2D matrix following column order required by LAPACK
  !
  k=1
  do j=1,n
  do i=j,n
     packvect(k)=real(a(i,j),kind=q)
     k=k+1
  enddo
  enddo
  !
  ! LAPACK subroutine
  !
  call dpptrf('L',n,packvect,locerr)
  if(locerr>0) then
    write(*,*) 'Error Code=',locerr
    write(*,*) 'Subfun_cholesky_lower_lapack: Errors in calling LAPACK subroutine-dpptrf.'
    ierr=1
    return
  endif
  !
  ! Extract the Cholesky lower matrix from the packed 1D array following column order
  !
  a(1:n,1:n)=real(0.0,kind=q)
  k=1
  do j=1,n
  do i=j,n
     a(i,j)=real(packvect(k))
     k=k+1
  enddo
  enddo
  !
  ierr=0
  deallocate(packvect)
  end subroutine Subfun_cholesky_lower_lapack

  subroutine Subfun_invert_pd_matrix_lapack(a,n,ierr)
  real, intent(inout)  :: a(:,:)
  integer, intent(in)  :: n
  integer, intent(out) :: ierr
  !
  ! Get inverse matrix using Cholesky decomposition and LAPACK subroutines.
  !
  external                  :: dpptrf, dpptri !From LAPACK
  integer                   :: i,j,k,m,locerr
  real(kind=q), allocatable :: packvect(:)
  !
  m=n
  allocate(packvect(m*(m+1)/2),STAT=locerr)
  if(locerr>0) then
    write(6,*) 'cholesky_lower_lapack: Errors in allocating memory.'
    ierr=1
    return
  endif
  !
  ! Store the lower part of a 2D matrix following column order required by LAPACK
  !
  k=1
  do j=1,m
  do i=j,m
     packvect(k)=real(a(i,j),kind=q)
     k=k+1
  enddo
  enddo
  !
  ! Factorize matrix using LAPACK subroutine
  !
  call dpptrf('L',m,packvect,locerr)
  if(locerr.ne.0) then
    do i=1,n
       write(32,*) 'i,value:',i,a(i,i)
    enddo
    write(*,*) 'Error Code=',locerr
    write(*,*) 'Subfun_invert_pd_matrix_lapack: Errors in calling LAPACK subroutine-dpptrf.'
    ierr=1
    return
  else
    !
    ! Compute inverse of the matrix 
    !
    call dpptri('L',m,packvect,locerr)
    !
    ! Assign lower triangle first and then upper triangle
    !
    k=1
    do j=1,m
    do i=j,m
     a(i,j)=real(packvect(k))
     k=k+1
    enddo
    enddo
    !
    ! Upper triangle
    !
    do i=1,m
    do j=i+1,m
       a(i,j)=a(j,i)
    enddo
    enddo
  endif
  ierr=0
  deallocate(packvect)
  end subroutine Subfun_invert_pd_matrix_lapack

  subroutine Subfun_invert_pd_matrix_lapack2(a,n,ierr)
  real(kind=q), intent(inout) :: a(:,:)
  integer, intent(in)         :: n
  integer, intent(out)        :: ierr
  !
  ! Get inverse matrix using Cholesky decomposition and LAPACK subroutines
  ! using sgetrf, and sgetri (multi-threads)
  !
  external                  :: sgetrf,sgetri
  external                  :: dgetrf,dgetri
  integer                   :: info,lda,lwork,locerr(2)
  real, allocatable         :: swork(:)
  real(kind=q), allocatable :: dwork(:)
  integer, allocatable      :: ipiv(:)
  !
  lda=n
  lwork=64*n
  if(q==kind(1.0)) then
    allocate(swork(lwork),STAT=locerr(1))
  else
    allocate(dwork(lwork),STAT=locerr(1))
  endif
  allocate(ipiv(n),STAT=locerr(2))
  if(any(locerr(1:2)>0)) then
    write(6,*) 'cholesky_lower_lapack: Errors in allocating memory.'
    ierr=1
    return
  endif
  !
  ! Set multiple threads
  !
  !call omp_set_num_threads(24)
  !
  ! Factorization
  !
  if(q==kind(1.0)) then
    call sgetrf(n,n,a,n,ipiv,info)
  else
    call dgetrf(n,n,a,n,ipiv,info)
  endif
  if(info.ne.0) then
    ierr=1
    if(allocated(ipiv)) deallocate(ipiv)
    return
  endif
  !
  ! Compute inverse of matrix a
  !
  if(q==kind(1.0)) then
    call sgetri(n,a,lda,ipiv,swork,lwork,info)
  else 
    call dgetri(n,a,lda,ipiv,dwork,lwork,info)
  endif
  if(info.ne.0) then
    ierr=1
    if(allocated(swork)) deallocate(swork)
    if(allocated(dwork)) deallocate(dwork)
    if(allocated(ipiv)) deallocate(ipiv)
    return
  endif
  !
  ierr=0
  if(allocated(swork)) deallocate(swork)
  if(allocated(dwork)) deallocate(dwork)
  if(allocated(ipiv)) deallocate(ipiv)
  end subroutine Subfun_invert_pd_matrix_lapack2
end module class_Subfun
