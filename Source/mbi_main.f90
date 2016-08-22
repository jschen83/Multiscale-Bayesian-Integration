!****************************************************************************************
! Purpose:
!   Multiscale Bayesian integeration, transformed from the MATLAB code developed by
! Haruko Wainwright.
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
program main_samp
  use class_Global
  use class_MISC
  use class_Post
  !  
  implicit none
  !
  type(io_files)           :: iofiles
  character(len=maxstrlen) :: fn0,fn1,fn2,fn3,runfile
  !
  type(gridDomain) :: grid     !Hold general grid information
  type(geoModel)   :: gmod(3)  !Variogram model information
  type(surveyData) :: walkdata !Walking survey data
  type(surveyData) :: cardata  !Car survey data
  type(surveyData) :: airdata  !Air survey data
  !
  real, dimension(:,:), allocatable    :: estmean,eststd
  real, dimension(:,:), allocatable    :: smoothmat
  integer, dimension(:,:), allocatable :: lusemat
  !  
  real               :: start_time,end_time
  integer            :: i,j,k,ntotal,ierr(6)
  integer, parameter :: locDebug=1
  !
  ! Read runfile 
  !  
  if(iargc().lt.1)then
    write(*,*) 'Enter parameter file name:'
    read(*,'(a)') runfile
  else
    call getarg(1,runfile)
  endif
  !
  ! Read runfile containing all input file names and get dimension size.
  ! It must be called before all other subroutines in MISC module.
  !
  call CPU_Time(start_time)
  call MISC_ReadRunfile(runfile,iofiles)
  !
  ! Open files for saving debuging and log information
  !
  fn0(iofiles%is8:iofiles%ie8+1)=iofiles%output_dir(iofiles%is8:&
                                 iofiles%ie8)//'/'
  fn1=fn0(iofiles%is8:iofiles%ie8+1)//iofiles%output_file(iofiles%is9:&
                                 iofiles%ie9)//'.dbg'
  fn2=fn0(iofiles%is8:iofiles%ie8+1)//iofiles%output_file(iofiles%is9:&
                                   iofiles%ie9)//'.log'
  fn3=fn0(iofiles%is8:iofiles%ie8+1)//iofiles%output_file(iofiles%is9:&
                                   iofiles%ie9)//'.est'
  open(idbg,file=fn1,action='Write',iostat=ierr(1))
  open(ilog,file=fn2,action='Write',iostat=ierr(2))
  open(30,file=fn3,action='Write',iostat=ierr(3))
  if(any(ierr(1:3)>0)) then
    write(6,*) 'Main Program: Error opening log and/or debug file.'
    stop
  endif  
  !
  ! Read model setting
  !
  call MISC_ReadModelSetting(grid)
  !
  ! Allocate memory for smoothmat and lusemat
  !
  ntotal=grid%nx*grid%ny
  allocate(smoothmat(grid%nsmooth,2),STAT=ierr(1))
  allocate(lusemat(grid%nx*grid%ny,3),STAT=ierr(2))
  if(any(ierr(1:2)>0)) then
    write(*,*) 'MAIN: Errors in allocating smoothmat and lusemat.'
    stop
  endif
  !
  ! Read data 
  !
  call MISC_ReadGeoModel(gmod)
  call MISC_ReadSmoothMatrix(smoothmat)
  call MISC_ReadLanduseMatrix(lusemat)
  call MISC_ReadWalkingData(walkdata)
  call MISC_ReadCarData(cardata)
  call MISC_ReadAirborneData(airdata)
  !
  ! Call posterior module
  !
  call Post_InputPara(grid,gmod,smoothmat,lusemat,walkdata)
  if(allocated(smoothmat)) deallocate(smoothmat)
  if(allocated(lusemat)) deallocate(lusemat)
  !
  ntotal=grid%nx*grid%ny
  allocate(estmean(ntotal,2),STAT=ierr(1))
  allocate(eststd(ntotal,2),STAT=ierr(2))
  if(any(ierr(1:2)>0)) then
    write(*,*) 'MAIN: Errors in allocating postmean and poststd.'
    stop
  endif
  call CPU_Time(end_time)
  write(*,*) 'CPU Time for preparing input parameters and data:',(end_time-start_time)
  !
  call CPU_Time(start_time)
  call Post_CalcMeanCova(cardata,airdata,estmean,eststd)
  call CPU_Time(end_time)
  write(*,*) 'Total CPU time for computing prior and posterior estimates:',(end_time-start_time)
  !
  ! Output estimates
  !
  write(30, *) 'nx and ny:'
  write(30, *) grid%nx,grid%ny
  write(30,*) 'xloc  yloc  PriorMean  PostMean  PriorStd  PostStd'
  do j=1,grid%ny
  do i=1,grid%nx
     k=(j-1)*grid%nx+i
     write(30,100) (i-1)*grid%dx,(j-1)*grid%dy,estmean(k,1:2),eststd(k,1:2)
  enddo
  enddo
  100 format(1x,F10.3,2x,F10.3,4x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6)
  close(30)
  close(idbg)
  close(ilog)
  if(allocated(estmean)) deallocate(estmean)
  if(allocated(eststd)) deallocate(eststd)
  call MISC_Free(walkdata,cardata,airdata)
  call Post_Free
end program main_samp
