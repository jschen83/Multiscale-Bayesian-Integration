!****************************************************************************************
! Purpose:
!   Provide various functions or subroutines for reading data from files.
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
module class_MISC
use class_Global
!
implicit none
!
  private
  public :: MISC_ReadRunfile, MISC_Free, &
            MISC_ReadModelSetting, MISC_ReadGeoModel, &
            MISC_ReadSmoothMatrix, MISC_ReadLanduseMatrix, & 
            MISC_ReadWalkingData, MISC_ReadAirborneData, MISC_ReadCarData 
  !
  type(gridDomain) :: grid
  !
  ! I/O filename 
  !
  type(io_files) :: misc_iofiles
  !
Contains

  subroutine MISC_Free(walkdata,cardata,airdata) 
  type(surveyData), intent(inout) :: walkdata
  type(surveyData), intent(inout) :: cardata
  type(surveyData), intent(inout) :: airdata
  !
  ! Deallocate memory
  !
  if(allocated(walkdata%landuse)) deallocate(walkdata%landuse)
  if(allocated(walkdata%datamat)) deallocate(walkdata%datamat)
  !
  if(allocated(cardata%landuse)) deallocate(cardata%landuse)
  if(allocated(cardata%datamat)) deallocate(cardata%datamat)
  !
  if(allocated(airdata%landuse)) deallocate(airdata%landuse)
  if(allocated(airdata%datamat)) deallocate(airdata%datamat)
  end subroutine MISC_Free

  subroutine MISC_ReadModelSetting(mygrid) 
  type(gridDomain), intent(inout) :: mygrid
  !
  ! Read grid setting from a file
  !
  character(len=maxstrlen) :: title,tmpfn,fname
  integer                  :: ierr
  integer, parameter       :: locDebug=0
  !
  tmpfn=trim(misc_iofiles%input_dir(misc_iofiles%is0:misc_iofiles%ie0))//'/'
  fname=trim(tmpfn)//adjustl(misc_iofiles%model_file(misc_iofiles%is1:misc_iofiles%ie1))
  open(20,file=fname,action='READ',iostat=ierr)
  if(ierr>0) then
    write(ilog,*) fname
    write(ilog,*) 'MISC_ReadModelSetting: Errors in open the model setting file.'
    stop
  endif
  !
  read(20,'(a)') title
  read(20,*) mygrid%nx
  read(20,*) mygrid%ny
  read(20,*) mygrid%nsmooth
  read(20,*) mygrid%dx 
  read(20,*) mygrid%dy 
  read(20,*) mygrid%gmean
  if(locDebug==1) then
    write(*,*) trim(title)
    write(*,*) 'nx=',mygrid%nx
    write(*,*) 'ny=',mygrid%ny
    write(*,*) 'nsmooth=',mygrid%nsmooth
    write(*,*) 'dx=',mygrid%dx
    write(*,*) 'dy=',mygrid%dy
    write(*,*) 'nsmooth=',mygrid%nsmooth
    write(*,*) 'gmean=',mygrid%gmean
  endif
  close(20)
  grid=mygrid
  end subroutine MISC_ReadModelSetting 

  subroutine MISC_ReadGeoModel(mygmod) 
  type(geoModel), intent(inout) :: mygmod(:)
  !
  ! Read geostatistical models from a file
  !
  character(len=maxstrlen) :: title,tmpfn,fname
  integer                  :: i,ierr
  real                     :: tmpvec(2)
  integer, parameter       :: locDebug=0
  !
  tmpfn=trim(misc_iofiles%input_dir(misc_iofiles%is0:misc_iofiles%ie0))//'/'
  fname=trim(tmpfn)//adjustl(misc_iofiles%geostat_file(misc_iofiles%is2:misc_iofiles%ie2))
  open(20,file=fname,action='READ',iostat=ierr)
  if(ierr>0) then
    write(ilog,*) fname
    write(ilog,*) 'MISC_ReadGeoModel: Errors in open the geostatistical file.'
    stop
  endif
  !
  read(20,'(a)') title
  if(locDebug==1) write(*,*) trim(title)
  !
  do i=1,3  !Assuming only three types of geostatistical models
     read(20,*) mygmod(i)%variogram
     read(20,*) tmpvec(1),tmpvec(2)
     mygmod(i)%crlengthx=tmpvec(1)
     mygmod(i)%crlengthy=tmpvec(2)
     read(20,*) mygmod(i)%variance
     read(20,*) mygmod(i)%nugget
     read(20,*) mygmod(i)%relnugget
  enddo
  close(20)
  !
  if(locDebug==1) then
     do i=1,3  !Assuming only three types of geostatistical models
        write(*,*) 'Variogram ID=',mygmod(i)%variogram
        write(*,*) 'Crlengthx=',mygmod(i)%crlengthx,'Crlengthy=',mygmod(i)%crlengthy
        write(*,*) 'Variance=',mygmod(i)%variance
        write(*,*) 'Nugget=',mygmod(i)%nugget
        write(*,*) 'Rel-nugget=',mygmod(i)%relnugget
     enddo
  endif
  end subroutine MISC_ReadGeoModel 

  subroutine MISC_ReadSmoothMatrix(smoothmat) 
  real, dimension(:,:), intent(inout) :: smoothmat
  !
  ! Read smoothing matrix
  !
  character(len=maxstrlen) :: title,tmpfn,fname
  integer                  :: i,ierr,ntotal
  integer, parameter       :: locDebug=0
  !
  tmpfn=trim(misc_iofiles%input_dir(misc_iofiles%is0:misc_iofiles%ie0))//'/'
  fname=trim(tmpfn)//adjustl(misc_iofiles%smoothing_file(misc_iofiles%is3:misc_iofiles%ie3))
  open(20,file=fname,action='READ',iostat=ierr)
  if(ierr>0) then
    write(ilog,*) fname
    write(ilog,*) 'MISC_ReadSmoothMatrix: Errors in open the smoothing file.'
    stop
  endif
  !
  read(20,'(a)') title
  if(locDebug==1) write(idbg,*) trim(title)
  !
  read(20,*) ntotal
  if(locDebug==1) write(idbg,*) 'ntotal=',ntotal
  !
  if(ntotal>grid%nsmooth) then
    write(*,*) 'Number of nsmooth should be larger. Check the input domain parameter file.'
    stop
  endif
  !
  read(20,'(a)') title
  if(locDebug==1) write(idbg,*) trim(title)
  do i=1,ntotal
     read(20,*) smoothmat(i,1),smoothmat(i,2)
     if(locDebug==1) write(idbg,*) smoothmat(i,1),smoothmat(i,2)
  enddo
  close(20)
  end subroutine MISC_ReadSmoothMatrix 

  subroutine MISC_ReadLanduseMatrix(lusemat) 
  integer, dimension(:,:), intent(inout) :: lusemat
  !
  ! Read land-use matrix
  !
  character(len=maxstrlen) :: title,tmpfn,fname
  integer                  :: i,ierr,tmpnx,tmpny,ntotal
  integer                  :: tmpint(3)
  integer, parameter       :: locDebug=0
  !
  tmpfn=trim(misc_iofiles%input_dir(misc_iofiles%is0:misc_iofiles%ie0))//'/'
  fname=trim(tmpfn)//adjustl(misc_iofiles%landuse_file(misc_iofiles%is4:misc_iofiles%ie4))
  open(20,file=fname,action='READ',iostat=ierr)
  if(ierr>0) then
    write(ilog,*) fname
    write(ilog,*) 'MISC_ReadLanduseMatrix: Errors in open the landuse file.'
    stop
  endif
  !
  read(20,'(a)') title
  if(locDebug==1) write(idbg,*) trim(title)
  !
  read(20,*) tmpnx,tmpny
  ntotal=tmpnx*tmpny
  if(locDebug==1) write(idbg,*) 'nx=',tmpnx,'ny=',tmpny
  !
  if(ntotal>grid%nx*grid%ny) then
    write(*,*) 'Number of nx*ny should be larger. Check the input domain parameter file.'
    stop
  endif
  !
  read(20,'(a)') title
  if(locDebug==1) write(idbg,*) trim(title)
  do i=1,ntotal
     read(20,*) tmpint(1:3)
     lusemat(i,1:3)=tmpint(1:3) 
     if(locDebug==1) write(idbg,*) lusemat(i,1:3)
  enddo
  close(20)
  end subroutine MISC_ReadLanduseMatrix 

  subroutine MISC_ReadWalkingData(walkdata) 
  type(surveyData), intent(inout) :: walkdata
  !
  ! Read car survey data
  !
  character(len=maxstrlen) :: title,tmpfn,fname
  integer                  :: i,ntotal,ierr(2)
  real                     :: tmpreal(4)
  integer, parameter       :: locDebug=0
  !
  tmpfn=trim(misc_iofiles%input_dir(misc_iofiles%is0:misc_iofiles%ie0))//'/'
  fname=trim(tmpfn)//adjustl(misc_iofiles%walking_data_file(misc_iofiles%is5:misc_iofiles%ie5))
  open(20,file=fname,action='READ',iostat=ierr(1))
  if(ierr(1)>0) then
    write(ilog,*) fname
    write(ilog,*) 'MISC_ReadWalkingData: Errors in open the walking data file.'
    stop
  endif
  !
  read(20,'(a)') title 
  if(locDebug==1) write(idbg,*) trim(title)
  !
  ! Initial those with ZERO (not use them)
  !
  walkdata%muvec(1:3)=0.0
  walkdata%stdvec(1:3)=0.0
  walkdata%varvec(1:3)=0.0
  !
  read(20,*) ntotal
  if(locDebug==1) write(idbg,*) 'ntotal=',ntotal
  walkdata%num=ntotal
  !
  allocate(walkdata%landuse(ntotal),STAT=ierr(1))
  allocate(walkdata%datamat(ntotal,3),STAT=ierr(2))
  if(any(ierr(1:2)>0)) then
    write(ilog,*) 'MISC_ReadWalkingData: Errors in allocating.'
    stop
  endif
  !
  read(20,'(a)') title
  if(locDebug==1) write(idbg,*) trim(title)
  do i=1,ntotal
     read(20,*) tmpreal(1:4)
     walkdata%datamat(i,1:3)=tmpreal(1:3)
     walkdata%landuse(i)=int(tmpreal(4))
     if(locDebug==1) then
       write(idbg,*) walkdata%datamat(i,1:3),walkdata%landuse(i)
     endif
  enddo
  close(20)
  end subroutine MISC_ReadWalkingData 

  subroutine MISC_ReadCarData(cardata) 
  type(surveyData), intent(inout) :: cardata
  !
  ! Read car survey data
  !
  character(len=maxstrlen) :: title,tmpfn,fname
  integer                  :: i,ntotal,ierr(2)
  real                     :: tmpreal(4)
  integer, parameter       :: locDebug=0
  !
  tmpfn=trim(misc_iofiles%input_dir(misc_iofiles%is0:misc_iofiles%ie0))//'/'
  fname=trim(tmpfn)//adjustl(misc_iofiles%car_data_file(misc_iofiles%is6:misc_iofiles%ie6))
  open(20,file=fname,action='READ',iostat=ierr(1))
  if(ierr(1)>0) then
    write(ilog,*) fname
    write(ilog,*) 'MISC_ReadCarData: Errors in open the car data file.'
    stop
  endif
  !
  read(20,'(a)') title  !1st line
  if(locDebug==1) write(idbg,*) trim(title)
  read(20,'(a)') title  !2nd line
  if(locDebug==1) write(idbg,*) trim(title)
  !
  do i=1,3 !assuming only three land-use type
     read(20,*) cardata%stdvec(i),cardata%muvec(i) !std and mean
     read(20,*) cardata%varvec(i)                  !variance of measurement errors
  enddo
  !
  read(20,*) ntotal
  if(locDebug==1) write(idbg,*) 'ntotal=',ntotal
  cardata%num=ntotal
  !
  allocate(cardata%landuse(ntotal),STAT=ierr(1))
  allocate(cardata%datamat(ntotal,3),STAT=ierr(2))
  if(any(ierr(1:2)>0)) then
    write(ilog,*) 'MISC_ReadCarData: Errors in allocating.'
    stop
  endif
  !
  read(20,'(a)') title
  if(locDebug==1) write(idbg,*) trim(title)
  do i=1,ntotal
     read(20,*) tmpreal(1:4)
     cardata%datamat(i,1:3)=tmpreal(1:3)
     cardata%landuse(i)=int(tmpreal(4))
     if(locDebug==1) then
       write(idbg,*) cardata%datamat(i,1:3),cardata%landuse(i)
     endif
  enddo
  close(20)
  end subroutine MISC_ReadCarData 

  subroutine MISC_ReadAirborneData(airdata) 
  type(surveyData), intent(inout) :: airdata
  !
  ! Read airborne survey data
  !
  character(len=maxstrlen) :: title,tmpfn,fname
  integer                  :: i,ntotal,ierr(2)
  real                     :: tmpreal(4)
  integer, parameter       :: locDebug=0
  !
  tmpfn=trim(misc_iofiles%input_dir(misc_iofiles%is0:misc_iofiles%ie0))//'/'
  fname=trim(tmpfn)//adjustl(misc_iofiles%airborne_data_file(misc_iofiles%is7:misc_iofiles%ie7))
  open(20,file=fname,action='READ',iostat=ierr(1))
  if(ierr(1)>0) then
    write(ilog,*) fname
    write(ilog,*) 'MISC_ReadAirborneData: Errors in open the airborne data file.'
    stop
  endif
  !
  read(20,'(a)') title  !1st line
  if(locDebug==1) write(idbg,*) trim(title)
  read(20,'(a)') title  !2nd line
  if(locDebug==1) write(idbg,*) trim(title)
  !
  do i=1,3 !assuming only three land-use type
     read(20,*) airdata%stdvec(i),airdata%muvec(i) !std and mean
     if(locDebug==1) write(idbg,*) airdata%stdvec(i),airdata%muvec(i)
     read(20,*) airdata%varvec(i) !Variance of measurement errors
     if(locDebug==1) write(idbg,*) airdata%varvec(i)
  enddo
  !
  read(20,*) ntotal
  airdata%num=ntotal
  if(locDebug==1) write(idbg,*) 'ntotal=',ntotal
  !
  allocate(airdata%landuse(ntotal),STAT=ierr(1))
  allocate(airdata%datamat(ntotal,3),STAT=ierr(2))
  if(any(ierr(1:2)>0)) then
    write(ilog,*) 'MISC_ReadAirborneData: Errors in allocating.'
    stop
  endif
  !
  read(20,'(a)') title
  if(locDebug==1) write(idbg,*) trim(title)
  do i=1,ntotal
     read(20,*) tmpreal(1:4)
     airdata%datamat(i,1:3)=tmpreal(1:3)
     airdata%landuse(i)=int(tmpreal(4))
     if(locDebug==1) then
       write(idbg,*) airdata%datamat(i,1:3),airdata%landuse(i)
     endif
  enddo
  close(20)
  end subroutine MISC_ReadAirborneData 

  subroutine MISC_ReadRunfile(runfile,iofiles)
  character(len=*), intent(in) :: runfile
  type(io_files), intent(out)  :: iofiles
  !
  character(len=maxstrlen) :: temp,title
  integer                  :: i,n,ix,is,ie,ierr
  integer, parameter       :: locDebug=0
  !
  open(20,file=runfile(1:len_trim(runfile)),action='READ',iostat=ierr)
  if(ierr >0) then
    write(*,*) 'MISC_ReadRunfile: Errors in openning'
    write(*,'(a)') runfile(1:len_trim(runfile))
    stop
  endif
  !  
  read(20,'(a)') title
  read(20,*) globDebug
  write(*,*) trim(title)
  write(*,*) 'Global Debug level=',globDebug
  !
  ! Input directory
  !
  read(20,'(a)') temp
  is=begwrd(temp,1)
  ie=endwrd(temp,1)
  iofiles%is0=is
  iofiles%ie0=ie
  n=ie-is+1
  if(n>maxstrlen) then
     write(*,100) temp(is:ie),maxstrlen
     stop
  endif
  iofiles%input_dir(1:n)=temp(is:ie)  
  if(locDebug==1) write(*,*) trim(iofiles%input_dir(1:n))
  ! 
  ! Model setting file
  !
  read(20,'(a)') temp
  is=begwrd(temp,1)
  ie=endwrd(temp,1)
  n=ie-is+1
  iofiles%is1=1
  iofiles%ie1=n
  if(n>maxstrlen) then
    write(*,100) temp(is:ie),maxstrlen
    stop
  endif
  iofiles%model_file(1:n)=temp(is:ie)
  if(locDebug==1) write(*,*) trim(iofiles%model_file(1:n))
  !
  ! Geostatistical model parameter file
  !
  read(20,'(a)') temp
  is=begwrd(temp,1)
  ie=endwrd(temp,1)
  n=ie-is+1
  iofiles%is2=1
  iofiles%ie2=n
  if(n>maxstrlen) then
    write(*,100) temp(is:ie),maxstrlen
    stop
  endif
  iofiles%geostat_file(1:n)=temp(is:ie)
  if(locDebug==1) write(*,*) trim(iofiles%geostat_file(1:n))
  !  
  ! Smoothing weight file
  !
  read(20,'(a)') temp
  is=begwrd(temp,1)
  ie=endwrd(temp,1)
  n=ie-is+1
  iofiles%is3=1
  iofiles%ie3=n
  if(n>maxstrlen) then
    write(*,100) temp(is:ie),maxstrlen
    stop
  endif
  iofiles%smoothing_file(1:n)=temp(is:ie)
  if(locDebug==1) write(*,*) trim(iofiles%smoothing_file(1:n))
  !  
  ! Landuse file
  !
  read(20,'(a)') temp
  is=begwrd(temp,1)
  ie=endwrd(temp,1)
  n=ie-is+1
  iofiles%is4=1
  iofiles%ie4=n
  if(n>maxstrlen) then
    write(*,100) temp(is:ie),maxstrlen
    stop
  endif
  iofiles%landuse_file(1:n)=temp(is:ie)
  if(locDebug==1) write(*,*) trim(iofiles%landuse_file(1:n))
  !  
  ! Walking data file
  !
  read(20,'(a)') temp
  is=begwrd(temp,1)
  ie=endwrd(temp,1)
  n=ie-is+1
  iofiles%is5=1
  iofiles%ie5=n
  if(n>maxstrlen) then
    write(*,100) temp(is:ie),maxstrlen
    stop
  endif
  iofiles%walking_data_file(1:n)=temp(is:ie)
  if(locDebug==1) write(*,*) trim(iofiles%walking_data_file(1:n))
  !  
  ! Car data file
  !
  read(20,'(a)') temp
  is=begwrd(temp,1)
  ie=endwrd(temp,1)
  n=ie-is+1
  iofiles%is6 = 1
  iofiles%ie6 = n
  if(n>maxstrlen) then
    write(*,100) temp(is:ie),maxstrlen
    stop
  endif
  iofiles%car_data_file(1:n)=temp(is:ie)
  if(locDebug==1) write(*,*) trim(iofiles%car_data_file(1:n))
  !  
  ! Airborne data file
  !
  read(20,'(a)') temp
  is=begwrd(temp,1)
  ie=endwrd(temp,1)
  n=ie-is+1
  iofiles%is7 = 1
  iofiles%ie7 = n
  if(n>maxstrlen) then
    write(*,100) temp(is:ie),maxstrlen
    stop
  endif
  iofiles%airborne_data_file(1:n)=temp(is:ie)
  if(locDebug==1) write(*,*) trim(iofiles%airborne_data_file(1:n))
  !
  ! Output directory name 
  !
  read(20,'(a)') temp
  is=begwrd(temp,1)
  ie=endwrd(temp,1)
  iofiles%is8=is
  iofiles%ie8=ie
  n=ie-is+1
  if(n>maxstrlen) then
    write(*,100) temp(is:ie),maxstrlen
  stop
  endif
  iofiles%output_dir(1:n)=temp(is:ie)
  if(locDebug==1) write(*,*) trim(iofiles%output_dir(1:n))
  ! 
  ! Output file name 
  !
  read(20,'(a)') temp
  is=begwrd(temp,1)
  ie=endwrd(temp,1)
  iofiles%is9 = is
  iofiles%ie9 = ie
  n=ie-is+1
  if(n>maxstrlen) then
    write(*,100) temp(is:ie), maxstrlen
  stop
  endif
  iofiles%output_file(1:n) = temp(is:ie)
  if(locDebug==1) write(*,*) trim(iofiles%output_file(1:n))
  close(20)
  !
  misc_iofiles=iofiles
  100 format('File: ',a,' length > maxstrlen(',i3,'), stopping.') 
  end subroutine MISC_ReadRunfile

  function begwrd(string,iwrd) result (rvint)
  character(len=*), intent(in) :: string
  integer, intent(in)          :: iwrd
  integer                      :: rvint
  !
  ! Returns the index of the first non-blank character in the iwrd'th
  ! non-blank word (word are seperated by spaces, tabs or commas).
  ! Returns len if iwrd'th word is not found.
  !
  integer   :: i, nword
  logical*2 :: wasblk
  intrinsic :: len
  !
  wasblk=.true.
  nword= 0
  do i=1,len(string)
     if(string(i:i).eq.' ' .or. string(i:i).eq.',' .or. string(i:i).eq.'  ') then
       !
       ! Current character is blank
       !
       wasblk=.true.
     else
       if(wasblk)  nword=nword+1
       wasblk=.false.
       if(nword.eq.iwrd) then
         rvint=i
         return
       endif
     endif
  enddo
  rvint= len(string)
  return
  end function begwrd

  function endwrd(string,iwrd) result (rvint)
  character(len=*), intent(in) :: string
  integer, intent(in)          :: iwrd
  integer                      :: rvint
  !
  ! Returns the index of the last non-blank character in the iwrd'th
  ! non-blank word (word are seperated by spaces, tabs or commas).
  ! Returns len if iwrd'th word is not found.
  !
  integer   :: i, nword
  logical*2 :: wasblk
  intrinsic :: len
  !
  wasblk=.true.
  nword= 0
  do i=1,len(string)
     if(string(i:i).eq.' ' .or. string(i:i).eq.',' .or. string(i:i).eq.'  ') then
       !
       ! Current character is blank
       !
       wasblk=.true.
       if(nword.eq.iwrd) RETURN
     else
       if(wasblk) nword=nword+1
       wasblk= .false.
       if(nword.eq.iwrd) rvint=i
     endif
  enddo
  rvint= len(string)
  return
  end function endwrd
end module class_MISC
