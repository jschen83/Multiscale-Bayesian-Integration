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

module class_Global
  !
  implicit none
  !
  ! Global parameters
  !
  integer, parameter      :: q=kind(1.0d0)
  integer, parameter      :: maxstrlen=256
  integer, parameter      :: idbg=40      !Output debug info
  integer, parameter      :: ilog=42      !Output log of running 
  real(kind=q), parameter :: pi=3.1415926_q
  !
  ! 1-Output debug information
  ! 0-Not output debug information
  !
  integer :: globDebug=0  !Default not output debug information
  !
  ! Grid Setting 
  !
  type gridDomain
    integer :: nx,ny
    integer :: nsmooth
    real    :: dx,dy
    real    :: gmean
  end type gridDomain
  !
  ! Geostatistical Model 
  !
  type geoModel
    integer :: variogram
    real    :: variance
    real    :: crlengthx
    real    :: crlengthy
    real    :: nugget
    real    :: relnugget
  end type geoModel
  !
  ! Survey data
  !
  type surveyData
    integer              :: num               !Total # of data points 
    real                 :: muvec(3)          !mean dose
    real                 :: stdvec(3)         !std of dose
    real                 :: varvec(3)         !measurements of doses
    integer, allocatable :: landuse(:)        !Landuse type
    real, allocatable    :: datamat(:,:)      !Data matrix
  end type surveyData
  ! 
  ! Get file names
  !
  type io_files
    character(len=maxstrlen) :: input_dir
    integer :: is0
    integer :: ie0  
    character(len=maxstrlen) :: model_file
    integer :: is1
    integer :: ie1
    character(len=maxstrlen) :: geostat_file
    integer :: is2
    integer :: ie2
    character(len=maxstrlen) :: smoothing_file
    integer :: is3
    integer :: ie3
    character(len=maxstrlen) :: landuse_file
    integer :: is4
    integer :: ie4
    character(len=maxstrlen) :: walking_data_file
    integer :: is5
    integer :: ie5
    character(len=maxstrlen) :: car_data_file
    integer :: is6
    integer :: ie6
    character(len=maxstrlen) :: airborne_data_file
    integer :: is7
    integer :: ie7
    character(len=maxstrlen) :: output_dir
    integer :: is8
    integer :: ie8  
    character(len=maxstrlen) :: output_file
    integer :: is9
    integer :: ie9   
  end type io_files
end module class_Global
