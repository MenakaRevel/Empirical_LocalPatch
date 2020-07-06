      program print_grid
! ===============================================
      implicit none
! index
      integer            ::  ix, iy
      integer            ::  nx, ny                        !! grid numbers
      real               ::  west0, north0                 !! map west and north edge
      real               ::  gsize                         !! map grid size
! vars
      real               ::  west, east, north, south      !! plot domain

      real,allocatable   ::  lon(:), lat(:)                !! lon lat
      integer,allocatable::  nextx(:,:)                    !! downstream x
! file
      character*128      ::  params
      character*128      ::  rfile1
      parameter             (params='./map/params.txt')
      parameter             (rfile1='./map/nextxy.bin')

      character*64       ::  buf
! ===============================================
      call getarg(1,buf)
      read(buf,*) west
      call getarg(2,buf)
      read(buf,*) east
      call getarg(3,buf)
      read(buf,*) north
      call getarg(4,buf)
      read(buf,*) south

      open(10,file=params,form='formatted')
      read(10,*) west0
      read(10,*) north0
      read(10,*) nx
      read(10,*) ny
      read(10,*) gsize
      close(10)

      allocate(lon(nx),lat(ny),nextx(nx,ny))

      open(11, file=rfile1, form='unformatted', access='direct', recl=4*nx*ny)
      read(11,rec=1) nextx
      close(11)
      
      do ix=1, nx
        lon(ix)=west0 + (real(ix)-1)*gsize
      end do
      do iy=1, ny
        lat(iy)=north0- (real(iy)-1)*gsize
      end do

! ===============================================
      do iy=1, ny
        do ix=1, nx
          if( nextx(ix,iy)/=-9999 )then 
            if( lon(ix)>=west .and. lon(ix)<=east .and. lat(iy)>=south .and. lat(iy)<=north )then
              print *, lon(ix), lat(iy), 1
            endif
          endif
        end do
      end do

      end program print_grid
