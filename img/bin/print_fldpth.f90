      program print_fd
! ===============================================
! Resolution
      integer             ::  ix, iy, jx, jy
      integer             ::  nx, ny
      real                ::  west0, north0                 !! map west and north edge
      real                ::  gsize                         !! map grid size
! 
      real                ::  west, east, north, south      !! plot domain

      real,allocatable    ::  lon(:,:),   lat(:,:)
      real                ::  lon1, lon2, lat1, lat2
      real                ::  len, elv, wth

      character*128       ::  params
      parameter              (params='./map/params.txt')
      character*256       ::  ffldpath, flonlat
      character*64        ::  buf
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

      allocate(lon(nx,ny),lat(nx,ny))

      flonlat='./map/lonlat.bin'
      open(11,file=flonlat,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) lon
      read(11,rec=2) lat
      close(11)

      ffldpath='./map/fldpth.txt'

      open(10,file=ffldpath,form='formatted')
      read(10,*) 
 1000 continue
        read(10,*,end=1090) ix, iy, jx, jy, len, elv, wth
        if( ix>0 .and. ix<=nx .and. iy>0 .and. iy<=ny .and. &
            jx>0 .and. jx<=nx .and. jy>0 .and. jy<=ny )then
          lon1=lon(ix,iy)
          lon2=lon(jx,jy)
          lat1=lat(ix,iy)
          lat2=lat(jx,jy)

          if( lon1>=west .and. lon1<=east .and. lat1<=north .and. lat1>=south )then
            write(6,'(4f10.3,3f10.1)') lon1, lat1, lon2, lat2, len, elv, wth
          elseif( lon2>=west .and. lon2<=east .and. lat2<=north .and. lat2>=south )then
            write(6,'(4f10.3,3f10.1)') lon2, lat2, lon1, lat1, len, elv, wth
          endif
        endif
      goto 1000
 1090 continue

      close(10)

      end program print_fd
