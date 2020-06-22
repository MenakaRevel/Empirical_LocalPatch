      program txt_vector
! ================================================
      implicit none
! Resolution
      integer             ::  ix, iy, jx, jy
      integer             ::  nx, ny                        !! grid numbers
      real                ::  west0, north0, east0, south0  !! map edge
      real                ::  gsize                         !! map grid size
      integer             ::  flood                         !! floodpalin layer
! domain
      real                ::  west, east, north, south      !! plot domain
! vars
      integer,allocatable ::  nextx(:,:), nexty(:,:)        !! downstream xy
      real,allocatable    ::  lon(:,:), lat(:,:)
      real,allocatable    ::  uparea(:,:)                   !! dorainage area [m2]
      real                ::  lon1, lon2, lat1, lat2
! files
      character*540       ::  params
      character*240       ::  camadir
      character*240       ::  mapname
      character*240       ::  rfile1, rfile2, rfile3
!      parameter              (params='./map/params.txt')
!      parameter              (rfile1='./map/nextxy.bin')
!      parameter              (rfile2='./map/lonlat.bin')
!      parameter              (rfile3='./map/uparea.bin')

      character*64        ::  buf
! ================================================
      call getarg(1,buf)
      read(buf,*) west
      call getarg(2,buf)
      read(buf,*) east
      call getarg(3,buf)
      read(buf,*) north
      call getarg(4,buf)
      read(buf,*) south
      call getarg(5,camadir)
      call getarg(6,mapname)
      !read(buf,"(A)") camadir
      !write(*,*) south, camadir
      ! updated accoring to params.txt of CaMa v396
      ! updated to add map floder name
      ! Menaka@IIS 2019/10/26
      if (mapname==" ") mapname="glb_06min"
      params=trim(camadir)//"/map/"//trim(mapname)//"/params.txt"
      print* , trim(params)
      open(10,file=trim(params),form='formatted')
      read(10,*) nx
      read(10,*) ny
      read(10,*) gsize
      read(10,*) flood
      read(10,*) west0
      read(10,*) east0
      read(10,*) south0
      read(10,*) north0
      
      close(10)
      
      allocate(nextx(nx,ny),nexty(nx,ny))
      allocate(lon(nx,ny),lat(nx,ny))
      allocate(uparea(nx,ny))

      rfile1=trim(camadir)//"/map/"//trim(mapname)//"/nextxy.bin"
      open(11,file=trim(rfile1),form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) nextx
      read(11,rec=2) nexty
      close(11)

      rfile2=trim(camadir)//"/map/"//trim(mapname)//"/lonlat.bin"
      open(11,file=trim(rfile2),form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) lon
      read(11,rec=2) lat
      close(11)

      rfile3=trim(camadir)//"/map/"//trim(mapname)//"/uparea.bin"
      open(11,file=trim(rfile3),form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) uparea
      close(11)

      do iy=1, ny
        do ix=1, nx
          if( nextx(ix,iy)>0 )then
            jx=nextx(ix,iy)
            jy=nexty(ix,iy)
            lon1=lon(ix,iy)
            lon2=lon(jx,jy)
            lat1=lat(ix,iy)
            lat2=lat(jx,jy)

            if( lon1>=west .and. lon1<=east .and. lat1<=north .and. lat1>=south )then
              print '(4f12.5,f12.1)', lon1, lat1, lon2, lat2, uparea(ix,iy)/1000.**2
            elseif( lon2>=west .and. lon2<=east .and. lat2<=north .and. lat2>=south )then
              print '(4f12.5,f12.1)', lon2, lat2, lon1, lat1, uparea(ix,iy)/1000.**2
            endif
          elseif( nextx(ix,iy)==-9 .or. nextx(ix,iy)==-10 )then
            print '(4f12.5,f12.1)', lon(ix,iy), lat(ix,iy), -999.,      -999.,      &
                                                   uparea(ix,iy)/1000.**2
          endif
        end do
      end do
      
! === Finish =====================================
      end program txt_vector
