program lpara_sfcelv
!=========================================================================
! get weight and gaussian weight and arange it into a easy acess test file;
! contains
! lon   lat    gaussian wgt
! Menaka@IIS 2020/06/02
!=========================================================================
!$ use omp_lib
implicit none
character(len=128)                    :: fname,buf,camadir,outdir
character(len=32)                     :: outname,varname,mapname,inname
integer                               :: syear,eyear
!-map variables
real                                  :: gsize,west, north, east, south ! map boundries
integer                               :: latpx,lonpx,nflp    ! pixel size, calculated
real,allocatable                      :: rivwth(:,:),rivlen(:,:),nextdst(:,:)
integer,allocatable                   :: nextX(:,:),nextY(:,:),ocean(:,:)
integer                               :: i,j,ios,N
real                                  :: lat,lon
real,allocatable                      :: weightage(:,:),gauss_weight(:,:)
integer,allocatable                   :: targetp(:,:),countp(:,:)
real                                  :: lag_dist,conflag
integer                               :: ix,iy,nx,ny
integer                               :: patch_size,patch_side,countnum,patch_nums
integer(kind=4)                       :: i_m,j_m
integer(kind=4)                       :: target_pixel,fn
character(len=8)                      :: llon,llat
real                                  :: threshold
!====================================================
call getarg(1,buf)
read(buf,*) N ! length of time series

call getarg(2,buf)
read(buf,*) syear ! start year

call getarg(3,buf)
read(buf,*) eyear ! end year

call getarg(4,buf)
read(buf,'(A)') outname ! variable name [e.g. sfcelv,rivout,outflw]

call getarg(5,buf)
read(buf,"(A)") mapname ! map name

call getarg(6,buf)
read(buf,'(A)') inname ! input runoff forcing name

call getarg(7,buf)
read(buf,"(A)") camadir
write(*,*) camadir

call getarg(8,buf)
read(buf,"(A)") outdir
write(*,*) outdir

call getarg(9,buf)
read(buf,*) threshold ! threshold for defining the local patch
!-
varname=outname
!==
fname=trim(camadir)//"/map/"//trim(mapname)//"/params.txt"
print *, fname
open(11,file=fname,form='formatted')
read(11,*) lonpx
read(11,*) latpx
read(11,*) nflp
read(11,*) gsize
read(11,*) west
read(11,*) east
read(11,*) south
read(11,*) north
close(11)
!-------
!======
patch_size=1000
patch_side=patch_size*2+1
patch_nums=patch_side**2


fname=trim(adjustl(outdir))//"/local_patch/lonlat.txt"
open(78,file=fname,status='replace')

21 format(i4.4,2x,i4.4,2x,f10.7)
22 format(a4,2x,a4,2x,a8)
23 format(i4.4,2x,i4.4,2x,i4.4,2x,i4.4)

!$ write(*,*)"omp threads",omp_get_num_threads()

!--
allocate(rivwth(lonpx,latpx),rivlen(lonpx,latpx),nextdst(lonpx,latpx))
allocate(nextX(lonpx,latpx),nextY(lonpx,latpx),ocean(lonpx,latpx))
allocate(weightage(lonpx,latpx),gauss_weight(lonpx,latpx),targetp(lonpx,latpx),countp(lonpx,latpx))
! read river width
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivwth.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivwth
    ! ocean is -9999
else
    write(*,*) "no file rivwth"
end if
close(34)

! read next grid information
! read nextX and nextY
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/nextxy.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextX
    read(34,rec=2) nextY
else
    write(*,*) "no file nextXY at:",fname
end if
close(34)

! make ocean mask from nextX data (1 is ocean; 0 is not ocean)
ocean = (nextX==-9999) * (-1)

! read river length
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivlen.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivlen
    ! ocean is -9999
else
    write(*,*) "no file rivlen",fname
end if
close(34)
! read distance to next grid
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/nxtdst.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextdst
else
    write(*,*) "no file nextdst",fname
end if
close(34)
!--
countp=0
targetp=0
!===
! for global map do not write below 60S
south=max(south,-60.0)
nx=lonpx
!ny=latpx-30.0/dble(gsize) ! writed only up -60S latitude
ny=(north-south)/dble(gsize)
!--
! do parallel
!$omp parallel default(none)&
!$omp& private(iy,lat_cent,lat,lon,llon,llat,fname,fn,lag,countnum,i,j,i_m,j_m,threshold)&
!$omp& shared(ocean,rivwth,targetp,countp,patch_size,weightage,gauss_weight)
!$omp do
!--
do ix = 1,nx ! pixels along longtitude direction
    do iy = 1,ny ! pixel along latitude direction
        lat = 90.0-(iy-1.0)/4.0
        lon = (ix-1.0)/4.0-180.0
        !remove ocean
        if (ocean(ix,iy)==1) then
            cycle
        end if
        write(*,*)"===",lat,lon,"===",ix,iy,rivwth(ix,iy)!$,omp_get_num_threads(),omp_get_thread_num(),"==="
        ! open emperical weightage
        write(llon,'(i4.4)') ix
        write(llat,'(i4.4)') iy
        ! read weightage
        fname=trim(adjustl(outdir))//"/weightage/"//trim(llon)//trim(llat)//".bin"
        print*, "read weightage",fname
        fn = 34
        call read_wgt(fname,lonpx,latpx,weightage)
        ! read gausssian weight
        print*, "read gausssian weight"
        fname=trim(adjustl(outdir))//"/gaussian_weight/"//trim(llon)//trim(llat)//".bin"
        fn = 34
        call read_wgt(fname,lonpx,latpx,gauss_weight)
        countnum=1
        ! file to save
        fn=72
        fname=trim(adjustl(outdir))//"/local_patch/patch"//trim(llon)//trim(llat)//".txt"
        open(fn,file=fname,status='replace')
        ! patch size should be >=1000
        do j=iy-patch_size,iy+patch_size
            do i=iy-patch_size,iy+patch_size
                i_m = i
                j_m = j
                if (mapname(1:3)=="glb") then
                    !print*, "global map"
                    call ixy2iixy(i,j,lonpx,latpx,i_m,j_m)
                else
                    !print*, "regional map"
                    i_m=min(i_m,lonpx)
                    i_m=max(i_m,1)
                    j_m=min(j_m,latpx)
                    j_m=max(j_m,1)
                end if
                ! weitage >= threshold is considered
                if (weightage(i_m,j_m) < threshold) then
                    cycle
                end if
                print*, i_m,j_m,weightage(i_m,j_m),threshold
                !write(*,*)i_m,j_m
                ! ocean removed
                if (ocean(i_m,j_m)==1) then
                    cycle
                    !continue
                end if
                ! non river pixels removed
                if (rivwth(i_m,j_m)<=0.0) then
                    cycle
                end if
                ! calculate lag distance
                call lag_distance(i_m,j_m,ix,iy,lonpx,latpx,nextX,nextY,nextdst,lag_dist)
                !---
                write(*,*) lag_dist,i,j,i_m,j_m,ix,iy
                !---
                ! only river pixels which connects to target pixel is considered
                if (lag_dist == -9999.0) then
                    cycle
                    !continue
                end if
                ! consistancy of weigtage
                call wgt_consistancy(i_m,j_m,ix,iy,lonpx,latpx,nextX,nextY,weightage,threshold,conflag)
                if (conflag == 0.0) then
                    cycle
                    !continue
                end if
                !---
                print*, i_m,j_m
                write(fn,21)i_m,j_m,gauss_weight(i_m,j_m)
                write(*,21)i_m,j_m,gauss_weight(i_m,j_m)
                ! find the target pixel
                target_pixel=1
                if (i_m==ix .and. j_m==iy) then
                  target_pixel=countnum
                end if
                countnum=countnum+1
            end do
        end do
        !--
        close(fn)
        countnum=countnum-1
        !$omp critical
        write(78,23) ix,iy,countnum,target_pixel
        write(*,*) ix,iy,countnum,target_pixel,threshold
        !========
        targetp(ix,iy)=target_pixel
        countp(ix,iy)=countnum
    end do
end do
!$omp end do
!$omp end parallel
!---
fname=trim(adjustl(outdir))//"/local_patch/countnum.bin"
open(84,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="replace",iostat=ios)
if(ios==0)then
    write(84,rec=1) countp
    write(84,rec=2) targetp
else
    write(*,*) "no weightage", fname
end if
close(84)
close(78)
deallocate(rivwth,rivlen,nextdst,nextX,nextY,ocean)
deallocate(weightage,gauss_weight,targetp,countp)
end program lpara_sfcelv
!*****************************************************************
subroutine lag_distance(i,j,x,y,nx,ny,nextX,nextY,nextdst,lag_dist)
implicit none 
!--
integer                             :: i,j,x,y,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: nextdst
!--
real                                :: lag_dist
integer                             :: ix,iy,iix,iiy,tx,ty,ud
real                                :: length,rl
!--
if (i==x .and. j==y) then
  ud=0
else
  ud=-1
end if
!--
if (ud==-1) then
  tx=x
  ty=y
  ix=i
  iy=j
  length=0.0
  lag_dist=0.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy 
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix==-9 .or. iy==-9) then
      ud=+1
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=+1
      exit
    end if
    !-- half of the present grid
    rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
    length=length+rl!/2.0
  end do
end if
!--
if (ud==+1) then
  tx=i
  ty=j
  ix=x
  iy=y
  length=0.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    !---
    if (ix==-9 .or. iy==-9) then
      ud=-9999
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=-9999
      exit
    end if
    !-- half of the present grid
    rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
    length=length+rl!/2.0
  end do
end if
!-- 
if (ud==-9999) then
  lag_dist=-9999
elseif (ud==0) then
  lag_dist=0.0
else
  lag_dist=length
end if
!---
return
!---
end subroutine lag_distance
!*****************************************************************
subroutine wgt_consistancy(i,j,x,y,nx,ny,nextX,nextY,weightage,thersold,conflag)
implicit none 
!--
integer                             :: i,j,x,y,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: weightage
!--
real                                :: conflag,thersold
integer                             :: ix,iy,iix,iiy,tx,ty,ud
real                                :: length
!--
if (i==x .and. j==y) then
  ud=0
else
  ud=-1
end if
!--
if (ud==-1) then
  tx=x
  ty=y
  ix=i
  iy=j
  length=0.0
  conflag=1.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy 
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix==-9 .or. iy==-9) then
      ud=+1
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=+1
      exit
    end if
    !--compare the weightage and thersold
    if (weightage(ix,iy) < thersold) then
      conflag=0.0
    end if 
  end do
end if
!--
if (ud==+1) then
  tx=i
  ty=j
  ix=x
  iy=y
  length=0.0
  conflag=1.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    !---
    if (ix==-9 .or. iy==-9) then
      ud=-9999
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=-9999
      exit
    end if
    !--compare the weightage and thersold
    if (weightage(ix,iy) < thersold) then
      conflag=0.0
    end if 
  end do
end if
!-- 
if (ud==-9999) then
  conflag=0.0
elseif (ud==0) then
  conflag=1.0
else
  conflag=conflag
end if
!---
return
!---
end subroutine wgt_consistancy
!**************************************************
subroutine read_wgt(fname,nx,ny,weightage)
!$ use omp_lib    
implicit none
character*128                      :: fname 
integer                            :: fn,nx,ny,ios
real,dimension(nx,ny)              :: weightage
fn=34
!$ fn= fn + omp_get_thread_num()
!!$ write(*,*) fn
open(fn,file=fname,form="unformatted",access="direct",recl=4*ny*nx,status="old",iostat=ios)
if(ios==0)then
    read(fn,rec=1) weightage
else
    write(*,*) "no weightage", fname
end if
close(fn)
!--
return
!---
end subroutine read_wgt
!**************************************************
function Gauss_wt(lag)
implicit none
real                                :: lag,Gauss_wt
real,parameter                      :: sigma=1000.0 !1000 km 
!---
Gauss_wt=exp(-(lag**2.0/(2.0*sigma**2.0)))  
!---
return
!---
end function Gauss_wt
!***************************************************   
function roundx(ix, nx)
implicit none
!-- for input -----------
integer                     ix, nx
!-- for output ----------
integer                     roundx
!------------------------
if (ix .ge. 1) then
  roundx = ix - int((ix -1)/nx)*nx
else
  roundx = nx - abs(mod(ix,nx))
end if 
return
end function roundx
!*****************************************************************
subroutine ixy2iixy(ix,iy, nx, ny, iix, iiy)
implicit none
!- for input -----------------
integer                   ix, iy, nx, ny
!- for output ----------------
integer                   iix, iiy,roundx
!-----------------------------
if (iy .lt. 1) then
  iiy = 2 - iy
  iix = ix + int(nx/2.0)
  iix = roundx(iix, nx)
else if (iy .gt. ny) then
  iiy = 2*ny -iy
  iix = ix + int(nx/2.0)
  iix = roundx(iix, nx)
else
  iiy = iy
  iix = roundx(ix, nx)
end if
return
end subroutine ixy2iixy
!*****************************************************************        
