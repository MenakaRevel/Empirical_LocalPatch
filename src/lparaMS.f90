program lparaMS_sfcelv
!=========================================================================
! select only the mainstream of the empirical local patch
! get weight and gaussian weight and arange it into a easy acess test file;
! contains
! lon   lat    gaussian wgt
! Menaka@IIS 2020/06/02
!=========================================================================
!$ use omp_lib
implicit none
!character*128                   :: fname,buf,camadir,outdir
!character*8                     :: yyyymmdd,nxtyyyymmdd
!real                            :: assimN,assimS,assimW,assimE,lat,lon
!character*2                     :: swot_day
!real,allocatable                :: swot_obs(:,:),global_xa(:,:,:)
!integer*4                       :: lon_cent,lat_cent,patch_size,patch_side,i,j,k,countnum,patch_nums
!integer,parameter               :: latpx=720,lonpx=1440
!real,dimension(lonpx,latpx)     :: rivwth,rivlen,nextdst,rivseq,uparea
!real,dimension(lonpx,latpx)     :: weightage,gauss_weight
!integer,dimension(lonpx,latpx)  :: nextX,nextY,ocean,targetp,countp
!
!real,allocatable                :: Wvec(:),lag(:),local_lag(:)!global_sum_xam(:,:),
!real,allocatable                :: wgt(:),local_wgt(:)
!real                            :: wt,lag_dist,conflag!lag_distance,Gauss_wt,
!integer*4                       :: i_m,j_m,ix,iy
!integer,allocatable             :: xlist(:),ylist(:)
!integer*4                       :: target_pixel,fn
!character*8                     :: llon,llat
!real                            :: threshold
!integer                         :: ios,info
!!for writing
!integer,dimension(lonpx,latpx)  :: upx,upy,lps ! upstrea x,y and length
!
!write(*,*) "lparaMS_sfcelv"
!call getarg(1,buf)
!read(buf,*) patch_size ! radius
!
!call getarg(2,buf)
!read(buf,"(A)") camadir
!write(*,*) camadir
!
!call getarg(3,buf)
!read(buf,"(A)") outdir
!
!call getarg(4,buf)
!read(buf,*) threshold

character(len=128)                    :: fname,buf,camadir,outdir
character(len=32)                     :: outname,varname,mapname,inname
integer                               :: syear,eyear
!-map variables
real                                  :: gsize,west, north, east, south ! map boundries
integer                               :: latpx,lonpx,nflp    ! pixel size, calculated
real,allocatable                      :: rivwth(:,:),rivlen(:,:),nextdst(:,:),uparea(:,:)
integer,allocatable                   :: nextX(:,:),nextY(:,:),rivseq(:,:),ocean(:,:)
integer                               :: i,j,ios,N
real                                  :: lat,lon
real,allocatable                      :: weightage(:,:),gauss_weight(:,:)
integer,allocatable                   :: targetp(:,:),countp(:,:)
real                                  :: lag_dist,conflag
integer                               :: ix,iy,nx,ny
integer                               :: patch_size,patch_side,countnum,patch_nums
integer(kind=4)                       :: i_m,j_m,iix,iiy,k
integer,allocatable                   :: xlist(:),ylist(:)
integer(kind=4)                       :: target_pixel,fn
character(len=8)                      :: llon,llat
real                                  :: threshold
integer                               :: info
!for writing
integer,allocatable                   :: upx(:,:),upy(:,:),lps(:,:) ! upstrea x,y and length
character(len=2)                      :: thrname
!====================================================
call getarg(1,buf)
read(buf,*) N ! length of time series

call getarg(2,buf)
read(buf,*) syear ! start year

call getarg(3,buf)
read(buf,*) eyear ! end year

! call getarg(4,buf)
! read(buf,'(A)') outname ! variable name [e.g. sfcelv,rivout,outflw]

call getarg(4,buf)
read(buf,"(A)") mapname ! map name

call getarg(5,buf)
read(buf,'(A)') inname ! input runoff forcing name

call getarg(6,buf)
read(buf,"(A)") camadir
write(*,*) camadir

call getarg(7,buf)
read(buf,"(A)") outdir
write(*,*) outdir

call getarg(8,buf)
read(buf,*) threshold ! threshold for defining the local patch
!-
! varname=outname
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

print*, threshold*100
write(thrname,'(i2.0)') int(threshold*100)
print*, thrname

fname=trim(adjustl(outdir))//"/local_patchMS/"//trim(mapname)//"_"//trim(inname)//"_"//trim(thrname)//"/lonlat.txt"
open(78,file=fname,status='replace')

21 format(i4.4,2x,i4.4,2x,f10.7)
22 format(a4,2x,a4,2x,a8)
23 format(i4.4,2x,i4.4,2x,i4.4,2x,i4.4)

!$ write(*,*)"omp threads",omp_get_num_threads()
!--
allocate(rivwth(lonpx,latpx),rivlen(lonpx,latpx),nextdst(lonpx,latpx),rivseq(lonpx,latpx),uparea(lonpx,latpx))
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
! read rivseq
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivseq.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivseq
    ! ocean is -9999
else
    write(*,*) "no file rivseq",fname
end if
close(34)
! read uparea
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/uparea.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) uparea
    ! ocean is -9999
else
    write(*,*) "no file uparea",fname
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
ny=(north-south)/gsize !dble(gsize)
!--
allocate(upx(lonpx,latpx),upy(lonpx,latpx),lps(lonpx,latpx))
!===========
! do parallel
!$omp parallel default(none)&
!$omp& private(ix,iy,lat,lon,llon,llat,fname,fn,lag,countnum,i,j,i_m,j_m,thresold)&
!$omp& shared(ocean,rivwth,targetp,countp,patch_size,weightage,gauss_weight)
!$omp do
!--
do ix = 1, nx !int((assimW+180)*4+1),int((assimE+180)*4+1),1
  do iy = 1, ny !int((90-assimN)*4+1),int((90-assimS)*4+1),1
        lat = 90.0-(real(iy)-1.0)/4.0
        lon = (real(ix)-1.0)/4.0-180.0
        !remove ocean
        if (ocean(ix,iy)==1) then
            cycle
            !continue
        end if
        ! ! remove rivwth <= 0m
        ! if (rivwth(ix,iy) <=0.0) then
        !     cycle
        ! end if
        write(*,*)"===",lat,lon,"===",ix,iy,rivwth(ix,iy)
        ! open emperical weightage
        write(llon,'(i4.4)') ix
        write(llat,'(i4.4)') iy
        fname=trim(adjustl(outdir))//"/weightage/"//trim(mapname)//"_"//trim(inname)//"/"//trim(llon)//trim(llat)//".bin"
        fn = 34
        call read_wgt(fname,lonpx,latpx,weightage)
        ! read gausssian weight
        fname=trim(adjustl(outdir))//"/gaussian_weight/"//trim(mapname)//"_"//trim(inname)//"_"//trim(thrname)//"/"//trim(llon)//trim(llat)//".bin"
        fn = 34
        call read_wgt(fname,lonpx,latpx,gauss_weight)
        ! get the mainstream pixels
        ! most upstream pixels and number of downstream pixels
        iix=0
        iiy=0
        k=0
        call patch_pixels(ix,iy,lonpx,latpx,threshold,weightage,nextX,nextY,nextdst,uparea,rivseq,iix,iiy,k)
        upx(ix,iy)=iix
        upy(ix,iy)=iiy
        lps(ix,iy)=k
        write(*,*)"===",ix,iy,"===",iix,iiy,k
        allocate(xlist(k),ylist(k))
        call downstream_pixels(iix,iiy,k,lonpx,latpx,weightage,nextX,nextY,xlist,ylist,info)
        if (info /= 0) then
            write(*,*) "info:",info
            deallocate(xlist,ylist)
            cycle
        end if
        ! file to save
        fn=72
        fname=trim(adjustl(outdir))//"/local_patchMS/"//trim(mapname)//"_"//trim(inname)//"_"//trim(thrname)//"/patch"//trim(llon)//trim(llat)//".txt"
        open(fn,file=fname,status='replace')
        !write(fn,22)"lon","lat","","forcast","assim"
        do i=1,k
            i_m=xlist(i)
            j_m=ylist(i)
            write(fn,21)i_m,j_m,gauss_weight(i_m,j_m)
            write(*,21)i_m,j_m,gauss_weight(i_m,j_m)
            ! find the target pixel
            if (i_m==ix .and. j_m==iy) then
                target_pixel=i
            end if
        end do
        !-------
        close(fn)
        countnum=k
        !write(*,*) "size",countnum 
        !$omp critical
        write(78,23) ix,iy,countnum,target_pixel
        write(*,*) ix,iy,countnum,target_pixel,threshold
        write(*,*) "************************************************"
!        !$omp end critical
        targetp(ix,iy)=target_pixel
        countp(ix,iy)=countnum
        !--
        deallocate(xlist,ylist)
    end do
end do
!$omp end do
!$omp end parallel
!---
fname=trim(adjustl(outdir))//"/local_patchMS/"//trim(mapname)//"_"//trim(inname)//"_"//trim(thrname)//"/countnum.bin"
open(84,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="replace",iostat=ios)
if(ios==0)then
    write(84,rec=1) countp
    write(84,rec=2) targetp
else
    write(*,*) "no countnum.bin", fname
end if
close(84)
!--
!---
fname=trim(adjustl(outdir))//"/local_patchMS/"//trim(mapname)//"_"//trim(inname)//"_"//trim(thrname)//"/lpara_patch.bin"
open(85,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="replace",iostat=ios)
if(ios==0)then
    write(85,rec=1) upx
    write(85,rec=2) upy
    write(85,rec=3) lps
else
    write(*,*) "no lpara_patch.bin", fname
end if
close(85)

close(78)

deallocate(rivwth,rivlen,nextdst,nextX,nextY,rivseq,uparea,ocean)
deallocate(weightage,gauss_weight,targetp,countp,upx,upy,lps)
end program lparaMS_sfcelv
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
    if (ix==-10 .or. iy==-10) then
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
    if (ix==-10 .or. iy==-10) then
      ud=+1
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
subroutine upstream(i,j,nx,ny,nextX,nextY,rivseq,uparea,x,y)
implicit none 
! find the upstream pixel with closest upstream area to the i,j
!--
integer                             :: i,j,nx,ny,x,y
integer,dimension(nx,ny)            :: nextX,nextY,rivseq
real,dimension(nx,ny)               :: nextdst, uparea
!--
real                                :: dA ! area differnce
integer                             :: ix,iy,iix,iiy,tx,ty,ud,d
real                                :: length,rl
!--
x=-9999
y=-9999
d=10 ! look at 100*100 box
dA=1.0e20 ! area differnce
!--
!write(*,*)i,j
do tx=i-d,i+d
  do ty=j-d,j+d
    !write(*,*)tx,ty
    call ixy2iixy(tx,ty, nx, ny, ix, iy)
    !write(*,*)nextX(ix,iy),nextY(ix,iy),ix,iy,uparea(ix,iy),rivseq(ix,iy)
    if (nextX(ix,iy) == i .and. nextY(ix,iy) == j) then
        !write(*,*)ix,iy
        if (abs(uparea(i,j)-uparea(ix,iy)) < dA) then
            dA=abs(uparea(i,j)-uparea(ix,iy))
            x=ix
            y=iy
            !write(*,*)x,y
        end if
    end if
  end do
end do
return
end subroutine upstream
!******************************************************************
subroutine patch_pixels(i,j,nx,ny,threshold,weight,nextX,nextY,nextdst,uparea,rivseq,x,y,k)
implicit none
! find the upstream pixel with closest upstream area to the i,j
!--
integer                             :: i,j,x,y,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY,rivseq
real,dimension(nx,ny)               :: weight,nextdst,uparea
!--
real                                :: threshold,dA ! area differnce
integer                             :: ix,iy,iix,iiy,tx,ty,ud,d,k
real                                :: length,rl
!integer,dimension(1000)             :: lx,ly
!--
k=0
!write(*,*)i,j
x=i
y=j
if (rivseq(i,j)>1) then ! no upstream
    ix=i
    iy=j
    do while (weight(ix,iy)>=threshold)
        call upstream(ix,iy,nx,ny,nextX,nextY,rivseq,uparea,iix,iiy)
        !write(*,*)iix,iiy,rivseq(ix,iy)
        if (weight(iix,iiy)<threshold) exit
        !write(*,*)iix,iiy,weight(iix,iiy),weight(iix,iiy)>=threshold
        if (rivseq(iix,iiy)==1) then
            x=iix
            y=iiy
            exit
        end if
        ix=iix
        iy=iiy
        k=k+1
    end do
    x=ix
    y=iy
end if
write(*,*)"upstream:",k!,weight(x,y)
ud=k
k=0
call downstream(x,y,nx,ny,threshold,weight,nextX,nextY,k)
write(*,*)"downstream:",k-ud !,max(k-ud,0)
return
end subroutine patch_pixels
!*****************************************************************
subroutine downstream(i,j,nx,ny,threshold,weight,nextX,nextY,k)
implicit none
!--
integer                             :: i,j,x,y,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: weight !nextdst
!--
real                                :: threshold,lag_dist
integer                             :: ix,iy,iix,iiy,tx,ty,ud,k
real                                :: length,rl
!--
k=1
!--
ix=i
iy=j
do while (weight(ix,iy)>=threshold)
  iix=ix
  iiy=iy 
  ix=nextX(iix,iiy)
  iy=nextY(iix,iiy)
  !write(*,*)ix,iy
  if (ix==-9 .or. iy==-9) then
    k=k+1
    exit
  end if
  if (ix==-10 .or. iy==-10) then
    k=k+1
    exit
  end if
  if (ix==-9999 .or. iy==-9999) then
    k=k+1
    exit
  end if
  k=k+1
end do
!---
k=k-1
if (k==0) then
    if (weight(i,j)>=threshold) then
        k=1
    end if
end if
!write(*,*)k,weight(i,j)
!---
return
!---
end subroutine downstream
!*****************************************************************
subroutine downstream_pixels(i,j,k,nx,ny,weight,nextX,nextY,lx,ly,info)
implicit none
!--
integer                             :: i,j,k,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: weight !nextdst
!--
real                                :: threshold,lag_dist
integer                             :: num,iix,iiy,ix,iy,info
real                                :: length,rl
integer,dimension(k)                :: lx,ly
!--
num=1
ix=i
iy=j
lx(num)=ix
ly(num)=iy
num=num+1
!--
do while (num<=k)
  iix=ix
  iiy=iy
  ix=nextX(iix,iiy)
  iy=nextY(iix,iiy)
  ! write(*,*)ix,iy
  if (ix==-9 .or. iy==-9) then
    num=num+1
    exit
  end if
  if (ix==-10 .or. iy==-10) then
    num=num+1
    exit
  end if
  if (ix==-9999 .or. iy==-9999) then
    num=num+1
    exit
  end if
  lx(num)=ix
  ly(num)=iy
  num=num+1
end do
!---
num=num-1
info=num-k
!---
return
!---
end subroutine downstream_pixels
!*****************************************************************
