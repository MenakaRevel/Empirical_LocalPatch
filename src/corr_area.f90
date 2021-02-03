program corr_area 
!$ use omp_lib    
implicit none
character*128                   :: fname,buf,camadir,outdir
character*8                     :: yyyymmdd,nxtyyyymmdd
real                            :: assimN,assimS,assimW,assimE,lat,lon
character*2                     :: swot_day
real,allocatable                :: swot_obs(:,:),global_xa(:,:,:)
integer*4                       :: lon_cent,lat_cent,patch_size,patch_side,i,j,k,countnum,patch_nums
integer,parameter               :: latpx=720,lonpx=1440
real,dimension(lonpx,latpx)     :: rivwth,rivlen,nextdst,weightage,gauss_weight,grarea,area
integer,dimension(lonpx,latpx)  :: nextX,nextY,ocean,targetp,countp,pixel

real,allocatable                :: Wvec(:),lag(:),local_lag(:)!global_sum_xam(:,:),
real,allocatable                :: wgt(:),local_wgt(:)
real                            :: wt,lag_dist!lag_distance,Gauss_wt,
integer*4                       :: i_m,j_m
integer,allocatable             :: xlist(:),ylist(:)
integer*4                       :: target_pixel,fn
character*8                     :: llon,llat
real                            :: thresold
integer                         :: ios

write(*,*) "corr_area"
call getarg(1,buf)
read(buf,*) patch_size ! radius

call getarg(2,buf)
read(buf,"(A)") camadir
write(*,*) camadir

call getarg(3,buf)
read(buf,"(A)") outdir

call getarg(4,buf)
read(buf,*) thresold


patch_side=patch_size*2+1
patch_nums=patch_side**2


!fname=trim(adjustl(outdir))//"local_patch/lonlat.txt"
!open(78,file=fname,status='replace')

21 format(i4.4,2x,i4.4,2x,f10.7)
22 format(a4,2x,a4,2x,a8)
23 format(i4.4,2x,i4.4,2x,i4.4,2x,i4.4)

!$ write(*,*)"omp threads",omp_get_num_threads()

! read river width
fname=trim(adjustl(camadir))//"map/glb_15min/rivwth.bin"
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
fname=trim(adjustl(camadir))//"map/glb_15min/nextxy.bin"
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
fname=trim(adjustl(camadir))//"map/glb_15min/rivlen.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivlen
    ! ocean is -9999
else
    write(*,*) "no file rivlen",fname
end if
close(34)

! read distance to next grid
fname=trim(adjustl(camadir))//"map/glb_15min/nxtdst.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextdst
else
    write(*,*) "no file nextdst",fname
end if
close(34)

! read distance to uparea
fname=trim(adjustl(camadir))//"map/glb_15min/grdare.bin"
write(*,*)"grarea", fname
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) grarea
else
    write(*,*) "no file grarea",fname
end if
close(34)
!--
countp=0
targetp=0
area=0
pixel=0
!--
! do parallel
!$omp parallel default(none)&
!$omp& private(lon_cent,lat_cent,lat,lon,llon,llat,fname,fn,lag,countnum,i,j,i_m,j_m,thresold)
!$omp& shared(ocean,rivwth,targetp,countp,patch_size,weightage,gauss_weight)
!$omp do
!--
do lon_cent = 1,lonpx !int((assimW+180)*4+1),int((assimE+180)*4+1),1
  do lat_cent = 1, latpx !int((90-assimN)*4+1),int((90-assimS)*4+1),1
        lat = 90.0-(lat_cent-1.0)/4.0
        lon = (lon_cent-1.0)/4.0-180.0
    
        !remove ocean
        if (ocean(lon_cent,lat_cent)==1) then
            cycle
            !continue
        end if  
        ! remove rivwth <= 0m
        if (rivwth(lon_cent,lat_cent) <=0.0) then
            cycle
        end if

        ! open emperical weightage
        write(llon,'(i4.4)') lon_cent
        write(llat,'(i4.4)') lat_cent
        fname=trim(adjustl(outdir))//"weightage/"//trim(llon)//trim(llat)//".bin"
        fn = 34
        call read_wgt(fname,lonpx,latpx,weightage)
        !-----------
        gauss_weight=(weightage(lon_cent,lat_cent)>=thresold)*(-1.0)
        area(lon_cent,lat_cent)=sum(((weightage>=thresold)*(-1.0))*grarea) ! km^2
        pixel(lon_cent,lat_cent)=sum(((weightage>=thresold)*(-1.0))) ! number of pixels 
        !write(*,*)lat,lon,sum(((weightage>=0.6)*(-1.0))*uparea*(10**(-6))),sum(((weightage>=0.6)*(-1.0)))
        write(*,*)lat,lon,area(lon_cent,lat_cent),pixel(lon_cent,lat_cent)
        !$omp critical
        !--
        !deallocate(lag,xlist,ylist,wgt)
    end do
end do
!$omp end do
!$omp end parallel
!---
fname=trim(adjustl(outdir))//"Corr_Area/area.bin"
open(84,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="replace",iostat=ios)
if(ios==0)then
    write(84,rec=1) area
else
    write(*,*) "no area", fname
end if
close(84)
!---
fname=trim(adjustl(outdir))//"Corr_Area/pixel.bin"
open(84,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="replace",iostat=ios)
if(ios==0)then
    write(84,rec=1) pixel
else
    write(*,*) "no pixel", fname
end if
close(84)

!close(78)



end program corr_area
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
