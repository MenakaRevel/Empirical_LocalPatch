program semivariance
!===================================================
! calculates the semivariances along the river
! modeifed to use rivseq
! Menaka@IIS
! files read  using netCDF4 format
! experimental semivarinces are strored in text files
! 2020/05/30
! Reference:
! Revel et al,. (2018), Revel et al,. (2019)
!====================================================
!$ use omp_lib
use netcdf
implicit none
character(len=128)                    :: fname,buf,camadir,outdir,tag
character(len=32)                     :: outname,varname,mapname,inname
integer                               :: syear,eyear
!-map variables
real                                  :: gsize,west, north, east, south ! map boundries
integer                               :: latpx,lonpx,nflp    ! pixel size, calculated
real,allocatable                      :: globalarray(:,:,:),globaltrue(:),offset(:)
real,allocatable                      :: nextdst(:,:)
integer,allocatable                   :: rivseq(:,:),nextX(:,:),nextY(:,:),ocean(:,:)
integer                               :: i,j,ios,N
integer                               :: ncidin,varidin
integer                               :: ix,iy,nx,ny
integer,dimension(3)                  :: start,count
real,parameter                        :: rmis=1.0e20
real,allocatable                      :: xf(:)
real                                  :: semivar,std
real,allocatable                      :: rlen(:)
integer                               :: seqnum,patch_nums,countnum
integer,allocatable                   :: ux(:),uy(:),xt(:),yt(:),xlist(:),ylist(:)
integer                               :: un,k
character(len=8)                      :: lon,lat,u
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
print*,"allocate"
allocate(nextX(lonpx,latpx),nextY(lonpx,latpx),ocean(lonpx,latpx),rivseq(lonpx,latpx),nextdst(lonpx,latpx))
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
! make ocean mask from storage data (1 is ocean; 0 is not ocean)
ocean = (nextX == -9999) * (-1)

! read nextdst
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/nxtdst.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextdst
    ! ocean is -9999
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
!--

patch_nums=100000
allocate(xt(patch_nums),yt(patch_nums),xf(patch_nums))
allocate(ux(patch_nums),uy(patch_nums),rlen(patch_nums))

fname=trim(adjustl(outdir))//"/semivar/"//trim(mapname)//"_"//trim(inname)//"/lonlat_list.txt"
open(79,file=fname,form="formatted",status='replace',iostat=ios)
write(79,'(a4,4x,a4,4x,a3,4x,a3)')"lon","lat","up","dn"

21  format(i4.4,2x,i4.4,2x,f8.2,2x,e14.5,2x,e14.5)
22  format(a4,2x,a4,2x,a8,2x,a14,2x,a14)
23  format(a4,4x,a4,4x,i5.5,4x,i5.5)
!!$write(*,*) omp_get_num_threads()
! for global map do not write below 60S
south=max(south,-60.0)
nx=lonpx
!ny=latpx-30.0/dble(gsize) ! writed only up -60S latitude
ny=(north-south)/gsize !dble(gsize)
print*,nx,ny
! find all the upstream , rivseq=1
allocate(xlist(1000000),ylist(1000000))
xlist=-9999
ylist=-9999
seqnum=0
i=1
!print*,"rivseq"
do ix = 1,nx !
    do iy = 1,ny
        if (rivseq(ix,iy) == 1) then
            !print*,ix,iy
            xlist(i)=ix
            ylist(i)=iy
            i=i+1
        end if
    end do
end do
seqnum=i-1
write(*,*) "seqnum",seqnum
write(*,*) "start calculation"
!---------------------
write(tag,'(i4.0,a,i4.0)')syear,"-",eyear
! allocate varibles
!=====================
allocate(globalarray(nx,ny,N),globaltrue(N),offset(N))
!===============
!*********************************
!============READING==============
!*********************************
!  ncidin
!  varidin
!---------------------------------
! read netCDF file
fname=trim(outdir)//"/CaMa_out/"//trim(mapname)//"_"//trim(inname)//"/"//trim(varname)//trim(tag)//".nc"
print*, "open ",trim(fname)
call nccheck( nf90_open(trim(fname), nf90_nowrite, ncidin) )
!call nccheck( nf90_open_par(fname, IOR(NF90_NETCDF4,NF90_MPIIO), MPI_COMM_WORLD, MPI_INFO_NULL, ncidin) )
print*, "inqure variable id", " ", trim(varname)
!call nccheck( nf90_inq_varid(ncidin, trim(varname),varidin) )
call nccheck( nf90_inq_varid(ncidin, "standardize",varidin) )
!===read the traget pixel==
start=(/1,1,1/)
count=(/nx,ny,N/)
print*, "read variable"
! get variable
call nccheck( nf90_get_var(ncidin,varidin,globalarray,start=start,count=count) )
!====close netcdf=====
call nccheck( nf90_close(ncidin ) )
! parallel calculation

!$omp parallel default(private) shared(ocean,rivwth,globalarray,nextX,nextY,outdir,patch_nums,nextdst,N,ix)
!!!default(shared) private(lat_cent,xt,xf,i,j,i_m,j_m,rlen,pixel, cov,corr,countnum)
!$omp do
do ix = 1,nx !
    do iy = 1,ny
        !remove ocean
        if (ocean(ix,iy) == 0) then
            countnum=1
            write(lon,'(i4.4)')ix
            write(lat,'(i4.4)')iy
            !!!===read the traget pixel==
            !!start=(/ix,iy,1/)
            !!count=(/1,1,N/)
            !!! get variable subset
            !!call nccheck( nf90_get_var(ncidin,varidin,globaltrue,start=start,count=count) )
            !===========
            globaltrue=globalarray(ix,iy,:)
            !===========
            write(*,*)"================================="
            write(*,*) lon,lat!,results
            write(*,*)"*********************************"
            !***********************************************
            !===============upstram pixels =================
            !===calculate number of pixels where rivseq=1===
            call up_river_pixel(ix,iy,xlist,ylist,lonpx,latpx,seqnum,patch_nums,nextX,nextY,ux,uy,un)
            write(*,*) "un",un!ux,uy
            do i=1,un
                !===wiriting file===
                write(u,'(a2, i5.5)') "up",i
                fname=trim(adjustl(outdir))//"/semivar/"//trim(mapname)//"_"//trim(inname)//"/"//trim(lon)//trim(lat)//"/"//trim(u)//".svg"
                write(*,*)fname
                open(34,file=fname,form="formatted",status='replace',iostat=ios)
                if (ios /= 0) then
                    write(*,*) "***cannot make file***"
                endif
                write(34,22)"lon","lat","dis","gamma","sig"
                write(*,*) "upstream pixel:",ux(i),uy(i)
                xt=-9999
                yt=-9999
                !===got through each upstream===
                call river_up(ux(i),uy(i),ix,iy,patch_nums,lonpx,latpx,nextX,nextY,nextdst,xt,yt,rlen,k)
                write(*,*)"-------------",k,"upstreams pixels"!,xt,yt
                do j=1,k
                    !!!!===offset pixel===
                    !!!start=(/xt(j),yt(j),1/)
                    !!!count=(/1,1,N/)
                    !!!! get variable subset
                    !!!call nccheck( nf90_get_var(ncidin,varidin,offset,start=start,count=count) )
                    !===========
                    offset=globalarray(xt(j),yt(j),:)
                    !===========
                    ! calculate experimental semivariance
                    call semi_var(globaltrue,offset,N,semivar,std)
                    write(34,21)xt(j),yt(j),rlen(j),semivar,std
                    write(*,21)xt(j),yt(j),rlen(j),semivar,std
                end do
                close(34)
            end do
            !***********************************************
            !===downstream pixels===
            xt=-9999
            yt=-9999
            call river_dn(ix,iy,patch_nums,lonpx,latpx,nextX,nextY,nextdst,xt,yt,rlen,k) 
            write(*,*)"-------------",k,"downstream pixels"
            write(u,'(a2, i5.5)') "dn",0
            fname=trim(adjustl(outdir))//"/semivar/"//trim(mapname)//"_"//trim(inname)//"/"//trim(lon)//trim(lat)//"/"//trim(u)//".svg"
            write(*,*)fname
            open(34,file=fname,form="formatted",status='replace',iostat=ios)
            if (ios /= 0) then
                write(*,*) "cannot make file"
            endif
            write(34,22)"lon","lat","dis","gamma","sig"
            write(*,*) "downstream pixel:" ,xt(k),yt(k)
            do j=1,k
                !!!start=(/xt(j),yt(j),1/)
                !!!count=(/1,1,N/)
                !!!! get variable subset
                !!!call nccheck( nf90_get_var(ncidin,varidin,offset,start=start,count=count) )
                !===========
                offset=globalarray(xt(j),yt(j),:)
                !===========
                ! calculate experimantal semivariance
                call semi_var(globaltrue,offset,N,semivar,std)
                write(34,21)xt(j),yt(j),rlen(j),semivar,std
                write(*,21)xt(j),yt(j),rlen(j),semivar,std
            end do
            close(34)
            if (k>0) then 
                k=1
            else
                k=0
            end if
            write(79,23)lon,lat,un,k
            write(*,*)"##################################"
        end if
     end do
end do
!$omp end do
!$omp end parallel
!====close netcdf=====
!call nccheck( nf90_close(ncidin ) )

!---
deallocate(nextX,nextY,ocean,rivseq,nextdst)
deallocate(globalarray,globaltrue,offset)
deallocate(xt,yt,xf,ux,uy,rlen)
deallocate(xlist,ylist)
close(79)
end program semivariance
!*****************************************************************
function roundx(ix, nx)
implicit none
!-- for input -----------
integer                     :: ix, nx
!-- for output ----------
integer                     :: roundx
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
integer                   :: ix, iy, nx, ny
!- for output ----------------
integer                   :: iix, iiy,roundx
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
subroutine uord(i,j,x,y,nx,ny,nextX,nextY,ud)
implicit none 
!--
! find the (i,j) upstream/downstream of (x,y) 
! upstream ud = -1  downstream ud = +1
! ud = -9999 : another river 
integer                     :: i,j,x,y,nx,ny
integer,dimension(nx,ny)    :: nextX,nextY
!--
integer                     :: ix,iy,iix,iiy,tx,ty,ud
!--
tx=x
ty=y
ix=i
iy=j
ud=-1
do while (ix/=tx .or. iy/=ty) 
  iix=ix
  iiy=iy
  ix=nextX(iix,iiy)
  iy=nextY(iix,iiy)
  if (ix==-9 .or. iy==-9) then
    ud=+1 
    exit
  end if
  if (ix == -10 .or. iy == -10) then  ! inland termination
    ud=+1 
    exit
  end if
end do
!---
if (ud==+1) then
  !-
  tx=i
  ty=j
  ix=x
  iy=y
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix ==-9 .or. iy==-9) then
      ud=-9999
      exit
    end if 
    if (ix == -10 .or. iy == -10) then ! inland termination
      ud=-9999
      exit
    end if 
  end do
end if
return
!---
end subroutine uord
!*****************************************************************
subroutine river_up(i,j,x,y,patch,nx,ny,nextX,nextY,nextdst,xpixel,ypixel,rlen,k)
implicit none
!--
integer                             :: i,j,x,y,patch,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: nextdst
!--
integer,dimension(patch)            :: xpixel,ypixel,lx,ly
real,dimension(patch)               :: rlen,rlen1
integer                             :: ix,iy,iix,iiy,tx,ty,k,l
real                                :: length,rl
!--
xpixel = -9999
ypixel = -9999
!---

rlen=-9999.0
rlen1=-9999.0
rl=0.0
length=0.0
!--
xpixel(1) =i
ypixel(1) =j
!---
rl=anint((nextdst(i,j)/1000.0)*100)/100.0
rlen1(1)=0.0
!-
k=2
!-
tx=x
ty=y
ix=i
iy=j
do while (ix/=tx .or. iy/=ty) 
  iix=ix
  iiy=iy
  ix=nextX(iix,iiy)
  iy=nextY(iix,iiy)
  if (ix==-9 .or. iy==-9) then
    exit
  end if
  if (ix==-9999 .or. iy==-9999) then
    exit
  end if  
  xpixel(k) =ix
  ypixel(k) =iy 
  rl=anint((nextdst(iix,iiy)/1000.0)*100)/100.0
  length=length+rl
  !write(*,*)"+++",ix,iy,length
  rlen1(k)=length
  k=k+1
end do
!---- 
k=k-1
rl=0.0
do l=1,k 
  !write(*,*)rl!en(k-l+1),k,l
  rlen(l)=rlen1(k)-rlen1(k-l+1)
  lx(l)=xpixel(k-l+1)
  ly(l)=ypixel(k-l+1)
  end do 
xpixel=lx 
ypixel=ly
!rlen=rlen1
!write(*,*)"--",xpixel(1),ypixel(1),rlen(1)
return
!---
end subroutine river_up
!*****************************************************************
subroutine river_dn(x,y,patch,nx,ny,nextX,nextY,nextdst,xpixel,ypixel,rlen,k)
implicit none 
!--
integer                             :: x,y,patch,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: nextdst
!--
integer,dimension(patch)            :: xpixel,ypixel
real,dimension(patch)               :: rlen
integer                             :: ix,iy,iix,iiy,tx,ty,k
real                                :: length,rl
!--
xpixel = -9999
ypixel = -9999
!---
xpixel(1) =x
ypixel(1) =y
rlen= 0.0
rl=0.0
length=0.0
!--
rl=anint((nextdst(x,y)/1000.0)*100)/100.0
rlen(1)=0.0
!-
k=2
ix=x
iy=y
do while (ix>0 .or. iy>0) 
  iix=ix
  iiy=iy
  !write(*,*)iix,iiy
  ix=nextX(iix,iiy)
  iy=nextY(iix,iiy)
  if (ix==-9 .or. iy==-9) then
    exit
  end if
  if (ix==-10 .or. iy==-10) then ! inland termination
    exit
  end if
  if (ix==-9999 .or. iy==-9999) then
    exit
  end if 
  xpixel(k) =ix
  ypixel(k) =iy 
  !--next distance => nextdst 
  rl=anint((nextdst(iix,iiy)/1000.0)*100)/100.0
  length=length+rl!/2.0
  rlen(k)=length
  k=k+1
  if (k>patch) then
    exit
  end if  
end do
!-- 
k=k-1
return
!---
end subroutine river_dn
!***************************************************************
subroutine up_river_pixel(x,y,xlist,ylist,nx,ny,seqnum,patch,nextX,nextY,ux,uy,un)
implicit none
!--
integer                             :: x,y,nx,ny,seqnum,patch
integer,dimension(seqnum)           :: xlist,ylist
integer,dimension(nx,ny)            :: nextX,nextY
integer,dimension(patch)            :: ux,uy
!--
integer                             :: ix,iy,k,un,ud
!--
un=1
do k=1, seqnum
    ix=xlist(k)
    iy=ylist(k)
    !-- find up or down of x,y
    call uord(ix,iy,x,y,nx,ny,nextX,nextY,ud)
    if (ud/=-1) then
        cycle
    end if
    ux(un)=ix
    uy(un)=iy
    un=un+1
end do
!--
un=un-1
return
!---
end subroutine up_river_pixel
!*****************************************************************
subroutine semi_var(t,h,n,semivar,std)
implicit none
!--
integer                   :: n,i
real,dimension(n)         :: t,h
!--
real                      :: semivar,std
!-
real                      :: p,v
!------------------------------------
p=0.0
v=0.0
semivar=0.0
do i=1,n
  p=p+(h(i)-t(i))**2.0
  v=v+(((h(i)-t(i))**2.0)/2.0)**2.0
end do
semivar=p/(2.0*real(n))
std=sqrt(abs((v/(real(n-1)+1.0e-20))-((semivar**2.0)*(real(n)/(real(n-1)+1.0e-20)))))
return
end subroutine semi_var
!******************************
subroutine nccheck(status)
use netcdf
implicit none
integer, intent ( in)          :: status
!=============================
if(status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  stop "Stopped"
end if
end subroutine nccheck
!***************************************************