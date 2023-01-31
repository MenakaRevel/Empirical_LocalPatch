program standardized
!===================================================
! standardized the trend and seasonality removed
! (data/mean) - std
! Menaka@IIS
! files read and write using netCDF4 format
! 2020/05/30
! Reference:
! Revel et al,. (2018), Revel et al,. (2019)
!====================================================
!$ use omp_lib
use netcdf
implicit none
character(len=128)                    :: fname,buf,camadir,outdir
character(len=32)                     :: outname,varname,mapname,inname
integer                               :: syear,eyear
integer                               :: smonth,sday,shour,smin
character(len=128)                    :: longname,tag,ctime
character(len=8)                      :: units
!-map variables
real                                  :: gsize,west, north, east, south ! map boundries
integer                               :: latpx,lonpx,nflp    ! pixel size, calculated
real,allocatable                      :: globaltrue(:),rmdtrnd(:)
integer,allocatable                   :: nextX(:,:),nextY(:,:),ocean(:,:)
integer                               :: i,ios,N,slice,chunk,Noff
real                                  :: a,b
integer                               :: timeid,varid,latid,lonid 
integer                               :: ncidin,ncidout,varidin,varidout
integer                               :: ix,iy,nx,ny
integer,allocatable                   :: dt(:)
real,allocatable                      :: d1lat(:), d1lon(:)
integer,dimension(3)                  :: start,count
real,parameter                        :: rmis=1.0e20
integer                               :: num_threads
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
read(buf,*) num_threads
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
allocate(nextX(lonpx,latpx),nextY(lonpx,latpx),ocean(lonpx,latpx))
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

!*********************************
!============WRITING==============
!*********************************
!  ncidout
!  varidout
!---------------------------------
! write netCDF file
!******
!https://www.unidata.ucar.edu/software/netcdf/examples/programs/ 
!******
! Create the netCDF file
! for global map do not write below 60S
south=max(south,-60.0)
nx=lonpx
!ny=latpx-30.0/dble(gsize) ! writed only up -60S latitude
ny=(north-south)/gsize !dble(gsize)
write(tag,'(i4.0,a,i4.0)')syear,"-",eyear
! edited the file name as /CaMa_out/{mapname}_{inputname}/{var}{syear}-{eyear}.nc
fname=trim(adjustl(outdir))//"/CaMa_out/"//trim(mapname)//"_"//trim(inname)//"/standardized"//trim(tag)//".nc"
print*, "create",fname
call nccheck( nf90_create(fname, NF90_NETCDF4, ncidout) )
!call nccheck( nf90_create_par(fname, IOR(NF90_NETCDF4,NF90_MPIIO), MPI_COMM_WORLD, MPI_INFO_NULL, ncidout) )
!=== set dimension ===
print*, "set dimension"
call nccheck( nf90_def_dim(ncidout, 'time', N, timeid) )
call nccheck( nf90_def_dim(ncidout, 'lat', ny, latid) )
call nccheck( nf90_def_dim(ncidout, 'lon', nx, lonid) )
!=== define variables ===
print*, "define variables"
call nccheck( nf90_def_var(ncidout, 'lat', nf90_float, (/latid/), varid) )
call nccheck( nf90_put_att(ncidout, varid, 'long_name','latitude') )
call nccheck( nf90_put_att(ncidout, varid, 'units','degrees_north') )

call nccheck( nf90_def_var(ncidout, 'lon', nf90_float, (/lonid/), varid) )
call nccheck( nf90_put_att(ncidout, varid, 'long_name','longitude') )
call nccheck( nf90_put_att(ncidout, varid, 'units','degrees_east') )
! Create time unit
print*, "create variables metadata"
smonth=1
sday=1
shour=0
smin=0
write(ctime,'(a14,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') 'days since ',syear,'-',smonth,'-',sday,' ',shour,":",smin
print*, "time"
call nccheck( nf90_def_var(ncidout, 'time', nf90_int, (/timeid/), varid) ) 
call nccheck( nf90_put_att(ncidout, varid, 'long_name','time') )
call nccheck( nf90_put_att(ncidout, varid, 'units',ctime) )

!===
print*, "variable"
call nccheck( nf90_def_var(ncidout, "standardize", nf90_float, &
         (/lonid,latid,timeid/), varidout) )
print*, "put attribute"
call nccheck( nf90_put_att(ncidout, varidout, 'long_name', "trend, sesonality removed and standardized water suface elevation") )
call nccheck( nf90_put_att(ncidout, varidout, 'units',     "none") )
call nccheck( nf90_put_att(ncidout, varidout, '_fillvalue',rmis) )
print*, "chunking : Write"
chunk=100
slice=100
Noff=int((real(N)/real(slice)))*slice
call nccheck( nf90_def_var_chunking(ncidout, varidout, 0, (/ chunk /)))

!===set collective I/O globally===
!call nccheck( nf90_var_par_access(ncidout, nf90_global, nf90_collective) )

!===end header===
print*, "end def"
call nccheck( nf90_enddef(ncidout) )

!=== put lon lat info ===
print*, "put lon, lat info"
allocate(d1lat(ny),d1lon(nx))
do ix=1,nx
  d1lon(ix)=west +(dble(ix)-0.5d0)*gsize !(east-west)  /dble(nx)
enddo
do iy=1,ny
  d1lat(iy)=north-(dble(iy)-0.5d0)*gsize !(north-(-60.0))/dble(ny)
enddo

call nccheck( nf90_inq_varid(ncidout,'lon',varid) )
call nccheck( nf90_put_var(ncidout,varid,d1lon))

call nccheck( nf90_inq_varid(ncidout,'lat',varid) )
call nccheck( nf90_put_var(ncidout,varid,d1lat))

!===put time===
print*, "put time"
allocate(dt(N))
dt = (/(i, i=1,N,1)/)! days since syear, smonth, sday

call nccheck( nf90_inq_varid(ncidout,'time',varid) )
call nccheck( nf90_put_var(ncidout,varid,dt) )

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
call nccheck( nf90_open(fname, nf90_nowrite, ncidin) )
!call nccheck( nf90_open_par(fname, IOR(NF90_NETCDF4,NF90_MPIIO), MPI_COMM_WORLD, MPI_INFO_NULL, ncidin) )
call nccheck( nf90_inq_varid(ncidin, trim(varname),varidin) )
!===set collective I/O globally===
!call nccheck( nf90_var_par_access(ncidin, nf90_global, nf90_collective) )
!===define chnuking ====
! print*, "chunking : Read"
! call nccheck( nf90_def_var_chunking(ncidin, varidin, 0, (/ ix,iy,1 /)))
!------------
allocate(globaltrue(N))

print*,N
! parallel calculation
print*, "define open mp threads", num_threads
!$ call omp_set_num_threads(num_threads)
!$ print*, omp_get_num_threads()
! print*, "open loop"
! !$omp parallel default(private) shared(ocean,N,iy,ncidin,varidin,ncidout,varidout)
! !$omp critical
! !$omp do
do ix = 1,nx ! pixels along longtitude direction
    do iy = 1,ny ! pixel along latitude direction
        start=(/ix,iy,1/)
        count=(/1,1,N/)
        !remove ocean
        if (ocean(ix,iy) == 0) then
            print*, ix, iy
            ! get variable subset
            call nccheck( nf90_get_var(ncidin,varidin,globaltrue,start=start,count=count) )
            !globaltrue = 100
            !write(*,*) globaltrue
            
            ! standadize
            !$omp parallel sections
            !$omp section
            call standadize(globaltrue,N)
            !$omp section
            !$omp end parallel sections
            !$print*, omp_get_num_threads()
            ! print*, globaltrue
            !--write variable--
            ! call nccheck( nf90_put_var(ncidout,varidout,globaltrue,start=start,count=count) )
            ! write variable as chuck and slice
            !$omp parallel sections
            !$omp section
            !$omp do
            do i=1,Noff,slice
                ! print*,i, N, Noff
                start=(/ix,iy,i/)
                count=(/1,1,slice/)
                ! print*, "L241:",start, count
                call nccheck( nf90_put_var(ncidout,varidout,globaltrue(i:i+slice),start=start,count=count) )
            end do
            !$omp end do
            start=(/ix,iy,Noff+1/)
            count=(/1,1,N-Noff/)
            ! print*, "last",start, count
            call nccheck( nf90_put_var(ncidout,varidout,globaltrue(Noff+1:N),start=start,count=count) )
            !$omp section
            !$omp end parallel sections
        end if
    end do
end do
! !$omp end do
! !$omp end critical
! !$omp end parallel
!---
!====close netcdf=====
call nccheck( nf90_close(ncidin ) )
call nccheck( nf90_close(ncidout) )
!---
deallocate(nextX,nextY,ocean,globaltrue) !,rmdtrnd)
end program standardized
!****************
SUBROUTINE standadize(data,N) 
! standadize the data
implicit none
integer                   :: N,i
real,dimension(N)         :: data,data2
real                      :: mean,std,var
! calculate mean 
mean = sum(data)/N
data2= data*data
var=0.0
do i=1,N
  var=var+(data(i)-mean)**2.0
end do
!----
std=sqrt(var/(N-1))
do i=1,N
  data(i) = (data(i) - mean)/(std+1.0e-20)
end do
return
END SUBROUTINE standadize
!***************
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
!*********************************
function units(varname)
implicit none
character(len=32)                 :: varname
character(len=8)                  :: units
if (trim(varname)=="sfcelv") then
    units="m"
end if
if (trim(varname)=="outflw" .or. trim(varname)=="rivout") then
    units="m3/s"
end if
return
end function units
!*********************************
function longname(varname)
implicit none
character(len=32)                 :: varname
character(len=128)                :: longname
if (trim(varname)=="sfcelv") then
    longname="water surface elevation"
end if
if (trim(varname)=="outflw") then
    longname="total discharge (river+floodplain)"
end if
if (trim(varname)=="rivout") then
    longname="river discharge"
end if

return
end function longname
!*********************************