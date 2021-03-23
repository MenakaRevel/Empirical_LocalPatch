program bin2nc
!$ use omp_lib
use netcdf
implicit none
character(len=128)                    :: fname,buf,camadir,outdir
character(len=32)                     :: outname,varname,mapname,inname
character(len=4)                      :: yyyy
integer                               :: year,syear,eyear,dyear
integer                               :: smonth,sday,shour,smin
character(len=128)                    :: longname,tag,ctime
character(len=8)                      :: units
!-map variables
real                                  :: gsize,west, north, east, south ! map boundries
integer                               :: latpx,lonpx,nflp    ! pixel size, calculated
real,allocatable                      :: globaltrue(:,:,:)
integer*4                             :: i,j!,i_m,j_m,pixel
integer                               :: days,ios,N!,countnum
integer                               :: ncid,timeid,varid,latid,lonid
integer                               :: ix,iy,nx,ny
integer,allocatable                   :: dt(:)
real,allocatable                      :: d1lat(:), d1lon(:)
integer,dimension(3)                  :: start,count
real,parameter                        :: rmis=1.0e20
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
!print*, north,south,gsize
!print*, nx,ny
write(tag,'(i4.0,a,i4.0)')syear,"-",eyear
! edited the file name as /CaMa_out/{mapname}_{inputname}/{var}{syear}-{eyear}.nc
fname=trim(adjustl(outdir))//"/CaMa_out/"//trim(mapname)//"_"//trim(inname)//"/"//trim(varname)//trim(tag)//".nc"
print*, "create",fname
call nccheck( nf90_create(fname, NF90_NETCDF4, ncid) )
!=== set dimension ===
print*, "set dimension"
call nccheck( nf90_def_dim(ncid, 'time', N, timeid) )
call nccheck( nf90_def_dim(ncid, 'lat', ny, latid) )
call nccheck( nf90_def_dim(ncid, 'lon', nx, lonid) )
!=== define variables ===
print*, "define variables"
call nccheck( nf90_def_var(ncid, 'lat', nf90_float, (/latid/), varid) )
call nccheck( nf90_put_att(ncid, varid, 'long_name','latitude') )
call nccheck( nf90_put_att(ncid, varid, 'units','degrees_north') )

call nccheck( nf90_def_var(ncid, 'lon', nf90_float, (/lonid/), varid) )
call nccheck( nf90_put_att(ncid, varid, 'long_name','longitude') )
call nccheck( nf90_put_att(ncid, varid, 'units','degrees_east') )
! Create time unit
print*, "create variables metadata"
smonth=1
sday=1
shour=0
smin=0
write(ctime,'(a14,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') 'days since ',syear,'-',smonth,'-',sday,' ',shour,":",smin
print*, "time"
call nccheck( nf90_def_var(ncid, 'time', nf90_int, (/timeid/), varid) ) 
call nccheck( nf90_put_att(ncid, varid, 'long_name','time') )
call nccheck( nf90_put_att(ncid, varid, 'units',ctime) )

!===
print*, "variable"
call nccheck( nf90_def_var(ncid, trim(varname), nf90_float, &
         &(/lonid,latid,timeid/), varid) )
print*, "put attribute"
call nccheck( nf90_put_att(ncid, varid, 'long_name', trim(longname(varname))) )
call nccheck( nf90_put_att(ncid, varid, 'units',     trim(units(varname))) )
call nccheck( nf90_put_att(ncid, varid, '_fillvalue',rmis) )

!===end header===
print*, "end def"
call nccheck( nf90_enddef(ncid) )

!=== put lon lat info ===
print*, "put lon, lat info"
allocate(d1lat(ny),d1lon(nx))
do ix=1,nx
  d1lon(ix)=west +(dble(ix)-0.5d0)*gsize !(east-west)  /dble(nx)
enddo
do iy=1,ny
  d1lat(iy)=north-(dble(iy)-0.5d0)*gsize !(north-(-60.0))/dble(ny)
enddo

call nccheck( nf90_inq_varid(ncid,'lon',varid) )
call nccheck( nf90_put_var(ncid,varid,d1lon))

call nccheck( nf90_inq_varid(ncid,'lat',varid) )
call nccheck( nf90_put_var(ncid,varid,d1lat))

!===put time===
print*, "put time"
allocate(dt(N))
do i=1,N
    dt(i)=i ! days since syear, smonth, sday
end do
call nccheck( nf90_inq_varid(ncid,'time',varid) )
call nccheck( nf90_put_var(ncid,varid,dt) )

!===put variable===
print*, "put variables"
count=(/nx,ny,1/)
start=(/1,1,1/)
call nccheck( nf90_inq_varid(ncid,trim(varname),varid) )
! read variable
i=1
do year=syear,eyear
    write(yyyy,'(i4.0)') year
    write(*,*) yyyy !day,, days(day) 
    days=dyear(year)
    allocate(globaltrue(lonpx,latpx,days))
    varname=outname
    fname=trim(adjustl(outdir))//"/CaMa_out/"//trim(mapname)//"_"//trim(inname)//"/"//trim(varname)//yyyy//".bin"
    !print *,"L153",trim(varname),days!fname
    open(34,file=fname,form="unformatted",access="direct",recl=4*days*lonpx*latpx,status="old",iostat=ios)
    if(ios==0)then
        read(34,rec=1) globaltrue(:,:,:)
    else
        write(*,*) "no",trim(varname),trim(fname)
    end if
    close(34)
    !print*, "L161",trim(varname)
    start(3)=i
    count(3)=days
    call nccheck( nf90_put_var(ncid,varid,globaltrue(1:nx,1:ny,1:days),start=start,count=count) )
    !do j=1,days
    !    start(3)=j+i
    !    call nccheck( nf90_put_var(ncid,varid,globaltrue(1:nx,1:ny,i+j),start=start,count=count) )
    !    !print*,"L165", trim(varname),j
    !end do
    deallocate(globaltrue)
    i=i+days
    !print*,"L169",trim(varname)
end do
!===close netCDF4===
call nccheck( nf90_close(ncid))

print* , "***", trim(varname)//trim(tag), "sucessfully created***"
deallocate(d1lat,d1lon,dt)
end program bin2nc
!*********************************
function dyear(year)
implicit none
integer                        :: year,dyear
!real                           :: mod
!--
! days of the year
! normal : 365
! leap   : 366
dyear=365
if (mod(dble(year),4.0)   == 0) dyear=366
if (mod(dble(year),100.0) == 0) dyear=365
if (mod(dble(year),400.0) == 0) dyear=366
return
end function dyear
!*********************************
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