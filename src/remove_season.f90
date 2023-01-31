program remove_season
!===================================================
! remove seanality of 3 month, 6 month and integer multiples of year
! seasonality was found using fft
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
real,allocatable                      :: globaltrue(:),rmdsesn(:),data1(:)
integer,allocatable                   :: nextX(:,:),nextY(:,:),ocean(:,:)
integer                               :: i,ios,N,NN,slice,chunk,Noff
integer                               :: timeid,varid,latid,lonid
integer                               :: ncidin,ncidout,varidin,varidout
integer                               :: ix,iy,nx,ny
integer,allocatable                   :: dt(:)
real,allocatable                      :: d1lat(:), d1lon(:)
integer,dimension(3)                  :: start,count
real,parameter                        :: rmis=1.0e20,coeff=0.1,P=365.0
!==========================================================
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
NN=2**(int(log(real(N))/log(2.0))+1)
!
write(*,*) N, NN
!
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
fname=trim(adjustl(outdir))//"/CaMa_out/"//trim(mapname)//"_"//trim(inname)//"/rmdsesn"//trim(tag)//".nc"
print*, "create",fname
call nccheck( nf90_create(fname, NF90_NETCDF4, ncidout) )
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
call nccheck( nf90_def_var(ncidout, "rmdsesn", nf90_float, &
         (/lonid,latid,timeid/), varidout) )
print*, "put attribute"
call nccheck( nf90_put_att(ncidout, varidout, 'long_name', 'trend & seasonlaty removed water surface elevation') ) !trim(longname(varname))
call nccheck( nf90_put_att(ncidout, varidout, 'units',     trim(units(varname))) )
call nccheck( nf90_put_att(ncidout, varidout, '_fillvalue',rmis) )
print*, "chunking : Write"
chunk=100
slice=100
Noff=int((real(N)/real(slice)))*slice
call nccheck( nf90_def_var_chunking(ncidout, varidout, 0, (/ chunk /)))

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
call nccheck( nf90_inq_varid(ncidin, trim(varname),varidin) )
allocate(globaltrue(N),rmdsesn(N))
!--
print*,N
!--
allocate(data1(NN))!,anomaly(lonpx,latpx,N))
!allocate(xt(N),xf(N),data1(NN),data11(NN),anomaly(lonpx,latpx,N))
!xt = (/(real(i), i=1,N,1)/)
! zero padding
data1=0.0
!data11=0.0

! parallel calculation
!$omp parallel default(private) shared(ocean,globaltrue,N,NN) 
!$omp do
!$write(*,*) omp_get_num_threads()
do ix = 1,nx ! pixels along longtitude direction
    do iy = 1,ny ! pixel along latitude direction
        start=(/ix,iy,1/)
        count=(/1,1,N/)
        !remove ocean
        if (ocean(ix,iy) /= 0) then
            cycle
        end if
        ! get variable subset
        call nccheck( nf90_get_var(ncidin,varidin,globaltrue,start=start,count=count) )
        ! zero padding
        data1=0.0
        data1(1:N)=globaltrue(:)
        call REALFT(data1,NN,+1)
        !remove frequency other than multiple of p days
        !covarite of 6
        !frequency of 90 and 180 also kept
        print*, ix,iy,NN,shape(data1)
        call iff_p(data1,NN,P)
        call REALFT(data1,NN,-1)
        !call FOUR1(data11,NN/2,-1)
        rmdsesn(:)=globaltrue(:)-(2.0/real(NN)) *data1(1:N) 
        !--write variable--
        ! call nccheck( nf90_put_var(ncidout,varidout,rmdsesn,start=start,count=count) )
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
    end do
end do
!$omp end do
!$omp end parallel
!----
print*, "save netCDF"
print*, trim(adjustl(outdir))//"/CaMa_out/"//trim(mapname)//"_"//trim(inname)//"/rmdsesn"//trim(tag)//".nc"
!====close netcdf=====
call nccheck( nf90_close(ncidin ) )
call nccheck( nf90_close(ncidout) )
!---
print*, "deallocate"
deallocate(nextX,nextY,ocean,d1lat,d1lon,dt,globaltrue,data1,rmdsesn)
print*,"end program"
end program remove_season
!**********************************
  SUBROUTINE REALFT (DATA, N, ISIGN)
!
! Numerical recipes FFT of a real function routine.
!
! Calculates the Fourier transform of a set of n real-valued data points.
! Replaces this data by the positive frequency half of its complex Fourier
! transform.  
! The real valued first and last components of the complex transform are returned as data(1) and data(2). 
! n must be a power of 2.  
! For isign=-1 the inverse transform is returned and must be multiplied by 2/n.
!
      REAL(8) WR, WI, WPR, WPI, WTEMP, THETA
      REAL(4) DATA ( * )
      THETA = 3.141592653589793D0 / DBLE (N / 2)
      C1 = 0.5
      IF (ISIGN.EQ.1) THEN
         C2 = - 0.5
         CALL FOUR1 (DATA, N / 2, + 1)
      ELSE
         C2 = 0.5
         THETA = - THETA
      ENDIF
      WPR = - 2.0D0 * DSIN (0.5D0 * THETA) **2
      WPI = DSIN (THETA)
      WR = 1.0D0 + WPR
      WI = WPI
      N2P3 = N + 3
      DO 11 I = 2, N / 4
         I1 = 2 * I - 1
         I2 = I1 + 1
         I3 = N2P3 - I2
         I4 = I3 + 1
         WRS = SNGL (WR)
         WIS = SNGL (WI)
         H1R = C1 * (DATA (I1) + DATA (I3) )
         H1I = C1 * (DATA (I2) - DATA (I4) )
         H2R = - C2 * (DATA (I2) + DATA (I4) )
         H2I = C2 * (DATA (I1) - DATA (I3) )
         DATA (I1) = H1R + WRS * H2R - WIS * H2I
         DATA (I2) = H1I + WRS * H2I + WIS * H2R
         DATA (I3) = H1R - WRS * H2R + WIS * H2I
         DATA (I4) = - H1I + WRS * H2I + WIS * H2R
         WTEMP = WR
         WR = WR * WPR - WI * WPI + WR
         WI = WI * WPR + WTEMP * WPI + WI
11    END DO
      IF (ISIGN.EQ.1) THEN
        H1R = DATA (1)
        DATA (1) = H1R + DATA (2)
        DATA (2) = H1R - DATA (2)
      ELSE
        H1R = DATA (1)
        DATA (1) = C1 * (H1R + DATA (2) )
        DATA (2) = C1 * (H1R - DATA (2) )
        CALL FOUR1 (DATA, N / 2, - 1)
      ENDIF
      RETURN
END SUBROUTINE REALFT
!********************
!  four1 -- this is the four1 routine from numerical recipes
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END SUBROUTINE FOUR1
!*************************
SUBROUTINE iff_p(data,N,P)
!remove frequency other than multiple of p days
!covarite of 6
! real and imaginary part saved in adjecent places 
implicit none
integer               :: N
real                  :: P
real,dimension(N)     :: data,data_temp
!--
integer               :: i,j,k,c
!---
data_temp=0.0
!--
!write(*,*)shape(data_temp)
! integer multiples of P
c=3
do i=1,int(real(N)/(2.0*P))
  j=int(real(N)/(real(i)*P))
  !write(*,*) i,j 
  do k=j-c,j+c 
    if (k<=2) then
      cycle
    end if
    data_temp(2*k+1)=data(2*k+1)
    data_temp(2*k+2)=data(2*k+2)
  end do
end do
! 180, 90 day cycles
! (1/2) * P and (1/4) * P cycles
do i=1,2
  j=int(real(N)/((1.0/2.0*real(i))*P))
  !write(*,*) i,j 
  do k=j-c,j+c 
    if (k<=2) then
      cycle
    end if
    data_temp(2*k+1)=data(2*k+1)
    data_temp(2*k+2)=data(2*k+2)
  end do
end do
! data(1) and data(2) is saved seperately
! first two indices are saved
data_temp(1)=data(1)
data_temp(2)=data(2)
!--
data=data_temp
return 
END SUBROUTINE iff_p
!***************************
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
