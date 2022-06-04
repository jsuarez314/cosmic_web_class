
! grid tools for quick location, search
module LSS_grid_tools

use LSS_tools

implicit none

	type :: gridcell
		integer :: np=0 ! how many points
		integer, allocatable :: iplist(:) ! list of index of particles
	end type

	type :: gridtype
		real(dl), allocatable :: xs(:), ys(:), zs(:), ws(:) ! x,y,z and weight
		integer :: np=0, nx=0, ny=0, nz=0 ! #-of-elemetns; #-of-cells in x,y,z
		real(dl) :: xmin,xmax,ymin,ymax,zmin,zmax,wmin,wmax, & ! range of xs, ys, zs, ws
				gridxmin,gridxmax,gridymin,gridymax,gridzmin,gridzmax, & ! range of 3D grid
				deltax,deltay,deltaz ! cell width
		type(gridcell), allocatable :: cellmat(:,:,:)
	end type gridtype

	integer, parameter :: maxnumip = 100000 ! maximal number used in array


contains

  !------------------------------------------
  ! calculate rho at some position
	subroutine grid_fixr_rho_drho(grid,x,y,z,fixr,numNB,touchboundary, &
			rhomean, rho, drhodx, drhody, drhodz, &
			nbarmean, nbar, dnbardx, dnbardy, dnbardz)
		! Dummy 
		type(gridtype), intent(in) :: grid
		real(dl),intent(in) :: x,y,z,fixr
		integer, intent(out) :: numNB
		logical, intent(out) :: touchboundary
		real(dl),intent(out) :: rhomean, rho, drhodx, drhody, drhodz
		real(dl),intent(out), optional :: nbarmean, nbar, dnbardx, dnbardy, dnbardz
		! Local
		integer :: iNB, ip
		real(dl) :: rlist(maxnumip), h, x1,y1,z1, nowr,weight, dweight
		integer :: iplist(maxnumip)
		
		! find out list of NBs
		call grid_fixr_search(grid,x,y,z,fixr,numNB,touchboundary,iplist,rlist)

		! compute rhomean (sum of weight / volume), rho, drho (density, gradient estimated using spline kernel); similarly, nbar means number density rather than weighted density
		rhomean=0; rho=0; drhodx=0;drhody=0;drhodz=0;	
		if(present(nbarmean)) then
		  if(.not.present(nbar).or..not.present(dnbardx) &
			.or..not.present(dnbardy).or..not.present(dnbardz)) then
		    print *, ' (grid_fixr_rho_drho) ERROR! '//&
			'nbar, dnbardx, dnbardy, dnbardz must appear together!!'
		    stop
		  endif
		  nbarmean=0; nbar=0; dnbardx=0; dnbardy=0; dnbardz=0
		endif
  		h = fixr / 2.0
		do iNB = 1, numNB
			nowr = rlist(iNB)
			ip = iplist(iNB); 
			x1=grid%xs(ip);y1=grid%ys(ip);z1=grid%zs(ip); weight = grid%ws(ip)
			dweight = der_w_kernel(nowr,h)
			! rho and drho
			rhomean = rhomean + weight
			rho = rho + weight*w_kernel(nowr, h)
			drhodx = drhodx + weight*(x-x1) / nowr * dweight
			drhody = drhody + weight*(y-y1) / nowr * dweight
			drhodz = drhodz + weight*(z-z1) / nowr * dweight
			! nbar and dnbar
			if(present(nbarmean)) then
				nbarmean = nbarmean + 1.0
				nbar = nbar + w_kernel(nowr, h)
				dnbardx = dnbardx + (x-x1) / nowr * dweight
				dnbardy = dnbardy + (y-y1) / nowr * dweight
				dnbardz = dnbardz + (z-z1) / nowr * dweight
			endif
		enddo
		rhomean = rhomean  / sphere_vol(fixr)
		if(present(nbarmean)) nbarmean = nbarmean  / sphere_vol(fixr)
	end subroutine grid_fixr_rho_drho
  !------------------------------------------
  ! check whether a point touches boundary
	logical function check_touchbd(grid,x,y,z,dr)
		! Dummy
		type(gridtype), intent(in) :: grid
		real(dl) :: x,y,z,dr
		if(abs(x-grid%xmin)<dr.or.abs(y-grid%ymin)<dr.or.abs(z-grid%zmin)<dr &
		.or.abs(x-grid%xmax)<dr.or.abs(y-grid%ymax)<dr.or.abs(z-grid%zmax)<dr &
		.or.abs(x-grid%gridxmin)<dr.or.abs(y-grid%gridymin)<dr.or.abs(z-grid%gridzmin)<dr &
		.or.abs(x-grid%gridxmax)<dr.or.abs(y-grid%gridymax)<dr.or.abs(z-grid%gridzmax)<dr) then
			check_touchbd = .true.
		else
			check_touchbd = .false.
		endif
	end function check_touchbd


  !------------------------------------------
  ! find ip of list, fixed position, fixed radius
	subroutine grid_fixr_search(grid,x,y,z,fixr, numNB,touchboundary,iplist,rlist)
		! Dummy 
		type(gridtype), intent(in) :: grid
		real(dl),intent(in) :: x,y,z,fixr
		integer, intent(out) :: numNB, iplist(maxnumip)
		logical, intent(out) :: touchboundary
		real(dl),optional,intent(out) :: rlist(maxnumip)
		! Local
		integer :: imin,imax,jmin,jmax,kmin,kmax, i,j,k,l,ip
		real(dl) :: nowr

  		imin = int((x-fixr-grid%gridxmin)/grid%deltax +1.0_dl)
  		imax = int((x+fixr-grid%gridxmin)/grid%deltax +1.0_dl)
  		jmin = int((y-fixr-grid%gridymin)/grid%deltay +1.0_dl)
  		jmax = int((y+fixr-grid%gridymin)/grid%deltay +1.0_dl)
  		kmin = int((z-fixr-grid%gridzmin)/grid%deltaz +1.0_dl)
  		kmax = int((z+fixr-grid%gridzmin)/grid%deltaz +1.0_dl)

		touchboundary = check_touchbd(grid,x,y,z,fixr)

		numNB = 0
		do i = max(1,imin), min(grid%nx,imax)
		do j = max(1,jmin), min(grid%ny,jmax)
		do k = max(1,kmin), min(grid%nz,kmax)
			do l = 1, grid%cellmat(i,j,k)%np
				ip = grid%cellmat(i,j,k)%iplist(l)
				nowr = distancexyz(grid%xs(ip),grid%ys(ip),grid%zs(ip),x,y,z)
				if(nowr < fixr) then
					numNB = numNB+1
					if(numNB > maxnumip) then
						print *, '(grid_fixr_search) Warning: outflow! '//&
						'increase size of maxnumip:', maxnumip, numNB
					endif
					iplist(numNB) = ip
					if(present(rlist)) rlist(numNB) = nowr
				endif
			enddo
		enddo
		enddo
		enddo
	end subroutine grid_fixr_search

  !------------------------------------------
  ! build up a grid; reading xyzw from file; automatic nx set as np^(1.0/3.0) (optional)
	subroutine grid_init_fromfile(grid,inputfile,printinfo,nx)
		! Dummy
		character(len=char_len) :: inputfile
		type(gridtype), intent(out) :: grid
		logical, intent(in) :: printinfo
		integer, optional, intent(in) :: nx
		! Local
		integer :: np
		real(dl), allocatable :: xyzw(:,:)

		if(printinfo) then
			write(*,'(A,A)') '   (grid_init) initializing 3d grid of cells from file: ', &
				trim(adjustl(inputfile))
		endif

		call read_in(inputfile, 4, np, xyzw)
		if(present(nx)) then
			call grid_init(grid, xyzw(1:np,1),xyzw(1:np,2),xyzw(1:np,3),xyzw(1:np,4), np, nx, printinfo)
			
		else
			call grid_init(grid, xyzw(1:np,1),xyzw(1:np,2),xyzw(1:np,3),xyzw(1:np,4), np, int(real(np)**(1.0/3.0)), printinfo)
		endif
	end subroutine
  !------------------------------------------
  ! build up a grid
	subroutine grid_init(grid,xs,ys,zs,ws,np,nx,printinfo)
		! Dummy
		integer :: np,nx
		real(dl) :: xs(np), ys(np), zs(np), ws(np)
		type(gridtype), intent(out) :: grid
		logical, intent(in) :: printinfo
		! Local
		integer :: i, ix,iy,iz
		integer, allocatable :: ipmat(:,:,:)
		real :: numhasdata, maxdatanum, avgdatanum, numpixel
		real(dl) :: delta, extrarat = 0.001 ! the spatial extension of grid is slightly larger than 0.1% of cell width, making sure all points are included into the grid

		! list of x,y,z
		if(printinfo) print *, '  (grid_init) initializing 3d grid of cells.'
		grid%np=np;
		if(allocated(grid%xs)) deallocate(grid%xs); if(allocated(grid%ys)) deallocate(grid%ys);
		if(allocated(grid%zs)) deallocate(grid%zs); if(allocated(grid%ws)) deallocate(grid%ws);
		allocate(grid%xs(np),grid%ys(np),grid%zs(np),grid%ws(np))
		grid%xs=xs; grid%ys=ys; grid%zs=zs; grid%ws=ws
		call grid_findrange(grid,.false.)
		if(printinfo) call grid_printrange(grid)

!		determine the range of 3d grid; keep deltax=deltay=deltaz
		grid%nx = nx
		! preliminary ranges
		delta = (grid%xmax-grid%xmin) / dble(grid%nx) 
		grid%gridxmin = grid%xmin - delta * extrarat
		grid%gridxmax = grid%xmax + delta * extrarat
 		grid%gridymin = grid%ymin - delta * extrarat 
		grid%gridymax = grid%ymax + delta * extrarat
		grid%gridzmin = grid%zmin - delta * extrarat 
		grid%gridzmax = grid%zmax + delta * extrarat
		! delta
		grid%deltax = (grid%gridxmax-grid%gridxmin) / dble(nx)
		grid%deltay=grid%deltax; grid%deltaz=grid%deltax
		! range, #-cell of y, z 
		grid%ny = ceiling( (grid%gridymax-grid%gridymin)/grid%deltay)
		grid%nz = ceiling( (grid%gridzmax-grid%gridzmin)/grid%deltaz)
		! re-calcluate maximal range as finial ranges
		grid%gridxmax = grid%gridxmin + grid%deltax*grid%nx
		grid%gridymax = grid%gridymin + grid%deltay*grid%ny
		grid%gridzmax = grid%gridzmin + grid%deltaz*grid%nz
		if(printinfo) call grid_printgrid(grid)

		! Build up the grid...

		if(allocated(grid%cellmat)) deallocate(grid%cellmat)
		allocate(grid%cellmat(grid%nx,grid%ny,grid%nz))

		allocate(ipmat(grid%nx,grid%ny,grid%nz))
		ipmat = 0
		do i=1,np
			call cellindex(grid,ix,iy,iz,xs(i),ys(i),zs(i))
			ipmat(ix,iy,iz) = ipmat(ix,iy,iz)+1
		enddo

		numhasdata=0; maxdatanum=0; avgdatanum=0
		do ix = 1, grid%nx
		do iy = 1, grid%ny
		do iz = 1, grid%nz
			if(ipmat(ix,iy,iz)>0) then
				grid%cellmat(ix,iy,iz)%np = ipmat(ix,iy,iz)
				allocate(grid%cellmat(ix,iy,iz)%iplist(grid%cellmat(ix,iy,iz)%np))
				! for stat
				numhasdata = 	numhasdata +1
				maxdatanum = max(maxdatanum, real(ipmat(ix,iy,iz)))
				avgdatanum = avgdatanum+ipmat(ix,iy,iz)
			endif
		enddo
		enddo
		enddo


		ipmat = 0
		do i=1,np
			call cellindex(grid,ix,iy,iz,xs(i),ys(i),zs(i))
			ipmat(ix,iy,iz) = ipmat(ix,iy,iz)+1
			grid%cellmat(ix,iy,iz)%iplist(ipmat(ix,iy,iz)) = i
		enddo

		numpixel = grid%nx*grid%ny*grid%nz
		avgdatanum=avgdatanum/real(numpixel)
		if(printinfo) then
		  write(*,'(18x,i7,A,f6.2,A,i6,A,f10.3)') int(numhasdata+0.5),'(',&
			real(numhasdata/real(numpixel))*100.0,'%) cells have point inside; maximal # is ',&
			 int(maxdatanum+0.5), ', mean num is ', real(avgdatanum)
		  print *, '  (grid_init done)'
		endif
	
	end subroutine grid_init


  !------------------------------------------
  ! find out range of xs,ys,zs,ws
	subroutine grid_findrange(grid,printinfo)
		! Dummy
		type(gridtype), intent(in) :: grid
		logical, intent(in) :: printinfo

		! find out ranges of x,y,z,w.		
		call find_min_max(grid%xs,grid%np,grid%xmin,grid%xmax)
		call find_min_max(grid%ys,grid%np,grid%ymin,grid%ymax)
		call find_min_max(grid%zs,grid%np,grid%zmin,grid%zmax)
		call find_min_max(grid%ws,grid%np,grid%wmin,grid%wmax)

		if(printinfo) call grid_printrange(grid)
	end subroutine grid_findrange


  !------------------------------------------
  ! print range of xs,ys,zs,ws
	subroutine grid_printrange(grid)
		type(gridtype), intent(in) :: grid

		write(*,'(20x,A,i12)'), '#-of points:        ', grid%np
		write(*,'(20x,A,2E14.5,A)'), 'Range of weight:    (', grid%wmin,grid%wmax,')'
		write(*,'(20x,A,2f9.3,A)'), 'Range of x:         (', grid%xmin,grid%xmax,')'
		write(*,'(20x,A,2f9.3,A)'), 'Range of y:         (', grid%ymin,grid%ymax,')'
		write(*,'(20x,A,2f9.3,A)'), 'Range of z:         (', grid%zmin,grid%zmax,')'

	end subroutine grid_printrange
  !------------------------------------------
  ! print range of xs,ys,zs,ws
	subroutine grid_printgrid(grid)
		type(gridtype), intent(in) :: grid

		write(*,'(21x,A,2f9.3,A)'), 'Grid range of x:   (', grid%gridxmin,grid%gridxmax,')'
		write(*,'(21x,A,2f9.3,A)'), 'Grid range of y:   (', grid%gridymin,grid%gridymax,')'
		write(*,'(21x,A,2f9.3,A)'), 'Grid range of z:   (', grid%gridzmin,grid%gridzmax,')'

		write(*,'(18x,A,i4,A,i4,A,i4,A,i11)')  ' Cell-#    =', &
				grid%nx,'*',grid%ny,'*',grid%nz,' = ',grid%nx*grid%ny*grid%nz
		write(*,'(18x,A,3(f10.5,A),e15.7)') ' Cell-size =', real(grid%deltax),' *',real(grid%deltay), &
				' *',real(grid%deltaz), ' =', grid%deltax*grid%deltay*grid%deltaz
	end subroutine grid_printgrid


  !------------------------------------------
  ! cell center position
  	subroutine cellpos(grid,ix,iy,iz,x,y,z)
		type(gridtype), intent(in) :: grid
		integer, intent(in) :: ix, iy, iz
  		real(dl), intent(out) :: x,y,z
  		x = grid%gridxmin + (ix-0.5)*grid%deltax
  		y = grid%gridymin + (iy-0.5)*grid%deltay
  		z = grid%gridzmin + (iz-0.5)*grid%deltaz
	end subroutine cellpos
  ! index of dewelling cell 
  	subroutine cellindex(grid,ix,iy,iz,x,y,z)
		type(gridtype), intent(in) :: grid
		integer, intent(out) :: ix, iy, iz
  		real(dl), intent(in) :: x,y,z
  		ix = int((x-grid%gridxmin) / grid%deltax +1.0)
  		iy = int((y-grid%gridymin) / grid%deltay +1.0)
  		iz = int((z-grid%gridzmin) / grid%deltaz +1.0)
	end subroutine cellindex
	



end module	

