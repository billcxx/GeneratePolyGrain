!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!	Task:	Generate Polycrystal grains for ferro code
! 	Author: Xiaoxing Cheng
! 	Date:	23/10/2013
! 	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

program GeneratePolyGrain

	implicit none

	type point
		real *8 x,y,z
	end type point
	
	type, extends(point) grainCorePoint
		integer index
		real *8 theta,phi,psi
	end type grainCorePoint

	type gridPoint
		integer grain
		logical isBorder
		logical isCorner
		integer otherSideOfBorder
		integer otherGrainAtCorner(3)
		real *8 theta,phi,psi
	end type gridPoint


	integer grainNum
	integer gridNumX,gridNumY,gridNumZ
	integer i,j,k,l
	integer grainIndex
	type(grainCorePoint), allocatable :: grainCore(:)
	type(gridPoint), allocatable :: grid(:,:,:)
	real *8	distance
	real *8 small
	real *8	wallSharpness;
	real *8 allocatable::PolyTR(:,:,:,:,:)


	allocate(grid(gridNumX+1,gridNumY+1,gridNumZ+1))
	allocate(grainCore(grainNum))

	small = Distance2PointsPeriodic(grid(1,1,1),grainCore(1))
! 	Find the nearest grain core to a grid point
! 	thus assign the grain index to the grid point
! 	taken into consideration of periodic boundary conditions


	do i = 1, gridNumX+1, 1
		do j = 1, gridNumY+1, 1
			do k = 1, gridNumZ+1, 1
! 	Setup the cordinate for the current grid point, if the point
! 	is out of gridNum then set its codinate to be the corresponding one
! 	in the box
				SetCordinateForPoint(grid(i,j,k),i,j,k)
				grid(i,j,k)=FindCorrespondPointInRange(grid(i,j,k),gridNumX,gridNumY,gridNumZ)

				do l = 1, grainNum, 1
					distance=Distance2PointsPeriodic(grid(i,j,k),grainCore(l),gridNumX,gridNumY,gridNumZ)
					if ( distance < small ) then
						small = distance
						grainIndex = l
					end if
				end do
				grid(i,j,k)%grain = grainIndex

				grid(i,j,k)%isBorder = .false.

			end do
		end do
	end do

! 	Find the boundary of grains by searching through each grid point and 
! 	mark those whose neighbour's belongs to different grains
! 	Taken into consideration of periodic boundary conditions

	grid(:,:,gridNumZ+1)=grid(:,:,1)
	grid(:,gridNumY+1,:)=grid(:,1,:)
	grid(gridNumX+1,:,:)=grid(1,:,:)

	do i = 1, gridNumX, 1
		do j = 1, gridNumY, 1
			do k = 1, gridNumZ, 1
				if ( grid(i,j,k)%grain /= grid(i,j,k+1)%grain ) then
					grid(i,j,k)%isBorder = .true.
					grid(i,j,k+1)%isBorder = .true.
					grid(i,j,k)%otherSideGrain = grid(i,j,k+1)%grain
					grid(i,j,k+1)%otherSideGrain = grid(i,j,k)%grain
				end if
			end do
		end do
	end do




end program GeneratePolyGrain

subroutine IsBoderInDirection(point1,axis)
	type(point), intent(in) :: point1
	integer, intent(in) :: axis

	select case (axis)
		case (1)
			i=i+1
		case (2)
			j=j+1
		case (3)
			k=k+1
		case default

				
			
			
	end select

	
end subroutine IsBoderInDirection(point,axis)


subroutine SetPointOrientation(point1,angle1,angle2,angle3)
	type(point), intent(in) :: point1 
	real *8, intent(in) :: angle1,angle2,angle3

	point1%theta = angle1
	point1%phi = angle2
	point1%psi = angle3
	
end subroutine SetPointOrientation


subroutine SetCordinateForPoint(point1,x1,y1,z1)
	type(point), intent(in) :: point1
	integer, intent(in) :: x1,y1,z1

	point1%x=x1
	point1%x=y1
	point1%x=z1
	
end subroutine SetCordinateForPoint



subroutine FindCorrespondPointInRange(point1,nx,ny,nz)
	type(point), intent(in) :: point1
	type(point), intent(put) :: point2
	integer , intent(in) :: nx
	integer , intent(in) :: ny
	integer , intent(in) :: nz

	point2 = point1
	if ( poin1%x>nx .or. point1%y>ny .or. point1%z>nz) then
		point2%x = point1%x-nx				
		point2%y = point1%y-ny				
		point2%z = point1%z-nz				
	end if

	return
end subroutine FindCorrespondPointInRange



subroutine Distance2PointsPeriodic(point1,point2,nx,ny,nz)
	type(point), intent(in) :: point1
	type(point), intent(in) :: point2
	type(point) pointA
	type(point) pointB
	integer *8, intent(in) :: nx
	integer *8, intent(in) :: ny
	integer *8, intent(in) :: nz
	integer i
	real *8 , intent(out):: distance
	real *8 , distancePeriod(7)

	pointA=point1
	pointB=point2
	distancePeriod(1)=GetDistanceBetween2Points(pointA,pointB)
	do i = 2, 3, 1
		pointB%x=point2%x+(-1)**i*nx
		distancePeriod(i)=GetDistanceBetween2Points(pointA,pointB)
	end do
	
	do i = 4, 5, 1
		pointB%y=point2%y+(-1)**i*ny
		distancePeriod(i)=GetDistanceBetween2Points(pointA,pointB)
	end do

	do i = 6, 7, 1
		pointB%z=point2%z+(-1)**i*nz
		distancePeriod(i)=GetDistanceBetween2Points(pointA,pointB)
	end do

	distance=minval(distancePeriod, dim=1)	
	
	return
end subroutine Distance2PointsPeriodic



subroutine GetDistanceBetween2Points(point1,point2)
	type(point), intent(in) :: point1
	type(point), intent(in) :: point2
	real *8 , intent(out):: distance

	distance=sqrt((point1%x-point2%x)**2+(point1%y-point2%y)**2+(point1%z-point2%z)**2)
	return
end subroutine GetDistanceBetween2Points