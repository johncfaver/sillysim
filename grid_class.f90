module grid_class
    
    implicit none
    private

    public :: integerGrid, realGrid, doubleGrid
    public :: createIntegerGrid
    
    type, abstract :: grid
        double precision                                :: resolution
        double precision, dimension(3)                  :: center
        double precision, dimension(3)                  :: dimensions
        double precision, dimension(3)                  :: origin
        integer, dimension(3)                           :: numPoints
        contains
            procedure :: coordinate
            procedure :: pointIsInGrid
            procedure :: coordToPoint
            procedure :: pointToCoord
            procedure :: writeXYZ
    end type grid

    type, extends(grid) :: integerGrid
        integer, dimension(:,:,:), allocatable          :: values
    end type integerGrid
     
    type, extends(grid) :: realGrid
        real, dimension(:,:,:), allocatable             :: values
    end type realGrid
    
    type, extends(grid) :: doubleGrid
        double precision, dimension(:,:,:), allocatable :: values
    end type doubleGrid
    
    contains
        function createIntegerGrid(center,dimensions,resolution)
            double precision, dimension(3)  :: center
            double precision, dimension(3)  :: dimensions
            double precision                :: resolution
            double precision, dimension(3)  :: origin
            type(integerGrid), target       :: createIntegerGrid
            integer                         :: n1,n2,n3
            createIntegerGrid%resolution = resolution
            createIntegerGrid%center     = center
            createIntegerGrid%dimensions = dimensions
            createIntegerGrid%numPoints  = dimensions / resolution + 1
            createIntegerGrid%origin     = center - dimensions/2.0d0
            n1 = createIntegerGrid%numPoints(1)
            n2 = createIntegerGrid%numPoints(2)
            n3 = createIntegerGrid%numPoints(3)
            allocate(createIntegerGrid%values(n1,n2,n3))
            createIntegerGrid%values = 0
        end function createIntegerGrid
    
        function coordinate(self, gridpoint)
            !From grid point, return realspace 3D coordinates
            class(grid)                     :: self
            integer, dimension(3)           :: gridpoint
            double precision, dimension(3)  :: coordinate
            coordinate = self%origin + (gridpoint-1)*self%resolution
        end function coordinate
        
        function pointIsInGrid(self, x,y,z)
            class(grid)             :: self
            double precision        :: x,y,z
            logical                 :: pointIsInGrid
            pointIsInGrid = .true.
            if(x < self%origin(1)) pointIsInGrid = .false.
            if(y < self%origin(2)) pointIsInGrid = .false.
            if(z < self%origin(3)) pointIsInGrid = .false.
            if(x > self%origin(1) + self%dimensions(1)) pointIsInGrid = .false.
            if(y > self%origin(2) + self%dimensions(2)) pointIsInGrid = .false.
            if(z > self%origin(3) + self%dimensions(3)) pointIsInGrid = .false.
            pointIsInGrid = pointIsInGrid
        end function pointIsInGrid

        function pointToCoord(self, a,b,c)
            class(grid)                     :: self
            double precision, dimension(3)  :: pointToCoord
            integer                         :: a,b,c
            pointToCoord(1) = self%origin(1) + (a-1)*self%resolution
            pointToCoord(2) = self%origin(2) + (b-1)*self%resolution
            pointToCoord(3) = self%origin(3) + (c-1)*self%resolution
        end function pointToCoord

        function coordToPoint(self, x,y,z)
            class(grid)             :: self
            double precision        :: x, y, z
            integer, dimension(3)   :: coordToPoint
            if(self%pointIsInGrid(x,y,z))then
                coordToPoint(1) = (x - self%origin(1))/self%resolution + 1
                coordToPoint(2) = (y - self%origin(2))/self%resolution + 1
                coordToPoint(3) = (z - self%origin(3))/self%resolution + 1
            else
                print *, 'WARNING: Point (',x,',',y,',',z,') is outside grid.'
                coordToPoint = -1
            endif
        end function coordToPoint

        subroutine writeXYZ(self,gridfile)
            class(grid)                 :: self
            integer                     :: i,j,k
            character(len=*), optional  :: gridfile
            character(len=50)           :: gridfile_ ="grid.xyz"
            if(present(gridfile)) gridfile_ = gridfile
            open(23,file=gridfile_,action="WRITE",status="REPLACE")
            write(23,*) "GRID"
            write(23,*) self%numPoints(1)*self%numPoints(2)*self%numPoints(3)
            do i=1,self%numPoints(1)
                do j=1,self%numPoints(2)
                    do k=1,self%numPoints(3)
                        select type(self)
                        class is(integerGrid)
                            write(23,*) self%values(i,j,k), self%coordinate([i,j,k])
                        class is(realGrid)
                            write(23,*) self%values(i,j,k), self%coordinate([i,j,k])
                        class is(doubleGrid)
                            write(23,*) self%values(i,j,k), self%coordinate([i,j,k])
                        end select
                    enddo
                enddo
            enddo
            close(23)
        end subroutine writeXYZ

end module grid_class
