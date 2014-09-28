! Various geometry-related functions
module geometry
    use atom_class
    
    implicit none
    private
   
    public :: deg2rad, rad2deg
    public :: r12, sumSquares, norm, normalized 
    public :: dotProduct, crossProduct
    public :: getBondLength, getBondLengthSquared
    public :: getAngle, getDihedral

    double precision, parameter :: pi = dacos(-1.0d0)
    double precision, parameter :: tau = 2.0d0*dacos(-1.0d0)
    double precision, parameter :: d2r = pi/180.0d0
    double precision, parameter :: r2d = 180.0d0/pi

contains
        function r12(a1,a2)
            ! Return vector from a1 to a2
            type(atom), intent(in)         :: a1,a2
            double precision, dimension(3) :: r12
            r12 = a2%coords - a1%coords 
        end function r12

        function norm(vector)
            ! Return magnitude (L-2 norm) of vector
            double precision, dimension(3), intent(in) :: vector
            double precision                           :: norm
            norm = dsqrt(sumSquares(vector)) 
        end function norm

        function normalized(vector)
            ! Return normalized vector
            double precision, dimension(3) :: vector, normalized
            normalized = vector / norm(vector)
        end function normalized

        function sumSquares(vector)
            ! Return Pythagorean sum (sum of squared elements)
            double precision, dimension(3), intent(in) :: vector
            double precision                           :: sumSquares
            integer                                    :: i
            sumSquares = 0.0d0
            do i=1,3
                sumSquares = sumSquares + vector(i)**2.0
            enddo
        end function sumSquares
        
        function dotProduct(v1,v2)
            double precision, dimension(3), intent(in) :: v1, v2
            double precision                           :: dotProduct
            integer                                    :: i
            dotProduct = 0.0d0
            do i=1,3
                dotProduct = dotProduct + v1(i) * v2(i)
            enddo
        end function dotProduct   

        function crossProduct(v1,v2)
            double precision, dimension(3), intent(in)  :: v1, v2
            double precision, dimension(3)              :: crossProduct
            crossProduct(1) = v1(2) * v2(3) - v1(3) * v2(2)
            crossProduct(2) = v1(3) * v2(1) - v1(1) * v2(3)
            crossProduct(3) = v1(1) * v2(2) - v1(2) * v2(1)
        end function crossProduct   
 
        function getBondLength(a1,a2)
            type(atom), intent(in)      :: a1, a2
            double precision            :: getBondLength
            getBondLength = norm(r12(a1,a2))
        end function getBondLength

        function getBondLengthSquared(a1,a2)
            type(atom), intent(in)      :: a1, a2
            double precision            :: getBondLengthSquared
            getBondLengthSquared = sumSquares(r12(a1,a2))
        end function getBondLengthSquared

        function getAngle(a1,a2,a3)
            ! Return angle of atoms in radians
            type(atom), intent(in)          :: a1, a2, a3
            double precision                :: getAngle, icos, isin
            double precision, dimension(3)  :: v21, v23 
            v21 = normalized(r12(a2,a1))
            v23 = normalized(r12(a2,a3))
            icos = dotProduct(v21,v23)
            isin = norm(crossProduct(v21,v23))
            getAngle = datan2(isin,icos)
            if(getAngle < 0.0d0) getAngle = getAngle + tau 
        end function getAngle
        
        function getDihedral(a1,a2,a3,a4)
            ! Return dihedral angle of atoms in radians
            type(atom), intent(in)          :: a1, a2, a3, a4
            double precision                :: getDihedral
            double precision, dimension(3)  :: n1, n2 
            n1 = crossProduct(r12(a1,a2),r12(a1,a3))
            n2 = crossProduct(r12(a2,a4),r12(a3,a4))
            getDihedral = dacos(dotProduct(n1,n2)/(norm(n1)*norm(n2)))
            if(dotProduct(r12(a3,a4),n1) < 0.0d0)then
                getDihedral = -getDihedral
            endif
        end function getDihedral

        function deg2rad(indeg)
            ! Convert degrees to radians
            double precision, intent(in)    :: indeg
            double precision                :: deg2rad
            deg2rad = indeg*d2r
        end function deg2rad
        
        function rad2deg(inrad)
            ! Convert radians to degrees
            double precision, intent(in)    :: inrad
            double precision                :: rad2deg
            rad2deg = inrad*r2d
        end function rad2deg

end module geometry
