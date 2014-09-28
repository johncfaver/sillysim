! A simulation contains an instance of a potential.
!
! A potential contains a set of potential terms:
!   twoBodyTerms:   an array of type (twoBodyPotentialTerm) for atom-atom two body terms
!   bond:           an instance of type (bondPotentialTerm) for bonds
!   angle:          an instance of type (anglePotentialTerm) for angles
!   dihedral:       an instance of type (dihedralPotentialTerm) for dihedrals
!
! Each potential term contains:
!   energy:         a pointer to function for calculating the energy for the term
!   gradient:       a pointer to function for calculating the 3D derivative dE/dr
!   gradient2:      a pointer to function for calculating the 3D 2nd derivative d2E/dr2
!   inUse:          a boolean value for switching the term on or off
!

module potential_class
    use molecule_class
    use atom_class
    use geometry

    implicit none
    private
    
    public :: potential, defaultPotentialFunction
    
    integer(1), parameter       :: maxTwoBodyTerms = 3
    double precision, parameter :: gconstant = 0.05d0

    interface 
        function twoBodyTerm(a1,a2)
            import                          :: atom
            type(atom), intent(in)          :: a1, a2
            double precision                :: twoBodyTerm
        end function twoBodyTerm
        function twoBodyTerm3d(a1,a2)
             import                         :: atom
             type(atom), intent(in)         :: a1, a2
             double precision, dimension(3) :: twoBodyTerm3d
        end function twoBodyTerm3d
        function bondTerm(b)
            import                          :: bond
            type(bond), intent(in)          :: b
            double precision                :: bondTerm
        end function bondTerm
        function bondTerm3d(b,i)
            import                          :: bond
            type(bond), intent(in)          :: b
            integer, intent(in)             :: i
            double precision, dimension(3)  :: bondTerm3d
        end function bondTerm3d
        function angleTerm(a)
            import                          :: angle
            type(angle), intent(in)         :: a
            double precision                :: angleTerm
        end function angleTerm
        function angleTerm3d(a,i)
            import                          :: angle
            type(angle), intent(in)         :: a
            integer, intent(in)             :: i
            double precision, dimension(3)  :: angleTerm3d
        end function angleTerm3d
        function dihedralTerm(d)
            import                          :: dihedral
            type(dihedral), intent(in)      :: d
            double precision                :: dihedralTerm
        end function dihedralTerm
        function dihedralTerm3d(d,i)
            import                          :: dihedral
            type(dihedral), intent(in)      :: d
            integer, intent(in)             :: i
            double precision, dimension(3)  :: dihedralTerm3d
        end function dihedralTerm3d
    end interface    

    type genericPotentialTerm
        logical :: inUse = .false.
        integer :: potentialType
    end type genericPotentialTerm

    type, extends(genericPotentialTerm) :: twoBodyPotentialTerm
        procedure (twoBodyTerm), pointer, nopass    :: energy => null()
        procedure (twoBodyTerm3d), pointer, nopass  :: gradient => null(), gradient2 => null()
    end type twoBodyPotentialTerm

    type, extends(genericPotentialTerm) :: bondPotentialTerm
        procedure (bondTerm), pointer, nopass       :: energy => null()
        procedure (bondTerm3d), pointer, nopass     :: gradient => null(), gradient2 => null()
    end type bondPotentialTerm

    type, extends(genericPotentialTerm) :: anglePotentialTerm
        procedure (angleTerm), pointer, nopass      :: energy => null()
        procedure (angleTerm3d), pointer, nopass    :: gradient => null(), gradient2 => null()
    end type anglePotentialTerm

    type, extends(genericPotentialTerm) :: dihedralPotentialTerm
        procedure (dihedralTerm), pointer, nopass   :: energy => null()
        procedure (dihedralTerm3d), pointer, nopass :: gradient => null(), gradient2 => null()
    end type dihedralPotentialTerm

    type potential
        type(twoBodyPotentialTerm), dimension(maxTwoBodyTerms) :: twoBodyTerms
        type(bondPotentialTerm)                                :: bond
        type(anglePotentialTerm)                               :: angle 
        type(dihedralPotentialTerm)                            :: dihedral
        contains
           procedure :: getMaxTwoBodyTerms
           procedure :: addTerm  
           procedure :: totalEnergy
           procedure :: printInfo
    end type potential

    contains
        function defaultPotentialFunction()
            ! Returns a default potential function
            ! Add bond, angle, dihedral, LJ, Coulomb
            type(potential) :: defaultPotentialFunction
            call defaultPotentialFunction%addTerm(1)
            call defaultPotentialFunction%addTerm(2)
            call defaultPotentialFunction%addTerm(3) 
            call defaultPotentialFunction%addTerm(4)
            call defaultPotentialFunction%addTerm(5)
        end function defaultPotentialFunction

        function getMaxTwoBodyTerms(self)
            class(potential) :: self
            integer          :: getMaxTwoBodyTerms
            getMaxTwoBodyTerms = maxTwoBodyTerms
        end function getMaxTwoBodyTerms

        subroutine printInfo(self)
            ! Prints information about this potential function.
            class(potential), intent(in) :: self
            integer                      :: i

            do i=1,self%getMaxTwoBodyTerms()
                if(self%twoBodyTerms(i)%inUse)then
                    write(*,'(A25,I2,A9,A30)') 'Two-body potential term #', i, ' in use: ', &
                        potentialTypeName(self%twoBodyTerms(i)%potentialType)
                else
                    exit
                endif
            enddo
            if(self%bond%inUse)     write(*,'(A22)') 'Bond potential in use.'
            if(self%angle%inUse)    write(*,'(A23)') 'Angle potential in use.'
            if(self%dihedral%inUse) write(*,'(A26)') 'Dihedral potential in use.'
        end subroutine printInfo 

        function totalEnergy(self,mol)
            !Compute and return total potential energy of molecule mol.
            class(potential)            :: self
            class(molecule), intent(in) :: mol
            integer                     :: i, j, k
            double precision            :: totalEnergy
            
            totalEnergy = 0.0
            do i=1,self%getMaxTwoBodyTerms()
                if(self%twoBodyTerms(i)%inUse)then
                    do j=1,mol%natoms
                        do k=j+1,mol%natoms
                            totalEnergy = totalEnergy + &
                            self%twoBodyTerms(i)%energy(mol%atoms(j),mol%atoms(k))
                        enddo
                    enddo 
                endif
            enddo 
            if(self%bond%inUse)then
                do i=1,mol%nbonds
                    totalEnergy = totalEnergy + &
                    self%bond%energy(mol%bonds(i))
                enddo
            endif
            if(self%angle%inUse)then
                do i=1,mol%nangles
                    totalEnergy = totalEnergy + &
                    self%angle%energy(mol%angles(i))
                enddo
            endif
            if(self%dihedral%inUse)then
                do i=1,mol%ndihedrals
                    totalEnergy = totalEnergy + &
                    self%dihedral%energy(mol%dihedrals(i))
                enddo
            endif
        end function totalEnergy

        subroutine addTerm(self,term)
            ! term is an integer code for a potential term type.
            class(potential)    :: self
            integer, intent(in) :: term
            integer             :: i
           
            if(term == 1)then
                self%bond%energy => bondEnergy
                self%bond%gradient => bondGradient3d
                self%bond%gradient2 => bondGradient3d2
                self%bond%inUse = .true.
                self%bond%potentialType = 1
            elseif(term == 2) then
                self%angle%energy => angleEnergy
                self%angle%gradient => angleGradient3d
                self%angle%gradient2 => angleGradient3d2
                self%angle%inUse = .true.
                self%angle%potentialType = 2
            elseif(term == 3) then
                self%dihedral%energy => dihedralEnergy
                self%dihedral%gradient => dihedralGradient3d
                self%dihedral%inUse = .true.
                self%dihedral%potentialType = 3
            endif
            if( term >= 4 )then
                do i=1,self%getMaxTwoBodyTerms()
                    if(.not. self%twoBodyTerms(i)%inUse) then
                        select case (term)
                            case(4)
                                self%twoBodyTerms(i)%energy => lennardJones
                                self%twoBodyTerms(i)%gradient => lennardJonesGradient3d
                            case(5)
                                self%twoBodyTerms(i)%energy => coulomb
                                self%twoBodyTerms(i)%gradient => coulombGradient3d
                            case(6)
                                self%twoBodyTerms(i)%energy => gravity
                                self%twoBodyTerms(i)%gradient => gravityGradient3d
                            case default
                                print *, 'Error: Bad two-body potential term requested.'
                                stop
                        end select
                        self%twoBodyTerms(i)%inUse = .true.
                        self%twoBodyTerms(i)%potentialType = term
                        exit
                    endif
                enddo
            endif
        end subroutine addTerm

        function potentialTypeName(i)
            ! Return the name of the potential given the integer code.
            integer, intent(in) :: i
            character(len=30)   :: potentialTypeName
       
            select case(i)              !123456789012345678901234567890
                case(1)
                    potentialTypeName = "Bond"
                case(2)
                    potentialTypeName = "Angle"
                case(3)
                    potentialTypeName = "Dihedral"
                case(4)
                    potentialTypeName = "LJ12-6"
                case(5)
                    potentialTypeName = "Point charge Electrostatic"
                case(6)
                    potentialTypeName = "Gravity"
            end select
        end function potentialTypeName

        function bondEnergy(b)
            ! U = k/2*[r-r0]^2
            type(bond), intent(in)  :: b
            double precision        :: bondEnergy
            bondEnergy = b%k * 0.5d0 * (norm(r12(b%a1,b%a2))-b%r0)**2
        end function bondEnergy
      
        function bondGradient(b)
            ! du/dr = k*[r-r0]
            type(bond), intent(in)  :: b
            double precision        :: bondGradient
            bondGradient = b%k * (norm(r12(b%a1,b%a2))-b%r0)
        end function bondGradient

        function bondGradient2(b)
            ! d2u/dr2 = k
            type(bond), intent(in)  :: b
            double precision        :: bondGradient2
            bondGradient2 = b%k 
        end function bondGradient2

        function bondGradient3d(b,i)
            ! return gradient on atom i (a1 or a2) due to bond b
            type(bond), intent(in)         :: b
            integer, intent(in)            :: i
            double precision, dimension(3) :: bondGradient3d
            bondGradient3d = bondGradient(b) * normalized(r12(b%a2,b%a1))
            if(i == 2) then
                bondGradient3d = -bondGradient3d
            endif
        end function bondGradient3d
  
        function bondGradient3d2(b,i)
            ! d2u/dr2 = k
            type(bond), intent(in)         :: b
            integer, intent(in)            :: i
            double precision, dimension(3) :: bondGradient3d2
            bondGradient3d2 = b%k * normalized(r12(b%a2,b%a1))
            if(i == 2) then
                bondGradient3d2 = -bondGradient3d2
            endif
        end function bondGradient3d2 

        function angleEnergy(a)
            ! U = k/2*[r-r0]^2
            type(angle), intent(in) :: a
            double precision        :: angleEnergy
            angleEnergy = a%k * 0.5d0 * (getAngle(a%a1,a%a2,a%a3)-a%r0)**2
        end function angleEnergy
      
        function angleGradient(a)
            ! du/dr = k*[r-r0]
            type(angle), intent(in) :: a
            double precision        :: angleGradient
            angleGradient = a%k * (getAngle(a%a1,a%a2,a%a3)-a%r0)
        end function angleGradient

        function angleGradient2(a)
            ! d2u/dr2 = k
            type(angle), intent(in) :: a
            double precision        :: angleGradient2
            angleGradient2 = a%k 
        end function angleGradient2

        function angleGradient3d(a,i)
            ! return gradient on atom i (a1, a2, or a3) due to angle a
            ! define orthonormal basis b1, b2:
            !      b1 is the a1->a2 vector
            !      b2 is the orthogonal vector in the plane of the angle. 
            type(angle), intent(in)        :: a
            integer, intent(in)            :: i
            double precision               :: magnitude, theta
            double precision, dimension(3) :: angleGradient3d
            double precision, dimension(3) :: b1,b2,u23
            theta = getAngle(a%a1,a%a2,a%a3)
            magnitude = angleGradient(a)
            b1 = normalized(r12(a%a1,a%a2))
            u23 = normalized(r12(a%a2,a%a3))
            b2 = u23 - dotProduct(u23,b1)*b1
            if(i == 1)then
                angleGradient3d = magnitude*(-b2)
            elseif(i == 2)then
                angleGradient3d = magnitude*(-dsin(theta)*b1+(b2-dcos(theta)*b2))
            elseif(i == 3)then      
                angleGradient3d = magnitude*(dsin(theta)*b1+dcos(theta)*b2)
            endif
        end function angleGradient3d

        function angleGradient3d2(a,i)
            type(angle), intent(in)        :: a
            integer, intent(in)            :: i
            double precision               :: magnitude, theta
            double precision, dimension(3) :: angleGradient3d2
            double precision, dimension(3) :: b1,b2,u23
            theta = getAngle(a%a1,a%a2,a%a3)
            magnitude = a%k
            b1 = normalized(r12(a%a1,a%a2))
            u23 = normalized(r12(a%a2,a%a3))
            b2 = u23 - dotProduct(u23,b1)*b1
            if(i == 1)then
                angleGradient3d2 = magnitude*(-b2)
            elseif(i == 2)then
                angleGradient3d2 = magnitude*(-dsin(theta)*b1+(b2-dcos(theta)*b2))
            elseif(i == 3)then      
                angleGradient3d2 = magnitude*(dsin(theta)*b1+dcos(theta)*b2)
            endif
        end function angleGradient3d2

        function dihedralEnergy(d)
            ! U = h/2 * [ 1 + cos(t*f-p) ]
            type(dihedral), intent(in)  :: d
            double precision            :: dihedralEnergy
            dihedralEnergy = d%height*0.5d0*(1.0d0+dcos(getDihedral(d%a1,d%a2,d%a3,d%a4)*d%frequency - d%phase))
        end function dihedralEnergy
      
        function dihedralGradient(d)
            ! dU/dx = h/2 * [ -sin(t*f-p)*f ]
            type(dihedral), intent(in)  :: d
            double precision            :: dihedralGradient
            dihedralGradient = d%height*0.5d0*(-dsin(getDihedral(d%a1,d%a2,d%a3,d%a4)*d%frequency - d%phase)*d%frequency)
        end function dihedralGradient

        function dihedralGradient3d(d,i)
            ! return gradient on atom i (a1, a2, a3, or a4) due to dihedral d
            type(dihedral), intent(in)      :: d
            integer, intent(in)             :: i
            double precision, dimension(3)  :: dihedralGradient3d,n1,n2,u12,u23,u34
            double precision                :: l12, l23, l34, a123, a234
            double precision, dimension(3)  :: a,b, f1, f2, f3, f4
            double precision                :: rs2j,rs2k, rrcj, rrck

            l12 = norm(r12(d%a1,d%a2))
            l23 = norm(r12(d%a2,d%a3))
            l34 = norm(r12(d%a3,d%a4))
            
            u12 = normalized(r12(d%a1,d%a2))
            u23 = normalized(r12(d%a2,d%a3))
            u34 = normalized(r12(d%a3,d%a4))

            a123 = getAngle(d%a1,d%a2,d%a3)
            a234 = getAngle(d%a2,d%a3,d%a4)
            
            rs2j = 1.0d0/(l12*dsin(a123)**2)
            rs2k = 1.0d0/(l34*dsin(a234)**2)
            rrcj = -l12/l23*dcos(a123)
            rrck = -l34/l23*dcos(a234)

            a = crossProduct(u12,u23)
            b = crossProduct(u23,u34)
           
            f1 = a * (-rs2j)
            f4 = b * rs2k
            
            a = f1 * (rrcj-1.0d0)
            b = f4 * rrck
            
            f2 = a - b

            f4 = -(f1+f2+f3)

            select case(i)
                case(1)
                    dihedralGradient3d = f1
                case(2)
                    dihedralGradient3d = f2
                case(3)
                    dihedralGradient3d = f3
                case(4)
                    dihedralGradient3d = f4
                case default
                    print *, 'ERROR IN DIHEDRAL GRADIENT'
            end select
        end function dihedralGradient3d

        function lennardJones(a1,a2)
            ! U =  e [ (s/r)**12 - (s/r)**6 ] 
            type(atom), intent(in)   :: a1, a2
            double precision         :: term, e, s, lennardJones
            s = (a1%ljs + a2%ljs)*0.5d0
            e = sqrt(a1%lje*a2%lje)
            term = (s/norm(r12(a1,a2)))**6
            lennardJones = e*(term**2 - term) 
        end function lennardJones

        function lennardJonesGradient(a1,a2)
            ! dU/dr = -6*e/r [ 2*(s/r)^12 - (s/r)^6 ]
            type(atom), intent(in)  :: a1, a2
            double precision        :: term, e, s, r, lennardJonesGradient
            s = (a1%ljs + a2%ljs)*0.5d0
            e = sqrt(a1%lje*a2%lje)
            r = norm(r12(a1,a2))
            term = (s/r)**6
            lennardJonesGradient = -6.0d0*e/r*(2.0d0*term**2 - term)
        end function lennardJonesGradient

        function lennardJonesGradient3d(a1,a2)
            ! return gradient in direction to atom 2 from atom 1
            type(atom), intent(in)          :: a1, a2
            double precision                :: e, s
            double precision, dimension(3)  :: lennardJonesGradient3d
            s = (a1%ljs + a2%ljs)*0.5
            e = sqrt(a1%lje * a2%lje)
            lennardJonesGradient3d = normalized(r12(a2,a1)) * lennardJonesGradient(a1,a2)
        end function lennardJonesGradient3d

        function coulomb(a1,a2)
            ! E = q1 * q2 / r
            type(atom), intent(in)  :: a1, a2
            double precision        :: coulomb
            coulomb = a1%charge * a2%charge / norm(r12(a1,a2))
        end function coulomb

        function coulombGradient(a1,a2)
            ! dE/dr = -(q1 * q2 / r^2)
            type(atom), intent(in)  :: a1, a2
            double precision        :: coulombGradient
            coulombGradient = -(a1%charge * a2%charge) / sumSquares(r12(a1,a2))
        end function coulombGradient

        function coulombGradient3d(a1,a2)
            ! return gradient in direction to atom 2 from atom 1
            type(atom), intent(in)          :: a1, a2
            double precision, dimension(3)  :: coulombGradient3d
            coulombGradient3d = normalized(r12(a2,a1)) * coulombGradient(a1,a2)
        end function coulombGradient3d

        function gravity(a1,a2)
            ! E = q1 * q2 / r
            type(atom), intent(in)  :: a1, a2
            double precision        :: gravity
            gravity = - gconstant * (a1%mass * a2%mass) / norm(r12(a1,a2))
        end function gravity

        function gravityGradient(a1,a2)
            ! dE/dr = -(q1 * q2 / r^2)
            type(atom), intent(in)  :: a1, a2
            double precision        :: gravityGradient
            gravityGradient = gconstant * a1%mass * a2%mass / sumSquares(r12(a1,a2))
        end function gravityGradient

        function gravityGradient3d(a1,a2)
            ! return gradient in direction to atom 2 from atom 1
            type(atom), intent(in)          :: a1, a2
            double precision, dimension(3)  :: gravityGradient3d
            gravityGradient3d = normalized(r12(a2,a1)) * gravityGradient(a1,a2)
        end function gravityGradient3d

end module potential_class
