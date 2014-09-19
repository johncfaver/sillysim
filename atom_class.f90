module atom_class
    implicit none
    private
    
    public ::  atom, createAtom

    character(len=2), dimension(92) :: elements = [ &
                      " H","He","Li","Be"," B"," C"," N"," O"," F","Ne", &
                      "Na","Mg","Al","Si"," P"," S","Cl","Ar"," K","Ca", &
                      "Sc","Ti"," V","Cr","Mn","Fe","Co","Ni","Cu","Zn", &
                      "Ga","Ge","As","Se","Br","Kr","Rb","Sr"," Y","Zr", &
                      "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn", &
                      "Sb","Te"," I","Xe","Cs","Ba","La","Ce","Pr","Nd", &
                      "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb", &
                      "Lu","Hf","Ta"," W","Re","Os","Ir","Pt","Au","Hg", &
                      "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th", &
                      "Pa"," U"]
 

    type atom
        character(len=2)                :: element
        integer                         :: aid
        double precision                :: x, y, z
        double precision, dimension(3)  :: coords
        real                            :: mass, LJe, LJs, charge
        contains
            procedure :: printInfo
            procedure :: translate
    end type atom

    contains
        function createAtom(e,x,y,z,LJe,LJs,charge,mass)
            ! Return atom object of element e placed at x,y,z
            type(atom)                      :: createAtom
            character(len=*), intent(in)    :: e 
            double precision, intent(in)    :: x, y, z
            real, intent(in), optional      :: LJe,LJs,charge,mass
            createAtom%element = e
            createAtom%x = x
            createAtom%y = y
            createAtom%z = z
            createAtom%coords = [x, y, z]
            createAtom%aid = getAid(e)
            if(present(mass)) then
                createAtom%mass = mass
            else 
                createAtom%mass = getmass(e)
            endif
            if(present(LJe)) then
                createAtom%LJe = LJe
            else 
                createAtom%LJe = getLJe(e)
            endif
            if(present(LJs)) then
                createAtom%LJs = LJs
            else
                createAtom%LJs = getLJs(e)
            endif
            if(present(charge)) then
                createAtom%charge = charge
            else
                createAtom%charge = getcharge(e)
            endif
        end function createAtom
       
        function getMass(e)
            ! Return mass of atom of element e
            real             :: getMass
            character(len=*) :: e
            select case (e)
                case("H")
                    getMass=1.00794
                case("C")
                    getMass=12.0107
                case("O")
                    getMass=15.9994
                case("Cl")
                    getMass=35.4527
                case default
                    getMass=200000.0
            end select
        end function getMass 

        function getLJe(e)
            ! Lennard-Jones epsilon values
            real             :: getLJe
            character(len=*) :: e
            select case (e)
                case("H")
                    getLJe=1.0
                case("C")
                    getLJe=2.0
                case("O")
                    getLJe=2.5
                case("Cl")
                    getLJe=3.0
                case default
                    getLJe=2.0
            end select
        end function getLJe
        
        function getLJs(e)
            ! Lennard-Jones sigma values
            real             :: getLJs
            character(len=*) :: e
            select case (e)
                case("H")
                    getLJs=1.0
                case("C")
                    getLJs=3.0
                case("O")
                    getLJs=3.5
                case("Cl")
                    getLJs=4.0
                case default
                    getLJs=3.0
            end select
        end function getLJs
  
        function getCharge(e)
            ! Atomic charge
            real             :: getCharge
            character(len=*) :: e
            select case (e)
                case("H")
                    getcharge=0.3
                case("C")
                    getcharge=-0.3
                case("O")
                    getcharge=-0.5
                case("Cl")
                    getcharge=-0.5
                case default
                    getcharge=0.0
            end select
        end function getcharge
 
        function getAid(e)
            ! Atomic number
            real             :: getAid
            character(len=*) :: e
            character(len=2) :: e2
            integer          :: i
            write(e2,'(A2)') e
            do i=1,size(elements)
                if(e2 == elements(i)) then
                    getAid = i
                    exit
                endif
            enddo
        end function getAid

        subroutine printInfo(self)
            class(atom) :: self
            write(*,'(A2,3(1X,F15.7))') &
                self%element, self%x, self%y, self%z
        end subroutine printInfo

        subroutine translate(self,v)
            ! Translate this atom by vector v
            class(atom)                                :: self
            double precision, intent(in), dimension(3) :: v
            self%x = self%x + v(1)
            self%y = self%y + v(2)
            self%z = self%z + v(3)
            self%coords = [self%x, self%y, self%z]
        end subroutine translate        

end module atom_class
