! A molecule contains an array of atoms.
! Molecules also contain definitions of bonds, angles, and dihedrals.
! A connectivity definition (e.g. bond) contains
!   a) pointers to the atoms involved (a1,a2,a3)
!   b) parameters for each definition (e.g. bond constant: k)
!   c) integer array ids which holds atom numbers in the molecule

module molecule_class
    use atom_class
    use geometry

    implicit none
    private 
    public :: molecule, createmolecule, createRandomMolecule
    public :: bond, angle, dihedral

    type bond
        type(atom), pointer   :: a1 => null()
        type(atom), pointer   :: a2 => null()
        double precision      :: k=250.0d0, r0=1.5d0
        integer, dimension(2) :: ids
    end type bond
    type angle
        type(atom), pointer   :: a1 => null()
        type(atom), pointer   :: a2 => null()
        type(atom), pointer   :: a3 => null()
        double precision      :: k=150.0d0, r0=1.911d0
        integer, dimension(3) :: ids
    end type angle
    type dihedral
        type(atom), pointer   :: a1 => null()
        type(atom), pointer   :: a2 => null()
        type(atom), pointer   :: a3 => null()
        type(atom), pointer   :: a4 => null()
        double precision      :: height=2.0d0, frequency=3.0d0, phase=1.047d0
        integer, dimension(4) :: ids
    end type dihedral

    type molecule
        character(len=50)   :: molname
        integer             :: natoms=0
        integer             :: nbonds=0, nangles=0, ndihedrals=0 
        type(atom),     dimension(:), allocatable :: atoms
        type(bond),     dimension(:), allocatable :: bonds
        type(angle),    dimension(:), allocatable :: angles
        type(dihedral), dimension(:), allocatable :: dihedrals
        contains
            procedure :: createBond
            procedure :: createAngle
            procedure :: createDihedral
            procedure :: guessTopology
            procedure :: printTopMap
            procedure :: printInfo
            procedure :: centerOfMolecule
    end type molecule

    contains
        function createMolecule(molname,atoms)
            ! Return a molecule object given a name and array of atoms.
            character(len=*), intent(in)     :: molname
            type(atom), dimension(:), target :: atoms
            integer                          :: maxBonds, maxAngles, maxDihedrals, n
            type(molecule)                   :: createMolecule

            n = size(atoms) 
            createMolecule%molname  = molname
            createMolecule%natoms   = n
            createMolecule%atoms    = atoms
            maxBonds     = n*4
            maxAngles    = n*6
            maxdihedrals = n*4
            allocate(createMolecule%bonds(maxBonds), &
                     createMolecule%angles(maxAngles), &
                     createMolecule%dihedrals(maxDihedrals))
        end function createMolecule

        function createRandomMolecule(molname,natoms,cubeSize)
            ! Return a randomly generated molecule given name, number of atoms, and box size.
            character(len=*), intent(in)                :: molname
            integer, optional                           :: natoms
            integer                                     :: natoms_=25
            real, optional                              :: cubeSize
            real                                        :: cubeSize_=20.0
            integer                                     :: i, j            
            double precision, dimension(:), allocatable :: rands
            type(atom), dimension(:), allocatable       :: atoms
            type(molecule)                              :: createRandomMolecule

            !rands will hold random numbers for atomic coordinates
            !cubesize_ is a length of side of enclosing cube
            !rands(1:3) will be x,y,z which range from 0 to cubesize_/2
            !rands(4:6) will determine the sign of rands(1:3)
            !rands(7)   will determine the atom type
           
            if(present(cubeSize)) cubeSize_ = cubeSize
            if(present(natoms)) natoms_ = natoms
            
            allocate(atoms(natoms_),rands(7)) 
            call RANDOM_SEED() 
            do i=1,natoms_
                call RANDOM_NUMBER(rands)
                do j=1,3
                    if (rands(3+j) > 0.5d0)then
                        rands(j) = -rands(j)
                    endif
                enddo
                ! rands(1:3) are range [-1,1]
                rands = rands*cubeSize_/2.0
                ! rands(1:3) are range [-cubesize_/2 : cubesize_/2 ]
                ! rands(7) is range  [ 0 : cubesize_/2 ]
                if (rands(7) < cubeSize_*0.125d0)then
                    atoms(i) = createAtom("H",rands(1),rands(2),rands(3))
                elseif (rands(7) < cubeSize_*0.25d0)then
                    atoms(i) = createAtom("C",rands(1),rands(2),rands(3))
                elseif (rands(7) < cubeSize_*0.375d0)then
                    atoms(i) = createAtom("O",rands(1),rands(2),rands(3))
                else
                    atoms(i) = createAtom("Cl",rands(1),rands(2),rands(3))
                endif
            enddo
            createRandomMolecule = createMolecule(molname,atoms)
            deallocate(atoms,rands)
        end function createRandomMolecule

        subroutine createBond(self,a,b,k,r0)
            ! Assigns a new bond to this molecule.
            class(molecule), target    :: self
            type(bond)                 :: newbond
            integer                    :: a, b
            double precision, optional :: k, r0
            
            if(a < 0 .or. b < 0 .or. a > self%natoms .or. b > self%natoms)then
                print *, 'Bad bond definition.'
                stop
            endif
            self%nbonds = self%nbonds + 1
            newbond%a1 => self%atoms(a)
            newbond%a2 => self%atoms(b)
            newbond%ids(1) = a
            newbond%ids(2) = b
            if(present(k)) newbond%k = k
            if(present(r0)) newbond%r0 = r0
            self%bonds(self%nbonds) = newbond
        end subroutine createBond
 
        subroutine createAngle(self,a,b,c,k,r0)
            ! Assigns a new angle to this molecule.
            class(molecule), target    :: self
            type(angle)                :: newangle
            integer                    :: a, b, c
            double precision, optional :: k, r0

            self%nangles = self%nangles + 1
            newangle%a1 => self%atoms(a)
            newangle%a2 => self%atoms(b)
            newangle%a3 => self%atoms(c)
            newangle%ids(1) = a
            newangle%ids(2) = b
            newangle%ids(3) = c
            if(present(k)) newangle%k = k
            if(present(r0)) newangle%r0 = deg2rad(r0)
            self%angles(self%nangles) = newangle
        end subroutine createangle
  
        subroutine createdihedral(self,a,b,c,d,height,frequency,phase)
            ! Assigns a new dihedral to this molecule.
            class(molecule), target    :: self
            type(dihedral)             :: newdihedral
            integer                    :: a, b, c, d
            double precision, optional :: height, frequency, phase

            self%ndihedrals = self%ndihedrals + 1
            newdihedral%a1 => self%atoms(a)
            newdihedral%a2 => self%atoms(b)
            newdihedral%a3 => self%atoms(c)
            newdihedral%a4 => self%atoms(d)
            newdihedral%ids(1) = a
            newdihedral%ids(2) = b
            newdihedral%ids(3) = c
            newdihedral%ids(4) = d
            if(present(height)) newdihedral%height = height
            if(present(frequency)) newdihedral%frequency = frequency
            if(present(phase)) newdihedral%phase = phase
            self%dihedrals(self%ndihedrals) = newdihedral
        end subroutine createdihedral
 
        subroutine guessTopology(self, guessBonds, guessAngles, guessDihedrals)
            !Try to automatically assign bonds/angles/dihedrals based on atomic coordinates.
            !This is very basic at the moment and will need to be improved.
            class(molecule)                      :: self
            integer                              :: i, j, nbound, c
            logical, dimension(:,:), allocatable :: bondMatrix
            logical, intent(in), optional        :: guessBonds, guessAngles, guessDihedrals
            logical                              :: guessBonds_, guessAngles_, guessDihedrals_
            integer, dimension(:), allocatable   :: tbonds
            double precision                     :: td

            guessBonds_     = .false.
            guessAngles_    = .false.
            guessDihedrals_ = .false.
            if(present(guessBonds))     guessBonds_     = guessBonds
            if(present(guessAngles))    guessAngles_    = guessAngles
            if(present(guessDihedrals)) guessDihedrals_ = guessDihedrals

            allocate(bondMatrix(self%natoms,self%natoms))
            bondMatrix = .false.

            do i=1,self%natoms
                do j=i+1,self%natoms
                    if(self%atoms(i)%aid == 1 .and. self%atoms(j)%aid == 1) cycle !Ignore H-H
                    
                    td = getBondLengthSquared(self%atoms(i),self%atoms(j))
                    if(td > 4.25d0) cycle !Maximum cutoff for bonds

                    if(self%atoms(i)%aid == 1 .or. self%atoms(j)%aid == 1)then       ! Bonds to H
                        if(self%atoms(i)%aid >= 15 .or. self%atoms(j)%aid >= 15)then ! i.e. S-H
                            if(td <= 2.0d0) then
                                if(guessBonds_) call self%createBond(i,j,r0=1.9d0) 
                                bondMatrix(i,j) = .true. 
                                bondMatrix(j,i) = .true.
                            endif
                        else
                            if(td <= 1.6d0) then
                                if(guessBonds_) call self%createBond(i,j,r0=1.1d0) ! i.e. C-H, O-H, N-H
                                bondMatrix(i,j) = .true. 
                                bondMatrix(j,i) = .true. 
                            endif
                        endif

                    elseif(self%atoms(i)%aid == 8 .or. self%atoms(j)%aid == 8)then ! i.e. C-O, N-O
                        if(td <= 2.6d0) then
                            if(guessBonds_) call self%createBond(i,j,r0=1.4d0)
                            bondMatrix(i,j) = .true.
                            bondMatrix(j,i) = .true.
                        endif
                    elseif(self%atoms(i)%aid >= 15 .or. self%atoms(j)%aid >= 15)then ! i.e. C-S, C-P
                        if(td <= 4.0d0) then
                            if(guessBonds_) call self%createBond(i,j,r0=1.9d0)
                            bondMatrix(i,j) = .true.
                            bondMatrix(j,i) = .true.
                        endif
                    elseif(td <= 3.25d0) then ! i.e. C-C, C-N
                        if(guessBonds_) call self%createBond(i,j,r0=1.4d0)
                        bondMatrix(i,j) = .true. 
                        bondMatrix(j,i) = .true.
                    endif
                enddo
            enddo 

            if(guessAngles_) then 
                allocate(tbonds(10))
                do i=1,self%natoms
                    tbonds = 0 
                    nbound = 0
                    do j=1,self%natoms
                        if(i==j) cycle
                        if(bondMatrix(i,j))then
                            nbound = nbound + 1
                            tbonds(nbound) = j
                        endif 
                    enddo
                    if(self%atoms(i)%aid == 6)then
                        if(nbound==4)then !sp3-like
                            call self%createAngle(tbonds(1),i,tbonds(2),r0=109.5d0)
                            call self%createAngle(tbonds(1),i,tbonds(3),r0=109.5d0)
                            call self%createAngle(tbonds(1),i,tbonds(4),r0=109.5d0)
                            call self%createAngle(tbonds(2),i,tbonds(3),r0=109.5d0)
                            call self%createAngle(tbonds(2),i,tbonds(4),r0=109.5d0)
                            call self%createAngle(tbonds(3),i,tbonds(4),r0=109.5d0)
                        elseif(nbound==3)then !sp2-like
                            call self%createAngle(tbonds(1),i,tbonds(2),r0=120.0d0)
                            call self%createAngle(tbonds(1),i,tbonds(3),r0=120.0d0)
                            call self%createAngle(tbonds(2),i,tbonds(3),r0=120.0d0)
                        elseif(nbound==2)then !sp-like
                            call self%createAngle(tbonds(1),i,tbonds(2),r0=180.0d0)
                        endif
                    elseif(self%atoms(i)%aid == 7) then
                        if(nbound==4)then 
                            call self%createAngle(tbonds(1),i,tbonds(2),r0=109.5d0)
                            call self%createAngle(tbonds(1),i,tbonds(3),r0=109.5d0)
                            call self%createAngle(tbonds(1),i,tbonds(4),r0=109.5d0)
                            call self%createAngle(tbonds(2),i,tbonds(3),r0=109.5d0)
                            call self%createAngle(tbonds(2),i,tbonds(4),r0=109.5d0)
                            call self%createAngle(tbonds(3),i,tbonds(4),r0=109.5d0)
                        elseif(nbound==3)then 
                            call self%createAngle(tbonds(1),i,tbonds(2),r0=120.0d0)
                            call self%createAngle(tbonds(1),i,tbonds(3),r0=120.0d0)
                            call self%createAngle(tbonds(2),i,tbonds(3),r0=120.0d0)
                        elseif(nbound==2)then 
                            call self%createAngle(tbonds(1),i,tbonds(2),r0=120.0d0)
                        endif
                    elseif(self%atoms(i)%aid == 8) then
                        if(nbound==2)then 
                            call self%createAngle(tbonds(1),i,tbonds(2),r0=120.0d0)
                        endif
                    else
                        if(nbound==2)then 
                            call self%createAngle(tbonds(1),i,tbonds(2),r0=120.0d0)
                        endif
                    endif    
                enddo
            endif

            if(allocated(bondMatrix)) deallocate(bondMatrix)
            if(allocated(tbonds)) deallocate(tbonds)

        end subroutine guessTopology

        subroutine printTopMap(self)
            !Print xyz file with H at midpoints of bonds (for visualization/debugging)
            class(molecule)                         :: self
            double precision, dimension(3)          :: tpoint
            integer                                 :: i
            
            open(22,file="topmap.xyz",action="WRITE",status="REPLACE")
            write(22,*) "TOPMAP"
            write(22,*) self%nbonds
            do i=1,self%nbonds
                tpoint = self%bonds(i)%a1%coords + r12(self%bonds(i)%a1,self%bonds(i)%a2)/2.0d0
                write(22,*) "H ", tpoint
            enddo
            close(22)
        end subroutine printTopMap

        subroutine printInfo(self)
            ! Prints information about this molecule.
            class(molecule) :: self
            integer         :: i

            write(*,'(A50)') self%molname
            write(*,'(A6,1X,I5)') 'Atoms:', self%natoms
            do i = 1,self%natoms
                call self%atoms(i)%printInfo()
            enddo
            do i=1,self%nbonds
                write(*,'(A5,I4,A9,I4,A5,I4,2(A4,F8.3))') 'Bond ', &
                         i,' between ',                 &
                         self%bonds(i)%ids(1),' and ',  &
                         self%bonds(i)%ids(2),          &
                         '  k=',self%bonds(i)%k,        &
                         ' r0=',self%bonds(i)%r0
            enddo
            do i=1,self%nangles
                write(*,'(A6,I4,A9,2(I4,A5),I4,2(A4,F8.3))') 'Angle ', &
                         i,' between ',  &
                         self%angles(i)%ids(1),' and ', &
                         self%angles(i)%ids(2),' and ', &
                         self%angles(i)%ids(3),         &
                         '  k=', self%angles(i)%k,      &
                         ' r0=', rad2deg(self%angles(i)%r0)
            enddo
            do i=1,self%ndihedrals
                write(*,'(A9,I4,A9,3(I4,A5),I4,3(A4,F8.3))') 'Dihedral ', &
                         i,' between ', &
                         self%dihedrals(i)%ids(1),' and ',      &
                         self%dihedrals(i)%ids(2),' and ',      &
                         self%dihedrals(i)%ids(3),' and ',      &
                         self%dihedrals(i)%ids(4),              &
                         '  h=', self%dihedrals(i)%height,      &
                         '  f=', self%dihedrals(i)%frequency,   &
                         '  p=', rad2deg(self%dihedrals(i)%phase)       
            enddo

        end subroutine printInfo

        function centerOfMolecule(self)
            ! Return center of molcule as vector
            class(molecule)                 :: self
            double precision, dimension(3)  :: centerOfMolecule 
            integer                         :: i

            centerOfMolecule = 0.0d0
            do i=1,self%natoms
                centerOfMolecule = centerOfMolecule + self%atoms(i)%coords
            enddo  
            centerOfMolecule = centerOfMolecule / self%natoms
        end function centerOfMolecule

    end module molecule_class
