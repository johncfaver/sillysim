! Simulator object
! Given a molecule, the simulator object can:
!       Calculate potential energy
!       Optimize the structure
!       Perform dynamics simulations
! 
module simulator_class
    use atom_class
    use molecule_class
    use potential_class
    use file_handler_class
    use geometry
    implicit none
    private

    public :: simulator, createSimulator
    
    type simulator
        type(molecule), pointer     :: simulatedMolecule => null()
        type(potential)             :: potentialFunction          
        character(len=80)           :: outputfile="out.pdb"
        integer                     :: nsteps=1, outputFrequency=0, nthreads=1
        real                        :: dt=0.01
        double precision, dimension(:,:), allocatable   :: forces, velocities, accelerations
        double precision, dimension(:,:,:), allocatable :: forceMatrix
        contains
            procedure :: simulate
            procedure :: minimize
            procedure :: potential_energy
            procedure :: kinetic_energy
            procedure :: calculate_forces
            procedure :: calculate_accelerations
            procedure :: update_velocities
            procedure :: update_positions
    end type simulator 

    contains
        function createSimulator(simulatedMolecule, potentialFunction, &
                                 outputfile, nthreads, nsteps, &
                                 outputFrequency, dt)
            !Simulator constructor. 
            !Requires molecule object.
            !Other variables are optional.
            type(simulator)                      :: createSimulator
            type(molecule), target               :: simulatedMolecule
            type(potential), optional            :: potentialFunction
            character(len=*), optional           :: outputfile
            integer, optional                    :: nsteps, outputFrequency, nthreads
            real, optional                       :: dt
            integer                              :: ios
            
            createSimulator%simulatedMolecule => simulatedMolecule
            if(present(potentialfunction)) then
                createSimulator%potentialFunction = potentialFunction
            else
                createSimulator%potentialFunction = defaultPotentialFunction()
            endif
            if(present(nsteps)) createSimulator%nsteps = nsteps
            if(present(dt)) createSimulator%dt = dt
            if(present(nthreads)) createSimulator%nthreads = nthreads
            if(present(outputfile)) createSimulator%outputfile = outputfile
            if(present(outputFrequency)) createSimulator%outputFrequency = outputFrequency
            if(createSimulator%outputFrequency == 0) createSimulator%outputFrequency = nsteps 

            open(20, file=createSimulator%outputfile, status="REPLACE", iostat=ios)
            if(ios < 0) print *, 'Error opening output file: ', createSimulator%outputfile
            close(20, iostat=ios)
        end function createSimulator

        subroutine minimize(self, nsteps, threshold)
            !Minimize molecule by stepping in the 
            !direction of the negative energy gradient by dt
            !for nsteps steps. Stop if energy difference < threshold
            class(simulator)                :: self
            integer                         :: natoms, istep, i
            integer, optional               :: nsteps
            double precision, dimension(3)  :: ivec
            double precision                :: le, pe, threshold_=1.0d-5
            double precision, optional      :: threshold
            real                            :: t_start, t_finish 

            natoms = self%simulatedMolecule%natoms
            print *, "Starting minimization of ", natoms , &
                     " atoms using ", self%nthreads, " threads."
            
            if(.not. allocated(self%forces)) allocate(self%forces(natoms,3))
            if(.not. allocated(self%forceMatrix)) allocate(self%forceMatrix(natoms,natoms,3))
          
            if(present(threshold)) threshold_ = threshold
            le = 1.0d100 
            if(present(nsteps)) self%nsteps = nsteps

            print *, 'INITIAL POTENTIAL ENERGY = ',self%potential_energy()
            call writeMolecule(self%outputfile, self%simulatedMolecule)
            call cpu_time(t_start)
            do istep = 1, self%nsteps
                call self%calculate_forces()
                !$OMP  PARALLEL DO SCHEDULE(STATIC) &
                !$OMP& DEFAULT(SHARED) PRIVATE(i,ivec) &
                !$OMP& NUM_THREADS(self%nthreads)
                do i=1,self%simulatedMolecule%natoms
                    ivec = normalized(self%forces(i,:))*self%dt
                    call self%simulatedMolecule%atoms(i)%translate(ivec)
                enddo
                !$OMP END PARALLEL DO 
                pe = self%potential_energy()
                if(dabs(le - pe) <= threshold_) then
                    exit
                else
                    if (mod(istep,self%outputFrequency) == 0)then
                        call writeMolecule(self%outputfile,self%simulatedMolecule,append=.true.)
                        print *, 'E= ',pe," at step ",istep, ' dE=',(le-pe)
                    endif 
                    le = pe
                endif        
            enddo
            call cpu_time(t_finish)
            call writeMolecule(self%outputfile,self%simulatedMolecule,append=.true.)
            print *, 'Minimization finished in ',istep,' steps (', (t_finish-t_start)/self%nthreads,'s).'
        end subroutine minimize

        subroutine simulate(self, nsteps)
            class(simulator)    :: self
            double precision    :: pe, ke, te
            integer             :: istep, natoms
            integer, optional   :: nsteps
            real                :: t_start, t_finish
            
            natoms = self%simulatedMolecule%natoms
            print *, "Starting simulation of ", natoms , &
                     " atoms using ", self%nthreads, " threads."
            
            if(.not. allocated(self%forces)) allocate(self%forces(natoms,3))
            if(.not. allocated(self%velocities)) allocate(self%velocities(natoms,3))
            if(.not. allocated(self%accelerations)) allocate(self%accelerations(natoms,3))
            if(.not. allocated(self%forceMatrix)) allocate(self%forceMatrix(natoms,natoms,3))

            if(present(nsteps)) self%nsteps = nsteps
            self%velocities = 0.0d0
            
            print *, 'INITIAL POTENTIAL ENERGY = ',self%potential_energy()

            call cpu_time(t_start)
            do istep = 1,self%nsteps
                call self%calculate_forces()
                call self%calculate_accelerations()
                call self%update_velocities()
                call self%update_positions()

                if (mod(istep,self%outputFrequency) == 0)then
                    call writeMolecule(self%outputfile,self%simulatedMolecule,append=.true.)
                    pe = self%potential_energy()
                    ke = self%kinetic_energy()
                    te = pe+ke
                    print *, 'PE= ',pe, &
                             "KE= ",ke, &
                             "TE= ",te, &
                             "AT STEP ", istep
                endif
            enddo
            call cpu_time(t_finish)
            print *, 'Simulation finished in ', (t_finish-t_start)/self%nthreads,'s.'
        end subroutine simulate 

        function potential_energy(self)
            class(simulator), target    :: self
            double precision            :: potential_energy
            type(molecule), pointer     :: imol => null()
            type(potential), pointer    :: ipot => null()
            integer                     :: i,j,k,x,ncells,nterms

            imol => self%simulatedMolecule
            ipot => self%potentialFunction
            ncells = (imol%natoms**2 - imol%natoms)*0.5
            nterms = ipot%getMaxTwoBodyTerms()
            potential_energy = 0.0d0

            !$OMP  PARALLEL DO SCHEDULE(STATIC) &
            !$OMP& DEFAULT(SHARED) PRIVATE(x,i,j,k) &
            !$OMP& NUM_THREADS(self%nthreads) REDUCTION(+:potential_energy)
            do x=1,ncells
                i = imol%natoms - nint(sqrt(2.0*(1+ncells-x)))
                j = mod(x+i*(i+1)/2-1,imol%natoms)+1
                do k=1,nterms
                    if(ipot%twoBodyTerms(k)%inUse)then
                        potential_energy = potential_energy + &
                            ipot%twoBodyTerms(k)%energy(imol%atoms(i),imol%atoms(j))
                    else
                        exit
                    endif
                enddo
            enddo
           !$OMP END PARALLEL DO
            if(ipot%bond%inUse)then
                do i=1,imol%nbonds
                    potential_energy = potential_energy + &
                        ipot%bond%energy(imol%bonds(i))
                enddo
            endif
            if(ipot%angle%inUse)then
                do i=1,imol%nangles
                    potential_energy = potential_energy + &
                        ipot%angle%energy(imol%angles(i))
                enddo
            endif
            if(ipot%dihedral%inUse)then
                do i=1,imol%ndihedrals
                    potential_energy = potential_energy + &
                        ipot%dihedral%energy(imol%dihedrals(i))
                enddo
            endif
        end function potential_energy 

        function kinetic_energy(self)
            class(simulator)  :: self
            double precision :: kinetic_energy
            integer          :: i

            kinetic_energy = 0.0d0
            !$OMP PARALLEL DO SCHEDULE(STATIC) &
            !$OMP DEFAULT(SHARED) PRIVATE(i) &
            !$OMP NUM_THREADS(self%nthreads) REDUCTION(+:kinetic_energy)
            do i=1,self%simulatedMolecule%natoms
                kinetic_energy = kinetic_energy + sumSquares(self%velocities(i,:)) &
                            * self%simulatedMolecule%atoms(i)%mass
            enddo
            !$OMP END PARALLEL DO
            kinetic_energy = kinetic_energy * 0.5
        end function kinetic_energy
 
        subroutine calculate_forces(self)
            class(simulator), target        :: self
            type(molecule),   pointer       :: imol => null()
            type(potential),  pointer       :: ipot => null()
            double precision, dimension(3)  :: f12
            integer                         :: i,j,k,x,ncells,id

            imol => self%simulatedMolecule 
            ipot => self%potentialFunction

            !FIRST ADD TWO BODY TERMS
            ncells = (imol%natoms**2 - imol%natoms)*0.5
            self%forces = 0.0d0
            self%forceMatrix = 0.0d0
            !uses https://stackoverflow.com/a/8135316 to vectorize F elements
            ! Force matrix F (Fij is force on atom i due to atom j)
            ! The total force on atom i is the sum over j
            !
            !$OMP  PARALLEL DO SCHEDULE(STATIC) & 
            !$OMP& DEFAULT(SHARED) PRIVATE(x,i,j,k,f12) &
            !$OMP& NUM_THREADS(self%nthreads)
            do x=1,ncells
                i = imol%natoms - nint(sqrt(2.0*(1+ncells-x)))
                j = mod(x+i*(i+1)/2-1,imol%natoms)+1
                do k=1,ipot%getMaxTwoBodyTerms()
                    if(ipot%twoBodyTerms(k)%inUse)then
                        f12 = -ipot%twoBodyTerms(k)%gradient(imol%atoms(i),imol%atoms(j))
                        self%forceMatrix(i,j,:) = self%forceMatrix(i,j,:) + f12
                        self%forceMatrix(j,i,:) = self%forceMatrix(j,i,:) - f12
                    else
                        exit
                    endif
                enddo
            enddo
            !$OMP END PARALLEL DO 
            
            !NOW COMBINE TWO BODY TERMS
            !$OMP  PARALLEL DO SCHEDULE(STATIC) &
            !$OMP& DEFAULT(SHARED) PRIVATE(i,j) &
            !$OMP& NUM_THREADS(self%nthreads)
            do i=1,imol%natoms
                do j=1,imol%natoms
                    if (j==i) continue
                    self%forces(i,:) = self%forces(i,:) + self%forceMatrix(i,j,:)
                enddo
            enddo
            !$OMP END PARALLEL DO 

            !NOW DO BOND/ANGLE/DIHEDRAL TERMS
            if(ipot%bond%inUse)then
                !$OMP  PARALLEL DO SCHEDULE(STATIC) & 
                !$OMP& DEFAULT(SHARED) PRIVATE(i,j,id,f12) &
                !$OMP& NUM_THREADS(self%nthreads)
                do i=1,imol%nbonds
                    do j=1,2
                        f12 = -ipot%bond%gradient(imol%bonds(i),j)
                        id  = imol%bonds(i)%ids(j)
                        !$OMP CRITICAL
                        self%forces(id,:) = self%forces(id,:) + f12
                        !$OMP END CRITICAL
                    enddo
                enddo
                !$OMP END PARALLEL DO 
            endif

            if(ipot%angle%inUse)then
                !$OMP  PARALLEL DO SCHEDULE(STATIC) & 
                !$OMP& DEFAULT(SHARED) PRIVATE(i,j,id,f12) &
                !$OMP& NUM_THREADS(self%nthreads)
                do i=1,imol%nangles
                    do j=1,3
                        f12 = -ipot%angle%gradient(imol%angles(i),j)
                        id  = imol%angles(i)%ids(j)
                        !$OMP CRITICAL
                        self%forces(id,:) = self%forces(id,:) + f12
                        !$OMP END CRITICAL
                    enddo
                enddo
                !$OMP END PARALLEL DO 
            endif
           
            if(ipot%dihedral%inUse)then
                !$OMP  PARALLEL DO SCHEDULE(STATIC) & 
                !$OMP& DEFAULT(SHARED) PRIVATE(i,j,id,f12) &
                !$OMP& NUM_THREADS(self%nthreads)
                do i=1,imol%ndihedrals
                    do j=1,4
                        f12 = -ipot%dihedral%gradient(imol%dihedrals(i),j)
                        id  = imol%dihedrals(i)%ids(j)
                        !$OMP CRITICAL
                        self%forces(id,:) = self%forces(id,:) + f12
                        !$OMP END CRITICAL
                    enddo
                enddo
                !$OMP END PARALLEL DO 
            endif
        end subroutine calculate_forces

        subroutine calculate_accelerations(self)
            class(simulator)                :: self
            integer                         :: i

            !$OMP  PARALLEL DO SCHEDULE(STATIC) &
            !$OMP& DEFAULT(SHARED) PRIVATE(i) &
            !$OMP& NUM_THREADS(self%nthreads)
            do i=1,self%simulatedMolecule%natoms
                self%accelerations(i,:) = self%forces(i,:)/self%simulatedMolecule%atoms(i)%mass
            enddo
            !$OMP END PARALLEL DO
        end subroutine calculate_accelerations

        subroutine update_velocities(self)
            class(simulator) :: self
            integer          :: i

            !$OMP  PARALLEL DO SCHEDULE(STATIC) &
            !$OMP& DEFAULT(SHARED) PRIVATE(i) &
            !$OMP& NUM_THREADS(self%nthreads)
            do i=1,self%simulatedMolecule%natoms
                self%velocities(i,:) = self%velocities(i,:) + self%accelerations(i,:)*self%dt
            enddo 
            !$OMP END PARALLEL DO
        end subroutine update_velocities
 
       subroutine update_positions(self)
            class(simulator)                :: self
            double precision, dimension(3)  :: ivec
            integer                         :: i

            !$OMP  PARALLEL DO SCHEDULE(STATIC) &
            !$OMP& DEFAULT(SHARED) PRIVATE(i,ivec) &
            !$OMP& NUM_THREADS(self%nthreads)
            do i=1,self%simulatedMolecule%natoms
                ivec = (self%velocities(i,:) + self%accelerations(i,:)*self%dt*0.5)*self%dt
                call self%simulatedMolecule%atoms(i)%translate(ivec)
            enddo
            !$OMP END PARALLEL DO
        end subroutine update_positions

end module simulator_class
