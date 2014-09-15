program atomtest
    use molecule_class
    use potential_class
    use simulator_class

    implicit none

    type(molecule)              :: mymol
    type(potential)             :: mypotential
    type(simulator)             :: mysim
    integer                     :: nthreads = 4

    !Make a random molecule with 7 atoms, placed in a cube with sides of 10 Å.
    mymol = createRandomMolecule("rando", 7, 10.0)

    !Make mymol octahedral by assigning bonds/angles
    !Atom 1 at the center, 2-5 in a plane, 6-7 in the perpendicular plane
    !                           6  3
    !                           | /
    !                           |/
    !                      2----1----4
    !                          /|
    !                         / |
    !                        5  7
    call mymol%createBond(1,2)
    call mymol%createBond(1,3)
    call mymol%createBond(1,4)
    call mymol%createBond(1,5)
    call mymol%createBond(1,6)
    call mymol%createBond(1,7)
    call mymol%createAngle(2,1,3,r0=90.d0)
    call mymol%createAngle(2,1,4,r0=180.d0)
    call mymol%createAngle(2,1,5,r0=90.d0)
    call mymol%createAngle(2,1,6,r0=90.d0)
    call mymol%createAngle(2,1,7,r0=90.d0)
    call mymol%createAngle(3,1,4,r0=90.d0)
    call mymol%createAngle(3,1,5,r0=180.d0)
    call mymol%createAngle(3,1,6,r0=90.d0)
    call mymol%createAngle(3,1,7,r0=90.d0)
    call mymol%createAngle(4,1,5,r0=90.d0)
    call mymol%createAngle(4,1,6,r0=90.d0)
    call mymol%createAngle(4,1,7,r0=90.d0)
    call mymol%createAngle(5,1,6,r0=90.d0)
    call mymol%createAngle(5,1,7,r0=90.d0)
    call mymol%createAngle(6,1,7,r0=180.d0)
    
    !Print description of molecule
    call mymol%printInfo()

    !Set up a potential function object
    call mypotential%addTerm(1) !Add bond terms
    call mypotential%addTerm(2) !Add angle terms

    !Print description of potential function
    call mypotential%printInfo()

    !Create simulator object
    mysim = createSimulator(mymol, mypotential, outputfile="out.pdb", &
                        nsteps=50000, outputFrequency=100, dt=0.001, nthreads=nthreads)                  

    !Print current potential energy of mymol evaluated with mypotential 
    print *, 'E=', mysim%potential_energy()

    !Minimize mymol using mypotential
    call mysim%minimize()

end program atomtest
