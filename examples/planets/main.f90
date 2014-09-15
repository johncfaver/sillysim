program atomtest
    use molecule_class
    use potential_class
    use simulator_class

    implicit none

    type(molecule)              :: mymol
    type(potential)             :: mypotential
    type(simulator)             :: mysim
    integer                     :: nthreads = 4

    !Make a random system with 25 atoms, placed in a cube with sides of 30 
    mymol = createRandomMolecule("galaxy", 20, 30.0)

    !Set up a potential function object
    call mypotential%addTerm(6) !Add gravity

    !Print description of potential function
    call mypotential%printInfo()

    !Create simulator 
    mysim = createSimulator(mymol, mypotential, outputfile="out.pdb", &
                        nsteps=100000, outputFrequency=100, dt=0.001, nthreads=nthreads)                  

    call mysim%simulate()

end program atomtest
