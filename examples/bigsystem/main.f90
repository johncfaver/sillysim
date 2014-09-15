program atomtest
    use molecule_class
    use potential_class
    use simulator_class

    implicit none

    type(molecule)              :: mymol
    type(potential)             :: mypotential
    type(simulator)             :: mysim
    integer                     :: nthreads = 6

        
    mymol = createRandomMolecule("example", 150, 30.0)

    call mymol%printInfo()

    !Create simulator with default potential function
    mysim = createSimulator(mymol, outputfile="out.pdb", &
                        nsteps=150000, outputFrequency=100, dt=0.001, nthreads=nthreads)                  

    !Minimize a bit until dE <= 0.001
    call mysim%minimize(threshold=1.0d-3)

    !Simulate
    call mysim%simulate()

end program atomtest
