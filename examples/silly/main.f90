program atomtest
    use molecule_class
    use potential_class
    use simulator_class

    implicit none

    type(molecule)              :: mymol
    type(potential)             :: mypotential
    type(simulator)             :: mysim
    integer                     :: nthreads = 6

    mymol = createRandomMolecule("example", 5, 10.0)

    !I don't really know what I'm making    
    call mymol%createBond(1,2,r0=2.0d0)
    call mymol%createBond(1,3,r0=2.5d0,k=100.0d0)
    call mymol%createBond(2,4,k=120.0d0,r0=5.0d0)
    call mymol%createBond(4,5)
    call mymol%createBond(5,1)
    call mymol%createAngle(1,2,3,r0=109.5d0)
    call mymol%createAngle(2,3,4,r0=180.d0)
    call mymol%createAngle(4,5,1,r0=90.d0)
    call mymol%createAngle(3,4,5,r0=90.d0)
    call mymol%createAngle(2,5,1,r0=120.d0)

    call mymol%printInfo()

    !Create simulator with default potential function
    mysim = createSimulator(mymol, outputfile="out.pdb", &
                        nsteps=100000, outputFrequency=100, dt=0.001, nthreads=nthreads)                  

    !Minimize a bit until dE <= 0.001
    call mysim%minimize(threshold=1.0d-3)

    print *, 'Minimized to E=', mysim%potential_energy()

    !Simulate
    call mysim%simulate()

end program atomtest
