program atomtest
    use molecule_class
    use potential_class
    use simulator_class

    implicit none

    type(molecule)              :: mymol
    type(potential)             :: mypotential
    type(simulator)             :: mysim
    integer                     :: nthreads = 6

    !Make a molecule with atoms placed randomly in a box
    mymol = createRandomMolecule("Test", natoms=10)

    !Define some covalent bonds
    call mymol%createBond(1,2)
    call mymol%createBond(1,3, r0=2.5d0, k=50.0d0)
    call mymol%createBond(5,3, r0=1.5d0)

    !Print information for mymol
    call mymol%printInfo()
   
    !Create potential energy function
    call mypotential%addTerm(1) !Add covalent bond terms
    call mypotential%addTerm(4) !Add Lennard Jones terms
   
    !Print information for mypotential 
    call mypotential%printInfo()

    !Create simulator object
    mysim = createSimulator(mymol, mypotential, outputfile="out.pdb", &
                        outputFrequency=1000, dt=0.001, nthreads=nthreads)                  

    !Relax system for 1000 steps
    call mysim%minimize(nsteps=1000)
   
    !Start dynamics simulation
    call mysim%simulate(nsteps=200000)

end program atomtest
