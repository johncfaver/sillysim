program atomtest
    use atom_class
    use molecule_class
    use potential_class
    use simulator_class
    use geometry

    implicit none

    type(atom), dimension(401)  :: myatoms
    type(molecule)              :: mymol
    type(potential)             :: mypotential
    type(simulator)             :: mysim
    integer                     :: nthreads = 6, i, j, nrow=20

    double precision :: space = 3.0d0, midpoint
    midpoint = space*nrow/2.0d0+1.5d0

    do i=1,size(myatoms)-1
        myatoms(i) = createAtom("CL",x=dble(mod(i,nrow))*space,y=dble((i-1)/nrow)*space, &
            z=0.0d0,charge=-0.5,mass=200.0,ljs=1.0,lje=0.5)
    enddo
    myatoms(size(myatoms)) = createAtom("O",x=midpoint,y=-5.0d0,z=0.0d0,charge=3.0,mass=10.0,ljs=0.2,lje=0.2)

    mymol = createMolecule("grid", myatoms)

    do i=1,100
        do j=i+1,100
            if(norm(r12(mymol%atoms(i),mymol%atoms(j))) - 1.0d0 < 1.0d-3 )then 
                call mymol%createBond(i,j,r0=1.0d0,k=1.0d7)
            endif
        enddo
    enddo

    call mymol%printInfo()
    call mypotential%addTerm(1)
    call mypotential%addTerm(4)


    call mypotential%addTerm(5)


!    Create simulator with default potential function
    mysim = createSimulator(mymol, mypotential, outputfile="out.pdb", &
                        outputFrequency=1000, dt=0.001, nthreads=nthreads)                  

    !Minimize a bit until dE <= 0.001
    !call mysim%minimize(threshold=1.0d-3)
    !Simulate
    call mysim%simulate(nsteps=200000)

end program atomtest
