program atomtest
    use atom_class
    use molecule_class
    use potential_class
    use simulator_class
    use file_handler_class
    use geometry

    implicit none

    type(molecule)              :: mymol
    type(potential)             :: mypot
    type(simulator)             :: mysim
    type(atom), dimension(:), allocatable :: myatoms
    integer                     :: nthreads = 6, n,nn,i,j,nwat, wcount
    double precision, dimension(3) :: f1,f2,f3,f4,a
    logical     :: clash

    !load a molecule
    mymol = readMolecule("a.pdb")
    !we will add nwat waters 
    nwat = 300

    n = mymol%natoms
    nn = n + 3*nwat
    allocate(myatoms(nn))
    do i=1,n
        myatoms(i) = mymol%atoms(i)
    enddo
    
    f1 = mymol%centerOfMolecule()
    f2 = [5.0d0,0.0d0,0.0d0]
    f3 = [0.0d0,5.0d0,0.0d0]
    f4 = [0.0d0,0.0d0,5.0d0]

    f1 = f1 - 5.0*f2 - 5.0*f3 - 5.0*f4
    ! 10 x 10 x 10 water box
    ! don't overlap with molecule
    wcount = 0
    do i=1,3*nwat, 3
        a = f1 + mod(i,10)*f2 + mod((i-1)/5,10)*f3 + (i-1)/100*f4 
        clash = .false.
        do j=1,n
            if(norm(mymol%atoms(j)%coords-a) < 2.0d0)then
                clash = .true.
                exit
            endif
        enddo
        if(clash) cycle
        myatoms(n+(wcount)*3+1) = createAtom("O",a(1),a(2),a(3))
        myatoms(n+(wcount)*3+2) = createAtom("H",a(1),a(2),a(3)+1.0d0)
        myatoms(n+(wcount)*3+3) = createAtom("H",a(1),a(2)+0.4d0,a(3)-0.7d0)
        wcount = wcount + 1
    enddo

    mymol = createMolecule("boxy",myatoms(1:n+wcount*3))
    call mymol%guessTopology(guessBonds=.true.,guessAngles=.true.)

    call mymol%printInfo()
    call mypot%addTerm(1)
    call mypot%addTerm(2)
    call mypot%addTerm(4)
    call mypot%addTerm(5)

    call mypot%printInfo()
    mysim = createSimulator(mymol,                   &
                            potentialFunction=mypot, &
                            nthreads=nthreads,       &
                            outputFrequency=100,     &
                            nbcutoff=15.0d0,        &
                            nbUpdateFrequency=50,   &
                            dt=0.01)

    !call mysim%minimize(nsteps=2)
    call mysim%simulate(nsteps=50000)

end program atomtest
