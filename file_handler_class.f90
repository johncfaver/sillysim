! Handles file input/output
! readMolecule will read a molecule from a file and return a molecule object
! writeMolecule will write the molecule to a file

module file_handler_class
    use molecule_class
    use atom_class
    implicit none
    private
    
    public :: readMolecule, writeMolecule

    character(len=*), parameter :: PDBFormat  = "(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,10X,2A2)"

    contains
        function readMolecule(filename)
            ! Return molecule object read from file.
            ! Figure out file type, and use appropriate reading function.
            type(molecule)   :: readMolecule
            character(len=*) :: filename
            if(index(filename,".pdb") > 0) then
                readMolecule = readPDB(filename)
            else
                print *, 'Unnkown file format. Could not read ', filename
            endif
        end function readMolecule

        subroutine writeMolecule(filename,molobj,append)
            ! Write molecule object to file.
            ! Figure out desired file type from name, then use appropriate writer function.
            ! Append to file if append is given and .true.
            character(len=*)  :: filename
            type(molecule)    :: molobj
            logical, optional :: append
            logical           :: shouldAppend = .FALSE.

            if(present(append)) shouldAppend = append
            if(index(filename,".pdb") > 0) then
                call writePDB(filename,molobj,shouldAppend)
            else
                print *, 'Unknown file format; Could not write ', filename
            endif
        end subroutine writeMolecule

        function readPDB(filename)
            character(len=*)                    :: filename
            type(molecule)                      :: readPDB
            type(atom),dimension(:),allocatable :: atoms
            character(len=50)                   :: title="Default Title"
            character(len=80)                   :: templine
            character(len=6)                    :: record
            character(len=4)                    :: atomName
            character(len=3)                    :: resName
            character(len=2)                    :: element
            character                           :: chain
            integer                             :: atomNumber, resNumber
            integer                             :: natoms=0, ios=1, iatom=0
            double precision                    :: x, y, z, occ, beta

            open(21, file=filename, status="old", iostat=ios, action="READ")
            do
                read(21,'(A80)',iostat=ios) templine
                if(ios /= 0) exit
                if(index(templine,"ATOM  ") == 1) natoms = natoms + 1
                if(index(templine,"HETATM") == 1) natoms = natoms + 1
                if(index(templine,"TITLE ") == 1) title  = templine(7:50)
            enddo
            allocate(atoms(natoms))
            rewind(21)
            do
                read(21,fmt='(A80)',iostat=ios) templine
                if(ios /= 0) exit
                if(index(templine,"ATOM  ") == 1 .or. &
                   index(templine,"HETATM") == 1) then
                    read(templine,PDBFormat) record, atomNumber, atomName, resName, &
                                    chain, resNumber, x, y, z , occ, beta, element
                    iatom = iatom + 1
                    atoms(iatom) = createAtom(element,x,y,z)
                endif
            enddo
            close(21,iostat=ios)
            readPDB = createMolecule(title, atoms)
            deallocate(atoms)    
        end function readPDB

        subroutine writePDB(filename,molobj,append)
            character(len=*)                             :: filename
            character(len=80), dimension(:), allocatable :: header, footer
            logical                                      :: append
            type(molecule)                               :: molobj
            integer                                      :: i
            
            allocate(header(2),footer(1))
            write(header(1),'(A6,A74)') "TITLE ", molobj%molname
            write(header(2),'(A5)') "MODEL"
            write(footer,'(A6)') "ENDMDL" 
            if(append) then
                open(20,file=filename,action="WRITE",position="APPEND")
            else
                open(20,file=filename,action="WRITE",status="REPLACE")
            endif
            do i= 1,size(header)
                write(20,'(A80)') header(i)
            enddo
            do i = 1,molobj%natoms
                write(20,PDBFormat) "ATOM  ",i,molobj%atoms(i)%element,"UNK","A", &
                   1,molobj%atoms(i)%coords,1.0,0.0,molobj%atoms(i)%element
            enddo
            write(20,'(A6)') footer
            close(20) 
            deallocate(header,footer)
        end subroutine writePDB

end module file_handler_class
