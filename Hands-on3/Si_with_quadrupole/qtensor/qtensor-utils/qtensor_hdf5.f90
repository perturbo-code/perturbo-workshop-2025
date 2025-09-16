! The easiest way is to compile with hdf5 wrapper
!
program main
    !
    use hdf5  
    !
    implicit none
    !
    !character(len=256) :: filename = "tlw_3.abo"
    character(len=256) :: filename
    !
    integer :: iunit, count, status, i, j, num_row, num_col, na, nat
    integer :: pert, pert_tmp
    !
    real(kind = 8), allocatable, dimension(:,:)     :: qtensor_tmp
    real(kind = 8), allocatable, dimension(:,:,:,:) :: qtensor
    !
    character(len=128) :: line, dir
    
    ! HDF5 variables
    integer(HID_T) :: file_id, dset_id, dspace_id
    integer(HSIZE_T), dimension(4) :: dims
    integer :: error
    !
    ! initial num_row, num_col
    num_row = 0
    num_col = 9
    !
    iunit = 37
    !
    write(*,"(A)", advance="no") "Please give the anaddb output file name: "
    read(*,"(a)") filename
    open(file=filename, unit = iunit, status = "old", action = "read")
    !
        do
            read(iunit, "(A)", iostat = status) line
            !
            if (status /= 0) then
                write(*,"(3x,A)") "Please check your inputfile"
                stop
            end if
            !
            if (index(line,'atom   dir       Qxx         Qyy') > 0) then
                exit
            end if
            !
        end do
        !
        ! we have found the quadrupole tensor
        !
        do while (.true.)
            !
            read(iunit, "(A)", iostat = status) line
            !
            if (status == -1 .or. index(trim(line), "===") > 0) then
                exit
            end if
            !
            num_row = num_row + 1 ! number of row for qtensor
        end do
        !
        write(*,"(2x,A,2x,I5)") "The number of atom in your system is: ", (num_row - 1) / 3
        !
        !num_row = num_row + 1
        !
        do i = 1, num_row + 1
            backspace(iunit)
        end do
        !
        num_row = num_row - 1
        nat = num_row / 3
        !
        allocate( qtensor(3,3,3,nat), qtensor_tmp(num_row,6))
        !
        write(*,"(2x,A)") "The Quadrupole tensor has read:"
        !
        do i = 1, num_row
            read(iunit,*) na, dir, (qtensor_tmp(i,j), j = 1, 6)
            write(*,"(2x,I5,4x,A3,x,6F15.7,2x)") na, dir, (qtensor_tmp(i,j), j = 1, 6)
        end do
        !
    close(iunit)
    !
    qtensor = 0.d0
    !
    pert_tmp = 0
    !
    do na = 1, nat
        do pert = 1, 3
            !
            pert_tmp = pert_tmp + 1
            !
            do i = 1, 3
                qtensor(i,i,pert,na) = qtensor_tmp(pert_tmp,i)
            end do
            !
            qtensor(2,3,pert,na) = qtensor_tmp(pert_tmp,4)
            qtensor(1,3,pert,na) = qtensor_tmp(pert_tmp,5)
            qtensor(1,2,pert,na) = qtensor_tmp(pert_tmp,6)
            qtensor(2,1,pert,na) = qtensor(1,2,pert,na)
            qtensor(3,1,pert,na) = qtensor(1,3,pert,na)
            qtensor(3,2,pert,na) = qtensor(2,3,pert,na)
        end do
        !pert_tmp = pert_tmp + 3
    end do
    !
    ! Write qtensor as h5
    call h5open_f(error)
    
    call h5fcreate_f("qtensor.h5", H5F_ACC_TRUNC_F, file_id, error)
    
    ! Set dimensions for the dataset
    dims = [3, 3, 3, nat]
    
    ! Create dataspace
    call h5screate_simple_f(4, dims, dspace_id, error)
    
    ! Create dataset
    call h5dcreate_f(file_id, "/qtensor", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    
    ! Write data to dataset
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qtensor, dims, error)
    
    ! Close everything
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
    !
    write(*,"(2x,A)") "The Quadrupole tensor has been written to qtensor.h5! Please rename it for your own system!!!"
    !
    !
end program main
