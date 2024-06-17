MODULE binding
CONTAINS
    SUBROUTINE config_filename(filename) BIND(C, name="binding_config_filename")
        CHARACTER(LEN=:), POINTER, INTENT(OUT) :: filename
        INTEGER :: arg_len

        CALL get_command_argument(1, LENGTH=arg_len)
        allocate(CHARACTER(LEN=arg_len) :: filename)
        CALL get_command_argument(1, filename)
    END

    SUBROUTINE free_config_filename(filename) BIND(C, name="binding_free_config_filename")
        CHARACTER(LEN=:), POINTER, INTENT(INOUT) :: filename
        deallocate(filename)
    END

    SUBROUTINE allocate_double_array(n, descriptor) BIND(C, NAME="binding_allocate_double_array")
        INTEGER(8), VALUE, INTENT(IN) :: n
        REAL(8), POINTER, INTENT(INOUT) :: descriptor(:)
        ALLOCATE(descriptor(n))
    END
END
