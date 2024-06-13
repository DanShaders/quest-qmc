MODULE dqmc_util

USE detour_dqmc_util_defs

IMPLICIT NONE

INTERFACE
    FUNCTION dqmc_matdiff(n, a, b) BIND(C, NAME="fortran_dqmc_util_dqmc_matdiff")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(IN) :: a(n, n)
        REAL(8), INTENT(IN) :: b(n, n)
        REAL(8) :: dqmc_matdiff
    END

    FUNCTION dqmc_matnorm(n, a) BIND(C, NAME="fortran_dqmc_util_dqmc_matnorm")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(IN) :: a(n, n)
        REAL(8) :: dqmc_matnorm
    END

    SUBROUTINE dqmc_eye(n, a) BIND(C, NAME="fortran_dqmc_util_dqmc_eye")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(INOUT) :: a(n, n)
    END

    SUBROUTINE dqmc_trans(n, at, a) BIND(C, NAME="fortran_dqmc_util_dqmc_trans")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(INOUT) :: at(n, n)
        REAL(8), INTENT(IN) :: a(n, n)
    END

    SUBROUTINE dqmc_scalecol(n, a, d) BIND(C, NAME="fortran_dqmc_util_dqmc_scalecol")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(INOUT) :: a(n, n)
        REAL(8), INTENT(IN) :: d(n)
    END

    SUBROUTINE dqmc_scalerow(n, a, d) BIND(C, NAME="fortran_dqmc_util_dqmc_scalerow")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(INOUT) :: a(n, n)
        REAL(8), INTENT(IN) :: d(n)
    END

    SUBROUTINE dqmc_scalecolinv(n, a, d) BIND(C, NAME="fortran_dqmc_util_dqmc_scalecolinv")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(INOUT) :: a(n, n)
        REAL(8), INTENT(IN) :: d(n)
    END

    SUBROUTINE dqmc_scalerowinv(n, a, d) BIND(C, NAME="fortran_dqmc_util_dqmc_scalerowinv")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(INOUT) :: a(n, n)
        REAL(8), INTENT(IN) :: d(n)
    END

    SUBROUTINE dqmc_signjackknife_real(n, avg, err, x, y, sgn, sum_sgn) BIND(C, NAME="fortran_dqmc_util_dqmc_signjackknife_real")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(OUT) :: avg
        REAL(8), INTENT(OUT) :: err
        REAL(8), INTENT(IN) :: x(:)
        REAL(8), INTENT(INOUT) :: y(:)
        REAL(8), INTENT(IN) :: sgn(:)
        REAL(8), INTENT(IN) :: sum_sgn
    END

    SUBROUTINE dqmc_signjackknife_complex(n, avg, err, x, y, sgn, sum_sgn) BIND(C, NAME="fortran_dqmc_util_dqmc_signjackknife_complex")
        INTEGER(4), INTENT(IN) :: n
        COMPLEX(8), INTENT(OUT) :: avg
        COMPLEX(8), INTENT(OUT) :: err
        COMPLEX(8), INTENT(IN) :: x(:)
        COMPLEX(8), INTENT(INOUT) :: y(:)
        COMPLEX(8), INTENT(IN) :: sgn(:)
        COMPLEX(8), INTENT(IN) :: sum_sgn
    END

    SUBROUTINE dqmc_jackknife_real(n, avg, err, x, y, sgn, sum_sgn) BIND(C, NAME="fortran_dqmc_util_dqmc_jackknife_real")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(OUT) :: avg
        REAL(8), INTENT(OUT) :: err
        REAL(8), INTENT(IN) :: x(n)
        REAL(8), INTENT(OUT) :: y(n)
        REAL(8), INTENT(OUT) :: sgn(n)
        REAL(8), INTENT(OUT) :: sum_sgn
    END

    SUBROUTINE dqmc_jackknife_complex(n, avg, err, x, y, sgn, sum_sgn) BIND(C, NAME="fortran_dqmc_util_dqmc_jackknife_complex")
        INTEGER(4), INTENT(IN) :: n
        COMPLEX(8), INTENT(OUT) :: avg
        COMPLEX(8), INTENT(OUT) :: err
        COMPLEX(8), INTENT(IN) :: x(n)
        COMPLEX(8), INTENT(OUT) :: y(n)
        COMPLEX(8), INTENT(OUT) :: sgn(n)
        COMPLEX(8), INTENT(OUT) :: sum_sgn
    END

    SUBROUTINE dqmc_geterr(n, err, avg, list) BIND(C, NAME="fortran_dqmc_util_dqmc_geterr")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(OUT) :: err
        REAL(8), INTENT(IN) :: avg
        REAL(8), INTENT(INOUT) :: list(n)
    END

    SUBROUTINE dqmc_geterr1(n, data, avg, err) BIND(C, NAME="fortran_dqmc_util_dqmc_geterr1")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(IN) :: data(n)
        REAL(8), INTENT(OUT) :: avg
        REAL(8), INTENT(OUT) :: err
    END

    SUBROUTINE dqmc_geterr2(n, sm, ssq, avg, err) BIND(C, NAME="fortran_dqmc_util_dqmc_geterr2")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(IN) :: sm
        REAL(8), INTENT(IN) :: ssq
        REAL(8), INTENT(OUT) :: avg
        REAL(8), INTENT(OUT) :: err
    END

    SUBROUTINE dqmc_error(message, no) BIND(C, NAME="fortran_dqmc_util_dqmc_error")
        CHARACTER(LEN=*), INTENT(IN) :: message
        INTEGER(4), INTENT(IN) :: no
    END

    SUBROUTINE dqmc_warning(message, no) BIND(C, NAME="fortran_dqmc_util_dqmc_warning")
        CHARACTER(LEN=*), INTENT(IN) :: message
        INTEGER(4), INTENT(IN) :: no
    END

    SUBROUTINE ran0(n, var, seed) BIND(C, NAME="fortran_dqmc_util_ran0")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(OUT) :: var(n)
        INTEGER(4), INTENT(INOUT) :: seed(4)
    END

    SUBROUTINE ran1(n, var, seed) BIND(C, NAME="fortran_dqmc_util_ran1")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(OUT) :: var(n)
        INTEGER(4), INTENT(INOUT) :: seed(4)
    END

    FUNCTION intran(l, seed) BIND(C, NAME="fortran_dqmc_util_intran")
        INTEGER(4), INTENT(IN) :: l
        INTEGER(4), INTENT(INOUT) :: seed(4)
        INTEGER(4) :: intran
    END

    SUBROUTINE rann(n, var, seed) BIND(C, NAME="fortran_dqmc_util_rann")
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(OUT) :: var(n)
        INTEGER(4), INTENT(INOUT) :: seed(4)
    END

    SUBROUTINE dumpa(m, n, a, opt) BIND(C, NAME="fortran_dqmc_util_dumpa")
        INTEGER(4), INTENT(IN) :: m
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(IN) :: a(m, n)
        INTEGER(4), INTENT(IN) :: opt
    END

    SUBROUTINE dqmc_print_realarray(n, m, title, label, avg, err, opt) BIND(C, NAME="fortran_dqmc_util_dqmc_print_realarray")
        INTEGER(4), INTENT(IN) :: n
        INTEGER(4), INTENT(IN) :: m
        CHARACTER(LEN=*), INTENT(IN) :: title
        CHARACTER(LEN=*), INTENT(IN) :: label(:)
        REAL(8), INTENT(IN) :: avg(:, :)
        REAL(8), INTENT(IN) :: err(:, :)
        INTEGER(4), INTENT(IN) :: opt
    END

    SUBROUTINE dqmc_print_complexarray(n, m, title, label, avg, err, opt) BIND(C, NAME="fortran_dqmc_util_dqmc_print_complexarray")
        INTEGER(4), INTENT(IN) :: n
        INTEGER(4), INTENT(IN) :: m
        CHARACTER(LEN=*), INTENT(IN) :: title
        CHARACTER(LEN=*), INTENT(IN) :: label(:)
        COMPLEX(8), INTENT(IN) :: avg(:, :)
        COMPLEX(8), INTENT(IN) :: err(:, :)
        INTEGER(4), INTENT(IN) :: opt
    END

    SUBROUTINE dqmc_print_eigenmode(n, m, title, value, opt) BIND(C, NAME="fortran_dqmc_util_dqmc_print_eigenmode")
        INTEGER(4), INTENT(IN) :: n
        INTEGER(4), INTENT(IN) :: m
        CHARACTER(LEN=*), INTENT(IN) :: title
        COMPLEX(8), INTENT(IN) :: value(:, :, :)
        INTEGER(4), INTENT(IN) :: opt
    END

    SUBROUTINE dqmc_getftk(value, n, nclass, class, na, nk, ft_wgt, phase, valuek) BIND(C, NAME="fortran_dqmc_util_dqmc_getftk")
        REAL(8), INTENT(IN) :: value(nclass)
        INTEGER(4), INTENT(IN) :: n
        INTEGER(4), INTENT(IN) :: nclass
        INTEGER(4), INTENT(IN) :: class(n, n)
        INTEGER(4), INTENT(IN) :: na
        INTEGER(4), INTENT(IN) :: nk
        COMPLEX(8), INTENT(IN) :: ft_wgt(n / na, nk)
        INTEGER(4), INTENT(IN) :: phase(n, n)
        COMPLEX(8), INTENT(OUT) :: valuek(nk * na * (na + 1) / 2)
    END

    SUBROUTINE dqmc_io_open(fname, inp_unit, out_unit) BIND(C, NAME="fortran_dqmc_util_dqmc_io_open")
        CHARACTER(LEN=*), INTENT(OUT) :: fname
        INTEGER(4), INTENT(OUT) :: inp_unit
        INTEGER(4), INTENT(OUT) :: out_unit
    END

    SUBROUTINE dqmc_open_file(fname, fstatus, file_unit) BIND(C, NAME="fortran_dqmc_util_dqmc_open_file")
        CHARACTER(LEN=*), INTENT(IN) :: fname
        CHARACTER(LEN=*), INTENT(IN) :: fstatus
        INTEGER(4), INTENT(OUT) :: file_unit
    END

    SUBROUTINE dqmc_count_records(n, file_unit) BIND(C, NAME="fortran_dqmc_util_dqmc_count_records")
        INTEGER(4), INTENT(OUT) :: n
        INTEGER(4), INTENT(IN) :: file_unit
    END

    FUNCTION move_to_record(string, iunit) BIND(C, NAME="fortran_dqmc_util_move_to_record")
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER(4), INTENT(IN) :: iunit
        LOGICAL(4) :: move_to_record
    END

    FUNCTION get_det(a) BIND(C, NAME="fortran_dqmc_util_get_det")
        REAL(8), INTENT(IN) :: a(3, 3)
        REAL(8) :: get_det
    END

    SUBROUTINE get_inverse(a, inv) BIND(C, NAME="fortran_dqmc_util_get_inverse")
        REAL(8), INTENT(IN) :: a(3, 3)
        REAL(8), INTENT(OUT) :: inv(3, 3)
    END
END INTERFACE

INTERFACE DQMC_JackKnife
 PROCEDURE DQMC_JackKnife_Real, DQMC_JackKnife_Complex
END INTERFACE DQMC_JackKnife

INTERFACE DQMC_SignJackKnife
 PROCEDURE DQMC_SignJackKnife_Real, DQMC_SignJackKnife_Complex
END INTERFACE DQMC_SignJackKnife

INTERFACE DQMC_Print_Array
 PROCEDURE DQMC_Print_RealArray, DQMC_Print_ComplexArray
END INTERFACE DQMC_Print_Array

CONTAINS
FUNCTION detour_dqmc_matdiff(n, a, b) BIND(C, NAME="detour_dqmc_util_dqmc_matdiff")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(IN) :: a(n, n)
    REAL(8), INTENT(IN) :: b(n, n)
    REAL(8) :: detour_dqmc_matdiff

    detour_dqmc_matdiff = dqmc_matdiff(n, a, b)
END

FUNCTION detour_dqmc_matnorm(n, a) BIND(C, NAME="detour_dqmc_util_dqmc_matnorm")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(IN) :: a(n, n)
    REAL(8) :: detour_dqmc_matnorm

    detour_dqmc_matnorm = dqmc_matnorm(n, a)
END

SUBROUTINE detour_dqmc_eye(n, a) BIND(C, NAME="detour_dqmc_util_dqmc_eye")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(INOUT) :: a(n, n)

    CALL dqmc_eye(n, a)
END

SUBROUTINE detour_dqmc_trans(n, at, a) BIND(C, NAME="detour_dqmc_util_dqmc_trans")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(INOUT) :: at(n, n)
    REAL(8), INTENT(IN) :: a(n, n)

    CALL dqmc_trans(n, at, a)
END

SUBROUTINE detour_dqmc_scalecol(n, a, d) BIND(C, NAME="detour_dqmc_util_dqmc_scalecol")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(INOUT) :: a(n, n)
    REAL(8), INTENT(IN) :: d(n)

    CALL dqmc_scalecol(n, a, d)
END

SUBROUTINE detour_dqmc_scalerow(n, a, d) BIND(C, NAME="detour_dqmc_util_dqmc_scalerow")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(INOUT) :: a(n, n)
    REAL(8), INTENT(IN) :: d(n)

    CALL dqmc_scalerow(n, a, d)
END

SUBROUTINE detour_dqmc_scalecolinv(n, a, d) BIND(C, NAME="detour_dqmc_util_dqmc_scalecolinv")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(INOUT) :: a(n, n)
    REAL(8), INTENT(IN) :: d(n)

    CALL dqmc_scalecolinv(n, a, d)
END

SUBROUTINE detour_dqmc_scalerowinv(n, a, d) BIND(C, NAME="detour_dqmc_util_dqmc_scalerowinv")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(INOUT) :: a(n, n)
    REAL(8), INTENT(IN) :: d(n)

    CALL dqmc_scalerowinv(n, a, d)
END

SUBROUTINE detour_dqmc_signjackknife_real(n, avg, err, x, y, sgn, sum_sgn) BIND(C, NAME="detour_dqmc_util_dqmc_signjackknife_real")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(OUT) :: avg
    REAL(8), INTENT(OUT) :: err
    REAL(8), INTENT(IN) :: x(:)
    REAL(8), INTENT(INOUT) :: y(:)
    REAL(8), INTENT(IN) :: sgn(:)
    REAL(8), INTENT(IN) :: sum_sgn

    CALL dqmc_signjackknife_real(n, avg, err, x, y, sgn, sum_sgn)
END

SUBROUTINE detour_dqmc_signjackknife_complex(n, avg, err, x, y, sgn, sum_sgn) BIND(C, NAME="detour_dqmc_util_dqmc_signjackknife_complex")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    COMPLEX(8), INTENT(OUT) :: avg
    COMPLEX(8), INTENT(OUT) :: err
    COMPLEX(8), INTENT(IN) :: x(:)
    COMPLEX(8), INTENT(INOUT) :: y(:)
    COMPLEX(8), INTENT(IN) :: sgn(:)
    COMPLEX(8), INTENT(IN) :: sum_sgn

    CALL dqmc_signjackknife_complex(n, avg, err, x, y, sgn, sum_sgn)
END

SUBROUTINE detour_dqmc_jackknife_real(n, avg, err, x, y, sgn, sum_sgn) BIND(C, NAME="detour_dqmc_util_dqmc_jackknife_real")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(OUT) :: avg
    REAL(8), INTENT(OUT) :: err
    REAL(8), INTENT(IN) :: x(n)
    REAL(8), INTENT(OUT) :: y(n)
    REAL(8), INTENT(OUT) :: sgn(n)
    REAL(8), INTENT(OUT) :: sum_sgn

    CALL dqmc_jackknife_real(n, avg, err, x, y, sgn, sum_sgn)
END

SUBROUTINE detour_dqmc_jackknife_complex(n, avg, err, x, y, sgn, sum_sgn) BIND(C, NAME="detour_dqmc_util_dqmc_jackknife_complex")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    COMPLEX(8), INTENT(OUT) :: avg
    COMPLEX(8), INTENT(OUT) :: err
    COMPLEX(8), INTENT(IN) :: x(n)
    COMPLEX(8), INTENT(OUT) :: y(n)
    COMPLEX(8), INTENT(OUT) :: sgn(n)
    COMPLEX(8), INTENT(OUT) :: sum_sgn

    CALL dqmc_jackknife_complex(n, avg, err, x, y, sgn, sum_sgn)
END

SUBROUTINE detour_dqmc_geterr(n, err, avg, list) BIND(C, NAME="detour_dqmc_util_dqmc_geterr")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(OUT) :: err
    REAL(8), INTENT(IN) :: avg
    REAL(8), INTENT(INOUT) :: list(n)

    CALL dqmc_geterr(n, err, avg, list)
END

SUBROUTINE detour_dqmc_geterr1(n, data, avg, err) BIND(C, NAME="detour_dqmc_util_dqmc_geterr1")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(IN) :: data(n)
    REAL(8), INTENT(OUT) :: avg
    REAL(8), INTENT(OUT) :: err

    CALL dqmc_geterr1(n, data, avg, err)
END

SUBROUTINE detour_dqmc_geterr2(n, sm, ssq, avg, err) BIND(C, NAME="detour_dqmc_util_dqmc_geterr2")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(IN) :: sm
    REAL(8), INTENT(IN) :: ssq
    REAL(8), INTENT(OUT) :: avg
    REAL(8), INTENT(OUT) :: err

    CALL dqmc_geterr2(n, sm, ssq, avg, err)
END

SUBROUTINE detour_dqmc_error(message, no) BIND(C, NAME="detour_dqmc_util_dqmc_error")
    USE detour_dqmc_util
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER(4), INTENT(IN) :: no

    CALL dqmc_error(message, no)
END

SUBROUTINE detour_dqmc_warning(message, no) BIND(C, NAME="detour_dqmc_util_dqmc_warning")
    USE detour_dqmc_util
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER(4), INTENT(IN) :: no

    CALL dqmc_warning(message, no)
END

SUBROUTINE detour_ran0(n, var, seed) BIND(C, NAME="detour_dqmc_util_ran0")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(OUT) :: var(n)
    INTEGER(4), INTENT(INOUT) :: seed(4)

    CALL ran0(n, var, seed)
END

SUBROUTINE detour_ran1(n, var, seed) BIND(C, NAME="detour_dqmc_util_ran1")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(OUT) :: var(n)
    INTEGER(4), INTENT(INOUT) :: seed(4)

    CALL ran1(n, var, seed)
END

FUNCTION detour_intran(l, seed) BIND(C, NAME="detour_dqmc_util_intran")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: l
    INTEGER(4), INTENT(INOUT) :: seed(4)
    INTEGER(4) :: detour_intran

    detour_intran = intran(l, seed)
END

SUBROUTINE detour_rann(n, var, seed) BIND(C, NAME="detour_dqmc_util_rann")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(OUT) :: var(n)
    INTEGER(4), INTENT(INOUT) :: seed(4)

    CALL rann(n, var, seed)
END

SUBROUTINE detour_dumpa(m, n, a, opt) BIND(C, NAME="detour_dqmc_util_dumpa")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: m
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(IN) :: a(m, n)
    INTEGER(4), INTENT(IN) :: opt

    CALL dumpa(m, n, a, opt)
END

SUBROUTINE detour_dqmc_print_realarray(n, m, title, label, avg, err, opt) BIND(C, NAME="detour_dqmc_util_dqmc_print_realarray")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    INTEGER(4), INTENT(IN) :: m
    CHARACTER(LEN=*), INTENT(IN) :: title
    CHARACTER(LEN=*), INTENT(IN) :: label(:)
    REAL(8), INTENT(IN) :: avg(:, :)
    REAL(8), INTENT(IN) :: err(:, :)
    INTEGER(4), INTENT(IN) :: opt

    CALL dqmc_print_realarray(n, m, title, label, avg, err, opt)
END

SUBROUTINE detour_dqmc_print_complexarray(n, m, title, label, avg, err, opt) BIND(C, NAME="detour_dqmc_util_dqmc_print_complexarray")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    INTEGER(4), INTENT(IN) :: m
    CHARACTER(LEN=*), INTENT(IN) :: title
    CHARACTER(LEN=*), INTENT(IN) :: label(:)
    COMPLEX(8), INTENT(IN) :: avg(:, :)
    COMPLEX(8), INTENT(IN) :: err(:, :)
    INTEGER(4), INTENT(IN) :: opt

    CALL dqmc_print_complexarray(n, m, title, label, avg, err, opt)
END

SUBROUTINE detour_dqmc_print_eigenmode(n, m, title, value, opt) BIND(C, NAME="detour_dqmc_util_dqmc_print_eigenmode")
    USE detour_dqmc_util
    INTEGER(4), INTENT(IN) :: n
    INTEGER(4), INTENT(IN) :: m
    CHARACTER(LEN=*), INTENT(IN) :: title
    COMPLEX(8), INTENT(IN) :: value(:, :, :)
    INTEGER(4), INTENT(IN) :: opt

    CALL dqmc_print_eigenmode(n, m, title, value, opt)
END

SUBROUTINE detour_dqmc_getftk(value, n, nclass, class, na, nk, ft_wgt, phase, valuek) BIND(C, NAME="detour_dqmc_util_dqmc_getftk")
    USE detour_dqmc_util
    REAL(8), INTENT(IN) :: value(nclass)
    INTEGER(4), INTENT(IN) :: n
    INTEGER(4), INTENT(IN) :: nclass
    INTEGER(4), INTENT(IN) :: class(n, n)
    INTEGER(4), INTENT(IN) :: na
    INTEGER(4), INTENT(IN) :: nk
    COMPLEX(8), INTENT(IN) :: ft_wgt(n / na, nk)
    INTEGER(4), INTENT(IN) :: phase(n, n)
    COMPLEX(8), INTENT(OUT) :: valuek(nk * na * (na + 1) / 2)

    CALL dqmc_getftk(value, n, nclass, class, na, nk, ft_wgt, phase, valuek)
END

SUBROUTINE detour_dqmc_io_open(fname, inp_unit, out_unit) BIND(C, NAME="detour_dqmc_util_dqmc_io_open")
    USE detour_dqmc_util
    CHARACTER(LEN=*), INTENT(OUT) :: fname
    INTEGER(4), INTENT(OUT) :: inp_unit
    INTEGER(4), INTENT(OUT) :: out_unit

    CALL dqmc_io_open(fname, inp_unit, out_unit)
END

SUBROUTINE detour_dqmc_open_file(fname, fstatus, file_unit) BIND(C, NAME="detour_dqmc_util_dqmc_open_file")
    USE detour_dqmc_util
    CHARACTER(LEN=*), INTENT(IN) :: fname
    CHARACTER(LEN=*), INTENT(IN) :: fstatus
    INTEGER(4), INTENT(OUT) :: file_unit

    CALL dqmc_open_file(fname, fstatus, file_unit)
END

SUBROUTINE detour_dqmc_count_records(n, file_unit) BIND(C, NAME="detour_dqmc_util_dqmc_count_records")
    USE detour_dqmc_util
    INTEGER(4), INTENT(OUT) :: n
    INTEGER(4), INTENT(IN) :: file_unit

    CALL dqmc_count_records(n, file_unit)
END

FUNCTION detour_move_to_record(string, iunit) BIND(C, NAME="detour_dqmc_util_move_to_record")
    USE detour_dqmc_util
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER(4), INTENT(IN) :: iunit
    LOGICAL(4) :: detour_move_to_record

    detour_move_to_record = move_to_record(string, iunit)
END

FUNCTION detour_get_det(a) BIND(C, NAME="detour_dqmc_util_get_det")
    USE detour_dqmc_util
    REAL(8), INTENT(IN) :: a(3, 3)
    REAL(8) :: detour_get_det

    detour_get_det = get_det(a)
END

SUBROUTINE detour_get_inverse(a, inv) BIND(C, NAME="detour_dqmc_util_get_inverse")
    USE detour_dqmc_util
    REAL(8), INTENT(IN) :: a(3, 3)
    REAL(8), INTENT(OUT) :: inv(3, 3)

    CALL get_inverse(a, inv)
END
END
