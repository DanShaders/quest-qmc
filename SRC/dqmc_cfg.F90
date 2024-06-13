MODULE dqmc_cfg

USE dqmc_util
USE detour_dqmc_cfg_defs

IMPLICIT NONE

INTERFACE
    SUBROUTINE dqmc_readln(str, ipt, status) BIND(C, NAME="fortran_dqmc_cfg_dqmc_readln")
        USE detour_dqmc_cfg_defs
        CHARACTER(LEN=*), INTENT(INOUT) :: str
        INTEGER(4), INTENT(IN) :: ipt
        INTEGER(4), INTENT(INOUT) :: status
    END

    SUBROUTINE dqmc_default_def(cfg) BIND(C, NAME="fortran_dqmc_cfg_dqmc_default_def")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(INOUT) :: cfg
    END

    SUBROUTINE dqmc_config_free(cfg) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_free")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(INOUT) :: cfg
    END

    SUBROUTINE dqmc_read_def(cfg, ipt) BIND(C, NAME="fortran_dqmc_cfg_dqmc_read_def")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(INOUT) :: cfg
        INTEGER(4), INTENT(IN) :: ipt
    END

    SUBROUTINE dqmc_print_def(cfg, opt) BIND(C, NAME="fortran_dqmc_cfg_dqmc_print_def")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(INOUT) :: cfg
        INTEGER(4), INTENT(IN) :: opt
    END

    FUNCTION dqmc_find_param(cfg, pname) BIND(C, NAME="fortran_dqmc_cfg_dqmc_find_param")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(IN) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: pname
        INTEGER(4) :: dqmc_find_param
    END

    SUBROUTINE dqmc_read_config(cfg) BIND(C, NAME="fortran_dqmc_cfg_dqmc_read_config")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(INOUT) :: cfg
    END

    SUBROUTINE dqmc_config_seti(cfg, name, value) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_seti")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(INOUT) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(4), INTENT(IN) :: value
    END

    SUBROUTINE dqmc_config_setr(cfg, name, value) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_setr")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(INOUT) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        REAL(8), INTENT(IN) :: value
    END

    SUBROUTINE dqmc_config_setpr(cfg, name, n, value) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_setpr")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(INOUT) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(4), INTENT(IN) :: n
        REAL(8), INTENT(IN) :: value(n)
    END

    SUBROUTINE dqmc_config_setpi(cfg, name, n, value) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_setpi")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(INOUT) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(4), INTENT(IN) :: n
        INTEGER(4), INTENT(IN) :: value(n)
    END

    SUBROUTINE dqmc_config_sets(cfg, name, value) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_sets")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(INOUT) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        CHARACTER(LEN=*), INTENT(IN) :: value
    END

    FUNCTION dqmc_config_isset(cfg, name) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_isset")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(IN) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        LOGICAL(4) :: dqmc_config_isset
    END

    SUBROUTINE dqmc_config_geti(cfg, name, value) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_geti")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(IN) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(4), INTENT(OUT) :: value
    END

    SUBROUTINE dqmc_config_getr(cfg, name, value) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_getr")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(IN) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        REAL(8), INTENT(OUT) :: value
    END

    SUBROUTINE dqmc_config_getpr(cfg, name, n, value) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_getpr")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(IN) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(4), INTENT(OUT) :: n
        REAL(8), INTENT(INOUT), POINTER :: value(:)
    END

    SUBROUTINE dqmc_config_getpi(cfg, name, n, value) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_getpi")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(IN) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(4), INTENT(OUT) :: n
        INTEGER(4), INTENT(INOUT), POINTER :: value(:)
    END

    SUBROUTINE dqmc_config_gets(cfg, name, value) BIND(C, NAME="fortran_dqmc_cfg_dqmc_config_gets")
        USE detour_dqmc_cfg_defs
        TYPE(config), INTENT(IN) :: cfg
        CHARACTER(LEN=*), INTENT(IN) :: name
        CHARACTER(LEN=*) :: value
    END
END INTERFACE

INTERFACE CFG_Set
    PROCEDURE DQMC_Config_SetI, DQMC_Config_SetR
    PROCEDURE DQMC_Config_SetS, DQMC_Config_SetPR
    PROCEDURE DQMC_Config_SetPI
END INTERFACE

INTERFACE CFG_Get
    PROCEDURE DQMC_Config_GetI, DQMC_Config_GetR
    PROCEDURE DQMC_Config_GetS, DQMC_Config_GetPR
    PROCEDURE DQMC_Config_GetPI
END INTERFACE

CONTAINS
SUBROUTINE detour_dqmc_readln(str, ipt, status) BIND(C, NAME="detour_dqmc_cfg_dqmc_readln")
    USE detour_dqmc_cfg
    CHARACTER(LEN=*), INTENT(INOUT) :: str
    INTEGER(4), INTENT(IN) :: ipt
    INTEGER(4), INTENT(INOUT) :: status

    CALL dqmc_readln(str, ipt, status)
END

SUBROUTINE detour_dqmc_default_def(cfg) BIND(C, NAME="detour_dqmc_cfg_dqmc_default_def")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(INOUT) :: cfg

    CALL dqmc_default_def(cfg)
END

SUBROUTINE detour_dqmc_config_free(cfg) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_free")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(INOUT) :: cfg

    CALL dqmc_config_free(cfg)
END

SUBROUTINE detour_dqmc_read_def(cfg, ipt) BIND(C, NAME="detour_dqmc_cfg_dqmc_read_def")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(INOUT) :: cfg
    INTEGER(4), INTENT(IN) :: ipt

    CALL dqmc_read_def(cfg, ipt)
END

SUBROUTINE detour_dqmc_print_def(cfg, opt) BIND(C, NAME="detour_dqmc_cfg_dqmc_print_def")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(INOUT) :: cfg
    INTEGER(4), INTENT(IN) :: opt

    CALL dqmc_print_def(cfg, opt)
END

FUNCTION detour_dqmc_find_param(cfg, pname) BIND(C, NAME="detour_dqmc_cfg_dqmc_find_param")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(IN) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: pname
    INTEGER(4) :: detour_dqmc_find_param

    detour_dqmc_find_param = dqmc_find_param(cfg, pname)
END

SUBROUTINE detour_dqmc_read_config(cfg) BIND(C, NAME="detour_dqmc_cfg_dqmc_read_config")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(INOUT) :: cfg

    CALL dqmc_read_config(cfg)
END

SUBROUTINE detour_dqmc_config_seti(cfg, name, value) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_seti")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(INOUT) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(4), INTENT(IN) :: value

    CALL dqmc_config_seti(cfg, name, value)
END

SUBROUTINE detour_dqmc_config_setr(cfg, name, value) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_setr")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(INOUT) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(8), INTENT(IN) :: value

    CALL dqmc_config_setr(cfg, name, value)
END

SUBROUTINE detour_dqmc_config_setpr(cfg, name, n, value) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_setpr")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(INOUT) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(4), INTENT(IN) :: n
    REAL(8), INTENT(IN) :: value(n)

    CALL dqmc_config_setpr(cfg, name, n, value)
END

SUBROUTINE detour_dqmc_config_setpi(cfg, name, n, value) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_setpi")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(INOUT) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(4), INTENT(IN) :: n
    INTEGER(4), INTENT(IN) :: value(n)

    CALL dqmc_config_setpi(cfg, name, n, value)
END

SUBROUTINE detour_dqmc_config_sets(cfg, name, value) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_sets")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(INOUT) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*), INTENT(IN) :: value

    CALL dqmc_config_sets(cfg, name, value)
END

FUNCTION detour_dqmc_config_isset(cfg, name) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_isset")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(IN) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL(4) :: detour_dqmc_config_isset

    detour_dqmc_config_isset = dqmc_config_isset(cfg, name)
END

SUBROUTINE detour_dqmc_config_geti(cfg, name, value) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_geti")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(IN) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(4), INTENT(OUT) :: value

    CALL dqmc_config_geti(cfg, name, value)
END

SUBROUTINE detour_dqmc_config_getr(cfg, name, value) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_getr")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(IN) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(8), INTENT(OUT) :: value

    CALL dqmc_config_getr(cfg, name, value)
END

SUBROUTINE detour_dqmc_config_getpr(cfg, name, n, value) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_getpr")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(IN) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(4), INTENT(OUT) :: n
    REAL(8), INTENT(INOUT), POINTER :: value(:)

    CALL dqmc_config_getpr(cfg, name, n, value)
END

SUBROUTINE detour_dqmc_config_getpi(cfg, name, n, value) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_getpi")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(IN) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(4), INTENT(OUT) :: n
    INTEGER(4), INTENT(INOUT), POINTER :: value(:)

    CALL dqmc_config_getpi(cfg, name, n, value)
END

SUBROUTINE detour_dqmc_config_gets(cfg, name, value) BIND(C, NAME="detour_dqmc_cfg_dqmc_config_gets")
    USE detour_dqmc_cfg
    TYPE(config), INTENT(IN) :: cfg
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*) :: value

    CALL dqmc_config_gets(cfg, name, value)
END
END
