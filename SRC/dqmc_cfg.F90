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
END INTERFACE

INTERFACE CFG_Get
    PROCEDURE DQMC_Config_GetI, DQMC_Config_GetR
    PROCEDURE DQMC_Config_GetS, DQMC_Config_GetPR
END INTERFACE

CONTAINS

SUBROUTINE detour_dqmc_readln(str, ipt, status) BIND(C, NAME="detour_dqmc_cfg_dqmc_readln")
    USE detour_dqmc_cfg
    CHARACTER(LEN=*), INTENT(INOUT) :: str
    INTEGER(4), INTENT(IN) :: ipt
    INTEGER(4), INTENT(INOUT) :: status

    CALL dqmc_readln(str, ipt, status)
END

END
