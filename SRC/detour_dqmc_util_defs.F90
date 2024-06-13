MODULE detour_dqmc_util_defs

IMPLICIT NONE

integer,  parameter :: WP = kind(1.0d0)  ! work precision
real(WP), parameter :: ZERO = 0.0D0      ! constant 0
real(WP), parameter :: ONE  = 1.0D0      ! constant 1
real(WP), parameter :: TWO  = 2.0D0      ! constant 1
real(WP), parameter :: HALF = 0.5D0      ! constant 1

integer,  parameter :: STDERR = 0        ! standard error output
integer,  parameter :: STDOUT = 6        ! standard output
integer,  parameter :: STDIN  = 5        ! standardinput

character(*), parameter :: FMT_STRINT  = "(a30, i12)"
character(*), parameter :: FMT_STRDBL  = "(a30, f19.6)"
character(*), parameter :: FMT_STR2BL  = "(a30, '(', f11.6, ',', f11.6, ')')"
character(*), parameter :: FMT_VALERR  = "(a30, f12.6,' +- ',f12.6)"
character(*), parameter :: FMT_INTPAR  = "(i3,i3)"
character(*), parameter :: FMT_DBLINE  = "(76('='))"
character(*), parameter :: FMT_SGLINE  = "(76('-'))"
character(*), parameter :: FMT_POINT   = "('point ; dx=', i3, ' ; dy=', i3, ' :')"

! Preset parameters for dlarnv() call. These parameters were defined previously in lapack_mod.F90
! which is no longer used. 
integer, parameter :: DLARNV_UNI_0_1  = 1
integer, parameter :: DLARNV_UNI_N1_1 = 2
integer, parameter :: DLARNV_NORMAL   = 3

END MODULE
