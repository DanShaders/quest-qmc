MODULE detour_dqmc_cfg_defs
USE dqmc_util
IMPLICIT NONE

  ! string parameters
  character, parameter :: COMMENT = "#"
  character, parameter :: SEPARAT = "="
  character, parameter :: COMMA   = ","

  ! string length
  integer, parameter :: slen = 60
  integer, parameter :: llen = 256
  integer, parameter :: alen = 10    ! array limit

  ! status param
  integer, parameter :: STAT_EOF     = -1
  integer, parameter :: STAT_COMMENT = 0
  integer, parameter :: STAT_NORMAL  = 1

  ! type def
  character(8), parameter :: TYPE_STR(3) =(/"Real*8  ","Integer ","Char(60)"/)

  integer, parameter :: TYPE_REAL    = 1 
  integer, parameter :: TYPE_INTEGER = 2 
  integer, parameter :: TYPE_STRING  = 3 

  type config
     integer(8) :: cpp_data
  end type config

END MODULE
