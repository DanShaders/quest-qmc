module detour_DQMC_Cfg

  use DQMC_util
  use detour_dqmc_cfg_defs
  implicit none 

  ! 
  ! This module contains subroutines to read input parameters.
  !

contains

  !---------------------------------------------------------------------!
  subroutine DQMC_ReadLn(str, IPT, status)
    !
    ! Purpose
    ! =======
    !    This subrotine reads in a line from file.
    !    It will get rid of #.
    !
    ! Arguments
    ! =========
    !
    character(len=llen)  :: str
    integer              :: status
    intent(inout)        :: str, status
    integer, intent(in)  :: IPT

    ! ... Local Variables ...
    integer                :: ios, pos

    ! ... Executable ...

    read (unit=IPT, FMT="(a)", iostat=ios)  str
    status = STAT_COMMENT

    ! end of file
    if (ios .ne. 0) then
       status = STAT_EOF
       return
    end if

    ! find comment # and get rid of the tailing part
    pos = scan(str, COMMENT, .false.)
    if (pos .ne. 0) then
       ! find the comment sign
       if (pos .ge. 2) then
          str = str(1:pos-1)
       else
          str = ""
       end if
    end if

    if (len_trim(str) .gt. 0) then
       status = STAT_NORMAL
    end if

  end subroutine DQMC_ReadLn

  !---------------------------------------------------------------------!

end module detour_DQMC_Cfg
