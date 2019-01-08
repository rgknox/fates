module FatesGlobals
  ! NOTE(bja, 201608) This is a temporary hack module to store global
  ! data used inside fates. It's main use it to explicitly call out
  ! global data that needs to be dealt with, but doesn't have an
  ! immediately obvious home.

  use FatesConstantsMod         , only : r8 => fates_r8
   
  implicit none

  public :: FatesGlobalsInit
  public :: fates_log
  public :: fates_global_verbose

  integer, private :: fates_log_
  logical, private :: fates_global_verbose_

contains



  ! =====================================================================================

  subroutine FatesGlobalsInit(log_unit,global_verbose)

    implicit none

    integer, intent(in) :: log_unit
    logical, intent(in) :: global_verbose

    fates_log_ = log_unit
    fates_global_verbose_ = global_verbose

  end subroutine FatesGlobalsInit

  ! =====================================================================================

  integer function fates_log()
    fates_log = fates_log_
  end function fates_log

  logical function fates_global_verbose()
    fates_global_verbose = fates_global_verbose_
  end function fates_global_verbose

  subroutine fates_endrun(msg) 

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    ! This subroutine was derived from CLM's
    ! endrun_vanilla() in abortutils.F90
    !
    use shr_sys_mod , only: shr_sys_abort
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: msg    ! string to be printed
    !-----------------------------------------------------------------------

    write(fates_log(),*)'ENDRUN:', msg
    call shr_sys_abort()

  end subroutine fates_endrun

  ! =====================================================================================

  subroutine err_msg_supp1()

    write(fates_log(),*)'                                                                             '
    write(fates_log(),*)'                                                                             '
    write(fates_log(),*)'      ___  ___    __    _____          ____     __    _____  _____           '
    write(fates_log(),*)'     /     |  \  /  \   |     |   |    |   \   /  \   |      |       |       '
    write(fates_log(),*)'     |     |--| /----\  ----| |---|    ----   /____\  ----|  ----|   |       '
    write(fates_log(),*)'     \___  |  |/      \ ____| |   |    |___/ /      \ ____|  ____|           '
    write(fates_log(),*)'                                                                     O       '
    write(fates_log(),*)''         
    write(fates_log(),*)''   
    write(fates_log(),*)''   
    write(fates_log(),*)''   
    write(fates_log(),*)''
    write(fates_log(),*)''
    write(fates_log(),*)''
    write(fates_log(),*)'                                       /..:-                                     '
    write(fates_log(),*)'                                 -.-o+++o+++           ``.--::////////::--.--.   '
    write(fates_log(),*)'                               .-+ooo+oosyos. `-:/osyhddhyyysssyssosssyyhs++o/+. '
    write(fates_log(),*)'                              `++oooosssyyhhhhhysssoooo++++++os+soy+++ss++o:---/ '
    write(fates_log(),*)'                             -++ooosyhhdhysso++++++++++++++++os+y+s++so++o:----+ '
    write(fates_log(),*)'                            .+oosyhddhs+++++++++++++++++++++++ooosooso++s:-----/ '
    write(fates_log(),*)'                           `+oyyddhso+++++++++++++++++++++++++++++os+++s/-----:. '
    write(fates_log(),*)'                           :shdyoo++++++++++++++++++++++++++++++oso+++os------:  '
    write(fates_log(),*)'                 `.-:/++osyddyo+++++++++++++++++++++++o++++++++oo+++++h:---::/   '
    write(fates_log(),*)'              `:ossossssyddyo+++++++++++++++++++++++++s+++++++oo+++++so-----::-  '
    write(fates_log(),*)'             -soooosyyyddhs+++++++++++++++++++++++++++oo++++++so++ooosyo+:-----/`'
    write(fates_log(),*)'            -oooosyyyhddyo+++++++++++oooooo++++++o+++++ooooo+++o+++++++++oso+//+o'
    write(fates_log(),*)'            ooosssyshdhso++++++++++++s+++osssso++so++++++++++++++++++oo++++++++o/'
    write(fates_log(),*)'           :ooooyyyddho++++++++++++++s++oossyyyhyoso+++++++++++++++++++//+++/-.  '
    write(fates_log(),*)'           /osoyysydho+++++++++++++++sooossssshhds+++oo+++++++++++++o:--:-.`     '
    write(fates_log(),*)'           .oos+-/dho+++++++++++++++++ssosssyyyhys/////++oo++++++////+-`         '
    write(fates_log(),*)'             `  :dho+++++++++++++++////+ooooo++::----------:://:----.            '
    write(fates_log(),*)'               `dhs+++++///:::::::::---------------:oo++//-..                    '
    write(fates_log(),*)'               shy++++//:/ossyyysssosoo+---------:oyyyshys                       '
    write(fates_log(),*)'              +dyo++//::ossssssoooooos/         `ssyysohyo                       '
    write(fates_log(),*)'            -shy++//:-: `oosssoosoos/`          :ooo++ysoo                       '
    write(fates_log(),*)'        -/+ooso+++/--/    -+ossso/-             .++++.`:/:                       '
    write(fates_log(),*)'     `:oooooooo++++/:/:                          .-`                             '
    write(fates_log(),*)'    :oo+ooosooooooooooy/                                                         '
    write(fates_log(),*)'  `+o+oo+/oooso+sssossyy/                                                        '
    write(fates_log(),*)' .+o+o+++o+s+sssoysoyyyss:                                                       '
    write(fates_log(),*)' +s+o+/o+soyoyoyoyssyssyys`                                                      '
    write(fates_log(),*)'  `.--///+::.  -osssssyssy/                                                      '
    write(fates_log(),*)'                `:/+sysyy/                                                      '
    write(fates_log(),*)''
    write(fates_log(),*)''
    
    return
  end subroutine err_msg_supp1




end module FatesGlobals
