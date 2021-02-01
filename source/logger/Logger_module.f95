MODULE Logger_module
    private
    public :: finalizeLogger
    public :: info, error, fatal
    public :: realToChar, intToChar
    integer, parameter ::  io = 2326544

    type, public :: Logger
    end type Logger

    type (Logger), pointer :: logger_single => null()

CONTAINS
    SUBROUTINE initializeLogger()
        if (.not. associated(logger_single)) then
            allocate(logger_single)
            open(io, file='source\resource\outputres\logs_logger.txt', status = 'REPLACE')
            write(io,*) 'initializeLogger'
        end if
    END SUBROUTINE initializeLogger

    SUBROUTINE info(str)
        character(*) :: str
        call writeLogs('info', str)
    END SUBROUTINE info

    SUBROUTINE error(str)
        character(*) :: str
        call writeLogs('ERROR', str)
    END SUBROUTINE error

    SUBROUTINE fatal(str)
        character(*) :: str
        call writeLogs('!!!!!!!!!FATAL!!!!!!!!!!', str)
    END SUBROUTINE fatal

    SUBROUTINE writeLogs(status, str)
        implicit none
        character(*) :: str
        character(*) :: status
        if(.not. associated(logger_single)) then
            call initializeLogger()
        end if
        write(io,*) status, ' : ', str
    END SUBROUTINE writeLogs

    SUBROUTINE finalizeLogger()
        if (associated(logger_single)) then
            write(io,*) 'finalizeLogger'
            close(io)
            deallocate(logger_single)
        end if
    END SUBROUTINE finalizeLogger

    CHARACTER(24) FUNCTION realToChar(r)
        implicit none
        real :: r
        real :: drob
        integer(8) :: drobInt, wholeInt
        drob = r - int(r)
        drobInt = int(drob * 1e8)
        wholeInt = int(r)
        realToChar = ''
        do
            realToChar = achar(mod(drobInt, 10) + 48) // realToChar
            if(drobInt/10 .eq. 0) then
                exit
            end if
            drobInt = drobInt/10
        end do
        realToChar = '.' // realToChar
        do
            realToChar = achar(mod(wholeInt, 10) + 48) // realToChar
            if(wholeInt/10 .eq. 0) then
                exit
            end if
            wholeInt = wholeInt/10
        end do
        return
    END FUNCTION realToChar

    CHARACTER(24) FUNCTION intToChar(ins)
        implicit none
        integer :: ins
        integer :: pin
        pin = ins
        intToChar = ''
        do
            intToChar = achar(mod(pin, 10) + 48) // intToChar
            if(pin/10 .eq. 0) then
                exit
            end if
            pin = int(pin/10)
        end do
        return
    END FUNCTION intToChar

END MODULE Logger_module
