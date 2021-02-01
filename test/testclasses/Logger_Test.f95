module Logger_Test
    use Logger_module
    real, parameter :: eps = 0.000001
    integer, parameter :: logs = 1
    CONTAINS

    SUBROUTINE realToChar_Test()
        real, parameter :: a = 92121.32
        character(*), parameter ::  strans = '92121.32'
        character(24) ::  str
        open(logs, file = "test/testlogs.txt", status = "old")
        str = trim(realToChar(a))
        if(str(1:8) .eq. strans) then
            write(*,*) ' Logger_Test:realToChar_Test::complete'
        else
            write(logs,*) ' Logger_Test:realToChar_Test::ошибка :: ' // str
            write(*,*) ' Logger_Test:realToChar_Test::failed'
        end if
        close(logs)
    END SUBROUTINE realToChar_Test

    SUBROUTINE realToChar_Test2()
        real, parameter :: a = 0.0000001
        character(*), parameter ::  strans = '0.0000001'
        character(24) ::  str
        open(logs, file = "test/testlogs.txt", status = "old")
        str = trim(realToChar(a))
        if(str(1:9) .eq. strans) then
            write(*,*) ' Logger_Test:realToChar_Test2::complete'
        else
            write(logs,*) ' Logger_Test:realToChar_Test2::ошибка :: ' // str
            write(*,*) ' Logger_Test:realToChar_Test2::failed'
        end if
        close(logs)
    END SUBROUTINE realToChar_Test2

    SUBROUTINE intToChar_Test()
        integer, parameter :: a = 9212132
        character(*), parameter ::  strans = '9212132'
        character(24) ::  str
        open(logs, file = "test/testlogs.txt", status = "old")
        str = trim(intToChar(a))
        if(str(1:7) .eq. strans(1:7)) then
            write(*,*) ' Logger_Test:intToChar_Test::complete'
        else
            write(logs,*) ' Logger_Test:intToChar_Test::ошибка :: ' // str
            write(*,*) ' Logger_Test:intToChar_Test::failed'
        end if
        close(logs)
    END SUBROUTINE intToChar_Test

end module Logger_Test
