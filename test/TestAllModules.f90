MODULE TestAllModules
    use SweepMethods_Test
    use Logger_Test
    implicit none

    private :: refresh_file
    public :: test_all

    CONTAINS

    SUBROUTINE test_all()
        call refresh_file()
        write(*,*) 'Тесты начаты'
        call progonka_test()
        call realToChar_Test()
        call intToChar_Test()
        call realToChar_Test2()
        write(*,*) 'Тесты завершены'
        return
    END SUBROUTINE

    SUBROUTINE refresh_file()
        open(1, file = "test/testlogs.txt", status = "replace")
        close(1)
    END SUBROUTINE


END MODULE TestAllModules
