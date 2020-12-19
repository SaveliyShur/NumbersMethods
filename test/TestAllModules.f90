MODULE TestAllModules
    use SweepMethods_Test
    implicit none

    private :: refresh_file
    public :: test_all

    CONTAINS

    SUBROUTINE test_all()
        call refresh_file()
        write(*,*) '����� ������'
        call progonka_test()
        write(*,*) '����� ���������'
        return
    END SUBROUTINE

    SUBROUTINE refresh_file()
        open(1, file = "test/testlogs.txt", status = "replace")
        close(1)
    END SUBROUTINE


END MODULE TestAllModules
