module TestAllModules
    use SweepMethods_Test
    implicit none

    contains

    subroutine test_all()
        write(*,*) '����� ������'
        call progonka_test()
        Write(*,*) '����� ���������'
        return
    end subroutine

end module TestAllModules
