MODULE  SweepMethods_Test
    use SweepMethods
    implicit none
    real(8), parameter :: eps = 0.000001
    integer, parameter :: logs = 1

    CONTAINS

    SUBROUTINE progonka_test()
        integer :: NI
        real(8),allocatable :: A(:),B(:),C(:),D(:), U(:), ANSWER(:)
        open(logs, file = "test/testlogs.txt", status = "old")
        NI = 4
        allocate(A(1:NI))
        allocate(B(1:NI))
        allocate(C(1:NI))
        allocate(D(1:NI))
        allocate(U(1:NI))
        allocate(ANSWER(1:NI))

        ANSWER = reshape((/ 1.0263797, -0.13189857, -0.11029192, 0.45514596 /), shape(ANSWER))
        A = reshape((/ 0, 3, 6, 5 /), shape(A))
        B = reshape((/ 10, 20, 100, 10 /), shape(B))
        C = reshape((/ 2, 4, 4, 0 /), shape(C))
        D = reshape((/ 10, 0, -10, 4 /), shape(D))

        call progonka(A,B,C,D,NI,U(1:NI))

        if(abs(U(4) - ANSWER(4))>eps) then
            write(logs,*) ' SweepMethod_Test:progonka_test::Ошибка в последнем слагаемом'
            write(*,*) ' SweepMethod_Test:progonka_test::failed'
        else if(maxval(abs(U(1:4) - ANSWER(1:4))) > eps) then
            write(logs,*) 'SweepMethod_Test:progonka_test::Ошибка в ответе'
            write(*,*) ' SweepMethod_Test:progonka_test::failed'
        else
            write(*,*) ' SweepMethod_Test:progonka_test::complete'
        end if
        close(logs)
        return
    END SUBROUTINE

END MODULE
