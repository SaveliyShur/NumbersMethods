!Модуль для решения задачи обтекания пластины
!При решении применяется численное решение уравнений Навье-Стокса методом установления
!Входные параметры задаются в resource/Input.txt
!Решение выводиться в файл resource/data_ns.plt
MODULE MethodOfEstablishing_Plate
    implicit none

    CONTAINS

    SUBROUTINE MethodOfEstablishinglSolve_Plate()
        integer, parameter:: IO = 1, IO_Residuals = 2, loggers = 3 ! input-output unit
        real, parameter :: Eps = 3e-5
        integer NI, NJ, NITER
        integer I,J, N
        real L,H,dx,dy, visk, U0, CFL
        real dt, A, U_Residuals, V_Residuals, P_Residuals
        real,allocatable :: X_Node(:,:),Y_Node(:,:)
        real,allocatable :: X_Cell(:,:),Y_Cell(:,:)
        real,allocatable :: U(:,:),V(:,:),P(:,:)
        real,allocatable :: U_n(:,:),V_n(:,:),P_n(:,:)
        real :: Ui_p, Vi_p, Ui_m, Vi_m, Uj_p, Vj_p, Vj_m, Uj_m
        real :: U_half_ip, U_half_im, U_half_jp, U_half_jm
        real :: V_half_ip, V_half_im, V_half_jp, V_half_jm
        real :: P_half_ip, P_half_im, P_half_jp, P_half_jm

        write(*,*) 'Read input file'
        open(IO,FILE='source\resource\Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) U0
        read(IO,*) visk
        read(IO,*) NITER
        read(IO,*) CFL
        close(IO)

        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates
        allocate(X_Cell(0:NI,0:NJ)) ! Cell Centers
        allocate(Y_Cell(0:NI,0:NJ)) ! Cell Centers

        !Cell-centered variables
        allocate(U(0:NI,0:NJ))   ! Velocity U
        allocate(V(0:NI,0:NJ))   ! Velocity V
        allocate(P(0:NI,0:NJ))   ! Pressure

        !Node variables
        allocate(U_n(0:NI,0:NJ))   ! Velocity U
        allocate(V_n(0:NI,0:NJ))   ! Velocity V
        allocate(P_n(0:NI,0:NJ))   ! Pressure

        dx=L/(NI-1)
        dy=H/(NJ-1)
        A = 1/(U0**2)
        dt = CFL * min(dx/U0, min(0.5*dx*dx/visk, 0.5*dy*dy/visk))

        do I=1,NI
          do J=1,NJ
            X_Node(I,J)=(I-1)*dx
            Y_Node(I,J)=(J-1)*dy
          end do
        end do

        X_Cell(0,1:NJ)=-dx/2
        Y_Cell(0,1:NJ)=Y_Node(1,1:NJ)+dy/2
        X_Cell(1:NI,0)=X_Node(1:NI,1)+dx/2
        Y_Cell(1:NI,0)=-dy/2
        do I=1,NI
          do J=1,NJ
            X_Cell(I,J)=X_Node(I,J)+dx/2
            Y_Cell(I,J)=Y_Node(I,J)+dy/2
          end do
        end do

    !Initial field
    call InitValue(U, V, P, NI, NJ, U0)
    call BoundValue(U, V, P, NI, NJ, U0)

    !Solve equation
    open(IO_Residuals,FILE='source/resource/residuals.dat', status = "replace")
    open(loggers, FILE='source/resource/logs.txt', status = "replace")
    do N = 1, NITER
        do J = 1, NJ-1                                                          !Мэйби коректирование индексов или разбивание циклов
            do I = 1, NI-1
                !Вычисление значений в полуцелых индексах
                Ui_p = 0.5 * (U(I,J) + U(I+1,J));
                Vi_p = 0.5 * (V(I,J) + V(I+1,J));
                Ui_m = 0.5 * (U(I-1,J) + U(I,J));
                Vi_m = 0.5 * (V(I-1,J) + V(I,J));
                Uj_p = 0.5 * (U(I,J) + U(I,J+1));
                Vj_p = 0.5 * (V(I,J) + V(I,J+1));
                Uj_m = 0.5 * (U(I,J-1) + U(I,J));
                Vj_m = 0.5 * (V(I,J-1) + V(I,J));

                if(Ui_p .ge. 0) then
                    U_half_ip = U(I,J)
                    V_half_ip = V(I,J)
                    P_half_ip = P(I+1,J)
                else
                    U_half_ip = U(I+1,J)
                    V_half_ip = V(I+1,J)
                    P_half_ip = P(I,J)
                endif

                if(Ui_m .ge. 0) then
                    U_half_im = U(I-1,J)
                    V_half_im = V(I-1,J)
                    P_half_im = P(I,J)
                else
                    U_half_im = U(I,J)
                    V_half_im = V(I,J)
                    P_half_im = P(I-1,J)
                endif

                if(Vi_p .ge. 0) then
                    U_half_jp = U(I,J)
                    V_half_jp = V(I,J)
                    P_half_ip = P(I,J+1)
                else
                    U_half_jp = U(I,J+1)
                    V_half_jp = V(I,J+1)
                    P_half_jp = P(I,J)
                endif

                if(Vi_m .ge. 0) then
                    U_half_jm = U(I,J-1)
                    V_half_jm = V(I,J-1)
                    P_half_jm = P(I,J)
                else
                    U_half_jm = U(I,J)
                    V_half_jm = V(I,J)
                    P_half_jm = P(I,J-1)
                endif

                !На границах
                if(J .eq. 1) then
                     V_half_jm = Vj_m
                     U_half_jm = Uj_m
                endif

!                if(( J .eq. 1) .and. (I .eq. 20)) then
!                    write(loggers, *) V_half_jm, U_half_jm
!                endif

                !Вычисление давления
                P_n(I,J) = P(I,J) - (dt/A)*( (U_half_ip - U_half_im)/dx + (V_half_jp - V_half_jm)/dy )

                !Высление продольной компоненты скорости
                U_n(I,J) = U(I,J) - dt*( (Ui_p *U_half_ip - Ui_m*U_half_ip)/dx &
                & + (Vj_p * U_half_jp - Vj_m * U_half_jm)/dy &
                & + (P_half_ip - P_half_im)/dx &
                & - (visk*(U(I+1,J) - U(I,J))/dx - visk*(U(I,J) - U(I-1,J))/dx)/dx &
                & - (visk*(U(I,J+1) - U(I,J))/dy - visk*(U(I,J) - U(I,J-1))/dy)/dy )

                !Вычисление поперечной компоненты скорости
                V_n(I,J) = V(I,J) - dt*( (Ui_p * V_half_ip - Ui_m * V_half_im)/dx &
                & + (Vj_p * V_half_jp - Vj_m * V_half_jm)/dy &
                & + (P_half_jp - P_half_jm)/dy &
                & - (visk*(V(I+1,J) - V(I,J))/dx - visk*(V(I,J) - V(I-1,J))/dx )/dx &
                & - (visk*(V(I,J+1) - V(I,J))/dy - visk*(V(I,J) - V(I,J-1))/dy )/dy )

                !Пересчет в граничных ячейках
                call BoundValue(U_n, V_n, P_n, NI, NJ, U0)
            end do
        end do

        !Проверяем сходимость и выводим невязки
        U_Residuals = maxval(abs(U_n-U))/maxval(abs(U_n))
        V_Residuals = maxval(abs(V_n-V))/maxval(abs(V_n))
        P_Residuals = maxval(abs(P_n-P))/maxval(abs(P_n))
        if( (U_Residuals.le.Eps ) .and. (V_Residuals.le.Eps ) .and. (P_Residuals.le.Eps ) ) then
            write(*,*) "MethodOfEstablishinglSolve_Plate:Complete"
            exit
        endif
        write(*,*) "N=", N
        write(IO_Residuals, *) dt*N, U_Residuals, V_Residuals, P_Residuals

        if(N .eq. 64) then
            call writeAnswer(IO,NI,NJ,X_Cell,Y_Cell,U,V,P)
            exit
        endif

        !Переопределяем для следующего шага
        U = U_n
        V = V_n
        P = P_n

        write(loggers, *) U(0,0), U(1,0), U(0,1), U(1,1)
    end do
    close(IO_Residuals)
    close(loggers)

    !call writeAnswer(IO,NI,NJ,X_Cell,Y_Cell,U_n,V_n,P_n)
    return
    END SUBROUTINE MethodOfEstablishinglSolve_Plate

    !Функция пишет ответ в файл
    SUBROUTINE writeAnswer(IO,NI,NJ,X,Y,U,V,P)
       implicit none

       integer NI,NJ,IO
       real, dimension(NI,NJ):: X,Y
       real, dimension(0:NI,0:NJ)::U,V,P

       write(*,*) 'Output data cell (Navier - Stokes) '
       open(IO,FILE='source/resource/data_ns.plt', status = "replace")
       call OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
       close(IO)
    END SUBROUTINE writeAnswer

    !Функция для вывода в формате техплот
    SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
        implicit none

        integer NI,NJ,IO
        real, dimension(NI,NJ):: X,Y
        real, dimension(0:NI,0:NJ)::U,V,P

        write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
        write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-5]=CELLCENTERED)'
        write(IO,'(100E25.16)') X(1:NI,1:NJ)
        write(IO,'(100E25.16)') Y(1:NI,1:NJ)
        write(IO,'(100E25.16)') U(1:NI-1,1:NJ-1)
        write(IO,'(100E25.16)') V(1:NI-1,1:NJ-1)
        write(IO,'(100E25.16)') P(1:NI-1,1:NJ-1)

    END SUBROUTINE OutputFields_Cell

    SUBROUTINE BoundValue(U,V,P,NI,NJ,U0)
        implicit none

        integer NI,NJ
        real :: U(0:NI,0:NJ), V(0:NI,0:NJ), P(0:NI,0:NJ)
        real U0
        !нижняя непроницаемая граница
        U(1:NI-1,0) = - U(1:NI-1,1)
        V(1:NI-1,0) = - V(1:NI-1,1)
        P(1:NI-1,0) = P(1:NI-1,1)

        !левая граница, вход
        U(0,0:NJ) =  U0
        V(0,0:NJ) = 0
        P(0,0:NJ) = P(1,0:NJ)

        !правая граница, выход
        U(NI,0:NJ) = U(NI-1,0:NJ)
        V(NI,0:NJ) = V(NI-1,0:NJ)
        P(NI,0:NJ) = 0

        !верхняя граница, выход
        U(1:NI-1,NJ) = U(1:NI-1,NJ-1)
        V(1:NI-1,NJ) = V(1:NI-1,NJ-1)
        P(1:NI-1,NJ) = 0
    END SUBROUTINE BoundValue

    SUBROUTINE InitValue(U,V,P,NI,NJ,U0)
        implicit none

        real :: U(0:NI,0:NJ), V(0:NI,0:NJ),P(0:NI,0:NJ)
        real U0
        integer NI,NJ

        U(0:NI,0:NJ) = U0
        V(0:NI,0:NJ) = 0
        P(0:NI,0:NJ) = 0

    END SUBROUTINE InitValue

END MODULE MethodOfEstablishing_Plate
