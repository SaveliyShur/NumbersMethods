module NavieStocksGas
    use Logger_module
    implicit none
    CONTAINS

    SUBROUTINE NavieStocksGas_Solver()
        use omp_lib
        implicit none

        integer, parameter:: IO = 565, IO_Residuals = 1132, loggers = 498
        real(8) :: Eps
        integer NI, NJ, NITER, ios
        integer I,J, N, num
        real(8) L,H,dx,dy, visk, U0, CFL, Diametr, Cc, gamm, ro0, Uoe, Ue
        real(8) dt, A, U_Residuals, V_Residuals, P_Residuals
        real(8),allocatable :: X_Cell(:),Y_Cell(:)
        real(8),allocatable :: U(:,:),V(:,:),P(:,:), ro(:,:)
        real(8),allocatable :: U_n(:,:),V_n(:,:),P_n(:,:), ro_n(:,:)
        real(8), allocatable:: U_cap(:,:), V_cap(:,:)
        real(8), allocatable:: U_i_half(:,:), V_i_half(:,:), P_i_half(:,:), ro_i_half(:,:)
        real(8), allocatable:: U_j_half(:,:), V_j_half(:,:), P_j_half(:,:), ro_j_half(:,:)

        call info('Read projects settings')
        open(IO,file='projectsettings.txt', STATUS='OLD', IOSTAT=ios)
        if(ios/=0) then
            call fatal('projects settings no found')
            stop 1
        end if
        read(IO,*) Eps
        read(IO,*) NITER
        read(IO,*) num
        close(IO)
        call info('Read projects settings :: complete')

        call info('Read input file')
        open(IO,file='resource\inputres\Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) U0
        read(IO,*) visk
        read(IO,*) CFL
        read(IO,*) Diametr
        read(IO,*) Cc
        read(IO,*) gamm
        read(IO,*) ro0
        read(IO,*) Uoe
        close(IO)
        call info('Read input file :: complete')

        allocate(X_Cell(0:NI))
        allocate(Y_Cell(0:NJ))

        allocate(U(0:NI,0:NJ))
        allocate(V(0:NI,0:NJ))
        allocate(P(0:NI,0:NJ))
        allocate(ro(0:NI,0:NJ))

        allocate(U_n(0:NI,0:NJ))
        allocate(V_n(0:NI,0:NJ))
        allocate(P_n(0:NI,0:NJ))
        allocate(ro_n(0:NI,0:NJ))

        allocate(U_cap(0:NI,0:NJ))
        allocate(V_cap(0:NI,0:NJ))
        allocate(U_i_half(0:NI,0:NJ), V_i_half(0:NI,0:NJ), P_i_half(0:NI,0:NJ), ro_i_half(0:NI,0:NJ))
        allocate(U_j_half(0:NI,0:NJ), V_j_half(0:NI,0:NJ), P_j_half(0:NI,0:NJ), ro_j_half(0:NI,0:NJ))

        dx=L/(NI-1)
        dy=H/(NJ-1)
        A = 1/(U0**2)
        dt = CFL * min(dx/U0, min(0.5*dx*dx/visk, 0.5*dy*dy/visk))
        Ue = U0/Uoe
        call info('dt=' // trim(realToChar(dt)))
        call omp_set_num_threads(num)

        do I=0,NI
            X_Cell(I)=(I-0.5)*dx
        end do
        do J=0,NJ
            Y_Cell(J)=(J-0.5)*dy
        end do


    call info('start NavieStocksGas_Solver solver for gas, parameters: eps=' // trim(realToChar(Eps)) // ' NITER='&
    & // trim(intToChar(NITER)) // ' numbers thread=' // trim(intToChar(num))  )

    !Initial field
    call InitValue(U, V, P, ro, NI, NJ, U0, Ue, ro0, Cc, gamm)
    call BoundValue(U, V, P, ro, NI, NJ, U0, Cc, ro0, gamm, H, Diametr, Ue)
    U_n = U
    V_n = V
    P_n = P

    call info('Write init data in data_ns_init.tec')
    open(IO,FILE='resource/outputres/data_ns_init.tec')
    call OutputFields_Cell(IO,NI,NJ,X_Cell,Y_Cell,U_n,V_n,P_n)
    close(IO)
    call info('Write data init :: Complete')

    !Solve equation
    open(IO_Residuals,FILE='resource/outputres/residuals.dat', status = "replace")

    do N = 1, NITER
        !Вычисление значений в полуцелых индексах
        !$omp parallel
        !$omp do
        do I = 0, NI-1
            do J = 1, NJ
                U_cap(I,J) = 0.5 * (U(I+1,J) + U(I,J))
                if (U_cap(I,J) .ge. 0.0) then
                    U_i_half(I,J) = U(I,J)
                    V_i_half(I,J) = V(I,J)
                    P_i_half(I,J) = P(I+1,J)
                    ro_i_half(I,J) = ro(I+1,J)
                else
                    U_i_half(I,J) = U(I+1,J)
                    V_i_half(I,J) = V(I+1,J)
                    P_i_half(I,J) = P(I,J)
                    ro_i_half(I,J) = ro(I,J)
                end if
            end do
        end do
        !$omp end do
        !$omp do
        do I = 1, NI
            do J = 0, NJ-1
                V_cap(I,J) = 0.5 * (V(I,J) + V(I,J+1))
                if (V_cap(I,J) .ge. 0.0) then
                    U_j_half(I,J) = U(I,J)
                    V_j_half(I,J) = V(I,J)
                    P_j_half(I,J) = P(I,J+1)
                    ro_j_half(I,J) = ro(I,J+1)
                else
                    U_j_half(I,J) = U(I,J+1)
                    V_j_half(I,J) = V(I,J+1)
                    P_j_half(I,J) = P(I,J)
                    ro_j_half(I,J) = ro(I,J)
                end if
            end do
        end do
        !$omp end do
        !$omp barrier
        !$omp do
        !Вычисление компонент решения
        do J = 1, NJ-1
            do I = 1, NI-1

                !Вычисление давления
                P_n(I,J) = P(I,J) - dt * (U0**2) *( (U_i_half(I,J)*ro_i_half(I,J) - U_i_half(I-1,J)*ro_i_half(I-1,J))/dx &
                & + (V_j_half(I,J)*ro_j_half(I,J) - V_j_half(I,J-1)*ro_j_half(I,J-1))/dy )

                !Высление продольной компоненты скорости
                U_n(I,J) = U(I,J) - dt/ro(I,J) &
                & *( (U_cap(I,J)*U_i_half(I,J)*ro_i_half(I,J) - U_cap(I-1,J)*U_i_half(I-1,J)*ro_i_half(I-1,J))/dx &
                & + (V_cap(I,J)*U_j_half(I,J)*ro_j_half(I,J) - V_cap(I,J-1)*U_j_half(I,J-1)*ro_j_half(I,J-1))/dy &
                & + (P_i_half(I,J) - P_i_half(I-1,J))/dx &
                & - 4.0/3.0 *visk*(U(I+1,J) - 2*U(I,J) + U(I-1,J))/(dx**2) &
                & - visk*(U(I,J+1) - 2*U(I,J) + U(I,J-1))/(dy**2) )

                !Вычисление поперечной компоненты скорости
                V_n(I,J) = V(I,J) - dt/ro(I,J) &
                & *( (V_cap(I,J)*V_j_half(I,J)*ro_j_half(I,J) - V_cap(I,J-1)*V_j_half(I,J-1)*ro_j_half(I,J-1))/dy &
                & + (U_cap(I,J)*V_i_half(I,J)*ro_i_half(I,J) - U_cap(I-1,J)*V_i_half(I-1,J)*ro_i_half(I-1,J))/dx &
                & + (P_j_half(I,J) - P_j_half(I,J-1))/dy &
                & - 4.0/3.0* visk*(V(I+1,J) - 2*V(I,J) + V(I-1,J))/(dx**2) &
                & - visk*(V(I,J+1) - 2*V(I,J) + V(I,J-1))/(dy**2) )

                ro_n(I,J) = (P(I,J)/Cc)**(1.0/gamm)

            end do
        end do
        !$omp end do
        !$omp end parallel

        !Пересчет в граничных ячейках
        call BoundValue(U_n, V_n, P_n, ro_n, NI, NJ, U0, Cc, ro0, gamm, H, Diametr, Ue)

        !Проверяем сходимость и выводим невязки
        U_Residuals = maxval(abs(U_n-U))/maxval(abs(U_n))
        V_Residuals = maxval(abs(V_n-V))/maxval(abs(V_n))
        P_Residuals = maxval(abs(P_n-P))/maxval(abs(P_n))
        if( (U_Residuals.le.Eps ) .and. (V_Residuals.le.Eps ) .and. (P_Residuals.le.Eps ) ) then
            write(*,*) "NavieStocksGas_Solver : Complete"
            call info('NavieStocksGas_Solver solver for liquid :: Complete')
            exit
        endif
        if(MOD(N,100) .eq. 0) then
            write(*,*) "N=", N, "eps=", max(U_Residuals, V_Residuals, P_Residuals)
        end if

        write(IO_Residuals, *) dt*N, U_Residuals, V_Residuals, P_Residuals

        !Переопределяем для следующего шага
        U = U_n
        V = V_n
        P = P_n
        ro = ro_n

        if(N .eq. NITER) then
            call error('NavieStocksGas_Solver solver for gas :: Solution underreported, eps=' //&
            & realToChar(max(U_Residuals, V_Residuals, P_Residuals)))
        end if
    end do

    close(IO_Residuals)

    call writeAnswer(IO,NI,NJ,X_Cell,Y_Cell,U_n,V_n,P_n)
    return
    END SUBROUTINE NavieStocksGas_Solver

    !Функция пишет ответ в файл
    SUBROUTINE writeAnswer(IO,NI,NJ,X,Y,U,V,P)
       implicit none

       integer NI,NJ,IO
       real(8), dimension(NI):: X
       real(8), dimension(NJ):: Y
       real(8), dimension(0:NI,0:NJ)::U,V,P
       call info('Write answer')
       write(*,*) 'Output data cell (Navier - Stokes) '
       open(IO,FILE='resource/outputres/data_ns.tec', status = "replace")
       call OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
       close(IO)
       call info('Write answer :: Complete')
    END SUBROUTINE writeAnswer

    !Функция для вывода в формате техплот
    SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
        implicit none

        integer NI,NJ,IO, i, j
        real(8), dimension(NI):: X
        real(8), dimension(NJ):: Y
        real(8), dimension(0:NI,0:NJ)::U,V,P

        U(1,:) = (U(1,:) + U(0,:))/2.0
        V(1,:) = (V(1,:) + V(0,:))/2.0
        P(1,:) = (P(1,:) + P(0,:))/2.0

        U(NI,:) = (U(NI,:) + U(NI-1,:))/2.0
        V(NI,:) = (V(NI,:) + V(NI-1,:))/2.0
        P(NI,:) = (P(NI,:) + P(NI-1,:))/2.0

        U(:,NJ) = (U(:,NJ-1) + U(:,NJ))/2.0
        V(:,NJ) = (V(:,NJ-1) + V(:,NJ))/2.0
        P(:,NJ) = (P(:,NJ-1) + P(:,NJ))/2.0

        U(:,1) = (U(:,1) + U(:,0))/2.0
        V(:,1) = (V(:,1) + V(:,0))/2.0
        P(:,1) = (P(:,1) + P(:,0))/2.0
        write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
        write(IO,*) 'ZONE I=', NI, ', J=', NJ
        do i = 1, NI
            do j = 1, NJ
                write(IO,*) X(I), Y(J), U(I,J), V(I,J), P(I,J)
            end do
        end do

    END SUBROUTINE OutputFields_Cell

    SUBROUTINE BoundValue(U,V,P,ro,NI,NJ,U0, C, ro0, gamm, H, Diametr, Ue)
        implicit none

        integer NI,NJ
        real(8) :: U(0:NI,0:NJ), V(0:NI,0:NJ), P(0:NI,0:NJ), ro(0:NI,0:NJ)
        real(8) :: U0, Ue, C, ro0, gamm, H, Diametr

        !левая граница
        U(0,1:nint(Diametr/H*NJ)) = U0
        U(0,(nint(Diametr/H*NJ)+1):NJ) = Ue
        V(0,1:NJ) = 0.0
        P(0,1:NJ) = P(1,1:NJ)
        ro(0,1:NJ) = ro(1,1:NJ)

        !правая граница
        U(NI,1:NJ) = U(NI-1,1:NJ)
        V(NI,1:NJ) = V(NI-1,1:NJ)
        P(NI,1:NJ) = C*(ro0**gamm)
        ro(NI,1:NJ) = ro0

        !нижняя граница
        U(1:NI,0) = - U(1:NI,1)
        V(1:NI,0) = - V(1:NI,1)
        P(1:NI,0) = P(1:NI,1)
        ro(1:NI,0) = ro(1:NI,1)

        !верхняя граница
        U(1:NI,NJ) = U(1:NI,NJ-1)
        V(1:NI,NJ) = V(1:NI,NJ-1)
        P(1:NI,NJ) = P(1:NI,NJ-1)
        ro(1:NI,NJ) = ro(1:NI,NJ-1)
    END SUBROUTINE BoundValue

    SUBROUTINE InitValue(U,V,P,ro,NI,NJ,U0,Ue, ro0, Cc, gamm)
        implicit none

        real(8) :: U(0:NI,0:NJ), V(0:NI,0:NJ),P(0:NI,0:NJ), ro(0:NI,0:NJ)
        real(8) U0, Ue, Cc, gamm, ro0
        integer NI,NJ

        U = Ue
        V = 0.0
        P = Cc*(ro0**gamm)
        ro = ro0
    END SUBROUTINE InitValue
end module NavieStocksGas
