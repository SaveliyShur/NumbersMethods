module NavieStocksGas
    use Logger_module
    use Listener
    use Writer
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
        real(8),allocatable :: X_Cell(:,:),Y_Cell(:,:)
        real(8),allocatable :: X_Node(:,:),Y_Node(:,:)
        real(8),allocatable :: U(:,:),V(:,:),P(:,:), ro(:,:)
        real(8),allocatable :: U_n(:,:),V_n(:,:),P_n(:,:), ro_n(:,:)
        real(8),allocatable :: U_Node(:,:),V_Node(:,:),P_Node(:,:), ro_Node(:,:)
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

        allocate(X_Cell(0:NI,0:NJ))
        allocate(Y_Cell(0:NI,0:NJ))
        allocate(X_Node(1:NI,1:NJ))
        allocate(Y_Node(1:NI,1:NJ))

        allocate(U(0:NI,0:NJ))
        allocate(V(0:NI,0:NJ))
        allocate(P(0:NI,0:NJ))
        allocate(ro(0:NI,0:NJ))

        allocate(U_n(0:NI,0:NJ))
        allocate(V_n(0:NI,0:NJ))
        allocate(P_n(0:NI,0:NJ))
        allocate(ro_n(0:NI,0:NJ))

        allocate(U_Node(1:NI,1:NJ))
        allocate(V_Node(1:NI,1:NJ))
        allocate(P_Node(1:NI,1:NJ))
        allocate(ro_Node(1:NI,1:NJ))

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

        do J = 1, NJ
            do I=1,NI
                X_Cell(I,J)=(I-0.5)*dx
                Y_Cell(I,J)=(J-0.5)*dy
                X_Node(I,J) = (I-1.0)*dx
                Y_Node(I,J) = (J-1.0)*dy
            end do
        end do

        do J=1,NJ
            X_Cell(0,J) = -dx/2.0
        end do
        X_Cell(0,0) = -dx/2.0
        Y_Cell(0,0) = -dy/2.0
        do I = 1, NI
            Y_Cell(I,0) = -dy/2.0
        end do

        call info('start NavieStocksGas_Solver solver for gas, parameters: eps=' // trim(realToChar(Eps)) // ' NITER='&
        & // trim(intToChar(NITER)) // ' numbers thread=' // trim(intToChar(num)) //' dt=' // trim(realToChar(dt))  )

        !Initial field
        call InitValue(U, V, P, ro, NI, NJ, U0, Ue, ro0, Cc, gamm)
        call BoundValue(U, V, P, ro, NI, NJ, U0, Cc, ro0, gamm, H, Diametr, Ue)
        U_n = U
        V_n = V
        P_n = P

        call info('Write init data in data_ns_init.tec')
        open(IO,FILE='resource/outputres/data_ns_init.tec')
        call OutputFields_Cell(IO,NI,NJ,X_Cell,Y_Cell,U_n,V_n,P_n, ro)
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
            if(N >=3) then
                write(IO_Residuals, *) dt*N, U_Residuals, V_Residuals, P_Residuals
            end if
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
        open(IO, file='resource/outputres/data_ns.dat')
        call OutputFields_Cell(IO,NI,NJ,X_Cell,Y_Cell,U_n,V_n,P_n, ro_n)
        close(IO)
        do I = 1, NI
            do J = 1, NJ
                U_Node(I,J) = (U_n(I,J) + U_n(I-1,J) + U_n(I,J-1) + U_n(I-1,J-1))/4.0
                V_Node(I,J) = (V_n(I,J) + V_n(I-1,J) + V_n(I,J-1) + V_n(I-1,J-1))/4.0
                P_Node(I,J) = (P_n(I,J) + P_n(I-1,J) + P_n(I,J-1) + P_n(I-1,J-1))/4.0
            end do
        end do
        call initializeListener()
        call writeListenXComponent(U_Node, X_Node, Y_Node, NJ, NI, 'U_NavStocks')
        call writeListenXComponent(V_Node, X_Node, Y_Node, NJ, NI, 'V_NavStocks')
        return
    END SUBROUTINE NavieStocksGas_Solver

    SUBROUTINE BoundValue(U,V,P,ro,NI,NJ,U0, C, ro0, gamm, H, Diametr, Ue)
        implicit none

        integer NI,NJ
        real(8) :: U(0:NI,0:NJ), V(0:NI,0:NJ), P(0:NI,0:NJ), ro(0:NI,0:NJ)
        real(8) :: U0, Ue, C, ro0, gamm, H, Diametr

        !левая граница
        U(0,1:NJ) = Ue
        U(0,10:nint(Diametr/H*NJ)+10) = U0
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
