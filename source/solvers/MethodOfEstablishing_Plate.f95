!������ ��� ������� ������ ��������� ��������
!��� ������� ����������� ��������� ������� ��������� �����-������ ������� ������������
!������� ��������� �������� � resource/Input.txt
!������� ���������� � ���� resource/data_ns.plt
MODULE MethodOfEstablishing_Plate
    use Logger_module
    use Listener
    implicit none
    CONTAINS

    SUBROUTINE MethodOfEstablishinglSolve_Plate()
        use omp_lib
        implicit none
        integer, parameter:: IO = 1145661, IO_Residuals = 25451, loggers = 331231, IO_Cf=10121989 ! input-output unit
        real(8) :: Eps
        integer :: NI, NJ, NITER, ios
        integer :: I,J, N, num
        real(8) :: L,H,dx,dy, visk, U0, CFL
        real(8) :: dt, A, U_Residuals, V_Residuals, P_Residuals
        real(8),allocatable :: X_Cell(:, :),Y_Cell(:,:)
        real(8),allocatable :: X_Node(:, :),Y_Node(:,:)
        real(8),allocatable :: U(:,:),V(:,:),P(:,:)
        real(8),allocatable :: U_n(:,:),V_n(:,:),P_n(:,:)
        real(8),allocatable :: U_Node(:,:),V_Node(:,:),P_Node(:,:)
        real(8), allocatable:: U_cap(:,:), V_cap(:,:)
        real(8), allocatable:: U_i_half(:,:), V_i_half(:,:), P_i_half(:,:)
        real(8), allocatable:: U_j_half(:,:), V_j_half(:,:), P_j_half(:,:)

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
        close(IO)
        call info('Read input file :: complete')

        allocate(X_Cell(0:NI,0:NJ))
        allocate(Y_Cell(0:NI,0:NJ))
        allocate(X_Node(1:NI,1:NJ))
        allocate(Y_Node(1:NI,1:NJ))

        allocate(U(0:NI,0:NJ))
        allocate(V(0:NI,0:NJ))
        allocate(P(0:NI,0:NJ))

        allocate(U_n(0:NI,0:NJ))
        allocate(V_n(0:NI,0:NJ))
        allocate(P_n(0:NI,0:NJ))

        allocate(U_Node(1:NI,1:NJ))
        allocate(V_Node(1:NI,1:NJ))
        allocate(P_Node(1:NI,1:NJ))

        allocate(U_cap(0:NI,0:NJ))
        allocate(V_cap(0:NI,0:NJ))
        allocate(U_i_half(0:NI,0:NJ), V_i_half(0:NI,0:NJ), P_i_half(0:NI,0:NJ))
        allocate(U_j_half(0:NI,0:NJ), V_j_half(0:NI,0:NJ), P_j_half(0:NI,0:NJ))

        dx=L/(NI-1)
        dy=H/(NJ-1)
        A = 1/(U0**2)
        dt = CFL * min(dx/U0, min(0.5*dx*dx/visk, 0.5*dy*dy/visk))
        call omp_set_num_threads(num);

        do J = 1, NJ
            do I=1,NI
                X_Cell(I,J)=(I-0.5)*dx
                Y_Cell(I,J)=(J-0.5)*dy
                X_Node(I,J) = (I-1)*dx
                Y_Node(I,J) = (J-1)*dy
            end do
        end do

        do J=1,NJ
            X_Cell(0,J) = -dx/2.0
        end do

        do I = 1, NI
            Y_Cell(I,0) = -dy/2.0
        end do

    call info('start MethodOfEstablishingl solver for liquid, parameters: eps=' // trim(realToChar(Eps)) // ' NITER='&
    & // trim(intToChar(NITER)) // ' numbers thread=' // trim(intToChar(num))  )

    !Initial field
    call InitValue(U, V, P, NI, NJ, U0)
    call BoundValue(U, V, P, NI, NJ, U0)
    U_n = U
    V_n = V
    P_n = P

    !Solve equation
    open(IO_Residuals,FILE='resource/outputres/residuals.dat', status = "replace")
    write(IO_Residuals,*) 'VARIABLES = "N", "U_res", "V_res", "P_res"'
    do N = 1, NITER
        !���������� �������� � ��������� ��������
        !$omp parallel
        !$omp do
        do I = 0, NI-1
            do J = 1, NJ
                U_cap(I,J) = 0.5 * (U(I+1,J) + U(I,J))
                if (U_cap(I,J) .ge. 0.0) then
                    U_i_half(I,J) = U(I,J)
                    V_i_half(I,J) = V(I,J)
                    P_i_half(I,J) = P(I+1,J)
                else
                    U_i_half(I,J) = U(I+1,J)
                    V_i_half(I,J) = V(I+1,J)
                    P_i_half(I,J) = P(I,J)
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
                else
                    U_j_half(I,J) = U(I,J+1)
                    V_j_half(I,J) = V(I,J+1)
                    P_j_half(I,J) = P(I,J)
                end if
            end do
        end do
        !$omp end do
        !$omp barrier
        !$omp do
        !���������� ��������� �������
        do J = 1, NJ-1
            do I = 1, NI-1
                !���������� ��������
                P_n(I,J) = P(I,J) - dt * (U0**2) *( (U_i_half(I,J) - U_i_half(I-1,J))/dx &
                & + (V_j_half(I,J) - V_j_half(I,J-1))/dy )

                !�������� ���������� ���������� ��������
                U_n(I,J) = U(I,J) - dt*( (U_cap(I,J)*U_i_half(I,J) - U_cap(I-1,J)*U_i_half(I-1,J))/dx &
                & + (V_cap(I,J)*U_j_half(I,J) - V_cap(I,J-1)*U_j_half(I,J-1))/dy &
                & + (P_i_half(I,J) - P_i_half(I-1,J))/dx &
                & - visk*(U(I+1,J) - 2*U(I,J) + U(I-1,J))/(dx**2) &
                & - visk*(U(I,J+1) - 2*U(I,J) + U(I,J-1))/(dy**2) )

                !���������� ���������� ���������� ��������
                V_n(I,J) = V(I,J) - dt*( (V_cap(I,J)*V_j_half(I,J) - V_cap(I,J-1)*V_j_half(I,J-1))/dy &
                & + (U_cap(I,J)*V_i_half(I,J) - U_cap(I-1,J)*V_i_half(I-1,J))/dx &
                & + (P_j_half(I,J) - P_j_half(I,J-1))/dy &
                & - visk*(V(I+1,J) - 2*v(I,J) + V(I-1,J))/(dx**2) &
                & - visk*(V(I,J+1) - 2*v(I,J) + V(I,J-1))/(dy**2) )


            end do
        end do
        !$omp end do
        !$omp end parallel

        !�������� � ��������� �������
        call BoundValue(U_n, V_n, P_n, NI, NJ, U0)

        !��������� ���������� � ������� �������
        U_Residuals = maxval(abs(U_n-U))/maxval(abs(U_n))
        V_Residuals = maxval(abs(V_n-V))/maxval(abs(V_n))
        P_Residuals = maxval(abs(P_n-P))/maxval(abs(P_n))
        if( (U_Residuals.le.Eps ) .and. (V_Residuals.le.Eps ) .and. (P_Residuals.le.Eps ) ) then
            write(*,*) "MethodOfEstablishinglSolve_Plate : Complete"
            call info('MethodOfEstablishingl solver for liquid :: Complete, N=' //trim(intToChar(N)))
            exit
        endif
        if(MOD(N,100) .eq. 0) then
            write(*,*) "N=", N, "eps=", max(U_Residuals, V_Residuals, P_Residuals)
        end if

        if(N >= 3) then
            write(IO_Residuals, *) N, U_Residuals, V_Residuals, P_Residuals
        end if

        !�������������� ��� ���������� ����
        U = U_n
        V = V_n
        P = P_n

        if(N .eq. NITER) then
            call error('MethodOfEstablishingl solver for liquid :: Solution underreported, eps=' //&
            & realToChar(max(U_Residuals, V_Residuals, P_Residuals)))
        end if
    end do

    close(IO_Residuals)
    call writeAnswer(IO,NI,NJ,X_Cell,Y_Cell,U_n,V_n,P_n)

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
    END SUBROUTINE MethodOfEstablishinglSolve_Plate

    !������� ����� ����� � ����
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

    !������� ��� ������ � ������� �������
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

    SUBROUTINE BoundValue(U,V,P,NI,NJ,U0)
        implicit none

        integer NI,NJ
        real(8) :: U(0:NI,0:NJ), V(0:NI,0:NJ), P(0:NI,0:NJ)
        real(8) U0

        !����� �������, ����
        U(0,1:NJ) =  U0
        V(0,1:NJ) = 0.0
        P(0,1:NJ) = P(1,1:NJ)

        !������ �������, �����
        U(NI,1:NJ) = U(NI-1,1:NJ)
        V(NI,1:NJ) = V(NI-1,1:NJ)
        P(NI,1:NJ) = 0.0

        !������ ������������� �������
        U(1:NI,0) = - U(1:NI,1)
        V(1:NI,0) = - V(1:NI,1)
        P(1:NI,0) = P(1:NI,1)

        !������� �������, �����
        U(1:NI,NJ) = U(1:NI,NJ-1)
        V(1:NI,NJ) = V(1:NI,NJ-1)
        P(1:NI,NJ) = 0.0
    END SUBROUTINE BoundValue

    SUBROUTINE InitValue(U,V,P,NI,NJ,U0)
        implicit none

        real(8) :: U(0:NI,0:NJ), V(0:NI,0:NJ),P(0:NI,0:NJ)
        real(8) U0
        integer NI,NJ

        U = U0
        V = 0.0
        P = 0.0
        U(:,0) = - U(:,1)

    END SUBROUTINE InitValue

    SUBROUTINE CfSolver(U0, visk, U, X_Node, Y_Node,  NI, NJ, IO_Cf)
        real(8) :: U0, visk, delta_max, delta_medium
        integer :: IO_Cf, NI, NJ, I
        real(8) :: X_Node(1:NI,1:NJ), U(1:NI,1:NJ), Y_Node(1:NI,1:NJ)
        real(8), allocatable :: Cf(:), Cf_Blazius(:),  delta1(:)
        call info('������ � ����� ������������ ������������� � ���� Cf_x.dat')
        allocate (delta1(NI-2))
        allocate(Cf(NI))
        allocate (Cf_Blazius(NI))
        open(IO_Cf, FILE='resource/outputres/Cf_x.dat', status = "replace")
         write(IO_Cf,*) 'VARIABLES = "X", "Cf", "Cf_Blazius"'
        do I = 3, NI
            Cf(I-2) = (2*visk/(U0**2))*(U(I,2)-U(I,1))/(Y_Node(I,2) - Y_Node(I,1))
            Cf_Blazius(I-2) = 0.664/sqrt(U0*X_Node(I,1)/visk)
            if(Cf(I-2) > Cf_Blazius(I-2)) then
                delta1(I-2) = (Cf(I-2) - Cf_Blazius(I-2))/Cf(I-2)
            else
                delta1(I-2) = (Cf_Blazius(I-2)-Cf(I-2))/Cf_Blazius(I-2)
            end if
            write(IO_Cf, *) X_Node(I-2,1), Cf(I-2), Cf_Blazius(I-2)
        end do
        close(IO_Cf)
        delta_max = MAXVAL(delta1(1:NI-2)) * 100
        delta_medium = SUM(delta1(1:NI-2))/(NI-2) *100
        open(IO_Cf, FILE='resource/outputres/Cf_delta.dat', status = "replace")
        write(IO_Cf, *) 'delta_max=', delta_max
        write(IO_Cf, *) 'delta_medium=', delta_medium
        close(IO_Cf)
        deallocate (Cf)
        deallocate (Cf_Blazius)
        deallocate (delta1)
        call info('������ � ����� ������������ ������������� � ���� Cf_x.dat :: Complete')
    END SUBROUTINE CfSolver

END MODULE MethodOfEstablishing_Plate
