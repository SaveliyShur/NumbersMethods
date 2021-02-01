!������ ��� ������� ������ ��������� ��������
!��� ������� ����������� ��������� ������� ��������� �����-������ ������� ������������
!������� ��������� �������� � resource/Input.txt
!������� ���������� � ���� resource/data_ns.plt
MODULE MethodOfEstablishing_Plate
    implicit none
    CONTAINS

    SUBROUTINE MethodOfEstablishinglSolve_Plate()
        use omp_lib
        implicit none
        integer, parameter:: IO = 1, IO_Residuals = 2, loggers = 3 ! input-output unit
        real, parameter :: Eps = 1e-7
        integer NI, NJ, NITER
        integer I,J, N, num
        real L,H,dx,dy, visk, U0, CFL
        real dt, A, U_Residuals, V_Residuals, P_Residuals
        real,allocatable :: X_Cell(:),Y_Cell(:)
        real,allocatable :: U(:,:),V(:,:),P(:,:)
        real,allocatable :: U_n(:,:),V_n(:,:),P_n(:,:)
        real, allocatable:: U_cap(:,:), V_cap(:,:)
        real, allocatable:: U_i_half(:,:), V_i_half(:,:), P_i_half(:,:)
        real, allocatable:: U_j_half(:,:), V_j_half(:,:), P_j_half(:,:)

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

        allocate(X_Cell(0:NI)) ! Cell Centers
        allocate(Y_Cell(0:NJ)) ! Cell Centers

        !Cell-centered variables
        allocate(U(0:NI,0:NJ))   ! Velocity U
        allocate(V(0:NI,0:NJ))   ! Velocity V
        allocate(P(0:NI,0:NJ))   ! Pressure

        !Node variables
        allocate(U_n(0:NI,0:NJ))   ! Velocity U
        allocate(V_n(0:NI,0:NJ))   ! Velocity V
        allocate(P_n(0:NI,0:NJ))   ! Pressure

        allocate(U_cap(0:NI,0:NJ))
        allocate(V_cap(0:NI,0:NJ))
        allocate(U_i_half(0:NI,0:NJ), V_i_half(0:NI,0:NJ), P_i_half(0:NI,0:NJ))
        allocate(U_j_half(0:NI,0:NJ), V_j_half(0:NI,0:NJ), P_j_half(0:NI,0:NJ))

        dx=L/(NI-1)
        dy=H/(NJ-1)
        A = 1/(U0**2)
        dt = CFL * min(dx/U0, min(0.5*dx*dx/visk, 0.5*dy*dy/visk))
        num = 4
        call omp_set_num_threads(num);
        do I=1,NI
            X_Cell(I)=(I-0.5)*dx
        end do
        do J=1,NJ
            Y_Cell(J)=(J-0.5)*dy
        end do
        X_Cell(0) = -dx/2.0
        Y_Cell(0) = -dy/2.0



    !Initial field
    call InitValue(U, V, P, NI, NJ, U0)
    call BoundValue(U, V, P, NI, NJ, U0)
    U_n = U
    V_n = V
    P_n = P
call writeAnswer(IO,NI,NJ,X_Cell,Y_Cell,U,V,P)

    !Solve equation
    open(IO_Residuals,FILE='source/resource/residuals.dat', status = "replace")
    open(loggers, FILE='source/resource/logs.txt', status = "replace")

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
            write(*,*) "MethodOfEstablishinglSolve_Plate:Complete"
            exit
        endif
        if(MOD(N,100) .eq. 0) then
            write(*,*) "N=", N, "eps=", max(U_Residuals, V_Residuals, P_Residuals)
        end if

        write(IO_Residuals, *) dt*N, U_Residuals, V_Residuals, P_Residuals

        !�������������� ��� ���������� ����
        U = U_n
        V = V_n
        P = P_n

        write(loggers, *) U(0,0), U(1,0), U(0,1), U(1,1)
    end do

    close(IO_Residuals)
    close(loggers)

    call writeAnswer(IO,NI,NJ,X_Cell,Y_Cell,U_n,V_n,P_n)
    return
    END SUBROUTINE MethodOfEstablishinglSolve_Plate

    !������� ����� ����� � ����
    SUBROUTINE writeAnswer(IO,NI,NJ,X,Y,U,V,P)
       implicit none

       integer NI,NJ,IO
       real, dimension(NI):: X
       real, dimension(NJ):: Y
       real, dimension(0:NI,0:NJ)::U,V,P

       write(*,*) 'Output data cell (Navier - Stokes) '
       open(IO,FILE='source/resource/data_ns.tec', status = "replace")
       call OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
       close(IO)
    END SUBROUTINE writeAnswer

    !������� ��� ������ � ������� �������
    SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
        implicit none

        integer NI,NJ,IO, i, j
        real, dimension(NI):: X
        real, dimension(NJ):: Y
        real, dimension(0:NI,0:NJ)::U,V,P

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
        real :: U(0:NI,0:NJ), V(0:NI,0:NJ), P(0:NI,0:NJ)
        real U0

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

        real :: U(0:NI,0:NJ), V(0:NI,0:NJ),P(0:NI,0:NJ)
        real U0
        integer NI,NJ

        U = U0
        V = 0.0
        P = 0.0
        U(:,0) = - U(:,1)

    END SUBROUTINE InitValue

END MODULE MethodOfEstablishing_Plate
