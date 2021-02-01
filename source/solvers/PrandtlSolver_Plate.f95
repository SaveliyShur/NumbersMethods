!Модуль для решения задачи обтекания пластины
!При решении применяется численное решение уравнений Прандтля
!Входные параметры задаются в resource/Input.txt
!Решение выводиться в файл resource/data_pt.plt

MODULE  PrandtlSolver_Plate
    use SweepMethods
    use Logger_module
    implicit none
    CONTAINS

    SUBROUTINE PrandtlSolve_Plate()
        integer, parameter:: IO = 5165146 ! input-output unit
        real :: Eps
        integer NI, NJ
        integer I,J, NITER, ios, s
        real L,H,dx,dy, visk, U0
        real,allocatable :: X_Node(:,:),Y_Node(:,:)
        real,allocatable :: X_Cell(:,:),Y_Cell(:,:)
        real,allocatable :: U_c(:,:),V_c(:,:),P_c(:,:)
        real,allocatable :: U_n(:,:),V_n(:,:),P_n(:,:)
        real, allocatable :: A(:), B(:), C(:), D(:)

        call info('Read projects settings')
        open(IO,file='projectsettings.txt', STATUS='OLD', IOSTAT=ios)
        if(ios/=0) then
            call fatal('projects settings no found')
            stop 1
        end if
        read(IO,*) Eps
        read(IO,*) NITER
        close(IO)
        call info('Read projects settings :: complete')

        call info('Read input file')
        write(*,*) 'Read input file'
        open(IO,FILE='source\resource\inputres\Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) U0
        read(IO,*) visk
        close(IO)
        call info('Read input file : Complete')

        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates
        allocate(X_Cell(0:NI,0:NJ)) ! Cell Centers
        allocate(Y_Cell(0:NI,0:NJ)) ! Cell Centers

        allocate(A(1:NJ))
        allocate(B(1:NJ))
        allocate(C(1:NJ))
        allocate(D(1:NJ))

    !*******************  Cell-centered variables **********
        allocate(U_c(0:NI,0:NJ))   ! Velocity U
        allocate(V_c(0:NI,0:NJ))   ! Velocity V
        allocate(P_c(0:NI,0:NJ))   ! Pressure

    !*******************  Node variables ******************
        allocate(U_n(NI,NJ))   ! Velocity U
        allocate(V_n(NI,NJ))   ! Velocity V
        allocate(P_n(NI,NJ))   ! Pressure

        dx=L/(NI-1)
        dy=H/(NJ-1)

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

    !************************* INITIAL FIELD *********************
        call info('Start Prandtl solver for liquid, parameters: eps=' // realToChar(Eps) // 'NITER='&
    & // intToChar(NITER)   )
        U_c = 0.
        V_c = 0.
        U_n = 0.
        V_n = 0.
        call InitValue(U_c,NI,NJ,U0)
        call BoundValue(U_c,V_c,NI,NJ,U0)
        call InitValue(U_n,NI,NJ,U0)
        call BoundValue(U_n,V_c,NI,NJ,U0)

    !****************** Solve equation ********************
        do I = 2, NI
            U_c(I,1:NJ) = U_c(I-1,1:NJ)
            V_c(I,1:NJ) = V_c(I-1,1:NJ)
            s = 0
            do
                s = s + 1
                A(1) = 0
                B(1) = 1
                C(1) = 0
                D(1) = 0

                A(NJ) = 0
                B(NJ) = 1
                C(NJ) = 0
                D(NJ) = U0

                do J=2,NJ-1
                        A(J) = -V_c(I,J-1)/(2*dy) - visk/(dy**2)
                        B(J) = U_c(I,J)/dx + 2*visk/(dy**2)
                        C(J) = V_c(I,J+1)/(2*dy) - visk/(dy**2)
                        D(J) = U_n(I-1,J)**2 / dx
                enddo

                call progonka(A,B,C,D,NJ,U_n(I,1:NJ))

                V_n(I,1) = 0.
                do J=2,NJ
                        V_n(I,J) = V_n(I,J-1) - (dy/(2.0 * dx)) * (U_n(I,J) - U_n(I-1,J) + U_n(I,J-1) - U_n(I-1,J-1))
                enddo

                If (((maxval(abs(U_n(I,1:NJ)-U_c(I,1:NJ)))/maxval(abs(U_n(I,1:NJ)))).LE.Eps).and.&
                    &((maxval(abs(V_n(I,1:NJ)-V_c(I,1:NJ)))/maxval(abs(V_n(I,1:NJ)))).LE.Eps)) then
                        write(*,*) "s = ", s,  "Prandtl solver for liquid :: Complete"
                        call info('Prandtl solver for liquid :: Complete')
                        exit
                endif

                If (s > NITER) then
                    call error('Prandtl solver for liquid :: Error, errotU=' &
                    & // realToChar(maxval(abs(U_n(I,1:NJ)-U_c(I,1:NJ)))/maxval(abs(U_n(I,1:NJ)))) &
                    & // 'errorV' // realToChar(maxval(abs(V_n(I,1:NJ)-V_c(I,1:NJ)))/maxval(abs(V_n(I,1:NJ)))))
                    write(*,*) 'Stop method for iter'
                    stop 2
                endif

                U_c=U_n
                V_c=V_n

            enddo
        enddo

    !****************** Output Results ********************

        write(*,*) 'Output data node (Prandtl)'
        call info('Output data node (Prandtl)')
        open(IO,FILE='source/resource/outputres/data_pr.tec')
        call OutputFields_Node(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P_n)
        close(IO)
        call info('Output data node (Prandtl) :: Complete')
        return

    END SUBROUTINE

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

    END SUBROUTINE

    SUBROUTINE OutputFields_Node(IO,NI,NJ,X,Y,U,V,P)
        implicit none

        integer NI,NJ,IO
        real, dimension(NI,NJ):: X,Y
        real, dimension(NI,NJ):: U,V,P

        write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
        write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
        write(IO,'(100E25.16)') X(1:NI,1:NJ)
        write(IO,'(100E25.16)') Y(1:NI,1:NJ)
        write(IO,'(100E25.16)') U(1:NI,1:NJ)
        write(IO,'(100E25.16)') V(1:NI,1:NJ)
        write(IO,'(100E25.16)') P(1:NI,1:NJ)

    END  SUBROUTINE

    !Set Boundary Conditio
    SUBROUTINE BoundValue(U,V,NI,NJ,U0)
        implicit none

        integer NI,NJ
        real :: U(1:NI,1:NJ), V(1:NI,1:NJ)
        real U0

        U(1:NI,1) = 0.
        V(1:NI,1) = 0.
        U(1:NI,NJ) = U0

    END SUBROUTINE

    !Set Initial Values
    SUBROUTINE InitValue(U,NI,NJ,U0)
        implicit none

        real :: U(1:NI,1:NJ)
        real U0
        integer NI,NJ

        U(1,1:NJ) = U0

    END SUBROUTINE
END MODULE
