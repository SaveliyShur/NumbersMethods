MODULE PrandtlSolverGas
    use Logger_module
    use SweepMethods
    implicit none
    CONTAINS

    SUBROUTINE PrandtlSolve_Gas()
        integer, parameter:: IO = 54156 ! input-output unit
        real(8) :: Eps
        integer NI, NJ
        integer I,J, NITER, ios, s
        real(8) L,H,dx,dy, visk, U0, Diametr, Cc, gamm, ro0, Uoe, Ue
        real(8),allocatable :: X_Node(:,:),Y_Node(:,:)
        real(8),allocatable :: U(:,:),V(:,:),P(:,:)
        real(8),allocatable :: U_n(:,:),V_n(:,:)
        real(8), allocatable :: A(:), B(:), C(:), D(:)
        real(8), allocatable :: ro(:,:)

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
        open(IO,FILE='resource\inputres\Input.txt', IOSTAT=ios)
        if(ios/=0) then
            call fatal('input settings no found')
            stop 1
        end if
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) U0
        read(IO,*) visk
        read(IO,*) Diametr
        read(IO,*) Cc
        read(IO,*) gamm
        read(IO,*) ro0
        read(IO,*) Uoe
        close(IO)
        call info('Read input file : Complete')

        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates

        allocate(A(1:NJ))
        allocate(B(1:NJ))
        allocate(C(1:NJ))
        allocate(D(1:NJ))

        allocate(U(NI,NJ))
        allocate(V(NI,NJ))
        allocate(P(NI,NJ))
        allocate(ro(NI,NJ))

        allocate(U_n(NI,NJ))
        allocate(V_n(NI,NJ))

        call info('Create mesh')
        dx=L/(NI-1)
        dy=H/(NJ-1)
        Ue = U0/Uoe

        do I=1,NI
          do J=1,NJ
            X_Node(I,J)=(I-1)*dx
            Y_Node(I,J)=(J-1)*dy
          end do
        end do
        call info('Create mesh :: Complete')

        call info('Start Prandtl solver for gas, parameters: eps=' // trim(realToChar(Eps)) // ' NITER='&
    & // trim(intToChar(NITER))   )
        call InitValue(U, V, P, ro, NI, NJ, U0, ro0, Cc, gamm, Ue)
        call BoundValue(U, V, NI, NJ, U0, H, Diametr, Ue)
        U_n = U
        V_n = V

        call info('Write init data in data_pr_init.tec')
        open(IO,FILE='resource/outputres/data_pr_init.tec')
        call OutputFields_Node(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P)
        close(IO)
        call info('Write data init :: Complete')

    !****************** Solve equation ********************
        do I = 2, NI
            U(I,1:NJ) = U(I-1,1:NJ)
            V(I,1:NJ) = V(I-1,1:NJ)
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
                D(NJ) = Ue

                do J=2,NJ-1
                        A(J) = -ro(I,J-1)*V(I,J-1)/(2*dy) - visk/(dy**2)
                        B(J) = ro(I,J)*U(I,J)/dx + 2*visk/(dy**2)
                        C(J) = ro(I,J+1)*V(I,J+1)/(2*dy) - visk/(dy**2)
                        D(J) = ro(I-1,J)*U(I-1,J)**2 / dx - (P(I,J) - P(I-1,J))/dx
                enddo

                call progonka(A,B,C,D,NJ,U_n(I,1:NJ))

                V_n(I,1) = 0.
                do J=2,NJ
                        V_n(I,J) = V_n(I,J-1)*ro(I,J-1)/ro(I,J) &
                        & - (dy/(2.0*dx*ro(I,J))) * (U_n(I,J)*ro(I,J) - U_n(I-1,J)*ro(I-1,J)&
                        & + U_n(I,J-1)*ro(I,J-1) - U_n(I-1,J-1)*ro(I-1,J-1))
                enddo

                If (((maxval(abs(U_n(I,1:NJ)-U(I,1:NJ)))/maxval(abs(U_n(I,1:NJ)))).LE.Eps).and.&
                    &((maxval(abs(V_n(I,1:NJ)-V(I,1:NJ)))/maxval(abs(V_n(I,1:NJ)))).LE.Eps)) then
                        write(*,*) "I=", I, " s=", s
                        exit
                endif

                If (s > NITER) then
                    call error('Prandtl solver for gas :: Error, errotU=' &
                    & // trim(realToChar(maxval(abs(U_n(I,1:NJ)-U(I,1:NJ)))/maxval(abs(U_n(I,1:NJ))))) &
                    & // ' errorV=' // trim(realToChar(maxval(abs(V_n(I,1:NJ)-V(I,1:NJ)))/maxval(abs(V_n(I,1:NJ))))) &
                    & // ' I=' // trim(intToChar(I)))
                    exit
                endif

                U = U_n
                V = V_n

            end do
        end do
        write(*,*) "Prandtl solver for gas :: Complete"
        call info('Prandtl solver for gas :: Complete')

    !****************** Output Results ********************

        write(*,*) 'Output data node (Prandtl)'
        call info('Output data node (Prandtl)')
        open(IO,FILE='resource/outputres/data_pr.tec')
        call OutputFields_Node(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P)
        close(IO)
        call info('Output data node (Prandtl) :: Complete')
        return

    END SUBROUTINE

    SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
        implicit none

        integer NI,NJ,IO
        real(8), dimension(NI,NJ):: X,Y
        real(8), dimension(0:NI,0:NJ)::U,V,P

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
        real(8), dimension(NI,NJ):: X,Y
        real(8), dimension(NI,NJ):: U,V,P

        write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
        write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
        write(IO,'(100E25.16)') X(1:NI,1:NJ)
        write(IO,'(100E25.16)') Y(1:NI,1:NJ)
        write(IO,'(100E25.16)') U(1:NI,1:NJ)
        write(IO,'(100E25.16)') V(1:NI,1:NJ)
        write(IO,'(100E25.16)') P(1:NI,1:NJ)

    END  SUBROUTINE

    !Set Boundary Conditio
    SUBROUTINE BoundValue(U,V,NI,NJ,U0, H, Diametr, Ue)
        implicit none

        integer NI,NJ
        real(8) :: U(1:NI,1:NJ), V(1:NI,1:NJ)
        real(8) :: U0, H, Diametr, Ue

        U(1:NI,1) = 0.
        V(1:NI,1) = 0.
        U(1,1:nint(Diametr/H*NJ)) = U0
        U(1,(nint(Diametr/H*NJ)+1):NJ) = Ue
        U(1:NI,NJ) = Ue

    END SUBROUTINE

    !Set Initial Values
    SUBROUTINE InitValue(U,V,P,ro,NI,NJ,U0, ro0, C, gamm, Ue)
        implicit none

        real(8) :: U(1:NI,1:NJ), P(1:NI,1:NJ), ro(1:NI,1:NJ)
        real(8) :: V(1:NI,1:NJ)
        real(8) U0, ro0, C, gamm, Ue
        integer NI,NJ
        ro = ro0
        P = C*(ro0**gamm)
        U = 0.0
        V = 0.0
    END SUBROUTINE


END MODULE
