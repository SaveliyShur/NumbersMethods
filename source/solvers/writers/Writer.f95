module Writer
    use Logger_module
    implicit none
    contains

    !Функция для вывода в формате техплот
    SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P, ro)
        implicit none

        integer NI,NJ,IO, i, j
        real(8), allocatable :: X(:,:)
        real(8), allocatable :: Y(:,:)
        real(8), allocatable::U(:,:),V(:,:),P(:,:), ro(:,:)
        call checkNan(X, NI, NJ)
        call checkNan(Y, NI, NJ)
        call checkNan(U, NI, NJ)
        call info('Write')
        write(*,*) 'Output data cell (Navier - Stokes) '
        write(IO, *) 'VARIABLES = "X", "Y", "U", "V", "P", "RO"'
        write(IO, *) 'ZONE I=', NI, ', J=', NJ, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
        write(IO, '(100E25.16)') X(1:NI, 1:NJ)
        write(IO, '(100E25.16)') Y(1:NI, 1:NJ)
        write(IO, '(100E25.16)') U(1:NI-1, 1:NJ-1)
        write(IO, '(100E25.16)') V(1:NI-1, 1:NJ-1)
        write(IO, '(100E25.16)') P(1:NI-1, 1:NJ-1)
        write(IO, '(100E25.16)') ro(1:NI-1, 1:NJ-1)
        call info('Write :: Complete')
    END SUBROUTINE OutputFields_Cell

    subroutine checkNan(Massiv, NI, NJ)
        real(8), allocatable :: Massiv(:, :)
        integer NI, NJ, I, J

        do I = 1, NI
            do J = 1, Nj
                if(isnan(Massiv(I,J))) then
                    call error('В ответ передано NaN в (' // trim(intToChar(I)) // ',' // trim(intToChar(J)) // ')' )
                end if
            end do
        end do
    end subroutine checkNan
end module Writer
