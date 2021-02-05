module Listener
    use Logger_module
    implicit none
    public :: initializeListener, writeListenXComponent

    private
    logical :: checkListener, isListen
    integer, parameter :: IOListener = 323480, IOListenerWriter = 232102
    real(8), parameter :: eps = 1e-6
    integer :: NLocal
    integer :: N
    real(8), allocatable :: X_Listeners(:)

    CONTAINS

    SUBROUTINE initializeListener()
        call info('Initialize Listener')
        checkListener = .false.
        isListen = .false.
    END SUBROUTINE initializeListener

    SUBROUTINE writeListenXComponent(U, X, Y, NJ, NI, fileName)
        character(100) :: fileNameConcat
        character(*) :: fileName
        integer :: NI, NJ, I, J
        real(8), dimension(1:NI, 1:NJ) :: U
        real(8), dimension(1:NI, 1:NJ) :: X, Y
        if(isListenerReal() .eqv. .false.) then
            return
        end if
        do I = 1, N
            if(giveX(X, NI, NJ, X_Listeners(I)) .eqv. .true.) then
                fileNameConcat = 'resource/outputres/' //trim(fileName) // trim(intToChar(NLocal)) // '.dat'
                call info('Запись локальных значений в файл ' // trim(fileNameConcat) &
                & // ' Для Х = ' // trim(realToChar(X_Listeners(I))))
                open(IOListenerWriter, file=trim(fileNameConcat), STATUS='replace')
                write(IOListenerWriter,*) 'VARIABLES = "Y", "U"'
                do J = 1, NJ
                    write(IOListenerWriter,*) U(NLocal, J), Y(NLocal,J)
                end do
                close(IOListenerWriter)
            end if
        end do
    END SUBROUTINE writeListenXComponent

    logical function isListenerReal()
        integer :: I, ios
        if(checkListener .eqv. .true.) then
            isListenerReal = isListen
            return
        end if
        call info('Read InputListenerReal.txt')
        open(IOListener,file='resource/inputres/InputListenerReal.txt', STATUS='OLD', IOSTAT=ios)
        if(ios/=0) then
            call info('Listener settings no found')
            checkListener = .true.
            isListen = .false.
            isListenerReal = .false.
            return
        end if

        read(IOListener, *) N
        allocate(X_Listeners(N))
        do I = 1, N
            read(IOListener, *) X_Listeners(I)
        end do
        call info('Read InputListenerReal.txt :: Complete')

        checkListener = .true.
        isListen = .true.
        isListenerReal = .true.
        return
    end function isListenerReal

    logical function giveX(X, NI, NJ, XNode)
        integer :: NI, I, NJ
        real(8), dimension(1:NI, 1:NJ) :: X
        real(8) :: XNode
        do I=1, NI
            if(abs(X(I,1) - XNode) < eps) then
                NLocal = I
                giveX = .true.
                return
            end if
        end do
        call error('Ячейка не досупна для вывода локальных значений X=' &
        & // trim(realToChar(XNode)))
        giveX = .false.
        return
    end function giveX
end module Listener
