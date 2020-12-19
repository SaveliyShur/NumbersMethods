PROGRAM Main
    use PrandtlSolver_Plate
    use TestAllModules
    use MethodOfEstablishing_Plate
    implicit none
    character(100) :: commandArgument

    if(COMMAND_ARGUMENT_COUNT() > 0) then
        call GET_COMMAND_ARGUMENT(1,commandArgument)
        if(commandArgument .EQ. "test") then
            call test_all()
        else
            write(*,*) "Нераспознанный аргумент"
        end if
    end if

    call MethodOfEstablishinglSolve_Plate()

end PROGRAM Main
