
program test
integer:: e

    e = 67
    call file('test2.dat', e)
  
end program test

subroutine file(nome_file,a)
    character nome_file*(*)
    integer:: a
  
    open(24, file=nome_file)
    write(24,*) 'test', a
    write(24,*) 'dio cane se hfveihfgviu'
    close(24)
end subroutine file