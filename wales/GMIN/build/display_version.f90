subroutine display_version(file)

integer file

write(file,10) "==========================================="
write(file,10)
write(file,10) " - "
write(file,10)
write(file,10) "Copyright (C) 1999-2010 David J. Wales"
write(file,10)
write(file,10) "SVN revision: " 
write(file,10) 
write(file,10) "Compilation time: Wed 10 Jan 16:17:43 GMT 2018"
write(file,10) "Compiled by bruce@atocha on GNU/Linux x86_64"
write(file,10)
write(file,10) "Compiler name:  "
write(file,10) "Compiler executable:  "
write(file,10) "Compiler flags: "
write(file,10) "Command-line options passed to makefile:  "
write(file,10)
write(file,10) "==========================================="

10 format(a)

end
