# define variables
HDRDIR  = ./

# set options for this machine
# specify which compilers to use for c and linking
CC	= g++
LD	= g++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I$(HDRDIR) -g

# link flags to be used 
LDFLAGS	= -g

# libraries to be linked in
LIBS	=  -llapack -lblas -lgfortran

# types of files we are going to construct rules for
.SUFFIXES: .cpp

# rule for .c files
.cpp.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.cpp

# list of objects to be compiled
OBJS    = main.o

main:$(OBJS) 
	$(LD)  $(LDFLAGS) -o main $(OBJS) $(LIBS)

# make it clear that main.o depends on matrix.hpp
main.o: matrix.hpp

# what to do if user types "make clean"
clean :
	rm -r $(OBJS)
