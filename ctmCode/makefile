#
# This is a Makefile for the nbyn interpolation program.  This
# program uses three classes MIN, AutosubShift and OniomXs.
# Typing 'make' or 'make count' will create the executable file.
#

# define some Makefile variables for the compiler and compiler flags
# to use Makefile variables later in the Makefile: $()
#
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#
# for C++ define  CC = g++
CC = g++
CFLAGS  = -g -Wall

# typing 'make' will invoke the first target entry in the file 
# (in this case the default target entry)
# you can name this target entry anything, but "default" or "all"
# are the most commonly used names by convention
#
default: count

# To create the executable file count we need the object files
# MIN.o, AutosubShift.o,OniomXs.o and the main program file:
#
count:  MIN.o AutosubShift.o OniomXs.o 
	$(CC) $(CFLAGS) -o count main.cpp  MIN.o AutosubShift.o OniomXs.o

# To create the object file MIN.o, we need the source
# files MIN.cpp, MIN.h:
#
MIN.o:  MIN.cpp MIN.h 
	$(CC) $(CFLAGS) -c MIN.cpp

# To create the object file AutosubShift.o, we need the source files
# AutosubShift.cpp and AutosubShift.h:
#
AutosubShift.o: AutosubShift.cpp AutosubShift.h
	$(CC) $(CFLAGS) -c AutosubShift.cpp

# To create the object file OniomXs.o, we need the source files
# OniomXs.cpp and OniomXs.h:
#
OniomXs.o:  OniomXs.cpp OniomXs.h 
	$(CC) $(CFLAGS) -c OniomXs.cpp

# To start over from scratch, type 'make clean'.  This
# removes the executable file, as well as old .o object
# files and *~ backup files:
#
clean: 
	$(RM) count *.o *~


